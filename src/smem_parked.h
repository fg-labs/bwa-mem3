#ifndef SMEM_PARKED_H
#define SMEM_PARKED_H

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "fmi_seed_api.h"    // SMEM
#include "grow_capacity.h"
#include "lockstep_width.h"  // SMEM_LOCKSTEP_N_MAX
#include "utils.h"           // xassert

/* In-order emission for the lockstep SMEM driver when reads finish out of order
 * (getSMEMsOnePosOneThread_lockstep in FMI_search.cpp).
 *
 * The driver refills a slot as soon as its read finishes, so reads can finish
 * ahead of an earlier, slower read. finish() emits a read straight to the output
 * when it is the next one due (`cursor`); otherwise it copies the read's matches
 * into `smems` and records them in `ring` until every earlier read is out. Reads
 * leave strictly in input order and each read's matches keep their order, so
 * the output is exactly what an in-order driver writes.
 *
 * The driver keeps every in-flight or parked read in [cursor, cursor + window)
 * (admits()), so input i has its own ring entry at i & (window - 1). An entry's
 * len is -1 when empty; 0 is a parked read with no matches. The window is what
 * lets the ring be a fixed size, and it caps how far one slow read lets the
 * others run ahead of it.
 *
 * The window bounds how many reads are parked, not how many SMEMs: `smems` is
 * appended to and rewound only when no read is parked, so within a call it can
 * reach about that call's whole output. Per thread and reused across calls:
 * `smems` only grows, like the driver's prev[]/match_buf[] cache. */
/* The smallest power of two >= v (v >= 1). Outside the struct so the ring size
 * below can use it in a constant expression. */
static inline constexpr int32_t smem_parked_pow2_at_least(int32_t v)
{
    return v <= 1 ? 1 : 2 * smem_parked_pow2_at_least((v + 1) / 2);
}

struct LockstepParked {
    /* The driver's window: 8 reads per slot, rounded up to a power of two. */
    static constexpr int32_t kReadsPerSlot = 8;
    static constexpr int32_t window_for(int32_t n_slots) { return smem_parked_pow2_at_least(kReadsPerSlot * n_slots); }

    /* Start a call: nothing emitted or parked, inputs from 0. */
    void reset(int32_t window_) {
        xassert(window_ >= 1 && window_ <= kRingMax && (window_ & (window_ - 1)) == 0,
                "LockstepParked: window must be a power of two <= the ring size");
        window  = window_;
        cursor  = 0;
        pending = 0;
        n       = 0;
        for (int32_t i = 0; i < window; i++) len[i] = -1;
    }

    /* Whether input `next` may start: it must fit in the window ahead of `cursor`. */
    bool admits(int32_t next) const { return next - cursor < window; }

    /* Input `idx` finished with `count` matches in `buf`. Writes every read that is
     * now due to out[*n_out...], in input order. */
    void finish(int32_t idx, const SMEM *buf, int32_t count, SMEM *out, int64_t *n_out) {
        if (idx != cursor) {
            const int32_t e = idx & (window - 1);
            xassert(idx > cursor && idx - cursor < window && len[e] < 0,
                    "LockstepParked: finished read outside the window or already parked");
            if (count > 0) {
                if (n + (size_t)count > cap) grow(n + (size_t)count);
                memcpy(smems + n, buf, (size_t)count * sizeof(SMEM));
            }
            start[e] = n;
            len[e]   = count;
            n += (size_t)count;
            pending++;
            return;
        }
        if (count > 0) memcpy(out + *n_out, buf, (size_t)count * sizeof(SMEM));
        *n_out += count;
        cursor++;
        while (pending > 0) {
            const int32_t e = cursor & (window - 1);
            if (len[e] < 0) break;
            if (len[e] > 0) memcpy(out + *n_out, smems + start[e], (size_t)len[e] * sizeof(SMEM));
            *n_out += len[e];
            len[e] = -1;
            pending--;
            cursor++;
        }
        if (pending == 0) n = 0;
    }

    ~LockstepParked() { free(smems); }

    int32_t cursor  = 0;   // next input to emit
    int32_t window  = 1;
    int32_t pending = 0;   // parked reads

private:
    static constexpr int32_t kRingMax = smem_parked_pow2_at_least(kReadsPerSlot * SMEM_LOCKSTEP_N_MAX);

    void grow(size_t need) {
        xassert(need <= SIZE_MAX / sizeof(SMEM), "LockstepParked: parked SMEM count overflows");
        const size_t want = grow_capacity(cap, need, sizeof(SMEM));
        SMEM *p = (SMEM *)realloc(smems, want * sizeof(SMEM));
        xassert(p != NULL, "LockstepParked: out of memory");
        smems = p;
        cap   = want;
    }

    SMEM   *smems = nullptr;   // parked matches, back to back
    size_t  n     = 0;         // SMEMs in use
    size_t  cap   = 0;
    size_t  start[kRingMax];
    int32_t len[kRingMax];
};

#endif /* SMEM_PARKED_H */
