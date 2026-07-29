/* SPDX-License-Identifier: MIT */

/* src/bam_rec_scratch.h — per-thread scratch for building one BAM record.
 *
 * Shared by the two BAM emitters, src/bam_writer.cpp and src/meth_bam.cpp, so
 * that a buffer optimization or a hardening fix lands on both. They emit the
 * same record shape from the same mem_aln_t, and the meth writer previously
 * carried its own per-record malloc/free of every one of these buffers.
 *
 * Each translation unit holds its own `static thread_local` instance. That is
 * one unused instance per thread in either mode (a run is --meth or it is not,
 * never both), which costs nothing until the buffers are first grown.
 */

#ifndef BWAMEM3_BAM_REC_SCRATCH_H
#define BWAMEM3_BAM_REC_SCRATCH_H

#include "htslib/kstring.h"

#include <stdint.h>

#include <cstdlib>

namespace bwamem3 {

/* Per-thread scratch for the transient bam_cigar / seq_text / qual_bin buffers
 * built when emitting a record. These were three malloc/free per emitted record;
 * bam_set1 copies them into b->data before returning, so a per-thread grow-only
 * buffer reused across records is safe and removes the allocator round-trips
 * (real pressure at tens of millions of records × many threads). Freed on
 * thread exit. */
struct BamRecScratch {
    uint32_t *cigar = nullptr; size_t cigar_cap = 0;   // in uint32 ops
    char     *seq   = nullptr; size_t seq_cap   = 0;   // in bytes (incl NUL)
    char     *qual  = nullptr; size_t qual_cap  = 0;   // in bytes
    /* Grow-only kstrings for the MC:Z (mate CIGAR, ~every paired record) and
     * SA:Z (other primary hits, multi-mapping records) aux tags. Reset .l = 0
     * per record and reused; the buffer persists so building them is no longer
     * a malloc/free per record (freed on thread exit like the buffers above). */
    kstring_t mc = {0, 0, nullptr};
    kstring_t sa = {0, 0, nullptr};
    uint32_t *ensure_cigar(size_t n) { return grow(cigar, cigar_cap, n); }
    char     *ensure_seq(size_t n)   { return grow(seq,   seq_cap,   n); }
    char     *ensure_qual(size_t n)  { return grow(qual,  qual_cap,  n); }
    ~BamRecScratch() { free(cigar); free(seq); free(qual); free(mc.s); free(sa.s); }

    /* Manages raw memory with a custom destructor: delete copies so an accidental
     * copy can't double-free. Only ever used as a static thread_local (never
     * copied), so this is future-proofing, not a live fix. Declaring the copy
     * ctor suppresses the implicit default ctor, so default it explicitly. */
    BamRecScratch() = default;
    BamRecScratch(const BamRecScratch &) = delete;
    BamRecScratch &operator=(const BamRecScratch &) = delete;

private:
    /* Grow `buf` to hold at least `n` elements, returning it (or nullptr).
     *
     * One implementation for all three buffers so a hardening fix cannot land on
     * some of them and not the others. Two properties the previous per-buffer
     * one-liners lacked:
     *
     *   - the size is checked for overflow before `malloc`. `n * sizeof(T)` is
     *     computed in size_t, so on a 32-bit build a large element count could
     *     wrap and under-allocate, after which the caller writes `n` elements
     *     into a shorter buffer. Not reachable on 64-bit from the current callers
     *     (n comes from an int CIGAR-op count or read length), so this is
     *     hardening against a future caller and a 32-bit target, not a live bug.
     *
     *   - the new block is allocated BEFORE the old one is freed. The old form
     *     freed first, so a failed malloc destroyed a buffer that was merely too
     *     small; on failure now, `buf` and `cap` are both left untouched and a
     *     later, smaller record can still use the existing buffer.
     *
     * Callers all check the returned pointer and fail the record, so neither of
     * these was a live NULL-dereference. */
    template <typename T>
    static T *grow(T *&buf, size_t &cap, size_t n) {
        if (n <= cap) return buf;
        if (n > SIZE_MAX / sizeof(T)) return nullptr;
        T *next = static_cast<T *>(malloc(n * sizeof(T)));
        if (next == nullptr) return nullptr;
        free(buf);
        buf = next;
        cap = n;
        return buf;
    }
};

}  // namespace bwamem3

#endif /* BWAMEM3_BAM_REC_SCRATCH_H */
