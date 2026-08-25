#include "lockstep_width.h"

#include <stdlib.h>
#include <stdint.h>
#include <stddef.h>
#include <stdio.h>
#include <time.h>
#include <errno.h>
#include <pthread.h>

/* Runtime phase-2 lockstep width. Defaults to the compile-time floor so a
 * binary that never calls the startup probe is identical to the old constant. */
int32_t g_smem_lockstep_n = SMEM_LOCKSTEP_N;

int32_t bwa3_lockstep_width_from_probe(int32_t raw_mlp) {
    if (raw_mlp < SMEM_LOCKSTEP_N)     return SMEM_LOCKSTEP_N;      /* floor (also covers <=0) */
    if (raw_mlp > SMEM_LOCKSTEP_N_MAX) return SMEM_LOCKSTEP_N_MAX;  /* ceiling */
    return raw_mlp;
}

int32_t bwa3_lockstep_width_parse_env(const char *env) {
    if (env == NULL || env[0] == '\0') return 0;  /* unset: run the probe */

    errno = 0;
    char *end = NULL;
    long v = strtol(env, &end, 10);
    if (end == env || *end != '\0') return -1;  /* not a clean integer */
    if (errno == ERANGE)            return -1;  /* overflowed the parse */
    if (v < 1)                      return -1;  /* non-positive */
    if (v > SMEM_LOCKSTEP_N_MAX)    return SMEM_LOCKSTEP_N_MAX;  /* clamp large to ceiling */
    return (int32_t)v;
}

/* ---- startup memory-level-parallelism probe --------------------------------
 *
 * A pointer chase whose next address is derived from the value just loaded is
 * a pure chain of dependent loads: the prefetcher cannot run ahead and every
 * cache-missing hop is a serialized memory stall. Running k such chains
 * interleaved lets the core keep up to k misses outstanding at once; per-access
 * latency falls as k rises while the core still has miss slots to spare and
 * flattens once it does not. The smallest k at (near) that floor is the core's
 * usable memory-level parallelism -- exactly how many reads' FM-index walks the
 * lockstep seeding driver can profitably keep in flight.
 *
 * The chase walks the FM-index's own cp_occ checkpoint array, NOT a scratch
 * buffer. A fresh allocation measures whatever cache level it happens to fit
 * in -- on a large system-level cache it reports near-cache latency and picks
 * the widest candidate for the wrong reason. cp_occ is gigabytes, already
 * resident, carries the real TLB behaviour, and is the array seeding actually
 * hammers, so probing it measures the thing that matters and costs no memory.
 */

static int64_t bwa3_now_ns(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (int64_t)ts.tv_sec * 1000000000LL + ts.tv_nsec;
}

/* Candidates span the plausible MLP range; sub-floor values still get
 * measured (a low knee just floors to the default in from_probe). Shared by
 * bwa3_measure_mlp and the per-candidate worker below. */
static const int32_t bwa3_mlp_sweep[] = {8, 12, 16, 24, 32, 48, 64};
static const size_t  bwa3_mlp_n_sweep = sizeof(bwa3_mlp_sweep) / sizeof(bwa3_mlp_sweep[0]);

/* Per-candidate probe: everything bwa3_measure_mlp's loop body used to do for
 * one sweep entry, unchanged (same best-of-2 timing, same warm-up, same hop
 * pattern) -- only extracted so it can run on its own thread. `cur`/`sink`
 * that were function-local before are now per-call stack locals, so each
 * concurrently-running candidate has its own, exactly as the sequential
 * version implicitly did between loop iterations. */
typedef struct {
    const char *bp;
    int64_t     n_blocks;
    size_t      stride, word_off;
    int32_t     k;
    double      ns;  /* out: ns/access, lower better */
} bwa3_mlp_candidate_t;

static void *bwa3_mlp_probe_candidate(void *arg_) {
    bwa3_mlp_candidate_t *arg = (bwa3_mlp_candidate_t *)arg_;
    const char    *bp = arg->bp;
    const int64_t  n_blocks = arg->n_blocks;
    const size_t   stride = arg->stride, word_off = arg->word_off;
    const int32_t  k = arg->k;
    const uint64_t MIX = 0x9E3779B97F4A7C15ULL;  /* golden-ratio odd multiplier */
    const int64_t  accesses = 200000;            /* equal work per candidate */
    const int64_t  steps = accesses / k;

    uint64_t cur[SMEM_LOCKSTEP_N_MAX];
    volatile uint64_t sink = 0;

    /* One hop: read a 64-bit word from block cur, mix with the index so the
     * value alone chooses the next block -- a genuine data dependency. */
    #define BWA3_HOP(idx) do { \
        uint64_t v = *(const uint64_t *)(bp + ((idx) % (uint64_t)n_blocks) * stride + word_off) ^ (idx); \
        (idx) = (v * MIX >> 24) % (uint64_t)n_blocks; \
    } while (0)

    /* Best-of-2: the faster timing rejects a scheduling hiccup or migration
     * that would otherwise inflate a single point. */
    int64_t best_dt = 0;
    for (int rep = 0; rep < 2; rep++) {
        for (int32_t j = 0; j < k; j++)
            cur[j] = (uint64_t)((int64_t)j * (n_blocks / k));
        /* Untimed warm-up so every width starts from the same TLB/cache
         * state and does not eat first-touch cost. */
        for (int w = 0; w < 64; w++)
            for (int32_t j = 0; j < k; j++) BWA3_HOP(cur[j]);
        const int64_t t0 = bwa3_now_ns();
        for (int64_t s = 0; s < steps; s++)
            for (int32_t j = 0; j < k; j++) BWA3_HOP(cur[j]);
        const int64_t t1 = bwa3_now_ns();
        for (int32_t j = 0; j < k; j++) sink += cur[j];  /* defeat DCE */
        const int64_t dt = (t1 - t0) > 0 ? (t1 - t0) : 1;
        if (best_dt == 0 || dt < best_dt) best_dt = dt;
    }
    (void)sink;
    #undef BWA3_HOP

    arg->ns = (double)best_dt / (double)(steps * k);
    return NULL;
}

int32_t bwa3_measure_mlp(const void *base, int64_t n_blocks,
                         size_t stride, size_t word_off) {
    /* Too small to measure DRAM behaviour (e.g. a tiny test index): let the
     * caller fall back to the compiled default. */
    if (base == NULL || n_blocks < 4096) return 0;

    const int32_t max_c = SMEM_LOCKSTEP_N_MAX;

    /* n_run: how many of the (ascending, so this is a simple prefix count)
     * sweep candidates fit under the compiled ceiling -- identical to the
     * old loop's `if (k > max_c) break`. Always all 7 at the shipped
     * SMEM_LOCKSTEP_N_MAX=64, but a smaller -D override can still cap it. */
    size_t n_run = 0;
    while (n_run < bwa3_mlp_n_sweep && bwa3_mlp_sweep[n_run] <= max_c) n_run++;
    if (n_run == 0) return 0;

    /* Run every candidate CONCURRENTLY, attempting one thread each, instead
     * of one after another: the seeding worker threads have not spawned yet
     * at this point in startup, so -t's worth of cores are otherwise idle for
     * the whole probe. Any candidate whose pthread_create fails is measured
     * inline in this calling thread instead, so a thread-count cap degrades to
     * a sequential sweep rather than losing candidates.
     * Same measurements, same knee-selection logic below
     * (byte-identical to the sequential version modulo ordinary
     * measurement noise) -- this only changes WHEN each candidate's work
     * runs, not what is measured or how the result is chosen. On an otherwise
     * idle host wall time trends from the SUM of all candidates toward the
     * SLOWEST one, but the size of that gain is hardware-dependent: the
     * candidates now contend for shared CPU, LLC, and memory bandwidth (see the
     * trade-off below), so it need not approach a single candidate or scale with
     * the candidate count -- on a quiet 16-core host the measured startup saving
     * was only ~3%.
     *
     * Trade-off, noted rather than hidden: candidates now share memory
     * bandwidth and LLC while probing (they didn't before), which could
     * shift the absolute ns/access values. Since actual seeding likewise
     * runs all -t threads concurrently and shares that same bandwidth, this
     * arguably measures the more representative condition, not a less
     * accurate one -- but it hasn't been validated against the sequential
     * baseline on real hardware, so treat that as an open question, not a
     * settled one. */
    bwa3_mlp_candidate_t cand[bwa3_mlp_n_sweep];
    pthread_t             tid[bwa3_mlp_n_sweep];
    /* pthread_t is opaque, so no value (0 included) is a portable
     * invalid-thread sentinel: a successful create can legitimately yield a
     * handle that compares equal to 0. Track creation explicitly so every
     * thread that was actually started is joined -- skipping a join would let
     * a still-running worker race on cand[si].ns and touch this frame's stack
     * after bwa3_measure_mlp has returned. */
    int                   created[bwa3_mlp_n_sweep];
    for (size_t si = 0; si < n_run; si++) {
        cand[si].bp = (const char *)base;
        cand[si].n_blocks = n_blocks;
        cand[si].stride = stride;
        cand[si].word_off = word_off;
        cand[si].k = bwa3_mlp_sweep[si];
        cand[si].ns = 0.0;
        created[si] = (pthread_create(&tid[si], NULL, bwa3_mlp_probe_candidate, &cand[si]) == 0);
        if (!created[si]) {
            /* Thread creation failure (e.g. a process-level thread-count
             * cap) -- fall back to running this one candidate inline rather
             * than losing the measurement entirely. */
            bwa3_mlp_probe_candidate(&cand[si]);
        }
    }
    for (size_t si = 0; si < n_run; si++)
        if (created[si]) pthread_join(tid[si], NULL);

    double ns[bwa3_mlp_n_sweep];
    for (size_t si = 0; si < n_run; si++) ns[si] = cand[si].ns;

    /* Best (lowest) per-access latency, then the knee = smallest candidate
     * within 5% of it. A flat top thus picks the cheaper (fewer-slot) width. */
    double best_ns = ns[0];
    for (size_t si = 1; si < n_run; si++)
        if (ns[si] < best_ns) best_ns = ns[si];
    int32_t knee = bwa3_mlp_sweep[n_run - 1];
    for (size_t si = 0; si < n_run; si++)
        if (ns[si] <= best_ns * 1.05) { knee = bwa3_mlp_sweep[si]; break; }

    /* BWA3_MLP_DEBUG: dump the sweep so the knee choice is auditable. */
    if (getenv("BWA3_MLP_DEBUG") != NULL) {
        fprintf(stderr, "[mlp] ns/access chasing cp_occ:");
        for (size_t si = 0; si < n_run; si++)
            fprintf(stderr, " %d=%.2f", bwa3_mlp_sweep[si], ns[si]);
        fprintf(stderr, "  knee=%d\n", knee);
    }
    return knee;
}

void bwa3_init_smem_lockstep_width(const void *base, int64_t n_blocks,
                                   size_t stride, size_t word_off) {
    static int done = 0;
    if (done) return;

    const char *env = getenv("BWA3_SMEM_LOCKSTEP_N");
    const int32_t pinned = bwa3_lockstep_width_parse_env(env);
    if (pinned > 0) {
        /* A valid override pins the width and disables the probe, so gated/CI
         * runs pay no measurement cost and stay deterministic. */
        g_smem_lockstep_n = pinned;
    } else {
        if (pinned < 0)
            /* Set but invalid: do NOT silently fall through to the floor -- that
             * would disable auto-tuning without a trace. Diagnose and probe. */
            fprintf(stderr,
                    "ERROR: BWA3_SMEM_LOCKSTEP_N=\"%s\" is not a positive integer "
                    "(<= %d); ignoring it and auto-tuning the lockstep width.\n",
                    env, SMEM_LOCKSTEP_N_MAX);
        g_smem_lockstep_n = bwa3_lockstep_width_from_probe(
            bwa3_measure_mlp(base, n_blocks, stride, word_off));
    }
    done = 1;
}
