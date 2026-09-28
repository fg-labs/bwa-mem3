#include "lockstep_width.h"

#include <stdlib.h>
#include <stdint.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <errno.h>
#include <mutex>

/* Runtime phase-2 lockstep width. Defaults to the compile-time floor so a
 * binary that never calls the startup probe is identical to the old constant. */
int32_t g_smem_lockstep_n = SMEM_LOCKSTEP_N;

int32_t bwa3_lockstep_width_from_probe(int32_t raw_mlp) {
    if (raw_mlp < SMEM_LOCKSTEP_N)     return SMEM_LOCKSTEP_N;      /* floor (also covers <=0) */
    if (raw_mlp > SMEM_LOCKSTEP_N_MAX) return SMEM_LOCKSTEP_N_MAX;  /* ceiling */
    return raw_mlp;
}

int32_t bwa3_lockstep_width_parse_env(const char *env) {
    if (env == NULL || env[0] == '\0') return 0;  /* unset/empty: no override */

    errno = 0;
    char *end = NULL;
    long v = strtol(env, &end, 10);
    if (end == env || *end != '\0') return -1;  /* not a clean integer */
    if (errno == ERANGE)            return -1;  /* overflowed the parse */
    if (v < 1)                      return -1;  /* non-positive */
    if (v > SMEM_LOCKSTEP_N_MAX)    return SMEM_LOCKSTEP_N_MAX;  /* clamp large to ceiling */
    return (int32_t)v;
}

int bwa3_lockstep_probe_enabled(const char *env) {
    /* Truthy opt-in: enable only for a present, non-empty value that is not "0".
     * Unset (NULL), empty, and "0" all leave the probe off, so neither an empty
     * value nor an explicit =0 can silently pay the startup sweep's cost. */
    return env != NULL && env[0] != '\0' && strcmp(env, "0") != 0;
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

int32_t bwa3_measure_mlp(const void *base, int64_t n_blocks,
                         size_t stride, size_t word_off) {
    /* Too small to measure DRAM behaviour (e.g. a tiny test index): let the
     * caller fall back to the compiled default. */
    if (base == NULL || n_blocks < 4096) return 0;

    const char    *bp  = (const char *)base;
    const uint64_t MIX = 0x9E3779B97F4A7C15ULL;  /* golden-ratio odd multiplier */
    const int32_t  max_c = SMEM_LOCKSTEP_N_MAX;
    const int64_t  accesses = 200000;            /* equal work per candidate */

    /* Candidates span the plausible MLP range; sub-floor values still get
     * measured (a low knee just floors to the default in from_probe). */
    static const int32_t sweep[] = {8, 12, 16, 24, 32, 48, 64};
    const size_t n_sweep = sizeof(sweep) / sizeof(sweep[0]);
    double  ns[sizeof(sweep) / sizeof(sweep[0])] = {0.0};  /* ns/access; lower better */

    uint64_t cur[SMEM_LOCKSTEP_N_MAX];
    volatile uint64_t sink = 0;

    /* One hop: read a 64-bit word from block cur, mix with the index so the
     * value alone chooses the next block -- a genuine data dependency. */
    #define BWA3_HOP(idx) do { \
        uint64_t v = *(const uint64_t *)(bp + ((idx) % (uint64_t)n_blocks) * stride + word_off) ^ (idx); \
        (idx) = (v * MIX >> 24) % (uint64_t)n_blocks; \
    } while (0)

    size_t n_run = 0;
    for (size_t si = 0; si < n_sweep; si++) {
        int32_t k = sweep[si];
        if (k > max_c) break;
        const int64_t steps = accesses / k;
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
        ns[si] = (double)best_dt / (double)(steps * k);
        n_run = si + 1;
    }
    (void)sink;
    #undef BWA3_HOP

    if (n_run == 0) return 0;

    /* Best (lowest) per-access latency, then the knee = smallest candidate
     * within 5% of it. A flat top thus picks the cheaper (fewer-slot) width. */
    double best_ns = ns[0];
    for (size_t si = 1; si < n_run; si++)
        if (ns[si] < best_ns) best_ns = ns[si];
    int32_t knee = sweep[n_run - 1];
    for (size_t si = 0; si < n_run; si++)
        if (ns[si] <= best_ns * 1.05) { knee = sweep[si]; break; }

    /* BWA3_MLP_DEBUG: dump the sweep so the knee choice is auditable. */
    if (getenv("BWA3_MLP_DEBUG") != NULL) {
        fprintf(stderr, "[mlp] ns/access chasing cp_occ:");
        for (size_t si = 0; si < n_run; si++)
            fprintf(stderr, " %d=%.2f", sweep[si], ns[si]);
        fprintf(stderr, "  knee=%d\n", knee);
    }
    return knee;
}

void bwa3_init_smem_lockstep_width(const void *base, int64_t n_blocks,
                                   size_t stride, size_t word_off) {
    /* One-time resolution. The sole intended caller (fastmap.cpp) runs this
     * once before the seeding workers spawn, so today it is single-threaded;
     * std::call_once matches the repo's other once-init sites
     * (simd_dispatch.cpp, read_memo.cpp) and keeps it correct if that ever
     * changes. */
    static std::once_flag once;
    std::call_once(once, [&]() {
    const char *env = getenv("BWA3_SMEM_LOCKSTEP_N");
    const int32_t pinned = bwa3_lockstep_width_parse_env(env);
    if (pinned > 0) {
        /* A valid override pins the width explicitly, so gated/CI runs pay no
         * measurement cost and stay deterministic. */
        g_smem_lockstep_n = pinned;
    } else {
        if (pinned < 0)
            /* Set but invalid: do NOT silently accept it. Diagnose, then fall
             * through to the same resolution as an unset value below. */
            fprintf(stderr,
                    "ERROR: BWA3_SMEM_LOCKSTEP_N=\"%s\" is not a positive integer "
                    "(<= %d); ignoring it (resolving as if unset).\n",
                    env, SMEM_LOCKSTEP_N_MAX);
        /* No explicit pin. The startup MLP probe is OPT-IN, not the default:
         * measured across architectures (x86/arm), core counts, and page sizes,
         * the lockstep width is a flat ~1% knob end-to-end, while the probe
         * itself costs a non-trivial single-threaded startup pass. So by default
         * we keep the compile-time SMEM_LOCKSTEP_N (the value the probe converges
         * to on all benchmarked hardware -- output is unchanged, the probe cost
         * is not paid). BWA3_SMEM_LOCKSTEP_PROBE=1 opts back into the probe to
         * self-calibrate on new/untested hardware. Truthy opt-in
         * (bwa3_lockstep_probe_enabled): any value other than unset, empty, or
         * "0" enables it, so unset it (or set it to 0) to disable. */
        if (bwa3_lockstep_probe_enabled(getenv("BWA3_SMEM_LOCKSTEP_PROBE")))
            g_smem_lockstep_n = bwa3_lockstep_width_from_probe(
                bwa3_measure_mlp(base, n_blocks, stride, word_off));
        /* else: g_smem_lockstep_n keeps its compile-time SMEM_LOCKSTEP_N init. */
    }
    });
}

/* ---- third-pass bwtseed lockstep: on/off policy (see lockstep_width.h) ----- */

/* The compile-time default for the third-pass bwtseed lockstep, and the single
 * place it is written down: on, on every platform and at every thread count
 * (measurements, and the range they cover, in lockstep_width.h). Used both to
 * initialize g_bwtseed_lockstep and as the fallback in bwa3_init_bwtseed_lockstep
 * -- so a later unset run always falls back here, never to a leaked prior pin. */
#define BWA3_BWTSEED_LOCKSTEP_DEFAULT 1

int32_t g_bwtseed_lockstep = BWA3_BWTSEED_LOCKSTEP_DEFAULT;

/* ---- third-pass bwtseed lockstep WIDTH (see lockstep_width.h) -------------- */

/* Runtime third-pass lockstep width. Defaults to the compile-time BWTSEED_LOCKSTEP_N
 * so a binary that never calls the resolver is identical to the constant. */
int32_t g_bwtseed_lockstep_n = BWTSEED_LOCKSTEP_N;

int32_t bwa3_bwtseed_width_parse_env(const char *env) {
    if (env == NULL || env[0] == '\0') return 0;  /* unset/empty: keep the default */

    errno = 0;
    char *end = NULL;
    long v = strtol(env, &end, 10);
    if (end == env || *end != '\0') return -1;              /* not a clean integer */
    if (errno == ERANGE)            return -1;              /* overflowed the parse */
    if (v < 1)                      return -1;              /* non-positive */
    if (v > BWTSEED_LOCKSTEP_N_MAX) return BWTSEED_LOCKSTEP_N_MAX;  /* clamp to ceiling */
    return (int32_t)v;
}

void bwa3_init_bwtseed_lockstep_width(void) {
    /* One-time resolution; the sole caller (fastmap.cpp) runs it once before the
     * seeding workers spawn. std::call_once matches the other once-init sites. */
    static std::once_flag once;
    std::call_once(once, [] {
        const char *env = getenv("BWA3_BWTSEED_LOCKSTEP_N");
        const int32_t pinned = bwa3_bwtseed_width_parse_env(env);
        if (pinned > 0) {
            g_bwtseed_lockstep_n = pinned;   /* explicit pin (also the gate/CI path) */
        } else if (pinned < 0) {
            fprintf(stderr,
                    "ERROR: BWA3_BWTSEED_LOCKSTEP_N=\"%s\" is not a positive integer "
                    "(<= %d); ignoring it (resolving as if unset).\n",
                    env, BWTSEED_LOCKSTEP_N_MAX);
        }
        /* else: g_bwtseed_lockstep_n keeps its compile-time default. */
    });
}

int32_t bwa3_bwtseed_lockstep_parse_env(const char *env) {
    if (env == NULL || env[0] == '\0') return -1;  /* unset/empty: the compiled default */
    if (strcmp(env, "0") == 0) return 0;
    if (strcmp(env, "1") == 0) return 1;
    return -2;                                     /* invalid: report, then the default */
}

void bwa3_init_bwtseed_lockstep(void) {
    /* Resolve on every run, not once process-wide, so a library caller that
     * invokes main_mem more than once gets each run's own pin rather than whichever
     * came first. The sole caller runs this on the main thread before the seeding
     * workers spawn, so the plain write to g_bwtseed_lockstep needs no
     * synchronization. A 0/1 pin is taken alone; unset and invalid both fall back
     * to the compile-time default -- never the mutable g_bwtseed_lockstep -- so an
     * earlier explicit pin does not leak into a later run whose environment is
     * unset. */
    const char *env = getenv("BWA3_BWTSEED_LOCKSTEP");
    const int32_t pinned = bwa3_bwtseed_lockstep_parse_env(env);
    if (pinned == -2)
        fprintf(stderr,
                "ERROR: BWA3_BWTSEED_LOCKSTEP=\"%s\" is not 0 or 1; ignoring it "
                "(resolving as if unset).\n", env);
    g_bwtseed_lockstep = (pinned == 0 || pinned == 1) ? pinned : BWA3_BWTSEED_LOCKSTEP_DEFAULT;
}
