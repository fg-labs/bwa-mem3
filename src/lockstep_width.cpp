#include "lockstep_width.h"

#include <stdlib.h>
#include <stdint.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <errno.h>
#include <mutex>
#include "lockstep_width_internal.h"
#include "simd_dispatch.h"  /* bwa3_host_cpu_identity */

/* ---- host-tuned default widths (policy and measurements: lockstep_width.h) - */

static int s_smem_lockstep_n_source    = BWA3_LOCKSTEP_WIDTH_SRC_DEFAULT;
static int s_bwtseed_lockstep_n_source = BWA3_LOCKSTEP_WIDTH_SRC_DEFAULT;

int bwa3_smem_lockstep_n_source(void)    { return s_smem_lockstep_n_source; }
int bwa3_bwtseed_lockstep_n_source(void) { return s_bwtseed_lockstep_n_source; }

int bwa3_amd_family_1ah_model_is_zen5(uint32_t model) {
    /* Linux bsp_init_amd, family 0x1a: Zen 5 = 0x00-0x2f, 0x40-0x4f, 0x60-0x7f,
     * 0xd0-0xd7 (Zen 6 = 0x50-0x5f, 0x80-0xaf, 0xc0-0xcf, 0xd8-0xef). */
    return model <= 0x2fu
        || (model >= 0x40u && model <= 0x4fu)
        || (model >= 0x60u && model <= 0x7fu)
        || (model >= 0xd0u && model <= 0xd7u);
}

int bwa3_cpu_is_amd_zen5(int is_amd, uint32_t family, uint32_t model) {
    return is_amd && family == 0x1au && bwa3_amd_family_1ah_model_is_zen5(model);
}

/* Clamp a width into [lo, hi]. */
static int32_t bwa3_clamp_width(int32_t w, int32_t lo, int32_t hi) {
    return w < lo ? lo : (w > hi ? hi : w);
}

int32_t bwa3_smem_host_default_width(int is_amd, uint32_t family, uint32_t model) {
    if (bwa3_cpu_is_amd_zen5(is_amd, family, model))
        return bwa3_clamp_width(BWA3_ZEN5_SMEM_LOCKSTEP_N,
                                SMEM_LOCKSTEP_N, SMEM_LOCKSTEP_N_MAX);
    return SMEM_LOCKSTEP_N;
}

int32_t bwa3_bwtseed_host_default_width(int is_amd, uint32_t family, uint32_t model) {
    if (bwa3_cpu_is_amd_zen5(is_amd, family, model))
        return bwa3_clamp_width(BWA3_ZEN5_BWTSEED_LOCKSTEP_N,
                                BWTSEED_LOCKSTEP_N, BWTSEED_LOCKSTEP_N_MAX);
    return BWTSEED_LOCKSTEP_N;
}

/* The host default for a width whose compile-time default is `dflt`: the
 * clamped host_width, labelled HOST only when it differs from `dflt`. */
static bwa3_lockstep_width_res_t bwa3_host_res(int32_t host_width, int32_t dflt,
                                               int32_t cap) {
    bwa3_lockstep_width_res_t r;
    r.width = bwa3_clamp_width(host_width, dflt, cap);
    r.src   = r.width != dflt ? BWA3_LOCKSTEP_WIDTH_SRC_HOST
                              : BWA3_LOCKSTEP_WIDTH_SRC_DEFAULT;
    return r;
}

bwa3_lockstep_width_res_t bwa3_resolve_smem_width(int32_t pinned, int probe_on,
                                                  int32_t host_width) {
    if (pinned > 0) {
        bwa3_lockstep_width_res_t r = {pinned, BWA3_LOCKSTEP_WIDTH_SRC_ENV};
        return r;
    }
    if (probe_on) {
        bwa3_lockstep_width_res_t r = {0, BWA3_LOCKSTEP_WIDTH_SRC_PROBE};
        return r;
    }
    return bwa3_host_res(host_width, SMEM_LOCKSTEP_N, SMEM_LOCKSTEP_N_MAX);
}

bwa3_lockstep_width_res_t bwa3_resolve_smem_width_from_probe(int32_t raw_mlp,
                                                             int32_t host_width) {
    if (raw_mlp > 0) {
        bwa3_lockstep_width_res_t r = {bwa3_lockstep_width_from_probe(raw_mlp),
                                       BWA3_LOCKSTEP_WIDTH_SRC_PROBE};
        return r;
    }
    return bwa3_host_res(host_width, SMEM_LOCKSTEP_N, SMEM_LOCKSTEP_N_MAX);
}

bwa3_lockstep_width_res_t bwa3_resolve_bwtseed_width(int32_t pinned,
                                                     int32_t host_width) {
    if (pinned > 0) {
        bwa3_lockstep_width_res_t r = {pinned, BWA3_LOCKSTEP_WIDTH_SRC_ENV};
        return r;
    }
    return bwa3_host_res(host_width, BWTSEED_LOCKSTEP_N, BWTSEED_LOCKSTEP_N_MAX);
}

const char *bwa3_lockstep_width_source_name(int src) {
    switch (src) {
    case BWA3_LOCKSTEP_WIDTH_SRC_ENV:   return "env pin";
    case BWA3_LOCKSTEP_WIDTH_SRC_PROBE: return "startup probe";
    case BWA3_LOCKSTEP_WIDTH_SRC_HOST:  return "host default: AMD Zen 5";
    default:                            return "default";
    }
}

/* Runtime phase-2 lockstep width. Defaults to the compile-time floor so a
 * binary that never calls the initializer is identical to the old constant. */
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
     * caller fall back to the host default. */
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
    if (pinned < 0)
        /* Set but invalid: do NOT silently accept it. Diagnose, then resolve
         * exactly as an unset value. */
        fprintf(stderr,
                "ERROR: BWA3_SMEM_LOCKSTEP_N=\"%s\" is not a positive integer "
                "(<= %d); ignoring it (resolving as if unset).\n",
                env, SMEM_LOCKSTEP_N_MAX);
    /* The startup MLP probe is OPT-IN, not the default: it costs a
     * single-threaded startup pass and picks the per-access latency knee, which
     * sits below the whole-aligner wall optimum on hosts that keep gaining from
     * a wider width (Zen 5: knee 16-32 across startups, wall optimum 32). So by
     * default the width is the host-tuned value where one was measured, else the
     * compile-time SMEM_LOCKSTEP_N. BWA3_SMEM_LOCKSTEP_PROBE=1 opts back into
     * the probe to self-calibrate on new/untested hardware. Truthy opt-in
     * (bwa3_lockstep_probe_enabled): any value other than unset, empty, or "0"
     * enables it, so unset it (or set it to 0) to disable. */
    const int probe_on = bwa3_lockstep_probe_enabled(getenv("BWA3_SMEM_LOCKSTEP_PROBE"));
    int is_amd = 0;
    uint32_t family = 0, model = 0;
    bwa3_host_cpu_identity(&is_amd, &family, &model);
    const int32_t host_width = bwa3_smem_host_default_width(is_amd, family, model);
    /* A valid pin is taken alone, so gated/CI runs pay no measurement cost and
     * stay deterministic; the probe runs only when the resolver asks for it. */
    bwa3_lockstep_width_res_t r = bwa3_resolve_smem_width(pinned, probe_on, host_width);
    if (r.src == BWA3_LOCKSTEP_WIDTH_SRC_PROBE)
        r = bwa3_resolve_smem_width_from_probe(
            bwa3_measure_mlp(base, n_blocks, stride, word_off), host_width);
    g_smem_lockstep_n = r.width;
    s_smem_lockstep_n_source = r.src;
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
 * so a binary that never calls the initializer is identical to the constant. */
int32_t g_bwtseed_lockstep_n = BWTSEED_LOCKSTEP_N;

int32_t bwa3_bwtseed_width_parse_env(const char *env) {
    if (env == NULL || env[0] == '\0') return 0;  /* unset/empty: use the host default */

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
        if (pinned < 0)
            fprintf(stderr,
                    "ERROR: BWA3_BWTSEED_LOCKSTEP_N=\"%s\" is not a positive integer "
                    "(<= %d); ignoring it (resolving as if unset).\n",
                    env, BWTSEED_LOCKSTEP_N_MAX);
        int is_amd = 0;
        uint32_t family = 0, model = 0;
        bwa3_host_cpu_identity(&is_amd, &family, &model);
        const bwa3_lockstep_width_res_t r = bwa3_resolve_bwtseed_width(
            pinned, bwa3_bwtseed_host_default_width(is_amd, family, model));
        g_bwtseed_lockstep_n = r.width;   /* a pin is also the gate/CI path */
        s_bwtseed_lockstep_n_source = r.src;
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
