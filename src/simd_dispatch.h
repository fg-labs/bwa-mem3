// src/simd_dispatch.h
#ifndef BWAMEM3_SIMD_DISPATCH_H
#define BWAMEM3_SIMD_DISPATCH_H

#include <stdio.h>          /* for FILE * in bwamem3_print_version_simd */
#include <stddef.h>         /* for size_t in bwamem3_format_host_floor_error */

#ifdef __cplusplus
extern "C" {
#endif

/* SIMD tier enum. Higher numeric values are strictly more capable: any
 * code path that runs at tier T also runs at tier > T. The dispatcher
 * uses this ordering for the BWAMEM3_FORCE_TIER downgrade check. */
enum bwamem3_tier {
    BWAMEM3_TIER_NONE     = 0,   /* scalar-only fallback */
    BWAMEM3_TIER_SSE41    = 1,
    BWAMEM3_TIER_SSE42    = 2,
    BWAMEM3_TIER_AVX      = 3,
    BWAMEM3_TIER_AVX2     = 4,
    BWAMEM3_TIER_AVX512BW = 5,
    BWAMEM3_TIER_NEON     = 6    /* arm64; not ordered against x86 tiers */
};

/* Initialize the SIMD dispatcher. Idempotent; safe to call multiple times.
 * Must be called before any kernel factory (make_banded_pair_wise_sw,
 * make_kswv) or sam_encode_* function-pointer call.
 *
 * Thread-safety: the implementation uses std::call_once, so concurrent calls
 * from multiple threads are safe and at most one thread runs the body. The
 * dispatch wrappers in simd_dispatch.cpp also defensively call this on first
 * use, so worker threads that race ahead of main() are still correct. The
 * intended usage pattern is still to call this once from the main thread
 * before spawning workers, to keep startup-time error messages (FORCE_TIER
 * warnings, debug-tier banner) on the main thread.
 */
void bwamem3_simd_init(void);

/* Returns the active tier. Valid only after bwamem3_simd_init() has run; if
 * called before init, returns BWAMEM3_TIER_NONE. */
int bwamem3_simd_tier(void);

/* Returns a stable string for a tier value: "sse41", "sse42", "avx",
 * "avx2", "avx512bw", "neon", "scalar", or "unknown". */
const char *bwamem3_simd_tier_name(int tier);

/* Pure comparison helper: returns 1 if host_tier can run a binary built
 * at build_tier, 0 otherwise. Handles NEON/x86 family orthogonality
 * defensively. Declared in the header so unit tests don't need to
 * forward-declare it (which would let signature drift go undetected). */
int bwamem3_check_host_floor(int host_tier, int build_tier);

/* Pure formatter: writes the host-floor refusal message into a caller-
 * provided buffer. Returns the byte count written (excluding terminator),
 * or -1 on truncation/null-buffer. Buffer is NUL-terminated within its
 * bounds when buf != NULL && bufsz > 0. */
int bwamem3_format_host_floor_error(char *buf, size_t bufsz,
                                    int host_tier, int build_tier);

/* Returns 1 if the running host can execute this binary's compiled-in
 * instructions, 0 otherwise. Reads the raw host capability (independent
 * of BWAMEM3_FORCE_TIER) so the query reflects the SIGILL risk for the
 * compiler-emitted non-kernel TU instructions, not the kernel dispatch
 * tier. Calls bwamem3_simd_init() to populate the cache. Pure query —
 * never exits. Used by the version subcommand to print a warning. */
int bwamem3_host_meets_floor(void);

/* Refuses to run if the host doesn't meet the build's floor: writes a clear
 * error to stderr and calls exit(2). On success, returns. Called from
 * main() after the early -h/--help/version short-circuits so diagnostic
 * commands (`bwa-mem3 version`, `bwa-mem3 <cmd> --help`) remain usable
 * on too-old hosts. */
void bwamem3_enforce_host_floor(void);

/* Prints the SIMD floor / runtime / optional warning banner to the given
 * stream. Calls bwamem3_simd_init() to ensure tiers are populated. Two
 * lines unconditionally (floor + runtime); a third warning line iff
 * bwamem3_host_meets_floor() returns 0. Never exits. */
void bwamem3_print_version_simd(FILE *f);

#ifdef __cplusplus
}
#endif

#endif /* BWAMEM3_SIMD_DISPATCH_H */
