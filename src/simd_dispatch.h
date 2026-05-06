// src/simd_dispatch.h
#ifndef BWAMEM3_SIMD_DISPATCH_H
#define BWAMEM3_SIMD_DISPATCH_H

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

#ifdef __cplusplus
}
#endif

#endif /* BWAMEM3_SIMD_DISPATCH_H */
