/* ext_reverse_copy.h -- the reversed left-extension windows of stage_seed_extension
 * (bwamem.cpp): dst[i] = src[n - 1 - i] for i in [0, n); n <= 0 copies nothing.
 *
 * n is a by-value parameter, so the loop bound is invariant and the compiler can vectorize the
 * loop (clang 19 on aarch64 reverses 32 bytes per iteration with rev64 + ext); __restrict spares
 * it a runtime overlap check. Read through the window struct instead, the bound must be reloaded
 * after every byte store, which may alias it, and the loop stays scalar. dst and src must not
 * overlap (the staging buffers vs the read or the fetched reference).
 *
 * Under AVX2 the copy is written out in 32-byte blocks (vpshufb within each 128-bit lane, then
 * vpermq to swap the lanes) and ends with one overlapping block instead of a byte loop: clang 19
 * vectorizes the plain loop at 16 bytes per iteration with a byte tail of up to 15, and the tail
 * held most of the copy's samples on x86. Lengths below 32 take two overlapping 16-, 8- or 4-byte
 * reversals. Every store writes the final value of each byte it covers, so the overlap is exact.
 * ext_reverse_inplace_avx2 is the in-place form (bwa_gen_cigar4's reverse-strand query and
 * reference, bwa.cpp; mate rescue's revseq uses seq_reverse.h). Other targets keep the byte loop. */
#ifndef BWA_MEM3_EXT_REVERSE_COPY_H
#define BWA_MEM3_EXT_REVERSE_COPY_H

#include <stdint.h>
#include <string.h>

#if defined(__AVX2__)
#include <immintrin.h>

static inline __m256i ext_rev32_avx2(__m256i v)
{
    const __m256i lane_rev = _mm256_setr_epi8(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0,
                                              15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
    return _mm256_permute4x64_epi64(_mm256_shuffle_epi8(v, lane_rev), 0x4e);
}

static inline __m128i ext_rev16_avx2(__m128i v)
{
    return _mm_shuffle_epi8(v, _mm_setr_epi8(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0));
}

/* [dlo, dhi) = [slo, shi) reversed, for a span of 0 < m < 32 bytes; the spans are either the same
 * (in place) or disjoint. Both end loads are made before either store, and each reversed load is
 * stored at the other end, so the overlap of the two stores in the middle is exact. */
static inline void ext_reverse_short_avx2(uint8_t *dlo, uint8_t *dhi, const uint8_t *slo,
                                          const uint8_t *shi, int64_t m)
{
    if (m >= 16) {
        const __m128i a = _mm_loadu_si128((const __m128i *)slo);
        const __m128i b = _mm_loadu_si128((const __m128i *)(shi - 16));
        _mm_storeu_si128((__m128i *)dlo, ext_rev16_avx2(b));
        _mm_storeu_si128((__m128i *)(dhi - 16), ext_rev16_avx2(a));
    } else if (m >= 8) {
        uint64_t a, b;
        memcpy(&a, slo, 8); memcpy(&b, shi - 8, 8);
        a = __builtin_bswap64(a); b = __builtin_bswap64(b);
        memcpy(dlo, &b, 8); memcpy(dhi - 8, &a, 8);
    } else if (m >= 4) {
        uint32_t a, b;
        memcpy(&a, slo, 4); memcpy(&b, shi - 4, 4);
        a = __builtin_bswap32(a); b = __builtin_bswap32(b);
        memcpy(dlo, &b, 4); memcpy(dhi - 4, &a, 4);
    } else if (m >= 2) {
        const uint8_t a = slo[0], b = shi[-1], c = slo[1];
        dlo[0] = b; dhi[-1] = a;
        if (m == 3) dlo[1] = c;
    } else if (m == 1) {
        dlo[0] = slo[0];
    }
}

/* s[i] <-> s[n - 1 - i] for i < n / 2; n <= 1 changes nothing. One 32-byte block from each end
 * per step toward the middle, then the middle (< 64 bytes) by two overlapping end loads. */
static inline void ext_reverse_inplace_avx2(uint8_t *s, int64_t n)
{
    uint8_t *lo = s, *hi = s + n;
    while (hi - lo >= 64) {
        const __m256i a = _mm256_loadu_si256((const __m256i *)lo);
        const __m256i b = _mm256_loadu_si256((const __m256i *)(hi - 32));
        _mm256_storeu_si256((__m256i *)lo, ext_rev32_avx2(b));
        _mm256_storeu_si256((__m256i *)(hi - 32), ext_rev32_avx2(a));
        lo += 32; hi -= 32;
    }
    const int64_t m = hi - lo;
    if (m >= 32) {
        const __m256i a = _mm256_loadu_si256((const __m256i *)lo);
        const __m256i b = _mm256_loadu_si256((const __m256i *)(hi - 32));
        _mm256_storeu_si256((__m256i *)lo, ext_rev32_avx2(b));
        _mm256_storeu_si256((__m256i *)(hi - 32), ext_rev32_avx2(a));
    } else if (m >= 2) {
        ext_reverse_short_avx2(lo, hi, lo, hi, m);
    }
}
#endif

static inline void ext_reverse_copy(uint8_t *__restrict dst, const uint8_t *__restrict src, int64_t n)
{
#if defined(__AVX2__)
    if (n >= 32) {
        int64_t i = 0;
        for (; i + 32 <= n; i += 32)
            _mm256_storeu_si256((__m256i *)(dst + i),
                                ext_rev32_avx2(_mm256_loadu_si256((const __m256i *)(src + n - 32 - i))));
        if (i < n)
            _mm256_storeu_si256((__m256i *)(dst + n - 32),
                                ext_rev32_avx2(_mm256_loadu_si256((const __m256i *)src)));
    } else if (n > 0) {
        ext_reverse_short_avx2(dst, dst + n, src, src + n, n);
    }
#else
    for (int64_t i = 0; i < n; ++i) dst[i] = src[n - 1 - i];
#endif
}

#endif
