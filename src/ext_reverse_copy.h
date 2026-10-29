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
 * held most of the copy's samples on x86. Lengths 4 to 31 take two overlapping 16-, 8- or 4-byte
 * reversals, and lengths 1 to 3 move single bytes. Every store writes the final value of each
 * byte it covers, so the overlap is exact. The block shuffles and the unaligned loads and stores
 * are seq_reverse.h's. Other targets keep the byte loop. */
#ifndef BWA_MEM3_EXT_REVERSE_COPY_H
#define BWA_MEM3_EXT_REVERSE_COPY_H

#include <stdint.h>
#include <string.h>

#if defined(__AVX2__)
#include "seq_reverse.h"

/* dst[i] = src[m - 1 - i] for 0 < m < 32 (ext_reverse_copy's short case; the spans are disjoint).
 * Both end loads are made before either store, and each reversed load is stored at the other end,
 * so the overlap of the two stores in the middle is exact. */
static inline void ext_reverse_short_avx2(uint8_t *__restrict dst, const uint8_t *__restrict src,
                                          int64_t m)
{
    if (m >= 16) {
        const seq_rev_v16 a = seq_rev_ld16(src), b = seq_rev_ld16(src + m - 16);
        seq_rev_st16(dst, seq_rev_rev16(b));
        seq_rev_st16(dst + m - 16, seq_rev_rev16(a));
    } else if (m >= 8) {
        const uint64_t a = seq_rev_ld8(src), b = seq_rev_ld8(src + m - 8);
        seq_rev_st8(dst, __builtin_bswap64(b));
        seq_rev_st8(dst + m - 8, __builtin_bswap64(a));
    } else if (m >= 4) {
        const uint32_t a = seq_rev_ld4(src), b = seq_rev_ld4(src + m - 4);
        seq_rev_st4(dst, __builtin_bswap32(b));
        seq_rev_st4(dst + m - 4, __builtin_bswap32(a));
    } else if (m >= 2) {
        const uint8_t a = src[0], b = src[m - 1], c = src[1];
        dst[0] = b; dst[m - 1] = a;
        if (m == 3) dst[1] = c;
    } else if (m == 1) {
        dst[0] = src[0];
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
                                seq_rev_rev32(_mm256_loadu_si256((const __m256i *)(src + n - 32 - i))));
        if (i < n)
            _mm256_storeu_si256((__m256i *)(dst + n - 32),
                                seq_rev_rev32(_mm256_loadu_si256((const __m256i *)src)));
    } else if (n > 0) {
        ext_reverse_short_avx2(dst, src, n);
    }
#else
    for (int64_t i = 0; i < n; ++i) dst[i] = src[n - 1 - i];
#endif
}

#endif
