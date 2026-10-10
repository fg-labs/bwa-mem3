/* ext_reverse_copy.h -- the reversed left-extension windows of stage_seed_extension
 * (bwamem.cpp): dst[i] = src[n - 1 - i] for i in [0, n).
 *
 * The old in-place byte loops ran one byte per iteration: their bound was a field of the
 * window struct reached through a reference, and a byte store may alias it, so the bound was
 * reloaded every byte and the loop could not vectorize. Here the bound is a local, and on
 * aarch64 the copy runs 16 bytes at a time (rev64 reverses each 8-byte half, ext swaps the
 * halves); the rest, and other targets, use the byte loop, which the compiler may now
 * vectorize. The result is the same bytes either way. dst and src must not overlap (the
 * staging buffers vs the read or the fetched reference). */
#ifndef BWA_MEM3_EXT_REVERSE_COPY_H
#define BWA_MEM3_EXT_REVERSE_COPY_H

#include <stdint.h>
#if defined(__aarch64__)
#include <arm_neon.h>
#endif

static inline void ext_reverse_copy(uint8_t *__restrict dst, const uint8_t *__restrict src, int64_t n)
{
    int64_t i = 0;
#if defined(__aarch64__)
    for (; i + 16 <= n; i += 16) {
        const uint8x16_t v = vrev64q_u8(vld1q_u8(src + n - 16 - i));
        vst1q_u8(dst + i, vextq_u8(v, v, 8));
    }
#endif
    for (; i < n; ++i) dst[i] = src[n - 1 - i];
}

#endif
