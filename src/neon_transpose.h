/* SPDX-License-Identifier: MIT */
/* neon_transpose.h -- in-register 16x16 byte and 8x8 halfword transposes, shared by the NEON SoA
 * packing (neon_soa_pack.h), the banded mate-rescue kernel (rescue_band.cpp) and the NEON rescue
 * filter's segment scans (rescue_prune_neon.h).
 *
 * Self-contained, unlike neon_soa_pack.h (which needs the kernels' SIMD_WIDTH8 / AMBIG_ / AMBQ
 * macros), so a translation unit that only transposes can include it on its own. NEON-only. */
#ifndef BWAMEM3_NEON_TRANSPOSE_H
#define BWAMEM3_NEON_TRANSPOSE_H

#if defined(__ARM_NEON) || defined(__aarch64__)
#include <arm_neon.h>

/* In-register 16x16 byte transpose: on entry r[j] holds 16 consecutive bytes
 * of lane j; on return r[k] holds byte k of every lane (the SoA row k). Four
 * zip stages -- bytes, halfwords, words, doublewords -- of 16 ops each. */
static inline void neon_transpose16x16_u8(uint8x16_t r[16])
{
    uint8x16_t a[16];
    for (int p = 0; p < 8; p++) {                 /* lanes (2p, 2p+1): k 0-7 | 8-15 */
        a[2 * p]     = vzip1q_u8(r[2 * p], r[2 * p + 1]);
        a[2 * p + 1] = vzip2q_u8(r[2 * p], r[2 * p + 1]);
    }
    uint16x8_t b[16];
    for (int g = 0; g < 4; g++) {                 /* lanes 4g..4g+3: k 0-3 | 4-7 | 8-11 | 12-15 */
        const uint16x8_t lo0 = vreinterpretq_u16_u8(a[4 * g]),     lo1 = vreinterpretq_u16_u8(a[4 * g + 2]);
        const uint16x8_t hi0 = vreinterpretq_u16_u8(a[4 * g + 1]), hi1 = vreinterpretq_u16_u8(a[4 * g + 3]);
        b[4 * g]     = vzip1q_u16(lo0, lo1);
        b[4 * g + 1] = vzip2q_u16(lo0, lo1);
        b[4 * g + 2] = vzip1q_u16(hi0, hi1);
        b[4 * g + 3] = vzip2q_u16(hi0, hi1);
    }
    uint32x4_t c[16];
    for (int h = 0; h < 2; h++) {                 /* lanes 8h..8h+7: k pairs (0,1) .. (14,15) */
        for (int q = 0; q < 4; q++) {
            const uint32x4_t x = vreinterpretq_u32_u16(b[8 * h + q]);
            const uint32x4_t y = vreinterpretq_u32_u16(b[8 * h + 4 + q]);
            c[8 * h + 2 * q]     = vzip1q_u32(x, y);
            c[8 * h + 2 * q + 1] = vzip2q_u32(x, y);
        }
    }
    for (int m = 0; m < 8; m++) {                 /* lanes 0-7 | 8-15 -> rows 2m, 2m+1 */
        const uint64x2_t x = vreinterpretq_u64_u32(c[m]);
        const uint64x2_t y = vreinterpretq_u64_u32(c[8 + m]);
        r[2 * m]     = vreinterpretq_u8_u64(vzip1q_u64(x, y));
        r[2 * m + 1] = vreinterpretq_u8_u64(vzip2q_u64(x, y));
    }
}

/* In-register 8x8 halfword transpose: on entry r[j] holds 8 consecutive
 * positions of lane j; on return r[k] holds position k of every lane (the
 * 16-bit SoA row k). Three zip stages -- halfwords, words, doublewords. */
static inline void neon_transpose8x8_u16(uint16x8_t r[8])
{
    uint16x8_t a[8];
    for (int p = 0; p < 4; p++) {                 /* lanes (2p, 2p+1): k 0-3 | 4-7 */
        a[2 * p]     = vzip1q_u16(r[2 * p], r[2 * p + 1]);
        a[2 * p + 1] = vzip2q_u16(r[2 * p], r[2 * p + 1]);
    }
    uint32x4_t b[8];
    for (int h = 0; h < 2; h++) {                 /* lanes 4h..4h+3: k pairs (0,1) .. (6,7) */
        const uint32x4_t lo0 = vreinterpretq_u32_u16(a[4 * h]),     lo1 = vreinterpretq_u32_u16(a[4 * h + 2]);
        const uint32x4_t hi0 = vreinterpretq_u32_u16(a[4 * h + 1]), hi1 = vreinterpretq_u32_u16(a[4 * h + 3]);
        b[4 * h]     = vzip1q_u32(lo0, lo1);
        b[4 * h + 1] = vzip2q_u32(lo0, lo1);
        b[4 * h + 2] = vzip1q_u32(hi0, hi1);
        b[4 * h + 3] = vzip2q_u32(hi0, hi1);
    }
    for (int m = 0; m < 4; m++) {                 /* lanes 0-3 | 4-7 -> rows 2m, 2m+1 */
        const uint64x2_t x = vreinterpretq_u64_u32(b[m]);
        const uint64x2_t y = vreinterpretq_u64_u32(b[4 + m]);
        r[2 * m]     = vreinterpretq_u16_u64(vzip1q_u64(x, y));
        r[2 * m + 1] = vreinterpretq_u16_u64(vzip2q_u64(x, y));
    }
}

/* The same permutation on signed halfwords (the rescue filter's int16 scan segments). */
static inline void neon_transpose8x8_s16(int16x8_t r[8])
{
    uint16x8_t u[8];
    for (int k = 0; k < 8; k++) u[k] = vreinterpretq_u16_s16(r[k]);
    neon_transpose8x8_u16(u);
    for (int k = 0; k < 8; k++) r[k] = vreinterpretq_s16_u16(u[k]);
}

#endif  /* __ARM_NEON || __aarch64__ */

#endif  /* BWAMEM3_NEON_TRANSPOSE_H */
