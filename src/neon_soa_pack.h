/* SPDX-License-Identifier: MIT */
/* neon_soa_pack.h -- tiled structure-of-arrays packing for the 16-lane 8-bit
 * and 8-lane 16-bit NEON Smith-Waterman kernels.
 *
 * The 128-bit kernels read their sequences as SoA: row k holds position k of
 * every lane -- 16 lanes of one byte for the 8-bit kernels, 8 lanes of one
 * halfword for the 16-bit kernels. Filling that layout one base at a time is a
 * strided scatter the compiler cannot vectorize (a load/compare/select/store
 * per base). This header packs a whole tile of positions of all lanes at once
 * instead: one contiguous load per lane, an in-register transpose (zip
 * stages), and one store per SoA row.
 *
 * NEON-only; x86_soa_pack.h is the x86 counterpart. */
#ifndef BWAMEM3_NEON_SOA_PACK_H
#define BWAMEM3_NEON_SOA_PACK_H

#if defined(__ARM_NEON) || defined(__aarch64__)
#include <arm_neon.h>
#include <stddef.h>
#include <stdint.h>
#ifndef __cplusplus
#include <stdbool.h>   /* self-contained bool for a hypothetical C includer */
#endif
#include <assert.h>
#include <string.h>
#include "neon_transpose.h"

static const uint8_t neon_soa_iota16[16] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};

/* The first r = len - kb bytes of a boundary tile (positions [kb, len) of a sequence, 1 <= r <= 15)
 * in bytes [0, r) of a vector, read without touching memory outside [seq, seq + len); bytes
 * [r, 16) are unspecified. With len >= 16 the tile's bases are the last r bytes of the 16 that
 * end at seq + len: one load and one table lookup. A shorter sequence (then kb == 0) is read as
 * its first and last a bytes, a the largest power of two <= len, which together cover it. */
static inline uint8x16_t neon_soa_tail16(const uint8_t *seq, int len, int kb)
{
    const uint8x16_t iota = vld1q_u8(neon_soa_iota16);
    const int r = len - kb;
    if (len >= 16)
        return vqtbl1q_u8(vld1q_u8(seq + len - 16), vaddq_u8(iota, vdupq_n_u8((uint8_t) (16 - r))));
    /* len < 16 here, so kb == 0 (the caller's r > 0 means kb < len); no assert: on this hot path
     * it changes the pack's inlining and measured 3% slower on AVX-512BW */
    uint64_t lo, hi;
    int a;
    if (len >= 8) {
        a = 8; memcpy(&lo, seq, 8); memcpy(&hi, seq + len - 8, 8);
    } else if (len >= 4) {
        uint32_t x, y; a = 4; memcpy(&x, seq, 4); memcpy(&y, seq + len - 4, 4); lo = x; hi = y;
    } else if (len >= 2) {
        uint16_t x, y; a = 2; memcpy(&x, seq, 2); memcpy(&y, seq + len - 2, 2); lo = x; hi = y;
    } else {
        a = 1; lo = seq[0]; hi = seq[0];
    }
    /* byte t < a is byte t of lo; byte a <= t < len is byte t - (len - a) of hi, at 8 + t - len + a */
    const uint8x16_t v = vcombine_u8(vcreate_u8(lo), vcreate_u8(hi));
    const uint8x16_t idx = vbslq_u8(vcltq_u8(iota, vdupq_n_u8((uint8_t) a)), iota,
                                    vaddq_u8(iota, vdupq_n_u8((uint8_t) (8 - len + a))));
    return vqtbl1q_u8(v, idx);
}

/* Pack 16 lanes' byte sequences into the 16-lane SoA layout the u8 kernels
 * read (row k = position k of every lane), 16 positions per transposed tile
 * instead of one strided byte store per base. Lane j's row k is
 *     seq[j][k]                     for k <  len[j]        (4 -> 8 when remap4to8)
 *     padA                          for len[j] <= k < padStart[j]
 *     padB                          for k >= padStart[j]
 * and rows [0, nrows) are written; the ranges are ordered, so len[j] <=
 * padStart[j] is required (a tile past padStart is a padB broadcast). Full 16-byte tiles inside a lane are one
 * vector load; a lane's boundary tile is one or two loads that stay inside its
 * sequence (neon_soa_tail16) blended with the pad bytes, and tiles entirely past
 * padStart are a broadcast. The remap is applied only to real bases (the
 * boundary tile's pad bytes are blended in after it), so pad bytes are never
 * rewritten even if padA or padB were 4. Byte-for-byte the same SoA the scalar
 * fill produced. */
static inline void neon_soa_pack16(uint8_t *soa,
                                   const uint8_t *const seq[SIMD_WIDTH8],
                                   const int len[SIMD_WIDTH8],
                                   const int padStart[SIMD_WIDTH8],
                                   int nrows, uint8_t padA, uint8_t padB,
                                   bool remap4to8)
{
    static_assert(SIMD_WIDTH8 == 16,
                  "neon_soa_pack16 hardcodes 16-wide tiles (t[16]/kb += 16); "
                  "SIMD_WIDTH8 must be 16 on the NEON tier");
    const uint8x16_t four  = vdupq_n_u8(AMBIG_);
    const uint8x16_t eight = vdupq_n_u8(AMBQ);
    const uint8x16_t padAv = vdupq_n_u8(padA);
    const uint8x16_t padBv = vdupq_n_u8(padB);
    const uint8x16_t iota = vld1q_u8(neon_soa_iota16);
    for (int j = 0; j < SIMD_WIDTH8; j++) assert(len[j] <= padStart[j]);   /* see the contract above */
    for (int kb = 0; kb < nrows; kb += 16) {
        uint8x16_t t[16];
        for (int j = 0; j < SIMD_WIDTH8; j++) {
            uint8x16_t v;
            if (kb + 16 <= len[j]) {
                v = vld1q_u8(seq[j] + kb);
                if (remap4to8) v = vbslq_u8(vceqq_u8(v, four), eight, v);
            } else if (kb >= padStart[j]) {
                v = padBv;
            } else {
                /* boundary tile: r real bases (0 <= r <= 15), then padA up to padStart, then padB */
                const int r = len[j] > kb ? len[j] - kb : 0;
                const int pa = padStart[j] - kb < 16 ? padStart[j] - kb : 16;   /* > 0 here */
                uint8x16_t real = r > 0 ? neon_soa_tail16(seq[j], len[j], kb) : padBv;
                if (remap4to8) real = vbslq_u8(vceqq_u8(real, four), eight, real);
                const uint8x16_t padv = vbslq_u8(vcltq_u8(iota, vdupq_n_u8((uint8_t) pa)), padAv, padBv);
                v = vbslq_u8(vcltq_u8(iota, vdupq_n_u8((uint8_t) r)), real, padv);
            }
            t[j] = v;
        }
        neon_transpose16x16_u8(t);
        const int kend = (kb + 16 < nrows) ? kb + 16 : nrows;
        for (int k = kb; k < kend; k++)
            vst1q_u8(soa + (size_t) k * SIMD_WIDTH8, t[k - kb]);
    }
}

/* 16-bit twin of neon_soa_pack16 for the 8-lane int16 SoA: 8 positions per
 * tile, each lane's 8 bytes loaded and widened to halfwords, the ambiguity
 * code remapped (AMBIG_ -> ambCode) on whole tile vectors (no pad code equals
 * 4), transposed, and stored as 8 SoA rows. Lane j's row k is
 *     seq[j][k]  (4 -> ambCode)   for k <  len[j]
 *     padA                        for len[j] <= k < padStart[j]
 *     padB                        for k >= padStart[j]
 * and rows [0, nrows) are written. Byte-for-byte the scalar fill's output. */
static inline void neon_soa_pack8_u16(int16_t *soa,
                                      const uint8_t *const seq[SIMD_WIDTH16],
                                      const int len[SIMD_WIDTH16],
                                      const int padStart[SIMD_WIDTH16],
                                      int nrows, uint16_t padA, uint16_t padB,
                                      uint16_t ambCode)
{
    static_assert(SIMD_WIDTH16 == 8,
                  "neon_soa_pack8_u16 hardcodes 8-wide tiles (t[8]/tmp[8]/kb += 8); "
                  "SIMD_WIDTH16 must be 8 on the NEON tier");
    const uint16x8_t four  = vdupq_n_u16(AMBIG_);
    const uint16x8_t ambv  = vdupq_n_u16(ambCode);
    const uint16x8_t padBv = vdupq_n_u16(padB);
    for (int kb = 0; kb < nrows; kb += 8) {
        uint16x8_t t[8];
        for (int j = 0; j < SIMD_WIDTH16; j++) {
            uint16x8_t v;
            if (kb + 8 <= len[j]) {
                v = vmovl_u8(vld1_u8(seq[j] + kb));
                v = vbslq_u16(vceqq_u16(v, four), ambv, v);
            } else if (kb >= padStart[j]) {
                v = padBv;
            } else {
                uint16_t tmp[8];
                for (int t2 = 0; t2 < 8; t2++) {
                    const int k = kb + t2;
                    tmp[t2] = k < len[j] ? (seq[j][k] == AMBIG_ ? ambCode : (uint16_t) seq[j][k])
                            : (k < padStart[j] ? padA : padB);
                }
                v = vld1q_u16(tmp);
            }
            t[j] = v;
        }
        neon_transpose8x8_u16(t);
        const int kend = (kb + 8 < nrows) ? kb + 8 : nrows;
        for (int k = kb; k < kend; k++)
            vst1q_s16(soa + (size_t) k * SIMD_WIDTH16, vreinterpretq_s16_u16(t[k - kb]));
    }
}

#endif  /* __ARM_NEON || __aarch64__ */

#endif  /* BWAMEM3_NEON_SOA_PACK_H */
