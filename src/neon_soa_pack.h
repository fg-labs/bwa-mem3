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
 * NEON-only; the x86 tiers of the same source keep their scalar loops. */
#ifndef BWAMEM3_NEON_SOA_PACK_H
#define BWAMEM3_NEON_SOA_PACK_H

#if defined(__ARM_NEON) || defined(__aarch64__)
#include <arm_neon.h>
#include <stddef.h>
#include <stdint.h>
#ifndef __cplusplus
#include <stdbool.h>   /* self-contained bool for a hypothetical C includer */
#endif
#include "neon_transpose.h"

/* Pack 16 lanes' byte sequences into the 16-lane SoA layout the u8 kernels
 * read (row k = position k of every lane), 16 positions per transposed tile
 * instead of one strided byte store per base. Lane j's row k is
 *     seq[j][k]                     for k <  len[j]        (4 -> 8 when remap4to8)
 *     padA                          for len[j] <= k < padStart[j]
 *     padB                          for k >= padStart[j]
 * and rows [0, nrows) are written. Full 16-byte tiles inside a lane are one
 * vector load; a lane's boundary tile is assembled byte-wise (no read past
 * its sequence), and tiles entirely past padStart are a broadcast. The remap
 * is applied only to real bases -- on the whole-tile vector (all bases) and, in
 * the boundary tile, per byte for k < len[j] -- so pad bytes are never rewritten
 * even if padA or padB were 4. Byte-for-byte the same SoA the scalar fill
 * produced. */
static inline void neon_soa_pack16(uint8_t *soa,
                                   const uint8_t *const seq[SIMD_WIDTH8],
                                   const int len[SIMD_WIDTH8],
                                   const int padStart[SIMD_WIDTH8],
                                   int nrows, uint8_t padA, uint8_t padB,
                                   bool remap4to8)
{
    static_assert(SIMD_WIDTH8 == 16,
                  "neon_soa_pack16 hardcodes 16-wide tiles (t[16]/tmp[16]/kb += 16); "
                  "SIMD_WIDTH8 must be 16 on the NEON tier");
    const uint8x16_t four  = vdupq_n_u8(AMBIG_);
    const uint8x16_t eight = vdupq_n_u8(AMBQ);
    const uint8x16_t padBv = vdupq_n_u8(padB);
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
                uint8_t tmp[16];
                for (int t2 = 0; t2 < 16; t2++) {
                    const int k = kb + t2;
                    if (k < len[j])
                        tmp[t2] = (remap4to8 && seq[j][k] == AMBIG_) ? AMBQ : seq[j][k];
                    else
                        tmp[t2] = (k < padStart[j]) ? padA : padB;
                }
                v = vld1q_u8(tmp);
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
