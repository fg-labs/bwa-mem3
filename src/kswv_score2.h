/* score2/te2 plateau scan of the 8-bit batched mate-rescue kernels (kswv).
 *
 * After the forward DP, every u8 kernel (NEON, AVX2, AVX-512BW) emulates scalar
 * ksw_u8's b[] build and score2 scan per lane over the row maxima it stored in
 * rowMax (row-major, W lanes per row): consecutive rows at or above minsc collapse
 * into one entry anchored at the run's first strictly-greatest row, and an entry
 * anchored inside the primary region [low, high] is excluded.
 *
 * Almost every row is below minsc and takes the scan's `continue`, which changes
 * no scan state; a dense walk would pay a strided load, a compare and a branch
 * for it in every lane. So the kernels first build a per-lane bitmap of the
 * rows that clear minsc with whole-row compares -- eight rows per byte, all lanes
 * at once -- and each lane then visits only its set bits, in ascending row order.
 * Every visited row is re-tested against minsc, so the bitmap only has to be a
 * superset of the qualifying rows, and the visited rows run the exact sequence of
 * state updates the dense walk does: the same score2/te2.
 *
 * Its own header so the unit tests can check the bitmap scan against the dense
 * one directly. */
#ifndef BWA_MEM3_KSWV_SCORE2_H
#define BWA_MEM3_KSWV_SCORE2_H

#include <stddef.h>
#include <stdint.h>
#include <string.h>

/* Bitmap bytes per lane (one per block of 8 rows) for nrows rows. */
static inline int kswv_score2_bits_blocks(int nrows) { return (nrows + 7) >> 3; }

/* bits[b * W + l] bit k is set iff row 8b + k (< nrows) of lane l has
 * rowMax[row * W + l] >= minsc[l]. bits must hold kswv_score2_bits_blocks(nrows) * W
 * bytes. The lane loop has a compile-time trip count, so it compiles to one
 * vector compare, mask and OR per row. */
template <int W>
static inline void kswv_score2_row_bits_u8(const uint8_t *rowMax, int nrows,
                                           const uint8_t *minsc, uint8_t *bits)
{
    const int nblk = kswv_score2_bits_blocks(nrows);
    for (int b = 0; b < nblk; b++) {
        uint8_t acc[W];
        for (int l = 0; l < W; l++) acc[l] = 0;
        const int kend = nrows - 8 * b < 8 ? nrows - 8 * b : 8;
        for (int k = 0; k < kend; k++) {
            const uint8_t *row = rowMax + (size_t)(8 * b + k) * W;
            const uint8_t bit = (uint8_t)(1u << k);
            for (int l = 0; l < W; l++)
                acc[l] |= (uint8_t)(row[l] >= minsc[l] ? bit : 0);
        }
        for (int l = 0; l < W; l++) bits[(size_t)b * W + l] = acc[l];
    }
}

/* Every row of every lane marked (BWA3_RESCUE_SCORE2_BITS=0): the scan below then visits
 * all of a lane's rows, the dense walk. */
template <int W>
static inline void kswv_score2_row_bits_all(int nrows, uint8_t *bits)
{
    memset(bits, 0xFF, (size_t)kswv_score2_bits_blocks(nrows) * W);
}

/* The per-lane scan over rows [0, nrows) of lane `lane`, visiting only the rows
 * set in bits. bits must cover at least nrows rows and include every row of the
 * lane with rowMax >= minsc (kswv_score2_row_bits_u8 at the same minsc, or any
 * superset: extra rows are re-tested and skipped). Writes -1/-1 when no entry
 * outside [low, high] qualifies. */
template <int W>
static inline void kswv_score2_scan_bits_u8(const uint8_t *rowMax, const uint8_t *bits,
                                            int nrows, int lane, int minsc, int low,
                                            int high, int *score2_out, int *te2_out)
{
    int score2 = -1, te2 = -1;
    /* b_pos = -2 sentinel: (-2 + 1) != any real row, forces an "append" on the
     * first qualifying row. */
    int b_score = -1, b_pos = -2;
    const int nblk = kswv_score2_bits_blocks(nrows);
    for (int b = 0; b < nblk; b++) {
        unsigned m = bits[(size_t)b * W + lane];
        while (m) {
            const int i2 = 8 * b + __builtin_ctz(m);
            m &= m - 1;
            if (i2 >= nrows) break;
            const int imax = (int)rowMax[(size_t)i2 * W + lane];
            if (imax < minsc) continue;

            if (b_pos + 1 != i2) {
                /* APPEND: flush the outgoing b[] entry to score2 first. */
                if (b_pos >= 0 && (b_pos < low || b_pos > high) && b_score > score2) {
                    score2 = b_score;
                    te2    = b_pos;
                }
                b_score = imax;
                b_pos   = i2;
            } else if (b_score < imax) {
                /* UPDATE: strict greater extends the run's anchor. Equal imax
                 * leaves b_pos in place (matches scalar). */
                b_score = imax;
                b_pos   = i2;
            }
        }
    }

    /* Flush trailing b[] entry. */
    if (b_pos >= 0 && (b_pos < low || b_pos > high) && b_score > score2) {
        score2 = b_score;
        te2    = b_pos;
    }
    *score2_out = score2;
    *te2_out    = te2;
}

#endif
