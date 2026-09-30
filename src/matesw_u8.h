/* u8 mate-rescue admission, shared by the aligner (bwamem_pair.cpp) and the banded-rescue harness
 * (test/rescue_band_harness.cpp), which must pick the 8-bit jobs exactly as the aligner does. */
#ifndef BWA_MEM3_MATESW_U8_H
#define BWA_MEM3_MATESW_U8_H

#include <stdint.h>
#include "bwamem.h"

/* ------------------------------------------------------------------------- *
 * u8 mate-rescue admission (tier selection)
 *
 * Mate rescue routes each pair to the 8-bit (u8) or 16-bit SW kernel. The u8
 * kernel biases scores by `shift` and stores them in a uint8_t, so it is only
 * byte-identical to the 16-bit kernel while the max attainable biased score
 * cannot saturate/early-break the u8 DP. The binding limit is the u8-only
 * early break `gmax + q->shift >= 255` in ksw.cpp (ksw_u8), i.e. the biased
 * max score must be <= 254; the batched kswv u8 kernel saturates one unit
 * later, so 254 is the tighter, shared bound.
 *
 * `shift` is the kernel's score bias: ksw derives it as `256 - min(matrix)`
 * taken as a uint8_t, which equals `-(min matrix entry)` for every matrix
 * bwa-mem3 builds (>= 1 in practice; the floor guards a degenerate matrix).
 *
 * There is a second, geometric provability condition. The u8 and 16-bit paths
 * pad the query to different widths (16*ceil(l/16) vs 8*ceil(l/8)); a pair may
 * flip tier only when those pad widths are equal, otherwise the two kernels run
 * different padded column counts and their tie-break/XS scan can diverge. That
 * equality holds exactly when l % 16 == 0 || l % 16 >= 9. Pairs already routed
 * to u8 under the historical `l_ms*a < 250` bound keep their routing regardless
 * (they were u8 before and stay u8), so the geometry clause only gates the new
 * admissions the tightened bound opens up.
 * ------------------------------------------------------------------------- */

/* Kernel score bias for a 5x5 scoring matrix, matching ksw_qinit's
 * `256 - (uint8_t)min` reduction. >= 1 for any real matrix (min <= 0). */
static inline int matesw_u8_shift(const int8_t *mat) {
    int mn = 0;
    for (int i = 0; i < 25; ++i) if (mat[i] < mn) mn = mat[i];
    int s = -mn;
    return s < 1 ? 1 : s;
}

/* Worst-case (max) bias across every matrix the SW at a site could use. Taking
 * the max only ever routes MORE pairs to the safe 16-bit tier, so it can never
 * over-admit. Under --meth the actual matrix is one of mat_ot/mat_ob (chosen by
 * the rescued mate's read#); folding both in keeps the batched pre/post split
 * consistent no matter which hypothesis a pair ends up scored under. */
static inline int matesw_u8_shift_for(const mem_opt_t *opt, const int8_t *base_mat) {
    int s = matesw_u8_shift(base_mat);
    if (opt->meth_mode) {
        int so = matesw_u8_shift(opt->mat_ot);
        int sb = matesw_u8_shift(opt->mat_ob);
        if (so > s) s = so;
        if (sb > s) s = sb;
    }
    return s;
}

/* True iff a mate of query length l_ms may be scored on the u8 tier
 * byte-identically. `a` is the match reward (opt->a); `shift` the kernel bias
 * from matesw_u8_shift_for(). */
static inline int matesw_use_u8(int l_ms, int a, int shift) {
    int prod = l_ms * a;                              /* max attainable SW score */
    return (prod + shift <= 254)                      /* no u8 saturation / early break */
        && (prod < 250                                /* historical admissions keep routing */
            || l_ms % 16 == 0 || (l_ms % 16) >= 9);   /* new admissions: pad geometry equal */
}

#endif
