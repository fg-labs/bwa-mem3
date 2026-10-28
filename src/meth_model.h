/* SPDX-License-Identifier: MIT */
/*
 * The methylation substitution model for --meth-context cpg.
 *
 * A 7-code reference alphabet: the usual A C G T N (0..4), plus C5, a C whose
 * next base in the window's frame is G, and G6, a G whose previous base in
 * the frame is C. Read codes stay 0..4. Under cpg a C->T (OT) or G->A (OB)
 * conversion is free only against C5 (OT) or G6 (OB); against a plain C or G
 * it is an ordinary mismatch.
 *
 * The frame is the window's fetch frame: a reverse-strand window (doubled
 * coordinate >= l_pac) comes out of bns_get_seq_into already
 * reverse-complemented, and is labelled in that orientation. CpG is its own
 * reverse complement, so a CpG gets the same two labels on either strand.
 */
#ifndef BWA_MEM3_METH_MODEL_H
#define BWA_MEM3_METH_MODEL_H

#include <stdint.h>

#include "bntseq.h"
#include "bwamem.h"
#include "meth_opts.h"
#include "utils.h"

enum { METH_A = 0, METH_C = 1, METH_G = 2, METH_T = 3, METH_N = 4,
       METH_C5 = 5, METH_G6 = 6, METH_M = 7 };

/* The plain base of a reference code: C5 -> C, G6 -> G, others unchanged. */
static inline uint8_t meth_code_base(uint8_t c) { return c >= METH_C5 ? c - 4 : c; }

/* Build o->mat7_ot / o->mat7_ob (row = reference code, column = read code,
 * stride 7) from o->mat and o->meth_scoring. Rows and columns 0..4 copy
 * o->mat; rows 5/6 and columns 5/6 copy those of C/G. Only C5 x read T (OT)
 * and G6 x read A (OB) hold the conversion value (0 under NEUTRAL, +a
 * otherwise); a plain C/T or G/A stays a mismatch. COLLAPSED is not combined
 * with cpg, so its mirror cell is not modelled. Called by
 * mem_opt_fill_meth_mat. */
void meth_model_fill_mats(mem_opt_t *o);

/* Relabel win[0..len) in place: C -> C5 where the next frame base is G, and
 * G -> G6 where the previous frame base is C. win holds 0..3 codes in fetch
 * frame starting at doubled coordinate rb, as bns_get_seq_into produces it.
 * The neighbours beyond either end of the window are read from pac in the
 * same frame. A neighbour outside contig rid, or a position or neighbour
 * inside a bns->ambs run (an original-reference N), never forms a CpG.
 * Window positions outside contig rid are not labelled.
 *
 * Caller contract:
 *   - win is a writable, thread-private copy, never a pointer into a shared
 *     reference image (bns_fetch_seq_v2_rid returns ref_string + beg when one
 *     is loaded);
 *   - bns/pac are the ORIGINAL-reference handles (meth_orig_bns /
 *     meth_orig_pac), not the converted seed index;
 *   - rb is the start the fetch actually used (its post-clamp *beg);
 *   - win holds codes 0..3, or labels from an earlier call; code 4 is not
 *     expected. */
void meth_label_cpg(uint8_t *win, int64_t len, int64_t rb, int rid,
                    const bntseq_t *bns, const uint8_t *pac);

/* Map C5 -> C and G6 -> G in place, for a consumer that needs plain codes. */
void meth_unlabel(uint8_t *win, int64_t len);

/* Whether this run scores under --meth-context cpg. */
static inline int meth_cpg_on(const mem_opt_t *opt)
{
    return opt->meth_mode && opt->meth_context == METH_CTX_CPG;
}

/* The substitution model for a (strand-adjusted) hypothesis: under cpg with a
 * hypothesis, the 7-code matrix (mat7_ot / mat7_ob) and *m = 7, so the
 * reference window must be labelled (meth_label_cpg) before it is scored;
 * otherwise mem_opt_meth_mat's 5-code matrix and *m = 5, on plain codes. */
static inline const int8_t *meth_model_mat(const mem_opt_t *opt, int hyp, int *m)
{
    if (meth_cpg_on(opt) && hyp >= 0) {
        *m = METH_M;
        return (hyp & 1) ? opt->mat7_ot : opt->mat7_ob;
    }
    *m = 5;
    return mem_opt_meth_mat(opt, hyp);
}

/* The alphabet size of a matrix meth_model_mat or mem_opt_meth_mat returned:
 * 7 for mat7_ot / mat7_ob, 5 for every other. For a consumer handed only the
 * matrix pointer (mem_seed_sw, mem_patch_reg).
 *
 * m is recovered by POINTER IDENTITY with opt's own arrays: a copy of a mat7
 * matrix (or one from another mem_opt_t) reads as m = 5 and would be indexed
 * at stride 5, silently. Pass the pointer meth_model_mat returned, unchanged. */
static inline int meth_mat_m(const mem_opt_t *opt, const int8_t *mat)
{
    return (mat == opt->mat7_ot || mat == opt->mat7_ob) ? METH_M : 5;
}

/* Label a fetched window for an m = 7 consumer: meth_label_cpg with the
 * caller contract checked (rid valid, the window not inside a loaded reference
 * image: under --meth no .0123 is loaded, so `ref_string` is NULL and every
 * window is a private decode). bns/pac must be the original-reference handles. */
static inline void meth_label_window(uint8_t *win, int64_t len, int64_t rb, int rid,
                                     const bntseq_t *bns, const uint8_t *pac,
                                     const uint8_t *ref_string)
{
    xassert(ref_string == NULL, "meth_label_window: a loaded reference image must never be labelled in place");
    xassert(bns != NULL && pac != NULL, "meth_label_window: no reference handles");
    xassert(rid >= 0 && rid < bns->n_seqs, "meth_label_window: window has no valid contig");
    meth_label_cpg(win, len, rb, rid, bns, pac);
}

#endif
