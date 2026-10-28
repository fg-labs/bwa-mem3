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
 * Window positions outside contig rid are not labelled. */
void meth_label_cpg(uint8_t *win, int64_t len, int64_t rb, int rid,
                    const bntseq_t *bns, const uint8_t *pac);

/* Map C5 -> C and G6 -> G in place, for a consumer that needs plain codes. */
void meth_unlabel(uint8_t *win, int64_t len);

#endif
