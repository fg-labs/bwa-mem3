/* SPDX-License-Identifier: MIT */
/* New code for bwa-mem3 (not from upstream bwa): exact refutations for the `aln --meth` search.
 *
 * The --meth search runs on the converted index. A search path is an alignment of the complemented
 * converted read (the job's seq, "S" below) to the index text T (the converted reference, both
 * copies, then its reverse complement), and under BWA_MODE_GAPE each of its differences (a
 * mismatch, a gap open or a gap extension) is one unit edit: a mismatch, or one inserted or
 * deleted base. So a lower bound on the unit edit distance between S and T at a locus bounds every
 * path to that locus, whatever the search's own limits (seed, gap opens, end skips), which only
 * remove paths. The functions here compute such bounds against T itself; a bound above a path's
 * remaining differences proves that no path reaches a hit there. Nothing here changes what a
 * search that is not cut short would report.
 *
 * T is read from the original reference (aln_meth_ref_t::pac) as base masks: a position holds a
 * set of bases, and a read base matches when it is in the set. The converted copy of an original C
 * (C->T copy) or G (G->A copy) may be either base, since the converted FASTA converts upper-case
 * bases only; an ambiguous original base (an .amb run, which the converted index fills with bases
 * of its own) may be any base. A read N matches nothing, as in the search. So a bound computed on
 * the masks is a bound on T as indexed. */
#ifndef BWA_MEM3_ALN_METH_PROVE_H
#define BWA_MEM3_ALN_METH_PROVE_H

#include <stdint.h>
#include <vector>

#include "aln_bwt.h"
#include "aln_meth.h"

namespace bwa_aln {

/* The base masks of T[y0, y1) (bit b: base b may be there) into out[0, y1 - y0). T has length
 * bwt seq_len = 2 * conv l_pac; positions in [conv l_pac, 2 * conv l_pac) are the reverse
 * complement of the converted forward text. Requires 0 <= y0 <= y1 <= 2 * conv l_pac. */
void aln_meth_text_masks(const aln_meth_ref_t &ref, int64_t y0, int64_t y1, uint8_t *out);

/* Banded unit edit distance of query q[0, n) (codes 0-3, or 4 for N, which matches nothing) against
 * text masks t[0, m), over the cells D[a][b] with dlo <= b - a <= dhi. Row 0 is D[0][b] = b
 * (anchored start), or 0 with free_start (the alignment may start anywhere in t); the result is
 * min over b of D[n][b] (the end in t is free), or kmax + 1 as soon as a whole row exceeds kmax
 * (row minima never decrease). Any alignment that stays in the band and costs at most kmax is
 * counted, so the result never exceeds the true minimum over such alignments. */
int aln_meth_band_ed(const ubyte_t *q, int n, const uint8_t *t, int m, int dlo, int dhi, bool free_start, int kmax,
                     std::vector<int16_t> *scratch);

/* Scratch for the window refutations (per search slot; reused). */
struct aln_meth_verify_scratch_t {
    std::vector<uint8_t> text;
    std::vector<int16_t> dp;
};

/* Whether S[0, len) has no alignment of at most kmax unit edits starting at T position
 * s in [y0, y0 + 2 * kmax] (the text after it is free): the semi-global distance on the band of
 * those starts. Returns 1 when none exists. */
int aln_meth_refute_window(const aln_meth_ref_t &ref, int64_t seq_len, const ubyte_t *seq, int len, int64_t y0, int kmax,
                           aln_meth_verify_scratch_t *s);
/* The same for every start in [y_lo, y_hi] (y_lo <= y_hi): one banded pass over their union. */
int aln_meth_refute_starts(const aln_meth_ref_t &ref, int64_t seq_len, const ubyte_t *seq, int len, int64_t y_lo,
                           int64_t y_hi, int kmax, aln_meth_verify_scratch_t *s);

/* The pigeonhole proof that a read has no hit (no search path within kmax differences reaches
 * i == 0). S[0, len) is cut into p consecutive pieces, piece j with an edit allowance e_j of 0 or 1,
 * such that sum(e_j + 1) > kmax. Any alignment of S with at most kmax unit edits then aligns some
 * piece within its allowance (otherwise its edits would number at least sum(e_j + 1)). Each piece's
 * occurrences within its allowance are enumerated completely by a backtracking search on the index
 * (mismatches, inserted and deleted bases, pruned only by the piece's own absent-substring bound),
 * every occurrence row gives the window of starts it allows for S (aln_meth_refute_window), and
 * the read has no hit when every window is refuted.
 *
 * width: the search's widths of S (bwt_cal_width order, so width[t].bid bounds S[0, t]), or null.
 * Returns 1 when the proof is complete, 0 when it is not (a window that may hold an alignment, more
 * than ALN_METH_PH_ROWS candidate rows, more than ALN_METH_PH_NODES search nodes, or a read too
 * short for its pieces): the caller then searches as usual. *steps gets the Occ lookups made,
 * *rows the candidate rows resolved. */
#ifndef ALN_METH_PH_ROWS
#define ALN_METH_PH_ROWS 4096
#endif
#ifndef ALN_METH_PH_NODES
#define ALN_METH_PH_NODES (1 << 18)
#endif
enum { ALN_METH_PH_MIN_PIECE = 12 }; // shorter pieces have too many occurrences to be worth it
struct aln_meth_ph_node_t { bwtint_t k, l; int32_t j; uint8_t piece, used; };
struct aln_meth_ph_scratch_t {
    std::vector<aln_meth_ph_node_t> cur, next;
    std::vector<int> bid;    // per piece: bid[off + t] bounds S[a, a + t]
    std::vector<bwtint_t> occ_k, occ_l;
    std::vector<int> occ_a;
    std::vector<bwtint_t> rows;
    std::vector<int64_t> starts;
    aln_meth_verify_scratch_t vs;
    int why = 0; // the last call's outcome: ALN_METH_PH_*
};
enum { ALN_METH_PH_PROVEN = 0, ALN_METH_PH_SHORT = 1, ALN_METH_PH_ROWS_CAP = 2, ALN_METH_PH_NODES_CAP = 3, ALN_METH_PH_WINDOW = 4,
       ALN_METH_PH_NOT_TEXT = 5 };
int aln_meth_pigeonhole(const bwt_t *bwt, const aln_meth_ref_t &ref, const ubyte_t *seq, int len, int kmax,
                        const bwt_width_t *width, aln_meth_ph_scratch_t *s, long *steps, long *rows,
                        int64_t max_rows = ALN_METH_PH_ROWS, long max_nodes = ALN_METH_PH_NODES);

} // namespace bwa_aln

#endif
