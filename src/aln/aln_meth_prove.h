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

} // namespace bwa_aln

#endif
