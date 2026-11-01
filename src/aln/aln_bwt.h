/* bwa's FM-index API (bwt_occ, bwt_2occ, bwt_2occ4, bwt_sa, bwt_match_exact_alt)
 * on top of bwa-mem3's FMI_search index.
 *
 * Stock bwa and bwa-mem3 build the same BWT (identical strings on every row,
 * with the `$` on the same row) and the same suffix array, except that bwa-mem3
 * keeps the `$` in its text and so counts one more row than bwa. The relation:
 *
 *   bwa seq_len   = bwa-mem3 reference_seq_len - 1  (= 2 * l_pac)
 *   bwa primary   = bwa-mem3 sentinel_index
 *   bwa L2[c]     = bwa-mem3 count[c] - 1
 *   bwa occ(k, c) = number of c in bwa-mem3 BWT rows [0, k]
 *   bwa sa(0)     = (bwtint_t)-1   (bwa-mem3 stores seq_len for the `$` suffix)
 *
 * test/aln_index_parity_test checks this row for row against a stock `bwa index`.
 *
 * bwt_match_exact_alt is from bwa's Apache2 branch (bwt.c; see README.md). This file
 * is MODIFIED from the original: the declarations are re-expressed in namespace
 * bwa_aln over bwa-mem3's index, and bwt_t is a different struct.
 *
 * The MIT License
 *
 * Copyright (c) 2008 Genome Research Ltd (GRL).
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#ifndef BWA_MEM3_ALN_BWT_H
#define BWA_MEM3_ALN_BWT_H

#include <stdint.h>

class FMI_search;

namespace bwa_aln {

typedef uint64_t bwtint_t;
typedef uint8_t ubyte_t;

struct bwt_t {
    bwtint_t primary;   /* row whose BWT character is `$` */
    bwtint_t L2[5];     /* cumulative base counts, `$` excluded */
    bwtint_t seq_len;   /* 2 * l_pac */
    FMI_search *fmi;    /* not owned; non-const because the SA accessor is non-const */
};

/* `fmi` must already be loaded (load_index) with the compressed SA. */
bwt_t *bwt_from_fmi(FMI_search *fmi);
void bwt_destroy(bwt_t *bwt); /* does not free the FMI */

/* Number of c in the `$`-removed BWT up to and including row k; 0 for k == (bwtint_t)-1. */
bwtint_t bwt_occ(const bwt_t *bwt, bwtint_t k, ubyte_t c);
void bwt_2occ(const bwt_t *bwt, bwtint_t k, bwtint_t l, ubyte_t c, bwtint_t *ok, bwtint_t *ol);
void bwt_2occ4(const bwt_t *bwt, bwtint_t k, bwtint_t l, bwtint_t cntk[4], bwtint_t cntl[4]);

/* Prefetches the Occ block that bwt_occ(bwt, p - 1, c) reads (any c). */
void bwt_prefetch_occ(const bwt_t *bwt, bwtint_t p);

/* SA[k] in bwa's row numbering (row 0 is the `$` suffix and returns (bwtint_t)-1). */
bwtint_t bwt_sa(const bwt_t *bwt, bwtint_t k);

int bwt_match_exact_alt(const bwt_t *bwt, int len, const ubyte_t *str, bwtint_t *k0, bwtint_t *l0);

} // namespace bwa_aln

#endif
