/* Ported from bwa (https://github.com/lh3/bwa), Apache2 branch, commit a08c7c0.
 * Copyright (c) Heng Li and the bwa contributors; licensed under the Apache License,
 * Version 2.0 (see LICENSE in this directory).
 *
 * MODIFIED from the original: moved into namespace bwa_aln; alignments come from each
 * read's own aln array instead of being read from a .sai file. See README.md. */
#ifndef BWA_MEM3_ALN_BWAPE_H
#define BWA_MEM3_ALN_BWAPE_H

#include "bwtaln.h"
#include "../bntseq.h"

namespace bwa_aln {

typedef struct {
	double avg, std, ap_prior;
	bwtint_t low, high, high_bayesian;
} isize_info_t;

/* sampe's defaults (caller frees with free()). */
pe_opt_t *bwa_init_pe_opt(void);

/* Estimates the insert-size distribution from the pairs where both ends are unique enough.
 * Returns 0, or -1 (ii->avg = -1) when it cannot; prints its findings to stderr. */
int infer_isize(int n_seqs, bwa_seq_t *seqs[2], isize_info_t *ii, double ap_prior, int64_t L);

/* Single-end positions, then pairing, for one batch. Each read's alignments are its aln /
 * n_aln. Consumes drand48(): call serially, in read order. When the batch gives no estimate
 * and last_ii->avg > 0, last_ii is used. Returns the number of alignments whose position
 * changed. */
int bwa_cal_pac_pos_pe(const bntseq_t *bns, const bwt_t *bwt, int n_seqs, bwa_seq_t *seqs[2],
					   isize_info_t *ii, const pe_opt_t *opt, const gap_opt_t *gopt,
					   const isize_info_t *last_ii);

/* Mate Smith-Waterman for unpaired reads. `pacseq` is the packed forward reference; it is
 * returned (stock may load it itself, here it is always given). */
ubyte_t *bwa_paired_sw(const bntseq_t *bns, const ubyte_t *pacseq, int n_seqs, bwa_seq_t *seqs[2],
					   const pe_opt_t *popt, const isize_info_t *ii);

} // namespace bwa_aln

#endif
