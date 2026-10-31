/* Ported from bwa (https://github.com/lh3/bwa), Apache2 branch, commit a08c7c0.
 * Copyright (c) Heng Li and the bwa contributors; licensed under the Apache License,
 * Version 2.0 (see LICENSE in this directory).
 *
 * MODIFIED from the original: moved into namespace bwa_aln; alignments come from each
 * read's own aln array instead of being read from a .sai file; bwa_cal_pac_pos_pe is cut into
 * four stages and bwa_paired_sw into a counting core and a report; bwa_pe_pairing and
 * bwa_pe_pair_grouped are new (for the unit tests). See README.md. */
#ifndef BWA_MEM3_ALN_BWAPE_H
#define BWA_MEM3_ALN_BWAPE_H

#include "bwtaln.h"
#include "../bntseq.h"
#include "../utils.h"

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
/* The driver (aln_run_pe) runs the four stages below itself, so it can run stages 2 and 4 in
 * parallel; this composition is kept as upstream's function for the unit tests, which check the
 * stages against it. aln_run_pe's own insert-size carry-over is tested separately. */
int bwa_cal_pac_pos_pe(const bntseq_t *bns, const bwt_t *bwt, int n_seqs, bwa_seq_t *seqs[2],
					   isize_info_t *ii, const pe_opt_t *opt, const gap_opt_t *gopt,
					   const isize_info_t *last_ii);

/* bwa_cal_pac_pos_pe's stages, in its order; it calls them one after another. Pass slices
 * (seqs[0] + i, seqs[1] + i) to run a stage on part of a batch.
 *  1. bwa_pe_aln2seq: flags and best hit of each read. Consumes drand48(): call serially over
 *     the whole batch, in read order (pair by pair, read 1 then read 2), before stage 2.
 *  2. bwa_pe_se_pos: single-end mapQ and position. Independent per pair; no RNG.
 *  3. bwa_pe_isize: infer_isize over the whole batch, with the last_ii fallback and -A.
 *  4. bwa_pe_pair: pairing and the extra hits. Independent per pair given ii; no RNG (its
 *     bwa_aln2seq_core calls use set_main = 0 and never reach the sampling branch). Returns the
 *     number of alignments whose position changed. */
void bwa_pe_aln2seq(int n_seqs, bwa_seq_t *seqs[2]);
void bwa_pe_se_pos(const bntseq_t *bns, const bwt_t *bwt, int n_seqs, bwa_seq_t *seqs[2], const gap_opt_t *gopt);
void bwa_pe_isize(const bwt_t *bwt, int n_seqs, bwa_seq_t *seqs[2], isize_info_t *ii, const pe_opt_t *opt,
				  const isize_info_t *last_ii);
int bwa_pe_pair(const bntseq_t *bns, const bwt_t *bwt, int n_seqs, bwa_seq_t *seqs[2],
				const isize_info_t *ii, const pe_opt_t *opt, const gap_opt_t *gopt);
/* sampe's pairing of one pair, for the unit tests' reference: `arr` holds upstream's
 * (position, k<<2 | strand<<1 | end) array for the pair (sorted in place); the alignments are
 * p[0]->aln and p[1]->aln. Returns the number of alignments whose position changed. */
int bwa_pe_pairing(bwa_seq_t *p[2], pair64_v *arr, const pe_opt_t *opt, int s_mm, const isize_info_t *ii);
/* The largest weight of a group of pairs in bwa_pe_pair (see bwa_pe_pair_grouped). */
constexpr long long ALN_PE_GROUP_ROWS = 8192;
/* bwa_pe_pair with the group size as a parameter (bwa_pe_pair passes ALN_PE_GROUP_ROWS); the unit
 * tests call it with other sizes. It makes its SA lookups in batches over groups of consecutive
 * pairs of weight up to group_rows: pairing rows, plus a bound on extra hits, plus one per pair;
 * a heavier pair is a group alone. The output does not depend on group_rows. */
int bwa_pe_pair_grouped(const bntseq_t *bns, const bwt_t *bwt, int n_seqs, bwa_seq_t *seqs[2],
						const isize_info_t *ii, const pe_opt_t *opt, const gap_opt_t *gopt, long long group_rows);

/* Mate Smith-Waterman for unpaired reads: upstream's bwa_paired_sw without its stderr report
 * (the driver calls these two directly; upstream's wrapper is not kept). Does nothing unless popt->is_sw and ii->avg >= 0;
 * otherwise adds its singleton (index 1) and discordant (index 0) counts to n_tot and
 * n_mapped. Independent per pair; no RNG. bwa_paired_sw_report prints the counts. */
void bwa_paired_sw_core(const bntseq_t *bns, const ubyte_t *pacseq, int n_seqs, bwa_seq_t *seqs[2], const pe_opt_t *popt,
						const isize_info_t *ii, uint64_t n_tot[2], uint64_t n_mapped[2]);
void bwa_paired_sw_report(const uint64_t n_tot[2], const uint64_t n_mapped[2]);

} // namespace bwa_aln

#endif
