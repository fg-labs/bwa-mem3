/* Ported from bwa (https://github.com/lh3/bwa), Apache2 branch, commit a08c7c0.
 * Copyright (c) Heng Li and the bwa contributors; licensed under the Apache License,
 * Version 2.0 (see LICENSE in this directory).
 *
 * MODIFIED from the original: the code is moved into namespace bwa_aln, the
 * extern "C" block is dropped, the bwa_aln2seq declaration is left out (it is ported
 * with the SAM code), and the declarations of the search entry points from upstream's
 * bwtaln.c (gap_init_opt, bwa_cal_maxdiff, bwt_cal_width) are added together with
 * aln_search_batch, aln_batch_max_gapo and the gap_job_t interleaved search, which are new. */
#ifndef BWA_MEM3_ALN_BWTGAP_H
#define BWA_MEM3_ALN_BWTGAP_H

#include "aln_bwt.h"
#include "bwtaln.h"

namespace bwa_aln {

typedef struct { // recursion stack
	uint32_t info; // score<<21 | i
	uint32_t n_mm:8, n_gapo:8, n_gape:8, state:2, n_seed_mm:6;
	uint32_t n_ins:16, n_del:16;
	int last_diff_pos;
	bwtint_t k, l; // (k,l) is the SA region of [i,n-1]
} gap_entry_t;

typedef struct {
	int n_entries, m_entries;
	gap_entry_t *stack;
} gap_stack1_t;

typedef struct {
	int n_stacks, best, n_entries;
	gap_stack1_t *stacks;
} gap_stack_t;

gap_stack_t *gap_init_stack2(int max_score);
gap_stack_t *gap_init_stack(int max_mm, int max_gapo, int max_gape, const gap_opt_t *opt);
void gap_destroy_stack(gap_stack_t *stack);
bwt_aln1_t *bwt_match_gap(bwt_t *const bwt, int len, const ubyte_t *seq, bwt_width_t *w,
						  bwt_width_t *seed_w, const gap_opt_t *opt, int *_n_aln, gap_stack_t *stack);

/* bwt_match_gap for one read as a resumable job (new code; see bwtgap.cpp). The caller sets
 * the inputs (len, seq, width, seed_width, opt, stack: what bwt_match_gap takes), then calls
 * gap_job_start, which returns 1 if the read is already done; otherwise gap_job_run until it
 * returns 1. Each gap_job_run call that returns 0 has prefetched the Occ blocks its next call
 * reads first: the next bwt_2occ4, or the next base of an exact-match tail (bwt_match_exact_alt,
 * which the job runs one base per call). When done, aln/n_aln hold what bwt_match_gap returns
 * (the caller owns aln). */
typedef struct {
	// inputs
	int len;
	const ubyte_t *seq;      // the complemented read, as bwt_match_gap takes it
	bwt_width_t *width, *seed_width;
	gap_opt_t opt;
	gap_stack_t *stack;
	// results
	bwt_aln1_t *aln;
	int n_aln;
	// state
	int best_score, best_diff, max_diff, best_cnt, max_entries, m_aln;
	int pending, i, m, m_seed; // pending: 0, or which lookup the next call resumes at (bwtgap.cpp)
	bwtint_t k, l;
	gap_entry_t e;
} gap_job_t;
int gap_job_start(gap_job_t *j, const bwt_t *bwt);
int gap_job_run(gap_job_t *j, const bwt_t *bwt);

/* bwa's stock option defaults (caller frees with free()). */
gap_opt_t *gap_init_opt(void);
/* Largest number of differences worth searching for an l-bp read at error rate err
 * so that a read with more is missed with probability below thres. */
int bwa_cal_maxdiff(int l, double err, double thres);
/* `width` must hold len+1 zeroed entries. */
int bwt_cal_width(const bwt_t *bwt, int len, const ubyte_t *str, bwt_width_t *width);

/* bwt_cal_width as a resumable walk (new code; see bwtaln.cpp), so that a read's widths can be
 * computed interleaved with other reads' work. width_walk_init, then width_walk_run until it
 * returns 1; width[0..len] then holds what bwt_cal_width writes (and `bid` its return value).
 * Exposed for the unit tests. */
typedef struct {
	const ubyte_t *str;
	bwt_width_t *width;
	int len, i, bid, pending;
	bwtint_t k, l;
} width_walk_t;
void width_walk_init(width_walk_t *w, const bwt_t *bwt, int len, const ubyte_t *str, bwt_width_t *width);
int width_walk_run(width_walk_t *w, const bwt_t *bwt);

/* max_gapo as the search applies it to this batch: opt->max_gapo capped at the number of
 * differences allowed for the longest read of the whole batch (bwa_cal_maxdiff of that
 * length when opt->fnr > 0, else opt->max_diff). Computed once per batch, never per read. */
int aln_batch_max_gapo(int n_seqs, const bwa_seq_t *seqs, const gap_opt_t *opt);

/* Searches every read of the batch on n_threads threads (reads are dealt round-robin).
 * On return each seqs[i] has aln/n_aln set and sa = 0, type = BWA_TYPE_NO_MATCH,
 * c1 = c2 = 0. name, seq, rseq and qual are left intact (stock bwa frees them and
 * complements seq in place; here the search complements a per-thread copy). */
void aln_search_batch(const bwt_t *bwt, int n_seqs, bwa_seq_t *seqs, const gap_opt_t *opt, int n_threads);

} // namespace bwa_aln

#endif
