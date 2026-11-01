/* Ported from bwa (https://github.com/lh3/bwa), Apache2 branch, commit a08c7c0.
 * Copyright (c) Heng Li and the bwa contributors; licensed under the Apache License,
 * Version 2.0 (see LICENSE in this directory).
 *
 * MODIFIED from the original: the code is moved into namespace bwa_aln, the
 * extern "C" block is dropped, the bwa_aln2seq declaration is left out (it is ported
 * with the SAM code), and the declarations of the search entry points from upstream's
 * bwtaln.c (gap_init_opt, bwa_cal_maxdiff, bwt_cal_width) are added together with
 * aln_search_batch, aln_batch_max_gapo, the gap_job_t interleaved search and its chunked
 * stack (gap_cstack_t, gap_pool_t), the width walks and gap_slots_run, which are new. */
#ifndef BWA_MEM3_ALN_BWTGAP_H
#define BWA_MEM3_ALN_BWTGAP_H

#include "aln_bwt.h"
#include "bwtaln.h"

namespace bwa_aln {

typedef struct { // recursion stack
	uint32_t info; // score<<21 | i
	// n_seed_mm is never read: gap_cpush stores 0 and gap_push leaves the slot's old bits. Do not
	// start reading it without setting it in both.
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
/* Upstream's search for one read, unchanged. The shipped search runs gap_job_* (below), a copy of
 * this loop as a resumable job; bwt_match_gap stays as the in-tree reference that the unit tests
 * check gap_job_run against (test_aln_search.cpp). No production code calls it. */
bwt_aln1_t *bwt_match_gap(bwt_t *const bwt, int len, const ubyte_t *seq, bwt_width_t *w,
						  bwt_width_t *seed_w, const gap_opt_t *opt, int *_n_aln, gap_stack_t *stack);

/* The stack gap_job_t uses (new code; see bwtgap.cpp). It holds the same entries as gap_stack_t
 * and pops them in exactly the same order (the lowest-score non-empty bucket, last in first
 * out within a bucket), but each bucket is a chain of fixed-size chunks taken from, and given
 * back to, a gap_pool_t, so the memory it holds follows the entries in use instead of each
 * bucket's high-water mark. A pool serves one thread's stacks (it has no locking); it keeps at
 * most max_free free chunks and frees the rest. */
constexpr int GAP_CHUNK_BYTES = 16384;
constexpr int GAP_CHUNK_ENTRIES = (int)((GAP_CHUNK_BYTES - sizeof(void*)) / sizeof(gap_entry_t));
typedef struct gap_chunk_s gap_chunk_t;
struct gap_chunk_s {
	gap_chunk_t *prev; // the chunk below this one in its bucket (or in the free list)
	gap_entry_t e[GAP_CHUNK_ENTRIES];
};
static_assert(sizeof(gap_chunk_t) <= GAP_CHUNK_BYTES, "a stack chunk must fit in GAP_CHUNK_BYTES");
static_assert(GAP_CHUNK_ENTRIES == 511, "src/aln/README.md states 511 entries per chunk");

typedef struct {
	gap_chunk_t *free_list;
	int n_free, max_free;
	// chunks allocated through this pool and not yet freed, in use or free: leak accounting,
	// checked by gap_destroy_pool and read by the unit tests
	long n_alloc;
} gap_pool_t;

/* Free stack chunks a thread's pool keeps for reuse (GAP_CHUNK_BYTES each, so 4 MiB); chunks
 * given back beyond this are freed. A read's stack gives its chunks back as its buckets drain and
 * when it ends, so the memory held follows the entries in use; this bounds what an idle thread
 * keeps. Output does not depend on it. */
constexpr int ALN_POOL_KEEP = 256;

typedef struct {
	gap_chunk_t *top; // the chunk holding the last entries pushed
	int top_n;        // entries in `top`
	int n_entries;    // entries in the bucket
} gap_bucket_t;

typedef struct {
	int n_stacks, best, n_entries;
	gap_bucket_t *stacks;
	gap_pool_t *pool;
} gap_cstack_t;

void gap_init_pool(gap_pool_t *pool, int max_free);
void gap_destroy_pool(gap_pool_t *pool); // every stack using it must have been released (checked)
gap_cstack_t *gap_init_cstack(int max_mm, int max_gapo, int max_gape, const gap_opt_t *opt, gap_pool_t *pool);
void gap_release_cstack(gap_cstack_t *stack); // empties the stack, giving its chunks back to the pool
void gap_destroy_cstack(gap_cstack_t *stack);

/* Out-of-line push, pop and reset of both stacks (gap_stack_t has a reset; gap_cstack_t's is
 * gap_release_cstack), for the unit tests that compare their pop orders. */
void gap_stack_push_ext(gap_stack_t *stack, int i, bwtint_t k, bwtint_t l, int n_mm, int n_gapo, int n_gape,
						int n_ins, int n_del, int state, int is_diff, const gap_opt_t *opt);
void gap_stack_pop_ext(gap_stack_t *stack, gap_entry_t *e);
void gap_stack_reset_ext(gap_stack_t *stack);
void gap_cstack_push_ext(gap_cstack_t *stack, int i, bwtint_t k, bwtint_t l, int n_mm, int n_gapo, int n_gape,
						 int n_ins, int n_del, int state, int is_diff, const gap_opt_t *opt);
void gap_cstack_pop_ext(gap_cstack_t *stack, gap_entry_t *e);

/* bwt_match_gap for one read as a resumable job (new code; see bwtgap.cpp). The caller sets
 * the inputs (len, seq, width, seed_width, opt: what bwt_match_gap takes, and a gap_cstack_t), then calls
 * gap_job_start, which returns 1 if the read is already done; otherwise gap_job_run until it
 * returns 1 (the driver steps jobs through gap_slots_run, which runs the same step, gap_job_step,
 * inlined; gap_job_run is that step as a call, for the unit tests). Each call that returns 0 has
 * prefetched the Occ blocks its next call reads first: the next bwt_2occ4, or the next base of an
 * exact-match tail (bwt_match_exact_alt, which the job runs one base per call). When done, aln/n_aln hold what bwt_match_gap returns
 * (the caller owns aln).
 *
 * The fields the job step reads on every call come first, so that they share as few cache lines as
 * possible: the resume state, the inputs it uses, and `c`, the options it reads, copied from opt by
 * gap_job_start (with seed_off = len - opt.seed_len precomputed).
 *
 * Invariants: `c` and `score_cut` are valid only after gap_job_start, and neither opt nor len may
 * change after it (the hit path reads opt, the step reads c). score_cut tracks best_score: every
 * write of best_score (gap_job_start, gap_job_hit's first hit) updates it too. */
enum { JOB_OCC4 = 1, JOB_EXACT = 2 }; // gap_job_t::pending: the lookup the next call resumes at
typedef struct { // the job step's copy of the options it reads (set by gap_job_start)
	int s_mm, s_gapo, s_gape, mode;
	int max_entries, max_gapo, max_gape, max_seed_diff;
	int indel_end_skip, max_del_occ, seed_off;
} gap_job_const_t;

typedef struct {
	// state: what the next step resumes from
	int pending, i, m, m_seed; // pending: 0, JOB_OCC4 or JOB_EXACT
	bwtint_t k, l;
	gap_entry_t e;
	// inputs
	gap_cstack_t *stack;     // emptied (chunks back to its pool) when the job ends
	const ubyte_t *seq;      // the complemented read, as bwt_match_gap takes it
	bwt_width_t *width, *seed_width;
	int len;
	// state read on every call
	int max_diff, score_cut; // score_cut: the score above which the search stops
	int max_entries;         // the most entries the stack has held (a statistic; the -m limit is c.max_entries)
	gap_job_const_t c;
	// inputs and state the hits use
	gap_opt_t opt;
	int best_score, best_diff, best_cnt, m_aln;
	// results
	bwt_aln1_t *aln;
	int n_aln;
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

/* bwt_cal_width as a resumable walk (new code; see bwtgap.cpp), so that a read's widths can be
 * computed interleaved with other reads' work. width_walk_init, then width_walk_run until it
 * returns 1; width[0..len] then holds what bwt_cal_width writes (and `bid` its return value).
 * Do not call width_walk_run again after it has returned 1 (it would bump `bid` and rewrite
 * width[len]); width_walk_init starts a new walk. Exposed for the unit tests. */
typedef struct {
	const ubyte_t *str;
	bwt_width_t *width;
	int len, i, bid, pending;
	bwtint_t k, l;
} width_walk_t;
void width_walk_init(width_walk_t *w, const bwt_t *bwt, int len, const ubyte_t *str, bwt_width_t *width);
int width_walk_run(width_walk_t *w, const bwt_t *bwt);

/* One read's place in the interleaved search (new code; see bwtgap.cpp): its two width walks
 * (whole read and seed), then its search job. bwa_cal_sa_reg_gap_group (bwtaln.cpp) sets them
 * up; gap_slots_run steps them. */
enum {
	GAP_SLOT_EMPTY = 0,  // no read
	GAP_SLOT_WIDTHS = 1, // the width walks run until widths_done and seed_done are both set
	GAP_SLOT_SEARCH = 2, // the search is started (set just before gap_job_start; stepped once it returns 0)
};
typedef struct {
	gap_job_t job;
	width_walk_t ww, sw;
	int state, widths_done, seed_done;
} gap_slot_core_t;
/* Steps the n slots (an array, so a step reads no pointer to reach its slot) in turn, starting with
 * the one after *pos (one width_walk_run of each unfinished walk, or one job step), until a slot's
 * widths are done (it is left in GAP_SLOT_WIDTHS for the caller to start its search) or its search
 * is done (GAP_SLOT_SEARCH, with job.aln/n_aln set). Sets *pos to that slot and returns it. Some
 * slot must be non-empty. */
int gap_slots_run(gap_slot_core_t *slots, int n, int *pos, const bwt_t *bwt);

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
