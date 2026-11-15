/* Ported from bwa (https://github.com/lh3/bwa), Apache2 branch, commit a08c7c0.
 * Copyright (c) Heng Li and the bwa contributors; licensed under the Apache License,
 * Version 2.0 (see LICENSE in this directory).
 *
 * MODIFIED from the original: the code is moved into namespace bwa_aln, the
 * extern "C" block is dropped, the bwa_aln2seq declaration is left out (it is ported
 * with the SAM code), and the declarations of the search entry points from upstream's
 * bwtaln.c (gap_init_opt, bwa_cal_maxdiff, bwt_cal_width) are added together with
 * aln_search_batch, aln_batch_max_gapo, the gap_job_t interleaved search and its chunked
 * stack (gap_cstack_t, gap_pool_t), the width walks, gap_slots_run and the search threads' read
 * counter (read_claim_t, ALN_CLAIM_CHUNK), the thread-creation hook (aln_thread_create) and the
 * duplicate-read handling (aln_dup_state_t, ALN_DUP_*, detail::find_dup_reads and the
 * detail::dup_sample_* checks), which are new. */
#ifndef BWA_MEM3_ALN_BWTGAP_H
#define BWA_MEM3_ALN_BWTGAP_H

#include "aln_bwt.h"
#include "bwtaln.h"

#include <pthread.h>

#include <atomic>
#include <cstdint>

namespace bwa_aln {

// info packs score into its upper 11 bits. Only actual scores must fit: a large
// optimistic bound alone must not reject otherwise representable searches.
static const int GAP_SCORE_LIMIT = 1 << 11;

// The stock search without NONSTOP stops at an entry whose score exceeds best_score + s_mm
// (upstream's test, made in 64 bits). Entry scores are packed in 11 bits (0..GAP_SCORE_LIMIT-1), so
// the cut is clamped to [-1, GAP_SCORE_LIMIT-1] and stored as an int: the comparison with any
// packable score is the same. The --meth search and NONSTOP do not apply this cut (their score_cut
// stays GAP_SCORE_LIMIT - 1; --meth rescores its hits instead, aln_meth.cpp).
static inline int gap_score_cut(int64_t best_score, int s_mm)
{
	const int64_t cut = best_score + s_mm;
	return cut < -1 ? -1 : cut > GAP_SCORE_LIMIT - 1 ? GAP_SCORE_LIMIT - 1 : (int)cut;
}

typedef struct { // recursion stack
	uint32_t info; // score<<21 | i
	// n_seed_mm is never read: gap_cpush stores 0 and gap_push leaves the slot's old bits. Do not
	// start reading it without setting it in both.
	uint32_t n_mm:8, n_gapo:8, n_gape:8, state:2, n_seed_mm:6;
	uint32_t n_ins:16, n_del:16;
	int last_diff_pos;
	bwtint_t k, l; // (k,l) is the SA region of [i,n-1]
} gap_entry_t;

/* The interleaved search's policies (bwtgap.cpp): AlnStock is the stock search; AlnMeth is
 * `aln --meth` (aln_meth.h). The job functions, the stack helpers, gap_slots_run and the batch
 * driver's group loop are templates over the policy; AlnStock's instantiation is the stock search
 * exactly (test/tools/disasm_identity.sh, a manual pre-release check, compares its instructions
 * with the build before the templates). AlnMeth's per-slot state lives in its ext_t, one per slot
 * beside the slots, never in the stock structs, so their layouts do not change. */
struct aln_meth_slot_t; // aln_meth.h
struct AlnStock {
	typedef gap_entry_t entry_t;
	typedef void ext_t;
	static constexpr bool meth = false;
};
struct AlnMeth {
	typedef gap_entry_t entry_t;
	typedef aln_meth_slot_t ext_t;
	static constexpr bool meth = true;
};

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
/* A chunk of another entry type: `prev` first, as in gap_chunk_t, and a
 * capacity derived from what the pool allocates, sizeof(gap_chunk_t), so the pool and
 * gap_release_cstack serve both and a chunk of either type fits the same allocation. */
template <class E> constexpr int gap_chunk_entries() { return (int)((sizeof(gap_chunk_t) - sizeof(void*)) / sizeof(E)); }
template <class E> struct gap_chunk_of {
	gap_chunk_t *prev;
	E e[gap_chunk_entries<E>()];
};
static_assert(gap_chunk_entries<gap_entry_t>() == GAP_CHUNK_ENTRIES, "the stock chunk's capacity");
static_assert(sizeof(gap_chunk_of<gap_entry_t>) == sizeof(gap_chunk_t), "the stock chunk is gap_chunk_t");

typedef struct {
	gap_chunk_t *free_list;
	int n_free, max_free;
	// chunks allocated through this pool and not yet freed, in use or free: leak accounting,
	// checked by gap_destroy_pool and read by the unit tests
	long n_alloc;
} gap_pool_t;

/* Free stack chunks a thread's pool keeps for reuse (sizeof(gap_chunk_t) each, about 4 MiB); chunks
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
enum { JOB_OCC4 = 1, JOB_EXACT = 2, JOB_FLUSH = 3 }; // gap_job_t::pending: the lookup the next call resumes at
// (JOB_FLUSH: AlnMeth only, the job waits for its hits to be resolved and rescored)
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
	int max_diff, score_cut; // score_cut: the score above which the search stops (clamped to the packed range)
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
	GAP_SLOT_PARKED = 3, // AlnMeth: the search waits for a flush (gap_slots_run_meth skips it)
};
typedef struct {
	gap_job_t job;
	width_walk_t ww, sw;
	int state, widths_done, seed_done;
#ifdef ALN_DISASM_POSITIVE_CONTROL
	int64_t positive_control_pad[2]; // test/tools/disasm_identity.sh's positive control: must fail the gate
#endif
} gap_slot_core_t;
/* Steps the n slots (an array, so a step reads no pointer to reach its slot) in turn, starting with
 * the one after *pos (one width_walk_run of each unfinished walk, or one job step), until a slot's
 * widths are done (it is left in GAP_SLOT_WIDTHS for the caller to start its search) or its search
 * is done (GAP_SLOT_SEARCH, with job.aln/n_aln set). Sets *pos to that slot and returns it. Some
 * slot must be non-empty. */
int gap_slots_run(gap_slot_core_t *slots, int n, int *pos, const bwt_t *bwt);

/* The AlnMeth instantiations (aln_meth.h): gap_slots_run with the slots' --meth state ext[0..n); it
 * also returns when a slot parks and the group should be flushed (half the slots, or every active
 * one, parked), with that slot's job pending JOB_FLUSH and its state GAP_SLOT_PARKED. gap_job_start
 * for a --meth job. Meth candidate discovery applies no rank-dependent width shadow. */
int gap_slots_run_meth(gap_slot_core_t *slots, int n, int *pos, const bwt_t *bwt, aln_meth_slot_t *ext);
int gap_job_start_meth(gap_job_t *j, const bwt_t *bwt);

/* max_gapo as the search applies it to this batch: opt->max_gapo capped at the number of
 * differences allowed for the longest read of the whole batch (bwa_cal_maxdiff of that
 * length when opt->fnr > 0, else opt->max_diff). Computed once per batch, never per read. */
int aln_batch_max_gapo(int n_seqs, const bwa_seq_t *seqs, const gap_opt_t *opt);

/* The duplicate check's state for one input stream (one per run for single-end, one per mate for
 * paired-end), owned by the caller: after a batch with fewer than one duplicate in
 * ALN_DUP_MIN_RATIO reads, the stream's next ALN_DUP_RECHECK - 1 batches skip the check. A
 * stream's batches are alike, which two mates' batches need not be. */
struct aln_dup_state_t {
	int skip_left = 0; // batches left before the next duplicate check
};

/* Searches every read of the batch on n_threads threads (each claims chunks of consecutive reads
 * from a shared counter, read_claim_t; the hits do not depend on the thread count). When `dup` is
 * not null and the batch has duplicates (detail::dup_sample_check), a read identical to an
 * earlier read of the batch (detail::find_dup_reads) is not searched: it gets its own copy of
 * that read's aln, which is what its search would return. With `dup` null every read is searched.
 * Returns the number of reads given a copy instead of a search (0 when every read was searched).
 * On return each seqs[i] has aln/n_aln set and sa = 0, type = BWA_TYPE_NO_MATCH,
 * c1 = c2 = 0. name, seq, rseq and qual are left intact (stock bwa frees them and
 * complements seq in place; here the search complements a per-thread copy). */
int aln_search_batch(const bwt_t *bwt, int n_seqs, bwa_seq_t *seqs, const gap_opt_t *opt, int n_threads,
					 aln_dup_state_t *dup = nullptr);

/* A batch with fewer than one duplicate in ALN_DUP_MIN_RATIO reads is searched read by read, and
 * the stream's next ALN_DUP_RECHECK - 1 batches skip the check (aln_dup_state_t). The check
 * samples about 1 read in 2^ALN_DUP_SAMPLE_SHIFT; a batch under ALN_DUP_SAMPLE_MIN reads gets the
 * full pass instead. */
constexpr int ALN_DUP_MIN_RATIO = 200;
constexpr int ALN_DUP_RECHECK = 8;
constexpr int ALN_DUP_SAMPLE_SHIFT = 5;
constexpr int ALN_DUP_SAMPLE_MIN = 8192;

/* Internal to aln_search_batch; declared here for the unit tests. */
namespace detail {

/* Sets rep[i] (rep has n_seqs entries) to the index of the first read in the batch with read i's
 * length and bases, so rep[i] == i for the first copy. Returns the number of reads with
 * rep[i] != i. Duplicates are found within one batch only. The fingerprints are computed on up to
 * n_threads threads (kt_for_nt), the table pass on the calling thread; the result does not depend
 * on n_threads. */
int find_dup_reads(int n_seqs, const bwa_seq_t *seqs, int *rep, int n_threads);

/* Whether read p is in the key sample (a key of its length and three 8-base words, so every copy
 * of a read is in or out together; about 1 read in 2^ALN_DUP_SAMPLE_SHIFT). */
bool dup_key_sampled(const bwa_seq_t *p);

/* The key sample's verdict: 1 if the reads in the key sample have at least one duplicate per
 * ALN_DUP_MIN_RATIO sampled reads, or if that sample is far from its expected size (few distinct
 * keys, so it cannot speak for the batch); else 0. It catches duplicates of any multiplicity
 * whose key is sampled. */
int dup_sample_by_key(int n_seqs, const bwa_seq_t *seqs);

/* The position sample's verdict: 1 if every 2^ALN_DUP_SAMPLE_SHIFT-th read (or a wider stride, to
 * cap the sample) has at least one duplicate per ALN_DUP_MIN_RATIO sampled reads. It catches a
 * read with many copies whatever its key, which the key sample misses when that key is not
 * sampled. */
int dup_sample_by_position(int n_seqs, const bwa_seq_t *seqs);

/* Whether a batch is worth the full duplicate pass: 1 if it is smaller than ALN_DUP_SAMPLE_MIN, or
 * if either sample (dup_sample_by_key, dup_sample_by_position) says so. */
int dup_sample_check(int n_seqs, const bwa_seq_t *seqs);

} // namespace detail

/* Reads a search thread claims at a time from read_claim_t. A read's hits do not depend on which
 * thread searches it; the chunk size only trades the counter's traffic against the imbalance at
 * the end of the batch. */
constexpr int ALN_CLAIM_CHUNK = 64;

/* Reads bwa_cal_sa_reg_gap_group (bwtaln.cpp) searches at once per thread. Only the scheduling
 * depends on it; each read's hits do not.
 * Measured at -t 16 (and -t 8 where noted) on hg38 (1M aDNA reads at default options and at
 * -l 1024 -n 0.04 -o 2, 1M merged reads, 2M 20-40 bp reads, 2M 2x36 bp pairs, 1.09M real pairs)
 * with the chunked stacks, built with clang 19:
 * - Apple M3 Ultra: 12 was 4-11% faster in wall time than 8; 16 was faster than 12 only on the
 *   20-40 bp reads and the pairs (1.5-3.6%) and slower on the other three.
 * - Intel Sapphire Rapids (16 vCPU, 8 cores with SMT, AVX-512 and AVX2 tiers): 8 was 2.5-3.5%
 *   faster than 12 on four of the six workloads and within 2.1% on the other two.
 * - AWS Graviton4 and Graviton3 (Linux arm64, 16 cores, no SMT): 12 was faster in wall time and
 *   CPU time than 8 on all six workloads at -t 16 and -t 8 (Graviton4 wall 5-14%, Graviton3
 *   0.4-9%). Peak memory at -t 16 rose by up to about 0.1 GiB over 8 (on the -l 1024 run: more
 *   stack entries in flight), and grows with the thread count.
 * So every aarch64 target uses 12. Other arm64 cores (Neoverse N1 / Graviton2, Ampere, Cortex-A,
 * Grace) were not measured and inherit 12 from these three. x86-64, where 12 measured slower, and
 * any other target keep 8, the depth before the chunked stacks. Exported for the unit tests, which
 * pin the value per target. */
#if defined(__aarch64__)
constexpr int ALN_GROUP = 12;
#else
constexpr int ALN_GROUP = 8;
#endif

/* The next unclaimed read of a batch of n_seqs reads, shared by its search threads (new code;
 * bwtaln.cpp). claim() returns the first read of the next ALN_CLAIM_CHUNK reads (fewer at the
 * end) and sets *end one past the last, or returns -1 once [0, n_seqs) is used up. The chunks of
 * concurrent claims are disjoint and cover [0, n_seqs). Aligned to its own cache line (128 bytes
 * covers Apple arm64): the threads write it, and on the stack (search_reads' frame in bwtaln.cpp,
 * next to its caller's batch options) it would otherwise share a line with the options every
 * search step reads. Exported for the unit tests. */
struct alignas(128) read_claim_t {
	// `first` (default 0) is where the claims start: a unit-test hook to reach the end of a batch
	// of up to INT_MAX reads without claiming every chunk before it.
	explicit read_claim_t(int n_seqs, int first = 0) : next_(first), n_seqs_(n_seqs) {}
	int claim(int *end);
private:
	// 64-bit so the claims past the end (one per thread) cannot wrap it for any int n_seqs
	std::atomic<int64_t> next_;
	const int n_seqs_;
};

/* aln_search_batch creates its search threads through this (pthread_create by default). A unit-
 * test hook: the tests swap in a function that fails, to cover the caller's fallback worker. */
extern int (*aln_thread_create)(pthread_t *, const pthread_attr_t *, void *(*)(void *), void *);

} // namespace bwa_aln

#endif
