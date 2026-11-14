/* Ported from bwa (https://github.com/lh3/bwa), Apache2 branch, commit a08c7c0.
 * Copyright (c) Heng Li and the bwa contributors; licensed under the Apache License,
 * Version 2.0 (see LICENSE in this directory).
 *
 * MODIFIED from the original: this file holds only the search half of upstream's
 * bwtaln.c (gap_init_opt, bwa_cal_maxdiff, bwt_cal_width and the per-thread search
 * bwa_cal_sa_reg_gap), in namespace bwa_aln. The search no longer complements the read
 * in place or frees name/seq/rseq/qual: it complements a per-thread scratch copy, so
 * the same reads can go on to the SAM stage. The batch-wide max_gapo cap that upstream
 * computes inside every thread is computed once in aln_search_batch before the fan-out
 * and passed down. bwa_cal_sa_reg_gap's one-read-at-a-time loop is replaced by
 * bwa_cal_sa_reg_gap_group, which prepares each read the same way and interleaves several
 * reads' searches (gap_job_t in bwtgap.cpp). Its widths are computed by width_walk_run
 * (bwtgap.cpp), bwt_cal_width's loop as a resumable walk; gap_slots_run (bwtgap.cpp) steps the
 * slots; and each slot's search stack is a gap_cstack_t drawing chunks from a per-thread pool.
 * The threads claim chunks of reads from a shared counter (read_claim_t), not i % n_threads == tid.
 * Identical reads in a batch share one search (detail::find_dup_reads, the detail::dup_sample_*
 * checks and the compaction in aln_search_batch, all new). The BAM input, .sai writer and
 * command-line parser are not ported. */
#include "bwtgap.h"
#include "aln_alloc.h"
#include "../kthread.h"

#include <math.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>
#include <atomic>
#include <memory>
#include <vector>

namespace bwa_aln {

gap_opt_t *gap_init_opt(void)
{
	gap_opt_t *o;
	o = (gap_opt_t*)aln_calloc(1, sizeof(gap_opt_t));
	/* IMPORTANT: s_mm*10 should be about the average base error
	   rate. Voilating this requirement will break pairing! */
	o->s_mm = 3; o->s_gapo = 11; o->s_gape = 4;
	o->max_diff = -1; o->max_gapo = 1; o->max_gape = 6;
	o->indel_end_skip = 5; o->max_del_occ = 10; o->max_entries = 2000000;
	o->mode = BWA_MODE_GAPE | BWA_MODE_COMPREAD;
	o->seed_len = 32; o->max_seed_diff = 2;
	o->fnr = 0.04;
	o->n_threads = 1;
	o->max_top2 = 30;
	o->trim_qual = 0;
	return o;
}

int bwa_cal_maxdiff(int l, double err, double thres)
{
	double elambda = exp(-l * err);
	double sum, y = 1.0;
	int k, x = 1;
	for (k = 1, sum = elambda; k < 1000; ++k) {
		y *= l * err;
		x *= k;
		sum += elambda * y / x;
		if (1.0 - sum < thres) return k;
	}
	return 2;
}

// width must be filled as zero
int bwt_cal_width(const bwt_t *bwt, int len, const ubyte_t *str, bwt_width_t *width)
{
	bwtint_t k, l, ok, ol;
	int i, bid;
	bid = 0;
	k = 0; l = bwt->seq_len;
	for (i = 0; i < len; ++i) {
		ubyte_t c = str[i];
		if (c < 4) {
			bwt_2occ(bwt, k - 1, l, c, &ok, &ol);
			k = bwt->L2[c] + ok + 1;
			l = bwt->L2[c] + ol;
		}
		if (k > l || c > 3) { // then restart
			k = 0;
			l = bwt->seq_len;
			++bid;
		}
		width[i].w = l - k + 1;
		width[i].bid = bid;
	}
	width[len].w = 0;
	width[len].bid = ++bid;
	return bid;
}

// Upstream does this at the top of every thread's bwa_cal_sa_reg_gap; the result is
// the same in every thread, so it is computed once per batch.
static gap_opt_t batch_opt(int n_seqs, const bwa_seq_t *seqs, const gap_opt_t *opt)
{
	int i, max_len;
	gap_opt_t local_opt = *opt;
	for (i = max_len = 0; i != n_seqs; ++i)
		if (seqs[i].len > max_len) max_len = seqs[i].len;
	if (opt->fnr > 0.0) local_opt.max_diff = bwa_cal_maxdiff(max_len, BWA_AVG_ERR, opt->fnr);
	if (local_opt.max_diff < local_opt.max_gapo) local_opt.max_gapo = local_opt.max_diff;
	return local_opt;
}

int aln_batch_max_gapo(int n_seqs, const bwa_seq_t *seqs, const gap_opt_t *opt)
{
	return batch_opt(n_seqs, seqs, opt).max_gapo;
}

/* New code: upstream's bwa_cal_sa_reg_gap with several of the thread's reads searched at once
 * (gap_job_*, bwtgap.cpp). Per read, everything bwa_cal_sa_reg_gap does is done in the same
 * way: the same widths, the same per-read max_diff and seed_len, the same complemented read
 * and a stack of its own (chunks from the thread's pool); only the reads' steps are interleaved,
 * so each read gets the hits bwt_match_gap gives it. A read's two width computations (whole read
 * and seed) are run as width walks, interleaved with the other slots' work, before its search
 * starts. `opt` is batch_opt()'s result: max_diff and max_gapo are the batch-wide values. */
namespace {
struct gap_slot_t {
	gap_slot_core_t *c; // the read's width walks and search (in `cores`), which gap_slots_run steps
	std::vector<bwt_width_t> w, seed_w;
	std::vector<ubyte_t> cseq;
	bwa_seq_t *p;
};
}

int read_claim_t::claim(int *end)
{
	const int64_t st = next_.fetch_add(ALN_CLAIM_CHUNK, std::memory_order_relaxed);
	if (st >= n_seqs_) return -1;
	*end = n_seqs_ - st > ALN_CLAIM_CHUNK? (int)st + ALN_CLAIM_CHUNK : n_seqs_; // st < n_seqs_: no overflow
	return (int)st;
}

int (*aln_thread_create)(pthread_t *, const pthread_attr_t *, void *(*)(void *), void *) = pthread_create;

/* The group driver's per-policy steps. AlnStock's are the stock driver's own expressions, so its
 * instantiation is the stock driver (test/tools/disasm_identity.sh checks it); another policy (aln
 * --meth, added later) specialises this with its own. ctx_t is what the batch passes each
 * thread beyond the reads (nothing for the stock search); the thread hands it to set_ctx before
 * the group runs. scope_t is constructed for the group's lifetime, after the slots and before the
 * reads: a policy with per-group state (aln --meth) keeps it there, its steps reach it through the
 * thread, and its destructor releases what the policy holds when the group ends. So the driver's
 * signature, its statements and its lambdas' captures are the stock driver's (every hook is an
 * expression the stock driver had, or nothing). A slot is named by its core (c - cores is its
 * index). */
template <class P> struct gap_group;
template <> struct gap_group<AlnStock> {
	struct ctx_t {};
	static inline void set_ctx(const ctx_t *) {}
	struct scope_t {
		scope_t(const bwt_t *, const gap_opt_t *, gap_slot_core_t *, int) {}
	};
	// read r (p) taken into slot c (--meth converts it here)
	static inline void slot_start(gap_slot_core_t *, bwa_seq_t *, int) {}
	// the bases slot c's width walks and search run on, read at each use as the stock driver does
	static inline const ubyte_t *seq(const gap_slot_core_t *, const bwa_seq_t *p) { return p->seq; }
	static inline int job_start(gap_job_t *j, const bwt_t *bwt) { return gap_job_start(j, bwt); }
	// steps the slots to the next one whose widths or search are done (gap_slots_run); a policy
	// that parks slots (aln --meth) also flushes them in here, so the driver never sees a parked one
	static inline int run(gap_slot_core_t *cores, int n, int *pos, const bwt_t *bwt)
	{
		return gap_slots_run(cores, n, pos, bwt);
	}
	// slot c's read is done and its hits were handed to the read (p->aln): a policy that keeps its
	// results elsewhere (aln --meth) takes them back here
	static inline void finish(gap_slot_core_t *, bwa_seq_t *) {}
};

template <class P>
static void bwa_cal_sa_reg_gap_group(read_claim_t *claim, const bwt_t *bwt, int n_seqs, bwa_seq_t *seqs,
									 const gap_opt_t *opt)
{
	typedef gap_group<P> H;
	const int G = ALN_GROUP;
	std::vector<gap_slot_t> slots(G);
	std::vector<gap_slot_core_t> cores(G);
	gap_pool_t pool;
	gap_init_pool(&pool, ALN_POOL_KEEP);
	typename H::scope_t scope(bwt, opt, cores.data(), G);
	for (int i = 0; i < G; ++i) {
		gap_slot_t &s = slots[i];
		s.c = &cores[i];
		s.c->job.stack = gap_init_cstack(opt->max_diff, opt->max_gapo, opt->max_gape, opt, &pool);
		s.c->state = GAP_SLOT_EMPTY;
		s.p = 0; // seed_w is sized on first use: only reads longer than -l use it
	}
	int next = 0, chunk_end = 0, n_active = 0; // the reads [next, chunk_end) are claimed and not yet taken
	bool claimed_all = false; // the counter is used up: no more claims (and no more traffic on it)
	// The thread's next read (claiming a new chunk when its own is used up), or -1 when none is left.
	auto next_read = [&]() -> int {
		if (next == chunk_end) {
			if (claimed_all || (next = claim->claim(&chunk_end)) < 0) {
				claimed_all = true;
				next = chunk_end = 0;
				return -1;
			}
		}
		return next++;
	};
	// Takes the thread's next read into slot s and starts its width walks; returns 0 when
	// none is left.
	auto fill = [&](gap_slot_t &s) -> int {
		const int r = next_read();
		if (r < 0) {
			s.p = 0;
			s.c->state = GAP_SLOT_EMPTY;
			return 0;
		}
		bwa_seq_t *p = seqs + r;
		p->sa = 0; p->type = BWA_TYPE_NO_MATCH; p->c1 = p->c2 = 0; p->n_aln = 0; p->aln = 0;
		H::slot_start(s.c, p, r);
		if ((int)s.w.size() < p->len + 1) s.w.resize(p->len + 1);
		width_walk_init(&s.c->ww, bwt, p->len, H::seq(s.c, p), s.w.data());
		s.c->widths_done = 0;
		s.c->seed_done = 1;
		if (p->len > opt->seed_len) {
			// sized here, not per slot up front: a large -l (seeding off) would otherwise cost
			// (-l + 1) x 16 B per slot even though no read is longer; here -l < p->len
			if (s.seed_w.size() < (size_t)opt->seed_len + 1) s.seed_w.resize((size_t)opt->seed_len + 1);
			width_walk_init(&s.c->sw, bwt, opt->seed_len, H::seq(s.c, p) + (p->len - opt->seed_len), s.seed_w.data());
			s.c->seed_done = 0;
		}
		s.c->state = GAP_SLOT_WIDTHS;
		s.p = p;
		return 1;
	};
	// Slot s's widths are done: sets up its search as bwa_cal_sa_reg_gap does and starts it.
	// Returns 1 if the read is already done (aln/n_aln set), else 0 (the search runs next).
	auto start = [&](gap_slot_t &s) -> int {
		bwa_seq_t *p = s.p;
		gap_job_t *j = &s.c->job;
		j->opt = *opt;
		if (opt->fnr > 0.0) j->opt.max_diff = bwa_cal_maxdiff(p->len, BWA_AVG_ERR, opt->fnr);
		j->opt.seed_len = opt->seed_len < p->len? opt->seed_len : 0x7fffffff;
		s.cseq.resize(p->len + 1);
		for (int k = 0; k < p->len; ++k) // we need to complement
			s.cseq[k] = H::seq(s.c, p)[k] > 3? 4 : 3 - H::seq(s.c, p)[k];
		j->len = p->len;
		j->seq = s.cseq.data();
		j->width = s.w.data();
		j->seed_width = p->len <= opt->seed_len? 0 : s.seed_w.data();
		s.c->state = GAP_SLOT_SEARCH;
		return H::job_start(j, bwt);
	};
	for (gap_slot_t &s : slots)
		if (fill(s)) ++n_active;
	int pos = G - 1;
	while (n_active) {
		gap_slot_t &s = slots[H::run(cores.data(), G, &pos, bwt)];
		if (s.c->state == GAP_SLOT_WIDTHS && !start(s)) continue; // its search runs in the next turns
		s.p->aln = s.c->job.aln; s.p->n_aln = s.c->job.n_aln; // the read is done
		H::finish(s.c, s.p);
		if (!fill(s)) --n_active;
	}
	for (gap_slot_t &s : slots) gap_destroy_cstack(s.c->job.stack);
	gap_destroy_pool(&pool);
}

/* What each search thread gets: the reads and options, and the policy's context when it has one
 * (the stock search's has none, so its threads' struct is as before). */
template <class C> struct thread_ctx_t {
	const C *ctx;
	void set_ctx(const C *c) { ctx = c; }
	const C *get_ctx() const { return ctx; }
};
template <> struct thread_ctx_t<gap_group<AlnStock>::ctx_t> {
	void set_ctx(const gap_group<AlnStock>::ctx_t *) {}
	const gap_group<AlnStock>::ctx_t *get_ctx() const { return nullptr; }
};
template <class P> struct thread_aux_t : thread_ctx_t<typename gap_group<P>::ctx_t> {
	read_claim_t *claim;
	const bwt_t *bwt;
	int n_seqs;
	bwa_seq_t *seqs;
	const gap_opt_t *opt;
};

template <class P>
static void *worker(void *data)
{
	thread_aux_t<P> *d = (thread_aux_t<P>*)data;
	gap_group<P>::set_ctx(d->get_ctx());
	bwa_cal_sa_reg_gap_group<P>(d->claim, d->bwt, d->n_seqs, d->seqs, d->opt);
	return 0;
}

/* The search keeps upstream's fan-out of fresh joinable pthreads per batch (the thread count is a
 * parameter), but the reads are not dealt round-robin as upstream does (thread t searching reads
 * t, t + n_threads, ...): each thread claims ALN_CLAIM_CHUNK consecutive reads at a time from a
 * shared counter, so a thread that drew cheap reads takes more of them and the threads finish
 * together. Each thread keeps its ALN_GROUP slots fed across chunk boundaries. Hits do not depend
 * on which thread searches a read, so the output is the same. The post-search stages use kt_for_nt;
 * the search is not moved onto it: kt_for_nt hands out item ranges of its own layout, so the slots
 * (stacks, width buffers) would have to become per-tid state and the interleave would drain at
 * every range boundary. One thread spawn per thread per batch of 0x40000 reads is small next to
 * the search. */
template <class P>
static void search_reads(const bwt_t *bwt, int n_seqs, bwa_seq_t *seqs, const gap_opt_t *bopt, int n_threads,
						 const typename gap_group<P>::ctx_t *ctx = nullptr)
{
	read_claim_t claim(n_seqs);
	if (n_threads <= 1) { // no multi-threading at all (bwa aln's -t 0 too)
		gap_group<P>::set_ctx(ctx);
		bwa_cal_sa_reg_gap_group<P>(&claim, bwt, n_seqs, seqs, bopt);
	} else {
		pthread_t *tid;
		pthread_attr_t attr;
		thread_aux_t<P> *data;
		int j;
		pthread_attr_init(&attr);
		pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
		data = (thread_aux_t<P>*)aln_calloc(n_threads, sizeof(thread_aux_t<P>));
		tid = (pthread_t*)aln_calloc(n_threads, sizeof(pthread_t));
		std::vector<char> started(n_threads, 0);
		bool any_failed = false;
		for (j = 0; j < n_threads; ++j) {
			data[j].claim = &claim; data[j].bwt = bwt;
			data[j].n_seqs = n_seqs; data[j].seqs = seqs; data[j].opt = bopt;
			data[j].set_ctx(ctx);
			started[j] = aln_thread_create(&tid[j], &attr, worker<P>, data + j) == 0;
			any_failed = any_failed || !started[j];
		}
		// A thread that could not be created must not be joined. The caller runs one worker in
		// place of the missing ones: it claims reads like the others, so it and the threads that
		// started search every read (it alone does if none started).
		if (any_failed) worker<P>(data);
		for (j = 0; j < n_threads; ++j)
			if (started[j]) pthread_join(tid[j], 0);
		free(data); free(tid);
		pthread_attr_destroy(&attr);
	}
}

/* New code: identical reads in a batch share one search. A read's hits depend only on its bases,
 * its length and the batch's options (the per-read max_diff and seed_len follow the length; the
 * batch-wide max_gapo cap is computed from every read before the search), so a read whose bases
 * and length equal an earlier read's in the same batch gets a copy of that read's hits instead of
 * a search of its own. Everything after the search (the drand48 pass included) still runs per
 * read, so the output is unchanged.
 *
 * The search itself is not changed: when a batch has duplicates, aln_search_batch moves its
 * distinct reads to the front of the batch, searches them there with search_reads and moves
 * every read back, so a batch without duplicates runs exactly the code it ran before. Finding
 * duplicates costs a pass over the batch, so samples decide first whether the pass is worth
 * running (detail::dup_sample_check). */
namespace {

static inline uint64_t load_u64(const ubyte_t *p)
{
	uint64_t x;
	memcpy(&x, p, 8);
	return x;
}

static inline uint64_t mix_u64(uint64_t h, uint64_t x)
{
	h ^= x;
	h *= 0xbf58476d1ce4e5b9ULL;
	return h ^ (h >> 31);
}

/* A fingerprint of the read's length and all its bases, read 8 at a time as 64-bit words (the
 * last, partial word zero-padded). It covers the whole read, not only its ends: reads of one
 * length that share their ends and differ in the middle (an amplicon with errors, constant flanks
 * around a variable insert) would otherwise all share a fingerprint and one probe chain, and the
 * table would compare each new one against every earlier one. A match is still confirmed with a
 * full compare. */
static inline uint64_t read_fingerprint(const bwa_seq_t *p)
{
	const int len = p->len;
	uint64_t h = mix_u64(0x9e3779b97f4a7c15ULL, (uint64_t)len);
	int i = 0;
	for (; i + 8 <= len; i += 8) h = mix_u64(h, load_u64(p->seq + i));
	if (i < len) {
		uint64_t t = 0;
		memcpy(&t, p->seq + i, (size_t)(len - i));
		h = mix_u64(h, t);
	}
	return mix_u64(h, 0);
}

/* A cheap key of the read's length and three 8-base words (start, middle, end; the whole read
 * when it is shorter than 8). Identical reads have equal keys, so sampling reads by key keeps or
 * drops every copy of a read together: the sample's share of duplicates estimates the batch's. */
static inline uint64_t read_sample_key(const bwa_seq_t *p)
{
	const int len = p->len;
	uint64_t h = mix_u64(0x2545f4914f6cdd1dULL, (uint64_t)len);
	if (len >= 8) {
		h = mix_u64(h, load_u64(p->seq));
		h = mix_u64(h, load_u64(p->seq + (len - 8) / 2));
		h = mix_u64(h, load_u64(p->seq + len - 8));
	} else if (len > 0) {
		uint64_t t = 0;
		memcpy(&t, p->seq, (size_t)len);
		h = mix_u64(h, t);
	}
	return mix_u64(h, 0);
}

struct dup_slot_t {
	uint32_t tag; // the fingerprint's low 32 bits
	int32_t idx;  // the first item with this slot's fingerprint, or -1 if the slot is empty
};

// The smallest power of two of at least 64 and at least 2 * n: a table of that many slots is at
// most half full with n items, so every probe ends.
static inline size_t dup_table_slots(int64_t n)
{
	size_t slots = 64;
	while ((int64_t)slots < 2 * n) slots <<= 1;
	return slots;
}

/* The table pass shared by the full pass and both samples: items k in [0, n), item k being read
 * seq_at(k) with fingerprint fp[k], in order. Sets rep[k] (if rep is not null) to the first item
 * with the same length and bases, and returns the number of items that are not the first copy.
 * `table` has dup_table_slots(n) slots, so it is at most half full and every probe ends. */
template <class SeqAt>
static inline int dup_table_pass(int n, SeqAt seq_at, const uint64_t *fp, int *rep, dup_slot_t *table)
{
	const size_t mask = dup_table_slots(n) - 1;
	for (size_t s = 0; s <= mask; ++s) table[s].tag = 0, table[s].idx = -1;
	int n_dups = 0;
	for (int k = 0; k < n; ++k) {
		const uint64_t h = fp[k];
		const uint32_t tag = (uint32_t)h;
		const bwa_seq_t *b = seq_at(k);
		int first = k;
		for (size_t at = (size_t)(h ^ (h >> 47)) & mask;; at = (at + 1) & mask) {
			dup_slot_t &e = table[at];
			if (e.idx < 0) { e.tag = tag, e.idx = k; break; }
			if (e.tag != tag || fp[e.idx] != h) continue;
			const bwa_seq_t *a = seq_at(e.idx);
			if (a->len == b->len && memcmp(a->seq, b->seq, b->len) == 0) {
				first = e.idx;
				++n_dups;
				break;
			}
		}
		if (rep) rep[k] = first;
	}
	return n_dups;
}

struct fp_ctx_t {
	const bwa_seq_t *seqs;
	uint64_t *fp;
};

static void fp_worker(void *data, int st, int n, int)
{
	const fp_ctx_t *c = (const fp_ctx_t*)data;
	for (int i = st; i < st + n; ++i) c->fp[i] = read_fingerprint(c->seqs + i);
}

/* The samples' buffers: fixed arrays allocated once per calling thread (under 1 MB), sized for a
 * sample of up to DUP_SAMPLE_CAP reads, 4x the expected sample of a 0x40000-read batch; a larger
 * key sample is out of range anyway (it sends the batch to the full pass), and the position
 * sample widens its stride to stay under the cap. */
constexpr int DUP_SAMPLE_CAP = 1 << 15;
struct dup_sample_buffers_t {
	uint64_t fp[DUP_SAMPLE_CAP];
	int idx[DUP_SAMPLE_CAP];
	dup_slot_t table[2 * DUP_SAMPLE_CAP];
};
static thread_local std::unique_ptr<dup_sample_buffers_t> dup_sample_buf;

static dup_sample_buffers_t &sample_buffers()
{
	if (!dup_sample_buf) dup_sample_buf.reset(new dup_sample_buffers_t);
	return *dup_sample_buf;
}

// Whether the n_s reads seqs[idx[0..n_s)] have at least one duplicate per ALN_DUP_MIN_RATIO reads.
static int sample_has_dups(const bwa_seq_t *seqs, dup_sample_buffers_t &sb, int n_s)
{
	for (int k = 0; k < n_s; ++k) sb.fp[k] = read_fingerprint(seqs + sb.idx[k]);
	auto seq_at = [&](int k) { return seqs + sb.idx[k]; };
	const int n_dups = dup_table_pass(n_s, seq_at, sb.fp, nullptr, sb.table);
	return (int64_t)n_dups * ALN_DUP_MIN_RATIO >= n_s;
}

} // namespace

namespace detail {

int find_dup_reads(int n_seqs, const bwa_seq_t *seqs, int *rep, int n_threads)
{
	if (n_seqs <= 1) {
		for (int i = 0; i < n_seqs; ++i) rep[i] = i;
		return 0;
	}
	std::vector<uint64_t> fp(n_seqs);
	fp_ctx_t c = {seqs, fp.data()};
	if (n_threads > 1 && n_seqs >= 4096) kt_for_nt(fp_worker, &c, n_seqs, n_threads);
	else fp_worker(&c, 0, n_seqs, 0);
	std::vector<dup_slot_t> table(dup_table_slots(n_seqs));
	return dup_table_pass(n_seqs, [seqs](int k) { return seqs + k; }, fp.data(), rep, table.data());
}

bool dup_key_sampled(const bwa_seq_t *p)
{
	return read_sample_key(p) >> (64 - ALN_DUP_SAMPLE_SHIFT) == 0;
}

int dup_sample_by_key(int n_seqs, const bwa_seq_t *seqs)
{
	dup_sample_buffers_t &sb = sample_buffers();
	// identical reads share a key, so a batch of few distinct reads samples all or none of a
	// read's copies; when the sample is far from its expected size, it cannot speak for the batch
	const int64_t expect = (int64_t)n_seqs >> ALN_DUP_SAMPLE_SHIFT;
	const int64_t most = expect * 4 < DUP_SAMPLE_CAP? expect * 4 : DUP_SAMPLE_CAP;
	int n_s = 0;
	for (int i = 0; i < n_seqs; ++i) {
		if (!dup_key_sampled(seqs + i)) continue;
		if (n_s == most) return 1; // over-full: also keeps idx[] within DUP_SAMPLE_CAP
		sb.idx[n_s++] = i;
	}
	if (n_s < expect / 4) return 1;
	return sample_has_dups(seqs, sb, n_s);
}

int dup_sample_by_position(int n_seqs, const bwa_seq_t *seqs)
{
	dup_sample_buffers_t &sb = sample_buffers();
	int64_t stride = (int64_t)1 << ALN_DUP_SAMPLE_SHIFT;
	while ((n_seqs + stride - 1) / stride > DUP_SAMPLE_CAP) stride <<= 1;
	int n_s = 0;
	for (int64_t i = 0; i < n_seqs; i += stride) sb.idx[n_s++] = (int)i;
	return sample_has_dups(seqs, sb, n_s);
}

int dup_sample_check(int n_seqs, const bwa_seq_t *seqs)
{
	if (n_seqs < ALN_DUP_SAMPLE_MIN) return 1; // too few reads to sample; the full pass is small
	return dup_sample_by_key(n_seqs, seqs) || dup_sample_by_position(n_seqs, seqs);
}

} // namespace detail

/* When aln_search_batch looks for duplicates. Samples (detail::dup_sample_check) decide first,
 * so a batch without duplicates costs one key per read and a small sample, not the full pass. A
 * batch the samples or the full pass find with fewer than one duplicate in ALN_DUP_MIN_RATIO reads
 * turns the check off for the stream's next ALN_DUP_RECHECK - 1 batches (aln_dup_state_t; the
 * check itself is a serial pass over the batch while the search threads wait). Only the work
 * changes: a read searched instead of copied gets the same hits. */
/* The duplicate gate and compaction every batch search shares (aln_search_batch; aln --meth's, when
 * it is added, is the second), so each caller gets the same skip window, one-read batches included.
 * search_all() searches the whole batch; search_uniq(n_u, uniq, rep) searches its compacted front,
 * the distinct reads (uniq[k] is the batch index of seqs[k]; rep[i] is read i's first copy); after
 * the batch is restored, copy_dup(i, rep[i]) fills duplicate i from its first copy. Returns the
 * number of duplicates, 0 when the batch was searched whole. */
template <class SearchAll, class SearchUniq, class CopyDup>
static int search_batch_dedup(int n_seqs, bwa_seq_t *seqs, int n_threads, aln_dup_state_t *dup,
							  SearchAll search_all, SearchUniq search_uniq, CopyDup copy_dup)
{
	if (!dup || n_seqs < 2) {
		if (dup && n_seqs == 1 && dup->skip_left > 0) --dup->skip_left; // a batch of the stream too
		search_all();
		return 0;
	}
	if (dup->skip_left > 0) {
		--dup->skip_left;
		search_all();
		return 0;
	}
	if (!detail::dup_sample_check(n_seqs, seqs)) {
		dup->skip_left = ALN_DUP_RECHECK - 1;
		search_all();
		return 0;
	}
	std::vector<int> rep(n_seqs);
	const int n_dups = detail::find_dup_reads(n_seqs, seqs, rep.data(), n_threads);
	if ((int64_t)n_dups * ALN_DUP_MIN_RATIO < n_seqs) {
		// too few to be worth the compaction (the samples were wrong for this batch)
		dup->skip_left = ALN_DUP_RECHECK - 1;
		search_all();
		return 0;
	}
	// Compact the batch in place: the distinct reads move to its front, in order, and the
	// duplicates' entries are set aside (only they take extra memory). The front is searched, then
	// every entry goes back to its own index, so the batch's order and its pointers are as before.
	std::vector<int> uniq;
	std::vector<bwa_seq_t> dups;
	uniq.reserve(n_seqs - n_dups);
	dups.reserve(n_dups);
	for (int i = 0; i < n_seqs; ++i) {
		if (rep[i] != i) {
			dups.push_back(seqs[i]);
			continue;
		}
		// forward: never past an unread entry
		if (i != (int)uniq.size()) seqs[uniq.size()] = seqs[i];
		uniq.push_back(i);
	}
	const int n_u = (int)uniq.size();
	search_uniq(n_u, uniq.data(), rep.data());
	for (int k = n_u - 1; k >= 0; --k) // backward: uniq[k] >= k, past every entry still to move
		if (uniq[k] != k) seqs[uniq[k]] = seqs[k];
	for (int i = 0, j = 0; i < n_seqs; ++i) {
		if (rep[i] == i) continue;
		seqs[i] = dups[j++];
		copy_dup(i, rep[i]);
	}
	return n_dups;
}

int aln_search_batch(const bwt_t *bwt, int n_seqs, bwa_seq_t *seqs, const gap_opt_t *opt, int n_threads,
					 aln_dup_state_t *dup)
{
	gap_opt_t bopt = batch_opt(n_seqs, seqs, opt);
	return search_batch_dedup(n_seqs, seqs, n_threads, dup,
		[&]() { search_reads<AlnStock>(bwt, n_seqs, seqs, &bopt, n_threads); },
		[&](int n_u, const int *, const int *) { search_reads<AlnStock>(bwt, n_u, seqs, &bopt, n_threads); },
		// Every duplicate gets an aln array of its own (its entry held none before the search), a
		// copy of its first copy's, so each read's aln can be freed with the read as before.
		[&](int i, int r) {
			bwa_seq_t *d = seqs + i;
			const bwa_seq_t *q = seqs + r;
			d->sa = 0; d->type = BWA_TYPE_NO_MATCH; d->c1 = d->c2 = 0;
			const size_t n = q->n_aln > 0? q->n_aln : 1;
			d->aln = (bwt_aln1_t*)aln_malloc(n * sizeof(bwt_aln1_t));
			if (q->n_aln > 0) memcpy(d->aln, q->aln, q->n_aln * sizeof(bwt_aln1_t));
			d->n_aln = q->n_aln;
		});
}

} // namespace bwa_aln
