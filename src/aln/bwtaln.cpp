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
 * The BAM input, .sai writer and command-line parser are not ported. */
#include "bwtgap.h"
#include "aln_alloc.h"

#include <math.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>
#include <atomic>
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

/* Reads searched at once per thread. Only the scheduling depends on it; each read's hits do not.
 * Measured at -t 16 on hg38 (1M aDNA reads at default options and at -l 1024 -n 0.04 -o 2, 1M
 * merged reads, 2M 20-40 bp reads, 2M 2x36 bp pairs, 1.09M real pairs) with the chunked stacks:
 * - Apple M3 Ultra: 12 was 4-11% faster in wall time than 8; 16 was faster than 12 only on the
 *   20-40 bp reads and the pairs (1.5-3.6%) and slower on the other three.
 * - Intel Sapphire Rapids (16 vCPU, 8 cores with SMT, AVX-512 and AVX2 tiers): 8 was 2.5-3.5%
 *   faster than 12 on four of the six workloads and within 2.1% on the other two.
 * No other host was swept. So 12 is used only on Apple arm64, and every other target, Linux
 * arm64 included, keeps 8, the depth before the chunked stacks. */
#if defined(__APPLE__) && defined(__aarch64__)
#define ALN_GROUP 12
#else
#define ALN_GROUP 8
#endif

int read_claim_t::claim(int *end)
{
	const int st = next_.fetch_add(ALN_CLAIM_CHUNK, std::memory_order_relaxed);
	if (st >= n_seqs_) return -1;
	*end = st + ALN_CLAIM_CHUNK < n_seqs_? st + ALN_CLAIM_CHUNK : n_seqs_;
	return st;
}

static void bwa_cal_sa_reg_gap_group(read_claim_t *claim, const bwt_t *bwt, int n_seqs, bwa_seq_t *seqs,
									 const gap_opt_t *opt)
{
	const int G = ALN_GROUP;
	std::vector<gap_slot_t> slots(G);
	std::vector<gap_slot_core_t> cores(G);
	gap_pool_t pool;
	gap_init_pool(&pool, ALN_POOL_KEEP);
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
		if ((int)s.w.size() < p->len + 1) s.w.resize(p->len + 1);
		width_walk_init(&s.c->ww, bwt, p->len, p->seq, s.w.data());
		s.c->widths_done = 0;
		s.c->seed_done = 1;
		if (p->len > opt->seed_len) {
			// sized here, not per slot up front: a large -l (seeding off) would otherwise cost
			// (-l + 1) x 16 B per slot even though no read is longer; here -l < p->len
			if (s.seed_w.size() < (size_t)opt->seed_len + 1) s.seed_w.resize((size_t)opt->seed_len + 1);
			width_walk_init(&s.c->sw, bwt, opt->seed_len, p->seq + (p->len - opt->seed_len), s.seed_w.data());
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
			s.cseq[k] = p->seq[k] > 3? 4 : 3 - p->seq[k];
		j->len = p->len;
		j->seq = s.cseq.data();
		j->width = s.w.data();
		j->seed_width = p->len <= opt->seed_len? 0 : s.seed_w.data();
		s.c->state = GAP_SLOT_SEARCH;
		return gap_job_start(j, bwt);
	};
	for (gap_slot_t &s : slots)
		if (fill(s)) ++n_active;
	int pos = G - 1;
	while (n_active) {
		gap_slot_t &s = slots[gap_slots_run(cores.data(), G, &pos, bwt)];
		if (s.c->state == GAP_SLOT_WIDTHS && !start(s)) continue; // its search runs in the next turns
		s.p->aln = s.c->job.aln; s.p->n_aln = s.c->job.n_aln; // the read is done
		if (!fill(s)) --n_active;
	}
	for (gap_slot_t &s : slots) gap_destroy_cstack(s.c->job.stack);
	gap_destroy_pool(&pool);
}

typedef struct {
	read_claim_t *claim;
	const bwt_t *bwt;
	int n_seqs;
	bwa_seq_t *seqs;
	const gap_opt_t *opt;
} thread_aux_t;

static void *worker(void *data)
{
	thread_aux_t *d = (thread_aux_t*)data;
	bwa_cal_sa_reg_gap_group(d->claim, d->bwt, d->n_seqs, d->seqs, d->opt);
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
void aln_search_batch(const bwt_t *bwt, int n_seqs, bwa_seq_t *seqs, const gap_opt_t *opt, int n_threads)
{
	gap_opt_t bopt = batch_opt(n_seqs, seqs, opt);
	read_claim_t claim(n_seqs);
	if (n_threads <= 1) { // no multi-threading at all (bwa aln's -t 0 too)
		bwa_cal_sa_reg_gap_group(&claim, bwt, n_seqs, seqs, &bopt);
	} else {
		pthread_t *tid;
		pthread_attr_t attr;
		thread_aux_t *data;
		int j;
		pthread_attr_init(&attr);
		pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
		data = (thread_aux_t*)aln_calloc(n_threads, sizeof(thread_aux_t));
		tid = (pthread_t*)aln_calloc(n_threads, sizeof(pthread_t));
		std::vector<char> started(n_threads, 0);
		bool any_failed = false;
		for (j = 0; j < n_threads; ++j) {
			data[j].claim = &claim; data[j].bwt = bwt;
			data[j].n_seqs = n_seqs; data[j].seqs = seqs; data[j].opt = &bopt;
			started[j] = pthread_create(&tid[j], &attr, worker, data + j) == 0;
			any_failed = any_failed || !started[j];
		}
		// A thread that could not be created must not be joined. The caller runs one worker in
		// place of the missing ones: it claims reads like the others, so it and the threads that
		// started search every read (it alone does if none started).
		if (any_failed) worker(data);
		for (j = 0; j < n_threads; ++j)
			if (started[j]) pthread_join(tid[j], 0);
		free(data); free(tid);
		pthread_attr_destroy(&attr);
	}
}

} // namespace bwa_aln
