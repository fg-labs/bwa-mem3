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
 * reads' searches (gap_job_t in bwtgap.cpp); its widths are computed by width_walk_run, which
 * is bwt_cal_width's loop as a resumable walk. The BAM input, .sai writer and command-line
 * parser are not ported. */
#include "bwtgap.h"

#include <math.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

namespace bwa_aln {

gap_opt_t *gap_init_opt(void)
{
	gap_opt_t *o;
	o = (gap_opt_t*)calloc(1, sizeof(gap_opt_t));
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

/* New code: bwt_cal_width above as a resumable walk, so that a read's widths can be computed
 * interleaved with other reads' work. width_walk_run runs bwt_cal_width's loop until the next
 * base that needs an Occ lookup, prefetches the two Occ blocks that lookup reads (rows k and
 * l + 1, as in gap_job_run) and returns 0; the next call does the lookup and carries on. It
 * returns 1 once width[0..len] is written, holding exactly what bwt_cal_width writes. */
void width_walk_init(width_walk_t *w, const bwt_t *bwt, int len, const ubyte_t *str, bwt_width_t *width)
{
	w->str = str; w->width = width; w->len = len;
	w->i = 0; w->bid = 0; w->pending = 0;
	w->k = 0; w->l = bwt->seq_len;
}

int width_walk_run(width_walk_t *w, const bwt_t *bwt)
{
	bwtint_t k = w->k, l = w->l, ok, ol;
	int i = w->i, bid = w->bid;
	for (; i < w->len; ++i) {
		ubyte_t c = w->str[i];
		if (c < 4) {
			if (!w->pending) { // yield before the lookup
				w->pending = 1;
				w->i = i; w->bid = bid; w->k = k; w->l = l;
				bwt_prefetch_occ(bwt, k);
				bwt_prefetch_occ(bwt, l + 1);
				return 0;
			}
			w->pending = 0;
			bwt_2occ(bwt, k - 1, l, c, &ok, &ol);
			k = bwt->L2[c] + ok + 1;
			l = bwt->L2[c] + ol;
		}
		if (k > l || c > 3) { // then restart
			k = 0;
			l = bwt->seq_len;
			++bid;
		}
		w->width[i].w = l - k + 1;
		w->width[i].bid = bid;
	}
	w->width[w->len].w = 0;
	w->width[w->len].bid = ++bid;
	w->i = i; w->bid = bid; w->k = k; w->l = l;
	return 1;
}

/* New code: upstream's bwa_cal_sa_reg_gap with several of the thread's reads searched at once
 * (gap_job_*, bwtgap.cpp). Per read, everything bwa_cal_sa_reg_gap does is done in the same
 * way: the same widths, the same per-read max_diff and seed_len, the same complemented read
 * and a stack of its own; only the reads' steps are interleaved, so each read gets the hits
 * bwt_match_gap gives it. A read's two width computations (whole read and seed) are run as
 * width walks, interleaved with the other slots' work, before its search starts. `opt` is
 * batch_opt()'s result: max_diff and max_gapo are the batch-wide values. */
namespace {
struct gap_slot_t {
	gap_job_t job;
	std::vector<bwt_width_t> w, seed_w;
	std::vector<ubyte_t> cseq;
	width_walk_t ww, sw;
	int widths_done, seed_done; // the read's search starts once both walks are done
	bwa_seq_t *p;
};
}

/* Reads searched at once per thread. Measured on Apple M2 Max (hg38, 262,144 aDNA reads,
 * -t 1, default options): 4 reads -36% wall against one at a time, 8 -43%, 12 -43%, 16 -41%,
 * 32 -40%; 8 was also the best at -t 12. */
#define ALN_GROUP 8

/* A slot's stack keeps the capacity its largest read needed; with ALN_GROUP stacks per thread
 * instead of one, that would multiply the stacks' peak memory. After a read whose stack grew
 * past this many entries, the slot's stack buffers are freed (they regrow on demand), so the
 * memory held is close to the one-stack-per-thread search's. Output does not depend on it. */
#define ALN_STACK_KEEP 65536

static void gap_stack_trim(gap_stack_t *stack)
{
	long cap = 0;
	for (int i = 0; i < stack->n_stacks; ++i) cap += stack->stacks[i].m_entries;
	if (cap <= ALN_STACK_KEEP) return;
	for (int i = 0; i < stack->n_stacks; ++i) {
		free(stack->stacks[i].stack);
		stack->stacks[i].stack = 0;
		stack->stacks[i].n_entries = stack->stacks[i].m_entries = 0;
	}
	stack->n_entries = 0;
	stack->best = stack->n_stacks;
}

static void bwa_cal_sa_reg_gap_group(int tid, int n_threads, bwt_t *const bwt, int n_seqs, bwa_seq_t *seqs,
									 const gap_opt_t *opt)
{
	const int G = ALN_GROUP;
	std::vector<gap_slot_t> slots(G);
	for (gap_slot_t &s : slots) {
		s.job.stack = gap_init_stack(opt->max_diff, opt->max_gapo, opt->max_gape, opt);
		s.seed_w.resize(opt->seed_len + 1);
		s.p = 0;
	}
	int next = tid, n_active = 0;
	// Takes the thread's next read into slot s and starts its width walks; returns 0 when
	// none is left.
	auto fill = [&](gap_slot_t &s) -> int {
		if (next >= n_seqs) {
			s.p = 0;
			return 0;
		}
		bwa_seq_t *p = seqs + next;
		next += n_threads;
		p->sa = 0; p->type = BWA_TYPE_NO_MATCH; p->c1 = p->c2 = 0; p->n_aln = 0; p->aln = 0;
		if ((int)s.w.size() < p->len + 1) s.w.resize(p->len + 1);
		width_walk_init(&s.ww, bwt, p->len, p->seq, s.w.data());
		s.widths_done = 0;
		s.seed_done = 1;
		if (p->len > opt->seed_len) {
			width_walk_init(&s.sw, bwt, opt->seed_len, p->seq + (p->len - opt->seed_len), s.seed_w.data());
			s.seed_done = 0;
		}
		s.p = p;
		return 1;
	};
	// Advances slot s by one turn: a step of each unfinished width walk, then the search.
	// Returns 1 when its read is done (aln/n_aln stored).
	auto step = [&](gap_slot_t &s) -> int {
		if (!s.widths_done || !s.seed_done) {
			if (!s.widths_done) s.widths_done = width_walk_run(&s.ww, bwt);
			if (!s.seed_done) s.seed_done = width_walk_run(&s.sw, bwt);
			if (!s.widths_done || !s.seed_done) return 0;
			bwa_seq_t *p = s.p;
			gap_job_t *j = &s.job;
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
			if (!gap_job_start(j, bwt) && !gap_job_run(j, bwt)) return 0;
		} else if (!gap_job_run(&s.job, bwt)) {
			return 0;
		}
		s.p->aln = s.job.aln; s.p->n_aln = s.job.n_aln;
		gap_stack_trim(s.job.stack);
		return 1;
	};
	for (gap_slot_t &s : slots)
		if (fill(s)) ++n_active;
	while (n_active) {
		for (gap_slot_t &s : slots) {
			if (!s.p) continue;
			if (step(s) && !fill(s)) --n_active;
		}
	}
	for (gap_slot_t &s : slots) gap_destroy_stack(s.job.stack);
}

typedef struct {
	int tid, n_threads;
	bwt_t *bwt;
	int n_seqs;
	bwa_seq_t *seqs;
	const gap_opt_t *opt;
} thread_aux_t;

static void *worker(void *data)
{
	thread_aux_t *d = (thread_aux_t*)data;
	bwa_cal_sa_reg_gap_group(d->tid, d->n_threads, d->bwt, d->n_seqs, d->seqs, d->opt);
	return 0;
}

void aln_search_batch(const bwt_t *bwt_, int n_seqs, bwa_seq_t *seqs, const gap_opt_t *opt, int n_threads)
{
	bwt_t *bwt = const_cast<bwt_t*>(bwt_); // the search only reads it
	gap_opt_t bopt = batch_opt(n_seqs, seqs, opt);
	if (n_threads <= 1) { // no multi-threading at all
		bwa_cal_sa_reg_gap_group(0, 1, bwt, n_seqs, seqs, &bopt);
	} else {
		pthread_t *tid;
		pthread_attr_t attr;
		thread_aux_t *data;
		int j;
		pthread_attr_init(&attr);
		pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
		data = (thread_aux_t*)calloc(n_threads, sizeof(thread_aux_t));
		tid = (pthread_t*)calloc(n_threads, sizeof(pthread_t));
		for (j = 0; j < n_threads; ++j) {
			data[j].tid = j; data[j].n_threads = n_threads; data[j].bwt = bwt;
			data[j].n_seqs = n_seqs; data[j].seqs = seqs; data[j].opt = &bopt;
			pthread_create(&tid[j], &attr, worker, data + j);
		}
		for (j = 0; j < n_threads; ++j) pthread_join(tid[j], 0);
		free(data); free(tid);
		pthread_attr_destroy(&attr);
	}
}

} // namespace bwa_aln
