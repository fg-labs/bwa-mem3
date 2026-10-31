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
 * and passed down. The BAM input, .sai writer and command-line parser are not ported. */
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

// `opt` is batch_opt()'s result: max_diff and max_gapo are the batch-wide values.
static void bwa_cal_sa_reg_gap(int tid, int n_threads, bwt_t *const bwt, int n_seqs, bwa_seq_t *seqs,
							   const gap_opt_t *opt)
{
	int i, j, max_l = 0;
	gap_stack_t *stack;
	bwt_width_t *w, *seed_w;
	gap_opt_t local_opt = *opt;
	std::vector<ubyte_t> cseq; // the complemented read; seq itself is left alone

	// initiate priority stack
	stack = gap_init_stack(local_opt.max_diff, local_opt.max_gapo, local_opt.max_gape, &local_opt);

	seed_w = (bwt_width_t*)calloc(opt->seed_len+1, sizeof(bwt_width_t));
	w = 0;
	for (i = 0; i != n_seqs; ++i) {
		bwa_seq_t *p = seqs + i;
		if (i % n_threads != tid) continue;
		p->sa = 0; p->type = BWA_TYPE_NO_MATCH; p->c1 = p->c2 = 0; p->n_aln = 0; p->aln = 0;
		if (max_l < p->len) {
			max_l = p->len;
			w = (bwt_width_t*)realloc(w, (max_l + 1) * sizeof(bwt_width_t));
			memset(w, 0, (max_l + 1) * sizeof(bwt_width_t));
		}
		bwt_cal_width(bwt, p->len, p->seq, w);
		if (opt->fnr > 0.0) local_opt.max_diff = bwa_cal_maxdiff(p->len, BWA_AVG_ERR, opt->fnr);
		local_opt.seed_len = opt->seed_len < p->len? opt->seed_len : 0x7fffffff;
		if (p->len > opt->seed_len)
			bwt_cal_width(bwt, opt->seed_len, p->seq + (p->len - opt->seed_len), seed_w);
		// core function
		cseq.resize(p->len + 1);
		for (j = 0; j < p->len; ++j) // we need to complement
			cseq[j] = p->seq[j] > 3? 4 : 3 - p->seq[j];
		p->aln = bwt_match_gap(bwt, p->len, cseq.data(), w, p->len <= opt->seed_len? 0 : seed_w, &local_opt, &p->n_aln, stack);
	}
	free(seed_w); free(w);
	gap_destroy_stack(stack);
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
	bwa_cal_sa_reg_gap(d->tid, d->n_threads, d->bwt, d->n_seqs, d->seqs, d->opt);
	return 0;
}

void aln_search_batch(const bwt_t *bwt_, int n_seqs, bwa_seq_t *seqs, const gap_opt_t *opt, int n_threads)
{
	bwt_t *bwt = const_cast<bwt_t*>(bwt_); // the search only reads it
	gap_opt_t bopt = batch_opt(n_seqs, seqs, opt);
	if (n_threads <= 1) { // no multi-threading at all
		bwa_cal_sa_reg_gap(0, 1, bwt, n_seqs, seqs, &bopt);
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
