/* Ported from bwa (https://github.com/lh3/bwa), Apache2 branch, commit a08c7c0.
 * Copyright (c) Heng Li and the bwa contributors; licensed under the Apache License,
 * Version 2.0 (see LICENSE in this directory).
 *
 * MODIFIED from the original: this file holds bwape.c's pairing half (sampe without the
 * .sai reading and command line) in namespace bwa_aln. Each read's alignments are taken from
 * its own aln array; the dead hash-table branch is dropped; the reference is always
 * the packed in-memory copy. The mate alignment uses bwa-mem3's ksw_align and ksw_global.
 * bwa_cal_pac_pos_pe is cut into four stages and bwa_paired_sw into a counting core and a report
 * (the two stderr lines), so the driver can run them in parallel; each holds upstream's loop body.
 * See README.md. */
#include "bwape.h"
#include "aln_alloc.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../bwa.h"
#include "../kvec.h"
#include "../ksw.h"
#include "../utils.h"
#include "bwase.h"
#include "bwaseqio.h"
#include "bwtgap.h"

namespace bwa_aln {

typedef struct {
	pair64_v arr;
	const bwt_aln1_t *aln[2];
} pe_data_t;

pe_opt_t *bwa_init_pe_opt(void)
{
	pe_opt_t *po;
	po = (pe_opt_t*)aln_calloc(1, sizeof(pe_opt_t));
	po->max_isize = 500;
	po->force_isize = 0;
	po->max_occ = 100000;
	po->n_multi = 3;
	po->N_multi = 10;
	po->type = BWA_PET_STD;
	po->is_sw = 1;
	po->ap_prior = 1e-5;
	return po;
}

// for normal distribution, this is about 3std
#define OUTLIER_BOUND 2.0

int infer_isize(int n_seqs, bwa_seq_t *seqs[2], isize_info_t *ii, double ap_prior, int64_t L)
{
	uint64_t x, *isizes, n_ap = 0;
	int n, i, tot, p25, p75, p50, max_len = 1, tmp;
	double skewness = 0.0, kurtosis = 0.0, y;

	ii->avg = ii->std = -1.0;
	ii->low = ii->high = ii->high_bayesian = 0;
	isizes = (uint64_t*)aln_calloc(n_seqs, 8);
	for (i = 0, tot = 0; i != n_seqs; ++i) {
		bwa_seq_t *p[2];
		p[0] = seqs[0] + i; p[1] = seqs[1] + i;
		if (p[0]->mapQ >= 20 && p[1]->mapQ >= 20) {
			x = (p[0]->pos < p[1]->pos)? p[1]->pos + p[1]->len - p[0]->pos : p[0]->pos + p[0]->len - p[1]->pos;
			if (x < 100000) isizes[tot++] = x;
		}
		if (p[0]->len > max_len) max_len = p[0]->len;
		if (p[1]->len > max_len) max_len = p[1]->len;
	}
	if (tot < 20) {
		fprintf(stderr, "[infer_isize] fail to infer insert size: too few good pairs\n");
		free(isizes);
		return -1;
	}
	ks_introsort_64(tot, isizes);
	p25 = isizes[(int)(tot*0.25 + 0.5)];
	p50 = isizes[(int)(tot*0.50 + 0.5)];
	p75 = isizes[(int)(tot*0.75 + 0.5)];
	tmp  = (int)(p25 - OUTLIER_BOUND * (p75 - p25) + .499);
	ii->low = tmp > max_len? tmp : max_len; // ii->low is unsigned
	ii->high = (int)(p75 + OUTLIER_BOUND * (p75 - p25) + .499);
	if (ii->low > ii->high) {
		fprintf(stderr, "[infer_isize] fail to infer insert size: upper bound is smaller than read length\n");
		free(isizes);
		return -1;
	}
	for (i = 0, x = n = 0; i < tot; ++i)
		if (isizes[i] >= ii->low && isizes[i] <= ii->high)
			++n, x += isizes[i];
	ii->avg = (double)x / n;
	for (i = 0; i < tot; ++i) {
		if (isizes[i] >= ii->low && isizes[i] <= ii->high) {
			double tmp = (isizes[i] - ii->avg) * (isizes[i] - ii->avg);
			ii->std += tmp;
			skewness += tmp * (isizes[i] - ii->avg);
			kurtosis += tmp * tmp;
		}
	}
	kurtosis = kurtosis/n / (ii->std / n * ii->std / n) - 3;
	ii->std = sqrt(ii->std / n); // it would be better as n-1, but n is usually very large
	skewness = skewness / n / (ii->std * ii->std * ii->std);
	for (y = 1.0; y < 10.0; y += 0.01)
		if (.5 * erfc(y / M_SQRT2) < ap_prior / L * (y * ii->std + ii->avg)) break;
	ii->high_bayesian = (bwtint_t)(y * ii->std + ii->avg + .499);
	for (i = 0; i < tot; ++i)
		if (isizes[i] > ii->high_bayesian) ++n_ap;
	ii->ap_prior = .01 * (n_ap + .01) / tot;
	if (ii->ap_prior < ap_prior) ii->ap_prior = ap_prior;
	free(isizes);
	fprintf(stderr, "[infer_isize] (25, 50, 75) percentile: (%d, %d, %d)\n", p25, p50, p75);
	if (isnan(ii->std) || p75 > 100000) {
		ii->low = ii->high = ii->high_bayesian = 0; ii->avg = ii->std = -1.0;
		fprintf(stderr, "[infer_isize] fail to infer insert size: weird pairing\n");
		return -1;
	}
	for (y = 1.0; y < 10.0; y += 0.01)
		if (.5 * erfc(y / M_SQRT2) < ap_prior / L * (y * ii->std + ii->avg)) break;
	ii->high_bayesian = (bwtint_t)(y * ii->std + ii->avg + .499);
	fprintf(stderr, "[infer_isize] low and high boundaries: %ld and %ld for estimating avg and std\n", (long)ii->low, (long)ii->high);
	fprintf(stderr, "[infer_isize] inferred external isize from %d pairs: %.3lf +/- %.3lf\n", n, ii->avg, ii->std);
	fprintf(stderr, "[infer_isize] skewness: %.3lf; kurtosis: %.3lf; ap_prior: %.2e\n", skewness, kurtosis, ii->ap_prior);
	fprintf(stderr, "[infer_isize] inferred maximum insert size: %ld (%.2lf sigma)\n", (long)ii->high_bayesian, y);
	return 0;
}

static int pairing(bwa_seq_t *p[2], pe_data_t *d, const pe_opt_t *opt, int s_mm, const isize_info_t *ii)
{
	int i, j, o_n, subo_n, cnt_chg = 0, low_bound = ii->low, max_len;
	uint64_t o_score, subo_score;
	pair64_t last_pos[2][2], o_pos[2];
	max_len = p[0]->full_len;
	if (max_len < p[1]->full_len) max_len = p[1]->full_len;
	if (low_bound < max_len) low_bound = max_len;

	// here v>=u. When ii is set, we check insert size with ii; otherwise with opt->max_isize
#define __pairing_aux(u,v) do { \
		bwtint_t l = (v).x + p[(v).y&1]->len - ((u).x); \
		if ((u).x != (uint64_t)-1 && (v).x > (u).x && l >= max_len \
			&& ((ii->high && l <= ii->high_bayesian) || (ii->high == 0 && l <= opt->max_isize))) \
		{ \
			uint64_t s = d->aln[(v).y&1][(v).y>>2].score + d->aln[(u).y&1][(u).y>>2].score; \
			s *= 10; \
			if (ii->high) s += (int)(-4.343 * log(.5 * erfc(M_SQRT1_2 * fabs(l - ii->avg) / ii->std)) + .499); \
			s = s<<32 | (uint32_t)hash_64((u).x<<32 | (v).x); \
			if (s>>32 == o_score>>32) ++o_n; \
			else if (s>>32 < o_score>>32) { subo_n += o_n; o_n = 1; } \
			else ++subo_n; \
			if (s < o_score) subo_score = o_score, o_score = s, o_pos[(u).y&1] = (u), o_pos[(v).y&1] = (v); \
			else if (s < subo_score) subo_score = s; \
		} \
	} while (0)

#define __pairing_aux2(q, w) do { \
		const bwt_aln1_t *r = d->aln[(w).y&1] + ((w).y>>2); \
		(q)->extra_flag |= SAM_FPP; \
		if ((q)->pos != (w).x || (q)->strand != ((w).y>>1&1)) { \
			(q)->n_mm = r->n_mm; (q)->n_gapo = r->n_gapo; (q)->n_gape = r->n_gape; (q)->strand = (w).y>>1&1; \
			(q)->score = r->score; \
			(q)->pos = (w).x; \
			if ((q)->mapQ > 0) ++cnt_chg; \
		} \
	} while (0)

	o_score = subo_score = (uint64_t)-1;
	o_n = subo_n = 0;
	ks_introsort_128(d->arr.n, d->arr.a);
	for (j = 0; j < 2; ++j) last_pos[j][0].x = last_pos[j][0].y = last_pos[j][1].x = last_pos[j][1].y = (uint64_t)-1;
	if (opt->type == BWA_PET_STD) {
		for (i = 0; i < d->arr.n; ++i) {
			pair64_t x = d->arr.a[i];
			int strand = x.y>>1&1;
			if (strand == 1) { // reverse strand, then check
				int y = 1 - (x.y&1);
				__pairing_aux(last_pos[y][1], x);
				__pairing_aux(last_pos[y][0], x);
			} else { // forward strand, then push
				last_pos[x.y&1][0] = last_pos[x.y&1][1];
				last_pos[x.y&1][1] = x;
			}
		}
	} else {
		fprintf(stderr, "[paring] not implemented yet!\n");
		exit(1);
	}
	// set pairing
	//fprintf(stderr, "[%ld, %d, %d, %d]\n", d->arr.n, (int)(o_score>>32), (int)(subo_score>>32), o_n);
	if (o_score != (uint64_t)-1) {
		int mapQ_p = 0; // this is the maximum mapping quality when one end is moved
		//fprintf(stderr, "%d, %d\n", o_n, subo_n);
		if (o_n == 1) {
			if (subo_score == (uint64_t)-1) mapQ_p = 29; // no sub-optimal pair
			else if ((subo_score>>32) - (o_score>>32) > s_mm * 10) mapQ_p = 23; // poor sub-optimal pair
			else {
				int n = subo_n > 255? 255 : subo_n;
				mapQ_p = ((subo_score>>32) - (o_score>>32)) / 2 - g_log_n[n];
				if (mapQ_p < 0) mapQ_p = 0;
			}
		}
		if ((p[0]->pos == o_pos[0].x && p[0]->strand == (o_pos[0].y>>1&1)) && (p[1]->pos == o_pos[1].x && p[1]->strand == (o_pos[1].y>>1&1))) { // both ends not moved
			if (p[0]->mapQ > 0 && p[1]->mapQ > 0) {
				int mapQ = p[0]->mapQ + p[1]->mapQ;
				if (mapQ > 60) mapQ = 60;
				p[0]->mapQ = p[1]->mapQ = mapQ;
			} else {
				if (p[0]->mapQ == 0) p[0]->mapQ = (mapQ_p + 7 < p[1]->mapQ)? mapQ_p + 7 : p[1]->mapQ;
				if (p[1]->mapQ == 0) p[1]->mapQ = (mapQ_p + 7 < p[0]->mapQ)? mapQ_p + 7 : p[0]->mapQ;
			}
		} else if (p[0]->pos == o_pos[0].x && p[0]->strand == (o_pos[0].y>>1&1)) { // [1] moved
			p[1]->seQ = 0; p[1]->mapQ = p[0]->mapQ;
			if (p[1]->mapQ > mapQ_p) p[1]->mapQ = mapQ_p;
		} else if (p[1]->pos == o_pos[1].x && p[1]->strand == (o_pos[1].y>>1&1)) { // [0] moved
			p[0]->seQ = 0; p[0]->mapQ = p[1]->mapQ;
			if (p[0]->mapQ > mapQ_p) p[0]->mapQ = mapQ_p;
		} else { // both ends moved
			p[0]->seQ = p[1]->seQ = 0;
			mapQ_p -= 20;
			if (mapQ_p < 0) mapQ_p = 0;
			p[0]->mapQ = p[1]->mapQ = mapQ_p;
		}
		__pairing_aux2(p[0], o_pos[0]);
		__pairing_aux2(p[1], o_pos[1]);
	}
	return cnt_chg;
}

/* bwa_cal_pac_pos_pe is upstream's single function cut into its stages, so that the driver can
 * run the RNG-free stages in parallel (see bwape.h). Stage bodies are upstream's, in order. */

void bwa_pe_aln2seq(int n_seqs, bwa_seq_t *seqs[2])
{
	int i, j;
	for (i = 0; i != n_seqs; ++i) {
		for (j = 0; j < 2; ++j) {
			bwa_seq_t *p = seqs[j] + i;
			p->n_multi = 0;
			p->extra_flag |= SAM_FPD | (j == 0? SAM_FR1 : SAM_FR2);
			// generate SE alignment and mapping quality
			bwa_aln2seq(p->n_aln, p->aln, p);
		}
	}
}

void bwa_pe_se_pos(const bntseq_t *bns, const bwt_t *bwt, int n_seqs, bwa_seq_t *seqs[2], const gap_opt_t *gopt)
{
	int i, j;
	for (i = 0; i != n_seqs; ++i) {
		bwa_seq_t *p[2];
		for (j = 0; j < 2; ++j) {
			p[j] = seqs[j] + i;
			if (p[j]->type == BWA_TYPE_UNIQUE || p[j]->type == BWA_TYPE_REPEAT) {
				int strand;
				int max_diff = gopt->fnr > 0.0? bwa_cal_maxdiff(p[j]->len, BWA_AVG_ERR, gopt->fnr) : gopt->max_diff;
				p[j]->seQ = p[j]->mapQ = bwa_approx_mapQ(p[j], max_diff);
				p[j]->pos = bwa_sa2pos(bns, bwt, p[j]->sa, p[j]->len + p[j]->ref_shift, &strand);
				p[j]->strand = strand;
				if (p[j]->pos == (bwtint_t)-1) p[j]->type = BWA_TYPE_NO_MATCH;
			}
		}
	}
}

void bwa_pe_isize(const bwt_t *bwt, int n_seqs, bwa_seq_t *seqs[2], isize_info_t *ii, const pe_opt_t *opt,
				  const isize_info_t *last_ii)
{
	// infer isize
	infer_isize(n_seqs, seqs, ii, opt->ap_prior, bwt->seq_len/2);
	if (ii->avg < 0.0 && last_ii->avg > 0.0) *ii = *last_ii;
	if (opt->force_isize) {
		fprintf(stderr, "[%s] discard insert size estimate as user's request.\n", "bwa_cal_pac_pos_pe");
		ii->low = ii->high = 0; ii->avg = ii->std = -1.0;
	}
}

int bwa_pe_pair(const bntseq_t *bns, const bwt_t *bwt, int n_seqs, bwa_seq_t *seqs[2],
				const isize_info_t *ii, const pe_opt_t *opt, const gap_opt_t *gopt)
{
	int i, j, cnt_chg = 0;
	pe_data_t *d;

	d = (pe_data_t*)aln_calloc(1, sizeof(pe_data_t));

	// PE
	for (i = 0; i != n_seqs; ++i) {
		bwa_seq_t *p[2];
		for (j = 0; j < 2; ++j) {
			p[j] = seqs[j] + i;
			d->aln[j] = p[j]->aln;
		}
		if ((p[0]->type == BWA_TYPE_UNIQUE || p[0]->type == BWA_TYPE_REPEAT)
			&& (p[1]->type == BWA_TYPE_UNIQUE || p[1]->type == BWA_TYPE_REPEAT))
		{ // only when both ends mapped
			pair64_t x;
			int j, k;
			long long n_occ[2];
			for (j = 0; j < 2; ++j) {
				n_occ[j] = 0;
				for (k = 0; k < p[j]->n_aln; ++k)
					n_occ[j] += d->aln[j][k].l - d->aln[j][k].k + 1;
			}
			if (n_occ[0] > opt->max_occ || n_occ[1] > opt->max_occ) continue;
			d->arr.n = 0;
			for (j = 0; j < 2; ++j) {
				for (k = 0; k < p[j]->n_aln; ++k) {
					const bwt_aln1_t *r = d->aln[j] + k;
					bwtint_t l;
					for (l = r->k; l <= r->l; ++l) { // calculate on the fly
						int strand;
						x.x = bwa_sa2pos(bns, bwt, l, p[j]->len + p[j]->ref_shift, &strand);
						x.y = k<<2 | strand<<1 | j;
						kv_push(pair64_t, d->arr, x);
					}
				}
			}
			cnt_chg += pairing(p, d, opt, gopt->s_mm, ii);
		}

		if (opt->N_multi || opt->n_multi) {
			for (j = 0; j < 2; ++j) {
				if (p[j]->type != BWA_TYPE_NO_MATCH) {
					int k, n_multi;
					if (!(p[j]->extra_flag&SAM_FPP) && p[1-j]->type != BWA_TYPE_NO_MATCH) {
						bwa_aln2seq_core(p[j]->n_aln, d->aln[j], p[j], 0, p[j]->c1+p[j]->c2-1 > opt->N_multi? opt->n_multi : opt->N_multi);
					} else bwa_aln2seq_core(p[j]->n_aln, d->aln[j], p[j], 0, opt->n_multi);
					for (k = 0, n_multi = 0; k < p[j]->n_multi; ++k) {
						int strand;
						bwt_multi1_t *q = p[j]->multi + k;
						q->pos = bwa_sa2pos(bns, bwt, q->pos, p[j]->len + q->ref_shift, &strand);
						q->strand = strand;
						if (q->pos != p[j]->pos && q->pos != (bwtint_t)-1)
							p[j]->multi[n_multi++] = *q;
					}
					p[j]->n_multi = n_multi;
				}
			}
		}
	}

	kv_destroy(d->arr);
	free(d);
	return cnt_chg;
}

int bwa_cal_pac_pos_pe(const bntseq_t *bns, const bwt_t *bwt, int n_seqs, bwa_seq_t *seqs[2],
					   isize_info_t *ii, const pe_opt_t *opt, const gap_opt_t *gopt,
					   const isize_info_t *last_ii)
{
	// SE
	bwa_pe_aln2seq(n_seqs, seqs);
	bwa_pe_se_pos(bns, bwt, n_seqs, seqs, gopt);
	bwa_pe_isize(bwt, n_seqs, seqs, ii, opt, last_ii);
	return bwa_pe_pair(bns, bwt, n_seqs, seqs, ii, opt, gopt);
}

#define SW_MIN_MATCH_LEN 20
#define SW_MIN_MAPQ 17

// cnt = n_mm<<16 | n_gapo<<8 | n_gape
static bwa_cigar_t *bwa_sw_core(bwtint_t l_pac, const ubyte_t *pacseq, int len, const ubyte_t *seq, int64_t *beg, int reglen, int *n_cigar, uint32_t *_cnt)
{
	kswr_t r;
	uint32_t *cigar32 = 0;
	bwa_cigar_t *cigar = 0;
	ubyte_t *ref_seq;
	bwtint_t k, x, y, l;
	int xtra, gscore;
	int8_t mat[25];

	bwa_fill_scmat(1, 3, mat);
	// check whether there are too many N's
	if (reglen < SW_MIN_MATCH_LEN || (int64_t)l_pac - *beg < len) return 0;
	for (k = 0, x = 0; k < len; ++k)
		if (seq[k] >= 4) ++x;
	if ((float)x/len >= 0.25 || len - x < SW_MIN_MATCH_LEN) return 0;

	// get reference subsequence
	ref_seq = (ubyte_t*)aln_calloc(reglen, 1);
	for (k = *beg, l = 0; l < reglen && k < l_pac; ++k)
		ref_seq[l++] = pacseq[k>>2] >> ((~k&3)<<1) & 3;

	// do alignment
	xtra = KSW_XSUBO | KSW_XSTART | (len < 250? KSW_XBYTE : 0);
	r = ksw_align(len, (uint8_t*)seq, l, ref_seq, 5, mat, 5, 1, xtra, 0);
	gscore = ksw_global(r.qe - r.qb + 1, &seq[r.qb], r.te - r.tb + 1, &ref_seq[r.tb], 5, mat, 5, 1, 50, n_cigar, &cigar32);
	cigar = (bwa_cigar_t*)cigar32;
	for (k = 0; k < *n_cigar; ++k)
		cigar[k] = __cigar_create((cigar32[k]&0xf), (cigar32[k]>>4));

	if (r.score < SW_MIN_MATCH_LEN || r.score2 == r.score || gscore != r.score) { // poor hit or tandem hits or weird alignment
		free(cigar); free(ref_seq); *n_cigar = 0;
		return 0;
	}

	// check whether the alignment is good enough
	for (k = 0, x = y = 0; k < *n_cigar; ++k) {
		bwa_cigar_t c = cigar[k];
		if (__cigar_op(c) == FROM_M) x += __cigar_len(c), y += __cigar_len(c);
		else if (__cigar_op(c) == FROM_D) x += __cigar_len(c);
		else y += __cigar_len(c);
	}
	if (x < SW_MIN_MATCH_LEN || y < SW_MIN_MATCH_LEN) { // not good enough
		free(cigar); free(ref_seq);
		*n_cigar = 0;
		return 0;
	}

	{ // update cigar and coordinate;
		int start = r.qb, end = r.qe + 1;
		*beg += r.tb;
		cigar = (bwa_cigar_t*)aln_realloc(cigar, sizeof(bwa_cigar_t) * (*n_cigar + 2));
		if (start) {
			memmove(cigar + 1, cigar, sizeof(bwa_cigar_t) * (*n_cigar));
			cigar[0] = __cigar_create(3, start);
			++(*n_cigar);
		}
		if (end < len) {
			/*cigar[*n_cigar] = 3<<14 | (len - end);*/
			cigar[*n_cigar] = __cigar_create(3, (len - end));
			++(*n_cigar);
		}
	}

	{ // set *cnt
		int n_mm, n_gapo, n_gape;
		n_mm = n_gapo = n_gape = 0;
		x = r.tb; y = r.qb;
		for (k = 0; k < *n_cigar; ++k) {
			bwa_cigar_t c = cigar[k];
			if (__cigar_op(c) == FROM_M) {
				for (l = 0; l < (__cigar_len(c)); ++l)
					if (ref_seq[x+l] < 4 && seq[y+l] < 4 && ref_seq[x+l] != seq[y+l]) ++n_mm;
				x += __cigar_len(c), y += __cigar_len(c);
			} else if (__cigar_op(c) == FROM_D) {
				x += __cigar_len(c), ++n_gapo, n_gape += (__cigar_len(c)) - 1;
			} else if (__cigar_op(c) == FROM_I) {
				y += __cigar_len(c), ++n_gapo, n_gape += (__cigar_len(c)) - 1;
			}
		}
		*_cnt = (uint32_t)n_mm<<16 | n_gapo<<8 | n_gape;
	}
	
	free(ref_seq);
	return cigar;
}

void bwa_paired_sw_core(const bntseq_t *bns, const ubyte_t *_pacseq, int n_seqs, bwa_seq_t *seqs[2], const pe_opt_t *popt,
						const isize_info_t *ii, uint64_t n_tot[2], uint64_t n_mapped[2])
{
	ubyte_t *pacseq = (ubyte_t*)_pacseq;
	int i;

	if (!popt->is_sw || ii->avg < 0.0) return;

	// perform mate alignment
	for (i = 0; i != n_seqs; ++i) {
		bwa_seq_t *p[2];
		p[0] = seqs[0] + i; p[1] = seqs[1] + i;
		if ((p[0]->mapQ >= SW_MIN_MAPQ || p[1]->mapQ >= SW_MIN_MAPQ) && (p[0]->extra_flag&SAM_FPP) == 0) { // unpaired and one read has high mapQ
			int k, n_cigar[2], is_singleton, mapQ = 0, mq_adjust[2];
			int64_t beg[2], end[2];
			bwa_cigar_t *cigar[2];
			uint32_t cnt[2];

			/* In the following, _pref points to the reference read
			 * which must be aligned; _pmate points to its mate which is
			 * considered to be modified. */

#define __set_rght_coor(_a, _b, _pref, _pmate) do {						\
				(_a) = (int64_t)_pref->pos + ii->avg - 3 * ii->std - _pmate->len * 1.5; \
				(_b) = (_a) + 6 * ii->std + 2 * _pmate->len;			\
				if ((_a) < (int64_t)_pref->pos + _pref->len) (_a) = _pref->pos + _pref->len; \
				if ((_b) > bns->l_pac) (_b) = bns->l_pac;				\
			} while (0)

#define __set_left_coor(_a, _b, _pref, _pmate) do {						\
				(_a) = (int64_t)_pref->pos + _pref->len - ii->avg - 3 * ii->std - _pmate->len * 0.5; \
				(_b) = (_a) + 6 * ii->std + 2 * _pmate->len;			\
				if ((_a) < 0) (_a) = 0;									\
				if ((_b) > _pref->pos) (_b) = _pref->pos;				\
			} while (0)
			
#define __set_fixed(_pref, _pmate, _beg, _cnt) do {						\
				_pmate->type = BWA_TYPE_MATESW;							\
				_pmate->pos = _beg;										\
				_pmate->seQ = _pref->seQ;								\
				_pmate->strand = (popt->type == BWA_PET_STD)? 1 - _pref->strand : _pref->strand; \
				_pmate->n_mm = _cnt>>16; _pmate->n_gapo = _cnt>>8&0xff; _pmate->n_gape = _cnt&0xff; \
				_pmate->extra_flag |= SAM_FPP;							\
				_pref->extra_flag |= SAM_FPP;							\
			} while (0)

			mq_adjust[0] = mq_adjust[1] = 255; // not effective
			is_singleton = (p[0]->type == BWA_TYPE_NO_MATCH || p[1]->type == BWA_TYPE_NO_MATCH)? 1 : 0;

			++n_tot[is_singleton];
			cigar[0] = cigar[1] = 0;
			n_cigar[0] = n_cigar[1] = 0;
			if (popt->type != BWA_PET_STD) continue; // other types of pairing is not considered
			for (k = 0; k < 2; ++k) { // p[1-k] is the reference read and p[k] is the read considered to be modified
				ubyte_t *seq;
				if (p[1-k]->type == BWA_TYPE_NO_MATCH) continue; // if p[1-k] is unmapped, skip
				{ // note that popt->type == BWA_PET_STD always true; in older versions, there was a branch for color-space FF/RR reads
					if (p[1-k]->strand == 0) { // then the mate is on the reverse strand and has larger coordinate
						__set_rght_coor(beg[k], end[k], p[1-k], p[k]);
						seq = p[k]->rseq;
					} else { // then the mate is on forward stand and has smaller coordinate
						__set_left_coor(beg[k], end[k], p[1-k], p[k]);
						seq = p[k]->seq;
						seq_reverse(p[k]->len, seq, 0); // because ->seq is reversed; this will reversed back shortly
					}
				}
				// perform SW alignment
				cigar[k] = bwa_sw_core(bns->l_pac, pacseq, p[k]->len, seq, &beg[k], end[k] - beg[k], &n_cigar[k], &cnt[k]);
				if (cigar[k] && p[k]->type != BWA_TYPE_NO_MATCH) { // re-evaluate cigar[k]
					int s_old, clip = 0, s_new;
					if (__cigar_op(cigar[k][0]) == 3) clip += __cigar_len(cigar[k][0]);
					if (__cigar_op(cigar[k][n_cigar[k]-1]) == 3) clip += __cigar_len(cigar[k][n_cigar[k]-1]);
					s_old = (int)((p[k]->n_mm * 9 + p[k]->n_gapo * 13 + p[k]->n_gape * 2) / 3. * 8. + .499);
					s_new = (int)(((cnt[k]>>16) * 9 + (cnt[k]>>8&0xff) * 13 + (cnt[k]&0xff) * 2 + clip * 3) / 3. * 8. + .499);
					s_old += -4.343 * log(ii->ap_prior / bns->l_pac);
					s_new += (int)(-4.343 * log(.5 * erfc(M_SQRT1_2 * 1.5) + .499)); // assume the mapped isize is 1.5\sigma
					if (s_old < s_new) { // reject SW alignment
						mq_adjust[k] = s_new - s_old;
						free(cigar[k]); cigar[k] = 0; n_cigar[k] = 0;
					} else mq_adjust[k] = s_old - s_new;
				}
				// now revserse sequence back such that p[*]->seq looks untouched
				if (popt->type == BWA_PET_STD) {
					if (p[1-k]->strand == 1) seq_reverse(p[k]->len, seq, 0);
				} else {
					if (p[1-k]->strand == 0) seq_reverse(p[k]->len, seq, 0);
				}
			}
			k = -1; // no read to be changed
			if (cigar[0] && cigar[1]) {
				k = p[0]->mapQ < p[1]->mapQ? 0 : 1; // p[k] to be fixed
				mapQ = abs(p[1]->mapQ - p[0]->mapQ);
			} else if (cigar[0]) k = 0, mapQ = p[1]->mapQ;
			else if (cigar[1]) k = 1, mapQ = p[0]->mapQ;
			if (k >= 0 && p[k]->pos != beg[k]) {
				++n_mapped[is_singleton];
				{ // recalculate mapping quality
					int tmp = (int)p[1-k]->mapQ - p[k]->mapQ/2 - 8;
					if (tmp <= 0) tmp = 1;
					if (mapQ > tmp) mapQ = tmp;
					p[k]->mapQ = p[1-k]->mapQ = mapQ;
					p[k]->seQ = p[1-k]->seQ = p[1-k]->seQ < mapQ? p[1-k]->seQ : mapQ;
					if (p[k]->mapQ > mq_adjust[k]) p[k]->mapQ = mq_adjust[k];
					if (p[k]->seQ > mq_adjust[k]) p[k]->seQ = mq_adjust[k];
				}
				// update CIGAR
				free(p[k]->cigar); p[k]->cigar = cigar[k]; cigar[k] = 0;
				p[k]->n_cigar = n_cigar[k];
				// update the rest of information
				__set_fixed(p[1-k], p[k], beg[k], cnt[k]);
			}
			free(cigar[0]); free(cigar[1]);
		}
	}
}

ubyte_t *bwa_paired_sw(const bntseq_t *bns, const ubyte_t *_pacseq, int n_seqs, bwa_seq_t *seqs[2], const pe_opt_t *popt, const isize_info_t *ii)
{
	ubyte_t *pacseq = (ubyte_t*)_pacseq;
	uint64_t n_tot[2], n_mapped[2];

	if (!popt->is_sw || ii->avg < 0.0) return pacseq;
	n_tot[0] = n_tot[1] = n_mapped[0] = n_mapped[1] = 0;
	bwa_paired_sw_core(bns, _pacseq, n_seqs, seqs, popt, ii, n_tot, n_mapped);
	bwa_paired_sw_report(n_tot, n_mapped);
	return pacseq;
}

void bwa_paired_sw_report(const uint64_t n_tot[2], const uint64_t n_mapped[2])
{
	fprintf(stderr, "[bwa_paired_sw] %lld out of %lld Q%d singletons are mated.\n",
			(long long)n_mapped[1], (long long)n_tot[1], SW_MIN_MAPQ);
	fprintf(stderr, "[bwa_paired_sw] %lld out of %lld Q%d discordant pairs are fixed.\n",
			(long long)n_mapped[0], (long long)n_tot[0], SW_MIN_MAPQ);
}

} // namespace bwa_aln
