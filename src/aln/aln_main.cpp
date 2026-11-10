/* SPDX-License-Identifier: MIT */
/* New code for bwa-mem3 (not from upstream bwa): the `aln` subcommand driver. The option
 * set is stock `bwa aln`'s plus the samse/sampe -n (as --max-xa-hits), -f and -r, and sampe's
 * paired-end options; see README.md. */
#include "aln_main.h"
#include "aln_known_snvs.h"
#include "aln_meth_realign.h"

#include <assert.h>
#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <limits.h>
#include <sys/stat.h>

#include <algorithm>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include <pthread.h>

#include "../FMI_search.h"
#include "../bwa.h"
#include "../bwa_shm.h"
#include "../kthread.h"
#include "../macro.h"
#include "../read_index_ele.h"
#include "aln_bwt.h"
#include "bwape.h"
#include "bwase.h"
#include "bwaseqio.h"
#include "aln_alloc.h"
#include "bwtgap.h"
#include "aln_meth.h"
#include "../meth_opts.h"

extern char *bwa_pg;

namespace bwa_aln {

void aln_opt_init(aln_opt_t *o)
{
	gap_opt_t *g = gap_init_opt();
	o->gap = *g;
	free(g);
	o->max_xa_hits = 3;
	pe_opt_t *pe = bwa_init_pe_opt();
	o->pe = *pe;
	free(pe);
	o->batch_size = 0x40000;
	o->share_dup_search = 1;
	o->meth = nullptr;
	o->meth_ref = nullptr;
}

/* The next batch is read on a helper thread while the current one is processed (background_task_t). */

/* The stages after the search run in parallel on kt_for (stage order, as follows): every drand48() draw
 * stays in one serial pass, in stock order, before them; infer_isize is a barrier between the
 * paired-end passes; the records are formatted in parallel and written in read order. */

namespace {

/* Joins read 2's barcode onto read 1's and gives both mates the result, keeping at most
 * BWA_MAX_BCLEN bases: main_aln limits -B to half that with two files, but the paired drivers
 * are public entry points. */
void join_pair_barcodes(bwa_seq_t *p0, bwa_seq_t *p1)
{
	if (!p0->bc[0] && !p1->bc[0]) return;
	const size_t l0 = strlen(p0->bc);
	snprintf(p0->bc + l0, sizeof p0->bc - l0, "%s", p1->bc);
	memcpy(p1->bc, p0->bc, sizeof p0->bc);
}

/* One formatted run of consecutive reads (or pairs) [st, st + n): `len` bytes at `off` in
 * thread `tid`'s buffer. `stop` >= 0: the pair at index `stop` has mismatched names, and the
 * run ends with that pair's records. */
struct fmt_seg_t {
	int st, tid, stop;
	size_t off, len;
};

/* Per-thread output buffers and segment lists for one pass, reused across batches. The
 * records are written after the whole pass rather than through kt_for_prefix as they finish:
 * kt_for_prefix hands on_prefix only (st, ed), so each item would still need its own buffer,
 * the items it never offers need a tail write after it returns, and the paired-end stop at
 * mismatched names would have to be applied in both places. Formatting is a small share of
 * the run next to the search, so the simpler write-after-pass is kept. */
struct fmt_bufs_t {
	std::vector<kstring_t> buf;
	std::vector<std::vector<fmt_seg_t>> segs;
	explicit fmt_bufs_t(int n_threads) : buf(n_threads), segs(n_threads)
	{
		for (auto &k : buf) k.l = k.m = 0, k.s = 0;
	}
	~fmt_bufs_t()
	{
		for (auto &k : buf) free(k.s);
	}
	void reset()
	{
		for (auto &k : buf) k.l = 0;
		for (auto &v : segs) v.clear();
	}
	/* Writes the segments in read order, up to and including the first stopped one. Sets
	 * *stop to that segment's stopping index (or -1). Returns 0 or -1 on a write error. */
	int write(aln_out_t *out, int *stop)
	{
		std::vector<fmt_seg_t> all;
		for (auto &v : segs) all.insert(all.end(), v.begin(), v.end());
		std::sort(all.begin(), all.end(), [](const fmt_seg_t &x, const fmt_seg_t &y) { return x.st < y.st; });
		*stop = -1;
		for (const fmt_seg_t &g : all) {
			if (aln_out_write_buf(out, buf[g.tid].s + g.off, g.len) != 0) return -1;
			if (g.stop >= 0) {
				*stop = g.stop;
				break;
			}
		}
		return 0;
	}
};

/* Runs one function on a helper thread while the caller works; wait() joins it. If the thread
 * cannot be created, wait() runs the function itself, so the result is the same either way. The
 * batch loops use it to read and parse the next batch while the current one is aligned: the
 * reader is then the only user of the input streams, and the batches are still aligned and
 * written in input order. The reader is one thread on top of the -t threads. A read error does not
 * exit on it (bwa_seq_defer_read_errors): the loop reports it after the current batch is written,
 * where stock's read of the next batch would have failed, and the run returns 1 (stock exits
 * there; returning lets the caller close the output, as for a failed close of the input). */
class background_task_t {
public:
	background_task_t() = default;
	background_task_t(const background_task_t &) = delete;
	background_task_t &operator=(const background_task_t &) = delete;
	~background_task_t() // a task still running (a path that skipped wait) is joined, never left behind
	{
		if (started_) pthread_join(th_, 0);
	}
	void start(std::function<void()> fn)
	{
		assert(!pending_); // one task at a time: wait() before the next start()
		fn_ = std::move(fn);
		pending_ = true;
		started_ = pthread_create(&th_, 0, run, this) == 0;
	}
	void wait()
	{
		if (!pending_) return; // nothing started
		if (started_) pthread_join(th_, 0);
		else fn_(); // the thread could not be created: run the task here
		started_ = pending_ = false;
	}
	// The task's result is not needed (the run stops early): a task that has not run, because the
	// thread could not be created, is dropped without running; a running one is joined, as a
	// blocked read cannot be interrupted safely.
	void abandon()
	{
		if (started_) pthread_join(th_, 0);
		started_ = pending_ = false;
	}
private:
	static void *run(void *self)
	{
		((background_task_t*)self)->fn_();
		return 0;
	}
	std::function<void()> fn_;
	pthread_t th_;
	bool started_ = false; // a thread is running fn_ (joined by wait or the destructor)
	bool pending_ = false; // start() was called and wait() not yet
};

struct se_ctx_t {
	const aln_opt_t *o;
	const bwt_t *bwt;
	const bntseq_t *bns;
	ubyte_t *pac;
	bwa_seq_t *seqs;
	fmt_bufs_t *fb;
	bwt_aln1_t **ks_aln = nullptr;         // --known-snvs: each read's window hits (null without)
	const int *ks_n = nullptr;
	std::vector<int64_t> *ks_count = nullptr; // per thread, KS_N_OUTCOMES counters each
	struct ks_rec_t *ks_rec = nullptr;     // per read: its sites, twin and XV (with ks_aln)
	int ks_phase = 0;                      // 0: place and print; 1: place only; 2: print only
};

/* --known-snvs, per read: the sites under its final placement and its twin's outcome. */
struct ks_rec_t {
	std::vector<ks_vote_t> votes;
	int twin_of = -1;            // the vote that gets a twin, or -1
	std::vector<uint8_t> tseq;   // the twin, in read order
	int64_t pos = 0;             // the read's final placement
	int strand = 0;
	int twin_same = -1, twin_mapq = 0;
};

/* Appends a record's XV tag before the newline ending k (--known-snvs). */
void ks_append_xv(kstring_t *k, const bntseq_t *bns, const ks_rec_t &r, std::string *xv)
{
	if (r.votes.empty()) return;
	ks_format_xv(bns, r.votes, r.twin_of, r.twin_same, r.twin_mapq, xv);
	--k->l;
	kputs("\tXV:Z:", k);
	kputs(xv->c_str(), k);
	kputc('\n', k);
}

/* Single end, reads [st, st + n): SA to position, refinement, MD, trim correction and SAM text;
 * bwa_cal_pac_pos (with its SA lookups batched) and bwa_refine_gapped are upstream's batch loops
 * run on the slice. */
void se_post_worker(void *data, int st, int n, int tid)
{
	se_ctx_t *c = (se_ctx_t*)data;
	bwa_seq_t *seqs = c->seqs + st;
	if (c->ks_phase != 2) {
		bwa_cal_pac_pos(c->bns, c->bwt, n, seqs, c->o->gap.max_diff, c->o->gap.fnr);
		std::vector<int64_t> alt_site(c->ks_aln ? n : 0, -1);
		if (c->ks_aln) {
			int64_t *cnt = c->ks_count->data() + (size_t)tid * KS_N_OUTCOMES;
			for (int i = 0; i < n; ++i)
				++cnt[ks_override(c->o->ks_windows, c->bns, c->bwt, &c->o->gap, c->o->max_xa_hits, seqs + i,
								  c->ks_aln[st + i], c->ks_n[st + i], &alt_site[i])];
		}
		bwa_refine_gapped(c->bns, n, seqs, c->pac);
		if (c->ks_rec) {
			for (int i = 0; i < n; ++i) {
				ks_rec_t &r = c->ks_rec[st + i];
				r.twin_of = ks_site_votes(c->o->known_snvs, c->bns, &c->o->gap, seqs + i, alt_site[i], -1, &r.votes);
				if (r.twin_of < 0) continue;
				ks_twin_seq(seqs + i, r.votes[r.twin_of], -1, &r.tseq);
				r.pos = (int64_t)seqs[i].pos;
				r.strand = seqs[i].strand;
			}
		}
		if (c->ks_phase == 1) return;
	}
	kstring_t *k = &c->fb->buf[tid];
	fmt_seg_t g;
	g.st = st, g.tid = tid, g.stop = -1, g.off = k->l;
	std::string xv;
	for (int i = 0; i < n; ++i) {
		bwa_print_sam1(k, c->bns, seqs + i, 0, c->o->gap.mode, c->o->gap.max_top2);
		if (c->ks_rec) ks_append_xv(k, c->bns, c->ks_rec[st + i], &xv); // before the record's newline
	}
	g.len = k->l - g.off;
	c->fb->segs[tid].push_back(g);
}

struct pe_ctx_t {
	const aln_opt_t *o;
	const bwt_t *bwt;
	const bntseq_t *bns;
	ubyte_t *pac;
	bwa_seq_t *seqs[2];
	const isize_info_t *ii;
	fmt_bufs_t *fb;
	std::vector<long long> cnt_chg;       // per thread
	std::vector<uint64_t> sw_cnt;         // per thread: n_tot[0..1], n_mapped[0..1]
	// --known-snvs (all null without): per end, per pair
	bwt_aln1_t **ks_aln[2] = {nullptr, nullptr}; // the window hits
	const int *ks_n[2] = {nullptr, nullptr};
	std::vector<pe_hit_t> *ks_hits[2] = {nullptr, nullptr}; // a credited end's hits, for the pairing
	std::vector<int64_t> *ks_halt[2] = {nullptr, nullptr};  // each of those hits' ALT site, or -1
	int64_t *ks_alt[2] = {nullptr, nullptr};                // the final placement's ALT site, or -1
	ks_rec_t *ks_rec[2] = {nullptr, nullptr};
	std::vector<int64_t> *ks_count = nullptr;               // per thread, KS_N_OUTCOMES counters each
	const pe_alt_sites_t *ks_rescue = nullptr;
	int ks_phase = 0;                                       // as se_ctx_t
};

/* Paired end, pass 1 (pairs [st, st + n)): single-end mapQ and position. */
void pe_pos_worker(void *data, int st, int n, int tid)
{
	(void)tid;
	pe_ctx_t *c = (pe_ctx_t*)data;
	bwa_seq_t *s[2] = {c->seqs[0] + st, c->seqs[1] + st};
	bwa_pe_se_pos(c->bns, c->bwt, n, s, &c->o->gap);
}

/* Paired end, pass 2 (after infer_isize): pairing and extra hits, mate Smith-Waterman,
 * refinement of both reads, then each pair's records, read 1 then read 2. */
void pe_post_worker(void *data, int st, int n, int tid)
{
	pe_ctx_t *c = (pe_ctx_t*)data;
	bwa_seq_t *s[2] = {c->seqs[0] + st, c->seqs[1] + st};
	const bool ks = c->ks_aln[0] != nullptr;
	if (c->ks_phase != 2) {
		// --known-snvs: each end credited after the insert size is estimated from stock MAPQs
		if (ks) {
			int64_t *cnt = c->ks_count->data() + (size_t)tid * KS_N_OUTCOMES;
			for (int j = 0; j < 2; ++j)
				for (int i = 0; i < n; ++i) {
					int64_t a;
					++cnt[ks_override(c->o->ks_windows, c->bns, c->bwt, &c->o->gap, 0, s[j] + i,
									  c->ks_aln[j][st + i], c->ks_n[j][st + i], &a, &c->ks_hits[j][st + i],
									  &c->ks_halt[j][st + i])];
				}
		}
		const std::vector<pe_hit_t> *h[2] = {ks ? c->ks_hits[0] + st : nullptr, ks ? c->ks_hits[1] + st : nullptr};
		c->cnt_chg[tid] += bwa_pe_pair(c->bns, c->bwt, n, s, c->ii, &c->o->pe, &c->o->gap, ks ? h : nullptr);
		int64_t *alt[2] = {nullptr, nullptr};
		if (ks)
			for (int j = 0; j < 2; ++j) {
				alt[j] = c->ks_alt[j] + st;
				for (int i = 0; i < n; ++i) { // the ALT site of the hit the pairing settled on
					const std::vector<pe_hit_t> &hj = c->ks_hits[j][st + i];
					alt[j][i] = -1;
					for (size_t k = 0; k < hj.size(); ++k)
						if ((bwtint_t)hj[k].pos == s[j][i].pos && hj[k].strand == (int)s[j][i].strand) {
							alt[j][i] = c->ks_halt[j][st + i][k];
							break;
						}
				}
			}
		uint64_t *sw = &c->sw_cnt[4 * tid];
		bwa_paired_sw_core(c->bns, c->pac, n, s, &c->o->pe, c->ii, sw, sw + 2, ks ? c->ks_rescue : nullptr,
						   ks ? alt : nullptr);
		for (int j = 0; j < 2; ++j) bwa_refine_gapped(c->bns, n, s[j], c->pac);
		if (c->ks_rec[0])
			for (int j = 0; j < 2; ++j)
				for (int i = 0; i < n; ++i) {
					ks_rec_t &r = c->ks_rec[j][st + i];
					r.twin_of = ks_site_votes(c->o->known_snvs, c->bns, &c->o->gap, s[j] + i, alt[j][i], -1, &r.votes);
					if (r.twin_of < 0) continue;
					ks_twin_seq(s[j] + i, r.votes[r.twin_of], -1, &r.tseq);
					r.pos = (int64_t)s[j][i].pos;
					r.strand = s[j][i].strand;
				}
		if (c->ks_phase == 1) return;
	}
	std::string xv;
	kstring_t *k = &c->fb->buf[tid];
	fmt_seg_t g;
	g.st = st, g.tid = tid, g.stop = -1, g.off = k->l;
	for (int i = 0; i < n; ++i) {
		bwa_seq_t *p0 = s[0] + i, *p1 = s[1] + i;
		join_pair_barcodes(p0, p1);
		bwa_print_sam1(k, c->bns, p0, p1, c->o->gap.mode, c->o->gap.max_top2);
		if (c->ks_rec[0]) ks_append_xv(k, c->bns, c->ks_rec[0][st + i], &xv);
		bwa_print_sam1(k, c->bns, p1, p0, c->o->gap.mode, c->o->gap.max_top2);
		if (c->ks_rec[1]) ks_append_xv(k, c->bns, c->ks_rec[1][st + i], &xv);
		if (strcmp(p0->name, p1->name) != 0) { // as sampe: stop after this pair
			g.stop = st + i;
			break;
		}
	}
	g.len = k->l - g.off;
	c->fb->segs[tid].push_back(g);
}

} // namespace

namespace {

/* The run's --known-snvs summary line (per-thread KS_N_OUTCOMES counters). */
void ks_report(const std::vector<int64_t> &count, const char *what)
{
	int64_t t[KS_N_OUTCOMES] = {0};
	for (size_t i = 0; i < count.size(); ++i) t[i % KS_N_OUTCOMES] += count[i];
	fprintf(stderr, "[M::aln] --known-snvs: %lld %s credited, %lld at a site with no change, "
			"%lld kept stock (too many hits or too long)\n", (long long)t[KS_CREDITED], what,
			(long long)t[KS_SAME], (long long)t[KS_INCOMPLETE]);
}

/* A read for the search from its bases in read order (as bwa_read_seq leaves one). */
void ks_make_read(bwa_seq_t *p, const uint8_t *read_order, int len, int is_comp)
{
	p->tid = -1;
	p->full_len = p->clip_len = p->len = (uint32_t)len;
	p->seq = (ubyte_t*)aln_calloc(len, 1);
	p->rseq = (ubyte_t*)aln_calloc(len, 1);
	memcpy(p->seq, read_order, len);
	memcpy(p->rseq, read_order, len);
	seq_reverse(p->len, p->seq, 0);
	seq_reverse(p->len, p->rseq, is_comp);
	p->name = aln_strdup("twin");
}

/* --known-snvs: every twin of a batch searched as one more batch (main and window index), then
 * placed without a draw, so the stock drand48 stream never sees them. */
void ks_place_twins(const aln_opt_t *o, const bwt_t *bwt, const bntseq_t *bns, std::vector<ks_rec_t> &rec)
{
	std::vector<int> idx;
	for (int i = 0; i < (int)rec.size(); ++i)
		if (rec[i].twin_of >= 0) idx.push_back(i);
	const int nt = (int)idx.size();
	if (nt == 0) return;
	const int is_comp = o->gap.mode & BWA_MODE_COMPREAD;
	bwa_seq_t *tw = (bwa_seq_t*)aln_calloc(nt, sizeof(bwa_seq_t));
	for (int j = 0; j < nt; ++j) ks_make_read(tw + j, rec[idx[j]].tseq.data(), (int)rec[idx[j]].tseq.size(), is_comp);
	const int n_threads = o->gap.n_threads < 1 ? 1 : o->gap.n_threads;
	aln_search_batch(o->ks_windows->bwt, nt, tw, &o->gap, n_threads, nullptr);
	std::vector<bwt_aln1_t*> win_aln(nt);
	std::vector<int> win_n(nt);
	for (int j = 0; j < nt; ++j) {
		win_aln[j] = tw[j].aln; win_n[j] = tw[j].n_aln;
		tw[j].aln = 0; tw[j].n_aln = 0;
	}
	aln_search_batch(bwt, nt, tw, &o->gap, n_threads, nullptr);
	for (int j = 0; j < nt; ++j) {
		ks_rec_t &r = rec[idx[j]];
		r.twin_same = ks_twin_place(o->ks_windows, bns, bwt, &o->gap, tw + j, win_aln[j], win_n[j], r.pos,
									r.strand, &r.twin_mapq);
		free(win_aln[j]);
	}
	bwa_free_read_seq(nt, tw);
}

/* --known-snvs, paired end: each twin travels with a copy of its mate (with the same site swapped
 * too when the mate reads REF or ALT there: the other allele's fragment), and the pairs go through the
 * batch's pipeline (both searches, crediting, pairing with the batch's insert size, ALT-aware
 * rescue) with the primary drawn per read, never from drand48. A twin is "same" when its end ends
 * on the read's strand within KS_TWIN_SLOP of its start; its MAPQ is the paired one. */
struct ks_twin_pe_ctx_t {
	const aln_opt_t *o;
	const bwt_t *bwt;
	const bntseq_t *bns;
	ubyte_t *pac;
	const isize_info_t *ii;
	bwa_seq_t *tw[2];
	bwt_aln1_t **win_aln[2];
	const int *win_n[2];
	const int *side;          // per twin pair: the end that holds the twin
	int *incomplete;          // per twin pair: crediting the twin was incomplete
	const pe_alt_sites_t *rescue;
};

void ks_twin_pe_worker(void *data, int st, int n, int tid)
{
	(void)tid;
	ks_twin_pe_ctx_t *c = (ks_twin_pe_ctx_t*)data;
	bwa_seq_t *s[2] = {c->tw[0] + st, c->tw[1] + st};
	for (int j = 0; j < 2; ++j)
		for (int i = 0; i < n; ++i) {
			bwa_seq_t *p = s[j] + i;
			p->n_multi = 0;
			p->extra_flag |= SAM_FPD | (j == 0 ? SAM_FR1 : SAM_FR2);
			ks_aln2seq_nodraw(c->bns, p);
		}
	bwa_pe_se_pos(c->bns, c->bwt, n, s, &c->o->gap);
	std::vector<std::vector<pe_hit_t>> hits[2];
	std::vector<int64_t> halt;
	for (int j = 0; j < 2; ++j) {
		hits[j].resize(n);
		for (int i = 0; i < n; ++i) {
			int64_t a;
			const ks_outcome_t r = ks_override(c->o->ks_windows, c->bns, c->bwt, &c->o->gap, 0, s[j] + i,
											   c->win_aln[j][st + i], c->win_n[j][st + i], &a, &hits[j][i], &halt);
			if (j == c->side[st + i]) c->incomplete[st + i] = r == KS_INCOMPLETE;
		}
	}
	const std::vector<pe_hit_t> *h[2] = {hits[0].data(), hits[1].data()};
	bwa_pe_pair(c->bns, c->bwt, n, s, c->ii, &c->o->pe, &c->o->gap, h);
	uint64_t sw[4] = {0, 0, 0, 0};
	bwa_paired_sw_core(c->bns, c->pac, n, s, &c->o->pe, c->ii, sw, sw + 2, c->rescue, nullptr);
}

void ks_place_twins_pe(const aln_opt_t *o, const bwt_t *bwt, const bntseq_t *bns, ubyte_t *pac,
					   const isize_info_t *ii, bwa_seq_t *const seqs[2], std::vector<ks_rec_t> rec[2],
					   const pe_alt_sites_t *rescue)
{
	std::vector<std::pair<int, int>> idx; // (pair, end holding the twin)
	for (int i = 0; i < (int)rec[0].size(); ++i)
		for (int j = 0; j < 2; ++j)
			if (rec[j][i].twin_of >= 0) idx.push_back({i, j});
	const int nt = (int)idx.size();
	if (nt == 0) return;
	const int is_comp = o->gap.mode & BWA_MODE_COMPREAD;
	std::vector<int> side(nt), incomplete(nt, 0);
	bwa_seq_t *tw[2];
	for (int j = 0; j < 2; ++j) tw[j] = (bwa_seq_t*)aln_calloc(nt, sizeof(bwa_seq_t));
	for (int t = 0; t < nt; ++t) {
		const int i = idx[t].first, j = idx[t].second;
		side[t] = j;
		const std::vector<uint8_t> &ts = rec[j][i].tseq;
		ks_make_read(tw[j] + t, ts.data(), (int)ts.size(), is_comp);
		const bwa_seq_t *m = seqs[1 - j] + i; // the mate, its bases in read order after refinement
		const ks_rec_t &mr = rec[1 - j][i];
		const ks_site_t *site = rec[j][i].votes[rec[j][i].twin_of].site;
		std::vector<uint8_t> ms;
		for (const ks_vote_t &v : mr.votes)
			if (v.site == site && (v.allele == 'R' || v.allele == 'A')) ks_twin_seq(m, v, -1, &ms);
		if (ms.empty()) ms.assign(m->seq, m->seq + m->len);
		ks_make_read(tw[1 - j] + t, ms.data(), (int)ms.size(), is_comp);
		tw[1 - j][t].full_len = m->full_len;
	}
	const int n_threads = o->gap.n_threads < 1 ? 1 : o->gap.n_threads;
	std::vector<bwt_aln1_t*> win_aln[2];
	std::vector<int> win_n[2];
	for (int j = 0; j < 2; ++j) {
		aln_search_batch(o->ks_windows->bwt, nt, tw[j], &o->gap, n_threads, nullptr);
		win_aln[j].resize(nt);
		win_n[j].resize(nt);
		for (int t = 0; t < nt; ++t) {
			win_aln[j][t] = tw[j][t].aln; win_n[j][t] = tw[j][t].n_aln;
			tw[j][t].aln = 0; tw[j][t].n_aln = 0;
		}
		aln_search_batch(bwt, nt, tw[j], &o->gap, n_threads, nullptr);
	}
	ks_twin_pe_ctx_t c;
	c.o = o, c.bwt = bwt, c.bns = bns, c.pac = pac, c.ii = ii;
	for (int j = 0; j < 2; ++j) c.tw[j] = tw[j], c.win_aln[j] = win_aln[j].data(), c.win_n[j] = win_n[j].data();
	c.side = side.data();
	c.incomplete = incomplete.data();
	c.rescue = rescue;
	kt_for_nt(ks_twin_pe_worker, &c, nt, n_threads);
	for (int t = 0; t < nt; ++t) {
		const int i = idx[t].first, j = idx[t].second;
		ks_rec_t &r = rec[j][i];
		const bwa_seq_t *q = tw[j] + t;
		r.twin_mapq = q->mapQ;
		if (incomplete[t]) r.twin_same = -1;
		else r.twin_same = q->type != BWA_TYPE_NO_MATCH && (int)q->strand == r.strand
			&& std::llabs((int64_t)q->pos - r.pos) <= KS_TWIN_SLOP;
	}
	for (int j = 0; j < 2; ++j) {
		for (bwt_aln1_t *a : win_aln[j]) free(a);
		bwa_free_read_seq(nt, tw[j]);
	}
}

} // namespace

int aln_run_se(const aln_opt_t *o, const bwt_t *bwt, const bntseq_t *bns, ubyte_t *pac,
			   const char *fn_fq, aln_out_t *out)
{
	if (o->batch_size <= 0) { // a batch of no reads would end the run at once, reporting success
		fprintf(stderr, "[E::aln] the batch size must be positive, not %d\n", o->batch_size);
		return 1;
	}
	bwa_seqio_t *ks = bwa_seq_open(fn_fq);
	if (!ks) {
		return 1;
	}
	bwa_seq_defer_read_errors(ks); // reported below, on this thread (see background_task_t)
	bwase_initialize();
	srand48(bns->seed); // once per run, never per batch
	const int n_threads = o->gap.n_threads < 1 ? 1 : o->gap.n_threads;
	fmt_bufs_t fb(n_threads);
	std::vector<int64_t> ks_count((size_t)n_threads * KS_N_OUTCOMES, 0);
	long long tot = 0;
	int ret = 0;
	struct batch_t { bwa_seq_t *seqs; int n; };
	auto read_batch = [&](batch_t *b) { b->seqs = bwa_read_seq(ks, o->batch_size, &b->n, o->gap.mode, o->gap.trim_qual); };
	batch_t cur = {}, next = {};
	read_batch(&cur);
	if (bwa_seq_report_read_error(ks)) ret = 1; // cur.seqs is then 0: no batch is processed
	background_task_t ahead;
	aln_dup_state_t dup_state; // this run's duplicate check (see aln_search_batch)
	aln_dup_state_t *dup = o->share_dup_search? &dup_state : nullptr;
	while (cur.seqs) {
		next = {}; // stays empty if the read is abandoned before it runs
		ahead.start([&] { read_batch(&next); }); // the next batch is read while this one is aligned
		bwa_seq_t *seqs = cur.seqs;
		const int n_seqs = cur.n;
		tot += n_seqs;
		const int n_shared = aln_search_batch(bwt, n_seqs, seqs, &o->gap, o->gap.n_threads, dup);
		if (bwa_verbose >= 4) fprintf(stderr, "[M::aln] %d of %d reads shared an identical read's search.\n", n_shared, n_seqs);
		// --known-snvs: the same search on the ALT windows, kept beside the main search's hits
		std::vector<bwt_aln1_t*> ks_aln;
		std::vector<int> ks_n;
		if (o->ks_windows) {
			ks_aln.resize(n_seqs);
			ks_n.resize(n_seqs);
			for (int i = 0; i < n_seqs; ++i) {
				ks_aln[i] = seqs[i].aln; ks_n[i] = seqs[i].n_aln;
				seqs[i].aln = 0; seqs[i].n_aln = 0;
			}
			aln_search_batch(o->ks_windows->bwt, n_seqs, seqs, &o->gap, o->gap.n_threads, nullptr);
			for (int i = 0; i < n_seqs; ++i) {
				std::swap(ks_aln[i], seqs[i].aln);
				std::swap(ks_n[i], seqs[i].n_aln);
			}
		}
		// serial and in read order: bwa_aln2seq_core draws from drand48
		for (int i = 0; i < n_seqs; ++i)
			bwa_aln2seq_core(seqs[i].n_aln, seqs[i].aln, &seqs[i], 1, o->max_xa_hits);
		se_ctx_t c = {o, bwt, bns, pac, seqs, &fb};
		std::vector<ks_rec_t> ks_rec;
		if (o->ks_windows) {
			c.ks_aln = ks_aln.data();
			c.ks_n = ks_n.data();
			c.ks_count = &ks_count;
			ks_rec.resize(n_seqs);
			c.ks_rec = ks_rec.data();
			c.ks_phase = 1;
			kt_for_nt(se_post_worker, &c, n_seqs, n_threads);
			ks_place_twins(o, bwt, bns, ks_rec);
			c.ks_phase = 2;
		}
		fb.reset();
		kt_for_nt(se_post_worker, &c, n_seqs, n_threads);
		int stop;
		if (fb.write(out, &stop) != 0) ret = 1;
		bwa_free_read_seq(n_seqs, seqs);
		for (bwt_aln1_t *a : ks_aln) free(a);
		if (ret) { // the batch read ahead is never aligned
			ahead.abandon();
			if (next.seqs) bwa_free_read_seq(next.n, next.seqs);
			break;
		}
		ahead.wait();
		if (bwa_verbose >= 3) fprintf(stderr, "[M::aln] %lld sequences have been processed.\n", tot);
		if (bwa_seq_report_read_error(ks)) { // where stock's read of the next batch would have exited
			ret = 1;
			break;
		}
		cur = next;
	}
	if (o->ks_windows && bwa_verbose >= 3) ks_report(ks_count, "reads");
	if (bwa_seq_n_bad(ks)) ret = 1; // malformed records were skipped (and named on stderr)
	// Not fatal (stock's err_gzclose exits): the caller still closes the output, so a BAM is
	// flushed and terminated before the run exits 1.
	if (bwa_seq_close(ks) != 0) ret = 1;
	return ret;
}

// The paired-end batch loop shared by plain and --meth runs: it reads each pair of batches (the
// next one while this one is aligned), reconciles their counts, writes the records and stops at
// mismatched names; `process` aligns one pair of batches into fb.
template <class Process>
static int run_pe_batches(const aln_opt_t *o, const char *fn_fq1, const char *fn_fq2,
						  aln_out_t *out, Process process)
{
	if (o->batch_size <= 0) { // as aln_run_se
		fprintf(stderr, "[E::aln] the batch size must be positive, not %d\n", o->batch_size);
		return 1;
	}
	bwa_seqio_t *ks[2];
	ks[0] = bwa_seq_open(fn_fq1);
	if (!ks[0]) {
		return 1;
	}
	ks[1] = bwa_seq_open(fn_fq2);
	if (!ks[1]) {
		bwa_seq_close(ks[0]);
		return 1;
	}
	bwase_initialize();
	const int n_threads = o->gap.n_threads < 1 ? 1 : o->gap.n_threads;
	fmt_bufs_t fb(n_threads);
	isize_info_t last_ii;
	last_ii.avg = -1.0;
	long long tot = 0;
	int ret = 0;
	bwa_seq_t *seqs[2];
	int n_in[2];
	bool r1_done = false;
	// A batch pair: R1's batch, then R2's, read only if R1's is non-empty (as stock reads them: a
	// read error or the end of R1 leaves R2 unread).
	struct pe_batch_t { bwa_seq_t *seqs[2]; int n[2]; };
	auto read_batch = [&](pe_batch_t *b) {
		b->seqs[1] = 0, b->n[1] = 0;
		if ((b->seqs[0] = bwa_read_seq(ks[0], o->batch_size, &b->n[0], o->gap.mode, o->gap.trim_qual)) != 0)
			b->seqs[1] = bwa_read_seq(ks[1], o->batch_size, &b->n[1], o->gap.mode, o->gap.trim_qual);
	};
	pe_batch_t cur = {}, next = {};
	bwa_seq_defer_read_errors(ks[0]); // reported below, on this thread (see background_task_t)
	bwa_seq_defer_read_errors(ks[1]);
	read_batch(&cur);
	if (bwa_seq_report_read_error(ks[0]) || bwa_seq_report_read_error(ks[1])) { // if R2's read failed, R1's batch was read
		for (int j = 0; j < 2; ++j)
			if (cur.seqs[j]) bwa_free_read_seq(cur.n[j], cur.seqs[j]);
		ret = 1;
	}
	background_task_t ahead;
	while (!ret) {
		seqs[0] = cur.seqs[0], seqs[1] = cur.seqs[1], n_in[0] = cur.n[0], n_in[1] = cur.n[1];
		if (seqs[0] == 0) {
			r1_done = true;
			break;
		}
		if (!seqs[1]) {
			// n_in[0] is only the R1 batch just read; later R1 batches are not read, so the
			// count is a lower bound
			fprintf(stderr, "[W::aln] read 2 file ended before read 1 file; ignoring the rest of read 1 (%d or more read(s))\n", n_in[0]);
			bwa_free_read_seq(n_in[0], seqs[0]);
			break;
		}
		// the next pair of batches is read while this one is aligned; started only past the end
		// checks above, so after R2 ends no further R1 batch is read, as before
		next = {}; // stays empty if the read is abandoned before it runs
		ahead.start([&] { read_batch(&next); });
		int n_seqs = n_in[1]; // as sampe: the second file's count
		if (n_in[0] != n_in[1]) {
			if (n_seqs > n_in[0]) n_seqs = n_in[0]; // sampe would read past the first file's reads
			fprintf(stderr, "[W::aln] read counts differ in this batch: R1 %d, R2 %d; using %d%s\n",
					n_in[0], n_in[1], n_seqs, n_seqs == n_in[1]? " (as bwa sampe)" : " (the smaller count; bwa sampe would read past R1)");
		}
		tot += n_seqs;
		process(n_seqs, n_in, seqs, fb, last_ii, tot - n_seqs);
		int stop;
		if (fb.write(out, &stop) != 0) ret = 1;
		if (stop >= 0) {
			fprintf(stderr, "[E::aln] paired reads have different names: \"%s\", \"%s\"\n", seqs[0][stop].name, seqs[1][stop].name);
			ret = 1;
		}
		bwa_free_read_seq(n_in[0], seqs[0]);
		bwa_free_read_seq(n_in[1], seqs[1]);
		if (ret) { // the batches read ahead are never aligned
			ahead.abandon();
			for (int j = 0; j < 2; ++j)
				if (next.seqs[j]) bwa_free_read_seq(next.n[j], next.seqs[j]);
			break;
		}
		ahead.wait();
		if (bwa_verbose >= 3) fprintf(stderr, "[M::aln] %lld read pairs have been processed.\n", tot);
		if (bwa_seq_report_read_error(ks[0]) || bwa_seq_report_read_error(ks[1])) { // where stock's read would have exited
			for (int j = 0; j < 2; ++j)
				if (next.seqs[j]) bwa_free_read_seq(next.n[j], next.seqs[j]);
			ret = 1;
			break;
		}
		cur = next;
	}
	if (r1_done && !ret) { // R1 ran out first: say so if R2 still has reads
		int n_extra = 0;
		bwa_seq_t *extra = bwa_read_seq(ks[1], 1, &n_extra, o->gap.mode, o->gap.trim_qual);
		if (bwa_seq_report_read_error(ks[1])) ret = 1;
		else if (extra) {
			fprintf(stderr, "[W::aln] read 1 file ended before read 2 file; ignoring the remaining read(s) of read 2\n");
			bwa_free_read_seq(n_extra, extra);
		}
	}
	if (bwa_seq_n_bad(ks[0]) || bwa_seq_n_bad(ks[1])) ret = 1; // malformed records, named on stderr
	if (bwa_seq_close(ks[0]) != 0) ret = 1;
	if (bwa_seq_close(ks[1]) != 0) ret = 1;
	return ret;
}

int aln_run_pe(const aln_opt_t *o, const bwt_t *bwt, const bntseq_t *bns, ubyte_t *pac,
			   const char *fn_fq1, const char *fn_fq2, aln_out_t *out)
{
	srand48(bns->seed); // once per run, never per batch
	// each file's own duplicate check: one mate's batches need not look like the other's
	aln_dup_state_t dup_state[2];
	aln_dup_state_t *dup[2] = {nullptr, nullptr};
	if (o->share_dup_search) dup[0] = &dup_state[0], dup[1] = &dup_state[1];
	const int n_threads = o->gap.n_threads < 1 ? 1 : o->gap.n_threads;
	std::vector<int64_t> ks_count((size_t)n_threads * KS_N_OUTCOMES, 0);
	ks_rescue_ctx_t ks_rctx = {o->known_snvs, bns};
	const pe_alt_sites_t ks_rescue = {ks_rescue_sites, &ks_rctx};
	const int ret = run_pe_batches(o, fn_fq1, fn_fq2, out,
		[&](int n_seqs, const int *n_in, bwa_seq_t **seqs, fmt_bufs_t &fb,
			isize_info_t &last_ii, int64_t) {
		// each file is its own batch: the gap-open cap follows that file's longest read
		std::vector<bwt_aln1_t*> ks_aln[2];
		std::vector<int> ks_n[2];
		for (int m = 0; m < 2; ++m) {
			const int n_shared = aln_search_batch(bwt, n_in[m], seqs[m], &o->gap, o->gap.n_threads, dup[m]);
			if (bwa_verbose >= 4)
				fprintf(stderr, "[M::aln] read %d: %d of %d reads shared an identical read's search.\n", m + 1, n_shared, n_in[m]);
			if (!o->ks_windows) continue;
			// --known-snvs: the same search on the ALT windows, kept beside the main search's hits
			ks_aln[m].resize(n_in[m]);
			ks_n[m].resize(n_in[m]);
			for (int i = 0; i < n_in[m]; ++i) {
				ks_aln[m][i] = seqs[m][i].aln; ks_n[m][i] = seqs[m][i].n_aln;
				seqs[m][i].aln = 0; seqs[m][i].n_aln = 0;
			}
			aln_search_batch(o->ks_windows->bwt, n_in[m], seqs[m], &o->gap, o->gap.n_threads, nullptr);
			for (int i = 0; i < n_in[m]; ++i) {
				std::swap(ks_aln[m][i], seqs[m][i].aln);
				std::swap(ks_n[m][i], seqs[m][i].n_aln);
			}
		}
		// bwa_cal_pac_pos_pe, stage by stage: the drand48 pass serially in stock order, then
		// positions in parallel, the insert-size barrier, then the per-pair work in parallel
		bwa_pe_aln2seq(n_seqs, seqs);
		pe_ctx_t c;
		c.o = o, c.bwt = bwt, c.bns = bns, c.pac = pac, c.seqs[0] = seqs[0], c.seqs[1] = seqs[1];
		c.fb = &fb;
		c.cnt_chg.assign(n_threads, 0);
		c.sw_cnt.assign(4 * n_threads, 0);
		kt_for_nt(pe_pos_worker, &c, n_seqs, n_threads);
		isize_info_t ii;
		// The insert-size and rescue reports follow -v, as on the --meth path (3 by default,
		// where they match stock sampe's stderr).
		bwa_pe_isize_ref(bwt->seq_len / 2, n_seqs, seqs, &ii, &o->pe, &last_ii, nullptr, bwa_verbose >= 3);
		c.ii = &ii;
		std::vector<std::vector<pe_hit_t>> ks_hits[2];
		std::vector<std::vector<int64_t>> ks_halt[2];
		std::vector<int64_t> ks_alt[2];
		std::vector<ks_rec_t> ks_rec[2];
		if (o->ks_windows) {
			for (int j = 0; j < 2; ++j) {
				ks_hits[j].resize(n_seqs);
				ks_halt[j].resize(n_seqs);
				ks_alt[j].assign(n_seqs, -1);
				ks_rec[j].resize(n_seqs);
				c.ks_aln[j] = ks_aln[j].data();
				c.ks_n[j] = ks_n[j].data();
				c.ks_hits[j] = ks_hits[j].data();
				c.ks_halt[j] = ks_halt[j].data();
				c.ks_alt[j] = ks_alt[j].data();
				c.ks_rec[j] = ks_rec[j].data();
			}
			c.ks_count = &ks_count;
			c.ks_rescue = &ks_rescue;
			c.ks_phase = 1;
			kt_for_nt(pe_post_worker, &c, n_seqs, n_threads);
			ks_place_twins_pe(o, bwt, bns, pac, &ii, seqs, ks_rec, &ks_rescue);
			c.ks_phase = 2;
		}
		fb.reset();
		kt_for_nt(pe_post_worker, &c, n_seqs, n_threads);
		for (int j = 0; j < 2; ++j)
			for (bwt_aln1_t *a : ks_aln[j]) free(a);
		long long cnt_chg = 0;
		uint64_t n_tot[2] = {0, 0}, n_mapped[2] = {0, 0};
		for (int t = 0; t < n_threads; ++t) {
			cnt_chg += c.cnt_chg[t];
			n_tot[0] += c.sw_cnt[4 * t], n_tot[1] += c.sw_cnt[4 * t + 1];
			n_mapped[0] += c.sw_cnt[4 * t + 2], n_mapped[1] += c.sw_cnt[4 * t + 3];
		}
		if (bwa_verbose >= 3) fprintf(stderr, "[M::aln] changed the coordinates of %lld alignments.\n", cnt_chg);
		if (bwa_verbose >= 3 && o->pe.is_sw && ii.avg >= 0.0) bwa_paired_sw_report(n_tot, n_mapped);
		last_ii = ii;
	});
	if (o->ks_windows && bwa_verbose >= 3) ks_report(ks_count, "ends");
	return ret;
}

namespace {

struct se_meth_ctx_t {
	const aln_opt_t *o;
	bwa_seq_t *seqs;
	aln_meth_read_t *res;
	fmt_bufs_t *fb;
	ks_rec_t *ks_rec = nullptr; // --known-snvs (as se_ctx_t)
	int ks_phase = 0;
};

/* --known-snvs under --meth: the twins of a batch searched (window index, then main with the
 * windows' candidates) and finished as single-end reads; their draws are seeded from ordinals past
 * any read's, so no read's draw changes. */
void ks_place_twins_meth(const aln_opt_t *o, const bwt_t *bwt, std::vector<ks_rec_t> &rec, int64_t ord0)
{
	std::vector<int> idx;
	for (int i = 0; i < (int)rec.size(); ++i)
		if (rec[i].twin_of >= 0) idx.push_back(i);
	const int nt = (int)idx.size();
	if (nt == 0) return;
	const int is_comp = o->gap.mode & BWA_MODE_COMPREAD;
	bwa_seq_t *tw = (bwa_seq_t*)aln_calloc(nt, sizeof(bwa_seq_t));
	for (int j = 0; j < nt; ++j) {
		const std::vector<uint8_t> &t = rec[idx[j]].tseq;
		bwa_seq_t *p = tw + j;
		p->tid = -1;
		p->full_len = p->clip_len = p->len = (uint32_t)t.size();
		p->seq = (ubyte_t*)aln_calloc(t.size(), 1);
		p->rseq = (ubyte_t*)aln_calloc(t.size(), 1);
		memcpy(p->seq, t.data(), t.size());
		memcpy(p->rseq, t.data(), t.size());
		seq_reverse(p->len, p->seq, 0);
		seq_reverse(p->len, p->rseq, is_comp);
		p->name = aln_strdup("twin");
	}
	const int n_threads = o->gap.n_threads < 1 ? 1 : o->gap.n_threads;
	const int64_t tord = ord0 + ((int64_t)1 << 48); // past any read's ordinal
	std::vector<aln_meth_read_t> wres((size_t)nt), tres((size_t)nt);
	aln_meth_search_batch(o->ks_windows->bwt, o->ks_windows->mref, o->meth, nt, tw, &o->gap, n_threads, 1, tord, 1,
						  o->max_xa_hits, nullptr, wres.data());
	std::vector<std::vector<aln_meth_candidate_t>> extra(nt);
	for (int j = 0; j < nt; ++j) ks_meth_candidates(o->ks_windows, o->meth_ref->bns, wres[j], &extra[j]);
	wres.clear();
	aln_meth_search_batch(bwt, o->meth_ref, o->meth, nt, tw, &o->gap, n_threads, 1, tord, 0, o->max_xa_hits, nullptr,
						  tres.data(), extra.data());
	aln_meth_finish_se(o->meth_ref, o->meth, &o->gap, o->max_xa_hits, nt, tw, tres.data());
	for (int j = 0; j < nt; ++j) {
		ks_rec_t &r = rec[idx[j]];
		const bwa_seq_t &t = tw[j];
		r.twin_mapq = t.type == BWA_TYPE_NO_MATCH ? 0 : t.mapQ;
		r.twin_same = t.type == BWA_TYPE_UNIQUE && (int)t.strand == r.strand &&
					  std::llabs((int64_t)t.pos - r.pos) <= KS_TWIN_SLOP;
	}
	bwa_free_read_seq(nt, tw);
}

/* Single end under --meth, reads [st, st + n): strata, main hit, refinement, MD/NM and SAM text. */
void se_meth_post_worker(void *data, int st, int n, int tid)
{
	se_meth_ctx_t *c = (se_meth_ctx_t*)data;
	const aln_opt_t *o = c->o;
	if (c->ks_phase != 2) {
		aln_meth_finish_se(o->meth_ref, o->meth, &o->gap, o->max_xa_hits, n, c->seqs + st, c->res + st);
		for (int i = 0; c->ks_rec && i < n; ++i) {
			const aln_meth_read_t &res = c->res[st + i];
			if (res.main < 0) continue;
			ks_rec_t &r = c->ks_rec[st + i];
			const aln_meth_hit_t &h = res.hits[res.main];
			r.twin_of = ks_site_votes(o->known_snvs, o->meth_ref->bns, &o->gap, c->seqs + st + i,
									  h.alt ? KS_ALT_ALL : -1, res.hyp, &r.votes);
			if (r.twin_of < 0) continue;
			ks_twin_seq(c->seqs + st + i, r.votes[r.twin_of], res.hyp, &r.tseq);
			r.pos = (int64_t)c->seqs[st + i].pos;
			r.strand = c->seqs[st + i].strand;
		}
		if (c->ks_phase == 1) return;
	}
	kstring_t *k = &c->fb->buf[tid];
	fmt_seg_t g;
	g.st = st, g.tid = tid, g.stop = -1, g.off = k->l;
	std::string xv;
	for (int i = 0; i < n; ++i) {
		aln_meth_print_sam1(k, o->meth_ref, o->meth, c->seqs + st + i, c->res + st + i, 1, 0, o->gap.mode, o->gap.max_top2);
		if (c->ks_rec && !c->ks_rec[st + i].votes.empty()) {
			const ks_rec_t &r = c->ks_rec[st + i];
			ks_format_xv(o->meth_ref->bns, r.votes, r.twin_of, r.twin_same, r.twin_mapq, &xv);
			--k->l;
			kputs("\tXV:Z:", k);
			kputs(xv.c_str(), k);
			kputc('\n', k);
		}
	}
	g.len = k->l - g.off;
	c->fb->segs[tid].push_back(g);
}

} // namespace

int aln_run_se_meth(const aln_opt_t *o, const bwt_t *bwt, const char *fn_fq, aln_out_t *out)
{
	if (o->batch_size <= 0) { // as aln_run_se
		fprintf(stderr, "[E::aln] the batch size must be positive, not %d\n", o->batch_size);
		return 1;
	}
	bwa_seqio_t *ks = bwa_seq_open(fn_fq);
	if (!ks) {
		return 1;
	}
	bwa_seq_defer_read_errors(ks);
	bwase_initialize();
	const int n_threads = o->gap.n_threads < 1 ? 1 : o->gap.n_threads;
	fmt_bufs_t fb(n_threads);
	long long tot = 0;
	int ret = 0;
	struct batch_t { bwa_seq_t *seqs; int n; };
	auto read_batch = [&](batch_t *b) { b->seqs = bwa_read_seq(ks, o->batch_size, &b->n, o->gap.mode, o->gap.trim_qual); };
	batch_t cur = {}, next = {};
	read_batch(&cur);
	if (bwa_seq_report_read_error(ks)) ret = 1;
	background_task_t ahead;
	aln_dup_state_t dup_state; // this run's duplicate check (keyed on the original bases)
	aln_dup_state_t *dup = o->share_dup_search? &dup_state : nullptr;
	std::vector<aln_meth_read_t> res;
	while (cur.seqs) {
		next = {}; // stays empty if the read is abandoned before it runs
		ahead.start([&] { read_batch(&next); }); // the next batch is read while this one is aligned
		bwa_seq_t *seqs = cur.seqs;
		const int n_seqs = cur.n;
		res.assign((size_t)n_seqs, aln_meth_read_t());
		// --known-snvs: search the ALT windows first; the windows they find covering their site are
		// resolved by each read's main search as candidates it may not have proposed itself
		std::vector<std::vector<aln_meth_candidate_t>> ks_extra;
		if (o->ks_windows) {
			std::vector<aln_meth_read_t> wres((size_t)n_seqs);
			aln_meth_search_batch(o->ks_windows->bwt, o->ks_windows->mref, o->meth, n_seqs, seqs, &o->gap,
								  o->gap.n_threads, 1, tot, 1, o->max_xa_hits, nullptr, wres.data());
			ks_extra.resize(n_seqs);
			for (int i = 0; i < n_seqs; ++i)
				ks_meth_candidates(o->ks_windows, o->meth_ref->bns, wres[i], &ks_extra[i]);
		}
		const int n_shared = aln_meth_search_batch(bwt, o->meth_ref, o->meth, n_seqs, seqs, &o->gap, o->gap.n_threads,
												   1, tot, 0, o->max_xa_hits, dup, res.data(),
												   ks_extra.empty() ? nullptr : ks_extra.data());
		if (bwa_verbose >= 4) fprintf(stderr, "[M::aln] %d of %d reads shared an identical read's search.\n", n_shared, n_seqs);
		tot += n_seqs;
		se_meth_ctx_t c = {o, seqs, res.data(), &fb};
		std::vector<ks_rec_t> ks_rec;
		if (o->ks_windows) {
			ks_rec.resize(n_seqs);
			c.ks_rec = ks_rec.data();
			c.ks_phase = 1;
			kt_for_nt(se_meth_post_worker, &c, n_seqs, n_threads);
			ks_place_twins_meth(o, bwt, ks_rec, tot);
			c.ks_phase = 2;
		}
		fb.reset();
		kt_for_nt(se_meth_post_worker, &c, n_seqs, n_threads);
		int stop;
		if (fb.write(out, &stop) != 0) ret = 1;
		bwa_free_read_seq(n_seqs, seqs);
		if (ret) { // the batch read ahead is never aligned
			ahead.abandon();
			if (next.seqs) bwa_free_read_seq(next.n, next.seqs);
			break;
		}
		ahead.wait();
		if (bwa_verbose >= 3) fprintf(stderr, "[M::aln] %lld sequences have been processed.\n", tot);
		if (bwa_seq_report_read_error(ks)) {
			ret = 1;
			break;
		}
		cur = next;
	}
	if (bwa_seq_n_bad(ks)) ret = 1;
	if (bwa_seq_close(ks) != 0) ret = 1;
	return ret;
}

namespace {

struct pe_meth_ctx_t {
	const aln_opt_t *o;
	bwa_seq_t *seqs[2];
	aln_meth_read_t *res[2];
	const isize_info_t *ii;
	fmt_bufs_t *fb;
	std::vector<long long> cnt_chg;
	std::vector<uint64_t> sw_cnt;
	ks_rec_t *ks_rec[2] = {nullptr, nullptr}; // --known-snvs (as se_ctx_t)
	int ks_phase = 0;
};

/* --known-snvs under --meth: whether a finished end's placement assumes its sites' ALTs. */
int64_t ks_meth_alt(const bwa_seq_t *s, const aln_meth_read_t &r)
{
	if (s->type == BWA_TYPE_MATESW) return r.ks_alt ? KS_ALT_ALL : -1;
	return r.main >= 0 && r.hits[r.main].alt ? KS_ALT_ALL : -1;
}

void pe_meth_post_worker(void *data, int st, int n, int tid)
{
	pe_meth_ctx_t *c = (pe_meth_ctx_t*)data;
	const aln_opt_t *o = c->o;
	bwa_seq_t *s[2] = {c->seqs[0] + st, c->seqs[1] + st};
	aln_meth_read_t *r[2] = {c->res[0] + st, c->res[1] + st};
	if (c->ks_phase != 2) {
		uint64_t *sw = &c->sw_cnt[4 * tid];
		c->cnt_chg[tid] += aln_meth_pe_stage2(o->meth_ref, o->meth, &o->gap, &o->pe, c->ii, n, s, r, sw, sw + 2);
		if (c->ks_rec[0])
			for (int j = 0; j < 2; ++j)
				for (int i = 0; i < n; ++i) {
					if (s[j][i].type == BWA_TYPE_NO_MATCH) continue;
					ks_rec_t &kr = c->ks_rec[j][st + i];
					kr.twin_of = ks_site_votes(o->known_snvs, o->meth_ref->bns, &o->gap, s[j] + i,
											   ks_meth_alt(s[j] + i, r[j][i]), r[j][i].hyp, &kr.votes);
					if (kr.twin_of < 0) continue;
					ks_twin_seq(s[j] + i, kr.votes[kr.twin_of], r[j][i].hyp, &kr.tseq);
					kr.pos = (int64_t)s[j][i].pos;
					kr.strand = s[j][i].strand;
				}
		if (c->ks_phase == 1) return;
	}
	std::string xv;
	kstring_t *k = &c->fb->buf[tid];
	fmt_seg_t g;
	g.st = st, g.tid = tid, g.stop = -1, g.off = k->l;
	for (int i = 0; i < n; ++i) {
		bwa_seq_t *p0 = s[0] + i, *p1 = s[1] + i;
		join_pair_barcodes(p0, p1);
		aln_meth_print_sam1(k, o->meth_ref, o->meth, p0, r[0] + i, 1, p1, o->gap.mode, o->gap.max_top2);
		if (c->ks_rec[0]) ks_append_xv(k, o->meth_ref->bns, c->ks_rec[0][st + i], &xv);
		aln_meth_print_sam1(k, o->meth_ref, o->meth, p1, r[1] + i, 0, p0, o->gap.mode, o->gap.max_top2);
		if (c->ks_rec[1]) ks_append_xv(k, o->meth_ref->bns, c->ks_rec[1][st + i], &xv);
		if (strcmp(p0->name, p1->name) != 0) {
			g.stop = st + i;
			break;
		}
	}
	g.len = k->l - g.off;
	c->fb->segs[tid].push_back(g);
}

/* --known-snvs under --meth, paired end: as ks_place_twins_pe, through the --meth pipeline (window
 * search, main search with the windows it finds, stages 1 and 2 with the batch's insert size). The
 * twins' ordinals start past any read's, so no read's draw changes. */
struct ks_twin_pe_meth_ctx_t {
	const aln_opt_t *o;
	const isize_info_t *ii;
	bwa_seq_t *tw[2];
	aln_meth_read_t *res[2];
};

void ks_twin_pe_meth_worker(void *data, int st, int n, int tid)
{
	(void)tid;
	ks_twin_pe_meth_ctx_t *c = (ks_twin_pe_meth_ctx_t*)data;
	const aln_opt_t *o = c->o;
	bwa_seq_t *s[2] = {c->tw[0] + st, c->tw[1] + st};
	aln_meth_read_t *r[2] = {c->res[0] + st, c->res[1] + st};
	aln_meth_pe_stage1(&o->gap, n, s, r);
	uint64_t sw[4] = {0, 0, 0, 0};
	aln_meth_pe_stage2(o->meth_ref, o->meth, &o->gap, &o->pe, c->ii, n, s, r, sw, sw + 2);
}

void ks_place_twins_pe_meth(const aln_opt_t *o, const bwt_t *bwt, const isize_info_t *ii, bwa_seq_t *const seqs[2],
							aln_meth_read_t *const res[2], std::vector<ks_rec_t> rec[2], int64_t ord0)
{
	std::vector<std::pair<int, int>> idx; // (pair, end holding the twin)
	for (int i = 0; i < (int)rec[0].size(); ++i)
		for (int j = 0; j < 2; ++j)
			if (rec[j][i].twin_of >= 0) idx.push_back({i, j});
	const int nt = (int)idx.size();
	if (nt == 0) return;
	const int is_comp = o->gap.mode & BWA_MODE_COMPREAD;
	bwa_seq_t *tw[2];
	for (int j = 0; j < 2; ++j) tw[j] = (bwa_seq_t*)aln_calloc(nt, sizeof(bwa_seq_t));
	for (int t = 0; t < nt; ++t) {
		const int i = idx[t].first, j = idx[t].second;
		const std::vector<uint8_t> &ts = rec[j][i].tseq;
		ks_make_read(tw[j] + t, ts.data(), (int)ts.size(), is_comp);
		const bwa_seq_t *m = seqs[1 - j] + i;
		const ks_site_t *site = rec[j][i].votes[rec[j][i].twin_of].site;
		std::vector<uint8_t> ms;
		for (const ks_vote_t &v : rec[1 - j][i].votes)
			if (v.site == site && (v.allele == 'R' || v.allele == 'A')) ks_twin_seq(m, v, res[1 - j][i].hyp, &ms);
		if (ms.empty()) ms.assign(m->seq, m->seq + m->len);
		ks_make_read(tw[1 - j] + t, ms.data(), (int)ms.size(), is_comp);
		tw[1 - j][t].full_len = m->full_len;
	}
	const int n_threads = o->gap.n_threads < 1 ? 1 : o->gap.n_threads;
	const int64_t tord = ord0 + ((int64_t)1 << 48); // past any read's ordinal
	std::vector<aln_meth_read_t> tres[2];
	for (int j = 0; j < 2; ++j) {
		std::vector<aln_meth_read_t> wres((size_t)nt);
		aln_meth_search_batch(o->ks_windows->bwt, o->ks_windows->mref, o->meth, nt, tw[j], &o->gap, n_threads,
							  j == 0 ? 1 : 0, tord, 1, o->max_xa_hits, nullptr, wres.data());
		std::vector<std::vector<aln_meth_candidate_t>> extra(nt);
		for (int t = 0; t < nt; ++t) ks_meth_candidates(o->ks_windows, o->meth_ref->bns, wres[t], &extra[t]);
		wres.clear();
		tres[j].assign((size_t)nt, aln_meth_read_t());
		aln_meth_search_batch(bwt, o->meth_ref, o->meth, nt, tw[j], &o->gap, n_threads, j == 0 ? 1 : 0, tord, 1,
							  o->max_xa_hits, nullptr, tres[j].data(), extra.data());
	}
	ks_twin_pe_meth_ctx_t c;
	c.o = o, c.ii = ii;
	for (int j = 0; j < 2; ++j) c.tw[j] = tw[j], c.res[j] = tres[j].data();
	kt_for_nt(ks_twin_pe_meth_worker, &c, nt, n_threads);
	for (int t = 0; t < nt; ++t) {
		const int i = idx[t].first, j = idx[t].second;
		ks_rec_t &r = rec[j][i];
		const bwa_seq_t *q = tw[j] + t;
		r.twin_mapq = q->type == BWA_TYPE_NO_MATCH ? 0 : q->mapQ;
		r.twin_same = q->type != BWA_TYPE_NO_MATCH && (int)q->strand == r.strand
			&& std::llabs((int64_t)q->pos - r.pos) <= KS_TWIN_SLOP;
	}
	for (int j = 0; j < 2; ++j) bwa_free_read_seq(nt, tw[j]);
}

} // namespace

int aln_run_pe_meth(const aln_opt_t *o, const bwt_t *bwt, const char *fn_fq1, const char *fn_fq2, aln_out_t *out)
{
	std::vector<aln_meth_read_t> res[2];
	aln_dup_state_t dup_state[2];
	aln_dup_state_t *dup[2] = {nullptr, nullptr};
	if (o->share_dup_search) dup[0] = &dup_state[0], dup[1] = &dup_state[1];
	const int n_threads = o->gap.n_threads < 1 ? 1 : o->gap.n_threads;
	return run_pe_batches(o, fn_fq1, fn_fq2, out,
		[&](int n_seqs, const int *n_in, bwa_seq_t **seqs, fmt_bufs_t &fb,
			isize_info_t &last_ii, int64_t tot) {
		for (int m = 0; m < 2; ++m) { // each file is its own batch and stream (read 1 C->T, read 2 G->A)
			res[m].assign((size_t)n_in[m], aln_meth_read_t());
			// --known-snvs: the ALT windows' search first, its windows resolved by the main search
			std::vector<std::vector<aln_meth_candidate_t>> ks_extra;
			if (o->ks_windows) {
				std::vector<aln_meth_read_t> wres((size_t)n_in[m]);
				aln_meth_search_batch(o->ks_windows->bwt, o->ks_windows->mref, o->meth, n_in[m], seqs[m], &o->gap,
									  o->gap.n_threads, m == 0 ? 1 : 0, tot, 1, o->max_xa_hits, nullptr, wres.data());
				ks_extra.resize(n_in[m]);
				for (int i = 0; i < n_in[m]; ++i) ks_meth_candidates(o->ks_windows, o->meth_ref->bns, wres[i], &ks_extra[i]);
			}
			const int n_shared = aln_meth_search_batch(bwt, o->meth_ref, o->meth, n_in[m], seqs[m], &o->gap, o->gap.n_threads,
													   m == 0 ? 1 : 0, tot, 1, o->max_xa_hits, dup[m], res[m].data(),
													   ks_extra.empty() ? nullptr : ks_extra.data());
			if (bwa_verbose >= 4)
				fprintf(stderr, "[M::aln] read %d: %d of %d reads shared an identical read's search.\n", m + 1, n_shared, n_in[m]);
		}
		aln_meth_read_t *r[2] = {res[0].data(), res[1].data()};
		aln_meth_pe_stage1(&o->gap, n_seqs, seqs, r);
		isize_info_t ii;
		// --known-snvs: a pair whose either end placed with an ALT credited is left out of the
		// estimate (its MAPQ 0 for the estimate only), which then rests on stock placements
		std::vector<std::pair<int, int>> ks_mapq; // (pair, its ends' MAPQs packed)
		if (o->ks_windows)
			for (int i = 0; i < n_seqs; ++i)
				if (ks_meth_alt(seqs[0] + i, r[0][i]) != -1 || ks_meth_alt(seqs[1] + i, r[1][i]) != -1) {
					ks_mapq.push_back({i, (int)(seqs[0][i].mapQ << 8 | seqs[1][i].mapQ)});
					seqs[0][i].mapQ = seqs[1][i].mapQ = 0;
				}
		// pairing is same-contig only, so the estimate counts only same-contig pairs
		bwa_pe_isize_ref(o->meth_ref->bns->l_pac, n_seqs, seqs, &ii, &o->pe, &last_ii, o->meth_ref->bns, bwa_verbose >= 3);
		for (const auto &x : ks_mapq) seqs[0][x.first].mapQ = x.second >> 8, seqs[1][x.first].mapQ = x.second & 0xff;
		pe_meth_ctx_t c;
		c.o = o, c.seqs[0] = seqs[0], c.seqs[1] = seqs[1], c.res[0] = r[0], c.res[1] = r[1], c.ii = &ii, c.fb = &fb;
		c.cnt_chg.assign(n_threads, 0);
		c.sw_cnt.assign(4 * n_threads, 0);
		std::vector<ks_rec_t> ks_rec[2];
		if (o->ks_windows) {
			for (int j = 0; j < 2; ++j) {
				ks_rec[j].resize(n_seqs);
				c.ks_rec[j] = ks_rec[j].data();
			}
			c.ks_phase = 1;
			kt_for_nt(pe_meth_post_worker, &c, n_seqs, n_threads);
			ks_place_twins_pe_meth(o, bwt, &ii, seqs, r, ks_rec, tot);
			c.ks_phase = 2;
		}
		fb.reset();
		kt_for_nt(pe_meth_post_worker, &c, n_seqs, n_threads);
		long long cnt_chg = 0;
		uint64_t n_tot[2] = {0, 0}, n_mapped[2] = {0, 0};
		for (int t = 0; t < n_threads; ++t) {
			cnt_chg += c.cnt_chg[t];
			n_tot[0] += c.sw_cnt[4 * t], n_tot[1] += c.sw_cnt[4 * t + 1];
			n_mapped[0] += c.sw_cnt[4 * t + 2], n_mapped[1] += c.sw_cnt[4 * t + 3];
		}
		if (bwa_verbose >= 3) fprintf(stderr, "[M::aln] changed the coordinates of %lld alignments.\n", cnt_chg);
		if (bwa_verbose >= 3 && o->pe.is_sw && ii.avg >= 0.0) bwa_paired_sw_report(n_tot, n_mapped);
		last_ii = ii;
	});
}

} // namespace bwa_aln

using namespace bwa_aln;

static int aln_usage(const aln_opt_t *o, FILE *fp)
{
	const gap_opt_t *g = &o->gap;
	fprintf(fp, "\n");
	fprintf(fp, "Usage:   bwa-mem3 aln [options] <prefix> <in.fq> [<in2.fq>]\n\n");
	fprintf(fp, "Aligns short reads with bwa's backtracking algorithm and writes SAM (or BAM) in one\n");
	fprintf(fp, "step: bwa aln followed by bwa samse (one FASTQ/FASTA, or '-' for stdin) or bwa sampe\n");
	fprintf(fp, "(two files), with no .sai files. On the tested inputs its records match bwa 0.7.19's\n");
	fprintf(fp, "(@PG differs) unless --meth, --short-insert, --ref-span-insert or --exact-rescue-window\n");
	fprintf(fp, "is given; `-n 1e-2` is a probability here and the integer 1 in stock.\n");
	fprintf(fp, "<prefix> is an index made by `bwa-mem3 index` (also one made with --meth, or one staged\n");
	fprintf(fp, "by `bwa-mem3 shm`).\n");
	fprintf(fp, "Each option shows its default and the stock command it comes from.\n\n");
	fprintf(fp, "Search options (from bwa aln):\n");
	fprintf(fp, "  -n NUM    max #diff (int) or missing prob under %.2f err rate (float) [%.2f]\n", BWA_AVG_ERR, g->fnr);
	fprintf(fp, "  -o INT    maximum number of gap opens [%d]\n", g->max_gapo);
	fprintf(fp, "  -e INT    maximum number of gap extensions, -1 for disabling long gaps [-1]\n");
	fprintf(fp, "  -i INT    do not put an indel within INT bp towards the ends [%d]\n", g->indel_end_skip);
	fprintf(fp, "  -d INT    maximum occurrences for extending a long deletion [%d]\n", g->max_del_occ);
	fprintf(fp, "  -l INT    seed length [%d]\n", g->seed_len);
	fprintf(fp, "  -k INT    maximum differences in the seed [%d]\n", g->max_seed_diff);
	fprintf(fp, "  -m INT    maximum entries in the queue [%d]\n", g->max_entries);
	fprintf(fp, "  -t INT    number of threads, at most %d; below 1 uses 1 [%d]\n", MAX_THREADS, g->n_threads);
	fprintf(fp, "  -M INT    mismatch penalty [%d]\n", g->s_mm);
	fprintf(fp, "  -O INT    gap open penalty [%d]\n", g->s_gapo);
	fprintf(fp, "  -E INT    gap extension penalty [%d]\n", g->s_gape);
	fprintf(fp, "  -R INT    stop searching when there are >INT equally best hits [%d]\n", g->max_top2);
	fprintf(fp, "  -q INT    quality threshold for read trimming down to %dbp [%d]\n", BWA_MIN_RDLEN, g->trim_qual);
	fprintf(fp, "  -B INT    length of barcode [0]\n");
	fprintf(fp, "  -L        add a log-scaled gap-open plus gap-extension count to the read-end gap exclusion (-i)\n");
	fprintf(fp, "  -N        non-iterative mode: search for all n-difference hits (slooow)\n");
	fprintf(fp, "  -I        the input is in the Illumina 1.3+ FASTQ-like format\n");
	fprintf(fp, "  -Y        filter Casava-filtered sequences\n");
	fprintf(fp, "\nOutput options (from bwa samse and sampe; they apply to single-end and paired-end output):\n");
	fprintf(fp, "  --max-xa-hits INT   maximum number of alternative hits in the XA tag [%d] (samse and sampe -n)\n", o->max_xa_hits);
	fprintf(fp, "  -r STR    read group header line such as '@RG\\tID:foo\\tSM:bar' [null] (samse and sampe -r)\n");
	fprintf(fp, "  -f FILE   file to write output to instead of stdout (samse and sampe -f)\n");
	fprintf(fp, "\nPaired-end options (two input files; from bwa sampe):\n");
	fprintf(fp, "  -a INT    maximum insert size [%d]\n", o->pe.max_isize);
	fprintf(fp, "  --max-occ INT   maximum occurrences for one end [%d] (sampe -o)\n", o->pe.max_occ);
	fprintf(fp, "  --max-xa-hits-discordant INT   maximum hits to output for discordant pairs [%d] (sampe -N)\n", o->pe.N_multi);
	fprintf(fp, "  -c FLOAT  prior of chimeric rate (lower bound) [%.1le]\n", o->pe.ap_prior);
	fprintf(fp, "  -s        disable Smith-Waterman for the unmapped mate\n");
	fprintf(fp, "  -A        disable insert size estimate (force -s)\n");
	fprintf(fp, "  --short-insert   estimate short inserts without a batch-wide read-length floor, pair mates\n"
	            "                   that start at the same base, and let mate rescue overlap the anchor [off]\n");
	fprintf(fp, "                   paired only, including --meth; -A overrides\n");
	fprintf(fp, "  --ref-span-insert\n"
	            "                   measure inserts from each alignment's reference span (read length\n"
	            "                   plus deletions minus insertions), not its read length, in insert-size\n"
	            "                   estimation, pairing (and -a) and mate rescue windows [off]\n");
	fprintf(fp, "                   paired only, including --meth\n");
	fprintf(fp, "  --exact-rescue-window\n"
	            "                   round the mate-rescue window outward from the insert range (avg +/-\n"
	            "                   3 std), instead of truncating its start and then its end as sampe\n"
	            "                   does, which drops up to two bases and moves the end of every window\n"
	            "                   whenever 6 std crosses an integer [off; on under --meth]\n");
	fprintf(fp, "                   paired only; --no-exact-rescue-window keeps sampe's window under --meth\n");
	fprintf(fp, "  -P        accepted and ignored (sampe -P loads the index into memory; always so here)\n");
	fprintf(fp, "\nMethylation (shared with bwa-mem3 mem):\n");
	meth_opts_usage_shared(fp);
	fprintf(fp, "   Under --meth, converted search proposes windows; original-reference global realignment\n");
	fprintf(fp, "   refines each (match +1, mismatch -3, gaps 5+length), then -M/-O/-E cost-ranks its CIGAR.\n");
	fprintf(fp, "   Strata and MAPQ use stock aln accounting; capped results make no exact strata claims.\n");
	fprintf(fp, "   TAPS/5base neutral breaks equal costs by fewest conversions. When this displaces an\n");
	fprintf(fp, "   equal-cost locus, MAPQ <= min(60,6*delta); delta is the smallest conversion-count gap.\n");
	fprintf(fp, "   Pair all retained loci up to --max-occ only without lookup/work cutoffs on either end.\n");
	fprintf(fp, "   Lookup-cap rescue uses <=255 resolved mate loci; rescued MAPQ stays 0, not proper-pair.\n");
	fprintf(fp, "   aln -B is barcode length, not mem's mismatch penalty. Records carry\n");
	fprintf(fp, "   XR/XG/XM:Z instead of XM:i/XO:i/XG:i.\n");
	fprintf(fp, "\nKnown SNVs:\n");
	fprintf(fp, "  --known-snvs FILE   VCF/BCF of SNVs (FILTER ignored): a read carrying a site's ALT is\n"
	            "                   placed as on a genome with both alleles, and every read over a site\n"
	            "                   gets XV:Z (allele, differences elsewhere, budget, allele-swapped twin's\n"
	            "                   placement) [off]; single-end, paired-end and --meth\n");
	fprintf(fp, "\nNew in bwa-mem3 aln:\n");
	fprintf(fp, "  -v INT    verbosity: 0-2 no progress reports, 3 progress, 4 diagnostics [3]\n");
	fprintf(fp, "            (some errors and warnings still print at 0)\n");
	fprintf(fp, "  --bam     write BAM instead of SAM\n");
	fprintf(fp, "  --bam-level INT     BAM compression level 0-9 [1]\n");
	fprintf(fp, "  --help    show this message\n");
	fprintf(fp, "\nThe long names stand in for stock samse/sampe -n and -o and sampe -N, whose letters\n");
	fprintf(fp, "are taken by the search options above.\n");
	fprintf(fp, "\n");
	return fp == stdout ? 0 : 1;
}

/* An input that stock would fail to open must fail before any output is created or the index
 * is loaded. "-" is standard input. */
static int aln_input_ok(const char *fn)
{
	struct stat st;
	memset(&st, 0, sizeof st);
	if (strcmp(fn, "-") == 0) return 1;
	if (stat(fn, &st) != 0) {
		fprintf(stderr, "[E::aln] cannot open the input '%s': %s\n", fn, strerror(errno));
		return 0;
	}
	if (S_ISDIR(st.st_mode)) {
		fprintf(stderr, "[E::aln] cannot open the input '%s': is a directory\n", fn);
		return 0;
	}
	if (access(fn, R_OK) != 0) {
		fprintf(stderr, "[E::aln] cannot open the input '%s': %s\n", fn, strerror(errno));
		return 0;
	}
	return 1;
}

/* Preflight existing output paths without opening them: a probe can create a file owned by
 * another process or give a FIFO reader an early EOF. Absent paths (including dangling symlink
 * targets) are left to aln_out_open, which reports failures after the index is loaded. */
static int aln_output_ok(const char *path)
{
	if (!path || strcmp(path, "-") == 0) return 1;
	struct stat st;
	const int existed = stat(path, &st) == 0;
	if (existed && S_ISDIR(st.st_mode)) {
		fprintf(stderr, "[E::aln] cannot open '%s' for writing: is a directory\n", path);
		return 0;
	}
	if (existed && access(path, W_OK) != 0) {
		fprintf(stderr, "[E::aln] cannot open '%s' for writing: %s\n", path, strerror(errno));
		return 0;
	}
	return 1;
}

/* 1 if -f names the regular file on standard input, which is read when an input is "-": opening
 * -f for writing would truncate it before it is read. */
static int aln_out_is_stdin(const char *out)
{
	struct stat so, si;
	if (!out || strcmp(out, "-") == 0) return 0;
	if (stat(out, &so) != 0 || fstat(STDIN_FILENO, &si) != 0 || !S_ISREG(si.st_mode)) return 0;
	return so.st_dev == si.st_dev && so.st_ino == si.st_ino;
}

/* 1 if `out` names an existing file that is also the input `in`: opening -f for writing would
 * truncate the input before it is read. "-" on either side is never the same file. */
static int aln_same_file(const char *out, const char *in)
{
	struct stat so, si;
	if (!out || strcmp(out, "-") == 0 || strcmp(in, "-") == 0) return 0;
	if (stat(out, &so) != 0 || stat(in, &si) != 0) return 0;
	return so.st_dev == si.st_dev && so.st_ino == si.st_ino;
}

/* With --bam, every @SQ line of an -r text (which replaces the index's @SQ lines, as in stock
 * samse) that names an index contig must give that contig's length: htslib would otherwise
 * write records past the declared end of the reference. A name the index lacks is left to the
 * record writer, which refuses a record on it. Returns 0 after naming the mismatch. */
static int aln_bam_sq_lengths_ok(const char *hdr, const bntseq_t *bns)
{
	std::vector<std::pair<std::string, int64_t>> idx;
	idx.reserve(bns->n_seqs);
	for (int i = 0; i < bns->n_seqs; ++i) idx.emplace_back(bns->anns[i].name, bns->anns[i].len);
	std::sort(idx.begin(), idx.end());
	for (const char *p = hdr; p && *p;) {
		const char *nl = strchr(p, '\n');
		const size_t n = nl ? (size_t)(nl - p) : strlen(p);
		if (n > 4 && strncmp(p, "@SQ\t", 4) == 0) {
			std::string sn, ln;
			for (const char *f = p + 4, *end = p + n; f < end;) {
				const char *t = (const char*)memchr(f, '\t', end - f);
				const char *fe = t ? t : end;
				if (fe - f > 3 && strncmp(f, "SN:", 3) == 0) sn.assign(f + 3, fe);
				else if (fe - f >= 3 && strncmp(f, "LN:", 3) == 0) ln.assign(f + 3, fe);
				f = fe + 1;
			}
			auto it = std::lower_bound(idx.begin(), idx.end(), std::make_pair(sn, (int64_t)INT64_MIN));
			if (it == idx.end() || it->first != sn) { p = nl ? nl + 1 : p + n; continue; }
			char *lend = NULL;
			errno = 0;
			const long long len = ln.empty() ? -1 : strtoll(ln.c_str(), &lend, 10);
			if (ln.empty() || errno != 0 || *lend != '\0' || len != (long long)it->second) {
				fprintf(stderr, "[E::aln] --bam: the -r header gives @SQ %s length %s, but the index contig is %lld bp\n",
						sn.c_str(), ln.empty() ? "(none)" : ln.c_str(), (long long)it->second);
				return 0;
			}
		}
		p = nl ? nl + 1 : p + n;
	}
	return 1;
}

/* Strict number parsing. Stock parses with atoi/atof, which turn "abc" into 0 and "20x" into 20;
 * here a malformed value, or one outside [lo, hi], is an error that names the option. */
static int aln_parse_int(const char *flag, const char *s, long lo, long hi, int *out)
{
	char *end = NULL;
	errno = 0;
	long v = strtol(s, &end, 10);
	if (*s == '\0' || *end != '\0' || errno == ERANGE || v < lo || v > hi) {
		if (lo == 0 && hi == INT_MAX)
			fprintf(stderr, "[E::aln] %s wants a non-negative integer, not '%s'\n", flag, s);
		else if (lo == INT_MIN && hi == INT_MAX)
			fprintf(stderr, "[E::aln] %s wants an integer, not '%s'\n", flag, s);
		else
			fprintf(stderr, "[E::aln] %s wants an integer from %ld to %ld, not '%s'\n", flag, lo, hi, s);
		return 0;
	}
	*out = (int)v;
	return 1;
}

static int aln_parse_double(const char *flag, const char *s, double lo, double *out)
{
	char *end = NULL;
	errno = 0;
	double v = strtod(s, &end);
	if (*s == '\0' || *end != '\0' || errno == ERANGE || !isfinite(v) || !(v >= lo)) {
		fprintf(stderr, "[E::aln] %s wants a finite number%s, not '%s'\n", flag, lo > -HUGE_VAL ? " of at least 0" : "", s);
		return 0;
	}
	*out = v;
	return 1;
}

enum { OPT_MAX_XA = 256, OPT_MAX_XA_DISC, OPT_MAX_OCC, OPT_BAM, OPT_BAM_LEVEL, OPT_HELP, OPT_SHORT_INSERT, OPT_REF_SPAN_INSERT, OPT_KNOWN_SNVS,
	   OPT_EXACT_RESCUE_WINDOW, OPT_NO_EXACT_RESCUE_WINDOW };

int main_aln(int argc, char *argv[])
{
	aln_opt_t o;
	aln_opt_init(&o);
	gap_opt_t *opt = &o.gap;
	int c, opte = -1, bam = 0, bam_level = 1, pe_flag = 0, bc_len = 0;
	int exact_rescue_window = -1; // -1: not given (on under --meth, off otherwise); else the last given form
	const char *out_path = NULL;
	const char *ks_path = NULL;
	ks_sites_t ks_sites;
	ks_windows_t ks_windows;
	std::unique_ptr<char, decltype(&free)> rg_line(NULL, &free);
	aln_meth_opt_t meth_opt;
	aln_meth_ref_t meth_ref;
	meth_cli_t meth_cli;
	memset(&meth_cli, 0, sizeof meth_cli);
	static const struct option long_opts[] = {
		{"max-xa-hits", required_argument, 0, OPT_MAX_XA},
		{"max-xa-hits-discordant", required_argument, 0, OPT_MAX_XA_DISC},
		{"max-occ", required_argument, 0, OPT_MAX_OCC},
		{"short-insert", no_argument, 0, OPT_SHORT_INSERT},
		{"ref-span-insert", no_argument, 0, OPT_REF_SPAN_INSERT},
		{"exact-rescue-window", no_argument, 0, OPT_EXACT_RESCUE_WINDOW},
		{"no-exact-rescue-window", no_argument, 0, OPT_NO_EXACT_RESCUE_WINDOW},
		{"known-snvs", required_argument, 0, OPT_KNOWN_SNVS},
		{"bam", no_argument, 0, OPT_BAM},
		{"bam-level", required_argument, 0, OPT_BAM_LEVEL},
		{"help", no_argument, 0, OPT_HELP},
		METH_LONG_OPTIONS_SHARED
		{0, 0, 0, 0}
	};

	// Negative values are refused where they are meaningless and stock crashes, aborts or drops
	// every read on them (-l, -m, -M, -O, -E, -B, an integer -n); other options take any integer, as
	// stock does (-e -1 is the documented "no long gaps", for example).
#define ALN_INT(flag, lo, dst) do { if (!aln_parse_int(flag, optarg, lo, INT_MAX, dst)) return 1; } while (0)
	bwa_rg_id[0] = '\0'; // bwa_set_rg's global: a run without -r must not tag reads with an earlier run's ID
	// A clean reparse (the unit tests call main_aln more than once). getopt keeps a pointer into the
	// previous argv when that parse ended on a flag such as -N, and the caller may have freed it.
#if defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || defined(__NetBSD__)
	optind = 1;
	optreset = 1; // BSD getopt keeps that pointer as a static `place`
#else
	optind = 0; // glibc and musl reinitialize on 0, dropping the stale pointer; optind=1 keeps it
#endif
	while ((c = getopt_long(argc, argv, "n:o:e:i:d:l:k:LR:m:t:v:NM:O:E:q:f:r:IYB:a:c:sAP", long_opts, 0)) >= 0) {
		switch (c) {
		case 'n':
			if (strpbrk(optarg, ".eE")) { // stock tests only for '.'; it reads `1e-2` as the integer 1
				double fnr;
				if (!aln_parse_double("-n", optarg, 0.0, &fnr)) return 1;
				if (!((float)fnr > 0.0f)) { // stock runs fnr <= 0 as -n -1 and leaves every read unmapped
					fprintf(stderr, "[E::aln] -n wants a probability above 0, not '%s'\n", optarg);
					return 1;
				}
				opt->fnr = fnr; // a float, as stock's `opt->fnr = atof(optarg)`
				opt->max_diff = -1;
			} else {
				ALN_INT("-n", 0, &opt->max_diff);
				opt->fnr = -1.0;
			}
			break;
		case 'o': ALN_INT("-o", INT_MIN, &opt->max_gapo); break;
		case 'e': ALN_INT("-e", INT_MIN, &opte); break;
		case 'M': ALN_INT("-M", 0, &opt->s_mm); break;
		case 'O': ALN_INT("-O", 0, &opt->s_gapo); break;
		case 'E': ALN_INT("-E", 0, &opt->s_gape); break;
		case 'd': ALN_INT("-d", INT_MIN, &opt->max_del_occ); break;
		case 'i': ALN_INT("-i", INT_MIN, &opt->indel_end_skip); break;
		case 'l': ALN_INT("-l", 0, &opt->seed_len); break;
		case 'k': ALN_INT("-k", INT_MIN, &opt->max_seed_diff); break;
		case 'm': ALN_INT("-m", 0, &opt->max_entries); break;
		case 't':
			if (!aln_parse_int("-t", optarg, INT_MIN, MAX_THREADS, &opt->n_threads)) return 1;
			break;
		case 'v': if (!aln_parse_int("-v", optarg, 0, 4, &bwa_verbose)) return 1; break;
		case 'L': opt->mode |= BWA_MODE_LOGGAP; break;
		case 'R': ALN_INT("-R", INT_MIN, &opt->max_top2); break;
		case 'q': ALN_INT("-q", INT_MIN, &opt->trim_qual); break;
		case 'N': opt->mode |= BWA_MODE_NONSTOP; opt->max_top2 = 0x7fffffff; break;
		case 'f': out_path = optarg; break;
		case 'r':
			// stock bwa_set_rg rejects a literal tab; bwa-mem3's bwa_set_rg does not
			if (strchr(optarg, '\t') != NULL) {
				fprintf(stderr, "[E::aln] the read group line contained literal <tab> characters -- replace with escaped tabs: \\t\n");
				return 1;
			}
			rg_line.reset(bwa_set_rg(optarg));
			if (!rg_line) return 1;
			{ // Escaping can introduce interior blank lines; the writer adds the final newline.
				const size_t rl = strlen(rg_line.get());
				if ((rl && rg_line.get()[rl - 1] == '\n') || strstr(rg_line.get(), "\n\n")) {
					fprintf(stderr, "[E::aln] the read group text would write an empty SAM header line\n");
					return 1;
				}
			}
			break;
		case 'I': opt->mode |= BWA_MODE_IL13; break;
		case 'Y': opt->mode |= BWA_MODE_CFY; break;
		case 'B': {
			// A negative length would set every mode bit and make every read "shorter than its
			// barcode" (stock drops them all and exits 0). Repeated -B values are OR-ed into the
			// mode, as stock does; over BWA_MAX_BCLEN is refused below, where the limit for one
			// or two read files is known (the shift is skipped then so it stays defined).
			int b;
			ALN_INT("-B", 0, &b);
			if (b > bc_len) bc_len = b;
			if (b <= BWA_MAX_BCLEN) opt->mode |= b << 24;
			break;
		}
		case 'a': ALN_INT("-a", INT_MIN, &o.pe.max_isize); pe_flag = 1; break;
		case 'c': if (!aln_parse_double("-c", optarg, -HUGE_VAL, &o.pe.ap_prior)) return 1; pe_flag = 1; break;
		case 's': o.pe.is_sw = 0; pe_flag = 1; break;
		case 'A': o.pe.force_isize = 1; pe_flag = 1; break;
		case 'P': pe_flag = 1; break; // sampe -P preloads the index; it is always in memory here
		// Alternative-hit generation adds one to these limits using int arithmetic.
		case OPT_MAX_XA:
			if (!aln_parse_int("--max-xa-hits", optarg, INT_MIN, INT_MAX - 1, &o.max_xa_hits)) return 1;
			o.pe.n_multi = o.max_xa_hits;
			break;
		case OPT_MAX_XA_DISC:
			if (!aln_parse_int("--max-xa-hits-discordant", optarg, INT_MIN, INT_MAX - 1, &o.pe.N_multi)) return 1;
			pe_flag = 1;
			break;
		case OPT_MAX_OCC: ALN_INT("--max-occ", INT_MIN, &o.pe.max_occ); pe_flag = 1; break;
		case OPT_SHORT_INSERT: o.pe.short_insert = 1; break;
		case OPT_REF_SPAN_INSERT: o.pe.ref_span_insert = 1; break;
		case OPT_EXACT_RESCUE_WINDOW: exact_rescue_window = 1; break;
		case OPT_NO_EXACT_RESCUE_WINDOW: exact_rescue_window = 0; break;
		case OPT_KNOWN_SNVS: ks_path = optarg; break;
		case OPT_BAM: bam = 1; break;
		case OPT_BAM_LEVEL:
			if (!aln_parse_int("--bam-level", optarg, 0, 9, &bam_level)) return 1;
			break;
		case OPT_HELP: return aln_usage(&o, stdout);
		default:
			if (meth_opts_parse_shared(c, optarg, &meth_opt.m, &meth_cli)) break;
			return 1;
		}
	}
#undef ALN_INT
	// getopt binds optional arguments only with '='; do not reinterpret a real path, nor an index
	// prefix whose FASTA is gone (its plain or converted index files are there).
	for (int pos = optind; meth_cli.set_meth && pos < argc; ++pos) {
		const char *canon = nullptr;
		const char *flag = meth_opts_stray_flag(argv[pos], &canon);
		if (!flag || access(argv[pos], F_OK) == 0) continue;
		if (pos == optind && (access((std::string(argv[pos]) + ".ann").c_str(), F_OK) == 0 ||
							  access((std::string(argv[pos]) + ".meth.bwt.2bit.64").c_str(), F_OK) == 0))
			continue;
		fprintf(stderr, "ERROR: '%s' looks like a value for %s, not a positional path; "
				"use %s=%s (comma-separated for --meth-tags)\n",
				argv[pos], flag, flag, canon ? canon : argv[pos]);
		return 1;
	}
	meth_opts_finalize_shared(&meth_opt.m, &meth_cli); // the shared checks and chemistry defaults
	if (meth_opt.m.meth_mode) {
		aln_meth_opt_finish(&meth_opt);
		o.meth = &meth_opt;
		o.meth_ref = &meth_ref;
		if (aln_meth_header_check(rg_line.get()) != 0) return 1;
	} else if (meth_opt.m.meth_set_as_failed) {
		fprintf(stderr, "[E::aln] --set-as-failed requires --meth\n");
		return 1;
	}
	// --meth turns the exact rescue window on unless --no-exact-rescue-window is given; plain aln
	// keeps sampe's truncated window unless --exact-rescue-window is given.
	o.pe.exact_rescue_window = exact_rescue_window >= 0 ? exact_rescue_window : meth_opt.m.meth_mode ? 1 : 0;
	if (opte > 0) {
		opt->max_gape = opte;
		opt->mode &= ~BWA_MODE_GAPE;
	}
	if (opt->n_threads < 1) {
		fprintf(stderr, "[W::aln] -t %d is below 1; using 1 thread\n", opt->n_threads);
		opt->n_threads = 1;
	}
	{
		// The search keeps one priority stack per score below
		// score(max_diff + 1, max_gapo + 1, max_gape + 1) and indexes it by an entry's score. With
		// non-negative penalties every entry fits when one more difference of each kind still adds
		// a positive score: all three penalties 0, or a -o below -1 that cancels the others, leaves
		// too few stacks, and stock writes past the end of them.
		const long long room = opt->max_gapo >= 0
			? (long long)opt->s_mm + opt->s_gapo + opt->s_gape
			: (long long)opt->s_mm + (long long)(opt->max_gapo + 1) * opt->s_gapo + ((long long)opt->max_gape + 1) * opt->s_gape;
		if (room <= 0) {
			fprintf(stderr, "[E::aln] -M %d -O %d -E %d with -o %d leaves the search no room for a difference; "
					"use a positive penalty%s\n", opt->s_mm, opt->s_gapo, opt->s_gape, opt->max_gapo,
					opt->max_gapo < -1 ? " or -o -1 or more" : "");
			return 1;
		}
	}

	if (optind + 2 > argc) return aln_usage(&o, stderr);
	if (optind + 4 <= argc) return aln_usage(&o, stderr);
	const int paired = optind + 3 == argc;
	if (o.pe.short_insert && !paired) {
		fprintf(stderr, "[E::aln] --short-insert requires paired-end input\n");
		return 1;
	}
	if (o.pe.ref_span_insert && !paired) {
		fprintf(stderr, "[E::aln] --ref-span-insert requires paired-end input\n");
		return 1;
	}
	if (exact_rescue_window == 1 && !paired) { // the --meth default and --no-... are fine with one file
		fprintf(stderr, "[E::aln] --exact-rescue-window requires paired-end input\n");
		return 1;
	}
	if (pe_flag && !paired)
		fprintf(stderr, "[W::aln] paired-end options (-a, -c, -s, -A, -P, --max-occ, --max-xa-hits-discordant) are ignored with one FASTQ\n");

	{
		const int max_bc = paired ? BWA_MAX_BCLEN / 2 : BWA_MAX_BCLEN; // paired: both barcodes are joined into one buffer
		if (bc_len > max_bc) {
			fprintf(stderr, "[E::aln] -B %d is too long: the maximum barcode length is %d%s\n", bc_len, max_bc,
					paired ? " with two read files (the two barcodes are joined)" : "");
			return 1;
		}
	}

	// the positional arguments: index prefix, read 1, and read 2 when paired
	const char *prefix = argv[optind], *fn_fq = argv[optind + 1], *fn_fq2 = paired ? argv[optind + 2] : nullptr;
	// --meth: the search runs on the converted index <prefix>.meth; the original index at <prefix>
	// supplies the reference (both written by `bwa-mem3 index --meth`). <prefix>.meth is accepted too,
	// unless <prefix>.meth is itself a reference indexed with --meth (it has <prefix>.meth.meth).
	std::string meth_orig, meth_conv;
	if (o.meth) {
		const std::string a = prefix, sfx = ".meth";
		const bool ends = a.size() >= sfx.size() && a.compare(a.size() - sfx.size(), sfx.size(), sfx) == 0;
		const std::string stripped = ends ? a.substr(0, a.size() - sfx.size()) : std::string();
		// A converted index counts whether its FM-index is on disk or staged with `bwa-mem3 shm`
		// (attached without its disk files).
		auto conv_present = [](const std::string &p) {
			return access((p + ".bwt.2bit.64").c_str(), F_OK) == 0 || bwa_shm_test(p.c_str()) == 1;
		};
		const bool as_orig = conv_present(a + sfx);
		const bool as_conv = ends && access((stripped + ".ann").c_str(), F_OK) == 0 && conv_present(a);
		if (as_orig && as_conv) {
			fprintf(stderr, "ERROR: --meth index '%s' is ambiguous: it names either a reference indexed with --meth "
					"or the converted index of '%s'; rename one of them\n", prefix, stripped.c_str());
			return 1;
		}
		if (ends && !as_orig) {
			meth_conv = a;
			meth_orig = a.substr(0, a.size() - sfx.size());
		} else {
			meth_conv = a + sfx;
			meth_orig = a;
		}
		prefix = meth_conv.c_str();
		// Probe the small original-reference files before loading the converted FM-index.
		for (const char *suffix : {".ann", ".amb", ".pac"}) {
			const std::string path = meth_orig + suffix;
			struct stat st;
			if (stat(path.c_str(), &st) != 0 || !S_ISREG(st.st_mode) || access(path.c_str(), R_OK) != 0) {
				fprintf(stderr, "ERROR: --meth cannot read the original reference index file '%s'; "
						"run `bwa-mem3 index --meth` on the reference\n", path.c_str());
				return 1;
			}
		}
		// A converted index staged with `bwa-mem3 shm` is attached without its disk files.
		const bool conv_staged = bwa_shm_test(meth_conv.c_str()) == 1;
		for (const char *suffix : {".bwt.2bit.64", ".ann", ".amb", ".pac"}) {
			if (conv_staged) break;
			const std::string path = meth_conv + suffix;
			struct stat st;
			if (stat(path.c_str(), &st) != 0 || !S_ISREG(st.st_mode) || access(path.c_str(), R_OK) != 0) {
				fprintf(stderr, "ERROR: --meth cannot read the converted index file '%s'; "
						"run `bwa-mem3 index --meth` on the reference\n", path.c_str());
				return 1;
			}
		}
	}
	{
		// An index staged with `bwa-mem3 shm` is attached by load_index without reading the
		// disk files, so only an index that is not staged must be on disk.
		std::string probe = std::string(prefix) + ".bwt.2bit.64";
		if (bwa_shm_test(prefix) != 1 && access(probe.c_str(), F_OK) != 0) {
			if (o.meth)
				fprintf(stderr, "ERROR: --meth cannot locate the converted index '%s' (no %s); "
						"run `bwa-mem3 index --meth` on the reference\n", prefix, probe.c_str());
			else
				fprintf(stderr, "[E::aln] fail to locate the index '%s' (no %s); run `bwa-mem3 index` first\n",
						prefix, probe.c_str());
			return 1;
		}
	}
	if (paired && strcmp(fn_fq, "-") == 0 && strcmp(argv[optind + 2], "-") == 0) {
		// the first '-' would take all of standard input and leave the second empty
		fprintf(stderr, "[E::aln] standard input ('-') can stand for only one of the two read files\n");
		return 1;
	}
	if (!aln_input_ok(fn_fq) || (paired && !aln_input_ok(argv[optind + 2])))
		return 1;
	for (int k = 0; k < (paired ? 2 : 1); ++k) {
		const char *in = argv[optind + 1 + k];
		if (aln_same_file(out_path, in) || (strcmp(in, "-") == 0 && aln_out_is_stdin(out_path))) {
			fprintf(stderr, "[E::aln] -f '%s' is also an input; writing it would truncate the reads before they are read\n", out_path);
			return 1;
		}
	}
	{ // -f naming an index file would truncate it before (or while) the index is read; the
	  // `.meth.*` files are an `index --meth` converted index, which `aln --meth` and `mem --meth` read
		static const char *const idx_sfx[] = {"", ".bwt.2bit.64", ".pac", ".ann", ".amb", ".0123", ".alt",
											  ".meth.bwt.2bit.64", ".meth.pac", ".meth.ann", ".meth.amb"};
		for (const std::string &base : {std::string(prefix), meth_orig}) { // --meth also reads the original index
			if (base.empty()) continue;
			for (const char *sfx : idx_sfx) {
				std::string comp = base + sfx;
				if (aln_same_file(out_path, comp.c_str())) {
					fprintf(stderr, "[E::aln] -f '%s' is the index file '%s'; writing it would overwrite the index\n",
							out_path, comp.c_str());
					return 1;
				}
			}
		}
		if (!meth_orig.empty()) { // --meth reads the original reference's header sidecar
			const std::string hdr_path = meth_orig + ".hdr";
			char *dict_path = bwa_index_dict_path(meth_orig.c_str());
			if (dict_path == nullptr) {
				fprintf(stderr, "[E::aln] out of memory checking -f '%s' against the header sidecars\n", out_path);
				return 1;
			}
			const char *side = aln_same_file(out_path, hdr_path.c_str()) ? hdr_path.c_str()
							 : aln_same_file(out_path, dict_path) ? dict_path : nullptr;
			const int refused = side != nullptr;
			if (refused) fprintf(stderr, "[E::aln] -f '%s' is the header sidecar '%s'; writing it would overwrite "
								 "the reference's header\n", out_path, side);
			free(dict_path);
			if (refused) return 1;
		}
	}
	if (!aln_output_ok(out_path)) return 1;
	FMI_search *fmi = new FMI_search(prefix);
	// aln searches with bwt_2occ4/call_one_step and never reads mem's seed jump table, so
	// skip building it (saves its memory and build time at start-up; output unchanged).
	fmi->load_index(/*load_pac=*/!o.meth, opt->n_threads, /*build_seed_jump=*/false); // under --meth the converted pac is not used
	const bntseq_t *bns = fmi->idx->bns;
	bwa_aln::bwt_t *bwt = bwt_from_fmi(fmi);
	char *meth_hdr_lines = NULL;
	if (o.meth) {
		if (aln_meth_ref_load(&meth_ref, meth_orig.c_str(), bns, opt->n_threads) != 0) {
			bwt_destroy(bwt);
			delete fmi;
			return 1;
		}
		meth_hdr_lines = bwa_load_hdr_from_index(meth_orig.c_str()); // the original reference's sidecar
	}
	if (ks_path) {
		// Sites are checked against the original reference (under --meth, not the converted one).
		const bntseq_t *ks_bns = o.meth ? meth_ref.bns : bns;
		const uint8_t *ks_pac = o.meth ? meth_ref.pac : fmi->idx->pac;
		if (ks_load(ks_path, ks_bns, ks_pac, &ks_sites) != 0) {
			free(meth_hdr_lines);
			if (o.meth) aln_meth_ref_free(&meth_ref);
			bwt_destroy(bwt);
			delete fmi;
			return 1;
		}
		if (bwa_verbose >= 3)
			fprintf(stderr, "[M::aln] --known-snvs: %lld SNV alleles from %lld records (%lld alleles not SNVs, "
					"%lld duplicates skipped)\n", (long long)ks_sites.n_sites, (long long)ks_sites.n_records,
					(long long)ks_sites.n_not_snv, (long long)ks_sites.n_duplicate);
		if (!ks_sites.empty()) {
			if (o.meth) { // --meth realigns its candidate windows with the ALTs in place too
				meth_ref.known_snvs = &ks_sites;
				o.known_snvs = &ks_sites;
			}
			const double t0 = realtime();
			if (ks_windows_build(&ks_sites, ks_bns, ks_pac, opt->n_threads, o.meth != nullptr, &ks_windows) != 0) {
				free(meth_hdr_lines);
				if (o.meth) aln_meth_ref_free(&meth_ref);
				bwt_destroy(bwt);
				delete fmi;
				return 1;
			}
			if (bwa_verbose >= 3)
				fprintf(stderr, "[M::aln] --known-snvs: indexed %zu ALT windows in %.1f s\n",
						ks_windows.win.size(), realtime() - t0);
			o.known_snvs = &ks_sites;
			o.ks_windows = &ks_windows;
		}
	}
	kstring_t hdr = {0, 0, 0};
	if (o.meth) {
		if (aln_meth_header(&hdr, &meth_ref, rg_line.get(), bwa_pg, meth_hdr_lines) != 0) {
			free(hdr.s);
			free(meth_hdr_lines);
			aln_meth_ref_free(&meth_ref);
			bwt_destroy(bwt);
			delete fmi;
			return 1;
		}
	}
	else aln_print_sam_hdr(&hdr, bns, rg_line.get(), bwa_pg);
	free(meth_hdr_lines);
	aln_out_t *out = (!bam || aln_bam_sq_lengths_ok(rg_line.get(), bns))
		? aln_out_open(out_path, bam, bam_level, hdr.s, opt->n_threads) : NULL;
	free(hdr.s);
	int ret = 1;
	if (out) {
		if (o.meth) ret = paired ? aln_run_pe_meth(&o, bwt, fn_fq, fn_fq2, out)
								 : aln_run_se_meth(&o, bwt, fn_fq, out);
		else ret = paired ? aln_run_pe(&o, bwt, bns, fmi->idx->pac, fn_fq, fn_fq2, out)
						  : aln_run_se(&o, bwt, bns, fmi->idx->pac, fn_fq, out);
	}
	if (out && aln_out_close(out) != 0) ret = 1;

	ks_windows_free(&ks_windows);
	bwt_destroy(bwt);
	if (o.meth) aln_meth_ref_free(&meth_ref);
	delete fmi;
	return ret;
}
