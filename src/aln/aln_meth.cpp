/* SPDX-License-Identifier: MIT */
/* New code for bwa-mem3 (not from upstream bwa): `aln --meth`. See aln_meth.h and README.md. */
#include "aln_meth.h"
#include "aln_meth_realign.h"
#include "aln_meth_search.h"

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <algorithm>
#include <string>
#include <map>

#include "../bwa.h"
#include "../compat_target.h"
#include "../ksw.h"
#include "../kthread.h"
#include "../meth_bam.h"
#include "../meth_model.h"
#include "../meth_opts.h"
#include "../meth_ref.h"
#include "../meth_xm.h"
#include "../utils.h"
#include "aln_alloc.h"
#include "bwape.h"
#include "bwase.h"
#include "bwaseqio.h"
#include "bwtgap.h"

namespace bwa_aln {

/* ---------------------------------------------------------------------------------------------
 * Options, references, header, conversion
 * ------------------------------------------------------------------------------------------- */

void aln_meth_opt_finish(aln_meth_opt_t *mo)
{
	int8_t base[25];
	bwa_fill_scmat(1, 3, base); // aln's scale: stock's refinement and mate rescue use bwa_fill_scmat(1, 3)
	meth_fill_mats(base, 1, &mo->m, mo->mat_ot, mo->mat_ob, mo->mat7_ot, mo->mat7_ob);
}

int aln_meth_ref_load(aln_meth_ref_t *r, const char *orig_prefix, const bntseq_t *conv_bns, int n_threads)
{
	r->conv_bns = conv_bns;
	if (meth_orig_ref_load(orig_prefix, &r->bns, &r->pac, n_threads < 1 ? 1 : n_threads) != 0) return -1;
	// Pairing is structural: r<name>,f<name> in original order, with identical lengths and
	// contiguous offsets. Metadata cannot detect two different sequences with the same layout.
	bool matches = conv_bns && (int64_t)conv_bns->n_seqs == 2LL * r->bns->n_seqs
		&& conv_bns->l_pac == 2 * r->bns->l_pac;
	for (int i = 0; matches && i < r->bns->n_seqs; ++i) {
		const bntann1_t &orig = r->bns->anns[i];
		for (int h = 0; h < 2; ++h) {
			const bntann1_t &conv = conv_bns->anns[2 * i + h];
			const std::string name = std::string(h ? "f" : "r") + orig.name;
			if (name != conv.name || conv.len != orig.len || conv.offset != 2 * orig.offset + h * (int64_t)orig.len)
				matches = false;
		}
	}
	if (!matches) {
		fprintf(stderr, "ERROR: --meth converted index does not match the original reference '%s' "
				"(contig names, lengths or offsets); rebuild with `bwa-mem3 index --meth`\n", orig_prefix);
		aln_meth_ref_free(r);
		return -1;
	}
	return 0;
}

void aln_meth_ref_free(aln_meth_ref_t *r)
{
	meth_orig_ref_free(&r->bns, &r->pac);
}

int aln_meth_header_check(const char *rg_line)
{
	// The original reference supplies the authoritative sequence dictionary.
	for (const char *line = rg_line; line && *line;) {
		if (strncmp(line, "@SQ", 3) == 0 && (line[3] == '\t' || line[3] == '\n' || line[3] == '\0')) {
			fprintf(stderr, "ERROR: --meth does not accept @SQ records in -r; the original index supplies them\n");
			return -1;
		}
		const char *next = strchr(line, '\n');
		line = next ? next + 1 : nullptr;
	}
	return 0;
}

int aln_meth_header(kstring_t *out, const aln_meth_ref_t *ref, const char *rg_line, const char *pg_line,
					const char *orig_hdr_lines)
{
	if (aln_meth_header_check(rg_line) != 0) return -1;
	return meth_hdr_text(out, ref->bns, BWAMEM3_DEFAULT_HD_LINE, rg_line, orig_hdr_lines, pg_line,
						 NULL, "bwa-mem3 aln --meth");
}

void aln_meth_convert(int read_no, const ubyte_t *in, int len, ubyte_t *out)
{
	const ubyte_t from = read_no ? 1 : 2, to = read_no ? 3 : 0; // read 1: C->T, read 2: G->A
	for (int k = 0; k < len; ++k) out[k] = in[k] == from ? to : in[k];
}

/* Original-reference model helpers for tags on the retained alignment. */
static inline const int8_t *meth_mat(const aln_meth_opt_t *mo, int hyp, int *m)
{
    *m = mo->m.meth_context == METH_CTX_CPG ? METH_M : 5;
    return *m == 5 ? (hyp ? mo->mat_ot : mo->mat_ob) : (hyp ? mo->mat7_ot : mo->mat7_ob);
}
static void fetch_window(const aln_meth_ref_t *ref, const aln_meth_opt_t *mo, int64_t rb, int64_t len,
                         int rid, std::vector<uint8_t> &win)
{
    aln_meth_fetch_window(*ref, *mo, rb, (int)len, rid, win);
}
static inline uint64_t splitmix64(uint64_t *x)
{
    uint64_t z = (*x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}
uint64_t aln_meth_draw_seed(const aln_meth_ref_t *ref, int read_no, int64_t ord)
{
    uint64_t r = (uint64_t)ref->bns->seed * 0x9e3779b97f4a7c15ULL ^ ((uint64_t)ord << 1 | (uint64_t)read_no);
    splitmix64(&r);
    return r;
}
/* For current best B, retained cost T=B+M excludes gaps when T<O and E>=0.
 * Its raw CIGAR then has <=floor(T/M) strict mismatches. Every such equal-span
 * window has an ungapped converted proposal with no more mismatches: genuine
 * matches/conversions project identically and either N counts strictly regardless
 * of PAC fill. Unseeded discovery admitting that many differences proposes every
 * surviving window by score M*floor(T/M). Refine the entire stratum before stopping.
 * B can only decrease, so later flushes can only narrow the required stratum.
 * This bounds equivalent proposals, NOT the cost of arbitrary search histories. */
bool aln_meth_can_stop(const aln_meth_slot_t *x, int next_score)
{
    const gap_job_t *job=x->job;
    const gap_opt_t &opt=job->opt;
    if (!aln_meth_stop_eligible(x) || x->res->occ.empty()) return false;
    int64_t best=INT64_MAX;
    for (const auto &o : x->res->occ) best=std::min(best,o.cost);
    // Stored costs fit INT_MAX; widen before adding penalties or multiplying.
    if (best<0 || best>INT_MAX) return false;
    const int64_t retained=best+(int64_t)opt.s_mm;
    if (retained>=opt.s_gapo) return false;
    const int64_t differences=retained/opt.s_mm;
    return opt.max_diff>=differences && (int64_t)next_score>differences*opt.s_mm;
}

void aln_meth_hit(aln_meth_slot_t *x, const gap_entry_t *e, bwtint_t k, bwtint_t l)
{
    // Reaching the same interval with the same span adds no new candidate window.
    // Other spans must still be evaluated, even if their normalized starts later coincide.
    const int span = x->len + (int)e->n_del - (int)e->n_ins;
    const auto key = std::make_tuple(k, l, span);
    if (!x->intervals.insert(key).second) return;
    x->pend.push_back({k, l, *e});
}
const ubyte_t *aln_meth_slot_start(aln_meth_slot_t *x, bwa_seq_t *p, aln_meth_read_t *res, int64_t ord, int, int has_copy)
{
    x->p = p; x->res = res; x->len = p->len;
    x->keep = has_copy ? ALN_METH_KEEP_COPIES : x->grp->keep_strata ? ALN_METH_KEEP_STRATA : 0;
    *res = aln_meth_read_t();
    x->rng = aln_meth_draw_seed(x->grp->ref, x->grp->read_no, ord);
    x->vseq.resize(p->len); x->rd_fwd.resize(p->len);
    aln_meth_convert(x->grp->read_no, p->seq, p->len, x->vseq.data());
    for (int k = 0; k < (int)p->len; ++k) x->rd_fwd[k] = p->seq[p->len - 1 - k];
    x->pend.clear(); x->windows.clear(); x->intervals.clear(); x->loci.clear();
    x->alignment_bytes = 0; x->unrep_cost = INT64_MAX;
    x->final_flush = x->stop = 0;
    x->flushed_score = -1;
    return x->vseq.data();
}
static bool conv_tiebreak(const aln_meth_opt_t *mo)
{
    return mo && mo->m.meth_chem == METH_CHEM_TAPS &&
        mo->m.meth_scoring == MEM_METH_SCORING_NEUTRAL;
}

static int conv_mapq_cap(int gap)
{
    return gap < 10 ? 6 * gap : 60;
}

// Rank the fixed-refinement CIGAR by stock cost, not by its refinement reward.
static bool alignment_better(const aln_meth_occ_t &a, const aln_meth_occ_t &b, bool conversions)
{
    if (a.cost != b.cost) return a.cost < b.cost;
    if (conversions && a.n_conv != b.n_conv) return a.n_conv < b.n_conv;
    if (a.alignment.cigar != b.alignment.cigar) return a.alignment.cigar < b.alignment.cigar;
    if (a.alignment.re != b.alignment.re) return a.alignment.re < b.alignment.re;
    return a.alignment.raw_cigar < b.alignment.raw_cigar;
}

/* Account on the raw global alignment: terminal gaps still cost O+(k-1)E even
 * when display normalization clips/removes them. Neutral conversions are free.
 * Returns false when the counts cannot be stored in a record (8-bit fields). A refined
 * window can drift far from its search hit (a low-complexity read whose converted form
 * matches a homopolymer, say); callers drop such a candidate, and a read whose dropped
 * candidate could have competed (cost within the retained window) is marked incomplete. */
static bool account_alignment(const aln_meth_ref_t *ref, const aln_meth_opt_t *mo,
                              const gap_opt_t *opt, int64_t rb, int64_t re,
                              const ubyte_t *query, int len, aln_meth_occ_t *o)
{
    int m;
    const int8_t *mat = meth_mat(mo, o->hyp, &m);
    std::vector<uint8_t> win;
    fetch_window(ref, mo, rb, re - rb, o->rid, win);
    int q = 0, t = 0;
    for (uint32_t item : o->alignment.raw_cigar) {
        const int op = item & 15, n = item >> 4;
        if (n <= 0 || op > 2 || (op != 2 && n > len - q) ||
            (op != 1 && n > (int)win.size() - t)) {
            fprintf(stderr, "ERROR: --meth invalid raw CIGAR accounting\n"); exit(1);
        }
        if (op == 0) {
            for (int k = 0; k < n; ++k) {
                const int read = query[q++], base = win[t++];
                if (read > 3 || mat[base * m + read] < 0) ++o->n_mm;
                else if (read != meth_code_base(base)) ++o->n_conv;
            }
        } else {
            ++o->n_gapo; o->n_gape += n - 1;
            if (op == 1) { o->n_ins += n; q += n; }
            else { o->n_del += n; t += n; }
        }
    }
    if (q != len || t != (int)win.size() || opt->s_mm < 0 || opt->s_gapo < 0 || opt->s_gape < 0) {
        fprintf(stderr, "ERROR: --meth invalid alignment dimensions or penalties\n"); exit(1);
    }
    // Counts are bounded by query+window lengths (<=131070); int penalties fit int64 products.
    o->cost = (int64_t)opt->s_mm * o->n_mm + (int64_t)opt->s_gapo * o->n_gapo +
              (int64_t)opt->s_gape * o->n_gape;
    return o->cost <= INT_MAX && o->n_mm <= 255 && o->n_gapo <= 255 && o->n_gape <= 255;
}
static void resolve_window(aln_meth_slot_t *x, const aln_meth_pend_t &h, bwtint_t sa)
{
    aln_meth_group_t *g = x->grp;
    const int span = x->len + (int)h.e.n_del - (int)h.e.n_ins;
    if (span <= 0) return;
    int strand;
    const bwtint_t pos = bwa_pos_from_sa(g->ref->conv_bns, sa, span, &strand);
    if (pos == (bwtint_t)-1) return;
    const int crid = bns_pos2rid(g->ref->conv_bns, pos);
    if (crid < 0) return;
    const bntann1_t &ca = g->ref->conv_bns->anns[crid];
    if (pos < (bwtint_t)ca.offset || pos - ca.offset > (bwtint_t)ca.len ||
        (bwtint_t)span > (bwtint_t)ca.len - (pos - ca.offset)) return;
    const int hyp = crid & 1, rid = crid >> 1;
    if ((hyp ^ strand) != g->read_no) { ++x->res->n_wrong_copy; return; }
    const int64_t rb = g->ref->bns->anns[rid].offset + (pos - ca.offset);
    const auto key = std::make_tuple(rid, hyp, strand, rb, rb + span);
    if (!x->windows.insert(key).second) return;
    aln_meth_occ_t o;
    o.rid = rid; o.hyp = hyp; o.strand = strand;
    const aln_meth_candidate_t candidate{rid, strand, hyp, rb, rb + span};
    const ubyte_t *query = strand ? x->p->rseq : x->rd_fwd.data();
    const auto status = aln_meth_realign_window(*g->ref, *g->mo, candidate, query, x->len,
                                               h.e.n_ins, h.e.n_del, &o.alignment);
    // Invalid/resource/empty normalized alignments must never be silently omitted competitors.
    if (status != aln_meth_eval_status::ok) {
        fprintf(stderr, "ERROR: --meth cannot evaluate a generated candidate window; stopping\n");
        exit(1);
    }
    if (!account_alignment(g->ref, g->mo, g->opt, rb, rb + span, query, x->len, &o)) {
        ++x->res->n_unrepresentable;
        x->unrep_cost = std::min(x->unrep_cost, o.cost);
        return;
    }
    ++x->res->n_rescored;
    const uint64_t bytes = (o.alignment.raw_cigar.size() + o.alignment.cigar.size()) * sizeof(uint32_t);
    // Cumulative CIGAR bytes conservatively bound ownership even when a locus is replaced.
    // This is separate from the wrapper's active-DP budget.
    if (bytes > (64ULL << 20) - x->alignment_bytes) {
        x->res->flags |= ALN_METH_BIND_CIGAR;
        x->stop = 1;
        return;
    }
    x->alignment_bytes += bytes;
    o.pos = o.alignment.rb;
    const auto locus = std::make_tuple(rid, hyp, strand, o.alignment.rb);
    const auto it = x->loci.find(locus);
    if (it == x->loci.end()) {
        x->loci.emplace(locus, x->res->occ.size());
        x->res->occ.push_back(std::move(o));
    } else if (alignment_better(o, x->res->occ[it->second], conv_tiebreak(g->mo)))
        x->res->occ[it->second] = std::move(o);
}
void aln_meth_flush_parked(aln_meth_group_t *g)
{
    const int64_t budget = std::max<int64_t>(0, (int64_t)g->mo->p_best + g->mo->p_second);
    if (budget > 65536) { fprintf(stderr, "ERROR: --meth lookup budget exceeds 65536 rows per read\n"); exit(1); }
    for (int i = 0; i < g->n_slots; ++i) {
        if (g->cores[i].state != GAP_SLOT_PARKED) continue;
        aln_meth_slot_t &x = g->slot[i];
        ++x.res->n_flushes;
        for (const auto &h : x.pend) {
            if (x.stop) break;
            const bwtint_t width = h.l - h.k + 1;
            const int64_t remaining = budget - x.res->n_lookups;
            const int n = (int)std::min<uint64_t>(width, (uint64_t)std::max<int64_t>(0, remaining));
            std::vector<bwtint_t> rows(n);
            for (int k = 0; k < n; ++k) rows[k] = h.k + k;
            bwt_sa_batch(g->bwt, rows.data(), n, rows.data());
            for (bwtint_t sa : rows) {
                if (x.stop) break;
                ++x.res->n_lookups;
                resolve_window(&x, h, sa);
            }
            if ((bwtint_t)n < width) { x.res->flags |= ALN_METH_BIND_BEST; x.stop = 1; }
        }
        x.pend.clear();
        g->cores[i].state = GAP_SLOT_SEARCH;
    }
    g->n_parked = 0;
}
// Kept out of line: its one caller is slot_finish, which the search kernel inlines, and LTO
// inlining it there changed the kernel's code (clang 19). aln_meth_redraw_hits must draw the same.
__attribute__((noinline)) void aln_meth_redraw(aln_meth_read_t *res, uint64_t rng)
{
    if (res->main < 0) return;
    int64_t best = INT64_MAX;
    for (const auto &o : res->occ) best = std::min(best, o.cost);
    int fewest = INT_MAX;
    if (res->n_conv_moved)
        for (const auto &o : res->occ) if (o.cost == best) fewest = std::min(fewest, o.n_conv);
    int held = -1, count = 0;
    for (int i = 0; i < (int)res->occ.size(); ++i) {
        if (res->occ[i].cost != best ||
            (res->n_conv_moved && res->occ[i].n_conv != fewest)) continue;
        ++count;
        if (count == 1 || (splitmix64(&rng) % (uint64_t)count) == 0) held = i;
    }
    res->main = held;
}
// Packs a locus into a hit; a CIGAR of more than one operation goes to the cig table.
static aln_meth_hit_t meth_hit(const aln_meth_occ_t &o, std::vector<uint32_t> &cig)
{
    const std::vector<uint32_t> &c = o.alignment.cigar;
    const int64_t span = o.alignment.re - o.alignment.rb;
    // pos is the alignment's start, so pos + span is its end (cap rescue's anchor bounds).
    if ((!c.empty() && o.alignment.rb != (int64_t)o.pos) ||
        span < INT32_MIN || span > INT32_MAX || o.cost < 0 || o.cost > INT_MAX ||
        o.n_mm < 0 || o.n_mm > 255 || o.n_gapo < 0 || o.n_gapo > 255 || o.n_gape < 0 || o.n_gape > 255 ||
        cig.size() > UINT32_MAX - c.size()) {
        fprintf(stderr, "ERROR: --meth locus does not fit its stored record\n"); exit(1);
    }
    aln_meth_hit_t h;
    h.pos = o.pos; h.span = (int32_t)span; h.cost = (int32_t)o.cost; h.rid = o.rid; h.n_conv = o.n_conv;
    h.n_cig = (uint32_t)c.size();
    if (c.size() == 1) h.cig = c[0];
    else if (c.size() > 1) { h.cig = (uint32_t)cig.size(); cig.insert(cig.end(), c.begin(), c.end()); }
    else h.cig = 0; // no CIGAR: no table offset to go stale
    h.strand = (uint8_t)o.strand; h.hyp = (uint8_t)o.hyp;
    h.n_mm = (uint8_t)o.n_mm; h.n_gapo = (uint8_t)o.n_gapo; h.n_gape = (uint8_t)o.n_gape;
    return h;
}
void aln_meth_store_loci(aln_meth_read_t *res)
{
    if (res->occ.empty()) return;
    xassert(res->hits.empty() && res->cig.empty(), "--meth loci stored twice"); // main indexes occ, then hits
    size_t n_cig = 0;
    for (const auto &o : res->occ) if (o.alignment.cigar.size() > 1) n_cig += o.alignment.cigar.size();
    res->hits.reserve(res->occ.size());
    res->cig.reserve(n_cig);
    for (const auto &o : res->occ) res->hits.push_back(meth_hit(o, res->cig));
    std::vector<aln_meth_occ_t>().swap(res->occ);
}
// Finishes a read whose main hit is drawn: stores its loci as hits and frees what only the search
// used (occ: raw_cigar, the counts and alignment fields no reader of a finished read uses). An
// incomplete read (flags != 0) is never paired or mate-rescued, never anchors cap rescue (which
// reads only its main hit's cost, as a bound) and reports no X0/X1/XA, so it keeps only its
// best-cost loci, in their order: a redraw of a duplicate copy ranges over the same loci, and
// c1/c2 still count every retained locus. Repeating it is a no-op. Kept out of line: one body
// shared by both callers.
static __attribute__((noinline)) void meth_release_finished(aln_meth_read_t *res)
{
    if (res->main < 0 || res->occ.empty()) return;
    if (res->flags) {
        const int64_t best = res->occ[res->main].cost;
        size_t n = 0;
        for (size_t i = 0; i < res->occ.size(); ++i) {
            if (res->occ[i].cost > best) continue;
            if ((int64_t)i == res->main) res->main = (int)n;
            if (n != i) res->occ[n] = std::move(res->occ[i]);
            ++n;
        }
        res->occ.resize(n);
    }
    aln_meth_store_loci(res);
}
void aln_meth_trim(aln_meth_read_t *res, int max_xa_hits)
{
    if (res->main < 0) return;
    // A read not yet stored stores only its main locus. That locus has the best cost, so an
    // incomplete read's best-cost filter keeps it.
    if (!res->occ.empty() && (res->flags || (int64_t)res->occ.size() > (int64_t)max_xa_hits + 1)) {
        aln_meth_occ_t main = std::move(res->occ[res->main]);
        res->occ.clear(); res->occ.push_back(std::move(main)); res->main = 0;
    }
    meth_release_finished(res);
    // A stored read (a duplicate copy, after its draw) trims its hits.
    if (res->hits.size() > 1 && (res->flags || (int64_t)res->hits.size() > (int64_t)max_xa_hits + 1)) {
        aln_meth_hit_t main = res->hits[res->main];
        std::vector<uint32_t> cig;
        if (main.n_cig > 1) {
            const uint32_t *c = aln_meth_hit_cigar(*res, main);
            cig.assign(c, c + main.n_cig);
            main.cig = 0;
        }
        res->cig.swap(cig);
        std::vector<aln_meth_hit_t>(1, main).swap(res->hits); res->main = 0;
    }
}
// aln_meth_redraw's draw over hits: keep the two in step (a unit test compares them seed by seed).
void aln_meth_redraw_hits(aln_meth_read_t *res, uint64_t rng)
{
    if (res->main < 0) return;
    int64_t best = INT64_MAX;
    for (const auto &o : res->hits) best = std::min<int64_t>(best, o.cost);
    int fewest = INT_MAX;
    if (res->n_conv_moved)
        for (const auto &o : res->hits) if (o.cost == best) fewest = std::min(fewest, (int)o.n_conv);
    int held = -1, count = 0;
    for (int i = 0; i < (int)res->hits.size(); ++i) {
        if (res->hits[i].cost != best ||
            (res->n_conv_moved && res->hits[i].n_conv != fewest)) continue;
        ++count;
        if (count == 1 || (splitmix64(&rng) % (uint64_t)count) == 0) held = i;
    }
    res->main = held;
}
// The out-of-line tail of aln_meth_slot_finish. keep 0, and STRATA on an incomplete read (pairing
// never reads it past its main hit): trim. STRATA on a complete read, and COPIES (copies redraw
// over the ties; aln_meth_search_batch trims after drawing them): release only.
// Keep this work here, not in slot_finish or aln_meth_redraw: slot_finish is inlined into the
// search kernel, and under LTO a change to either altered the kernel's search-loop code
// (clang 19, x86 and arm64).
static __attribute__((noinline)) void meth_slot_finish_tail(aln_meth_read_t *res, int keep, int max_xa_hits)
{
    if (!keep || (keep == ALN_METH_KEEP_STRATA && res->flags)) aln_meth_trim(res, max_xa_hits);
    else meth_release_finished(res);
}
void aln_meth_slot_finish(aln_meth_slot_t *x)
{
    aln_meth_read_t *r = x->res;
    if (r->occ.empty()) return;
    int64_t best = INT64_MAX;
    for (const auto &o : r->occ) best = std::min(best, o.cost);
    const int64_t cut = (x->grp->opt->mode & BWA_MODE_NONSTOP) ? INT64_MAX
        : best + (int64_t)x->grp->opt->s_mm;
    // A dropped candidate that could be a best or retained worse locus leaves X0/X1 unknown.
    if (r->n_unrepresentable && x->unrep_cost <= cut) r->flags |= ALN_METH_BIND_UNREP;
    std::vector<aln_meth_occ_t> keep;
    for (auto &o : r->occ) {
        if (o.cost > cut) continue;
        if (o.cost == best) ++r->c1;
        else ++r->c2; // stock c2 counts every retained worse cost, not only the next distinct one
        keep.push_back(std::move(o));
    }
    std::sort(keep.begin(), keep.end(), [](const aln_meth_occ_t &a, const aln_meth_occ_t &b) {
        return std::make_tuple(a.rid, a.hyp, a.strand, a.pos) < std::make_tuple(b.rid, b.hyp, b.strand, b.pos);
    });
    r->occ.swap(keep); r->main = 0;
    if (conv_tiebreak(x->grp->mo) && !r->flags) {
        int fewest = INT_MAX, gap = INT_MAX;
        for (const auto &o : r->occ) if (o.cost == best) fewest = std::min(fewest, o.n_conv);
        for (const auto &o : r->occ) if (o.cost == best && o.n_conv > fewest) {
            ++r->n_conv_moved;
            gap = std::min(gap, o.n_conv - fewest);
        }
        if (r->n_conv_moved) {
            r->conv_gap = gap;
            r->c1 -= r->n_conv_moved;
            r->c2 += r->n_conv_moved;
        }
    }
    aln_meth_redraw(r, x->rng);
    meth_slot_finish_tail(r, x->keep, x->grp->max_xa_hits);
}
void aln_meth_report(int read_no, int n, const aln_meth_read_t *res)
{
    if (bwa_verbose < 4) return;
    long lookups = 0, wrong = 0, scored = 0, flushes = 0, incomplete = 0, dropped = 0;
    for (int i = 0; i < n; ++i) {
        lookups += res[i].n_lookups; wrong += res[i].n_wrong_copy; scored += res[i].n_rescored;
        dropped += res[i].n_unrepresentable;
        flushes += res[i].n_flushes; incomplete += res[i].flags != 0;
    }
    fprintf(stderr, "[M::aln] --meth read %d: %d reads, %ld SA lookups (%ld wrong-copy), %ld unique windows scored (%ld dropped as unrepresentable), %ld flushes; %ld incomplete; stock cost/count MAPQ\n",
            read_no ? 1 : 2, n, lookups, wrong, scored, dropped, flushes, incomplete);
}

/* bwa_cal_md1 with NM counting only the mismatches the mode does not free (and gap lengths); MD
 * stays literal (conversions are listed, as mem does). */
/* Matrix-parameterized twin of bwase.cpp's bwa_cal_md1; the stock twin stays unchanged. */
static char *meth_cal_md(const aln_meth_ref_t *ref, const aln_meth_opt_t *mo, int hyp, int n_cigar,
						 bwa_cigar_t *cigar, int len, bwtint_t pos, const ubyte_t *seq, kstring_t *str, int *_nm)
{
	const int64_t l_pac = ref->bns->l_pac;
	int m;
	const int8_t *mat = meth_mat(mo, hyp, &m);
	// reference span, labelled for cpg
	int64_t span = 0;
	if (cigar) {
		for (int k = 0; k < n_cigar; ++k)
			if (__cigar_op(cigar[k]) == FROM_M || __cigar_op(cigar[k]) == FROM_D) span += __cigar_len(cigar[k]);
	} else span = len;
	if ((int64_t)pos + span > l_pac) span = l_pac - (int64_t)pos;
	std::vector<uint8_t> win;
	fetch_window(ref, mo, (int64_t)pos, span, bns_pos2rid(ref->bns, (int64_t)pos), win);
	int64_t x = 0, y = 0;
	int u = 0, nm = 0;
	str->l = 0;
	auto put_mm = [&](int64_t xi, int64_t yi) {
		const uint8_t fb = win[xi], c = meth_code_base(fb), rb = seq[yi];
		if (c > 3 || rb > 3 || c != rb) {
			kputl(u, str);
			kputc("ACGTN"[c], str);
			u = 0;
			if (rb > 3 || mat[fb * m + rb] < 0) ++nm;
		} else ++u;
	};
	if (cigar) {
		for (int k = 0; k < n_cigar; ++k) {
			const int l = __cigar_len(cigar[k]), op = __cigar_op(cigar[k]);
			if (op == FROM_M) {
				for (int z = 0; z < l && x + z < span; ++z) put_mm(x + z, y + z);
				x += l; y += l;
			} else if (op == FROM_I || op == FROM_S) {
				y += l;
				if (op == FROM_I) nm += l;
			} else if (op == FROM_D) {
				kputl(u, str);
				kputc('^', str);
				for (int z = 0; z < l && x + z < span; ++z) kputc("ACGTN"[meth_code_base(win[x + z])], str);
				u = 0;
				x += l; nm += l;
			}
		}
	} else {
		for (int z = 0; z < len && z < span; ++z) put_mm(z, z);
	}
	kputl(u, str);
	*_nm = nm;
	return aln_strdup(str->s);
}

/* Restore the known read tail in wide BAM operations before legacy 14-bit splitting.
 * The query is still the trimmed, oriented sequence when MD/NM is computed. */
static bwa_cigar_t *legacy_cigar(const uint32_t *cigar, size_t n_cigar, int *n,
                               int len, int full_len, int strand, int max_ops = INT_MAX)
{
    if (len <= 0 || full_len < len || strand < 0 || strand > 1 || n_cigar == 0) {
        fprintf(stderr, "ERROR: --meth invalid CIGAR/trim dimensions\n"); exit(1);
    }
    std::vector<uint32_t> wide(cigar, cigar + n_cigar);
    const uint64_t trim = (uint64_t)full_len - len;
    if (trim) {
        const size_t end = strand ? 0 : wide.size() - 1;
        const bool merge = (wide[end] & 15) == 4;
        const uint64_t clip = trim + (merge ? wide[end] >> 4 : 0);
        if (clip > (UINT32_MAX >> 4)) {
            fprintf(stderr, "ERROR: --meth restored soft clip exceeds BAM capacity\n"); exit(1);
        }
        const uint32_t item = (uint32_t)(clip << 4) | 4;
        if (merge) wide[end] = item;
        else if (strand) wide.insert(wide.begin(), item);
        else wide.push_back(item);
    }
    uint64_t count = 0, query_span = 0;
    for (uint32_t item : wide) {
        const unsigned code = item & 15, length = item >> 4;
        if (!length || (code != 0 && code != 1 && code != 2 && code != 4)) {
            fprintf(stderr, "ERROR: --meth invalid stored CIGAR\n"); exit(1);
        }
        count += (length + (uint64_t)CIGAR_LN_MASK - 1) / CIGAR_LN_MASK;
        if (code != 2) query_span += length;
    }
    if (query_span != (uint64_t)full_len || count > (uint64_t)max_ops ||
        count > SIZE_MAX / sizeof(bwa_cigar_t)) {
        fprintf(stderr, "ERROR: --meth restored CIGAR exceeds record capacity or query span\n"); exit(1);
    }
    auto *out = (bwa_cigar_t*)aln_malloc((size_t)count * sizeof(bwa_cigar_t));
    int k = 0;
    for (uint32_t item : wide) {
        int length = item >> 4, code = item & 15;
        if (code == 4) code = FROM_S;
        while (length) {
            const int take = std::min(length, CIGAR_LN_MASK);
            out[k++] = __cigar_create(code, take); length -= take;
        }
    }
    *n = k;
    return out;
}
static bwa_cigar_t *legacy_cigar(const std::vector<uint32_t> &cigar, int *n,
                               int len, int full_len, int strand, int max_ops = INT_MAX)
{
    return legacy_cigar(cigar.data(), cigar.size(), n, len, full_len, strand, max_ops);
}
static bwa_cigar_t *legacy_cigar(const aln_meth_read_t *r, const aln_meth_hit_t &h, int *n,
                               int len, int full_len, int strand, int max_ops = INT_MAX)
{
    return legacy_cigar(aln_meth_hit_cigar(*r, h), h.n_cig, n, len, full_len, strand, max_ops);
}
// Select/account first; emission can then keep a pairing decision's MAPQ unchanged.
static void meth_set_main(bwa_seq_t *s, aln_meth_read_t *r, const gap_opt_t *opt)
{
    s->mapQ = s->seQ = 0;
    if (r->main < 0) { s->type = BWA_TYPE_NO_MATCH; r->hyp = -1; return; }
    const auto &o = r->hits[r->main];
    s->pos = o.pos; s->strand = o.strand; s->score = (int)o.cost; r->hyp = o.hyp;
    s->type = r->flags || r->c1 > 1 ? BWA_TYPE_REPEAT : BWA_TYPE_UNIQUE;
    s->c1 = r->flags ? 0 : r->c1; s->c2 = r->flags ? 0 : r->c2;
    s->n_mm = o.n_mm; s->n_gapo = o.n_gapo; s->n_gape = o.n_gape;
    s->ref_shift = (int)o.span - (int)s->len;
    const int max_diff = opt->fnr > 0.0 ? bwa_cal_maxdiff(s->len, BWA_AVG_ERR, opt->fnr) : opt->max_diff;
    s->mapQ = s->seQ = r->flags ? 0 : bwa_approx_mapQ(s, max_diff);
    if (r->n_conv_moved && s->mapQ > conv_mapq_cap(r->conv_gap))
        s->mapQ = s->seQ = conv_mapq_cap(r->conv_gap);
}

static void meth_emit(const aln_meth_ref_t *ref, const aln_meth_opt_t *mo,
                      int max_xa_hits, int n, bwa_seq_t *seqs, aln_meth_read_t *res)
{
    kstring_t str = {0, 0, 0};
    for (int i = 0; i < n; ++i) {
        bwa_seq_t *s = seqs + i; aln_meth_read_t *r = res + i;
        seq_reverse(s->len, s->seq, 0);
        if (s->type == BWA_TYPE_NO_MATCH) { s->len = s->full_len; continue; }
        const auto &o = r->hits[r->main];
        s->cigar = legacy_cigar(r, o, &s->n_cigar, s->len, s->len, s->strand);
        int nm;
        s->md = meth_cal_md(ref, mo, o.hyp, s->n_cigar, s->cigar, s->len, s->pos,
                            s->strand ? s->rseq : s->seq, &str, &nm);
        if (nm > 4095) { fprintf(stderr, "ERROR: --meth NM exceeds record capacity\n"); exit(1); }
        s->nm = nm;
        s->n_multi = 0; s->multi = nullptr;
        if (!r->flags && r->hits.size() > 1 && (int64_t)r->hits.size() <= (int64_t)max_xa_hits + 1) {
            s->multi = (bwt_multi1_t*)aln_calloc(r->hits.size()-1, sizeof(bwt_multi1_t));
            for (int k = 0; k < (int)r->hits.size(); ++k) {
                if (k == r->main) continue;
                const auto &alt = r->hits[k]; bwt_multi1_t *q = &s->multi[s->n_multi++];
                q->pos = alt.pos; q->strand = alt.strand;
                int nc, edits;
                q->cigar = legacy_cigar(r, alt, &nc, s->len, s->len, q->strand, 32767);
                if (nc > 32767) { fprintf(stderr, "ERROR: --meth XA CIGAR exceeds record capacity\n"); exit(1); }
                q->n_cigar = nc;
                char *md = meth_cal_md(ref, mo, alt.hyp, nc, q->cigar, s->len, q->pos,
                                      q->strand ? s->rseq : s->seq, &str, &edits);
                free(md);
                if (edits > 255) { // the XA edit field is 8 bits: omit XA for this read, keep the primary
                    for (int j = 0; j < s->n_multi; ++j) free(s->multi[j].cigar);
                    free(s->multi); s->multi = nullptr; s->n_multi = 0;
                    break;
                }
                q->mm = edits; q->gap = 0;
                if (s->len != s->full_len) {
                    free(q->cigar);
                    q->cigar = legacy_cigar(r, alt, &nc, s->len, s->full_len, q->strand, 32767);
                    q->n_cigar = nc;
                }
            }
        }
        if (s->len != s->full_len) {
            free(s->cigar);
            s->cigar = legacy_cigar(r, o, &s->n_cigar, s->len, s->full_len, s->strand);
        }
        s->len = s->full_len;
    }
    free(str.s);
}

void aln_meth_finish_se(const aln_meth_ref_t *ref, const aln_meth_opt_t *mo, const gap_opt_t *opt,
                        int max_xa_hits, int n, bwa_seq_t *seqs, aln_meth_read_t *res)
{
    for (int i = 0; i < n; ++i) meth_set_main(seqs + i, res + i, opt);
    meth_emit(ref, mo, max_xa_hits, n, seqs, res);
}

/* ---------------------------------------------------------------------------------------------
 * sampe under --meth
 * ------------------------------------------------------------------------------------------- */

void aln_meth_pe_stage1(const gap_opt_t *opt, int n, bwa_seq_t *seqs[2], aln_meth_read_t *res[2])
{
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < 2; ++j) {
			bwa_seq_t *p = seqs[j] + i;
			p->extra_flag |= SAM_FPD | (j == 0 ? SAM_FR1 : SAM_FR2);
			meth_set_main(p, res[j] + i, opt);
		}
}

#define SW_MIN_MATCH_LEN 20
#define SW_MIN_MAPQ 17

/* bwape.cpp's bwa_sw_core with the anchor's matrix (labelled window under cpg) on the original
 * reference; *_cnt's n_mm counts only the mismatches the mode does not free. */
static bwa_cigar_t *meth_sw_core(const aln_meth_ref_t *ref, const aln_meth_opt_t *mo, int hyp, int rid,
                                 int len, const ubyte_t *seq, int64_t *beg, int64_t end,
                                 int *n_cigar, uint32_t *_cnt,
                                 const gap_opt_t *account_opt = nullptr, aln_meth_occ_t *owned = nullptr)
{
    *n_cigar = 0;
    if (len <= 0 || len > 65535 || rid < 0 || rid >= ref->bns->n_seqs) return nullptr;
    const auto &ann = ref->bns->anns[rid];
    const int64_t contig_end = ann.offset + ann.len;
    // The anchor owns the window, even when its unbounded start falls in a neighbor.
    *beg = std::max<int64_t>(*beg, ann.offset);
    end = std::min(end, contig_end);
    const int64_t span = end - *beg;
    if (span < SW_MIN_MATCH_LEN) return nullptr;
    // The local-rescue profile and traceback must remain within a bounded evaluation window.
    if (span > 65535 || 2ULL * len * span + 128ULL * (len + span) > (64ULL << 20)) {
        fprintf(stderr, "ERROR: --meth rescue window exceeds allocation budget\n"); exit(1);
    }
    int ambiguous = 0;
    for (int k = 0; k < len; ++k) ambiguous += seq[k] >= 4;
    if (ambiguous * 4 >= len || len - ambiguous < SW_MIN_MATCH_LEN) return nullptr;
    int m;
    const int8_t *mat = meth_mat(mo, hyp, &m);
    std::vector<uint8_t> target;
    fetch_window(ref, mo, *beg, span, rid, target);
    const int xtra = KSW_XSUBO | KSW_XSTART | (len < 250 ? KSW_XBYTE : 0);
    const kswr_t r = ksw_align(len, const_cast<uint8_t*>(seq), (int)span, target.data(), m, mat, 5, 1, xtra, nullptr);
    if (r.score < SW_MIN_MATCH_LEN || r.score2 == r.score || r.qb < 0 || r.tb < 0 ||
        r.qe < r.qb || r.te < r.tb || r.qe >= len || r.te >= span ||
        r.qe-r.qb+1 < SW_MIN_MATCH_LEN || r.te-r.tb+1 < SW_MIN_MATCH_LEN) return nullptr;
    // Reuse the bounded, certified global refiner for the selected local rescue segment.
    const aln_meth_candidate_t candidate{rid, 0, hyp, *beg+r.tb, *beg+r.te+1};
    aln_meth_alignment_t alignment;
    if (aln_meth_realign_window(*ref, *mo, candidate, seq+r.qb, r.qe-r.qb+1, 0, 0, &alignment)
        != aln_meth_eval_status::ok) {
        fprintf(stderr, "ERROR: --meth cannot evaluate mate rescue\n"); exit(1);
    }
    if (alignment.reward != r.score) return nullptr; // stock rejects inconsistent local/global scores
    if (alignment.rb < ann.offset || alignment.re > contig_end || alignment.rb >= alignment.re) {
        fprintf(stderr, "ERROR: --meth rescued interval leaves anchor contig\n"); exit(1);
    }
    int q = r.qb, t = r.tb, mm = 0, opens = 0, extensions = 0;
    for (uint32_t item : alignment.raw_cigar) {
        const int op=item&15, length=item>>4;
        if (op==0) {
            for (int k=0;k<length;++k) {
                const int base=target[t++], read=seq[q++];
                if (read>3 || mat[base*m+read]<0) ++mm;
            }
        } else {
            ++opens; extensions+=length-1;
            if(op==1) q+=length; else t+=length;
        }
    }
    if (mm>255 || opens>255 || extensions>255) return nullptr; // not representable: no rescue
    *_cnt=(uint32_t)mm<<16 | (uint32_t)opens<<8 | extensions;
    // Local-query clipping and global terminal normalization are merged before legacy splitting.
    auto cigar=alignment.cigar;
    if (r.qb) {
        if ((cigar.front()&15)==4) cigar.front() += (uint32_t)r.qb<<4;
        else cigar.insert(cigar.begin(), (uint32_t)r.qb<<4|4);
    }
    const int tail=len-r.qe-1;
    if (tail) {
        if ((cigar.back()&15)==4) cigar.back() += (uint32_t)tail<<4;
        else cigar.push_back((uint32_t)tail<<4|4);
    }
    if (owned) {
        *owned = aln_meth_occ_t();
        owned->rid=rid;owned->hyp=hyp;owned->pos=alignment.rb;
        owned->alignment=alignment;
        if (!account_alignment(ref,mo,account_opt,candidate.rb,candidate.re,seq+r.qb,r.qe-r.qb+1,owned))
            return nullptr;
        owned->alignment.cigar=cigar;
    }
    *beg=alignment.rb;
    return legacy_cigar(cigar,n_cigar,len,len,0);
}

/* Bounded local evidence can improve placement, but cannot resolve omitted genomic loci. */
static int meth_cap_rescue(const aln_meth_ref_t *ref, const aln_meth_opt_t *mo,
                          const gap_opt_t *opt, const isize_info_t *ii,
                          bwa_seq_t *p[2], aln_meth_read_t *r[2], int end, bool short_insert,
                          bool ref_span) try
{
    const int mate=1-end;
    auto *pm=p[end];auto *pa=p[mate];auto *ra=r[mate];auto *rm=r[end];
    if (rm->flags & ALN_METH_BIND_CIGAR) return 0; // retain the capped placement without more DP
    if (ra->flags || ra->main<0 || ra->hits.empty()) return 0;
    if (!isfinite(ii->avg) || !isfinite(ii->std) || ii->avg<0 || ii->std<0 ||
        ii->avg>INT_MAX || ii->std>INT_MAX) {
        fprintf(stderr,"ERROR: --meth invalid cap-rescue insert-size range\n");exit(1);
    }
    const int64_t anchor_cost=ra->hits[ra->main].cost;
    int anchors=0;
    for(const auto &a:ra->hits) anchors+=a.cost==anchor_cost;
    if(anchors>255) return 0; // Refuse instead of truncating the anchor competitors.
    if(pm->len<=0 || pm->len>65535) {
        fprintf(stderr,"ERROR: --meth cap-rescue query length out of range\n");exit(1);
    }
    const int64_t bound=rm->main<0?INT64_MAX:rm->hits[rm->main].cost;
    std::map<std::tuple<int,int,int,int64_t,int64_t>,aln_meth_occ_t> windows;
    uint64_t bytes=0;
    int best_anchor=-1;
    int64_t best_score=INT64_MAX;
    std::pair<int64_t,int> best_rank={INT64_MAX,INT_MAX},runner_rank=best_rank;
    int tied_best=0,valid_pairs=0;
    int best_conv=INT_MAX;
    aln_meth_occ_t best;
    auto canonical=[&](const aln_meth_occ_t &o,const aln_meth_hit_t &a) {
        return std::make_tuple(std::string(ref->bns->anns[o.rid].name),
            o.pos-ref->bns->anns[o.rid].offset,o.strand,o.hyp,
            o.alignment.cigar,o.alignment.re-o.pos,a.pos-ref->bns->anns[a.rid].offset,a.strand);
    };
    for(int index=0;index<(int)ra->hits.size();++index) {
        const auto &a=ra->hits[index];if(a.cost!=anchor_cost) continue;
        const int64_t a_re=(int64_t)a.pos+a.span;
        if(a.rid<0 || a.rid>=ref->bns->n_seqs || a.span<=0) {
            fprintf(stderr,"ERROR: --meth invalid cap-rescue anchor\n");exit(1);
        }
        const auto &ann=ref->bns->anns[a.rid];
        if(a.pos<(uint64_t)ann.offset || a_re>ann.offset+ann.len) {
            fprintf(stderr,"ERROR: --meth cap-rescue anchor leaves contig\n");exit(1);
        }
        const double begin=a.strand ? a_re-ii->avg-3*ii->std-pm->len*.5
                                   : a.pos+ii->avg-3*ii->std-pm->len*1.5;
        const double finish=begin+6*ii->std+2*pm->len;
        if(!isfinite(begin) || !isfinite(finish) || begin<INT64_MIN || finish>=INT64_MAX) {
            fprintf(stderr,"ERROR: --meth cap-rescue window out of range\n");exit(1);
        }
        int64_t beg=(int64_t)begin,last=(int64_t)finish;
        // As in mate rescue, --short-insert lets the window cover the anchor's own span.
        if(a.strand) last=std::min<int64_t>(last,short_insert?a_re:(int64_t)a.pos);
        else beg=std::max<int64_t>(beg,short_insert?(int64_t)a.pos:a_re);
        beg=std::max<int64_t>(beg,ann.offset);last=std::min<int64_t>(last,ann.offset+ann.len);
        if(last<=beg) continue;
        const auto key=std::make_tuple(a.rid,a.hyp,1-a.strand,beg,last);
        aln_meth_occ_t candidate;
        const auto cached=windows.find(key);
        if(cached!=windows.end()) {
            if(cached->second.alignment.cigar.empty()) continue;
            candidate=cached->second;
        } else {
            ubyte_t *query=a.strand?pm->seq:pm->rseq;
            if(a.strand) seq_reverse(pm->len,query,0);
            int nc=0;uint32_t cnt=0;
            auto *cigar=meth_sw_core(ref,mo,a.hyp,a.rid,pm->len,query,&beg,last,&nc,&cnt,opt,&candidate);
            if(a.strand) seq_reverse(pm->len,query,0);
            if(!cigar) {windows.emplace(key,aln_meth_occ_t());continue;}
            free(cigar);
            candidate.strand=1-a.strand;
            windows.emplace(key,candidate);
        }
        const auto &display=candidate.alignment.cigar;
        const uint64_t used=(display.size()+candidate.alignment.raw_cigar.size())*sizeof(uint32_t);
        if(used>(64ULL<<20)-bytes) {
            fprintf(stderr,"ERROR: --meth cap-rescue CIGAR budget exceeded\n");exit(1);
        }
        bytes+=used;
        if(display.empty() || (display.front()&15)==4 || (display.back()&15)==4 || candidate.cost>bound) continue;
        // Keep the candidate inside the fragment: the forward read neither starts nor ends after the reverse.
        if(short_insert && (a.strand ? (candidate.alignment.rb>(int64_t)a.pos || candidate.alignment.re>a_re)
                                     : (candidate.alignment.rb<(int64_t)a.pos || candidate.alignment.re<a_re)))
            continue;
        const int64_t length=a.strand?a_re-candidate.alignment.rb:candidate.alignment.re-a.pos;
        // --ref-span-insert: the outer span must cover both alignments' reference spans, as in pairing.
        if(length<(ref_span?std::max<int64_t>(a.span,candidate.alignment.re-candidate.alignment.rb)
                            :(int64_t)std::max(p[0]->len,p[1]->len))) continue;
        int64_t score=(a.cost+candidate.cost)*10;
        if(ii->high) {
            if(ii->std<=0) {fprintf(stderr,"ERROR: --meth invalid cap-rescue insert deviation\n");exit(1);}
            const double penalty=-4.343*log(.5*erfc(M_SQRT1_2*fabs(length-ii->avg)/ii->std));
            // An underflowed Gaussian tail cannot support this candidate pair.
            if(!isfinite(penalty)) continue;
            // Match pairing's opt-in rejection of numerically unrepresentable likelihoods.
            if(short_insert && penalty+.499>INT_MAX) continue;
            if(penalty>INT_MAX) {
                fprintf(stderr,"ERROR: --meth cap-rescue pair score out of range\n");exit(1);
            }
            score+=(int)(penalty+.499);
        }
        const int conversions=conv_tiebreak(mo)?a.n_conv+candidate.n_conv:0;
        ++valid_pairs;
        const auto rank=std::make_pair(score,conversions);
        if(rank<best_rank) {runner_rank=best_rank;best_rank=rank;tied_best=1;}
        else {runner_rank=std::min(runner_rank,rank);if(rank==best_rank) ++tied_best;}
        const bool better=best_anchor<0 || score<best_score ||
            (score==best_score && (conversions<best_conv ||
             (conversions==best_conv && canonical(candidate,a)<canonical(best,ra->hits[best_anchor]))));
        if(better) {
            best_score=score;best_conv=conversions;
            best_anchor=index;best=std::move(candidate);
        }
    }
    if(best_anchor<0) return 0;
    if(bwa_verbose>=4) fprintf(stderr,"[M::aln] cap rescue: %d evaluated pairs, %d best ties; "
        "best cost %lld, runner-up cost %lld; incomplete MAPQ 0\n",valid_pairs,tied_best,
        (long long)best_rank.first,(long long)(runner_rank.first==INT64_MAX?-1:runner_rank.first));
    const auto &a=ra->hits[best_anchor];
    const bool unchanged=ra->main==best_anchor;
    ra->main=best_anchor;ra->hyp=a.hyp;
    const int saved_q=pa->mapQ,saved_se=pa->seQ;
    meth_set_main(pa,ra,opt);
    if(unchanged && ra->c1==1) {pa->mapQ=saved_q;pa->seQ=saved_se;}
    else pa->mapQ=pa->seQ=0;
    std::vector<uint32_t> cig;const aln_meth_hit_t stored=meth_hit(best,cig);
    rm->cig.swap(cig);std::vector<aln_meth_hit_t>(1,stored).swap(rm->hits);
    rm->main=0;rm->hyp=a.hyp;
    rm->c1=rm->c2=0;
    const auto &o=rm->hits[0];
    free(pm->cigar);pm->cigar=legacy_cigar(rm,o,&pm->n_cigar,pm->len,pm->len,o.strand);
    pm->pos=o.pos;pm->strand=o.strand;pm->type=BWA_TYPE_MATESW;
    pm->score=(int)o.cost;pm->n_mm=o.n_mm;pm->n_gapo=o.n_gapo;pm->n_gape=o.n_gape;
    pm->ref_shift=(int)o.span-(int)pm->len;
    pm->mapQ=pm->seQ=pm->c1=pm->c2=0;
    pm->extra_flag&=~SAM_FPP;pa->extra_flag&=~SAM_FPP;
    return 1;
} catch(const std::bad_alloc &) {
    fprintf(stderr,"ERROR: --meth cannot allocate cap-rescue candidates\n");exit(1);
}

int aln_meth_pe_stage2(const aln_meth_ref_t *ref, const aln_meth_opt_t *mo, const gap_opt_t *gopt, const pe_opt_t *popt,
					   const isize_info_t *ii, int n, bwa_seq_t *seqs[2], aln_meth_read_t *res[2],
					   uint64_t n_tot[2], uint64_t n_mapped[2])
{
	int cnt_chg = 0;
	// Complete means no lookup/work cutoff in configured discovery, not all genome hits.
	const int max_occ = popt->max_occ;
	std::vector<bwt_aln1_t> a[2];
	std::vector<int32_t> span[2]; // --ref-span-insert: each hit's reference span, for pairing
	std::vector<pair64_t> arr;
	kstring_t str = {0, 0, 0};
    std::vector<std::pair<int,int>> initial_selection(n);
    for(int i=0;i<n;++i) initial_selection[i]={res[0][i].main,res[1][i].main};
	for (int i = 0; i < n; ++i) {
		bwa_seq_t *p[2] = {seqs[0] + i, seqs[1] + i};
		aln_meth_read_t *r[2] = {res[0] + i, res[1] + i};
		// Explicitly incomplete ends cannot be promoted by pairing or rescue.
		const bool incomplete = r[0]->flags || r[1]->flags;
		// pairing on the retained stock-cost occurrences
		if ((p[0]->type == BWA_TYPE_UNIQUE || p[0]->type == BWA_TYPE_REPEAT)
			&& (p[1]->type == BWA_TYPE_UNIQUE || p[1]->type == BWA_TYPE_REPEAT)) {
			if (!incomplete && (int)r[0]->hits.size() <= max_occ && (int)r[1]->hits.size() <= max_occ) {
				arr.clear();
				for (int j = 0; j < 2; ++j) {
					a[j].assign(r[j]->hits.size(), bwt_aln1_t());
					if (popt->ref_span_insert) {
						span[j].resize(r[j]->hits.size());
						for (size_t x = 0; x < r[j]->hits.size(); ++x) span[j][x] = r[j]->hits[x].span;
					}
					for (size_t x = 0; x < r[j]->hits.size(); ++x) {
						const aln_meth_hit_t &o = r[j]->hits[x];
						a[j][x].n_mm = o.n_mm; a[j][x].n_gapo = o.n_gapo; a[j][x].n_gape = o.n_gape;
						if (o.cost > ((1<<20)-1)) {
                            fprintf(stderr, "ERROR: --meth pairing cost exceeds record capacity\n"); exit(1);
                        }
                        a[j][x].score = o.cost;
						pair64_t q;
						q.x = o.pos;
						q.y = (uint64_t)x << 2 | (uint64_t)o.strand << 1 | (uint64_t)j;
						arr.push_back(q);
					}
				}
				pair64_v v;
				v.a = arr.data(); v.n = v.m = arr.size();
				bwt_aln1_t *save[2] = {p[0]->aln, p[1]->aln};
				p[0]->aln = a[0].data(); p[1]->aln = a[1].data();
                // Stock multiplies s_mm by ten; clamp only above all representable pair-score gaps.
                const int32_t *spans[2] = {span[0].data(), span[1].data()};
                cnt_chg += bwa_pe_pairing_same_contig(ref->bns, p, &v, popt,
                                                     std::min(gopt->s_mm, INT_MAX/10), ii,
                                                     popt->ref_span_insert ? spans : nullptr);
				p[0]->aln = save[0]; p[1]->aln = save[1];
				// the ends' new loci (pairing may move either): their hypotheses and gaps
				for (int j = 0; j < 2; ++j) {
					for (size_t x = 0; x < r[j]->hits.size(); ++x) {
						const aln_meth_hit_t &o = r[j]->hits[x];
						if (o.pos == p[j]->pos && o.strand == (int)p[j]->strand) {
							r[j]->main = (int)x; r[j]->hyp = o.hyp;
                            p[j]->n_mm=o.n_mm; p[j]->n_gapo=o.n_gapo; p[j]->n_gape=o.n_gape;
                            p[j]->score=(int)o.cost;
                            p[j]->ref_shift=(int)o.span-(int)p[j]->len;
							break;
						}
					}
				}
                // Summing unchanged ends' qualities is not independent evidence for a conversion choice.
                // Compare owned selections: an unchanged end can also start with seQ=0.
                for (int j = 0; j < 2; ++j)
                    if (r[j]->n_conv_moved && r[j]->main == (j ? initial_selection[i].second : initial_selection[i].first) && p[j]->mapQ > conv_mapq_cap(r[j]->conv_gap))
                        p[j]->mapQ = conv_mapq_cap(r[j]->conv_gap);
				xassert(!((p[0]->extra_flag & SAM_FPP) && r[0]->hyp != r[1]->hyp), "a proper pair's ends must share a hypothesis");
			}
		}

	}
    if(popt->is_sw && ii->avg>=0 && popt->type==BWA_PET_STD)
        for(int i=0;i<n;++i) {
            bwa_seq_t *p[2]={seqs[0]+i,seqs[1]+i};
            aln_meth_read_t *r[2]={res[0]+i,res[1]+i};
            const bool cap0=r[0]->flags&ALN_METH_BIND_BEST,cap1=r[1]->flags&ALN_METH_BIND_BEST;
            if(cap0!=cap1) cnt_chg+=meth_cap_rescue(ref,mo,gopt,ii,p,r,cap0?0:1,popt->short_insert,
                                                popt->ref_span_insert);
        }
	// mate rescue (bwa_paired_sw_core with the anchor's matrix; the rescued mate takes its hypothesis)
	if (popt->is_sw && ii->avg >= 0.0) {
		const int64_t l_pac = ref->bns->l_pac;
		for (int i = 0; i < n; ++i) {
			bwa_seq_t *p[2] = {seqs[0] + i, seqs[1] + i};
			aln_meth_read_t *r[2] = {res[0] + i, res[1] + i};
			if (r[0]->flags || r[1]->flags) continue;
			if (!((p[0]->mapQ >= SW_MIN_MAPQ || p[1]->mapQ >= SW_MIN_MAPQ) && (p[0]->extra_flag & SAM_FPP) == 0)) continue;
			int n_cigar[2] = {0, 0}, mapQ = 0, mq_adjust[2] = {255, 255};
			int64_t beg[2], end[2];
			bwa_cigar_t *cigar[2] = {0, 0};
			uint32_t cnt[2];
			const int is_singleton = (p[0]->type == BWA_TYPE_NO_MATCH || p[1]->type == BWA_TYPE_NO_MATCH) ? 1 : 0;
			++n_tot[is_singleton];
			if (popt->type != BWA_PET_STD) continue;
			for (int k = 0; k < 2; ++k) {
				bwa_seq_t *pref = p[1 - k], *pmate = p[k];
				if (pref->type == BWA_TYPE_NO_MATCH) continue;
				ubyte_t *seq;
				if (pref->strand == 0) {
					beg[k] = (int64_t)pref->pos + ii->avg - 3 * ii->std - pmate->len * 1.5;
					end[k] = beg[k] + 6 * ii->std + 2 * pmate->len;
					if (beg[k] < aln_rescue_window_floor(pref, popt)) beg[k] = aln_rescue_window_floor(pref, popt);
					if (end[k] > l_pac) end[k] = l_pac;
					seq = pmate->rseq;
				} else {
					beg[k] = aln_rescue_anchor_end(pref, popt) - ii->avg - 3 * ii->std - pmate->len * 0.5;
					end[k] = beg[k] + 6 * ii->std + 2 * pmate->len;
					if (beg[k] < 0) beg[k] = 0;
					if (end[k] > aln_rescue_window_ceiling(pref, popt)) end[k] = aln_rescue_window_ceiling(pref, popt);
					seq = pmate->seq;
					seq_reverse(pmate->len, seq, 0);
				}
                const int anchor_rid = bns_pos2rid(ref->bns, pref->pos);
                cigar[k] = meth_sw_core(ref, mo, r[1-k]->hyp, anchor_rid, pmate->len, seq,
                                        &beg[k], end[k], &n_cigar[k], &cnt[k]);
				if (cigar[k] && popt->short_insert && !aln_rescue_within_fragment(pref, beg[k], cigar[k], n_cigar[k])) {
					free(cigar[k]); cigar[k] = 0; n_cigar[k] = 0;
				}
				if (cigar[k] && pmate->type != BWA_TYPE_NO_MATCH) {
					int64_t s_old, s_new; int clip = 0;
					if (__cigar_op(cigar[k][0]) == 3) clip += __cigar_len(cigar[k][0]);
					if (__cigar_op(cigar[k][n_cigar[k] - 1]) == 3) clip += __cigar_len(cigar[k][n_cigar[k] - 1]);
					s_old = (int64_t)(((int64_t)pmate->n_mm * gopt->s_mm + (int64_t)pmate->n_gapo * gopt->s_gapo + (int64_t)pmate->n_gape * gopt->s_gape) / 3. * 8. + .499);
					s_new = (int64_t)(((int64_t)(cnt[k] >> 16) * gopt->s_mm + (int64_t)(cnt[k] >> 8 & 0xff) * gopt->s_gapo + (int64_t)(cnt[k] & 0xff) * gopt->s_gape + (int64_t)clip * gopt->s_mm) / 3. * 8. + .499);
					s_old += -4.343 * log(ii->ap_prior / l_pac);
					s_new += (int)(-4.343 * log(.5 * erfc(M_SQRT1_2 * 1.5) + .499));
					if (s_old < s_new) {
						mq_adjust[k] = (int)std::min<int64_t>(255, s_new - s_old);
						free(cigar[k]); cigar[k] = 0; n_cigar[k] = 0;
					} else mq_adjust[k] = (int)std::min<int64_t>(255, s_old - s_new);
				}
				if (pref->strand == 1) seq_reverse(pmate->len, seq, 0);
			}
			int k = -1;
			if (cigar[0] && cigar[1]) {
				k = p[0]->mapQ < p[1]->mapQ ? 0 : 1;
				mapQ = abs(p[1]->mapQ - p[0]->mapQ);
			} else if (cigar[0]) k = 0, mapQ = p[1]->mapQ;
			else if (cigar[1]) k = 1, mapQ = p[0]->mapQ;
			// an unmapped mate's pos is 0, which is not a placement: a rescue there still installs;
			// so does one at the mate's own start on the other strand (the rescue takes the strand
			// opposite the anchor's, so a mate on the anchor's strand is a different placement)
			if (k >= 0 && (p[k]->type == BWA_TYPE_NO_MATCH || p[k]->pos != (bwtint_t)beg[k] ||
						   p[k]->strand == p[1 - k]->strand)) {
				++n_mapped[is_singleton];
				int tmp = (int)p[1 - k]->mapQ - p[k]->mapQ / 2 - 8;
				if (tmp <= 0) tmp = 1;
				if (mapQ > tmp) mapQ = tmp;
				p[k]->mapQ = p[1 - k]->mapQ = mapQ;
				p[k]->seQ = p[1 - k]->seQ = p[1 - k]->seQ < mapQ ? p[1 - k]->seQ : mapQ;
				if (p[k]->mapQ > mq_adjust[k]) p[k]->mapQ = mq_adjust[k];
				if (p[k]->seQ > mq_adjust[k]) p[k]->seQ = mq_adjust[k];
				free(p[k]->cigar); p[k]->cigar = cigar[k]; cigar[k] = 0;
				p[k]->n_cigar = n_cigar[k];
				p[k]->type = BWA_TYPE_MATESW;
				p[k]->pos = beg[k];
				p[k]->seQ = p[1 - k]->seQ;
				p[k]->strand = 1 - p[1 - k]->strand;
				p[k]->n_mm = cnt[k] >> 16; p[k]->n_gapo = cnt[k] >> 8 & 0xff; p[k]->n_gape = cnt[k] & 0xff;
                    const int64_t cost=(int64_t)p[k]->n_mm*gopt->s_mm+(int64_t)p[k]->n_gapo*gopt->s_gapo+(int64_t)p[k]->n_gape*gopt->s_gape;
                    if(cost>INT_MAX) { fprintf(stderr,"ERROR: --meth rescued cost exceeds record capacity\n"); exit(1); }
                    p[k]->score=(int)cost;
				p[k]->extra_flag |= SAM_FPP;
				p[1 - k]->extra_flag |= SAM_FPP;
				r[k]->hyp = r[1 - k]->hyp; // the rescued mate is on the anchor's genome strand
			}
			free(cigar[0]); free(cigar[1]);
		}
	}
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < 2; ++j) {
            bwa_seq_t *s = seqs[j] + i;
            aln_meth_read_t *r = res[j] + i;
            // Rescue must not inflate an unchanged conversion-selected anchor either.
            if (r->n_conv_moved && s->type != BWA_TYPE_MATESW &&
                r->main == (j ? initial_selection[i].second : initial_selection[i].first) &&
                s->mapQ > conv_mapq_cap(r->conv_gap)) s->mapQ = conv_mapq_cap(r->conv_gap);
            int lim = popt->n_multi;
            if (!(s->extra_flag & SAM_FPP) && seqs[1-j][i].type != BWA_TYPE_NO_MATCH)
                lim = (int64_t)s->c1 + s->c2 - 1 > popt->N_multi ? popt->n_multi : popt->N_multi;
            if (s->type == BWA_TYPE_MATESW) {
                seq_reverse(s->len, s->seq, 0);
                int nm;
                s->md = meth_cal_md(ref, mo, r->hyp, s->n_cigar, s->cigar, s->len, s->pos,
                                    s->strand ? s->rseq : s->seq, &str, &nm);
                if (nm > 4095) { fprintf(stderr, "ERROR: --meth rescued NM exceeds capacity\n"); exit(1); }
                s->nm = nm;
                std::vector<uint32_t> cigar;
                for (int k = 0; k < s->n_cigar; ++k) {
                    const int op = __cigar_op(s->cigar[k]);
                    cigar.push_back((uint32_t)__cigar_len(s->cigar[k]) << 4 | (op == FROM_S ? 4 : op));
                }
                free(s->cigar);
                s->cigar = legacy_cigar(cigar, &s->n_cigar, s->len, s->full_len, s->strand);
                s->len = s->full_len;
            } else meth_emit(ref, mo, lim, 1, s, r);
        }
    }
	free(str.s);
	return cnt_chg;
}

/* ---------------------------------------------------------------------------------------------
 * Output
 * ------------------------------------------------------------------------------------------- */

void aln_meth_print_sam1(kstring_t *out, const aln_meth_ref_t *ref, const aln_meth_opt_t *mo, bwa_seq_t *p,
						 const aln_meth_read_t *r, int read_no, const bwa_seq_t *mate, int mode, int max_top2)
{
	const bntseq_t *bns = ref->bns;
	const char *xr = read_no ? "CT" : "GA";
	int j;
	if (p->type != BWA_TYPE_NO_MATCH || (mate && mate->type != BWA_TYPE_NO_MATCH)) {
		int seqid, nn, am = 0, flag = p->extra_flag;
		if (p->type == BWA_TYPE_NO_MATCH) {
			p->pos = mate->pos;
			p->strand = mate->strand;
			flag |= SAM_FSU;
			j = 1;
		} else j = (int)(pos_end(p) - p->pos);
		nn = aln_bns_cnt_ambi(bns, p->pos, j, &seqid);
		const int rec_seqid = seqid; // the XA loop below reuses seqid for the other hits' contigs
		if (p->type != BWA_TYPE_NO_MATCH && (uint64_t)p->pos + (uint64_t)j - (uint64_t)bns->anns[seqid].offset > (uint64_t)bns->anns[seqid].len)
			flag |= SAM_FSU;
		if (p->strand) flag |= SAM_FSR;
		if (mate) {
			if (mate->type != BWA_TYPE_NO_MATCH) {
				if (mate->strand) flag |= SAM_FMR;
			} else flag |= SAM_FMU;
		}
		const int mapped = p->type != BWA_TYPE_NO_MATCH && !(flag & SAM_FSU);
		int hyp = -1;
		if (p->type != BWA_TYPE_NO_MATCH) hyp = r->hyp;
		if (mapped && mo->m.meth_set_as_failed && ((mo->m.meth_set_as_failed == 'f') == (hyp == 1))) flag |= 0x200;
		kputs(p->name, out); kputc('\t', out); kputl(flag, out); kputc('\t', out);
		kputs(bns->anns[seqid].name, out); kputc('\t', out);
		kputl((long)(p->pos - bns->anns[seqid].offset + 1), out); kputc('\t', out);
		kputl(p->mapQ, out); kputc('\t', out);
		if (p->cigar) {
			for (j = 0; j != p->n_cigar; ++j) { kputl(__cigar_len(p->cigar[j]), out); kputc("MIDS"[__cigar_op(p->cigar[j])], out); }
		} else if (p->type == BWA_TYPE_NO_MATCH) kputc('*', out);
		else { kputl(p->len, out); kputc('M', out); }
		if (mate && mate->type != BWA_TYPE_NO_MATCH) {
			int m_seqid;
			long long isize;
			am = mate->seQ < p->seQ ? mate->seQ : p->seQ;
			aln_bns_cnt_ambi(bns, mate->pos, mate->len, &m_seqid);
			kputc('\t', out); kputs((seqid == m_seqid) ? "=" : bns->anns[m_seqid].name, out); kputc('\t', out);
			isize = (seqid == m_seqid) ? pos_5(mate) - pos_5(p) : 0;
			if (p->type == BWA_TYPE_NO_MATCH) isize = 0;
			kputl((long)(mate->pos - bns->anns[m_seqid].offset + 1), out); kputc('\t', out);
			kputl((long)isize, out); kputc('\t', out);
		} else if (mate) {
			kputs("\t=\t", out); kputl((long)(p->pos - bns->anns[seqid].offset + 1), out); kputs("\t0\t", out);
		} else kputs("\t*\t0\t0\t", out);
		bwa_print_seq(out, p);
		kputc('\t', out);
		if (p->qual) {
			if (p->strand) seq_reverse(p->len, p->qual, 0);
			kputs((const char*)p->qual, out);
		} else kputc('*', out);
		if (bwa_rg_id[0]) { kputs("\tRG:Z:", out); kputs(bwa_rg_id, out); }
		if (p->bc[0]) { kputs("\tBC:Z:", out); kputs(p->bc, out); }
		if (p->clip_len < (int)p->full_len) { kputs("\tXC:i:", out); kputl(p->clip_len, out); }
		if (p->type != BWA_TYPE_NO_MATCH) {
			char XT = "NURM"[p->type];
			if (nn > 10) XT = 'N';
			kputs("\tXT:A:", out); kputc(XT, out); kputc('\t', out);
			kputs((mode & BWA_MODE_COMPREAD) ? "NM" : "CM", out); kputs(":i:", out); kputl(p->nm, out);
			if (nn) { kputs("\tXN:i:", out); kputl(nn, out); }
			if (mate) { kputs("\tSM:i:", out); kputl((int)p->seQ, out); kputs("\tAM:i:", out); kputl(am, out); }
			if (p->type != BWA_TYPE_MATESW && !r->flags) {
				kputs("\tX0:i:", out); kputl((int)p->c1, out);
				if ((int)p->c1 <= max_top2) { kputs("\tX1:i:", out); kputl((int)p->c2, out); }
			}
			if (p->md) { kputs("\tMD:Z:", out); kputs(p->md, out); }
			if (p->n_multi) {
				kputs("\tXA:Z:", out);
				for (int i = 0; i < p->n_multi; ++i) {
					bwt_multi1_t *q = p->multi + i;
					j = (int)(pos_end_multi(q, p->len) - q->pos);
					nn = aln_bns_cnt_ambi(bns, q->pos, j, &seqid);
					kputs(bns->anns[seqid].name, out); kputc(',', out); kputc(q->strand ? '-' : '+', out);
					kputl((long)(q->pos - bns->anns[seqid].offset + 1), out); kputc(',', out);
					if (q->cigar) {
						for (int k = 0; k < (int)q->n_cigar; ++k) { kputl(__cigar_len(q->cigar[k]), out); kputc("MIDS"[__cigar_op(q->cigar[k])], out); }
					} else { kputl(p->len, out); kputc('M', out); }
					kputc(',', out); kputl(q->gap + q->mm, out); kputc(';', out);
				}
			}
		}
		if (mo->m.meth_tags & MEM_METH_TAG_XR) { kputs("\tXR:Z:", out); kputs(xr, out); }
		if (mapped && hyp >= 0) {
			if (mo->m.meth_tags & MEM_METH_TAG_XG) { kputs("\tXG:Z:", out); kputs(hyp ? "CT" : "GA", out); }
			if (mo->m.meth_tags & MEM_METH_TAG_XM) {
				// BAM CIGAR ops (S is 4 in BAM, 3 in aln's FROM_S)
				std::vector<uint32_t> bc;
				if (p->cigar) {
					for (int k = 0; k < p->n_cigar; ++k) {
						const int op = __cigar_op(p->cigar[k]);
						bc.push_back((uint32_t)__cigar_len(p->cigar[k]) << 4 | (uint32_t)(op == FROM_S ? 4 : op));
					}
				} else bc.push_back((uint32_t)p->len << 4);
				std::string seq_text;
				seq_text.reserve(p->full_len);
				if (p->strand == 0) for (int i = 0; i < (int)p->full_len; ++i) seq_text.push_back("ACGTN"[p->seq[i]]);
				else for (int i = (int)p->full_len - 1; i >= 0; --i) seq_text.push_back("TGCAN"[p->seq[i]]);
				char *xm = meth_build_xm(bns, ref->pac, rec_seqid, (int64_t)(p->pos - bns->anns[rec_seqid].offset), hyp == 1,
										 bc.data(), (int)bc.size(), seq_text.c_str(), (int)p->full_len, (meth_chem_t)mo->m.meth_chem);
				if (xm) { kputs("\tXM:Z:", out); kputs(xm, out); }
			}
		}
		kputc('\n', out);
	} else {
		int flag = p->extra_flag | SAM_FSU;
		if (mate && mate->type == BWA_TYPE_NO_MATCH) flag |= SAM_FMU;
		kputs(p->name, out); kputc('\t', out); kputl(flag, out); kputs("\t*\t0\t0\t*\t*\t0\t0\t", out);
		bwa_print_seq(out, p);
		kputc('\t', out);
		if (p->qual) {
			if (p->strand) seq_reverse(p->len, p->qual, 0);
			kputs((const char*)p->qual, out);
		} else kputc('*', out);
		if (bwa_rg_id[0]) { kputs("\tRG:Z:", out); kputs(bwa_rg_id, out); }
		if (p->bc[0]) { kputs("\tBC:Z:", out); kputs(p->bc, out); }
		if (p->clip_len < (int)p->full_len) { kputs("\tXC:i:", out); kputl(p->clip_len, out); }
		if (mo->m.meth_tags & MEM_METH_TAG_XR) { kputs("\tXR:Z:", out); kputs(xr, out); }
		kputc('\n', out);
	}
}

} // namespace bwa_aln
