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
static inline uint8_t pac_base(const uint8_t *pac, int64_t x) { return pac[x >> 2] >> ((~x & 3) << 1) & 3; }
static inline const int8_t *meth_mat(const aln_meth_opt_t *mo, int hyp, int *m)
{
    *m = mo->m.meth_context == METH_CTX_CPG ? METH_M : 5;
    return *m == 5 ? (hyp ? mo->mat_ot : mo->mat_ob) : (hyp ? mo->mat7_ot : mo->mat7_ob);
}
static void fetch_window(const aln_meth_ref_t *ref, const aln_meth_opt_t *mo, int64_t rb, int64_t len,
                         int rid, std::vector<uint8_t> &win)
{
    win.resize((size_t)len);
    for (int64_t i = 0; i < len; ++i) win[i] = pac_base(ref->pac, rb + i);
    if (mo->m.meth_context == METH_CTX_CPG && len) meth_label_cpg(win.data(), len, rb, rid, ref->bns, ref->pac);
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
void aln_meth_hit(aln_meth_slot_t *x, const gap_entry_t *e, bwtint_t k, bwtint_t l)
{
    // Reaching the same interval with the same span adds no new candidate window.
    // Other spans must still be evaluated, even if their normalized starts later coincide.
    const int span = x->len + (int)e->n_del - (int)e->n_ins;
    const auto key = std::make_tuple(k, l, span);
    if (!x->intervals.insert(key).second) return;
    x->pend.push_back({k, l, *e});
}
const ubyte_t *aln_meth_slot_start(aln_meth_slot_t *x, bwa_seq_t *p, aln_meth_read_t *res, int64_t ord, int, int keep)
{
    x->p = p; x->res = res; x->len = p->len; x->keep = keep || x->grp->keep_strata;
    *res = aln_meth_read_t();
    x->rng = aln_meth_draw_seed(x->grp->ref, x->grp->read_no, ord);
    x->vseq.resize(p->len); x->rd_fwd.resize(p->len);
    aln_meth_convert(x->grp->read_no, p->seq, p->len, x->vseq.data());
    for (int k = 0; k < (int)p->len; ++k) x->rd_fwd[k] = p->seq[p->len - 1 - k];
    x->pend.clear(); x->windows.clear(); x->intervals.clear(); x->loci.clear();
    x->alignment_bytes = 0; x->unrep_cost = INT64_MAX;
    x->final_flush = x->stop = 0;
    return x->vseq.data();
}
// Rank the fixed-refinement CIGAR by stock cost, not by its refinement reward.
static bool alignment_better(const aln_meth_occ_t &a, const aln_meth_occ_t &b)
{
    if (a.cost != b.cost) return a.cost < b.cost;
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
    } else if (alignment_better(o, x->res->occ[it->second]))
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
void aln_meth_redraw(aln_meth_read_t *res, uint64_t rng)
{
    if (res->main < 0) return;
    int64_t best = INT64_MAX;
    for (const auto &o : res->occ) best = std::min(best, o.cost);
    int held = -1, count = 0;
    for (int i = 0; i < (int)res->occ.size(); ++i) {
        if (res->occ[i].cost != best) continue;
        ++count;
        if (count == 1 || (splitmix64(&rng) % (uint64_t)count) == 0) held = i;
    }
    res->main = held;
}
void aln_meth_trim(aln_meth_read_t *res, int max_xa_hits)
{
    if (res->main < 0) return;
    if (res->flags || (int64_t)res->occ.size() > (int64_t)max_xa_hits + 1) {
        aln_meth_occ_t main = std::move(res->occ[res->main]);
        res->occ.clear(); res->occ.push_back(std::move(main)); res->main = 0;
    }
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
    aln_meth_redraw(r, x->rng);
    if (!x->keep) aln_meth_trim(r, x->grp->max_xa_hits);
}
void aln_meth_report(int read_no, int n, const aln_meth_read_t *res)
{
    if (bwa_verbose < 4 && !getenv("BWA_ALN_METH_STATS")) return;
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
				for (int z = 0; z < l && x + z < span; ++z) kputc("ACGT"[meth_code_base(win[x + z])], str);
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
static bwa_cigar_t *legacy_cigar(const std::vector<uint32_t> &cigar, int *n,
                               int len, int full_len, int strand, int max_ops = INT_MAX)
{
    if (len <= 0 || full_len < len || strand < 0 || strand > 1 || cigar.empty()) {
        fprintf(stderr, "ERROR: --meth invalid CIGAR/trim dimensions\n"); exit(1);
    }
    std::vector<uint32_t> wide = cigar;
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
void aln_meth_finish_se(const aln_meth_ref_t *ref, const aln_meth_opt_t *mo, const gap_opt_t *opt,
                        int max_xa_hits, int n, bwa_seq_t *seqs, aln_meth_read_t *res)
{
    kstring_t str = {0, 0, 0};
    for (int i = 0; i < n; ++i) {
        bwa_seq_t *s = seqs + i; aln_meth_read_t *r = res + i;
        seq_reverse(s->len, s->seq, 0);
        s->mapQ = s->seQ = 0;
        if (r->main < 0) { s->type = BWA_TYPE_NO_MATCH; r->hyp = -1; s->len = s->full_len; continue; }
        const auto &o = r->occ[r->main];
        s->pos = o.pos; s->strand = o.strand; s->score = (int)o.cost; r->hyp = o.hyp;
        s->type = r->flags || r->c1 > 1 ? BWA_TYPE_REPEAT : BWA_TYPE_UNIQUE;
        s->c1 = r->flags ? 0 : r->c1; s->c2 = r->flags ? 0 : r->c2;
        s->cigar = legacy_cigar(o.alignment.cigar, &s->n_cigar, s->len, s->len, s->strand);
        int nm;
        s->md = meth_cal_md(ref, mo, o.hyp, s->n_cigar, s->cigar, s->len, s->pos,
                            s->strand ? s->rseq : s->seq, &str, &nm);
        if (nm > 4095) { fprintf(stderr, "ERROR: --meth NM exceeds record capacity\n"); exit(1); }
        s->nm = nm;
        s->n_mm = o.n_mm; s->n_gapo = o.n_gapo; s->n_gape = o.n_gape;
        const int max_diff = opt->fnr > 0.0 ? bwa_cal_maxdiff(s->len, BWA_AVG_ERR, opt->fnr) : opt->max_diff;
        s->mapQ = s->seQ = r->flags ? 0 : bwa_approx_mapQ(s, max_diff);
        s->n_multi = 0; s->multi = nullptr;
        if (!r->flags && r->occ.size() > 1 && (int64_t)r->occ.size() <= (int64_t)max_xa_hits + 1) {
            s->multi = (bwt_multi1_t*)aln_calloc(r->occ.size()-1, sizeof(bwt_multi1_t));
            for (int k = 0; k < (int)r->occ.size(); ++k) {
                if (k == r->main) continue;
                const auto &alt = r->occ[k]; bwt_multi1_t *q = &s->multi[s->n_multi++];
                q->pos = alt.pos; q->strand = alt.strand;
                int nc, edits;
                q->cigar = legacy_cigar(alt.alignment.cigar, &nc, s->len, s->len, q->strand, 32767);
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
                    q->cigar = legacy_cigar(alt.alignment.cigar, &nc, s->len, s->full_len, q->strand, 32767);
                    q->n_cigar = nc;
                }
            }
        }
        if (s->len != s->full_len) {
            free(s->cigar);
            s->cigar = legacy_cigar(o.alignment.cigar, &s->n_cigar, s->len, s->full_len, s->strand);
        }
        s->len = s->full_len;
    }
    free(str.s);
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
