/* SPDX-License-Identifier: MIT */
/* New code for bwa-mem3 (not from upstream bwa): `aln --meth`. See aln_meth.h and README.md. */
#include "aln_meth.h"
#include "aln_meth_realign.h"
#include "aln_meth_search.h"
#include "aln_known_snvs.h"

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
	aln_meth_fast_init(*mo, &mo->fast);
}

// Pairing is structural: r<name>,f<name> in original order, with identical lengths and
// contiguous offsets. Metadata cannot detect two different sequences with the same layout.
static int meth_ref_check(aln_meth_ref_t *r, const char *orig_desc, const bntseq_t *conv_bns, const char *hint)
{
	r->conv_bns = conv_bns;
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
				"(contig names, lengths or offsets)%s\n", orig_desc, hint);
		aln_meth_ref_free(r);
		return -1;
	}
	return 0;
}

int aln_meth_ref_load(aln_meth_ref_t *r, const char *orig_prefix, const bntseq_t *conv_bns, int n_threads)
{
	r->conv_bns = conv_bns;
	if (meth_orig_ref_load(orig_prefix, &r->bns, &r->pac, n_threads < 1 ? 1 : n_threads) != 0) return -1;
	return meth_ref_check(r, orig_prefix, conv_bns, "; rebuild with `bwa-mem3 index --meth`");
}

int aln_meth_ref_adopt(aln_meth_ref_t *r, bntseq_t *orig_bns, uint8_t *orig_pac, const bntseq_t *conv_bns)
{
	r->bns = orig_bns;
	r->pac = orig_pac;
	return meth_ref_check(r, "(in-memory windows)", conv_bns, " (internal error: the in-memory window references are inconsistent)");
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
// A read's mismatch budget: stock's -n, a fraction (fnr) or a count.
static inline int read_max_diff(const gap_opt_t *opt, int len)
{
    return opt->fnr > 0.0 ? bwa_cal_maxdiff(len, BWA_AVG_ERR, opt->fnr) : opt->max_diff;
}
// --meth-mapq conversion-blind: a locus's cost with its mirror-cell mismatches free.
static inline int64_t mirror_free_cost(int64_t cost, int s_mm, int n_mirror)
{
    return cost - (int64_t)s_mm * n_mirror;
}
static inline const int8_t *meth_mat(const aln_meth_opt_t *mo, int hyp, int *m)
{
    *m = mo->m.meth_context == METH_CTX_CPG ? METH_M : 5;
    return *m == 5 ? (hyp ? mo->mat_ot : mo->mat_ob) : (hyp ? mo->mat7_ot : mo->mat7_ob);
}
static void fetch_window(const aln_meth_ref_t *ref, const aln_meth_opt_t *mo, int64_t rb, int64_t len,
                         int rid, std::vector<uint8_t> &win, int64_t site = -1)
{
    aln_meth_fetch_window(*ref, *mo, rb, (int)len, rid, win, site);
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
 * Seeding (reads longer than -l) admits at most -k differences in the seed, so a seeded
 * search proposes every such window only when floor(T/M) <= -k: no proposal with that
 * many mismatches in all can have more in its seed.
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
    if (job->seed_width && differences>opt.max_seed_diff) return false;
    // Every competitive raw alignment is ungapped. A same-read ALT bound (aln_meth_alt_certificate)
    // strictly above its mismatch allowance rules out all omitted single-allele improvements.
    if (x->grp->ref && x->grp->ref->known_snvs && x->alt_bound.value <= differences) return false;
    return opt.max_diff>=differences && (int64_t)next_score>differences*opt.s_mm;
}

// e's gaps (see ALN_METH_GAP) into gap[0, e.n_gapo)
static void entry_gaps(const aln_meth_slot_t *x, const gap_entry_t &e, uint32_t *gap)
{
    if (e.n_gapo == 1) gap[0] = (uint32_t)e.last_diff_pos;
    else if (e.n_gapo > 1)
        for (int g = 0; g < (int)e.n_gapo; ++g) gap[g] = x->garena[(size_t)e.last_diff_pos * ALN_METH_MAX_GAPS + g];
}
__attribute__((noinline, cold)) int aln_meth_gap_arena(aln_meth_slot_t *x, const gap_entry_t &p, bool ext, uint32_t w)
{
    uint32_t gap[ALN_METH_MAX_GAPS] = {0, 0, 0, 0};
    entry_gaps(x, p, gap);
    int n = p.n_gapo;
    if (ext) gap[n - 1] += 1u << 16;
    else if (n < ALN_METH_MAX_GAPS) gap[n++] = w;
    else n = ALN_METH_MAX_GAPS + 1;
    if (n > ALN_METH_MAX_GAPS || x->garena.size() / ALN_METH_MAX_GAPS >= (size_t)INT_MAX) {
        fprintf(stderr, "ERROR: --meth search path has more gaps than it can record\n"); exit(1);
    }
    const int idx = (int)(x->garena.size() / ALN_METH_MAX_GAPS);
    x->garena.insert(x->garena.end(), gap, gap + ALN_METH_MAX_GAPS);
    return idx;
}
/* Stock aln's cuts, driven by the verified best B (the lowest cost of a stored row; every row's
 * cost is at least its path's converted score, see score_path). next is the score of the entry
 * the search pops next: every entry below it has been expanded and every hit it made resolved.
 * - Score cut: no row of a path above B + M can be retained, and a path's rows cost at least its
 *   score, so the search stops past B + M (job->score_cut), at any time.
 * - Difference ceiling: once B < next, every row that can cost B has been resolved, so the stored
 *   rows at B are final and stock's top-2 ceiling (their fewest differences + 1, at most -n)
 *   applies. It is taken from the stored rows, as slot_finish takes it: a row a same-locus
 *   replacement discarded does not lower it. While B == next, a later row at B could still
 *   replace the stored row with the fewest differences and raise the finish's ceiling, so the
 *   ceiling waits. A row's converted differences are at most its strict ones, so the ceiling
 *   never cuts a path to a retained row.
 * - -R: with the best stratum complete, more than -R loci at B within the ceiling (under the
 *   TAPS/5base tie-break, the fewest-conversion ones of them) end the search, as stock's top2b
 *   stop does.
 * --meth-mapq conversion-blind only observes the search (blind_note_need): the cuts, and every
 * limit, are those of --meth-mapq scoring, so the genomic result is the same. */
static bool conv_tiebreak(const aln_meth_opt_t *mo);
static bool mirror_free_better(const aln_meth_occ_t &a, const aln_meth_occ_t &b, int s_mm);
// --meth-mapq conversion-blind: the read's loci so far for the mirror-free strata, each with its
// rows (its occ row and its mirror table rows, rows[first, first + n)), its occ row (null for a
// locus only the mirror table holds) and its lowest mirror-free row (mirror_free_better). occ loci
// first, in occ order, then the table's other loci in table order.
struct blind_locus_t { const aln_meth_occ_t *occ, *best; uint32_t first, n; };
static void blind_loci(const aln_meth_read_t *r, int s_mm, std::vector<blind_locus_t> &loci,
                       std::vector<const aln_meth_occ_t *> &rows)
{
    loci.clear(); rows.clear();
    if (r->mirror_rows.empty()) { // the common case: each occ row is its locus's only row
        for (const auto &o : r->occ) {
            loci.push_back(blind_locus_t{&o, &o, (uint32_t)rows.size(), 1});
            rows.push_back(&o);
        }
        return;
    }
    static thread_local aln_flat_index_t<aln_meth_locus_key_t> ix;
    static thread_local std::vector<std::pair<uint32_t, const aln_meth_occ_t *>> tmp;
    ix.clear(); tmp.clear();
    for (const auto &o : r->occ) {
        ix.insert(aln_meth_locus_key_t{o.rid, o.hyp, o.strand, o.alignment.rb}, loci.size());
        tmp.emplace_back((uint32_t)loci.size(), &o);
        loci.push_back(blind_locus_t{&o, &o, 0, 0});
    }
    for (const auto &o : r->mirror_rows) {
        const auto ins = ix.insert(aln_meth_locus_key_t{o.rid, o.hyp, o.strand, o.alignment.rb}, loci.size());
        if (ins.second) loci.push_back(blind_locus_t{nullptr, &o, 0, 0});
        const uint32_t k = (uint32_t)ix.val[ins.first];
        tmp.emplace_back(k, &o);
        if (mirror_free_better(o, *loci[k].best, s_mm)) loci[k].best = &o;
    }
    for (const auto &e : tmp) ++loci[e.first].n;
    uint32_t at = 0;
    for (auto &l : loci) { l.first = at; at += l.n; l.n = 0; }
    rows.resize(tmp.size());
    for (const auto &e : tmp) rows[loci[e.first].first + loci[e.first].n++] = e.second;
}
/* --meth-mapq conversion-blind: whether the mirror-free strata (meth_blind_strata, stock's selection
 * on mirror-free costs) can miss a row because of the cuts, decided where they first apply: the
 * search's next entry has score next, the best stratum is complete, the search has just lowered
 * its difference ceiling to c, and stopping says whether -R ends it here. It reads the search's
 * state and writes only blind_need and need_next, so the search is that of --meth-mapq scoring.
 * A row's converted score is at most its mirror-free cost, which is at most its cost, and the same
 * for its differences. Every row below next is resolved, under the read's own -n (the ceiling was
 * not lowered before), so the mirror-free best B_mf (at most B), the fewest mirror-free
 * differences D_mf among the loci at B_mf (each at its lowest mirror-free row, which the mirror
 * table holds) and the count of those loci are final. A row the strata still need costs at most
 * B_mf + M mirror-free, within the score cut (B + M), with at most min(-n, D_mf + 1) differences:
 * - next past B_mf + M: no such row is left, the strata are complete;
 * - -R stops the search here: a mirror-free run stops too when more than -R loci at B_mf are
 *   within its ceiling (its strata past that are partial, as stock's are), else they are unknown;
 * - otherwise the ceiling keeps every such row when min(-n, D_mf + 1) <= c, else it may not. */
static __attribute__((noinline)) void blind_note_need(aln_meth_slot_t *x, int next, int64_t c, bool stopping)
{
    const gap_opt_t &opt = x->job->opt;
    const bool gape = opt.mode & BWA_MODE_GAPE;
    auto mf_diffs = [&](const aln_meth_occ_t &o) { return (int64_t)o.n_mm - o.n_mirror + o.n_gapo + (gape ? o.n_gape : 0); };
    static thread_local std::vector<blind_locus_t> bl;
    static thread_local std::vector<const aln_meth_occ_t *> rows;
    static thread_local std::vector<std::pair<int64_t, int64_t>> loci; // (mirror-free cost, differences)
    blind_loci(x->res, opt.s_mm, bl, rows);
    loci.clear();
    int64_t best = INT64_MAX;
    for (const auto &l : bl) {
        loci.emplace_back(mirror_free_cost(l.best->cost, opt.s_mm, l.best->n_mirror), mf_diffs(*l.best));
        best = std::min(best, loci.back().first);
    }
    x->need_next = next;
    if (best == INT64_MAX || (int64_t)next > best + opt.s_mm) return;
    int64_t fewest = INT64_MAX;
    for (const auto &l : loci) if (l.first == best) fewest = std::min(fewest, l.second);
    if (stopping) {
        long n = 0;
        for (const auto &l : loci) n += l.first == best && l.second <= fewest + 1;
        x->blind_need = n <= opt.max_top2;
    } else x->blind_need = std::min<int64_t>(x->max_diff, fewest + 1) > c;
}
int aln_meth_boundary(aln_meth_slot_t *x, int next)
{
    gap_job_t *job = x->job;
    const gap_opt_t &opt = job->opt;
    job->score_cut = std::min(job->score_cut, gap_score_cut(x->best, opt.s_mm));
    if (x->best >= next) return 0;
    const bool conv = conv_tiebreak(x->grp->mo), gape = opt.mode & BWA_MODE_GAPE;
    auto diffs = [&](const aln_meth_occ_t &o) { return (int64_t)o.n_mm + o.n_gapo + (gape ? o.n_gape : 0); };
    int64_t fewest_d = INT64_MAX;
    for (const auto &o : x->res->occ) if (o.cost == x->best) fewest_d = std::min(fewest_d, diffs(o));
    const int64_t ceiling = fewest_d >= INT_MAX - 1 ? INT_MAX : fewest_d + 1;
    if (ceiling < job->max_diff) job->max_diff = (int)ceiling;
    // Count the best loci the finish keeps: those within the ceiling, and under the TAPS/5base
    // tie-break the fewest-conversion ones of them.
    auto kept = [&](const aln_meth_occ_t &o) { return o.cost == x->best && diffs(o) <= ceiling; };
    int fewest = INT_MAX;
    if (conv) for (const auto &o : x->res->occ) if (kept(o)) fewest = std::min(fewest, o.n_conv);
    long n = 0;
    for (const auto &o : x->res->occ) n += kept(o) && (!conv || o.n_conv == fewest);
    const bool stop = n > opt.max_top2;
    if (x->blind && x->need_next < 0) blind_note_need(x, next, job->max_diff, stop);
    return stop;
}
void aln_meth_hit(aln_meth_slot_t *x, const gap_entry_t *e, bwtint_t k, bwtint_t l)
{
    // Reaching the same interval along a path with the same gaps adds no new candidate rows; a path
    // with other gaps (another placement in a tandem repeat, say) can cost differently.
    aln_meth_pend_t h;
    h.k = k; h.l = l; h.e = *e;
    for (int g = 0; g < ALN_METH_MAX_GAPS; ++g) h.gap[g] = 0;
    entry_gaps(x, *e, h.gap);
    aln_meth_iv_key_t key{k, l, (uint32_t)e->n_gapo, {h.gap[0], h.gap[1], h.gap[2], h.gap[3]}};
    if (!x->intervals.insert(key, 0).second) return;
    x->pend.push_back(h);
}
const ubyte_t *aln_meth_slot_start(aln_meth_slot_t *x, bwa_seq_t *p, aln_meth_read_t *res, int64_t ord, int max_diff, int has_copy)
{
    x->p = p; x->res = res; x->len = p->len; x->max_diff = max_diff;
    x->keep = has_copy ? ALN_METH_KEEP_COPIES : x->grp->keep_strata ? ALN_METH_KEEP_STRATA : 0;
    *res = aln_meth_read_t();
    x->alt_bound = aln_meth_alt_bound_t();
    x->rng = aln_meth_draw_seed(x->grp->ref, x->grp->read_no, ord);
    x->vseq.resize(p->len); x->rd_fwd.resize(p->len); x->rd_rev.resize(p->len);
    aln_meth_convert(x->grp->read_no, p->seq, p->len, x->vseq.data());
    // The search takes any code above 3 as ambiguous, but the realigner accepts only 0-4: aln
    // reads '-' as 5, so the refinement queries carry it as 4 (N). p->seq and p->rseq keep it.
    for (int k = 0; k < (int)p->len; ++k) {
        x->rd_fwd[k] = std::min<ubyte_t>(p->seq[p->len - 1 - k], 4);
        x->rd_rev[k] = std::min<ubyte_t>(p->rseq[k], 4);
    }
    x->pend.clear(); x->windows.clear(); x->intervals.clear(); x->loci.clear(); x->garena.clear();
    x->mirror_loci.clear(); x->mirror_bytes = 0;
    x->sa_ix.clear(); x->sa_rows.clear();
    x->bits_ready = 0;
    x->alignment_bytes = 0; x->unrep_cost = INT64_MAX;
    x->final_flush = x->stop = 0;
    x->flushed_score = -1;
    const gap_opt_t *o = x->grp->opt;
    x->cut = !(o->mode & BWA_MODE_NONSTOP) && o->s_mm >= 0 && o->s_gapo >= 0 && o->s_gape >= 0 && !x->grp->ref->known_snvs &&
        !x->grp->ref->no_cut;
    x->best = INT64_MAX; x->flushed = 0;
    x->blind = x->grp->mo && meth_mapq_conversion_blind(&x->grp->mo->m);
    x->need_next = -1; x->blind_need = 0;
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
 * candidate could have competed (cost within the retained window) is marked incomplete.
 * win is the alignment's reference window as aln_meth_fetch_window returns it. */
static bool account_alignment(const aln_meth_opt_t *mo, const gap_opt_t *opt, const std::vector<uint8_t> &win,
                              const ubyte_t *query, int len, aln_meth_occ_t *o)
{
    int m;
    const int8_t *mat = meth_mat(mo, o->hyp, &m);
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
                if (read > 3 || mat[base * m + read] < 0) {
                    ++o->n_mm;
                    // the mirror cell: the converted base read over a reference base that is
                    // already it (OT: read C over T; OB: read G over A)
                    if (read == (o->hyp ? 1 : 2) && meth_code_base(base) == (o->hyp ? 3 : 0)) ++o->n_mirror;
                }
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
// account_alignment on the window [rb, re) of o's contig, fetched here.
static bool account_alignment(const aln_meth_ref_t *ref, const aln_meth_opt_t *mo,
                              const gap_opt_t *opt, int64_t rb, int64_t re,
                              const ubyte_t *query, int len, aln_meth_occ_t *o, int64_t site = -1)
{
    std::vector<uint8_t> win;
    fetch_window(ref, mo, rb, re - rb, o->rid, win, site);
    return account_alignment(mo, opt, win, query, len, o);
}
static void resolve_one(aln_meth_slot_t *x, const aln_meth_candidate_t &candidate, const ubyte_t *query,
                        int64_t site, int n_ins, int n_del, bool bounded = false);

/* ---------------------------------------------------------------------------------------------
 * Path cost: a hit's rows are scored on the original reference along the search path that found
 * them, as stock aln scores a hit by its path (no window DP in the search). The converted search
 * sees a converted mismatch wherever the original reference has one, and the gap terms are the
 * same, so a row's strict cost is never below its path's converted score.
 * ------------------------------------------------------------------------------------------- */

// The path's gaps in query orientation (aln_meth_layout_t) from its gap events (ALN_METH_GAP),
// whose indexes are into the complemented read the job steps over, S. S is the query of strand 1
// (rd_rev, from p->rseq) base for base; strand 0's query (rd_fwd) is S reversed.
static void meth_layout(const uint32_t *gap, int n_gapo, int len, int strand, aln_meth_layout_t *L)
{
    L->n = n_gapo;
    for (int g = 0; g < n_gapo; ++g) {
        const int i = (int)(gap[g] & 0xffff), gl = (int)(gap[g] >> 16 & 0x7fff), del = (int)(gap[g] >> 31);
        // S index: an insertion's first base in S order; for a deletion, the S bases before it
        const int s0 = del ? i + 1 : i - gl + 1;
        L->qpos[g] = strand ? s0 : del ? len - s0 : len - s0 - gl;
        L->len[g] = gl; L->del[g] = del;
    }
    for (int a = 1; a < L->n; ++a) // ascending qpos (at most ALN_METH_MAX_GAPS)
        for (int b = a; b > 0 && L->qpos[b - 1] > L->qpos[b]; --b) {
            std::swap(L->qpos[b - 1], L->qpos[b]); std::swap(L->len[b - 1], L->len[b]); std::swap(L->del[b - 1], L->del[b]);
        }
}
static int layout_span(const aln_meth_layout_t &L, int len)
{
    int span = len;
    for (int g = 0; g < L.n; ++g) span += L.del[g] ? L.len[g] : -L.len[g];
    return span;
}
// The layout as a raw CIGAR (BAM op codes); false when it is not a valid alignment of len bases.
static bool layout_cigar(const aln_meth_layout_t &L, int len, std::vector<uint32_t> *c)
{
    c->clear();
    int q = 0;
    for (int g = 0; g < L.n; ++g) {
        if (L.qpos[g] < q || L.len[g] <= 0) return false;
        if (L.qpos[g] > q) c->push_back((uint32_t)(L.qpos[g] - q) << 4);
        c->push_back((uint32_t)L.len[g] << 4 | (L.del[g] ? 2u : 1u));
        q = L.qpos[g] + (L.del[g] ? 0 : L.len[g]);
    }
    if (q > len) return false;
    if (q < len) c->push_back((uint32_t)(len - q) << 4);
    return true;
}
// The layout of a raw CIGAR of M/I/D operations (false: another operation, or too many gaps).
static bool cigar_layout(const uint32_t *c, size_t n, aln_meth_layout_t *L)
{
    L->n = 0;
    int q = 0;
    for (size_t k = 0; k < n; ++k) {
        const int op = c[k] & 15, l = (int)(c[k] >> 4);
        if (op == 0) { q += l; continue; }
        if ((op != 1 && op != 2) || L->n == ALN_METH_MAX_GAPS) return false;
        L->qpos[L->n] = q; L->len[L->n] = l; L->del[L->n] = op == 2; ++L->n;
        if (op == 1) q += l;
    }
    return true;
}
// Display normalization of a raw global CIGAR, as aln_meth_realign_window does it: a terminal
// insertion is soft-clipped, a terminal deletion is dropped from the reference span.
static void normalize_display(aln_meth_alignment_t *a)
{
    a->cigar = a->raw_cigar;
    if ((a->cigar.back() & 15) == 1) a->cigar.back() = (a->cigar.back() & ~15u) | 4;
    if ((a->cigar.front() & 15) == 1) a->cigar.front() = (a->cigar.front() & ~15u) | 4;
    if ((a->cigar.back() & 15) == 2) { a->re -= a->cigar.back() >> 4; a->cigar.pop_back(); }
    if (!a->cigar.empty() && (a->cigar.front() & 15) == 2) { a->rb += a->cigar.front() >> 4; a->cigar.erase(a->cigar.begin()); }
}
// The difference counts of the path's gap events made inside the seed, as the search counts them
// against -k: a difference is constrained when the step that makes it is at an index above the
// seed's offset (bwt_match_gap's seed test runs for ii = i - seed_off > 0).
static int seed_gap_diffs(const uint32_t *gap, int n_gapo, int seed_off, bool gape)
{
    int d = 0;
    for (int g = 0; g < n_gapo; ++g) {
        const int i = (int)(gap[g] & 0xffff), gl = (int)(gap[g] >> 16 & 0x7fff), del = (int)(gap[g] >> 31);
        if (i > seed_off) ++d; // the open
        if (!gape) continue;
        if (del) d += i > seed_off ? gl - 1 : 0; // a deletion's extensions step at the same index
        else for (int s = i - 1; s > i - gl; --s) d += s > seed_off; // an insertion's at i-1, i-2, ...
    }
    return d;
}
// Whether query position q (orientation ori: 0 rd_fwd, 1 rd_rev) is inside the seed's constrained
// part (S index above seed_off).
static inline bool seed_pos(int q, int ori, int len, int seed_off) { return ori ? q > seed_off : q < len - 1 - seed_off; }
// The read's bit planes for orientation ori, with the seed mask when the read is seeded.
static const aln_meth_read_bits_t &slot_bits(aln_meth_slot_t *x, int ori)
{
    if (!(x->bits_ready >> ori & 1)) {
        const ubyte_t *query = ori ? x->rd_rev.data() : x->rd_fwd.data();
        aln_meth_read_bits_t &b = x->bits[ori];
        aln_meth_read_bits(query, x->len, &b);
        b.seed.clear();
        if (x->job && x->job->seed_width) {
            b.seed.assign(b.code.size(), 0);
            for (int q = 0; q < x->len; ++q)
                if (seed_pos(q, ori, x->len, x->job->c.seed_off)) b.seed[q >> 5] |= 1ULL << (62 - 2 * (q & 31));
        }
        x->bits_ready |= 1 << ori;
    }
    return x->bits[ori];
}
/* Scores path L of the read at window [rb, rb + span) of contig rid under hypothesis hyp, strand
 * strand, with site allele site alone in place (--known-snvs; -1: the reference), into o (counts,
 * cost, CIGARs, start); *seed_mm gets its strict mismatches in the seed's constrained part (0 for an
 * unseeded read), and *seed_mirror the mirror-cell mismatches among them. The bit-parallel blocks serve the 5-letter model on the reference; otherwise the
 * window is fetched (with that site's ALT, CpG labels) and counted base by base as
 * account_alignment counts. */
static bool score_path(aln_meth_slot_t *x, const aln_meth_layout_t &L, int rid, int hyp, int strand, int64_t rb,
                       int span, int64_t site, aln_meth_occ_t *o, int *seed_mm, int *seed_mirror)
{
    aln_meth_group_t *g = x->grp;
    const int len = x->len, ori = strand ? 1 : 0;
    o->rid = rid; o->hyp = hyp; o->strand = strand; o->ks_site = site;
    o->alignment.rb = rb; o->alignment.re = rb + span;
    if (!layout_cigar(L, len, &o->alignment.raw_cigar)) return false;
    *seed_mm = *seed_mirror = 0;
    const bool seeded = x->job && x->job->seed_width;
    if (site < 0 && g->mo->fast.ok) {
        const aln_meth_read_bits_t &rd = slot_bits(x, ori);
        aln_meth_amb_bits(*g->ref, rb, span, &x->amb);
        aln_meth_block_counts_t c;
        int q = 0, w = 0;
        for (int k = 0; k <= L.n; ++k) {
            const int end = k < L.n ? L.qpos[k] : len;
            if (end > q) aln_meth_score_block(*g->ref, g->mo->fast, hyp, rb + w, q, end - q, rd, x->amb.data(), w, &c);
            w += end - q; q = end;
            if (k == L.n) break;
            ++o->n_gapo; o->n_gape += L.len[k] - 1;
            if (L.del[k]) { o->n_del += L.len[k]; w += L.len[k]; }
            else { o->n_ins += L.len[k]; q += L.len[k]; }
        }
        o->n_mm = c.n_mm; o->n_conv = c.n_conv; o->n_mirror = c.n_mirror;
        *seed_mm = c.n_seed_mm; *seed_mirror = c.n_seed_mirror;
    } else {
        std::vector<uint8_t> &win = x->win;
        fetch_window(g->ref, g->mo, rb, span, rid, win, site);
        const ubyte_t *query = ori ? x->rd_rev.data() : x->rd_fwd.data();
        int m;
        const int8_t *mat = meth_mat(g->mo, hyp, &m);
        int q = 0, t = 0;
        for (uint32_t item : o->alignment.raw_cigar) {
            const int op = item & 15, n = (int)(item >> 4);
            if (op == 0) {
                for (int k = 0; k < n; ++k, ++q, ++t) {
                    const int read = query[q], base = win[t];
                    if (read > 3 || mat[base * m + read] < 0) {
                        ++o->n_mm;
                        // the mirror cell, as account_alignment counts it
                        const bool mirror = read == (hyp ? 1 : 2) && meth_code_base(base) == (hyp ? 3 : 0);
                        o->n_mirror += mirror;
                        if (seeded && seed_pos(q, ori, len, x->job->c.seed_off)) {
                            ++*seed_mm;
                            *seed_mirror += mirror;
                        }
                    } else if (read != meth_code_base(base)) ++o->n_conv;
                }
            } else {
                ++o->n_gapo; o->n_gape += n - 1;
                if (op == 1) { o->n_ins += n; q += n; } else { o->n_del += n; t += n; }
            }
        }
    }
    o->cost = (int64_t)g->opt->s_mm * o->n_mm + (int64_t)g->opt->s_gapo * o->n_gapo + (int64_t)g->opt->s_gape * o->n_gape;
    normalize_display(&o->alignment);
    o->pos = o->alignment.rb;
    return true;
}
// bounded (--known-snvs, a window the ALT-window search proposed): the read's own search never
// proposes an alignment over its difference limit, so neither may this path. The bound is the
// gap search's: n_mm + n_gapo (+ n_gape under the default gap-extension mode) <= the read's
// max_diff (x->max_diff, the limit store_row also applies), at most max_gapo gap opens and
// max_gape gap extensions.
// mirror_free: the mirror mismatches count as conversions, as --meth-mapq conversion-blind takes them.
static bool within_read_limits(const aln_meth_slot_t *x, const aln_meth_occ_t &o, bool mirror_free = false)
{
    const gap_opt_t *op = x->grp->opt;
    const int diffs = o.n_mm - (mirror_free ? o.n_mirror : 0) + o.n_gapo + ((op->mode & BWA_MODE_GAPE) ? o.n_gape : 0);
    return diffs <= x->max_diff && o.n_gapo <= std::max(0, op->max_gapo) && o.n_gape <= std::max(0, op->max_gape);
}
// --meth-mapq conversion-blind: whether row a is a better row of its locus than b with their
// mirror mismatches free. Under collapsed scoring a mirror mismatch is a conversion and a row's
// CIGARs (its search path's) are the same, so this is store_row's order on mirror-free costs:
// the cost, the site allele, then the CIGARs. The TAPS conversion count is left out
// (conversion-blind is EM-seq only).
static bool mirror_free_better(const aln_meth_occ_t &a, const aln_meth_occ_t &b, int s_mm)
{
    const int64_t ca = mirror_free_cost(a.cost, s_mm, a.n_mirror), cb = mirror_free_cost(b.cost, s_mm, b.n_mirror);
    if (ca != cb) return ca < cb;
    if (a.ks_site != b.ks_site) return (uint64_t)(a.ks_site + 1) < (uint64_t)(b.ks_site + 1);
    if (a.alignment.cigar != b.alignment.cigar) return a.alignment.cigar < b.alignment.cigar;
    if (a.alignment.re != b.alignment.re) return a.alignment.re < b.alignment.re;
    return a.alignment.raw_cigar < b.alignment.raw_cigar;
}
// --meth-mapq conversion-blind: notes row o in the mirror table (aln_meth_read_t::mirror_rows)
// when it is within the read's limits with its mirror mismatches free, as collapsed scoring
// applies them: -n, -k, a record's fields and, for a bounded window, within_read_limits. A copy is
// kept; o is not changed, so store_row goes on as under scoring. The table keeps, per locus, the
// best row by mirror_free_better. Without the stock cuts only rows with a mirror mismatch need an
// entry: a row without one costs its genomic cost mirror-free, so the row that beats it in occ is
// no worse mirror-free (meth_blind_strata takes the better of the two). Under the cuts the
// mirror-free strata apply a difference ceiling, and a mirror-free run keeps a locus's best row
// among those it reaches (meth_blind_strata), which can be a row without a mirror mismatch that
// occ does not hold: store_row then also passes the row its locus's occ row beats or replaces,
// and the table keeps the best per locus and number of mirror-free differences. A row store_row
// drops otherwise is not needed: without a mirror mismatch it is over the read's limits or past
// the best + M on both costs.
// The table has a CIGAR budget of its own; past it the strata are unknown (blind_partial).
static __attribute__((noinline)) void store_mirror_row(aln_meth_slot_t *x, const aln_meth_occ_t &o, int seed_diffs,
                                                       bool bounded, int seed_mirror)
{
    const gap_opt_t *op = x->grp->opt;
    aln_meth_read_t *r = x->res;
    const int64_t differences = (int64_t)o.n_mm - o.n_mirror + o.n_gapo + ((op->mode & BWA_MODE_GAPE) ? o.n_gape : 0);
    const bool seeded = x->job && x->job->seed_width;
    if (differences > x->max_diff || (seeded && seed_diffs - seed_mirror > op->max_seed_diff)) return;
    if (mirror_free_cost(o.cost, op->s_mm, o.n_mirror) > INT_MAX || o.n_mm - o.n_mirror > 255
        || o.n_gapo > 255 || o.n_gape > 255) return;
    if (bounded && !within_read_limits(x, o, true)) return;
    if (r->blind_partial) return;
    // Under the cut the mirror-free strata keep rows within M of the mirror-free best, which is at
    // most the best so far (store_row's drop, on mirror-free costs).
    if (x->cut && x->best != INT64_MAX && mirror_free_cost(o.cost, op->s_mm, o.n_mirror) > x->best + op->s_mm) return;
    const aln_meth_mirror_key_t key{aln_meth_locus_key_t{o.rid, o.hyp, o.strand, o.alignment.rb},
                                    x->cut ? (int)differences : 0};
    size_t pos;
    const int32_t at = x->mirror_loci.find(key, &pos);
    if (at >= 0 && !mirror_free_better(o, r->mirror_rows[x->mirror_loci.val[at]], op->s_mm)) return;
    const uint64_t bytes = (o.alignment.raw_cigar.size() + o.alignment.cigar.size()) * sizeof(uint32_t);
    if (bytes > (64ULL << 20) - x->mirror_bytes) {
        r->blind_partial = 1;
        std::vector<aln_meth_occ_t>().swap(r->mirror_rows);
        return;
    }
    x->mirror_bytes += bytes;
    if (at >= 0) r->mirror_rows[x->mirror_loci.val[at]] = o;
    else {
        x->mirror_loci.insert(key, r->mirror_rows.size());
        r->mirror_rows.push_back(o);
    }
}
// Applies the strict -n/-k filters and stores o as its locus's candidate (the cheaper wins).
// seed_mirror: the mirror-cell mismatches among seed_diffs.
static void store_row(aln_meth_slot_t *x, aln_meth_occ_t &o, int seed_diffs, bool bounded = false, int seed_mirror = 0)
{
    aln_meth_group_t *g = x->grp;
    // --meth-mapq conversion-blind observes the rows with a mirror mismatch before the filters
    // below (and, under the cuts, the rows a locus's occ row beats or replaces).
    if (o.n_mirror && x->blind) store_mirror_row(x, o, seed_diffs, bounded, seed_mirror);
    // Stock aln never reports a hit with more differences than the read's -n allows, or more in its
    // seed than -k allows. The converted search cannot see a read C over a reference T, so it can
    // reach a row within -n (-k) converted differences whose original-reference path has more; drop
    // it before it is a locus (no X0/X1, MAPQ or pairing). Differences are counted as stock aln
    // counts them: mismatches plus gap opens, plus gap extensions unless -e is positive.
    const bool gape = g->opt->mode & BWA_MODE_GAPE;
    const int64_t differences = (int64_t)o.n_mm + o.n_gapo + (gape ? o.n_gape : 0);
    const bool seeded = x->job && x->job->seed_width;
    const bool over_diff = differences > x->max_diff;
    const bool over_seed = !over_diff && seed_diffs > g->opt->max_seed_diff && seeded;
    if (over_diff || over_seed) {
        ++(over_diff ? x->res->n_over_diff : x->res->n_over_seed);
        return;
    }
    if (!(o.cost <= INT_MAX && o.n_mm <= 255 && o.n_gapo <= 255 && o.n_gape <= 255)) {
        ++x->res->n_unrepresentable;
        x->unrep_cost = std::min(x->unrep_cost, o.cost);
        return;
    }
    if (bounded && !within_read_limits(x, o)) return;
    ++x->res->n_rescored;
    // Under the cut a row above the current best + M can never be retained (the best only falls).
    // The mirror table saw the row first (a row without a mirror mismatch costs as much mirror-free).
    if (x->cut && x->best != INT64_MAX && o.cost > x->best + g->opt->s_mm) return;
    x->best = std::min(x->best, o.cost);
    const uint64_t bytes = (o.alignment.raw_cigar.size() + o.alignment.cigar.size()) * sizeof(uint32_t);
    // Cumulative CIGAR bytes conservatively bound ownership even when a locus is replaced.
    if (bytes > (64ULL << 20) - x->alignment_bytes) {
        x->res->flags |= ALN_METH_BIND_CIGAR;
        x->stop = ALN_METH_STOP_CIGAR;
        return;
    }
    x->alignment_bytes += bytes;
    const auto ins = x->loci.insert(aln_meth_locus_key_t{o.rid, o.hyp, o.strand, o.alignment.rb}, x->res->occ.size());
    if (ins.second) x->res->occ.push_back(std::move(o));
    else {
        const size_t at = x->loci.val[ins.first];
        const aln_meth_occ_t &cur = x->res->occ[at];
        // --known-snvs: the cheaper wins; on equal cost the reference, then the lower
        // (position, ALT) site allele (aln_site_key order, -1 first). Two rows of the same allele
        // (the reference, or one site allele) compare as reference rows do.
        const bool better = o.ks_site == cur.ks_site ? alignment_better(o, cur, conv_tiebreak(g->mo))
                          : o.cost != cur.cost ? o.cost < cur.cost
                          : (uint64_t)(o.ks_site + 1) < (uint64_t)(cur.ks_site + 1);
        // --meth-mapq conversion-blind under the cuts: the losing row of a locus, if it has no
        // mirror mismatch (the table already holds those), for the mirror-free strata
        // (store_mirror_row). It passed the read's limits, so it passes them mirror-free too.
        const aln_meth_occ_t &loser = better ? cur : o;
        if (x->cut && x->blind && !loser.n_mirror) store_mirror_row(x, loser, 0, false, 0);
        if (better) x->res->occ[at] = std::move(o);
    }
}
// Scores a row along path L as the reference and, under --known-snvs, once per site allele in its
// window with that ALT alone in place (as resolve_with_sites does for a refined window); a site the
// row does not cover leaves it as the reference's. false: L is not an alignment of the read.
static bool score_row_sites(aln_meth_slot_t *x, const aln_meth_layout_t &L, int rid, int hyp, int strand,
                            int64_t rb, int span, int seed_gaps, bool bounded)
{
    aln_meth_group_t *g = x->grp;
    auto one = [&](int64_t site) {
        aln_meth_occ_t o;
        int seed_mm = 0, seed_mirror = 0;
        if (!score_path(x, L, rid, hyp, strand, rb, span, site, &o, &seed_mm, &seed_mirror)) return false;
        if (site >= 0) {
            o.alt = ks_site_rel(g->ref->known_snvs, g->ref->bns, rid, o.alignment.rb, site);
            if (aln_site_key_pos(site) < o.alignment.rb || aln_site_key_pos(site) >= o.alignment.re || !o.alt) return true;
        }
        store_row(x, o, seed_mm + seed_gaps, bounded, seed_mirror);
        return true;
    };
    if (!one(-1)) return false;
    if ((x->res->flags & ALN_METH_BIND_CIGAR) || !g->ref->known_snvs) return true;
    std::vector<int64_t> keys;
    aln_meth_window_sites(*g->ref, rb, span, &keys);
    for (const int64_t key : keys) {
        if (x->res->flags & ALN_METH_BIND_CIGAR) break; // the read's CIGAR storage is full
        one(key);
    }
    return true;
}
/* --known-snvs: resolves a candidate window against the reference, then once per site allele in it
 * (aln_meth_window_sites) with that ALT alone in place; its locus keeps the cheapest (see
 * resolve_one). */
static void resolve_with_sites(aln_meth_slot_t *x, const aln_meth_candidate_t &candidate, const ubyte_t *query,
                               int n_ins, int n_del, bool bounded = false)
{
    resolve_one(x, candidate, query, -1, n_ins, n_del, bounded);
    if ((x->res->flags & ALN_METH_BIND_CIGAR) || !x->grp->ref->known_snvs) return;
    std::vector<int64_t> keys;
    aln_meth_window_sites(*x->grp->ref, candidate.rb, (int)(candidate.re - candidate.rb), &keys);
    for (const int64_t key : keys) {
        if (x->res->flags & ALN_METH_BIND_CIGAR) return; // the read's CIGAR storage is full
        resolve_one(x, candidate, query, key, n_ins, n_del, bounded);
    }
}
// Scores one row of hit h (SA value sa) along h's path, as the reference and each site allele.
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
    // --known-snvs, the ALT-window search: every distinct read-length window reaches here once or
    // more, before exact-window dedup, refinement and locus collapse can drop or move it. One
    // covering its own window's site is an ungapped ALT placement with the history's mismatches
    // (a gapped history of equal span does not give the window's mismatch count: bound it by 0).
    if (x->res->site_bound.valid && span == x->len) {
        xassert(rid >= 0 && (size_t)rid < g->ref->window_site.size(), "--known-snvs window without a recorded site");
        const int64_t off = (int64_t)(pos - ca.offset), site = g->ref->window_site[rid];
        if (site >= off && site < off + span)
            x->res->site_bound.value = std::min(x->res->site_bound.value, h.e.n_gapo ? 0 : (int)h.e.n_mm);
    }
    const int64_t rb = g->ref->bns->anns[rid].offset + (pos - ca.offset);
    aln_meth_win_key_t key{rid, hyp, strand, rb, rb + span, aln_meth_layout_t()};
    meth_layout(h.gap, h.e.n_gapo, x->len, strand, &key.lay);
    if (layout_span(key.lay, x->len) != span) {
        fprintf(stderr, "ERROR: --meth search path's gaps (%d) do not give its reference span\n", (int)h.e.n_gapo); exit(1);
    }
    if (!x->windows.insert(key, 0).second) return;
    // A 65535-base read with a deletion spans more than a window holds. Drop the window as
    // unrepresentable, with cost 0 as its bound: any retained placement is then incomplete.
    if (span > ALN_MAX_READ_LEN) {
        ++x->res->n_unrepresentable;
        x->unrep_cost = 0;
        return;
    }
    const int seed_gaps = x->job && x->job->seed_width
        ? seed_gap_diffs(h.gap, h.e.n_gapo, x->job->c.seed_off, g->opt->mode & BWA_MODE_GAPE) : 0;
    if (!score_row_sites(x, key.lay, rid, hyp, strand, rb, span, seed_gaps, false)) {
        fprintf(stderr, "ERROR: --meth search path (%d gaps: %#x %#x %#x %#x) is not an alignment of the %d-base read\n",
                (int)h.e.n_gapo, h.gap[0], h.gap[1], h.gap[2], h.gap[3], x->len);
        exit(1);
    }
}
// A --known-snvs candidate window without a path (its CIGAR is not one): refined by the window DP,
// as before path costs, and filtered by store_row without -k (it has no seed).
static void resolve_one(aln_meth_slot_t *x, const aln_meth_candidate_t &candidate, const ubyte_t *query,
                        int64_t site, int n_ins, int n_del, bool bounded)
{
    aln_meth_group_t *g = x->grp;
    // Realigns the window (band_limit as aln_meth_realign_window takes it) into a fresh o, with
    // fits whether its counts fit a record; false when the result misses its site.
    aln_meth_occ_t o;
    bool fits = false;
    // The accounting reads the same window the refinement fetched (same window, contig and site).
    std::vector<uint8_t> &win = x->win;
    auto evaluate = [&](int band_limit) {
        o = aln_meth_occ_t();
        o.rid = candidate.rid; o.hyp = candidate.hyp; o.strand = candidate.strand; o.ks_site = site;
        const auto status = aln_meth_realign_window(*g->ref, *g->mo, candidate, query, x->len,
                                                   n_ins, n_del, &o.alignment, site, &win, band_limit);
        // Invalid/resource/empty normalized alignments must never be silently omitted competitors.
        if (status != aln_meth_eval_status::ok) {
            fprintf(stderr, "ERROR: --meth cannot evaluate a generated candidate window; stopping\n");
            exit(1);
        }
        if (site >= 0) { // --known-snvs: a site the alignment does not cover leaves it as the reference's
            o.alt = ks_site_rel(g->ref->known_snvs, g->ref->bns, candidate.rid, o.alignment.rb, site);
            if (aln_site_key_pos(site) < o.alignment.rb || aln_site_key_pos(site) >= o.alignment.re || !o.alt) return false;
        }
        fits = account_alignment(g->mo, g->opt, win, query, x->len, &o);
        return true;
    };
    if (!evaluate(-1)) return;
    // A bounded window whose best-scoring alignment breaks the read's limits (within_read_limits)
    // is realigned once more with its indels held to the window's own extent (n_ins, n_del:
    // ungapped for a window of the read's length); store_row keeps it only if that alignment is
    // within the limits.
    if (bounded && fits && !within_read_limits(x, o)) {
        const int narrow = std::max(n_ins, n_del);
        if (o.alignment.band > narrow && !evaluate(narrow)) return;
    }
    if (!fits) { // counts past a record's fields
        ++x->res->n_unrepresentable;
        x->unrep_cost = std::min(x->unrep_cost, o.cost);
        return;
    }
    o.pos = o.alignment.rb;
    store_row(x, o, 0, bounded);
}
void aln_meth_flush_parked(aln_meth_group_t *g)
{
    const int64_t budget = std::max<int64_t>(0, (int64_t)g->mo->p_best + g->mo->p_second);
    if (budget > 65536) { fprintf(stderr, "ERROR: --meth lookup budget exceeds 65536 rows per read\n"); exit(1); }
    for (int i = 0; i < g->n_slots; ++i) {
        if (g->cores[i].state != GAP_SLOT_PARKED) continue;
        aln_meth_slot_t &x = g->slot[i];
        ++x.res->n_flushes;
        x.flushed = 1;
        // A flush whose pending rows exceed what is left of the budget binds it: the read ends
        // incomplete (ALN_METH_BIND_BEST) whatever those rows resolve to, keeping only a MAPQ-0
        // best-known placement. Resolve just ALN_METH_BIND_ROWS more rows for that placement
        // rather than the rest of the budget. Earlier flushes are unchanged, so the same reads
        // end incomplete.
        // An interval already resolved (reached again along a path with other gaps) reuses its SA
        // rows: it costs no lookups, so it does not count toward the budget.
        auto cached = [&](const aln_meth_pend_t &h) {
            size_t pos;
            return x.sa_ix.find(aln_meth_kl_key_t{h.k, h.l}, &pos);
        };
        // Two paths can reach one uncached interval in the same flush: its rows are looked up once.
        static thread_local std::vector<std::pair<bwtint_t, bwtint_t>> fresh;
        fresh.clear();
        for (const auto &h : x.pend) if (cached(h) < 0) fresh.emplace_back(h.k, h.l);
        std::sort(fresh.begin(), fresh.end());
        fresh.erase(std::unique(fresh.begin(), fresh.end()), fresh.end());
        int64_t limit = budget, pending = 0;
        for (const auto &kl : fresh) pending += (int64_t)(kl.second - kl.first + 1);
        if (pending > budget - x.res->n_lookups)
            limit = std::min<int64_t>(budget, x.res->n_lookups + ALN_METH_BIND_ROWS);
        for (const auto &h : x.pend) {
            if (x.stop) break;
            const bwtint_t width = h.l - h.k + 1;
            const int32_t c = cached(h);
            if (c >= 0) {
                const size_t off = x.sa_ix.val[c];
                for (bwtint_t r = 0; r < width && !x.stop; ++r) resolve_window(&x, h, x.sa_rows[off + r]);
                continue;
            }
            const int64_t remaining = limit - x.res->n_lookups;
            const int n = (int)std::min<uint64_t>(width, (uint64_t)std::max<int64_t>(0, remaining));
            const size_t off = x.sa_rows.size();
            x.sa_rows.resize(off + n);
            bwtint_t *rows = x.sa_rows.data() + off;
            for (int k = 0; k < n; ++k) rows[k] = h.k + k;
            bwt_sa_batch(g->bwt, rows, n, rows);
            if ((bwtint_t)n == width) x.sa_ix.insert(aln_meth_kl_key_t{h.k, h.l}, off);
            for (int k = 0; k < n && !x.stop; ++k) { ++x.res->n_lookups; resolve_window(&x, h, x.sa_rows[off + k]); }
            // A CIGAR stop on this interval's looked-up rows is final: keep it over the lookup
            // stop, which aln_meth_slot_finish lifts to resolve the ALT-window candidates.
            if ((bwtint_t)n < width) {
                x.res->flags |= ALN_METH_BIND_BEST;
                if (!x.stop) x.stop = ALN_METH_STOP_LOOKUP;
            }
        }
        x.pend.clear();
        g->cores[i].state = GAP_SLOT_SEARCH;
    }
    g->n_parked = 0;
}
// The reservoir draw among the best-cost loci (the fewest-conversion ones once conversions
// displaced an equal-cost locus); -1 when there are none. One body for the search's occ and a
// finished read's hits, so a duplicate copy's draw over hits matches the search's draw.
template <class Locus>
static int draw_main(const std::vector<Locus> &loci, int n_conv_moved, uint64_t rng)
{
    int64_t best = INT64_MAX;
    for (const auto &o : loci) best = std::min<int64_t>(best, o.cost);
    int fewest = INT_MAX;
    if (n_conv_moved)
        for (const auto &o : loci) if (o.cost == best) fewest = std::min(fewest, (int)o.n_conv);
    int held = -1, count = 0;
    for (int i = 0; i < (int)loci.size(); ++i) {
        if (loci[i].cost != best || (n_conv_moved && loci[i].n_conv != fewest)) continue;
        ++count;
        if (count == 1 || (splitmix64(&rng) % (uint64_t)count) == 0) held = i;
    }
    return held;
}
// --meth-mapq conversion-blind: collapsed scoring's order of a read's loci (slot_finish's sort).
static bool blind_collapsed_order(const aln_meth_blind_row_t &a, const aln_meth_blind_row_t &b)
{
    return std::make_tuple(a.rid, a.hyp, a.strand, a.pos) < std::make_tuple(b.rid, b.hyp, b.strand, b.pos);
}
// --meth-mapq conversion-blind: the index of the main hit collapsed scoring draws among the n
// mirror-free best loci (in blind_collapsed_order): draw_main over them with the same seed.
// collapsed draws over all its loci in that order, but only a best-cost locus takes a draw, and
// conversion-blind is EM-seq only (no TAPS conversion tie-break). The draw reads only how many
// loci there are, not the loci. -1 when there are none.
static int blind_draw(size_t n, uint64_t rng)
{
    int held = -1;
    for (int i = 0; i < (int)n; ++i) // draw_main's loop over loci that all take a draw
        if (i == 0 || (splitmix64(&rng) % (uint64_t)(i + 1)) == 0) held = i;
    return held;
}
// Kept out of line: its one caller is slot_finish, which the search kernel inlines, and LTO
// inlining it there changed the kernel's code (clang 19).
__attribute__((noinline)) void aln_meth_redraw(aln_meth_read_t *res, uint64_t rng)
{
    if (res->main < 0) return;
    res->main = draw_main(res->occ, res->n_conv_moved, rng);
}
// Whether aln_meth_trim cuts a read with n loci to its main hit.
static bool trims(const aln_meth_read_t *res, size_t n, int max_xa_hits)
{
    return res->flags || (int64_t)n > (int64_t)max_xa_hits + 1;
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
    h.alt = (uint8_t)o.alt; // ks_site_rel's index, relative to h.pos == o.alignment.rb
    h.n_mirror = (uint8_t)o.n_mirror; // <= n_mm <= 255
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
    // Free the search-only raw CIGARs first, so they are not held alongside the packed hits.
    for (auto &o : res->occ) std::vector<uint32_t>().swap(o.alignment.raw_cigar);
    aln_meth_store_loci(res);
}
void aln_meth_trim(aln_meth_read_t *res, int max_xa_hits)
{
    if (res->main < 0) return;
    // A read not yet stored stores only its main locus. That locus has the best cost, so an
    // incomplete read's best-cost filter keeps it.
    if (!res->occ.empty() && trims(res, res->occ.size(), max_xa_hits)) {
        aln_meth_occ_t main = std::move(res->occ[res->main]);
        res->occ.clear(); res->occ.push_back(std::move(main)); res->main = 0;
    }
    meth_release_finished(res);
    // A stored read (a duplicate copy, after its draw) trims its hits.
    if (res->hits.size() > 1 && trims(res, res->hits.size(), max_xa_hits)) {
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
// aln_meth_redraw's draw over a finished read's hits (draw_main is shared).
void aln_meth_redraw_hits(aln_meth_read_t *res, uint64_t rng)
{
    if (res->main < 0) return;
    res->main = draw_main(res->hits, res->n_conv_moved, rng);
}
void aln_meth_copy_dup(aln_meth_read_t *dst, aln_meth_read_t *src, uint64_t rng, int keep_all, int max_xa_hits)
{
    // src is finished: its loci are stored as hits (has_copy kept them all for these draws).
    xassert(src->main < 0 || src->occ.empty(), "--meth duplicate copied before its loci were stored");
    // Copy every field except the loci: src's are set aside for the assignment, then restored.
    std::vector<aln_meth_hit_t> hits, copied;
    std::vector<uint32_t> cig, copied_cig;
    std::vector<aln_meth_blind_row_t> ties; // a copy keeps only the locus it draws (blind_pick)
    hits.swap(src->hits); cig.swap(src->cig); ties.swap(src->blind_ties);
    *dst = *src;
    src->hits.swap(hits); src->cig.swap(cig); src->blind_ties.swap(ties);
    // a copy ran no search: the work counters (aln_meth_report) stay with its first copy
    dst->n_lookups = dst->n_wrong_copy = dst->n_rescored = dst->n_flushes = dst->n_unrepresentable = 0;
    dst->n_over_diff = dst->n_over_seed = 0;
    if (src->main >= 0) dst->main = draw_main(src->hits, src->n_conv_moved, rng);
    dst->blind_main = blind_draw(src->blind_ties.size(), rng); // its own draw, as the main hit's
    if (dst->blind_main >= 0) dst->blind_pick = src->blind_ties[dst->blind_main];
    if (src->main >= 0 && !keep_all && trims(src, src->hits.size(), max_xa_hits)) {
        // Copy only the drawn primary, not the representative's whole locus list.
        aln_meth_hit_t main = src->hits[dst->main];
        if (main.n_cig > 1) {
            const uint32_t *c = aln_meth_hit_cigar(*src, main);
            copied_cig.assign(c, c + main.n_cig);
            main.cig = 0;
        }
        copied.push_back(main);
        dst->main = 0;
    } else { copied = src->hits; copied_cig = src->cig; }
    // Fresh vectors: dst's old buffers are freed.
    dst->hits.swap(copied); dst->cig.swap(copied_cig);
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
// --meth-mapq conversion-blind: locus o at row o (its lowest mirror-free row), named pos.
static aln_meth_blind_row_t blind_row(const aln_meth_occ_t &o, bwtint_t pos, int s_mm)
{
    aln_meth_blind_row_t row;
    xassert(o.alignment.re - o.alignment.rb <= ALN_MAX_READ_LEN, "--meth-mapq locus spans more than a row holds");
    row.pos = pos; row.span = (uint16_t)(o.alignment.re - o.alignment.rb);
    row.cost = (int32_t)mirror_free_cost(o.cost, s_mm, o.n_mirror);
    row.rid = o.rid; row.strand = (uint8_t)o.strand; row.hyp = (uint8_t)o.hyp; row.alt = (uint8_t)o.alt;
    row.n_mm = (uint8_t)(o.n_mm - o.n_mirror);
    row.n_gapo = (uint8_t)o.n_gapo; row.n_gape = (uint8_t)o.n_gape;
    return row;
}
// --meth-mapq conversion-blind: the strata with the mirror cell freed, over every candidate
// before the stock cut (a locus the cut drops can still tie or beat the main hit once its mirror
// mismatches are free), each locus at its lowest mirror-free row among its occ row and its mirror
// table rows (blind_loci), which also hold the loci whose every row the -n (or -k) limit dropped
// but would keep with the mirror mismatches free. That is the row collapsed scoring keeps there.
// The search proposed every window within the genomic best cost plus s_mm (aln_meth_can_stop); a
// mirror-free cost never exceeds the genomic cost, so every locus within s_mm of the mirror-free
// best was proposed too. Under stock's cuts (cuts) the strata also apply stock's top-2 ceiling on
// mirror-free differences, with each locus at the row a mirror-free run that takes the cuts keeps
// there, so they are stock's selection on mirror-free costs, as collapsed scoring computes it,
// unless the cuts may have pruned a row they need (blind_note_need; blind_incomplete). The loci
// in the strata that are not hits (the genomic finish drops them, on its cost cut, genomic_cut,
// or its difference ceiling, genomic_ceiling, or they are not in occ) are kept as rows for the
// paired-end shadow pairing, and a hit whose row is not its own as blind_alt. Out of line, as
// meth_slot_finish_tail is.
static __attribute__((noinline)) void meth_blind_strata(const aln_meth_slot_t *x, aln_meth_read_t *r,
                                                        int64_t genomic_cut, int64_t genomic_ceiling, bool cuts)
{
    const gap_opt_t *opt = x->grp->opt;
    const bool gape = opt->mode & BWA_MODE_GAPE;
    if (r->blind_partial) { std::vector<aln_meth_occ_t>().swap(r->mirror_rows); return; }
    // the cuts may have pruned a row these strata need (blind_note_need): the cap is 0
    r->blind_incomplete = x->blind_need;
    // Each locus once: its row, its occ row (null when only the mirror table holds it), and
    // whether it is a hit.
    static thread_local std::vector<blind_locus_t> bl;
    static thread_local std::vector<const aln_meth_occ_t *> rows;
    blind_loci(r, opt->s_mm, bl, rows);
    struct cand_t { const aln_meth_occ_t *row, *occ; bool hit; };
    static thread_local std::vector<cand_t> cand;
    cand.clear();
    for (const auto &l : bl) {
        const bool hit = l.occ && l.occ->cost <= genomic_cut
            && (int64_t)l.occ->n_mm + l.occ->n_gapo + (gape ? l.occ->n_gape : 0) <= genomic_ceiling;
        cand.push_back(cand_t{l.best, l.occ, hit});
    }
    auto mf_diffs = [&](const aln_meth_occ_t &o) {
        return (int64_t)o.n_mm - o.n_mirror + o.n_gapo + (gape ? o.n_gape : 0);
    };
    auto mf = [&](const aln_meth_occ_t &o) { return mirror_free_cost(o.cost, opt->s_mm, o.n_mirror); };
    int64_t best = INT64_MAX;
    for (const auto &e : cand) best = std::min(best, mf(*e.row));
    const int64_t cut = (opt->mode & BWA_MODE_NONSTOP) ? INT64_MAX : best + (int64_t)opt->s_mm;
    // stock's top-2 ceiling on the mirror-free strata: the mirror-free best rows' fewest
    // mirror-free differences + 1, as slot_finish takes it on the genomic rows
    int64_t ceiling = INT64_MAX;
    if (cuts)
        for (const auto &e : cand)
            if (mf(*e.row) == best) ceiling = std::min(ceiling, mf_diffs(*e.row) + 1);
    // A mirror-free run that takes the cuts reaches every row at its best cost before it applies
    // the ceiling, and a row above that cost only within the ceiling: a locus whose lowest
    // mirror-free row is above the best and over the ceiling holds its best row within the
    // ceiling instead (if it has none, the ceiling drops the locus below).
    if (cuts)
        for (size_t k = 0; k < cand.size(); ++k) {
            cand_t &e = cand[k];
            if (mf(*e.row) == best || mf_diffs(*e.row) <= ceiling) continue;
            const blind_locus_t &l = bl[k];
            for (uint32_t j = l.first; j < l.first + l.n; ++j)
                if (mf_diffs(*rows[j]) <= ceiling
                    && (mf_diffs(*e.row) > ceiling || mirror_free_better(*rows[j], *e.row, opt->s_mm)))
                    e.row = rows[j];
        }
    r->blind_best = best;
    r->blind_cut = cut; r->blind_ceiling = ceiling;
    // The best-cost loci, for the draw of collapsed's main hit (a read that pairs or has
    // duplicate copies). Only the drawn locus is kept, unless the read's copies draw again over
    // them all (blind_ties); the rest are gathered in a per-thread buffer.
    static thread_local std::vector<aln_meth_blind_row_t> ties_buf;
    std::vector<aln_meth_blind_row_t> &ties = x->keep == ALN_METH_KEEP_COPIES ? r->blind_ties : ties_buf;
    ties.clear();
    for (const auto &e : cand) {
        const aln_meth_occ_t &o = *e.row;
        // a hit whose row is not its own (the table can hold a copy of the occ row), keyed as the
        // hit is (aln_meth_blind_view)
        if (e.hit && e.row != e.occ
            && (mirror_free_better(*e.row, *e.occ, opt->s_mm) || mirror_free_better(*e.occ, *e.row, opt->s_mm)))
            r->blind_alt.push_back(blind_row(o, e.occ->pos, opt->s_mm));
        const int64_t c = mf(o);
        if (c > cut || mf_diffs(o) > ceiling) continue;
        if (c == best) {
            ++r->blind_c1;
            if (x->keep) ties.push_back(blind_row(o, o.alignment.rb, opt->s_mm));
        } else ++r->blind_c2;
        // store_mirror_row and account_alignment bound the cost and counts
        if (!e.hit) r->blind_rows.push_back(blind_row(o, o.alignment.rb, opt->s_mm));
    }
    std::sort(r->blind_rows.begin(), r->blind_rows.end(), blind_collapsed_order);
    r->blind_main = blind_draw(ties.size(), x->rng);
    if (r->blind_main >= 0) {
        // The loci are distinct in blind_collapsed_order (one row per locus), so the drawn index
        // names the same locus whether the ties are sorted or only partitioned around it.
        if (&ties == &r->blind_ties) std::sort(ties.begin(), ties.end(), blind_collapsed_order);
        else std::nth_element(ties.begin(), ties.begin() + r->blind_main, ties.end(), blind_collapsed_order);
        r->blind_pick = ties[r->blind_main];
    }
    // sorted for aln_meth_blind_view's binary search (a pair can hold many hits)
    std::sort(r->blind_alt.begin(), r->blind_alt.end(), [](const aln_meth_blind_row_t &a, const aln_meth_blind_row_t &b) {
        return a.pos != b.pos ? a.pos < b.pos : a.strand < b.strand;
    });
    std::vector<aln_meth_occ_t>().swap(r->mirror_rows);
}
// The --meth-mapq gate for meth_blind_strata (x->blind, as store_row's), out of line so
// slot_finish keeps a single call here.
static __attribute__((noinline)) void meth_blind_cap(const aln_meth_slot_t *x, aln_meth_read_t *r, int64_t cut,
                                                     int64_t ceiling)
{
    if (x->blind) meth_blind_strata(x, r, cut, ceiling, x->cut);
}
void aln_meth_slot_finish(aln_meth_slot_t *x)
{
    aln_meth_read_t *r = x->res;
    if (x->extra) { // --known-snvs: windows the ALT-window search found, resolved as any other
        // The lookup budget ends only the genomic lookups: these windows are already found, so
        // they are still resolved, under the CIGAR budget. ALN_METH_BIND_BEST keeps the read
        // incomplete (MAPQ 0, no pairing promotion).
        if (x->stop == ALN_METH_STOP_LOOKUP) x->stop = 0;
        for (const aln_meth_candidate_t &c : *x->extra) {
            // a lookup-budget, certificate or search-cut stop still scores them; a CIGAR cap does not
            if (x->res->flags & ALN_METH_BIND_CIGAR) break;
            if ((c.hyp ^ c.strand) != x->grp->read_no) continue;
            const int64_t span = c.re - c.rb;
            // The window search's path (its hit's CIGAR) when it is one: scored as a search row.
            aln_meth_win_key_t key{c.rid, c.hyp, c.strand, c.rb, c.re, aln_meth_layout_t()};
            const bool path = !c.cig.empty() && cigar_layout(c.cig.data(), c.cig.size(), &key.lay) &&
                layout_span(key.lay, x->len) == span;
            if (!path) key.lay.n = -1;
            if (!x->windows.insert(key, 0).second) continue;
            // As in resolve_window: a span past a window's limit is dropped as unrepresentable.
            if (span > ALN_MAX_READ_LEN) {
                ++r->n_unrepresentable;
                x->unrep_cost = 0;
                continue;
            }
            if (path) { score_row_sites(x, key.lay, c.rid, c.hyp, c.strand, c.rb, (int)span, 0, true); continue; }
            const ubyte_t *query = c.strand ? x->rd_rev.data() : x->rd_fwd.data();
            resolve_with_sites(x, c, query, std::max(0, x->len - (int)span), std::max(0, (int)span - x->len), true);
        }
    }
    if (r->occ.empty()) {
        // No record to cap, but collapsed can still place the read (its mate's shadow pairing)
        meth_blind_cap(x, r, -1, -1);
        std::vector<aln_meth_occ_t>().swap(r->mirror_rows);
        // Every candidate was dropped: any of them could have been the best locus.
        if (r->n_unrepresentable) r->flags |= ALN_METH_BIND_UNREP;
        return;
    }
    int64_t best = INT64_MAX;
    for (const auto &o : r->occ) best = std::min(best, o.cost);
    const bool nonstop = x->grp->opt->mode & BWA_MODE_NONSTOP;
    const int64_t cut = nonstop ? INT64_MAX : best + (int64_t)x->grp->opt->s_mm;
    // Stock's top-2 ceiling, for a search that takes the cuts: a retained row has at most the best
    // rows' fewest differences + 1 (the search applies the same ceiling once the best stratum is
    // complete, aln_meth_boundary). A search without the cuts keeps the cost window alone.
    const bool gape = x->grp->opt->mode & BWA_MODE_GAPE;
    auto diffs = [&](const aln_meth_occ_t &o) { return (int64_t)o.n_mm + o.n_gapo + (gape ? o.n_gape : 0); };
    int64_t ceiling = INT64_MAX;
    if (x->cut) for (const auto &o : r->occ) if (o.cost == best) ceiling = std::min(ceiling, diffs(o) + 1);
    // A dropped candidate that could be a best or retained worse locus leaves X0/X1 unknown.
    if (r->n_unrepresentable && x->unrep_cost <= cut) r->flags |= ALN_METH_BIND_UNREP;
    meth_blind_cap(x, r, cut, ceiling);
    std::vector<aln_meth_occ_t> keep;
    for (auto &o : r->occ) {
        if (o.cost > cut || diffs(o) > ceiling) continue;
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
    long lookups = 0, wrong = 0, scored = 0, flushes = 0, incomplete = 0, dropped = 0, over = 0, blind = 0;
    for (int i = 0; i < n; ++i) {
        lookups += res[i].n_lookups; wrong += res[i].n_wrong_copy; scored += res[i].n_rescored;
        dropped += res[i].n_unrepresentable; over += res[i].n_over_diff + res[i].n_over_seed;
        flushes += res[i].n_flushes; incomplete += res[i].flags != 0; blind += res[i].blind_incomplete;
    }
    fprintf(stderr, "[M::aln] --meth read %d: %d reads, %ld SA lookups (%ld wrong-copy), %ld unique windows scored (%ld dropped as unrepresentable), %ld flushes; %ld incomplete; %ld rows over -n/-k dropped; %ld with mirror-free strata the cuts may have pruned; stock cost/count MAPQ\n",
            read_no ? 1 : 2, n, lookups, wrong, scored, dropped, flushes, incomplete, over, blind);
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
	int64_t span = cigar ? bwa_cigar_ref_span(cigar, n_cigar) : len;
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
    s->mapQ = s->seQ = r->flags ? 0 : bwa_approx_mapQ(s, read_max_diff(opt, s->len));
    if (r->n_conv_moved && s->mapQ > conv_mapq_cap(r->conv_gap))
        s->mapQ = s->seQ = conv_mapq_cap(r->conv_gap);
}

// --meth-mapq conversion-blind: collapsed scoring's single-end MAPQ for end s of read r, at its
// main hit (blind_main): the stock MAPQ on the mirror-free strata's counts and that row's mismatches.
static int meth_blind_main_mapq(const aln_meth_read_t *r, const bwa_seq_t *s, const gap_opt_t *opt)
{
    bwa_seq_t t = *s;
    t.c1 = r->blind_c1; t.c2 = r->blind_c2; t.n_mm = r->blind_pick.n_mm;
    return bwa_approx_mapQ(&t, read_max_diff(opt, s->len));
}
// --meth-mapq conversion-blind: the stock single-end MAPQ the main hit would get with the mirror
// cell freed (meth_blind_strata; the main hit's locus at its mirror-free row), or 0 when a
// mirror-free locus beats it, the mirror-free ceiling drops it, or the strata are unknown
// (blind_partial, blind_incomplete); 255 (no cap) when the read has no strata (an incomplete read
// already has MAPQ 0). Read after meth_set_main, on the same max_diff.
static int meth_blind_se_cap(const aln_meth_read_t *r, const bwa_seq_t *s, const gap_opt_t *opt)
{
    if (r->main < 0 || r->flags) return 255;
    if (r->blind_partial || r->blind_incomplete) return 0;
    if (r->blind_best == INT64_MAX) return 255;
    // the main hit is not a mirror-free best locus: past the best, or over the mirror-free ceiling
    if (!aln_meth_blind_selected(*r, r->hits[r->main], opt->s_mm, opt->mode & BWA_MODE_GAPE)) return 0;
    const aln_meth_blind_row_t m = aln_meth_blind_view(*r, r->hits[r->main], opt->s_mm);
    if (m.cost > r->blind_best) return 0;
    bwa_seq_t t = *s;
    t.c1 = r->blind_c1; t.c2 = r->blind_c2; t.n_mm = m.n_mm;
    return bwa_approx_mapQ(&t, read_max_diff(opt, s->len));
}
// Lowers a record's MAPQ to its end's cap, keeping the replaced MAPQ for YQ:i.
static void meth_apply_blind_cap(bwa_seq_t *s, aln_meth_read_t *r, int cap)
{
    if (s->type == BWA_TYPE_NO_MATCH || s->mapQ <= cap) return;
    r->yq = s->mapQ;
    s->mapQ = cap;
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
    const bool blind = meth_mapq_conversion_blind(&mo->m);
    for (int i = 0; i < n; ++i) {
        meth_set_main(seqs + i, res + i, opt);
        if (blind) meth_apply_blind_cap(seqs + i, res + i, meth_blind_se_cap(res + i, seqs + i, opt));
    }
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
 * reference; *_cnt's n_mm counts only the mismatches the mode does not free. --known-snvs: with site
 * set the window carries that site allele's ALT (aln_meth_fetch_window), and *score gets the local score.
 * With fail set, a window past the allocation budget or a rescue it cannot evaluate sets *fail and
 * gives no alignment, where it otherwise stops the run. */
static bwa_cigar_t *meth_sw_core(const aln_meth_ref_t *ref, const aln_meth_opt_t *mo, int hyp, int rid,
                                 int len, const ubyte_t *raw_seq, int64_t *beg, int64_t end,
                                 int *n_cigar, uint32_t *_cnt,
                                 const gap_opt_t *account_opt = nullptr, aln_meth_occ_t *owned = nullptr,
                                 int64_t site = -1, int *score = nullptr, int *tie = nullptr,
                                 bool *fail = nullptr)
{
    *n_cigar = 0;
    if (len <= 0 || len > 65535 || rid < 0 || rid >= ref->bns->n_seqs) return nullptr;
    // aln reads a FASTQ '-' as code 5; local rescue and the realigner take codes 0-4, so the
    // rescue scores it as N. The mate's own buffers keep it.
    std::vector<uint8_t> query(raw_seq, raw_seq + len);
    for (uint8_t &c : query) c = std::min<uint8_t>(c, 4);
    const uint8_t *seq = query.data();
    const auto &ann = ref->bns->anns[rid];
    const int64_t contig_end = ann.offset + ann.len;
    // The anchor owns the window, even when its unbounded start falls in a neighbor.
    *beg = std::max<int64_t>(*beg, ann.offset);
    end = std::min(end, contig_end);
    const int64_t span = end - *beg;
    if (span < SW_MIN_MATCH_LEN) return nullptr;
    // The local-rescue profile and traceback must remain within a bounded evaluation window.
    if (span > 65535 || 2ULL * len * span + 128ULL * (len + span) > (64ULL << 20)) {
        if (fail) { *fail = true; return nullptr; }
        fprintf(stderr, "ERROR: --meth rescue window exceeds allocation budget\n"); exit(1);
    }
    int ambiguous = 0;
    for (int k = 0; k < len; ++k) ambiguous += seq[k] >= 4;
    if (ambiguous * 4 >= len || len - ambiguous < SW_MIN_MATCH_LEN) return nullptr;
    int m;
    const int8_t *mat = meth_mat(mo, hyp, &m);
    std::vector<uint8_t> target;
    fetch_window(ref, mo, *beg, span, rid, target, site);
    const int xtra = KSW_XSUBO | KSW_XSTART | (len < 250 ? KSW_XBYTE : 0);
    const kswr_t r = ksw_align(len, const_cast<uint8_t*>(seq), (int)span, target.data(), m, mat, 5, 1, xtra, nullptr);
    if (tie && r.score >= SW_MIN_MATCH_LEN && r.score2 == r.score) *tie = r.score; // tied at two loci
    if (r.score < SW_MIN_MATCH_LEN || r.score2 == r.score || r.qb < 0 || r.tb < 0 ||
        r.qe < r.qb || r.te < r.tb || r.qe >= len || r.te >= span ||
        r.qe-r.qb+1 < SW_MIN_MATCH_LEN || r.te-r.tb+1 < SW_MIN_MATCH_LEN) return nullptr;
    // Reuse the bounded, certified global refiner for the selected local rescue segment.
    const aln_meth_candidate_t candidate{rid, 0, hyp, *beg+r.tb, *beg+r.te+1};
    aln_meth_alignment_t alignment;
    if (aln_meth_realign_window(*ref, *mo, candidate, seq+r.qb, r.qe-r.qb+1, 0, 0, &alignment, site)
        != aln_meth_eval_status::ok) {
        if (fail) { *fail = true; return nullptr; }
        fprintf(stderr, "ERROR: --meth cannot evaluate mate rescue\n"); exit(1);
    }
    if (alignment.reward != r.score) return nullptr; // stock rejects inconsistent local/global scores
    if (alignment.rb < ann.offset || alignment.re > contig_end || alignment.rb >= alignment.re) {
        if (fail) { *fail = true; return nullptr; }
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
    if (score) *score = r.score;
    return legacy_cigar(cigar,n_cigar,len,len,0);
}

/* Bounded local evidence can improve placement, but cannot resolve omitted genomic loci. */
static int meth_cap_rescue(const aln_meth_ref_t *ref, const aln_meth_opt_t *mo,
                          const gap_opt_t *opt, const isize_info_t *ii,
                          bwa_seq_t *p[2], aln_meth_read_t *r[2], int end, bool short_insert,
                          bool ref_span, bool exact_window) try
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
        // The anchor's end as mate rescue takes it (aln_rescue_anchor_end): its read length, or
        // its reference span under --ref-span-insert. It places a reverse anchor's window and,
        // without --short-insert, floors a forward anchor's.
        const int64_t a_end=ref_span?a_re:(int64_t)a.pos+pa->len;
        const double begin=a.strand ? a_end-ii->avg-3*ii->std-pm->len*.5
                                   : a.pos+ii->avg-3*ii->std-pm->len*1.5;
        const double finish=begin+6*ii->std+2*pm->len;
        if(!isfinite(begin) || !isfinite(finish) || begin<INT64_MIN || finish>=INT64_MAX) {
            fprintf(stderr,"ERROR: --meth cap-rescue window out of range\n");exit(1);
        }
        int64_t beg=(int64_t)begin,last=(int64_t)finish;
        // --exact-rescue-window: the same window as mate rescue's, rounded outward from the anchor.
        if(exact_window) aln_rescue_window(a.strand?a_end:(int64_t)a.pos,!a.strand,pm->len,ii,true,&beg,&last);
        // As in mate rescue, --short-insert lets the window cover the anchor's own span; the
        // reverse anchor's ceiling is then its reference end either way.
        if(a.strand) last=std::min<int64_t>(last,short_insert?a_re:(int64_t)a.pos);
        else beg=std::max<int64_t>(beg,short_insert?(int64_t)a.pos:a_end);
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

/* Mate rescue for one pair whose ends are complete and that qualifies (an end of MAPQ at least
 * SW_MIN_MAPQ, not a proper pair; standard pair type): bwa_paired_sw_core with the anchor's matrix
 * (hyp[] is each end's hypothesis), on p[] in place. Returns the rescued end, or -1, and sets *site
 * to the site allele the rescue substituted (--known-snvs), or -1. The caller gives the rescued end
 * its anchor's hypothesis. With fail set, what would stop the run (meth_sw_core's, or a rescued
 * cost past the record) sets *fail instead, and the result is not meaningful. */
static int meth_mate_rescue(const aln_meth_ref_t *ref, const aln_meth_opt_t *mo, const gap_opt_t *gopt,
                            const pe_opt_t *popt, const isize_info_t *ii, bwa_seq_t *const p[2],
                            const int hyp[2], int64_t *site, bool *fail = nullptr)
{
    const int64_t l_pac = ref->bns->l_pac;
    int rescued = -1;
    *site = -1;
    int n_cigar[2] = {0, 0}, mapQ = 0, mq_adjust[2] = {255, 255};
    int64_t site_used[2] = {-1, -1};
    int64_t beg[2], end[2];
    bwa_cigar_t *cigar[2] = {0, 0};
    uint32_t cnt[2];
    for (int k = 0; k < 2; ++k) {
        bwa_seq_t *pref = p[1 - k], *pmate = p[k];
        if (pref->type == BWA_TYPE_NO_MATCH) continue;
        ubyte_t *seq;
        if (pref->strand == 0) {
            aln_rescue_window((int64_t)pref->pos, true, pmate->len, ii, popt->exact_rescue_window, &beg[k], &end[k]);
            if (beg[k] < aln_rescue_window_floor(pref, popt)) beg[k] = aln_rescue_window_floor(pref, popt);
            if (end[k] > l_pac) end[k] = l_pac;
            seq = pmate->rseq;
        } else {
            aln_rescue_window(aln_rescue_anchor_end(pref, popt), false, pmate->len, ii, popt->exact_rescue_window,
                              &beg[k], &end[k]);
            if (beg[k] < 0) beg[k] = 0;
            if (end[k] > aln_rescue_window_ceiling(pref, popt)) end[k] = aln_rescue_window_ceiling(pref, popt);
            seq = pmate->seq;
            seq_reverse(pmate->len, seq, 0);
        }
        const int anchor_rid = bns_pos2rid(ref->bns, pref->pos);
        const int64_t beg0 = beg[k];
        std::vector<int64_t> keys; // --known-snvs: the window's site alleles
        if (ref->known_snvs && end[k] > beg0)
            aln_meth_window_sites(*ref, beg0, (int)std::min<int64_t>(end[k] - beg0, INT_MAX), &keys,
                                  anchor_rid); // meth_sw_core clips each trial to this contig too
        if (keys.empty()) {
            cigar[k] = meth_sw_core(ref, mo, hyp[1-k], anchor_rid, pmate->len, seq,
                                    &beg[k], end[k], &n_cigar[k], &cnt[k], nullptr, nullptr, -1, nullptr,
                                    nullptr, fail);
        } else { // the window as the reference, then with each site allele alone (aln_rescue_pick)
            std::vector<pe_rescue_trial_t> trials(keys.size() + 1, pe_rescue_trial_t{nullptr, 0, 0, beg0, 0, -1, 0});
            for (size_t t = 0; t < trials.size(); ++t) {
                pe_rescue_trial_t &x = trials[t];
                x.site = t ? keys[t - 1] : -1;
                x.cigar = meth_sw_core(ref, mo, hyp[1-k], anchor_rid, pmate->len, seq,
                                       &x.beg, end[k], &x.n_cigar, &x.cnt, nullptr, nullptr, x.site, &x.score,
                                       &x.tie, fail);
            }
            const int pick = aln_rescue_pick(trials.data(), (int)trials.size(), pref, popt);
            for (int t = 0; t < (int)trials.size(); ++t) {
                if (t != pick) { free(trials[t].cigar); continue; }
                cigar[k] = trials[t].cigar, n_cigar[k] = trials[t].n_cigar, cnt[k] = trials[t].cnt;
                beg[k] = trials[t].beg, site_used[k] = trials[t].site;
            }
        }
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
        if(cost>INT_MAX && !fail) { fprintf(stderr,"ERROR: --meth rescued cost exceeds record capacity\n"); exit(1); }
        if(cost>INT_MAX) *fail = true;
        p[k]->score=(int)std::min<int64_t>(cost,INT_MAX);
        p[k]->extra_flag |= SAM_FPP;
        p[1 - k]->extra_flag |= SAM_FPP;
        *site = site_used[k], rescued = k;
    }
    free(cigar[0]); free(cigar[1]);
    return rescued;
}

// --meth-mapq conversion-blind: what meth_blind_shadow knows of collapsed scoring's pairing.
enum meth_blind_shadow_t {
    METH_BLIND_NONE = 0, // nothing: the ends' caps are already set
    METH_BLIND_PAIRED,   // the collapsed state of the ends, after its pairing
    METH_BLIND_SAME      // collapsed's pairing is the genomic one: its state is the genomic pairing's
};
// The shadow pairing's buffers, reused across pairs.
struct meth_blind_scratch_t {
    std::vector<aln_meth_blind_row_t> rows[2];
    std::vector<bwt_aln1_t> a[2];
    std::vector<int32_t> span[2];
    std::vector<pair64_t> arr;
};
// The fields of a collapsed end at row m that pairing and mate rescue read (meth_set_main's).
static void meth_blind_set_row(bwa_seq_t *s, int *hyp, const aln_meth_blind_row_t &m)
{
    s->pos = m.pos; s->strand = m.strand; *hyp = m.hyp;
    s->n_mm = m.n_mm; s->n_gapo = m.n_gapo; s->n_gape = m.n_gape; s->score = m.cost;
    s->ref_shift = m.span - (int)s->len;
}
// --meth-mapq conversion-blind, a pair whose ends are complete and of which the genomic search maps
// at least one: the state a collapsed run gives the ends before its mate rescue (cs, each end's
// hypothesis in hyp), which the conversion-blind MAPQ caps each end by once that rescue is rerun.
// An end collapsed places (blind_main) starts at its main hit with collapsed's single-end MAPQ; an
// end it does not place is unmapped (cs NO_MATCH). With both ends placed it pairs collapsed's
// candidates: the mirror-free strata (each end's hits within the mirror-free selection, the cost
// cut blind_best + s_mm and the ceiling, at their mirror-free rows, aln_meth_blind_view, and
// blind_rows), in collapsed's order, with collapsed's
// insert size ii (same_ii: whether it is the genomic pairing's); with more than --max-occ loci on
// an end collapsed does not pair, and each end stays at its main hit. Pairing then reselects each
// end's row as the genomic pairing does (the first at its start and strand). Returns
// METH_BLIND_SAME when nothing differs from the genomic run (the same insert size, no mirror
// mismatch on any hit, no hit outside the mirror-free selection, no row only the mirror-free
// strata have, no moved row, neither single-end cap below its MAPQ), where collapsed's state is
// the genomic pairing's. Returns METH_BLIND_NONE, with both caps set to 0, when an end's
// mirror-free strata are unknown (blind_partial, blind_incomplete) or a
// row's cost is past the pairing record's 20 bits. Otherwise returns METH_BLIND_PAIRED, with
// cs[j].cigar null.
static meth_blind_shadow_t meth_blind_shadow(const aln_meth_ref_t *ref, const gap_opt_t *gopt, const pe_opt_t *popt,
                                             const isize_info_t *ii, bool same_ii, bwa_seq_t *const p[2],
                                             aln_meth_read_t *const r[2], bwa_seq_t cs[2], int hyp[2],
                                             meth_blind_scratch_t &sc)
{
    if (r[0]->blind_partial || r[1]->blind_partial || r[0]->blind_incomplete || r[1]->blind_incomplete) {
        r[0]->blind_cap = r[1]->blind_cap = 0;
        return METH_BLIND_NONE;
    }
    const bool gape = gopt->mode & BWA_MODE_GAPE;
    bool differs = !same_ii;
    for (int j = 0; j < 2 && !differs; ++j) {
        differs = !r[j]->blind_rows.empty() || !r[j]->blind_alt.empty() || r[j]->blind_cap < p[j]->mapQ;
        // a hit outside the mirror-free selection (past its cut or over its ceiling)
        for (size_t x = 0; x < r[j]->hits.size() && !differs; ++x)
            differs = r[j]->hits[x].n_mirror > 0 || !aln_meth_blind_selected(*r[j], r[j]->hits[x], gopt->s_mm, gape);
    }
    if (!differs) return METH_BLIND_SAME;
    bool fits = true;
    for (int j = 0; j < 2; ++j) {
        cs[j] = *p[j];
        cs[j].cigar = nullptr; cs[j].n_cigar = 0;
        cs[j].extra_flag &= ~SAM_FPP;
        if (r[j]->blind_main < 0) { // collapsed leaves it unmapped
            cs[j].type = BWA_TYPE_NO_MATCH;
            cs[j].mapQ = cs[j].seQ = 0;
            hyp[j] = -1;
            continue;
        }
        cs[j].type = r[j]->blind_c1 > 1 ? BWA_TYPE_REPEAT : BWA_TYPE_UNIQUE;
        cs[j].mapQ = cs[j].seQ = meth_blind_main_mapq(r[j], p[j], gopt);
        meth_blind_set_row(&cs[j], &hyp[j], r[j]->blind_pick);
        fits = fits && (int64_t)r[j]->blind_c1 + r[j]->blind_c2 <= popt->max_occ;
    }
    if (r[0]->blind_main < 0 || r[1]->blind_main < 0 || !fits) return METH_BLIND_PAIRED; // no pairing
    sc.arr.clear();
    for (int j = 0; j < 2; ++j) {
        // the hits within the mirror-free selection (aln_meth_blind_selected: its cost cut and
        // ceiling), in collapsed's order as finishing sorts them, merged with blind_rows (sorted)
        std::vector<aln_meth_blind_row_t> &rows = sc.rows[j];
        rows.clear();
        for (const aln_meth_hit_t &h : r[j]->hits)
            if (aln_meth_blind_selected(*r[j], h, gopt->s_mm, gape)) rows.push_back(aln_meth_blind_view(*r[j], h, gopt->s_mm));
        if (!std::is_sorted(rows.begin(), rows.end(), blind_collapsed_order))
            std::sort(rows.begin(), rows.end(), blind_collapsed_order);
        const size_t n_hits = rows.size();
        rows.insert(rows.end(), r[j]->blind_rows.begin(), r[j]->blind_rows.end());
        std::inplace_merge(rows.begin(), rows.begin() + n_hits, rows.end(), blind_collapsed_order);
        // a cost past the pairing record's 20 bits stops a collapsed run (as below): unknown
        for (const aln_meth_blind_row_t &o : rows)
            if (o.cost > ((1 << 20) - 1)) {
                r[0]->blind_cap = r[1]->blind_cap = 0;
                return METH_BLIND_NONE;
            }
        sc.a[j].assign(rows.size(), bwt_aln1_t());
        sc.span[j].resize(rows.size());
        for (size_t x = 0; x < rows.size(); ++x) {
            const aln_meth_blind_row_t &o = rows[x];
            bwt_aln1_t &b = sc.a[j][x];
            b.n_mm = o.n_mm; b.n_gapo = o.n_gapo; b.n_gape = o.n_gape; b.score = o.cost;
            sc.span[j][x] = o.span;
            pair64_t q;
            q.x = o.pos;
            q.y = (uint64_t)x << 2 | (uint64_t)o.strand << 1 | (uint64_t)j;
            sc.arr.push_back(q);
        }
        cs[j].aln = sc.a[j].data();
    }
    const int32_t *spans[2] = {sc.span[0].data(), sc.span[1].data()};
    bwa_seq_t *sp[2] = {&cs[0], &cs[1]};
    pair64_v v;
    v.a = sc.arr.data(); v.n = v.m = sc.arr.size();
    bwa_pe_pairing_same_contig(ref->bns, sp, &v, popt, std::min(gopt->s_mm, INT_MAX/10), ii,
                               popt->ref_span_insert ? spans : nullptr);
    // each end's row as aln_meth_pe_stage2 reselects it after pairing: the first at its start and
    // strand (with --ref-span-insert, also of its span and cost)
    for (int j = 0; j < 2; ++j) {
        cs[j].aln = p[j]->aln;
        for (const aln_meth_blind_row_t &o : sc.rows[j])
            if (o.pos == cs[j].pos && o.strand == (int)cs[j].strand
                && (!popt->ref_span_insert || (o.span - (int)cs[j].len == cs[j].ref_shift && o.cost == cs[j].score))) {
                meth_blind_set_row(&cs[j], &hyp[j], o);
                break;
            }
    }
    return METH_BLIND_PAIRED;
}

// --meth-mapq conversion-blind: collapsed's mate rescue of a pair, rerun with collapsed's matrices
// (cmo) and insert size (ii) from cs, its state before that rescue (with same, the genomic pairing's
// state, p, which this copies), as aln_meth_pe_stage2 runs it. Each end's blind_place_* get where
// collapsed then leaves it, and its MAPQ there (0 where it leaves the end unmapped). What would stop
// a collapsed run's rescue (a window past meth_sw_core's allocation bound, a rescue it cannot
// evaluate, a rescued cost past the record) leaves both at MAPQ 0 instead.
static void meth_blind_place(const aln_meth_ref_t *ref, const aln_meth_opt_t *cmo, const gap_opt_t *gopt,
                             const pe_opt_t *popt, const isize_info_t *ii, bwa_seq_t *const p[2],
                             aln_meth_read_t *const r[2], bwa_seq_t cs[2], int hyp[2], bool same)
{
    if (same)
        for (int j = 0; j < 2; ++j) {
            cs[j] = *p[j];
            cs[j].cigar = nullptr; cs[j].n_cigar = 0;
            hyp[j] = r[j]->hyp;
        }
    bool unknown = false;
    if (popt->is_sw && ii->avg >= 0.0 && popt->type == BWA_PET_STD
        && (cs[0].mapQ >= SW_MIN_MAPQ || cs[1].mapQ >= SW_MIN_MAPQ) && (cs[0].extra_flag & SAM_FPP) == 0) {
        bwa_seq_t *cp[2] = {&cs[0], &cs[1]};
        int64_t site;
        const int k = meth_mate_rescue(ref, cmo, gopt, popt, ii, cp, hyp, &site, &unknown);
        if (k >= 0) { free(cs[k].cigar); cs[k].cigar = nullptr; }
    }
    for (int j = 0; j < 2; ++j) {
        r[j]->blind_place_pos = cs[j].pos; r[j]->blind_place_strand = cs[j].strand;
        r[j]->blind_place_mapq = unknown || cs[j].type == BWA_TYPE_NO_MATCH ? 0 : cs[j].mapQ;
    }
}

void aln_meth_blind_isize(const aln_meth_ref_t *ref, const gap_opt_t *gopt, const pe_opt_t *popt, int n,
                          bwa_seq_t *seqs[2], aln_meth_read_t *const res[2], const uint8_t *excluded,
                          const isize_info_t *last_ii, isize_info_t *ii)
{
    // Each end's fields infer_isize reads, set to collapsed's single-end selection and then restored.
    struct saved_t { bwtint_t pos; int strand, mapq, ref_shift; };
    std::vector<saved_t> saved[2];
    for (int j = 0; j < 2; ++j) {
        saved[j].resize(n);
        for (int i = 0; i < n; ++i) {
            bwa_seq_t *s = seqs[j] + i;
            const aln_meth_read_t &r = res[j][i];
            saved[j][i] = saved_t{s->pos, (int)s->strand, (int)s->mapQ, s->ref_shift};
            const bool placed = s->type == BWA_TYPE_UNIQUE || s->type == BWA_TYPE_REPEAT
                || (s->type == BWA_TYPE_NO_MATCH && r.blind_main >= 0); // collapsed places it
            if (!placed || r.flags) continue;
            // collapsed's selection is unknown: not an observation
            if (r.blind_partial || r.blind_incomplete) s->mapQ = 0;
            else if (r.blind_main >= 0) {
                const aln_meth_blind_row_t &m = r.blind_pick;
                s->mapQ = meth_blind_main_mapq(&r, s, gopt);
                s->pos = m.pos; s->strand = m.strand; s->ref_shift = m.span - (int)s->len;
            }
        }
    }
    bwa_pe_isize_ref(ref->bns->l_pac, n, seqs, ii, popt, last_ii, ref->bns, 0, excluded);
    for (int j = 0; j < 2; ++j)
        for (int i = 0; i < n; ++i) {
            bwa_seq_t *s = seqs[j] + i;
            s->pos = saved[j][i].pos; s->strand = saved[j][i].strand;
            s->mapQ = saved[j][i].mapq; s->ref_shift = saved[j][i].ref_shift;
        }
}

int aln_meth_pe_stage2(const aln_meth_ref_t *ref, const aln_meth_opt_t *mo, const gap_opt_t *gopt, const pe_opt_t *popt,
					   const isize_info_t *ii, int n, bwa_seq_t *seqs[2], aln_meth_read_t *res[2],
					   uint64_t n_tot[2], uint64_t n_mapped[2], const isize_info_t *blind_ii)
{
	int cnt_chg = 0;
	// Complete means no lookup/work cutoff in configured discovery, not all genome hits.
	const int max_occ = popt->max_occ;
	std::vector<bwt_aln1_t> a[2];
	std::vector<int32_t> span[2]; // --ref-span-insert: each hit's reference span, for pairing
	std::vector<pair64_t> arr;
	meth_blind_scratch_t shadow; // --meth-mapq conversion-blind: the shadow pairing's buffers
	long n_blind_unknown = 0; // pairs with an end whose mirror-free strata are unknown
	kstring_t str = {0, 0, 0};
    std::vector<std::pair<int,int>> initial_selection(n);
    for(int i=0;i<n;++i) initial_selection[i]={res[0][i].main,res[1][i].main};
    // --meth-mapq conversion-blind: each end starts capped at its single-end conversion-blind
    // MAPQ (stage 1's selection); for a complete pair with a mapped end, where collapsed's pairing
    // and mate rescue leave each end (meth_blind_shadow, meth_blind_place) replaces that.
    const bool blind = meth_mapq_conversion_blind(&mo->m);
    // whether collapsed's insert size is the genomic pairing's (pairing reads every field)
    // (pairing reads low and high, and with high set avg, std and high_bayesian)
    const bool same_ii = !blind_ii || (blind_ii->low == ii->low && blind_ii->high == ii->high
        && (ii->high == 0 || (blind_ii->avg == ii->avg && blind_ii->std == ii->std
                              && blind_ii->high_bayesian == ii->high_bayesian)));
    // An end whose mirror-free strata are unknown (blind_partial, blind_incomplete) leaves the
    // pair's unknown, whether or not the pairing below runs: both ends get 0, without the shadow.
    if (blind)
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < 2; ++j) res[j][i].blind_cap = meth_blind_se_cap(res[j] + i, seqs[j] + i, gopt);
            const bool unknown = res[0][i].blind_partial || res[1][i].blind_partial || res[0][i].blind_incomplete
                || res[1][i].blind_incomplete;
            if (unknown) res[0][i].blind_cap = res[1][i].blind_cap = 0;
            n_blind_unknown += res[0][i].blind_incomplete || res[1][i].blind_incomplete;
        }
    // collapsed scoring's options, for its mate rescue (meth_blind_place): its matrices
    aln_meth_opt_t cmo;
    if (blind) {
        cmo = *mo;
        cmo.m.meth_scoring = MEM_METH_SCORING_COLLAPSED;
        cmo.m.meth_mapq = MEM_METH_MAPQ_SCORING;
        aln_meth_opt_finish(&cmo);
    }
	for (int i = 0; i < n; ++i) {
		bwa_seq_t *p[2] = {seqs[0] + i, seqs[1] + i};
		aln_meth_read_t *r[2] = {res[0] + i, res[1] + i};
		// Explicitly incomplete ends cannot be promoted by pairing or rescue.
		const bool incomplete = r[0]->flags || r[1]->flags;
		const bool mapped = (p[0]->type == BWA_TYPE_UNIQUE || p[0]->type == BWA_TYPE_REPEAT)
			&& (p[1]->type == BWA_TYPE_UNIQUE || p[1]->type == BWA_TYPE_REPEAT);
		// --meth-mapq conversion-blind: collapsed scoring's state of the ends before its mate rescue
		// (meth_blind_shadow), from where they start, for a complete pair the genomic search maps an
		// end of. An end the genomic search left unmapped is placed by collapsed when it has mirror-free
		// strata (rows -n or -k dropped only for their mirror mismatches).
		bwa_seq_t cs[2];
		int cs_hyp[2] = {-1, -1};
		const meth_blind_shadow_t shadow_st = blind && !incomplete
			&& (p[0]->type != BWA_TYPE_NO_MATCH || p[1]->type != BWA_TYPE_NO_MATCH)
			? meth_blind_shadow(ref, gopt, popt, blind_ii ? blind_ii : ii, same_ii, p, r, cs, cs_hyp, shadow)
			: METH_BLIND_NONE;
		// pairing on the retained stock-cost occurrences
		if (mapped) {
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
				const int32_t *spans[2] = {span[0].data(), span[1].data()};
				v.a = arr.data(); v.n = v.m = arr.size();
				bwt_aln1_t *save[2] = {p[0]->aln, p[1]->aln};
				p[0]->aln = a[0].data(); p[1]->aln = a[1].data();
                // Stock multiplies s_mm by ten; clamp only above all representable pair-score gaps.
                cnt_chg += bwa_pe_pairing_same_contig(ref->bns, p, &v, popt,
                                                     std::min(gopt->s_mm, INT_MAX/10), ii,
                                                     popt->ref_span_insert ? spans : nullptr);
				p[0]->aln = save[0]; p[1]->aln = save[1];
				// the ends' new loci (pairing may move either): their hypotheses and gaps. With
				// --ref-span-insert pairing installs the selected hit's cost and span even when the end
				// keeps its start and strand, so match those too rather than the first hit there.
				for (int j = 0; j < 2; ++j) {
					for (size_t x = 0; x < r[j]->hits.size(); ++x) {
						const aln_meth_hit_t &o = r[j]->hits[x];
						if (o.pos == p[j]->pos && o.strand == (int)p[j]->strand
						    && (!popt->ref_span_insert
						        || ((int)o.span - (int)p[j]->len == p[j]->ref_shift && (int)o.cost == p[j]->score))) {
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
		// collapsed's mate rescue, rerun from that state with collapsed's matrices and insert size,
		// and where it leaves each end (meth_blind_place): an end is capped at collapsed's MAPQ where
		// collapsed places it as the genomic run finally does, else 0
		if (shadow_st != METH_BLIND_NONE)
			meth_blind_place(ref, &cmo, gopt, popt, blind_ii ? blind_ii : ii, p, r, cs, cs_hyp,
			                 shadow_st == METH_BLIND_SAME);
	}
    if(popt->is_sw && ii->avg>=0 && popt->type==BWA_PET_STD)
        for(int i=0;i<n;++i) {
            bwa_seq_t *p[2]={seqs[0]+i,seqs[1]+i};
            aln_meth_read_t *r[2]={res[0]+i,res[1]+i};
            const bool cap0=r[0]->flags&ALN_METH_BIND_BEST,cap1=r[1]->flags&ALN_METH_BIND_BEST;
            if(cap0!=cap1) cnt_chg+=meth_cap_rescue(ref,mo,gopt,ii,p,r,cap0?0:1,popt->short_insert,
                                                popt->ref_span_insert,popt->exact_rescue_window);
        }
	// mate rescue (bwa_paired_sw_core with the anchor's matrix; the rescued mate takes its hypothesis)
	if (popt->is_sw && ii->avg >= 0.0) {
		for (int i = 0; i < n; ++i) {
			bwa_seq_t *p[2] = {seqs[0] + i, seqs[1] + i};
			aln_meth_read_t *r[2] = {res[0] + i, res[1] + i};
			if (r[0]->flags || r[1]->flags) continue;
			if (!((p[0]->mapQ >= SW_MIN_MAPQ || p[1]->mapQ >= SW_MIN_MAPQ) && (p[0]->extra_flag & SAM_FPP) == 0)) continue;
			const int is_singleton = (p[0]->type == BWA_TYPE_NO_MATCH || p[1]->type == BWA_TYPE_NO_MATCH) ? 1 : 0;
			++n_tot[is_singleton];
			if (popt->type != BWA_PET_STD) continue;
			const int hyp[2] = {r[0]->hyp, r[1]->hyp};
			int64_t site;
			const int k = meth_mate_rescue(ref, mo, gopt, popt, ii, p, hyp, &site);
			if (k >= 0) {
				++n_mapped[is_singleton];
				r[k]->hyp = r[1 - k]->hyp; // the rescued mate is on the anchor's genome strand
				r[k]->ks_site = site;
			}
		}
	}
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < 2; ++j) {
            bwa_seq_t *s = seqs[j] + i;
            aln_meth_read_t *r = res[j] + i;
            // Rescue must not inflate an unchanged conversion-selected anchor either. A rescued
            // (MATESW) end has left its conversion-selected locus (r->main is stale for it): its
            // MAPQ comes from the pair and is already bounded by the capped anchor's.
            if (r->n_conv_moved && s->type != BWA_TYPE_MATESW &&
                r->main == (j ? initial_selection[i].second : initial_selection[i].first) &&
                s->mapQ > conv_mapq_cap(r->conv_gap)) s->mapQ = conv_mapq_cap(r->conv_gap);
            // Where collapsed's state is known (meth_blind_place), an end is capped at collapsed's
            // MAPQ where collapsed leaves it at this end's start and strand, else 0. Otherwise
            // (no shadow: an incomplete end, or unknown strata or costs, which set the caps) a
            // rescued (MATESW) end, which has no strata of its own, takes its anchor's cap.
            if (blind) {
                int cap = r->blind_cap;
                if (r->blind_place_mapq >= 0)
                    cap = r->blind_place_pos == s->pos && r->blind_place_strand == (int)s->strand ? r->blind_place_mapq : 0;
                else if (s->type == BWA_TYPE_MATESW) cap = res[1 - j][i].blind_cap;
                meth_apply_blind_cap(s, r, cap);
            }
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
	if (blind && bwa_verbose >= 4)
		fprintf(stderr, "[M::aln] --meth-mapq conversion-blind: %ld of %d pairs have an end whose mirror-free strata the cuts may have pruned (both ends capped to 0)\n",
				n_blind_unknown, n);
	free(str.s);
	return cnt_chg;
}

/* ---------------------------------------------------------------------------------------------
 * Output
 * ------------------------------------------------------------------------------------------- */

// SEQ text for a --meth record. aln reads a FASTQ '-' as code 5 (upstream's table, kept for
// stock parity); stock prints it as a NUL byte. --meth records have no stock counterpart, so
// any code above 3 is printed as N.
static void meth_seq_text(std::string *out, const bwa_seq_t *p)
{
	out->clear();
	out->reserve(p->full_len);
	if (p->strand == 0) for (int i = 0; i < (int)p->full_len; ++i) out->push_back("ACGTN"[std::min<ubyte_t>(p->seq[i], 4)]);
	else for (int i = (int)p->full_len - 1; i >= 0; --i) out->push_back("TGCAN"[std::min<ubyte_t>(p->seq[i], 4)]);
}
static void meth_print_seq(kstring_t *out, const bwa_seq_t *p)
{
	std::string text;
	meth_seq_text(&text, p);
	kputsn(text.data(), (int)text.size(), out);
}

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
		meth_print_seq(out, p);
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
			if (r->yq >= 0) { kputs("\tYQ:i:", out); kputl(r->yq, out); }
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
				meth_seq_text(&seq_text, p);
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
		meth_print_seq(out, p);
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
