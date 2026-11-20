/* SPDX-License-Identifier: MIT */
/* Internal candidate generation hooks; no search ancestry is retained. */
#ifndef BWA_MEM3_ALN_METH_SEARCH_H
#define BWA_MEM3_ALN_METH_SEARCH_H
#include <limits.h>
#include <algorithm>
#include <map>
#include <utility>
#include <set>
#include <tuple>
#include "aln_meth.h"
#include "aln_meth_realign.h"
#include "aln_meth_prove.h"
#include <vector>
#include "bwtgap.h"
namespace bwa_aln {
// SA rows a read resolves at the flush that binds its lookup budget (at most what is left of it),
// to pick its best-known placement: that flush makes the read incomplete either way.
enum { ALN_METH_BIND_ROWS = 64 };
// A hit waiting for the next flush: its interval, entry and the path's gaps (aln_meth_hit unpacks them).
struct aln_meth_pend_t { bwtint_t k, l; gap_entry_t e; uint32_t gap[ALN_METH_MAX_GAPS]; };
/* A flat, open-addressed index of distinct keys (insertion order kept in key, a value per key),
 * cleared per read. Only membership and lookup are used, as the std::set / std::map it replaces
 * were, so the results do not depend on the container. K needs operator== and a hash(). */
template <class K> struct aln_flat_index_t {
    std::vector<K> key;
    std::vector<size_t> val;
    std::vector<int32_t> tab;
    size_t size() const { return key.size(); }
    void clear() {
        key.clear(); val.clear();
        if (tab.size() != 64) tab.assign(64, -1);
        else std::fill(tab.begin(), tab.end(), -1);
    }
    // the slot of k: its index in key, or -1 (and *pos is where to insert it)
    int32_t find(const K &k, size_t *pos) const {
        const size_t mask = tab.size() - 1;
        for (size_t i = k.hash() & mask;; i = (i + 1) & mask) {
            const int32_t t = tab[i];
            if (t < 0) { *pos = i; return -1; }
            if (key[t] == k) { *pos = i; return t; }
        }
    }
    void grow() {
        tab.assign(tab.size() * 2, -1);
        const size_t mask = tab.size() - 1;
        for (size_t j = 0; j < key.size(); ++j) {
            size_t i = key[j].hash() & mask;
            while (tab[i] >= 0) i = (i + 1) & mask;
            tab[i] = (int32_t)j;
        }
    }
    // inserts k with value v unless present; returns its index and whether it was inserted
    std::pair<size_t, bool> insert(const K &k, size_t v) {
        if (tab.empty()) tab.assign(64, -1);
        size_t pos;
        const int32_t t = find(k, &pos);
        if (t >= 0) return {(size_t)t, false};
        tab[pos] = (int32_t)key.size();
        key.push_back(k); val.push_back(v);
        if (key.size() * 2 > tab.size()) grow();
        return {key.size() - 1, true};
    }
};
static inline uint64_t aln_flat_mix(uint64_t h) {
    h ^= h >> 33; h *= 0xff51afd7ed558ccdULL; h ^= h >> 33; h *= 0xc4ceb9fe1a85ec53ULL; return h ^ (h >> 33);
}
// A hit interval and the gap events of the path that reached it: two paths with different gap
// placements can reach the same interval (a gap in a tandem repeat) and score differently on the
// original reference, so each is resolved.
struct aln_meth_iv_key_t {
    bwtint_t k, l; uint32_t n; uint32_t gap[ALN_METH_MAX_GAPS];
    bool operator==(const aln_meth_iv_key_t &o) const {
        if (k != o.k || l != o.l || n != o.n) return false;
        for (uint32_t g = 0; g < n; ++g) if (gap[g] != o.gap[g]) return false;
        return true;
    }
    uint64_t hash() const {
        uint64_t h = (uint64_t)k * 0x9e3779b97f4a7c15ULL ^ (uint64_t)l ^ (uint64_t)n << 56;
        for (uint32_t g = 0; g < n; ++g) h = h * 0x100000001b3ULL ^ gap[g];
        return aln_flat_mix(h);
    }
};
struct aln_meth_kl_key_t {
    bwtint_t k, l;
    bool operator==(const aln_meth_kl_key_t &o) const { return k == o.k && l == o.l; }
    uint64_t hash() const { return aln_flat_mix((uint64_t)k * 0x9e3779b97f4a7c15ULL ^ (uint64_t)l); }
};
/* A path's alignment in the orientation the rows are scored in (the query as resolve passes it,
 * against the forward reference window): n gaps, ascending by qpos. An insertion covers query
 * [qpos, qpos + len); a deletion of len reference bases follows query base qpos - 1. */
struct aln_meth_layout_t {
    int n = 0;
    int qpos[ALN_METH_MAX_GAPS], len[ALN_METH_MAX_GAPS], del[ALN_METH_MAX_GAPS];
    bool operator==(const aln_meth_layout_t &o) const {
        if (n != o.n) return false;
        for (int g = 0; g < n; ++g)
            if (qpos[g] != o.qpos[g] || len[g] != o.len[g] || del[g] != o.del[g]) return false;
        return true;
    }
    uint64_t hash() const {
        uint64_t h = (uint64_t)n;
        for (int g = 0; g < n; ++g) h = h * 0x100000001b3ULL ^ ((uint64_t)qpos[g] << 32 | (uint64_t)len[g] << 1 | (uint64_t)del[g]);
        return h;
    }
};
struct aln_meth_win_key_t {
    int rid, hyp, strand; int64_t rb, re;
    aln_meth_layout_t lay; // n = -1: a window without a path (refined by DP)
    bool operator==(const aln_meth_win_key_t &o) const {
        return rid == o.rid && hyp == o.hyp && strand == o.strand && rb == o.rb && re == o.re && lay == o.lay;
    }
    uint64_t hash() const {
        return aln_flat_mix((uint64_t)rb * 0x9e3779b97f4a7c15ULL ^ (uint64_t)re ^ ((uint64_t)(hyp * 2 + strand) << 60) ^ (uint64_t)rid << 40 ^ lay.hash());
    }
};
struct aln_meth_locus_key_t {
    int rid, hyp, strand; int64_t rb;
    bool operator==(const aln_meth_locus_key_t &o) const {
        return rid == o.rid && hyp == o.hyp && strand == o.strand && rb == o.rb;
    }
    uint64_t hash() const { return aln_flat_mix((uint64_t)rb ^ ((uint64_t)(hyp * 2 + strand) << 60) ^ (uint64_t)rid << 40); }
};
struct aln_meth_group_t;
enum { ALN_METH_KEEP_STRATA = 1, ALN_METH_KEEP_COPIES = 2 }; // aln_meth_slot_t::keep
// aln_meth_slot_t::stop: why a read stopped resolving windows (0 while it has not); TOP2 is
// stock's -R stop under the cut (aln_meth_boundary), CERT a saturation certificate (a read with
// aln_meth_slot_t::cert; one that gives up incomplete instead stops with LOOKUP)
enum { ALN_METH_STOP_CIGAR = 1, ALN_METH_STOP_LOOKUP = 2, ALN_METH_STOP_TOP2 = 3, ALN_METH_STOP_CERT = 4 };
struct aln_meth_slot_t {
    aln_meth_group_t *grp = nullptr;
    gap_job_t *job = nullptr;
    bwa_seq_t *p = nullptr;
    aln_meth_read_t *res = nullptr;
    // keep (ALN_METH_KEEP_*): 0 trims a finished read to its main hit (plus XA loci); STRATA keeps
    // a complete read's loci for pairing; COPIES keeps them for duplicate copies to redraw over.
    // max_diff: the read's -n limit (bwa_cal_maxdiff on its length for a fractional -n), set by
    // aln_meth_slot_start; a refined window with more differences is dropped (resolve_one).
    int len = 0, keep = 0, final_flush = 0, stop = 0, flushed_score = -1, max_diff = INT_MAX;
    // cut: the search takes stock aln's cuts from the verified best (aln_meth_boundary): set unless
    // -N or a negative penalty; under --known-snvs set in the ALT-window search, and in the main
    // search only for a read that fits an ALT window (ks_fits_window). best: the lowest cost of a
    // stored row. flushed: a flush has run. gen: the --known-snvs ALT-window search, a candidate
    // generator: it takes only the cuts that later rows of the main search cannot invalidate
    // (aln_meth_boundary_gen) and keeps every row it found.
    int cut = 0, flushed = 0, gen = 0;
    // cert: rows are verified stratum by stratum and certificates end the read once its best
    // stratum is complete (single-end, under cut, not under --known-snvs). ivs: every hit flushed,
    // for the draw.
    int cert = 0;
    // extras_done: --known-snvs, the ALT-window candidates (extra) are scored (aln_meth_score_extras).
    // ceil_applied: the tightest top-2 difference ceiling aln_meth_boundary applied (INT_MAX: none);
    // r_stop: it ended the search on -R. Both for the finish-time audit.
    int extras_done = 0, ceil_applied = INT_MAX, r_stop = 0;
    // prove: the search tries the pigeonhole proof once, at its first entry past the exact stratum
    // (aln_meth_prove_nohit): set under BWA_MODE_GAPE outside the --known-snvs ALT-window search, for
    // reads allowed at least ALN_METH_PH_MIN_DIFF differences.
    int prove = 0;
    // prove2_at: once the search has made more Occ4 expansions (n_steps) than this with no hit yet,
    // the proof is tried again with larger caps (aln_meth_prove_nohit2); INT64_MAX when it is not.
    int64_t n_steps = 0, prove2_at = INT64_MAX;
    aln_meth_ph_scratch_t ph;
    std::vector<aln_meth_pend_t> ivs;
    uint64_t perm_rng = 0;
    int64_t best = INT64_MAX;
    // cap_best: the --known-snvs ALT-window generator's seed, an upper bound on the read's final best
    // known before its search (INT64_MAX: none). The generator's best starts there, so its score
    // cut, K(B), -R and horizon (aln_meth_boundary_gen, meth_gen_finish) use min(cap_best, its own
    // best) from the first entry, before it has a row. Unset (INT64_MAX) on any other search.
    int64_t cap_best = INT64_MAX;
    uint64_t rng = 0, alignment_bytes = 0;
    int64_t unrep_cost = INT64_MAX; // the lowest cost of a refined window dropped as unrepresentable
    // Gap lists of the search's paths with more than one gap, ALN_METH_MAX_GAPS words each (see
    // ALN_METH_GAP; an entry's last_diff_pos indexes its path's list). Cleared per read.
    std::vector<uint32_t> garena;
    // The SA rows of each interval resolved in full (sa_ix: its offset in sa_rows), for an interval
    // that a path with other gaps reaches again. Cleared per read.
    aln_flat_index_t<aln_meth_kl_key_t> sa_ix;
    std::vector<bwtint_t> sa_rows;
    // Refinement queries for each strand; any code above 4 (aln reads '-' as 5) is stored as 4.
    std::vector<ubyte_t> vseq, rd_fwd, rd_rev;
    std::vector<aln_meth_pend_t> pend;
    aln_flat_index_t<aln_meth_iv_key_t> intervals;
    // Exact windows dedup before DP; every distinct span is evaluated before locus collapse.
    aln_flat_index_t<aln_meth_win_key_t> windows;
    aln_flat_index_t<aln_meth_locus_key_t> loci; // value: the locus's index in res->occ
    // The read's two query orientations as bit planes (rd_fwd, then p->rseq), built on first use.
    aln_meth_read_bits_t bits[2];
    int bits_ready = 0;
    // --known-snvs: candidate windows from the ALT-window search, resolved before the read finishes
    const aln_meth_extra_t *extra = nullptr;
    // the ordinals of the read's duplicate copies (their draws are made here, aln_meth_slot_finish)
    const std::vector<int64_t> *copy_ords = nullptr;
    // Scratch: the reference window of the candidate being resolved (refined, then accounted), and
    // the ambiguity bits of a row's window (score_path).
    std::vector<uint8_t> win;
    std::vector<uint64_t> amb;
};
struct aln_meth_group_t {
    const bwt_t *bwt = nullptr;
    const aln_meth_ref_t *ref = nullptr;
    const aln_meth_opt_t *mo = nullptr;
    const gap_opt_t *opt = nullptr;
    int read_no = 1, keep_strata = 0, max_xa_hits = 3;
    // probe: the --known-snvs probe of the main search (aln_meth_probe_batch): its exact stratum only,
    // at most ALN_METH_PROBE_ROWS rows, giving each read an upper bound on its final best.
    int probe = 0;
    int n_slots = 0, n_active = 0, n_parked = 0, k_park = 1;
    std::vector<aln_meth_slot_t> slot;
    gap_slot_core_t *cores = nullptr;
};
// Applies stock aln's cuts at a score boundary, for a search whose next entry has score next:
// see aln_meth.cpp. Returns 1 when the read is done (stock's -R stop, or for a certificate read a
// certificate, aln_meth_certify).
int aln_meth_boundary(aln_meth_slot_t *, int next);
// Keep batches bounded. Under the cut, flush at the first hit and at every boundary; without it
// (-N, a negative penalty, or a --known-snvs read that does not fit an ALT window: the search runs
// to its -n) flush only for the bound and to find the lookup budget's binding flush at its stratum.
static inline int aln_meth_needs_flush(aln_meth_slot_t *x, int score) {
    if (x->pend.size() >= 32) return 1;
    if (x->cut) { // flush at the first hit and at every score boundary; cut at boundaries
        if (!x->pend.empty()) return !x->flushed || score > x->flushed_score;
        if (score <= x->flushed_score) return 0;
        x->flushed_score = score;
        // a certificate read's boundary sets its own stop (aln_meth_certify)
        return x->best != INT64_MAX && aln_meth_boundary(x, score) && (x->stop || (x->stop = ALN_METH_STOP_TOP2));
    }
    if (score <= x->flushed_score) return 0;
    // A skipped boundary replaces a flush/resume: do not re-examine this stratum.
    x->flushed_score=score;
    if (x->pend.empty()) return 0;
    // Keep lookup-cap detection ahead of later work-cap failures, as with boundary flushing.
    int64_t remaining=(int64_t)x->grp->mo->p_best+x->grp->mo->p_second;
    if (remaining<0) remaining=0;
    remaining-=x->res->n_lookups;
    for (const auto &h : x->pend) {
        const bwtint_t width=h.l-h.k+1;
        if (remaining<0 || width>(uint64_t)remaining) return 1;
        remaining-=(int64_t)width;
    }
    return 0;
}
static inline int aln_meth_needs_flush(const void *, int) { return 0; }
static inline int aln_meth_has_pending(const aln_meth_slot_t *x) { return !x->pend.empty(); }
static inline int aln_meth_has_pending(const void *) { return 0; }
static inline void aln_meth_set_final(aln_meth_slot_t *x) { x->final_flush = 1; }
static inline void aln_meth_set_final(void *) {}
static inline void aln_meth_flag(aln_meth_slot_t *x, int f) { x->res->flags |= f; }
static inline void aln_meth_flag(void *, int) {}
static inline int aln_meth_resume(aln_meth_slot_t *x, gap_job_t *j, gap_entry_t *e) {
    if (x->final_flush || x->stop) return 1;
    *e = j->e;
    const int score=(int)(e->info>>21);
    x->flushed_score=score;
    return x->cut && x->best != INT64_MAX && aln_meth_boundary(x, score);
}
static inline int aln_meth_resume(void *, gap_job_t *, gap_entry_t *) { return 1; }
void aln_meth_hit(aln_meth_slot_t *, const gap_entry_t *, bwtint_t, bwtint_t);
static inline void aln_meth_hit(void *, const gap_entry_t *, bwtint_t, bwtint_t) {}
// A child of entry p (n_gapo >= 1) whose gaps need an arena list: p's gaps with the last one a base
// longer (ext), or plus the new gap w. Returns the list's index.
int aln_meth_gap_arena(aln_meth_slot_t *, const gap_entry_t &p, bool ext, uint32_t w);
// has_copy: the read has duplicate copies that will redraw over its loci.
// cap_best: the ALT-window generator's seed (aln_meth_slot_t::cap_best); ignored by other searches.
const ubyte_t *aln_meth_slot_start(aln_meth_slot_t *, bwa_seq_t *, aln_meth_read_t *, int64_t, int, int has_copy,
                                   int64_t cap_best = INT64_MAX);
void aln_meth_slot_finish(aln_meth_slot_t *);
/* Whether row a is a better row for its locus than b: lower stock cost; then fewer stock
 * differences (mismatches and gap opens, and gap extensions under gape), so that the locus keeps,
 * among its rows at the best cost, the fewest differences any of them has, whichever row was found
 * first, and the top-2 ceiling taken from the stored rows (best_ceiling) is the same in any order;
 * then, under the TAPS/5base tie-break (conversions), fewer conversions; then (--known-snvs) the
 * reference before a site allele, and the lower aln_site_key; then the normalized and raw CIGARs and
 * the end, for a deterministic order. */
bool aln_meth_row_better(const aln_meth_occ_t &a, const aln_meth_occ_t &b, bool conversions, bool gape);
/* A search path's gap differences inside the seed (above seed_off), counted from its gap events
 * (ALN_METH_GAP), and the same count from the path as a layout (meth_layout), as an ALT-window
 * candidate's CIGAR gives it. The two must agree: -k treats both kinds of row alike. */
int aln_meth_seed_gap_diffs(const uint32_t *gap, int n_gapo, int seed_off, bool gape);
int aln_meth_layout_seed_diffs(const uint32_t *gap, int n_gapo, int len, int strand, int seed_off, bool gape);
/* --known-snvs: scores the read's ALT-window candidates (x->extra) as rows of the main search. The
 * driver calls it once the job is started (its seed is then known) and before the search's first
 * step, so the read's best and its cuts start from them; aln_meth_slot_finish calls it for a slot
 * that has not (no-op once done). */
void aln_meth_score_extras(aln_meth_slot_t *);
/* The pigeonhole proof for the read of slot x, whose search j has found no hit in its exact stratum
 * (no interval yet): returns 1 when no search path within the read's -n differences reaches a hit
 * (aln_meth_pigeonhole), so ending the search there leaves every result as the full search would
 * (no hit, row, lookup or draw). A read whose root bound leaves no spare difference is not tried:
 * its search is cheap. */
enum { ALN_METH_PH_MIN_DIFF = 3 };
int aln_meth_prove_nohit(aln_meth_slot_t *x, const gap_job_t *j);
int aln_meth_prove_nohit2(aln_meth_slot_t *x, const gap_job_t *j);
void aln_meth_report(int, int, const aln_meth_read_t *);
// The rows (SA lookups) a --known-snvs probe resolves per read (aln_meth_probe_batch): a few verified
// rows suffice for an upper bound, and a wide interval would cost its whole width.
enum { ALN_METH_PROBE_ROWS = 32 };
// Finishes a read (stores its occ as hits; see aln_meth_store_loci), then trims it to its main hit
// (an incomplete read always; a complete one when it has more than max_xa_hits + 1 loci).
void aln_meth_trim(aln_meth_read_t *, int);
// Moves every occ locus to hits, in order (main then indexes hits), and frees occ. No-op when occ
// is empty, as it is once a read is finished.
void aln_meth_store_loci(aln_meth_read_t *);
void aln_meth_redraw(aln_meth_read_t *, uint64_t);
// aln_meth_redraw over a finished read's hits: the same draw (one shared body), the same rng.
void aln_meth_redraw_hits(aln_meth_read_t *, uint64_t);
/* Fills a duplicate read's result from its representative's finished, not-yet-trimmed src (its
 * loci stored as hits): every field but the work counters is copied and the main hit is redrawn
 * with rng among src's best hits. Unless keep_all, a result aln_meth_trim(max_xa_hits) would cut
 * to its main hit gets only that hit (any max_xa_hits < 0 prints no XA, so it always does). src is
 * unchanged on return. ord: the copy's ordinal; a copy src's search drew for (a repeat
 * certificate's copy_draw) takes that draw instead. */
void aln_meth_copy_dup(aln_meth_read_t *dst, aln_meth_read_t *src, uint64_t rng, int keep_all, int max_xa_hits,
                       int64_t ord = -1);
uint64_t aln_meth_draw_seed(const aln_meth_ref_t *, int, int64_t);
void aln_meth_flush_parked(aln_meth_group_t *);
}
#endif
