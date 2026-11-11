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
#include <vector>
#include "bwtgap.h"
namespace bwa_aln {
struct aln_meth_pend_t { bwtint_t k, l; gap_entry_t e; };
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
struct aln_meth_iv_key_t {
    bwtint_t k, l; int span;
    bool operator==(const aln_meth_iv_key_t &o) const { return k == o.k && l == o.l && span == o.span; }
    uint64_t hash() const { return aln_flat_mix((uint64_t)k * 0x9e3779b97f4a7c15ULL ^ (uint64_t)l ^ ((uint64_t)span << 48)); }
};
struct aln_meth_win_key_t {
    int rid, hyp, strand; int64_t rb, re;
    bool operator==(const aln_meth_win_key_t &o) const {
        return rid == o.rid && hyp == o.hyp && strand == o.strand && rb == o.rb && re == o.re;
    }
    uint64_t hash() const {
        return aln_flat_mix((uint64_t)rb * 0x9e3779b97f4a7c15ULL ^ (uint64_t)re ^ ((uint64_t)(hyp * 2 + strand) << 60) ^ (uint64_t)rid << 40);
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
    uint64_t rng = 0, alignment_bytes = 0;
    int64_t unrep_cost = INT64_MAX; // the lowest cost of a refined window dropped as unrepresentable
    std::vector<ubyte_t> vseq, rd_fwd;
    std::vector<aln_meth_pend_t> pend;
    aln_flat_index_t<aln_meth_iv_key_t> intervals;
    // Exact windows dedup before DP; every distinct span is evaluated before locus collapse.
    aln_flat_index_t<aln_meth_win_key_t> windows;
    aln_flat_index_t<aln_meth_locus_key_t> loci; // value: the locus's index in res->occ
    // The read's two query orientations as bit planes (rd_fwd, then p->rseq), built on first use.
    aln_meth_read_bits_t bits[2];
    int bits_ready = 0;
    // --known-snvs: candidate windows from the ALT-window search, resolved before the read finishes
    const std::vector<aln_meth_candidate_t> *extra = nullptr;
    // Scratch: the reference window of the candidate being resolved (refined, then accounted).
    std::vector<uint8_t> win;
};
struct aln_meth_group_t {
    const bwt_t *bwt = nullptr;
    const aln_meth_ref_t *ref = nullptr;
    const aln_meth_opt_t *mo = nullptr;
    const gap_opt_t *opt = nullptr;
    int read_no = 1, keep_strata = 0, max_xa_hits = 3;
    int n_slots = 0, n_active = 0, n_parked = 0, k_park = 1;
    std::vector<aln_meth_slot_t> slot;
    gap_slot_core_t *cores = nullptr;
};
// These guards do not depend on candidates that the next flush might discover. A seeded read
// (longer than -l) can stop only once its stratum fits the seed's difference limit (-k), which
// aln_meth_can_stop checks against the best cost; -k 0 never can.
static inline bool aln_meth_stop_eligible(const aln_meth_slot_t *x) {
    const gap_opt_t &opt=x->job->opt;
    return opt.s_mm>0 && opt.s_gape>=0 && opt.s_gapo>opt.s_mm &&
        !(opt.mode&BWA_MODE_NONSTOP) && (!x->job->seed_width || opt.max_seed_diff>0) && !x->res->flags;
}
bool aln_meth_can_stop(const aln_meth_slot_t *, int next_score);
// Keep batches bounded, but only refine score boundaries when stopping is possible.
static inline int aln_meth_needs_flush(aln_meth_slot_t *x, int score) {
    if (x->pend.size() >= 32) return 1;
    if (score <= x->flushed_score) return 0;
    if (aln_meth_stop_eligible(x)) return 1;
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
    if (aln_meth_can_stop(x,score)) return 1;
    return 0;
}
static inline int aln_meth_resume(void *, gap_job_t *, gap_entry_t *) { return 1; }
void aln_meth_hit(aln_meth_slot_t *, const gap_entry_t *, bwtint_t, bwtint_t);
static inline void aln_meth_hit(void *, const gap_entry_t *, bwtint_t, bwtint_t) {}
// has_copy: the read has duplicate copies that will redraw over its loci.
const ubyte_t *aln_meth_slot_start(aln_meth_slot_t *, bwa_seq_t *, aln_meth_read_t *, int64_t, int, int has_copy);
void aln_meth_slot_finish(aln_meth_slot_t *);
void aln_meth_report(int, int, const aln_meth_read_t *);
// Finishes a read (stores its occ as hits; see aln_meth_store_loci), then trims it to its main hit
// (an incomplete read always; a complete one when it has more than max_xa_hits + 1 loci).
void aln_meth_trim(aln_meth_read_t *, int);
// Moves every occ locus to hits, in order (main then indexes hits), and frees occ. No-op when occ
// is empty, as it is once a read is finished.
void aln_meth_store_loci(aln_meth_read_t *);
void aln_meth_redraw(aln_meth_read_t *, uint64_t);
// aln_meth_redraw over a finished read's hits (a duplicate copy's draw): the same draw, the same
// rng. A copy of aln_meth_redraw's logic, whose body the search kernel's code depends on: change both.
void aln_meth_redraw_hits(aln_meth_read_t *, uint64_t);
uint64_t aln_meth_draw_seed(const aln_meth_ref_t *, int, int64_t);
void aln_meth_flush_parked(aln_meth_group_t *);
}
#endif
