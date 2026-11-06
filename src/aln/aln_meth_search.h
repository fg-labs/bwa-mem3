/* SPDX-License-Identifier: MIT */
/* Internal candidate generation hooks; no search ancestry is retained. */
#ifndef BWA_MEM3_ALN_METH_SEARCH_H
#define BWA_MEM3_ALN_METH_SEARCH_H
#include <map>
#include <set>
#include <tuple>
#include "aln_meth.h"
#include "bwtgap.h"
namespace bwa_aln {
struct aln_meth_pend_t { bwtint_t k, l; gap_entry_t e; };
struct aln_meth_group_t;
struct aln_meth_slot_t {
    aln_meth_group_t *grp = nullptr;
    gap_job_t *job = nullptr;
    bwa_seq_t *p = nullptr;
    aln_meth_read_t *res = nullptr;
    int len = 0, keep = 0, final_flush = 0, stop = 0;
    uint64_t rng = 0, alignment_bytes = 0;
    std::vector<ubyte_t> vseq, rd_fwd;
    std::vector<aln_meth_pend_t> pend;
    std::set<std::tuple<bwtint_t,bwtint_t,int>> intervals;
    // Exact windows dedup before DP; every distinct span is evaluated before locus collapse.
    std::set<std::tuple<int,int,int,int64_t,int64_t>> windows;
    std::map<std::tuple<int,int,int,int64_t>, size_t> loci;
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
// Flush bounded pending intervals before popping the next search entry; no rank-based pruning.
static inline int aln_meth_needs_flush(const aln_meth_slot_t *x, int) { return x->pend.size() >= 32; }
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
    return 0;
}
static inline int aln_meth_resume(void *, gap_job_t *, gap_entry_t *) { return 1; }
void aln_meth_hit(aln_meth_slot_t *, const gap_entry_t *, bwtint_t, bwtint_t);
static inline void aln_meth_hit(void *, const gap_entry_t *, bwtint_t, bwtint_t) {}
const ubyte_t *aln_meth_slot_start(aln_meth_slot_t *, bwa_seq_t *, aln_meth_read_t *, int64_t, int, int);
void aln_meth_slot_finish(aln_meth_slot_t *);
void aln_meth_report(int, int, const aln_meth_read_t *);
void aln_meth_trim(aln_meth_read_t *, int);
void aln_meth_redraw(aln_meth_read_t *, uint64_t);
uint64_t aln_meth_draw_seed(const aln_meth_ref_t *, int, int64_t);
void aln_meth_flush_parked(aln_meth_group_t *);
}
#endif
