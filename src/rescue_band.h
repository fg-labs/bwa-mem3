/* Exact banded mate-rescue DP (passes 0 and 1 of the 8-bit kswv rescue path).
 *
 * rescue_prune.h narrows a rescue window to the hull of the 5-mer diagonal components whose bound
 * reaches minsc. This file goes further: it runs the pass-0 DP only inside the diagonal BANDS of
 * those components, and reassembles exactly the kswv outputs {score, te, qe, score2} of the hull
 * job from the bands' per-row maxima; pass 1 (start recovery, tb / qb) then runs in a band derived
 * from the pass-0 score (below). mem_matesw_batch_post is unchanged: the parent SeqPair is still
 * the hull window.
 *
 * Why it is exact (test/rescue_band_harness.cpp checks every output field against kswv on the full
 * window, on generated adversarial jobs and on real rescue jobs):
 *  - Every alignment scoring >= tau lies inside ONE diagonal component at threshold tau (the
 *    K-mer bound decays by c per diagonal while a gap excursion costs more per diagonal), and
 *    starts at or after the component's first hit row - (K - 1), and ends by row
 *    i1(tau) = dmaxhit + quanta - 1 + tail(ub, tau) (rescue_prune_params: max(0, ub - tau - 2) at
 *    the default scoring).
 *  - A zero-state DP restricted to any cell set (a band) never exceeds the full DP, and is exact
 *    on every cell of every alignment it contains. So per-row max over bands is <= the true row
 *    max, and equal wherever an alignment scoring >= tau ends. Widening a band (more diagonals,
 *    as grouping lanes into one vector does) keeps both properties, and so does skipping cells that do not
 *    exist in the full DP (query column j < 0 or j >= quanta: the kernel computes, per row, only
 *    the union of its lanes' live columns). The fused (G-based) cells drop an in-row gap run
 *    followed directly by a vertical one; the band then need not match the original cell's
 *    values, but it still holds the full-DP value at the end of every alignment above, because
 *    such an alignment has an equal-scoring twin with the two runs swapped that lies in the same
 *    band (rescue_band.cpp, the note above RB_CELL1).
 *  - Round 1 computes the components at T1 = max(minsc, ub2 / 2) (ub2 = second-largest component
 *    bound at minsc), or at ub1 - delta for a lone near-perfect primary (plan()). The merged
 *    S / te / qe / score2 are exact iff T1 == minsc, or S >= T1 and either score2 >= T1 -- the
 *    termination test --
 *    or the whole hull lies inside the zone [te - Z, te + Z], Z = ceil(S / a) (every row with
 *    R >= minsc is in the hull, so no b[] anchor is out of zone and score2 is -1). Otherwise round 2
 *    computes every component at minsc (skipping those confined to (te - Z, te) with bound < S) and
 *    is exact unconditionally.
 *  - score2 is reproduced from the merged row maxima with kswv's own emulation (lagged rising-row
 *    zeroing, kswv.cpp:~932-936, then the b[] scan, kswv.cpp:~1436-1519); the result
 *    depends only on rows the rule above keeps exact.
 *
 * Pass 1 (tb, qb) of every 8-bit job is banded too (take_pass1 / run_pass1),
 * whether its pass 0 was banded or ran through kswv: the argument below needs only the exact
 * pass-0 (S, te, qe) and the scoring, not the pass-0 components. kswv phase 1 runs the
 * DP on the reversed prefixes ref[te..0] x q[qe..0] and stops at the first row whose max reaches S;
 * tb = te - that row, qb = qe - the first column holding S in it. Every alignment scoring S inside
 * that rectangle ends at (te, qe) (an earlier end row would contradict te being the first row of S
 * in pass 0, an earlier end column in row te the same for qe), so the answer is the start of the
 * score-S alignment A* ending at (te, qe) with the latest start row, then the latest start column.
 * A* with M matches, X mismatches, I inserted and D deleted bases has M + X + I <= qe + 1 (its
 * query span) and scores S <= a M - b X - (gap costs), so a deletion run gives
 * S <= a (qe + 1) - o_del - e_del D and an insertion run S <= a (qe + 1 - I) - o_ins - e_ins I;
 * hence D <= Dmax = max(0, (a (qe + 1) - S - o_del) / e_del), I <= Imax =
 * max(0, (a (qe + 1) - S - o_ins) / (a + e_ins)), and it spans at most qe + 1 + Dmax rows (at the
 * default scoring Dmax = max(0, qe - S - 5) and Imax = Dmax / 2). In reversed
 * coordinates A* starts at (0, 0) and every cell lies on a diagonal r - c in [-Imax, Dmax]. A
 * zero-state DP restricted to that band (cells outside nonexistent: H = E = F = 0) never exceeds
 * kswv's reversed DP and holds S at A*'s end, so it has the same first row reaching S and, in that
 * row, the same first column holding S: tb / qb are byte-identical, ties included. The pad columns
 * [qe + 1, quanta) are in both DPs; they cannot hold S on the first row reaching S (a pad cell's
 * value comes from a row above or, reduced by a gap, from its own row), so they never decide
 * either argmax. The band stops a lane group once every lane has reached its S (gmax cannot exceed
 * S, so nothing after that row changes te or qe), like kswv's KSW_XSTOP freeze. Jobs whose band
 * is not cheaper than kswv, or whose banded max is not S (cannot happen; kept as a guard), run
 * kswv phase 1.
 *
 * Scope: the NEON kernel (aarch64, 16 lanes) and the AVX2 kernel (x86, 32 lanes,
 * rescue_band_kernel_x86.h), for any scoring the kernels take (rb_scoring: the score table
 * {a, -b, -1}, separate deletion and insertion gap costs), 8-bit and non-meth: pass 0 for the jobs
 * rescue_prune_applies() prunes (minsc in [5, 255], plan()), pass 1 for every such job. Without a
 * SIMD kernel the hull path runs.
 * Env: the BWA3_RESCUE_BAND* knobs and BWA3_RESCUE_PRUNE_STATS, listed with their defaults in
 * rescue_env.h.
 *
 * Overview and gates: docs/src/developer-guide/rescue-banding.md. */
#ifndef BWA_MEM3_RESCUE_BAND_H
#define BWA_MEM3_RESCUE_BAND_H

#include <stdint.h>
#include <vector>
#include "kswv.h"
#include "rescue_prune.h"

/* The scoring the band kernels run with (per batch: every pair of a batch shares mem_opt_t): match
 * a, mismatch b (N scores -1), gap of length L costs o + e L per type (deletion: reference advances,
 * the kernels' vertical E; insertion: query advances, their in-row F). shift = kswv's 8-bit bias,
 * -min(a, -b, -1) = max(1, b): a pass-0 score with S + shift >= 255 is kswv's saturated 255. */
struct rb_scoring {
    int a = 1, b = 4, o_del = 6, e_del = 1, o_ins = 6, e_ins = 1;
    int shift() const { return b > 1 ? b : 1; }
    /* Deletion and insertion alike: one gap subtraction serves both (the kernels' Sym form). */
    bool sym_gaps() const { return o_del + e_del == o_ins + e_ins && e_del == e_ins; }
    /* score2's zone half-width around te: kswv's ceil(S / qmax), qmax = a (kswv.cpp). */
    int zone(int S) const { return (S + a - 1) / a; }
    bool operator==(const rb_scoring &o) const
    {
        return a == o.a && b == o.b && o_del == o.o_del && e_del == o.e_del && o_ins == o.o_ins && e_ins == o.e_ins;
    }
    /* What the kernels need: positive a, b and extends, non-negative opens, and kswv's 8-bit score
     * table (+a, -b, a + shift <= 255): the first check of rescue_prune_params::from. */
    bool valid() const
    {
        return a >= 1 && b >= 1 && e_del >= 1 && e_ins >= 1 && o_del >= 0 && o_ins >= 0 && a <= 127 && b <= 128
               && a + shift() <= 255;
    }
    static rb_scoring from(const rescue_prune_params &p)
    {
        rb_scoring s;
        s.a = p.a; s.b = p.b; s.o_del = p.o_del; s.e_del = p.e_del; s.o_ins = p.o_ins; s.e_ins = p.e_ins;
        return s;
    }
};

/* A diagonal component at some threshold, in full-window coordinates (unshifted diagonals
 * d = i - j). ub = max interval bound; i0 = min over hit diagonals of (first hit end row - (K - 1));
 * dmaxhit = the component's highest diagonal holding a hit. */
struct rb_comp { int ub, i0, dlo, dhi, dmaxhit; };

/* A band job in parent (hull-window) coordinates: rows [r0, r1], diagonals [dlo, dhi] (d = r - j). */
struct rb_band { int32_t r0, r1, dlo, dhi, ub; };

/* Enumerate the components of a filter view at threshold tau within diagonal indices [xa, xb)
 * (shifted, x = d + off). Appends to out; returns false if more than cap would be produced.
 * Exported for the exactness harness, which checks that the NEON and scalar views of one job give
 * the same components. */
bool rescue_band_components(const rescue_prune_view &v, int tau, int xa, int xb,
                            std::vector<rb_comp> &out, int cap);

struct rescue_band_stats {
    uint64_t parents = 0, bands1 = 0, bands2 = 0, r2_band = 0, r2_kswv = 0, single = 0;
    uint64_t cells_req = 0, cells_pad = 0, hull_cells = 0, groups = 0, planned_no = 0, comp_cap = 0, tight = 0;
    uint64_t p1_band = 0, p1_kswv = 0, p1_guard = 0;
};

/* Per-batch banding state, one per tid (mem_cache::rescue_band). mem_matesw_batch_pre plans (plan +
 * commit per enqueued pair), mem_sam_pe_batch runs pass 0 for the banded parents, then reset(). Every buffer is
 * grow-only and reused across batches: nothing is allocated per job. */
class RescueBandBatch {
public:
    /* Plan the bands of a B2 hull [hb, he] of a full window of len1 rows against a query of len2.
     * v must be the view rescue_prune_window returned with that decision under p (the scoring and
     * threshold), which must be the batch's scoring (set_scoring) and without --meth
     * (rescue_prune_params::band_ok); the caller guarantees the 8-bit path (rescue_prune_applies).
     * Returns true when banding is chosen (the cost model says it beats kswv on the hull);
     * commit() then binds it to the pair's regid. */
    bool plan(const rescue_prune_view &v, const rescue_prune_params &p, int len1, int len2, int hb, int he);
    /* The scoring of this batch's pass-0 and pass-1 band jobs (valid()); set before plan() /
     * take_pass1() and kept across reset(). */
    void set_scoring(const rb_scoring &sc) { sc_ = sc; }
    /* Bind the pending plan (or "not banded") to regid. Call for every enqueued pair. */
    void commit(int regid);
    bool banded(int regid) const
    {
        return regid >= 0 && regid < (int)regid2parent_.size() && regid2parent_[regid] >= 0;
    }
    /* Stable-partition pairs[0, n) so the banded pairs come first; returns their count. */
    int partition(SeqPair *pairs, int n);
    /* Pass 0 for the banded parents pairs[0, nb): fills aln[regid].{score, te, qe, score2, te2}
     * exactly as kswv phase 0 on the hull would. Parents that need a round 2 which is not banded
     * run through kswv (phase 0) on the hull. */
    void run_pass0(const SeqPair *pairs, int nb, const uint8_t *seqBufRef,
                   const uint8_t *seqBufQer, kswr_t *aln, Ikswv *kswv);
    /* Pass 1 of an 8-bit rescue job. sp is the pair as prepared for kswv phase 1 (reversed
     * prefixes of lengths te + 1 and qe + 1 in the sequence buffers, len2 = qe + 1, h0 =
     * KSW_XSTOP | S) and r its pass-0 result, scored with set_scoring()'s scoring. The CALLER
     * guarantees that scoring and non-meth (a banded parent implies both). Queues the job for the banded pass 1 and returns true when
     * BWA3_RESCUE_BAND_P1 admits it and the band is cheaper than kswv; false means the caller
     * runs it through kswv phase 1 as before. */
    bool take_pass1(const SeqPair &sp, const kswr_t &r, bool banded_parent);
    /* Run the queued pass-1 jobs: fills aln[regid].{tb, qb} exactly as kswv phase 1 would. */
    void run_pass1(const uint8_t *seqBufRef, const uint8_t *seqBufQer, kswr_t *aln, Ikswv *kswv);
    void reset();
    /* Counters of the current batch (folded into the process totals and cleared by reset() when
     * BWA3_RESCUE_PRUNE_STATS=1). */
    const rescue_band_stats &stats() const { return stats_; }

private:
    rescue_band_stats stats_;
    struct parent_rec { int32_t b0, nb, c0, nc, T1, minsc; };   // the batch's scoring is sc_
    struct job {
        const uint8_t *ref, *qry;
        int32_t len2, quanta, r0, nrows, dlo, w, parent, key, target;
    };
    struct pstate {
        int32_t regid, L, S, te, qe, s2, te2, rbuf, rec, round;
    };
    /* A queued pass-1 job: the prepared pair, its reversed row count and band [-imax, dmax]. */
    struct p1job {
        SeqPair sp;
        int32_t nrows, imax, dmax;
    };
    /* Per-lane result of one band job: gmax, its first row (-1 if gmax == 0) and first column. */
    struct lane_res { int32_t g, te, qe; };
    void run_jobs(bool pass1);
    void finish_parent(pstate &p);
    std::vector<int32_t> regid2parent_;
    std::vector<parent_rec> recs_;
    std::vector<rb_band> bands_;
    std::vector<rb_comp> c19_, cT_;
    std::vector<job> jobs_;
    std::vector<int32_t> order_;
    std::vector<pstate> ps_;
    std::vector<uint8_t> rpool_;
    std::vector<SeqPair> spscratch_;
    std::vector<int32_t> kswv_list_, r2buf_;
    std::vector<rb_band> kept_;
    std::vector<p1job> p1_;
    std::vector<lane_res> p1res_;
    int pending_ = -1;
    rb_scoring sc_;
};

/* Free a batch made by rescue_band_batch_new (NULL is a no-op). The aligner keeps one per tid in
 * mem_cache (rescue_band, bwamem.h), next to the batch's other per-tid state. */
RescueBandBatch *rescue_band_batch_new();
void rescue_band_batch_free(RescueBandBatch *b);
bool rescue_band_enabled();
/* The default of the pass-0 cost gate (BWA3_RESCUE_BAND_COST) when kswv runs at SIMD tier `tier`
 * (simd_dispatch.h): 0 (no pass-0 banding) at BWAMEM3_TIER_AVX512BW, where kswv sweeps 64 lanes
 * against the band kernel's 32, and 85 at every other tier. See rb_cost_pct in rescue_band.cpp. */
int rescue_band_cost_pct_default(int tier);

#endif
