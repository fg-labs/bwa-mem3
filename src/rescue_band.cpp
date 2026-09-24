/* Exact banded mate-rescue DP: component planning, the 16-lane NEON banded kernel, grouping and
 * the per-parent merge. See rescue_band.h for the design and the exactness argument. */
#include "rescue_band.h"
#include "rescue_env.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#if defined(__aarch64__)
#include <arm_neon.h>
#endif

/* ------------------------------------------------------------------------------------------ */
/* Toggles and tuning (read once; the full BWA3_RESCUE_* list is in rescue_env.h).             */
/* ------------------------------------------------------------------------------------------ */

bool rescue_band_enabled()
{
#if defined(__aarch64__)
    static const bool on = rescue_env_flag("BWA3_RESCUE_BAND");
    return on;
#else
    return false;
#endif
}

/* Band a parent only when  100 * sum_bands rows * (min(width, quanta) + RB_OVH)  <  pct * hull_rows *
 * quanta, i.e. when the band cells (per row the kernel computes at most ~quanta of them, see
 * rb_dp_core) undercut the hull. The banded DP and kswv cost about the same per cell, so pct < 100
 * covers the band path's fixed costs (SoA build, merge, score2). BWA3_RESCUE_BAND_COST overrides
 * pct; 85 was the best of {40, 55, 70, 85, 100, 130} locally (M-series, wgs-like HG002). */
static int rb_cost_pct()
{
    static const int v = rescue_env_int("BWA3_RESCUE_BAND_COST", 85);
    return v;
}
static const int RB_OVH = 8;          // per-row fixed cost of a band lane, in cell units
static const int RB_COMP_CAP = 64;    // more components than this: keep the hull (kswv)

/* Round 2 (the few parents whose round-1 result is not provably final): banded at minsc (default)
 * or kswv on the hull (BWA3_RESCUE_BAND_R2=0). Both are exact. */
static bool rb_r2_banded()
{
    static const bool on = rescue_env_flag("BWA3_RESCUE_BAND_R2");
    return on;
}

/* Tight top band (see plan()): T1 = ub1 - delta when the hull is predicted to fit inside the
 * zone [te - S, te + S]. BWA3_RESCUE_BAND_TIGHT = delta (default below); 0 disables. */
static int rb_tight_delta()
{
    static const int v = rescue_env_int("BWA3_RESCUE_BAND_TIGHT", 8);
    return v;
}

/* Which pass-1 jobs run banded (BWA3_RESCUE_BAND_P1): 0 none (kswv phase 1 for all), 1 banded
 * parents only, 2 (default) every eligible 8-bit job. All are exact; see rescue_band.h. */
static int rb_p1_mode()
{
    static const int v = rescue_env_int("BWA3_RESCUE_BAND_P1", 2);
    return v;
}

/* Pass-1 cost model (take_pass1): band iff 100 * (min(width, quanta) + RB_OVH) < pct * quanta.
 * Separate from pass 0's pct because the two compare different things (a band against kswv on the
 * same reversed rows, versus bands against the whole hull). The default 130 bands every job with
 * quanta >= 32 even at full width: with the early exit the band kernel beat kswv phase 1 at every
 * width locally (M-series, r set: pass-1 thread-s 0.28-0.32 at 130 vs 0.35-0.45 at 85 and 0.80-0.93
 * with kswv only), presumably because its groups are sorted by width and rows while a kswv group
 * runs to its slowest lane. Retune on Graviton; BWA3_RESCUE_BAND_P1_COST overrides. */
static int rb_p1_cost_pct()
{
    static const int v = rescue_env_int("BWA3_RESCUE_BAND_P1_COST", 130);
    return v;
}

/* BWA3_RESCUE_FSCAN selects the kernel: 0 = the original cell (rb_dp_core<false>), 1 = the fused
 * G-based cell of rb_dp_core<true> (see RB_CELL1), 2 (default, or unset) = rb_dp_wave2, the fused
 * cell on two rows per step with the direct qe scan. All three are byte-identical in every value the
 * caller reads. The band path runs only at default scoring, so the insertion and deletion gap costs
 * always agree, as the fused cell requires. */
static int rb_kernel()
{
    static const int v = rescue_env_int("BWA3_RESCUE_FSCAN", 2);
    return v;
}

/* BWA3_RESCUE_BAND_SHIFT=0 disables the per-lane band shift of run_jobs (see there); for A/B. */
static bool rb_shift_on()
{
    static const bool on = rescue_env_flag("BWA3_RESCUE_BAND_SHIFT");
    return on;
}

static bool rb_stats_on()
{
    static const bool on = rescue_env_opt_in("BWA3_RESCUE_PRUNE_STATS");
    return on;
}

/* Process-wide totals, folded in by reset() when stats are on and printed at exit. */
namespace {
struct rb_global_stats {
    std::atomic<uint64_t> parents{0}, bands1{0}, bands2{0}, r2_band{0}, r2_kswv{0}, single{0};
    std::atomic<uint64_t> cells_req{0}, cells_pad{0}, hull_cells{0}, groups{0}, planned_no{0}, comp_cap{0}, tight{0};
    std::atomic<uint64_t> p1_band{0}, p1_kswv{0}, p1_guard{0};
    std::atomic<uint64_t> t_plan_ns{0}, t_band_ns{0};
    ~rb_global_stats()
    {
        if (!rb_stats_on()) return;
        fprintf(stderr, "[RESCUE_BAND] banded_parents=%llu tight_T1_planned=%llu declined_by_cost=%llu comp_cap=%llu "
                        "bands_r1=%llu single_band=%llu round2_banded=%llu round2_kswv=%llu bands_r2=%llu "
                        "groups=%llu cells_req=%llu cells_padded=%llu hull_cells_replaced=%llu "
                        "pass1_banded=%llu pass1_kswv=%llu pass1_guard_fallback=%llu "
                        "plan_s=%.3f band_s=%.3f\n",
                (unsigned long long)parents, (unsigned long long)tight, (unsigned long long)planned_no,
                (unsigned long long)comp_cap, (unsigned long long)bands1, (unsigned long long)single,
                (unsigned long long)r2_band, (unsigned long long)r2_kswv, (unsigned long long)bands2,
                (unsigned long long)groups, (unsigned long long)cells_req,
                (unsigned long long)cells_pad, (unsigned long long)hull_cells,
                (unsigned long long)p1_band, (unsigned long long)p1_kswv, (unsigned long long)p1_guard,
                t_plan_ns * 1e-9, t_band_ns * 1e-9);
    }
};
rb_global_stats g_rb_stats;
}  // namespace

RescueBandBatch &rescue_band_batch()
{
    static thread_local RescueBandBatch b;
    return b;
}

/* ------------------------------------------------------------------------------------------ */
/* Components                                                                                  */
/* ------------------------------------------------------------------------------------------ */

template <class BND>
static inline bool rb_scan_run(const rescue_prune_view &v, BND bnd, int a, int b, std::vector<rb_comp> &out, int cap)
{
    rb_comp k{0, INT_MAX, a - v.off, b - 1 - v.off, INT_MIN};
    for (int x = a; x < b; x++) {
        const int bx = bnd(x);
        if (bx > k.ub) k.ub = bx;
        if (v.cnt[x]) {
            if (v.minrow[x] < k.i0) k.i0 = v.minrow[x];
            k.dmaxhit = x - v.off;
        }
    }
    if (k.dmaxhit == INT_MIN) return true;   // no hit: dropped (optsim comps_tau)
    if ((int)out.size() >= cap) return false;
    out.push_back(k);
    return true;
}

template <class BND>
static bool rb_components_scalar(const rescue_prune_view &v, BND bnd, int tau, int xa, int xb,
                                 std::vector<rb_comp> &out, int cap)
{
    for (int x = xa; x < xb;) {
        if (bnd(x) < tau) { x++; continue; }
        const int a = x;
        while (x < xb && bnd(x) >= tau) x++;
        if (!rb_scan_run(v, bnd, a, x, out, cap)) return false;
    }
    return true;
}

#if defined(__aarch64__)
/* NEON view (bnd16 + the filter's bitsets). mw = {bnd >= 19}, hw = {bnd >= 19 and cnt > 0}.
 * Every component at tau >= 19 lies inside one run of mw, and inside such a run hw is exactly the
 * "has a hit" set, so dmaxhit is the last hw bit of the component. ub and i0 are vector max / masked
 * min over the component's diagonals (bnd16, cnt and minrow are readable up to the filter's
 * rounded-up diagonal count, so 8-wide loads never leave the arrays; lanes past b are masked). */
static inline void rb_run_ub_i0(const rescue_prune_view &v, int a, int b, int &ub, int &i0)
{
    static const int16_t iota[8] = {0, 1, 2, 3, 4, 5, 6, 7};
    const int16x8_t IO = vld1q_s16(iota);
    const int16x8_t big = vdupq_n_s16(32767), small = vdupq_n_s16(-32768);
    int16x8_t mx = small, mn = big, mx2 = small, mn2 = big;
    /* one block of 8 diagonals; in = lanes inside [a, b) */
    auto edge = [&](int x) {
        const int16x8_t pos = vaddq_s16(IO, vdupq_n_s16((int16_t)x));
        const uint16x8_t in = vandq_u16(vcgeq_s16(pos, vdupq_n_s16((int16_t)a)), vcltq_s16(pos, vdupq_n_s16((int16_t)b)));
        const uint16x8_t c = vld1q_u16(v.cnt + x);
        mx = vmaxq_s16(mx, vbslq_s16(in, vld1q_s16(v.bnd16 + x), small));
        mn = vminq_s16(mn, vbslq_s16(vandq_u16(in, vtstq_u16(c, c)), vld1q_s16(v.minrow + x), big));
    };
    /* inner blocks lie wholly inside [a, b): two independent accumulator pairs */
    auto inner = [&](int x, int16x8_t &MX, int16x8_t &MN) {
        const uint16x8_t c = vld1q_u16(v.cnt + x);
        MX = vmaxq_s16(MX, vld1q_s16(v.bnd16 + x));
        MN = vminq_s16(MN, vbslq_s16(vtstq_u16(c, c), vld1q_s16(v.minrow + x), big));
    };
    const int x0 = a & ~7, x1 = b & ~7;
    edge(x0);
    int x = x0 + 8;
    for (; x + 8 < x1; x += 16) { inner(x, mx, mn); inner(x + 8, mx2, mn2); }
    if (x < x1) inner(x, mx, mn);
    if (x1 > x0 && x1 < b) edge(x1);
    mx = vmaxq_s16(mx, mx2);
    mn = vminq_s16(mn, mn2);
    ub = vmaxvq_s16(mx);
    i0 = vminvq_s16(mn);
}

static bool rb_components_neon(const rescue_prune_view &v, int tau, int xa, int xb,
                               std::vector<rb_comp> &out, int cap)
{
    using rescue_prune_neon::neon_next;
    using rescue_prune_neon::neon_last;
    auto emit = [&](int a, int b) -> bool {
        const int dm = neon_last(v.hw, a, b);
        if (dm < 0) return true;   // no hit: dropped
        if ((int)out.size() >= cap) return false;
        rb_comp k;
        rb_run_ub_i0(v, a, b, k.ub, k.i0);
        k.dlo = a - v.off; k.dhi = b - 1 - v.off; k.dmaxhit = dm - v.off;
        out.push_back(k);
        return true;
    };
    if (tau == rescue_prune_neon::MINSC) {
        /* The whole view at MINSC: the filter already listed these components (same runs of mw,
         * same ub / i0 / dmaxhit computation, rescue_prune_neon.h step 7). Reproduce emit()'s
         * result: append until out holds cap, false iff any component did not fit. */
        if (xa == 0 && xb == v.nd && v.ncomp >= 0 && (v.ncomp == v.ncomp_stored || (int)out.size() + v.ncomp_stored >= cap)) {
            for (int k = 0; k < v.ncomp; k++) {
                if ((int)out.size() >= cap) return false;
                const rescue_prune_neon::Comp &c = v.comps[k];
                rb_comp K;
                K.ub = c.ub; K.i0 = c.i0; K.dlo = c.a - v.off; K.dhi = c.b - 1 - v.off; K.dmaxhit = c.dmax - v.off;
                out.push_back(K);
            }
            return true;
        }
        for (int d = neon_next(v.mw, xa, xb, 0); d < xb;) {
            const int b = neon_next(v.mw, d, xb, ~0ull);
            if (!emit(d, b)) return false;
            d = neon_next(v.mw, b, xb, 0);
        }
        return true;
    }
    /* tau > 19: runs of bnd >= tau inside [xa, xb) (a run of mw), 64 diagonals per bitmask word;
     * run starts and ends are the bit transitions of the word (branch per run, not per diagonal). */
    static const uint8_t bw[16] = {1, 2, 4, 8, 16, 32, 64, 128, 1, 2, 4, 8, 16, 32, 64, 128};
    const uint8x16_t BW = vld1q_u8(bw);
    const int16x8_t T = vdupq_n_s16((int16_t)tau);
    int run = -1;   // start of the open run, or -1
    for (int x = xa; x < xb; x += 64) {
        /* bnd16 holds nd + 64 readable entries (rescue_prune_neon.h), so the word's loads stay in
         * the array; bits at or past xb are masked off below */
        const int16_t *src = v.bnd16 + x;
        uint8x16_t g[4];
        for (int q = 0; q < 4; q++) {
            const uint16x8_t c0 = vcgeq_s16(vld1q_s16(src + 16 * q), T), c1 = vcgeq_s16(vld1q_s16(src + 16 * q + 8), T);
            g[q] = vandq_u8(vcombine_u8(vmovn_u16(c0), vmovn_u16(c1)), BW);
        }
        const uint8x16_t p = vpaddq_u8(vpaddq_u8(g[0], g[1]), vpaddq_u8(g[2], g[3]));
        uint64_t m = vgetq_lane_u64(vreinterpretq_u64_u8(vpaddq_u8(p, p)), 0);   // bit k: x + k
        if (xb - x < 64) m &= (1ull << (xb - x)) - 1;
        uint64_t t = m ^ ((m << 1) | (run >= 0 ? 1u : 0u));   // bit k: the state changes at x + k
        while (t) {
            const int k = __builtin_ctzll(t);
            t &= t - 1;
            if (run < 0) run = x + k;
            else { if (!emit(run, x + k)) return false; run = -1; }
        }
    }
    if (run >= 0 && !emit(run, xb)) return false;
    return true;
}
#endif

bool rescue_band_components(const rescue_prune_view &v, int tau, int xa, int xb,
                            std::vector<rb_comp> &out, int cap)
{
#if defined(__aarch64__)
    if (v.bnd16 && v.mw && v.hw && tau >= rescue_prune_neon::MINSC) return rb_components_neon(v, tau, xa, xb, out, cap);
#endif
    if (v.bnd16) {
        const int16_t *b16 = v.bnd16;
        return rb_components_scalar(v, [b16](int x) { return (int)b16[x]; }, tau, xa, xb, out, cap);
    }
    const int32_t *fw = v.fwd, *bw = v.bwd;
    const uint16_t *cn = v.cnt;
    return rb_components_scalar(v, [fw, bw, cn](int x) { return 5 + fw[x] + bw[x] - ((int)cn[x] - 1); },
                                tau, xa, xb, out, cap);
}

/* ------------------------------------------------------------------------------------------ */
/* score2 (kswv emulation)                                                                     */
/* ------------------------------------------------------------------------------------------ */

void rescue_band_score2(const uint8_t *R, int stride, int n, int row0, int S, int te, int minsc,
                        int *score2, int *te2)
{
    const int low = te - S, high = te + S;
    int s2 = -1, t2 = -1, bs = -1, bp = -2;
    int cur = n > 0 ? R[0] : 0;
    for (int k = 0; k < n; k++) {
        const int nxt = k + 1 < n ? R[(size_t)(k + 1) * stride] : 0;
        const int v = nxt > cur ? 0 : cur;   // kswv.cpp lagged store: rising row -> 0
        cur = nxt;
        if (v < minsc) continue;
        const int i = row0 + k;
        if (bp + 1 != i) {
            if (bp >= 0 && (bp < low || bp > high) && bs > s2) { s2 = bs; t2 = bp; }
            bs = v; bp = i;
        } else if (bs < v) {
            bs = v; bp = i;
        }
    }
    if (bp >= 0 && (bp < low || bp > high) && bs > s2) { s2 = bs; t2 = bp; }
    *score2 = s2;
    *te2 = t2;
}

/* ------------------------------------------------------------------------------------------ */
/* Planning (mem_matesw_batch_pre)                                                             */
/* ------------------------------------------------------------------------------------------ */

bool RescueBandBatch::plan(const rescue_prune_view &v, int len1, int len2, int hb, int he, int minsc)
{
    pending_ = -1;
    if (v.nd < 0 || minsc < 19 || minsc > 255) return false;
    const auto t0 = rb_stats_on() ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point();
    const int quanta = ((len2 + 15) / 16) * 16;
    const int H = he - hb + 1;
    bool ok = v.off == quanta && H > 0;
    c19_.clear();
    if (ok && !rescue_band_components(v, minsc, 0, v.nd, c19_, RB_COMP_CAP)) { ok = false; stats.comp_cap++; }
    if (ok && c19_.empty()) ok = false;
    int T1 = minsc;
    if (ok) {
        int ub1 = 0, ub2 = 0;
        for (const rb_comp &K : c19_) {
            if (K.ub > ub1) { ub2 = ub1; ub1 = K.ub; }
            else if (K.ub > ub2) ub2 = K.ub;
        }
        T1 = std::max(minsc, ub2 / 2);
        /* Tight top band. When the whole hull lies inside the zone [te - S, te + S], the true
         * score2 is exactly -1 (every row with R >= minsc is in the hull, so no b[] anchor can be
         * out of zone), and round 1 is final as soon as S >= T1 (run_pass0's termination test).
         * Containment needs S >= len2 - 1 in practice (the hull starts at the primary's first hit
         * row - 4 ~ te - len2 + 1, and a mismatched primary's out-of-zone rising edge is a genuine
         * kswv score2 contributor), so predict with S = ub1 (a mismatch costs the score and the
         * bound 5 each) and only when ub1 <= len2 (above that the bound is repeat noise). A wrong
         * prediction only costs a round 2 (exact either way). */
        const int delta = rb_tight_delta();
        if (delta > 0) {
            const rb_comp *K1 = &c19_[0];
            for (const rb_comp &K : c19_) if (K.ub > K1->ub) K1 = &K;
            const int Tt = std::max(minsc, ub1 - delta);
            const int te_pred = K1->dmaxhit + len2 - 1 - hb;
            if (ub1 <= len2 && Tt > T1 && te_pred - ub1 <= 0 && te_pred + ub1 >= H - 1) { T1 = Tt; stats.tight++; }
        }
    }
    const std::vector<rb_comp> *cT = &c19_;
    if (ok && T1 > minsc) {
        cT_.clear();
        for (const rb_comp &K : c19_)
            /* K.ub is the max bound over K's diagonals, so K holds no diagonal with bound >= T1
             * when K.ub < T1: the scan would find nothing. */
            if (K.ub >= T1 && !rescue_band_components(v, T1, K.dlo + v.off, K.dhi + v.off + 1, cT_, RB_COMP_CAP)) {
                ok = false; stats.comp_cap++; break;
            }
        cT = &cT_;
    }
    /* Band of component K at threshold tau, in parent (hull) coordinates; false if it would leave
     * the hull (cannot happen: the hull is the union of the minsc extents, which contain every
     * extent at tau >= minsc -- checked anyway, and the parent then keeps the hull). */
    auto to_band = [&](const rb_comp &K, int tau, rb_band &b) -> bool {
        const int i0 = std::max(0, K.i0);
        const int i1 = std::min(len1 - 1, K.dmaxhit + quanta - 1 + std::max(0, K.ub - tau - 2));
        b.r0 = i0 - hb; b.r1 = i1 - hb;
        if (b.r0 < 0 || b.r1 >= H || b.r0 > b.r1) return false;
        b.dlo = std::max(K.dlo - hb, b.r0 - quanta + 1);   // cells below: j >= quanta, nonexistent
        b.dhi = std::min(K.dhi - hb, b.r1);                // cells above: j < 0, nonexistent
        b.ub = K.ub;
        if (b.dlo > b.dhi) return false;
        /* Rows with no live cell: r < dlo (every j < 0: H = 0, so starting at dlo is the same zero
         * state) and r > dhi + quanta - 1 (every j >= quanta). */
        b.r0 = std::max(b.r0, b.dlo);
        b.r1 = std::min(b.r1, b.dhi + quanta - 1);
        return b.r0 <= b.r1;
    };
    const size_t b0 = bands_.size();
    if (ok) {
        long cost = 0;
        for (const rb_comp &K : *cT) {
            rb_band b;
            if (!to_band(K, T1, b)) { ok = false; break; }
            bands_.push_back(b);
            /* per row the kernel computes at most min(width, quanta) cells (live-cell clipping) */
            cost += (long)(b.r1 - b.r0 + 1) * (std::min(b.dhi - b.dlo + 1, quanta) + RB_OVH);
        }
        if (ok && cost * 100 >= (long)H * quanta * rb_cost_pct()) { ok = false; stats.planned_no++; }
    }
    const size_t c0 = bands_.size();
    if (ok && T1 > minsc) {
        for (const rb_comp &K : c19_) {
            rb_band b;
            if (!to_band(K, minsc, b)) { ok = false; break; }
            bands_.push_back(b);
        }
    }
    if (!ok) {
        bands_.resize(b0);
    } else {
        parent_rec r;
        r.b0 = (int32_t)b0; r.nb = (int32_t)(c0 - b0);
        r.c0 = (int32_t)c0; r.nc = (int32_t)(bands_.size() - c0);
        r.T1 = T1; r.minsc = minsc;
        pending_ = (int)recs_.size();
        recs_.push_back(r);
    }
    if (rb_stats_on())
        stats.t_plan += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return ok;
}

void RescueBandBatch::commit(int regid)
{
    if ((int)regid2parent_.size() <= regid) regid2parent_.resize((size_t)regid + 1024, -1);
    regid2parent_[regid] = pending_;
    pending_ = -1;
}

int RescueBandBatch::partition(SeqPair *pairs, int n)
{
    if (recs_.empty()) return 0;
    if ((int)spscratch_.size() < n) spscratch_.resize(n);
    int nb = 0;
    for (int i = 0; i < n; i++) if (banded(pairs[i].regid)) nb++;
    if (nb == 0) return 0;
    int a = 0, o = nb;
    for (int i = 0; i < n; i++) spscratch_[banded(pairs[i].regid) ? a++ : o++] = pairs[i];
    memcpy(pairs, spscratch_.data(), (size_t)n * sizeof(SeqPair));
    return nb;
}

void RescueBandBatch::reset()
{
    if (rb_stats_on()) {
        g_rb_stats.parents += stats.parents; g_rb_stats.bands1 += stats.bands1;
        g_rb_stats.bands2 += stats.bands2; g_rb_stats.r2_band += stats.r2_band;
        g_rb_stats.r2_kswv += stats.r2_kswv; g_rb_stats.single += stats.single;
        g_rb_stats.cells_req += stats.cells_req; g_rb_stats.cells_pad += stats.cells_pad;
        g_rb_stats.hull_cells += stats.hull_cells; g_rb_stats.groups += stats.groups;
        g_rb_stats.planned_no += stats.planned_no; g_rb_stats.comp_cap += stats.comp_cap;
        g_rb_stats.tight += stats.tight;
        g_rb_stats.p1_band += stats.p1_band; g_rb_stats.p1_kswv += stats.p1_kswv;
        g_rb_stats.p1_guard += stats.p1_guard;
        g_rb_stats.t_plan_ns += (uint64_t)(stats.t_plan * 1e9);
        g_rb_stats.t_band_ns += (uint64_t)(stats.t_band * 1e9);
        stats = rescue_band_stats();
    }
    /* regid2parent_ entries are rewritten by commit() for every regid of the next batch. */
    recs_.clear();
    bands_.clear();
    p1_.clear();
    pending_ = -1;
}

/* ------------------------------------------------------------------------------------------ */
/* NEON 16-lane banded kernel                                                                  */
/* ------------------------------------------------------------------------------------------ */
/* Band coordinates: cell (row r, band k) of lane l is (i, d) = (r0_l + r, dlo_l + k) in parent
 * coordinates, query column j = r - k + o_l with o_l = r0_l - dlo_l. Diagonal predecessor = same
 * k previous row, E predecessor = k - 1 previous row, F predecessor = k + 1 same row (so k runs
 * W-1 .. 0). The query SoA is indexed by p = r - k + C (C = W - 1, lane-uniform):
 * A[p][l] = code(q_l[p - C + o_l]), so every cell reads one lane-uniform address (no gather).
 * Every lane of a group runs the group's width W: a lane's band is WIDENED to W (more diagonals
 * above its own), which rescue_band.h shows is exact, so no per-lane width mask is needed.
 *
 * Query codes: base 0..3; N -> 8; pad [len2, quanta) -> 0x40; nonexistent (j < 0 or
 * j >= quanta) -> 0xC0. 0x40 and 0xC0 index outside the 16-entry table, so vqtbl1q gives 0: the
 * pad scores 0 exactly like kswv's NEON_QPAD8 column (m11 = h00). Ref codes: base 0..3; N -> 4;
 * rows past the lane's range -> 0x80 (also the row-inactive flag). Score index = q ^ r:
 * 0 match (+1), 1..3 mismatch (-4), 4..15 N (-1) -- the kswv table at default scoring.
 * Nonexistent cells: j < 0 cells have only j < 0 (or row -1) predecessors and score 0, so they
 * stay H = E = F = 0 by induction; j >= quanta cells never feed a j < quanta cell (diagonal, E and
 * F all move to larger or equal j), so they are only excluded from the row max and the qe
 * snapshot (QL mask). The F chain is the 2-cell prefix-scan form: h0 (everything but F) does not
 * depend on F, so f_in(k-1) = max(h0(k) - 7, f_in(k) - 1) and f_in(k-2) = max(h0(k-1) - 7,
 * h0(k) - 8, f_in(k) - 2): one qsub + max per two cells on the loop-carried path (the same values
 * as the 1-cell recurrence).
 * qe = min j with H == gmax in row te == max k; lanes that improved snapshot their row's H at the
 * end of each run of consecutive improving rows (double-buffered H keeps row r-1 available). */
#if defined(__aarch64__)
namespace {

struct rb_work {
    std::vector<uint8_t> A, QL, REF, H, E, R, SNAP, ST, RL;
    alignas(16) uint16_t te[16];
    alignas(16) uint8_t gmax[16];
    /* Per-lane early-exit target (pass 1: the pass-0 score S; 0 for empty lanes). */
    alignas(16) uint8_t target[16];
    /* rb_dp_wave2's qe: band index of the lane's qe cell, k = kb_chunk * 255 + kb_off - 1
     * (rb_qe_scan); the other kernels leave it unset and snapshot the H row instead (SNAP). */
    alignas(16) uint8_t kb_chunk[16], kb_off[16];
    template <class V> static void fit(V &v, size_t n) { if (v.size() < n) v.resize(n); }
};

static inline uint64_t rb_mask16(uint8x16_t m)
{
    return vget_lane_u64(vreinterpret_u64_u8(vshrn_n_u16(vreinterpretq_u16_u8(m), 4)), 0);
}

/* In-register 16x16 byte transpose: r[l] byte t -> r[t] byte l. */
static inline void rb_transpose16(uint8x16_t r[16])
{
    for (int i = 0; i < 16; i += 2) {
        const uint8x16_t a = r[i], b = r[i + 1];
        r[i] = vtrn1q_u8(a, b); r[i + 1] = vtrn2q_u8(a, b);
    }
    for (int base = 0; base < 16; base += 4)
        for (int o = 0; o < 2; o++) {
            const uint16x8_t a = vreinterpretq_u16_u8(r[base + o]), b = vreinterpretq_u16_u8(r[base + o + 2]);
            r[base + o] = vreinterpretq_u8_u16(vtrn1q_u16(a, b));
            r[base + o + 2] = vreinterpretq_u8_u16(vtrn2q_u16(a, b));
        }
    for (int base = 0; base < 16; base += 8)
        for (int o = 0; o < 4; o++) {
            const uint32x4_t a = vreinterpretq_u32_u8(r[base + o]), b = vreinterpretq_u32_u8(r[base + o + 4]);
            r[base + o] = vreinterpretq_u8_u32(vtrn1q_u32(a, b));
            r[base + o + 4] = vreinterpretq_u8_u32(vtrn2q_u32(a, b));
        }
    for (int o = 0; o < 8; o++) {
        const uint64x2_t a = vreinterpretq_u64_u8(r[o]), b = vreinterpretq_u64_u8(r[o + 8]);
        r[o] = vreinterpretq_u8_u64(vtrn1q_u64(a, b));
        r[o + 8] = vreinterpretq_u8_u64(vtrn2q_u64(a, b));
    }
}

/* ST holds 16 lane-major rows of stride `stride` (a multiple of 16); write n * 16 interleaved
 * bytes (position-major, lane-minor) to dst. */
static inline void rb_transpose_to(const uint8_t *ST, int stride, int n, uint8_t *dst)
{
    uint8x16_t r[16];
    for (int b = 0; b < n; b += 16) {
        for (int l = 0; l < 16; l++) r[l] = vld1q_u8(ST + (size_t)l * stride + b);
        rb_transpose16(r);
        const int m = std::min(16, n - b);
        for (int t = 0; t < m; t++) vst1q_u8(dst + (size_t)(b + t) * 16, r[t]);
    }
}

/* Copy row `row`'s H into SNAP for the lanes in msk (live cells only; dead cells 0). */
static inline void rb_snapshot(rb_work &w, const uint8_t *Hrow, int row, int W, uint8x16_t msk)
{
    uint8_t *SNAP = w.SNAP.data();
    const uint8_t *qp = w.QL.data() + (size_t)row * 16;
    for (int k = W - 1; k >= 0; k--, qp += 16) {
        const uint8x16_t m = vandq_u8(msk, vld1q_u8(qp));
        vst1q_u8(SNAP + k * 16, vbslq_u8(m, vld1q_u8(Hrow + k * 16), vbicq_u8(vld1q_u8(SNAP + k * 16), msk)));
    }
}

/* omax = max over lanes of o_l, ominq = min over lanes of (o_l - quanta_l + 1): row r only computes
 * k in [max(0, r + ominq), min(W - 1, r + omax)], the union of the lanes' live cells
 * (0 <= j < quanta_l). Cells above that range are j < 0 in every lane: never written, so they stay 0
 * in both H buffers and in E, which is exactly their value. Cells below it are j >= quanta in every
 * lane: dead, and they only ever feed dead cells, so skipping them changes nothing live. Both bounds
 * advance by one per row, so every computed cell's predecessors were computed (or are those zeros). */
/* early: stop after the first row in which every lane's gmax has reached w.target. Pass 1 only:
 * there gmax can never exceed the target S, so once a lane reaches S no later row improves it, and
 * its te (last strict improvement) and qe snapshot (taken for that row at the latest when the loop
 * ends) are final -- exactly kswv's freeze at KSW_XSTOP. */
/* FScan selects the G-based cell (see RB_CELL1); both instantiations are byte-identical in every
 * value the caller reads (gmax, te, SNAP, R). */
template <bool FScan>
static long rb_dp_core(rb_work &w, int W, int NR, int omax, int ominq, int omaskq, bool early)
{
    long computed = 0;
    alignas(16) static const int8_t tblv[16] = {1, -4, -4, -4, -1, -1, -1, -1,
                                                -1, -1, -1, -1, -1, -1, -1, -1};
    const uint8x16_t tbl = vld1q_u8((const uint8_t *)tblv);
    const uint8x16_t v7 = vdupq_n_u8(7), v1 = vdupq_n_u8(1), v2 = vdupq_n_u8(2), v80 = vdupq_n_u8(0x80);
    uint8_t *Hc = w.H.data(), *Hp = Hc + (size_t)W * 16, *E = w.E.data();
    const uint8_t *A = w.A.data(), *QL = w.QL.data(), *REF = w.REF.data();
    uint8_t *Rout = w.R.data();
    memset(Hc, 0, (size_t)W * 32);
    memset(E, 0, (size_t)(W + 1) * 16);
    uint8x16_t gmax = vdupq_n_u8(0), pend = vdupq_n_u8(0);
    const uint8x16_t tgt = vld1q_u8(w.target);
    int rlast = NR - 1;
    uint16x8_t te_lo = vdupq_n_u16(0), te_hi = vdupq_n_u16(0);
    for (int r = 0; r < NR; r++) {
        const uint8x16_t rref = vld1q_u8(REF + (size_t)r * 16);
        uint8x16_t f = vdupq_n_u8(0), rmax = vdupq_n_u8(0);
        const int khi = std::min(W - 1, r + omax), klo = std::max(0, r + ominq);
        computed += khi >= klo ? khi - klo + 1 : 0;
        const size_t p0 = (size_t)(r + (W - 1 - khi)) * 16;                    // p = r + W-1-k
        const uint8_t *ap = A + p0, *qp = QL + p0;
        int k = khi;
        /* k >= kun: every real lane has j < quanta (or j < 0, where H is 0), so the live mask is
         * not needed; below kun some lane may sit on a dead (j >= quanta) cell. */
        const int kun = std::max(klo, r + omaskq);
#define RB_CELL2(MASK)                                                                          \
        {                                                                                      \
            const uint8x16_t qa = vld1q_u8(ap), qb = vld1q_u8(ap + 16);   /* cells k, k-1 */  \
            const uint8x16_t ea = vld1q_u8(E + k * 16), eb = vld1q_u8(E + (k - 1) * 16);      \
            const int8x16_t sa = vreinterpretq_s8_u8(vqtbl1q_u8(tbl, veorq_u8(qa, rref)));    \
            const int8x16_t sb = vreinterpretq_s8_u8(vqtbl1q_u8(tbl, veorq_u8(qb, rref)));    \
            const uint8x16_t h0a = vmaxq_u8(vsqaddq_u8(vld1q_u8(Hp + k * 16), sa), ea);       \
            const uint8x16_t h0b = vmaxq_u8(vsqaddq_u8(vld1q_u8(Hp + (k - 1) * 16), sb), eb); \
            const uint8x16_t h07a = vqsubq_u8(h0a, v7), h07b = vqsubq_u8(h0b, v7);            \
            const uint8x16_t fb = vmaxq_u8(h07a, vqsubq_u8(f, v1));        /* f_in(k-1) */    \
            const uint8x16_t ha = vmaxq_u8(h0a, f), hb = vmaxq_u8(h0b, fb);                    \
            f = vmaxq_u8(vmaxq_u8(h07b, vqsubq_u8(h07a, v1)), vqsubq_u8(f, v2)); /* f_in(k-2) */ \
            vst1q_u8(Hc + k * 16, ha);                                                         \
            vst1q_u8(Hc + (k - 1) * 16, hb);                                                   \
            /* FScan: row max over G (h0), E opened from G - 7 (see RB_CELL1). */            \
            const uint8x16_t ra = FScan ? h0a : ha, rb = FScan ? h0b : hb;                    \
            if (MASK) rmax = vmaxq_u8(rmax, vmaxq_u8(vandq_u8(ra, vld1q_u8(qp)), vandq_u8(rb, vld1q_u8(qp + 16)))); \
            else rmax = vmaxq_u8(rmax, vmaxq_u8(ra, rb));                                      \
            /* Cell k writes E slot k+1 and cell k-1 writes slot k, which cell k already read. */ \
            vst1q_u8(E + (k + 1) * 16, vmaxq_u8(FScan ? h07a : vqsubq_u8(ha, v7), vqsubq_u8(ea, v1))); \
            vst1q_u8(E + k * 16, vmaxq_u8(FScan ? h07b : vqsubq_u8(hb, v7), vqsubq_u8(eb, v1))); \
        }
/* FScan, the kswv G-based cell in band coordinates. Here the gap carried in a register is F, the
 * in-row gap (query advances, k - 1, same row), and E is the one carried through memory, the
 * vertical gap (reference advances: slot k + 1 read by row r + 1 at k + 1, same query column).
 * h0 = max(diagonal, E) is G, and F already opens from it, so F -> E (vertical opened straight
 * out of an in-row gap) is the transition E carries and F does not have the mirror of. FScan
 * opens E from G - 7 as well: one sat(G - 7) serves both gaps and E no longer waits on F.
 *  - Values: the zero-state DP keeps every alignment path except those with an in-row gap run
 *    immediately followed by a vertical one. Such a path's twin, the same two runs in the other
 *    order between the same two cells, scores the same (6 + a + 6 + b either way round) and F
 *    admits it. In the FULL DP that makes every H equal (kswv's point 1 at
 *    KSWV_NEON_U8_CELL_PAIR_FS). In a band the twin can leave it -- vertical first climbs to
 *    band index k + b, possibly past W - 1 -- so the band values need not match the original
 *    kernel's. What the band path relies on still holds: (i) any restricted DP never exceeds the
 *    full DP; (ii) the band holds the full-DP value at the end of every alignment scoring >= tau
 *    (pass 0) or of A* (pass 1) -- rescue_band.h -- because that alignment, if optimal with an
 *    in-row run next to a vertical one, has an equally scoring twin (swapping runs never merges
 *    or splits a run of an OPTIMAL alignment, so repeated swaps remove every such adjacency), and
 *    the twin lies in the band too: pass 0's containment covers EVERY alignment scoring >= tau,
 *    and pass 1's diagonal range [-Imax, Dmax] bounds every prefix of any alignment with A*'s
 *    I and D counts. gmax, te, the qe column (cells equal to gmax: exact by (i) + (ii)), the
 *    early exit and the rows R the merge reads are therefore those of the original kernel.
 *  - Row max over G, not H: within a row F only carries earlier cells' G decayed (F = 0 at khi,
 *    sat(x - 7) <= x), and those cells are live or j < 0 (G = 0) -- dead j >= quanta cells lie
 *    after every live one -- so the row max over H equals the row max over G, masked or not.
 *    SNAP still copies H (rb_snapshot reads Hc/Hp, which hold H). */
#define RB_CELL1(MASK)                                                                          \
        {                                                                                      \
            const uint8x16_t q = vld1q_u8(ap);                                                 \
            const uint8x16_t e = vld1q_u8(E + k * 16);                                         \
            const int8x16_t sc = vreinterpretq_s8_u8(vqtbl1q_u8(tbl, veorq_u8(q, rref)));     \
            const uint8x16_t h0 = vmaxq_u8(vsqaddq_u8(vld1q_u8(Hp + k * 16), sc), e);          \
            const uint8x16_t h = vmaxq_u8(h0, f);                                              \
            vst1q_u8(Hc + k * 16, h);                                                          \
            const uint8x16_t r_ = FScan ? h0 : h;                                              \
            rmax = vmaxq_u8(rmax, (MASK) ? vandq_u8(r_, vld1q_u8(qp)) : r_);                   \
            const uint8x16_t h07 = vqsubq_u8(h0, v7);                                          \
            vst1q_u8(E + (k + 1) * 16, vmaxq_u8(FScan ? h07 : vqsubq_u8(h, v7), vqsubq_u8(e, v1))); \
            f = vmaxq_u8(h07, vqsubq_u8(f, v1));                                               \
        }
        for (; k >= kun + 1; k -= 2, ap += 32, qp += 32) RB_CELL2(false)
        if (k == kun) { RB_CELL1(false) k--; ap += 16; qp += 16; }
        for (; k >= klo + 1; k -= 2, ap += 32, qp += 32) RB_CELL2(true)
        for (; k >= klo; k--, ap += 16, qp += 16) RB_CELL1(true)
#undef RB_CELL2
#undef RB_CELL1
        vst1q_u8(Rout + (size_t)r * 16, rmax);
        const uint8x16_t imp = vandq_u8(vcgtq_u8(rmax, gmax), vcltq_u8(rref, v80));
        const uint8x16_t flush = vbicq_u8(pend, imp);
        if (rb_mask16(flush)) rb_snapshot(w, Hp, r - 1, W, flush);
        gmax = vbslq_u8(imp, rmax, gmax);
        const uint16x8_t rv = vdupq_n_u16((uint16_t)r);
        te_lo = vbslq_u16(vreinterpretq_u16_u8(vzip1q_u8(imp, imp)), rv, te_lo);
        te_hi = vbslq_u16(vreinterpretq_u16_u8(vzip2q_u8(imp, imp)), rv, te_hi);
        pend = imp;
        std::swap(Hc, Hp);
        rlast = r;
        if (early && vminvq_u8(vcgeq_u8(gmax, tgt)) == 0xFF) break;
    }
    if (rb_mask16(pend)) rb_snapshot(w, Hp, rlast, W, pend);
    vst1q_u16(w.te, te_lo);
    vst1q_u16(w.te + 8, te_hi);
    vst1q_u8(w.gmax, gmax);
    return computed;
}

/* qe of the lanes in msk, read off row `row` (H in Hrow; gmax is those lanes' gmax, which this row
 * reached): the largest band index k in the row's computed range [klo, khi] whose LIVE cell holds
 * gmax -- the same cell rb_snapshot + lane_result find, without copying the row. Cells outside the
 * range are either never written (above khi: zero, and gmax > 0) or stale from earlier rows but dead
 * for every lane (below klo: j >= quanta, live mask 0), exactly the cells the snapshot masks out.
 * The index is kept as two bytes, k = chunk * 255 + off - 1: within a chunk of 255 indices the
 * match position is (match & kv) with kv = 1..255, so a plain max finds the largest matching k and
 * 0 means no match; the last matching chunk wins. A flushed lane always matches: its row max over
 * live cells was gmax in this row (FScan: the row max over H equals the row max over G, RB_CELL1).
 * Lanes outside msk keep their values. Cost: 2 loads + ~4 ops per cell-vector over the row's range,
 * against 3 loads + 4 ops + 1 store over all W slots for the row copy. */
static inline void rb_qe_scan(rb_work &w, const uint8_t *Hrow, int row, int W, int omax, int ominq,
                              uint8x16_t msk, uint8x16_t gmax)
{
    const int khi = std::min(W - 1, row + omax), klo = std::max(0, row + ominq);
    const uint8x16_t v1 = vdupq_n_u8(1), v2 = vdupq_n_u8(2);
    uint8x16_t chunk = vld1q_u8(w.kb_chunk), off = vld1q_u8(w.kb_off);
    const uint8_t *QL = w.QL.data();
    for (int k0 = klo, c = klo / 255; k0 <= khi; k0 = (c + 1) * 255, c++) {
        const int k1 = std::min(khi, c * 255 + 254);
        uint8x16_t m0 = vdupq_n_u8(0), m1 = vdupq_n_u8(0);
        uint8x16_t kv = vdupq_n_u8((uint8_t)(k0 - c * 255 + 1));
        const uint8_t *hp = Hrow + (size_t)k0 * 16, *qp = QL + (size_t)(row + W - 1 - k0) * 16;   // p = row + W-1-k
        int k = k0;
        for (; k + 1 <= k1; k += 2, hp += 32, qp -= 32, kv = vaddq_u8(kv, v2)) {
            const uint8x16_t e0 = vandq_u8(vceqq_u8(vld1q_u8(hp), gmax), vld1q_u8(qp));
            const uint8x16_t e1 = vandq_u8(vceqq_u8(vld1q_u8(hp + 16), gmax), vld1q_u8(qp - 16));
            m0 = vmaxq_u8(m0, vandq_u8(e0, kv));
            m1 = vmaxq_u8(m1, vandq_u8(e1, vaddq_u8(kv, v1)));
        }
        if (k <= k1) m0 = vmaxq_u8(m0, vandq_u8(vandq_u8(vceqq_u8(vld1q_u8(hp), gmax), vld1q_u8(qp)), kv));
        const uint8x16_t mm = vmaxq_u8(m0, m1);
        const uint8x16_t any = vtstq_u8(mm, mm);
        off = vbslq_u8(any, mm, off);
        chunk = vbslq_u8(any, vdupq_n_u8((uint8_t)c), chunk);
    }
    vst1q_u8(w.kb_chunk, vbslq_u8(msk, chunk, vld1q_u8(w.kb_chunk)));
    vst1q_u8(w.kb_off, vbslq_u8(msk, off, vld1q_u8(w.kb_off)));
}

/* The default kernel: the fused (FScan) cell of RB_CELL1 on TWO rows per step, (r, k) paired with
 * (r + 1, k + 1). The pair shares the query slot p = r + W - 1 - k (one A / QL load for both), and
 * row r + 1's cell takes its E predecessor (E_out of (r, k)) and its diagonal predecessor (H of
 * (r, k + 1), the previous step's ha) from registers, so per two cells the step loads A, E and Hp
 * once and stores H twice (row r's, for the qe scan) and E once (row r + 1's, for the next pair):
 * 22 ops + 3 loads + 3 stores against 24 + 6 + 4 for RB_CELL2, with two independent F chains.
 * H rotates through three buffers (rows r - 1, r, r + 1) so the qe of a run ending at row r - 1
 * (flushed at r) and one ending at r (flushed at r + 1) both still have their rows. Row r + 1's
 * range is row r's shifted by one, clipped to [0, W - 1]: the top cell (r, W - 1) has no partner
 * when row r already reaches W - 1 (its E_out would feed (r + 1, W), nonexistent), and the bottom
 * cell (r + 1, 0) has none when both rows start at 0 (its E predecessor is slot 0, which is never
 * written: 0, as in the one-row kernel). E slot k + 2, written by (r + 1, k + 1), was last read by
 * (r, k + 2) two steps earlier, so the in-place E buffer is safe. Every value each cell reads is the
 * one the one-row kernel reads (same buffers or the register copy of the same value), the row maxima
 * and the imp / te / flush bookkeeping are the same per row, the early exit is tested after each
 * row (row r + 1's updates are discarded when row r already qualifies), and qe comes from
 * rb_qe_scan on the same row and mask the snapshot would use. Pass 1 additionally flushes only
 * lanes whose gmax already equals the target S: gmax never exceeds S (rescue_band.h), so a run
 * ending below S is superseded by the run that reaches S, and a lane that never reaches it fails
 * run_pass1's guard and is rerun through kswv, as before. */
static long rb_dp_wave2(rb_work &w, int W, int NR, int omax, int ominq, int omaskq, bool early)
{
    long computed = 0;
    alignas(16) static const int8_t tblv[16] = {1, -4, -4, -4, -1, -1, -1, -1,
                                                -1, -1, -1, -1, -1, -1, -1, -1};
    const uint8x16_t tbl = vld1q_u8((const uint8_t *)tblv);
    const uint8x16_t v7 = vdupq_n_u8(7), v1 = vdupq_n_u8(1), v80 = vdupq_n_u8(0x80);
    uint8_t *Hp = w.H.data(), *Ha = Hp + (size_t)W * 16, *Hb = Hp + (size_t)W * 32, *E = w.E.data();
    const uint8_t *A = w.A.data(), *QL = w.QL.data(), *REF = w.REF.data();
    uint8_t *Rout = w.R.data();
    memset(Hp, 0, (size_t)W * 48);
    memset(E, 0, (size_t)(W + 2) * 16);
    uint8x16_t gmax = vdupq_n_u8(0), pend = vdupq_n_u8(0);
    const uint8x16_t tgt = vld1q_u8(w.target);
    int rlast = NR - 1;
    uint16x8_t te_lo = vdupq_n_u16(0), te_hi = vdupq_n_u16(0);
    /* Row-level update for `row` with row max rmax: R, imp, gmax, te, pend; returns the flush mask
     * (lanes whose improving run ended at row - 1; pass 1: only those already at their target). */
    auto row_update = [&](int row, uint8x16_t rmax, uint8x16_t rref) -> uint8x16_t {
        vst1q_u8(Rout + (size_t)row * 16, rmax);
        const uint8x16_t imp = vandq_u8(vcgtq_u8(rmax, gmax), vcltq_u8(rref, v80));
        uint8x16_t flush = vbicq_u8(pend, imp);
        if (early) flush = vandq_u8(flush, vceqq_u8(gmax, tgt));
        gmax = vbslq_u8(imp, rmax, gmax);
        const uint16x8_t rv = vdupq_n_u16((uint16_t)row);
        te_lo = vbslq_u16(vreinterpretq_u16_u8(vzip1q_u8(imp, imp)), rv, te_lo);
        te_hi = vbslq_u16(vreinterpretq_u16_u8(vzip2q_u8(imp, imp)), rv, te_hi);
        pend = imp;
        rlast = row;
        return flush;
    };
    /* One row on its own (the last row of an odd count, or a row whose pair has an empty range):
     * RB_CELL1 with Hp -> Ha and E in place. */
    auto single_row = [&](int r) -> bool {
        const uint8x16_t rref = vld1q_u8(REF + (size_t)r * 16);
        uint8x16_t f = vdupq_n_u8(0), rmax = vdupq_n_u8(0);
        const int khi = std::min(W - 1, r + omax), klo = std::max(0, r + ominq);
        computed += khi >= klo ? khi - klo + 1 : 0;
        const uint8_t *ap = A + (size_t)(r + (W - 1 - khi)) * 16, *qp = QL + (ap - A);
        const int kun = std::max(klo, r + omaskq);
        for (int k = khi; k >= klo; k--, ap += 16, qp += 16) {
            const uint8x16_t q = vld1q_u8(ap), e = vld1q_u8(E + k * 16);
            const int8x16_t sc = vreinterpretq_s8_u8(vqtbl1q_u8(tbl, veorq_u8(q, rref)));
            const uint8x16_t h0 = vmaxq_u8(vsqaddq_u8(vld1q_u8(Hp + k * 16), sc), e);
            vst1q_u8(Ha + k * 16, vmaxq_u8(h0, f));
            rmax = vmaxq_u8(rmax, k >= kun ? h0 : vandq_u8(h0, vld1q_u8(qp)));
            const uint8x16_t h07 = vqsubq_u8(h0, v7);
            vst1q_u8(E + (k + 1) * 16, vmaxq_u8(h07, vqsubq_u8(e, v1)));
            f = vmaxq_u8(h07, vqsubq_u8(f, v1));
        }
        const uint8x16_t flush = row_update(r, rmax, rref);
        if (rb_mask16(flush)) rb_qe_scan(w, Hp, r - 1, W, omax, ominq, flush, gmax);
        std::swap(Hp, Ha);
        return early && vminvq_u8(vcgeq_u8(gmax, tgt)) == 0xFF;
    };
    int r = 0;
    bool stop = false;
    for (; !stop && r + 1 < NR; r += 2) {
        const int khi_a = std::min(W - 1, r + omax), klo_a = std::max(0, r + ominq);
        const int khi_b = std::min(W - 1, r + 1 + omax), klo_b = std::max(0, r + 1 + ominq);
        if (khi_a < klo_a || khi_b < klo_b) { stop = single_row(r) || single_row(r + 1); continue; }
        const uint8x16_t rref_a = vld1q_u8(REF + (size_t)r * 16), rref_b = vld1q_u8(REF + (size_t)(r + 1) * 16);
        uint8x16_t fa = vdupq_n_u8(0), fb = vdupq_n_u8(0), rmax_a = vdupq_n_u8(0), rmax_b = vdupq_n_u8(0);
        uint8x16_t diag = vdupq_n_u8(0);   // H of (r, k + 1); zero above row r's range, never written
        const int kun_a = std::max(klo_a, r + omaskq), kun_b = std::max(klo_b, r + 1 + omaskq);
        const int kun = std::max(kun_a, kun_b - 1);   // both cells of a step unmasked iff k >= kun
        const uint8_t *ap = A + (size_t)(r + (W - 1 - khi_a)) * 16, *qp = QL + (ap - A);
        int k = khi_a;
        if (khi_a == W - 1) {   // (r, W - 1) alone: its partner (r + 1, W) does not exist
            const uint8x16_t q = vld1q_u8(ap), e = vld1q_u8(E + k * 16);
            const int8x16_t sc = vreinterpretq_s8_u8(vqtbl1q_u8(tbl, veorq_u8(q, rref_a)));
            const uint8x16_t h0 = vmaxq_u8(vsqaddq_u8(vld1q_u8(Hp + k * 16), sc), e);
            diag = vmaxq_u8(h0, fa);
            vst1q_u8(Ha + k * 16, diag);
            rmax_a = k >= kun_a ? h0 : vandq_u8(h0, vld1q_u8(qp));
            fa = vqsubq_u8(h0, v7);   // max(h07, sat(0 - 1))
            k--; ap += 16; qp += 16;
        }
#define RB_W2_STEP(MASK)                                                                        \
        {                                                                                      \
            const uint8x16_t q = vld1q_u8(ap), e = vld1q_u8(E + k * 16);                       \
            const int8x16_t sca = vreinterpretq_s8_u8(vqtbl1q_u8(tbl, veorq_u8(q, rref_a)));   \
            const uint8x16_t h0a = vmaxq_u8(vsqaddq_u8(vld1q_u8(Hp + k * 16), sca), e);        \
            const uint8x16_t h07a = vqsubq_u8(h0a, v7);                                        \
            const uint8x16_t ha = vmaxq_u8(h0a, fa);                                           \
            vst1q_u8(Ha + k * 16, ha);                                                         \
            fa = vmaxq_u8(h07a, vqsubq_u8(fa, v1));                                            \
            const uint8x16_t ea = vmaxq_u8(h07a, vqsubq_u8(e, v1));   /* E in of (r+1, k+1) */ \
            const int8x16_t scb = vreinterpretq_s8_u8(vqtbl1q_u8(tbl, veorq_u8(q, rref_b)));   \
            const uint8x16_t h0b = vmaxq_u8(vsqaddq_u8(diag, scb), ea);                        \
            const uint8x16_t h07b = vqsubq_u8(h0b, v7);                                        \
            vst1q_u8(Hb + (k + 1) * 16, vmaxq_u8(h0b, fb));                                    \
            fb = vmaxq_u8(h07b, vqsubq_u8(fb, v1));                                            \
            vst1q_u8(E + (k + 2) * 16, vmaxq_u8(h07b, vqsubq_u8(ea, v1)));                     \
            diag = ha;                                                                         \
            if (MASK) {                                                                        \
                const uint8x16_t m = vld1q_u8(qp);                                             \
                rmax_a = vmaxq_u8(rmax_a, vandq_u8(h0a, m));                                   \
                rmax_b = vmaxq_u8(rmax_b, vandq_u8(h0b, m));                                   \
            } else {                                                                           \
                rmax_a = vmaxq_u8(rmax_a, h0a);                                                \
                rmax_b = vmaxq_u8(rmax_b, h0b);                                                \
            }                                                                                  \
        }
        for (; k >= kun; k--, ap += 16, qp += 16) RB_W2_STEP(false)
        for (; k >= klo_a; k--, ap += 16, qp += 16) RB_W2_STEP(true)
#undef RB_W2_STEP
        if (klo_b == klo_a) {   // (r + 1, 0) alone: its partner (r, -1) does not exist; E in = slot 0 = 0
            const uint8x16_t q = vld1q_u8(ap);
            const int8x16_t scb = vreinterpretq_s8_u8(vqtbl1q_u8(tbl, veorq_u8(q, rref_b)));
            const uint8x16_t h0b = vsqaddq_u8(diag, scb);
            vst1q_u8(Hb + klo_a * 16, vmaxq_u8(h0b, fb));
            vst1q_u8(E + (klo_a + 1) * 16, vqsubq_u8(h0b, v7));
            rmax_b = vmaxq_u8(rmax_b, vandq_u8(h0b, vld1q_u8(qp)));
        }
        computed += khi_a - klo_a + 1;
        /* Row r, then row r + 1; one cross-domain test covers both flush masks (and, in pass 1, one
         * the early exit, with row r re-tested only when the pair qualifies). */
        const uint8x16_t flush_a = row_update(r, rmax_a, rref_a);
        const uint8x16_t gmax_a = gmax;
        if (early && vminvq_u8(vcgeq_u8(gmax_a, tgt)) == 0xFF) {
            if (rb_mask16(flush_a)) rb_qe_scan(w, Hp, r - 1, W, omax, ominq, flush_a, gmax_a);
            std::swap(Hp, Ha);
            stop = true;
            break;
        }
        computed += khi_b - klo_b + 1;
        const uint8x16_t flush_b = row_update(r + 1, rmax_b, rref_b);
        if (rb_mask16(vorrq_u8(flush_a, flush_b))) {
            if (rb_mask16(flush_a)) rb_qe_scan(w, Hp, r - 1, W, omax, ominq, flush_a, gmax_a);
            if (rb_mask16(flush_b)) rb_qe_scan(w, Ha, r, W, omax, ominq, flush_b, gmax);
        }
        uint8_t *t = Hp; Hp = Hb; Hb = Ha; Ha = t;
        stop = early && vminvq_u8(vcgeq_u8(gmax, tgt)) == 0xFF;
    }
    if (!stop && r < NR) single_row(r);
    uint8x16_t fin = pend;
    if (early) fin = vandq_u8(fin, vceqq_u8(gmax, tgt));
    if (rb_mask16(fin)) rb_qe_scan(w, Hp, rlast, W, omax, ominq, fin, gmax);
    vst1q_u16(w.te, te_lo);
    vst1q_u16(w.te + 8, te_hi);
    vst1q_u8(w.gmax, gmax);
    return computed;
}

}  // namespace
#endif

/* ------------------------------------------------------------------------------------------ */
/* Group driver + merge                                                                        */
/* ------------------------------------------------------------------------------------------ */

#ifdef RB_PROFILE
struct rb_prof_t {
    double build = 0, dp = 0, post = 0, sort = 0;
    long whist[9] = {0}, rhist[9] = {0}, groups = 0; double cells = 0, lanes = 0;
    ~rb_prof_t()
    {
        fprintf(stderr, "[RB_PROFILE] sort %.3f build %.3f dp %.3f post %.3f s; groups %ld lanes %.0f padded cells %.3g\n",
                sort, build, dp, post, groups, lanes, cells);
        fprintf(stderr, "  W hist (<=8,16,32,64,128,256,512,1024,>):");
        for (long x : whist) fprintf(stderr, " %ld", x);
        fprintf(stderr, "\n  NR hist (<=64,128,256,512,1024,2048,4096,8192,>):");
        for (long x : rhist) fprintf(stderr, " %ld", x);
        fprintf(stderr, "\n");
    }
};
static rb_prof_t g_rb_prof;
static inline int rb_log_bucket(int v, int base)
{
    int b = 0;
    while (b < 8 && v > base << b) b++;
    return b;
}
#define RB_T(x) const auto x = std::chrono::steady_clock::now()
#define RB_ACC(field, a, b) g_rb_prof.field += std::chrono::duration<double>(b - a).count()
#else
#define RB_T(x)
#define RB_ACC(field, a, b)
#endif

#if defined(__aarch64__)
/* rescue_band_score2 for a CONTIGUOUS row-max array R[0, n) whose slack R[n, n + 16] is zero
 * (so row n - 1 sees a 0 successor, as it does in the merged hull view). The lagged zeroing and the
 * >= minsc test run 16 rows per vector; only qualifying rows reach the scalar b[] emulation. */
static void rb_score2_vec(const uint8_t *R, int n, int row0, int S, int te, int minsc, int *score2, int *te2)
{
    const int low = te - S, high = te + S;
    int s2 = -1, t2 = -1, bs = -1, bp = -2;
    const uint8x16_t ms = vdupq_n_u8((uint8_t)minsc);
    alignas(16) uint8_t vb[16];
    for (int k0 = 0; k0 < n; k0 += 16) {
        const uint8x16_t cur = vld1q_u8(R + k0), nxt = vld1q_u8(R + k0 + 1);
        const uint8x16_t v = vbicq_u8(cur, vcgtq_u8(nxt, cur));   // rising row -> 0
        uint64_t m = rb_mask16(vcgeq_u8(v, ms));
        if (!m) continue;
        vst1q_u8(vb, v);
        do {
            const int t = __builtin_ctzll(m) >> 2;
            m &= ~(0xFull << (t * 4));
            const int val = vb[t], i = row0 + k0 + t;
            if (bp + 1 != i) {
                if (bp >= 0 && (bp < low || bp > high) && bs > s2) { s2 = bs; t2 = bp; }
                bs = val; bp = i;
            } else if (bs < val) {
                bs = val; bp = i;
            }
        } while (m);
    }
    if (bp >= 0 && (bp < low || bp > high) && bs > s2) { s2 = bs; t2 = bp; }
    *score2 = s2;
    *te2 = t2;
}
#endif

static inline int rb_width_bucket(int w)
{
    return w <= 16 ? (w + 3) >> 2 : w <= 64 ? 4 + ((w - 16 + 7) >> 3) : w <= 128 ? 10 + ((w - 64 + 15) >> 4)
                                                                                 : 14 + ((w - 128 + 31) >> 5);
}

/* Run jobs_ in groups of 16 lanes. Pass 0 (pass1 == false) merges each lane into its parent
 * (ps_[parent]); pass 1 only needs each lane's gmax / first row / first column, stored in
 * p1res_[parent], and skips the row-max transpose. Stage counters cover pass 0 only. */
void RescueBandBatch::run_jobs(bool pass1)
{
#if defined(__aarch64__)
    static thread_local rb_work w;
    const int n = (int)jobs_.size();
    if (n == 0) return;
    RB_T(ts0);
    order_.resize(n);
    for (int i = 0; i < n; i++) {
        order_[i] = i;
        /* Group by width class, then rows: lanes of one group then share their live-cell clipping
         * (rb_dp_core), which measured better than grouping by the clipped width. */
        jobs_[i].key = (rb_width_bucket(jobs_[i].w) << 16) | std::min(jobs_[i].nrows, 65535);
    }
    std::sort(order_.begin(), order_.end(), [this](int a, int b) { return jobs_[a].key < jobs_[b].key; });
    RB_T(ts1);
    RB_ACC(sort, ts0, ts1);
    const bool st = rb_stats_on() && !pass1;
    for (int g = 0; g < n; g += 16) {
        const int nl = std::min(16, n - g);
        const job *L[16];
        int W = 1, NR = 1;
        for (int l = 0; l < nl; l++) {
            L[l] = &jobs_[order_[g + l]];
            W = std::max(W, L[l]->w);
            NR = std::max(NR, L[l]->nrows);
        }
        if (st) {
            stats.groups++;
            for (int l = 0; l < nl; l++) stats.cells_req += (uint64_t)std::min(L[l]->w, L[l]->quanta) * L[l]->nrows;
        }
        /* ---- SoA build ---- */
        RB_T(tb0);
#ifdef RB_PROFILE
        g_rb_prof.groups++; g_rb_prof.lanes += nl;
        g_rb_prof.whist[rb_log_bucket(W, 8)]++; g_rb_prof.rhist[rb_log_bucket(NR, 64)]++;
#endif
        const int C = W - 1, P = NR + W - 1;
        const int Pp = (P + 15) & ~15, NRp = (NR + 15) & ~15;
        const int stride = std::max(Pp, NRp);
        rb_work::fit(w.A, (size_t)Pp * 16); rb_work::fit(w.QL, (size_t)Pp * 16);
        rb_work::fit(w.REF, (size_t)NRp * 16); rb_work::fit(w.H, (size_t)W * 48);
        rb_work::fit(w.E, (size_t)(W + 2) * 16); rb_work::fit(w.R, (size_t)NRp * 16);
        rb_work::fit(w.SNAP, (size_t)W * 16); rb_work::fit(w.ST, (size_t)stride * 16);
        /* Band shift: a lane narrower than the group is widened to W anyway; put the spare
         * diagonals BELOW its band (dlo - delta) rather than above whenever that moves its query
         * offset o = r0 - dlo up toward the group's largest, so the lanes' live column ranges
         * line up and the per-row union the kernel computes (omax - ominq wide) shrinks: -9.5% of
         * pass-0 cells on real jobs. Exact: the lane's cells are a superset of its band's, which
         * rescue_band.h shows keeps every property the merge relies on, and the pass-1 band
         * [-Imax, Dmax] only gets more insertion diagonals. oe[l] replaces o everywhere below. */
        int oe[16];
        {
            int om = INT_MIN;
            for (int l = 0; l < nl; l++) om = std::max(om, L[l]->r0 - L[l]->dlo);
            for (int l = 0; l < nl; l++) {
                const int o = L[l]->r0 - L[l]->dlo;
                oe[l] = rb_shift_on() ? o + std::min(om - o, W - L[l]->w) : o;
            }
        }
        uint8_t *ST = w.ST.data();
        for (int l = 0; l < 16; l++) {
            uint8_t *row = ST + (size_t)l * stride;
            memset(row, 0xC0, Pp);
            if (l >= nl) continue;
            const job &J = *L[l];
            const int pq0 = C - oe[l];   // p where j = 0
            const int a = std::max(0, pq0), b = std::min(P, pq0 + J.len2);
            const uint8_t *q = J.qry - pq0;
            for (int p = a; p < b; p++) { const uint8_t c = q[p]; row[p] = c < 4 ? c : 8; }
            const int c0 = std::max(0, pq0 + J.len2), c1 = std::min(P, pq0 + J.quanta);
            if (c1 > c0) memset(row + c0, 0x40, c1 - c0);
        }
        rb_transpose_to(ST, stride, P, w.A.data());
        {
            const uint8x16_t vc0 = vdupq_n_u8(0xC0);
            const uint8_t *Ap = w.A.data();
            uint8_t *Qp = w.QL.data();
            for (int p = 0; p < P; p++) vst1q_u8(Qp + p * 16, vcltq_u8(vld1q_u8(Ap + p * 16), vc0));
        }
        for (int l = 0; l < 16; l++) {
            uint8_t *row = ST + (size_t)l * stride;
            int rows = 0;
            if (l < nl) {
                const job &J = *L[l];
                rows = J.nrows;
                const uint8_t *src = J.ref + J.r0;
                for (int r = 0; r < rows; r++) row[r] = src[r] < 4 ? src[r] : 4;
            }
            if (NRp > rows) memset(row + rows, 0x80, NRp - rows);
        }
        rb_transpose_to(ST, stride, NR, w.REF.data());
        /* Prefetch the next group's reference windows (staged by _pre, possibly out of L1/L2). */
        for (int l = 0; l < 16 && g + 16 + l < n; l++) {
            const job &J = jobs_[order_[g + 16 + l]];
            for (int r = 0; r < J.nrows; r += 64) __builtin_prefetch(J.ref + J.r0 + r, 0, 1);
            __builtin_prefetch(J.qry, 0, 1);
            __builtin_prefetch(J.qry + 64, 0, 1);
        }
        /* ---- DP ---- */
        RB_T(tb1);
        int omax = INT_MIN, ominq = INT_MAX, omaskq = INT_MIN;
        for (int l = 0; l < nl; l++) {
            const int o = oe[l];
            omax = std::max(omax, o);
            ominq = std::min(ominq, o - L[l]->quanta + 1);
            omaskq = std::max(omaskq, o - L[l]->quanta + 1);
        }
        for (int l = 0; l < 16; l++) w.target[l] = pass1 && l < nl ? (uint8_t)L[l]->target : 0;
        const int kern = rb_kernel();
        const long computed = kern >= 2 ? rb_dp_wave2(w, W, NR, omax, ominq, omaskq, pass1)
                            : kern == 1 ? rb_dp_core<true>(w, W, NR, omax, ominq, omaskq, pass1)
                                        : rb_dp_core<false>(w, W, NR, omax, ominq, omaskq, pass1);
        if (st) stats.cells_pad += (uint64_t)16 * computed;
#ifdef RB_PROFILE
        g_rb_prof.cells += 16.0 * computed;
#endif
        RB_T(tb2);
        RB_ACC(build, tb0, tb1);
        RB_ACC(dp, tb1, tb2);
        /* Lane l's gmax, the first row reaching it and the first column holding it there. */
        auto lane_result = [&](int l) {
            const job &J = *L[l];
            lane_res x{w.gmax[l], -1, 0};
            if (x.g > 0) {
                const int rt = w.te[l];
                x.te = J.r0 + rt;
                if (kern >= 2) x.qe = rt - ((int)w.kb_chunk[l] * 255 + (int)w.kb_off[l] - 1) + oe[l];
                else
                    for (int k = W - 1; k >= 0; k--)   // max k == min j
                        if (w.SNAP[k * 16 + l] == x.g) { x.qe = rt - k + oe[l]; break; }
            }
            return x;
        };
        if (pass1) {
            for (int l = 0; l < nl; l++) p1res_[L[l]->parent] = lane_result(l);
            continue;
        }
        /* ---- merge into the parents ---- */
        /* R back to lane-major (16x16 transposes), each lane's rows past its range zeroed, so
         * per-lane work below is contiguous vector code. */
        const int rls = NRp + 32;
        rb_work::fit(w.RL, (size_t)rls * 16);
        {
            uint8x16_t rr[16];
            uint8_t *RL = w.RL.data();
            for (int b = 0; b < NR; b += 16) {
                for (int t = 0; t < 16; t++) rr[t] = vld1q_u8(w.R.data() + (size_t)(b + t) * 16);
                rb_transpose16(rr);
                for (int l = 0; l < nl; l++) vst1q_u8(RL + (size_t)l * rls + b, rr[l]);
            }
            for (int l = 0; l < nl; l++) memset(RL + (size_t)l * rls + L[l]->nrows, 0, rls - L[l]->nrows);
        }
        for (int l = 0; l < nl; l++) {
            const job &J = *L[l];
            pstate &p = ps_[J.parent];
            const lane_res x = lane_result(l);
            const int g8 = x.g, te = x.te, qe = x.qe;
            const uint8_t *Rl = w.RL.data() + (size_t)l * rls;
            if (p.rbuf < 0) {   // the parent's only band: its outputs are the parent's
                p.S = g8; p.te = te; p.qe = qe;
                const int minsc = recs_[p.rec].minsc;
                if (g8 >= minsc) rb_score2_vec(Rl, J.nrows, J.r0, g8, te, minsc, &p.s2, &p.te2);
                else p.s2 = p.te2 = -1;
            } else {
                /* rbuf has >= 32 bytes of zeroed slack past the hull, and Rl is zero past nrows,
                 * so whole 16-row chunks can be max-merged. */
                uint8_t *buf = rpool_.data() + p.rbuf + J.r0;
                for (int r = 0; r < J.nrows; r += 16) vst1q_u8(buf + r, vmaxq_u8(vld1q_u8(buf + r), vld1q_u8(Rl + r)));
                if (g8 > p.S) { p.S = g8; p.te = te; p.qe = qe; }
                else if (g8 == p.S && g8 > 0) {
                    if (te < p.te) { p.te = te; p.qe = qe; }
                    else if (te == p.te && qe < p.qe) p.qe = qe;
                }
            }
        }
        RB_T(tb3);
        RB_ACC(post, tb2, tb3);
    }
#endif
}

void RescueBandBatch::finish_parent(pstate &p)
{
    if (p.rbuf >= 0) {
        const int minsc = recs_[p.rec].minsc;
        if (p.S >= minsc) {
#if defined(__aarch64__)
            rb_score2_vec(rpool_.data() + p.rbuf, p.L, 0, p.S, p.te, minsc, &p.s2, &p.te2);
#else
            rescue_band_score2(rpool_.data() + p.rbuf, 1, p.L, 0, p.S, p.te, minsc, &p.s2, &p.te2);
#endif
        } else
            p.s2 = p.te2 = -1;
    }
}

void RescueBandBatch::run_pass0(const SeqPair *pairs, int nb, const uint8_t *seqBufRef,
                                const uint8_t *seqBufQer, kswr_t *aln, Ikswv *kswv)
{
    if (nb <= 0) return;
    const bool st = rb_stats_on();
    const auto t0 = st ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point();
    ps_.clear();
    jobs_.clear();
    size_t rused = 0;
    auto add_jobs = [&](int pi, const SeqPair &sp, const rb_band *bb, int nbb, int round) {
        pstate &p = ps_[pi];
        p.S = 0; p.te = -1; p.qe = 0; p.s2 = p.te2 = -1; p.round = round;
        p.rbuf = -1;
        if (nbb > 1) {   // merged row-max buffer: the hull rows + 32 zero bytes of slack
            p.rbuf = (int32_t)rused;
            rused += (size_t)p.L + 32;
            if (rpool_.size() < rused) rpool_.resize(rused + 65536);
            memset(rpool_.data() + p.rbuf, 0, (size_t)p.L + 32);
        }
        const int quanta = ((sp.len2 + 15) / 16) * 16;
        for (int k = 0; k < nbb; k++) {
            const rb_band &b = bb[k];
            job J;
            J.ref = seqBufRef + sp.idr; J.qry = seqBufQer + sp.idq;
            J.len2 = sp.len2; J.quanta = quanta;
            J.r0 = b.r0; J.nrows = b.r1 - b.r0 + 1; J.dlo = b.dlo; J.w = b.dhi - b.dlo + 1;
            J.parent = pi; J.key = 0;
            jobs_.push_back(J);
        }
    };
    auto write = [&](const pstate &p) {
        kswr_t &a = aln[p.regid];
        const int shift = 4;   // kswv 8-bit bias at default scoring: 256 - min(1, -4, -1)
        a.score = p.S + shift < 255 ? p.S : 255;
        a.te = p.te;
        a.qe = p.qe;
        if (a.score == 255) { a.score2 = a.te2 = -1; }
        else { a.score2 = p.s2; a.te2 = p.te2; }
    };
    /* ---- round 1: every banded parent's components at T1 ---- */
    ps_.resize(nb);
    for (int i = 0; i < nb; i++) {
        const SeqPair &sp = pairs[i];
        const int ri = regid2parent_[sp.regid];
        const parent_rec &rec = recs_[ri];
        pstate &p = ps_[i];
        p.regid = sp.regid; p.L = sp.len1; p.rec = ri;
        add_jobs(i, sp, &bands_[rec.b0], rec.nb, 1);
        if (st) {
            stats.parents++; stats.bands1 += rec.nb; stats.single += rec.nb == 1;
            stats.hull_cells += (uint64_t)sp.len1 * (((sp.len2 + 15) / 16) * 16);
        }
    }
    run_jobs(false);
    /* ---- termination test; round 2 for the rest ---- */
    int n_r2 = 0, n_kswv = 0;
    jobs_.clear();
    rused = 0;
    static thread_local std::vector<int32_t> kswv_list, r2buf;   // grow-only scratch
    kswv_list.clear();
    r2buf.clear();
    for (int i = 0; i < nb; i++) {
        pstate &p = ps_[i];
        finish_parent(p);
        const parent_rec &rec = recs_[p.rec];
        /* Exact iff T1 == minsc, or S >= T1 (S / te / qe exact) and either score2 >= T1 (the
         * corrected rule) or the hull lies inside [te - S, te + S] (then score2 is -1, and the merged
         * value is -1 too: its anchors are rows with a true R >= minsc, all inside the hull). */
        const bool done = rec.T1 == rec.minsc
            || (p.S >= rec.T1 && (p.s2 >= rec.T1 || (p.te - p.S <= 0 && p.te + p.S >= p.L - 1)));
        if (done) { write(p); continue; }
        r2buf.push_back(i);
    }
    for (int32_t i : r2buf) {
        pstate &p = ps_[i];
        const parent_rec &rec = recs_[p.rec];
        if (!rb_r2_banded() || rec.nc == 0) { kswv_list.push_back(i); n_kswv++; continue; }
        /* Components at minsc; skip those confined to (te - S, te) whose bound is < S when the
         * round-1 S (>= T1) is already exact (optsim v4). */
        const int Sh = p.S >= rec.T1 ? p.S : 0, teh = p.te;
        const rb_band *cb = &bands_[rec.c0];
        int keep = 0;
        static thread_local std::vector<rb_band> kept;
        kept.clear();
        for (int k = 0; k < rec.nc; k++) {
            const rb_band &b = cb[k];
            if (Sh >= rec.minsc && b.ub < Sh && b.r0 > teh - Sh && b.r1 < teh) continue;
            kept.push_back(b); keep++;
        }
        add_jobs(i, pairs[i], kept.data(), keep, 2);
        n_r2++;
        if (st) stats.bands2 += keep;
    }
    if (n_r2) {
        run_jobs(false);
        for (int32_t i : r2buf) {
            pstate &p = ps_[i];
            if (p.round != 2) continue;
            finish_parent(p);
            write(p);
        }
    }
    if (st) { stats.r2_band += n_r2; stats.r2_kswv += n_kswv; }
    if (n_kswv) {
        if ((int)spscratch_.size() < n_kswv + 64) spscratch_.resize(n_kswv + 64);
        int m = 0;
        for (int32_t i : kswv_list) spscratch_[m++] = pairs[i];
        kswv->getScores8(spscratch_.data(), (uint8_t *)seqBufRef, (uint8_t *)seqBufQer, aln, m, 1, 0);
    }
    if (st) stats.t_band += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

/* ------------------------------------------------------------------------------------------ */
/* Pass 1 (start recovery)                                                                     */
/* ------------------------------------------------------------------------------------------ */

bool RescueBandBatch::take_pass1(const SeqPair &sp, const kswr_t &r, bool banded_parent)
{
#if defined(__aarch64__)
    if (rb_p1_mode() < (banded_parent ? 1 : 2)) return false;
    const int S = r.score, te = r.te, qe = r.qe;
    /* 8-bit and unsaturated (kswv's 255 sentinel is S + shift >= 255, shift 4), and a real end. */
    if (S <= 0 || S + 4 >= 255 || te < 0 || qe < 0 || sp.len2 != qe + 1) { stats.p1_kswv++; return false; }
    /* The band of rescue_band.h: A* has D <= dmax deleted and I <= imax inserted bases, and spans at
     * most qe + 1 + dmax rows of the te + 1 reversed ones. Diagonals past the last row (d > nrows - 1)
     * or past the query (d < -(quanta - 1)) hold no cell, so the band is clipped to them. */
    const int quanta = ((qe + 1 + 15) / 16) * 16;
    const int dall = std::max(0, qe - S - 5);
    const int nrows = std::min(te + 1, qe + 1 + dall);
    const int dmax = std::min(dall, nrows - 1), imax = std::min(dall / 2, quanta - 1);
    /* Both kswv and the band stop at the first row reaching S (the band via rb_dp_core's early
     * exit), so they run about the same rows: compare the per-row cells. */
    const int w = imax + dmax + 1;
    if ((long)(std::min(w, quanta) + RB_OVH) * 100 >= (long)quanta * rb_p1_cost_pct()) { stats.p1_kswv++; return false; }
    p1_.push_back(p1job{sp, nrows, imax, dmax});
    return true;
#else
    (void)sp; (void)r; (void)banded_parent;
    return false;
#endif
}

void RescueBandBatch::run_pass1(const uint8_t *seqBufRef, const uint8_t *seqBufQer, kswr_t *aln, Ikswv *kswv)
{
    const int n = (int)p1_.size();
    if (n == 0) return;
    jobs_.clear();
    for (int i = 0; i < n; i++) {
        const p1job &P = p1_[i];
        job J;
        J.ref = seqBufRef + P.sp.idr; J.qry = seqBufQer + P.sp.idq;
        J.len2 = P.sp.len2; J.quanta = ((P.sp.len2 + 15) / 16) * 16;
        J.r0 = 0; J.nrows = P.nrows; J.dlo = -P.imax; J.w = P.imax + P.dmax + 1;
        J.parent = i; J.key = 0; J.target = P.sp.h0 & 0xffff;   // KSW_XSTOP | S
        jobs_.push_back(J);
    }
    if ((int)p1res_.size() < n) p1res_.resize(n);
    run_jobs(true);
    /* kswv writes tb / qb only when its phase-1 score equals the pass-0 score; for the band that
     * always holds (rescue_band.h), so a lane that misses S is a broken invariant: rerun it
     * through kswv rather than write anything the kernel did not prove. */
    if ((int)spscratch_.size() < n + 64) spscratch_.resize(n + 64);   // kswv reads whole lane groups
    int m = 0;
    for (int i = 0; i < n; i++) {
        const lane_res &x = p1res_[i];
        kswr_t &a = aln[p1_[i].sp.regid];
        if (x.g == a.score) {
            a.tb = a.te - x.te;
            a.qb = a.qe - x.qe;
            continue;
        }
        spscratch_[m++] = p1_[i].sp;
    }
    stats.p1_band += n - m;
    stats.p1_guard += m;
    if (m) kswv->getScores8(spscratch_.data(), (uint8_t *)seqBufRef, (uint8_t *)seqBufQer, aln, m, 1, 1);
    p1_.clear();
}
