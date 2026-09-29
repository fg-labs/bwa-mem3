/* Exact mate-rescue window pruning.
 *
 * A mate-rescue Smith-Waterman job (mate `q` against reference window `ref`) contributes only
 * score/te/qe (pass 0), tb/qb (pass 1) and score2 (-> csub) to the output, and nothing at all
 * when score < minsc (the post gate fails). This filter proves, from exact 5-mer matches, which
 * part of the window can hold any row whose value matters, so the DP can run on a sub-window
 * (or not at all) with byte-identical results.
 *
 * Lemma (a = 1, b = 4, gap open 6 / extend 1; windows and mates without N -- an N breaks the
 * 2-bit 5-mer encoding, so any N makes the filter return RESCUE_PRUNE_FULL -- and the kernel's
 * query-pad columns, which score <= +1, treated as matches): a local alignment scores at most
 *     4 + (#exact 5-mer hits on its cells) - sum over gaps of (2 + |diagonal change|),
 * because a run of r matches holds max(0, r - 4) 5-mers and every mismatch costs 4. Hence an
 * alignment whose diagonals span the interval [d1, d2] scores at most
 *     5 + sum_{d in [d1, d2]} (cnt_d - 1),
 * with cnt_d the number of 5-mer hits on diagonal d = i - j.
 *
 * Decisions:
 *   RESCUE_PRUNE_B1   no diagonal interval reaches minsc: score < minsc is proven, so the job's
 *                     outputs are dead. The caller does not enqueue it and takes the ordinary
 *                     failing-rescue path instead (MATESW_GAR_PROVEN_FAIL in bwamem_pair.cpp).
 *   RESCUE_PRUNE_B2   rows [*hb, *he] contain every row whose max can reach minsc. A zero-state
 *                     DP started at *hb reproduces all consumed outputs: every alignment scoring
 *                     >= minsc starts at or after the first 5-run of its component (the prefix
 *                     before it nets <= 0), and the tail after its last hit nets <= -(2 + D) per
 *                     gap, which bounds the extent to (last hit diagonal) + quanta - 1 +
 *                     (ub - minsc - 2).
 *   RESCUE_PRUNE_FULL no pruning (N present, lengths too short or too long, minsc < 5, or too many
 *                     hits to be worth it).
 *
 * The bound and extents were validated against the kswv kernel's own per-row maxima on
 * 31.6 M wgs/wes rescue jobs (zero violations) and by whole-run output identity.
 *
 * Overview and gates: docs/src/developer-guide/rescue-pruning.md. */
#ifndef BWA_MEM3_RESCUE_PRUNE_H
#define BWA_MEM3_RESCUE_PRUNE_H

#include <stdint.h>
#include <string.h>
#include <algorithm>
#include <vector>
#include "kswv_quantum.h"
#include "rescue_prune_neon.h"

enum { RESCUE_PRUNE_FULL = 0, RESCUE_PRUNE_B1 = 1, RESCUE_PRUNE_B2 = 2 };

/* The scorings the lemma above is derived for. The filter takes no scoring and would return
 * false B1 "proofs" under any other (b = 2, say), so every caller must check this first. */
static inline bool rescue_prune_scoring_ok(int a, int b, int o_del, int e_del, int o_ins, int e_ins)
{
    return a == 1 && b == 4 && o_del == 6 && e_del == 1 && o_ins == 6 && e_ins == 1;
}

/* Per-thread scratch for the scalar path. Fixed capacity (windows are capped at 30000 rows and
 * queries at 1024 bases), so nothing is allocated per call. The query tables are cached: they
 * depend only on the oriented mate, which repeats across the anchors rescued with it. */
struct rescue_prune_scratch {
    static const int QCAP = 1024, DCAP = 30000 + 1024 + 2;
    int16_t head[1024], nxt[QCAP];
    uint16_t qcnt[1024];
    uint8_t qcache[QCAP];
    int qlen_c = -1;
    uint16_t cnt[DCAP];
    int16_t minrow[DCAP];
    int32_t fwd[DCAP], bwd[DCAP];
};

/* Scalar implementation (off aarch64, a threshold other than 19, or beyond the NEON filter's
 * capacity or int16 range). Same decisions as the NEON one wherever both run. */
static inline int rescue_prune_window_scalar(const uint8_t *ref, int len1, const uint8_t *q, int len2,
                                             int minsc, int max_hits, rescue_prune_scratch &s,
                                             int *hb, int *he)
{
    uint8_t orv = 0;
    for (int i = 0; i < len1; i++) orv |= ref[i];
    for (int j = 0; j < len2; j++) orv |= q[j];
    if (orv & 0xFC) return RESCUE_PRUNE_FULL;  // N present: its score (-1) breaks the lemma's accounting

    if (len2 != s.qlen_c || memcmp(q, s.qcache, (size_t)len2) != 0) {
        memset(s.head, 0xFF, sizeof s.head);
        memset(s.qcnt, 0, sizeof s.qcnt);
        int c = (q[0] << 6) | (q[1] << 4) | (q[2] << 2) | q[3];
        for (int j = 4; j < len2; j++) {
            c = ((c << 2) | q[j]) & 1023;
            s.nxt[j] = s.head[c];
            s.head[c] = (int16_t)j;
            s.qcnt[c]++;
        }
        memcpy(s.qcache, q, (size_t)len2);
        s.qlen_c = len2;
    }
    /* Gate: total hits from per-code query counts, without enumerating them. */
    int c = (ref[0] << 6) | (ref[1] << 4) | (ref[2] << 2) | ref[3];
    int nhits = 0;
    for (int i = 4; i < len1; i++) { c = ((c << 2) | ref[i]) & 1023; nhits += s.qcnt[c]; }
    if (nhits > max_hits) return RESCUE_PRUNE_FULL;

    const int quanta = kswv_query_quantum8(len2), off = quanta, nd = len1 + quanta + 1;
    memset(s.cnt, 0, (size_t)nd * sizeof s.cnt[0]);
    c = (ref[0] << 6) | (ref[1] << 4) | (ref[2] << 2) | ref[3];
    for (int i = 4; i < len1; i++) {
        c = ((c << 2) | ref[i]) & 1023;
        for (int j = s.head[c]; j >= 0; j = s.nxt[j]) {
            const int d = i - j + off;
            if (!s.cnt[d]) s.minrow[d] = (int16_t)(i - 4);
            s.cnt[d]++;
        }
    }
    int best = -1, f = 0;
    for (int d = 0; d < nd; d++) {
        f = (int)s.cnt[d] - 1 + (f > 0 ? f : 0);
        s.fwd[d] = f;
        if (f > best) best = f;
    }
    if (5 + best < minsc) return RESCUE_PRUNE_B1;
    int b = 0;
    for (int d = nd - 1; d >= 0; d--) {
        b = (int)s.cnt[d] - 1 + (b > 0 ? b : 0);
        s.bwd[d] = b;
    }
    int lo = len1, hi = -1;
    for (int d = 0; d < nd;) {
        if (5 + s.fwd[d] + s.bwd[d] - ((int)s.cnt[d] - 1) < minsc) { d++; continue; }
        int ub = 0, i0 = len1, dmax = -1;
        while (d < nd) {
            const int bnd = 5 + s.fwd[d] + s.bwd[d] - ((int)s.cnt[d] - 1);
            if (bnd < minsc) break;
            ub = std::max(ub, bnd);
            if (s.cnt[d]) { i0 = std::min(i0, (int)s.minrow[d]); dmax = d; }
            d++;
        }
        if (dmax < 0) continue;
        lo = std::min(lo, i0);
        hi = std::max(hi, (dmax - off) + quanta - 1 + std::max(0, ub - minsc - 2));
    }
    if (hi < 0) return RESCUE_PRUNE_B1;
    *hb = std::max(0, lo);
    *he = std::min(len1 - 1, hi);
    return RESCUE_PRUNE_B2;
}

/* The per-diagonal arrays behind a RESCUE_PRUNE_B2 decision, for deriving diagonal components at
 * thresholds above minsc (rescue_band.h). Diagonal index x in [0, nd) is the unshifted diagonal
 * d = i - j = x - off, with nd = len1 + off + 1 and off the query's 8-bit quantum. Exactly one of
 * bnd16 (NEON: the bound precomputed) or fwd / bwd (scalar: bnd = 5 + fwd + bwd - (cnt - 1)) is
 * set; mw / hw are set with bnd16. The pointers alias the filter's per-thread scratch, so a view is
 * valid only until the next filter call on the same thread. Readable lengths, which band planning's
 * vector loads rely on: on NEON (the filter accepts nd + 32 <= NeonScratch::CAP) cnt and minrow at
 * least nd + 32 entries and bnd16 at least nd + 64 (NeonScratch::VIEW_PAD); on the scalar path nd,
 * which is all the scalar planner reads. */
struct rescue_prune_view {
    int nd = -1, off = 0;
    const uint16_t *cnt = nullptr;
    const int16_t *minrow = nullptr;
    const int16_t *bnd16 = nullptr;
    const int32_t *fwd = nullptr, *bwd = nullptr;
    const uint64_t *mw = nullptr;   // NEON: bitset of diagonals with bnd >= 19
    const uint64_t *hw = nullptr;   // NEON: bitset of diagonals with bnd >= 19 and a hit
#if defined(__aarch64__)
    /* NEON: the filter's own components at bnd >= 19 with a hit (all of [0, nd)), in diagonal
     * order; the first min(ncomp, ncomp_stored) are in comps. ncomp < 0: not available. */
    const rescue_prune_neon::Comp *comps = nullptr;
    int ncomp = -1, ncomp_stored = 0;
#endif
};

/* The view of a rescue_prune_window_scalar call on (len1, len2) that returned RESCUE_PRUNE_B2. */
static inline rescue_prune_view rescue_prune_scalar_view(const rescue_prune_scratch &s, int len1, int len2)
{
    rescue_prune_view v;
    v.off = kswv_query_quantum8(len2);
    v.nd = len1 + v.off + 1;
    v.cnt = s.cnt; v.minrow = s.minrow; v.fwd = s.fwd; v.bwd = s.bwd;
    return v;
}

#if defined(__aarch64__)
/* The calling thread's NEON filter scratch (tables, per-diagonal arrays, query cache, memo). */
static inline rescue_prune_neon::NeonScratch &rescue_prune_neon_scratch()
{
    static thread_local rescue_prune_neon::NeonScratch ns;
    return ns;
}
#endif

/* Calls on this thread that the NEON filter answered from its repeat memo (0 off aarch64), for
 * BWA3_RESCUE_PRUNE_STATS and the tests. */
static inline uint64_t rescue_prune_memo_hits()
{
#if defined(__aarch64__)
    return rescue_prune_neon_scratch().memo_hits;
#else
    return 0;
#endif
}

/* Decide how much of a rescue window must be computed. Valid only under a scoring
 * rescue_prune_scoring_ok() admits and on the 8-bit kernel path; the caller checks both.
 *   ref, len1   reference window, bases 0-3 (>= 4 is N)
 *   q, len2     oriented mate, bases 0-3 (>= 4 is N)
 *   minsc       the rescue score threshold (min_seed_len * a)
 *   max_hits    return FULL when the window and mate share more 5-mer hits than this: the filter
 *               would cost more than the DP rows it can save
 *   hb, he      inclusive sub-window rows, set for RESCUE_PRUNE_B2
 *   view        optional: for RESCUE_PRUNE_B2, the per-diagonal arrays the decision came from
 *               (rescue_prune_view); left empty (nd = -1) otherwise */
static inline int rescue_prune_window(const uint8_t *ref, int len1, const uint8_t *q, int len2,
                                      int minsc, int max_hits, int *hb, int *he,
                                      rescue_prune_view *view = nullptr)
{
    if (view) *view = rescue_prune_view();
    *hb = *he = -1;
    // minsc < 5: the lemma's base term (5) already exceeds the threshold, so nothing is provable.
    if (minsc < 5 || len1 < 5 || len2 < 5 || len2 > rescue_prune_scratch::QCAP || len1 > 30000)
        return RESCUE_PRUNE_FULL;
#if defined(__aarch64__)
    if (minsc == rescue_prune_neon::MINSC) {  // identical decisions, ~2x faster (rescue_prune_neon.h)
        rescue_prune_neon::NeonScratch &ns = rescue_prune_neon_scratch();
        const rescue_prune_neon::Job jb{len1, len2, 0, 0, -1, -1, ref, q};
        int h, e;
        const rescue_prune_neon::Kind k = rescue_prune_neon::lean_neon(jb, ns, h, e, max_hits);
        if (k == rescue_prune_neon::B1) return RESCUE_PRUNE_B1;
        if (k == rescue_prune_neon::FULL) return RESCUE_PRUNE_FULL;
        if (k == rescue_prune_neon::B2) {
            *hb = h; *he = e;
            if (view) {
                view->off = kswv_query_quantum8(len2);
                view->nd = len1 + view->off + 1;
                view->cnt = ns.cnt; view->minrow = ns.minrow; view->bnd16 = ns.bnd;
                view->mw = ns.mw; view->hw = ns.hw;
                view->comps = ns.comps; view->ncomp = ns.ncomp;
                view->ncomp_stored = std::min(ns.ncomp, (int)rescue_prune_neon::NeonScratch::COMP_CAP);
            }
            return RESCUE_PRUNE_B2;
        }
        // FALLBACK: beyond the NEON path's capacity or int16 range -> int32 scalar filter below.
    }
#endif
    static thread_local rescue_prune_scratch s;
    const int kind = rescue_prune_window_scalar(ref, len1, q, len2, minsc, max_hits, s, hb, he);
    if (view && kind == RESCUE_PRUNE_B2) *view = rescue_prune_scalar_view(s, len1, len2);
    return kind;
}

#endif
