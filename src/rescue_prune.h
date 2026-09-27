/* Exact mate-rescue window pruning.
 *
 * A mate-rescue Smith-Waterman job (mate `q` against reference window `ref`) contributes only
 * score/te/qe (pass 0), tb/qb (pass 1) and score2 (-> csub) to the output, and nothing at all
 * when score < minsc (the post gate fails). This filter proves, from exact 5-mer matches, which
 * part of the window can hold any row whose value matters, so the DP can run on a sub-window
 * (or not at all) with byte-identical results.
 *
 * Lemma (a = 1, b = 4, gap open 6 / extend 1; N and query pad score <= +1 so they are treated
 * as matches): a local alignment scores at most
 *     4 + (#exact 5-mer hits on its cells) - sum over gaps of (2 + |diagonal change|),
 * because a run of r matches holds max(0, r - 4) 5-mers and every mismatch costs 4. Hence an
 * alignment whose diagonals span the interval [d1, d2] scores at most
 *     5 + sum_{d in [d1, d2]} (cnt_d - 1),
 * with cnt_d the number of 5-mer hits on diagonal d = i - j.
 *
 * Decisions:
 *   RESCUE_PRUNE_B1   no diagonal interval reaches minsc: score < minsc is proven, so the job's
 *                     outputs are dead. The caller enqueues it with len1 = 0.
 *   RESCUE_PRUNE_B2   rows [*hb, *he] contain every row whose max can reach minsc. A zero-state
 *                     DP started at *hb reproduces all consumed outputs: every alignment scoring
 *                     >= minsc starts at or after the first 5-run of its component (the prefix
 *                     before it nets <= 0), and the tail after its last hit nets <= -(2 + D) per
 *                     gap, which bounds the extent to (last hit diagonal) + quanta - 1 +
 *                     (ub - minsc - 2).
 *   RESCUE_PRUNE_FULL no pruning (N present, lengths too short, or too many hits to be worth it).
 *
 * The bound and extents were validated against the kswv kernel's own per-row maxima on
 * 31.6 M wgs/wes rescue jobs (zero violations) and by whole-run output identity. */
#ifndef BWA_MEM3_RESCUE_PRUNE_H
#define BWA_MEM3_RESCUE_PRUNE_H

#include <stdint.h>
#include <string.h>
#include <algorithm>
#include <vector>
#include "rescue_prune_neon.h"
#include "rescue_prune_x86.h"

enum { RESCUE_PRUNE_FULL = 0, RESCUE_PRUNE_B1 = 1, RESCUE_PRUNE_B2 = 2 };

/* Scoring parameters of the bound (general-scoring derivation, section 2): match a, mismatch b,
 * gap of length L costs o + e*L per type. The lemma generalises with K - 1 <= min(b, o_del,
 * o_ins + e_ins) / a (a run of r >= K matches holds r - K + 1 K-mer hits, every separator costs at
 * least (K - 1) a), the per-diagonal charge c = min over gap types of min(e, o + e - (K - 1) a), the
 * interval bound a (K - 1) + c + sum_d (a cnt_d - c), and the tail budget (Dmax) of
 * (ub - tau - o_del + (K - 1) a) / e_del deletion columns past the last hit diagonal. At the
 * defaults (1, 4, 6/1, 6/1) these are K = 5, c = 1, 5 + sum (cnt_d - 1) and ub - tau - 2. */
struct rescue_prune_params {
    int K = 5, a = 1, b = 4, c = 1;
    int o_del = 6, e_del = 1, o_ins = 6, e_ins = 1;
    int minsc = 19;
    bool valid = false;
    /* --meth (set_meth): the filter runs on copies of the window and the mate with base conv_from
     * rewritten to conv_to (C -> T for OT, G -> A for OB). Exact matching of the converted copies
     * relates every pair the meth matrices score at most +a (the freed cell, its collapsed mirror,
     * the identity) and nothing that scores -b in all of them except the mirror under genomic /
     * neutral scoring: a superset of the match-like cells, and over-counting hits only loosens the
     * bound, so the lemma holds. */
    int conv_from = -1, conv_to = -1;
    /* --meth, relation-expanded (set_meth_rel): the filter matches the exact relation of the genomic /
     * neutral matrices instead of converted copies. A query base relx relates to the reference bases
     * relx and relx ^ 2 (OT: read T to reference T or C; OB: read A to reference A or G); every other
     * pair is the identity. The freed cell scores +a (genomic) or 0 (neutral) <= a and every
     * unrelated real cell -b, so the lemma holds; unlike the converted copies it does not relate the
     * mirror (reference T / read C), which on mostly unconverted reads (TAPS) is what saturates the
     * windows. Not for collapsed scoring, whose mirror cell scores +a: that needs set_meth. -1: off. */
    int relx = -1;

    int base() const { return a * (K - 1) + c; }
    int weight(int cnt) const { return a * cnt - c; }
    /* Rows past the last hit diagonal's end that an alignment with interval bound ub can still
     * reach at score >= tau. */
    int tail(int ub, int tau) const
    {
        const int n = ub - tau - o_del + (K - 1) * a;
        return n > 0 ? n / e_del : 0;
    }
    /* The SIMD filters take K = 5 and any weights (simd_wt); a <= 16 keeps a cnt within int16. */
    bool simd_ok() const { return valid && K == 5 && a <= 16; }
    rescue_prune_neon::Wt simd_wt() const
    {
        rescue_prune_neon::Wt w;
        w.base = base(); w.a = a; w.c = c; w.toff = o_del - (K - 1) * a; w.e = e_del;
        return w;
    }
    bool default_scoring() const
    {
        return conv_from < 0 && relx < 0 && a == 1 && b == 4 && o_del == 6 && e_del == 1 && o_ins == 6 && e_ins == 1;
    }
    /* --meth: the rescued mate is scored with mat_ot (hyp 1: reference C / read T freed) or mat_ob
     * (hyp 0: reference G / read A freed); collapsed scoring also frees the mirror cell. */
    void set_meth(int hyp)
    {
        conv_from = hyp ? 1 : 2;
        conv_to = hyp ? 3 : 0;
        relx = -1;
    }
    /* The relation-expanded alternative to set_meth (relx above): genomic and neutral scoring only. */
    void set_meth_rel(int hyp)
    {
        conv_from = conv_to = -1;
        relx = hyp ? 3 : 0;
    }

    /* Validity conditions V0-V5 in order (V0 first: it guards every division). k_min is the smallest
     * K accepted: K >= 3 is sound, but below 5 hits are so dense (4x per K step on random sequence)
     * that the bound rarely drops below minsc; production uses 5, the fuzz goes down to 3. */
    static rescue_prune_params from(int a, int b, int o_del, int e_del, int o_ins, int e_ins, int minsc,
                                    int k_min = 5)
    {
        rescue_prune_params p;
        p.a = a; p.b = b; p.o_del = o_del; p.e_del = e_del; p.o_ins = o_ins; p.e_ins = e_ins;
        p.minsc = minsc;
        // V0: positive match, mismatch and extends; int8 score tables (+a, -b, a + shift <= 255).
        if (a < 1 || b < 1 || e_del < 1 || e_ins < 1 || o_del < 0 || o_ins < 0 || a > 127 || b > 128
            || a + std::max(1, b) > 255)
            return p;
        // K: a separator (mismatch, deletion, insertion) must cost at least (K - 1) a (V1-V3).
        const int km1 = std::min(4, std::min(b, std::min(o_del, o_ins + e_ins)) / a);
        p.K = km1 + 1;
        if (p.K < std::max(3, k_min)) return p;
        const int s = p.a * (p.K - 1);
        if (b < s || o_del < s || o_del + e_del < s || o_ins + e_ins < s) return p;   // V1-V3
        p.c = std::min(std::min(e_del, e_ins), std::min(o_del + e_del - s, o_ins + e_ins - s));
        if (p.c <= 0) return p;          // V5: the bound must decay
        if (minsc <= s) return p;        // V4: a hit-free window cannot reach minsc
        p.valid = true;
        return p;
    }
    static rescue_prune_params defaults(int minsc) { return from(1, 4, 6, 1, 6, 1, minsc); }
};

/* Per-thread scratch for the scalar path. Fixed capacity (windows are capped at 30000 rows and
 * queries at 1024 bases), so nothing is allocated per call. The query tables are cached: they
 * depend only on the oriented mate, which repeats across the anchors rescued with it. */
struct rescue_prune_scratch {
    static const int QCAP = 1024, DCAP = 30000 + 1024 + 2;
    /* Under a relation (relx >= 0) one query K-mer relates to up to 2^K reference codes, so the
     * query table holds entries (query position ent_j, next entry of the same code ent_nxt) instead
     * of query positions; a mate needing more than ECAP entries is not pruned (FULL). The SIMD filters
     * use the same cap. */
    static const int ECAP = rescue_prune_neon::REL_ECAP;
    int16_t head[1024], nxt[QCAP];
    int16_t ent_j[ECAP], ent_nxt[ECAP];
    uint16_t qcnt[1024];
    uint8_t qcache[QCAP];
    int qlen_c = -1, qk_c = 0, qrel_c = -1;   // the cached query tables are for this length, K and relx
    bool q_over = false;                        // ... and exceeded ECAP
    uint16_t cnt[DCAP];
    int16_t minrow[DCAP];
    int32_t fwd[DCAP], bwd[DCAP];
    int view_nd = -1;   // nd of the last call that returned B2 (arrays valid), else -1
    int view_base = 5, view_a = 1, view_c = 1, view_minsc = 19;   // that call's bound constants
};

/* Scalar implementation: any valid scoring, and the SIMD filters' fallback. At the parameters the
 * SIMD filters accept (simd_ok) it makes the same decisions as they do. */
static inline int rescue_prune_window_scalar(const uint8_t *ref, int len1, const uint8_t *q, int len2,
                                             const rescue_prune_params &p, int max_hits,
                                             rescue_prune_scratch &s, int *hb, int *he)
{
    const int K = p.K, mask = (1 << (2 * K)) - 1, minsc = p.minsc, base = p.base();
    uint8_t orv = 0;
    for (int i = 0; i < len1; i++) orv |= ref[i];
    for (int j = 0; j < len2; j++) orv |= q[j];
    if (orv & 0xFC) return RESCUE_PRUNE_FULL;  // N present: its score (-1) breaks the lemma's accounting

    if (len2 < K || len1 < K) return RESCUE_PRUNE_FULL;
    const int relx = p.relx;
    if (len2 != s.qlen_c || K != s.qk_c || relx != s.qrel_c || memcmp(q, s.qcache, (size_t)len2) != 0) {
        memset(s.head, 0xFF, sizeof s.head);
        memset(s.qcnt, 0, sizeof s.qcnt);
        s.q_over = false;
        int c = 0;
        for (int j = 0; j < K - 1; j++) c = (c << 2) | q[j];
        if (relx < 0) {
            for (int j = K - 1; j < len2; j++) {
                c = ((c << 2) | q[j]) & mask;
                s.nxt[j] = s.head[c];
                s.head[c] = (int16_t)j;
                s.qcnt[c]++;
            }
        } else {
            /* The reference codes related to the K-mer ending at j: its code with any subset of the
             * relx positions flipped by ^ 2 (M = those positions' high bits), one entry each. A
             * reference K-mer has one code, so every (row, j) hit is counted once. Entries go in with
             * j ascending, so each code's chain is j-descending, as nxt is. */
            int ne = 0;
            for (int j = K - 1; j < len2 && !s.q_over; j++) {
                c = ((c << 2) | q[j]) & mask;
                int M = 0;
                for (int t = 0; t < K; t++)
                    if (q[j - t] == relx) M |= 2 << (2 * t);
                for (int sm = M;; sm = (sm - 1) & M) {
                    if (ne == rescue_prune_scratch::ECAP) { s.q_over = true; break; }
                    const int rc = c ^ sm;
                    s.ent_j[ne] = (int16_t)j;
                    s.ent_nxt[ne] = s.head[rc];
                    s.head[rc] = (int16_t)ne++;
                    s.qcnt[rc]++;
                    if (!sm) break;
                }
            }
        }
        memcpy(s.qcache, q, (size_t)len2);
        s.qlen_c = len2;
        s.qk_c = K;
        s.qrel_c = relx;
    }
    if (s.q_over) return RESCUE_PRUNE_FULL;
    int c0 = 0;
    for (int i = 0; i < K - 1; i++) c0 = (c0 << 2) | ref[i];
    /* Gate: total hits from per-code query counts, without enumerating them. */
    int c = c0, nhits = 0;
    for (int i = K - 1; i < len1; i++) { c = ((c << 2) | ref[i]) & mask; nhits += s.qcnt[c]; }
    if (nhits > max_hits) return RESCUE_PRUNE_FULL;

    const int quanta = ((len2 + 15) / 16) * 16, off = quanta, nd = len1 + quanta + 1;
    s.view_nd = -1;
    memset(s.cnt, 0, (size_t)nd * sizeof s.cnt[0]);
    c = c0;
    for (int i = K - 1; i < len1; i++) {
        c = ((c << 2) | ref[i]) & mask;
        if (relx < 0) {
            for (int j = s.head[c]; j >= 0; j = s.nxt[j]) {
                const int d = i - j + off;
                if (!s.cnt[d]) s.minrow[d] = (int16_t)(i - (K - 1));
                s.cnt[d]++;
            }
        } else {
            for (int e = s.head[c]; e >= 0; e = s.ent_nxt[e]) {
                const int d = i - s.ent_j[e] + off;
                if (!s.cnt[d]) s.minrow[d] = (int16_t)(i - (K - 1));
                s.cnt[d]++;
            }
        }
    }
    /* Kadane over the diagonal weights a cnt_d - c: fwd[d] / bwd[d] are the best sums of an
     * interval ending / starting at d. */
    int best = -p.c, f = 0;
    for (int d = 0; d < nd; d++) {
        f = p.weight(s.cnt[d]) + (f > 0 ? f : 0);
        s.fwd[d] = f;
        if (f > best) best = f;
    }
    if (base + best < minsc) return RESCUE_PRUNE_B1;
    int b = 0;
    for (int d = nd - 1; d >= 0; d--) {
        b = p.weight(s.cnt[d]) + (b > 0 ? b : 0);
        s.bwd[d] = b;
    }
    int lo = len1, hi = -1;
    for (int d = 0; d < nd;) {
        if (base + s.fwd[d] + s.bwd[d] - p.weight(s.cnt[d]) < minsc) { d++; continue; }
        int ub = 0, i0 = len1, dmax = -1;
        while (d < nd) {
            const int bnd = base + s.fwd[d] + s.bwd[d] - p.weight(s.cnt[d]);
            if (bnd < minsc) break;
            ub = std::max(ub, bnd);
            if (s.cnt[d]) { i0 = std::min(i0, (int)s.minrow[d]); dmax = d; }
            d++;
        }
        if (dmax < 0) continue;
        lo = std::min(lo, i0);
        hi = std::max(hi, (dmax - off) + quanta - 1 + p.tail(ub, minsc));
    }
    if (hi < 0) return RESCUE_PRUNE_B1;
    *hb = std::max(0, lo);
    *he = std::min(len1 - 1, hi);
    s.view_nd = nd;
    s.view_base = base; s.view_a = p.a; s.view_c = p.c; s.view_minsc = minsc;
    return RESCUE_PRUNE_B2;
}

/* Per-diagonal arrays of the last rescue_prune_window() call on this thread, for deriving
 * diagonal components at thresholds above minsc (rescue_band.h). Valid only right after a call
 * that returned RESCUE_PRUNE_B2, and only when nd >= 0 (the SIMD paths' rare scalar fallback
 * leaves nd = -1). Diagonal index x in [0, nd) is the unshifted diagonal d = i - j = x - off.
 * Exactly one of bnd16 (SIMD filter: bnd precomputed) or fwd/bwd (scalar: bnd = base + fwd + bwd -
 * (a cnt - c)) is set; base / a / c are the call's constants either way. mw, when set, is the SIMD filter's bitset of diagonals with bnd >= minsc
 * (the minsc of the filter call that produced the view). */
struct rescue_prune_view {
    int nd = -1, off = 0;
    int minsc = 0;   // the filter call's threshold (SIMD: of mw, hw and comps)
    int base = 5, a = 1, c = 1;   // bound constants: bnd = base + fwd + bwd - (a cnt - c)
    const uint16_t *cnt = nullptr;
    const int16_t *minrow = nullptr;
    const int16_t *bnd16 = nullptr;
    const int32_t *fwd = nullptr, *bwd = nullptr;
    const uint64_t *mw = nullptr;   // SIMD filter: bitset of diagonals with bnd >= minsc
    const uint64_t *hw = nullptr;   // SIMD filter: bitset of diagonals with bnd >= minsc and a hit
    /* SIMD filter: its own components at bnd >= minsc with a hit (all of [0, nd)), in diagonal
     * order; the first min(ncomp, ncomp_stored) are in comps. ncomp < 0: not available. */
    const rescue_prune_neon::Comp *comps = nullptr;
    int ncomp = -1, ncomp_stored = 0;
};

static inline rescue_prune_scratch &rescue_prune_scalar_scratch()
{
    static thread_local rescue_prune_scratch s;
    return s;
}
#if defined(__aarch64__)
static inline rescue_prune_neon::NeonScratch &rescue_prune_neon_scratch()
{
    static thread_local rescue_prune_neon::NeonScratch ns;
    return ns;
}
#elif defined(__AVX2__)
static inline rescue_prune_x86::X86Scratch &rescue_prune_x86_scratch()
{
    static thread_local rescue_prune_x86::X86Scratch xs;
    return xs;
}
#endif
/* 0 = none, 1 = scalar, 2 = NEON, 3 = x86: which scratch the last rescue_prune_window() call used. */
static inline int &rescue_prune_last_path()
{
    static thread_local int p = 0;
    return p;
}

static inline rescue_prune_view rescue_prune_last_view()
{
    rescue_prune_view v;
    const int path = rescue_prune_last_path();
#if defined(__aarch64__)
    if (path == 2) {
        const rescue_prune_neon::NeonScratch &ns = rescue_prune_neon_scratch();
        if (ns.view_nd < 0) return v;
        v.nd = ns.view_nd;
        v.off = ((ns.qlen_c + 15) / 16) * 16;   // off = quanta of the last (cached) query
        v.minsc = ns.view_minsc;
        v.base = ns.view_wt.base; v.a = ns.view_wt.a; v.c = ns.view_wt.c;
        v.cnt = ns.cnt; v.minrow = ns.minrow; v.bnd16 = ns.bnd; v.mw = ns.mw; v.hw = ns.hw;
        v.comps = ns.comps; v.ncomp = ns.ncomp;
        v.ncomp_stored = std::min(ns.ncomp, (int)rescue_prune_neon::NeonScratch::COMP_CAP);
        return v;
    }
#elif defined(__AVX2__)
    if (path == 3) {
        const rescue_prune_x86::X86Scratch &xs = rescue_prune_x86_scratch();
        if (xs.view_nd < 0) return v;
        v.nd = xs.view_nd;
        v.off = ((xs.qlen_c + 15) / 16) * 16;   // off = quanta of the last (cached) query
        v.minsc = xs.view_minsc;
        v.base = xs.view_wt.base; v.a = xs.view_wt.a; v.c = xs.view_wt.c;
        v.cnt = xs.cnt; v.minrow = xs.minrow; v.bnd16 = xs.bnd; v.mw = xs.mw; v.hw = xs.hw;
        v.comps = xs.comps; v.ncomp = xs.ncomp;
        v.ncomp_stored = std::min(xs.ncomp, (int)rescue_prune_x86::X86Scratch::COMP_CAP);
        return v;
    }
#endif
    if (path == 1) {
        const rescue_prune_scratch &s = rescue_prune_scalar_scratch();
        if (s.view_nd < 0) return v;
        v.nd = s.view_nd;
        v.off = ((s.qlen_c + 15) / 16) * 16;
        v.minsc = s.view_minsc; v.base = s.view_base; v.a = s.view_a; v.c = s.view_c;
        v.cnt = s.cnt; v.minrow = s.minrow; v.fwd = s.fwd; v.bwd = s.bwd;
    }
    return v;
}

/* Decide how much of a rescue window must be computed.
 *   ref, len1   reference window, bases 0-3 (>= 4 is N)
 *   q, len2     oriented mate, bases 0-3 (>= 4 is N)
 *   p           the scoring and the rescue score threshold p.minsc (min_seed_len * a); an invalid
 *               p (rescue_prune_params::from refused it) always gives FULL
 *   max_hits    return FULL when the window and mate share more K-mer hits than this: the filter
 *               would cost more than the DP rows it can save
 *   hb, he      inclusive sub-window rows, set for RESCUE_PRUNE_B2 */
static inline int rescue_prune_window(const uint8_t *ref, int len1, const uint8_t *q, int len2,
                                      const rescue_prune_params &p, int max_hits, int *hb, int *he)
{
    *hb = *he = -1;
    rescue_prune_last_path() = 0;
    if (!p.valid || len1 < 5 || len2 < 5 || len2 > rescue_prune_scratch::QCAP || len1 > 30000)
        return RESCUE_PRUNE_FULL;
    if (p.conv_from >= 0) {   // --meth: exact matching of converted copies (see conv_from)
        static thread_local uint8_t cref[30000], cq[rescue_prune_scratch::QCAP];
        const uint8_t f = (uint8_t)p.conv_from, t = (uint8_t)p.conv_to;
        for (int i = 0; i < len1; i++) cref[i] = ref[i] == f ? t : ref[i];
        for (int j = 0; j < len2; j++) cq[j] = q[j] == f ? t : q[j];
        rescue_prune_params pi = p;
        pi.conv_from = pi.conv_to = -1;
        return rescue_prune_window(cref, len1, cq, len2, pi, max_hits, hb, he);
    }
    const int minsc = p.minsc;
    (void)minsc;
    if (!p.simd_ok()) goto scalar;
#if defined(__aarch64__)
    {   // identical decisions to the scalar filter at any minsc, ~2x faster (rescue_prune_neon.h)
        rescue_prune_neon::NeonScratch &ns = rescue_prune_neon_scratch();
        rescue_prune_last_path() = 2;
        const rescue_prune_neon::Job jb{len1, len2, 0, 0, -1, -1, ref, q};
        int h, e;
        const rescue_prune_neon::Kind k = rescue_prune_neon::lean_neon(jb, ns, h, e, max_hits, minsc, p.simd_wt(), p.relx);
        if (k == rescue_prune_neon::SCALAR) goto scalar;   // int32 path: > 32000 hits, long windows
        if (k == rescue_prune_neon::B1) return RESCUE_PRUNE_B1;
        if (k == rescue_prune_neon::FULL) return RESCUE_PRUNE_FULL;
        *hb = h; *he = e;
        return RESCUE_PRUNE_B2;
    }
#elif defined(__AVX2__)
    {   // identical decisions to the scalar filter at any minsc (rescue_prune_x86.h)
        rescue_prune_x86::X86Scratch &xs = rescue_prune_x86_scratch();
        rescue_prune_last_path() = 3;
        const rescue_prune_neon::Job jb{len1, len2, 0, 0, -1, -1, ref, q};
        int h, e;
        const rescue_prune_neon::Kind k = rescue_prune_x86::lean_x86(jb, xs, h, e, max_hits, minsc, p.simd_wt(), p.relx);
        if (k == rescue_prune_neon::SCALAR) goto scalar;   // int32 path: > 32000 hits, long windows
        if (k == rescue_prune_neon::B1) return RESCUE_PRUNE_B1;
        if (k == rescue_prune_neon::FULL) return RESCUE_PRUNE_FULL;
        *hb = h; *he = e;
        return RESCUE_PRUNE_B2;
    }
#endif
scalar:
    rescue_prune_last_path() = 1;
    return rescue_prune_window_scalar(ref, len1, q, len2, p, max_hits, rescue_prune_scalar_scratch(),
                                      hb, he);
}

/* The default scoring at threshold minsc (the SIMD filters' domain). */
static inline int rescue_prune_window(const uint8_t *ref, int len1, const uint8_t *q, int len2,
                                      int minsc, int max_hits, int *hb, int *he)
{
    return rescue_prune_window(ref, len1, q, len2, rescue_prune_params::defaults(minsc), max_hits,
                               hb, he);
}
static inline int rescue_prune_window_scalar(const uint8_t *ref, int len1, const uint8_t *q, int len2,
                                             int minsc, int max_hits, rescue_prune_scratch &s,
                                             int *hb, int *he)
{
    return rescue_prune_window_scalar(ref, len1, q, len2, rescue_prune_params::defaults(minsc),
                                      max_hits, s, hb, he);
}

#endif
