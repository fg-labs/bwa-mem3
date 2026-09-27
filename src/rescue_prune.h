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
 * --meth (rescue_prune_params::set_meth): the filter matches C -> T (OT) or G -> A (OB) converted
 * copies of the window and the mate exactly, which counts every cell the meth matrices score +a or
 * more as a hit (and some more, which only loosens the bound).
 *
 * Other scorings (rescue_prune_params): the same argument with K-mers, K - 1 <= min(b, o_del,
 * o_ins + e_ins) / a (capped at 5), a per-diagonal charge c and the constant a (K - 1) + c in place
 * of 1 and 5, and the tail below scaled to the deletion costs; rescue_prune_params::from refuses a
 * scoring the argument does not hold for (then every window runs in full).
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
 *                     (ub - minsc - 2) at the defaults (rescue_prune_params::tail in general),
 *                     quanta being the 8-bit kernels' query padding, which is at least the
 *                     16-bit kernels' (kswv_quantum.h), so the bound holds on both widths.
 *   RESCUE_PRUNE_FULL no pruning (N present, lengths too short or too long, a scoring or threshold
 *                     the lemma does not hold for -- at the defaults minsc < 5 -- or too many hits
 *                     to be worth it).
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
#include "rescue_prune_x86.h"

enum { RESCUE_PRUNE_FULL = 0, RESCUE_PRUNE_B1 = 1, RESCUE_PRUNE_B2 = 2 };

/* kswv's 8-bit bias for the score table {+a, -b, -1}: -min(-b, -1) (a --meth matrix's minimum is
 * still -b, mem_opt_fill_meth_mat). */
static inline int kswv8_shift(int b) { return std::max(1, b); }

/* Whether kswv's 8-bit kernels run a scoring exactly: positive match, mismatch and extends,
 * non-negative opens, the 8-bit score table (a <= 127, b <= 128, a + shift <= 255), and gap
 * constants that fit the byte lanes kswv loads them into (o + e <= 255 per gap type: kswv.cpp
 * broadcasts o + e and e as u8, so a larger sum wraps there). Every exact rescue shortcut is argued
 * against those kernels, so both pruning (rescue_prune_params::from) and the band kernels
 * (rb_scoring::valid) require it; a scoring it refuses keeps the full-window kswv path. */
static inline bool kswv8_scoring_ok(int a, int b, int o_del, int e_del, int o_ins, int e_ins)
{
    return a >= 1 && b >= 1 && e_del >= 1 && e_ins >= 1 && o_del >= 0 && o_ins >= 0 && a <= 127 && b <= 128
           && a + kswv8_shift(b) <= 255 && o_del + e_del <= 255 && o_ins + e_ins <= 255;
}

/* Scoring parameters of the bound: match a, mismatch b, a gap of length L costs o + e L per type.
 * The lemma generalises with K - 1 <= min(b, o_del, o_ins + e_ins) / a (a run of r >= K matches
 * holds r - K + 1 K-mer hits, and every separator -- mismatch, deletion, insertion -- costs at least
 * (K - 1) a), the per-diagonal charge c = min over gap types of min(e, o + e - (K - 1) a), the
 * interval bound a (K - 1) + c + sum_d (a cnt_d - c), and a tail of
 * (ub - tau - o_del + (K - 1) a) / e_del deletion columns past the last hit diagonal. At the
 * defaults (1, 4, 6/1, 6/1) these are K = 5, c = 1, 5 + sum (cnt_d - 1) and ub - tau - 2. The
 * filter itself takes only these parameters, so a scoring from() refuses (valid false) must not be
 * filtered: every caller checks valid, and rescue_prune_window returns FULL on an invalid p. */
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
     * bound, so the lemma holds. Default collapsed scoring (b = 2a) forces K = 3 and is refused;
     * --meth -B 4 and the genomic / neutral scorings (b = 4a) are not. */
    int conv_from = -1, conv_to = -1;

    int base() const { return a * (K - 1) + c; }
    int weight(int cnt) const { return a * cnt - c; }
    /* Rows past the last hit diagonal's end that an alignment with interval bound ub can still
     * reach at score >= tau. */
    int tail(int ub, int tau) const
    {
        const int n = ub - tau - o_del + (K - 1) * a;
        return n <= 0 ? 0 : e_del == 1 ? n : n / e_del;   // no division at the default -E 1
    }
    /* The SIMD filters take K = 5 and any weights (simd_wt); a <= 16 keeps a cnt in int16 (they
     * fall back to the scalar filter where a sum would leave it). */
    bool simd_ok() const { return valid && K == 5 && a <= 16; }
    rescue_prune_neon::Wt simd_wt() const
    {
        rescue_prune_neon::Wt w;
        w.base = base(); w.a = a; w.c = c; w.toff = o_del - (K - 1) * a; w.e = e_del;
        return w;
    }
    /* --meth: the rescued mate is scored with mat_ot (hyp 1: reference C / read T freed) or mat_ob
     * (hyp 0: reference G / read A freed); collapsed scoring also frees the mirror cell. */
    void set_meth(int hyp)
    {
        conv_from = hyp ? 1 : 2;
        conv_to = hyp ? 3 : 0;
    }

    /* The validity conditions in order, kswv8_scoring_ok first so that no division sees a zero.
     * k_min is the smallest K accepted: K >= 3 is sound, but below 5 hits are so dense (4x per K
     * step on random sequence) that the bound rarely drops below minsc. Production takes the
     * default 5, and the SIMD filters take only K = 5 (simd_ok); k_min below 5, and with it the
     * scalar filter's K = 3 / 4 path, exist for the band harness's random scorings, which check
     * the lemma's generality against kswv. */
    static rescue_prune_params from(int a, int b, int o_del, int e_del, int o_ins, int e_ins, int minsc,
                                    int k_min = 5)
    {
        rescue_prune_params p;
        p.a = a; p.b = b; p.o_del = o_del; p.e_del = e_del; p.o_ins = o_ins; p.e_ins = e_ins;
        p.minsc = minsc;
        if (!kswv8_scoring_ok(a, b, o_del, e_del, o_ins, e_ins)) return p;
        // K: a separator (mismatch, deletion, insertion) must cost at least (K - 1) a.
        const int km1 = std::min(4, std::min(b, std::min(o_del, o_ins + e_ins)) / a);
        p.K = km1 + 1;
        if (p.K < std::max(3, k_min)) return p;
        const int s = p.a * (p.K - 1);
        if (b < s || o_del < s || o_del + e_del < s || o_ins + e_ins < s) return p;
        p.c = std::min(std::min(e_del, e_ins), std::min(o_del + e_del - s, o_ins + e_ins - s));
        if (p.c <= 0) return p;          // the bound must decay along a hit-free stretch
        if (minsc <= s) return p;        // a hit-free window must not reach minsc (minsc >= 5 at defaults)
        p.valid = true;
        return p;
    }
    static rescue_prune_params defaults(int minsc) { return from(1, 4, 6, 1, 6, 1, minsc); }
};

/* Per-thread scratch for the scalar path. Fixed capacity (windows are capped at WCAP rows and
 * queries at QCAP bases), so nothing is allocated per call. The query tables are cached: they
 * depend only on the oriented mate, which repeats across the anchors rescued with it. */
struct rescue_prune_scratch {
    static const int WCAP = 30000, QCAP = 1024, DCAP = WCAP + QCAP + 2;
    int16_t head[1024], nxt[QCAP];
    uint16_t qcnt[1024];
    uint8_t qcache[QCAP];
    int qlen_c = -1, qk_c = 0;   // the cached query tables are for this length and K
    uint16_t cnt[DCAP];
    int16_t minrow[DCAP];
    int32_t fwd[DCAP], bwd[DCAP];
};

/* Scalar implementation: any valid scoring (p.valid, which the caller checks), a build without a
 * SIMD filter, and beyond the SIMD filter's capacity or int16 range. At the parameters the SIMD
 * filters take (simd_ok) it makes the same decisions as they do. */
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
    if (len2 != s.qlen_c || K != s.qk_c || memcmp(q, s.qcache, (size_t)len2) != 0) {
        memset(s.head, 0xFF, sizeof s.head);
        memset(s.qcnt, 0, sizeof s.qcnt);
        int c = 0;
        for (int j = 0; j < K - 1; j++) c = (c << 2) | q[j];
        for (int j = K - 1; j < len2; j++) {
            c = ((c << 2) | q[j]) & mask;
            s.nxt[j] = s.head[c];
            s.head[c] = (int16_t)j;
            s.qcnt[c]++;
        }
        memcpy(s.qcache, q, (size_t)len2);
        s.qlen_c = len2;
        s.qk_c = K;
    }
    int c0 = 0;
    for (int i = 0; i < K - 1; i++) c0 = (c0 << 2) | ref[i];
    /* Gate: total hits from per-code query counts, without enumerating them. */
    int c = c0, nhits = 0;
    for (int i = K - 1; i < len1; i++) { c = ((c << 2) | ref[i]) & mask; nhits += s.qcnt[c]; }
    if (nhits > max_hits) return RESCUE_PRUNE_FULL;

    const int quanta = kswv_query_quantum8(len2), off = quanta, nd = len1 + quanta + 1;
    memset(s.cnt, 0, (size_t)nd * sizeof s.cnt[0]);
    c = c0;
    for (int i = K - 1; i < len1; i++) {
        c = ((c << 2) | ref[i]) & mask;
        for (int j = s.head[c]; j >= 0; j = s.nxt[j]) {
            const int d = i - j + off;
            if (!s.cnt[d]) s.minrow[d] = (int16_t)(i - (K - 1));
            s.cnt[d]++;
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
    return RESCUE_PRUNE_B2;
}

/* The per-diagonal arrays behind a RESCUE_PRUNE_B2 decision, for deriving diagonal components at
 * thresholds above minsc (rescue_band.h). Diagonal index x in [0, nd) is the unshifted diagonal
 * d = i - j = x - off, with nd = len1 + off + 1 and off the query's 8-bit quantum. Exactly one of
 * bnd16 (SIMD filter: the bound precomputed) or fwd / bwd (scalar: bnd = base + fwd + bwd -
 * (a cnt - c), with the call's constants base / a / c, which are 5 / 1 / 1 at the default scoring)
 * is set; mw / hw are set with bnd16. minsc is the threshold of the filter call that produced the view
 * (on the SIMD filter, the threshold of mw, hw and comps). The pointers alias the filter's per-thread scratch, so a view is
 * valid only until the next filter call on the same thread. Readable lengths, which band planning's
 * vector loads rely on: on the SIMD filter (it accepts nd + 32 <= NeonScratch::CAP / X86Scratch::CAP) cnt and
 * minrow at least nd + 32 entries and bnd16 at least nd + 64 (the scratch's VIEW_PAD); on the scalar path nd,
 * which is all the scalar planner reads. */
struct rescue_prune_view {
    int nd = -1, off = 0;
    /* The SIMD filter answered this call from its repeat memo: the window bytes, the mate bytes,
     * the lengths, the hit gate, the threshold and the bound weights equal the previous filter
     * call's on this thread. */
    bool repeat = false;
    /* keyed: the SIMD filter decided, and key summarizes its inputs (window, mate, lengths, hit
     * gate, threshold, bound weights). Equal inputs have equal keys; a caller finding an earlier job
     * with this key compares the bytes before treating it as the same job (rescue_prune_neon.h
     * whash / qhash). */
    bool keyed = false;
    uint64_t key = 0;
    int minsc = 0;   // the filter call's threshold
    int base = 5, a = 1, c = 1;   // the call's bound constants: bnd = base + fwd + bwd - (a cnt - c)
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

/* The view of a rescue_prune_window_scalar call on (len1, len2) under p that returned
 * RESCUE_PRUNE_B2. */
static inline rescue_prune_view rescue_prune_scalar_view(const rescue_prune_scratch &s, int len1, int len2,
                                                         const rescue_prune_params &p)
{
    rescue_prune_view v;
    v.minsc = p.minsc;
    v.base = p.base(); v.a = p.a; v.c = p.c;
    v.off = kswv_query_quantum8(len2);
    v.nd = len1 + v.off + 1;
    v.cnt = s.cnt; v.minrow = s.minrow; v.fwd = s.fwd; v.bwd = s.bwd;
    return v;
}

/* The SIMD filter this build dispatches to, at any threshold: NEON on aarch64
 * (rescue_prune_neon.h), the SSE4.1 / SSSE3 port on the AVX2-floor x86 builds (rescue_prune_x86.h),
 * none elsewhere. Both return the scalar filter's decisions and the same view. */
#if defined(__aarch64__)
#define RESCUE_PRUNE_HAVE_SIMD 1
typedef rescue_prune_neon::NeonScratch rescue_prune_simd_scratch_t;
static inline rescue_prune_neon::Kind rescue_prune_simd_lean(const rescue_prune_neon::Job &jb,
                                                             rescue_prune_simd_scratch_t &s, int &hb,
                                                             int &he, int max_hits, int minsc,
                                                             const rescue_prune_neon::Wt &wt)
{
    return rescue_prune_neon::lean_neon(jb, s, hb, he, max_hits, minsc, wt);
}
#elif defined(__AVX2__)
#define RESCUE_PRUNE_HAVE_SIMD 1
typedef rescue_prune_x86::X86Scratch rescue_prune_simd_scratch_t;
static inline rescue_prune_neon::Kind rescue_prune_simd_lean(const rescue_prune_neon::Job &jb,
                                                             rescue_prune_simd_scratch_t &s, int &hb,
                                                             int &he, int max_hits, int minsc,
                                                             const rescue_prune_neon::Wt &wt)
{
    return rescue_prune_x86::lean_x86(jb, s, hb, he, max_hits, minsc, wt);
}
#else
#define RESCUE_PRUNE_HAVE_SIMD 0
#endif

#if RESCUE_PRUNE_HAVE_SIMD
/* The calling thread's SIMD filter scratch (tables, per-diagonal arrays, query cache, memo). */
static inline rescue_prune_simd_scratch_t &rescue_prune_simd_scratch()
{
    static thread_local rescue_prune_simd_scratch_t ss;
    return ss;
}
#endif

/* Calls on this thread that the SIMD filter answered from its repeat memo (0 without one), for
 * BWA3_RESCUE_PRUNE_STATS and the tests. */
static inline uint64_t rescue_prune_memo_hits()
{
#if RESCUE_PRUNE_HAVE_SIMD
    return rescue_prune_simd_scratch().memo_hits;
#else
    return 0;
#endif
}

/* The default hit gate (BWA3_RESCUE_PRUNE_MAX_HITS unset), for a run with banding on or off, under
 * --meth (meth) or not: 400 for the hull path. On aarch64, 1000 when banding is on and the run's
 * pruned windows can be banded (no --meth): banding turns more of the pruned windows into savings,
 * which pays for the filter on the denser windows (best of {400, 1000, 3000, 10^4, 10^9} measured on
 * WGS-like data at minsc 19). The SIMD filter runs at every minsc, so the gate does not depend on
 * it; 16-bit jobs (pruned but never banded) and scorings the scalar filter decides (a > 16) take
 * 1000 on aarch64 too. The 16-bit pruning figures (Graviton 4, prune on vs off, -t 16, wgs-5M /
 * wes-5M: -A 2 -10.0 / -14.9 %, -A 3 -B 12 -O 18 -E 3 -9.3 / -13.3 %, every job 16-bit) are whole
 * runs at that gate; 400 for those jobs has not been measured. A --meth pruning run never bands its
 * pruned windows (rescue_band_meth_on), so it takes 400 (Graviton 4, EM-seq panel 5 M pairs, wall,
 * 400 vs 1000: genomic -3.3 %, collapsed -B 4 flat). On x86, 400 with banding too: the x86 kswv
 * kernels are cheap enough that the extra filter work on dense windows does not pay (prune + band,
 * wall, 1000 vs 400: Zen 3 AVX2 wgs-5M 74.52 vs 73.57 s, wes-5M 38.57 vs 37.65 s; Zen 5 AVX-512
 * wgs-5M 27.53 vs 26.73 s). Output is identical at every value. */
static inline int rescue_prune_max_hits_default(bool banding, bool meth)
{
#if defined(__aarch64__)
    return banding && !meth ? 1000 : 400;
#else
    (void)banding; (void)meth;
    return 400;
#endif
}

/* Whether pruning pays for the scoring and threshold p (the cost gate beside the exactness
 * preconditions, p.valid among them, which the caller checks), for a run under --meth (meth, with
 * EM-seq / bisulfite chemistry: emseq) or not, whose kswv runs at the AVX-512BW tier (avx512) or not.
 * On aarch64 everywhere but --meth with other chemistry, the scalar filter included (Graviton 4,
 * prune on vs off, wall, wes-5M / wgs-5M: -B 6 -13.9 / -5.7 %, -O 8 -E 2 -6.9 / -2.9 %,
 * -x intractg -6.4 / -1.1 %). The --meth filter matches converted copies (set_meth), which fits reads
 * whose unmethylated C's are converted; TAPS reads are mostly unconverted, so collapsing them to three
 * letters leaves little to prune (EM-seq panel, 5 M pairs: genomic -9.0 %, collapsed -B 4 -10.1 %;
 * TAPS -0.1 %). On x86 no --meth, where the cheaper kswv leaves nothing to win (Zen 5, EM-seq genomic
 * +0.1 %, collapsed -B 4 +1.6 %, TAPS +5.0 %), and only where the SIMD filter runs (simd_ok: K = 5,
 * a <= 16): the scalar filter costs more than it saves against the cheaper x86 kswv (Zen 5 AVX-512,
 * prune on vs off, scalar-filtered: -O 8 -E 2 +13.3 / +4.7 %, -x intractg +14.6 / +4.9 %), while the
 * SIMD filter at those scorings wins or breaks even (Zen 3 AVX2: -O 8 -E 2 -2.6 / -7.2 %,
 * -x intractg -0.3 / -6.2 %; Zen 5: all within 1 %), and not at the
 * AVX-512BW tier at minsc >= 25: there the 64-lane kswv is cheap and few rescues pass at a high
 * threshold, so the filter costs more than it saves (Zen 5 wgs-5M: +1.0 / +2.4 / +1.2 / +1.3 % at
 * -k 25 / 28 / 32 / 40; wes-5M -1.3 / +0.9 / -0.8 / -1.2 %). AVX2 still wins at -k 32 (Zen 3 wes-5M
 * -3.8 %). */
static inline bool rescue_prune_cost_ok(const rescue_prune_params &p, bool meth, bool emseq, bool avx512)
{
#if defined(__aarch64__)
    (void)p; (void)avx512;
    return !meth || emseq;
#else
    (void)emseq;
    return !meth && p.simd_ok() && !(avx512 && p.minsc >= 25);
#endif
}

/* Decide how much of a rescue window must be computed, for a job on either kswv width. The hull
 * bound counts the query-pad columns of the 8-bit kernels (kswv_query_quantum8) as matches; the
 * 16-bit kernels pad to kswv_query_quantum16, which is never more, and score their pad columns
 * no higher, so on a 16-bit job the bound only over-counts and both decisions stay exact.
 *   ref, len1   reference window, bases 0-3 (>= 4 is N)
 *   q, len2     oriented mate, bases 0-3 (>= 4 is N)
 *   p           the scoring and the rescue score threshold p.minsc (min_seed_len * a); an invalid
 *               p (rescue_prune_params::from refused it) always gives FULL
 *   max_hits    return FULL when the window and mate share more K-mer hits than this: the filter
 *               would cost more than the DP rows it can save
 *   hb, he      inclusive sub-window rows, set for RESCUE_PRUNE_B2
 *   view        optional: for RESCUE_PRUNE_B2, the per-diagonal arrays the decision came from
 *               (rescue_prune_view); left empty (nd = -1) otherwise */
/* With the SIMD filter inlined into it this is a few thousand instructions, and whether the aligner's
 * per-job caller (mem_matesw_batch_pre) should hold it or call it differs by architecture, each way
 * measured against the other in one interleaved round (-t 16, 12 reps, wgs-5M / wes-5M, output
 * identical): on aarch64 the call is faster (Graviton 4, inlined vs called: wall +0.12 % / +0.08 %,
 * user CPU +0.14 % / +0.18 %), on x86 the inlined body is (Zen 5 AVX-512BW, called vs inlined: wall
 * +0.70 % / +1.21 %, user CPU +0.77 % / +0.82 %). The attribute pins each. */
#if defined(__aarch64__)
#define RESCUE_PRUNE_WINDOW_INLINE __attribute__((noinline))
#else
#define RESCUE_PRUNE_WINDOW_INLINE inline
#endif
static RESCUE_PRUNE_WINDOW_INLINE int rescue_prune_window(const uint8_t *ref, int len1, const uint8_t *q, int len2,
                                                          const rescue_prune_params &p, int max_hits, int *hb, int *he,
                                                          rescue_prune_view *view = nullptr)
{
    if (view) *view = rescue_prune_view();
    *hb = *he = -1;
    // An invalid p (at the defaults, minsc < 5: the base term 5 already reaches the threshold, so
    // nothing is provable) keeps the full window.
    if (!p.valid || len1 < 5 || len2 < 5 || len2 > rescue_prune_scratch::QCAP
        || len1 > rescue_prune_scratch::WCAP)
        return RESCUE_PRUNE_FULL;
    if (p.conv_from >= 0) {   // --meth: exact matching of converted copies (see conv_from)
        static thread_local uint8_t cref[rescue_prune_scratch::WCAP], cq[rescue_prune_scratch::QCAP];
        const uint8_t f = (uint8_t)p.conv_from, t = (uint8_t)p.conv_to;
        /* The conversion maps a base to a base (f, t < 4), so an N survives it: check for one in
         * the same pass and skip the filter call that would return FULL for it. The filters below
         * then run on the copies; nothing they read of p depends on the conversion. (Not a
         * recursive call, so the function stays one straight path.) */
        uint8_t orv = 0;
        for (int i = 0; i < len1; i++) { orv |= ref[i]; cref[i] = ref[i] == f ? t : ref[i]; }
        for (int j = 0; j < len2; j++) { orv |= q[j]; cq[j] = q[j] == f ? t : q[j]; }
        if (orv & 0xFC) return RESCUE_PRUNE_FULL;   // N: as the filters decide it
        ref = cref; q = cq;
    }
    const int minsc = p.minsc;
#if RESCUE_PRUNE_HAVE_SIMD
    if (p.simd_ok()) {   // the scalar filter's decisions at any minsc and weights, faster
        rescue_prune_simd_scratch_t &ss = rescue_prune_simd_scratch();
        const rescue_prune_neon::Job jb{len1, len2, 0, 0, -1, -1, ref, q};
        int h, e;
        const uint64_t hits0 = ss.memo_hits;
        const rescue_prune_neon::Wt wt = p.simd_wt();
        const rescue_prune_neon::Kind k = rescue_prune_simd_lean(jb, ss, h, e, max_hits, minsc, wt);
        if (view) {
            view->repeat = ss.memo_hits != hits0;
            view->keyed = k != rescue_prune_neon::FALLBACK;
            const uint64_t wkey = (uint64_t)(uint32_t)wt.base ^ (uint64_t)(uint32_t)wt.a << 16
                                  ^ (uint64_t)(uint32_t)wt.c << 32 ^ (uint64_t)(uint32_t)wt.toff << 40
                                  ^ (uint64_t)(uint32_t)wt.e << 52;
            view->key = (ss.whash ^ (ss.qhash * 0x9E3779B97F4A7C15ULL)) + (uint64_t)max_hits * 0xD6E8FEB86659FD93ULL
                        + (uint64_t)minsc * 0xA24BAED4963EE407ULL + wkey * 0xC2B2AE3D27D4EB4FULL;
        }
        if (k == rescue_prune_neon::B1) return RESCUE_PRUNE_B1;
        if (k == rescue_prune_neon::FULL) return RESCUE_PRUNE_FULL;
        if (k == rescue_prune_neon::B2) {
            *hb = h; *he = e;
            if (view) {
                view->off = kswv_query_quantum8(len2);
                view->nd = len1 + view->off + 1;
                view->minsc = minsc;
                view->base = p.base(); view->a = p.a; view->c = p.c;
                view->cnt = ss.cnt; view->minrow = ss.minrow; view->bnd16 = ss.bnd;
                view->mw = ss.mw; view->hw = ss.hw;
                view->comps = ss.comps; view->ncomp = ss.ncomp;
                view->ncomp_stored = std::min(ss.ncomp, (int)rescue_prune_simd_scratch_t::COMP_CAP);
            }
            return RESCUE_PRUNE_B2;
        }
        // FALLBACK: beyond the SIMD filter's capacity or int16 range -> int32 scalar filter below.
    }
#endif
    static thread_local rescue_prune_scratch s;
    const int kind = rescue_prune_window_scalar(ref, len1, q, len2, p, max_hits, s, hb, he);
    if (view && kind == RESCUE_PRUNE_B2) *view = rescue_prune_scalar_view(s, len1, len2, p);
    return kind;
}

#endif
