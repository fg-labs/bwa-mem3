// Regression test for the 8-bit banded-SW z-drop drift in the HIGH-h0 /
// SMALL-zdrop regime.
//
// The 8-bit kernel's z-drop test (smithWaterman128_8 wide epilogue) mirrors
// scalarBandedSWA's drift term (i-max_i) - (mj-max_j). The vector reconstructs
// max_j from a per-lane best-column side channel that carries a +1 frame bias,
// which cancels against the row-max column's matching +1 bias -- but ONLY once a
// row has beaten the seed (xrow >= 1). While the running best is still the h0
// seed (xrow == 0, scalar's max_i == max_j == -1), the reconstruction
// yc = y128 + (xrow-1) produced -1 instead of the sentinel-correct max_j+1 == 0,
// under-counting the drift by one and firing the z-drop one row too early. That
// diverged from scalar precisely when the seed score h0 sits high relative to a
// small zdrop, so the row max never beats h0 before the (spurious) z-drop -- the
// exact region a relaxed 8-bit routing envelope would newly admit.
//
// This fixture drives that region directly: tandem-repeat pairs (the geometry
// most prone to premature z-drop) with h0 well above zdrop+1 at small zdrop, and
// requires getScores8 to be byte-identical to scalarBandedSWA on EVERY result
// field. The only precondition the byte-wide DP actually needs is that no cell
// overflow [0,255] (h0 + min(len1,len2)*maxStep < 255 - maxStep); the removed
// h0 <= zdrop+1 condition was a re-baseline artifact, so it is deliberately NOT
// imposed here.
//
// Scoring is restricted to the default matrix and moderate gap costs. Harsh
// gap/mismatch settings (e.g. -B9 -O16 -E3) additionally trip an unrelated
// int8-overflow in the wrapper's band-width clamp (tracked separately), which
// would mask this z-drop regression; they are intentionally excluded.
//
// Deterministic (fixed RNG seed); exits non-zero on any mismatch.

#include <cstdio>
#include <cstdint>
#include <vector>
#include <random>

#include "bandedSWA.h"

namespace {

struct Out { int score, tle, gtle, qle, gscore, max_off; };

// Byte-fit precondition of the 8-bit kernel (the only bound the DP needs) plus
// the length / target>=query / band conditions that remain load-bearing once
// the re-baseline-era h0 <= zdrop+1 condition is dropped. Note: h0 is allowed to
// exceed zdrop+1 here -- that is the whole point.
bool envelope_ok(int len1, int len2, int w, int h0, int maxStep) {
    int shorter = len1 < len2 ? len1 : len2;
    return len1 < MAX_SEQ_LEN8 && len2 < MAX_SEQ_LEN8 && len1 >= len2 &&
           w <= 127 && (h0 + shorter * maxStep) < 255 - maxStep;
}

struct Params { int a, b, o, e, zdrop, w; const char *label; };

// One (scoring parameters) x (high-h0 pairs) trial. Returns the number of
// getScores8-vs-scalar mismatches over the admitted pairs. `admitted` and
// `above` report coverage so a vacuous trial is visible rather than a silent pass.
long run_trial(const Params &P, long n, long *admitted_out, long *above_out) {
    const int maxStep   = P.a > 1 ? P.a : 1;
    const int end_bonus = 5;
    const int STRIDE    = 300;

    // Seeds strictly above the old h0 <= zdrop+1 gate, low enough that a useful
    // query still fits under the byte-fit bound.
    const int h0min = P.zdrop + 2;
    const int h0max = (254 - maxStep) / 2;
    if (h0max <= h0min) { *admitted_out = 0; *above_out = 0; return 0; }
    int maxlen = (254 - maxStep - h0max) / maxStep;
    if (maxlen > 120) maxlen = 120;
    if (maxlen < 12)  { *admitted_out = 0; *above_out = 0; return 0; }

    int8_t mat[25];
    { int k = 0;
      for (int i = 0; i < 4; ++i) { for (int j = 0; j < 4; ++j) mat[k++] = (i == j) ? (int8_t)P.a : (int8_t)-P.b; mat[k++] = -1; }
      for (int j = 0; j < 5; ++j) mat[k++] = -1; }

    BandedPairWiseSW bsw(P.o, P.e, P.o, P.e, P.zdrop, end_bonus, mat, (int8_t)P.a, (int8_t)P.b, 1);

    std::mt19937_64 rng(0xEC7C0DE4ull);
    std::vector<uint8_t> ref((size_t)STRIDE * n, 0), qer((size_t)STRIDE * n, 0);
    std::vector<SeqPair> pairs(n);
    std::vector<Out> oracle(n);
    std::uniform_int_distribution<int> lenD(10, maxlen);
    std::uniform_int_distribution<int> hD(h0min, h0max);
    std::uniform_int_distribution<int> unit(3, 12);

    for (long c = 0; c < n; ++c) {
        int aa = lenD(rng), bb = lenD(rng);
        int len1 = aa > bb ? aa : bb, len2 = aa > bb ? bb : aa;   // len1 >= len2
        int h0 = hD(rng);
        uint8_t *s1 = &ref[(size_t)c * STRIDE];
        uint8_t *s2 = &qer[(size_t)c * STRIDE];
        int u = unit(rng); uint8_t ubuf[16];
        for (int i = 0; i < u; i++) ubuf[i] = (uint8_t)(rng() % 4);
        for (int i = 0; i < len1; i++) s1[i] = ubuf[i % u];       // tandem repeat
        int ti = 0;
        for (int i = 0; i < len2; i++) {
            uint8_t base = (ti < len1) ? s1[ti] : (uint8_t)(rng() % 4);
            if ((int)(rng() % 100) < 5) base = (uint8_t)(rng() % 4);
            s2[i] = base;
            if ((rng() % 40) == 0) { if (rng() & 1) ti += 2; } else ti++;
        }
        SeqPair &p = pairs[c];
        p.id = c; p.len1 = len1; p.len2 = len2; p.h0 = h0;
        p.idr = (int)((size_t)c * STRIDE); p.idq = (int)((size_t)c * STRIDE);
        p.seqid = c; p.regid = c;
        p.score = p.tle = p.gtle = p.qle = p.gscore = p.max_off = -1;
        Out &O = oracle[c];
        O.score = bsw.scalarBandedSWA(len2, s2, len1, s1, P.w, h0,
                                      &O.qle, &O.tle, &O.gtle, &O.gscore, &O.max_off);
    }

    std::vector<long> idx; idx.reserve(n);
    long above = 0;
    for (long c = 0; c < n; ++c) {
        const SeqPair &p = pairs[c];
        if (!envelope_ok(p.len1, p.len2, P.w, p.h0, maxStep)) continue;
        idx.push_back(c);
        if (p.h0 > P.zdrop + 1) above++;
    }
    const long adm = (long)idx.size();
    *admitted_out = adm; *above_out = above;
    if (adm == 0) return 0;

    // getScores8 pads the batch to a whole number of SIMD lanes and writes the
    // trailing padding lanes; size to the rounded capacity.
    const long round_adm = ((adm + SIMD_WIDTH8 - 1) / SIMD_WIDTH8) * SIMD_WIDTH8;
    std::vector<SeqPair> ap(round_adm);
    for (long k = 0; k < adm; ++k) ap[k] = pairs[idx[k]];
    bsw.getScores8(ap.data(), ref.data(), qer.data(), (int32_t)adm, 1, P.w);

    long diffs = 0;
    for (long k = 0; k < adm; ++k) {
        const SeqPair &g = ap[k];
        const Out &O = oracle[idx[k]];
        if (g.score != O.score || g.tle != O.tle || g.gtle != O.gtle ||
            g.qle != O.qle || g.gscore != O.gscore || g.max_off != O.max_off) {
            if (diffs < 4)
                fprintf(stderr,
                    "  [%s] DIFF len1=%d len2=%d h0=%d | 8 %d/%d/%d/%d/%d/%d | scalar %d/%d/%d/%d/%d/%d\n",
                    P.label, g.len1, g.len2, g.h0,
                    g.score, g.tle, g.gtle, g.qle, g.gscore, g.max_off,
                    O.score, O.tle, O.gtle, O.qle, O.gscore, O.max_off);
            diffs++;
        }
    }
    return diffs;
}

} // namespace

int main() {
    // Small zdrop with h0 far above zdrop+1 is the failing regime; a few moderate
    // (non-harsh) scoring points guard against a parameter-specific reintroduction.
    const Params sweep[] = {
        {1, 4,  6, 1,   5, 100, "defaults zdrop=5"},
        {1, 4,  6, 1,  10, 100, "defaults zdrop=10"},
        {1, 4,  6, 1,  20, 100, "defaults zdrop=20"},
        {1, 4,  6, 1,  50, 100, "defaults zdrop=50"},
        {1, 4,  6, 1, 100, 100, "defaults zdrop=100 (h0>101)"},
        {2, 8, 12, 2,  20, 100, "-A2 -B8 -O12 -E2 zdrop=20"},
        {3, 4,  6, 1,  20, 100, "-A3 zdrop=20"},
        {1, 4,  6, 1,  10,  20, "defaults zdrop=10 narrow band"},
    };
    const long n = 20000;

    long total_diffs = 0, trials = 0, vacuous = 0, total_above = 0;
    for (const Params &P : sweep) {
        long adm = 0, above = 0;
        long d = run_trial(P, n, &adm, &above);
        total_diffs += d;
        total_above += above;
        trials++;
        if (!adm) vacuous++;
        fprintf(stderr, "  [%-30s] zdrop=%-3d w=%-3d admitted=%6ld h0>gate=%6ld diffs=%ld  %s\n",
                P.label, P.zdrop, P.w, adm, above, d,
                d ? "FAIL" : (adm ? "ok" : "VACUOUS"));
    }

    fprintf(stderr, "[high-h0-zdrop] %ld sets, %ld vacuous, h0>gate total=%ld, diffs=%ld\n",
            trials, vacuous, total_above, total_diffs);

    // The fixture is only meaningful if it actually exercised the newly-admitted
    // region (h0 > zdrop+1) on the small-zdrop sets.
    if (total_above == 0 || vacuous == trials) {
        fprintf(stderr, "bandedswa_high_h0_zdrop_test: FAIL — fixture admitted no h0>zdrop+1 pairs\n");
        return 2;
    }
    if (total_diffs != 0) {
        fprintf(stderr, "bandedswa_high_h0_zdrop_test: FAIL — %ld getScores8/scalar mismatches\n", total_diffs);
        return 1;
    }
    fprintf(stderr, "bandedswa_high_h0_zdrop_test: OK\n");
    return 0;
}
