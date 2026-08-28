// test/unit/test_bandedswa_zdrop_gate.cpp
//
// Regression test for the vector banded-SW z-drop gate at `-d 0` (z-drop disabled).
//
// scalarBandedSWA gates its z-drop early-exit on `zdrop > 0` (bandedSWA.cpp, the
// `else if (zdrop > 0)` branch), so with z-drop disabled it never truncates. The
// vector kernels (8-bit and 16-bit, every tier) applied the z-drop test
// UNCONDITIONALLY, so at `zdrop == 0` they compared the running drop against a
// zero threshold and killed a lane as soon as its score fell one point below the
// row max -- truncating alignments the scalar reference runs to completion. The
// two then disagreed on the query-end fields (gscore/gtle) and, once the vector
// stops early, on the local fields too.
//
// The fix gates every vector z-drop mask update on `zdrop > 0`, matching the
// scalar. On the DEFAULT `zdrop = 100` the gate is a no-op (the branch is taken),
// so this changes nothing on the shipping path; it only restores agreement at
// `zdrop == 0`. This case must therefore PASS on the fixed build and FAIL on the
// old one.
//
// Each pair is a matching prefix followed by a divergent tail, so the alignment
// score peaks early and then falls -- exactly the shape that makes an
// unconditional z-drop fire at `zdrop == 0`. Deterministic (fixed RNG seed).
#include <cstdint>
#include <random>
#include <vector>

#include "doctest/doctest.h"
#include "bandedSWA.h"

#if HAVE_BSW_VECTOR_8_16

namespace {

struct Out { int score, tle, gtle, qle, gscore, max_off; };

/// Compare getScores<width> to scalarBandedSWA at zdrop == 0 over n peak-then-drop
/// pairs. Returns the mismatch count (per the documented query-end contract) and
/// reports how many pairs were checked so a vacuous trial cannot pass.
long zdrop_trial(int width, int a, int b, int o, int e, int end_bonus, int w,
                 int minlen, int maxlen, long n, long *checked)
{
    const int STRIDE = maxlen + 256;
    const int zdrop = 0;   // z-drop DISABLED -- the case under test

    int8_t mat[25];
    { int k = 0;
      for (int i = 0; i < 4; ++i) { for (int j = 0; j < 4; ++j) mat[k++] = (i == j) ? (int8_t)a : (int8_t)-b; mat[k++] = -1; }
      for (int j = 0; j < 5; ++j) mat[k++] = -1; }

    BandedPairWiseSW bsw(o, e, o, e, zdrop, end_bonus, mat, (int8_t)a, (int8_t)b, 1);

    std::mt19937_64 rng(width == 8 ? 0x2D809A7Eull : 0x2D8016A7ull);
    std::vector<uint8_t> ref((size_t)STRIDE * n, 0), qer((size_t)STRIDE * n, 0);
    std::vector<SeqPair> pairs(((n + 63) / 64) * 64);   // SIMD-width round-up + slack
    std::vector<Out> oracle(n);
    std::uniform_int_distribution<int> lenD(minlen, maxlen);

    for (long c = 0; c < n; ++c) {
        int la = lenD(rng), lb = lenD(rng);
        int len1 = la > lb ? la : lb, len2 = la > lb ? lb : la;   /* len1 >= len2 */
        uint8_t *s1 = &ref[(size_t)c * STRIDE];   /* target (ref) */
        uint8_t *s2 = &qer[(size_t)c * STRIDE];   /* query  (read) */
        for (int i = 0; i < len1; ++i) s1[i] = (uint8_t)(rng() % 4);
        /* Full-length match with a short mismatch VALLEY near the middle. Without
         * z-drop the optimal alignment crosses the valley and runs to the read end
         * (high score, gscore set); an unconditional z-drop at zdrop==0 fires at the
         * valley dip and dies there, reporting only the pre-valley peak -- so score,
         * tle, qle, gscore and gtle all diverge. The suffix (len - valley) matches
         * outscore the 2-base valley (2a matches lost to 2b), so the global max is
         * genuinely past the valley. */
        for (int i = 0; i < len2; ++i) s2[i] = (i < len1) ? s1[i] : (uint8_t)(rng() % 4);
        int mid = len2 / 2;
        for (int d = 0; d < 2; ++d) {
            int j = mid + d;
            if (j < len2) { uint8_t rb = (j < len1) ? s1[j] : s2[j]; s2[j] = (uint8_t)((rb + 1) % 4); }
        }
        int h0 = a;   /* small seed */
        SeqPair &p = pairs[c];
        p.id = c; p.len1 = len1; p.len2 = len2; p.h0 = h0;
        p.idr = (int)((size_t)c * STRIDE); p.idq = (int)((size_t)c * STRIDE);
        p.seqid = c; p.regid = c;
        p.score = p.tle = p.gtle = p.qle = p.gscore = p.max_off = -1;
        Out &O = oracle[c];
        O.score = bsw.scalarBandedSWA(len2, s2, len1, s1, w, h0,
                                      &O.qle, &O.tle, &O.gtle, &O.gscore, &O.max_off);
    }

    *checked = n;
    if (width == 8) bsw.getScores8(pairs.data(), ref.data(), qer.data(), (int32_t)n, 1, w);
    else            bsw.getScores16(pairs.data(), ref.data(), qer.data(), (int32_t)n, 1, w);

    long diffs = 0;
    for (long c = 0; c < n; ++c) {
        const SeqPair &g = pairs[c];
        const Out &O = oracle[c];
        /* Same query-end contract as the other bandedSWA parity tests: local
         * fields unconditional; gscore/gtle only when a to-end alignment is
         * observable on either side. */
        const bool local_differs = g.score != O.score || g.tle != O.tle ||
                                   g.qle != O.qle || g.max_off != O.max_off;
        const bool toend_observable = g.gscore > 0 || O.gscore > 0;
        const bool toend_differs = toend_observable &&
                                   (g.gscore != O.gscore || g.gtle != O.gtle);
        if (local_differs || toend_differs) {
            if (diffs < 5) {
                MESSAGE("  w" << width << " len1=" << g.len1 << " len2=" << g.len2
                        << " | vec " << g.score << "/" << g.tle << "/" << g.gtle << "/"
                        << g.qle << "/" << g.gscore << "/" << g.max_off
                        << " | scalar " << O.score << "/" << O.tle << "/" << O.gtle << "/"
                        << O.qle << "/" << O.gscore << "/" << O.max_off);
            }
            diffs++;
        }
    }
    return diffs;
}

} // namespace

TEST_CASE("bandedSWA vector z-drop is disabled at zdrop==0 (vector == scalar)"
          * doctest::test_suite("unit/bandedswa-zdrop")) {
    long checked = 0;
    // Small deterministic corpus (256 pairs): every pair is a peak-then-valley
    // shape that the pre-fix unconditional z-drop truncates, so ~all pairs
    // diverge on the old kernel -- a few hundred is ample to reproduce the
    // failure while staying inside the 100 ms unit-test budget. The 8-bit and
    // 16-bit widths test the same invariant with different lengths, so each is a
    // SUBCASE of the one case.
    SUBCASE("getScores8") {
        const long diffs = zdrop_trial(8, 1, 4, 6, 1, 5, 100, 60, 200, 256, &checked);
        CHECK(checked > 0);
        CHECK(diffs == 0);
    }
    SUBCASE("getScores16") {
        const long diffs = zdrop_trial(16, 1, 4, 6, 1, 5, 100, 120, 400, 256, &checked);
        CHECK(checked > 0);
        CHECK(diffs == 0);
    }
}

#else  // !HAVE_BSW_VECTOR_8_16

// Vector 8/16-bit kernels are not built for this tier (scalar-only build), so the
// getScores<width>-vs-scalar parity above has nothing to exercise. Register the
// case with a doctest::skip() decorator so a scalar-only build reports it as an
// explicit SKIP -- not a passing test -- and cannot silently claim "no failures"
// for a gate that never ran. (getScores<width> uses the compile-time baseline
// kernel, so each CI matrix row -- NEON, AVX2, ... -- exercises its own tier; a
// tier is covered by building the suite for it.)
TEST_CASE("bandedSWA vector z-drop gate: SKIPPED (vector 8/16-bit kernels not built)"
          * doctest::test_suite("unit/bandedswa-zdrop")
          * doctest::skip()) {
    // Body is never run under doctest::skip(); the case exists only so the
    // unsupported build state is visible as a skip in the test report.
}

#endif // HAVE_BSW_VECTOR_8_16
