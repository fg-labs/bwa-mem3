// test/unit/test_bandedswa_asym.cpp
//
// γ (D3 asymmetric BS scoring) gate. The SIMD banded SW kernels
// (getScores8 / getScores16) must score an arbitrary, possibly *asymmetric*
// substitution matrix identically to the scalar oracle (scalarBandedSWA, which
// already consumes the full mat[25]). The generic-matrix seam landed (commit
// e62a253): the SIMD path now consumes the full mat[25] rather than a symmetric
// XOR-LUT, so on a bisulfite OT matrix (ref-C x read-T freed) it matches the
// scalar oracle on converted reads.
//
// Expected state (now GREEN): both the symmetric control and the asymmetric
// cases PASS. Convention: target-major, mat[ref*5 + read].

#include <cstdint>
#include <random>
#include <vector>

#include "doctest/doctest.h"
#include "bandedSWA.h"

#if HAVE_BSW_VECTOR_8_16

namespace {

// Standard symmetric nt4 matrix (diagonal = +a, off-diagonal = -b, N = ambig).
void build_sym_mat(int8_t mat[25], int a, int b, int ambig) {
    int k = 0;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) mat[k++] = (i == j) ? (int8_t)a : (int8_t)(-b);
        mat[k++] = (int8_t)ambig;
    }
    for (int j = 0; j < 5; ++j) mat[k++] = (int8_t)ambig;
}

// OT bisulfite matrix: symmetric, then free ref-C x read-T (a C->T conversion is
// not a mismatch). A,C,G,T = 0,1,2,3; target-major mat[ref*5 + read].
void build_ot_mat(int8_t mat[25], int a, int b, int ambig) {
    build_sym_mat(mat, a, b, ambig);
    mat[1 * 5 + 3] = (int8_t)a;   // ref C (1), read T (3) -> match (allowed conversion)
}

struct Out { int score, tle, gtle, qle, gscore, max_off; };

// Generate n converted-read pairs (read = C->T projection of the ref prefix,
// plus sparse real mismatches) and compare getScores<width> to scalarBandedSWA
// under the given matrix. Returns the number of pairs whose SIMD score differs
// from the scalar oracle.
int score_mismatches(int width, const int8_t mat[25], int a, int b,
                     int n, int maxlen, unsigned long seed) {
    const int ambig = -1, o = 6, e = 1, end_bonus = 5, w = 100, zdrop = 100;
    const int STRIDE = maxlen + 256;
    BandedPairWiseSW bsw(o, e, o, e, zdrop, end_bonus, const_cast<int8_t*>(mat), a, b, 1);

    std::vector<uint8_t> ref((size_t)STRIDE * n, 0), qer((size_t)STRIDE * n, 0);
    std::vector<SeqPair> pairs(((n + 63) / 64) * 64);   // SIMD-width round-up + slack
    std::vector<Out> oracle(n);
    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<int> lenD(40, maxlen);

    for (int c = 0; c < n; ++c) {
        int la = lenD(rng), lb = lenD(rng);
        int len1 = la > lb ? la : lb, len2 = la > lb ? lb : la;
        uint8_t *s1 = &ref[(size_t)c * STRIDE];   // target (ref), original bases
        uint8_t *s2 = &qer[(size_t)c * STRIDE];   // query (read)
        for (int i = 0; i < len1; ++i) s1[i] = (uint8_t)(rng() % 4);
        for (int i = 0; i < len2; ++i) {
            uint8_t rb = s1[i];
            uint8_t base = (rb == 1) ? 3 : rb;                 // C->T conversion in the read
            if ((int)(rng() % 100) < 3) base = (uint8_t)(rng() % 4);  // 3% real mismatch
            s2[i] = base;
        }
        SeqPair &sp = pairs[c];
        sp.id = c; sp.len1 = len1; sp.len2 = len2; sp.h0 = 40;
        sp.idr = (int)((size_t)c * STRIDE); sp.idq = (int)((size_t)c * STRIDE);
        sp.seqid = c; sp.regid = c;
        sp.score = sp.tle = sp.gtle = sp.qle = sp.gscore = sp.max_off = -1;
        Out &O = oracle[c];
        O.score = bsw.scalarBandedSWA(len2, s2, len1, s1, w, sp.h0,
                                      &O.qle, &O.tle, &O.gtle, &O.gscore, &O.max_off);
    }
    if (width == 8) bsw.getScores8(pairs.data(), ref.data(), qer.data(), (int32_t)n, 1, w);
    else            bsw.getScores16(pairs.data(), ref.data(), qer.data(), (int32_t)n, 1, w);

    int diff = 0;
    for (int c = 0; c < n; ++c) if (pairs[c].score != oracle[c].score) ++diff;
    return diff;
}

} // namespace

TEST_CASE("bandedSWA symmetric control: getScores8 == scalar (sanity)"
          * doctest::test_suite("unit/bandedswa-asym")) {
    int8_t mat[25]; build_sym_mat(mat, 1, 4, -1);
    CHECK(score_mismatches(8, mat, 1, 4, 4000, 200, 1) == 0);
}

TEST_CASE("bandedSWA getScores8 honors an asymmetric (OT bisulfite) matrix"
          * doctest::test_suite("unit/bandedswa-asym")) {
    int8_t mat[25]; build_ot_mat(mat, 1, 4, -1);
    // GREEN since the generic-matrix seam (e62a253): the SIMD path consumes the
    // full mat[25] and rewards the freed ref-C x read-T cell like the scalar oracle.
    CHECK(score_mismatches(8, mat, 1, 4, 4000, 200, 2) == 0);
}

TEST_CASE("bandedSWA getScores16 honors an asymmetric (OT bisulfite) matrix"
          * doctest::test_suite("unit/bandedswa-asym")) {
    int8_t mat[25]; build_ot_mat(mat, 1, 4, -1);
    CHECK(score_mismatches(16, mat, 1, 4, 4000, 400, 3) == 0);
}

#endif // HAVE_BSW_VECTOR_8_16
