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
#include <string>
#include <vector>

#include "doctest/doctest.h"
#include "bandedSWA.h"
#include "ext_result.h"
#include "simd_dispatch.h"

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

using bwa_tests::ExtResult;

// Generate n converted-read pairs (read = C->T projection of the ref prefix,
// plus sparse real mismatches) and compare getScores<width> to scalarBandedSWA
// under the given matrix. Returns the number of pairs whose SIMD result record
// differs from the scalar oracle on ANY field (score, tle, gtle, qle, gscore,
// max_off) -- an ambiguous-base regression can perturb a secondary field while
// leaving score intact, so score alone is not a sufficient parity gate. When
// `out` is given it receives the kernel's result records (pairs [0, n)); a
// `compact_groups` >= 0 sets the 8-bit lane compaction (0 = off) instead of the
// tier's default. With `compare_scalar` false the scalar oracle is skipped and 0 is
// returned (for an output-only call whose records are compared elsewhere).
int score_mismatches(int width, const int8_t mat[25], int a, int b,
                     int n, int maxlen, unsigned long seed, int n_pct = 0,
                     std::vector<SeqPair> *out = nullptr, int compact_groups = -1,
                     bool compare_scalar = true) {
    const int ambig = -1, o = 6, e = 1, end_bonus = 5, w = 100, zdrop = 100;
    const int STRIDE = maxlen + 256;
    BandedPairWiseSW bsw(o, e, o, e, zdrop, end_bonus, const_cast<int8_t*>(mat), a, b, 1);
    if (compact_groups >= 0) bsw.set_lane_compaction(compact_groups, BSW_COMPACT_MIN_W_DEFAULT);

    std::vector<uint8_t> ref((size_t)STRIDE * n, 0), qer((size_t)STRIDE * n, 0);
    std::vector<SeqPair> pairs(((n + 63) / 64) * 64);   // SIMD-width round-up + slack
    std::vector<ExtResult> oracle(compare_scalar ? n : 0);
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
        // Ambiguous-base (N = nt4 code 4) injection: exercises the query-side
        // 4->8 remap in neon_soa_pack16 (target N stays 4). The scalar oracle
        // scores N through the same mat[25] (row/col 4), so both paths must agree.
        if (n_pct > 0) {
            for (int i = 0; i < len1; ++i)
                if ((int)(rng() % 100) < n_pct) s1[i] = 4;     // target N
            for (int i = 0; i < len2; ++i)
                if ((int)(rng() % 100) < n_pct) s2[i] = 4;     // query N -> code 8 in the kernel
        }
        SeqPair &sp = pairs[c];
        sp.id = c; sp.len1 = len1; sp.len2 = len2; sp.h0 = 40;
        sp.idr = (int)((size_t)c * STRIDE); sp.idq = (int)((size_t)c * STRIDE);
        sp.seqid = c; sp.regid = c;
        sp.score = sp.tle = sp.gtle = sp.qle = sp.gscore = sp.max_off = -1;
        if (compare_scalar) {
            ExtResult &O = oracle[c];
            O.score = bsw.scalarBandedSWA(len2, s2, len1, s1, w, sp.h0,
                                          &O.qle, &O.tle, &O.gtle, &O.gscore, &O.max_off);
        }
    }
    if (width == 8) bsw.getScores8(pairs.data(), ref.data(), qer.data(), (int32_t)n, 1, w);
    else            bsw.getScores16(pairs.data(), ref.data(), qer.data(), (int32_t)n, 1, w);

    if (out) out->assign(pairs.begin(), pairs.begin() + n);
    if (!compare_scalar) return 0;
    int diff = 0;
    for (int c = 0; c < n; ++c) {
        const SeqPair &p = pairs[c];
        const ExtResult &o = oracle[c];
        if (p.score != o.score || p.tle != o.tle || p.gtle != o.gtle ||
            p.qle != o.qle || p.gscore != o.gscore || p.max_off != o.max_off)
            ++diff;
    }
    return diff;
}

} // namespace

// These parity cases deliberately run getScores8/16 on whatever tier the host
// (or BWAMEM3_FORCE_TIER) selects, NOT only ARM NEON: exercising every tier is
// what lets this suite catch a tier-specific kernel divergence -- e.g. a 256-bit
// AVX2 gap-open sign bug that leaked a negative F into the band's zero-scan and
// perturbed gtle, invisible to a NEON-only gate. The NEON path itself is covered
// by the arm64 CI matrix rows. This sentinel makes the active tier visible in the
// test log and asserts a real vector kernel (not the scalar fallback) is running,
// so a misconfigured runner that quietly dropped to scalar fails loudly here.
TEST_CASE("bandedSWA parity gate exposes and asserts the active SIMD tier"
          * doctest::test_suite("unit/bandedswa-asym")) {
    bwamem3_simd_init();
    const int tier = bwamem3_simd_tier();
    MESSAGE("bandedSWA parity gate active SIMD tier: " << std::string(bwamem3_simd_tier_name(tier)));
    CHECK(tier != BWAMEM3_TIER_NONE);
}

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

// Ambiguous-base coverage for the packed 8-bit path. Every other 8-bit
// extension parity case builds sequences from rng()%4, so the query-side N->8
// remap in neon_soa_pack16 (and the target-side "N stays 4" fill) went
// untested. With ~10% N in both target and query, getScores8 must still match
// the scalar oracle, which scores N through the same mat[25] (row/col 4).
// 512 pairs (32 full SIMD batches of 16 lanes) rather than 4000: the scalar
// oracle runs once per pair, so a large corpus can blow the unit-test time
// budget under the sanitizer CI configs. The fixed seed makes the corpus
// deterministic and reproducible, and at 10% N over the generated lengths every
// batch carries N at target and query boundary and interior positions.
TEST_CASE("bandedSWA getScores8 matches scalar with ambiguous (N) bases"
          * doctest::test_suite("unit/bandedswa-asym")) {
    int8_t mat[25]; build_sym_mat(mat, 1, 4, -1);
    CHECK(score_mismatches(8, mat, 1, 4, 512, 200, 4, /*n_pct=*/10) == 0);
}

TEST_CASE("bandedSWA getScores16 matches scalar with ambiguous (N) bases (control)"
          * doctest::test_suite("unit/bandedswa-asym")) {
    int8_t mat[25]; build_sym_mat(mat, 1, 4, -1);
    CHECK(score_mismatches(16, mat, 1, 4, 512, 400, 5, /*n_pct=*/10) == 0);
}

#if defined(__aarch64__)
// Byte-identity guard for the aarch64 diagonal, not a statement of the intended score. The
// 128-bit 8-bit kernel (SSE4.1 and NEON) has always added a cell's substitution score s as the
// SBT_SPLIT8 halves (max(s, 0), max(-s, 0)), which split s = -128 into (0, 0): a -128 score adds
// nothing. On aarch64 it now adds s with one USQADD instead, and maps -128 to 0 up front to keep
// that pre-existing result (BSW8_SCORES_NOMIN8 in bandedSWA.cpp, used by both 8-bit drivers).
// Pin it in each matrix arm, with lane compaction off and on: a matrix whose mismatch cells are
// -128 (-B 128) scores every pair as the one whose mismatch cells are 0.
// (The scalar oracle, the 16-bit kernel and the AVX2 / AVX-512 8-bit kernels subtract the full
// 128. That pre-existing divergence of the 128-bit 8-bit kernel at -B 128 is out of scope here;
// changing it is a deliberate output change that must update this test.)
TEST_CASE("bandedSWA getScores8 (NEON) keeps the pre-existing -128 substitution result (0) in every matrix arm"
          * doctest::test_suite("unit/bandedswa-asym")) {
    // the arms: symmetric (XOR LUT), one freed cell (RANK1), two freed cells (AMAT)
    for (int compact = 0; compact < 2; ++compact)
    for (int arm = 0; arm < 3; ++arm) {
        int8_t m128[25], m0[25];
        for (int k = 0; k < 2; ++k) {
            int8_t *mat = k ? m0 : m128;
            const int b = k ? 0 : 128;
            if (arm == 0) build_sym_mat(mat, 1, b, -1);
            else build_ot_mat(mat, 1, b, -1);
            if (arm == 2) mat[2 * 5 + 0] = 1;   // ref G (2), read A (0) freed as well (OB)
        }
        CAPTURE(arm);
        CAPTURE(compact);
        const int groups = compact ? BSW_COMPACT_GROUPS_NEON : 0;
        REQUIRE(m128[0 * 5 + 1] == -128);
        std::vector<SeqPair> r128, r0;
        // The -128 run is compared to the zero-penalty run, not to the scalar oracle (which
        // subtracts the full 128); the zero-penalty run is itself checked against the oracle,
        // so two identically wrong kernel runs cannot pass.
        (void) score_mismatches(8, m128, 1, 128, 256, 200, 6 + arm, /*n_pct=*/5, &r128, groups,
                                /*compare_scalar=*/false);
        CHECK(score_mismatches(8, m0, 1, 0, 256, 200, 6 + arm, /*n_pct=*/5, &r0, groups) == 0);
        int diff = 0, scored = 0;
        for (size_t c = 0; c < r0.size(); ++c) {
            const SeqPair &x = r128[c], &y = r0[c];
            scored += y.score > 0;
            diff += x.score != y.score || x.tle != y.tle || x.gtle != y.gtle ||
                    x.qle != y.qle || x.gscore != y.gscore || x.max_off != y.max_off;
        }
        CHECK(scored > 0);
        CHECK(diff == 0);
    }
}
#endif

#endif // HAVE_BSW_VECTOR_8_16
