// test/unit/test_kswv_correctness.cpp
//
// Self-consistency test for kswv::getScores8 versus scalar ksw_align2.
// Ported from test/kswv_selftest.cpp, now using the bwa_tests framework.
//
// Locks in bit-identical batched-SIMD output for 10,000 random pairs plus
// curated edge cases. Runs on every CI matrix row (SSE4.1, AVX2, AVX2
// clang, AVX2 no-mimalloc, multi-arch, ARM64 Linux, macOS ARM64).

#include <random>
#include <vector>

#include "doctest/doctest.h"

#include "kswr_cmp.h"
#include "ksw_runner.h"
#include "kswv_runner.h"
#include "scoring.h"
#include "seqpair_gen.h"

#if BWA_TESTS_HAVE_KSWV

namespace {

// Edge-case pair builder. Deterministic given seed.
std::vector<bwa_tests::TestPair> build_edge_cases(std::mt19937 &rng) {
    using bwa_tests::gen_exact_match_pair;
    using bwa_tests::gen_all_mismatch_pair;
    using bwa_tests::gen_homopolymer_pair;
    using bwa_tests::gen_sub_cluster_pair;
    using bwa_tests::gen_with_n_bases_pair;
    using bwa_tests::gen_random_pair;

    std::vector<bwa_tests::TestPair> pairs;
    pairs.push_back(gen_exact_match_pair(50));
    pairs.push_back(gen_exact_match_pair(100));
    pairs.push_back(gen_all_mismatch_pair(50));
    pairs.push_back(gen_all_mismatch_pair(100));
    pairs.push_back(gen_homopolymer_pair(50, 0));
    pairs.push_back(gen_homopolymer_pair(50, 1));
    pairs.push_back(gen_homopolymer_pair(50, 2));
    pairs.push_back(gen_homopolymer_pair(50, 3));
    pairs.push_back(gen_sub_cluster_pair(rng, 100, 150, 40, 3));
    pairs.push_back(gen_sub_cluster_pair(rng, 100, 150, 40, 10));
    pairs.push_back(gen_with_n_bases_pair(rng, 100, 150, 5));
    pairs.push_back(gen_with_n_bases_pair(rng, 100, 150, 20));
    pairs.push_back(gen_random_pair(rng, 20, 40));
    pairs.push_back(gen_random_pair(rng, 10, 30));
    return pairs;
}

std::vector<bwa_tests::TestPair> build_bulk_random(std::mt19937 &rng, int n) {
    std::vector<bwa_tests::TestPair> pairs;
    pairs.reserve(n);
    std::uniform_int_distribution<int> qlen_d(50, 128);
    std::uniform_int_distribution<int> rlen_d(100, 250);
    for (int i = 0; i < n; i++) {
        pairs.push_back(bwa_tests::gen_random_pair(rng, qlen_d(rng), rlen_d(rng)));
    }
    return pairs;
}

} // namespace

TEST_CASE("kswv::getScores8 matches scalar ksw_align2 on 10k random + curated edge pairs"
          * doctest::test_suite("unit/kswv")) {

    auto mat = bwa_tests::build_scoring_matrix(1, 4, 1);

    std::mt19937 rng(42);
    auto pairs = build_edge_cases(rng);
    auto bulk  = build_bulk_random(rng, 10000);
    pairs.insert(pairs.end(), bulk.begin(), bulk.end());

    std::vector<kswr_t> scalar_aln;
    scalar_aln.reserve(pairs.size());
    for (const auto &p : pairs) {
        scalar_aln.push_back(bwa_tests::run_scalar_ksw(p, mat));
    }

    auto batched = bwa_tests::run_kswv_batch(pairs, mat);
    REQUIRE(batched.size() == pairs.size());

    int score_mism = 0, coord_mism = 0, score2_mism = 0;
    for (size_t i = 0; i < pairs.size(); i++) {
        CAPTURE(i);
        CAPTURE(pairs[i].tag);
        CAPTURE(scalar_aln[i].score);
        CAPTURE(batched[i].score);
        const bool score_ok  = bwa_tests::kswr_score_eq(scalar_aln[i], batched[i]);
        const bool coord_ok  = bwa_tests::kswr_coords_eq(scalar_aln[i], batched[i]);
        const bool score2_ok = bwa_tests::kswr_score2_eq(scalar_aln[i], batched[i]);
        // Only emit per-pair CHECKs on mismatch — at 10k pairs the doctest
        // log would otherwise carry ~30k passing assertions per run, which
        // dwarfs the actual signal on a regression.
        if (!score_ok)  { ++score_mism;  CHECK(score_ok); }
        if (!coord_ok)  { ++coord_mism;  CHECK(coord_ok); }
        if (!score2_ok) { ++score2_mism; CHECK(score2_ok); }
    }

    MESSAGE("kswv vs scalar: score_mism=" << score_mism
            << " coord_mism=" << coord_mism
            << " score2_mism=" << score2_mism
            << " over " << pairs.size() << " pairs");
    // Aggregate gates so a regression with no per-pair CHECK still fails
    // the test (e.g. if all 10k pairs happen to mismatch in coords only,
    // those per-pair CHECKs cover it; this catches counting drift too).
    CHECK(score_mism == 0);
    CHECK(coord_mism == 0);
    CHECK(score2_mism == 0);
}

TEST_CASE("kswv handles every curated edge case identically to scalar"
          * doctest::test_suite("unit/kswv")) {

    auto mat = bwa_tests::build_scoring_matrix(1, 4, 1);
    std::mt19937 rng(42);

    // Helper lambda: run scalar+batched on a single pair and check all
    // three comparators (score, coords, score2). Mirrors the bulk test so a
    // regression in any of those dimensions is pinpointed by SUBCASE name.
    auto check_pair_parity = [&](const bwa_tests::TestPair &p) {
        auto s = bwa_tests::run_scalar_ksw(p, mat);
        auto b = bwa_tests::run_kswv_batch({p}, mat);
        CHECK(bwa_tests::kswr_score_eq(s, b[0]));
        CHECK(bwa_tests::kswr_coords_eq(s, b[0]));
        CHECK(bwa_tests::kswr_score2_eq(s, b[0]));
    };

    SUBCASE("exact match len 50")   { check_pair_parity(bwa_tests::gen_exact_match_pair(50)); }
    SUBCASE("exact match len 100")  { check_pair_parity(bwa_tests::gen_exact_match_pair(100)); }
    SUBCASE("all mismatch len 50")  { check_pair_parity(bwa_tests::gen_all_mismatch_pair(50)); }
    SUBCASE("homopolymer A")        { check_pair_parity(bwa_tests::gen_homopolymer_pair(50, 0)); }
    SUBCASE("sub cluster len 10")   { check_pair_parity(bwa_tests::gen_sub_cluster_pair(rng, 100, 150, 40, 10)); }
    SUBCASE("20% N bases")          { check_pair_parity(bwa_tests::gen_with_n_bases_pair(rng, 100, 150, 20)); }
}

#endif // BWA_TESTS_HAVE_KSWV
