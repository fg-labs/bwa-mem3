// test/unit/test_kswv_correctness.cpp
//
// Self-consistency test for kswv::getScores8 versus scalar ksw_align2.
// Ported from test/kswv_selftest.cpp, now using the bwa_tests framework.
//
// Locks in bit-identical batched-SIMD output for 10,000 random pairs plus
// curated edge cases. Runs on every CI matrix row (SSE4.1, AVX2, AVX2
// clang, AVX2 no-mimalloc, multi-arch, ARM64 Linux, macOS ARM64).

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <vector>

#include "doctest/doctest.h"

#include "kswr_cmp.h"
#include "ksw_runner.h"
#include "kswv.h"
#include "kswv_runner.h"
#include "scoring.h"
#include "seqpair_batch.h"
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
    using bwa_tests::gen_tandem_repeat_pair;

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
    // Odd query lengths so a column range ends on an odd count and the two-row
    // sweep hits its odd-tail branch (partial final block / biased-body jdummy).
    // 100/150 alone leave even partial blocks (100%16=4, 150%16=6).
    pairs.push_back(gen_random_pair(rng, 101, 151));
    pairs.push_back(gen_random_pair(rng, 151, 201));
    pairs.push_back(gen_random_pair(rng, 157, 200));
    // Tandem repeats give a genuine suboptimal alignment, so score2/te2 are
    // exercised (random pairs rarely have a meaningful second-best).
    pairs.push_back(gen_tandem_repeat_pair(rng, 100, 150));
    pairs.push_back(gen_tandem_repeat_pair(rng, 101, 157));
    pairs.push_back(gen_tandem_repeat_pair(rng, 128, 200));
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

// RAII override of an environment variable: set on construction, restore the
// prior value (or unset if it was absent) on destruction, so doctest run order
// cannot leak an override into a sibling case that assumes the defaults.
class ScopedEnv {
public:
    ScopedEnv(const char *name, const char *value) : name_(name) {
        const char *cur = getenv(name);
        had_ = cur != nullptr;
        if (had_) saved_ = cur;
        setenv(name, value, 1);
    }
    ~ScopedEnv() {
        if (had_) setenv(name_.c_str(), saved_.c_str(), 1);
        else      unsetenv(name_.c_str());
    }
    ScopedEnv(const ScopedEnv &) = delete;
    ScopedEnv &operator=(const ScopedEnv &) = delete;
private:
    std::string name_;
    std::string saved_;
    bool had_ = false;
};

// Strict all-field kswr_t equality (every observable field). For A/B arms that
// run the same recurrence and so must agree exactly -- not the tolerant
// scalar-vs-batched compare in kswr_cmp.h.
inline bool kswr_all_fields_eq(const kswr_t &a, const kswr_t &b) {
    return a.score == b.score && a.te == b.te && a.qe == b.qe
        && a.score2 == b.score2 && a.te2 == b.te2
        && a.tb == b.tb && a.qb == b.qb;
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

TEST_CASE("kswv::getScores16 matches scalar ksw_align2 on 10k random + curated edge pairs"
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

    // use16 = true drives the 16-bit kernel (kswv256_16 / kswv512_16 /
    // kswv_neon_16) regardless of the production l_ms*a routing.
    auto batched = bwa_tests::run_kswv_batch(pairs, mat,
                                             bwa_tests::DEFAULT_GAP_OPEN,
                                             bwa_tests::DEFAULT_GAP_EXTEND,
                                             0, /*use16=*/true);
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
        if (!score_ok)  { ++score_mism;  CHECK(score_ok); }
        if (!coord_ok)  { ++coord_mism;  CHECK(coord_ok); }
        if (!score2_ok) { ++score2_mism; CHECK(score2_ok); }
    }

    MESSAGE("kswv16 vs scalar: score_mism=" << score_mism
            << " coord_mism=" << coord_mism
            << " score2_mism=" << score2_mism
            << " over " << pairs.size() << " pairs");
    CHECK(score_mism == 0);
    CHECK(coord_mism == 0);
    CHECK(score2_mism == 0);
}

TEST_CASE("kswv::getScores16 matches scalar on high scores that overflow the 8-bit kernel"
          * doctest::test_suite("unit/kswv")) {
    // match=14 with long, near-exact pairs pushes alignment scores well
    // past 255, where the 8-bit kernel saturates. The 16-bit kernel must
    // still reproduce scalar ksw_align2 exactly. match=14 also makes
    // min_seed_len*match = 19*14 = 266 >= 250, so default_xtra_flags drops
    // KSW_XBYTE for BOTH the scalar reference and the batch — i.e. both run
    // the word path, exactly the case production routes to getScores16.
    auto mat = bwa_tests::build_scoring_matrix(14, 8, 1);

    std::mt19937 rng(1234);
    std::vector<bwa_tests::TestPair> pairs;
    std::uniform_int_distribution<int> qlen_d(80, 128);
    std::uniform_int_distribution<int> rlen_d(150, 250);
    for (int i = 0; i < 5000; i++) {
        pairs.push_back(bwa_tests::gen_sub_cluster_pair(
            rng, qlen_d(rng), rlen_d(rng), 40, 2));
    }

    std::vector<kswr_t> scalar_aln;
    scalar_aln.reserve(pairs.size());
    for (const auto &p : pairs) {
        scalar_aln.push_back(bwa_tests::run_scalar_ksw(p, mat));
    }

    auto batched = bwa_tests::run_kswv_batch(pairs, mat,
                                             bwa_tests::DEFAULT_GAP_OPEN,
                                             bwa_tests::DEFAULT_GAP_EXTEND,
                                             0, /*use16=*/true);
    REQUIRE(batched.size() == pairs.size());

    int score_mism = 0, coord_mism = 0, score2_mism = 0, over255 = 0;
    for (size_t i = 0; i < pairs.size(); i++) {
        CAPTURE(i);
        CAPTURE(pairs[i].tag);
        CAPTURE(scalar_aln[i].score);
        CAPTURE(batched[i].score);
        if (scalar_aln[i].score > 255) ++over255;
        const bool score_ok  = bwa_tests::kswr_score_eq(scalar_aln[i], batched[i]);
        const bool coord_ok  = bwa_tests::kswr_coords_eq(scalar_aln[i], batched[i]);
        const bool score2_ok = bwa_tests::kswr_score2_eq(scalar_aln[i], batched[i]);
        if (!score_ok)  { ++score_mism;  CHECK(score_ok); }
        if (!coord_ok)  { ++coord_mism;  CHECK(coord_ok); }
        if (!score2_ok) { ++score2_mism; CHECK(score2_ok); }
    }

    MESSAGE("kswv16 high-score: over255=" << over255
            << " score_mism=" << score_mism
            << " coord_mism=" << coord_mism
            << " score2_mism=" << score2_mism
            << " over " << pairs.size() << " pairs");
    // Guard the test's own premise: at least some pairs must exceed the
    // 8-bit ceiling, else this wouldn't be testing the 16-bit range.
    CHECK(over255 > 0);
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

// The two-row (KSWV_NEON_U8_CELL_PAIR) and one-row (KSWV_NEON_U8_CELL) macros
// restate the same u8-rescue recurrence, so a change applied to one but not the
// other would diverge silently. rescue_rowpair_enabled() reads BWA3_RESCUE_ROWPAIR
// on every call (the NEON u8 kernel only), so this drives one batch through
// getScores8 with pairing off then on. Two gates:
//   * strict all-field A/B (one-row == two-row) catches drift between the macros;
//   * an INDEPENDENT scalar oracle on the one-row arm — the shared KSWV_U8_EPILOGUE
//     means a bug there shifts BOTH arms together, so the A/B alone would stay
//     green; the existing scalar-oracle cases run only the default (pairing-on)
//     mode, leaving the one-row path without a reference otherwise.
TEST_CASE("kswv u8 rescue: BWA3_RESCUE_ROWPAIR off == on, and one-row matches scalar"
          * doctest::test_suite("unit/kswv")) {

#if !defined(__ARM_NEON) && !defined(__aarch64__)
    // BWA3_RESCUE_ROWPAIR gates only the NEON u8 rescue kernel; on x86 tiers
    // run_kswv_batch dispatches to kswv256/512_u8 where both arms run identical
    // code, so this case would pass without exercising the row-pair logic.
    MESSAGE("skipped: BWA3_RESCUE_ROWPAIR affects the NEON u8 rescue kernel only");
    return;
#else
    auto mat = bwa_tests::build_scoring_matrix(1, 4, 1);

    std::mt19937 rng(1234);
    auto pairs = build_edge_cases(rng);
    // A few hundred bulk pairs is plenty: macro drift surfaces on the first
    // affected cell, and each batch sweeps the full DP twice, so keep well under
    // the ~100 ms per-case budget (the scalar-oracle cases above use 10k).
    auto bulk  = build_bulk_random(rng, 300);
    pairs.insert(pairs.end(), bulk.begin(), bulk.end());

    // Independent oracle: scalar ksw_align2 for every pair.
    std::vector<kswr_t> scalar_aln;
    scalar_aln.reserve(pairs.size());
    for (const auto &p : pairs) scalar_aln.push_back(bwa_tests::run_scalar_ksw(p, mat));

    // ScopedEnv save/restores the caller's env so doctest run order can't leak
    // the override into a sibling case that assumes the default (pairing on).
    std::vector<kswr_t> one_row, two_row;
    {
        ScopedEnv rp("BWA3_RESCUE_ROWPAIR", "0");
        one_row = bwa_tests::run_kswv_batch(pairs, mat);
    }
    {
        ScopedEnv rp("BWA3_RESCUE_ROWPAIR", "1");
        two_row = bwa_tests::run_kswv_batch(pairs, mat);
    }

    REQUIRE(one_row.size() == pairs.size());
    REQUIRE(two_row.size() == pairs.size());

    int drift = 0, oracle_mism = 0;
    for (size_t i = 0; i < pairs.size(); i++) {
        const kswr_t &a = one_row[i];
        const kswr_t &b = two_row[i];
        // (1) Strict all-field A/B: rowpair is the same kernel row-blocked, so the
        // two arms must match exactly -- not the scalar-vs-batched tolerant compare.
        const bool eq = kswr_all_fields_eq(a, b);
        if (!eq) {
            ++drift;
            CAPTURE(i); CAPTURE(pairs[i].tag);
            CAPTURE(a.score); CAPTURE(b.score);
            CAPTURE(a.te); CAPTURE(b.te); CAPTURE(a.qe); CAPTURE(b.qe);
            CHECK(eq);
        }
        // (2) Independent oracle on the one-row arm: a bug in the shared epilogue
        // shifts both arms together, so check one_row against scalar ksw_align2
        // (tolerant compare — scalar vs batched has known-benign discrepancies).
        const bool o_score  = bwa_tests::kswr_score_eq(scalar_aln[i], a);
        const bool o_coord  = bwa_tests::kswr_coords_eq(scalar_aln[i], a);
        const bool o_score2 = bwa_tests::kswr_score2_eq(scalar_aln[i], a);
        if (!(o_score && o_coord && o_score2)) {
            ++oracle_mism;
            CAPTURE(i); CAPTURE(pairs[i].tag);
            CAPTURE(scalar_aln[i].score); CAPTURE(a.score);
            CHECK(o_score); CHECK(o_coord); CHECK(o_score2);
        }
    }
    MESSAGE("rowpair off-vs-on drift=" << drift << ", one-row vs scalar mism="
            << oracle_mism << " over " << pairs.size() << " pairs");
    CHECK(drift == 0);
    CHECK(oracle_mism == 0);
#endif
}

// The two-row sweep recovers the query-end column either inline (a per-cell
// strict-greater blend of a carried column counter) or lazily (row i's H written
// back into the dead diagonal slot H0[j], both rows' running maxima checkpointed
// per QE_BLK, and the shared KSWV_U8_EPILOGUE rescan run for both rows).
// rescue_lazyqe_enabled() reads BWA3_RESCUE_LAZYQE on every call, so this drives
// one batch through the paired kernel with the lazy path off then on. Both arms
// run the same recurrence, so the A/B is strict all-field; the inline arm is
// additionally checked against the scalar oracle so a shared-epilogue bug cannot
// hide by shifting both arms together. Pairing is forced ON for both arms
// (BWA3_RESCUE_ROWPAIR=1) so the lazy toggle is the only difference.
TEST_CASE("kswv u8 rescue: BWA3_RESCUE_LAZYQE off == on in the two-row sweep, and inline matches scalar"
          * doctest::test_suite("unit/kswv")) {

#if !defined(__ARM_NEON) && !defined(__aarch64__)
    MESSAGE("skipped: BWA3_RESCUE_LAZYQE affects the NEON u8 rescue kernel only");
    return;
#else
    auto mat = bwa_tests::build_scoring_matrix(1, 4, 1);

    std::mt19937 rng(4321);
    auto pairs = build_edge_cases(rng);
    auto bulk  = build_bulk_random(rng, 300);
    pairs.insert(pairs.end(), bulk.begin(), bulk.end());

    std::vector<kswr_t> scalar_aln;
    scalar_aln.reserve(pairs.size());
    for (const auto &p : pairs) scalar_aln.push_back(bwa_tests::run_scalar_ksw(p, mat));

    // ScopedEnv save/restores both vars so doctest run order cannot leak an
    // override into a sibling case that assumes the defaults (pairing on, lazy
    // on). Pairing is held ON so the lazy toggle is the only difference.
    std::vector<kswr_t> inline_qe, lazy_qe;
    {
        ScopedEnv rp("BWA3_RESCUE_ROWPAIR", "1");
        {
            ScopedEnv lq("BWA3_RESCUE_LAZYQE", "0");
            inline_qe = bwa_tests::run_kswv_batch(pairs, mat);
        }
        {
            ScopedEnv lq("BWA3_RESCUE_LAZYQE", "1");
            lazy_qe = bwa_tests::run_kswv_batch(pairs, mat);
        }
    }

    REQUIRE(inline_qe.size() == pairs.size());
    REQUIRE(lazy_qe.size() == pairs.size());

    int drift = 0, oracle_mism = 0;
    for (size_t i = 0; i < pairs.size(); i++) {
        const kswr_t &a = inline_qe[i];
        const kswr_t &b = lazy_qe[i];
        const bool eq = kswr_all_fields_eq(a, b);
        if (!eq) {
            ++drift;
            CAPTURE(i); CAPTURE(pairs[i].tag);
            CAPTURE(a.score); CAPTURE(b.score);
            CAPTURE(a.te); CAPTURE(b.te); CAPTURE(a.qe); CAPTURE(b.qe);
            CAPTURE(a.tb); CAPTURE(b.tb); CAPTURE(a.qb); CAPTURE(b.qb);
            CHECK(eq);
        }
        const bool o_score  = bwa_tests::kswr_score_eq(scalar_aln[i], a);
        const bool o_coord  = bwa_tests::kswr_coords_eq(scalar_aln[i], a);
        const bool o_score2 = bwa_tests::kswr_score2_eq(scalar_aln[i], a);
        if (!(o_score && o_coord && o_score2)) {
            ++oracle_mism;
            CAPTURE(i); CAPTURE(pairs[i].tag);
            CAPTURE(scalar_aln[i].score); CAPTURE(a.score);
            CHECK(o_score); CHECK(o_coord); CHECK(o_score2);
        }
    }
    MESSAGE("lazyqe off-vs-on drift=" << drift << ", inline vs scalar mism="
            << oracle_mism << " over " << pairs.size() << " pairs");
    CHECK(drift == 0);
    CHECK(oracle_mism == 0);
#endif
}

// The NEON 16-bit rescue kernel carries the same two-row sweep and lazy
// query-end recovery as the u8 kernel, behind the same env toggles. Drive one
// batch through getScores16 (use16 = true) in each of the three configurations
// -- one-row, two-row inline, two-row lazy -- and require strict all-field
// equality between them, with the one-row arm checked against the scalar
// oracle. High-scoring pairs (long exact matches) make the 16-bit tier the
// production route, but the kernel is exercised directly here regardless.
TEST_CASE("kswv u16 rescue: ROWPAIR/LAZYQE configurations agree, and one-row matches scalar"
          * doctest::test_suite("unit/kswv")) {

#if !defined(__ARM_NEON) && !defined(__aarch64__)
    MESSAGE("skipped: the toggles affect the NEON 16-bit rescue kernel only");
    return;
#else
    auto mat = bwa_tests::build_scoring_matrix(1, 4, 1);

    std::mt19937 rng(8642);
    auto pairs = build_edge_cases(rng);
    auto bulk  = build_bulk_random(rng, 300);
    pairs.insert(pairs.end(), bulk.begin(), bulk.end());

    std::vector<kswr_t> scalar_aln;
    scalar_aln.reserve(pairs.size());
    for (const auto &p : pairs) scalar_aln.push_back(bwa_tests::run_scalar_ksw(p, mat));

    auto run16 = [&]() {
        return bwa_tests::run_kswv_batch(pairs, mat,
                                         bwa_tests::DEFAULT_GAP_OPEN,
                                         bwa_tests::DEFAULT_GAP_EXTEND,
                                         0, /*use16=*/true);
    };
    // ScopedEnv save/restores both vars across the three configurations.
    std::vector<kswr_t> one_row, pair_inline, pair_lazy;
    {
        ScopedEnv rp("BWA3_RESCUE_ROWPAIR", "0");
        ScopedEnv lq("BWA3_RESCUE_LAZYQE", "1");
        one_row = run16();
    }
    {
        ScopedEnv rp("BWA3_RESCUE_ROWPAIR", "1");
        {
            ScopedEnv lq("BWA3_RESCUE_LAZYQE", "0");
            pair_inline = run16();
        }
        {
            ScopedEnv lq("BWA3_RESCUE_LAZYQE", "1");
            pair_lazy = run16();
        }
    }

    REQUIRE(one_row.size() == pairs.size());
    REQUIRE(pair_inline.size() == pairs.size());
    REQUIRE(pair_lazy.size() == pairs.size());

    auto same = [](const kswr_t &a, const kswr_t &b) {
        return kswr_all_fields_eq(a, b);
    };
    int drift_inline = 0, drift_lazy = 0, oracle_mism = 0;
    for (size_t i = 0; i < pairs.size(); i++) {
        const kswr_t &a = one_row[i];
        if (!same(a, pair_inline[i])) {
            ++drift_inline;
            CAPTURE(i); CAPTURE(pairs[i].tag);
            CAPTURE(a.score); CAPTURE(pair_inline[i].score);
            CAPTURE(a.te); CAPTURE(pair_inline[i].te); CAPTURE(a.qe); CAPTURE(pair_inline[i].qe);
            CHECK(same(a, pair_inline[i]));
        }
        if (!same(a, pair_lazy[i])) {
            ++drift_lazy;
            CAPTURE(i); CAPTURE(pairs[i].tag);
            CAPTURE(a.score); CAPTURE(pair_lazy[i].score);
            CAPTURE(a.te); CAPTURE(pair_lazy[i].te); CAPTURE(a.qe); CAPTURE(pair_lazy[i].qe);
            CHECK(same(a, pair_lazy[i]));
        }
        const bool o_score  = bwa_tests::kswr_score_eq(scalar_aln[i], a);
        const bool o_coord  = bwa_tests::kswr_coords_eq(scalar_aln[i], a);
        const bool o_score2 = bwa_tests::kswr_score2_eq(scalar_aln[i], a);
        if (!(o_score && o_coord && o_score2)) {
            ++oracle_mism;
            CAPTURE(i); CAPTURE(pairs[i].tag);
            CAPTURE(scalar_aln[i].score); CAPTURE(a.score);
            CHECK(o_score); CHECK(o_coord); CHECK(o_score2);
        }
    }
    MESSAGE("u16 one-row vs pair-inline drift=" << drift_inline
            << ", vs pair-lazy drift=" << drift_lazy
            << ", one-row vs scalar mism=" << oracle_mism
            << " over " << pairs.size() << " pairs");
    CHECK(drift_inline == 0);
    CHECK(drift_lazy == 0);
    CHECK(oracle_mism == 0);
#endif
}

// The biased u8 body (BWA3_RESCUE_USQADD=0) is a distinct instantiation: it
// replaces the saturating add with add+bias+saturating-subtract and splits the
// column range at a data-dependent jdummy, whose partial final block can be an
// odd column count -- the exact path the two-column unroll's odd-tail branch
// (`if (j < jend_) { ...; d1 = d1b; }`) handles. No other case drives USQADD=0.
// Pairing is held ON so it is the two-row biased body (and its odd tail) that
// USQADD toggles; the odd query lengths in build_edge_cases make the tail odd.
// Strict all-field A/B against the default (saturating) arm, plus the biased arm
// against the independent scalar oracle.
TEST_CASE("kswv u8 rescue: BWA3_RESCUE_USQADD off == on, and biased body matches scalar"
          * doctest::test_suite("unit/kswv")) {
#if !defined(__ARM_NEON) && !defined(__aarch64__)
    MESSAGE("skipped: BWA3_RESCUE_USQADD affects the NEON u8 rescue kernel only");
    return;
#else
    auto mat = bwa_tests::build_scoring_matrix(1, 4, 1);

    std::mt19937 rng(2468);
    auto pairs = build_edge_cases(rng);
    auto bulk  = build_bulk_random(rng, 300);
    pairs.insert(pairs.end(), bulk.begin(), bulk.end());

    std::vector<kswr_t> scalar_aln;
    scalar_aln.reserve(pairs.size());
    for (const auto &p : pairs) scalar_aln.push_back(bwa_tests::run_scalar_ksw(p, mat));

    std::vector<kswr_t> biased, saturating;
    {
        ScopedEnv rp("BWA3_RESCUE_ROWPAIR", "1");
        {
            ScopedEnv uq("BWA3_RESCUE_USQADD", "0");
            biased = bwa_tests::run_kswv_batch(pairs, mat);
        }
        {
            ScopedEnv uq("BWA3_RESCUE_USQADD", "1");
            saturating = bwa_tests::run_kswv_batch(pairs, mat);
        }
    }

    REQUIRE(biased.size() == pairs.size());
    REQUIRE(saturating.size() == pairs.size());

    int drift = 0, oracle_mism = 0;
    for (size_t i = 0; i < pairs.size(); i++) {
        const kswr_t &a = biased[i];
        const kswr_t &b = saturating[i];
        if (!kswr_all_fields_eq(a, b)) {
            ++drift;
            CAPTURE(i); CAPTURE(pairs[i].tag);
            CAPTURE(a.score); CAPTURE(b.score);
            CAPTURE(a.te); CAPTURE(b.te); CAPTURE(a.qe); CAPTURE(b.qe);
            CHECK(kswr_all_fields_eq(a, b));
        }
        const bool o_score  = bwa_tests::kswr_score_eq(scalar_aln[i], a);
        const bool o_coord  = bwa_tests::kswr_coords_eq(scalar_aln[i], a);
        const bool o_score2 = bwa_tests::kswr_score2_eq(scalar_aln[i], a);
        if (!(o_score && o_coord && o_score2)) {
            ++oracle_mism;
            CAPTURE(i); CAPTURE(pairs[i].tag);
            CAPTURE(scalar_aln[i].score); CAPTURE(a.score);
            CHECK(o_score); CHECK(o_coord); CHECK(o_score2);
        }
    }
    MESSAGE("usqadd off-vs-on drift=" << drift << ", biased vs scalar mism="
            << oracle_mism << " over " << pairs.size() << " pairs");
    CHECK(drift == 0);
    CHECK(oracle_mism == 0);
#endif
}

// BWA3_RESCUE_FSCAN selects a separate instantiation of every u8 kernel body
// (NEON, AVX2, AVX-512BW) -- the G-based cell (F opens from G only, one
// sat(G - oe) shared by both gaps, row max over G), the query-only boundary
// mask, and no reference-pad mask. Each is argued byte-identical in kswv.cpp;
// this pins it. For each body the flag reaches -- on NEON one-row, two-row
// inline argmax, two-row lazy, and the biased (USQADD=0) two-row body; on x86
// the one body per tier, which the other toggles do not reach -- one batch runs
// with FSCAN off then on and must agree on every field. The FSCAN arm of the
// default configuration is also checked against the independent scalar
// oracle. build_edge_cases' ragged lengths put pad rows and pad columns in most
// lane groups.
//
// Two scoring regimes. The FSCAN cell drops the E->F transition and relies on
// F->E to cover it (a gap pair's twin in the other order). bwa's default
// (mismatch 4, gaps 6+1) never puts an insertion directly beside a deletion on
// an optimal path -- one mismatch (4) beats two gap opens (14) -- so it would
// stay green even if the cell lost BOTH transitions. Mismatch 20 with gaps
// 1+1 makes the adjacent gap pair (4) the cheap move, so it would not.
TEST_CASE("kswv u8 rescue: BWA3_RESCUE_FSCAN off == on in every u8 body, and FSCAN matches scalar"
          * doctest::test_suite("unit/kswv")) {
    std::mt19937 rng(97531);
    auto pairs = build_edge_cases(rng);
    auto bulk  = build_bulk_random(rng, 300);
    pairs.insert(pairs.end(), bulk.begin(), bulk.end());

    struct Regime { std::string name; int mismatch, gap_open, gap_extend; };
    const Regime regimes[] = {
        {"default (B4 O6 E1)",       4, 6, 1},
        {"adjacent gaps (B20 O1 E1)", 20, 1, 1},
    };
    struct Config { std::string name; const char *rowpair, *lazyqe, *usqadd; bool oracle; };
    const Config configs[] = {
#if defined(__ARM_NEON) || defined(__aarch64__)
        {"one-row",         "0", "1", "1", false},
        {"two-row inline",  "1", "0", "1", false},
        {"two-row lazy",    "1", "1", "1", true},   // the production default
        {"two-row biased",  "1", "1", "0", false},
#else
        {"x86",             "1", "1", "1", true},   // toggles other than FSCAN are NEON-only
#endif
    };
    int oracle_mism = 0;
    for (const Regime &r : regimes) {
        auto mat = bwa_tests::build_scoring_matrix(1, r.mismatch, 1);
        std::vector<kswr_t> scalar_aln;
        scalar_aln.reserve(pairs.size());
        for (const auto &p : pairs)
            scalar_aln.push_back(bwa_tests::run_scalar_ksw(p, mat, r.gap_open, r.gap_extend));
        for (const Config &c : configs) {
            ScopedEnv rp("BWA3_RESCUE_ROWPAIR", c.rowpair);
            ScopedEnv lq("BWA3_RESCUE_LAZYQE", c.lazyqe);
            ScopedEnv uq("BWA3_RESCUE_USQADD", c.usqadd);
            std::vector<kswr_t> off, on;
            {
                ScopedEnv fs("BWA3_RESCUE_FSCAN", "0");
                off = bwa_tests::run_kswv_batch(pairs, mat, r.gap_open, r.gap_extend);
            }
            {
                ScopedEnv fs("BWA3_RESCUE_FSCAN", "1");
                on = bwa_tests::run_kswv_batch(pairs, mat, r.gap_open, r.gap_extend);
            }
            REQUIRE(off.size() == pairs.size());
            REQUIRE(on.size() == pairs.size());

            int drift = 0;
            for (size_t i = 0; i < pairs.size(); i++) {
                const kswr_t &a = off[i];
                const kswr_t &b = on[i];
                if (!kswr_all_fields_eq(a, b)) {
                    ++drift;
                    CAPTURE(r.name); CAPTURE(c.name); CAPTURE(i); CAPTURE(pairs[i].tag);
                    CAPTURE(a.score); CAPTURE(b.score);
                    CAPTURE(a.te); CAPTURE(b.te); CAPTURE(a.qe); CAPTURE(b.qe);
                    CAPTURE(a.score2); CAPTURE(b.score2); CAPTURE(a.te2); CAPTURE(b.te2);
                    CAPTURE(a.tb); CAPTURE(b.tb); CAPTURE(a.qb); CAPTURE(b.qb);
                    CHECK(kswr_all_fields_eq(a, b));
                }
                if (c.oracle) {
                    const bool o_score  = bwa_tests::kswr_score_eq(scalar_aln[i], b);
                    const bool o_coord  = bwa_tests::kswr_coords_eq(scalar_aln[i], b);
                    const bool o_score2 = bwa_tests::kswr_score2_eq(scalar_aln[i], b);
                    if (!(o_score && o_coord && o_score2)) {
                        ++oracle_mism;
                        CAPTURE(r.name); CAPTURE(i); CAPTURE(pairs[i].tag);
                        CAPTURE(scalar_aln[i].score); CAPTURE(b.score);
                        CHECK(o_score); CHECK(o_coord); CHECK(o_score2);
                    }
                }
            }
            MESSAGE("fscan off-vs-on [" << r.name << "] (" << c.name << ") drift=" << drift
                    << " over " << pairs.size() << " pairs");
            CHECK(drift == 0);
        }
    }
    CHECK(oracle_mism == 0);
}

// The 16-bit twin of the case above: BWA3_RESCUE_FSCAN selects the G-based
// int16 cell in the NEON (one-row, two-row inline, two-row lazy), AVX2 and
// AVX-512BW 16-bit bodies. Same two scoring regimes and why, plus the regime
// the 16-bit kernel exists for -- match 14, scores well past 255 -- whose
// sub-cluster pairs put the gaps next to high-scoring matches.
TEST_CASE("kswv u16 rescue: BWA3_RESCUE_FSCAN off == on in every u16 body, and FSCAN matches scalar"
          * doctest::test_suite("unit/kswv")) {
    std::mt19937 rng(86420);
    auto pairs = build_edge_cases(rng);
    auto bulk  = build_bulk_random(rng, 300);
    pairs.insert(pairs.end(), bulk.begin(), bulk.end());
    std::vector<bwa_tests::TestPair> high;
    {
        std::uniform_int_distribution<int> qlen_d(80, 128);
        std::uniform_int_distribution<int> rlen_d(150, 250);
        for (int i = 0; i < 300; i++)
            high.push_back(bwa_tests::gen_sub_cluster_pair(rng, qlen_d(rng), rlen_d(rng), 40, 2));
    }

    struct Regime {
        std::string name; int match, mismatch, gap_open, gap_extend;
        const std::vector<bwa_tests::TestPair> *pairs;
    };
    const Regime regimes[] = {
        {"default (A1 B4 O6 E1)",       1,  4, 6, 1, &pairs},
        {"adjacent gaps (A1 B20 O1 E1)", 1, 20, 1, 1, &pairs},
        {"high score (A14 B8 O6 E1)",   14,  8, 6, 1, &high},
    };
    struct Config { std::string name; const char *rowpair, *lazyqe; bool oracle; };
    const Config configs[] = {
#if defined(__ARM_NEON) || defined(__aarch64__)
        {"one-row",         "0", "1", false},
        {"two-row inline",  "1", "0", false},
        {"two-row lazy",    "1", "1", true},   // the production default
#else
        {"x86",             "1", "1", true},   // toggles other than FSCAN are NEON-only
#endif
    };
    int oracle_mism = 0;
    for (const Regime &r : regimes) {
        const auto &ps = *r.pairs;
        auto mat = bwa_tests::build_scoring_matrix(r.match, r.mismatch, 1);
        std::vector<kswr_t> scalar_aln;
        scalar_aln.reserve(ps.size());
        for (const auto &p : ps)
            scalar_aln.push_back(bwa_tests::run_scalar_ksw(p, mat, r.gap_open, r.gap_extend));
        for (const Config &c : configs) {
            ScopedEnv rp("BWA3_RESCUE_ROWPAIR", c.rowpair);
            ScopedEnv lq("BWA3_RESCUE_LAZYQE", c.lazyqe);
            std::vector<kswr_t> off, on;
            {
                ScopedEnv fs("BWA3_RESCUE_FSCAN", "0");
                off = bwa_tests::run_kswv_batch(ps, mat, r.gap_open, r.gap_extend,
                                                0, /*use16=*/true);
            }
            {
                ScopedEnv fs("BWA3_RESCUE_FSCAN", "1");
                on = bwa_tests::run_kswv_batch(ps, mat, r.gap_open, r.gap_extend,
                                               0, /*use16=*/true);
            }
            REQUIRE(off.size() == ps.size());
            REQUIRE(on.size() == ps.size());

            int drift = 0;
            for (size_t i = 0; i < ps.size(); i++) {
                const kswr_t &a = off[i];
                const kswr_t &b = on[i];
                if (!kswr_all_fields_eq(a, b)) {
                    ++drift;
                    CAPTURE(r.name); CAPTURE(c.name); CAPTURE(i); CAPTURE(ps[i].tag);
                    CAPTURE(a.score); CAPTURE(b.score);
                    CAPTURE(a.te); CAPTURE(b.te); CAPTURE(a.qe); CAPTURE(b.qe);
                    CAPTURE(a.score2); CAPTURE(b.score2); CAPTURE(a.te2); CAPTURE(b.te2);
                    CAPTURE(a.tb); CAPTURE(b.tb); CAPTURE(a.qb); CAPTURE(b.qb);
                    CHECK(kswr_all_fields_eq(a, b));
                }
                if (c.oracle) {
                    const bool o_score  = bwa_tests::kswr_score_eq(scalar_aln[i], b);
                    const bool o_coord  = bwa_tests::kswr_coords_eq(scalar_aln[i], b);
                    const bool o_score2 = bwa_tests::kswr_score2_eq(scalar_aln[i], b);
                    if (!(o_score && o_coord && o_score2)) {
                        ++oracle_mism;
                        CAPTURE(r.name); CAPTURE(i); CAPTURE(ps[i].tag);
                        CAPTURE(scalar_aln[i].score); CAPTURE(b.score);
                        CHECK(o_score); CHECK(o_coord); CHECK(o_score2);
                    }
                }
            }
            MESSAGE("u16 fscan off-vs-on [" << r.name << "] (" << c.name << ") drift=" << drift
                    << " over " << ps.size() << " pairs");
            CHECK(drift == 0);
        }
    }
    CHECK(oracle_mism == 0);
}

// ---------------------------------------------------------------------------
// BWA3_RESCUE_FSCAN fuzz: directed batches, each built to steer the 11-op cell
// into one corner of its byte-identity argument (kswv.cpp: the FScan notes,
// KSWV_NEON_U8_CELL_PAIR_FS and the pad-row note above kswv_neon_u8). The two
// FSCAN cases above cover it on random and curated pairs; these pin the corners
// random pairs reach rarely, if at all.
// ---------------------------------------------------------------------------

namespace {

// Gap costs with the insertion and deletion sides set independently.
// run_kswv_batch / run_scalar_ksw tie them together, which hides the one
// FSCAN precondition with a nontrivial case: o + e equal, e different.
struct GapCosts { int o_del, e_del, o_ins, e_ins; };

const GapCosts kDefaultGaps  = {6, 1, 6, 1};
const GapCosts kAdjacentGaps = {1, 1, 1, 1};

std::vector<uint8_t> random_bases(std::mt19937 &rng, int n) {
    std::uniform_int_distribution<int> base(0, 3);
    std::vector<uint8_t> s(static_cast<size_t>(std::max(0, n)));
    for (auto &b : s) b = static_cast<uint8_t>(base(rng));
    return s;
}

int uniform(std::mt19937 &rng, int lo, int hi) {
    return std::uniform_int_distribution<int>(lo, hi)(rng);
}

void append(std::vector<uint8_t> &dst, const std::vector<uint8_t> &src) {
    dst.insert(dst.end(), src.begin(), src.end());
}

// Scalar ksw_align2 with independent gap costs (see run_scalar_ksw for the
// const_cast: ksw_align2 restores the prefixes it reverses under KSW_XSTART).
// `use16` clears KSW_XBYTE so the oracle is ksw_i16, the 16-bit kernels' twin:
// ksw_u8 and ksw_i16 pad the query to different quanta (16 / 8), and a pad
// column carries its diagonal on at score 0, so an alignment ending on the last
// query column lifts a different number of later row maxima in each -- which
// moves score2 (score, te, qe, tb and qb never read a pad column's value).
kswr_t run_scalar_gaps(const bwa_tests::TestPair &p, const bwa_tests::ScoringMatrix &mat,
                       const GapCosts &g, bool use16) {
    int xtra = bwa_tests::default_xtra_flags(static_cast<int>(mat[0]));
    if (use16) xtra &= ~KSW_XBYTE;
    return ksw_align2(static_cast<int>(p.qry.size()), const_cast<uint8_t *>(p.qry.data()),
                      static_cast<int>(p.ref.size()), const_cast<uint8_t *>(p.ref.data()),
                      5, mat.data(), g.o_del, g.e_del, g.o_ins, g.e_ins, xtra, nullptr);
}

// run_kswv_batch with independent gap costs: phase 0, then phase 1 on the
// survivors, through the compile-time-tier kswv. Results in input order.
std::vector<kswr_t> run_kswv_gaps(const std::vector<bwa_tests::TestPair> &pairs,
                                  const bwa_tests::ScoringMatrix &mat, const GapCosts &g,
                                  bool use16) {
    const int xtra = bwa_tests::default_xtra_flags(static_cast<int>(mat[0]));
    bwa_tests::BatchBuffers bb(pairs, xtra);
    int32_t max_ref = 0, max_qry = 0;
    for (const auto &p : pairs) {
        max_ref = std::max(max_ref, static_cast<int32_t>(p.ref.size()));
        max_qry = std::max(max_qry, static_cast<int32_t>(p.qry.size()));
    }
    kswv k(g.o_del, g.e_del, g.o_ins, g.e_ins, mat[0], mat[1], 1, max_ref, max_qry);
    if (use16) k.getScores16(bb.pairs(), bb.ref_buf(), bb.qer_buf(), bb.aln(), bb.n(), 1, 0);
    else       k.getScores8(bb.pairs(), bb.ref_buf(), bb.qer_buf(), bb.aln(), bb.n(), 1, 0);
    const int survivors = bb.prepare_phase1();
    if (use16) k.getScores16(bb.pairs(), bb.ref_buf(), bb.qer_buf(), bb.aln(), survivors, 1, 1);
    else       k.getScores8(bb.pairs(), bb.ref_buf(), bb.qer_buf(), bb.aln(), survivors, 1, 1);
    // aln() is indexed by regid == input index; prepare_phase1 compacts only pairs().
    return std::vector<kswr_t>(bb.aln(), bb.aln() + bb.n());
}

std::vector<kswr_t> scalar_all(const std::vector<bwa_tests::TestPair> &pairs,
                               const bwa_tests::ScoringMatrix &mat, const GapCosts &g,
                               bool use16) {
    std::vector<kswr_t> out;
    out.reserve(pairs.size());
    for (const auto &p : pairs) out.push_back(run_scalar_gaps(p, mat, g, use16));
    return out;
}

// FSCAN off vs on over `pairs` in every kernel body the flag reaches, strict on
// every field, plus the FSCAN arm of the production configuration against the
// scalar oracle of the same width. Mirrors the two FSCAN cases above,
// generalised to independent gap costs.
void check_fscan_off_on(const std::string &what, const std::vector<bwa_tests::TestPair> &pairs,
                        const bwa_tests::ScoringMatrix &mat, const GapCosts &g, bool use16) {
    // The dispatcher takes the FSCAN cell only when both gaps open at the same
    // o + e; anything else would silently compare the original body to itself.
    REQUIRE(g.o_ins + g.e_ins == g.o_del + g.e_del);
    const auto scalar = scalar_all(pairs, mat, g, use16);
    struct Config { const char *name, *rowpair, *lazyqe, *usqadd; bool oracle; };
    const Config configs[] = {
#if defined(__ARM_NEON) || defined(__aarch64__)
        {"one-row",        "0", "1", "1", false},
        {"two-row inline", "1", "0", "1", false},
        {"two-row lazy",   "1", "1", "1", true},   // the production default
        {"two-row biased", "1", "1", "0", false},  // u8 only: USQADD has no 16-bit body
#else
        {"x86",            "1", "1", "1", true},   // toggles other than FSCAN are NEON-only
#endif
    };
    int oracle_mism = 0;
    for (const Config &c : configs) {
        if (use16 && std::strcmp(c.usqadd, "0") == 0) continue;
        ScopedEnv rp("BWA3_RESCUE_ROWPAIR", c.rowpair);
        ScopedEnv lq("BWA3_RESCUE_LAZYQE", c.lazyqe);
        ScopedEnv uq("BWA3_RESCUE_USQADD", c.usqadd);
        std::vector<kswr_t> off, on;
        {
            ScopedEnv fs("BWA3_RESCUE_FSCAN", "0");
            off = run_kswv_gaps(pairs, mat, g, use16);
        }
        {
            ScopedEnv fs("BWA3_RESCUE_FSCAN", "1");
            on = run_kswv_gaps(pairs, mat, g, use16);
        }
        REQUIRE(off.size() == pairs.size());
        REQUIRE(on.size() == pairs.size());
        int drift = 0;
        for (size_t i = 0; i < pairs.size(); i++) {
            const kswr_t &a = off[i];
            const kswr_t &b = on[i];
            if (!kswr_all_fields_eq(a, b)) {
                ++drift;
                CAPTURE(what); CAPTURE(c.name); CAPTURE(i); CAPTURE(pairs[i].tag);
                CAPTURE(a.score); CAPTURE(b.score);
                CAPTURE(a.te); CAPTURE(b.te); CAPTURE(a.qe); CAPTURE(b.qe);
                CAPTURE(a.score2); CAPTURE(b.score2); CAPTURE(a.te2); CAPTURE(b.te2);
                CAPTURE(a.tb); CAPTURE(b.tb); CAPTURE(a.qb); CAPTURE(b.qb);
                CHECK(kswr_all_fields_eq(a, b));
            }
            if (c.oracle) {
                const kswr_t &s = scalar[i];
                const bool ends_ok = s.score <= 0 || (s.te == b.te && s.qe == b.qe);
                const bool ok = bwa_tests::kswr_score_eq(s, b) && ends_ok
                    && bwa_tests::kswr_coords_eq(s, b) && bwa_tests::kswr_score2_eq(s, b);
                if (!ok) {
                    ++oracle_mism;
                    CAPTURE(what); CAPTURE(i); CAPTURE(pairs[i].tag);
                    CAPTURE(s.score); CAPTURE(b.score); CAPTURE(s.te); CAPTURE(b.te);
                    CAPTURE(s.qe); CAPTURE(b.qe); CAPTURE(s.score2); CAPTURE(b.score2);
                    CAPTURE(s.tb); CAPTURE(b.tb); CAPTURE(s.qb); CAPTURE(b.qb);
                    CHECK(ok);
                }
            }
        }
        MESSAGE(what << " (" << std::string(c.name) << (use16 ? ", u16" : ", u8") << "): fscan off-vs-on drift="
                << drift << " over " << pairs.size() << " pairs");
        CHECK(drift == 0);
    }
    CHECK(oracle_mism == 0);
}

struct Regime { const char *name; int mismatch; GapCosts gaps; };

// bwa's default scoring, and the regime where an insertion next to a deletion
// is the cheap move (see the first FSCAN case for why both are needed).
const Regime kRegimes[] = {
    {"default (B4 O6 E1)",        4, kDefaultGaps},
    {"adjacent gaps (B20 O1 E1)", 20, kAdjacentGaps},
};

} // namespace

// Reference pad rows (i >= len1 in a lane whose group runs longer) lose the
// per-cell m11 zeroing under FSCAN; the argument is that a pad row's max never
// EXCEEDS the previous row's, so every strict comparison that reads the row max
// (gmax/te/qe, the lagged rowMax zeroing feeding score2, freeze/done) is
// unchanged. The boundary of that argument is a TIE: an alignment ending on
// the last real row at a query column short of the last carries its score
// diagonally into the first pad row unchanged (0xFF scores 0), so that pad row
// holds exactly the previous row's max. Two layouts, ragged reference lengths
// so every lane group has pad rows:
//   tie-end:    the best alignment ends on the last row (te, qe must not move);
//   tie-score2: a secondary ends on the last row, outside the best one's zone,
//               so score2 reads that row's max through the lagged zeroing.
TEST_CASE("kswv FSCAN: a pad row tying the last real row's max changes nothing"
          * doctest::test_suite("unit/kswv")) {
    std::mt19937 rng(424242);
    std::vector<bwa_tests::TestPair> pairs;
    std::vector<int> want_s2;   // tie-score2 pairs: the secondary's score, else -1
    for (int i = 0; i < 400; i++) {
        const int qlen = uniform(rng, 60, 120);
        const auto q = random_bases(rng, qlen);
        bwa_tests::TestPair p;
        p.qry = q;
        if (i % 2 == 0) {
            const int k = uniform(rng, qlen / 2, qlen - 5);   // ends at query column k - 1 < qlen - 1
            p.ref = random_bases(rng, uniform(rng, 20, 200));
            p.ref.insert(p.ref.end(), q.begin(), q.begin() + k);
            p.tag = "tie-end";
            want_s2.push_back(-1);
        } else {
            const int k = uniform(rng, 25, qlen - 10);
            p.ref = random_bases(rng, uniform(rng, 5, 40));
            append(p.ref, q);
            append(p.ref, random_bases(rng, uniform(rng, qlen + 5, qlen + 60)));
            p.ref.insert(p.ref.end(), q.begin(), q.begin() + k);
            p.tag = "tie-score2";
            want_s2.push_back(k);
        }
        pairs.push_back(p);
    }
    for (const Regime &r : kRegimes) {
        const auto mat = bwa_tests::build_scoring_matrix(1, r.mismatch, 1);
        const auto scalar = scalar_all(pairs, mat, r.gaps, false);
        // Non-vacuity: the layouts must actually produce the tie they are for.
        int tie_end = 0, tie_s2 = 0;
        for (size_t i = 0; i < pairs.size(); i++) {
            const kswr_t &s = scalar[i];
            if (want_s2[i] < 0 && s.te == static_cast<int>(pairs[i].ref.size()) - 1
                && s.qe < static_cast<int>(pairs[i].qry.size()) - 1) ++tie_end;
            if (want_s2[i] >= 0 && s.score2 == want_s2[i]) ++tie_s2;
        }
        MESSAGE(std::string(r.name) << ": tie-end " << tie_end << ", tie-score2 " << tie_s2);
        REQUIRE(tie_end > 100);
        REQUIRE(tie_s2 > 100);
        for (bool use16 : {false, true})
            check_fscan_off_on(std::string("pad-row tie, ") + r.name, pairs, mat, r.gaps, use16);
    }
}

// Phase 1 runs on the reversed prefixes ref[te..0] x q[qe..0], so its len1 =
// te + 1 and len2 = qe + 1 are ragged within a group whatever the phase-0
// layout was: pad rows AND pad columns, both under the dropped reference mask.
// Alignments starting at reference position 0 reach S on phase 1's LAST real
// row (tb = 0), the row right before the first pad row; a start past query
// column 0 (qb > 0) makes the phase-1 query ragged too. Reference tails of
// 0-250 bases spread te, hence the phase-1 row counts, across each group.
TEST_CASE("kswv FSCAN: phase 1 with ragged pad rows reaching S on the last real row"
          * doctest::test_suite("unit/kswv")) {
    std::mt19937 rng(515151);
    std::vector<bwa_tests::TestPair> pairs;
    for (int i = 0; i < 400; i++) {
        const int qlen = uniform(rng, 40, 150);
        const auto q = random_bases(rng, qlen);
        const int qb = uniform(rng, 0, qlen / 3);
        bwa_tests::TestPair p;
        p.qry = q;
        p.ref.assign(q.begin() + qb, q.end());
        // A mismatch or two strictly inside, so the start is still column qb
        // and phase 1 is not a pure diagonal.
        for (int m = uniform(rng, 0, 2); m > 0; m--) {
            const int at = uniform(rng, 3, static_cast<int>(p.ref.size()) - 4);
            p.ref[at] = static_cast<uint8_t>((p.ref[at] + 1) & 3);
        }
        append(p.ref, random_bases(rng, uniform(rng, 0, 250)));
        p.tag = "p1-start-at-0";
        pairs.push_back(p);
    }
    for (const Regime &r : kRegimes) {
        const auto mat = bwa_tests::build_scoring_matrix(1, r.mismatch, 1);
        const auto scalar = scalar_all(pairs, mat, r.gaps, false);
        int tb0 = 0, tb0_qb = 0;
        for (size_t i = 0; i < pairs.size(); i++) {
            if (scalar[i].score >= bwa_tests::DEFAULT_MIN_SEED_LEN && scalar[i].tb == 0) {
                ++tb0;
                if (scalar[i].qb > 0) ++tb0_qb;
            }
        }
        MESSAGE(std::string(r.name) << ": tb == 0 on " << tb0 << " pairs (" << tb0_qb << " with qb > 0)");
        REQUIRE(tb0 > 200);
        REQUIRE(tb0_qb > 100);
        for (bool use16 : {false, true})
            check_fscan_off_on(std::string("phase-1 pad rows, ") + r.name, pairs, mat, r.gaps,
                               use16);
    }
}

// FSCAN drops the E -> F transition and shares one sat(G - oe) between the
// gaps; both arguments lean on the saturating clamp, so they are tested at the
// two ends of the u8 range:
//   near 254: a mate of the longest length the u8 tier admits (len * a + shift
//             <= 254, shift = the mismatch cost: 250 at B4, 234 at B20), exact
//             but for alternating insertions and deletions, some adjacent, in
//             its last 30 columns -- E and F opening from H in the 220s-240s;
//   near 0:   short matches (5-9 bases; 2-6 at 1+1 gaps) separated by
//             alternating 1-2 base insertions and deletions, and 20-base matches
//             across 25-40 base deletions, so H hovers just above 0 and E / F
//             clamp at 0 mid-gap.
TEST_CASE("kswv FSCAN: alternating E / F near 0 and near the u8 ceiling"
          * doctest::test_suite("unit/kswv")) {
    std::mt19937 rng(606060);
    for (const Regime &r : kRegimes) {
        const int max_qlen = 254 - r.mismatch;   // the u8 admission bound at match 1
        // Near-0 match runs: 1+1 gaps are cheap enough that 5-9 base runs climb
        // well clear of 0, so the adjacent regime uses 2-6.
        const int run_lo = r.gaps.o_del + r.gaps.e_del > 2 ? 5 : 2;
        std::vector<bwa_tests::TestPair> high, low;
        for (int i = 0; i < 200; i++) {
            const int qlen = uniform(rng, max_qlen - 6, max_qlen);
            const auto q = random_bases(rng, qlen);
            bwa_tests::TestPair p;
            p.qry = q;
            p.ref = random_bases(rng, uniform(rng, 0, 60));
            const int zone = qlen - 30;
            p.ref.insert(p.ref.end(), q.begin(), q.begin() + zone);
            bool drop_query = uniform(rng, 0, 1) != 0;
            for (int j = zone; j < qlen;) {
                const int run = uniform(rng, 0, 4);   // 0: the next gap is adjacent
                for (int t = 0; t < run && j < qlen; t++) p.ref.push_back(q[j++]);
                if (drop_query) j++;                  // query base with no reference base
                else append(p.ref, random_bases(rng, 1));
                drop_query = !drop_query;
            }
            append(p.ref, random_bases(rng, uniform(rng, 0, 60)));
            p.tag = "alt-gaps-high";
            high.push_back(p);
        }
        for (int i = 0; i < 300; i++) {
            bwa_tests::TestPair p;
            p.qry = random_bases(rng, uniform(rng, 30, 90));
            p.ref = random_bases(rng, uniform(rng, 0, 30));
            const int qlen = static_cast<int>(p.qry.size());
            if (i % 3 == 0 && qlen >= 40) {
                // 20 matched, a 25-40 base deletion, then the rest matched.
                p.ref.insert(p.ref.end(), p.qry.begin(), p.qry.begin() + 20);
                append(p.ref, random_bases(rng, uniform(rng, 25, 40)));
                p.ref.insert(p.ref.end(), p.qry.begin() + 20, p.qry.end());
                p.tag = "long-gap-low";
            } else {
                bool drop_query = uniform(rng, 0, 1) != 0;
                for (int j = 0; j < qlen;) {
                    const int run = uniform(rng, run_lo, run_lo + 4);
                    for (int t = 0; t < run && j < qlen; t++) p.ref.push_back(p.qry[j++]);
                    if (drop_query) j += uniform(rng, 1, 2);
                    else append(p.ref, random_bases(rng, uniform(rng, 1, 2)));
                    drop_query = !drop_query;
                    if (uniform(rng, 0, 5) == 0 && !p.ref.empty())
                        p.ref.back() = static_cast<uint8_t>((p.ref.back() + 1) & 3);
                }
                p.tag = "alt-gaps-low";
            }
            append(p.ref, random_bases(rng, uniform(rng, 0, 30)));
            low.push_back(p);
        }
        const auto mat = bwa_tests::build_scoring_matrix(1, r.mismatch, 1);
        const auto s_high = scalar_all(high, mat, r.gaps, false);
        const auto s_low = scalar_all(low, mat, r.gaps, false);
        int n_high = 0, n_low = 0;
        for (const kswr_t &s : s_high) n_high += s.score >= 200;
        for (const kswr_t &s : s_low) n_low += s.score > 0 && s.score < 20;
        MESSAGE(std::string(r.name) << ": score >= 200 on " << n_high << "/" << high.size()
                << ", 0 < score < 20 on " << n_low << "/" << low.size());
        REQUIRE(n_high > 100);
        REQUIRE(n_low > 100);
        for (bool use16 : {false, true}) {
            check_fscan_off_on(std::string("alternating gaps near 254, ") + r.name, high, mat,
                               r.gaps, use16);
            check_fscan_off_on(std::string("alternating gaps near 0, ") + r.name, low, mat, r.gaps,
                               use16);
        }
    }
}

// The dispatcher's FSCAN condition compares o + e, not (o, e): one T = sat(G -
// oe) serves both gaps, and each gap still extends at its own e. So a scoring
// with o_ins + e_ins == o_del + e_del but e_ins != e_del takes the FSCAN cell,
// and the twin-path argument (an E -> F path's vertical-first twin costs the
// same) must hold with the two extension costs different. Multi-base
// insertions and deletions make the extension cost matter; the non-vacuity
// check proves it does by scoring the batch with the two sides swapped.
TEST_CASE("kswv FSCAN: equal gap-open totals with different extension costs"
          * doctest::test_suite("unit/kswv")) {
    std::mt19937 rng(717171);
    std::vector<bwa_tests::TestPair> pairs;
    for (int i = 0; i < 400; i++) {
        const int qlen = uniform(rng, 60, 140);
        const auto q = random_bases(rng, qlen);
        bwa_tests::TestPair p;
        p.qry = q;
        p.ref = random_bases(rng, uniform(rng, 0, 80));
        int j = 0;
        for (int ev = uniform(rng, 2, 4); ev > 0 && j < qlen; ev--) {
            const int run = uniform(rng, 8, 30);
            for (int t = 0; t < run && j < qlen; t++) p.ref.push_back(q[j++]);
            const int len = uniform(rng, 2, 6);
            if (uniform(rng, 0, 1)) j += len;   // insertion: query bases with no reference
            else append(p.ref, random_bases(rng, len));
        }
        if (j < qlen) p.ref.insert(p.ref.end(), q.begin() + j, q.end());
        append(p.ref, random_bases(rng, uniform(rng, 0, 80)));
        p.tag = "multi-base indels";
        pairs.push_back(p);
    }
    struct Case { const char *name; int mismatch; GapCosts gaps, swapped; };
    const Case cases[] = {
        {"B4 del 6+1 / ins 5+2",  4, {6, 1, 5, 2}, {5, 2, 6, 1}},
        {"B4 del 5+2 / ins 6+1",  4, {5, 2, 6, 1}, {6, 1, 5, 2}},
        {"B4 del 4+3 / ins 6+1",  4, {4, 3, 6, 1}, {6, 1, 4, 3}},
        {"B20 del 2+1 / ins 1+2", 20, {2, 1, 1, 2}, {1, 2, 2, 1}},
    };
    for (const Case &c : cases) {
        const auto mat = bwa_tests::build_scoring_matrix(1, c.mismatch, 1);
        const auto scalar = scalar_all(pairs, mat, c.gaps, false);
        const auto swapped = scalar_all(pairs, mat, c.swapped, false);
        int differ = 0;
        for (size_t i = 0; i < pairs.size(); i++) differ += scalar[i].score != swapped[i].score;
        MESSAGE(std::string(c.name) << ": score changes with the gap sides swapped on " << differ << "/"
                << pairs.size() << " pairs");
        REQUIRE(differ > 50);
        for (bool use16 : {false, true})
            check_fscan_off_on(c.name, pairs, mat, c.gaps, use16);
    }
}

#endif // BWA_TESTS_HAVE_KSWV
