// test/unit/test_kswv_correctness.cpp
//
// Self-consistency test for kswv::getScores8 versus scalar ksw_align2.
// Ported from test/kswv_selftest.cpp, now using the bwa_tests framework.
//
// Locks in bit-identical batched-SIMD output for 10,000 random pairs plus
// curated edge cases. Runs on every CI matrix row (SSE4.1, AVX2, AVX2
// clang, AVX2 no-mimalloc, multi-arch, ARM64 Linux, macOS ARM64).

#include <cstdlib>
#include <random>
#include <string>
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

// The dispatchers route to FScan only when the insertion and deletion
// open-plus-extend sums (oe = o + e) are equal and fit the lane, and no o or e
// is negative (fscan_scoring_ok); otherwise the original bodies run. Pin both
// sides of that gate against the scalar oracle, in every u8 and u16 body the
// build reaches, with separate deletion and insertion costs, which the
// symmetric runs above never reach:
//   - unequal sums (6+1 vs 5+1, and 1+1 vs 2+1 with adjacent gaps cheap):
//     must go to the original bodies, since FScan's shared sat(G - oe) would
//     open one of the gaps at the wrong cost;
//   - equal sums with unequal extends and opens (5+2 vs 6+1): goes through
//     FScan, whose proof needs only oe >= e for each gap;
//   - equal sums with a negative open (-1+8 vs 6+1): oe < e, so the original
//     bodies again.
TEST_CASE("kswv rescue: the FSCAN gate keeps asymmetric gaps exact against scalar"
          * doctest::test_suite("unit/kswv")) {
    std::mt19937 rng(424242);
    auto pairs = build_edge_cases(rng);
    auto bulk  = build_bulk_random(rng, 300);
    pairs.insert(pairs.end(), bulk.begin(), bulk.end());

    struct Gaps { std::string name; int mismatch, o_del, e_del, o_ins, e_ins; };
    const Gaps regimes[] = {
        {"open sums differ (del 6+1, ins 5+1)",      4, 6, 1, 5, 1},
        {"open sums differ (del 5+1, ins 6+1)",      4, 5, 1, 6, 1},
        {"open sums differ, adjacent gaps cheap",   20, 1, 1, 2, 1},
        {"equal sums, extends differ (5+2 vs 6+1)",  4, 5, 2, 6, 1},
        // oe < e on the deletion side: FScan's proof fails, so the gate must
        // route to the original bodies despite the equal sums.
        {"equal sums, negative open (-1+8 vs 6+1)",  4, -1, 8, 6, 1},
    };
    for (const Gaps &g : regimes) {
        auto mat = bwa_tests::build_scoring_matrix(1, g.mismatch, 1);
        for (const bool use16 : {false, true}) {
            // The oracle of the same width: the 16-bit kernels need ksw_i16
            // (KSW_XBYTE clear), since ksw_u8 pads the query to a different
            // quantum and can report a different score2.
            int xtra = bwa_tests::default_xtra_flags(1);
            if (use16) xtra &= ~KSW_XBYTE;
            std::vector<kswr_t> ref;
            ref.reserve(pairs.size());
            for (const auto &p : pairs)
                ref.push_back(bwa_tests::run_scalar_ksw_gaps(p, mat, g.o_del, g.e_del, g.o_ins,
                                                             g.e_ins, xtra));
            for (const std::string fscan : {"0", "1"}) {
                ScopedEnv fs("BWA3_RESCUE_FSCAN", fscan.c_str());
                const auto batch = bwa_tests::run_kswv_batch_gaps(pairs, mat, g.o_del, g.e_del,
                                                                  g.o_ins, g.e_ins, 0, use16);
                REQUIRE(batch.size() == pairs.size());
                int mism = 0;
                for (size_t i = 0; i < pairs.size(); i++) {
                    const bool ok = bwa_tests::kswr_score_eq(ref[i], batch[i]) &&
                                    bwa_tests::kswr_coords_eq(ref[i], batch[i]) &&
                                    bwa_tests::kswr_score2_eq(ref[i], batch[i]);
                    if (!ok) {
                        ++mism;
                        CAPTURE(g.name); CAPTURE(use16); CAPTURE(fscan); CAPTURE(i);
                        CAPTURE(pairs[i].tag);
                        CAPTURE(ref[i].score); CAPTURE(batch[i].score);
                        CAPTURE(ref[i].te); CAPTURE(batch[i].te);
                        CAPTURE(ref[i].score2); CAPTURE(batch[i].score2);
                        CHECK(ok);
                    }
                }
                MESSAGE("fscan gate [" << g.name << "] " << std::string(use16 ? "u16" : "u8")
                        << " BWA3_RESCUE_FSCAN=" << fscan << " mismatches=" << mism
                        << " over " << pairs.size() << " pairs");
                CHECK(mism == 0);
            }
        }
    }
}

#endif // BWA_TESTS_HAVE_KSWV
