// test/unit/test_kswv_rowpair_u8.cpp
//
// The u8 kswv FScan body sweeps two target rows per column (BWA3_RESCUE_ROWPAIR,
// default on) on NEON and on AVX-512BW; BWA3_RESCUE_ROWPAIR=0 restores the
// one-row sweep. The two sweeps must give every output field bit for bit.
//
// The kernel is reached through make_kswv (bwa_tests::run_kswv_batch_runtime),
// which dispatches on the RUNTIME tier and honors BWAMEM3_FORCE_TIER, so this
// reaches the AVX-512BW kernel on a capable host (run_kswv_batch binds the tier this translation unit is compiled
// for, AVX2 on x86, and would not). On AVX2 both arms run the same one-row
// kernel, so the case would pass without exercising the pair: it is registered
// as skipped there (and on the SSE tiers, which have no batched kswv kernel).
//
// Shapes the two-row sweep treats specially:
//   * an odd last row, which takes the one-row sweep (odd and even maxLen1);
//   * partial lane groups (len1 = len2 = 0 dummy lanes, jsplit = 0);
//   * N (4) query and reference bases, and query code 5 (what a '-' read base
//     encoded to before it became N; the kernels still score it from the
//     table), next to the pad blend;
//   * the signed-domain body (BWA3_RESCUE_USQADD on, the x86 default) and the
//     biased one (off);
//   * the batch exit landing on either row of a pair. A batch exits only once
//     every lane of a full group has reached its KSW_XSTOP score, so
//     exit_batch() builds full groups of planted exact matches whose last stop
//     falls on an even row (row i of a pair) in one group and an odd row (row
//     i+1) in the other, so both of the pair's exit branches run (KSW_XSUBO is
//     set too, so the score2 scan runs over the rows processed up to the exit,
//     and KSW_XSTART, so phase 1 runs from the stopped ends).
//     Every lane has frozen by then, so the exit row itself never reaches the
//     output and an exit that never fired would give the same fields. Nor can
//     it observe the pair's exit bookkeeping (limit, the i += 1 / i += 2 step,
//     processed_rows): rows past a frozen batch's exit hold rowMax 0, so a slip
//     there changes no output; that bookkeeping is checked by reading the code.
//     The case checks that both exit branches give the one-row fields, and that
//     both arms match scalar ksw_u8 (ksw_align2 with KSW_XBYTE, through the
//     tolerant comparators of kswr_cmp.h) under the same stop flags, so an
//     epilogue regression shared by the two arms is caught too.
// Both phases run (phase 1 recovers tb/qb from the reversed prefixes).

#include <cstdlib>
#include <random>
#include <string>
#include <vector>

#include "doctest/doctest.h"

#include "ksw_runner.h"
#include "kswr_cmp.h"
#include "scoped_env.h"
#include "kswv_runner.h"  // BWA_TESTS_HAVE_KSWV
#include "scoring.h"
#include "seqpair_gen.h"

#if BWA_TESTS_HAVE_KSWV

#include "ksw.h"
#include "kswv.h"
#include "simd_dispatch.h"  // bwamem3_simd_tier

namespace {

constexpr int kMaxRefLen = 461;
constexpr int kMaxQryLen = 151;
constexpr uint8_t kDash  = 5;   // the former '-' code; still a table-scored query byte

// The tiers whose u8 FScan body sweeps two rows: NEON and AVX-512BW. AVX2 has
// a batched kswv kernel but only the one-row sweep; the SSE classes exit().
// Read when the test case registers (before main), so a BWAMEM3_FORCE_TIER
// override must already be in the environment.
bool rowpair_tier() {
    bwamem3_simd_init();
    const int tier = bwamem3_simd_tier();
    return tier == BWAMEM3_TIER_AVX512BW || tier == BWAMEM3_TIER_NEON;
}

using bwa_tests::ScopedEnv;   // scoped_env.h

// Phase 0 then phase 1 on the u8 kernel at the runtime tier; returns the
// kswr_t of every pair, indexed by pair.
std::vector<kswr_t> run_u8(const std::vector<bwa_tests::TestPair> &pairs,
                           const bwa_tests::ScoringMatrix &mat, int xtra) {
    auto out = bwa_tests::run_kswv_batch_runtime(pairs, mat, xtra, kMaxRefLen, kMaxQryLen);
    REQUIRE(out.size() == pairs.size());   // empty when the kernel needs scalar
    return out;
}

// Sprinkle N into both sequences and code 5 into the query, about 2% each.
void add_ambiguous(std::mt19937 &rng, bwa_tests::TestPair &p) {
    std::uniform_int_distribution<int> pct(0, 99);
    for (auto &b : p.ref) if (pct(rng) < 2) b = 4;
    for (auto &b : p.qry) {
        const int r = pct(rng);
        if (r < 2) b = 4;
        else if (r < 4) b = kDash;
    }
}

// Tandem-repeat pairs (a high primary plus near-threshold secondaries, so
// score2 is live) with ragged query lengths, ragged reference lengths whose
// longest is `max_rlen` (so maxLen1 is odd or even on demand), and an odd pair
// count so the last group is partial.
std::vector<bwa_tests::TestPair> ragged_batch(std::mt19937 &rng, int n_pairs, int max_rlen) {
    std::uniform_int_distribution<int> rlen_d(200, max_rlen - 1);
    std::vector<bwa_tests::TestPair> pairs;
    for (int i = 0; i < n_pairs; i++) {
        const int qlen = (i % 3 == 0) ? 75 : (i % 3 == 1) ? 143 : 151;
        const int rlen = (i == n_pairs / 2) ? max_rlen : rlen_d(rng);
        pairs.push_back(bwa_tests::gen_tandem_repeat_pair(rng, qlen, rlen));
        add_ambiguous(rng, pairs.back());
    }
    return pairs;
}

constexpr int kStopQlen = 100;
constexpr int kStopRlen = 300;

// Two full 64-lane groups. Each lane's query is planted exactly in a random
// reference at offset o, so it first reaches its stop score kStopQlen (match
// score 1) on row o + kStopQlen - 1. The last lane of each group sets that
// group's exit row: 280 (even, row i of a pair) in group 0 and 279 (odd, row
// i+1) in group 1, both before the last row; every other lane stops by row 249.
std::vector<bwa_tests::TestPair> exit_batch(std::mt19937 &rng) {
    std::uniform_int_distribution<int> base(0, 3), off(0, 150);
    std::vector<bwa_tests::TestPair> pairs;
    for (int g = 0; g < 2; g++) {
        for (int l = 0; l < 64; l++) {
            bwa_tests::TestPair p;
            p.ref.resize(kStopRlen);
            for (auto &b : p.ref) b = static_cast<uint8_t>(base(rng));
            p.qry.resize(kStopQlen);
            for (auto &b : p.qry) b = static_cast<uint8_t>(base(rng));
            const int o = (l == 63) ? (g == 0 ? 181 : 180) : off(rng);
            for (int k = 0; k < kStopQlen; k++) p.ref[o + k] = p.qry[k];
            p.tag = "exit g" + std::to_string(g) + " l" + std::to_string(l);
            pairs.push_back(std::move(p));
        }
    }
    return pairs;
}

// Compare `arm` with `ref` pair by pair, capturing every field on a mismatch.
// strict: every field bit for bit (kswr_all_fields_eq), for the one-row vs
// two-row A/B, whose arms run the same recurrence. Otherwise `ref` is scalar
// ksw_align2 and the tolerant scalar-vs-batched comparators of kswr_cmp.h apply
// (score, tb/qb, score2, te/qe), as in test_kswv_correctness.cpp's oracle.
int compare_results(const std::vector<kswr_t> &ref, const std::vector<kswr_t> &arm,
                    const std::vector<bwa_tests::TestPair> &pairs, const std::string &what,
                    bool strict) {
    REQUIRE(ref.size() == arm.size());
    int mism = 0;
    for (size_t i = 0; i < arm.size(); i++) {
        const bool eq = strict
            ? bwa_tests::kswr_all_fields_eq(ref[i], arm[i])
            : bwa_tests::kswr_score_eq(ref[i], arm[i]) && bwa_tests::kswr_coords_eq(ref[i], arm[i])
              && bwa_tests::kswr_score2_eq(ref[i], arm[i]) && bwa_tests::kswr_ends_eq(ref[i], arm[i]);
        if (!eq) {
            ++mism;
            CAPTURE(what); CAPTURE(i); CAPTURE(pairs[i].tag);
            CAPTURE(ref[i].score); CAPTURE(arm[i].score);
            CAPTURE(ref[i].te); CAPTURE(arm[i].te);
            CAPTURE(ref[i].qe); CAPTURE(arm[i].qe);
            CAPTURE(ref[i].score2); CAPTURE(arm[i].score2);
            CAPTURE(ref[i].te2); CAPTURE(arm[i].te2);
            CAPTURE(ref[i].tb); CAPTURE(arm[i].tb);
            CAPTURE(ref[i].qb); CAPTURE(arm[i].qb);
            CHECK(eq);
        }
    }
    MESSAGE(what << std::string(strict ? ": one-row vs two-row" : ": scalar vs kswv") << " mismatches="
            << mism << " over " << arm.size() << " pairs (tier " << bwamem3_simd_tier() << ")");
    return mism;
}

// One-row vs two-row on every field; with an oracle, each arm against it too.
void check_rowpair(const std::vector<bwa_tests::TestPair> &pairs, int xtra, const char *what,
                   const std::vector<kswr_t> *oracle = nullptr) {
    const auto mat = bwa_tests::default_scoring_matrix();
    // The AVX-512BW two-row sweep lives in the FScan body, so an inherited
    // BWA3_RESCUE_FSCAN=0 would run the one-row body in both arms.
    ScopedEnv fs("BWA3_RESCUE_FSCAN", "1");
    for (const char *usq : {"1", "0"}) {
        ScopedEnv uq("BWA3_RESCUE_USQADD", usq);
        std::vector<kswr_t> one_row, two_row;
        {
            ScopedEnv rp("BWA3_RESCUE_ROWPAIR", "0");
            one_row = run_u8(pairs, mat, xtra);
        }
        {
            ScopedEnv rp("BWA3_RESCUE_ROWPAIR", "1");
            two_row = run_u8(pairs, mat, xtra);
        }
        const std::string label = std::string(what) + " usqadd=" + usq;
        CHECK(compare_results(one_row, two_row, pairs, label, /*strict=*/true) == 0);
        if (oracle) {
            CHECK(compare_results(*oracle, one_row, pairs, label + " one-row", false) == 0);
            CHECK(compare_results(*oracle, two_row, pairs, label + " two-row", false) == 0);
        }
    }
}

} // namespace

TEST_CASE("kswv u8 rescue at the runtime tier: BWA3_RESCUE_ROWPAIR off == on"
          * doctest::test_suite("unit/kswv")
          * doctest::skip(!rowpair_tier())) {
    const int xtra = bwa_tests::default_xtra_flags(1);

    SUBCASE("ragged, even maxLen1") {
        std::mt19937 rng(20261026u);
        check_rowpair(ragged_batch(rng, 129, 460), xtra, "ragged even");
    }
    SUBCASE("ragged, odd maxLen1") {
        std::mt19937 rng(20261027u);
        check_rowpair(ragged_batch(rng, 129, 461), xtra, "ragged odd");
    }
    SUBCASE("batch exit on row i and on row i+1 of a pair") {
        std::mt19937 rng(20261028u);
        const auto pairs = exit_batch(rng);
        const int stop_xtra = KSW_XSTOP | KSW_XSUBO | KSW_XSTART | kStopQlen;
        // Check the planted stop rows: every lane reaches its stop score, and each
        // group's last lane peaks on the row built in (280 even, 279 odd), so the
        // batch can exit on either row of a pair. This does not prove the exit
        // fired: every lane has frozen by then, so the outputs are the same either
        // way. The subcase runs both exit branches for crash, indexing and
        // ordering safety; it cannot detect an exit that never fires.
        {
            const auto res = run_u8(pairs, bwa_tests::default_scoring_matrix(), stop_xtra);
            int max_te[2] = {-1, -1};
            for (size_t i = 0; i < pairs.size(); i++) {
                CAPTURE(pairs[i].tag);
                REQUIRE(res[i].score >= kStopQlen);
                const int g = static_cast<int>(i / 64);
                if (res[i].te > max_te[g]) max_te[g] = res[i].te;
            }
            CHECK(max_te[0] == 280);
            CHECK(max_te[1] == 279);
        }
        // Independent oracle: scalar ksw_align2 under the same stop flags, on
        // its u8 path (KSW_XBYTE) like the kernel under test.
        std::vector<kswr_t> scalar;
        for (const auto &p : pairs)
            scalar.push_back(bwa_tests::run_scalar_ksw(p, bwa_tests::default_scoring_matrix(),
                                                       bwa_tests::DEFAULT_GAP_OPEN,
                                                       bwa_tests::DEFAULT_GAP_EXTEND,
                                                       stop_xtra | KSW_XBYTE));
        check_rowpair(pairs, stop_xtra, "exit rows", &scalar);
    }
}

#endif // BWA_TESTS_HAVE_KSWV
