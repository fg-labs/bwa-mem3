// test/unit/test_kswv_dash_query.cpp
//
// A '-' query base must score the same in every batch composition and at every
// tier.
//
// nst_nt4_table maps FASTQ '-' to 5, and matesw_orient copies a forward-oriented
// mate's codes through unchanged, so a 5 can reach kswv as a query byte. The
// x86 8-bit kernels used to pad the query tail with DUMMY5 = 5 as well, and
// they elide the "is this the pad?" compare and blend below the group's jdummy
// (compute_jdummy: the first column any lane pads at, which depends on the
// OTHER lanes' len2). A '-' below jdummy was table-scored (w_ambig) and one at
// or after it was scored as pad (a zero delta), so the same job scored
// differently depending on which jobs shared its lane group. The NEON kernel
// pads with NEON_QPAD8 = 16, which no query byte can take, and always
// table-scored '-'.
//
// Separately, the AVX-512BW 16-bit kernel's score table held 0 at the xor
// indices only a '-' query base reaches (5 ^ {0..3} = 4..7, 5 ^ AMBR16 = 10),
// where the AVX2 and NEON 16-bit tables hold w_ambig.
//
// The oracle is scalar ksw_align2 over a 6-letter matrix whose '-' column holds
// what every tier's xor table gives a '-' against a base: w_ambig (xor 4..7).
// The test references contain no N, so the '-' x N entry (which the 8-bit and
// 16-bit tables score differently) is never read.
//
// The x86 groupings only matter on an AVX2 or AVX-512BW host (or under
// BWAMEM3_FORCE_TIER there); on NEON the test checks that kernel against the
// same oracle, which it already matched.

#include <random>
#include <vector>

#include "doctest/doctest.h"

#include "ksw_runner.h"
#include "kswr_cmp.h"
#include "kswv_runner.h"  // BWA_TESTS_HAVE_KSWV
#include "scoring.h"
#include "seqpair_batch.h"

#if BWA_TESTS_HAVE_KSWV

#include "ksw.h"
#include "kswv.h"
#include "simd_dispatch.h"  // make_kswv

namespace {

constexpr int kRefLen  = 460;
constexpr uint8_t kDash = 5;   // nst_nt4_table['-']

// Same guard as test_kswv_mixed_qlen.cpp: only AVX2, AVX-512BW and NEON have a
// batched kswv kernel; the SSE classes exit().
bool batched_kswv_available() {
    bwamem3_simd_init();
    const int tier = bwamem3_simd_tier();
    return tier == BWAMEM3_TIER_AVX2 || tier == BWAMEM3_TIER_AVX512BW ||
           tier == BWAMEM3_TIER_NEON;
}

// mat6[target * 6 + query]: the 5x5 matrix extended with a '-' query column.
// `dash_vs_base` is what a '-' scores against reference A/C/G/T; every other
// N or '-' entry is w_ambig (-1), as in build_scoring_matrix.
std::vector<int8_t> dash_matrix(int match, int mismatch, int dash_vs_base) {
    std::vector<int8_t> m(36, -1);
    for (int t = 0; t < 4; t++) {
        for (int q = 0; q < 4; q++) m[t * 6 + q] = (int8_t)(t == q ? match : -mismatch);
        m[t * 6 + 5] = (int8_t)dash_vs_base;
    }
    return m;
}

kswr_t scalar(const bwa_tests::TestPair &p, std::vector<int8_t> &mat6, int xtra) {
    return ksw_align2((int)p.qry.size(), const_cast<uint8_t *>(p.qry.data()),
                      (int)p.ref.size(), const_cast<uint8_t *>(p.ref.data()),
                      6, mat6.data(),
                      bwa_tests::DEFAULT_GAP_OPEN, bwa_tests::DEFAULT_GAP_EXTEND,
                      bwa_tests::DEFAULT_GAP_OPEN, bwa_tests::DEFAULT_GAP_EXTEND,
                      xtra, nullptr);
}

std::vector<uint8_t> random_bases(std::mt19937 &rng, int n) {
    std::uniform_int_distribution<int> base(0, 3);
    std::vector<uint8_t> s((size_t)n);
    for (auto &b : s) b = (uint8_t)base(rng);
    return s;
}

// A pair whose query is ref[off, off + len), so it aligns end to end.
bwa_tests::TestPair exact_pair(std::mt19937 &rng, int off, int len) {
    bwa_tests::TestPair p;
    p.ref = random_bases(rng, kRefLen);
    p.qry.assign(p.ref.begin() + off, p.ref.begin() + off + len);
    return p;
}

// The '-' job: len2 = 100 (quantum 112, so its own pad starts at 100) with a
// '-' at column 50.
constexpr int kDashLen = 100;
constexpr int kDashCol = 50;

bwa_tests::TestPair dash_pair(uint64_t seed) {
    std::mt19937 rng((unsigned)seed);
    bwa_tests::TestPair p = exact_pair(rng, 120, kDashLen);
    p.qry[kDashCol] = kDash;
    p.tag = "dash";
    return p;
}

// 64 pairs, a whole number of groups at every 8-bit width (16, 32, 64), so no
// group is partial: a partial group's dummy lanes have len2 = 0, which pins
// jsplit (and with it the FScan bodies' jdummy) to column 0, and would put the
// '-' in the pad-tested range in both compositions. Lane 0 is the '-' job and
// every other lane's len2 is a multiple of 16, so it pads at no column and the
// group's jdummy is the '-' job's own len2, above the '-' column.
// `with_short_lane` makes lane 1 a len2 = 20 job, which pads from column 20 and
// so pulls the jdummy of the group holding lanes 0 and 1 below the '-' column.
std::vector<bwa_tests::TestPair> batch(uint64_t seed, bool with_short_lane) {
    constexpr int kPairs = 64;
    const int lens[] = {96, 112, 128};
    std::mt19937 rng((unsigned)(seed * 7919 + 1));
    std::vector<bwa_tests::TestPair> pairs;
    pairs.push_back(dash_pair(seed));
    if (with_short_lane) pairs.push_back(exact_pair(rng, 200, 20));
    while ((int)pairs.size() < kPairs) {
        pairs.push_back(exact_pair(rng, 60, lens[pairs.size() % 3]));
    }
    return pairs;
}

// Run the batch through the runtime tier's kernel (phase 0) and return the
// kswr_t of every pair, in input order.
std::vector<kswr_t> run_batch(const std::vector<bwa_tests::TestPair> &pairs,
                              int match, int mismatch, bool use16, int xtra) {
    bwa_tests::BatchBuffers bb(pairs, xtra);
    auto pwsw = make_kswv(bwa_tests::DEFAULT_GAP_OPEN, bwa_tests::DEFAULT_GAP_EXTEND,
                          bwa_tests::DEFAULT_GAP_OPEN, bwa_tests::DEFAULT_GAP_EXTEND,
                          match, -mismatch, 1, kRefLen, 256);
    REQUIRE(pwsw->needsScalar() == false);
    if (use16) pwsw->getScores16(bb.pairs(), bb.ref_buf(), bb.qer_buf(), bb.aln(), bb.n(), 1, 0);
    else       pwsw->getScores8 (bb.pairs(), bb.ref_buf(), bb.qer_buf(), bb.aln(), bb.n(), 1, 0);
    std::vector<kswr_t> out;
    for (int i = 0; i < bb.n(); i++) out.push_back(bb.aln()[i]);
    return out;
}

void check_dash(uint64_t seed, bool use16) {
    if (!batched_kswv_available()) {
        MESSAGE("tier has no batched kswv kernel; skipping");
        return;
    }
    // match=1 keeps KSW_XBYTE (byte mode, both sides); match=14 drops it (word
    // mode, both sides), as in test_kswv_mixed_qlen.cpp.
    const int match = use16 ? 14 : 1, mismatch = use16 ? 8 : 4;
    const int xtra = bwa_tests::default_xtra_flags(match);
    std::vector<int8_t> mat6 = dash_matrix(match, mismatch, /*dash_vs_base=*/-1);

    // Not vacuous: scoring the '-' as pad (0) instead of w_ambig must change
    // the oracle's answer for the '-' job, or a pad-scored '-' would pass.
    const bwa_tests::TestPair dp = dash_pair(seed);
    std::vector<int8_t> mat6_pad = dash_matrix(match, mismatch, /*dash_vs_base=*/0);
    const kswr_t want_dash = scalar(dp, mat6, xtra);
    REQUIRE(want_dash.score != scalar(dp, mat6_pad, xtra).score);

    for (bool with_short_lane : {false, true}) {
        const auto pairs = batch(seed, with_short_lane);
        const auto got = run_batch(pairs, match, mismatch, use16, xtra);
        for (size_t i = 0; i < pairs.size(); i++) {
            const kswr_t want = scalar(pairs[i], mat6, xtra);
            CAPTURE(seed);
            CAPTURE(use16);
            CAPTURE(with_short_lane);
            CAPTURE(i);
            CHECK(got[i].score == want.score);
            CHECK(got[i].te == want.te);
            CHECK(got[i].qe == want.qe);
            CHECK(bwa_tests::kswr_score2_eq(want, got[i]));
            CHECK(got[i].te2 == want.te2);
        }
        CHECK(got[0].score == want_dash.score);
    }
}

} // namespace

TEST_CASE("kswv::getScores8 scores a '-' query base the same in every lane group"
          * doctest::test_suite("unit/kswv")) {
    for (uint64_t seed = 1; seed <= 8; ++seed) check_dash(seed, /*use16=*/false);
}

TEST_CASE("kswv::getScores16 scores a '-' query base as w_ambig at every tier"
          * doctest::test_suite("unit/kswv")) {
    for (uint64_t seed = 1; seed <= 8; ++seed) check_dash(seed, /*use16=*/true);
}

#endif  // BWA_TESTS_HAVE_KSWV
