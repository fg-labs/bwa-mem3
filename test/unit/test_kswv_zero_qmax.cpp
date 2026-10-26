// test/unit/test_kswv_zero_qmax.cpp
//
// The batched kswv kernels on a scoring whose largest substitution score
// (qmax) is 0.
//
// kswv stores the match and mismatch scores as int8_t, so a match score of 128
// or more wraps: -A 128 with the default mismatch scaling gives match -128 and
// mismatch 0, and -A 256 gives match 0. The score2 window divides each lane's
// score by qmax, so qmax == 0 raised SIGFPE on x86 (integer divide by zero)
// and was undefined behaviour on NEON. The AVX-512BW 16-bit kernel divided for
// every batch; the others only when a lane was live (8-bit: not saturated,
// 16-bit: scored above 0). The kernels now report no second-best
// (score2 = te2 = -1) when qmax <= 0, and the AVX-512BW 16-bit kernel returns
// early on a batch with no live lane, as the others do.
//
// A batch with no live lane also reports no second-best now on every kernel.
// Five of the six used to return before writing score2 and te2, leaving
// whatever the result buffer held, which the rescue reads as the hit's csub.
// And the 8-bit admission (matesw_use_u8) and the 8-bit kernels' saturation
// guard bound the score by the largest substitution score, not the match
// score, which a negative mismatch penalty exceeds.
//
// These tests call the kernels directly, so they hold whatever scoring the
// command line lets through. They run on the tier make_kswv picks (or
// BWAMEM3_FORCE_TIER selects), so the AVX-512BW kernels are only exercised on
// a host that has AVX-512BW.

#include "doctest/doctest.h"

#include "bwa.h"         // bwa_fill_scmat
#include "bwamem.h"      // mem_opt_init
#include "matesw_u8.h"   // matesw_u8_qmax(_for), matesw_u8_shift(_for), matesw_use_u8

#include <cstdlib>
#include <cstring>

#include "kswv_runner.h"  // BWA_TESTS_HAVE_KSWV
#include "seqpair_batch.h"
#include "seqpair_gen.h"

#if BWA_TESTS_HAVE_KSWV

#include "kswv.h"
#include "simd_dispatch.h"  // make_kswv, bwamem3_simd_tier

#include <random>
#include <vector>

namespace {

constexpr int kRefLen = 300;
constexpr int kQryLen = 151;

bool batched_kswv_available() {
    bwamem3_simd_init();
    const int tier = bwamem3_simd_tier();
    return tier == BWAMEM3_TIER_AVX2 || tier == BWAMEM3_TIER_AVX512BW ||
           tier == BWAMEM3_TIER_NEON;
}

// More pairs than one AVX-512BW 8-bit group (64 lanes), so every tier runs at
// least one full group and one partial group.
std::vector<bwa_tests::TestPair> batch(uint64_t seed) {
    std::mt19937 rng(static_cast<unsigned>(seed));
    std::vector<bwa_tests::TestPair> pairs;
    for (int i = 0; i < 70; i++)
        pairs.push_back(bwa_tests::gen_random_pair(rng, kQryLen, kRefLen));
    return pairs;
}

// Run phase 0 with match = mismatch = 0 (qmax 0) and the given gap penalties.
// Returns the number of pairs that scored above 0, and checks that every pair
// reports no second-best.
int run_zero_qmax(bool use16, int gap_open, int gap_extend) {
    const auto pairs = batch(use16 ? 2 : 1);
    // minsc 0: every lane is in the KSW_XSUBO mask, so a lane that is live
    // reaches the score2 window.
    bwa_tests::BatchBuffers bb(pairs, KSW_XSUBO | KSW_XSTART | (use16 ? 0 : KSW_XBYTE));
    auto pwsw = make_kswv(gap_open, gap_extend, gap_open, gap_extend,
                          /*w_match=*/0, /*w_mismatch=*/0, 1, kRefLen, kQryLen);
    REQUIRE(pwsw->needsScalar() == false);
    if (use16)
        pwsw->getScores16(bb.pairs(), bb.ref_buf(), bb.qer_buf(), bb.aln(), bb.n(), 1, 0);
    else
        pwsw->getScores8(bb.pairs(), bb.ref_buf(), bb.qer_buf(), bb.aln(), bb.n(), 1, 0);

    int scored = 0;
    for (int i = 0; i < bb.n(); i++) {
        const kswr_t &r = bb.aln()[bb.pairs()[i].regid];
        CAPTURE(use16);
        CAPTURE(i);
        CHECK(r.score2 == -1);
        CHECK(r.te2 == -1);
        if (r.score > 0) scored++;
    }
    return scored;
}

// A 16-bit batch at the default scoring in which no lane can score: every
// query is all C against an all-A reference. The result buffer starts with a
// score2 and te2 no kernel writes, so a kernel that returns without writing
// them leaves the marker behind.
void run_dead_group16() {
    std::vector<bwa_tests::TestPair> pairs;
    for (int i = 0; i < 70; i++) pairs.push_back(bwa_tests::gen_all_mismatch_pair(kQryLen));
    bwa_tests::BatchBuffers bb(pairs, KSW_XSUBO | KSW_XSTART | 30);
    for (int i = 0; i < bb.n() + bwa_tests::MAX_SIMD_WIDTH8; i++) {
        bb.aln()[i].score2 = 12345;
        bb.aln()[i].te2 = 12345;
    }
    auto pwsw = make_kswv(6, 1, 6, 1, /*w_match=*/1, /*w_mismatch=*/-4, 1,
                          2 * kQryLen, kQryLen);
    REQUIRE(pwsw->needsScalar() == false);
    pwsw->getScores16(bb.pairs(), bb.ref_buf(), bb.qer_buf(), bb.aln(), bb.n(), 1, 0);
    for (int i = 0; i < bb.n(); i++) {
        const kswr_t &r = bb.aln()[bb.pairs()[i].regid];
        CAPTURE(i);
        CHECK(r.score <= 0);
        CHECK(r.score2 == -1);
        CHECK(r.te2 == -1);
    }
}

} // namespace

TEST_CASE("kswv::getScores16 with no live lane reports no second-best"
          * doctest::test_suite("unit/kswv")) {
    if (!batched_kswv_available()) {
        MESSAGE("tier has no batched kswv kernel; skipping");
        return;
    }
    run_dead_group16();
}

TEST_CASE("kswv::getScores8 with a zero maximum substitution score does not divide by it"
          * doctest::test_suite("unit/kswv")) {
    if (!batched_kswv_available()) {
        MESSAGE("tier has no batched kswv kernel; skipping");
        return;
    }
    // Every lane scores 0, which the 8-bit kernels count as live (not
    // saturated), so each group reaches the score2 window.
    //
    // There is no 8-bit counterpart of run_dead_group16. An 8-bit lane is dead
    // only when it saturates (score 255), and kswv_u8_saturation_guard, which
    // is active in release builds, aborts on any group whose longest query
    // could reach that (len2 * qmax + shift > 254) before it is scored. With
    // qmax <= 0 and gap penalties of at least 0 (open) and 1 (extend), as the
    // parser requires, no lane scores above 0. So the 8-bit live == 0 return
    // cannot be reached without the guard firing first.
    CHECK(run_zero_qmax(/*use16=*/false, 6, 1) == 0);
}

TEST_CASE("kswv::getScores16 with a zero maximum substitution score and no live lane"
          * doctest::test_suite("unit/kswv")) {
    if (!batched_kswv_available()) {
        MESSAGE("tier has no batched kswv kernel; skipping");
        return;
    }
    // No lane scores above 0. The AVX-512BW 16-bit kernel used to divide here.
    CHECK(run_zero_qmax(/*use16=*/true, 6, 1) == 0);
}

TEST_CASE("kswv::getScores16 with a zero maximum substitution score and live lanes"
          * doctest::test_suite("unit/kswv")) {
    if (!batched_kswv_available()) {
        MESSAGE("tier has no batched kswv kernel; skipping");
        return;
    }
    // A negative gap extension rewards gaps, so lanes score above 0 with no
    // positive substitution score: the live lanes reach the score2 window on
    // every tier, and the test fails if none did.
    CHECK(run_zero_qmax(/*use16=*/true, 1, -2) > 0);
}

#endif  // BWA_TESTS_HAVE_KSWV

TEST_CASE("matesw_use_u8 bounds the score by the largest substitution score"
          * doctest::test_suite("unit/kswv")) {
    int8_t mat[25];
    // Default scoring: the match score is the largest entry, so the bound is
    // the one the match score gave.
    bwa_fill_scmat(1, 4, mat);
    CHECK(matesw_u8_qmax(mat) == 1);
    CHECK(matesw_use_u8(200, matesw_u8_qmax(mat), matesw_u8_shift(mat)));
    // A negative mismatch penalty scores a mismatch +4: a 200-base mate can
    // reach 800, which the 8-bit kernels cannot hold.
    bwa_fill_scmat(1, -4, mat);
    CHECK(matesw_u8_qmax(mat) == 4);
    CHECK_FALSE(matesw_use_u8(200, matesw_u8_qmax(mat), matesw_u8_shift(mat)));
    CHECK(matesw_use_u8(50, matesw_u8_qmax(mat), matesw_u8_shift(mat)));
}

TEST_CASE("matesw_u8_qmax_for folds the --meth matrices into the 8-bit admission"
          * doctest::test_suite("unit/kswv")) {
    mem_opt_t *o = mem_opt_init();
    REQUIRE(o != NULL);
    bwa_fill_scmat(1, 4, o->mat);
    // Synthetic per-hypothesis matrices with larger maxima than the base: OT
    // tops out at 3 and OB at 5, so each fold is visible on its own.
    memcpy(o->mat_ot, o->mat, sizeof(o->mat_ot));
    memcpy(o->mat_ob, o->mat, sizeof(o->mat_ob));
    o->mat_ot[1 * 5 + 3] = 3;
    o->mat_ob[2 * 5 + 0] = 5;
    const int shift = matesw_u8_shift_for(o, o->mat);

    // Outside --meth the per-hypothesis matrices are not used, so they do not count.
    o->meth_mode = 0;
    CHECK(matesw_u8_qmax_for(o, o->mat) == 1);
    CHECK(matesw_use_u8(200, matesw_u8_qmax_for(o, o->mat), shift));

    // Under --meth the largest entry of any of the three matrices is the bound.
    o->meth_mode = 1;
    CHECK(matesw_u8_qmax_for(o, o->mat) == 5);
    // A 200-base mate can reach 1000 under OB, which the 8-bit kernels cannot hold.
    CHECK_FALSE(matesw_use_u8(200, matesw_u8_qmax_for(o, o->mat), shift));
    CHECK(matesw_use_u8(40, matesw_u8_qmax_for(o, o->mat), shift));

    // OT alone above the base also moves the bound.
    o->mat_ob[2 * 5 + 0] = o->mat[2 * 5 + 0];
    CHECK(matesw_u8_qmax_for(o, o->mat) == 3);
    CHECK_FALSE(matesw_use_u8(200, matesw_u8_qmax_for(o, o->mat), shift));
    free(o);
}
