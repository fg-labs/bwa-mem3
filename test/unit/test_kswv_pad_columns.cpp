// test/unit/test_kswv_pad_columns.cpp
//
// Two kswv rescue-kernel cuts, each behind its own switch (NEON u8 and 16-bit,
// and the AVX-512BW u8 kernel):
//   BWA3_RESCUE_PADCOLS    the u8 FScan bodies run [jzero, jsplit), the columns
//                          past the group's longest query and short of its
//                          narrowest quantum, with no score gather (the score is
//                          0 in every lane there);
//   BWA3_RESCUE_QE_REVSCAN the lazy query-end rescan (u8 and 16-bit) scans the
//                          reached block last to first.
// Every on/off combination, in the two-row sweep and the one-row body, must give
// the same result fields as both switches off, and that reference must match
// scalar ksw_align2 (ksw_u8 for the 8-bit kernel, ksw_i16 for the 16-bit one).
// The batch is laid out as explicit lane groups of the kernel's width, 16 on
// NEON and 64 on AVX-512BW (kGroup; the runner keeps input
// order), so the zero range is non-empty (uniform lengths that are not a multiple
// of 16, and ragged lengths inside one quantum), empty (a lane whose length is a
// multiple of 16, mixed quanta), or absent (a partial final group). The rescan
// order only shows where the row that sets a lane's best score ties its max at
// several columns: a reference shorter than its query, with both homopolymer or
// periodic, puts the last row's max at every column from the first hit on,
// inside the first hit's QE_BLK block and in every block after it.

#include <random>
#include <string>
#include <vector>

#include "doctest/doctest.h"

#include "kswr_cmp.h"
#include "scoped_env.h"
#include "ksw_runner.h"
#include "kswv_runner.h"
#include "scoring.h"
#include "seqpair_gen.h"

#if BWA_TESTS_HAVE_KSWV

namespace {

using bwa_tests::ScopedEnv;   // scoped_env.h

// Lanes per u8 kernel group: a group's zero range is [max len2, min quantum), so
// the groups below are laid out at the width of the kernel that runs them.
#if defined(__AVX512BW__) && !defined(__ARM_NEON) && !defined(__aarch64__)
constexpr int kGroup = 64;
#else
constexpr int kGroup = 16;
#endif

// A query of period `period` (random bases) against a reference holding copies
// of it with a few substitutions, so a row's max recurs at columns `period` apart.
bwa_tests::TestPair periodic_pair(std::mt19937 &rng, int qlen, int rlen, int period)
{
    std::uniform_int_distribution<int> base(0, 3);
    std::vector<uint8_t> unit(period);
    for (auto &b : unit) b = static_cast<uint8_t>(base(rng));
    bwa_tests::TestPair p;
    p.qry.resize(qlen);
    for (int j = 0; j < qlen; j++) p.qry[j] = unit[j % period];
    p.ref.resize(rlen);
    std::uniform_int_distribution<int> pos(0, rlen - 1);
    for (int i = 0; i < rlen; i++) p.ref[i] = unit[i % period];
    for (int k = 0; k < rlen / 40; k++) p.ref[pos(rng)] = static_cast<uint8_t>(base(rng));
    p.tag = "periodic_" + std::to_string(qlen) + "_" + std::to_string(period);
    return p;
}

// A reference shorter than its query, both runs of one base (or of one short
// period): the last row, which sets the best score, reaches its max at the
// reference's last base and holds it at every later query column.
bwa_tests::TestPair tie_pair(std::mt19937 &rng, int qlen, int rlen, int period)
{
    std::uniform_int_distribution<int> base(0, 3);
    std::vector<uint8_t> unit(period);
    for (auto &b : unit) b = static_cast<uint8_t>(base(rng));
    bwa_tests::TestPair p;
    p.qry.resize(qlen);
    p.ref.resize(rlen);
    for (int j = 0; j < qlen; j++) p.qry[j] = unit[j % period];
    for (int i = 0; i < rlen; i++) p.ref[i] = unit[i % period];
    p.tag = "tie_" + std::to_string(qlen) + "_" + std::to_string(rlen) + "_" + std::to_string(period);
    return p;
}

// One group whose query lengths are drawn from [qlo, qhi] and whose row
// maxima tie (tie_pair), with references from a quarter to a half of the query.
void add_tie_group(std::vector<bwa_tests::TestPair> &pairs, std::mt19937 &rng, int qlo, int qhi)
{
    std::uniform_int_distribution<int> qd(qlo, qhi);
    for (int l = 0; l < kGroup; l++) {
        const int q = qd(rng);
        std::uniform_int_distribution<int> rd(q / 4, q / 2);
        pairs.push_back(tie_pair(rng, q, rd(rng), 1 + l % 3));
    }
}

// One group of 150 bp queries against a 40 bp run of one base, each query
// starting with a mismatching prefix of 7 * (lane mod 16) bases. Every lane improves on
// the same rows (the reference's), and on each of them lane l first reaches the
// row max at column 7 * l + row and ties it from there on, so the lanes' first
// hits fall in different QE_BLK blocks while the earlier lanes still tie in the
// later blocks. That is what tells a rescan that keeps a lane's first block from
// one that lets a later block's scan overwrite it.
void add_offset_tie_group(std::vector<bwa_tests::TestPair> &pairs)
{
    for (int l = 0; l < kGroup; l++) {
        bwa_tests::TestPair p;
        const int s = 7 * (l % 16);
        p.qry.assign(150, 0);
        for (int j = 0; j < s; j++) p.qry[j] = 3;   // never matches the reference's 0s
        p.ref.assign(40, 0);
        p.tag = "offset_tie_" + std::to_string(s);
        pairs.push_back(p);
    }
}

// One group whose query lengths are drawn from [qlo, qhi].
void add_group(std::vector<bwa_tests::TestPair> &pairs, std::mt19937 &rng, int qlo, int qhi,
               int nlanes = kGroup)
{
    std::uniform_int_distribution<int> qd(qlo, qhi);
    for (int l = 0; l < nlanes; l++) {
        const int q = qd(rng);
        switch (l % 4) {
        case 0: pairs.push_back(bwa_tests::gen_tandem_repeat_pair(rng, q, 2 * q)); break;
        case 1: pairs.push_back(periodic_pair(rng, q, 2 * q + 7, 3 + l % 5)); break;
        case 2: pairs.push_back(bwa_tests::gen_with_n_bases_pair(rng, q, 2 * q, 5)); break;
        default: pairs.push_back(bwa_tests::gen_random_pair(rng, q, 2 * q)); break;
        }
    }
}

std::vector<bwa_tests::TestPair> build_groups(bool short_only)
{
    std::mt19937 rng(20261023);
    std::vector<bwa_tests::TestPair> pairs;
    add_group(pairs, rng, 150, 150);   // zero range [150, 160)
    add_group(pairs, rng, 100, 100);   // [100, 112)
    add_group(pairs, rng, 145, 150);   // ragged, one quantum: [150, 160)
    add_group(pairs, rng, 23, 23);     // [23, 32)
    add_group(pairs, rng, 101, 101);   // odd range length
    add_group(pairs, rng, 140, 144);   // empty if a lane draws 144 (its quantum), else [max, 144)
    add_group(pairs, rng, 60, 128);    // mixed quanta: empty
    add_tie_group(pairs, rng, 150, 150);   // ties in [150, 160) and across QE_BLK blocks
    add_tie_group(pairs, rng, 60, 75);
    add_offset_tie_group(pairs);
    if (!short_only) add_group(pairs, rng, 200, 250);
    add_group(pairs, rng, 150, 150, 7);   // partial final group
    return pairs;
}

struct Arm { const char *padcols, *revscan, *rowpair; };

void check_arms(const std::vector<bwa_tests::TestPair> &pairs, bool use16)
{
    // The zero range needs the FScan body with USQADD on, and the two-row sweep
    // rescans lazily only with LAZYQE on: pin all three so the arms test what
    // they name whatever the environment sets.
    ScopedEnv fs("BWA3_RESCUE_FSCAN", "1"), uq("BWA3_RESCUE_USQADD", "1"),
              lq("BWA3_RESCUE_LAZYQE", "1");
    auto mat = bwa_tests::build_scoring_matrix(1, 4, 1);
    // The scalar oracle at the kernel's own width: ksw_u8 for the 8-bit kernel,
    // ksw_i16 (KSW_XBYTE cleared) for the 16-bit one, whose second-best
    // bookkeeping differs from ksw_u8's on tied rows.
    const int xtra = bwa_tests::default_xtra_flags(1) & (use16 ? ~KSW_XBYTE : ~0);
    std::vector<kswr_t> scalar;
    scalar.reserve(pairs.size());
    for (const auto &p : pairs)
        scalar.push_back(bwa_tests::run_scalar_ksw(p, mat, bwa_tests::DEFAULT_GAP_OPEN,
                                                   bwa_tests::DEFAULT_GAP_EXTEND, xtra));

    std::vector<kswr_t> ref;
    {
        ScopedEnv a("BWA3_RESCUE_PADCOLS", "0"), b("BWA3_RESCUE_QE_REVSCAN", "0"),
                  c("BWA3_RESCUE_ROWPAIR", "1");
        ref = bwa_tests::run_kswv_batch(pairs, mat, bwa_tests::DEFAULT_GAP_OPEN, bwa_tests::DEFAULT_GAP_EXTEND, 0, use16);
    }
    REQUIRE(ref.size() == pairs.size());

    // BWA3_RESCUE_PADCOLS reaches only the 8-bit kernel: the 16-bit runs vary
    // the rescan and the sweep only.
    const std::vector<Arm> arms = use16
        ? std::vector<Arm>{{"0", "1", "1"}, {"0", "0", "0"}, {"0", "1", "0"}}
        : std::vector<Arm>{{"1", "0", "1"}, {"0", "1", "1"}, {"1", "1", "1"},
                           {"0", "0", "0"}, {"1", "0", "0"}, {"0", "1", "0"}, {"1", "1", "0"}};
    for (const Arm &arm : arms) {
        ScopedEnv a("BWA3_RESCUE_PADCOLS", arm.padcols), b("BWA3_RESCUE_QE_REVSCAN", arm.revscan),
                  c("BWA3_RESCUE_ROWPAIR", arm.rowpair);
        const auto got = bwa_tests::run_kswv_batch(pairs, mat, bwa_tests::DEFAULT_GAP_OPEN, bwa_tests::DEFAULT_GAP_EXTEND, 0, use16);
        REQUIRE(got.size() == pairs.size());
        int drift = 0;
        for (size_t i = 0; i < pairs.size(); i++) {
            if (bwa_tests::kswr_all_fields_eq(ref[i], got[i])) continue;
            ++drift;
            CAPTURE(arm.padcols); CAPTURE(arm.revscan); CAPTURE(arm.rowpair);
            CAPTURE(i); CAPTURE(pairs[i].tag);
            CAPTURE(ref[i].score); CAPTURE(got[i].score);
            CAPTURE(ref[i].qe); CAPTURE(got[i].qe); CAPTURE(ref[i].te); CAPTURE(got[i].te);
            CAPTURE(ref[i].score2); CAPTURE(got[i].score2);
            CAPTURE(ref[i].qb); CAPTURE(got[i].qb); CAPTURE(ref[i].tb); CAPTURE(got[i].tb);
            CHECK(bwa_tests::kswr_all_fields_eq(ref[i], got[i]));
        }
        CHECK(drift == 0);
    }

    // Every arm equals the reference, so checking it against the scalar oracle
    // checks them all.
    int mism = 0;
    for (size_t i = 0; i < pairs.size(); i++) {
        // kswr_coords_eq compares tb/qb and only at scores of a seed length or
        // more; the end (te, qe) is what the rescan computes, so compare it at
        // every nonzero score.
        const bool end_ok = scalar[i].score == 0 ||
                            (scalar[i].te == ref[i].te && scalar[i].qe == ref[i].qe);
        const bool ok = bwa_tests::kswr_score_eq(scalar[i], ref[i]) && end_ok &&
                        bwa_tests::kswr_coords_eq(scalar[i], ref[i]) &&
                        bwa_tests::kswr_score2_eq(scalar[i], ref[i]);
        if (!ok) {
            ++mism;
            CAPTURE(i); CAPTURE(pairs[i].tag);
            CAPTURE(scalar[i].score); CAPTURE(ref[i].score);
            CAPTURE(scalar[i].qe); CAPTURE(ref[i].qe); CAPTURE(scalar[i].te); CAPTURE(ref[i].te);
            CAPTURE(scalar[i].score2); CAPTURE(ref[i].score2);
            CHECK(ok);
        }
    }
    CHECK(mism == 0);
}

} // namespace

TEST_CASE("kswv u8 rescue: BWA3_RESCUE_PADCOLS and BWA3_RESCUE_QE_REVSCAN on == off, in both sweeps"
          * doctest::test_suite("unit/kswv")) {
#if !defined(__ARM_NEON) && !defined(__aarch64__) && !defined(__AVX512BW__)
    MESSAGE("skipped: both switches affect the NEON and AVX-512BW u8 rescue kernels only");
#else
    check_arms(build_groups(true), false);
#endif
}

TEST_CASE("kswv u16 rescue: BWA3_RESCUE_QE_REVSCAN on == off, in both sweeps"
          * doctest::test_suite("unit/kswv")) {
#if !defined(__ARM_NEON) && !defined(__aarch64__)
    MESSAGE("skipped: the switch affects the NEON 16-bit rescue kernel only");
#else
    check_arms(build_groups(false), true);
#endif
}

#endif  // BWA_TESTS_HAVE_KSWV
