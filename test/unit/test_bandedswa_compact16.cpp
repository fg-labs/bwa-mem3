// test/unit/test_bandedswa_compact16.cpp
//
// Same-row lane compaction in getScores16 (src/bandedSWA_compact16.inc, NEON). The driver runs
// a batch in superblocks of K 8-lane groups in block lockstep, R rows per block, and at block
// ends moves the live lanes of the emptiest group into dead lanes of the others once they fit.
// It is a pure speed change: every pair must return the same six result fields (score, tle,
// qle, gscore, gtle, max_off) as the plain wrapper, the pair scored alone and the scalar kernel.
//
// Each batch is built to exercise the driver's moving parts, and the kernel's 16-bit
// compaction counters confirm that it did (lanes moved, groups retired by evacuation):
//   - staggered finishes: matches of mixed lengths with divergent and z-dropping pairs among
//     them, at 2, 3 and 8 groups per superblock and 1, 3, 8 and 16 rows per block, at wide,
//     medium and narrow bands, with the symmetric, rank-1 (--meth) and generic matrices, and
//     with seeds and lengths past the 8-bit range;
//   - stale cells: an incoming lane lands in the slot of a high-scoring lane that finished
//     left of the incoming lane's head, so the move must clear those cells;
//   - receiving groups with a narrower column range (short queries) and on a diagonal right of
//     the incoming lanes, so the move must widen the receiving group's columns;
//   - the setting's clamps, the one-group gate and the tier default.
// Removing the receiving slot's clear, or the move's nend, ncol or minq update, makes this file
// fail. (The nbeg update is kept for the same reason as the others, but no batch here, and no
// profile batch of the kernel goldens, has needed it: a group's nbeg trails its lanes' heads.)
// The reference is the same kernel copy with compaction off (set_lane_compaction16(0, 0)),
// batch for batch, each pair scored alone and the independent scalarBandedSWA recurrence.
// Ambiguous-base and tied-maximum batches also compare all six scalar fields after moves.
// Tiers without 16-bit compaction (x86) are
// checked to ignore the setting.

#include <algorithm>
#include <cstdint>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "doctest/doctest.h"
#include "bandedSWA.h"
#include "bsw_batch.h"
#include "ext_ladder.h"

#if HAVE_BSW_VECTOR_8_16

namespace {

using bwa_tests::ExtPair;
using bwa_tests::ExtResult;
using bwa_tests::ExtScoring;

ExtScoring scoring(int a, int b, int o_del, int e_del, int o_ins, int e_ins, int zdrop, int end_bonus)
{
    ExtScoring sc;
    sc.a = a; sc.b = b; sc.o_del = o_del; sc.e_del = e_del; sc.o_ins = o_ins; sc.e_ins = e_ins;
    sc.zdrop = zdrop; sc.pen_clip = end_bonus;
    sc.fill_mat();
    return sc;
}

/// `sc` with an asymmetric --meth-style matrix: kind 1 frees reference C against read T (the
/// rank-1 path); kind 2 also G against A (the generic matrix path). Kind 0 is symmetric.
ExtScoring with_matrix(ExtScoring sc, int kind)
{
    if (kind >= 1) sc.mat[1 * 5 + 3] = (int8_t)sc.a;
    if (kind >= 2) sc.mat[2 * 5 + 0] = (int8_t)sc.a;
    return sc;
}

std::unique_ptr<IBandedPairWiseSW> make(const ExtScoring &sc, int groups, int rows)
{
    std::unique_ptr<IBandedPairWiseSW> k(new BandedPairWiseSW(sc.o_del, sc.e_del, sc.o_ins, sc.e_ins, sc.zdrop,
                                                              sc.pen_clip, sc.mat, (int8_t)sc.a, (int8_t)sc.b, 1));
    k->set_lane_compaction16(groups, rows);
    return k;
}

std::string show(const ExtResult &o)
{
    return std::to_string(o.score) + "/" + std::to_string(o.tle) + "/" + std::to_string(o.qle) + "/" +
           std::to_string(o.gscore) + "/" + std::to_string(o.gtle) + "/" + std::to_string(o.max_off);
}

std::vector<uint8_t> random_bases(std::mt19937_64 &rng, int n)
{
    std::vector<uint8_t> v(n);
    for (auto &c : v) c = (uint8_t)(rng() % 4);
    return v;
}

/// A pair that matches perfectly for `len` bases; the target runs `extra` bases further.
ExtPair perfect(std::mt19937_64 &rng, int len, int extra, int h0)
{
    ExtPair p;
    p.target = random_bases(rng, len + extra);
    p.query.assign(p.target.begin(), p.target.begin() + len);
    p.h0 = h0;
    return p;
}

/// A pair with no useful alignment.
ExtPair divergent(std::mt19937_64 &rng, int len1, int len2, int h0)
{
    ExtPair p;
    p.target = random_bases(rng, len1);
    p.query = random_bases(rng, len2);
    p.h0 = h0;
    return p;
}

/// A pair aligned `d` columns right of the main diagonal (the first `d` query bases inserted).
ExtPair shifted(std::mt19937_64 &rng, int d, int len, int h0)
{
    ExtPair p;
    p.target = random_bases(rng, len + 40);
    p.query = random_bases(rng, d);
    p.query.insert(p.query.end(), p.target.begin(), p.target.begin() + len);
    p.h0 = h0;
    return p;
}

/// Staggered finishes: `groups` * 8 pairs, mostly matches of lengths spread over [8, 400] with
/// a few mismatches and 16-bit-range seeds, a quarter divergent or z-dropping.
std::vector<ExtPair> staggered_batch(std::mt19937_64 &rng, int groups)
{
    std::vector<ExtPair> v;
    for (int q = 0; q < 8 * groups; q++) {
        const int kind = (int)(rng() % 8);
        if (kind == 0) {
            const int len2 = 10 + (int)(rng() % 120);
            v.push_back(divergent(rng, len2 + (int)(rng() % 100), len2, 10 + (int)(rng() % 60)));
        } else if (kind == 1) {   // a match that z-drops part-way
            ExtPair p = perfect(rng, 40 + (int)(rng() % 200), (int)(rng() % 40), 20 + (int)(rng() % 200));
            for (size_t c = p.query.size() / 2; c < p.query.size(); c++) p.query[c] = (uint8_t)(rng() % 4);
            v.push_back(p);
        } else {
            const int len = rng() % 3 == 0 ? 200 + (int)(rng() % 200) : 8 + (int)(rng() % 120);
            ExtPair p = perfect(rng, len, (int)(rng() % 60), 5 + (int)(rng() % 300));
            for (int m = 0; m < len / 50; m++) p.query[rng() % p.query.size()] = (uint8_t)(rng() % 4);
            v.push_back(p);
        }
    }
    return v;
}

/// Score `batch` with compaction at `groups` / `rows` and with it off (the whole batch, then
/// each pair alone) and against the scalar recurrence; every field must agree. Returns the
/// number of disagreements; adds the
/// compacted run's counters to `moves` and `retired`.
int check_batch(const ExtScoring &sc, const std::vector<ExtPair> &batch, int w, int groups, int rows,
                const std::string &what, uint64_t &moves, uint64_t &retired)
{
    std::vector<const ExtPair *> ptrs;
    for (const ExtPair &p : batch) ptrs.push_back(&p);
    const std::unique_ptr<IBandedPairWiseSW> cmp = make(sc, groups, rows);
    const std::unique_ptr<IBandedPairWiseSW> plain = make(sc, 0, 0);
    const std::vector<ExtResult> got = bwa_tests::score_bsw_batch(*cmp, ptrs, 16, w);
    const std::vector<ExtResult> ref = bwa_tests::score_bsw_batch(*plain, ptrs, 16, w);
    const BswCompactCounts c = cmp->lane_compaction16_counts();
    moves += c.moves;
    retired += c.retired;
    int bad = 0;
    for (size_t q = 0; q < batch.size(); q++) {
        const ExtResult alone = bwa_tests::score_bsw_batch(*plain, {ptrs[q]}, 16, w)[0];
        const ExtResult scalar = bwa_tests::run_mem3_scalar(sc, batch[q], w);
        if (got[q] == ref[q] && got[q] == alone && got[q] == scalar) continue;
        if (bad++ < 4)
            MESSAGE(what << " K=" << groups << " R=" << rows << " w=" << w << " pair " << q << ": compacted "
                         << show(got[q]) << ", plain " << show(ref[q]) << ", alone " << show(alone)
                         << ", scalar " << show(scalar)
                         << " (len1=" << batch[q].target.size() << " len2=" << batch[q].query.size()
                         << " h0=" << batch[q].h0 << "; score/tle/qle/gscore/gtle/max_off)");
    }
    return bad;
}

} // namespace

#if defined(__aarch64__)

static void check_staggered_finishes(int selected_rep) {
    std::mt19937_64 rng(0xC0A1E516ull);
    uint64_t moves = 0, retired = 0;
    int bad = 0;
    const ExtScoring scs[] = {scoring(1, 4, 6, 1, 6, 1, 100, 5), scoring(1, 4, 6, 1, 6, 1, 20, 5),
                              scoring(2, 5, 7, 2, 5, 3, 150, 0), scoring(1, 9, 12, 1, 3, 2, 60, 20)};
    for (int groups : {2, 3, 8}) {
        for (int rows : {1, 3, 8, 16}) {
            for (int rep = 0; rep < 3; rep++) {
                // every scoring with each matrix kind (symmetric, rank-1, generic) across the reps
                const int s = (groups + rows + rep) % 4;
                const ExtScoring sc = with_matrix(scs[s], (rows + rep) % 3);
                const std::vector<ExtPair> batch = staggered_batch(rng, groups);
                if (rep != selected_rep) continue;   // keep the original RNG stream across split cases
                for (int w : {100, 40, 12})
                    bad += check_batch(sc, batch, w, groups, rows, "staggered", moves, retired);
            }
        }
    }
    MESSAGE("16-bit: " << moves << " moves, " << retired << " groups retired by evacuation");
    CHECK(moves > 100);
    CHECK(retired > 20);
    CHECK(bad == 0);
}

TEST_CASE("lane compaction 16-bit: staggered finishes, first scoring and matrix set"
          * doctest::test_suite("unit/bandedswa")) {
    check_staggered_finishes(0);
}

TEST_CASE("lane compaction 16-bit: staggered finishes, second scoring and matrix set"
          * doctest::test_suite("unit/bandedswa")) {
    check_staggered_finishes(1);
}

TEST_CASE("lane compaction 16-bit: staggered finishes, third scoring and matrix set"
          * doctest::test_suite("unit/bandedswa")) {
    check_staggered_finishes(2);
}

TEST_CASE("lane compaction 16-bit: ambiguous bases match the scalar oracle"
          * doctest::test_suite("unit/bandedswa")) {
    const ExtScoring base = scoring(1, 4, 6, 1, 6, 1, 100, 5);
    for (int matrix : {0, 1, 2}) {
        std::mt19937_64 rng(0xA6B16C7ull);
        const ExtScoring sc = with_matrix(base, matrix);
        std::vector<ExtPair> batch = staggered_batch(rng, 3);
        for (size_t q = 0; q < batch.size(); q++) {
            for (size_t j = q % 7; j < batch[q].query.size(); j += 17) batch[q].query[j] = 4;
            for (size_t j = q % 5; j < batch[q].target.size(); j += 19) batch[q].target[j] = 4;
        }
        for (int w : {12, 40, 100}) {
            uint64_t moves = 0, retired = 0;
            CHECK(check_batch(sc, batch, w, 3, 1, "ambiguous bases", moves, retired) == 0);
            CHECK(moves > 0);
        }
    }
}

TEST_CASE("lane compaction 16-bit: tied maxima match the scalar oracle"
          * doctest::test_suite("unit/bandedswa")) {
    // A match/mismatch cycle adds 1 then subtracts 1: later rows tie the first maximum.
    const ExtScoring sc = scoring(1, 1, 6, 1, 6, 1, 100, 5);
    std::vector<ExtPair> batch;
    for (int len : {12, 24, 36, 48, 60, 72, 84, 96, 108, 120, 132, 144, 156, 168, 180, 192}) {
        ExtPair tied;
        tied.query.assign(len, 0);
        tied.target.assign(len + 20, 0);
        for (size_t j = 1; j < tied.query.size(); j += 2) tied.query[j] = 1;
        tied.h0 = 20;
        batch.push_back(tied);
    }
    std::swap(batch.front(), batch.back());   // keep a long survivor among the short lanes
    for (int w : {12, 40, 100}) {
        uint64_t moves = 0, retired = 0;
        CHECK(check_batch(sc, batch, w, 2, 1, "tied maxima", moves, retired) == 0);
        CHECK(moves > 0);
    }
}

TEST_CASE("lane compaction 16-bit: a lane moved into a just-finished lane's slot ignores its stale cells"
          * doctest::test_suite("unit/bandedswa")) {
    // Superblocks of two groups. Group 0: two incoming lanes aligned d columns right of the
    // main diagonal, and short pairs that finish at once. Group 1: dead lanes that score high
    // on the main diagonal and finish at row r0, and long matches. When the dead lanes finish,
    // the live lanes fit in one group, so group 0 is evacuated and the incoming lanes land in
    // the dead lanes' slots, whose cells left of the incoming heads are large: only the move's
    // clear keeps them out of the incoming lanes' horizontal carry.
    std::mt19937_64 rng(0x5107C16Eull);
    uint64_t moves = 0, retired = 0;
    int bad = 0, batches = 0;
    for (int e : {2, 3})
        for (int w : {50, 100})
            for (int d : {4, 8})
                for (int r0 : {30, 50})
                    for (int rows : {1, 8}) {
                        const ExtScoring sc = scoring(1, 4, 6, e, 6, e, 100, 5);
                        std::vector<ExtPair> batch;
                        for (int q = 0; q < 2; q++) batch.push_back(shifted(rng, d, 110, 6 + e * d + 15));
                        while (batch.size() < 8) batch.push_back(perfect(rng, 4, 0, 10));
                        for (int q = 0; q < 2; q++) batch.push_back(perfect(rng, r0, 0, 120));
                        while (batch.size() < 16) batch.push_back(perfect(rng, 130, 10, batches % 2 ? 30 : 5));
                        const uint64_t m0 = moves;
                        bad += check_batch(sc, batch, w, 2, rows,
                                           "stale cells e=" + std::to_string(e) + " d=" + std::to_string(d) +
                                               " r0=" + std::to_string(r0),
                                           moves, retired);
                        CHECK(moves > m0);   // the incoming lanes moved
                        batches++;
                    }
    MESSAGE("16-bit: " << batches << " stale-cell batches, " << moves << " moves");
    CHECK(bad == 0);
}

TEST_CASE("lane compaction 16-bit: a lane moved into a group with a narrower column range keeps its own"
          * doctest::test_suite("unit/bandedswa")) {
    // The receiving group's live lanes have short queries (40) on long targets, so its column
    // range (nend, ncol) is far narrower than the incoming lanes'. The incoming lanes have
    // queries of 160-240 and align d columns right of the main diagonal, so right after the
    // move their best cells lie past the receiving group's nend. The move must widen nend and
    // ncol, and lower the group's minq, or the incoming lanes lose those columns.
    std::mt19937_64 rng(0x0C011A16ull);
    uint64_t moves = 0, retired = 0;
    int bad = 0;
    for (int len : {160, 240})
        for (int w : {60, 100})
            for (int rows : {1, 8})
                for (int d : {30, 50}) {
                    const ExtScoring sc = scoring(1, 4, 6, 1, 6, 1, 100, 5);
                    std::vector<ExtPair> batch;
                    for (int q = 0; q < 2; q++) batch.push_back(shifted(rng, d, len, 6 + d + 20));
                    while (batch.size() < 8) batch.push_back(perfect(rng, 4, 0, 10));
                    for (int q = 0; q < 2; q++) batch.push_back(perfect(rng, 20, 0, 60));
                    while (batch.size() < 16) batch.push_back(perfect(rng, 40, 300, 90));
                    const uint64_t m0 = moves;
                    bad += check_batch(sc, batch, w, 2, rows, "narrow receiver len=" + std::to_string(len), moves,
                                       retired);
                    CHECK(moves > m0);
                }
    CHECK(bad == 0);
}

TEST_CASE("lane compaction 16-bit: a lane moved into a group further right keeps its left columns"
          * doctest::test_suite("unit/bandedswa")) {
    // The receiving group's live lanes align d columns right of the main diagonal, while the
    // incoming lanes stay on the main diagonal. Check that these moves preserve every result
    // field at narrow bands; this batch does not isolate the move's nbeg update.
    std::mt19937_64 rng(0x1EF7C016ull);
    uint64_t moves = 0, retired = 0;
    int bad = 0;
    for (int d : {40, 80})
        for (int w : {20, 40})
            for (int rows : {1, 8}) {
                const ExtScoring sc = scoring(1, 4, 6, 1, 6, 1, 100, 5);
                std::vector<ExtPair> batch;
                for (int q = 0; q < 2; q++) batch.push_back(perfect(rng, 180, 20, 60));
                while (batch.size() < 8) batch.push_back(perfect(rng, 4, 0, 10));
                for (int q = 0; q < 2; q++) batch.push_back(perfect(rng, 24, 0, 30));
                while (batch.size() < 16) batch.push_back(shifted(rng, d, 160, 6 + d + 20));
                const uint64_t m0 = moves;
                bad += check_batch(sc, batch, w, 2, rows, "right receiver d=" + std::to_string(d), moves, retired);
                CHECK(moves > m0);
            }
    CHECK(bad == 0);
}

TEST_CASE("lane compaction 16-bit: the default setting, the clamps and the one-group gate"
          * doctest::test_suite("unit/bandedswa")) {
    std::mt19937_64 rng(0xDEFA1716ull);
    const ExtScoring sc = scoring(1, 4, 6, 1, 6, 1, 100, 5);
    const std::vector<ExtPair> batch = staggered_batch(rng, 4);
    std::vector<const ExtPair *> ptrs;
    for (const ExtPair &p : batch) ptrs.push_back(&p);
    uint64_t moves = 0, retired = 0;
    SUBCASE("the default compacts at the tier's groups and rows") {
        const int lane_groups = BSW_COMPACT_GROUPS_MAX + 1;
        const std::vector<ExtPair> default_batch = staggered_batch(rng, lane_groups);
        std::vector<const ExtPair *> default_ptrs;
        for (const ExtPair &p : default_batch) default_ptrs.push_back(&p);
        std::unique_ptr<IBandedPairWiseSW> k(new BandedPairWiseSW(sc.o_del, sc.e_del, sc.o_ins, sc.e_ins, sc.zdrop,
                                                                  sc.pen_clip, sc.mat, (int8_t)sc.a, (int8_t)sc.b, 1));
        const std::unique_ptr<IBandedPairWiseSW> plain = make(sc, 0, 0);
        const std::vector<ExtResult> got = bwa_tests::score_bsw_batch(*k, default_ptrs, 16, 100);
        const std::vector<ExtResult> ref = bwa_tests::score_bsw_batch(*plain, default_ptrs, 16, 100);
        CHECK(got == ref);
        const int groups = std::min(BSW16_COMPACT_GROUPS_NEON, BSW_COMPACT_GROUPS_MAX);
        const uint64_t expected = groups > 0 ? (uint64_t)((lane_groups + groups - 1) / groups) : 0;
        CHECK(k->lane_compaction16_counts().superblocks == expected);
        CHECK(plain->lane_compaction16_counts().superblocks == 0);
    }
    SUBCASE("one lane group runs the plain wrapper") {
        const std::unique_ptr<IBandedPairWiseSW> k = make(sc, 8, 8);
        const std::vector<const ExtPair *> one(ptrs.begin(), ptrs.begin() + 8);
        bwa_tests::score_bsw_batch(*k, one, 16, 100);
        CHECK(k->lane_compaction16_counts().superblocks == 0);
    }
    SUBCASE("groups past the cap use the capped number of groups per superblock") {
        const std::vector<ExtPair> big = staggered_batch(rng, BSW_COMPACT_GROUPS_MAX + 1);
        std::vector<const ExtPair *> big_ptrs;
        for (const ExtPair &p : big) big_ptrs.push_back(&p);
        const std::unique_ptr<IBandedPairWiseSW> capped = make(sc, 1000, 1);
        const std::unique_ptr<IBandedPairWiseSW> plain = make(sc, 0, 0);
        CHECK(bwa_tests::score_bsw_batch(*capped, big_ptrs, 16, 100) ==
              bwa_tests::score_bsw_batch(*plain, big_ptrs, 16, 100));
        CHECK(capped->lane_compaction16_counts().superblocks == 2);
    }
    SUBCASE("one configured group per superblock runs four superblocks without moves") {
        const std::unique_ptr<IBandedPairWiseSW> one = make(sc, 1, 1);
        const std::unique_ptr<IBandedPairWiseSW> plain = make(sc, 0, 0);
        CHECK(bwa_tests::score_bsw_batch(*one, ptrs, 16, 100) ==
              bwa_tests::score_bsw_batch(*plain, ptrs, 16, 100));
        CHECK(one->lane_compaction16_counts().superblocks == 4);
        CHECK(one->lane_compaction16_counts().moves == 0);
    }
    SUBCASE("nonpositive row settings preserve result records") {
        CHECK(check_batch(sc, batch, 100, 1000, 0, "zero rows", moves, retired) == 0);
        CHECK(check_batch(sc, batch, 100, 1, -5, "negative rows", moves, retired) == 0);
    }
}

#else   // tiers without 16-bit compaction

TEST_CASE("lane compaction 16-bit: tiers without it ignore the setting"
          * doctest::test_suite("unit/bandedswa")) {
    std::mt19937_64 rng(0x0FF16C0Dull);
    const ExtScoring sc = scoring(1, 4, 6, 1, 6, 1, 100, 5);
    const std::vector<ExtPair> batch = staggered_batch(rng, 4);
    uint64_t moves = 0, retired = 0;
    CHECK(check_batch(sc, batch, 100, 8, 8, "ignored", moves, retired) == 0);
    CHECK(moves == 0);
    CHECK(retired == 0);
}

#endif
#endif   // HAVE_BSW_VECTOR_8_16
