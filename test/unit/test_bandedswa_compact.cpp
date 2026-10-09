// test/unit/test_bandedswa_compact.cpp
//
// Same-row lane compaction in getScores8 (src/bandedSWA_compact.inc on NEON,
// src/bandedSWA_compact256.inc on AVX2, src/bandedSWA_compact512.inc on
// AVX-512BW). The driver runs a batch in superblocks of K lane groups in row
// lockstep and, when the live lanes fit in one group fewer, moves the live lanes
// of the emptiest group into dead lanes of the others and retires it. It is a
// pure speed change: every pair must return
// the same six result fields (score, tle, qle, gscore, gtle, max_off) as the
// plain wrapper, and as the pair scored alone.
//
// Each batch here is built to force the driver's moving parts, and the kernel's
// compaction counters confirm that it did (lanes moved, groups retired by
// evacuation):
//   - staggered finishes: perfect-match pairs of mixed lengths with divergent
//     and z-dropping pairs among them, so groups empty at different rows and are
//     evacuated in cascades, at 2, 3, 6 and 8 groups per superblock, at the shipped
//     w >= 31 gate and at every band, with the symmetric, rank-1 (--meth) and
//     generic scoring matrices;
//   - stale cells: an incoming lane lands in the slot of a high-scoring lane
//     that just finished left of the incoming lane's head. The kernel needs
//     zero there (those cells feed the horizontal carry into the lane's first
//     in-band column), so the move must clear them;
//   - early moves, while the incoming lanes still read their own h0 seed;
//   - receiving groups on another diagonal (right or left of the incoming
//     lane's, narrow bands, longer queries), so the moved lane's band, head,
//     max_off and query end differ from everything the group had.
// Also pinned: each tier's default setting, the band and one-group gates, the
// group-count clamp, and the plan (bsw_compact_plan) on hand-made live masks.
// Removing the receiving slot's clear, dropping any per-lane field from the
// move, or skipping any of its four column-range updates (nbeg, nend, ncol,
// minq) makes this file fail.
//
// The reference is the same kernel copy with compaction off (set_lane_compaction
// (0, 0)), batch for batch, and each pair scored alone. Tiers without compaction
// (the x86 128-bit tiers, SSE4.1 to AVX) are checked to ignore the setting.

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "doctest/doctest.h"
#include "bandedSWA.h"
#include "bsw_batch.h"
#include "bsw_compact.h"
#include "ext_ladder.h"

#if HAVE_BSW_VECTOR_8_16

namespace {

using bwa_tests::ExtPair;
using bwa_tests::ExtResult;
using bwa_tests::ExtScoring;

using Factory = IBandedPairWiseSW *(*)(int, int, int, int, int, int, const int8_t *, int8_t, int8_t, int);

/// A tier that compacts, with its getScores8 lane count.
struct Tier {
    std::string name;
    Factory make;
    int lanes;
};

#if defined(__aarch64__)
/* Defined by libbwa.a's bandedSWA.o: the BSW8_ROW_LEAN setting it was built with. */
extern "C" int bsw8_row_lean_enabled(void);

IBandedPairWiseSW *make_neon(int o_del, int e_del, int o_ins, int e_ins, int zdrop, int end_bonus,
                             const int8_t *mat, int8_t a, int8_t b, int nthreads)
{
    return new BandedPairWiseSW(o_del, e_del, o_ins, e_ins, zdrop, end_bonus, mat, a, b, nthreads);
}
#endif

std::vector<Tier> compacting_tiers()
{
    std::vector<Tier> ts;
#if defined(__aarch64__)
    ts.push_back({"neon", make_neon, 16});
#elif defined(__x86_64__) || defined(__i386__)
    __builtin_cpu_init();
    if (__builtin_cpu_supports("avx2")) ts.push_back({"avx2", make_bsw_kernel_avx2, 32});
    if (__builtin_cpu_supports("avx512bw")) ts.push_back({"avx512bw", make_bsw_kernel_avx512bw, 64});
#endif
    // A tier this host cannot run is skipped, and an empty list would make every
    // equivalence test loop over zero tiers and pass vacuously, so the skip must be
    // visible in the test outcome rather than a silent MESSAGE: an empty list records
    // a WARN (counted by doctest). BWA_REQUIRE_COMPACTION_TIER turns a missing tier
    // into a hard failure, so a driver's field-by-field equivalence cannot go
    // unchecked where it must run: set to 1 it requires some compacting tier; set
    // to a tier name (the AVX-512 CI job sets avx512bw) it requires that tier, so a
    // misspelled name fails rather than passing. Unset, empty or 0: no requirement.
    const char *req = std::getenv("BWA_REQUIRE_COMPACTION_TIER");
    if (req != nullptr && req[0] != '\0' && std::string(req) != "0") {
        const std::string want(req);
        const bool named = want != "1";
        const bool found = named ? std::any_of(ts.begin(), ts.end(), [&](const Tier &t) { return t.name == want; })
                                 : !ts.empty();
        if (!found)
            FAIL("lane compaction: BWA_REQUIRE_COMPACTION_TIER=" << want << " but this host has "
                 << (named ? "no such compacting SIMD tier" : "no compacting SIMD tier")
                 << "; the equivalence tests cannot run");
    } else if (ts.empty()) {
        WARN_MESSAGE(false, "lane compaction: this host has no compacting SIMD tier (needs avx2 or "
                            "avx512bw on x86); the compacting-tier equivalence tests are skipped here");
    }
    return ts;
}

std::unique_ptr<IBandedPairWiseSW> make(Factory f, const ExtScoring &sc, int groups, int min_w)
{
    std::unique_ptr<IBandedPairWiseSW> k(
        f(sc.o_del, sc.e_del, sc.o_ins, sc.e_ins, sc.zdrop, sc.pen_clip, sc.mat, (int8_t)sc.a, (int8_t)sc.b, 1));
    k->set_lane_compaction(groups, min_w);
    return k;
}

ExtScoring scoring(int a, int b, int o_del, int e_del, int o_ins, int e_ins, int zdrop, int end_bonus)
{
    ExtScoring sc;
    sc.a = a; sc.b = b; sc.o_del = o_del; sc.e_del = e_del; sc.o_ins = o_ins; sc.e_ins = e_ins;
    sc.zdrop = zdrop; sc.pen_clip = end_bonus;
    sc.fill_mat();
    return sc;
}

/// `sc` with an asymmetric --meth-style matrix: kind 1 scores reference C against read T
/// as a match (one freed cell, the rank-1 path); kind 2 also G against A (the generic
/// matrix path). Kind 0 leaves the symmetric matrix.
ExtScoring with_matrix(ExtScoring sc, int kind)
{
    if (kind >= 1) sc.mat[1 * 5 + 3] = (int8_t)sc.a;
    if (kind >= 2) sc.mat[2 * 5 + 0] = (int8_t)sc.a;
    return sc;
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

/// A pair that matches perfectly for `len` bases (query == target prefix); the
/// target runs `extra` bases further.
ExtPair perfect(std::mt19937_64 &rng, int len, int extra, int h0)
{
    ExtPair p;
    p.target = random_bases(rng, len + extra);
    p.query.assign(p.target.begin(), p.target.begin() + len);
    p.h0 = h0;
    return p;
}

/// A pair with no useful alignment: its rows score zero almost at once.
ExtPair divergent(std::mt19937_64 &rng, int len1, int len2, int h0)
{
    ExtPair p;
    p.target = random_bases(rng, len1);
    p.query = random_bases(rng, len2);
    p.h0 = h0;
    return p;
}

/// A pair whose alignment runs on the diagonal `d` columns right of the main
/// one: the first `d` query bases are inserted, and the h0 seed pays for them.
ExtPair shifted(std::mt19937_64 &rng, int d, int len, int h0)
{
    ExtPair p;
    p.target = random_bases(rng, len + 40);
    p.query = random_bases(rng, d);
    p.query.insert(p.query.end(), p.target.begin(), p.target.begin() + len);
    p.h0 = h0;
    return p;
}

/// Staggered finishes: `groups` * `lanes` pairs, mostly perfect matches of
/// lengths spread over [8, 150], with a quarter divergent or mismatch-laden.
std::vector<ExtPair> staggered_batch(std::mt19937_64 &rng, int lanes, int groups)
{
    std::vector<ExtPair> v;
    for (int q = 0; q < lanes * groups; q++) {
        const int kind = (int)(rng() % 8);
        if (kind == 0) {
            const int len2 = 10 + (int)(rng() % 60);
            v.push_back(divergent(rng, len2 + (int)(rng() % 80), len2, 10 + (int)(rng() % 30)));
        } else if (kind == 1) {   // a match that z-drops part-way
            ExtPair p = perfect(rng, 30 + (int)(rng() % 90), (int)(rng() % 40), 20 + (int)(rng() % 40));
            for (size_t c = p.query.size() / 2; c < p.query.size(); c++) p.query[c] = (uint8_t)(rng() % 4);
            v.push_back(p);
        } else {
            const int len = rng() % 3 == 0 ? 100 + (int)(rng() % 50) : 8 + (int)(rng() % 60);
            v.push_back(perfect(rng, len, (int)(rng() % 30), 5 + (int)(rng() % 60)));
        }
    }
    return v;
}

/// Slot reuse, in superblocks of two groups. Group 0: the incoming lanes and
/// short pairs that finish by row 4. Group 1: the `dead` lanes and long perfect
/// matches seeded with `h0_long`. Once the dead lanes finish, the live lanes fit
/// in one group, so group 0 is evacuated and each incoming lane lands in a
/// just-finished lane's slot.
std::vector<ExtPair> reuse_batch(std::mt19937_64 &rng, int lanes, const std::vector<ExtPair> &incoming,
                                 const std::vector<ExtPair> &dead, int h0_long)
{
    std::vector<ExtPair> g0 = incoming, g1 = dead;
    while ((int)g0.size() < lanes) g0.push_back(perfect(rng, 4, 0, 10));
    while ((int)g1.size() < lanes) g1.push_back(perfect(rng, 130, 10, h0_long));
    g0.insert(g0.end(), g1.begin(), g1.end());
    return g0;
}

/// Score `batch` on `tier` with compaction at `groups` / `min_w` and with it off
/// (the whole batch, then each pair alone); every field must agree. Returns the
/// number of disagreements; adds the compacted run's counters to `moves` and
/// `retired`.
int check_batch(const Tier &tier, const ExtScoring &sc, const std::vector<ExtPair> &batch, int w, int groups,
                int min_w, const std::string &what, uint64_t &moves, uint64_t &retired)
{
    std::vector<const ExtPair *> ptrs;
    for (const ExtPair &p : batch) ptrs.push_back(&p);
    const std::unique_ptr<IBandedPairWiseSW> cmp = make(tier.make, sc, groups, min_w);
    const std::unique_ptr<IBandedPairWiseSW> plain = make(tier.make, sc, 0, 0);
    const std::vector<ExtResult> got = bwa_tests::score_bsw_batch(*cmp, ptrs, 8, w);
    const std::vector<ExtResult> ref = bwa_tests::score_bsw_batch(*plain, ptrs, 8, w);
    const BswCompactCounts c = cmp->lane_compaction_counts();
    moves += c.moves;
    retired += c.retired;
    int bad = 0;
    for (size_t q = 0; q < batch.size(); q++) {
        const ExtResult alone = bwa_tests::score_bsw_batch(*plain, {ptrs[q]}, 8, w)[0];
        if (got[q] == ref[q] && got[q] == alone) continue;
        if (bad++ < 4)
            MESSAGE(tier.name << " " << what << " K=" << groups << " w=" << w << " pair " << q << ": compacted "
                              << show(got[q]) << ", plain " << show(ref[q]) << ", alone " << show(alone)
                              << " (len1=" << batch[q].target.size() << " len2=" << batch[q].query.size()
                              << " h0=" << batch[q].h0 << "; score/tle/qle/gscore/gtle/max_off)");
    }
    return bad;
}

} // namespace

TEST_CASE("lane compaction: staggered finishes, cascading retirements, every field equals plain and alone"
          * doctest::test_suite("unit/bandedswa")) {
    for (const Tier &tier : compacting_tiers()) {
        std::mt19937_64 rng(0xC0A1E5CEull);
        uint64_t moves = 0, retired = 0;
        int bad = 0;
        const ExtScoring scs[] = {scoring(1, 4, 6, 1, 6, 1, 100, 5), scoring(1, 4, 6, 1, 6, 1, 20, 5),
                                  scoring(2, 5, 7, 2, 5, 3, 150, 0), scoring(1, 9, 12, 1, 3, 2, 60, 20)};
        for (int groups : {2, 3, 6, 8}) {   // 6: the AVX2 default
            for (int rep = 0; rep < 6; rep++) {
                // every scoring with each matrix kind (symmetric, rank-1, generic) across the reps
                const ExtScoring sc = with_matrix(scs[rep % 4], rep % 3);
                const std::vector<ExtPair> batch = staggered_batch(rng, tier.lanes, groups);
                for (int w : {100, 40})
                    bad += check_batch(tier, sc, batch, w, groups, 31, "staggered", moves, retired);
                bad += check_batch(tier, sc, batch, 12, groups, 0, "staggered, every band", moves, retired);
            }
        }
        MESSAGE(tier.name << ": " << moves << " moves, " << retired << " groups retired by evacuation");
        CHECK(moves > 100);
        CHECK(retired > 20);
        CHECK(bad == 0);
    }
}

TEST_CASE("lane compaction: a lane moved into a just-finished lane's slot ignores its stale cells"
          * doctest::test_suite("unit/bandedswa")) {
    // The dead lanes score high on the main diagonal; the incoming lanes align d
    // columns right of it, so the dead lane's cells left of the incoming lane's
    // head are large. Gap extends of 2 to 4 and modest scores keep every lane of
    // the receiving group off its band clamps on the row after the move, so the
    // band trim (which zeroes cells left of a lane's head when some lane's band
    // edge was clamped) does not run there: only the move's own clear stands
    // between the stale cells and the incoming lane's horizontal carry.
    for (const Tier &tier : compacting_tiers()) {
        std::mt19937_64 rng(0x5107C1EAull);
        const int n_in = std::max(1, tier.lanes / 4);
        uint64_t moves = 0, retired = 0;
        int bad = 0, batches = 0;
        for (int e : {2, 3, 4})
            for (int w : {50, 100})
                for (int d : {4, 8})
                    for (int r0 : {30, 50}) {
                        const ExtScoring sc = scoring(1, 4, 6, e, 6, e, 100, 5);
                        std::vector<ExtPair> in, dead;
                        for (int q = 0; q < n_in; q++) in.push_back(shifted(rng, d, 110, 6 + e * d + 15));
                        for (int q = 0; q < n_in; q++) dead.push_back(perfect(rng, r0, 0, 120));
                        const std::vector<ExtPair> batch =
                            reuse_batch(rng, tier.lanes, in, dead, batches % 2 ? 30 : 5);
                        const uint64_t m0 = moves;
                        bad += check_batch(tier, sc, batch, w, 2, 31,
                                           "stale cells e=" + std::to_string(e) + " d=" + std::to_string(d) +
                                               " r0=" + std::to_string(r0),
                                           moves, retired);
                        CHECK(moves - m0 == (uint64_t)n_in);   // every incoming lane moved, once
                        batches++;
                    }
        MESSAGE(tier.name << ": " << batches << " stale-cell batches, " << moves << " moves");
        CHECK(bad == 0);
    }
}

TEST_CASE("lane compaction: a lane moved early keeps its own seed, head and band edge"
          * doctest::test_suite("unit/bandedswa")) {
    // Moves at rows 5 to 9, while the incoming lanes still read their h0-prefix
    // seed (column -1 and the seed row), into the slots of low-seed lanes whose
    // trackers differ from theirs.
    for (const Tier &tier : compacting_tiers()) {
        std::mt19937_64 rng(0xEA21C0DEull);
        const int n_in = std::max(1, tier.lanes / 4);
        uint64_t moves = 0, retired = 0;
        int bad = 0;
        const ExtScoring scs[] = {scoring(1, 4, 6, 1, 6, 1, 100, 5), scoring(1, 4, 6, 2, 4, 1, 100, 5),
                                  scoring(1, 4, 3, 1, 6, 2, 100, 0)};
        for (const ExtScoring &sc : scs)
            for (int r0 : {5, 7, 9})
                for (int h0_in : {60, 110, 150}) {
                    std::vector<ExtPair> in, dead;
                    for (int q = 0; q < n_in; q++)
                        in.push_back(q % 2 ? shifted(rng, 1 + (int)(rng() % 6), 70, h0_in)
                                           : perfect(rng, 60 + (int)(rng() % 30), (int)(rng() % 20), h0_in));
                    for (int q = 0; q < n_in; q++) dead.push_back(perfect(rng, r0, (int)(rng() % 3), 3 + (int)(rng() % 10)));
                    const std::vector<ExtPair> batch = reuse_batch(rng, tier.lanes, in, dead, 10);
                    for (int w : {100, 40})
                        bad += check_batch(tier, sc, batch, w, 2, 31, "early move r0=" + std::to_string(r0), moves,
                                           retired);
                }
        MESSAGE(tier.name << ": " << moves << " early moves");
        CHECK(moves >= (uint64_t)(2 * 27 * n_in));
        CHECK(bad == 0);
    }
}

TEST_CASE("lane compaction: a lane moved into a group on another diagonal keeps its band and trackers"
          * doctest::test_suite("unit/bandedswa")) {
    // The receiving group's lanes align on a diagonal s columns right (s > 0) or
    // left (s < 0) of the main one, with steep gap extends so their bands stay
    // narrow; its dead lanes have long queries on short targets. The incoming
    // lanes align on the main diagonal with short queries, or start with a
    // deletion of k target bases (read through the column -1 seed for k rows,
    // past the move). So the moved lane's band is not inside the receiving
    // group's column range, its head, max_off and seed differ from the dead
    // lane's, and its query end comes before any of the group's.
    for (const Tier &tier : compacting_tiers()) {
        std::mt19937_64 rng(0x6E0D1A60ull);
        const int n_in = std::max(1, tier.lanes / 4);
        uint64_t moves = 0, retired = 0;
        int bad = 0;
        for (int e : {3, 4})
            for (int s : {-8, 10})
                for (int r0 : {12, 24})
                    for (int w : {50, 100}) {
                        const ExtScoring sc = scoring(1, 4, 6, e, 6, e, 100, 5);
                        // a pair on diagonal s: s > 0 inserts s query bases first, s < 0 deletes -s target bases
                        auto on_diag = [&](int len1, int len2, int h0) {
                            ExtPair p;
                            p.target = random_bases(rng, len1);
                            if (s > 0) {
                                p.query = random_bases(rng, s);
                                p.query.insert(p.query.end(), p.target.begin(), p.target.begin() + std::min(len1, len2 - s));
                            } else {
                                p.query.assign(p.target.begin() - s, p.target.begin() + std::min(len1, len2 - s));
                            }
                            while ((int)p.query.size() < len2) p.query.push_back((uint8_t)(rng() % 4));
                            p.query.resize(len2);
                            p.h0 = h0;
                            return p;
                        };
                        const int hs = 6 + e * std::abs(s) + 20;
                        std::vector<ExtPair> g0, g1;
                        for (int q = 0; q < n_in; q++) {
                            if (q % 2 == 0) {
                                g0.push_back(perfect(rng, 50 + (int)(rng() % 20), 30, 30));
                            } else {   // deletion of k target bases first: column -1 seed for k rows
                                const int k = r0 + 4;
                                ExtPair p;
                                p.query = random_bases(rng, 50);
                                p.target = random_bases(rng, k);
                                p.target.insert(p.target.end(), p.query.begin(), p.query.end());
                                p.target.resize(p.target.size() + 20, 2);
                                p.h0 = 6 + e * k + 30;
                                g0.push_back(p);
                            }
                        }
                        while ((int)g0.size() < tier.lanes) g0.push_back(perfect(rng, 4, 0, 10));
                        for (int q = 0; q < n_in; q++) g1.push_back(on_diag(r0, 120, hs + 20));   // finish at row r0
                        while ((int)g1.size() < tier.lanes) g1.push_back(on_diag(130, 120, hs));
                        g0.insert(g0.end(), g1.begin(), g1.end());
                        const uint64_t m0 = moves;
                        bad += check_batch(tier, sc, g0, w, 2, 31,
                                           "diagonal s=" + std::to_string(s) + " e=" + std::to_string(e) +
                                               " r0=" + std::to_string(r0),
                                           moves, retired);
                        CHECK(moves - m0 == (uint64_t)n_in);
                    }
        MESSAGE(tier.name << ": " << moves << " moves across diagonals");
        CHECK(bad == 0);
    }
}

TEST_CASE("lane compaction: a lane moved into the slot of a lane on a right-shifted diagonal keeps its head and max_off"
          * doctest::test_suite("unit/bandedswa")) {
    // Every lane of the receiving group matches for 8 bases, then takes a k-base
    // insertion and keeps scoring k columns right of the main diagonal, past its
    // earlier best (so its max_off is about k); steep deletions keep the left
    // edge of its band, and of the group's column range, right of the main
    // diagonal. Its dead lanes finish at row r0. The incoming lanes align on the
    // main diagonal, so a stale head, or a group range not widened to the moved
    // lane's, would cut off their own diagonal, and a stale max_off would
    // outlast theirs.
    for (const Tier &tier : compacting_tiers()) {
        std::mt19937_64 rng(0x7EAD0FF5ull);
        const int n_in = std::max(1, tier.lanes / 4);
        uint64_t moves = 0, retired = 0;
        int bad = 0;
        for (int k : {8, 12})
            for (int r0 : {32, 40})
                for (int w : {50, 100}) {
                    const ExtScoring sc = scoring(1, 4, 6, 4, 6, 1, 100, 5);
                    std::vector<ExtPair> g0, g1;
                    for (int q = 0; q < n_in; q++) g0.push_back(perfect(rng, 80 + (int)(rng() % 20), 20, 30));
                    while ((int)g0.size() < tier.lanes) g0.push_back(perfect(rng, 4, 0, 10));
                    // 8 matches, a k-base insertion, then matches k columns right of the main diagonal
                    auto right_shifted = [&](int len1) {
                        ExtPair p;
                        p.target = random_bases(rng, len1);
                        p.query.assign(p.target.begin(), p.target.begin() + 8);
                        const std::vector<uint8_t> ins = random_bases(rng, k);
                        p.query.insert(p.query.end(), ins.begin(), ins.end());
                        p.query.insert(p.query.end(), p.target.begin() + 8, p.target.begin() + std::min(len1, 110));
                        p.h0 = 10;
                        return p;
                    };
                    for (int q = 0; q < n_in; q++) g1.push_back(right_shifted(r0));
                    while ((int)g1.size() < tier.lanes) g1.push_back(right_shifted(140));
                    g0.insert(g0.end(), g1.begin(), g1.end());
                    const uint64_t m0 = moves;
                    bad += check_batch(tier, sc, g0, w, 2, 31, "right-shifted dead k=" + std::to_string(k), moves,
                                       retired);
                    CHECK(moves - m0 == (uint64_t)n_in);
                }
        MESSAGE(tier.name << ": " << moves << " moves into right-shifted slots");
        CHECK(bad == 0);
    }
}

TEST_CASE("lane compaction: the default setting, the band and batch-size gates and the group clamp"
          * doctest::test_suite("unit/bandedswa")) {
    const ExtScoring sc = scoring(1, 4, 6, 1, 6, 1, 100, 5);
    for (const Tier &tier : compacting_tiers()) {
        INFO("tier " << tier.name);
        std::mt19937_64 rng(0xDEFA0175ull);
        const std::vector<ExtPair> batch = staggered_batch(rng, tier.lanes, 3);
        std::vector<const ExtPair *> all, one_group;
        for (const ExtPair &p : batch) all.push_back(&p);
        for (int q = 0; q < tier.lanes; q++) one_group.push_back(&batch[q]);
        // superblocks one fresh kernel runs on `pairs` at band w, with compaction set to
        // `groups` / `min_w` (groups < -100: the tier's default, untouched)
        auto superblocks = [&](const std::vector<const ExtPair *> &pairs, int w, int groups, int min_w) {
            std::unique_ptr<IBandedPairWiseSW> k(tier.make(sc.o_del, sc.e_del, sc.o_ins, sc.e_ins, sc.zdrop,
                                                           sc.pen_clip, sc.mat, (int8_t)sc.a, (int8_t)sc.b, 1));
            if (groups > -100) k->set_lane_compaction(groups, min_w);
            (void)bwa_tests::score_bsw_batch(*k, pairs, 8, w);
            return k->lane_compaction_counts();
        };
#if defined(__aarch64__)
        const bool on_by_default = bsw8_row_lean_enabled() != 0;   // off on Apple silicon
        const int default_groups = on_by_default ? BSW_COMPACT_GROUPS_NEON : 0;
#else
        const bool on_by_default = true;
        const int default_groups = tier.name == "avx2" ? BSW_COMPACT_GROUPS_AVX2 : BSW_COMPACT_GROUPS_AVX512;
#endif
        // the batch is 3 groups: one superblock at any default of 3 or more, none when off
        const uint64_t default_superblocks = default_groups > 0 ? (uint64_t)((3 + default_groups - 1) / default_groups) : 0;
        MESSAGE(tier.name << ": compaction " << std::string(on_by_default ? "on" : "off") << " by default");
        CHECK(superblocks(all, 31, -1000, 0).superblocks == default_superblocks);
        CHECK(superblocks(all, 30, -1000, 0).superblocks == 0);          // w <= 30 keeps the plain wrapper
        CHECK(superblocks(one_group, 100, 8, 0).superblocks == 0);       // one group: nothing to compact
        CHECK(superblocks(all, 100, 0, 0).superblocks == 0);
        CHECK(superblocks(all, 100, -5, 0).superblocks == 0);
        CHECK(superblocks(all, 40, 8, 41).superblocks == 0);
        CHECK(superblocks(all, 41, 8, 41).superblocks == 1);
        CHECK(superblocks(all, 100, 1, 0).superblocks == 3);             // one group per superblock
        CHECK(superblocks(all, 100, 1, 0).moves == 0);
        CHECK(superblocks(all, 100, 2, 0).superblocks == 2);
        CHECK(superblocks(all, 100, 1000, 0).superblocks == 1);          // capped at BSW_COMPACT_GROUPS_MAX
        // every field still equals the plain path at the cap and at one group
        uint64_t moves = 0, retired = 0;
        CHECK(check_batch(tier, sc, batch, 100, 1000, 0, "groups capped", moves, retired) == 0);
        CHECK(check_batch(tier, sc, batch, 100, 1, 0, "one group", moves, retired) == 0);
        // the default group count itself: a 13-group batch runs ceil(13 / K) superblocks,
        // which tells apart the shipped defaults (4, 6, 8) and their neighbours
        const std::vector<ExtPair> big = staggered_batch(rng, tier.lanes, 13);
        std::vector<const ExtPair *> all13;
        for (const ExtPair &p : big) all13.push_back(&p);
        const uint64_t default_superblocks13 =
            default_groups > 0 ? (uint64_t)((13 + default_groups - 1) / default_groups) : 0;
        CHECK(superblocks(all13, 31, -1000, 0).superblocks == default_superblocks13);
    }
}

TEST_CASE("lane compaction: the plan evacuates the emptiest group into the first holes"
          * doctest::test_suite("unit/bandedswa")) {
    BswMove mv[4 * 64];
    {   // 16 lanes: 2 + 3 + 16 live in 3 groups fit in 2; group 0 (fewest) empties into group 1
        uint64_t live[3] = {0x0011, 0x0007, 0xFFFF};
        uint8_t active[3] = {1, 1, 1};
        const int nm = bsw_compact_plan(live, active, 3, 16, mv);
        REQUIRE(nm == 2);
        CHECK((mv[0].sv == 0 && mv[0].sl == 0 && mv[0].dv == 1 && mv[0].dl == 3));
        CHECK((mv[1].sv == 0 && mv[1].sl == 4 && mv[1].dv == 1 && mv[1].dl == 4));
        CHECK((active[0] == 0 && active[1] == 1 && active[2] == 1));
        CHECK((live[0] == 0 && live[1] == 0x001F && live[2] == 0xFFFF));
    }
    {   // a group with no live lane retires without moves; ties go to the lower index
        uint64_t live[4] = {0, 0x0101, 0x0300, 0xFFFE};
        uint8_t active[4] = {1, 1, 1, 1};
        const int nm = bsw_compact_plan(live, active, 4, 16, mv);
        // 2 + 2 + 15 = 19 live need 2 groups of the 3 active: group 1 (2 live, tied with
        // group 2, lower index) empties into group 2's first holes, lanes 0 and 1
        REQUIRE(nm == 2);
        CHECK(active[0] == 0);
        CHECK(active[1] == 0);
        CHECK((mv[0].sv == 1 && mv[0].sl == 0 && mv[0].dv == 2 && mv[0].dl == 0));
        CHECK((mv[1].sv == 1 && mv[1].sl == 8 && mv[1].dv == 2 && mv[1].dl == 1));
        CHECK((live[2] == 0x0303 && live[3] == 0xFFFE));
    }
    {   // nothing to do when the live lanes need every active group (28, then 17 live)
        uint64_t live[2] = {0xFFFF, 0x0FFF};
        uint8_t active[2] = {1, 1};
        CHECK(bsw_compact_plan(live, active, 2, 16, mv) == 0);
        uint64_t live2[2] = {0x01FF, 0xFF00};
        CHECK(bsw_compact_plan(live2, active, 2, 16, mv) == 0);
        CHECK((active[0] == 1 && active[1] == 1));
        // 8 + 8 = 16 live fit in one group: group 0 (tie, lower index) moves into group 1
        uint64_t live3[2] = {0x00FF, 0xFF00};
        CHECK(bsw_compact_plan(live3, active, 2, 16, mv) == 8);
        CHECK((live3[0] == 0 && live3[1] == 0xFFFF && active[0] == 0));
    }
    {   // 32 lanes (AVX2): 31 + 1 + 32 = 64 live in 3 groups need 2: group 1 fills group 0's
        // one hole, lane 31 (group 2 is full: its 32 set bits are the whole 32-lane mask)
        uint64_t live[3] = {0x7FFFFFFFull, 0x00010000ull, 0xFFFFFFFFull};
        uint8_t active[3] = {1, 1, 1};
        const int nm = bsw_compact_plan(live, active, 3, 32, mv);
        REQUIRE(nm == 1);
        CHECK((mv[0].sv == 1 && mv[0].sl == 16 && mv[0].dv == 0 && mv[0].dl == 31));
        CHECK((live[0] == 0xFFFFFFFFull && live[1] == 0 && live[2] == 0xFFFFFFFFull));
        CHECK((active[0] == 1 && active[1] == 0 && active[2] == 1));
    }
    {   // 64 lanes, 63 + 1 + 1 + 1 = 66 live in 4 groups need 2: group 1 fills group 0's one
        // hole, then group 2 (group 0 now full) moves into group 3
        uint64_t live[4] = {~(uint64_t)0 >> 1, 1, (uint64_t)1 << 63, 2};
        uint8_t active[4] = {1, 1, 1, 1};
        const int nm = bsw_compact_plan(live, active, 4, 64, mv);
        REQUIRE(nm == 2);
        CHECK((mv[0].sv == 1 && mv[0].sl == 0 && mv[0].dv == 0 && mv[0].dl == 63));
        CHECK((mv[1].sv == 2 && mv[1].sl == 63 && mv[1].dv == 3 && mv[1].dl == 0));
        CHECK(live[0] == ~(uint64_t)0);
        CHECK((active[0] + active[1] + active[2] + active[3]) == 2);
        CHECK(__builtin_popcountll(live[0]) + __builtin_popcountll(live[3]) == 66);
    }
}

TEST_CASE("lane compaction: tiers without it ignore the setting"
          * doctest::test_suite("unit/bandedswa")) {
#if defined(__x86_64__) || defined(__i386__)
    __builtin_cpu_init();
    struct Plain { const char *name; Factory make; bool ok; };
    const Plain tiers[] = {{"sse41", make_bsw_kernel_sse41, __builtin_cpu_supports("sse4.1") != 0},
                           {"sse42", make_bsw_kernel_sse42, __builtin_cpu_supports("sse4.2") != 0},
                           {"avx", make_bsw_kernel_avx, __builtin_cpu_supports("avx") != 0}};
    std::mt19937_64 rng(0x0FFC0A1Eull);
    const ExtScoring sc = scoring(1, 4, 6, 1, 6, 1, 100, 5);
    for (const Plain &t : tiers) {
        if (!t.ok) continue;
        INFO("tier " << t.name);
        const Tier tier{t.name, t.make, 32};
        uint64_t moves = 0, retired = 0;
        const std::vector<ExtPair> batch = staggered_batch(rng, 32, 3);
        CHECK(check_batch(tier, sc, batch, 100, 3, 0, "no-op", moves, retired) == 0);
        std::vector<const ExtPair *> ptrs;
        for (const ExtPair &p : batch) ptrs.push_back(&p);
        const std::unique_ptr<IBandedPairWiseSW> k = make(t.make, sc, 3, 0);
        (void)bwa_tests::score_bsw_batch(*k, ptrs, 8, 100);
        CHECK(k->lane_compaction_counts().superblocks == 0);   // the compaction driver never ran
        CHECK(moves == 0);
    }
#else
    MESSAGE("lane compaction: every aarch64 build has the driver; nothing to check here");
#endif
}

#endif // HAVE_BSW_VECTOR_8_16

