// test/integration/test_bandedswa_compact16.cpp
//
// Same-row lane compaction in getScores16 (src/bandedSWA_compact16.inc on NEON,
// bandedSWA_compact16_256.inc on AVX2, bandedSWA_compact16_512.inc on AVX-512BW). The driver runs
// a batch in superblocks of K lane groups in block lockstep, R rows per block, and at block
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
// Fast clamps, one-group gates and tier defaults remain in test/unit/test_bandedswa_compact16.cpp.
// Removing the receiving slot's clear, or the move's nend, ncol or minq update, makes this file
// fail. (The nbeg update is kept for the same reason as the others, but no batch here, and no
// profile batch of the kernel goldens, has needed it: a group's nbeg trails its lanes' heads.)
// The reference is the same kernel copy with compaction off (set_lane_compaction16(0, 0)),
// batch for batch, each pair scored alone and the independent scalarBandedSWA recurrence.
// Ambiguous-base and tied-maximum batches also compare all six scalar fields after moves.
// Every compacting tier this host runs is checked (batches sized in its lane count). The unit
// tests separately check that noncompacting SSE tiers ignore the setting.

#include "bsw_compact16_test.h"

#if HAVE_BSW_VECTOR_8_16

static void check_staggered_finishes(int selected_rep) {
    for (const Tier &t : compacting_tiers()) {
        INFO("tier " << t.name);
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
                    const std::vector<ExtPair> batch = staggered_batch(rng, t.lanes, groups);
                    if (rep != selected_rep) continue;   // preserve the RNG stream
                    for (int w : {100, 40, 12})
                        bad += check_batch(t, sc, batch, w, groups, rows, "staggered", moves, retired);
                }
            }
        }
        MESSAGE(t.name << " 16-bit: " << moves << " moves, " << retired << " groups retired by evacuation");
        CHECK(moves > 100);
        CHECK(retired > 20);
        CHECK(bad == 0);
    }
}

TEST_CASE("lane compaction 16-bit: staggered finishes, first scoring and matrix set"
          * doctest::test_suite("integration/bandedswa")) {
    check_staggered_finishes(0);
}

TEST_CASE("lane compaction 16-bit: staggered finishes, second scoring and matrix set"
          * doctest::test_suite("integration/bandedswa")) {
    check_staggered_finishes(1);
}

TEST_CASE("lane compaction 16-bit: staggered finishes, third scoring and matrix set"
          * doctest::test_suite("integration/bandedswa")) {
    check_staggered_finishes(2);
}

TEST_CASE("lane compaction 16-bit: ambiguous bases match the scalar oracle"
          * doctest::test_suite("integration/bandedswa")) {
    for (const Tier &t : compacting_tiers()) {
        INFO("tier " << t.name);
        const ExtScoring base = scoring(1, 4, 6, 1, 6, 1, 100, 5);
        for (int matrix : {0, 1, 2}) {
            std::mt19937_64 rng(0xA6B16C7ull);
            const ExtScoring sc = with_matrix(base, matrix);
            std::vector<ExtPair> batch = staggered_batch(rng, t.lanes, 3);
            for (size_t q = 0; q < batch.size(); q++) {
                for (size_t j = q % 7; j < batch[q].query.size(); j += 17) batch[q].query[j] = 4;
                for (size_t j = q % 5; j < batch[q].target.size(); j += 19) batch[q].target[j] = 4;
            }
            for (int w : {12, 40, 100}) {
                uint64_t moves = 0, retired = 0;
                CHECK(check_batch(t, sc, batch, w, 3, 1, "ambiguous bases", moves, retired) == 0);
                CHECK(moves > 0);
            }
        }
    }
}

TEST_CASE("lane compaction 16-bit: tied maxima match the scalar oracle"
          * doctest::test_suite("integration/bandedswa")) {
    for (const Tier &t : compacting_tiers()) {
        INFO("tier " << t.name);
        // A match/mismatch cycle adds 1 then subtracts 1: later rows tie the first maximum.
        const ExtScoring sc = scoring(1, 1, 6, 1, 6, 1, 100, 5);
        std::vector<ExtPair> batch;
        for (int q = 0; q < 2 * t.lanes; q++) {
            const int len = 12 * (q % 16 + 1);
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
            CHECK(check_batch(t, sc, batch, w, 2, 1, "tied maxima", moves, retired) == 0);
            CHECK(moves > 0);
        }
    }
}

TEST_CASE("lane compaction 16-bit: a lane moved into a just-finished lane's slot ignores its stale cells"
          * doctest::test_suite("integration/bandedswa")) {
    // Superblocks of two groups. Group 0: two incoming lanes aligned d columns right of the
    // main diagonal, and short pairs that finish at once. Group 1: dead lanes that score high
    // on the main diagonal and finish at row r0, and long matches. When the dead lanes finish,
    // the live lanes fit in one group, so group 0 is evacuated and the incoming lanes land in
    // the dead lanes' slots, whose cells left of the incoming heads are large: only the move's
    // clear keeps them out of the incoming lanes' horizontal carry.
    for (const Tier &t : compacting_tiers()) {
        INFO("tier " << t.name);
        const size_t L = (size_t)t.lanes;
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
                            while (batch.size() < L) batch.push_back(perfect(rng, 4, 0, 10));
                            for (int q = 0; q < 2; q++) batch.push_back(perfect(rng, r0, 0, 120));
                            while (batch.size() < 2 * L) batch.push_back(perfect(rng, 130, 10, batches % 2 ? 30 : 5));
                            const uint64_t m0 = moves;
                            bad += check_batch(t, sc, batch, w, 2, rows,
                                               "stale cells e=" + std::to_string(e) + " d=" + std::to_string(d) +
                                                   " r0=" + std::to_string(r0),
                                               moves, retired);
                            CHECK(moves > m0);   // the incoming lanes moved
                            batches++;
                        }
        MESSAGE(t.name << " 16-bit: " << batches << " stale-cell batches, " << moves << " moves");
        CHECK(bad == 0);
    }
}

TEST_CASE("lane compaction 16-bit: a lane moved into a group with a narrower column range keeps its own"
          * doctest::test_suite("integration/bandedswa")) {
    // The receiving group's live lanes have short queries (40) on long targets, so its column
    // range (nend, ncol) is far narrower than the incoming lanes'. The incoming lanes have
    // queries of 160-240 and align d columns right of the main diagonal, so right after the
    // move their best cells lie past the receiving group's nend. The move must widen nend and
    // ncol, and lower the group's minq, or the incoming lanes lose those columns.
    for (const Tier &t : compacting_tiers()) {
        INFO("tier " << t.name);
        const size_t L = (size_t)t.lanes;
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
                        while (batch.size() < L) batch.push_back(perfect(rng, 4, 0, 10));
                        for (int q = 0; q < 2; q++) batch.push_back(perfect(rng, 20, 0, 60));
                        while (batch.size() < 2 * L) batch.push_back(perfect(rng, 40, 300, 90));
                        const uint64_t m0 = moves;
                        bad += check_batch(t, sc, batch, w, 2, rows, "narrow receiver len=" + std::to_string(len),
                                           moves, retired);
                        CHECK(moves > m0);
                    }
        CHECK(bad == 0);
    }
}

TEST_CASE("lane compaction 16-bit: a lane moved into a group further right keeps its left columns"
          * doctest::test_suite("integration/bandedswa")) {
    // The receiving group's live lanes align d columns right of the main diagonal, while the
    // incoming lanes stay on the main diagonal. Check that these moves preserve every result
    // field at narrow bands; this batch does not isolate the move's nbeg update.
    for (const Tier &t : compacting_tiers()) {
        INFO("tier " << t.name);
        const size_t L = (size_t)t.lanes;
        std::mt19937_64 rng(0x1EF7C016ull);
        uint64_t moves = 0, retired = 0;
        int bad = 0;
        for (int d : {40, 80})
            for (int w : {20, 40})
                for (int rows : {1, 8}) {
                    const ExtScoring sc = scoring(1, 4, 6, 1, 6, 1, 100, 5);
                    std::vector<ExtPair> batch;
                    for (int q = 0; q < 2; q++) batch.push_back(perfect(rng, 180, 20, 60));
                    while (batch.size() < L) batch.push_back(perfect(rng, 4, 0, 10));
                    for (int q = 0; q < 2; q++) batch.push_back(perfect(rng, 24, 0, 30));
                    while (batch.size() < 2 * L) batch.push_back(shifted(rng, d, 160, 6 + d + 20));
                    const uint64_t m0 = moves;
                    bad += check_batch(t, sc, batch, w, 2, rows, "right receiver d=" + std::to_string(d), moves,
                                       retired);
                    CHECK(moves > m0);
                }
        CHECK(bad == 0);
    }
}

#endif   // HAVE_BSW_VECTOR_8_16
