// test/unit/test_chain_cap.cpp — the score-gated chain-extension cap.
//
// Covers the two functions that PRODUCE the gate's pruning signal (test_extend_csub.cpp
// covers only the downstream MAPQ seed):
//
//   mem_chain_cap_extend  — the --max-extend-chains count cap + the --extend-tie-frac
//                           competitiveness gate + the --extend-tie-floor "always keep
//                           top-N" floor. Trims a[] in place, returns the kept count, and
//                           reports the MAX weight among dropped chains via capped_w_out
//                           (the seed --extend-csub feeds into MAPQ).
//   mem_capped_pair_subo  — the paired-end competitor-score floor: per-arm
//                           capped_w*score/span, summed, clamped strictly below o so a
//                           pruned pair is never forced to MAPQ 0.
//
// Invariants pinned here: the gate-OFF early return (tie_frac<=0, no count cap) leaves
// chains untouched with capped_w==0 (the byte-identity-off contract); the best chain
// always survives for tie_frac in [0,1] (gate = floor(tie_frac*w_best) <= w_best);
// tie_floor keeps exactly the top-N regardless of the fraction; capped_w is the MAX (not
// min/last) dropped weight; and the PE floor clamps to o-1. Every expectation is derived
// from the contract; fixtures are built in memory, no data files.

#include <cstring>
#include <vector>

#include "doctest/doctest.h"
#include "bwamem.h"

namespace {

// A chain with an explicit weight and a caller-chosen identity (stored in the otherwise-
// unused `pos` field). mem_chain_cap_extend uses a[i].w directly when > 0 (no seed access)
// and never reads `pos`, so `pos` is a clean marker that lets a test tell WHICH of two
// equal-weight chains survived a tie. It drops only free seeds when m > SEEDS_PER_CHAIN —
// m==0 here, so no allocation is needed.
mem_chain_t chain(uint32_t w, int64_t id = 0) {
    mem_chain_t c;
    memset(&c, 0, sizeof(c));
    c.w = w;
    c.pos = id;
    return c;
}

struct CapResult {
    int n;
    int capped_w;
    std::vector<uint32_t> weights;  // surviving weights, in SURVIVOR ORDER (not sorted):
                                    // mem_chain_cap_extend must preserve input order, and
                                    // a stable top-N tie must keep the earlier input chain.
    std::vector<int64_t> ids;       // survivors' identities (`pos`), in the same order.
};

CapResult run_cap(std::vector<mem_chain_t> chains, int max_n, float tie_frac, int tie_floor) {
    int capped_w = -999;  // sentinel; the function must overwrite it
    int n = mem_chain_cap_extend(chains.data(), (int)chains.size(), max_n, tie_frac, tie_floor,
                                 &capped_w);
    CapResult r;
    r.n = n;
    r.capped_w = capped_w;
    for (int i = 0; i < n; i++) { r.weights.push_back(chains[i].w); r.ids.push_back(chains[i].pos); }
    return r;
}

// One paired-end arm: a single chosen region [qb,qe) with score, plus the vector-level
// capped_w that mem_capped_pair_subo reads.
struct PairArm {
    mem_alnreg_t regs[1];
    mem_alnreg_v v;
    PairArm(int qb, int qe, int score, int capped_w) {
        memset(regs, 0, sizeof(regs));
        regs[0].qb = qb;
        regs[0].qe = qe;
        regs[0].score = score;
        memset(&v, 0, sizeof(v));
        v.a = regs;
        v.n = 1;
        v.m = 1;
        v.capped_w = capped_w;
    }
};

}  // namespace

TEST_CASE("cap: off-state (tie_frac=0, no count cap) keeps every chain, capped_w=0"
          * doctest::test_suite("unit/chain_cap")) {
    // This is the --extend-tie-frac 0 / default path: byte-identical to baseline, nothing dropped.
    auto r = run_cap({chain(100), chain(50), chain(10)}, /*max_n=*/0, /*tie_frac=*/0.0f, /*tie_floor=*/0);
    CHECK(r.n == 3);
    CHECK(r.capped_w == 0);
}

TEST_CASE("cap: the count cap keeps the top-N by weight and reports max dropped weight"
          * doctest::test_suite("unit/chain_cap")) {
    // max_n=2 keeps the two highest (100, 50); drops {30, 10}; capped_w = max dropped = 30.
    auto r = run_cap({chain(100), chain(30), chain(50), chain(10)}, /*max_n=*/2, 0.0f, 0);
    CHECK(r.n == 2);
    CHECK(r.weights == std::vector<uint32_t>{100, 50});  // survivor order = input order
    CHECK(r.capped_w == 30);
}

TEST_CASE("cap: a count-cap tie keeps the earlier input chain (stable top-N order)"
          * doctest::test_suite("unit/chain_cap")) {
    // Two equal-weight chains (50) below the best (100); max_n=2 keeps the best plus exactly
    // one of the tied pair. The stable (weight desc, index asc) selection must keep the
    // EARLIER input chain — id=2, not id=3. Sorting the survivors would hide a regression here.
    auto r = run_cap({chain(100, /*id=*/1), chain(50, /*id=*/2), chain(50, /*id=*/3)},
                     /*max_n=*/2, 0.0f, 0);
    CHECK(r.n == 2);
    CHECK(r.weights == std::vector<uint32_t>{100, 50});
    CHECK(r.ids == std::vector<int64_t>{1, 2});  // earlier tied chain survives, in input order
    CHECK(r.capped_w == 50);
}

TEST_CASE("cap: the competitiveness gate prunes non-competitive chains with no count cap"
          * doctest::test_suite("unit/chain_cap")) {
    // best=100, gate = floor(0.9*100) = 90. Keep w>=90 (100, 95); drop 80, 50.
    auto r = run_cap({chain(100), chain(95), chain(80), chain(50)}, /*max_n=*/0, 0.9f, /*tie_floor=*/0);
    CHECK(r.n == 2);
    CHECK(r.weights == std::vector<uint32_t>{100, 95});
    CHECK(r.capped_w == 80);  // MAX of {80, 50}
}

TEST_CASE("cap: a chain whose weight equals the gate threshold is kept (>=, not >)"
          * doctest::test_suite("unit/chain_cap")) {
    // best=100, gate = floor(0.5*100) = 50. Weight 50 must survive; 49 must drop.
    auto r = run_cap({chain(100), chain(50), chain(49)}, 0, 0.5f, 0);
    CHECK(r.weights == std::vector<uint32_t>{100, 50});
    CHECK(r.capped_w == 49);
}

TEST_CASE("cap: tie_floor keeps exactly the top-N; lower-ranked chains still face the gate"
          * doctest::test_suite("unit/chain_cap")) {
    // best=100, gate=90. tie_floor=2 forces ranks 0,1 (100, 80) to stay; rank 2 (50) fails the gate.
    auto kept2 = run_cap({chain(100), chain(80), chain(50)}, 0, 0.9f, /*tie_floor=*/2);
    CHECK(kept2.n == 2);
    CHECK(kept2.weights == std::vector<uint32_t>{100, 80});
    CHECK(kept2.capped_w == 50);
    // tie_floor=3 keeps all three regardless of the gate; nothing dropped.
    auto kept3 = run_cap({chain(100), chain(80), chain(50)}, 0, 0.9f, /*tie_floor=*/3);
    CHECK(kept3.n == 3);
    CHECK(kept3.capped_w == 0);
}

TEST_CASE("cap: the best chain always survives (n=1, and tie_frac=1 keeps only exact-weight ties)"
          * doctest::test_suite("unit/chain_cap")) {
    auto one = run_cap({chain(42)}, /*max_n=*/0, /*tie_frac=*/1.0f, /*tie_floor=*/0);
    CHECK(one.n == 1);
    CHECK(one.weights == std::vector<uint32_t>{42});
    CHECK(one.capped_w == 0);
    // tie_frac=1.0 -> gate = w_best; only chains tied at the best weight survive.
    auto ties = run_cap({chain(100), chain(99), chain(100)}, 0, 1.0f, 0);
    CHECK(ties.n == 2);
    CHECK(ties.weights == std::vector<uint32_t>{100, 100});
    CHECK(ties.capped_w == 99);
}

TEST_CASE("cap: capped_w is the MAX dropped weight, not the last or the min"
          * doctest::test_suite("unit/chain_cap")) {
    // max_n=1 keeps only the best (100); drops {70, 90, 20} in that order — max is 90, mid-list.
    auto r = run_cap({chain(100), chain(70), chain(90), chain(20)}, /*max_n=*/1, 0.0f, 0);
    CHECK(r.n == 1);
    CHECK(r.weights == std::vector<uint32_t>{100});
    CHECK(r.capped_w == 90);
}

TEST_CASE("pair_subo: no capped_w on either arm returns 0 (leaves subo unchanged)"
          * doctest::test_suite("unit/chain_cap")) {
    PairArm a0(0, 100, 80, /*capped_w=*/0), a1(0, 100, 80, 0);
    mem_alnreg_v arms[2] = {a0.v, a1.v};
    int z[2] = {0, 0};
    CHECK(mem_capped_pair_subo(arms, z, /*o=*/160) == 0);
}

TEST_CASE("pair_subo: sums per-arm capped_w*score/span across both mates"
          * doctest::test_suite("unit/chain_cap")) {
    // arm0: 40*80/100 = 32 ; arm1: 20*60/100 = 12 ; total 44 (o large, no clamp).
    PairArm a0(0, 100, 80, 40), a1(0, 100, 60, 20);
    mem_alnreg_v arms[2] = {a0.v, a1.v};
    int z[2] = {0, 0};
    CHECK(mem_capped_pair_subo(arms, z, /*o=*/1000) == 44);
}

TEST_CASE("pair_subo: clamps the floor strictly below o so q_pe is never forced to 0"
          * doctest::test_suite("unit/chain_cap")) {
    // arm0 alone estimates 200*80/100 = 160; with o=100 it must clamp to 99.
    PairArm a0(0, 100, 80, 200), a1(0, 100, 80, 0);
    mem_alnreg_v arms[2] = {a0.v, a1.v};
    int z[2] = {0, 0};
    CHECK(mem_capped_pair_subo(arms, z, /*o=*/100) == 99);
}

TEST_CASE("pair_subo: an out-of-bounds chosen index skips that arm"
          * doctest::test_suite("unit/chain_cap")) {
    PairArm a0(0, 100, 80, 40), a1(0, 100, 80, 40);
    mem_alnreg_v arms[2] = {a0.v, a1.v};
    int z_hi[2] = {0, 5};  // z[1]=5 >= n(1) -> arm1 skipped, only arm0's 32 counts
    CHECK(mem_capped_pair_subo(arms, z_hi, 1000) == 32);
    int z_neg[2] = {-1, 0};  // z[0]<0 -> arm0 skipped, only arm1's 32 counts
    CHECK(mem_capped_pair_subo(arms, z_neg, 1000) == 32);
}

TEST_CASE("pair_subo: an arm with zero span or non-positive score contributes nothing"
          * doctest::test_suite("unit/chain_cap")) {
    PairArm a0(50, 50, 80, 40) /*span 0*/, a1(0, 100, 0, 40) /*score 0*/;
    mem_alnreg_v arms[2] = {a0.v, a1.v};
    int z[2] = {0, 0};
    CHECK(mem_capped_pair_subo(arms, z, 1000) == 0);
}
