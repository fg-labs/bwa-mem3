// test/unit/test_extend_csub.cpp — mem_seed_capped_sub (--extend-csub).
//
// When --max-extend-chains / --extend-tie-frac prunes chains before banded-SW, the
// pruned competitor is never scored, so it never reaches `sub` in
// mem_approx_mapq_se and MAPQ comes out overstated. mem_seed_capped_sub seeds it
// back in from the max WEIGHT among the dropped chains (mem_alnreg_v::capped_w).
//
// It seeds the `sub` (2nd-best score) slot, NOT `csub`. `csub` is the tandem-hit
// slot, and the PE MAPQ recombination caps on it mate-BLIND *after* the +40 pair
// promotion -- routing a dropped (mate-non-concordant, under the mate-aware --fast
// cap) chain there over-demotes correct paired reads in paralog/pseudogene regions,
// where the mate could disambiguate but the csub cap ignores it. `sub` feeds the
// promotable q_se, so pairing rescues reads the mate resolves. SE MAPQ and XS both
// read max(sub,csub), so single-end output is byte-identical to seeding csub.
//
// The estimate is the all-match upper bound `capped_w * match_a` -- a dropped chain
// of weight w could not have scored more than w*a. That bound is deliberately
// generous, which is why the CLAMP is the load-bearing part of this function:
//
//   mem_approx_mapq_se does `sub = max(sub, csub); if (sub >= score) return 0;`
//
// so an unclamped estimate at/above the region's own score forces MAPQ 0. An
// earlier unclamped version did exactly that, and the PE recombination of two
// zeroed mates then produced spurious MAPQ60 *promotions* -- measured 7x worse
// than doing nothing. Clamping to score-1 turns the same signal into a
// proportional MAPQ reduction instead. The cases below pin that clamp, the seed
// targeting `sub` (never `csub`), the off-state no-op that keeps --extend-csub
// byte-identical when nothing was dropped, and the primary-only / max-not-overwrite
// semantics.
//
// Every expectation is derived from the function's contract, not recorded from a
// run. Fixtures are built in memory; no test data files are read.

#include <cstring>
#include <vector>

#include "doctest/doctest.h"
#include "bwamem.h"

namespace {

// A primary region spanning [qb,qe) with the given score and pre-existing sub/csub.
// secondary < 0 marks it primary (what mem_mark_primary_se leaves behind).
mem_alnreg_t reg(int qb, int qe, int score, int sub = 0, int csub = 0, int secondary = -1) {
    mem_alnreg_t a;
    memset(&a, 0, sizeof(a));
    a.qb = qb; a.qe = qe;
    a.score = a.truesc = score;
    a.sub = sub;
    a.csub = csub;
    a.secondary = secondary;
    a.meth_hypothesis = a.meth_strand_hyp = -1;
    return a;
}

// Wrap regions in the mem_alnreg_v the function takes, with capped_w set.
struct RegV {
    std::vector<mem_alnreg_t> storage;
    mem_alnreg_v v;
    RegV(std::vector<mem_alnreg_t> regs, int capped_w) : storage(std::move(regs)) {
        memset(&v, 0, sizeof(v));
        v.a = storage.data();
        v.n = storage.size();
        v.m = storage.size();
        v.capped_w = capped_w;
    }
};

const int MATCH_A = 1;  // opt->a default

}  // namespace

TEST_CASE("mem_seed_capped_sub: no-op when nothing was dropped"
          * doctest::test_suite("unit/extend_csub")) {
    // capped_w <= 0 is the state whenever the cap/gate is off, or is on but
    // dropped nothing. This is what keeps --extend-csub byte-identical on the
    // default path and on --fast with the gate off.
    SUBCASE("capped_w == 0 (cap/gate off, or dropped nothing)") {
        RegV r({reg(0, 100, 80, /*sub=*/7)}, /*capped_w=*/0);
        mem_seed_capped_sub(&r.v, MATCH_A);
        CHECK(r.v.a[0].sub == 7);
    }
    SUBCASE("capped_w < 0 (never-set sentinel)") {
        RegV r({reg(0, 100, 80, /*sub=*/7)}, /*capped_w=*/-1);
        mem_seed_capped_sub(&r.v, MATCH_A);
        CHECK(r.v.a[0].sub == 7);
    }
}

TEST_CASE("mem_seed_capped_sub: seeds sub, never csub (the mate-promotable slot)"
          * doctest::test_suite("unit/extend_csub")) {
    // The fix: the pruned-competitor estimate goes to `sub` (promotable by PE
    // pairing), and `csub` (the mate-blind tandem cap) is left untouched.
    RegV r({reg(0, 100, 80, /*sub=*/0, /*csub=*/7)}, /*capped_w=*/30);
    mem_seed_capped_sub(&r.v, MATCH_A);
    CHECK(r.v.a[0].sub == 30);   // estimate landed in sub
    CHECK(r.v.a[0].csub == 7);   // csub is NOT the target -- pre-existing value preserved
}

TEST_CASE("mem_seed_capped_sub: below the score, the estimate is used as-is"
          * doctest::test_suite("unit/extend_csub")) {
    // weight 30 * a 1 = 30 < score 80, so no clamping is needed and sub takes
    // the raw all-match bound.
    RegV r({reg(0, 100, 80)}, /*capped_w=*/30);
    mem_seed_capped_sub(&r.v, MATCH_A);
    CHECK(r.v.a[0].sub == 30);
}

TEST_CASE("mem_seed_capped_sub: clamps strictly below the region score"
          * doctest::test_suite("unit/extend_csub")) {
    // The regression guard. For every dropped weight at or above the region's
    // score, sub must land at exactly score-1 -- never >= score, which would
    // trip `sub >= score -> return 0` in mem_approx_mapq_se.
    const int score = 80;
    int w = 0;
    SUBCASE("weight == score")            { w = 80; }
    SUBCASE("weight one above score")     { w = 81; }
    SUBCASE("weight well above score")    { w = 200; }
    SUBCASE("weight enormously above")    { w = 10000; }
    RegV r({reg(0, 100, score)}, w);
    mem_seed_capped_sub(&r.v, MATCH_A);
    CHECK(r.v.a[0].sub == score - 1);
    CHECK(r.v.a[0].sub < score);
}

TEST_CASE("mem_seed_capped_sub: clamped estimate never forces MAPQ to zero"
          * doctest::test_suite("unit/extend_csub")) {
    // The clamp's actual purpose, stated in terms of the consumer. An enormous
    // dropped weight should REDUCE mapq relative to no pruning, not zero it --
    // zeroing is what produced the spurious PE promotions.
    mem_opt_t *opt = mem_opt_init();

    mem_alnreg_t base = reg(0, 100, 80);
    base.seedcov = 60;
    int mapq_untouched = mem_approx_mapq_se(opt, &base);

    RegV r({base}, /*capped_w=*/10000);
    mem_seed_capped_sub(&r.v, opt->a);
    int mapq_capped = mem_approx_mapq_se(opt, &r.v.a[0]);

    CHECK(mapq_untouched > 0);
    CHECK(mapq_capped > 0);                  // NOT zeroed -- the whole point
    CHECK(mapq_capped < mapq_untouched);     // but genuinely reduced

    free(opt);
}

TEST_CASE("mem_seed_capped_sub: secondary regions are left alone"
          * doctest::test_suite("unit/extend_csub")) {
    // Only primaries (secondary < 0) carry the pruned-competitor estimate;
    // mem_mark_primary_se runs before this, and secondaries already have a
    // meaningful sub relationship to their primary.
    RegV r({reg(0, 100, 80, /*sub=*/3, /*csub=*/0, /*secondary=*/0)}, /*capped_w=*/30);
    mem_seed_capped_sub(&r.v, MATCH_A);
    CHECK(r.v.a[0].sub == 3);
}

TEST_CASE("mem_seed_capped_sub: takes the max, never lowers an existing sub"
          * doctest::test_suite("unit/extend_csub")) {
    // sub may already hold a real scored competitor. The estimate must not erase it.
    SUBCASE("existing sub is larger -> preserved") {
        RegV r({reg(0, 100, 80, /*sub=*/50)}, /*capped_w=*/30);
        mem_seed_capped_sub(&r.v, MATCH_A);
        CHECK(r.v.a[0].sub == 50);
    }
    SUBCASE("estimate is larger -> raised") {
        RegV r({reg(0, 100, 80, /*sub=*/10)}, /*capped_w=*/30);
        mem_seed_capped_sub(&r.v, MATCH_A);
        CHECK(r.v.a[0].sub == 30);
    }
}

TEST_CASE("mem_seed_capped_sub: skips degenerate regions"
          * doctest::test_suite("unit/extend_csub")) {
    // A zero/negative query span or a non-positive score has no meaningful
    // per-region bound to compare against; both are skipped rather than clamped
    // (score-1 would be nonsense for score <= 0).
    SUBCASE("empty query span") {
        RegV r({reg(50, 50, 80, /*sub=*/1)}, /*capped_w=*/30);
        mem_seed_capped_sub(&r.v, MATCH_A);
        CHECK(r.v.a[0].sub == 1);
    }
    SUBCASE("non-positive score") {
        RegV r({reg(0, 100, 0, /*sub=*/1)}, /*capped_w=*/30);
        mem_seed_capped_sub(&r.v, MATCH_A);
        CHECK(r.v.a[0].sub == 1);
    }
}

TEST_CASE("mem_seed_capped_sub: each primary is clamped against its own score"
          * doctest::test_suite("unit/extend_csub")) {
    // capped_w is per-READ, but the clamp is per-REGION, so one dropped weight
    // lands differently on regions of differing score.
    RegV r({reg(0, 100, 80), reg(0, 100, 40), reg(0, 100, 200)}, /*capped_w=*/100);
    mem_seed_capped_sub(&r.v, MATCH_A);
    CHECK(r.v.a[0].sub == 79);   // clamped to its own score-1
    CHECK(r.v.a[1].sub == 39);   // clamped lower
    CHECK(r.v.a[2].sub == 100);  // under score, used as-is
}
