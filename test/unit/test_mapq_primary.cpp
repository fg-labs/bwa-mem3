// test/unit/test_mapq_primary.cpp -- single-end MAPQ (mem_approx_mapq_se) and
// primary/secondary marking (mem_mark_primary_se), worked by hand.
//
// Every record's MAPQ comes from mem_approx_mapq_se, and which of a read's hits
// is primary, secondary or a supplementary candidate comes from
// mem_mark_primary_se. The existing tests check one clamp of the MAPQ inputs
// (test_extend_csub.cpp) and the single-hit fast path of the marking
// (test_unique_mapper_fastpath.cpp); the formula itself and the general marking
// path (overlap rule, tie order, sub and sub_n, ALT hits) had no test.
//
// Expected values follow bwa's documented model (bwa 0.7.x bwamem.c,
// mem_approx_mapq_se and mem_mark_primary_se):
//   sub = max(sub or min_seed_len * a, csub); 0 when sub >= score; otherwise,
//   with l = max(query span, reference span) and identity
//   1 - (l * a - score) / (a + b) / l, at the default -Q 50:
//   t = (l < 50 ? 1 : int(ln 50) / ln l) * identity^2,
//   MAPQ = (int)(6.02 * (score - sub) / a * t^2 + .499); with -Q 0:
//   (int)(30 * (1 - sub / score) * ln(seedcov) + .499), times identity^2 below
//   0.95 identity; then less (int)(4.343 * ln(sub_n + 1) + .499), clamped to
//   [0, 60], times (1 - frac_rep). int(ln 50) is 3, not 3.912: bwa stores the
//   coefficient in an int field, and the output depends on that truncation.
//   Hits are ranked by score, then by hash_64(read id + input index) (bwa's
//   64-bit integer hash, written out in wang_hash64.h); a hit is
//   secondary to the first higher-ranked primary it overlaps on the read by at
//   least mask_level (0.5) of the shorter span; that primary's sub is the first
//   such hit's score and its sub_n counts those within max(a + b, o + e) of it.
// The arithmetic is written next to each check.

#include "doctest/doctest.h"
#include "bwamem.h"
#include "wang_hash64.h"

#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <vector>

namespace {

struct OptGuard {
    mem_opt_t *o;
    OptGuard() : o(mem_opt_init()) {}
    ~OptGuard() { free(o); }
};

// A region over read bases [qb, qe) and reference [rb, rb + rlen).
mem_alnreg_t region(int score, int qb, int qe, int64_t rb = 0, int rlen = -1) {
    mem_alnreg_t a;
    memset(&a, 0, sizeof(a));
    a.qb = qb;
    a.qe = qe;
    a.rb = rb;
    a.re = rb + (rlen < 0 ? qe - qb : rlen);
    a.score = a.truesc = score;
    a.seedcov = qe - qb;
    a.secondary = a.secondary_all = -1;
    return a;
}

}  // namespace

TEST_CASE("mem_approx_mapq_se: unique hit, sub-optimal second, csub"
          * doctest::test_suite("unit/mapq_primary")) {
    OptGuard g;
    // l = 100: t = 3 / ln 100 = .651442, t^2 = .424376.
    SUBCASE("a unique hit: sub = 19, (int)(6.02 * 81 * .424376 + .499) = 207 -> 60") {
        mem_alnreg_t a = region(100, 0, 100);
        CHECK(mem_approx_mapq_se(g.o, &a) == 60);
    }
    SUBCASE("a second hit at 90: (int)(6.02 * 10 * .424376 + .499) = 26") {
        // With an untruncated ln 50 (3.912) this would be 43.
        mem_alnreg_t a = region(100, 0, 100);
        a.sub = 90;
        CHECK(mem_approx_mapq_se(g.o, &a) == 26);
    }
    SUBCASE("csub above sub is used instead: csub 90, sub 80 -> 26") {
        mem_alnreg_t a = region(100, 0, 100);
        a.sub = 80;
        a.csub = 90;
        CHECK(mem_approx_mapq_se(g.o, &a) == 26);
    }
}

TEST_CASE("mem_approx_mapq_se: sub_n, frac_rep and the cap at 60"
          * doctest::test_suite("unit/mapq_primary")) {
    OptGuard g;
    mem_alnreg_t a = region(100, 0, 100);
    a.sub = 90;   // 26 before the adjustments below
    SUBCASE("sub_n 2: 26 - (int)(4.343 * ln 3 + .499) = 26 - (int)(5.27) = 21") {
        // Truncating 4.343 * ln 3 = 4.77 without the .499 would give 22.
        a.sub_n = 2;
        CHECK(mem_approx_mapq_se(g.o, &a) == 21);
    }
    SUBCASE("frac_rep .25: (int)(26 * .75 + .499) = (int)(19.999) = 19, not 20") {
        a.frac_rep = 0.25f;
        CHECK(mem_approx_mapq_se(g.o, &a) == 19);
    }
    SUBCASE("the cap applies before frac_rep: unique, frac_rep .5 -> (int)(60 * .5 + .499) = 30") {
        mem_alnreg_t u = region(100, 0, 100);
        u.frac_rep = 0.5f;
        CHECK(mem_approx_mapq_se(g.o, &u) == 30);
    }
}

TEST_CASE("mem_approx_mapq_se: identity and alignment length"
          * doctest::test_suite("unit/mapq_primary")) {
    OptGuard g;
    SUBCASE("a reference span longer than the query (110 vs 100) sets l") {
        // l = 110: 3 / ln 110 = .638228; identity 1 - (110 - 90) / 5 / 110 =
        // .963636, t = .638228 * .928595 = .592656, t^2 = .351241; sub 85:
        // (int)(6.02 * 5 * .351241 + .499) = (int)(11.07) = 11.
        mem_alnreg_t a = region(90, 0, 100, 0, 110);
        a.sub = 85;
        CHECK(mem_approx_mapq_se(g.o, &a) == 11);
    }
    SUBCASE("below -Q 50 bp the length factor is 1") {
        // l = 49: t = 1; score 49, sub 44: (int)(6.02 * 5 + .499) = 30.
        mem_alnreg_t a = region(49, 0, 49);
        a.sub = 44;
        CHECK(mem_approx_mapq_se(g.o, &a) == 30);
        // l = 50 is not below 50: t = 3 / ln 50 = .766870, t^2 = .588089;
        // score 50, sub 45: (int)(6.02 * 5 * .588089 + .499) = (int)(18.20) = 18.
        mem_alnreg_t b = region(50, 0, 50);
        b.sub = 45;
        CHECK(mem_approx_mapq_se(g.o, &b) == 18);
    }
}

TEST_CASE("mem_approx_mapq_se: the -Q 0 formula scales by identity^2 below 0.95"
          * doctest::test_suite("unit/mapq_primary")) {
    OptGuard g;
    g.o->mapQ_coef_len = 0;
    g.o->mapQ_coef_fac = 0;
    SUBCASE("identity .94 (score 70 over 100 bp) scales by identity^2") {
        // (int)(30 * (1 - 50/70) * ln 100 + .499) = (int)(39.97) = 39, then
        // (int)(39 * .8836 + .499) = (int)(34.96) = 34.
        mem_alnreg_t a = region(70, 0, 100);
        a.sub = 50;
        CHECK(mem_approx_mapq_se(g.o, &a) == 34);
    }
    SUBCASE("identity .984 (score 92 over 100 bp) is not scaled") {
        // (int)(30 * (1 - 70/92) * ln 100 + .499) = (int)(33.54) = 33; scaled by
        // .984^2 it would be (int)(32.45) = 32.
        mem_alnreg_t a = region(92, 0, 100);
        a.sub = 70;
        CHECK(mem_approx_mapq_se(g.o, &a) == 33);
    }
}

TEST_CASE("mem_mark_primary_se: overlap of at least half the shorter hit makes a secondary"
          * doctest::test_suite("unit/mapq_primary")) {
    // A covers read bases [0, 100) with 100 points, B [s, s + 100) with 90. The
    // overlap 100 - s is significant when >= .5 * 100.
    OptGuard g;
    const struct {
        int s;
        bool secondary;
    } cases[] = {{50, true}, {51, false}};
    for (const auto &c : cases) {
        CAPTURE(c.s);
        std::vector<mem_alnreg_t> a;
        a.push_back(region(90, c.s, c.s + 100, 5000));
        a.push_back(region(100, 0, 100, 100));
        const int n_pri = mem_mark_primary_se(g.o, 2, a.data(), 42);
        CHECK(n_pri == 2);   // the number of non-ALT hits, not of primaries
        CHECK(a[0].score == 100);
        CHECK(a[0].secondary == -1);
        CHECK(a[1].score == 90);
        CHECK(a[1].secondary == (c.secondary ? 0 : -1));
        CHECK(a[1].secondary_all == a[1].secondary);
        // A's sub is B's score only when B is its secondary; sub_n stays 0
        // because 100 - 90 > 7.
        CHECK(a[0].sub == (c.secondary ? 90 : 0));
        CHECK(a[0].sub_n == 0);
    }
}

TEST_CASE("mem_mark_primary_se: sub is the first overlapping hit, sub_n counts those within 7"
          * doctest::test_suite("unit/mapq_primary")) {
    // Four hits on the same read bases, scores 100, 95, 93 and 92 (input order
    // shuffled): 95 is the first one under 100, so sub = 95; within
    // max(a + b, o + e) = 7 of 100 are 95 and 93 (7 below, inclusive), not 92:
    // sub_n = 2. All three point at the primary.
    OptGuard g;
    std::vector<mem_alnreg_t> a;
    a.push_back(region(93, 0, 100, 3000));
    a.push_back(region(92, 0, 100, 4000));
    a.push_back(region(100, 0, 100, 1000));
    a.push_back(region(95, 0, 100, 2000));
    mem_mark_primary_se(g.o, 4, a.data(), 7);
    CHECK(a[0].score == 100);
    CHECK(a[1].score == 95);
    CHECK(a[2].score == 93);
    CHECK(a[3].score == 92);
    CHECK(a[0].secondary == -1);
    CHECK(a[0].sub == 95);
    CHECK(a[0].sub_n == 2);
    for (int i = 1; i < 4; ++i) {
        CAPTURE(i);
        CHECK(a[i].secondary == 0);
    }
}

TEST_CASE("mem_mark_primary_se: a chain of overlaps marks against primaries only"
          * doctest::test_suite("unit/mapq_primary")) {
    // A [0, 60) 60 points, B [60, 100) 40, C [0, 60) 58, D [40, 100) 39.
    // A is primary; B overlaps A by 0 -> primary (a supplementary candidate);
    // C overlaps A fully -> secondary of A; D overlaps A by 20 < .5 * 60 but B by
    // 40 >= .5 * 40 -> secondary of B. A.sub = 58 (within 7: sub_n 1);
    // B.sub = 39 (within 7: sub_n 1).
    OptGuard g;
    std::vector<mem_alnreg_t> a;
    a.push_back(region(40, 60, 100, 9000));
    a.push_back(region(39, 40, 100, 7000));
    a.push_back(region(58, 0, 60, 5000));
    a.push_back(region(60, 0, 60, 1000));
    mem_mark_primary_se(g.o, 4, a.data(), 3);
    REQUIRE(a[0].score == 60);
    REQUIRE(a[1].score == 58);
    REQUIRE(a[2].score == 40);
    REQUIRE(a[3].score == 39);
    CHECK(a[0].secondary == -1);
    CHECK(a[1].secondary == 0);
    CHECK(a[2].secondary == -1);
    CHECK(a[3].secondary == 2);
    CHECK(a[0].sub == 58);
    CHECK(a[0].sub_n == 1);
    CHECK(a[2].sub == 39);
    CHECK(a[2].sub_n == 1);
}

TEST_CASE("mem_mark_primary_se: a hit overlapping two primaries is the secondary of the higher-ranked one"
          * doctest::test_suite("unit/mapq_primary")) {
    // A [0, 60) 60 points and B [60, 100) 40 do not overlap: both primary. C
    // [30, 90) 35 overlaps A by 30 >= .5 * min(60, 60) and B by 30 >= .5 *
    // min(60, 40): significant for both, so it is the secondary of the first
    // primary in rank order, A. A.sub = 35 (60 - 35 > 7: sub_n 0); B keeps sub
    // 0, so B has MAPQ (int)(6.02 * (40 - 19) + .499) = 126 -> 60 (l = 40 < 50:
    // t = 1); as C's primary it would have had (int)(6.02 * 5 + .499) = 30. A
    // has (int)(6.02 * 25 * (3 / ln 60)^2 + .499) = (int)(81.30) -> 60.
    OptGuard g;
    std::vector<mem_alnreg_t> a;
    a.push_back(region(35, 30, 90, 5000));
    a.push_back(region(40, 60, 100, 3000));
    a.push_back(region(60, 0, 60, 1000));
    const int n_pri = mem_mark_primary_se(g.o, 3, a.data(), 11);
    CHECK(n_pri == 3);
    REQUIRE(a[0].score == 60);
    REQUIRE(a[1].score == 40);
    REQUIRE(a[2].score == 35);
    CHECK(a[0].secondary == -1);
    CHECK(a[1].secondary == -1);
    CHECK(a[2].secondary == 0);
    CHECK(a[0].sub == 35);
    CHECK(a[0].sub_n == 0);
    CHECK(a[1].sub == 0);
}

TEST_CASE("mem_mark_primary_se: equal scores are ranked by hash_64(read id + input index)"
          * doctest::test_suite("unit/mapq_primary")) {
    // Two hits with the same score on the same read bases. The primary is the one
    // with the smaller hash_64(id + i), i its input index; over many read ids
    // each input position must win sometimes, so a tie broken by input order (or
    // by a different hash) disagrees for some id.
    OptGuard g;
    int first_won = 0, second_won = 0;
    for (int64_t id = 0; id < 64; ++id) {
        CAPTURE(id);
        std::vector<mem_alnreg_t> a;
        a.push_back(region(100, 0, 100, 1000));
        a.push_back(region(100, 0, 100, 2000));
        mem_mark_primary_se(g.o, 2, a.data(), id);
        const bool first_wins = bwa_tests::wang_hash64((uint64_t)id) < bwa_tests::wang_hash64((uint64_t)id + 1);
        CHECK(a[0].rb == (first_wins ? 1000 : 2000));
        CHECK(a[0].secondary == -1);
        CHECK(a[1].secondary == 0);
        (first_wins ? first_won : second_won)++;
    }
    CHECK(first_won > 0);
    CHECK(second_won > 0);
}

TEST_CASE("mem_mark_primary_se: ALT hits do not compete with primary-assembly hits"
          * doctest::test_suite("unit/mapq_primary")) {
    // A primary-assembly hit A and an ALT hit on the same read bases. Hits are
    // first ranked all together (secondary_all records that ranking), then
    // re-sorted with primary-assembly hits first and marked among those alone.
    // A has no primary-assembly competitor, so it is primary with sub 0
    // (single-end MAPQ 60) whatever the ALT hit scores. The return value is the
    // number of primary-assembly hits.
    OptGuard g;
    SUBCASE("the ALT hit scores higher (100 vs 90)") {
        // In the all-hit ranking the ALT hit leads and A is its secondary: A keeps
        // the ALT score in alt_sc and points at it (index 1 after the re-sort);
        // the ALT hit, secondary to nothing, has secondary and secondary_all -1.
        std::vector<mem_alnreg_t> a;
        a.push_back(region(100, 0, 100, 2000));
        a[0].is_alt = 1;
        a.push_back(region(90, 0, 100, 1000));
        const int n_pri = mem_mark_primary_se(g.o, 2, a.data(), 5);
        CHECK(n_pri == 1);
        REQUIRE((int)a[0].is_alt == 0);
        CHECK(a[0].score == 90);
        CHECK(a[0].secondary == -1);
        CHECK(a[0].sub == 0);
        CHECK(a[0].alt_sc == 100);
        CHECK(a[0].secondary_all == 1);
        CHECK(mem_approx_mapq_se(g.o, &a[0]) == 60);
        CHECK((int)a[1].is_alt == 1);
        CHECK(a[1].secondary == -1);
        CHECK(a[1].secondary_all == -1);
    }
    SUBCASE("the ALT hit scores lower (95 vs 100)") {
        // The ALT hit is A's secondary in both rankings and is marked INT_MAX (an
        // ALT secondary: listed, never emitted on its own). A's sub, set to 95 in
        // the all-hit round, is cleared in the primary-assembly round. sub_n is
        // not cleared there, and an ALT hit within 7 of a primary-assembly hit
        // does not count toward it: 0.
        std::vector<mem_alnreg_t> a;
        a.push_back(region(95, 0, 100, 2000));
        a[0].is_alt = 1;
        a.push_back(region(100, 0, 100, 1000));
        const int n_pri = mem_mark_primary_se(g.o, 2, a.data(), 5);
        CHECK(n_pri == 1);
        REQUIRE((int)a[0].is_alt == 0);
        CHECK(a[0].score == 100);
        CHECK(a[0].secondary == -1);
        CHECK(a[0].sub == 0);
        CHECK(a[0].sub_n == 0);
        CHECK(a[0].alt_sc == 0);
        CHECK(a[0].secondary_all == -1);
        CHECK(mem_approx_mapq_se(g.o, &a[0]) == 60);
        CHECK((int)a[1].is_alt == 1);
        CHECK(a[1].secondary == INT_MAX);
        CHECK(a[1].secondary_all == 0);
    }
}
