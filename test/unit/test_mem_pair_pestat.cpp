// test/unit/test_mem_pair_pestat.cpp -- the insert-size estimator (mem_pestat) and
// the pair orientation (mem_infer_dir) against bwa's documented model.
//
// mem_pestat decides, once per batch, which orientations can pair and the insert
// window [low, high] inside which a pair is proper; every 0x2 bit, every pair
// score and every rescue window downstream reads it. It counts the candidate
// inserts into a histogram and walks the counts in ascending order, which its
// own comment says must give bwa's sorted-array estimator bit for bit. The
// existing tests compare its parallel and serial forms with each other, and
// check one symmetric distribution whose percentiles do not move when a rank
// is off by one; neither would catch a wrong percentile rank, a wrong outlier
// or mapping multiplier, or an inclusive bound made exclusive.
//
// The oracle below is bwa's estimator as published (bwa 0.7.x bwamem_pair.c,
// mem_pestat): collect the insert sizes of unique same-contig pairs per
// orientation into an array, sort it, read the 25th/50th/75th percentiles at
// index (int)(p * n + .499), take the mean and standard deviation of the values
// inside [p25 - 2 IQR, p75 + 2 IQR], then set the proper-pair window to
// [p25 - 3 IQR, p75 + 3 IQR] widened to cover mean -/+ 4 std, low at least 1;
// an orientation with fewer than 10 pairs, or fewer than 5% of the most common
// orientation's, fails. avg is a sum of integers divided once, which is exact on
// every target, so it is compared bitwise; std's sum of squares may be contracted
// to fused multiply-adds in one build and not the other (clang does so wherever
// the target has FMA), so it is compared to 1e-12. Hand-worked distributions pin
// each rule with the arithmetic written next to the check.

#include "doctest/doctest.h"
#include "bwamem.h"
#include "stderr_capture.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <random>
#include <vector>

namespace {

const int64_t L_PAC = 1000000;

struct OptGuard {
    mem_opt_t *o;
    OptGuard() : o(mem_opt_init()) {}
    ~OptGuard() { free(o); }
};

// Per-read region lists for a cohort of reads (2 reads per pair); owns storage.
struct Cohort {
    std::vector<std::vector<mem_alnreg_t> > store;
    std::vector<mem_alnreg_v> v;

    explicit Cohort(int n_reads) : store((size_t)n_reads), v((size_t)n_reads) {}

    // Appends one region to read `i`: reference [rb, rb + len), query [qb, qb + len).
    void add(int i, int rid, int64_t rb, int score, int qb = 0, int len = 100) {
        mem_alnreg_t a;
        memset(&a, 0, sizeof(a));
        a.rid = rid;
        a.rb = rb;
        a.re = rb + len;
        a.qb = qb;
        a.qe = qb + len;
        a.score = score;
        store[(size_t)i].push_back(a);
    }
    mem_alnreg_v *regs() {
        for (size_t i = 0; i < v.size(); ++i) {
            v[i].n = v[i].m = store[i].size();
            v[i].a = store[i].empty() ? NULL : store[i].data();
            v[i].capped_w = 0;
        }
        return v.data();
    }
    int n() const { return (int)v.size(); }
};

// Doubled-strand coordinate of a reverse-strand region whose forward-strand span
// ends (exclusive) at fwd_end: the reverse strand is the mirror about 2*l_pac.
int64_t rev_rb(int64_t fwd_end) { return 2 * L_PAC - fwd_end; }

// Pair p of the cohort as an FR pair: read 1 forward at [pos, pos + 100), read 2
// reverse ending at pos + dist (inclusive), so the insert mem_pestat sees is
// `dist` (the fragment length minus one: TLEN - 1).
void add_fr(Cohort &c, int p, int64_t pos, int64_t dist, int rid = 0) {
    c.add(2 * p, rid, pos, 100);
    c.add(2 * p + 1, rid, rev_rb(pos + dist + 1), 100);
}

// ---- the oracle: bwa's sorted-array estimator -------------------------------

struct OraclePes {
    int low, high, failed;
    double avg, std;
};

// The query-overlap rule of bwa's cal_sub: the score of the first later region
// that overlaps region 0 on the read by at least mask_level of the shorter
// span, else min_seed_len * a.
int oracle_cal_sub(const mem_opt_t *opt, const mem_alnreg_v &r) {
    for (size_t j = 1; j < r.n; ++j) {
        const int b = std::max(r.a[j].qb, r.a[0].qb), e = std::min(r.a[j].qe, r.a[0].qe);
        const int min_l = std::min(r.a[j].qe - r.a[j].qb, r.a[0].qe - r.a[0].qb);
        if (e > b && e - b >= min_l * opt->mask_level) return r.a[j].score;
    }
    return opt->min_seed_len * opt->a;
}

// Orientation by table: project read 2 onto read 1's strand (mirror it when the
// strands differ); FF/RR when the strands agree, FR/RF when they differ, the
// first letter of each pair naming the strand of the upstream end.
int oracle_dir(int64_t b1, int64_t b2, int64_t *dist) {
    const bool rev1 = b1 >= L_PAC, rev2 = b2 >= L_PAC;
    const int64_t p2 = rev1 == rev2 ? b2 : 2 * L_PAC - 1 - b2;
    const bool downstream = p2 > b1;
    *dist = downstream ? p2 - b1 : b1 - p2;
    if (rev1 == rev2) return downstream ? 0 : 3;
    return downstream ? 1 : 2;
}

void oracle_pestat(const mem_opt_t *opt, Cohort &c, OraclePes out[4]) {
    mem_alnreg_v *regs = c.regs();
    std::vector<int64_t> is[4];
    for (int i = 0; i < c.n() / 2; ++i) {
        const mem_alnreg_v &r0 = regs[2 * i], &r1 = regs[2 * i + 1];
        if (r0.n == 0 || r1.n == 0) continue;
        if (oracle_cal_sub(opt, r0) > 0.8 * r0.a[0].score) continue;
        if (oracle_cal_sub(opt, r1) > 0.8 * r1.a[0].score) continue;
        if (r0.a[0].rid != r1.a[0].rid) continue;
        int64_t d;
        const int dir = oracle_dir(r0.a[0].rb, r1.a[0].rb, &d);
        if (d > 0 && d <= opt->max_ins) is[dir].push_back(d);
    }
    size_t most = 0;
    for (int d = 0; d < 4; ++d) most = std::max(most, is[d].size());
    for (int d = 0; d < 4; ++d) {
        OraclePes &r = out[d];
        memset(&r, 0, sizeof(r));
        std::vector<int64_t> &q = is[d];
        const size_t n = q.size();
        if (n < 10) {
            r.failed = 1;
            continue;
        }
        std::sort(q.begin(), q.end());
        const int64_t p25 = q[(size_t)(.25 * n + .499)], p75 = q[(size_t)(.75 * n + .499)];
        r.low = (int)(p25 - 2.0 * (p75 - p25) + .499);
        if (r.low < 1) r.low = 1;
        r.high = (int)(p75 + 2.0 * (p75 - p25) + .499);
        int x = 0;
        r.avg = 0;
        for (size_t i = 0; i < n; ++i)
            if (q[i] >= r.low && q[i] <= r.high) r.avg += q[i], ++x;
        r.avg /= x;
        r.std = 0;
        for (size_t i = 0; i < n; ++i)
            if (q[i] >= r.low && q[i] <= r.high) r.std += (q[i] - r.avg) * (q[i] - r.avg);
        r.std = sqrt(r.std / x);
        r.low = (int)(p25 - 3.0 * (p75 - p25) + .499);
        r.high = (int)(p75 + 3.0 * (p75 - p25) + .499);
        if (r.low > r.avg - 4.0 * r.std) r.low = (int)(r.avg - 4.0 * r.std + .499);
        if (r.high < r.avg + 4.0 * r.std) r.high = (int)(r.avg + 4.0 * r.std + .499);
        if (r.low < 1) r.low = 1;
        if (n < most * 0.05) r.failed = 1;
    }
}

void run_pestat(const mem_opt_t *opt, Cohort &c, mem_pestat_t pes[4]) {
    memset(pes, 0x5a, 4 * sizeof(mem_pestat_t));  // nothing may survive from before
    bwa_tests::QuietStderr quiet;                  // it reports every orientation
    mem_pestat(opt, L_PAC, c.n(), c.regs(), pes);
}

// Returns whether mem_pestat fitted FR.
bool check_against_oracle(const mem_opt_t *opt, Cohort &c) {
    OraclePes want[4];
    oracle_pestat(opt, c, want);
    mem_pestat_t got[4];
    run_pestat(opt, c, got);
    for (int d = 0; d < 4; ++d) {
        CAPTURE(d);
        REQUIRE(got[d].failed == want[d].failed);
        if (want[d].failed) continue;
        CHECK(got[d].low == want[d].low);
        CHECK(got[d].high == want[d].high);
        CHECK(memcmp(&got[d].avg, &want[d].avg, sizeof(double)) == 0);
        CHECK(got[d].std == doctest::Approx(want[d].std).epsilon(1e-12));
    }
    return !got[1].failed;
}

// A cohort whose inserts come from a skewed, gapped distribution, so neighbouring
// ranks hold different values and an off-by-one percentile moves the window.
// Some pairs fall into every branch of the candidate filter, and some have a
// distance of 0, which no orientation counts.
Cohort random_cohort(std::mt19937_64 &rng, int n_pairs, int64_t center, int64_t spread) {
    Cohort c(2 * n_pairs);
    for (int p = 0; p < n_pairs; ++p) {
        const int kind = (int)(rng() % 40);
        const int64_t pos = 20000 + (int64_t)(rng() % 400000);
        // Right-skewed: the sum of two draws, one of them squared down.
        const int64_t u = (int64_t)(rng() % (uint64_t)spread);
        const int64_t w = (int64_t)(rng() % (uint64_t)spread);
        int64_t dist = center + u + (w * w) / (spread * 2);
        if (kind == 0) {                          // a far outlier inside max_ins
            dist = center + 20 * spread + (int64_t)(rng() % 1000);
        } else if (kind == 1) {                   // a short outlier
            dist = 1 + (int64_t)(rng() % 40);
        }
        if (kind == 2) {                          // read 2 unmapped
            c.add(2 * p, 0, pos, 100);
        } else if (kind == 7) {                   // read 1 unmapped
            c.add(2 * p + 1, 0, rev_rb(pos + dist + 1), 100);
        } else if (kind == 8) {                   // read 2 ambiguous: sub 81 > 0.8 * 100
            c.add(2 * p, 0, pos, 100);
            c.add(2 * p + 1, 0, rev_rb(pos + dist + 1), 100);
            c.add(2 * p + 1, 0, rev_rb(pos + dist + 5001), 81);
        } else if (kind == 9) {                   // both forward at one position: distance 0
            c.add(2 * p, 0, pos, 100);
            c.add(2 * p + 1, 0, pos, 100);
        } else if (kind == 3) {                   // different contigs
            c.add(2 * p, 0, pos, 100);
            c.add(2 * p + 1, 1, rev_rb(pos + dist + 1), 100);
        } else if (kind == 4) {                   // read 1 ambiguous: sub 81 > 0.8 * 100
            c.add(2 * p, 0, pos, 100);
            c.add(2 * p, 0, pos + 5000, 81);
            c.add(2 * p + 1, 0, rev_rb(pos + dist + 1), 100);
        } else if (kind == 5) {                   // RF: read 2 reverse, last base pos - dist
            c.add(2 * p, 0, pos, 100);
            c.add(2 * p + 1, 0, rev_rb(pos - dist + 1), 100);
        } else if (kind == 6) {                   // FF, read 2 downstream
            c.add(2 * p, 0, pos, 100);
            c.add(2 * p + 1, 0, pos + dist, 100);
        } else {
            add_fr(c, p, pos, dist);
        }
    }
    return c;
}

}  // namespace

TEST_CASE("mem_infer_dir: orientation and distance of hand-placed pairs"
          * doctest::test_suite("unit/pair")) {
    // l_pac = 1000000. A forward region starts at rb; a reverse region with
    // forward span [s, e) has rb = 2000000 - e, and its projection onto the
    // forward strand, 1999999 - rb = e - 1, is its last forward base.
    struct Case {
        const char *what;
        int64_t b1, b2;
        int dir;
        int64_t dist;
    } cases[] = {
        // read 1 F at [1000, 1100), read 2 R at [1200, 1300): FR, 1299 - 1000.
        {"FR, read 1 upstream", 1000, rev_rb(1300), 1, 299},
        // read 1 R at [1200, 1300), read 2 F at [1000, 1100): read 2 projected
        // onto read 1's (reverse) strand is 1999999 - 1000 = 1998999, and read 1
        // is at 2000000 - 1300 = 1998700: downstream by 299, still FR.
        {"FR, read 2 upstream", rev_rb(1300), 1000, 1, 299},
        // read 1 F at [1200, 1300), read 2 R at [1000, 1100): read 2's last base
        // 1099 lies upstream of read 1's first: RF, 1200 - 1099 = 101.
        {"RF", 1200, rev_rb(1100), 2, 101},
        // both forward, read 2 downstream: FF, 1500 - 1000.
        {"FF", 1000, 1500, 0, 500},
        // both forward, read 2 upstream: RR in mem_infer_dir's numbering, 500.
        {"FF with read 2 upstream", 1500, 1000, 3, 500},
        // both reverse, read 2 downstream on the forward strand, so upstream on
        // the reverse strand: 3, distance 1998500 - 1998300 = 200.
        {"RR", rev_rb(1500), rev_rb(1700), 3, 200},
        // the same start on the same strand: distance 0, orientation 3.
        {"same position", 4000, 4000, 3, 0},
    };
    for (const Case &k : cases) {
        CAPTURE(k.what);
        int64_t dist = -1;
        CHECK(mem_infer_dir(L_PAC, k.b1, k.b2, &dist) == k.dir);
        CHECK(dist == k.dist);
    }
}

TEST_CASE("mem_pestat: random skewed cohorts equal bwa's sorted-array estimator"
          * doctest::test_suite("unit/pair")) {
    OptGuard g;
    std::mt19937_64 rng(20261024);
    // Sizes straddle MIN_DIR_CNT for the minority orientations; spreads from a
    // tight library (IQR of a few bp) to a wide one.
    const struct {
        int n;
        int64_t spread;
    } cohorts[] = {{12, 3}, {77, 25}, {200, 160}, {641, 25}};
    int fitted = 0;
    for (const auto &k : cohorts) {
        CAPTURE(k.n);
        CAPTURE(k.spread);
        Cohort c = random_cohort(rng, k.n, 150 + (int64_t)(rng() % 400), k.spread);
        fitted += check_against_oracle(g.o, c);
    }
    // Guard against a vacuous pass: most cohorts must fit an FR distribution.
    CHECK(fitted >= 3);
}

TEST_CASE("mem_pestat: percentile ranks on twelve distinct inserts"
          * doctest::test_suite("unit/pair")) {
    // FR inserts 101..112, one each. Ranks: p25 = q[(int)(3 + .499)] = q[3] = 104,
    // p75 = q[(int)(9 + .499)] = q[9] = 110, IQR 6. Outlier bounds
    // [104 - 12, 110 + 12] = [92, 122] keep all twelve: avg = 1278 / 12 = 106.5,
    // std = sqrt(2 * (0.25 + 2.25 + 6.25 + 12.25 + 20.25 + 30.25) / 12)
    //     = sqrt(143 / 12) = 3.4521. Proper window [104 - 18, 110 + 18] = [86, 128];
    // avg -/+ 4 std = [92.69, 120.31] lies inside it, so it is not widened.
    // A rank off by one in either direction moves p25 or p75 and the window.
    OptGuard g;
    Cohort c(24);
    for (int p = 0; p < 12; ++p) add_fr(c, p, 5000 + 1000 * p, 101 + p);
    mem_pestat_t pes[4];
    run_pestat(g.o, c, pes);
    CHECK(pes[0].failed == 1);
    CHECK(pes[2].failed == 1);
    CHECK(pes[3].failed == 1);
    REQUIRE(pes[1].failed == 0);
    CHECK(pes[1].low == 86);
    CHECK(pes[1].high == 128);
    CHECK(pes[1].avg == 106.5);
    CHECK(pes[1].std == doctest::Approx(std::sqrt(143.0 / 12.0)).epsilon(1e-12));
}

TEST_CASE("mem_pestat: fractional percentile ranks round at .5"
          * doctest::test_suite("unit/pair")) {
    // The twelve-insert case's ranks .25n and .75n are whole, so they cannot
    // tell (int)(p * n + .499) from a truncation or a ceiling. FR inserts
    // 101..115 (n = 15): p25 = q[(int)(3.75 + .499)] = q[4] = 105 (truncated:
    // q[3]), p75 = q[(int)(11.25 + .499)] = q[11] = 112, IQR 7. Outlier bounds
    // [91, 126] keep all fifteen: avg 108, std = sqrt(2 * 140 / 15) = 4.3205.
    // Proper window [105 - 21, 112 + 21] = [84, 133]; avg -/+ 4 std = [90.72,
    // 125.28] lies inside it. n = 13 and 14 put the other ranks off the integer
    // (.75 * 13 = 9.75 rounds up; .25 * 14 = 3.5 does not) and are compared
    // with the oracle.
    OptGuard g;
    for (int n = 13; n <= 15; ++n) {
        CAPTURE(n);
        Cohort c(2 * n);
        for (int p = 0; p < n; ++p) add_fr(c, p, 5000 + 1000 * p, 101 + p);
        REQUIRE(check_against_oracle(g.o, c));
        if (n != 15) continue;
        mem_pestat_t pes[4];
        run_pestat(g.o, c, pes);
        CHECK(pes[1].low == 84);
        CHECK(pes[1].high == 133);
        CHECK(pes[1].avg == 108.0);
        CHECK(pes[1].std == doctest::Approx(std::sqrt(280.0 / 15.0)).epsilon(1e-12));
    }
}

TEST_CASE("mem_pestat: an insert past p75 + 2 IQR is left out of the mean"
          * doctest::test_suite("unit/pair")) {
    // Inserts 10 x 100, 10 x 110, 1 x 131 (n = 21): p25 = q[5] = 100, p75 =
    // q[16] = 110, IQR 10, outlier bounds [80, 130]. 131 is outside them (it
    // would be inside a 2.5 or 3 IQR bound): avg = 2100 / 20 = 105, std 5.
    // Proper window [70, 140] contains avg -/+ 4 std = [85, 125].
    OptGuard g;
    Cohort c(42);
    for (int p = 0; p < 21; ++p) add_fr(c, p, 4000 + 800 * p, p < 10 ? 100 : p < 20 ? 110 : 131);
    mem_pestat_t pes[4];
    run_pestat(g.o, c, pes);
    REQUIRE(pes[1].failed == 0);
    CHECK(pes[1].avg == 105.0);
    CHECK(pes[1].std == doctest::Approx(5.0).epsilon(1e-12));
    CHECK(pes[1].low == 70);
    CHECK(pes[1].high == 140);
}

TEST_CASE("mem_pestat: the window widens to mean -/+ 4 std"
          * doctest::test_suite("unit/pair")) {
    // Inserts 4 x 80, 6 x 100, 6 x 110, 4 x 130 (n = 20). p25 = q[5] = 100,
    // p75 = q[15] = 110, IQR 10; outlier bounds [80, 130] keep all twenty.
    // avg = (320 + 600 + 660 + 520) / 20 = 105; std = sqrt((4*625 + 6*25 + 6*25 +
    // 4*625) / 20) = sqrt(265) = 16.2788. The IQR window [70, 140] is narrower than
    // 105 -/+ 4 * 16.2788 = [39.885, 170.115], so both ends widen:
    // low = (int)(39.885 + .499) = 40, high = (int)(170.115 + .499) = 170.
    OptGuard g;
    Cohort c(40);
    const int64_t ins[4] = {80, 100, 110, 130};
    const int cnt[4] = {4, 6, 6, 4};
    int p = 0;
    for (int k = 0; k < 4; ++k)
        for (int j = 0; j < cnt[k]; ++j, ++p) add_fr(c, p, 3000 + 700 * p, ins[k]);
    mem_pestat_t pes[4];
    run_pestat(g.o, c, pes);
    REQUIRE(pes[1].failed == 0);
    CHECK(pes[1].avg == 105.0);
    CHECK(pes[1].std == doctest::Approx(std::sqrt(265.0)).epsilon(1e-12));
    CHECK(pes[1].low == 40);
    CHECK(pes[1].high == 170);
}

TEST_CASE("mem_pestat: a window reaching below 1 is clamped to 1"
          * doctest::test_suite("unit/pair")) {
    OptGuard g;
    // Inserts 1..12: p25 = q[3] = 4, p75 = q[9] = 10, IQR 6. The outlier low
    // (int)(4 - 12 + .499) = -7 clamps to 1, keeping all twelve: avg 6.5,
    // std = sqrt(143 / 12) = 3.452. Proper low (int)(4 - 18 + .499) = -13 is
    // below avg - 4 std = -7.31, so it is not raised, and clamps to 1; high =
    // (int)(10 + 18 + .499) = 28 > avg + 4 std = 20.31.
    Cohort c(24);
    for (int p = 0; p < 12; ++p) add_fr(c, p, 7000 + 600 * p, 1 + p);
    mem_pestat_t pes[4];
    run_pestat(g.o, c, pes);
    REQUIRE(pes[1].failed == 0);
    CHECK(pes[1].low == 1);
    CHECK(pes[1].high == 28);
    CHECK(pes[1].avg == 6.5);
}

TEST_CASE("mem_pestat: an orientation needs 10 pairs and 5% of the largest count"
          * doctest::test_suite("unit/pair")) {
    OptGuard g;
    mem_pestat_t pes[4];
    SUBCASE("9 FR pairs fail, 10 do not") {
        for (int n = 9; n <= 10; ++n) {
            CAPTURE(n);
            Cohort c(2 * n);
            for (int p = 0; p < n; ++p) add_fr(c, p, 1000 + 500 * p, 300);
            run_pestat(g.o, c, pes);
            CHECK(pes[1].failed == (n < 10 ? 1 : 0));
        }
    }
    SUBCASE("10 RF pairs survive beside 200 FR (10 >= 10) but not beside 201 (10 < 10.05)") {
        for (int n_fr = 200; n_fr <= 201; ++n_fr) {
            CAPTURE(n_fr);
            Cohort c(2 * (n_fr + 10));
            int p = 0;
            for (; p < n_fr; ++p) add_fr(c, p, 1000 + 1000 * p, 300);
            for (int j = 0; j < 10; ++j, ++p) {
                // RF: read 1 forward at pos, read 2 reverse with its last base at
                // pos - 201 (upstream), an RF insert of 201.
                const int64_t pos = 600000 + 1000 * j;
                c.add(2 * p, 0, pos, 100);
                c.add(2 * p + 1, 0, rev_rb(pos - 200), 100);
            }
            run_pestat(g.o, c, pes);
            CHECK(pes[1].failed == 0);
            CHECK(pes[2].failed == (n_fr == 201 ? 1 : 0));
        }
    }
}

TEST_CASE("mem_pestat: the candidate filter's edges"
          * doctest::test_suite("unit/pair")) {
    // Nine ordinary FR pairs at insert 300 plus one probe pair: the orientation
    // fits (10 pairs) exactly when the probe is counted.
    OptGuard g;
    struct Probe {
        const char *what;
        int sub_score;   // score of a second region overlapping read 1, or -1
        int64_t dist;
        bool counted;
    } probes[] = {
        {"read 1's second hit at 0.8 of its score (80 > 80 is false)", 80, 300, true},
        {"read 1's second hit above 0.8 of its score", 81, 300, false},
        {"insert equal to max_ins (10000)", -1, 10000, true},
        {"insert one past max_ins", -1, 10001, false},
    };
    for (const Probe &pr : probes) {
        CAPTURE(pr.what);
        Cohort c(20);
        for (int p = 0; p < 9; ++p) add_fr(c, p, 1000 + 500 * p, 300);
        const int64_t pos = 50000;
        c.add(18, 0, pos, 100);
        if (pr.sub_score >= 0) c.add(18, 0, pos + 30000, pr.sub_score);
        c.add(19, 0, rev_rb(pos + pr.dist + 1), 100);
        mem_pestat_t pes[4];
        run_pestat(g.o, c, pes);
        CHECK(pes[1].failed == (pr.counted ? 0 : 1));
    }
}
