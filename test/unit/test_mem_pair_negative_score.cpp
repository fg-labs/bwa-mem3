// test/unit/test_mem_pair_negative_score.cpp — mem_pair on a region whose
// score is negative, and on a pair score that is not finite.
//
// mem_pair packs each candidate's score into the high 32 bits of a uint64_t
// sort key. It used to read the score back with `y >> 32`, which zero-extends:
// a negative score came back as ~4.29e9. Upstream never sees one (every region
// score is at least the seed's len*a), but under --meth the seed is rescored in
// original space against the per-strand matrix, so a seed whose converted-space
// match hides penalized cells scores below zero and the extended region keeps
// it. The pair score then overflowed int in a double-to-int conversion, which
// is undefined: x86 (cvttsd2si) gave INT_MIN, clamped to 0 (no pair), while
// aarch64 (fcvtzs) saturated to INT_MAX, so arm64 builds paired the read, its
// q_se came out of raw_mapq(score - csub) negative and was written as MAPQ
// (uint8)(-65) = 191, with no AS (score < 0).
//
// mem_pair now reads the score back signed and scores any candidate with a
// negative-score end as q = 0 (what x86 builds produced), on every target. The
// candidate is still pushed, so it counts toward *sub / *n_sub as it did on x86.
//
// On x86 the old code also passes the negative-score cases: its out-of-range
// conversion gives INT_MIN, which clamps to the asserted q = 0. So those cases
// pin the bug on aarch64 (CI has arm64 jobs) and guard the semantics on x86.
//
// Expected values are literals worked by hand at the default -A 1 (asserted),
// not a copy of mem_pair's formula: q = (int)(s0 + s1 + .721 * ln(2 erfc(|ns| /
// sqrt 2)) + .499). At ns = 0 the log term is .721 * ln 2 = .4998.

#include "doctest/doctest.h"
#include "bwamem.h"
#include "bntseq.h"

#include <cstdlib>
#include <cstring>
#include <cstdint>

namespace {

struct OptGuard {
    mem_opt_t *o;
    OptGuard() : o(mem_opt_init()) {}
    ~OptGuard() { free(o); }
};

const int64_t L_PAC = 1000000;

// A 150 bp region. A forward region's sort-key position is rb; a reverse one
// sits in the second half of the doubled reference with key position `pos`.
mem_alnreg_t region(int64_t pos, int score, bool rev = false) {
    mem_alnreg_t r;
    memset(&r, 0, sizeof r);
    r.rb = rev ? (L_PAC << 1) - 1 - pos : pos;
    r.re = r.rb + 150; r.qb = 0; r.qe = 150;
    r.rid = 0; r.score = r.truesc = score;
    r.secondary = r.secondary_all = -1;
    r.meth_hypothesis = r.meth_strand_hyp = -1;
    return r;
}

struct PairResult {
    int o, sub, n_sub, z[2];
};

struct PairFixture {
    OptGuard g;
    bntann1_t ann;
    bntseq_t bns;
    mem_pestat_t pes[4];
    bseq1_t s[2];
    PairFixture() {
        memset(&ann, 0, sizeof ann);
        ann.offset = 0; ann.len = L_PAC;
        memset(&bns, 0, sizeof bns);
        bns.l_pac = L_PAC; bns.n_seqs = 1; bns.anns = &ann;
        for (int d = 0; d < 4; ++d) {
            pes[d].low = 0; pes[d].high = 1000; pes[d].failed = 0;
            pes[d].avg = 300.; pes[d].std = 50.;
        }
        memset(s, 0, sizeof s);
        REQUIRE(g.o->a == 1);   // the literals below are worked at -A 1
    }
    PairResult run(mem_alnreg_t *r0, int n0, mem_alnreg_t *r1, int n1) {
        mem_alnreg_v a[2];
        a[0].n = a[0].m = n0; a[0].a = r0;
        a[1].n = a[1].m = n1; a[1].a = r1;
        int n_pri[2] = { n0, n1 };
        PairResult res = { -1, -1, -1, { -1, -1 } };
        res.o = mem_pair(g.o, &bns, 0, pes, s, a, 7, &res.sub, &res.n_sub, res.z, n_pri);
        return res;
    }
};

} // namespace

TEST_SUITE("unit/pair") {

// When no pair scores above 0, mem_pair still writes z[] for the q = 0
// candidate; the caller ignores z[] whenever the returned score is 0, so these
// cases assert the score, *sub and *n_sub only.

TEST_CASE("mem_pair: a negative-score end forms no proper pair") {
    PairFixture f;
    mem_alnreg_t r0[1] = { region(1000, -11) };
    mem_alnreg_t r1[1] = { region(1300, 23) };
    const PairResult p = f.run(r0, 1, r1, 1);
    CHECK(p.o == 0);
    CHECK(p.sub == 0);
    CHECK(p.n_sub == 0);
}

TEST_CASE("mem_pair: both ends negative form no proper pair") {
    PairFixture f;
    mem_alnreg_t r0[1] = { region(1000, -11) };
    mem_alnreg_t r1[1] = { region(1300, -4) };
    const PairResult p = f.run(r0, 1, r1, 1);
    CHECK(p.o == 0);
    CHECK(p.sub == 0);
    CHECK(p.n_sub == 0);
}

TEST_CASE("mem_pair: a negative-score candidate scores 0 beside a real pair") {
    PairFixture f;
    // End 0: a real hit at 1000 and a negative-score one at 1100; end 1: one
    // hit at 1300 that pairs with both by distance (300 = avg, and 200).
    mem_alnreg_t r0[2] = { region(1000, 40), region(1100, -20) };
    mem_alnreg_t r1[1] = { region(1300, 30) };
    const PairResult p = f.run(r0, 2, r1, 1);
    CHECK(p.o == 70);          // (int)(40 + 30 + .4998 + .499)
    CHECK(p.z[0] == 0);
    CHECK(p.z[1] == 0);
    CHECK(p.sub == 0);         // the negative-score candidate's q is 0, not INT_MAX
    // ...and it is still pushed, so it counts as a near pair (0 - 0 <= 7), as
    // on x86 before the fix; skipping it would make this 0 and raise q_pe.
    CHECK(p.n_sub == 1);
}

TEST_CASE("mem_pair: a zero-score end still takes the pair-score formula") {
    PairFixture f;
    mem_alnreg_t r0[1] = { region(1000, 0) };
    mem_alnreg_t r1[1] = { region(1300, 30) };
    const PairResult p = f.run(r0, 1, r1, 1);
    CHECK(p.o == 30);          // (int)(0 + 30 + .4998 + .499)
    CHECK(p.z[0] == 0);
    CHECK(p.z[1] == 0);
    CHECK(p.sub == 0);
    CHECK(p.n_sub == 0);
}

TEST_CASE("mem_pair: non-negative scores are unchanged") {
    PairFixture f;
    SUBCASE("insert size at the mean") {
        mem_alnreg_t r0[1] = { region(1000, 40) };
        mem_alnreg_t r1[1] = { region(1300, 30) };
        const PairResult p = f.run(r0, 1, r1, 1);
        CHECK(p.o == 70);      // (int)(70 + .4998 + .499)
        CHECK(p.z[0] == 0);
        CHECK(p.z[1] == 0);
        CHECK(p.sub == 0);
        CHECK(p.n_sub == 0);
    }
    SUBCASE("insert size two standard deviations out") {
        // dist 400, avg 300, std 50: ns = 2, erfc(2 / sqrt 2) = .0455003,
        // .721 * ln(.0910006) = -1.7282, so (int)(70 - 1.7282 + .499) = 68.
        mem_alnreg_t r0[1] = { region(1000, 40) };
        mem_alnreg_t r1[1] = { region(1400, 30) };
        const PairResult p = f.run(r0, 1, r1, 1);
        CHECK(p.o == 68);
        CHECK(p.sub == 0);
        CHECK(p.n_sub == 0);
    }
    SUBCASE("forward/reverse orientation") {
        // Only one orientation is valid. End 1 on the reverse strand with key
        // position 1300 pairs with the forward end 0 at 1000 at dist 300.
        for (int d = 0; d < 4; ++d) f.pes[d].failed = 1;
        f.pes[1].failed = 0;
        mem_alnreg_t r0[1] = { region(1000, 40) };
        mem_alnreg_t r1[1] = { region(1300, 30, true) };
        const PairResult p = f.run(r0, 1, r1, 1);
        CHECK(p.o == 70);
        CHECK(p.z[0] == 0);
        CHECK(p.z[1] == 0);
        CHECK(p.sub == 0);
        CHECK(p.n_sub == 0);
    }
}

TEST_CASE("mem_pair: an insert size where erfc underflows scores 0") {
    // std 1 puts dist 1000 at ns = 700: erfc(495) underflows to 0 and the log
    // term is -inf, whose int conversion is undefined. mem_pair scores it 0
    // (x86 and aarch64 both produced INT_MIN there, clamped to 0).
    PairFixture f;
    for (int d = 0; d < 4; ++d) f.pes[d].std = 1.;
    mem_alnreg_t r0[1] = { region(1000, 40) };
    mem_alnreg_t r1[1] = { region(2000, 30) };
    const PairResult p = f.run(r0, 1, r1, 1);
    CHECK(p.o == 0);
    CHECK(p.sub == 0);
    CHECK(p.n_sub == 0);
}

TEST_CASE("mem_pair: finite scores saturate only after the insert-size adjustment") {
    // Catch an unchecked out-of-range cast, or clamping the sum before applying
    // the insert-size penalty. At the mean the fractional term is .9988;
    // two standard deviations out it is -1.2292.
    const struct {
        int s0, s1;
        int64_t pos1;
        int want;
    } cases[] = {
        {1073741823, 1073741823, 1300, 2147483646},
        {1073741823, 1073741824, 1300, INT32_MAX},
        {1073741824, 1073741824, 1300, INT32_MAX},
        {1100000000, 1100000000, 1300, INT32_MAX},
        {INT32_MAX, INT32_MAX, 1300, INT32_MAX},
        {1073741824, 1073741824, 1400, 2147483646},
    };
    for (const auto &c : cases) {
        CAPTURE(c.s0); CAPTURE(c.s1); CAPTURE(c.pos1);
        PairFixture f;
        mem_alnreg_t r0[1] = { region(1000, c.s0) };
        mem_alnreg_t r1[1] = { region(c.pos1, c.s1) };
        const PairResult p = f.run(r0, 1, r1, 1);
        CHECK(p.o == c.want);
        CHECK(p.z[0] == 0);
        CHECK(p.z[1] == 0);
        CHECK(p.sub == 0);
        CHECK(p.n_sub == 0);
    }
}

TEST_CASE("mem_pair: a saturated pair retains its runner-up and selected indices") {
    PairFixture f;
    mem_alnreg_t r0[2] = { region(1000, 0), region(1000, 1100000000) };
    mem_alnreg_t r1[1] = { region(1300, 1100000000) };
    const PairResult p = f.run(r0, 2, r1, 1);
    CHECK(p.o == INT32_MAX);
    CHECK(p.z[0] == 1);
    CHECK(p.z[1] == 0);
    CHECK(p.sub == 1100000000);
    CHECK(p.n_sub == 1);
}

TEST_CASE("mem_pair: saturated runner-ups count only nearby scores") {
    PairFixture f;
    mem_alnreg_t r0[3] = { region(1000, 0), region(1000, 1100000000), region(1000, 1100000001) };
    mem_alnreg_t r1[1] = { region(1300, 1100000000) };
    const PairResult p = f.run(r0, 3, r1, 1);
    CHECK(p.o == INT32_MAX);
    CHECK(p.sub == INT32_MAX);
    CHECK(p.n_sub == 1); // the 1.1-billion candidate is far below the runner-up
    CHECK(p.z[0] == 2); // larger low-32-bit hash: 2824057045 vs 1411446537
    CHECK(p.z[1] == 0);
}

TEST_CASE("mem_pair_resolve_batch_post: high pair scores have defined MAPQ arithmetic") {
    // Catch overflow in the unpaired score, meth margin and MAPQ conversion
    // after mem_pair returns INT32_MAX. The single-end competitor is ten
    // points behind, so its own MAPQ rounds to zero in the legacy model.
    const struct {
        int s0, s1, penalty, flag, legacy_mapq, default_mapq;
    } cases[] = {
        {1100000000, 1100000000, 17, 1, 0, 60},
        {1073741823, 1073741824, 17, 3, 40, 60},
        {INT32_MAX, 0, INT32_MAX, 3, 40, 60},
    };
    for (const auto &c : cases) for (int meth = 0; meth <= 1; ++meth)
    for (int legacy = 0; legacy <= 1; ++legacy) {
        CAPTURE(c.s0); CAPTURE(c.s1); CAPTURE(c.penalty); CAPTURE(meth); CAPTURE(legacy);
        PairFixture f;
        f.g.o->flag |= MEM_F_NO_RESCUE;
        f.g.o->meth_mode = meth;
        f.g.o->pen_unpaired = c.penalty;
        if (legacy) f.g.o->mapQ_coef_len = 0;
        mem_alnreg_t regs[2] = { region(1000, c.s0), region(1300, c.s1) };
        mem_alnreg_v a[2] = {};
        for (int i = 0; i < 2; ++i) {
            regs[i].sub = regs[i].csub = regs[i].score > 10 ? regs[i].score - 10 : 0;
            regs[i].seedcov = 150;
            a[i].n = a[i].m = 1;
            a[i].a = &regs[i];
        }
        mem_cache *mmc = (mem_cache *)calloc(1, sizeof(mem_cache));
        REQUIRE(mmc != NULL);
        int32_t gcnt = 0;
        int n_pri[2], z[2], q_se[2], extra_flag, paired;
        mem_pair_resolve_batch_post(f.g.o, &f.bns, NULL, f.pes, 7, f.s, a,
                                   NULL, mmc, gcnt, 0, n_pri, z, q_se, &extra_flag, &paired);
        free(mmc);
        CHECK(paired == 1);
        CHECK(extra_flag == c.flag);
        CHECK(z[0] == 0);
        CHECK(z[1] == 0);
        CHECK(n_pri[0] == 1);
        CHECK(n_pri[1] == 1);
        const int want = legacy ? c.legacy_mapq : c.default_mapq;
        CHECK(q_se[0] == want);
        CHECK(q_se[1] == (c.s1 == 0 ? 0 : want));
    }
}

TEST_CASE("mem_pair_resolve_batch_post: meth second-best sum stays wide") {
    PairFixture f;
    f.g.o->flag |= MEM_F_NO_RESCUE;
    f.g.o->meth_mode = 1;
    f.g.o->pen_unpaired = 100000000;
    f.g.o->mapQ_coef_len = 0;
    mem_alnreg_t r0[2] = { region(1000, 1100000000), region(10000, 1099999990) };
    mem_alnreg_t r1[2] = { region(1300, 1100000000), region(20000, 1099999990) };
    for (int i = 0; i < 2; ++i) {
        r0[i].seedcov = r1[i].seedcov = 150;
        r0[i].csub = r1[i].csub = 1099999990;
    }
    const PairResult p = f.run(r0, 2, r1, 2);
    CHECK(p.o == INT32_MAX);
    CHECK(p.sub == 0);
    CHECK(p.n_sub == 0);
    CHECK(p.z[0] == 0);
    CHECK(p.z[1] == 0);
    mem_alnreg_v a[2] = {};
    a[0].a = r0; a[1].a = r1;
    for (int i = 0; i < 2; ++i) a[i].n = a[i].m = 2;
    mem_cache *mmc = (mem_cache *)calloc(1, sizeof(mem_cache));
    REQUIRE(mmc != NULL);
    int32_t gcnt = 0;
    int n_pri[2], z[2], q_se[2], flag, paired;
    mem_pair_resolve_batch_post(f.g.o, &f.bns, NULL, f.pes, 7, f.s, a,
                               NULL, mmc, gcnt, 0, n_pri, z, q_se, &flag, &paired);
    free(mmc);
    CHECK(paired == 1);
    CHECK(flag == 3);
    CHECK(z[0] == 0);
    CHECK(z[1] == 0);
    CHECK(n_pri[0] == 2);
    CHECK(n_pri[1] == 2);
    CHECK(a[0].n == 2);
    CHECK(a[1].n == 2);
    CHECK(r0[1].score == 1099999990);
    CHECK(r1[1].score == 1099999990);
    // INT32_MAX + 4 - 2*1099999990 < 0; a narrow sum instead raises MAPQ to 40.
    CHECK(q_se[0] == 0);
    CHECK(q_se[1] == 0);
}

TEST_CASE("mem_approx_mapq_se: consistent long-span match product stays wide") {
    OptGuard g;
    g.o->a = 127;
    mem_alnreg_t r = region(1000, 2000000000);
    r.qe = 20000000;
    r.re = r.rb + 20000000;
    r.seedcov = 20000000;
    r.sub = r.csub = 1999000000;
    r.sub_n = 1000;
    // l*A = 2540000000; identity .7938931298 gives raw MAPQ 600.
    // The 30-point multiplicity penalty still leaves MAPQ above 60.
    CHECK(mem_approx_mapq_se(g.o, &r) == 60);
    g.o->mapQ_coef_len = 0;
    CHECK(mem_approx_mapq_se(g.o, &r) == 0);
}

TEST_CASE("mem_approx_mapq_se: suboptimal count correction precedes clipping") {
    OptGuard g;
    mem_alnreg_t r = region(1000, 100);
    r.qe = 100; r.re = r.rb + 100; r.seedcov = 100;
    r.sub = r.csub = 70;
    r.sub_n = 100;
    // Default raw MAPQ 77 minus rounded 4.343*ln(101) = 20.
    CHECK(mem_approx_mapq_se(g.o, &r) == 57);
    g.o->mapQ_coef_len = 0;
    // Legacy raw MAPQ rounds to 41; same correction leaves 21.
    CHECK(mem_approx_mapq_se(g.o, &r) == 21);
    r.sub_n = INT32_MAX;
    CHECK(mem_approx_mapq_se(g.o, &r) == 0);
    g.o->mapQ_coef_len = 50;
    CHECK(mem_approx_mapq_se(g.o, &r) == 0);
}

TEST_CASE("mem_approx_mapq_se: legacy MAPQ detects a wrapped identity numerator") {
    OptGuard g;
    g.o->a = 127;
    g.o->mapQ_coef_len = 0;
    mem_alnreg_t r = region(1000, 100000000);
    r.qe = 20000000;
    r.re = r.rb + 20000000;
    r.seedcov = 20000000;
    r.sub = r.csub = 50000000;
    // l*A - score = 2440000000 exceeds INT32_MAX. Wide identity is about
    // 0.0687, giving MAPQ 1; a wrapped numerator skips scaling and clips to 60.
    CHECK(mem_approx_mapq_se(g.o, &r) == 1);
}

TEST_CASE("mem_pair_resolve_batch_post: tandem caps cannot produce negative MAPQ") {
    for (int meth = 0; meth <= 1; ++meth) for (int capped = 0; capped < 2; ++capped)
    for (int legacy = 0; legacy <= 1; ++legacy) {
        CAPTURE(meth); CAPTURE(capped); CAPTURE(legacy);
        PairFixture f;
        f.g.o->flag |= MEM_F_NO_RESCUE;
        f.g.o->meth_mode = meth;
        if (legacy) f.g.o->mapQ_coef_len = 0;
        mem_alnreg_t regs[2] = { region(1000, 100), region(1300, 100) };
        mem_alnreg_v a[2] = {};
        for (int i = 0; i < 2; ++i) {
            regs[i].seedcov = 150;
            regs[i].csub = i == capped ? 110 : 90;
            a[i].n = a[i].m = 1;
            a[i].a = &regs[i];
        }
        mem_cache *mmc = (mem_cache *)calloc(1, sizeof(mem_cache));
        REQUIRE(mmc != NULL);
        int32_t gcnt = 0;
        int n_pri[2], z[2], q_se[2], flag, paired;
        mem_pair_resolve_batch_post(f.g.o, &f.bns, NULL, f.pes, 7, f.s, a,
                                   NULL, mmc, gcnt, 0, n_pri, z, q_se, &flag, &paired);
        free(mmc);
        CHECK(paired == 1);
        CHECK(flag == 3);
        CHECK(q_se[capped] == 0);
        CHECK(q_se[!capped] == (legacy ? 53 : 56)); // single-end MAPQ plus the boost of 40
    }
}

} // TEST_SUITE("unit/pair")
