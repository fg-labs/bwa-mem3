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
