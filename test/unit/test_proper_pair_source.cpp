// test/unit/test_proper_pair_source.cpp — which alignment the proper-pair bit
// (FLAG 0x2) is derived from.
//
// bwa and bwa-mem2 both derive 0x2 from the top-scoring region a[0], even when
// the record they emit is a[which] (bwa bwamem_pair.c:411; bwa-mem2 inherited it
// verbatim from the 0.7.17 port). fg-labs/bwa-mem3#17 changed bwa-mem3 to derive
// it from a[which] instead, so the bit describes the record it rides on.
//
// That is defensible, but 0x2 is aligner-defined ("properly aligned according to
// the aligner"), so it is a choice rather than a correction -- and it made
// bwa-mem3's DEFAULT output differ from both upstreams at once
// (fg-labs/bwa-mem3#362). This project prefers upstream output on the default
// path, so the #17 behavior is opt-in and the default matches bwa.
//
// The divergence is reachable only when a read has ALT hits (`which != 0`
// requires `n_pri < a.n`), so on any run without a `.alt` sidecar the two
// derivations are the same value and this option is inert. That is exactly why
// the behavioral cases below build the alnreg arrays by hand: no fixture in the
// suite can reach the branch, so nothing else here would catch a consumer that
// ignored the option and hardcoded one index.
//
// Both mem_sam_pe no-pairing blocks route through mem_proper_pair_extra_flag, so
// these cases cover the scalar and batched paths at their single shared seam;
// test/regression/compat_byte_identical.sh pins that neither block grew a
// private copy of the decision, and proper_pair_alt.sh drives the option through
// a real ALT-aware alignment end to end.

#include "doctest/doctest.h"
#include "bwamem.h"

#include <cstdlib>
#include <cstring>

namespace {
// mem_opt_init() allocates; free through a guard so a failing CHECK cannot leak.
struct OptGuard {
    mem_opt_t *o;
    OptGuard() : o(mem_opt_init()) {}
    ~OptGuard() { free(o); }
};

// Size of one strand of the packed concat-with-reverse-complement index. Every
// rb below is < l_pac, i.e. forward strand, so mem_infer_dir's r1 == r2 == 0.
const int64_t L_PAC = 1000000;

// The expected FLAG 0x2 for a forward/forward pair, computed WITHOUT calling
// mem_infer_dir or any other production code: the distance is the plain
// coordinate difference, and with both ends on the same strand the orientation
// is pes index 0. This is the independent oracle -- if the production selection
// reads the wrong region, the coordinate it picks is wrong, and this disagrees.
//
// Valid only for rb1 > rb0. mem_infer_dir returns
// `(r1 == r2 ? 0 : 1) ^ (p2 > b1 ? 0 : 3)`, so a same-strand pair maps to index
// 0 only when the second coordinate is the greater one; reversed, it maps to
// index 3. Every fixture passed to this oracle below satisfies rb1 > rb0. (The
// which == {1,0} subcase does reach index 3, which is why it asserts the
// literal 0 rather than calling this.)
int expected_flag_ff(int64_t rb0, int64_t rb1, const mem_pestat_t &ff) {
    const int64_t dist = rb1 > rb0 ? rb1 - rb0 : rb0 - rb1;
    const bool proper = !ff.failed && dist >= ff.low && dist <= ff.high;
    return proper ? 2 : 0;
}

// Two regions per mate: index 0 is the top-scoring PRIMARY region, index 1 is
// the first ALT region (what `which` points at when the primary scores below
// -T). Only `rb` is read by the code under test.
struct PairFixture {
    mem_alnreg_t regs[2][2];
    mem_alnreg_v a[2];
    mem_pestat_t pes[4];

    PairFixture(int64_t pri0, int64_t alt0, int64_t pri1, int64_t alt1) {
        std::memset(regs, 0, sizeof(regs));
        regs[0][0].rb = pri0; regs[0][1].rb = alt0;
        regs[1][0].rb = pri1; regs[1][1].rb = alt1;
        for (int i = 0; i < 2; ++i) { a[i].n = 2; a[i].m = 2; a[i].a = regs[i]; }
        std::memset(pes, 0, sizeof(pes));
        // Only pes index 0 carries a usable window; leave the other three
        // failed so a mis-inferred direction shows up as "not proper" rather
        // than silently passing against some other window. Same-strand pairs
        // land on index 0 when the second coordinate is the greater one and on
        // index 3 otherwise (see expected_flag_ff), so a subcase whose opt-in
        // side reverses the operands correctly reports 0 through pes[3].failed.
        pes[0].low = 100; pes[0].high = 500; pes[0].failed = 0;
        for (int d = 1; d < 4; ++d) pes[d].failed = 1;
    }
};
} // namespace

TEST_CASE("proper-pair bit derives from a[0] by default, as both upstreams do"
          * doctest::test_suite("unit/pair")) {
    OptGuard g;
    REQUIRE(g.o != nullptr);
    // 0 = derive from a[0] (bwa / bwa-mem2). This is the load-bearing assertion:
    // flipping this default silently re-breaks byte-identity with BOTH upstreams
    // on every ALT-aware run, and only on ALT-aware runs, so no non-ALT test
    // would catch it.
    CHECK(g.o->proper_pair_from_emitted == 0);
}

TEST_CASE("the #17 behavior is reachable, i.e. the field is a real switch"
          * doctest::test_suite("unit/pair")) {
    // Guards against "fixing" the divergence by deleting the option outright:
    // #17's reasoning still stands for anyone who wants the bit to describe the
    // emitted record, so it must remain selectable.
    OptGuard g;
    REQUIRE(g.o != nullptr);
    g.o->proper_pair_from_emitted = 1;
    CHECK(g.o->proper_pair_from_emitted == 1);
}

TEST_CASE("FLAG 0x2 tracks the selected region, not merely the stored option"
          * doctest::test_suite("unit/pair")) {
    OptGuard g;
    REQUIRE(g.o != nullptr);

    SUBCASE("a[0] in the window, a[which] outside: only the default pairs") {
        // Primaries 300 apart -> inside [100,500]. ALTs 20,000 apart -> outside.
        PairFixture f(/*pri0=*/10000, /*alt0=*/70000,
                      /*pri1=*/10300, /*alt1=*/90000);
        const int which[2] = { 1, 1 };

        g.o->proper_pair_from_emitted = 0;
        CHECK(mem_proper_pair_extra_flag(g.o, L_PAC, f.a, which, f.pes)
              == expected_flag_ff(f.regs[0][0].rb, f.regs[1][0].rb, f.pes[0]));
        g.o->proper_pair_from_emitted = 1;
        CHECK(mem_proper_pair_extra_flag(g.o, L_PAC, f.a, which, f.pes)
              == expected_flag_ff(f.regs[0][1].rb, f.regs[1][1].rb, f.pes[0]));
        // ...and the two oracles must actually disagree, or the case proves
        // nothing about which index was read.
        CHECK(expected_flag_ff(f.regs[0][0].rb, f.regs[1][0].rb, f.pes[0]) == 2);
        CHECK(expected_flag_ff(f.regs[0][1].rb, f.regs[1][1].rb, f.pes[0]) == 0);
    }

    SUBCASE("a[which] in the window, a[0] outside: only the opt-in pairs") {
        // The mirror image, so a consumer that hardcoded either index fails one
        // subcase or the other rather than passing both by luck.
        PairFixture f(/*pri0=*/10000, /*alt0=*/70000,
                      /*pri1=*/90000, /*alt1=*/70250);
        const int which[2] = { 1, 1 };

        g.o->proper_pair_from_emitted = 0;
        CHECK(mem_proper_pair_extra_flag(g.o, L_PAC, f.a, which, f.pes) == 0);
        g.o->proper_pair_from_emitted = 1;
        CHECK(mem_proper_pair_extra_flag(g.o, L_PAC, f.a, which, f.pes) == 2);
    }

    SUBCASE("each mate uses its OWN which[], not the other's") {
        // `which` is computed per mate, so one mate can emit its ALT region
        // while the other emits its primary. Every case above passes
        // which == {1,1}, where reading which[0] for both mates is
        // indistinguishable from correct -- so a transposed index would slip
        // through. These two cases are asymmetric in opposite directions, so
        // either transposition fails one of them.
        SUBCASE("which == {0,1}") {
            // Default: a[0][0] vs a[1][0] -> 300, inside. Opt-in: a[0][0] vs
            // a[1][1] -> 80,000, outside. Using which[0] for BOTH would read
            // a[1][0] on the opt-in side and wrongly report proper.
            PairFixture f(/*pri0=*/10000, /*alt0=*/70000,
                          /*pri1=*/10300, /*alt1=*/90000);
            const int which[2] = { 0, 1 };
            g.o->proper_pair_from_emitted = 0;
            CHECK(mem_proper_pair_extra_flag(g.o, L_PAC, f.a, which, f.pes) == 2);
            g.o->proper_pair_from_emitted = 1;
            CHECK(mem_proper_pair_extra_flag(g.o, L_PAC, f.a, which, f.pes) == 0);
        }
        SUBCASE("which == {1,0}") {
            // The mirror: using which[1] for BOTH would read a[0][0] on the
            // opt-in side and, again, wrongly report proper.
            PairFixture f(/*pri0=*/10000, /*alt0=*/70000,
                          /*pri1=*/10300, /*alt1=*/90000);
            const int which[2] = { 1, 0 };
            g.o->proper_pair_from_emitted = 0;
            CHECK(mem_proper_pair_extra_flag(g.o, L_PAC, f.a, which, f.pes) == 2);
            g.o->proper_pair_from_emitted = 1;
            CHECK(mem_proper_pair_extra_flag(g.o, L_PAC, f.a, which, f.pes) == 0);
        }
    }

    SUBCASE("which == 0 (no ALT hits): the option is inert") {
        // The structural claim the whole design rests on -- without a `.alt`
        // sidecar `which` is always 0, so both settings must agree on every
        // input. Checked on a properly paired and an improperly paired case.
        for (int64_t rb1 : { (int64_t)10300, (int64_t)90000 }) {
            PairFixture f(10000, 70000, rb1, 70250);
            const int which[2] = { 0, 0 };
            const int want =
                expected_flag_ff(f.regs[0][0].rb, f.regs[1][0].rb, f.pes[0]);
            g.o->proper_pair_from_emitted = 0;
            CHECK(mem_proper_pair_extra_flag(g.o, L_PAC, f.a, which, f.pes) == want);
            g.o->proper_pair_from_emitted = 1;
            CHECK(mem_proper_pair_extra_flag(g.o, L_PAC, f.a, which, f.pes) == want);
        }
    }

    SUBCASE("a failed orientation is never called proper") {
        // pes[d].failed short-circuits ahead of the window test; without this a
        // window that happens to contain the distance would mask the check.
        PairFixture f(10000, 70000, 10300, 70250);
        const int which[2] = { 1, 1 };
        f.pes[0].failed = 1;
        for (int emitted = 0; emitted < 2; ++emitted) {
            g.o->proper_pair_from_emitted = emitted;
            CHECK(mem_proper_pair_extra_flag(g.o, L_PAC, f.a, which, f.pes) == 0);
        }
    }
}
