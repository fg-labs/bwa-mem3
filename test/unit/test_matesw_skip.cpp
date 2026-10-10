// test/unit/test_matesw_skip.cpp — mem_matesw_skip, the mate-rescue orientation check.
//
// Before rescuing a mate around a region at rb, mem_matesw_batch_pre and mem_matesw_batch_post
// decide, per orientation r, whether to skip it: when pes[r] failed, or when some region of the
// mate already lies in orientation r at a distance inside [pes[r].low, pes[r].high]. Both used to
// run that as a loop over the mate's regions calling mem_infer_dir; mem_matesw_skip is the shared
// form, which off aarch64 is branch-free and stops once all four orientations are set (aarch64
// keeps the loop).
//
// The oracle below is the loop it replaced, verbatim, over mem_infer_dir. Every case compares
// the four flags and the return value against it: random region sets on both strands, insert
// windows that are empty, a single distance wide, or include 0, failed orientations, region
// counts on both sides of the 8-region block the early exit checks between, and regions at the
// strand boundary and at distance 0 (p2 == rb, which mem_infer_dir maps to orientation 3 or 2).
//
// Fixtures are built in-memory; no test data files are read.
#include "doctest/doctest.h"
#include "bwamem.h"
#include <cstdint>
#include <random>
#include <vector>

namespace {

// The orientation check as mem_matesw_batch_pre / _post wrote it before mem_matesw_skip.
bool oracle(int64_t l_pac, int64_t rb, const mem_alnreg_v *ma, const mem_pestat_t pes[4], int skip[4])
{
    for (int r = 0; r < 4; ++r) skip[r] = pes[r].failed ? 1 : 0;
    for (size_t i = 0; i < ma->n; ++i) {
        int64_t dist;
        const int r = mem_infer_dir(l_pac, rb, ma->a[i].rb, &dist);
        if (dist >= pes[r].low && dist <= pes[r].high) skip[r] = 1;
    }
    return skip[0] + skip[1] + skip[2] + skip[3] == 4;
}

void check_same(int64_t l_pac, int64_t rb, std::vector<mem_alnreg_t> &regs, const mem_pestat_t pes[4])
{
    mem_alnreg_v ma = {regs.size(), regs.size(), regs.empty() ? nullptr : regs.data(), 0};
    int want[4], got[4];
    const bool want_all = oracle(l_pac, rb, &ma, pes, want);
    const bool got_all = mem_matesw_skip(l_pac, rb, &ma, pes, got);
    CAPTURE(l_pac);
    CAPTURE(rb);
    CAPTURE(ma.n);
    CHECK(got_all == want_all);
    for (int r = 0; r < 4; ++r) {
        CAPTURE(r);
        CHECK(got[r] == want[r]);
    }
}

mem_alnreg_t region_at(int64_t rb)
{
    mem_alnreg_t a = {};
    a.rb = rb;
    return a;
}

} // namespace

TEST_CASE("mate rescue skip: matches the per-region mem_infer_dir loop on random mates"
          * doctest::test_suite("unit/pair")) {
    std::mt19937_64 rng(0x5eed5c4a);
    const int64_t L_PAC = 3100000000LL;   // hg38-sized, so 2 * l_pac exceeds 32 bits
    auto pick = [&](int64_t lo, int64_t hi) {
        return std::uniform_int_distribution<int64_t>(lo, hi)(rng);
    };
    for (int iter = 0; iter < 20000; ++iter) {
        mem_pestat_t pes[4] = {};
        for (int r = 0; r < 4; ++r) {
            pes[r].failed = pick(0, 3) == 0;
            pes[r].low = (int)pick(-50, 400);
            pes[r].high = pes[r].low + (int)pick(-20, 900);   // sometimes empty (high < low)
        }
        const int64_t rb = pick(0, 2 * L_PAC - 1);
        // the mate near the region on its strand, on the mirrored strand, and anywhere
        const int64_t mirror = 2 * L_PAC - 1 - rb;
        const size_t n = (size_t)pick(0, 40);
        std::vector<mem_alnreg_t> regs;
        for (size_t i = 0; i < n; ++i) {
            const int kind = (int)pick(0, 2);
            int64_t b2 = kind == 0 ? rb + pick(-1500, 1500) : kind == 1 ? mirror + pick(-1500, 1500)
                                                                         : pick(0, 2 * L_PAC - 1);
            if (b2 < 0) b2 = 0;
            if (b2 > 2 * L_PAC - 1) b2 = 2 * L_PAC - 1;
            regs.push_back(region_at(b2));
        }
        check_same(L_PAC, rb, regs, pes);
    }
}

TEST_CASE("mate rescue skip: window edges, distance 0, the strand boundary and the 8-region blocks"
          * doctest::test_suite("unit/pair")) {
    const int64_t L_PAC = 1000000;
    const int64_t RB = 400000;   // forward strand
    mem_pestat_t pes[4] = {};
    for (int r = 0; r < 4; ++r) pes[r].low = 100, pes[r].high = 500;

    SUBCASE("no regions: only failed orientations are set") {
        std::vector<mem_alnreg_t> regs;
        check_same(L_PAC, RB, regs, pes);
        pes[1].failed = pes[3].failed = 1;
        check_same(L_PAC, RB, regs, pes);
        pes[0].failed = pes[2].failed = 1;
        check_same(L_PAC, RB, regs, pes);   // all four failed: a consistent pair is reported
    }
    SUBCASE("each distance at and around both window edges, both strands, both sides") {
        for (const int64_t d : {0, 1, 99, 100, 101, 499, 500, 501}) {
            for (const int64_t b2 : {RB + d, RB - d, 2 * L_PAC - 1 - (RB + d), 2 * L_PAC - 1 - (RB - d)}) {
                std::vector<mem_alnreg_t> regs = {region_at(b2)};
                check_same(L_PAC, RB, regs, pes);
            }
        }
    }
    SUBCASE("a window that includes distance 0") {
        for (int r = 0; r < 4; ++r) pes[r].low = 0;
        for (const int64_t b2 : {RB, 2 * L_PAC - 1 - RB, RB + 1, RB - 1}) {
            std::vector<mem_alnreg_t> regs = {region_at(b2)};
            check_same(L_PAC, RB, regs, pes);
        }
    }
    SUBCASE("regions at the strand boundary, a on either strand") {
        for (const int64_t rb : {L_PAC - 1, L_PAC}) {
            std::vector<mem_alnreg_t> regs = {region_at(L_PAC - 1), region_at(L_PAC), region_at(0),
                                              region_at(2 * L_PAC - 1)};
            check_same(L_PAC, rb, regs, pes);
        }
        // A mate exactly at the first reverse-strand position (or the last forward one), in the
        // window of an a on either side of the boundary: which strand it is on decides both the
        // orientation and, through the mirror, the distance.
        // The same for an a exactly at either side of the boundary, with mates in its window.
        for (const int64_t rb : {L_PAC - 300, L_PAC - 1, L_PAC, L_PAC + 300}) {
            for (const int64_t b2 : {L_PAC - 301, L_PAC - 300, L_PAC - 1, L_PAC, L_PAC + 299, L_PAC + 300}) {
                std::vector<mem_alnreg_t> regs = {region_at(b2)};
                check_same(L_PAC, rb, regs, pes);
            }
        }
    }
    SUBCASE("the hit for each orientation at every position of 1 to 25 regions") {
        // In-window distances for each orientation of a forward rb: FF ahead (r = 0), the
        // mirrored strand ahead (1) and behind (2), FF behind (3); filler regions match nothing.
        const int64_t hit[4] = {RB + 300, 2 * L_PAC - 1 - (RB + 300), 2 * L_PAC - 1 - (RB - 300), RB - 300};
        for (int r = 0; r < 4; ++r) {
            int64_t d;
            REQUIRE(mem_infer_dir(L_PAC, RB, hit[r], &d) == r);
        }
        for (size_t n = 1; n <= 25; ++n) {
            for (size_t at = 0; at < n; ++at) {
                std::vector<mem_alnreg_t> regs(n, region_at(RB + 5000));
                regs[at] = region_at(hit[at % 4]);
                check_same(L_PAC, RB, regs, pes);
                // all four orientations present, the last one at `at`: the early exit point
                std::vector<mem_alnreg_t> all(n + 3, region_at(RB + 5000));
                for (int r = 0; r < 3; ++r) all[r] = region_at(hit[(r + 1) % 4]);
                all[at + 3] = region_at(hit[0]);
                check_same(L_PAC, RB, all, pes);
            }
        }
    }
}
