// test/unit/test_matesw_skip.cpp — mem_matesw_skip, the mate-rescue orientation check.
//
// Before rescuing a mate around a region at rb, mem_matesw_batch_pre and mem_matesw_batch_post
// decide, per orientation r, whether to skip it: when pes[r] failed, or when some region of the
// mate already lies in orientation r at a distance inside [pes[r].low, pes[r].high]. Both used to
// run that as a loop over the mate's regions calling mem_infer_dir; on every architecture the test
// is now four intervals of mate rb values, one per orientation, scanned with range compares that
// stop once all four orientations are set (mem_matesw_skip is that form as one call).
//
// The oracle below is the loop it replaced, verbatim, over mem_infer_dir. Every case compares
// the four flags and the return value against it: random region sets on both strands, insert
// windows that are empty, a single distance wide, or include 0, failed orientations, region
// counts on both sides of the 8-region block the early exit checks between, and regions at the
// strand boundary and at distance 0 (p2 == rb, which mem_infer_dir maps to orientation 3 or 2).
//
// mem_matesw_skip is the interval form (mem_matesw_skip_query_init + mem_matesw_skip_scan), and
// mem_sam_pe_batch_pre answers repeat mate lists with mem_matesw_skip_probe over the sorted rb
// values; both are compared against the same oracle, including regions outside [0, 2*l_pac),
// windows with low <= 0 or low > high, and tiny references where most positions are near the
// strand boundary.
//
// Fixtures are built in-memory; no test data files are read.
#include "doctest/doctest.h"
#include "bwamem.h"
#include <algorithm>
#include <cstdint>
#include <ostream>
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

unsigned oracle_bits(int64_t l_pac, int64_t rb, std::vector<mem_alnreg_t> &regs, const mem_pestat_t pes[4])
{
    mem_alnreg_v ma = {regs.size(), regs.size(), regs.empty() ? nullptr : regs.data(), 0};
    int want[4];
    oracle(l_pac, rb, &ma, pes, want);
    return (unsigned)(want[0] | want[1] << 1 | want[2] << 2 | want[3] << 3);
}

// The first case the interval forms got wrong, for the failure message.
struct Mismatch {
    int64_t l_pac = 0, rb = 0;
    size_t n = 0;
    unsigned want = 0, scan = 0, probe = 0;
    bool set = false;
};

std::ostream &operator<<(std::ostream &os, const Mismatch &m)
{
    if (!m.set) return os << "none";
    return os << "l_pac=" << m.l_pac << " rb=" << m.rb << " n=" << m.n << " want=" << m.want
              << " scan=" << m.scan << " probe=" << m.probe;
}

// The scan and the sorted probe of the interval form, each against the oracle; returns whether
// both matched (so a caller can count mismatches instead of asserting per case) and records the
// first mismatch in *first.
bool interval_forms_match(int64_t l_pac, int64_t rb, std::vector<mem_alnreg_t> &regs,
                          const mem_pestat_t pes[4], Mismatch *first)
{
    const unsigned want = oracle_bits(l_pac, rb, regs, pes);
    mem_matesw_skip_query_t q;
    mem_matesw_skip_query_init(l_pac, rb, pes, &q);
    const unsigned scan = mem_matesw_skip_scan(&q, regs.data(), regs.size());
    std::vector<int64_t> sorted;
    for (const mem_alnreg_t &m : regs) sorted.push_back(m.rb);
    std::sort(sorted.begin(), sorted.end());
    const unsigned probe = mem_matesw_skip_probe(&q, sorted.data(), sorted.size());
    const bool ok = scan == want && probe == want;
    if (!ok && !first->set) {
        first->l_pac = l_pac, first->rb = rb, first->n = regs.size();
        first->want = want, first->scan = scan, first->probe = probe, first->set = true;
    }
    return ok;
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

TEST_CASE("mate rescue skip: the interval scan and the sorted probe match the oracle"
          * doctest::test_suite("unit/pair")) {
    std::mt19937_64 rng(0x1a7e5a1);
    auto pick = [&](int64_t lo, int64_t hi) {
        return std::uniform_int_distribution<int64_t>(lo, hi)(rng);
    };
    int mismatches = 0, trials = 0;
    Mismatch first;
    for (int iter = 0; iter < 15000; ++iter) {
        // Tiny references put most positions next to the strand boundary or a reference end;
        // hg38-sized ones exercise coordinates past 32 bits.
        const int64_t l_pac = iter % 3 == 0 ? pick(1, 40) : iter % 3 == 1 ? pick(41, 5000) : 3100000000LL;
        const int64_t t = 2 * l_pac - 1;
        mem_pestat_t pes[4] = {};
        for (int r = 0; r < 4; ++r) {
            pes[r].failed = pick(0, 4) == 0;
            pes[r].low = (int)pick(-60, 300);
            pes[r].high = pes[r].low + (int)pick(-30, 600);   // sometimes empty (high < low)
        }
        const int n = (int)pick(0, 48);
        std::vector<mem_alnreg_t> regs;
        // Several anchors per mate list, like the repeat lists the sorted probe serves.
        for (int k = 0; k < 8; ++k) {
            const int64_t rb = pick(0, 9) == 0 ? pick(-100, t + 100) : pick(0, t);
            if (k == 0) {
                for (int i = 0; i < n; ++i) {
                    const int kind = (int)pick(0, 4);
                    const int64_t b2 = kind == 0 ? rb + pick(-700, 700)
                                     : kind == 1 ? t - rb + pick(-700, 700)
                                     : kind == 2 ? l_pac + pick(-3, 2)          // the strand boundary
                                     : kind == 3 ? pick(-200, t + 200)         // outside [0, 2*l_pac) too
                                                 : pick(0, t);
                    regs.push_back(region_at(b2));
                }
            }
            ++trials;
            if (!interval_forms_match(l_pac, rb, regs, pes, &first)) ++mismatches;
        }
    }
    CHECK(trials == 120000);
    CAPTURE(first);
    CHECK(mismatches == 0);
}

TEST_CASE("mate rescue skip: each orientation's interval ends, one region at a time"
          * doctest::test_suite("unit/pair")) {
    // A region exactly at, and one past, every end of every orientation's interval, for windows
    // that start below, at and above distance 0: the off-by-one forms (rb vs rb + 1 on the
    // ahead side, min(rb, rb - low) vs rb - low on the behind side, the l_pac clip) each move
    // one of these regions across the end.
    const int64_t L_PAC = 100000;
    const int64_t T = 2 * L_PAC - 1;
    int mismatches = 0, trials = 0;
    Mismatch first;
    for (const int64_t rb : {(int64_t)0, (int64_t)250, L_PAC - 300, L_PAC - 1, L_PAC, L_PAC + 1,
                             L_PAC + 299, T - 250, T}) {
        for (const int low : {-5, -1, 0, 1, 2, 100}) {
            for (const int width : {-1, 0, 1, 400}) {
                mem_pestat_t pes[4] = {};
                for (int r = 0; r < 4; ++r) pes[r].low = low, pes[r].high = low + width;
                mem_matesw_skip_query_t q;
                mem_matesw_skip_query_init(L_PAC, rb, pes, &q);
                for (int r = 0; r < 4; ++r) {
                    for (const int64_t end : {q.lo[r], q.hi[r]}) {
                        for (const int64_t d : {-1, 0, 1}) {
                            std::vector<mem_alnreg_t> regs = {region_at(end + d)};
                            ++trials;
                            if (!interval_forms_match(L_PAC, rb, regs, pes, &first)) ++mismatches;
                        }
                    }
                }
            }
        }
    }
    CHECK(trials == 9 * 6 * 4 * 4 * 2 * 3);
    CAPTURE(first);
    CHECK(mismatches == 0);
}
