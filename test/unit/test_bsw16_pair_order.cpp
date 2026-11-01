// test/unit/test_bsw16_pair_order.cpp — the 16-bit kernel's per-call pair order
// (bsw16_len1_major_order, sortPairsLen16).
//
// Pinned here, on synthetic pairs:
//   - the key choice: len1-major only when strictly more than half of the pairs have
//     len1 < len2; ties, all-equal-length pairs, n = 0 and n = 1 go the stated way;
//   - the order sortPairsLen16 produces matches an independent oracle, a std::stable_sort of
//     the input on (major key, minor key), for both choices and for the diverted call, whose
//     len1 key is max(len1, len2); equal pairs keep input order (stable) and no pair is
//     dropped or duplicated (every slot's full record is checked);
//   - n = 0 leaves the array untouched.

#include <algorithm>
#include <cstdint>
#include <random>
#include <vector>

#include "doctest/doctest.h"
#include "bandedSWA.h"
#include "bwamem.h"

namespace {

SeqPair make_pair(int len1, int len2, int id)
{
    SeqPair sp = {};
    sp.len1 = len1; sp.len2 = len2; sp.id = id; sp.seqid = id; sp.regid = id;
    return sp;
}

// The order sortPairsLen16 must produce: stable on (major, minor).
std::vector<SeqPair> oracle(std::vector<SeqPair> v, bool len1_major, bool diverted)
{
    auto l1 = [&](const SeqPair &p) { return diverted ? std::max(p.len1, p.len2) : p.len1; };
    std::stable_sort(v.begin(), v.end(), [&](const SeqPair &a, const SeqPair &b) {
        const int am = len1_major ? l1(a) : a.len2, bm = len1_major ? l1(b) : b.len2;
        if (am != bm) return am < bm;
        const int an = len1_major ? a.len2 : l1(a), bn = len1_major ? b.len2 : l1(b);
        return an < bn;
    });
    return v;
}

void check_sorted(std::vector<SeqPair> in, bool diverted)
{
    const bool len1_major = bsw16_len1_major_order(in.data(), (int32_t) in.size());
    const std::vector<SeqPair> want = oracle(in, len1_major, diverted);
    std::vector<SeqPair> tmp(in.size() + 1);
    std::vector<int32_t> hist(MAX_SEQ_LEN16 + 1, 0x5a5a5a5a);
    sortPairsLen16(in.data(), (int32_t) in.size(), tmp.data(), hist.data(), diverted);
    REQUIRE(in.size() == want.size());
    for (size_t i = 0; i < in.size(); i++) {
        CAPTURE(i);
        CHECK(in[i].id == want[i].id);
        CHECK(in[i].len1 == want[i].len1);
        CHECK(in[i].len2 == want[i].len2);
    }
}

}  // namespace

TEST_CASE("bsw16 pair order: the key choice follows the len1 < len2 majority"
          * doctest::test_suite("unit/bwamem"))
{
    SUBCASE("n = 0 goes len2-major") {
        CHECK_FALSE(bsw16_len1_major_order(nullptr, 0));
    }
    SUBCASE("n = 1") {
        SeqPair shorter_ref = make_pair(30, 80, 0), longer_ref = make_pair(80, 30, 0), equal = make_pair(50, 50, 0);
        CHECK(bsw16_len1_major_order(&shorter_ref, 1));
        CHECK_FALSE(bsw16_len1_major_order(&longer_ref, 1));
        CHECK_FALSE(bsw16_len1_major_order(&equal, 1));   // len1 == len2 is not len1 < len2
    }
    SUBCASE("a tie goes len2-major") {
        std::vector<SeqPair> v = {make_pair(30, 80, 0), make_pair(80, 30, 1)};
        CHECK_FALSE(bsw16_len1_major_order(v.data(), 2));
        v = {make_pair(10, 20, 0), make_pair(20, 10, 1), make_pair(15, 25, 2), make_pair(25, 25, 3)};
        CHECK_FALSE(bsw16_len1_major_order(v.data(), 4));
    }
    SUBCASE("a strict majority of len1 < len2 goes len1-major") {
        std::vector<SeqPair> v = {make_pair(10, 20, 0), make_pair(20, 10, 1), make_pair(15, 25, 2)};
        CHECK(bsw16_len1_major_order(v.data(), 3));
        v.push_back(make_pair(40, 40, 3));                 // 2 of 4: a tie again
        CHECK_FALSE(bsw16_len1_major_order(v.data(), 4));
    }
    SUBCASE("a strict majority of len1 >= len2 goes len2-major") {
        std::vector<SeqPair> v = {make_pair(300, 200, 0), make_pair(250, 250, 1), make_pair(20, 90, 2)};
        CHECK_FALSE(bsw16_len1_major_order(v.data(), 3));
    }
}

TEST_CASE("bsw16 pair order: sortPairsLen16 is the stable two-key order"
          * doctest::test_suite("unit/bwamem"))
{
    std::mt19937 rng(20261101);
    SUBCASE("n = 0 leaves the array untouched") {
        SeqPair guard = make_pair(7, 9, 42);
        std::vector<SeqPair> tmp(1);
        std::vector<int32_t> hist(MAX_SEQ_LEN16 + 1, 0);
        sortPairsLen16(&guard, 0, tmp.data(), hist.data(), false);
        CHECK(guard.id == 42);
        CHECK(guard.len1 == 7);
        CHECK(guard.len2 == 9);
    }
    SUBCASE("n = 1") {
        check_sorted({make_pair(30, 80, 0)}, false);
        check_sorted({make_pair(80, 30, 0)}, true);
    }
    SUBCASE("len1 < len2 majority (short-read shape), with repeated keys") {
        for (int trial = 0; trial < 20; trial++) {
            std::vector<SeqPair> v;
            const int n = 1 + (int) (rng() % 200);
            for (int k = 0; k < n; k++) {
                int len2 = 60 + (int) (rng() % 4) * 10;          // few distinct len2: ties
                int len1 = (rng() % 8 == 0) ? len2 + (int) (rng() % 20) : (int) (rng() % len2);
                v.push_back(make_pair(len1, len2, k));
            }
            check_sorted(v, false);
        }
    }
    SUBCASE("len1 >= len2 majority (long-extension shape), with repeated keys") {
        for (int trial = 0; trial < 20; trial++) {
            std::vector<SeqPair> v;
            const int n = 1 + (int) (rng() % 300);
            for (int k = 0; k < n; k++) {
                int len2 = 100 + (int) (rng() % 5) * 50;
                int len1 = (rng() % 8 == 0) ? (int) (rng() % len2) : len2 + (int) (rng() % 4) * 25;
                v.push_back(make_pair(len1, len2, k));
            }
            check_sorted(v, false);
            check_sorted(v, true);   // the diverted call's max(len1, len2) key
        }
    }
    SUBCASE("an exact tie sorts len2-major") {
        std::vector<SeqPair> v = {make_pair(90, 10, 0), make_pair(10, 90, 1), make_pair(50, 10, 2),
                                  make_pair(10, 50, 3), make_pair(50, 50, 4), make_pair(20, 40, 5)};
        // 3 of 6 have len1 < len2: not a strict majority
        REQUIRE_FALSE(bsw16_len1_major_order(v.data(), (int32_t) v.size()));
        check_sorted(v, false);
    }
}
