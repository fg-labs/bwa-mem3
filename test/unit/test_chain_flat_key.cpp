// test/unit/test_chain_flat_key.cpp — the flat chaining index's packed (pos, chain index) key
// and the parse of its chain-count cap.
//
// ChainFlat (bwamem.cpp) keeps one ascending uint64 array of chain_flat_pack(pos, ci) keys in
// place of a sorted position array plus a parallel chain-index array. Its byte-identity to the
// kbtree path rests on three properties pinned here against a plain reference model (a sorted
// vector of (pos, ci) pairs):
//
//   - packing round-trips pos and ci at the field bounds (pos up to CHAIN_FLAT_POS_LIMIT - 1,
//     ci up to CHAIN_FLAT_CI_LIMIT - 1), and the key order is the position order;
//   - chain_flat_rank_le(pos) is the number of chains whose position is <= pos, whatever the
//     chain indices (so the predecessor is the same chain the position-only index finds);
//   - the predecessor key carries both the position the tie test reads and the chain index.
//
// Only the position-range predicate (chain_flat_pos_ok) is tested here, not ChainFlat::add's
// decline on it. chain_flat_parse_cap is tested for the effective cap and its status: ok,
// malformed (the default is used) or clamped to CHAIN_FLAT_CI_LIMIT, so every packed chain
// index fits in its field.

#include <algorithm>
#include <random>
#include <utility>
#include <vector>

#include "doctest/doctest.h"
#include "chain_flat_key.h"

namespace {

// Reference: number of entries with position <= pos (positions are the sort key).
size_t ref_rank_le(const std::vector<std::pair<int64_t, uint32_t>> &v, int64_t pos) {
    size_t n = 0;
    for (const auto &e : v) n += e.first <= pos;
    return n;
}

}  // namespace

TEST_CASE("chain_flat_key: field bounds, round-trips and key order"
          * doctest::test_suite("unit/chain_flat_key")) {
    const uint32_t ci_max = (uint32_t)CHAIN_FLAT_CI_LIMIT - 1;  // 2^20 - 1
    const int64_t pos_max = CHAIN_FLAT_POS_LIMIT - 1;           // 2^43 - 1
    CHECK(ci_max == (1u << 20) - 1);
    CHECK(pos_max == ((int64_t)1 << 43) - 1);

    // Round-trips across the range, including pos = 0, an hg38-sized 2*l_pac and pos_max.
    const int64_t poss[] = {0, 1, 12345, (int64_t)6200000000LL, pos_max};
    const uint32_t cis[] = {0, 1, 511, ci_max};
    for (int64_t p : poss) {
        CHECK(chain_flat_pos_ok(p));
        for (uint32_t c : cis) {
            const uint64_t k = chain_flat_pack(p, c);
            CHECK(chain_flat_pos(k) == p);
            CHECK(chain_flat_ci(k) == c);
        }
    }
    CHECK(chain_flat_pack(0, 0) == 0);

    // The top key: pos = 2^43 - 1 with ci = 2^20 - 1 sets every bit below 2^63, sign bit clear.
    const uint64_t top = chain_flat_pack(pos_max, ci_max);
    CHECK(top == ((uint64_t)1 << 63) - 1);

    // Positions outside [0, 2^43) are rejected by the predicate (ChainFlat::add declines them).
    CHECK_FALSE(chain_flat_pos_ok(pos_max + 1));
    CHECK_FALSE(chain_flat_pos_ok(-1));
    CHECK_FALSE(chain_flat_pos_ok(INT64_MAX));
    CHECK_FALSE(chain_flat_pos_ok(INT64_MIN));

    // ci = 2^20 - 1 never spills into pos: the next position up is still a larger key.
    CHECK(chain_flat_pack(41, ci_max) < chain_flat_pack(42, 0));
    CHECK(chain_flat_pack(pos_max - 1, ci_max) < chain_flat_pack(pos_max, 0));
}

TEST_CASE("chain_flat_key: rank_le on empty, single-key and full-range arrays"
          * doctest::test_suite("unit/chain_flat_key")) {
    CHECK(chain_flat_rank_le(nullptr, 0, 5) == 0);
    const uint64_t one[] = {chain_flat_pack(10, 7)};
    CHECK(chain_flat_rank_le(one, 1, 9) == 0);
    CHECK(chain_flat_rank_le(one, 1, 10) == 1);
    CHECK(chain_flat_rank_le(one, 1, 11) == 1);

    const uint32_t ci_max = (uint32_t)CHAIN_FLAT_CI_MASK;
    const int64_t pos_max = CHAIN_FLAT_POS_LIMIT - 1;
    const uint64_t ends[] = {chain_flat_pack(0, ci_max), chain_flat_pack(pos_max, 0)};
    CHECK(chain_flat_rank_le(ends, 2, 0) == 1);
    CHECK(chain_flat_rank_le(ends, 2, pos_max - 1) == 1);
    CHECK(chain_flat_rank_le(ends, 2, pos_max) == 2);
}

TEST_CASE("chain_flat_key: rank_le counts every key at a tied position"
          * doctest::test_suite("unit/chain_flat_key")) {
    // Equal positions with different chain indices (ChainFlat itself declines a tie before
    // inserting one, but rank_le must still count them all and stop before the next position).
    const uint64_t k[] = {chain_flat_pack(5, 3), chain_flat_pack(7, 0), chain_flat_pack(7, 9),
                          chain_flat_pack(7, (uint32_t)CHAIN_FLAT_CI_MASK), chain_flat_pack(8, 1)};
    const size_t n = sizeof(k) / sizeof(k[0]);
    CHECK(chain_flat_rank_le(k, n, 4) == 0);
    CHECK(chain_flat_rank_le(k, n, 5) == 1);
    CHECK(chain_flat_rank_le(k, n, 6) == 1);
    CHECK(chain_flat_rank_le(k, n, 7) == 4);  // all three pos-7 keys, including ci = 2^20 - 1
    CHECK(chain_flat_rank_le(k, n, 8) == 5);
    // The predecessor of a probe at 7 is a pos-7 key, so the tie test sees the equal position.
    CHECK(chain_flat_pos(k[chain_flat_rank_le(k, n, 7) - 1]) == 7);
    // A probe one below the tie lands on the previous position, not on any pos-7 key.
    CHECK(chain_flat_pos(k[chain_flat_rank_le(k, n, 6) - 1]) == 5);
}

TEST_CASE("chain_flat_key: sorted insert + rank_le match a position-ordered reference"
          * doctest::test_suite("unit/chain_flat_key")) {
    // Replays ChainFlat's insert pattern: each new chain gets the next ci and is inserted at
    // rank_le(pos). Positions come from a small range (many ties), an hg38-sized range and a
    // range just below the limit; equal positions are allowed here (ChainFlat itself declines
    // them) to check that rank_le still counts every chain at that position. Sized to stay well
    // under 100 ms in an unoptimized build.
    std::mt19937_64 rng(20261025);
    const std::pair<int64_t, int64_t> ranges[] = {
        {0, 64}, {0, (int64_t)6200000000LL}, {CHAIN_FLAT_POS_LIMIT - 1000, CHAIN_FLAT_POS_LIMIT - 1}};
    for (const auto &rg : ranges) {
        for (int trial = 0; trial < 8; ++trial) {
            std::uniform_int_distribution<int64_t> d(rg.first, rg.second);
            std::vector<uint64_t> keys;
            std::vector<std::pair<int64_t, uint32_t>> ref;  // (pos, ci), stable by pos
            const int n = 1 + (int)(rng() % 200);
            for (int i = 0; i < n; ++i) {
                const int64_t pos = d(rng);
                const size_t r = chain_flat_rank_le(keys.data(), keys.size(), pos);
                REQUIRE(r == ref_rank_le(ref, pos));
                if (r > 0) {
                    // The predecessor: same chain and position as the reference's.
                    REQUIRE(chain_flat_ci(keys[r - 1]) == ref[r - 1].second);
                    REQUIRE(chain_flat_pos(keys[r - 1]) == ref[r - 1].first);
                }
                const uint32_t ci = (uint32_t)i;
                keys.insert(keys.begin() + (ptrdiff_t)r, chain_flat_pack(pos, ci));
                ref.insert(ref.begin() + (ptrdiff_t)r, {pos, ci});
            }
            REQUIRE(std::is_sorted(keys.begin(), keys.end()));
            for (size_t i = 0; i < keys.size(); ++i) {
                REQUIRE(chain_flat_pos(keys[i]) == ref[i].first);
                REQUIRE(chain_flat_ci(keys[i]) == ref[i].second);
            }
            // Probe positions that are not inserted as well.
            for (int q = 0; q < 20; ++q) {
                const int64_t pos = d(rng);
                REQUIRE(chain_flat_rank_le(keys.data(), keys.size(), pos) == ref_rank_le(ref, pos));
            }
        }
    }
}

TEST_CASE("chain_flat_key: parse_cap gives the effective cap and its status"
          * doctest::test_suite("unit/chain_flat_key")) {
    const int def = 512;
    const int lim = (int)CHAIN_FLAT_CI_LIMIT;  // 1048576
    struct Case {
        const char *in;
        int cap;
        chain_flat_cap_status_t status;
    };
    const Case cases[] = {
        {nullptr, def, CHAIN_FLAT_CAP_OK},          // unset: the default
        {"", def, CHAIN_FLAT_CAP_OK},               // empty: the default
        {"0", 0, CHAIN_FLAT_CAP_OK},                // every read to the kbtree
        {"512", 512, CHAIN_FLAT_CAP_OK},
        {"1048575", lim - 1, CHAIN_FLAT_CAP_OK},
        {"1048576", lim, CHAIN_FLAT_CAP_OK},        // exactly the ceiling: not clamped
        {"1048577", lim, CHAIN_FLAT_CAP_CLAMPED},
        {"2000000", lim, CHAIN_FLAT_CAP_CLAMPED},
        {"2147483648", lim, CHAIN_FLAT_CAP_CLAMPED},  // above INT_MAX
        {"99999999999999999999999999", lim, CHAIN_FLAT_CAP_CLAMPED},  // strtol ERANGE overflow
        {"-1", def, CHAIN_FLAT_CAP_MALFORMED},
        {"-99999999999999999999999999", def, CHAIN_FLAT_CAP_MALFORMED},  // negative overflow
        {"abc", def, CHAIN_FLAT_CAP_MALFORMED},
        {"12abc", def, CHAIN_FLAT_CAP_MALFORMED},
        {"not-a-number", def, CHAIN_FLAT_CAP_MALFORMED},
    };
    for (const Case &c : cases) {
        CAPTURE(c.in ? c.in : "(null)");
        const chain_flat_cap_t r = chain_flat_parse_cap(c.in, def);
        CHECK(r.cap == c.cap);
        CHECK(r.status == c.status);
        // Every effective cap keeps the largest packed chain index (cap - 1) in its field.
        CHECK(r.cap >= 0);
        CHECK((uint64_t)r.cap <= CHAIN_FLAT_CI_LIMIT);
    }
}
