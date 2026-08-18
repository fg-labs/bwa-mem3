// Unit tests for mem_chain_weight()'s saturating clamp (fg-labs/bwa-mem3#309).
//
// The clamp used to saturate at (1<<30)-1 while mem_chain_t::w is a 27-bit
// bitfield, so a weight in [2^27, 2^30) wrapped modulo 2^27 on store: a very
// heavy chain came back with a TINY w, which could lose it the
// `c->w < opt->min_chain_weight` gate or a `drop_ratio` shadowing comparison it
// should have won. The clamp now saturates at MEM_CHAIN_W_MAX, so an
// out-of-range weight stays large.
//
// This is unreachable on real data -- chain weight accumulates non-overlapping
// seed spans and takes min(query, reference) coverage, so it is bounded by the
// query length, and saturating needs a single ~134 Mbp read. That is exactly
// why it is worth a unit test: nothing else in the tree can exercise it, so
// without this the relationship between the clamp and the field is asserted
// only by a comment. mem_chain_weight() takes a plain mem_chain_t*, so a
// synthetic chain reaches it directly.

#include "doctest/doctest.h"
#include "bwamem.h"

#include <cstdint>
#include <vector>

namespace {

// Build a chain of non-overlapping seeds, each `len` long, laid end to end on
// both axes so query and reference coverage both equal n_seeds * len (the
// weight is the min of the two).
int weight_of(std::vector<mem_seed_t> &seeds, int n_seeds, int64_t len)
{
    seeds.clear();
    seeds.resize(static_cast<size_t>(n_seeds));
    for (int i = 0; i < n_seeds; ++i) {
        seeds[static_cast<size_t>(i)].qbeg = static_cast<int32_t>(i * len);
        seeds[static_cast<size_t>(i)].rbeg = i * len;
        seeds[static_cast<size_t>(i)].len  = static_cast<int32_t>(len);
    }

    mem_chain_t c{};
    c.n = static_cast<int32_t>(n_seeds);
    c.m = static_cast<int32_t>(n_seeds);
    c.seeds = seeds.data();
    return mem_chain_weight(&c);
}

}  // namespace

// An INDEPENDENT oracle for the field maximum: the width `27` and the value
// `(1u << 27) - 1` are spelled out here rather than read from MEM_CHAIN_W_MAX,
// so the tests below cannot silently track a wrong edit to that shared contract.
// The static_asserts pin MEM_CHAIN_W_BITS / MEM_CHAIN_W_MAX to this oracle at
// compile time; if the header width ever changes, this file fails to build
// rather than quietly re-deriving its own expectations from the bug.
constexpr uint32_t expected_max = (uint32_t{1} << 27) - 1;
static_assert(MEM_CHAIN_W_BITS == 27, "chain-weight field width drifted from 27 bits");
static_assert(MEM_CHAIN_W_MAX == expected_max, "MEM_CHAIN_W_MAX no longer equals (1<<27)-1");

TEST_CASE("mem_chain_weight: ordinary weights are returned unclamped"
          * doctest::test_suite("unit/chain_weight_clamp")) {
    std::vector<mem_seed_t> seeds;
    SUBCASE("three 100 bp seeds cover 300 bp on both axes") {
        CHECK(weight_of(seeds, 3, 100) == 300);
    }
    SUBCASE("a single seed the length of a long read") {
        CHECK(weight_of(seeds, 1, 25000) == 25000);
    }
}

TEST_CASE("mem_chain_weight: saturates at MEM_CHAIN_W_MAX rather than wrapping"
          * doctest::test_suite("unit/chain_weight_clamp")) {
    std::vector<mem_seed_t> seeds;

    SUBCASE("just under the field maximum: returned as-is") {
        const int64_t just_under = static_cast<int64_t>(expected_max) - 1;
        CHECK(weight_of(seeds, 1, just_under) == static_cast<int>(just_under));
    }

    SUBCASE("exactly the field maximum: still exact, not clamped off by one") {
        CHECK(weight_of(seeds, 1, static_cast<int64_t>(MEM_CHAIN_W_MAX))
              == static_cast<int>(MEM_CHAIN_W_MAX));
    }

    SUBCASE("over the maximum: saturates") {
        // Two seeds of 2^26 each sum to 2^27, which is MEM_CHAIN_W_MAX + 1 --
        // the smallest value the old (1<<30)-1 clamp let through to wrap, and it
        // wrapped to 0.
        const int over = weight_of(seeds, 2, static_cast<int64_t>(1) << 26);
        CHECK(over == static_cast<int>(MEM_CHAIN_W_MAX));
    }

    SUBCASE("far over the maximum: still saturates, not accidentally right") {
        // The old clamp returned 2^29 here because this value was below
        // (1<<30)-1. Storing 2^29 in the 27-bit field truncated it to 0.
        // Pin saturation directly.
        const int way_over = weight_of(seeds, 8, static_cast<int64_t>(1) << 26);
        CHECK(way_over == static_cast<int>(MEM_CHAIN_W_MAX));
    }
}

TEST_CASE("mem_chain_weight: a saturated weight survives the bitfield round trip"
          * doctest::test_suite("unit/chain_weight_clamp")) {
    // The actual defect was on STORE, not in the return value: the caller does
    // `c->w = mem_chain_weight(c)`. Reproduce that assignment and check the
    // value read back is still large, which is what the min_chain_weight gate
    // and the drop_ratio comparison depend on.
    std::vector<mem_seed_t> seeds;
    const int w = weight_of(seeds, 2, static_cast<int64_t>(1) << 26);

    mem_chain_t c{};
    c.w = static_cast<uint32_t>(w);
    // Copy out before comparing: doctest's CHECK binds its operands by
    // reference, and a reference cannot bind to a bitfield.
    const uint32_t stored = c.w;
    CHECK(stored == expected_max);

    // What the OLD clamp would have stored for the same chain. It returned the
    // true weight (2^27) unchanged, since 2^27 < (1<<30), and the narrowing
    // store then took it modulo 2^27 -- so the heaviest chain in this test
    // would have been recorded as the lightest possible one. Computed rather
    // than assigned out of range: `c.w = 1u << 27` is a constant that does not
    // fit the field and trips -Wbitfield-constant-conversion.
    const uint32_t old_stored = static_cast<uint32_t>(1u << 27) & expected_max;
    CHECK(old_stored == 0u);
    CHECK(stored != old_stored);
}
