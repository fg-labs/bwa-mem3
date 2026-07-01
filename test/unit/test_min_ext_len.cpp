// Unit tests for the --min-ext-len short-seed extension filter
// (mem_chain_drop_short_seeds). In a chain that still holds a seed >=
// min_ext_len (an anchor), seeds shorter than the threshold are dropped so they
// are never extended. A chain with NO such anchor -- every seed shorter than the
// threshold -- is left untouched, so the filter never empties a chain and never
// loses a read (recall-safety invariant). min_ext_len <= 0 is a no-op (default),
// which keeps output byte-identical to baseline.

#include "doctest/doctest.h"
#include "bwamem.h"

#include <vector>

namespace {

// Build a chain from a list of seed lengths, run the filter, and return the
// lengths of the surviving seeds in order. Backing storage is owned by the
// caller-supplied vector so the chain's seed pointer stays valid for the call.
std::vector<int> drop_short(std::vector<mem_seed_t> &seeds, int min_ext_len)
{
    mem_chain_t c{};
    c.n = static_cast<int32_t>(seeds.size());
    c.m = static_cast<int32_t>(seeds.size());
    c.seeds = seeds.empty() ? nullptr : seeds.data();

    const int new_n = mem_chain_drop_short_seeds(&c, min_ext_len);
    CHECK(new_n == c.n);

    std::vector<int> out;
    out.reserve(static_cast<size_t>(c.n));
    for (int i = 0; i < c.n; ++i) out.push_back(c.seeds[i].len);
    return out;
}

std::vector<mem_seed_t> make_seeds(const std::vector<int> &lens)
{
    std::vector<mem_seed_t> seeds(lens.size());  // mem_seed_t ctor sets defaults
    for (size_t i = 0; i < lens.size(); ++i) {
        seeds[i].len = lens[i];
        seeds[i].rbeg = static_cast<int64_t>(i) * 1000;  // distinct, to track order
        seeds[i].qbeg = static_cast<int32_t>(i);
    }
    return seeds;
}

}  // namespace

TEST_CASE("min_ext_len=0 is a no-op (byte-identical default)"
          * doctest::test_suite("unit/min_ext_len"))
{
    auto seeds = make_seeds({10, 50, 19, 100});
    CHECK(drop_short(seeds, 0) == std::vector<int>{10, 50, 19, 100});
}

TEST_CASE("negative threshold is also a no-op"
          * doctest::test_suite("unit/min_ext_len"))
{
    auto seeds = make_seeds({10, 50, 19, 100});
    CHECK(drop_short(seeds, -5) == std::vector<int>{10, 50, 19, 100});
}

TEST_CASE("drops seeds shorter than threshold and preserves order"
          * doctest::test_suite("unit/min_ext_len"))
{
    auto seeds = make_seeds({10, 50, 19, 100, 30});
    // 10 and 19 dropped; 50, 100, 30 kept in original order.
    CHECK(drop_short(seeds, 30) == std::vector<int>{50, 100, 30});
}

TEST_CASE("preserves full seed records, not just lengths"
          * doctest::test_suite("unit/min_ext_len"))
{
    auto seeds = make_seeds({10, 50, 19, 100});  // rbeg = 0, 1000, 2000, 3000
    mem_chain_t c{};
    c.n = c.m = static_cast<int32_t>(seeds.size());
    c.seeds = seeds.data();

    REQUIRE(mem_chain_drop_short_seeds(&c, 30) == 2);
    // Survivors are the original index-1 and index-3 seeds, copied whole:
    CHECK(c.seeds[0].len == 50);   CHECK(c.seeds[0].rbeg == 1000);
    CHECK(c.seeds[1].len == 100);  CHECK(c.seeds[1].rbeg == 3000);
}

TEST_CASE("boundary: seed length == threshold is kept (>=)"
          * doctest::test_suite("unit/min_ext_len"))
{
    auto seeds = make_seeds({29, 30, 31});
    CHECK(drop_short(seeds, 30) == std::vector<int>{30, 31});
}

TEST_CASE("all seeds shorter than threshold -> chain left intact (non-emptying)"
          * doctest::test_suite("unit/min_ext_len"))
{
    // No seed >= 30 means there is no anchor whose extension covers the short
    // seeds, so the filter must leave the chain untouched rather than empty it
    // (recall-safety: an all-short chain is a read's only evidence).
    auto seeds = make_seeds({10, 19, 25});
    CHECK(drop_short(seeds, 30) == std::vector<int>{10, 19, 25});
}

TEST_CASE("all seeds at least threshold -> all kept"
          * doctest::test_suite("unit/min_ext_len"))
{
    auto seeds = make_seeds({30, 40, 150});
    CHECK(drop_short(seeds, 30) == std::vector<int>{30, 40, 150});
}

TEST_CASE("empty chain stays empty"
          * doctest::test_suite("unit/min_ext_len"))
{
    std::vector<mem_seed_t> seeds;  // n == 0
    CHECK(drop_short(seeds, 30).empty());
}
