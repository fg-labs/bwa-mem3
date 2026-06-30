// Unit tests for smem_dedup_inplace: adjacent fully-identical SMEM deduplication.
// The function compacts a sorted SMEM array in-place, removing duplicates
// (entries identical on rid,m,n,k,l,s) from adjacent positions.
//
// Precondition tested here: identical SMEMs are adjacent (true after
// sortSMEMs + per-read intv_lt1 sort in mem_collect_smem).

#include "doctest/doctest.h"
#include "smem_dedup.h"

#include <cstring>
#include <vector>

namespace {

// Build a SMEM with all fields set explicitly.
SMEM make_smem(uint32_t rid, uint32_t m, uint32_t n,
               int64_t k, int64_t l, int64_t s)
{
    SMEM sm;
    memset(&sm, 0, sizeof(sm));
    sm.rid = rid;
    sm.m   = m;
    sm.n   = n;
    sm.k   = k;
    sm.l   = l;
    sm.s   = s;
    return sm;
}

// Run dedup and return the kept (rid,m,n,k,l,s) tuples as a vector of SMEMs.
std::vector<SMEM> run_dedup(std::vector<SMEM> arr)
{
    int64_t new_n = smem_dedup_inplace(arr.data(),
                                       static_cast<int64_t>(arr.size()));
    arr.resize(static_cast<size_t>(new_n));
    return arr;
}

// Compare two SMEMs on all six fields.
bool smem_eq(const SMEM &a, const SMEM &b)
{
    return a.rid == b.rid &&
           a.m   == b.m   &&
           a.n   == b.n   &&
           a.k   == b.k   &&
           a.l   == b.l   &&
           a.s   == b.s;
}

}  // namespace

// ---------------------------------------------------------------------------
// Edge cases: n=0 and n=1
// ---------------------------------------------------------------------------

TEST_CASE("n=0 returns 0"
          * doctest::test_suite("unit/smem_dedup"))
{
    // nullptr is safe at n=0: the `n <= 1` early-return fires before the array is
    // touched. The real caller (mem_collect_smem) can also pass an unpopulated array
    // when a read has zero SMEMs, so the empty case must short-circuit cleanly.
    int64_t result = smem_dedup_inplace(nullptr, 0);
    CHECK(result == 0);
}

TEST_CASE("n=1 returns 1 unchanged"
          * doctest::test_suite("unit/smem_dedup"))
{
    SMEM sm = make_smem(0, 5, 20, 100, 110, 2);
    int64_t result = smem_dedup_inplace(&sm, 1);
    REQUIRE(result == 1);
    CHECK(smem_eq(sm, make_smem(0, 5, 20, 100, 110, 2)));
}

// ---------------------------------------------------------------------------
// Adjacent duplicates are removed
// ---------------------------------------------------------------------------

TEST_CASE("two adjacent identical SMEMs -> one kept"
          * doctest::test_suite("unit/smem_dedup"))
{
    std::vector<SMEM> arr = {
        make_smem(0, 5, 20, 100, 110, 2),
        make_smem(0, 5, 20, 100, 110, 2),
    };
    auto result = run_dedup(arr);
    REQUIRE(result.size() == 1);
    CHECK(smem_eq(result[0], make_smem(0, 5, 20, 100, 110, 2)));
}

TEST_CASE("two identical then one distinct -> 2 kept"
          * doctest::test_suite("unit/smem_dedup"))
{
    SMEM dup  = make_smem(0, 5, 20, 100, 110, 2);
    SMEM uniq = make_smem(0, 7, 25, 200, 210, 1);
    std::vector<SMEM> arr = {dup, dup, uniq};
    auto result = run_dedup(arr);
    REQUIRE(result.size() == 2);
    CHECK(smem_eq(result[0], dup));
    CHECK(smem_eq(result[1], uniq));
}

TEST_CASE("triple duplicate -> one kept"
          * doctest::test_suite("unit/smem_dedup"))
{
    SMEM sm = make_smem(1, 3, 15, 50, 60, 5);
    std::vector<SMEM> arr = {sm, sm, sm};
    auto result = run_dedup(arr);
    REQUIRE(result.size() == 1);
    CHECK(smem_eq(result[0], sm));
}

// ---------------------------------------------------------------------------
// All distinct entries preserved
// ---------------------------------------------------------------------------

TEST_CASE("all distinct entries are all preserved"
          * doctest::test_suite("unit/smem_dedup"))
{
    std::vector<SMEM> arr = {
        make_smem(0, 1, 10, 10, 20, 3),
        make_smem(0, 2, 12, 30, 40, 2),
        make_smem(0, 4, 18, 80, 90, 1),
    };
    auto result = run_dedup(arr);
    REQUIRE(result.size() == 3);
    CHECK(smem_eq(result[0], arr[0]));
    CHECK(smem_eq(result[1], arr[1]));
    CHECK(smem_eq(result[2], arr[2]));
}

// ---------------------------------------------------------------------------
// rid boundary: same (m,n,k,l,s) but different rid -> NOT merged
// ---------------------------------------------------------------------------

TEST_CASE("same m,n,k,l,s but different rid -> both kept"
          * doctest::test_suite("unit/smem_dedup"))
{
    SMEM a = make_smem(0, 5, 20, 100, 110, 2);
    SMEM b = make_smem(1, 5, 20, 100, 110, 2);  // rid differs
    std::vector<SMEM> arr = {a, b};
    auto result = run_dedup(arr);
    REQUIRE(result.size() == 2);
    CHECK(smem_eq(result[0], a));
    CHECK(smem_eq(result[1], b));
}

// ---------------------------------------------------------------------------
// Non-adjacent duplicate (a different entry separates them): the precondition
// guarantees identical SMEMs are adjacent after sorting, so a non-adjacent
// re-occurrence is treated as a NEW distinct entry — dedup only removes
// adjacent duplicates.
// ---------------------------------------------------------------------------

TEST_CASE("non-adjacent duplicate treated as distinct (precondition: sorted-adjacent)"
          * doctest::test_suite("unit/smem_dedup"))
{
    SMEM sm_a = make_smem(0, 5, 20, 100, 110, 2);
    SMEM sm_b = make_smem(0, 7, 25, 200, 210, 1);  // different
    // sm_a appears at positions 0 and 2, separated by sm_b at 1
    std::vector<SMEM> arr = {sm_a, sm_b, sm_a};
    auto result = run_dedup(arr);
    // All three kept: the second sm_a at [2] differs from the last-kept sm_b at [1]
    REQUIRE(result.size() == 3);
    CHECK(smem_eq(result[0], sm_a));
    CHECK(smem_eq(result[1], sm_b));
    CHECK(smem_eq(result[2], sm_a));
}

// ---------------------------------------------------------------------------
// Partial equality: one field differs -> not a duplicate
// ---------------------------------------------------------------------------

TEST_CASE("entries differing only on k are distinct"
          * doctest::test_suite("unit/smem_dedup"))
{
    SMEM a = make_smem(0, 5, 20, 100, 110, 2);
    SMEM b = make_smem(0, 5, 20, 999, 110, 2);  // k differs
    std::vector<SMEM> arr = {a, b};
    auto result = run_dedup(arr);
    REQUIRE(result.size() == 2);
}

TEST_CASE("entries differing only on s are distinct"
          * doctest::test_suite("unit/smem_dedup"))
{
    SMEM a = make_smem(0, 5, 20, 100, 110, 2);
    SMEM b = make_smem(0, 5, 20, 100, 110, 3);  // s differs
    std::vector<SMEM> arr = {a, b};
    auto result = run_dedup(arr);
    REQUIRE(result.size() == 2);
}

// ---------------------------------------------------------------------------
// Mixed run: duplicates interspersed with unique entries
// ---------------------------------------------------------------------------

TEST_CASE("mixed: dup dup unique dup dup unique -> 4 kept"
          * doctest::test_suite("unit/smem_dedup"))
{
    SMEM x = make_smem(0, 1, 10, 10, 20, 5);
    SMEM y = make_smem(0, 2, 15, 50, 60, 3);
    SMEM z = make_smem(0, 3, 18, 90, 95, 1);
    // x x y z z y  -> after dedup: x y z y
    std::vector<SMEM> arr = {x, x, y, z, z, y};
    auto result = run_dedup(arr);
    REQUIRE(result.size() == 4);
    CHECK(smem_eq(result[0], x));
    CHECK(smem_eq(result[1], y));
    CHECK(smem_eq(result[2], z));
    CHECK(smem_eq(result[3], y));
}
