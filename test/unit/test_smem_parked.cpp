// test/unit/test_smem_parked.cpp
//
// Unit tests for LockstepParked (src/smem_parked.h), the in-order emitter the
// lockstep SMEM driver uses once it refills a slot as soon as its read finishes.
// Reads may finish in any order inside the window; the output must be every
// read's matches, in input order, each read's matches in their own order. The
// whole-driver parity against the scalar walker is
// test/smem_lockstep_parity_test.cpp (Case 15 drives out-of-order finishing).

#include "doctest/doctest.h"
#include "../../src/smem_parked.h"

#include <vector>

namespace {

// A read's matches: `count` SMEMs tagged with the read (rid) and their position
// within the read (m), so the output's order can be checked field by field.
std::vector<SMEM> matches_of(int32_t read, int32_t count)
{
    std::vector<SMEM> v((size_t)count);
    for (int32_t j = 0; j < count; j++) {
        v[(size_t)j] = SMEM();
        v[(size_t)j].rid = (uint32_t)read;
        v[(size_t)j].m = (uint32_t)j;
        v[(size_t)j].k = 1000 * read + j;
    }
    return v;
}

struct Emitter {
    LockstepParked parked;
    std::vector<SMEM> out = std::vector<SMEM>(4096);
    int64_t n_out = 0;

    void finish(int32_t read, int32_t count)
    {
        const std::vector<SMEM> m = matches_of(read, count);
        parked.finish(read, m.data(), count, out.data(), &n_out);
    }
};

// The output is reads 0..n-1 in order, read r with counts[r] matches.
void check_in_order(const Emitter &e, const std::vector<int32_t> &counts)
{
    int64_t at = 0;
    for (size_t r = 0; r < counts.size(); r++) {
        for (int32_t j = 0; j < counts[r]; j++, at++) {
            REQUIRE(at < e.n_out);
            CAPTURE(r);
            CAPTURE(j);
            CHECK(e.out[(size_t)at].rid == (uint32_t)r);
            CHECK(e.out[(size_t)at].m == (uint32_t)j);
            CHECK(e.out[(size_t)at].k == (int64_t)(1000 * r + j));
        }
    }
    CHECK(at == e.n_out);
}

}  // namespace

TEST_CASE("LockstepParked: the window is 8 reads per slot, rounded up to a power of two"
          * doctest::test_suite("unit/smem"))
{
    CHECK(LockstepParked::window_for(1) == 8);
    CHECK(LockstepParked::window_for(3) == 32);
    CHECK(LockstepParked::window_for(4) == 32);
    CHECK(LockstepParked::window_for(5) == 64);
    CHECK(LockstepParked::window_for(16) == 128);
    CHECK(LockstepParked::window_for(SMEM_LOCKSTEP_N_MAX) >= 8 * SMEM_LOCKSTEP_N_MAX);
}

TEST_CASE("LockstepParked: reads finishing out of order leave in input order"
          * doctest::test_suite("unit/smem"))
{
    Emitter e;
    e.parked.reset(16);
    // 2, 1, 3 park; 0 releases all four. Read 1 has no matches.
    e.finish(2, 3);
    e.finish(1, 0);
    e.finish(3, 2);
    CHECK(e.n_out == 0);
    CHECK(e.parked.pending == 3);
    CHECK(e.parked.cursor == 0);
    e.finish(0, 4);
    CHECK(e.parked.pending == 0);
    CHECK(e.parked.cursor == 4);
    check_in_order(e, {4, 0, 3, 2});
}

TEST_CASE("LockstepParked: a read that is due drains only the consecutive parked reads"
          * doctest::test_suite("unit/smem"))
{
    Emitter e;
    e.parked.reset(16);
    e.finish(1, 2);
    e.finish(3, 1);   // read 2 not done: 3 must wait behind it
    e.finish(0, 1);
    CHECK(e.parked.cursor == 2);
    CHECK(e.parked.pending == 1);
    check_in_order(e, {1, 2});
    e.finish(2, 0);
    CHECK(e.parked.cursor == 4);
    CHECK(e.parked.pending == 0);
    check_in_order(e, {1, 2, 0, 1});
}

TEST_CASE("LockstepParked: in-order reads go straight out, empty ones included"
          * doctest::test_suite("unit/smem"))
{
    Emitter e;
    e.parked.reset(4);
    const std::vector<int32_t> counts = {0, 3, 0, 0, 1, 2};
    for (size_t r = 0; r < counts.size(); r++) {
        e.finish((int32_t)r, counts[r]);
        CHECK(e.parked.pending == 0);
    }
    check_in_order(e, counts);
}

TEST_CASE("LockstepParked: ring entries are reused for later windows"
          * doctest::test_suite("unit/smem"))
{
    // Window 4, 40 reads: each window of 4 finishes in reverse (3, 2, 1, 0 of the
    // window), so every ring entry is filled and drained ten times, and parked
    // reads at i and i + 4 never coexist (the driver's admits() keeps them apart).
    Emitter e;
    e.parked.reset(4);
    std::vector<int32_t> counts;
    for (int32_t r = 0; r < 40; r++) counts.push_back(r % 5);   // 0..4, empties included
    for (int32_t base = 0; base < 40; base += 4) {
        for (int32_t k = 3; k >= 0; k--) {
            const int32_t r = base + k;
            REQUIRE(e.parked.admits(r));
            e.finish(r, counts[(size_t)r]);
        }
        CHECK(e.parked.cursor == base + 4);
        CHECK(e.parked.pending == 0);
    }
    check_in_order(e, counts);
}

TEST_CASE("LockstepParked: admits() bounds inputs to the window ahead of the cursor"
          * doctest::test_suite("unit/smem"))
{
    Emitter e;
    e.parked.reset(8);
    CHECK(e.parked.admits(7));
    CHECK_FALSE(e.parked.admits(8));
    e.finish(0, 1);
    CHECK(e.parked.admits(8));
    CHECK_FALSE(e.parked.admits(9));
    e.finish(5, 1);   // parked: the cursor stays at 1
    CHECK_FALSE(e.parked.admits(9));
}

TEST_CASE("LockstepParked: reset() starts a fresh call after parked reads"
          * doctest::test_suite("unit/smem"))
{
    // A previous call must leave nothing behind: its matches do not reappear and
    // its ring entries read as empty.
    Emitter e;
    e.parked.reset(8);
    e.finish(1, 2);
    e.finish(2, 1);
    e.finish(0, 1);
    REQUIRE(e.parked.pending == 0);
    e.parked.reset(8);
    e.n_out = 0;
    CHECK(e.parked.cursor == 0);
    e.finish(1, 3);
    CHECK(e.n_out == 0);
    e.finish(0, 2);
    check_in_order(e, {2, 3});
}

TEST_CASE("LockstepParked: reset() drops reads still parked"
          * doctest::test_suite("unit/smem"))
{
    // The driver always drains before it returns, but reset() must not depend on
    // that: reads parked when it runs are dropped, their ring entries read as
    // empty, and their matches are not emitted in the next call.
    Emitter e;
    e.parked.reset(8);
    e.finish(1, 2);
    e.finish(3, 1);
    REQUIRE(e.parked.pending == 2);
    e.parked.reset(8);
    e.n_out = 0;
    CHECK(e.parked.pending == 0);
    e.finish(2, 2);   // parked, so the next emit walks the ring
    e.finish(0, 1);   // read 1's old entry must read as empty: the cursor stops at 1
    CHECK(e.parked.cursor == 1);
    CHECK(e.n_out == 1);
    e.finish(1, 3);
    CHECK(e.parked.cursor == 3);
    e.finish(3, 1);
    check_in_order(e, {1, 3, 2, 1});
}

TEST_CASE("LockstepParked: the parked buffer grows past its first allocation"
          * doctest::test_suite("unit/smem"))
{
    // Read 0 is last: 63 reads of 40 matches each park first (2520 SMEMs).
    Emitter e;
    e.parked.reset(64);
    std::vector<int32_t> counts(64, 40);
    counts[0] = 5;
    for (int32_t r = 63; r >= 1; r--) e.finish(r, counts[(size_t)r]);
    CHECK(e.parked.pending == 63);
    e.finish(0, counts[0]);
    check_in_order(e, counts);
}
