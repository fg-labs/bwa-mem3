// test/unit/test_tie_break_id.cpp — the seed of the equal-score tie-breaks.
//
// mem_mark_primary_se (which equal-score region becomes primary) and mem_pair
// (which equal-score pairing is emitted) break ties with a hash of an id. bwa and
// bwa-mem2 pass the read's position in the input, so the same read can make a
// different choice when the input is reordered, subset or sharded.
// --tie-break name passes a hash of the read name instead.
//
// These cases pin the contract of mem_tie_break_id, the one place the id is
// chosen: the default returns the input-position id unchanged (the byte-identity
// guarantee), and name mode returns a value that depends on the name alone, fits
// the bound mem_pair's `int id << 8` needs, and is a fixed function of the name
// across builds (golden values), so shards aligned by different builds agree.
// test/regression/tie_break_name.sh drives the option end to end.

#include "doctest/doctest.h"
#include "bwamem.h"

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <set>

// mem_pair takes the id as an int and computes `id << 8`; a name-derived id must
// therefore stay below 2^23. Checked at compile time so widening the mask fails
// the build rather than overflowing at run time.
static_assert((MEM_TIE_BREAK_NAME_MASK << 8) <= (unsigned long long)INT32_MAX,
              "a name tie-break id must fit mem_pair's int id << 8");

namespace {
// mem_opt_init() allocates; free through a guard so a failing CHECK cannot leak.
struct OptGuard {
    mem_opt_t *o;
    OptGuard() : o(mem_opt_init()) {}
    ~OptGuard() { free(o); }
};
} // namespace

TEST_CASE("tie-break defaults to the input-position id"
          * doctest::test_suite("unit/tie_break")) {
    OptGuard g;
    CHECK(g.o->tie_break_mode == MEM_TIE_BREAK_INDEX);
    // Unchanged for any id, including ones past 2^23 and 2^32, so the default
    // path passes exactly what bwa passes.
    for (uint64_t id : {0ULL, 1ULL, 12345ULL, 0x7FFFFFULL, 0x800000ULL, 0x1234567890ULL}) {
        CHECK(mem_tie_break_id(g.o, "read/1", id) == id);
        CHECK(mem_tie_break_id(g.o, "another", id) == id);
    }
}

TEST_CASE("tie-break name: the id depends on the name, not the input position"
          * doctest::test_suite("unit/tie_break")) {
    OptGuard g;
    g.o->tie_break_mode = MEM_TIE_BREAK_NAME;
    const char *name = "A00744:46:HV3C3DSXX:2:2330:14281:12633";
    const uint64_t first = mem_tie_break_id(g.o, name, 0);
    for (uint64_t id : {1ULL, 2ULL, 999999ULL, 0x1234567890ULL})
        CHECK(mem_tie_break_id(g.o, name, id) == first);
    CHECK(first <= MEM_TIE_BREAK_NAME_MASK);
}

TEST_CASE("tie-break name: the id is a fixed function of the name"
          * doctest::test_suite("unit/tie_break")) {
    // Golden values: the hash is part of the output contract, since it decides
    // which tied placement name mode reports. Changing it would make two builds
    // place the same read differently, so any change must be deliberate and
    // documented. Regenerate only together with a release note.
    OptGuard g;
    g.o->tie_break_mode = MEM_TIE_BREAK_NAME;
    CHECK(mem_tie_break_id(g.o, "read1", 0) == 1317512ULL);
    CHECK(mem_tie_break_id(g.o, "A00744:46:HV3C3DSXX:2:2330:14281:12633", 0) == 5557930ULL);
    CHECK(mem_tie_break_id(g.o, "", 0) == 7727541ULL);
}

TEST_CASE("tie-break name: different names spread over the id range"
          * doctest::test_suite("unit/tie_break")) {
    // Names that differ only in their last field (Illumina y coordinate) must not
    // collapse onto a few ids, or reads sharing a tile would tie-break alike.
    OptGuard g;
    g.o->tie_break_mode = MEM_TIE_BREAK_NAME;
    std::set<uint64_t> ids;
    char name[64];
    const int n = 4096;
    for (int y = 0; y < n; ++y) {
        snprintf(name, sizeof(name), "A00744:46:HV3C3DSXX:2:2330:14281:%d", y);
        const uint64_t id = mem_tie_break_id(g.o, name, 0);
        CHECK(id <= MEM_TIE_BREAK_NAME_MASK);
        ids.insert(id);
    }
    // 4096 draws from 2^23 values: the expected number of collisions is about 1.
    CHECK(ids.size() >= (size_t)(n - 16));
}
