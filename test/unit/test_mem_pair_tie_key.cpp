// test/unit/test_mem_pair_tie_key.cpp — the key mem_pair hashes to break ties
// between equal-score pairings.
//
// bwa and bwa-mem2 compute `y ^ id<<8` with a signed int id (the pair's index in
// the input). From the 2^23-th pair on, `id<<8` overflows an int: undefined
// behaviour. The compilers this project ships with emit a 32-bit shift and then
// sign-extend the result to 64 bits, so released output already reflects those
// bits. mem_pair_tie_key spells them out in well-defined arithmetic.
//
// These cases pin the key below the overflow point (where it must equal the
// plain unsigned formula), at the boundary, past it, and for a negative id (a
// 64-bit index truncated to int). The expected values are written out by hand
// from the 32-bit-shift-then-sign-extend rule, so they are an oracle independent
// of the implementation.

#include "doctest/doctest.h"
#include "bwamem.h"

#include <cstdint>

TEST_CASE("mem_pair_tie_key: below 2^23 it is y ^ (id << 8)"
          * doctest::test_suite("unit/pair")) {
    const uint64_t ys[] = {0ULL, 1ULL, 0x0000000500000003ULL, 0xFFFFFFFFFFFFFFFFULL};
    const int ids[] = {0, 1, 2, 255, 12345, (1 << 23) - 1};
    for (uint64_t y : ys)
        for (int id : ids)
            CHECK(mem_pair_tie_key(y, id) == (y ^ ((uint64_t)id << 8)));
}

TEST_CASE("mem_pair_tie_key: past 2^23 the shift wraps in 32 bits and sign-extends"
          * doctest::test_suite("unit/pair")) {
    const uint64_t y = 0x0000000700000002ULL;
    // 2^23 << 8 = 0x80000000: the sign bit of a 32-bit int, so it sign-extends.
    CHECK(mem_pair_tie_key(y, 1 << 23) == (y ^ 0xFFFFFFFF80000000ULL));
    // 0x01234567 << 8 = 0x123456700; the 32-bit wrap drops the leading 1, and
    // 0x23456700 is positive, so nothing is sign-extended.
    CHECK(mem_pair_tie_key(y, 0x01234567) == (y ^ 0x0000000023456700ULL));
    // 0x00FFFFFF << 8 = 0xFFFFFF00: negative as an int32.
    CHECK(mem_pair_tie_key(y, 0x00FFFFFF) == (y ^ 0xFFFFFFFFFFFFFF00ULL));
    // INT32_MAX << 8 wraps to 0xFFFFFF00 as well.
    CHECK(mem_pair_tie_key(y, INT32_MAX) == (y ^ 0xFFFFFFFFFFFFFF00ULL));
}

TEST_CASE("mem_pair_tie_key: a negative id (truncated 64-bit index) is well defined"
          * doctest::test_suite("unit/pair")) {
    const uint64_t y = 0x0000000100000000ULL;
    // -1 << 8 = 0xFFFFFF00 = -256.
    CHECK(mem_pair_tie_key(y, -1) == (y ^ 0xFFFFFFFFFFFFFF00ULL));
    // INT32_MIN << 8 wraps to 0.
    CHECK(mem_pair_tie_key(y, INT32_MIN) == y);
}
