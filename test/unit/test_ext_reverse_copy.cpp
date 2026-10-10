// test/unit/test_ext_reverse_copy.cpp
//
// ext_reverse_copy (src/ext_reverse_copy.h) writes dst[i] = src[n - 1 - i]: the
// reversed left-extension windows of stage_seed_extension. The compiler may
// vectorize it in blocks with a byte-loop tail, so every length around the
// block sizes, at every destination offset mod 64, is checked against
// std::reverse_copy, with guard bytes on either side of the destination left
// untouched.

#include <algorithm>
#include <cstdint>
#include <vector>

#include "doctest/doctest.h"
#include "ext_reverse_copy.h"

TEST_CASE("ext_reverse_copy: equals std::reverse_copy at lengths around the block sizes" * doctest::test_suite("unit/extension")) {
    std::vector<int64_t> lens;
    for (int64_t n = 0; n <= 70; ++n) lens.push_back(n);
    for (int64_t n : {127, 128, 129, 255, 256, 1000, 4097}) lens.push_back(n);
    uint32_t x = 12345;
    for (int64_t n : lens) {
        // src is exactly n bytes, so an out-of-range read is out of bounds (a
        // sanitizer build reports it) rather than landing in readable padding.
        std::vector<uint8_t> src(n);
        for (auto &c : src) { x = x * 1103515245u + 12345u; c = (uint8_t)(x >> 24); }
        for (int64_t off = 0; off < 64; ++off) {
            // dst starts 16 + off bytes in, so there are always guard bytes before it.
            std::vector<uint8_t> got(n + 96, 0xAA), want(n + 96, 0xAA);
            std::reverse_copy(src.begin(), src.end(), want.begin() + 16 + off);
            ext_reverse_copy(got.data() + 16 + off, src.data(), n);
            CHECK_MESSAGE(got == want, "n=" << n << " off=" << off);
        }
    }
}

TEST_CASE("ext_reverse_copy: a non-positive length copies nothing" * doctest::test_suite("unit/extension")) {
    const std::vector<uint8_t> src(64, 0x55);
    for (int64_t n : {int64_t(0), int64_t(-1), int64_t(-64), INT64_MIN}) {
        std::vector<uint8_t> got(96, 0xAA);
        ext_reverse_copy(got.data() + 16, src.data(), n);
        CHECK_MESSAGE(got == std::vector<uint8_t>(96, 0xAA), "n=" << n);
    }
}
