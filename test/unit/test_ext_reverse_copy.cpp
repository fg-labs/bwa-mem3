// test/unit/test_ext_reverse_copy.cpp
//
// ext_reverse_copy (src/ext_reverse_copy.h) writes dst[i] = src[n - 1 - i]: the
// reversed left-extension windows of stage_seed_extension. On aarch64 it copies
// 16 bytes at a time with a byte-loop tail, so every length around the block
// size is checked against the byte loop, with the bytes on either side of the
// destination left untouched.

#include <cstdint>
#include <vector>

#include "doctest/doctest.h"
#include "ext_reverse_copy.h"

TEST_CASE("ext_reverse_copy: equals the byte loop for every length" * doctest::test_suite("unit/extension")) {
    std::vector<int64_t> lens;
    for (int64_t n = 0; n <= 70; ++n) lens.push_back(n);
    for (int64_t n : {127, 128, 129, 255, 256, 1000, 4097}) lens.push_back(n);
    uint32_t x = 12345;
    for (int64_t n : lens) {
        std::vector<uint8_t> src(n + 32), got(n + 32, 0xAA), want(n + 32, 0xAA);
        for (auto &c : src) { x = x * 1103515245u + 12345u; c = (uint8_t)(x >> 24); }
        const uint8_t *s = src.data() + 16;   // src[-16..n+16) readable, only [0, n) may be read
        for (int64_t i = 0; i < n; ++i) want[16 + i] = s[n - 1 - i];
        ext_reverse_copy(got.data() + 16, s, n);
        CHECK_MESSAGE(got == want, "n=" << n);
    }
}
