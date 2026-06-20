// test/unit/test_neon_movemask.cpp
//
// Scalar-vs-NEON parity for the two hand-written NEON movemask helpers:
//
//   - neon_movemask_u8       (neon_utils.h)   -- 16x uint8,  -> uint16 mask
//   - _mm_movemask_epi16     (simd_compat.h)  -- 8x uint16,  -> uint8  mask
//
// Both pack the per-lane high (sign) bit into an integer mask. They are used
// only on the native NEON kswv mate-rescue kernels (kswv_neon_u8/16), where the
// mask drives work-reduction / early-exit branches.
//
// WHY THIS TEST EXISTS. PR #160 once rewrote both helpers as
// `vand(vshr_logical(v), {1<<lane})`. A *logical* shift yields 0x01 per set
// lane, and 0x01 & (1<<lane) is zero for every lane but 0 (and 8), so the mask
// silently collapsed to {0,1,0x100,0x101}. End-to-end SAM output stayed
// bit-identical (the broken mask only defeats an optimization, the slow path
// still computes the right score), so every output-level test passed while the
// kernels burned ~+22% instructions. This test checks the helpers *directly*
// against a scalar reference so that class of regression cannot pass silently
// again. The signed (arithmetic) shift is what makes the vand correct.
//
// Runs on the NEON build only; a no-op placeholder keeps the file valid on x86.

#include <array>
#include <cstdint>
#include <random>
#include <vector>

#include "doctest/doctest.h"

#if defined(__ARM_NEON) || defined(__aarch64__)

#include <arm_neon.h>

#include "neon_utils.h"   // neon_movemask_u8
#include "simd_compat.h"  // _mm_movemask_epi16, __m128i, _mm_loadu_si128

namespace {

// Reference: bit i of the result is the MSB of lane i. Only the high bit of
// each lane may contribute; all lower bits are ignored.
uint16_t ref_movemask_u8(const uint8_t b[16]) {
    uint16_t m = 0;
    for (int i = 0; i < 16; ++i)
        if (b[i] & 0x80u) m |= static_cast<uint16_t>(1u << i);
    return m;
}

int ref_movemask_epi16(const uint16_t w[8]) {
    int m = 0;
    for (int i = 0; i < 8; ++i)
        if (w[i] & 0x8000u) m |= (1 << i);
    return m;
}

}  // namespace

TEST_SUITE("unit/neon_movemask") {

TEST_CASE("neon_movemask_u8 matches scalar reference (edge + sweep + random)") {
    std::vector<std::array<uint8_t, 16>> cases = {
        {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},                  // none -> 0x0000
        {{0x80,0x80,0x80,0x80,0x80,0x80,0x80,0x80,
          0x80,0x80,0x80,0x80,0x80,0x80,0x80,0x80}},                          // all  -> 0xffff
        {{0x80,0,0x80,0,0,0,0,0, 0,0,0x80,0,0,0,0,0x80}},                     // lanes 0,2,10,15 -> 0x8405
        {{0,0x80,0,0,0,0,0,0, 0,0,0,0,0,0,0,0}},                              // lane 1 only -> 0x0002
        // only the MSB matters: 0x7f (MSB clear, low bits set) must NOT set a bit
        {{0x7f,0xff,0x7f,0xff,0x7f,0xff,0x7f,0xff,
          0x7f,0xff,0x7f,0xff,0x7f,0xff,0x7f,0xff}},                          // -> 0xaaaa
    };
    for (const auto& b : cases) {
        uint8x16_t v = vld1q_u8(b.data());
        CHECK(neon_movemask_u8(v) == ref_movemask_u8(b.data()));
    }

    // Single-lane sweep: exactly one high bit set, in turn, in each of 16 lanes.
    for (int lane = 0; lane < 16; ++lane) {
        std::array<uint8_t, 16> b{};
        b[lane] = 0x80;
        uint8x16_t v = vld1q_u8(b.data());
        CHECK(neon_movemask_u8(v) == static_cast<uint16_t>(1u << lane));
    }

    // Randomized fuzz.
    std::mt19937 rng(0xC0FFEEu);
    std::uniform_int_distribution<int> byte(0, 255);
    int n = 0;
    for (int t = 0; t < 20000; ++t) {
        std::array<uint8_t, 16> b{};
        for (auto& x : b) x = static_cast<uint8_t>(byte(rng));
        uint8x16_t v = vld1q_u8(b.data());
        CHECK(neon_movemask_u8(v) == ref_movemask_u8(b.data()));
        ++n;
    }
    MESSAGE("neon_movemask_u8: " << cases.size() << " edge + 16 sweep + " << n << " random vectors");
}

TEST_CASE("_mm_movemask_epi16 matches scalar reference (edge + sweep + random)") {
    std::vector<std::array<uint16_t, 8>> cases = {
        {{0, 0, 0, 0, 0, 0, 0, 0}},                                          // none -> 0x00
        {{0x8000,0x8000,0x8000,0x8000,0x8000,0x8000,0x8000,0x8000}},          // all  -> 0xff
        {{0x8000,0,0x8000,0,0,0,0,0x8000}},                                   // lanes 0,2,7 -> 0x85
        {{0,0x8000,0,0,0,0,0,0}},                                            // lane 1 only -> 0x02
        {{0x7fff,0xffff,0x7fff,0xffff,0x7fff,0xffff,0x7fff,0xffff}},          // MSB-only -> 0xaa
    };
    for (const auto& w : cases) {
        __m128i v = _mm_loadu_si128(reinterpret_cast<const __m128i*>(w.data()));
        CHECK(_mm_movemask_epi16(v) == ref_movemask_epi16(w.data()));
    }

    for (int lane = 0; lane < 8; ++lane) {
        std::array<uint16_t, 8> w{};
        w[lane] = 0x8000;
        __m128i v = _mm_loadu_si128(reinterpret_cast<const __m128i*>(w.data()));
        CHECK(_mm_movemask_epi16(v) == (1 << lane));
    }

    std::mt19937 rng(0xBADC0DEu);
    std::uniform_int_distribution<int> half(0, 0xFFFF);
    int n = 0;
    for (int t = 0; t < 20000; ++t) {
        std::array<uint16_t, 8> w{};
        for (auto& x : w) x = static_cast<uint16_t>(half(rng));
        __m128i v = _mm_loadu_si128(reinterpret_cast<const __m128i*>(w.data()));
        CHECK(_mm_movemask_epi16(v) == ref_movemask_epi16(w.data()));
        ++n;
    }
    MESSAGE("_mm_movemask_epi16: " << cases.size() << " edge + 8 sweep + " << n << " random vectors");
}

}  // TEST_SUITE

#else  // non-NEON build: helpers don't exist; keep the file valid.

TEST_CASE("neon_movemask parity [skipped: non-NEON build]") {
    CHECK(true);
}

#endif
