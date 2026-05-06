// test/unit/test_simd_dispatch.cpp
//
// Unit tests for the SIMD dispatcher: init idempotency and tier-name table.

#include "doctest/doctest.h"
#include "../../src/simd_dispatch.h"
#include <stdlib.h>
#include <string.h>

TEST_CASE("simd: init is idempotent and yields a known tier") {
    bwamem3_simd_init();
    int t1 = bwamem3_simd_tier();
    bwamem3_simd_init();
    int t2 = bwamem3_simd_tier();
    CHECK(t1 == t2);
    const char *name = bwamem3_simd_tier_name(t1);
    CHECK(strcmp(name, "unknown") != 0);
}

TEST_CASE("simd: tier_name maps every known enum to a non-unknown string") {
    CHECK(strcmp(bwamem3_simd_tier_name(BWAMEM3_TIER_NONE),     "scalar")   == 0);
    CHECK(strcmp(bwamem3_simd_tier_name(BWAMEM3_TIER_SSE41),    "sse41")    == 0);
    CHECK(strcmp(bwamem3_simd_tier_name(BWAMEM3_TIER_SSE42),    "sse42")    == 0);
    CHECK(strcmp(bwamem3_simd_tier_name(BWAMEM3_TIER_AVX),      "avx")      == 0);
    CHECK(strcmp(bwamem3_simd_tier_name(BWAMEM3_TIER_AVX2),     "avx2")     == 0);
    CHECK(strcmp(bwamem3_simd_tier_name(BWAMEM3_TIER_AVX512BW), "avx512bw") == 0);
    CHECK(strcmp(bwamem3_simd_tier_name(BWAMEM3_TIER_NEON),     "neon")     == 0);
}
