// test/unit/test_simd_dispatch.cpp
//
// Unit tests for the SIMD dispatcher: init idempotency, the tier-name table,
// the host-floor check, and the x86 floor and kernel tier rankings.

#include "doctest/doctest.h"
#include "../../src/simd_dispatch.h"
#include <stdlib.h>
#include <string.h>

TEST_CASE("simd: init is idempotent and yields a known tier" * doctest::test_suite("unit/simd_dispatch")) {
    bwamem3_simd_init();
    int t1 = bwamem3_simd_tier();
    bwamem3_simd_init();
    int t2 = bwamem3_simd_tier();
    CHECK(t1 == t2);
    const char *name = bwamem3_simd_tier_name(t1);
    CHECK(strcmp(name, "unknown") != 0);
}

TEST_CASE("simd: tier_name maps every known enum to a non-unknown string" * doctest::test_suite("unit/simd_dispatch")) {
    CHECK(strcmp(bwamem3_simd_tier_name(BWAMEM3_TIER_NONE),     "scalar")   == 0);
    CHECK(strcmp(bwamem3_simd_tier_name(BWAMEM3_TIER_SSE41),    "sse41")    == 0);
    CHECK(strcmp(bwamem3_simd_tier_name(BWAMEM3_TIER_SSE42),    "sse42")    == 0);
    CHECK(strcmp(bwamem3_simd_tier_name(BWAMEM3_TIER_AVX),      "avx")      == 0);
    CHECK(strcmp(bwamem3_simd_tier_name(BWAMEM3_TIER_AVX2),     "avx2")     == 0);
    CHECK(strcmp(bwamem3_simd_tier_name(BWAMEM3_TIER_AVX512BW), "avx512bw") == 0);
    CHECK(strcmp(bwamem3_simd_tier_name(BWAMEM3_TIER_NEON),     "neon")     == 0);
}

TEST_CASE("check_host_floor: host above floor returns 1" * doctest::test_suite("unit/simd_dispatch")) {
    CHECK(bwamem3_check_host_floor(BWAMEM3_TIER_AVX512BW, BWAMEM3_TIER_AVX2) == 1);
    CHECK(bwamem3_check_host_floor(BWAMEM3_TIER_AVX2,     BWAMEM3_TIER_SSE41) == 1);
}

TEST_CASE("check_host_floor: host equal to floor returns 1" * doctest::test_suite("unit/simd_dispatch")) {
    CHECK(bwamem3_check_host_floor(BWAMEM3_TIER_AVX2,     BWAMEM3_TIER_AVX2) == 1);
    CHECK(bwamem3_check_host_floor(BWAMEM3_TIER_SSE41,    BWAMEM3_TIER_SSE41) == 1);
}

TEST_CASE("check_host_floor: host below floor returns 0" * doctest::test_suite("unit/simd_dispatch")) {
    CHECK(bwamem3_check_host_floor(BWAMEM3_TIER_SSE41,    BWAMEM3_TIER_AVX2) == 0);
    CHECK(bwamem3_check_host_floor(BWAMEM3_TIER_AVX,      BWAMEM3_TIER_AVX2) == 0);
}

TEST_CASE("check_host_floor: NEON-on-NEON is always 1" * doctest::test_suite("unit/simd_dispatch")) {
    CHECK(bwamem3_check_host_floor(BWAMEM3_TIER_NEON,     BWAMEM3_TIER_NEON) == 1);
}

TEST_CASE("check_host_floor: x86 vs NEON is rejected (orthogonal)" * doctest::test_suite("unit/simd_dispatch")) {
    CHECK(bwamem3_check_host_floor(BWAMEM3_TIER_AVX2,     BWAMEM3_TIER_NEON) == 0);
    CHECK(bwamem3_check_host_floor(BWAMEM3_TIER_NEON,     BWAMEM3_TIER_AVX2) == 0);
}

TEST_CASE("check_host_floor: TIER_NONE (scalar) boundary" * doctest::test_suite("unit/simd_dispatch")) {
    /* Scalar build accepts any host (build_tier == 0 is the lowest). */
    CHECK(bwamem3_check_host_floor(BWAMEM3_TIER_NONE,     BWAMEM3_TIER_NONE)  == 1);
    CHECK(bwamem3_check_host_floor(BWAMEM3_TIER_SSE41,    BWAMEM3_TIER_NONE)  == 1);
    /* Scalar host can't run SSE4.1 build. */
    CHECK(bwamem3_check_host_floor(BWAMEM3_TIER_NONE,     BWAMEM3_TIER_SSE41) == 0);
}

TEST_CASE("format_host_floor_error: message contains both tier names" * doctest::test_suite("unit/simd_dispatch")) {
    char buf[1024];
    int n = bwamem3_format_host_floor_error(buf, sizeof(buf),
                                            BWAMEM3_TIER_SSE41,
                                            BWAMEM3_TIER_AVX2);
    REQUIRE(n > 0);
    CHECK(strstr(buf, "avx2")  != NULL);  /* build floor */
    CHECK(strstr(buf, "sse41") != NULL);  /* detected host */
    CHECK(strstr(buf, "SIMD")  != NULL);  /* identifies the failure class */
}

TEST_CASE("format_host_floor_error: message mentions BASELINE_ARCH remediation" * doctest::test_suite("unit/simd_dispatch")) {
    char buf[1024];
    int n = bwamem3_format_host_floor_error(buf, sizeof(buf),
                                            BWAMEM3_TIER_SSE41,
                                            BWAMEM3_TIER_AVX2);
    REQUIRE(n > 0);
    CHECK(strstr(buf, "BASELINE_ARCH") != NULL);
}

TEST_CASE("format_host_floor_error: small buffer returns -1, leaves buf NUL-terminated" * doctest::test_suite("unit/simd_dispatch")) {
    char buf[16];
    int n = bwamem3_format_host_floor_error(buf, sizeof(buf),
                                            BWAMEM3_TIER_SSE41,
                                            BWAMEM3_TIER_AVX2);
    CHECK(n == -1);
    /* Caller may not depend on contents, but length is bounded — buffer
     * must still be valid C string for safe further use. */
    CHECK(buf[sizeof(buf) - 1] == '\0');
}

TEST_CASE("format_host_floor_error: bufsz == 0 returns -1 without dereferencing buf" * doctest::test_suite("unit/simd_dispatch")) {
    char dummy = 'X';
    int n = bwamem3_format_host_floor_error(&dummy, 0,
                                            BWAMEM3_TIER_SSE41,
                                            BWAMEM3_TIER_AVX2);
    CHECK(n == -1);
    /* buf was not touched (the early return must precede any write). */
    CHECK(dummy == 'X');
}

TEST_CASE("format_host_floor_error: bufsz == 1 returns -1 and writes only NUL" * doctest::test_suite("unit/simd_dispatch")) {
    char buf[1] = { 'X' };
    int n = bwamem3_format_host_floor_error(buf, sizeof(buf),
                                            BWAMEM3_TIER_SSE41,
                                            BWAMEM3_TIER_AVX2);
    CHECK(n == -1);
    CHECK(buf[0] == '\0');
}

TEST_CASE("host_meets_floor: returns 1 on the test runner (host meets build)" * doctest::test_suite("unit/simd_dispatch")) {
    /* The test runner host meets its own build's floor by construction —
     * otherwise the test binary itself would have SIGILL'd before we got
     * here. Just verifies the wrapper calls through correctly. */
    CHECK(bwamem3_host_meets_floor() == 1);
}

TEST_CASE("enforce_host_floor: returns without exiting on the test runner" * doctest::test_suite("unit/simd_dispatch")) {
    /* Same reasoning as host_meets_floor's smoke test: if the runner met
     * the floor (it must have, to get here), enforce returns without
     * exiting. The exit-2 path on too-old hosts is verified by the
     * integration test in test/regression/host_floor_enforce.sh, which
     * uses BWAMEM3_TESTING_HOST_TIER injection. */
    bwamem3_enforce_host_floor();
    /* Reaching this line means it didn't exit. */
    CHECK(true);
}

/* The avx512bw kernels are built with AVX-512 VL and DQ as well as F and BW, but an
 * arch=avx512bw build's non-kernel code uses F and BW only. So a host with BW but not
 * VL or DQ meets that build's floor and runs the avx2 kernels. Arguments: avx512bw,
 * avx512vl, avx512dq, avx2, avx, sse4.2, sse4.1. */
TEST_CASE("x86 tiers: AVX-512BW alone meets the avx512bw floor but gets the avx2 kernels" * doctest::test_suite("unit/simd_dispatch")) {
    CHECK(bwamem3_x86_floor_tier(1, 0, 0, 1, 1, 1, 1)  == BWAMEM3_TIER_AVX512BW);
    CHECK(bwamem3_x86_kernel_tier(1, 0, 0, 1, 1, 1, 1) == BWAMEM3_TIER_AVX2);
    CHECK(bwamem3_check_host_floor(bwamem3_x86_floor_tier(1, 0, 0, 1, 1, 1, 1),
                                   BWAMEM3_TIER_AVX512BW) == 1);
}

TEST_CASE("x86 tiers: AVX-512BW with only one of VL and DQ still gets the avx2 kernels" * doctest::test_suite("unit/simd_dispatch")) {
    CHECK(bwamem3_x86_floor_tier(1, 1, 0, 1, 1, 1, 1)  == BWAMEM3_TIER_AVX512BW);
    CHECK(bwamem3_x86_kernel_tier(1, 1, 0, 1, 1, 1, 1) == BWAMEM3_TIER_AVX2);
    CHECK(bwamem3_x86_floor_tier(1, 0, 1, 1, 1, 1, 1)  == BWAMEM3_TIER_AVX512BW);
    CHECK(bwamem3_x86_kernel_tier(1, 0, 1, 1, 1, 1, 1) == BWAMEM3_TIER_AVX2);
}

TEST_CASE("x86 tiers: AVX-512BW, VL and DQ select the avx512bw kernels" * doctest::test_suite("unit/simd_dispatch")) {
    CHECK(bwamem3_x86_floor_tier(1, 1, 1, 1, 1, 1, 1)  == BWAMEM3_TIER_AVX512BW);
    CHECK(bwamem3_x86_kernel_tier(1, 1, 1, 1, 1, 1, 1) == BWAMEM3_TIER_AVX512BW);
    CHECK(bwamem3_check_host_floor(bwamem3_x86_floor_tier(1, 1, 1, 1, 1, 1, 1),
                                   BWAMEM3_TIER_AVX512BW) == 1);
}

TEST_CASE("x86 tiers: VL and DQ without AVX-512BW count for nothing" * doctest::test_suite("unit/simd_dispatch")) {
    CHECK(bwamem3_x86_floor_tier(0, 1, 1, 1, 1, 1, 1)  == BWAMEM3_TIER_AVX2);
    CHECK(bwamem3_x86_kernel_tier(0, 1, 1, 1, 1, 1, 1) == BWAMEM3_TIER_AVX2);
    CHECK(bwamem3_check_host_floor(bwamem3_x86_floor_tier(0, 1, 1, 1, 1, 1, 1),
                                   BWAMEM3_TIER_AVX512BW) == 0);
}

TEST_CASE("x86 tiers: below avx512bw the floor and kernel tiers agree" * doctest::test_suite("unit/simd_dispatch")) {
    CHECK(bwamem3_x86_floor_tier(0, 0, 0, 0, 1, 1, 1)  == BWAMEM3_TIER_AVX);
    CHECK(bwamem3_x86_kernel_tier(0, 0, 0, 0, 1, 1, 1) == BWAMEM3_TIER_AVX);
    CHECK(bwamem3_x86_floor_tier(0, 0, 0, 0, 0, 1, 1)  == BWAMEM3_TIER_SSE42);
    CHECK(bwamem3_x86_kernel_tier(0, 0, 0, 0, 0, 1, 1) == BWAMEM3_TIER_SSE42);
    CHECK(bwamem3_x86_floor_tier(0, 0, 0, 0, 0, 0, 1)  == BWAMEM3_TIER_SSE41);
    CHECK(bwamem3_x86_kernel_tier(0, 0, 0, 0, 0, 0, 1) == BWAMEM3_TIER_SSE41);
    CHECK(bwamem3_x86_floor_tier(0, 0, 0, 0, 0, 0, 0)  == BWAMEM3_TIER_NONE);
    CHECK(bwamem3_x86_kernel_tier(0, 0, 0, 0, 0, 0, 0) == BWAMEM3_TIER_NONE);
}

/* The `SIMD floor:` parenthetical names the oldest CPU each floor runs on. It is not a
 * psABI level: the floor check verifies only the tier's own subsets (an avx512bw build's
 * non-kernel code needs AVX-512 F and BW), not a full x86-64-v2 / -v3 / -v4 set. */
TEST_CASE("simd_floor_label: names a CPU family per x86 floor, not a psABI level" * doctest::test_suite("unit/simd_dispatch")) {
    CHECK(strcmp(bwamem3_simd_floor_label(BWAMEM3_TIER_AVX512BW), "Skylake-X 2017+")    == 0);
    CHECK(strcmp(bwamem3_simd_floor_label(BWAMEM3_TIER_AVX2),     "Haswell 2013+")      == 0);
    CHECK(strcmp(bwamem3_simd_floor_label(BWAMEM3_TIER_AVX),      "Sandy Bridge 2011+") == 0);
    CHECK(strcmp(bwamem3_simd_floor_label(BWAMEM3_TIER_SSE42),    "Nehalem 2008+")      == 0);
    CHECK(strcmp(bwamem3_simd_floor_label(BWAMEM3_TIER_SSE41),    "Penryn 2007+")       == 0);
    CHECK(strcmp(bwamem3_simd_floor_label(BWAMEM3_TIER_NONE),     "scalar")             == 0);
}

TEST_CASE("print_version_simd: the floor line never names a psABI level" * doctest::test_suite("unit/simd_dispatch")) {
    FILE *f = tmpfile();
    REQUIRE(f != NULL);
    bwamem3_print_version_simd(f);
    rewind(f);
    char buf[512];
    const size_t n = fread(buf, 1, sizeof(buf) - 1, f);
    fclose(f);
    buf[n] = '\0';
    CHECK(strstr(buf, "SIMD floor: ") != NULL);
    CHECK(strstr(buf, "x86-64-v") == NULL);
}
