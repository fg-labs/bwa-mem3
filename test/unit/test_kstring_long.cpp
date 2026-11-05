// Long integer appenders must match decimal formatting without touching spare bytes.
#include "doctest/doctest.h"
#include "kstring.h"

#include <climits>
#include <cstdio>
#include <cstring>
#include <random>

namespace {

void check_unsigned_long(unsigned long value) {
    char expected[32];
    const int n = snprintf(expected, sizeof(expected), "%lu", value);
    REQUIRE(n > 0);
    char buffer[64];
    memset(buffer, '#', sizeof(buffer));
    buffer[0] = 'x';
    kstring_t out = {1, sizeof(buffer), buffer};
    kputul_u(value, &out);
    CHECK(out.l == (size_t)n + 1);
    CHECK(buffer[0] == 'x');
    CHECK(memcmp(buffer + 1, expected, (size_t)n) == 0);
    // Unchecked appenders must not write a trailing NUL or touch later bytes.
    for (size_t i = out.l; i < sizeof(buffer); ++i) CHECK(buffer[i] == '#');
}

void check_signed_long(long value) {
    char expected[32];
    const int n = snprintf(expected, sizeof(expected), "%ld", value);
    REQUIRE(n > 0);
    char buffer[64];
    memset(buffer, '#', sizeof(buffer));
    buffer[0] = 'x';
    kstring_t out = {1, sizeof(buffer), buffer};
    kputl_u(value, &out);
    CHECK(out.l == (size_t)n + 1);
    CHECK(buffer[0] == 'x');
    CHECK(memcmp(buffer + 1, expected, (size_t)n) == 0);
    for (size_t i = out.l; i < sizeof(buffer); ++i) CHECK(buffer[i] == '#');
}

} // namespace

TEST_CASE("kstring long appenders: decimal boundaries and full unsigned range"
          * doctest::test_suite("unit/kstring_long")) {
    for (unsigned long value : {0ul, 1ul, 9ul, 10ul, 11ul, 99ul, 100ul, 101ul,
                                (unsigned long)UINT_MAX, ULONG_MAX})
        check_unsigned_long(value);
    for (unsigned long power = 10; ; power *= 10) {
        check_unsigned_long(power - 1);
        check_unsigned_long(power);
        check_unsigned_long(power + 1);
        if (power > ULONG_MAX / 10) break;
    }
    std::mt19937_64 rng(91827);
    for (int i = 0; i < 2000; ++i) check_unsigned_long((unsigned long)rng());
}

TEST_CASE("kstring signed long appender: sign, LONG_MIN and negative template lengths"
          * doctest::test_suite("unit/kstring_long")) {
    for (long value : {LONG_MIN, LONG_MIN + 1, LONG_MAX, -1000000000l,
                       -10000l, -100l, -99l, -10l, -9l, -1l, 0l, 1l, 9l, 10l,
                       99l, 100l, 101l})
        check_signed_long(value);
#if LONG_MAX > 2147483647L
    check_signed_long(-4294967295l);
    check_signed_long(4294967295l);
#endif
    std::mt19937_64 rng(8127);
    for (int i = 0; i < 2000; ++i) {
        const long value = (long)(rng() & (unsigned long)LONG_MAX);
        check_signed_long(value);
        check_signed_long(-value);
    }
}
