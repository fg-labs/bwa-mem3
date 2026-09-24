// test/unit/test_rescue_env.cpp
//
// rescue_env_int (src/rescue_env.h) is the shared reader for the integer
// BWA3_RESCUE_* tuning knobs. Its contract: a knob is a non-negative decimal
// integer; anything else -- a sign, leading or trailing characters, overflow --
// is reported and the default used, so a typo cannot silently set a knob.

#include <climits>
#include <cstdlib>

#include "doctest/doctest.h"
#include "rescue_env.h"

namespace {

constexpr const char *kKnob = "BWA3_RESCUE_ENV_TEST_KNOB";
constexpr int kDefault = 12345;  // distinct from every value parsed below

int read_as(const char *value)
{
    setenv(kKnob, value, 1);
    const int v = rescue_env_int(kKnob, kDefault);
    unsetenv(kKnob);
    return v;
}

} // namespace

TEST_CASE("rescue_env_int: unset or empty knobs take the default"
          * doctest::test_suite("unit/util")) {
    unsetenv(kKnob);
    CHECK(rescue_env_int(kKnob, kDefault) == kDefault);
    CHECK(read_as("") == kDefault);
}

TEST_CASE("rescue_env_int: non-negative decimal integers are read as given"
          * doctest::test_suite("unit/util")) {
    CHECK(read_as("0") == 0);
    CHECK(read_as("7") == 7);
    CHECK(read_as("100000000") == 100000000);
    CHECK(read_as("2147483647") == INT_MAX);
}

TEST_CASE("rescue_env_int: a sign, stray characters or overflow take the default"
          * doctest::test_suite("unit/util")) {
    SUBCASE("leading plus sign") { CHECK(read_as("+85") == kDefault); }
    SUBCASE("negative zero") { CHECK(read_as("-0") == kDefault); }
    SUBCASE("negative value") { CHECK(read_as("-5") == kDefault); }
    SUBCASE("leading space") { CHECK(read_as(" 85") == kDefault); }
    SUBCASE("leading tab") { CHECK(read_as("\t85") == kDefault); }
    SUBCASE("trailing space") { CHECK(read_as("85 ") == kDefault); }
    SUBCASE("trailing letter") { CHECK(read_as("85x") == kDefault); }
    SUBCASE("hex prefix") { CHECK(read_as("0x10") == kDefault); }
    SUBCASE("no digits") { CHECK(read_as("abc") == kDefault); }
    SUBCASE("one past INT_MAX") { CHECK(read_as("2147483648") == kDefault); }
    SUBCASE("overflows long") { CHECK(read_as("99999999999999999999") == kDefault); }
}
