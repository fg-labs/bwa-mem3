// Unit tests for mem_opt_parse_meth_tags (--meth-tags).
//
// Grammar: `all` | `none` | a comma-separated list of tag names that is either
// all plain (an inclusion set) or all `^`-prefixed (subtracted from the full
// set). Mixing the two forms is rejected rather than given an arbitrary
// meaning. Tag names are case-insensitive.

#include "doctest/doctest.h"
#include "bwamem.h"

#include <string>

namespace {

// Parse `spec` and require success; returns the resulting bitmask.
int parse_ok(const char *spec)
{
    int bits = -1;
    const char *err = nullptr;
    REQUIRE(mem_opt_parse_meth_tags(spec, &bits, &err) == 0);
    return bits;
}

// Parse `spec` and require failure; returns the reason string.
const char *parse_err(const char *spec)
{
    int bits = -1;
    const char *err = nullptr;
    REQUIRE(mem_opt_parse_meth_tags(spec, &bits, &err) == -1);
    REQUIRE(err != nullptr);
    return err;
}

}  // namespace

TEST_CASE("all and none are the two extremes" * doctest::description("--meth-tags"))
{
    CHECK(parse_ok("all") == MEM_METH_TAGS_ALL);
    CHECK(parse_ok("none") == 0);
    // The whole point of `none` is that it is distinguishable from a parse
    // failure: it must succeed and yield an empty mask, not error.
    CHECK(parse_ok("ALL") == MEM_METH_TAGS_ALL);
    CHECK(parse_ok("None") == 0);
}

TEST_CASE("inclusion lists select exactly the named tags")
{
    CHECK(parse_ok("XR") == MEM_METH_TAG_XR);
    CHECK(parse_ok("XG") == MEM_METH_TAG_XG);
    CHECK(parse_ok("XM") == MEM_METH_TAG_XM);
    CHECK(parse_ok("XR,XG") == (MEM_METH_TAG_XR | MEM_METH_TAG_XG));
    CHECK(parse_ok("XM,XR,XG") == MEM_METH_TAGS_ALL);
    // Order and repetition are immaterial -- it is a set.
    CHECK(parse_ok("XG,XR") == parse_ok("XR,XG"));
    CHECK(parse_ok("XR,XR") == MEM_METH_TAG_XR);
}

TEST_CASE("exclusions subtract from the full set")
{
    // The issue's motivating case: keep the cheap strand labels, drop the
    // read-length call string.
    CHECK(parse_ok("^XM") == (MEM_METH_TAG_XR | MEM_METH_TAG_XG));
    CHECK(parse_ok("^XR") == (MEM_METH_TAG_XG | MEM_METH_TAG_XM));
    CHECK(parse_ok("^XR,^XG") == MEM_METH_TAG_XM);
    CHECK(parse_ok("^XR,^XG,^XM") == 0);
    // An exclusion spec and the equivalent inclusion spec agree.
    CHECK(parse_ok("^XM") == parse_ok("XR,XG"));
}

TEST_CASE("tag names are case-insensitive")
{
    CHECK(parse_ok("xr,xg") == parse_ok("XR,XG"));
    CHECK(parse_ok("^xm") == parse_ok("^XM"));
    CHECK(parse_ok("Xm") == MEM_METH_TAG_XM);
}

TEST_CASE("malformed specs are rejected with a reason")
{
    // Mixing forms has no obvious reading -- is "XR,^XM" only-XR or
    // everything-but-XM? Reject rather than pick one silently.
    CHECK(parse_err("XR,^XM") != nullptr);
    CHECK(parse_err("^XM,XR") != nullptr);
    // Unknown or malformed tag names.
    CHECK(parse_err("XZ") != nullptr);
    CHECK(parse_err("NM") != nullptr);
    CHECK(parse_err("XRX") != nullptr);
    CHECK(parse_err("X") != nullptr);
    CHECK(parse_err("^") != nullptr);
    // Empty and trailing-comma specs.
    CHECK(parse_err("") != nullptr);
    CHECK(parse_err("XR,") != nullptr);
    CHECK(parse_err(",XR") != nullptr);
    CHECK(parse_err("XR,,XG") != nullptr);
}

TEST_CASE("an empty list element is not diagnosed as an unknown tag")
{
    // A missing element and a misspelled one are different mistakes and must
    // read differently: there is no tag name in ",XR" to be unknown, and
    // saying so sends the user hunting for a typo that isn't there.
    using std::string;
    const string misspelled = parse_err("XZ");
    CHECK(misspelled.find("unknown") != string::npos);
    for (const char *empty_elem : {",XR", "XR,,XG", "XR,"}) {
        const string reason = parse_err(empty_elem);
        CHECK(reason != misspelled);
        CHECK(reason.find("unknown") == string::npos);
    }
}

TEST_CASE("a NULL spec is rejected rather than dereferenced")
{
    int bits = -1;
    const char *err = nullptr;
    CHECK(mem_opt_parse_meth_tags(nullptr, &bits, &err) == -1);
    CHECK(err != nullptr);
}

TEST_CASE("the default opt emits the full Bismark set")
{
    mem_opt_t *o = mem_opt_init();
    CHECK(o->meth_tags == MEM_METH_TAGS_ALL);
    free(o);
}
