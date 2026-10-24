// Unit tests for mem_skip_contained_ext_sound -- the proof envelope of the
// two-wave contained-seed extension skip (src/bwamem.cpp).
//
// The skip is byte-identical to the reference extension path exactly when
// cal_max_gap is non-decreasing in its distance argument, which holds for every
// non-negative match score. The extension driver reads this predicate, not the
// skip_contained_ext flag, so a scoring outside the envelope runs the reference
// path whatever the flag says.
//
// main_mem now refuses -A outside 1..127, so the CLI can no longer reach a
// scoring outside the envelope (the regression that used to run -A -1 and -A 0
// through compat_contained_skip_identity.sh is now a parse refusal in
// arg_range_validation.sh). The predicate stays as the driver's own guard, and
// this test keeps its boundary honest.

#include "doctest/doctest.h"
#include "bwamem.h"

#include <cstdlib>

namespace {

int sound_at(int a)
{
    mem_opt_t *o = mem_opt_init();
    o->a = a;
    const int r = mem_skip_contained_ext_sound(o);
    free(o);
    return r;
}

}  // namespace

TEST_CASE("contained-seed skip envelope: every match score the CLI accepts is inside it"
          * doctest::test_suite("unit/extension")) {
    CHECK(sound_at(1) == 1);     // the default
    CHECK(sound_at(2) == 1);
    CHECK(sound_at(127) == 1);   // the largest -A main_mem accepts
}

TEST_CASE("contained-seed skip envelope: the boundary is a >= 0"
          * doctest::test_suite("unit/extension")) {
    // a == 0 keeps cal_max_gap constant in d, so it is still non-decreasing.
    CHECK(sound_at(0) == 1);
    // A negative match score makes cal_max_gap decreasing: outside the envelope.
    CHECK(sound_at(-1) == 0);
    CHECK(sound_at(-128) == 0);
}

TEST_CASE("contained-seed skip envelope: it does not depend on the skip flag"
          * doctest::test_suite("unit/extension")) {
    mem_opt_t *o = mem_opt_init();
    o->a = -1;
    o->skip_contained_ext = 1;
    CHECK(mem_skip_contained_ext_sound(o) == 0);
    o->a = 1;
    o->skip_contained_ext = 0;
    CHECK(mem_skip_contained_ext_sound(o) == 1);
    free(o);
}
