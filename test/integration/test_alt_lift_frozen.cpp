// test/integration/test_alt_lift_frozen.cpp — the exclude-BED table gate on a real index.
//
// The same gate test/unit/test_alt_lift.cpp runs on generated indexes (restricted
// loading of a .alt with --alt-project-exclude-bed equals loading a FLAG-masked
// copy of it, field by field), on a caller-supplied index, .alt pair and BED.
// It reads those files, so it lives here rather than in the unit suite.

#include "../unit/alt_lift_gate.h"

// Release gate on a real index, configured only by environment: unset, the case
// is skipped; set, all four inputs are required. Only when no assertion in the
// case has failed does it print one ALT_GATE_FROZEN_COMPLETE line with the counts
// it checked, so a release run can tell a completed gate from a skipped or
// failed one.
TEST_CASE("exclude BED frozen table gate" * doctest::test_suite("integration/alt_lift_frozen") *
          doctest::skip(!alt_gate_any_env_set())) {
    const long failed_before = g_alt_gate_failed_asserts.load();
    std::string value[4];
    for (int i = 0; i < 4; ++i) {
        const char *v = getenv(kAltGateEnv[i]);
        INFO("partial configuration: " << std::string(kAltGateEnv[i]) << " is unset");
        REQUIRE((v != nullptr && *v != '\0'));
        value[i] = v;
    }
    const AltGateCounts c = run_alt_table_gate({value[0], value[1], value[2], value[3]});
    REQUIRE(c.rows > 0);
    REQUIRE(c.excluded > 0);
    REQUIRE(c.retained > 0);
    REQUIRE(c.input_intervals > 0);
    alt_gate_emit_marker(stdout, failed_before, c);
}
