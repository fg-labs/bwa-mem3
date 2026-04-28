// test/integration/placeholder.cpp
//
// Placeholder integration test that ships while the 5 legacy binaries
// (fmi_test, smem2_test, bwt_seed_strategy_test, sa2ref_test, xeonbsw)
// get migrated one-by-one in follow-up PRs. This case always passes — it
// exists only so the integration binary links and CI has something to
// invoke.
//
// Delete this file in the PR that migrates the final legacy binary.

#include "doctest/doctest.h"

TEST_CASE("integration binary is wired up"
          * doctest::test_suite("integration")) {
    CHECK(true);
}
