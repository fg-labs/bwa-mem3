#include "doctest/doctest.h"
#include "bwamem.h"

TEST_CASE("XA cutoff keeps exact 0.8 boundary" * doctest::test_suite("unit/xa")) {
    mem_opt_t *opt = mem_opt_init();
    REQUIRE(opt != nullptr);
    // get_pri_idx compares integer scores after promoting the stored ratio.
    // A float-backed 0.8 rounds upward and wrongly drops the score-32 hit.
    CHECK(32 >= 40 * static_cast<double>(opt->XA_drop_ratio));
    CHECK_FALSE(31 >= 40 * static_cast<double>(opt->XA_drop_ratio));

    // Preserve caller-supplied precision on either side of the same boundary.
    opt->XA_drop_ratio = 0.7999999999999999;
    CHECK(32 >= 40 * static_cast<double>(opt->XA_drop_ratio));
    opt->XA_drop_ratio = 0.8000000000000002;
    CHECK_FALSE(32 >= 40 * static_cast<double>(opt->XA_drop_ratio));
    free(opt);
}
