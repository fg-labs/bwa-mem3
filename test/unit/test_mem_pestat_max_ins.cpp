// test/unit/test_mem_pestat_max_ins.cpp — mem_pestat bounds opt->max_ins before
// it sizes the insert-size histogram.
//
// mem_pestat now counts insert sizes into a direct-address array rather than
// collecting and sorting them, so opt->max_ins went from being a pure FILTER on
// observed values -- where any value was harmless -- to being an ARRAY BOUND:
// the function allocates 4 * (max_ins + 1) * sizeof(int64_t) up front.
//
// mem_opt_init() sets max_ins to 10000 (320 kB) and no CLI option changes it,
// but mem_pestat reads the field from the caller's mem_opt_t and is a library
// entry point, so the value is not the aligner's to guarantee. Unbounded:
//
//   INT_MAX  -- `(size_t)4 * (max_ins + 1)` overflows int before the conversion
//               (undefined behaviour), and any surviving value asks for tens of
//               GB of zeroed counters.
//   negative -- `max_ins + 1` converts to an enormous size_t.
//
// Both are clamped instead: a negative value to 0, and anything past the cap to
// the cap. n == 0 here, so the candidate loop never runs and `regs` is never
// dereferenced -- the allocation is the whole point of the call, and the
// assertions below only need mem_pestat to return normally with every
// orientation marked failed for want of pairs.
//
// Without the clamp this test does not merely fail an assertion; it aborts on
// the xassert for a failed calloc, or takes the process down with it.
//
// Fixtures are built in-memory; no test data files are read.

#include <climits>
#include <cstdlib>

#include "doctest/doctest.h"
#include "bwamem.h"

namespace {

// Every orientation must report `failed` after a call with no pairs, whatever
// max_ins was: total[d] is 0, which is below MIN_DIR_CNT.
void expect_all_orientations_failed(const mem_pestat_t pes[4])
{
    for (int d = 0; d < 4; ++d) CHECK(pes[d].failed == 1);
}

// Run mem_pestat over an empty cohort with a given max_ins. Returns via `pes`.
void pestat_with_max_ins(int max_ins, mem_pestat_t pes[4])
{
    mem_opt_t *opt = mem_opt_init();
    REQUIRE(opt != NULL);
    opt->max_ins = max_ins;
    // l_pac is only read by mem_infer_dir inside the candidate loop, which n == 0
    // skips entirely; any value does.
    mem_pestat(opt, /*l_pac=*/1000000, /*n=*/0, /*regs=*/NULL, pes);
    free(opt);
}

} // namespace

TEST_CASE("mem_pestat: the default max_ins is accepted unchanged")
{
    mem_opt_t *opt = mem_opt_init();
    REQUIRE(opt != NULL);
    // Pins the premise the clamp is built around: the aligner's own value is
    // well inside the cap, so no production run takes the clamped path and this
    // guard cannot perturb output.
    CHECK(opt->max_ins == 10000);
    free(opt);

    mem_pestat_t pes[4];
    pestat_with_max_ins(10000, pes);
    expect_all_orientations_failed(pes);
}

TEST_CASE("mem_pestat: an oversized max_ins does not size the histogram")
{
    mem_pestat_t pes[4];
    pestat_with_max_ins(INT_MAX, pes);
    expect_all_orientations_failed(pes);
}

TEST_CASE("mem_pestat: a negative max_ins does not wrap to a huge allocation")
{
    mem_pestat_t pes[4];
    pestat_with_max_ins(-1, pes);
    expect_all_orientations_failed(pes);

    // -INT_MAX rather than INT_MIN: negating INT_MIN is itself UB, and the point
    // here is a large-magnitude negative, not that exact bit pattern.
    pestat_with_max_ins(-INT_MAX, pes);
    expect_all_orientations_failed(pes);
}
