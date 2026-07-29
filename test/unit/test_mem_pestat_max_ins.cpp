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
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <unistd.h>

#include "doctest/doctest.h"
#include "bwamem.h"

namespace {

/* The cap mem_pestat clamps to, mirrored from bwamem_pair.cpp. Duplicated on
 * purpose: the boundary test below is only meaningful if it names the value
 * independently of the code under test. */
const int kMaxInsCap = 1000000;

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

/* As above, but returns whatever mem_pestat wrote to stderr.
 *
 * Needed for the boundary case: with n == 0 every orientation ends up failed
 * whether or not the clamp fired, so the return value cannot distinguish "at the
 * cap, passed through" from "past the cap, clamped". The clamp's diagnostic is
 * the only externally visible difference, which makes it the thing to assert on. */
std::string pestat_stderr_with_max_ins(int max_ins, mem_pestat_t pes[4])
{
    fflush(stderr);
    const int saved_stderr = dup(fileno(stderr));
    REQUIRE(saved_stderr >= 0);

    /* Honour TMPDIR, as test/unit/test_fmi_pread_from_stream.cpp does: CI
     * runners point it at a writable scratch area, and /tmp is RAM-backed
     * tmpfs on several modern distributions. */
    const char *tmpdir = getenv("TMPDIR");
    if (tmpdir == NULL || tmpdir[0] == '\0') tmpdir = "/tmp";
    std::string tmpl_s = tmpdir;
    if (tmpl_s[tmpl_s.size() - 1] != '/') tmpl_s += '/';
    tmpl_s += "pestat_stderr_XXXXXX";
    std::vector<char> tmpl(tmpl_s.begin(), tmpl_s.end());
    tmpl.push_back('\0');
    const int fd = mkstemp(tmpl.data());
    REQUIRE(fd >= 0);
    REQUIRE(dup2(fd, fileno(stderr)) >= 0);

    pestat_with_max_ins(max_ins, pes);

    fflush(stderr);
    REQUIRE(dup2(saved_stderr, fileno(stderr)) >= 0);
    close(saved_stderr);

    std::string captured;
    REQUIRE(lseek(fd, 0, SEEK_SET) == 0);
    char buf[4096];
    ssize_t got;
    while ((got = read(fd, buf, sizeof(buf))) > 0) captured.append(buf, (size_t)got);
    close(fd);
    unlink(tmpl.data());
    return captured;
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

TEST_CASE("mem_pestat: max_ins at the cap passes through, one past it is clamped")
{
    // The clamp is `max_ins < 0 || max_ins > kMaxInsCap`, so the cap itself is
    // ACCEPTED and cap+1 is the first clamped value. An off-by-one either way --
    // `>=` here, or a cap of 999999 -- moves that edge, and only these two inputs
    // sit close enough to notice.
    mem_pestat_t pes[4];

    const std::string at_cap = pestat_stderr_with_max_ins(kMaxInsCap, pes);
    expect_all_orientations_failed(pes);
    CHECK(at_cap.find("out of range") == std::string::npos);

    const std::string past_cap = pestat_stderr_with_max_ins(kMaxInsCap + 1, pes);
    expect_all_orientations_failed(pes);
    CHECK(past_cap.find("out of range") != std::string::npos);
    // Names the value it fell back to, not just that it complained.
    CHECK(past_cap.find("using 1000000") != std::string::npos);
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
