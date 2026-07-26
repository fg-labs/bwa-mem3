// test/unit/test_index_memory_estimate.cpp — the index memory preflight arithmetic.
//
// `bwa-mem3 index` refuses a build whose estimated peak exceeds --max-memory.
// Two properties decide whether that refusal is correct, and both used to be
// inline arithmetic no test could reach:
//
//   libsais_sa_is_64bit(N)        — which SA width libsais will pick, since the
//                                   int64 path doubles the per-suffix cost.
//   libsais_estimate_peak_bytes(N)— the estimate compared against the budget.
//
// The regression these pin: a `min(50% of RAM, 32 GiB)` default meant hg38's
// estimate (12 B/base over a doubled 6.43 Gbase text = 71.9 GiB) lost to a
// 32 GiB budget on every host, including a 256 GiB server. The thresholds below
// are derived from the constants themselves — INT32_MAX - 10000 for the SA-width
// switch, 6/12 B/base for the cost — not recorded from a run.
//
// `N` throughout is the DOUBLED text length (2 * l_pac), which is what libsais
// is handed; a `--meth` seed index doubles it again (4 * l_pac).

#include "doctest/doctest.h"

#include "libsais_build.h"

#include <cstdint>

namespace {

constexpr int64_t kGiB          = 1LL << 30;
constexpr int64_t kSwitchPoint  = (int64_t)INT32_MAX - 10000;   // N + 1 >= this => int64 SA

// Reference lengths (l_pac) of the genomes this preflight actually gates.
constexpr int64_t kHg38Lpac   = 3217346917LL;
constexpr int64_t kChr17Lpac   = 83257441LL;

} // namespace

TEST_CASE("libsais_sa_is_64bit: switches exactly at INT32_MAX - 10000") {
    CHECK(libsais_sa_is_64bit(kSwitchPoint - 2) == false);
    CHECK(libsais_sa_is_64bit(kSwitchPoint - 1) == true);   // N + 1 == threshold
    CHECK(libsais_sa_is_64bit(kSwitchPoint) == true);
    CHECK(libsais_sa_is_64bit(0) == false);
}

TEST_CASE("libsais_estimate_peak_bytes: 6 B/base int32, 12 B/base int64") {
    CHECK(libsais_estimate_peak_bytes(1000) == 1001 * 6);
    CHECK(libsais_estimate_peak_bytes(kSwitchPoint) == (kSwitchPoint + 1) * 12);
    // Crossing the width switch makes the estimate jump, not grow smoothly.
    const int64_t below = libsais_estimate_peak_bytes(kSwitchPoint - 2);
    const int64_t above = libsais_estimate_peak_bytes(kSwitchPoint - 1);
    CHECK(above > below);
    CHECK(above >= below * 2);
}

TEST_CASE("libsais_estimate_peak_bytes: hg38 needs ~72 GiB, well past the old 32 GiB cap") {
    const int64_t hg38 = libsais_estimate_peak_bytes(2 * kHg38Lpac);
    CHECK(libsais_sa_is_64bit(2 * kHg38Lpac) == true);
    // ~71.9 GiB. Bracketed rather than pinned so a deliberate recalibration of
    // the per-base margin does not fail this test spuriously.
    CHECK(hg38 > 71 * kGiB);
    CHECK(hg38 < 73 * kGiB);
    // The defect: the removed default could never clear it, on any host.
    CHECK(hg38 > 32 * kGiB);
}

TEST_CASE("libsais_estimate_peak_bytes: a --meth seed costs ~2x its original") {
    // --meth builds the original over 2 * l_pac, then a seed over 4 * l_pac
    // (an f... C->T and an r... G->A contig per chromosome). Both hg38 halves
    // take the int64 path, so the ratio is exactly the text-length ratio.
    const int64_t orig = libsais_estimate_peak_bytes(2 * kHg38Lpac);
    const int64_t seed = libsais_estimate_peak_bytes(4 * kHg38Lpac);
    CHECK(libsais_sa_is_64bit(4 * kHg38Lpac) == true);
    CHECK(seed > orig);
    CHECK(seed == doctest::Approx((double)orig * 2).epsilon(0.001));
    CHECK(seed > 143 * kGiB);
    CHECK(seed < 145 * kGiB);

    // chr17 stays on the int32 path for BOTH halves, so the 2x holds there too
    // — this is the fixture the shell test exercises.
    const int64_t c_orig = libsais_estimate_peak_bytes(2 * kChr17Lpac);
    const int64_t c_seed = libsais_estimate_peak_bytes(4 * kChr17Lpac);
    CHECK(libsais_sa_is_64bit(4 * kChr17Lpac) == false);
    CHECK(c_seed == doctest::Approx((double)c_orig * 2).epsilon(0.001));

    // And the gap the up-front --meth preflight closes: a budget can clear the
    // original while the seed cannot fit, which used to be discovered only
    // after the original index had already been built and written.
    const int64_t budget = (c_orig + c_seed) / 2;
    CHECK(c_orig <= budget);
    CHECK(c_seed > budget);
}
