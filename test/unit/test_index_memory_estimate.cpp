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
// switch, 8/12 B/base for the cost — plus, in the last two cases, peak-RSS
// figures recorded from real builds.
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

TEST_CASE("libsais_estimate_peak_bytes: 8 B/base int32, 12 B/base int64") {
    CHECK(libsais_estimate_peak_bytes(1000) == 1001 * 8);
    CHECK(libsais_estimate_peak_bytes(kSwitchPoint) == (kSwitchPoint + 1) * 12);
    // Crossing the width switch makes the estimate jump, not grow smoothly:
    // 8 -> 12 B/base, so a one-base increase costs 1.5x.
    const int64_t below = libsais_estimate_peak_bytes(kSwitchPoint - 2);
    const int64_t above = libsais_estimate_peak_bytes(kSwitchPoint - 1);
    CHECK(above > below);
    CHECK(above == doctest::Approx((double)below * 1.5).epsilon(0.001));
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

// The estimate is a promise: the build it gates must not exceed it. These are
// peak-RSS figures MEASURED with /usr/bin/time -l and RSS sampling on one
// 64 GiB 12-core arm64 host, so they are a floor on the required margin rather
// than a universal bound -- another host's allocator or thread count can peak
// higher. They exist to stop the per-base constants being lowered back to a
// value that real builds overrun, which is exactly what the int32 constant did
// at 6 B/base.
namespace {

struct MeasuredPeak {
    const char* label;
    int64_t     N;            // doubled text length handed to libsais
    int64_t     peak_bytes;   // observed peak RSS
};

// Single builds: `index` on cumulative hg38 slices (and each --meth run's
// first build, which is an ordinary build of the original reference).
const MeasuredPeak kSingleBuildPeaks[] = {
    {"chr17       (int32)",   166514882,   902348800LL},
    {"chr17 seed  (int32)",   333029764,  1776664576LL},
    {"chr1        (int32)",   497912844,  2812051456LL},
    {"chr1-2      (int32)",   982299902,  5185159168LL},
    {"chr1-3      (int32)",  1378891020,  8129052672LL},
    {"chr1-4      (int32)",  1759320130, 10151133184LL},
    {"chr1-5      (int32)",  2122396648, 12057460736LL},
    {"chr1-6      (int64)",  2464008606, 24248434688LL},
    {"chr1-7      (int64)",  2782700552, 27249328128LL},
    {"chr1-8      (int64)",  3072977824, 29934387200LL},
    {"chr1-9      (int64)",  3349767258, 30999101440LL},
};

// `index --meth` whole-process peaks, keyed by the SEED build's N (= 4*l_pac),
// which is the dominant build and therefore what the preflight checks.
const MeasuredPeak kMethProcessPeaks[] = {
    {"chr17 --meth  (int32)",  333029764,  2311536640LL},
    {"chr1 --meth   (int32)",  995825688,  7044399104LL},
    {"chr1-2 --meth (int32)", 1964599804, 11325145088LL},
    {"chr1-3 --meth (int64)", 2757782040, 26843906048LL},
    {"chr1-4 --meth (int64)", 3518640260, 35315105792LL},
    {"chr1-5 --meth (int64)", 4244793296, 38338084864LL},
};

} // namespace

TEST_CASE("libsais_estimate_peak_bytes: covers every measured single-build peak") {
    for (const MeasuredPeak& m : kSingleBuildPeaks) {
        const int64_t est = libsais_estimate_peak_bytes(m.N);
        CAPTURE(m.label);
        CAPTURE(est);
        CAPTURE(m.peak_bytes);
        CHECK(est > m.peak_bytes);
    }
}

TEST_CASE("libsais_estimate_peak_bytes: covers every measured --meth process peak") {
    // This is the property the old int32 constant broke: the seed build's
    // estimate has to cover the whole invocation, because the original build's
    // memory is still resident when the seed build runs.
    for (const MeasuredPeak& m : kMethProcessPeaks) {
        const int64_t est = libsais_estimate_peak_bytes(m.N);
        CAPTURE(m.label);
        CAPTURE(est);
        CAPTURE(m.peak_bytes);
        CHECK(est > m.peak_bytes);
    }
}
