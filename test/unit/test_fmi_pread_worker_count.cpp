// test/unit/test_fmi_pread_worker_count.cpp
//
// Regression test for the parallel index load's chunk-size floor.
//
// FMI_search::load_index splits each big index array (cp_occ, the compressed SA
// samples) across pread workers. The split is only worth its create/join cost
// while every worker still gets a bandwidth-sized slice, so the worker count is
// capped such that no chunk falls below FMI_PREAD_MIN_CHUNK.
//
// The original cap computed `nbytes / min_chunk + 1`, which hands out one
// worker more than the array can cover at that size: an array exactly one
// min_chunk long was split across two workers, giving each half the floor.
// Reachable whenever a caller requests more workers than the array can feed --
// a test-sized reference, or a large BWA3_LOAD_THREADS override -- since on a
// GB-scale index the load's own 8-worker cap binds long before this one.

#include "doctest/doctest.h"
#include "../../src/FMI_search.h"

#include <cstddef>

namespace {

// Largest chunk any worker receives when `nbytes` is split `nthreads` ways the
// way parallel_pread does it (base size plus one byte for the first `rem`).
size_t max_chunk_bytes(size_t nbytes, int nthreads)
{
    const size_t n = (size_t)nthreads;
    return nbytes / n + (nbytes % n != 0 ? 1 : 0);
}

// Smallest chunk any worker receives under the same split.
size_t min_chunk_bytes(size_t nbytes, int nthreads)
{
    return nbytes / (size_t)nthreads;
}

}  // namespace

TEST_CASE("fmi_pread_worker_count: never splits below the chunk floor") {
    const size_t floor_bytes = FMI_PREAD_MIN_CHUNK;

    // The regression: an array exactly at the floor must stay single-worker,
    // however many workers the caller asks for. Pre-fix this returned 2.
    CHECK(fmi_pread_worker_count(floor_bytes, 8) == 1);
    CHECK(fmi_pread_worker_count(floor_bytes, 2) == 1);

    // Just under two floors is still one worker's worth of work.
    CHECK(fmi_pread_worker_count(2 * floor_bytes - 1, 8) == 1);

    // Two full floors can support two workers, not three.
    CHECK(fmi_pread_worker_count(2 * floor_bytes, 8) == 2);
    CHECK(fmi_pread_worker_count(3 * floor_bytes - 1, 8) == 2);
}

TEST_CASE("fmi_pread_worker_count: honors the caller's request when it fits") {
    const size_t floor_bytes = FMI_PREAD_MIN_CHUNK;

    // Plenty of bytes to go around: the caller's count is returned as-is.
    CHECK(fmi_pread_worker_count(64 * floor_bytes, 8) == 8);
    CHECK(fmi_pread_worker_count(64 * floor_bytes, 1) == 1);

    // hg38-scale cp_occ (~6 GB) against the load's own 8-worker cap.
    CHECK(fmi_pread_worker_count(6UL << 30, 8) == 8);
}

TEST_CASE("fmi_pread_worker_count: degenerate inputs yield one worker") {
    // Sub-floor arrays cannot honor the floor at any count; one worker is the
    // only correct answer (the caller short-circuits nbytes == 0 before this,
    // but the arithmetic must not divide by zero or return 0 regardless).
    CHECK(fmi_pread_worker_count(0, 8) == 1);
    CHECK(fmi_pread_worker_count(1, 8) == 1);
    CHECK(fmi_pread_worker_count(FMI_PREAD_MIN_CHUNK - 1, 8) == 1);

    // Non-positive requests clamp up, never to 0 (a 0 return would make
    // parallel_pread divide by zero when sizing its chunks).
    CHECK(fmi_pread_worker_count(64 * FMI_PREAD_MIN_CHUNK, 0) == 1);
    CHECK(fmi_pread_worker_count(64 * FMI_PREAD_MIN_CHUNK, -4) == 1);
}

TEST_CASE("fmi_pread_worker_count: the resulting split respects the floor") {
    // The property the count exists to guarantee, checked against the same
    // chunk arithmetic parallel_pread uses: every worker gets at least
    // FMI_PREAD_MIN_CHUNK, unless the whole array is smaller than that.
    const size_t floor_bytes = FMI_PREAD_MIN_CHUNK;
    const size_t sizes[] = {
        1, floor_bytes - 1, floor_bytes, floor_bytes + 1,
        2 * floor_bytes - 1, 2 * floor_bytes, 5 * floor_bytes + 12345,
        1UL << 30, 6UL << 30,
    };

    for (size_t i = 0; i < sizeof(sizes) / sizeof(sizes[0]); ++i) {
        const size_t nbytes = sizes[i];
        for (int requested = 1; requested <= 16; ++requested) {
            const int n = fmi_pread_worker_count(nbytes, requested);
            REQUIRE(n >= 1);
            CHECK(n <= requested);
            // Chunks stay balanced: at most one byte between largest and
            // smallest, so no worker is left holding the whole array.
            CHECK(max_chunk_bytes(nbytes, n) - min_chunk_bytes(nbytes, n) <= 1);
            if (nbytes >= floor_bytes) {
                CHECK(min_chunk_bytes(nbytes, n) >= floor_bytes);
            } else {
                CHECK(n == 1);
            }
        }
    }
}
