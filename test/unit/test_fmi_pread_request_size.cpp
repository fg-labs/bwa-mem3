// test/unit/test_fmi_pread_request_size.cpp
//
// Per-call size clamp for the parallel index read.
//
// parallel_pread splits an index array into index_load_threads() chunks (capped
// at 8) and each worker loops pread() until its chunk is consumed. Handing the
// whole remaining chunk to a single pread() is fine on Linux but fails on macOS,
// where a count greater than INT_MAX returns EINVAL. With the 8-worker cap that
// makes every index over ~17.2GB unloadable: the hg38 --meth FM-index is 20.9GB,
// so an 8-way split asks for 2.61GB per call and the load aborts with
//
//     ERROR: pread failed during index load: Invalid argument
//
// The failure is silent on Linux and total on macOS, so the invariant worth
// pinning is simply that the requested count is always a legal pread() count.

#include "doctest/doctest.h"
#include "../../src/FMI_search.h"

#include <cstddef>
#include <climits>

TEST_CASE("fmi_pread_request_size never asks pread() for more than INT_MAX")
{
    const size_t int_max = (size_t)INT_MAX;

    // The chunk sizes an 8-way split produces for real indexes, plus the
    // boundary itself. Every one must be a legal pread() count.
    const size_t sizes[] = {
        0,
        1,
        (size_t)1 << 20,                 // 1MiB
        int_max - 1,
        int_max,
        int_max + 1,                     // first illegal count on macOS
        2614094383u,                     // 20.9GB --meth index / 8 workers
        (size_t)20912755063u,            // the whole --meth index in one chunk
    };
    for (size_t n : sizes)
        CHECK(fmi_pread_request_size(n) <= int_max);
}

TEST_CASE("fmi_pread_request_size makes progress and never over-reads")
{
    // Must never return more than remains, or the worker would read past its
    // chunk and corrupt the neighbouring array.
    CHECK(fmi_pread_request_size(0) == 0);
    CHECK(fmi_pread_request_size(1) == 1);
    CHECK(fmi_pread_request_size(4096) == 4096);

    // Must always make progress on a huge chunk, otherwise the worker loop
    // spins forever instead of aborting or completing.
    CHECK(fmi_pread_request_size((size_t)20912755063u) > 0);

    // A chunk larger than the clamp is consumed in a finite number of calls.
    size_t remaining = 2614094383u, calls = 0;
    while (remaining > 0) {
        size_t want = fmi_pread_request_size(remaining);
        REQUIRE(want > 0);
        REQUIRE(want <= remaining);
        remaining -= want;
        if (++calls > 64) break;         // guard against a non-progressing clamp
    }
    CHECK(remaining == 0);
    CHECK(calls <= 64);
}
