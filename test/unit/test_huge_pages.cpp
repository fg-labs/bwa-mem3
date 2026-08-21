// Unit tests for bwamem_huge_pages_needed (the --huge-pages page-count math).
//
// The reservation itself touches the OS hugetlb pool and mimalloc, so it is
// exercised end-to-end on a Linux host; here we pin the pure sizing contract
// that decides how many 1 GB pages to reserve for a given index footprint:
//
//   * ceil((bwt_bytes + pac_bytes) / 1 GiB) + a fixed 2-page margin
//   * pac_bytes <= 0 is ignored (only the FM-index counts)
//   * bwt_bytes <= 0 returns 0 -> caller declines to reserve
//
// Inputs are synthetic multiples of 1 GiB; no real index or fixture is read.

#include "doctest/doctest.h"

#include "bwa_hugepages.h"

static const long long GB = 1024LL * 1024LL * 1024LL;

TEST_CASE("huge_pages_needed: bwt footprint rounds up to whole 1 GB pages + margin"
          * doctest::test_suite("unit/huge_pages")) {
    SUBCASE("unknown index size declines (bwt <= 0 returns 0)") {
        CHECK(bwamem_huge_pages_needed(0, 0) == 0);
        CHECK(bwamem_huge_pages_needed(-1, 5 * GB) == 0);
    }
    SUBCASE("exactly 1 GiB -> 1 page + 2 margin") {
        CHECK(bwamem_huge_pages_needed(GB, 0) == 3);
    }
    SUBCASE("1 byte over 1 GiB rounds up to 2 pages + margin") {
        CHECK(bwamem_huge_pages_needed(GB + 1, 0) == 4);
    }
    SUBCASE("just under 1 GiB still needs a whole page + margin") {
        CHECK(bwamem_huge_pages_needed(1, 0) == 3);
    }
}

TEST_CASE("huge_pages_needed: pac adds to the footprint, non-positive pac ignored"
          * doctest::test_suite("unit/huge_pages")) {
    SUBCASE("pac bytes add to the bwt footprint") {
        CHECK(bwamem_huge_pages_needed(2 * GB, GB) == 5);   // 3 pages + 2 margin
        CHECK(bwamem_huge_pages_needed(2 * GB, 0) == 4);    // 2 pages + 2 margin
    }
    SUBCASE("non-positive pac is ignored, not subtracted") {
        CHECK(bwamem_huge_pages_needed(2 * GB, -100) == bwamem_huge_pages_needed(2 * GB, 0));
    }
}

TEST_CASE("huge_pages_needed: large index sizes to 13 pages"
          * doctest::test_suite("unit/huge_pages")) {
    // A ~10.5 GiB index (an hg38-scale FM-index + packed reference) rounds to
    // 11 whole pages, + 2 margin = 13. Expressed synthetically, two ways.
    SUBCASE("bwt alone that rounds to 11 pages") {
        CHECK(bwamem_huge_pages_needed(11 * GB, 0) == 13);
        CHECK(bwamem_huge_pages_needed(10 * GB + 1, 0) == 13); // ceil(10 GiB + 1 B) = 11
    }
    SUBCASE("bwt + pac summing to 11 pages") {
        CHECK(bwamem_huge_pages_needed(10 * GB, GB) == 13);
    }
}

TEST_CASE("classify_huge_reservation: rc==0 alone is not full success"
          * doctest::test_suite("unit/huge_pages")) {
    // mimalloc returns rc==0 once at least one page per NUMA node is reserved,
    // which can be fewer than `need`. The free-count delta (before - after) is
    // the true count, so a partial reservation is detected instead of reported
    // as a clean win.
    long got = -12345;

    SUBCASE("full: rc==0 and the free count dropped by need") {
        // need=14, pool went 20 -> 6, so 14 pages actually landed.
        CHECK(bwamem_classify_huge_reservation(0, 20, 6, 14, &got) == BWAMEM_HUGE_FULL);
        CHECK(got == 14);
    }
    SUBCASE("partial: rc==0 but the free count dropped by fewer than need") {
        // need=14, pool went 20 -> 11, so only 9 pages landed — the bug this guards.
        CHECK(bwamem_classify_huge_reservation(0, 20, 11, 14, &got) == BWAMEM_HUGE_PARTIAL);
        CHECK(got == 9);
    }
    SUBCASE("full: exactly need pages landed (boundary)") {
        CHECK(bwamem_classify_huge_reservation(0, 14, 0, 14, &got) == BWAMEM_HUGE_FULL);
        CHECK(got == 14);
    }
    SUBCASE("failed: rc!=0 regardless of the count") {
        CHECK(bwamem_classify_huge_reservation(-1, 20, 6, 14, &got) == BWAMEM_HUGE_FAILED);
    }
    SUBCASE("unmeasurable free reads: rc==0 is UNVERIFIED, count reported -1") {
        // A -1 from either free-count read (pool vanished, /sys unreadable) means
        // the delta cannot be computed. mimalloc's rc==0 only proves >= 1 page per
        // NUMA node, so with no trustworthy per-call count we must not claim a full
        // reservation — classify as UNVERIFIED.
        CHECK(bwamem_classify_huge_reservation(0, -1, 6, 14, &got) == BWAMEM_HUGE_UNVERIFIED);
        CHECK(got == -1);
        CHECK(bwamem_classify_huge_reservation(0, 20, -1, 14, &got) == BWAMEM_HUGE_UNVERIFIED);
        CHECK(got == -1);
    }
    SUBCASE("free count rose (another process freed pages): UNVERIFIED, not full") {
        // The free_hugepages count is host-global; a rise means another process
        // perturbed the pool, so the delta is not ours to interpret.
        CHECK(bwamem_classify_huge_reservation(0, 6, 20, 14, &got) == BWAMEM_HUGE_UNVERIFIED);
        CHECK(got == -1);
    }
    SUBCASE("null reserved_out is allowed") {
        CHECK(bwamem_classify_huge_reservation(0, 20, 11, 14, nullptr) == BWAMEM_HUGE_PARTIAL);
    }
}
