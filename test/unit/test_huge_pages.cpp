// Unit tests for the --huge-pages helpers: the page-count math
// (bwamem_huge_pages_needed), the reservation classifier, the available-page
// count (free_hugepages minus resv_hugepages), and its sysfs reader, which is
// pointed at a temporary directory of synthetic counter files.
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

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <string>

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
    // which can be fewer than `need`. The available-count delta (before - after)
    // is the count reserved, so a partial reservation is detected instead of
    // reported as a clean success. The inputs are bwamem_huge_pages_available()
    // values (free minus reserved), not the raw free_hugepages counter.
    long got = -12345;

    SUBCASE("full: rc==0 and the available count dropped by need") {
        // need=14, pool went 20 -> 6, so 14 pages actually landed.
        CHECK(bwamem_classify_huge_reservation(0, 20, 6, 14, &got) == BWAMEM_HUGE_FULL);
        CHECK(got == 14);
    }
    SUBCASE("partial: rc==0 but the available count dropped by fewer than need") {
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
    SUBCASE("unmeasurable available reads: rc==0 is UNVERIFIED, count reported -1") {
        // A -1 from either available-count read (pool vanished, /sys unreadable) means
        // the delta cannot be computed. mimalloc's rc==0 only proves >= 1 page per
        // NUMA node, so with no trustworthy per-call count we must not claim a full
        // reservation — classify as UNVERIFIED.
        CHECK(bwamem_classify_huge_reservation(0, -1, 6, 14, &got) == BWAMEM_HUGE_UNVERIFIED);
        CHECK(got == -1);
        CHECK(bwamem_classify_huge_reservation(0, 20, -1, 14, &got) == BWAMEM_HUGE_UNVERIFIED);
        CHECK(got == -1);
    }
    SUBCASE("available count rose (another process freed pages): UNVERIFIED, not full") {
        // The available count is host-global; a rise means another process
        // perturbed the pool, so the delta is not ours to interpret.
        CHECK(bwamem_classify_huge_reservation(0, 6, 20, 14, &got) == BWAMEM_HUGE_UNVERIFIED);
        CHECK(got == -1);
    }
    SUBCASE("null reserved_out is allowed") {
        CHECK(bwamem_classify_huge_reservation(0, 20, 11, 14, nullptr) == BWAMEM_HUGE_PARTIAL);
    }
}

TEST_CASE("huge_pages_available: free minus reserved hugetlb pages"
          * doctest::test_suite("unit/huge_pages")) {
    // A hugetlb mapping reserves its pages up front but leaves free_hugepages
    // until each page is touched. Right after mimalloc reserves 13 of a 14-page
    // pool on a single-NUMA-node host it has touched one page (its arena
    // metadata): free=13, resv=12. Only one page is left for
    // anyone else, and the reservation took 13, not the 1 a free-only delta shows.
    SUBCASE("reserved-but-untouched pages are not available") {
        CHECK(bwamem_huge_pages_available(14, 0) == 14);
        CHECK(bwamem_huge_pages_available(13, 12) == 1);
        CHECK(bwamem_huge_pages_available(8, 7) == 1);
    }
    SUBCASE("the difference is a real subtraction, not just 0 or 1") {
        CHECK(bwamem_huge_pages_available(20, 5) == 15);
        CHECK(bwamem_huge_pages_available(14, 13) == 1);
        CHECK(bwamem_huge_pages_available(30, 11) == 19);
    }
    SUBCASE("a mimalloc reservation classifies as FULL on the available count") {
        long got = -1;
        const long before = bwamem_huge_pages_available(14, 0);
        const long after = bwamem_huge_pages_available(13, 12);
        CHECK(bwamem_classify_huge_reservation(0, before, after, 13, &got) == BWAMEM_HUGE_FULL);
        CHECK(got == 13);
    }
    SUBCASE("the classification does not depend on how many pages mimalloc touched") {
        // One metadata page per NUMA node may be touched, so the after-state
        // varies; the reserved count is the same.
        const long before = bwamem_huge_pages_available(14, 0);
        for (long touched = 1; touched <= 3; ++touched) {
            long got = -1;
            const long after = bwamem_huge_pages_available(14 - touched, 13 - touched);
            CHECK(bwamem_classify_huge_reservation(0, before, after, 13, &got) == BWAMEM_HUGE_FULL);
            CHECK(got == 13);
        }
        long got = -1;
        CHECK(bwamem_classify_huge_reservation(0, before, bwamem_huge_pages_available(1, 0), 13, &got)
              == BWAMEM_HUGE_FULL);
        CHECK(got == 13);
    }
    SUBCASE("reserve count at or above free floors at 0") {
        CHECK(bwamem_huge_pages_available(3, 3) == 0);
        CHECK(bwamem_huge_pages_available(3, 5) == 0);
    }
    SUBCASE("unreadable reserve count falls back to the free count") {
        CHECK(bwamem_huge_pages_available(9, -1) == 9);
        CHECK(bwamem_huge_pages_available(9, 0) == 9);   // zero is a valid reading
    }
    SUBCASE("unreadable free count means no pool") {
        CHECK(bwamem_huge_pages_available(-1, 0) == -1);
        CHECK(bwamem_huge_pages_available(-1, 4) == -1);
    }
}

// Make a temp pool directory holding the given counter files (a value of nullptr
// leaves the file out). Returns the directory path.
static std::string make_pool(const char *free_val, const char *resv_val) {
    char tmpl[] = "/tmp/bwamem_hp_testXXXXXX";
    REQUIRE(mkdtemp(tmpl) != nullptr);
    const std::string dir(tmpl);
    auto put = [&](const char *name, const char *val) {
        if (val == nullptr) return;
        FILE *f = fopen((dir + "/" + name).c_str(), "w");
        REQUIRE(f != nullptr);
        fprintf(f, "%s\n", val);
        fclose(f);
    };
    put("free_hugepages", free_val);
    put("resv_hugepages", resv_val);
    return dir;
}

static void remove_pool(const std::string &dir) {
    remove((dir + "/free_hugepages").c_str());
    remove((dir + "/resv_hugepages").c_str());
    rmdir(dir.c_str());
}

TEST_CASE("huge_pages_available_at: reads free_hugepages and resv_hugepages from a pool dir"
          * doctest::test_suite("unit/huge_pages")) {
    // Pins the counter names and the free-minus-resv wiring the reservation uses.
    SUBCASE("free minus reserved") {
        const std::string dir = make_pool("13", "12");
        CHECK(bwamem_avail_hugepages_at(dir.c_str()) == 1);
        CHECK(bwamem_read_hugepage_counter(dir.c_str(), "free_hugepages") == 13);
        CHECK(bwamem_read_hugepage_counter(dir.c_str(), "resv_hugepages") == 12);
        remove_pool(dir);
    }
    SUBCASE("a missing resv_hugepages falls back to the free count") {
        const std::string dir = make_pool("9", nullptr);
        CHECK(bwamem_avail_hugepages_at(dir.c_str()) == 9);
        remove_pool(dir);
    }
    SUBCASE("a missing free_hugepages means no pool") {
        const std::string dir = make_pool(nullptr, "4");
        CHECK(bwamem_avail_hugepages_at(dir.c_str()) == -1);
        remove_pool(dir);
    }
    SUBCASE("a missing directory means no pool") {
        CHECK(bwamem_avail_hugepages_at("/nonexistent/bwamem_hp_test") == -1);
        CHECK(bwamem_read_hugepage_counter("/nonexistent/bwamem_hp_test", "free_hugepages") == -1);
    }
}
