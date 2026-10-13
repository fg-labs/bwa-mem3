// Unit tests for the MADV_HUGEPAGE hint on the index arrays (src/bwa_madvise.h).
//
// madvise() fails with EINVAL when its start address is not page-aligned, and
// the index arrays come from _mm_malloc (64-byte alignment) and calloc.
// bwamem_madv_hugepage therefore advises the page-aligned interior of the
// buffer, computed by bwamem_page_align_inward. This pins:
//
//   1. the pure rounding (start rounded up, end rounded down, empty when no
//      whole page fits) at 4 KB, 16 KB and 64 KB pages, and its guards; and
//   2. on Linux with THP built in, that the wrapper actually marks the aligned
//      interior of an unaligned buffer MADV_HUGEPAGE (VmFlags "hg" in
//      /proc/self/smaps), leaves the partial pages at either end alone, and
//      advises nothing when no whole page fits.

#include "doctest/doctest.h"

#include "../../src/bwa_madvise.h"

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>

#if defined(__linux__)
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

namespace {

struct Range {
    size_t len;
    uintptr_t start;
};

Range inward(uintptr_t p, size_t sz, size_t page)
{
    uintptr_t start = 0;
    const size_t len = bwamem_page_align_inward(p, sz, page, &start);
    return Range{len, start};
}

}  // namespace

TEST_CASE("page_align_inward: rounds the start up and the end down"
          * doctest::test_suite("unit/madvise")) {
    const size_t pages[] = {4096, 16384, 65536};
    for (size_t pg : pages) {
        CAPTURE(pg);
        const uintptr_t base = 1000 * (uintptr_t)pg;   // page-aligned origin

        SUBCASE("already aligned start and end: unchanged") {
            Range r = inward(base, 8 * pg, pg);
            CHECK(r.start == base);
            CHECK(r.len == 8 * pg);
        }
        SUBCASE("64-byte-aligned start (posix_memalign): start rounds up a page") {
            Range r = inward(base + 64, 8 * pg, pg);
            CHECK(r.start == base + pg);
            CHECK(r.len == 7 * pg);   // end base+8pg+64 rounds down to base+8pg
        }
        SUBCASE("aligned start, unaligned end: end rounds down") {
            Range r = inward(base, 8 * pg + 100, pg);
            CHECK(r.start == base);
            CHECK(r.len == 8 * pg);
        }
        SUBCASE("unaligned both ends") {
            Range r = inward(base + 1, 8 * pg - 2, pg);
            CHECK(r.start == base + pg);
            CHECK(r.len == 6 * pg);
        }
        SUBCASE("start one byte below a page boundary") {
            Range r = inward(base + pg - 1, 2 * pg + 1, pg);
            CHECK(r.start == base + pg);
            CHECK(r.len == 2 * pg);
        }
        SUBCASE("exactly one whole page after rounding") {
            Range r = inward(base + 64, 2 * pg - 64, pg);
            CHECK(r.start == base + pg);
            CHECK(r.len == pg);
        }
    }
}

TEST_CASE("page_align_inward: empty when no whole page fits"
          * doctest::test_suite("unit/madvise")) {
    const size_t pages[] = {4096, 16384, 65536};
    for (size_t pg : pages) {
        CAPTURE(pg);
        const uintptr_t base = 1000 * (uintptr_t)pg;
        CHECK(inward(base, 0, pg).len == 0);              // empty input
        CHECK(inward(base, pg - 1, pg).len == 0);         // sub-page, aligned start
        CHECK(inward(base + 64, pg, pg).len == 0);        // one page long, straddles
        CHECK(inward(base + 64, 2 * pg - 65, pg).len == 0);  // ends 1 byte short
        CHECK(inward(base + 1, 1, pg).len == 0);
    }
}

TEST_CASE("page_align_inward: rejects bad page sizes and address overflow"
          * doctest::test_suite("unit/madvise")) {
    uintptr_t start = 42;
    CHECK(bwamem_page_align_inward(4096, 8192, 0, &start) == 0);
    CHECK(bwamem_page_align_inward(3 * 4096, 8192, 3000, &start) == 0);   // not 2^k
    CHECK(bwamem_page_align_inward(UINTPTR_MAX - 100, 4096, 4096, &start) == 0);
    CHECK(bwamem_page_align_inward(UINTPTR_MAX - 10, 5, 4096, &start) == 0);
    // p rounds up fine, but p + sz wraps: the wrapped end lies below the start.
    CHECK(bwamem_page_align_inward(UINTPTR_MAX - 2 * 4096 + 1, 4 * 4096, 4096, &start) == 0);
    CHECK(start == 42);   // untouched on every rejected input
}

#if defined(__linux__) && defined(MADV_HUGEPAGE)
namespace {

// The VmFlags line of the /proc/self/smaps VMA containing `addr`, or "" if absent.
std::string vmflags_at(uintptr_t addr)
{
    FILE *f = fopen("/proc/self/smaps", "r");
    if (f == NULL) return "";
    char line[4096];
    bool in = false;
    std::string flags;
    while (fgets(line, sizeof line, f) != NULL) {
        unsigned long lo = 0, hi = 0;
        if (sscanf(line, "%lx-%lx ", &lo, &hi) == 2) {   // a VMA header line
            in = (addr >= lo && addr < hi);
            continue;
        }
        if (in && strncmp(line, "VmFlags:", 8) == 0) {
            flags = std::string(line + 8);
            break;
        }
    }
    fclose(f);
    return flags;
}

bool has_flag(const std::string &flags, const char *f)
{
    const std::string padded = " " + flags;
    const std::string want = std::string(" ") + f;
    for (size_t i = padded.find(want); i != std::string::npos; i = padded.find(want, i + 1)) {
        const char next = padded[i + want.size()];
        if (next == ' ' || next == '\n' || next == '\0') return true;
    }
    return false;
}

}  // namespace

TEST_CASE("madv_hugepage: advises the aligned interior of an unaligned buffer"
          * doctest::test_suite("unit/madvise")) {
    struct stat st;
    if (stat("/sys/kernel/mm/transparent_hugepage", &st) != 0) {
        MESSAGE("SKIP: kernel has no transparent huge pages");
        return;
    }
    const size_t pg = (size_t)sysconf(_SC_PAGESIZE);
    const size_t npages = 16;
    void *m = mmap(NULL, npages * pg, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    REQUIRE(m != MAP_FAILED);
    const uintptr_t base = (uintptr_t)m;
    // Every live mapping has a VmFlags line; none means smaps cannot serve as
    // the oracle here, and the CHECK_FALSEs below would pass vacuously.
    if (vmflags_at(base).empty()) {
        munmap(m, npages * pg);
        MESSAGE("SKIP: /proc/self/smaps has no VmFlags for a live mapping");
        return;
    }

    // A 64-byte-aligned buffer like the index arrays: [base + 64, base + 14 pages
    // + 64), so pages 1..13 are whole and pages 0 and 14 are partial.
    bwamem_madv_hugepage((void *)(base + 64), (npages - 2) * pg);

    CHECK(has_flag(vmflags_at(base + pg), "hg"));                       // first whole page
    CHECK(has_flag(vmflags_at(base + (npages - 3) * pg), "hg"));        // last whole page
    CHECK_FALSE(has_flag(vmflags_at(base), "hg"));                      // partial head page
    CHECK_FALSE(has_flag(vmflags_at(base + (npages - 2) * pg), "hg"));  // partial tail page

    munmap(m, npages * pg);
}

TEST_CASE("madv_hugepage: advises nothing when no whole page fits"
          * doctest::test_suite("unit/madvise")) {
    struct stat st;
    if (stat("/sys/kernel/mm/transparent_hugepage", &st) != 0) {
        MESSAGE("SKIP: kernel has no transparent huge pages");
        return;
    }
    const size_t pg = (size_t)sysconf(_SC_PAGESIZE);
    void *m = mmap(NULL, 2 * pg, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    REQUIRE(m != MAP_FAILED);
    const uintptr_t base = (uintptr_t)m;
    if (vmflags_at(base).empty()) {
        munmap(m, 2 * pg);
        MESSAGE("SKIP: /proc/self/smaps has no VmFlags for a live mapping");
        return;
    }

    // One page long but 64 bytes in: it straddles two pages and covers neither,
    // so rounding outward instead of inward would advise both.
    bwamem_madv_hugepage((void *)(base + 64), pg);

    CHECK_FALSE(has_flag(vmflags_at(base), "hg"));
    CHECK_FALSE(has_flag(vmflags_at(base + pg), "hg"));

    munmap(m, 2 * pg);
}
#endif
