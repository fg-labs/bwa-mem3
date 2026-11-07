// Unit tests for --huge-pages: the 1 GB rounding, the mapping and pre-fault
// helpers' failure paths, and FMI_search::load_index's fallback, which must
// load the same arrays whether or not the 1 GB mapping works.
//
// The 1 GB mapping itself needs a Linux host with a reserved 1 GB hugetlb pool.
// By default the load test runs the forced fallback, the real ON mode (whose
// size rule keeps a test-sized index on default pages, so it takes no pool
// pages) and, on Linux, the mapped-table path on an ordinary anonymous mapping
// (BWAMEM_HUGETLB_TEST_PLAIN), so none of it depends on the host's pool. Set
// BWA_TEST_HUGETLB_POOL=1 on a host with at least 1 free 1 GB page to also map
// the test index's occurrence table on a real 1 GB page
// (BWAMEM_HUGETLB_ON_ANY_SIZE), check it, and check that freeing the index
// unmaps it. On Linux, with the variable unset or set on a host with no free
// 1 GB page, that case does not run and the test prints why. Elsewhere no
// mapping can succeed, so the case always runs and checks the fallback.

#include "doctest/doctest.h"

#include "bwa_hugepages.h"
#include "stderr_capture.h"
#include "temp_index.h"
#include "run_in_child.h"

#include <errno.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#if defined(__linux__)
#include <sys/mman.h>
#endif

#include <cstdio>
#include <fstream>
#include <memory>
#include <random>
#include <string>
#include <vector>

static const long long GB = 1024LL * 1024LL * 1024LL;

// Free pages in the host's 1 GB hugetlb pool, or -1 if the kernel exposes no
// such pool (not Linux, or a CPU without 1 GB page support).
static long hugetlb_1g_free_pages() {
    std::ifstream f("/sys/kernel/mm/hugepages/hugepages-1048576kB/free_hugepages");
    long n = -1;
    if (!(f >> n)) return -1;
    return n;
}

#if defined(__linux__)
// The KernelPageSize, in kB, of this process's mapping that contains `addr`, or
// -1 if no mapping contains it. Reads /proc/self/smaps, so it sees only this
// process's mappings, unlike the host-wide pool counts.
static long kernel_page_kb_at(const void *addr) {
    FILE *fp = std::fopen("/proc/self/smaps", "r");
    if (fp == NULL) return -1;
    const uintptr_t a = (uintptr_t)addr;
    bool in_vma = false;
    long kb = -1;
    char line[4096];
    while (std::fgets(line, sizeof(line), fp) != NULL) {
        // A mapping's header line starts "lo-hi "; its field lines never parse
        // as two hex numbers joined by '-'.
        unsigned long lo, hi;
        if (std::sscanf(line, "%lx-%lx ", &lo, &hi) == 2) {
            in_vma = (a >= lo && a < hi);
        } else if (in_vma && std::sscanf(line, "KernelPageSize: %ld kB", &kb) == 1) {
            break;
        }
    }
    std::fclose(fp);
    return in_vma ? kb : -1;
}
#endif
TEST_CASE("hugetlb_round_up: whole 1 GiB pages, 0 for empty or overflowing sizes"
          * doctest::test_suite("unit/huge_pages")) {
    const size_t G = (size_t)GB;
    CHECK(bwamem_hugetlb_round_up(0) == 0);
    CHECK(bwamem_hugetlb_round_up(1) == G);
    CHECK(bwamem_hugetlb_round_up(G) == G);
    CHECK(bwamem_hugetlb_round_up(G + 1) == 2 * G);
    CHECK(bwamem_hugetlb_round_up(SIZE_MAX) == 0);
    CHECK(bwamem_hugetlb_round_up(SIZE_MAX - G + 2) == 0);
    // The largest size that still rounds: 2^34 - 1 whole pages on a 64-bit
    // size_t (an off-by-one in the overflow test would return 0 here).
    if (sizeof(size_t) == 8) CHECK(bwamem_hugetlb_round_up(SIZE_MAX - G + 1) == SIZE_MAX - G + 1);
}

namespace {

std::string random_acgt(std::mt19937 &rng, size_t n) {
    std::string s(n, 'A');
    for (char &c : s) c = "ACGT"[rng() & 3];
    return s;
}

} // namespace

TEST_CASE("hugetlb_map: failures leave the buffer unmapped"
          * doctest::test_suite("unit/huge_pages")) {
    bwamem_hugetlb_buf_t b = {(void *)&b, 123};
    SUBCASE("a zero size is refused") {
        CHECK(bwamem_hugetlb_map(0, &b) == EINVAL);
        CHECK(b.base == NULL);
        CHECK(b.map_bytes == 0);
    }
    SUBCASE("an overflowing size is refused") {
        CHECK(bwamem_hugetlb_map(SIZE_MAX, &b) == EINVAL);
        CHECK(b.base == NULL);
        CHECK(b.map_bytes == 0);
    }
#if !defined(__linux__)
    SUBCASE("non-Linux: every mapping fails with ENOSYS") {
        CHECK(bwamem_hugetlb_map(4096, &b) == ENOSYS);
        CHECK(b.base == NULL);
    }
#endif
    SUBCASE("an unmapped buffer: prefault and unmap are no-ops") {
        bwamem_hugetlb_buf_t u = {NULL, 0};
        CHECK(bwamem_hugetlb_prefault(&u, 4) == 0);
        CHECK(bwamem_hugetlb_prefault(NULL, 4) == 0);
        bwamem_hugetlb_unmap(&u);
        bwamem_hugetlb_unmap(NULL);
        CHECK(u.base == NULL);
        CHECK(u.map_bytes == 0);
    }
#if !defined(__linux__)
    SUBCASE("non-Linux: the test mappings also fail with ENOSYS") {
        CHECK(bwamem_hugetlb_map_plain_for_test(4096, &b) == ENOSYS);
        CHECK(b.base == NULL);
        CHECK(b.map_bytes == 0);
        CHECK(bwamem_hugetlb_map_protected_for_test(4096, &b) == ENOSYS);
        CHECK(b.base == NULL);
        CHECK(b.map_bytes == 0);
    }
#endif
#if defined(__linux__)
    SUBCASE("a mapping whose pre-fault fails on worker threads is unmapped and reported") {
        // A 3 GiB PROT_NONE anonymous mapping commits no memory, and
        // MADV_POPULATE_WRITE cannot fault it in for writing (EFAULT; EINVAL on
        // a kernel without it). Its three pages with nthreads=4 make prefault
        // create two worker threads besides the caller, so this runs the
        // thread creation, the start barrier and the shared failure flag. The
        // test seam excludes the caller from page claiming, so only a spawned
        // worker can produce the reported failure.
        const size_t len = 3 * (size_t)GB;
        void *p = mmap(NULL, len, PROT_NONE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
        REQUIRE(p != MAP_FAILED);
        bwamem_hugetlb_buf_t m = {p, len};
        int workers = 0;
        CHECK(bwamem_hugetlb_prefault_workers_for_test(&m, 4, &workers) == 1);
        CHECK(workers > 0);
        CHECK(m.base == NULL);
        CHECK(m.map_bytes == 0);
        // The whole range was unmapped: mincore fails with ENOMEM on every page.
        for (size_t off = 0; off < len; off += (size_t)GB) {
            unsigned char vec = 0;
            errno = 0;
            const int rc = mincore((char *)p + off, 4096, &vec);
            const int err = errno;  // before any CHECK, whose reporting can set errno
            CHECK(rc == -1);
            CHECK(err == ENOMEM);
        }
    }
    SUBCASE("the plain test mapping rounds up to whole 1 GiB and unmaps") {
        CHECK(bwamem_hugetlb_map_plain_for_test(4096, &b) == 0);
        REQUIRE(b.base != NULL);
        CHECK(b.map_bytes == (size_t)GB);
        void *p = b.base;
        bwamem_hugetlb_unmap(&b);
        CHECK(b.base == NULL);
        CHECK(b.map_bytes == 0);
        unsigned char vec = 0;
        CHECK(mincore(p, 4096, &vec) == -1);
    }
#endif
}

#if defined(__linux__)
TEST_CASE("load_index with --huge-pages: mapped table with failed prefault falls back safely"
          * doctest::test_suite("unit/huge_pages")) {
    std::mt19937 rng(20261104);
    bwa_tests::TempIndex base({{"chrA", random_acgt(rng, 2200)}});
    const int status = bwa_tests::run_in_child([&base] {
        FMI_search *want = base.fmi();
        {
            FMI_search got(base.prefix().c_str());
            std::string log;
            {
                bwa_tests::StderrCapture cap;
                got.load_index(false, 2, false, BWAMEM_HUGETLB_TEST_PREFAULT_FAIL);
                log = cap.text();
            }
            const std::string tag = "[M::load_index] --huge-pages: ";
            const size_t at = log.find(tag);
            if (at == std::string::npos || log.find(tag, at + tag.size()) != std::string::npos ||
                log.find(tag + "could not fault in the 1 GB pages for the occurrence table") == std::string::npos)
                _exit(91);
            if (got.reference_seq_len != want->reference_seq_len || got.sa_compx != want->sa_compx ||
                got.cp_occ_size_bytes() != want->cp_occ_size_bytes() ||
                got.sa_ms_byte_size_bytes() != want->sa_ms_byte_size_bytes() ||
                got.sa_ls_word_size_bytes() != want->sa_ls_word_size_bytes()) _exit(92);
            if (memcmp(got.cp_occ_data(), want->cp_occ_data(), (size_t)want->cp_occ_size_bytes()) != 0 ||
                memcmp(got.sa_ms_byte_data(), want->sa_ms_byte_data(), (size_t)want->sa_ms_byte_size_bytes()) != 0 ||
                memcmp(got.sa_ls_word_data(), want->sa_ls_word_data(), (size_t)want->sa_ls_word_size_bytes()) != 0)
                _exit(93);
            for (int64_t r = 0; r < want->reference_seq_len; ++r)
                if (got.get_sa_entry_compressed(r) != want->get_sa_entry_compressed(r)) _exit(94);
        }
        _exit(0);
    });
    CHECK((WIFEXITED(status) && WEXITSTATUS(status) == 0));
}
#endif

TEST_CASE("load_index with --huge-pages: same arrays whether or not the 1 GB mapping works"
          * doctest::test_suite("unit/huge_pages")) {
    // The index is built and loaded once with the mode off; a second FMI_search
    // then loads the same files with each mode. FORCE_FAIL always takes the
    // fallback. ON keeps a test-sized index on default pages (too small for
    // 1 GB pages). TEST_PLAIN (Linux) maps the occurrence table on an ordinary
    // anonymous mapping, so the mapped path runs (the table filled into the
    // mapping, used, and unmapped by the destructor) without a pool; elsewhere
    // it reports that 1 GB pages are Linux-only, as ON_ANY_SIZE does. With
    // BWA_TEST_HUGETLB_POOL=1, ON_ANY_SIZE really maps the occurrence table on
    // a 1 GB page (on a host with a pool) and the mapping is checked. Every way,
    // the occurrence table and every SA lookup must match the default load.
    std::mt19937 rng(20261101);
    const std::vector<bwa_tests::Contig> ref = {
        {"chrA", random_acgt(rng, 1500)}, {"chrB", random_acgt(rng, 700)}};
    bwa_tests::TempIndex base(ref);
    FMI_search *want = base.fmi();

    const char *pool_env = getenv("BWA_TEST_HUGETLB_POOL");
    const bool want_pool = pool_env != NULL && strcmp(pool_env, "1") == 0;
    std::vector<bwamem_hugetlb_mode_t> modes = {BWAMEM_HUGETLB_FORCE_FAIL, BWAMEM_HUGETLB_ON,
                                                BWAMEM_HUGETLB_TEST_PLAIN};
    long pool_free_before = -1;
#if defined(__linux__)
    // The real 1 GB mapping runs only when asked for and a page is free, so a
    // pool-less host skips it with a reason instead of failing CHECK(mapped).
    pool_free_before = hugetlb_1g_free_pages();
    if (!want_pool)
        MESSAGE("SKIP: real 1 GB mapping not checked (set BWA_TEST_HUGETLB_POOL=1 on a host "
                "with a free 1 GB hugetlb page)");
    else if (pool_free_before < 1)
        MESSAGE("SKIP: real 1 GB mapping not checked: BWA_TEST_HUGETLB_POOL=1 but "
                << (pool_free_before < 0 ? "this host has no 1 GB hugetlb pool"
                                         : "the 1 GB hugetlb pool has no free page"));
    else
        modes.push_back(BWAMEM_HUGETLB_ON_ANY_SIZE);
#else
    // No mapping can succeed here: this runs the ENOSYS fallback and its line.
    modes.push_back(BWAMEM_HUGETLB_ON_ANY_SIZE);
#endif
    for (bwamem_hugetlb_mode_t m : modes) {
        CAPTURE((int)m);
        std::string log;
        std::unique_ptr<FMI_search> got_owner(new FMI_search(base.prefix().c_str()));
        FMI_search &got = *got_owner;
        {
            bwa_tests::StderrCapture cap;
            got.load_index(/*load_pac=*/false, /*n_threads=*/2, /*build_seed_jump=*/false, m);
            log = cap.text();
        }
        CAPTURE(log);
        // Exactly one outcome line, and it is one of the five outcomes.
        const std::string tag = "[M::load_index] --huge-pages: ";
        size_t n_lines = 0;
        for (size_t at = log.find(tag); at != std::string::npos; at = log.find(tag, at + 1))
            ++n_lines;
        CHECK(n_lines == 1);
        const auto has = [&log, &tag](const std::string &rest) {
            return log.find(tag + rest) != std::string::npos;
        };
        const bool mapped = has("occurrence table (") &&
                            log.find(" x 1 GB pages\n") != std::string::npos &&
                            log.find("is too small") == std::string::npos;
        const bool too_small = has("occurrence table (") &&
                               log.find("is too small for 1 GB pages") != std::string::npos;
        const bool not_mapped = has("could not map the occurrence table");
        const bool not_faulted = has("could not fault in the 1 GB pages for the occurrence table");
        const bool no_os = has("1 GB pages are supported only on Linux");
        CHECK((int)mapped + (int)too_small + (int)not_mapped + (int)not_faulted + (int)no_os == 1);
        if (m == BWAMEM_HUGETLB_FORCE_FAIL) CHECK(not_mapped);
        if (m == BWAMEM_HUGETLB_ON) CHECK(too_small);
#if defined(__linux__)
        if (m == BWAMEM_HUGETLB_TEST_PLAIN) {
            // The table was read into the mapping, not a fallback allocation.
            CHECK(mapped);
            const uintptr_t page = 4096;
            CHECK(((uintptr_t)got.cp_occ_data() % page) == 0);
        }
        if (m == BWAMEM_HUGETLB_ON_ANY_SIZE) {
            // The mapped occurrence table starts on a 1 GiB boundary; the SA
            // samples are on the default allocation, not 1 GB pages.
            const uintptr_t G1 = (uintptr_t)GB;
            CHECK(mapped);
            CHECK(((uintptr_t)got.cp_occ_data() % G1) == 0);
            CHECK(((uintptr_t)got.sa_ms_byte_data() % G1) != 0);
            CHECK(((uintptr_t)got.sa_ls_word_data() % G1) != 0);
        }
#else
        if (m == BWAMEM_HUGETLB_TEST_PLAIN || m == BWAMEM_HUGETLB_ON_ANY_SIZE) CHECK(no_os);
#endif
        REQUIRE(got.reference_seq_len == want->reference_seq_len);
        REQUIRE(got.sa_compx == want->sa_compx);
        REQUIRE(got.cp_occ_size_bytes() == want->cp_occ_size_bytes());
        CHECK(memcmp(got.cp_occ_data(), want->cp_occ_data(), (size_t)want->cp_occ_size_bytes()) == 0);
        CHECK(memcmp(got.sa_ms_byte_data(), want->sa_ms_byte_data(),
                     (size_t)want->sa_ms_byte_size_bytes()) == 0);
        CHECK(memcmp(got.sa_ls_word_data(), want->sa_ls_word_data(),
                     (size_t)want->sa_ls_word_size_bytes()) == 0);
        int64_t bad = 0;
        for (int64_t r = 0; r < want->reference_seq_len; ++r)
            bad += got.get_sa_entry_compressed(r) != want->get_sa_entry_compressed(r);
        CHECK(bad == 0);
#if defined(__linux__)
        // While the index is loaded its occurrence table sits on a 1 GB page of
        // this process; the destructor unmaps it, which returns the page to the
        // pool. Checked on this process's own mappings, not the host-wide free
        // count, which another process could change during the test.
        const void *occ = got.cp_occ_data();
        if (m == BWAMEM_HUGETLB_ON_ANY_SIZE) CHECK(kernel_page_kb_at(occ) == 1048576);
        got_owner.reset();
        if (m == BWAMEM_HUGETLB_ON_ANY_SIZE) CHECK(kernel_page_kb_at(occ) != 1048576);
#endif
    }
}
