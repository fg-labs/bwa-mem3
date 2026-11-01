// Unit tests for --huge-pages: the 1 GB rounding, the mapping and pre-fault
// helpers' failure paths, and FMI_search::load_index's fallback, which must
// load the same arrays whether or not the 1 GB mapping works.
//
// The 1 GB mapping itself needs a Linux host with a reserved 1 GB hugetlb pool.
// The load test runs the real mapping attempt (BWAMEM_HUGETLB_ON) as well as
// the forced fallback, so on such a host it also exercises the mapped path;
// elsewhere both take the fallback.

#include "doctest/doctest.h"

#include "bwa_hugepages.h"
#include "stderr_capture.h"
#include "temp_index.h"

#include <errno.h>
#include <stdint.h>
#include <string.h>

#include <random>
#include <string>
#include <vector>

static const long long GB = 1024LL * 1024LL * 1024LL;
TEST_CASE("hugetlb_round_up: whole 1 GiB pages, 0 for empty or overflowing sizes"
          * doctest::test_suite("unit/huge_pages")) {
    const size_t G = (size_t)GB;
    CHECK(bwamem_hugetlb_round_up(0) == 0);
    CHECK(bwamem_hugetlb_round_up(1) == G);
    CHECK(bwamem_hugetlb_round_up(G) == G);
    CHECK(bwamem_hugetlb_round_up(G + 1) == 2 * G);
    CHECK(bwamem_hugetlb_round_up(SIZE_MAX) == 0);
    CHECK(bwamem_hugetlb_round_up(SIZE_MAX - G + 2) == 0);
}

namespace {

// Sets the index mapping mode for one scope and restores OFF afterwards, so a
// failing CHECK cannot leak the mode into other tests.
struct HugetlbModeScope {
    explicit HugetlbModeScope(bwamem_hugetlb_mode_t m) { bwamem_hugetlb_set_index_mode(m); }
    ~HugetlbModeScope() { bwamem_hugetlb_set_index_mode(BWAMEM_HUGETLB_OFF); }
};

std::string random_acgt(std::mt19937 &rng, size_t n) {
    std::string s(n, 'A');
    for (char &c : s) c = "ACGT"[rng() & 3];
    return s;
}

} // namespace

TEST_CASE("hugetlb_mode_scope: sets the mode and restores the one it found"
          * doctest::test_suite("unit/huge_pages")) {
    HugetlbModeScope outer(BWAMEM_HUGETLB_OFF);
    {
        bwamem_hugetlb_mode_scope on(BWAMEM_HUGETLB_ON);
        CHECK(bwamem_hugetlb_index_mode() == BWAMEM_HUGETLB_ON);
        {
            bwamem_hugetlb_mode_scope fail(BWAMEM_HUGETLB_FORCE_FAIL);
            CHECK(bwamem_hugetlb_index_mode() == BWAMEM_HUGETLB_FORCE_FAIL);
        }
        // Restores ON, the mode it found, not OFF.
        CHECK(bwamem_hugetlb_index_mode() == BWAMEM_HUGETLB_ON);
    }
    CHECK(bwamem_hugetlb_index_mode() == BWAMEM_HUGETLB_OFF);
}

TEST_CASE("hugetlb_map: failures leave the buffer unmapped"
          * doctest::test_suite("unit/huge_pages")) {
    CHECK(bwamem_hugetlb_index_mode() == BWAMEM_HUGETLB_OFF);
    bwamem_hugetlb_buf_t b = {(void *)&b, 123};
    SUBCASE("a zero size is refused") {
        CHECK(bwamem_hugetlb_map(0, &b) != 0);
        CHECK(b.base == NULL);
        CHECK(b.map_bytes == 0);
    }
    SUBCASE("forced failure (the fallback the tests rely on)") {
        HugetlbModeScope mode(BWAMEM_HUGETLB_FORCE_FAIL);
        const int err = bwamem_hugetlb_map(4096, &b);
#if defined(__linux__)
        CHECK(err == ECANCELED);
#else
        CHECK(err == ENOSYS);
#endif
        CHECK(b.base == NULL);
        CHECK(b.map_bytes == 0);
    }
    SUBCASE("unmapped buffers: prefault and unmap are no-ops") {
        bwamem_hugetlb_buf_t bufs[2] = {{NULL, 0}, {NULL, 0}};
        CHECK(bwamem_hugetlb_prefault(bufs, 2, 4) == 0);
        CHECK(bwamem_hugetlb_prefault(NULL, 0, 4) == 0);
        bwamem_hugetlb_unmap(&bufs[0]);
        bwamem_hugetlb_unmap(NULL);
        CHECK(bufs[0].base == NULL);
    }
}

TEST_CASE("load_index with --huge-pages: same arrays whether or not the 1 GB mapping works"
          * doctest::test_suite("unit/huge_pages")) {
    // The index is built and loaded once with the mode off; a second FMI_search
    // then loads the same files with the mode on. FORCE_FAIL always takes the
    // fallback; ON really tries mmap(MAP_HUGETLB), which on a host without a
    // 1 GB pool also falls back, and on one with a pool maps the arrays. Either
    // way the occurrence table and every SA lookup must match the default load.
    std::mt19937 rng(20261101);
    const std::vector<bwa_tests::Contig> ref = {
        {"chrA", random_acgt(rng, 1500)}, {"chrB", random_acgt(rng, 700)}};
    bwa_tests::TempIndex base(ref);
    FMI_search *want = base.fmi();

    for (bwamem_hugetlb_mode_t m : {BWAMEM_HUGETLB_FORCE_FAIL, BWAMEM_HUGETLB_ON}) {
        CAPTURE((int)m);
        std::string log;
        FMI_search got(base.prefix().c_str());
        {
            HugetlbModeScope mode(m);
            bwa_tests::StderrCapture cap;
            got.load_index(/*load_pac=*/false, /*n_threads=*/2);
            log = cap.text();
        }
        CAPTURE(log);
        // One outcome line per buffer.
        for (const char *what : {"occurrence table", "suffix-array samples"}) {
            CAPTURE(what);
            const auto has = [&log](const std::string &needle) {
                return log.find(needle) != std::string::npos;
            };
            const bool mapped = has(std::string("--huge-pages: ") + what + " (");
            const bool not_mapped = has(std::string("--huge-pages: could not map the ") + what);
            const bool not_faulted = has(std::string("--huge-pages: could not fault in the 1 GB pages for the ") + what);
            CHECK((int)mapped + (int)not_mapped + (int)not_faulted == 1);
            if (m == BWAMEM_HUGETLB_FORCE_FAIL) CHECK(not_mapped);
        }
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
    }
    CHECK(bwamem_hugetlb_index_mode() == BWAMEM_HUGETLB_OFF);
}
