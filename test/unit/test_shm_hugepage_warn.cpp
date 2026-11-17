// Unit tests for the shm huge-page warning (src/bwa_shm.cpp): when it is
// printed (bwa_shm_hugepage_warn_wanted) and its text. The smaps parser whose
// verdict feeds it is tested in test_smaps_vma_parser.cpp. Nothing here
// depends on the host's /dev/shm mount.

#include "doctest/doctest.h"

#include "bwa_shm.h"

#include <string>

TEST_CASE("shm hugepage warning: one line that names the page sizes, the fix and the restage"
          * doctest::test_suite("unit/shm_hugepage_warn")) {
    const uint64_t KB4 = 4ULL << 10, MB2 = 2ULL << 20;
    char buf[1024];
    const int n = bwa_shm_format_hugepage_warning(buf, sizeof(buf), "bwa_shm_attach", "hg38", KB4, MB2);
    REQUIRE(n > 0);
    REQUIRE((size_t)n < sizeof(buf));
    const std::string s(buf);
    CHECK(s.rfind("[W::bwa_shm_attach] ", 0) == 0);
    CHECK(s.find('\n') == s.size() - 1);
    CHECK(s.find("'hg38'") != std::string::npos);
    CHECK(s.find("on 4 KB pages, not 2 MB huge pages") != std::string::npos);
    CHECK(s.find("uses more CPU than against a disk-loaded index") == std::string::npos);
    CHECK(s.find("a quarter to a third more") == std::string::npos);
    CHECK(s.find("sudo mount -o remount,huge=advise /dev/shm") != std::string::npos);
    CHECK(s.find("/etc/fstab") != std::string::npos);
    CHECK(s.find("bwa-mem3 shm -d") != std::string::npos);
    // A NULL buffer returns the length, which the caller sizes its buffer with.
    CHECK(bwa_shm_format_hugepage_warning(NULL, 0, "bwa_shm_attach", "hg38", KB4, MB2) == n);
    // An arm64 kernel with 64 KB base pages has 512 MB PMD pages.
    char b64[1024];
    bwa_shm_format_hugepage_warning(b64, sizeof(b64), "bwa_shm_stage", "hg38", 64ULL << 10, 512ULL << 20);
    CHECK(std::string(b64).find("on 64 KB pages, not 512 MB huge pages") != std::string::npos);
}

TEST_CASE("shm hugepage warning: only where the remount could help"
          * doctest::test_suite("unit/shm_hugepage_warn")) {
    const uint64_t MB2 = 2ULL << 20, GB21 = 21ULL << 30;
    const char *advise = "always within_size [advise] never deny force\n";
    const char *never = "always within_size advise [never] deny force\n";
    const char *deny = "always within_size advise never [deny] force\n";
    // A resident segment with no huge page, of one PMD page or more, on a
    // kernel with shmem THP not set to deny: the remount and restage can put
    // it on huge pages.
    CHECK(bwa_shm_hugepage_warn_wanted(0, MB2, MB2, advise) == 1);
    CHECK(bwa_shm_hugepage_warn_wanted(0, GB21, MB2, never) == 1);
    // On huge pages already, or no verdict from smaps.
    CHECK(bwa_shm_hugepage_warn_wanted(1, GB21, MB2, advise) == 0);
    CHECK(bwa_shm_hugepage_warn_wanted(-1, GB21, MB2, advise) == 0);
    // Under one PMD page (phiX, the CI fixtures): never PMD-mapped, whatever
    // the mount, so the advice cannot help. With 64 KB base pages the PMD page
    // is 512 MB.
    CHECK(bwa_shm_hugepage_warn_wanted(0, MB2 - 1, MB2, advise) == 0);
    CHECK(bwa_shm_hugepage_warn_wanted(0, 64 << 10, MB2, advise) == 0);
    CHECK(bwa_shm_hugepage_warn_wanted(0, 0, MB2, advise) == 0);
    CHECK(bwa_shm_hugepage_warn_wanted(0, 100ULL << 20, 512ULL << 20, advise) == 0);
    CHECK(bwa_shm_hugepage_warn_wanted(0, 512ULL << 20, 512ULL << 20, advise) == 1);
    // shmem_enabled = deny turns huge pages off for every tmpfs mount, whatever
    // its huge= option, so a remount cannot help either.
    CHECK(bwa_shm_hugepage_warn_wanted(0, GB21, MB2, deny) == 0);
    CHECK(bwa_shm_hugepage_warn_wanted(0, GB21, MB2, "[deny]") == 0);
    // No shmem THP in the kernel (the sysfs file is absent).
    CHECK(bwa_shm_hugepage_warn_wanted(0, GB21, MB2, NULL) == 0);
}
