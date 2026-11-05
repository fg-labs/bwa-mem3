// Unit tests for `bwa-mem3 shm --huge-pages`: the /proc/mounts parser that
// finds a hugetlbfs mount with 1 GB pages, the redirect stub that the POSIX
// segment holds in place of the index, and the hugetlbfs file helper's
// failure paths.
//
// Creating the file on 1 GB pages needs a Linux host with a hugetlbfs mount
// (pagesize=1G) and at least one free 1 GB page. Set BWA_TEST_HUGETLB_DIR to
// that mount to also stage a test index into a file there, compare it byte for
// byte with the packed index, and check that it is mapped on 1 GB pages.

#include "doctest/doctest.h"

#include "bwa_shm.h"
#include "temp_dir.h"
#include "temp_index.h"

#include <errno.h>
#include <fcntl.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <random>
#include <string>

namespace {

const uint64_t GIB = 1ull << 30;

int parse(const char *line, std::string *dir) {
    char buf[PATH_MAX];
    const int rc = bwa_shm_parse_hugetlbfs_1g_mount(line, buf, sizeof(buf));
    if (rc == 1) *dir = buf;
    return rc;
}

bwa_shm_hugetlb_stub_t good_stub() {
    bwa_shm_hugetlb_stub_t s;
    REQUIRE(bwa_shm_hugetlb_stub_encode("/mnt/huge-1g/bwaidx-hg38.fa", 11 * GIB - 5,
                                        11 * GIB, &s) == 0);
    return s;
}

} // namespace

TEST_CASE("parse_hugetlbfs_1g_mount: accepts hugetlbfs mounts with 1 GB pages only"
          * doctest::test_suite("unit/shm_hugetlb")) {
    std::string dir;
    // The kernel prints a 1 GiB page size as 1024M.
    CHECK(parse("none /mnt/huge-1g hugetlbfs rw,relatime,pagesize=1024M 0 0\n", &dir) == 1);
    CHECK(dir == "/mnt/huge-1g");
    // Spellings a mount command takes, pagesize not last, no trailing newline.
    CHECK(parse("hugetlbfs /a hugetlbfs rw,pagesize=1G,size=12G 0 0", &dir) == 1);
    CHECK(dir == "/a");
    CHECK(parse("x /b hugetlbfs pagesize=1048576K,rw 0 0", &dir) == 1);
    CHECK(dir == "/b");
    // Tabs as separators.
    CHECK(parse("x\t/c\thugetlbfs\trw,pagesize=1024M\t0\t0", &dir) == 1);
    CHECK(dir == "/c");

    // The default 2 MB mount, other page sizes and other file systems.
    CHECK(parse("hugetlbfs /dev/hugepages hugetlbfs rw,relatime,pagesize=2M 0 0", &dir) == 0);
    CHECK(parse("x /d hugetlbfs rw,pagesize=10240M 0 0", &dir) == 0);
    CHECK(parse("x /d hugetlbfs rw,pagesize=1024 0 0", &dir) == 0);
    CHECK(parse("x /d hugetlbfs rw,relatime 0 0", &dir) == 0);
    CHECK(parse("tmpfs /dev/shm tmpfs rw,pagesize=1024M 0 0", &dir) == 0);
    CHECK(parse("x /d hugetlbfsx rw,pagesize=1024M 0 0", &dir) == 0);
    // An option that only ends in pagesize=1024M is not that option.
    CHECK(parse("x /d hugetlbfs rw,xpagesize=1024M 0 0", &dir) == 0);

    // Malformed lines.
    CHECK(parse("", &dir) == 0);
    CHECK(parse("\n", &dir) == 0);
    CHECK(parse("x /d hugetlbfs", &dir) == 0);
    CHECK(parse("x relative hugetlbfs pagesize=1024M 0 0", &dir) == 0);
    char small[4];
    CHECK(bwa_shm_parse_hugetlbfs_1g_mount("x /abcd hugetlbfs pagesize=1G 0 0",
                                           small, sizeof(small)) == 0);
    CHECK(bwa_shm_parse_hugetlbfs_1g_mount("x /a hugetlbfs pagesize=1G 0 0",
                                           small, sizeof(small)) == 1);
    CHECK(bwa_shm_parse_hugetlbfs_1g_mount(NULL, small, sizeof(small)) == 0);
}

TEST_CASE("parse_hugetlbfs_1g_mount: decodes the octal escapes of /proc/mounts"
          * doctest::test_suite("unit/shm_hugetlb")) {
    std::string dir;
    CHECK(parse("none /mnt/huge\\0401g hugetlbfs rw,pagesize=1024M 0 0", &dir) == 1);
    CHECK(dir == "/mnt/huge 1g");
    CHECK(parse("none /a\\134b\\011c hugetlbfs rw,pagesize=1024M 0 0", &dir) == 1);
    CHECK(dir == "/a\\b\tc");
    // A backslash not followed by three octal digits is kept as is.
    CHECK(parse("none /a\\09x hugetlbfs rw,pagesize=1024M 0 0", &dir) == 1);
    CHECK(dir == "/a\\09x");
    CHECK(parse("none /a\\040 hugetlbfs rw,pagesize=1024M 0 0", &dir) == 1);
    CHECK(dir == "/a ");
    CHECK(parse("none /a\\04 hugetlbfs rw,pagesize=1024M 0 0", &dir) == 1);
    CHECK(dir == "/a\\04");
}

TEST_CASE("hugetlb_stub: encode then decode round-trips"
          * doctest::test_suite("unit/shm_hugetlb")) {
    const bwa_shm_hugetlb_stub_t s = good_stub();
    bwa_shm_hugetlb_stub_t d;
    REQUIRE(bwa_shm_hugetlb_stub_decode(&s, sizeof(s), &d) == 0);
    CHECK(d.magic == BWA_SHM_HUGETLB_STUB_MAGIC);
    CHECK(d.version == BWA_SHM_HUGETLB_STUB_VERSION);
    CHECK(d.total_size == 11 * GIB - 5);
    CHECK(d.file_size == 11 * GIB);
    CHECK(std::string(d.path) == "/mnt/huge-1g/bwaidx-hg38.fa");
    // Distinct from a packed index's header, which the attach tells apart by
    // its first eight bytes.
    CHECK(BWA_SHM_HUGETLB_STUB_MAGIC != BWA_SHM_MAGIC);
    CHECK(memcmp(&s.magic, "BWAHUGE", 8) == 0);
}

TEST_CASE("hugetlb_stub: encode refuses fields that would not decode"
          * doctest::test_suite("unit/shm_hugetlb")) {
    bwa_shm_hugetlb_stub_t s;
    CHECK(bwa_shm_hugetlb_stub_encode("relative/bwaidx-x", 1, GIB, &s) == -1);
    CHECK(bwa_shm_hugetlb_stub_encode("/mnt/x", 1, GIB, &s) == -1);          // not bwaidx-*
    CHECK(bwa_shm_hugetlb_stub_encode("/mnt/bwaidx-", 1, GIB, &s) == -1);    // empty name
    CHECK(bwa_shm_hugetlb_stub_encode("/mnt/bwaidx-x/y", 1, GIB, &s) == -1); // prefix not last
    CHECK(bwa_shm_hugetlb_stub_encode("/mnt/bwaidx-x", 0, GIB, &s) == -1);
    CHECK(bwa_shm_hugetlb_stub_encode("/mnt/bwaidx-x", GIB + 1, GIB, &s) == -1);
    CHECK(bwa_shm_hugetlb_stub_encode("/mnt/bwaidx-x", 1, GIB + 4096, &s) == -1);
    CHECK(bwa_shm_hugetlb_stub_encode(NULL, 1, GIB, &s) == -1);
    CHECK(bwa_shm_hugetlb_stub_encode("/mnt/bwaidx-x", 1, GIB, NULL) == -1);
    std::string long_path = "/" + std::string(PATH_MAX, 'a') + "/bwaidx-x";
    CHECK(bwa_shm_hugetlb_stub_encode(long_path.c_str(), 1, GIB, &s) == -1);
}

TEST_CASE("hugetlb_stub: decode rejects anything but a well-formed stub"
          * doctest::test_suite("unit/shm_hugetlb")) {
    bwa_shm_hugetlb_stub_t d;
    {
        const bwa_shm_hugetlb_stub_t s = good_stub();
        CHECK(bwa_shm_hugetlb_stub_decode(&s, sizeof(s) - 1, &d) == -1);
        CHECK(bwa_shm_hugetlb_stub_decode(NULL, sizeof(s), &d) == -1);
    }
    {
        bwa_shm_hugetlb_stub_t s = good_stub();
        s.magic = BWA_SHM_MAGIC;
        CHECK(bwa_shm_hugetlb_stub_decode(&s, sizeof(s), &d) == -1);
    }
    {
        bwa_shm_hugetlb_stub_t s = good_stub();
        s.version = BWA_SHM_HUGETLB_STUB_VERSION + 1;
        CHECK(bwa_shm_hugetlb_stub_decode(&s, sizeof(s), &d) == -1);
    }
    {
        // No NUL anywhere in the path.
        bwa_shm_hugetlb_stub_t s = good_stub();
        memset(s.path, 'a', sizeof(s.path));
        s.path[0] = '/';
        CHECK(bwa_shm_hugetlb_stub_decode(&s, sizeof(s), &d) == -1);
    }
    {
        // A stub must never point the destroy at an arbitrary file.
        bwa_shm_hugetlb_stub_t s = good_stub();
        strcpy(s.path, "/etc/passwd");
        CHECK(bwa_shm_hugetlb_stub_decode(&s, sizeof(s), &d) == -1);
    }
    {
        bwa_shm_hugetlb_stub_t s = good_stub();
        s.total_size = s.file_size + 1;
        CHECK(bwa_shm_hugetlb_stub_decode(&s, sizeof(s), &d) == -1);
    }
    {
        bwa_shm_hugetlb_stub_t s = good_stub();
        s.file_size -= 4096;
        s.total_size = 1;
        CHECK(bwa_shm_hugetlb_stub_decode(&s, sizeof(s), &d) == -1);
    }
}

TEST_CASE("hugetlb_close: removes only a file this process created"
          * doctest::test_suite("unit/shm_hugetlb")) {
    bwa_tests::TempDir dir;
    dir.write("bwaidx-x", "keep");
    const std::string path = dir.file("bwaidx-x");

    // As bwa_shm_hugetlb_create leaves it after O_EXCL found the file there:
    // the path is set for messages, but there is no descriptor.
    bwa_shm_hugetlb_file_t f = {-1, NULL, 0, {0}};
    snprintf(f.path, sizeof(f.path), "%s", path.c_str());
    bwa_shm_hugetlb_close(&f, /*unlink_file=*/true);
    CHECK(access(path.c_str(), F_OK) == 0);
    CHECK(f.path[0] == '\0');

    // With the descriptor that created it, the file is removed.
    f.fd = open(path.c_str(), O_RDONLY);
    REQUIRE(f.fd >= 0);
    snprintf(f.path, sizeof(f.path), "%s", path.c_str());
    bwa_shm_hugetlb_close(&f, /*unlink_file=*/true);
    CHECK(access(path.c_str(), F_OK) != 0);
    CHECK(f.fd == -1);

    bwa_shm_hugetlb_close(NULL, true);   // no-op
}

TEST_CASE("hugetlb_create: failures leave nothing behind"
          * doctest::test_suite("unit/shm_hugetlb")) {
    bwa_tests::TempDir dir;
    bwa_shm_hugetlb_file_t f;
#if defined(__linux__)
    // An ordinary directory is not a hugetlbfs mount with 1 GB pages.
    CHECK(bwa_shm_hugetlb_create(dir.path().c_str(), "idx", GIB, 1, &f) == ENOTSUP);
    CHECK(access(dir.file("bwaidx-idx").c_str(), F_OK) != 0);
    CHECK(bwa_shm_hugetlb_create((dir.path() + "/missing").c_str(), "idx", GIB, 1, &f) == ENOENT);
    CHECK(bwa_shm_hugetlb_create(dir.path().c_str(), "a/b", GIB, 1, &f) == EINVAL);
    CHECK(bwa_shm_hugetlb_create(dir.path().c_str(), "", GIB, 1, &f) == EINVAL);
    CHECK(bwa_shm_hugetlb_create(dir.path().c_str(), "idx", 0, 1, &f) == EINVAL);
#else
    CHECK(bwa_shm_hugetlb_create(dir.path().c_str(), "idx", GIB, 1, &f) == ENOSYS);
#endif
    CHECK(f.fd == -1);
    CHECK(f.base == NULL);
    CHECK(f.map_len == 0);
}

#if defined(__linux__)

namespace {

// KernelPageSize in kB of the /proc/self/smaps entry holding `addr`, or -1.
long kernel_page_kb(const void *addr) {
    FILE *fp = fopen("/proc/self/smaps", "r");
    if (fp == NULL) return -1;
    char line[512];
    bool in = false;
    long kb = -1;
    const uintptr_t a = (uintptr_t)addr;
    while (fgets(line, sizeof(line), fp) != NULL) {
        unsigned long lo, hi;
        if (sscanf(line, "%lx-%lx ", &lo, &hi) == 2) {   // a mapping's header line
            in = lo <= a && a < hi;
            continue;
        }
        if (in && sscanf(line, "KernelPageSize: %ld kB", &kb) == 1) break;
    }
    fclose(fp);
    return in ? kb : -1;
}

std::string random_acgt(size_t n, unsigned seed) {
    std::mt19937 rng(seed);
    std::string s(n, 'A');
    for (char &c : s) c = "ACGT"[rng() & 3];
    return s;
}

} // namespace

TEST_CASE("hugetlb_create: a packed index on 1 GB pages matches the packed bytes"
          * doctest::test_suite("unit/shm_hugetlb")) {
    const char *hdir = getenv("BWA_TEST_HUGETLB_DIR");
    if (hdir == NULL || hdir[0] == '\0') {
        MESSAGE("BWA_TEST_HUGETLB_DIR not set; skipping the 1 GB page staging check");
        return;
    }
    bwa_tests::TempIndex idx({{"chrA", random_acgt(5000, 7)}, {"chrB", random_acgt(3000, 8)}});
    const std::string name = idx.prefix().substr(idx.prefix().rfind('/') + 1);

    bwa_shm_layout_t layout;
    REQUIRE(bwa_shm_compute(idx.prefix().c_str(), &layout) == 0);
    bwa_shm_hugetlb_file_t f;
    REQUIRE(bwa_shm_hugetlb_create(hdir, name.c_str(), layout.total_size, 2, &f) == 0);
    CHECK(f.map_len == GIB);
    CHECK(kernel_page_kb(f.base) == 1048576);
    REQUIRE(bwa_shm_pack_into(&layout, (uint8_t *)f.base) == 0);
    bwa_shm_layout_free(&layout);

    uint8_t *packed = NULL;
    uint64_t packed_len = 0;
    REQUIRE(bwa_shm_pack_from_disk(idx.prefix().c_str(), &packed, &packed_len) == 0);
    CHECK(memcmp(f.base, packed, packed_len) == 0);
    free(packed);

    // A second create of the same name finds the file and leaves it alone.
    bwa_shm_hugetlb_file_t g;
    CHECK(bwa_shm_hugetlb_create(hdir, name.c_str(), GIB, 1, &g) == EEXIST);
    bwa_shm_hugetlb_close(&g, true);
    const std::string path = f.path;
    CHECK(access(path.c_str(), F_OK) == 0);

    bwa_shm_hugetlb_close(&f, /*unlink_file=*/true);
    CHECK(access(path.c_str(), F_OK) != 0);
}

#endif
