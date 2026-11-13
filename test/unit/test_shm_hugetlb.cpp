// Unit tests for `bwa-mem3 shm --huge-pages`: the /proc/mounts parser that
// finds a hugetlbfs mount with 1 GB pages, and the redirect stub that the POSIX
// segment holds in place of the index. The stub's first 64 bytes must read as
// a bwa_shm_header_t, so the staging state machine (stager PID, stage token,
// magic published last) works on a redirect unchanged.
//
// The stage, attach and removal paths are exercised by shm_stage_test.

#include "doctest/doctest.h"

#include "bwa_shm.h"

#include <stddef.h>
#include <stdint.h>
#include <string.h>

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

TEST_CASE("hugetlb_stub: its header fields sit where the index header keeps them"
          * doctest::test_suite("unit/shm_hugetlb")) {
    // The probe, the wait loop and the token check read a stub through
    // bwa_shm_header_t, so these must line up on every platform.
    CHECK(sizeof(bwa_shm_hugetlb_stub_t) == 64 + 4096);
    CHECK(offsetof(bwa_shm_hugetlb_stub_t, magic) == offsetof(bwa_shm_header_t, magic));
    CHECK(offsetof(bwa_shm_hugetlb_stub_t, kind) == offsetof(bwa_shm_header_t, n_sections));
    CHECK(offsetof(bwa_shm_hugetlb_stub_t, total_size) ==
          offsetof(bwa_shm_header_t, total_size));
    CHECK(offsetof(bwa_shm_hugetlb_stub_t, stager_pid) == offsetof(bwa_shm_header_t, reserved));
    CHECK(offsetof(bwa_shm_hugetlb_stub_t, stage_token) ==
          offsetof(bwa_shm_header_t, reserved) + sizeof(uint64_t));
    CHECK(offsetof(bwa_shm_hugetlb_stub_t, path) == sizeof(bwa_shm_header_t));
}

TEST_CASE("hugetlb_stub: encode then decode round-trips"
          * doctest::test_suite("unit/shm_hugetlb")) {
    bwa_shm_hugetlb_stub_t s = good_stub();
    s.stager_pid = 1234;
    s.stage_token = 0xabcdef01ull;
    bwa_shm_hugetlb_stub_t d;
    REQUIRE(bwa_shm_hugetlb_stub_decode(&s, sizeof(s), &d) == 0);
    CHECK(d.magic == BWA_SHM_HUGETLB_STUB_MAGIC);
    CHECK(d.version == BWA_SHM_HUGETLB_STUB_VERSION);
    CHECK(d.kind == BWA_SHM_HUGETLB_KIND_1G);
    CHECK(d.total_size == 11 * GIB - 5);
    CHECK(d.file_size == 11 * GIB);
    CHECK(d.stager_pid == 1234);
    CHECK(d.stage_token == 0xabcdef01ull);
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
    CHECK(bwa_shm_hugetlb_stub_encode("/mnt/bwaidx-x", 0, GIB, &s) == -1);
    CHECK(bwa_shm_hugetlb_stub_encode("/mnt/bwaidx-x", GIB + 1, GIB, &s) == -1);
    CHECK(bwa_shm_hugetlb_stub_encode("/mnt/bwaidx-x", 1, GIB + 4096, &s) == -1);
    CHECK(bwa_shm_hugetlb_stub_encode("/mnt/bwaidx-x", 1, 0, &s) == -1);
    CHECK(bwa_shm_hugetlb_stub_encode(NULL, 1, GIB, &s) == -1);
    CHECK(bwa_shm_hugetlb_stub_encode("/mnt/bwaidx-x", 1, GIB, NULL) == -1);
    std::string long_path = "/" + std::string(BWA_SHM_HUGETLB_PATH_MAX, 'a') + "/bwaidx-x";
    CHECK(bwa_shm_hugetlb_stub_encode(long_path.c_str(), 1, GIB, &s) == -1);
}

TEST_CASE("hugetlb_stub: decode rejects anything but a well-formed complete stub"
          * doctest::test_suite("unit/shm_hugetlb")) {
    bwa_shm_hugetlb_stub_t d;
    {
        const bwa_shm_hugetlb_stub_t s = good_stub();
        CHECK(bwa_shm_hugetlb_stub_decode(&s, sizeof(s) - 1, &d) == -1);
        CHECK(bwa_shm_hugetlb_stub_decode(&s, sizeof(s) + 1, &d) == -1);
        CHECK(bwa_shm_hugetlb_stub_decode(NULL, sizeof(s), &d) == -1);
        CHECK(bwa_shm_hugetlb_stub_decode(&s, sizeof(s), NULL) == -1);
    }
    {
        bwa_shm_hugetlb_stub_t s = good_stub();
        s.magic = BWA_SHM_MAGIC;
        CHECK(bwa_shm_hugetlb_stub_decode(&s, sizeof(s), &d) == -1);
    }
    {
        // Magic 0 is a stage still in progress, not a stub to follow.
        bwa_shm_hugetlb_stub_t s = good_stub();
        s.magic = 0;
        CHECK(bwa_shm_hugetlb_stub_decode(&s, sizeof(s), &d) == -1);
    }
    {
        bwa_shm_hugetlb_stub_t s = good_stub();
        s.version = BWA_SHM_HUGETLB_STUB_VERSION + 1;
        CHECK(bwa_shm_hugetlb_stub_decode(&s, sizeof(s), &d) == -1);
    }
    {
        bwa_shm_hugetlb_stub_t s = good_stub();
        s.kind = 0;
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
        bwa_shm_hugetlb_stub_t s = good_stub();
        strcpy(s.path, "relative/bwaidx-hg38.fa");
        CHECK(bwa_shm_hugetlb_stub_decode(&s, sizeof(s), &d) == -1);
    }
    {
        bwa_shm_hugetlb_stub_t s = good_stub();
        s.total_size = 0;
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
