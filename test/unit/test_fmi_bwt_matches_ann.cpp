// Structural index sanity checks and diagnostics, not content provenance.
// Disk failures advise a force rebuild; segment-internal failures advise dropping
// and restaging the segment, not rebuilding files a later attach would ignore.

#include "FMI_search.h"
#include "fm_index_writer.h"
#include "bntseq.h"
#include "temp_dir.h"
#include "temp_index.h"
#include "run_in_child.h"
#include "fmi_seed_api.h"
#include <cstdint>
#include <fcntl.h>
#include <unistd.h>

#include "doctest/doctest.h"
#include "stderr_capture.h"

#include <string>

TEST_SUITE_BEGIN("unit/fmi_index_consistency");

TEST_CASE("fmi_bwt_matches_ann: a consistent pair passes silently on both paths")
{
    for (bool shm : {false, true}) {
        bwa_tests::StderrCapture cap;
        CHECK(fmi_bwt_matches_ann("t", "idx", 2 * 5386 + 1, 5386, shm));
        CHECK(cap.text().empty());
    }
}

TEST_CASE("fmi_bwt_matches_ann: the shm path's advice is shm -d, the disk path's a rebuild")
{
    std::string disk, shm;
    {
        bwa_tests::StderrCapture cap;
        CHECK_FALSE(fmi_bwt_matches_ann("t", "idx", 2 * 5386 + 1, 5000));
        disk = cap.text();
    }
    {
        bwa_tests::StderrCapture cap;
        CHECK_FALSE(fmi_bwt_matches_ann("t", "idx", 2 * 5386 + 1, 5000, /*from_shm=*/true));
        shm = cap.text();
    }
    CHECK(disk.find("Rebuild the index") != std::string::npos);
    CHECK(disk.find("shm -d") == std::string::npos);
    CHECK(shm.find("shared memory") != std::string::npos);
    CHECK(shm.find("`bwa-mem3 shm -d`") != std::string::npos);
    CHECK(shm.find("Rebuild the index") == std::string::npos);
}

TEST_CASE("index consistency rejects impossible lengths before arithmetic")
{
    bwa_tests::StderrCapture cap;
    CHECK_FALSE(fmi_bwt_matches_ann("t", "idx", 100, INT64_MAX));
    CHECK_FALSE(fmi_bwt_matches_ann("t", "idx", 100, -1));
    CHECK_FALSE(fmi_bwt_matches_ann("t", "idx", INT64_MIN, 1, true));
}

TEST_CASE("index size checks accept legacy files and validate every writer tag")
{
    bwa_tests::TempDir dir("index_size");
    const std::string path = dir.file("idx.bwt.2bit.64");
    const int fd = open(path.c_str(), O_CREAT | O_RDWR, 0600);
    REQUIRE(fd >= 0);
    const int64_t rows = 4097;
    for (int c = 0; c <= CP_SHIFT; ++c) {
        CAPTURE(c);
        const Bwt2bitLayout layout = fmi_bwt2bit_layout(rows, c);
        // Independent byte accounting from the stock file layout, not the helper under test.
        const int64_t expected = 48 + (rows / 64 + 1) * 64 + (rows / (1LL << c) + 1) * 5 + 16;
        REQUIRE(layout.total_bytes == expected);
        REQUIRE(ftruncate(fd, layout.total_bytes) == 0);
        const int64_t tag = c;
        REQUIRE(pwrite(fd, &tag, sizeof tag, layout.off_sa_compx_tag) == sizeof tag);
        CHECK(detect_sa_compx(fd, layout.total_bytes, rows, SA_COMPX) == c);
        CHECK(fmi_bwt_size_ok("t", path.c_str(), layout.total_bytes, rows, c, fd));
    }
    const Bwt2bitLayout layout = fmi_bwt2bit_layout(rows, SA_COMPX);
    REQUIRE(ftruncate(fd, layout.off_sa_compx_tag) == 0);
    CHECK(detect_sa_compx(fd, layout.off_sa_compx_tag, rows, SA_COMPX) == SA_COMPX);
    CHECK(fmi_bwt_size_ok("t", path.c_str(), layout.off_sa_compx_tag, rows, SA_COMPX, fd));
    REQUIRE(ftruncate(fd, layout.total_bytes) == 0);
    for (int64_t tag : {int64_t(0), int64_t(6), int64_t(7), INT64_MAX}) {
        CAPTURE(tag);
        REQUIRE(pwrite(fd, &tag, sizeof tag, layout.off_sa_compx_tag) == sizeof tag);
        bwa_tests::StderrCapture cap;
        CHECK_FALSE(fmi_bwt_size_ok("t", path.c_str(), layout.total_bytes, rows, SA_COMPX, fd));
    }
    {
        bwa_tests::StderrCapture cap;
        CHECK_FALSE(fmi_bwt_size_ok("t", path.c_str(), 0, INT64_MAX, 3, fd));
        CHECK_FALSE(fmi_bwt_size_ok("t", path.c_str(), 0, rows, 64, fd));
    }
    REQUIRE(close(fd) == 0);
}

TEST_CASE("metadata header peek validates counts and retains unsigned seeds")
{
    bwa_tests::TempDir dir("ann_header");
    const std::string path = dir.file("ref.ann");
    int64_t len = -1;
    int32_t count = -1;
    uint32_t seed = 0;
    dir.write("ref.ann", "19 2 4294967295\n");
    REQUIRE(bns_read_ann_header(path.c_str(), &len, &count, &seed));
    CHECK(len == 19);
    CHECK(count == 2);
    CHECK(seed == UINT32_MAX);
    for (const char *header : {"-1 1 11\n", "19 -1 11\n", "19 2147483648 11\n", "19\n"}) {
        dir.write("ref.ann", header);
        bwa_tests::StderrCapture cap;
        CHECK_FALSE(bns_read_ann_header(path.c_str(), &len, &count, &seed));
    }
}

TEST_CASE("packed reference checks size and all four length remainders without moving stream")
{
    bwa_tests::TempDir dir("pac_length");
    const std::string path = dir.file("ref.pac");
    for (int64_t len = 16; len < 20; ++len) {
        std::string bytes(6, '\0');
        bytes.back() = (char)(len % 4);
        dir.write("ref.pac", bytes);
        FILE *fp = fopen(path.c_str(), "rb");
        REQUIRE(fp != nullptr);
        REQUIRE(fseek(fp, 2, SEEK_SET) == 0);
        CHECK(bns_pac_fp_size_ok("t", path.c_str(), fp, len));
        CHECK(ftell(fp) == 2);
        REQUIRE(fclose(fp) == 0);
        for (unsigned char bad : {4, 255}) {
            bytes.back() = (char)bad;
            dir.write("ref.pac", bytes);
            fp = fopen(path.c_str(), "rb");
            REQUIRE(fp != nullptr);
            bwa_tests::StderrCapture cap;
            CHECK_FALSE(bns_pac_fp_size_ok("t", path.c_str(), fp, len));
            REQUIRE(fclose(fp) == 0);
        }
    }
}

TEST_CASE("seed facade refuses inconsistent metadata before the record conversion")
{
    bwa_tests::TempIndex fixture({{"guard", std::string(200, 'A')}});
    const std::string ann = fixture.prefix() + ".ann";
    FILE *fp = fopen(ann.c_str(), "r+b");
    REQUIRE(fp != nullptr);
    REQUIRE(fwrite("201", 1, 3, fp) == 3);
    REQUIRE(fclose(fp) == 0);
    bwa_tests::TempDir dir("facade_bad");
    const std::string err = dir.file("err");
    const int status = bwa_tests::run_in_child([&] {
        FmiSeed *f = fmi_seed_open(fixture.prefix().c_str());
        fmi_seed_close(f);
    }, err.c_str());
    REQUIRE(WIFEXITED(status));
    CHECK(WEXITSTATUS(status) == EXIT_FAILURE);
    CHECK(dir.read("err").find("do not belong together") != std::string::npos);
    CHECK(dir.read("err").find("index --force") != std::string::npos);
}

TEST_CASE("disk BNS restore re-checks the restored metadata against the loaded FM rows")
{
    // 200 bases: l_pac 200, FM rows 2 * 200 + 1. The restore is given the rows of
    // a different FM index, as if .ann was replaced after the FM arrays loaded.
    bwa_tests::TempIndex fixture({{"guard", std::string(200, 'A')}});
    bwa_tests::TempDir dir("restore_recheck");
    const std::string err = dir.file("err");
    for (int which : {BWA_IDX_BNS, BWA_IDX_ALL}) {
        CAPTURE(which);
        const int status = bwa_tests::run_in_child([&] {
            indexEle idx;
            idx.bwa_idx_load_ele(fixture.prefix().c_str(), which, 1, 2 * 201 + 1);
        }, err.c_str());
        REQUIRE(WIFEXITED(status));
        CHECK(WEXITSTATUS(status) == EXIT_FAILURE);
        CHECK(dir.read("err").find("do not belong together") != std::string::npos);
    }
    // The matching row count, and the unchecked default, both load.
    for (int64_t rows : {int64_t(2 * 200 + 1), int64_t(-1)}) {
        CAPTURE(rows);
        const int status = bwa_tests::run_in_child([&] {
            indexEle idx;
            idx.bwa_idx_load_ele(fixture.prefix().c_str(), BWA_IDX_ALL, 1, rows);
            _exit(idx.idx->bns->l_pac == 200 && idx.idx->pac != NULL ? 0 : 3);
        }, err.c_str());
        REQUIRE(WIFEXITED(status));
        CHECK(WEXITSTATUS(status) == 0);
    }
}

TEST_SUITE_END();
