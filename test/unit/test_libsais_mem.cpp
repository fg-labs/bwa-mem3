// In-memory doubled-pac packer (src/libsais_build.cpp): must produce exactly the bytes
// write_doubled_pac writes to disk.

#include "doctest/doctest.h"

#include "../framework/temp_dir.h"
#include "../../src/fm_index_writer.h"
#include "../../src/bntseq.h"
#include "../../src/bwa.h"
#include "../../src/libsais_build.h"

#include <cstdint>
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

namespace {

std::vector<uint8_t> read_file(const std::string &path)
{
    std::vector<uint8_t> v;
    FILE *f = std::fopen(path.c_str(), "rb");
    REQUIRE(f != nullptr);
    uint8_t buf[4096];
    size_t n;
    while ((n = std::fread(buf, 1, sizeof buf, f)) > 0) v.insert(v.end(), buf, buf + n);
    std::fclose(f);
    return v;
}

// Locate a test/fixtures file whatever the working directory is (the tests run from
// the repo root or from test/).
std::string fixture(const char *name)
{
    const std::string src = __FILE__;   // .../test/unit/test_libsais_mem.cpp
    const std::string dir = src.substr(0, src.rfind('/'));
    for (const std::string &base : {dir + "/../fixtures/", std::string("test/fixtures/"),
                                    std::string("fixtures/"), std::string("../fixtures/")}) {
        FILE *f = std::fopen((base + name).c_str(), "rb");
        if (f) { std::fclose(f); return base + name; }
    }
    FAIL("fixture not found: " << name);
    return "";
}

} // namespace

TEST_CASE("libsais_mem: pack_doubled equals write_doubled_pac")
{
    bwa_tests::TempDir tmp("libsais_mem");
    const int64_t sizes[] = {1, 3, 4, 5, 1003};
    for (int64_t l_pac : sizes) {
        const size_t fwd_bytes = (size_t)((l_pac + 3) / 4);
        std::vector<uint8_t> fwd(fwd_bytes, 0);
        uint32_t st = 777u + (uint32_t)l_pac;
        for (int64_t l = 0; l < l_pac; ++l) {
            st = st * 1664525u + 1013904223u;
            fwd[(size_t)(l >> 2)] |= (uint8_t)(((st >> 24) & 3) << ((~l & 3) << 1));
        }
        const std::string in = tmp.path() + "/fwd.pac", out = tmp.path() + "/doubled.pac";
        FILE *f = std::fopen(in.c_str(), "wb");
        REQUIRE(f != nullptr);
        REQUIRE(std::fwrite(fwd.data(), 1, fwd_bytes, f) == fwd_bytes);
        std::fclose(f);
        for (int threads : {1, 4}) {
            CAPTURE(l_pac);
            CAPTURE(threads);
            write_doubled_pac(in.c_str(), l_pac, out.c_str(), threads);
            std::vector<uint8_t> disk = read_file(out);
            const size_t nbytes = (size_t)((2 * l_pac + 3) / 4);
            REQUIRE(disk.size() >= nbytes);
            PackedText view(fwd.data(), l_pac);
            std::unique_ptr<uint8_t[]> mem = pack_doubled(view, l_pac, threads);
            REQUIRE(mem != nullptr);
            for (size_t i = 0; i < nbytes; ++i) CHECK(mem[i] == disk[i]);
        }
    }
}

TEST_CASE("libsais_mem: MemSink holds the same bytes the file writer writes")
{
    bwa_tests::TempDir tmp("libsais_mem_sink");
    const int64_t N = 4098;   // doubled-text length; the index holds N + 1 = 4099 rows
    std::vector<uint8_t> buf((size_t)N + 1, 0);   // libsais alphabet: 1..4 = ACGT, buf[N] = 0
    uint32_t st = 12345u;
    int64_t count[5] = {0, 0, 0, 0, 0};
    int64_t freq[4] = {0, 0, 0, 0};
    for (int64_t i = 0; i < N; ++i) {
        st = st * 1664525u + 1013904223u;
        const uint8_t b = (uint8_t)((st >> 24) & 3);
        buf[(size_t)i] = (uint8_t)(b + 1);
        ++freq[b];
    }
    for (int c = 1; c <= 4; ++c) count[c] = count[c - 1] + freq[c - 1];

    std::vector<int32_t> sa((size_t)N + 1);
    for (int64_t i = 0; i <= N; ++i) sa[(size_t)i] = (int32_t)i;
    std::sort(sa.begin(), sa.end(), [&](int32_t a, int32_t b) {
        while (true) {   // the unique 0 terminator at buf[N] ends every comparison
            if (buf[(size_t)a] != buf[(size_t)b]) return buf[(size_t)a] < buf[(size_t)b];
            ++a; ++b;
        }
    });

    const int sa_compx = 3;
    for (int threads : {1, 4}) {
        CAPTURE(threads);
        const std::string path = tmp.path() + "/idx.bwt.2bit.64";
        int64_t sent_file = -1;
        write_fm_index_streaming(path.c_str(), buf.data(), sa.data(), false, N, count, &sent_file,
                                 sa_compx, threads);
        const std::vector<uint8_t> file = read_file(path);

        FmArrays arr;
        MemSink sink(N + 1, sa_compx, &arr);
        int64_t sent_mem = -1;
        write_fm_index_streaming(sink, buf.data(), sa.data(), false, N, count, &sent_mem,
                                 sa_compx, threads);

        const Bwt2bitLayout L = fmi_bwt2bit_layout(N + 1, sa_compx);
        REQUIRE(file.size() == (size_t)L.total_bytes);
        CHECK(sent_mem == sent_file);
        CHECK(arr.ref_seq_len == N + 1);
        CHECK(std::memcmp(arr.count, file.data() + 8, sizeof arr.count) == 0);
        CHECK(arr.sentinel_index == sent_file);
        CHECK(arr.sa_compx == sa_compx);
        CHECK(std::memcmp(arr.cp_occ, file.data() + L.off_cp_occ,
                          (size_t)L.cp_occ_count * sizeof(CP_OCC)) == 0);
        CHECK(std::memcmp(arr.sa_ms_byte, file.data() + L.off_ms_byte,
                          (size_t)L.sa_sample_count) == 0);
        CHECK(std::memcmp(arr.sa_ls_word, file.data() + L.off_ls_word,
                          (size_t)L.sa_sample_count * sizeof(uint32_t)) == 0);
    }
}

TEST_CASE("libsais_mem: memory build equals the index file")
{
    for (const char *fa_name : {"phix.fa", "synthetic_1mb.fa"}) {
        const std::string fa = fixture(fa_name);
        bwa_tests::TempDir tmp("libsais_mem_build");
        const std::string prefix = tmp.path() + "/idx";
        REQUIRE(bwa_idx_build(fa.c_str(), prefix.c_str(), 0, 3) == 0);
        bntseq_t *bns = bns_restore(prefix.c_str());
        REQUIRE(bns != nullptr);
        const int64_t l_pac = bns->l_pac;
        bns_destroy(bns);
        const std::vector<uint8_t> pac = read_file(prefix + ".pac");
        REQUIRE(pac.size() >= (size_t)((l_pac + 3) / 4));
        const std::vector<uint8_t> file = read_file(prefix + ".bwt.2bit.64");
        const int sa_compx = 3;
        const Bwt2bitLayout L = fmi_bwt2bit_layout(2 * l_pac + 1, sa_compx);
        REQUIRE(file.size() == (size_t)L.total_bytes);
        for (int threads : {1, 8}) {
            CAPTURE(fa_name);
            CAPTURE(threads);
            LibsaisBuildOpts opts;
            opts.num_threads = threads;
            opts.sa_compx = sa_compx;
            FmArrays arr;
            REQUIRE(libsais_build_fm_memory(pac.data(), l_pac, opts, &arr) == 0);
            CHECK(arr.ref_seq_len == 2 * l_pac + 1);
            CHECK(arr.sa_compx == sa_compx);
            CHECK(std::memcmp(&arr.ref_seq_len, file.data(), 8) == 0);
            CHECK(std::memcmp(arr.count, file.data() + 8, sizeof arr.count) == 0);
            int64_t sent;
            std::memcpy(&sent, file.data() + L.off_sentinel, sizeof sent);
            CHECK(arr.sentinel_index == sent);
            CHECK(std::memcmp(arr.cp_occ, file.data() + L.off_cp_occ,
                              (size_t)L.cp_occ_count * sizeof(CP_OCC)) == 0);
            CHECK(std::memcmp(arr.sa_ms_byte, file.data() + L.off_ms_byte,
                              (size_t)L.sa_sample_count) == 0);
            CHECK(std::memcmp(arr.sa_ls_word, file.data() + L.off_ls_word,
                              (size_t)L.sa_sample_count * sizeof(uint32_t)) == 0);
        }
    }
}

TEST_CASE("libsais_mem: a memory budget refusal allocates nothing")
{
    const int64_t l_pac = 1000;
    std::vector<uint8_t> pac((size_t)(l_pac + 3) / 4, 0x1b);
    LibsaisBuildOpts opts;
    opts.max_memory_bytes = 1;   // far below the estimate
    FmArrays arr;
    CHECK(libsais_build_fm_memory(pac.data(), l_pac, opts, &arr) == 3);
    CHECK(arr.cp_occ == nullptr);
    CHECK(arr.sa_ms_byte == nullptr);
    CHECK(arr.sa_ls_word == nullptr);
}
