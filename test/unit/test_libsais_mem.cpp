// In-memory doubled-pac packer (src/libsais_build.cpp): must produce exactly the bytes
// write_doubled_pac writes to disk.

#include "doctest/doctest.h"

#include "../framework/temp_dir.h"
#include "../../src/fm_index_writer.h"
#include "../../src/libsais_build.h"

#include <cstdint>
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

namespace {

// 2-bit pack bases (0-3), first base in the high bits, unused low bits of a partial final byte
// zero: the .pac byte layout, built here one base at a time as an oracle for pack_doubled.
std::vector<uint8_t> pack_bases(const std::vector<uint8_t> &bases)
{
    std::vector<uint8_t> p((bases.size() + 3) / 4, 0);
    for (size_t i = 0; i < bases.size(); ++i) p[i / 4] |= (uint8_t)(bases[i] << (2 * (3 - i % 4)));
    return p;
}

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

} // namespace

TEST_CASE("libsais_mem: pack_doubled packs forward then reverse complement"
          * doctest::test_suite("unit/libsais_mem"))
{
    // ACG -> doubled ACG CGT = 0 1 2 1 2 3 -> 00011001 1011(0000).
    const std::vector<uint8_t> fwd = pack_bases({0, 1, 2});
    REQUIRE(fwd.size() == 1);
    REQUIRE(fwd[0] == 0x18);
    PackedText view(fwd.data(), 3);
    std::unique_ptr<uint8_t[]> mem = pack_doubled(view, 3, 1);
    REQUIRE(mem != nullptr);
    CHECK(mem[0] == 0x19);
    CHECK(mem[1] == 0xB0);
}

TEST_CASE("libsais_mem: pack_doubled and write_doubled_pac equal an independent packing"
          * doctest::test_suite("unit/libsais_mem"))
{
    bwa_tests::TempDir tmp("libsais_mem");
    // Odd lengths end the doubled text in a partial byte; even ones fill the last byte and
    // write_doubled_pac then adds a pad byte before the count byte.
    const int64_t sizes[] = {1, 3, 4, 5, 1003, 1004};
    for (int64_t l_pac : sizes) {
        std::vector<uint8_t> bases((size_t)l_pac);
        uint32_t st = 777u + (uint32_t)l_pac;
        for (int64_t l = 0; l < l_pac; ++l) {
            st = st * 1664525u + 1013904223u;
            bases[(size_t)l] = (uint8_t)((st >> 24) & 3);
        }
        std::vector<uint8_t> doubled(bases);
        for (int64_t l = l_pac - 1; l >= 0; --l) doubled.push_back((uint8_t)(3 - bases[(size_t)l]));
        const std::vector<uint8_t> fwd = pack_bases(bases);
        const std::vector<uint8_t> want = pack_bases(doubled);
        const size_t nbytes = (size_t)((2 * l_pac + 3) / 4);
        REQUIRE(want.size() == nbytes);

        const std::string in = tmp.path() + "/fwd.pac", out = tmp.path() + "/doubled.pac";
        FILE *f = std::fopen(in.c_str(), "wb");
        REQUIRE(f != nullptr);
        REQUIRE(std::fwrite(fwd.data(), 1, fwd.size(), f) == fwd.size());
        std::fclose(f);
        for (int threads : {1, 4}) {
            CAPTURE(l_pac);
            CAPTURE(threads);
            PackedText view(fwd.data(), l_pac);
            std::unique_ptr<uint8_t[]> mem = pack_doubled(view, l_pac, threads);
            REQUIRE(mem != nullptr);
            CHECK(std::vector<uint8_t>(mem.get(), mem.get() + nbytes) == want);

            write_doubled_pac(in.c_str(), l_pac, out.c_str(), threads);
            const std::vector<uint8_t> disk = read_file(out);
            const bool pad = (2 * l_pac) % 4 == 0;
            REQUIRE(disk.size() == nbytes + (pad ? 2 : 1));
            CHECK(std::vector<uint8_t>(disk.data(), disk.data() + nbytes) == want);
            if (pad) CHECK(disk[nbytes] == 0);
            CHECK(disk.back() == (uint8_t)((2 * l_pac) % 4));
        }
    }
}

TEST_CASE("libsais_mem: MemSink holds the same bytes the file writer writes"
          * doctest::test_suite("unit/libsais_mem"))
{
    bwa_tests::TempDir tmp("libsais_mem_sink");
    // N is the doubled-text length; the index holds N + 1 rows. 4099 rows fill no block or
    // sample stride exactly; 4096 rows do, so the layout's last cp_occ entry and last SA sample
    // are reserved but never written, and must read back zero as they do from the file.
    for (const int64_t N : {(int64_t)4098, (int64_t)4095}) {
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
            if (a == b) return false;   // irreflexive, as std::sort requires
            while (true) {   // the unique 0 terminator at buf[N] ends every other comparison
                if (buf[(size_t)a] != buf[(size_t)b]) return buf[(size_t)a] < buf[(size_t)b];
                ++a; ++b;
            }
        });

        for (const int sa_compx : {0, 3})
            for (int threads : {1, 4}) {
                CAPTURE(N);
                CAPTURE(sa_compx);
                CAPTURE(threads);
                const std::string path = tmp.path() + "/idx.bwt.2bit.64";
                int64_t sent_file = -1;
                write_fm_index_streaming(path.c_str(), buf.data(), sa.data(), false, N, count,
                                         &sent_file, sa_compx, threads);
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
                if ((N + 1) % ((int64_t)1 << sa_compx) == 0) {
                    CHECK(arr.sa_ms_byte[L.sa_sample_count - 1] == 0);
                    CHECK(arr.sa_ls_word[L.sa_sample_count - 1] == 0);
                }
            }
    }
}
