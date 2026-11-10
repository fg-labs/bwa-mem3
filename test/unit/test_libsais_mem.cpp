// In-memory doubled-pac packer (src/libsais_build.cpp): must produce exactly the bytes
// write_doubled_pac writes to disk.

#include "doctest/doctest.h"

#include "../framework/temp_dir.h"
#include "../../src/libsais_build.h"

#include <cstdint>
#include <cstdio>
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
