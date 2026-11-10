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
