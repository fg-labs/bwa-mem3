// emit_0123 must write a fresh `.0123` file: a symlink or hard link planted at
// that path (after the index command's output preflight) must not redirect the
// write into the file the link names, such as the input FASTA.

#include "index_prelude.h"
#include "packed_text.h"
#include "temp_dir.h"

#include <cstdint>
#include <string>
#include <sys/stat.h>
#include <unistd.h>

#include "doctest/doctest.h"

TEST_SUITE_BEGIN("unit/index_prelude");

namespace {

constexpr int64_t kBases = 13;

// Writes a packed text of kBases bases (base i = i mod 4) as dir/name and returns its path.
std::string write_pac(const bwa_tests::TempDir &dir, const std::string &name)
{
    std::string pac((kBases + 3) / 4, '\0');
    for (int64_t l = 0; l < kBases; ++l)
        pac[(size_t)(l >> 2)] = (char)((uint8_t)pac[(size_t)(l >> 2)] | ((l & 3) << ((~l & 3) << 1)));
    dir.write(name, pac);
    return dir.file(name);
}

std::string expected_0123()
{
    std::string text;
    for (int64_t l = 0; l < kBases; ++l) text.push_back((char)(l & 3));
    return text;
}

} // namespace

TEST_CASE("emit_0123: a planted symlink is replaced, not followed")
{
    bwa_tests::TempDir dir;
    const std::string victim_text = ">victim\nACGTACGTACGT\n";
    dir.write("victim.fa", victim_text);
    PackedText pac(write_pac(dir, "ref.pac"), kBases);
    REQUIRE(symlink(dir.file("victim.fa").c_str(), dir.file("ref.0123").c_str()) == 0);

    REQUIRE(emit_0123(pac, dir.file("ref").c_str(), 1) == 0);

    CHECK(dir.read("victim.fa") == victim_text);
    struct stat st;
    REQUIRE(lstat(dir.file("ref.0123").c_str(), &st) == 0);
    CHECK(S_ISREG(st.st_mode));
    CHECK(dir.read("ref.0123") == expected_0123());
}

TEST_CASE("emit_0123: a planted hard link is replaced, not truncated through")
{
    bwa_tests::TempDir dir;
    const std::string victim_text = ">victim\nACGTACGTACGT\n";
    dir.write("victim.fa", victim_text);
    PackedText pac(write_pac(dir, "ref.pac"), kBases);
    REQUIRE(link(dir.file("victim.fa").c_str(), dir.file("ref.0123").c_str()) == 0);

    REQUIRE(emit_0123(pac, dir.file("ref").c_str(), 1) == 0);

    CHECK(dir.read("victim.fa") == victim_text);
    struct stat st;
    REQUIRE(lstat(dir.file("ref.0123").c_str(), &st) == 0);
    CHECK(st.st_nlink == 1);
    CHECK(dir.read("ref.0123") == expected_0123());
}

TEST_CASE("emit_0123: an existing regular .0123 is replaced with the new contents")
{
    bwa_tests::TempDir dir;
    dir.write("ref.0123", std::string(64, '\7'));
    PackedText pac(write_pac(dir, "ref.pac"), kBases);

    REQUIRE(emit_0123(pac, dir.file("ref").c_str(), 2) == 0);

    CHECK(dir.read("ref.0123") == expected_0123());
}

TEST_SUITE_END();
