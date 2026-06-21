// Unit tests for src/meth_xm.cpp — the Bismark XM:Z walk. Exercises
// CpG/CHG/CHH context, mismatches, indels, soft clips, unknown context,
// out-of-bounds context, and both top/bottom strand orientations.
//
// Each test builds a tiny in-memory ORIGINAL (un-converted) bns + 2-bit pac
// via OrigRefFixture (the same handles the D3 runtime loads from disk), then
// drives meth_build_xm against it. Same code path the production runtime uses.

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <string>

#include "doctest/doctest.h"

#include "meth_xm.h"
#include "meth_orig_ref_fixture.h"

namespace {

inline uint32_t bam_cigar_pack(uint32_t len, uint32_t op) {
    return (len << 4) | (op & 0xf);
}

}  // namespace

TEST_CASE("meth_build_xm: is_top_strand=1 simple matches; CpG/CHG/CHH; methylated/unmethylated") {
    // ref: A C G T C A G T C A T A   (positions 0..11)
    //   pos 1: C; ref[2]=G            -> CpG
    //   pos 4: C; ref[5]=A; ref[6]=G  -> CHG
    //   pos 8: C; ref[9]=A; ref[10]=T -> CHH
    meth_test::OrigRefFixture f("ACGTCAGTCATA");

    // pos 1 methylated (read=C), pos 4 unmethylated (read=T), pos 8 methylated (read=C).
    const std::string read = "ACGTTAGTCATA";
    uint32_t cigar[] = { bam_cigar_pack(12, /*M*/0) };
    char *xm = meth_build_xm(f.orig_bns, f.orig_pac, f.real_tid, /*pos=*/0, /*is_top_strand=*/1,
                             cigar, 1, read.c_str(), (int)read.size());
    REQUIRE(xm != nullptr);
    CHECK(std::string(xm) == ".Z..x...H...");
}

TEST_CASE("meth_build_xm: is_top_strand=0 mirrors (G-marker, upstream context)") {
    // forward = T A C T C G A T A T   (positions 0..9)
    //   pos 5 = G (bottom-strand C). Bottom-strand context-1 = forward[4] = C
    //   (bottom 'G' downstream) -> CpG.
    meth_test::OrigRefFixture f("TACTCGATAT");

    // perfect match -> read=G at pos 5 = methylated CpG -> Z
    const std::string read_meth = "TACTCGATAT";
    uint32_t cigar[] = { bam_cigar_pack(10, /*M*/0) };
    char *xm = meth_build_xm(f.orig_bns, f.orig_pac, f.real_tid, /*pos=*/0, /*is_top_strand=*/0,
                             cigar, 1, read_meth.c_str(), (int)read_meth.size());
    REQUIRE(xm != nullptr);
    CHECK(std::string(xm) == ".....Z....");

    // unmethylated -> read=A at pos 5 -> z
    const std::string read_unmeth = "TACTCAATAT";
    char *xm2 = meth_build_xm(f.orig_bns, f.orig_pac, f.real_tid, 0, /*is_top_strand=*/0,
                              cigar, 1, read_unmeth.c_str(), (int)read_unmeth.size());
    REQUIRE(xm2 != nullptr);
    CHECK(std::string(xm2) == ".....z....");
}

TEST_CASE("meth_build_xm: insertion emits '.' per inserted base; deletion no emit") {
    meth_test::OrigRefFixture f("ACGTACGT");

    // 2M, 2I (inserted "GG"), 6M.
    const std::string read = "ACGGGTACGT";
    uint32_t cigar[] = {
        bam_cigar_pack(2, /*M*/0),
        bam_cigar_pack(2, /*I*/1),
        bam_cigar_pack(6, /*M*/0),
    };
    char *xm = meth_build_xm(f.orig_bns, f.orig_pac, f.real_tid, /*pos=*/0, /*is_top_strand=*/1,
                             cigar, 3, read.c_str(), (int)read.size());
    REQUIRE(xm != nullptr);
    CHECK(std::string(xm) == ".Z.....Z..");

    // Deletion case.
    const std::string read_d = "ACACGT";
    uint32_t cigar_d[] = {
        bam_cigar_pack(2, /*M*/0),
        bam_cigar_pack(2, /*D*/2),
        bam_cigar_pack(4, /*M*/0),
    };
    char *xm_d = meth_build_xm(f.orig_bns, f.orig_pac, f.real_tid, 0, /*is_top_strand=*/1,
                               cigar_d, 3, read_d.c_str(), (int)read_d.size());
    REQUIRE(xm_d != nullptr);
    CHECK(std::string(xm_d) == ".Z.Z..");
}

TEST_CASE("meth_build_xm: soft-clip emits '.' per clipped base") {
    meth_test::OrigRefFixture f("ACGT");
    const std::string read = "NNACGT";
    uint32_t cigar[] = {
        bam_cigar_pack(2, /*S*/4),
        bam_cigar_pack(4, /*M*/0),
    };
    char *xm = meth_build_xm(f.orig_bns, f.orig_pac, f.real_tid, /*pos=*/0, /*is_top_strand=*/1,
                             cigar, 2, read.c_str(), (int)read.size());
    REQUIRE(xm != nullptr);
    CHECK(std::string(xm) == "...Z..");
}

TEST_CASE("meth_build_xm: read base != C and != T at ref-C emits '.'") {
    meth_test::OrigRefFixture f("ACGT");
    // ref pos 1 = C (CpG); read at pos 1 = 'A' (mismatch, neither C nor T).
    const std::string read = "AAGT";
    uint32_t cigar[] = { bam_cigar_pack(4, /*M*/0) };
    char *xm = meth_build_xm(f.orig_bns, f.orig_pac, f.real_tid, 0, /*is_top_strand=*/1,
                             cigar, 1, read.c_str(), (int)read.size());
    REQUIRE(xm != nullptr);
    CHECK(std::string(xm) == "....");
}

TEST_CASE("meth_build_xm: ref N at context emits 'u'/'U' (unknown context)") {
    meth_test::OrigRefFixture f("ACNTA");
    const std::string read = "ACNTA";
    uint32_t cigar[] = { bam_cigar_pack(5, /*M*/0) };
    char *xm = meth_build_xm(f.orig_bns, f.orig_pac, f.real_tid, 0, /*is_top_strand=*/1,
                             cigar, 1, read.c_str(), (int)read.size());
    REQUIRE(xm != nullptr);
    // pos 1: C, ref[2]=N -> unknown, read=C meth -> 'U'
    CHECK(std::string(xm) == ".U...");
}

TEST_CASE("meth_build_xm: out-of-bounds context (read aligned at end of contig)") {
    // ref: T A C  (3 bp). Read aligns at pos 0, len 3. ref[3] OOB -> 'N'.
    meth_test::OrigRefFixture f("TAC");
    const std::string read = "TAC";
    uint32_t cigar[] = { bam_cigar_pack(3, /*M*/0) };
    char *xm = meth_build_xm(f.orig_bns, f.orig_pac, f.real_tid, 0, /*is_top_strand=*/1,
                             cigar, 1, read.c_str(), (int)read.size());
    REQUIRE(xm != nullptr);
    CHECK(std::string(xm) == "..U");
}
