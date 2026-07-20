// Unit tests for the CHEMISTRY dimension of src/meth_xm.cpp — the
// methylated/unmethylated polarity of the Bismark XM:Z call string.
//
// Bisulfite/em-seq and TAPS produce the SAME aligner-visible base change
// (C->T on the converted strand) but assign it the OPPOSITE meaning:
//
//   em-seq : UNmethylated C converts.  read T at ref C => UNmethylated (z)
//                                      read C at ref C => methylated   (Z)
//   TAPS   : METHYLATED  C converts.   read T at ref C => methylated   (Z)
//                                      read C at ref C => UNmethylated (z)
//
// The context classification (CpG / CHG / CHH, and which base is the
// "C of interest" per strand) is chemistry-INDEPENDENT: it is a property of
// the reference, not of the chemistry. Only the meth/unmeth decision flips.
//
// These tests pin that contract for both chemistries. The em-seq cases
// duplicate the expectations in test_meth_xm_build.cpp on purpose: they are
// the control that proves a TAPS regression is a real polarity flip and not
// a broken context walk.

#include <cstdint>
#include <string>

#include "doctest/doctest.h"

#include "meth_xm.h"
#include "meth_orig_ref_fixture.h"

namespace {

inline uint32_t bam_cigar_pack(uint32_t len, uint32_t op) {
    return (len << 4) | (op & 0xf);
}

}  // namespace

TEST_CASE("meth_build_xm chemistry: top strand, TAPS inverts em-seq meth state") {
    // ref: A C G T C A G T C A T A   (positions 0..11)
    //   pos 1: C; ref[2]=G            -> CpG
    //   pos 4: C; ref[5]=A; ref[6]=G  -> CHG
    //   pos 8: C; ref[9]=A; ref[10]=T -> CHH
    meth_test::OrigRefFixture f("ACGTCAGTCATA");

    // read keeps C at pos 1 and 8, converts C->T at pos 4.
    const std::string read = "ACGTTAGTCATA";
    uint32_t cigar[] = { bam_cigar_pack(12, /*M*/0) };

    SUBCASE("em-seq: unconverted C = methylated") {
        char *xm = meth_build_xm(f.orig_bns, f.orig_pac, f.real_tid, /*pos=*/0,
                                 /*is_top_strand=*/1, cigar, 1,
                                 read.c_str(), (int)read.size(),
                                 METH_CHEM_EMSEQ);
        REQUIRE(xm != nullptr);
        CHECK(std::string(xm) == ".Z..x...H...");
    }

    SUBCASE("TAPS: unconverted C = UNmethylated (exact case inversion)") {
        char *xm = meth_build_xm(f.orig_bns, f.orig_pac, f.real_tid, /*pos=*/0,
                                 /*is_top_strand=*/1, cigar, 1,
                                 read.c_str(), (int)read.size(),
                                 METH_CHEM_TAPS);
        REQUIRE(xm != nullptr);
        // Same positions, same context letters, opposite case.
        CHECK(std::string(xm) == ".z..X...h...");
    }
}

TEST_CASE("meth_build_xm chemistry: bottom strand, TAPS inverts em-seq meth state") {
    // forward = T A C T C G A T A T (positions 0..9)
    //   pos 5 = G (bottom-strand C); bottom context-1 = forward[4] = C -> CpG
    meth_test::OrigRefFixture f("TACTCGATAT");
    uint32_t cigar[] = { bam_cigar_pack(10, /*M*/0) };

    const std::string read_g = "TACTCGATAT";  // G retained (unconverted)
    const std::string read_a = "TACTCAATAT";  // G->A (converted)

    SUBCASE("em-seq: retained G = methylated, converted A = unmethylated") {
        char *xm_g = meth_build_xm(f.orig_bns, f.orig_pac, f.real_tid, 0,
                                   /*is_top_strand=*/0, cigar, 1,
                                   read_g.c_str(), (int)read_g.size(),
                                   METH_CHEM_EMSEQ);
        REQUIRE(xm_g != nullptr);
        CHECK(std::string(xm_g) == ".....Z....");

        char *xm_a = meth_build_xm(f.orig_bns, f.orig_pac, f.real_tid, 0,
                                   /*is_top_strand=*/0, cigar, 1,
                                   read_a.c_str(), (int)read_a.size(),
                                   METH_CHEM_EMSEQ);
        REQUIRE(xm_a != nullptr);
        CHECK(std::string(xm_a) == ".....z....");
    }

    SUBCASE("TAPS: retained G = UNmethylated, converted A = methylated") {
        char *xm_g = meth_build_xm(f.orig_bns, f.orig_pac, f.real_tid, 0,
                                   /*is_top_strand=*/0, cigar, 1,
                                   read_g.c_str(), (int)read_g.size(),
                                   METH_CHEM_TAPS);
        REQUIRE(xm_g != nullptr);
        CHECK(std::string(xm_g) == ".....z....");

        char *xm_a = meth_build_xm(f.orig_bns, f.orig_pac, f.real_tid, 0,
                                   /*is_top_strand=*/0, cigar, 1,
                                   read_a.c_str(), (int)read_a.size(),
                                   METH_CHEM_TAPS);
        REQUIRE(xm_a != nullptr);
        CHECK(std::string(xm_a) == ".....Z....");
    }
}

TEST_CASE("meth_build_xm chemistry: non-C/T read base emits '.' under BOTH chemistries") {
    // A base that is neither the retained nor the converted base is not a
    // recognised methylation event in either chemistry -- it is a mismatch.
    // ref pos 1 = C, ref[2] = G -> CpG; read carries 'A' there.
    meth_test::OrigRefFixture f("ACGT");
    const std::string read = "AAGT";
    uint32_t cigar[] = { bam_cigar_pack(4, /*M*/0) };

    char *xm_e = meth_build_xm(f.orig_bns, f.orig_pac, f.real_tid, 0, 1, cigar, 1,
                               read.c_str(), (int)read.size(), METH_CHEM_EMSEQ);
    REQUIRE(xm_e != nullptr);
    CHECK(std::string(xm_e) == "....");

    char *xm_t = meth_build_xm(f.orig_bns, f.orig_pac, f.real_tid, 0, 1, cigar, 1,
                               read.c_str(), (int)read.size(), METH_CHEM_TAPS);
    REQUIRE(xm_t != nullptr);
    CHECK(std::string(xm_t) == "....");
}

TEST_CASE("meth_build_xm chemistry: unknown context stays 'u'/'U' with chemistry-correct case") {
    // ref[2] = N -> context unknown at the ref-C in pos 1.
    // em-seq: retained C => methylated => 'U'.  TAPS: retained C => 'u'.
    meth_test::OrigRefFixture f("ACNTA");
    const std::string read = "ACNTA";
    uint32_t cigar[] = { bam_cigar_pack(5, /*M*/0) };

    char *xm_e = meth_build_xm(f.orig_bns, f.orig_pac, f.real_tid, 0, 1, cigar, 1,
                               read.c_str(), (int)read.size(), METH_CHEM_EMSEQ);
    REQUIRE(xm_e != nullptr);
    CHECK(std::string(xm_e) == ".U...");

    char *xm_t = meth_build_xm(f.orig_bns, f.orig_pac, f.real_tid, 0, 1, cigar, 1,
                               read.c_str(), (int)read.size(), METH_CHEM_TAPS);
    REQUIRE(xm_t != nullptr);
    CHECK(std::string(xm_t) == ".u...");
}
