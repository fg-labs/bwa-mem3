// test/unit/test_mem_io_sam_record.cpp
//
// The SAM text record (mem_aln2sam) and the BAM record (mem_aln_to_bam) of
// hand-built alignments, field by field, against lines written out by hand.
//
// mem_aln2sam formats one mem_aln_t into a SAM line with unchecked writes into
// a buffer it sizes once, from an estimate of every field's length, so an
// estimate that falls short writes past the allocation. Every line here is
// also checked to fit the buffer's recorded capacity, and the long-record
// cases make each estimated field in turn long enough that dropping its term
// from the estimate would overrun it.
//
// The common records (forward and reverse single-end hits, a gapped CIGAR,
// proper pairs, --compat's tag set) are already pinned end to end by
// test/regression/mem_cli_basics.sh, chr22_parity.sh and
// compat_byte_identical.sh; the cases here are the ones those runs do not
// produce (FASTA input, unmapped and secondary records, split hits with their
// clips and SA, -Y, -M, ALT hits, mates on another contig or unmapped, TLEN
// edges, every optional field at once) and the BAM writer's copy of each. The
// expected lines follow the SAM specification and bwa's record model
// (bwamem.c mem_aln2sam, which bwa-mem3 reproduces):
//
//   - FLAG ORs in 0x1 with a mate, 0x4 when unmapped, 0x8 when the mate is
//     unmapped, 0x10/0x20 for the read's/mate's strand; the internal -M bit
//     0x10000 prints as 0x100.
//   - An unmapped read with a mapped mate takes the mate's contig, position
//     and strand and prints CIGAR "*"; a mapped read with an unmapped mate
//     gives the mate its own position (RNEXT "=", PNEXT = POS, TLEN 0).
//   - TLEN of a pair on one contig runs from this read's 5' end to its mate's
//     (the last mapped base of a reverse-strand hit, the first of a forward
//     one): p1 - p0 + 1 when the mate's 5' end is further right, p1 - p0 - 1
//     when it is further left, 0 when they meet. For an inward
//     forward/reverse pair that is the template span, positive on the
//     leftmost read.
//   - A reverse-strand record prints the reverse complement of SEQ and the
//     reversed QUAL; a secondary (0x100) record prints "*" for both.
//   - A supplementary record (which > 0) prints its clips as H and drops the
//     clipped bases from SEQ/QUAL, unless -Y (MEM_F_SOFTCLIP) is set or the hit
//     is on an ALT contig.
//   - Optional fields, in order: NM MD (aligned only), MC (mate CIGAR), MQ (mate
//     MAPQ), AS (score >= 0), XS (sub >= 0), RG (-R), SA (the other
//     non-secondary hits, as rname,pos,strand,CIGAR,MAPQ,NM; with S clips), pa
//     (score / ALT score, three decimals), XA, HN, the FASTQ comment (-C), XR
//     (-V: the contig's FASTA header comment with TABs as spaces).
//
// The BAM writer must carry the same record: each case is also converted with
// mem_aln_to_bam and printed back with htslib's sam_format1, which has to give
// the same line (htslib prints a float tag with %g, so pa:f:1.250 reads back as
// pa:f:1.25). Neither writer is used as the oracle for the other.

// htslib headers first: they share the KSTRING_H include guard with
// bwa-mem3's kstring.h, and the BAM writer is compiled against htslib's.
#include "htslib/kstring.h"
#include "htslib/sam.h"

#include "doctest/doctest.h"
#include "mem_record_fixture.h"

#include "bam_writer.h"
#include "bwa.h"
#include "bwamem.h"

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

namespace {

using bwa_tests::OptGuard;
typedef bwa_tests::RecordRef Ref;
typedef bwa_tests::RecordRead Read;
typedef bwa_tests::RecordAlns Alns;

// The unmapped record mem_reg2sam writes for a read with no alignment: the
// conversion of no region.
mem_aln_t unmapped(const mem_opt_t *opt, const Ref &ref, const Read &r) {
    return mem_reg2aln(opt, &ref.bns, NULL, r.s.l_seq, r.s.seq, NULL);
}

std::string sam_line(const mem_opt_t *opt, const Ref &ref, Read &r,
                     const std::vector<mem_aln_t> &list, int which, const mem_aln_t *mate) {
    __kstring_t str = {0, 0, NULL};
    mem_aln2sam(opt, &ref.bns, &str, &r.s, (int)list.size(), list.data(), which, mate);
    REQUIRE(str.s != NULL);
    // The writes are unchecked and do not grow m: a line (and its NUL) past
    // the capacity the estimate set has been written out of bounds.
    CHECK(str.l < str.m);
    REQUIRE(strlen(str.s) == str.l);  // NUL-terminated at the record's end
    std::string out(str.s, str.l);
    free(str.s);
    return out;
}

// The same record through the BAM writer, printed back as SAM text with a
// trailing newline (sam_format1 omits it).
std::string bam_line(const mem_opt_t *opt, const Ref &ref, Read &r,
                     const std::vector<mem_aln_t> &list, int which, const mem_aln_t *mate) {
    const std::string htext = ref.header_text();
    sam_hdr_t *h = sam_hdr_parse(htext.size(), htext.c_str());
    REQUIRE(h != NULL);
    bam1_t *b = bam_writer_alloc();
    REQUIRE(b != NULL);
    const int rc = mem_aln_to_bam(b, opt, &ref.bns, &r.s, (int)list.size(), list.data(), which, mate);
    kstring_t ks = {0, 0, NULL};
    std::string out;
    if (rc == 0 && sam_format1(h, b, &ks) >= 0) out = std::string(ks.s, ks.l) + "\n";
    free(ks.s);
    bam_writer_free(b);
    sam_hdr_destroy(h);
    CHECK(rc == 0);
    return out;
}

// Both writers against one expected line (BAM against `bam_expected` when the
// two differ only in htslib's float formatting).
void check_both(const mem_opt_t *opt, const Ref &ref, Read &r, const std::vector<mem_aln_t> &list,
                int which, const mem_aln_t *mate, const std::string &expected,
                const std::string &bam_expected = std::string()) {
    CAPTURE(which);
    CHECK(sam_line(opt, ref, r, list, which, mate) == expected);
    CHECK(bam_line(opt, ref, r, list, which, mate)
          == (bam_expected.empty() ? expected : bam_expected));
}

// Saves and restores the -R read-group ID (a process global).
struct RgIdGuard {
    char saved[sizeof(bwa_rg_id)];
    RgIdGuard() { memcpy(saved, bwa_rg_id, sizeof(saved)); }
    ~RgIdGuard() { memcpy(bwa_rg_id, saved, sizeof(saved)); }
};


// Read r1: ACGTNACGTT / ABCDEFGHIJ. Its reverse complement is AACGTNACGT and
// its reversed quality JIHGFEDCBA.
const char *kBases = "ACGTNACGTT";
const char *kQual = "ABCDEFGHIJ";

}  // namespace

TEST_CASE("mem_io record: FASTA input and an absent XS"
          * doctest::test_suite("unit/mem_io")) {
    OptGuard g;
    Ref ref;
    Alns alns;
    Read r("r1", kBases, kQual);

    SUBCASE("no quality (FASTA input): QUAL is *, on either strand") {
        Read fa("r1", kBases, "");
        std::vector<mem_aln_t> f{alns.make(0, 99, false, "10M", "10", 0, 60, 10, 0)};
        check_both(g.opt, ref, fa, f, 0, NULL,
                   "r1\t0\tchrA\t100\t60\t10M\t*\t0\t0\tACGTNACGTT\t*\t"
                   "NM:i:0\tMD:Z:10\tAS:i:10\tXS:i:0\n");
        std::vector<mem_aln_t> rv{alns.make(0, 99, true, "10M", "10", 0, 60, 10, 0)};
        check_both(g.opt, ref, fa, rv, 0, NULL,
                   "r1\t16\tchrA\t100\t60\t10M\t*\t0\t0\tAACGTNACGT\t*\t"
                   "NM:i:0\tMD:Z:10\tAS:i:10\tXS:i:0\n");
    }
    SUBCASE("a negative sub score: XS absent") {
        std::vector<mem_aln_t> l{alns.make(0, 9, false, "3M1I2M2D4M", "5^GA1C2", 4, 12, 0, -1)};
        check_both(g.opt, ref, r, l, 0, NULL,
                   "r1\t0\tchrA\t10\t12\t3M1I2M2D4M\t*\t0\t0\tACGTNACGTT\tABCDEFGHIJ\t"
                   "NM:i:4\tMD:Z:5^GA1C2\tAS:i:0\n");
    }
}

TEST_CASE("mem_io record: unmapped and secondary records"
          * doctest::test_suite("unit/mem_io")) {
    OptGuard g;
    Ref ref;
    Alns alns;
    Read r("r1", kBases, kQual);

    SUBCASE("unmapped: FLAG 4, RNAME *, POS 0, CIGAR *, SEQ/QUAL as read, AS/XS 0") {
        std::vector<mem_aln_t> l{unmapped(g.opt, ref, r)};
        check_both(g.opt, ref, r, l, 0, NULL,
                   "r1\t4\t*\t0\t0\t*\t*\t0\t0\tACGTNACGTT\tABCDEFGHIJ\tAS:i:0\tXS:i:0\n");
    }
    SUBCASE("unmapped with -C: the comment follows AS and XS, as on a mapped record") {
        Read rc("r1", kBases, kQual, "BC:Z:ACGT\tXN:i:7");
        std::vector<mem_aln_t> l{unmapped(g.opt, ref, rc)};
        check_both(g.opt, ref, rc, l, 0, NULL,
                   "r1\t4\t*\t0\t0\t*\t*\t0\t0\tACGTNACGTT\tABCDEFGHIJ\tAS:i:0\tXS:i:0\t"
                   "BC:Z:ACGT\tXN:i:7\n");
    }
    SUBCASE("secondary: SEQ and QUAL are *, no SA even beside a primary") {
        // A reverse-strand secondary: 0x100 | 0x10 = 272.
        std::vector<mem_aln_t> l{alns.make(0, 99, false, "10M", "10", 0, 60, 10, 9),
                                 alns.make(1, 49, true, "10M", "4A5", 1, 0, 9, -1)};
        l[1].flag = 0x100;
        check_both(g.opt, ref, r, l, 1, NULL,
                   "r1\t272\tchrB\t50\t0\t10M\t*\t0\t0\t*\t*\tNM:i:1\tMD:Z:4A5\tAS:i:9\n");
        // The primary's SA lists only other non-secondary hits: none here.
        check_both(g.opt, ref, r, l, 0, NULL,
                   "r1\t0\tchrA\t100\t60\t10M\t*\t0\t0\tACGTNACGTT\tABCDEFGHIJ\t"
                   "NM:i:0\tMD:Z:10\tAS:i:10\tXS:i:9\n");
    }
}

TEST_CASE("mem_io record: split alignments with hard and soft clips and SA and -M"
          * doctest::test_suite("unit/mem_io")) {
    OptGuard g;
    Ref ref;
    Alns alns;
    Read r("r1", kBases, kQual);
    // Read bases 0-5 (ACGTNA) on chrA:100-105, bases 6-9 (CGTT) on chrB:50-53.
    std::vector<mem_aln_t> l{alns.make(0, 99, false, "6M4S", "6", 0, 60, 6, 0),
                             alns.make(1, 49, false, "6S4M", "4", 0, 20, 4, 0)};
    l[0].HN = 0;
    l[1].HN = 0;
    l[1].flag = 0x800;

    SUBCASE("forward supplementary: H clips, SEQ/QUAL of the aligned part only") {
        // The primary keeps S; each record's SA names the other hit with S.
        check_both(g.opt, ref, r, l, 0, NULL,
                   "r1\t0\tchrA\t100\t60\t6M4S\t*\t0\t0\tACGTNACGTT\tABCDEFGHIJ\t"
                   "NM:i:0\tMD:Z:6\tAS:i:6\tXS:i:0\tSA:Z:chrB,50,+,6S4M,20,0;\tHN:i:0\n");
        check_both(g.opt, ref, r, l, 1, NULL,
                   "r1\t2048\tchrB\t50\t20\t6H4M\t*\t0\t0\tCGTT\tGHIJ\t"
                   "NM:i:0\tMD:Z:4\tAS:i:4\tXS:i:0\tSA:Z:chrA,100,+,6M4S,60,0;\tHN:i:0\n");
    }
    SUBCASE("reverse supplementary: the clip is at the CIGAR's end") {
        // On the reverse strand the read's first 6 bases are the last 6 of the
        // reference-oriented alignment: 4M6S. SEQ is revcomp(CGTT) = AACG,
        // QUAL reverse(GHIJ) = JIHG.
        l[1] = alns.make(1, 49, true, "4M6S", "4", 0, 20, 4, 0);
        l[1].HN = 0;
        l[1].flag = 0x800;
        check_both(g.opt, ref, r, l, 1, NULL,
                   "r1\t2064\tchrB\t50\t20\t4M6H\t*\t0\t0\tAACG\tJIHG\t"
                   "NM:i:0\tMD:Z:4\tAS:i:4\tXS:i:0\tSA:Z:chrA,100,+,6M4S,60,0;\tHN:i:0\n");
        check_both(g.opt, ref, r, l, 0, NULL,
                   "r1\t0\tchrA\t100\t60\t6M4S\t*\t0\t0\tACGTNACGTT\tABCDEFGHIJ\t"
                   "NM:i:0\tMD:Z:6\tAS:i:6\tXS:i:0\tSA:Z:chrB,50,-,4M6S,20,0;\tHN:i:0\n");
    }
    SUBCASE("clips at both ends, both strands") {
        // Forward 2S5M3S: read bases 2-6 GTNAC, quals CDEFG.
        l[1] = alns.make(1, 49, false, "2S5M3S", "5", 0, 20, 5, 0);
        l[1].flag = 0x800;
        check_both(g.opt, ref, r, l, 1, NULL,
                   "r1\t2048\tchrB\t50\t20\t2H5M3H\t*\t0\t0\tGTNAC\tCDEFG\t"
                   "NM:i:0\tMD:Z:5\tAS:i:5\tXS:i:0\tSA:Z:chrA,100,+,6M4S,60,0;\n");
        // Reverse 1S5M4S: the reference-oriented first clip is the read's 3'
        // end, so the aligned read bases are 4-8 (NACGT): revcomp ACGTN, quals
        // reverse(EFGHI) = IHGFE.
        l[1] = alns.make(1, 49, true, "1S5M4S", "5", 0, 20, 5, 0);
        l[1].flag = 0x800;
        check_both(g.opt, ref, r, l, 1, NULL,
                   "r1\t2064\tchrB\t50\t20\t1H5M4H\t*\t0\t0\tACGTN\tIHGFE\t"
                   "NM:i:0\tMD:Z:5\tAS:i:5\tXS:i:0\tSA:Z:chrA,100,+,6M4S,60,0;\n");
    }
    SUBCASE("-Y: supplementary keeps S clips and the whole SEQ/QUAL") {
        g.opt->flag |= MEM_F_SOFTCLIP;
        check_both(g.opt, ref, r, l, 1, NULL,
                   "r1\t2048\tchrB\t50\t20\t6S4M\t*\t0\t0\tACGTNACGTT\tABCDEFGHIJ\t"
                   "NM:i:0\tMD:Z:4\tAS:i:4\tXS:i:0\tSA:Z:chrA,100,+,6M4S,60,0;\tHN:i:0\n");
        l[1].is_rev = 1;
        check_both(g.opt, ref, r, l, 1, NULL,
                   "r1\t2064\tchrB\t50\t20\t6S4M\t*\t0\t0\tAACGTNACGT\tJIHGFEDCBA\t"
                   "NM:i:0\tMD:Z:4\tAS:i:4\tXS:i:0\tSA:Z:chrA,100,+,6M4S,60,0;\tHN:i:0\n");
    }
    SUBCASE("a supplementary hit on an ALT contig keeps S clips") {
        l[1].is_alt = 1;
        check_both(g.opt, ref, r, l, 1, NULL,
                   "r1\t2048\tchrB\t50\t20\t6S4M\t*\t0\t0\tACGTNACGTT\tABCDEFGHIJ\t"
                   "NM:i:0\tMD:Z:4\tAS:i:4\tXS:i:0\tSA:Z:chrA,100,+,6M4S,60,0;\tHN:i:0\n");
    }
    SUBCASE("-M: the shorter hit prints FLAG 256 and keeps its SEQ and SA") {
        // mem_reg2sam marks it 0x10000 instead of 0x800 under -M; the writers
        // print that bit as 0x100. It is not a 0x100 record internally, so SEQ
        // (hard-clipped) and both SA tags are still written, as in bwa.
        l[1].flag = 0x10000;
        check_both(g.opt, ref, r, l, 1, NULL,
                   "r1\t256\tchrB\t50\t20\t6H4M\t*\t0\t0\tCGTT\tGHIJ\t"
                   "NM:i:0\tMD:Z:4\tAS:i:4\tXS:i:0\tSA:Z:chrA,100,+,6M4S,60,0;\tHN:i:0\n");
        check_both(g.opt, ref, r, l, 0, NULL,
                   "r1\t0\tchrA\t100\t60\t6M4S\t*\t0\t0\tACGTNACGTT\tABCDEFGHIJ\t"
                   "NM:i:0\tMD:Z:6\tAS:i:6\tXS:i:0\tSA:Z:chrB,50,+,6S4M,20,0;\tHN:i:0\n");
    }
    SUBCASE("three hits: SA lists the other two in list order") {
        l.push_back(alns.make(0, 499, true, "3M7S", "3", 1, 5, 3, 1));
        l[2].flag = 0x800;
        check_both(g.opt, ref, r, l, 1, NULL,
                   "r1\t2048\tchrB\t50\t20\t6H4M\t*\t0\t0\tCGTT\tGHIJ\tNM:i:0\tMD:Z:4\tAS:i:4\t"
                   "XS:i:0\tSA:Z:chrA,100,+,6M4S,60,0;chrA,500,-,3M7S,5,1;\tHN:i:0\n");
        // Two records written into one buffer are the two lines in order. The
        // reverse 3M7S hit aligns read bases 7-9 (GTT): SEQ AAC, QUAL JIH.
        __kstring_t str = {0, 0, NULL};
        mem_aln2sam(g.opt, &ref.bns, &str, &r.s, 3, l.data(), 0, NULL);
        mem_aln2sam(g.opt, &ref.bns, &str, &r.s, 3, l.data(), 2, NULL);
        REQUIRE(str.s != NULL);
        CHECK(std::string(str.s, str.l)
              == "r1\t0\tchrA\t100\t60\t6M4S\t*\t0\t0\tACGTNACGTT\tABCDEFGHIJ\tNM:i:0\tMD:Z:6\t"
                 "AS:i:6\tXS:i:0\tSA:Z:chrB,50,+,6S4M,20,0;chrA,500,-,3M7S,5,1;\tHN:i:0\n"
                 "r1\t2064\tchrA\t500\t5\t3M7H\t*\t0\t0\tAAC\tJIH\tNM:i:1\tMD:Z:3\t"
                 "AS:i:3\tXS:i:1\tSA:Z:chrA,100,+,6M4S,60,0;chrB,50,+,6S4M,20,0;\n");
        free(str.s);
    }
}

TEST_CASE("mem_io record: mate fields and FLAG bits and TLEN and MC and MQ"
          * doctest::test_suite("unit/mem_io")) {
    OptGuard g;
    Ref ref;
    Alns alns;
    Read r("r1", kBases, kQual);
    // Read 1 forward on chrA:100-109, read 2 reverse on chrA:300-309.
    mem_aln_t r1 = alns.make(0, 99, false, "10M", "10", 0, 60, 10, 0);
    mem_aln_t r2 = alns.make(0, 299, true, "10M", "10", 0, 50, 10, 0);

    SUBCASE("a forward/reverse pair on one contig, from each end") {
        // 0x1 (added) | 0x2 | 0x40 (caller) | 0x20 (mate reverse) = 99.
        // TLEN: leftmost base 100, rightmost 309 -> 309 - 100 + 1 = 210.
        std::vector<mem_aln_t> l1{r1};
        l1[0].flag = 0x42;
        check_both(g.opt, ref, r, l1, 0, &r2,
                   "r1\t99\tchrA\t100\t60\t10M\t=\t300\t210\tACGTNACGTT\tABCDEFGHIJ\t"
                   "NM:i:0\tMD:Z:10\tMC:Z:10M\tMQ:i:50\tAS:i:10\tXS:i:0\n");
        // 0x1 | 0x2 | 0x80 | 0x10 (own strand) = 147; TLEN -210.
        std::vector<mem_aln_t> l2{r2};
        l2[0].flag = 0x82;
        check_both(g.opt, ref, r, l2, 0, &r1,
                   "r1\t147\tchrA\t300\t50\t10M\t=\t100\t-210\tAACGTNACGT\tJIHGFEDCBA\t"
                   "NM:i:0\tMD:Z:10\tMC:Z:10M\tMQ:i:60\tAS:i:10\tXS:i:0\n");
    }
    SUBCASE("TLEN runs 5' end to 5' end: a forward/forward and an outward pair") {
        // Both forward: 5' ends 100 and 300 -> 300 - 100 + 1 = 201, not the
        // span 210. FLAG 0x40 | 0x1 = 65 (mate forward).
        mem_aln_t ff = alns.make(0, 299, false, "10M", "10", 0, 50, 10, 0);
        std::vector<mem_aln_t> l1{r1};
        l1[0].flag = 0x40;
        check_both(g.opt, ref, r, l1, 0, &ff,
                   "r1\t65\tchrA\t100\t60\t10M\t=\t300\t201\tACGTNACGTT\tABCDEFGHIJ\t"
                   "NM:i:0\tMD:Z:10\tMC:Z:10M\tMQ:i:50\tAS:i:10\tXS:i:0\n");
        // This read reverse on 100-109 (5' end 109), mate forward at 300:
        // 300 - 109 + 1 = 192. FLAG 0x40 | 0x1 | 0x10 = 81.
        mem_aln_t out = alns.make(0, 99, true, "10M", "10", 0, 60, 10, 0);
        std::vector<mem_aln_t> l2{out};
        l2[0].flag = 0x40;
        check_both(g.opt, ref, r, l2, 0, &ff,
                   "r1\t81\tchrA\t100\t60\t10M\t=\t300\t192\tAACGTNACGT\tJIHGFEDCBA\t"
                   "NM:i:0\tMD:Z:10\tMC:Z:10M\tMQ:i:50\tAS:i:10\tXS:i:0\n");
    }
    SUBCASE("TLEN counts reference bases: a deletion lengthens, an insertion does not") {
        // Mate 4M2D3M3I covers 300..308 (9 reference bases): 308 - 100 + 1 = 209.
        mem_aln_t m = alns.make(0, 299, true, "4M2D3M3I", "4^AC3", 5, 50, 1, 0);
        std::vector<mem_aln_t> l1{r1};
        l1[0].flag = 0x40;
        check_both(g.opt, ref, r, l1, 0, &m,
                   "r1\t97\tchrA\t100\t60\t10M\t=\t300\t209\tACGTNACGTT\tABCDEFGHIJ\t"
                   "NM:i:0\tMD:Z:10\tMC:Z:4M2D3M3I\tMQ:i:50\tAS:i:10\tXS:i:0\n");
    }
    SUBCASE("mate on another contig: RNEXT is its name, TLEN 0") {
        mem_aln_t m = alns.make(1, 9, false, "10M", "10", 0, 7, 10, 0);
        std::vector<mem_aln_t> l1{r1};
        l1[0].flag = 0x40;
        check_both(g.opt, ref, r, l1, 0, &m,
                   "r1\t65\tchrA\t100\t60\t10M\tchrB\t10\t0\tACGTNACGTT\tABCDEFGHIJ\t"
                   "NM:i:0\tMD:Z:10\tMC:Z:10M\tMQ:i:7\tAS:i:10\tXS:i:0\n");
    }
    SUBCASE("mate unmapped: 0x8, the mate placed at this read, no MC") {
        Read r_other("r1", kBases, kQual);
        mem_aln_t m = unmapped(g.opt, ref, r_other);
        // 0x40 | 0x1 | 0x8 = 73; the copied mate takes this read's strand too.
        std::vector<mem_aln_t> l1{r1};
        l1[0].flag = 0x40;
        check_both(g.opt, ref, r, l1, 0, &m,
                   "r1\t73\tchrA\t100\t60\t10M\t=\t100\t0\tACGTNACGTT\tABCDEFGHIJ\t"
                   "NM:i:0\tMD:Z:10\tMQ:i:0\tAS:i:10\tXS:i:0\n");
        // Reverse read: 0x40 | 0x1 | 0x8 | 0x10 | 0x20 = 121.
        std::vector<mem_aln_t> l2{r2};
        l2[0].flag = 0x40;
        check_both(g.opt, ref, r, l2, 0, &m,
                   "r1\t121\tchrA\t300\t50\t10M\t=\t300\t0\tAACGTNACGT\tJIHGFEDCBA\t"
                   "NM:i:0\tMD:Z:10\tMQ:i:0\tAS:i:10\tXS:i:0\n");
    }
    SUBCASE("read unmapped, mate mapped: placed at the mate, on its strand") {
        // 0x80 | 0x1 | 0x4 | 0x10 (mate's strand copied) | 0x20 = 181; CIGAR *,
        // no NM/MD; SEQ reverse-complemented because the copied strand is reverse.
        std::vector<mem_aln_t> l{unmapped(g.opt, ref, r)};
        l[0].flag |= 0x80;
        check_both(g.opt, ref, r, l, 0, &r2,
                   "r1\t181\tchrA\t300\t0\t*\t=\t300\t0\tAACGTNACGT\tJIHGFEDCBA\t"
                   "MC:Z:10M\tMQ:i:50\tAS:i:0\tXS:i:0\n");
        // Forward mate: 0x80 | 0x1 | 0x4 = 133.
        check_both(g.opt, ref, r, l, 0, &r1,
                   "r1\t133\tchrA\t100\t0\t*\t=\t100\t0\tACGTNACGTT\tABCDEFGHIJ\t"
                   "MC:Z:10M\tMQ:i:60\tAS:i:0\tXS:i:0\n");
    }
}

TEST_CASE("mem_io record: every optional field in bwa's order"
          * doctest::test_suite("unit/mem_io")) {
    OptGuard g;
    RgIdGuard rg_guard;
    Ref ref;
    Alns alns;
    // -C: the comment is appended verbatim after a TAB.
    Read r("r1", kBases, kQual, "BC:Z:ACGT\tXN:i:7");
    mem_aln_t p = alns.make(0, 99, false, "10M", "10", 0, 60, 10, 3);
    p.flag = 0x42;
    p.alt_sc = 8;  // pa = 10 / 8 = 1.250
    p.HN = 2;
    std::string xa = "chrB,+11,10M,1;chrA,-801,5M5S,0;";
    p.XA = const_cast<char *>(xa.c_str());
    mem_aln_t supp = alns.make(1, 49, false, "6S4M", "4", 0, 20, 4, 0);
    supp.flag = 0x800;
    mem_aln_t mate = alns.make(0, 299, true, "10M", "10", 0, 50, 10, 0);
    std::vector<mem_aln_t> l{p, supp};

    strcpy(bwa_rg_id, "grp1");
    g.opt->flag |= MEM_F_REF_HDR;  // -V: XR from chrA's FASTA comment, TAB -> space
    check_both(g.opt, ref, r, l, 0, &mate,
               "r1\t99\tchrA\t100\t60\t10M\t=\t300\t210\tACGTNACGTT\tABCDEFGHIJ\tNM:i:0\tMD:Z:10\t"
               "MC:Z:10M\tMQ:i:50\tAS:i:10\tXS:i:3\tRG:Z:grp1\tSA:Z:chrB,50,+,6S4M,20,0;\t"
               "pa:f:1.250\tXA:Z:chrB,+11,10M,1;chrA,-801,5M5S,0;\tHN:i:2\tBC:Z:ACGT\tXN:i:7\t"
               "XR:Z:first contig x\n",
               "r1\t99\tchrA\t100\t60\t10M\t=\t300\t210\tACGTNACGTT\tABCDEFGHIJ\tNM:i:0\tMD:Z:10\t"
               "MC:Z:10M\tMQ:i:50\tAS:i:10\tXS:i:3\tRG:Z:grp1\tSA:Z:chrB,50,+,6S4M,20,0;\t"
               "pa:f:1.25\tXA:Z:chrB,+11,10M,1;chrA,-801,5M5S,0;\tHN:i:2\tBC:Z:ACGT\tXN:i:7\t"
               "XR:Z:first contig x\n");

    // chrB has an empty comment: no XR. Without -V: no XR either.
    std::vector<mem_aln_t> on_b{alns.make(1, 0, false, "10M", "10", 0, 60, 10, 0)};
    check_both(g.opt, ref, r, on_b, 0, NULL,
               "r1\t0\tchrB\t1\t60\t10M\t*\t0\t0\tACGTNACGTT\tABCDEFGHIJ\tNM:i:0\tMD:Z:10\t"
               "AS:i:10\tXS:i:0\tRG:Z:grp1\tBC:Z:ACGT\tXN:i:7\n");
    g.opt->flag &= ~MEM_F_REF_HDR;
    std::vector<mem_aln_t> on_a{alns.make(0, 0, false, "10M", "10", 0, 60, 10, 0)};
    check_both(g.opt, ref, r, on_a, 0, NULL,
               "r1\t0\tchrA\t1\t60\t10M\t*\t0\t0\tACGTNACGTT\tABCDEFGHIJ\tNM:i:0\tMD:Z:10\t"
               "AS:i:10\tXS:i:0\tRG:Z:grp1\tBC:Z:ACGT\tXN:i:7\n");
}

TEST_CASE("mem_io record: a record of long fields appended to a buffer that holds text"
          * doctest::test_suite("unit/mem_io")) {
    // mem_aln2sam sizes its buffer once from the field lengths and then writes
    // without bounds checks, so each field is made long enough to dominate that
    // estimate. The expected line is assembled field by field in spec order.
    OptGuard g;
    RgIdGuard rg_guard;
    Ref ref;
    ref.names[1] = std::string(3000, 'c') + "_alt";
    ref.annos[0] = std::string(2000, 'd');
    ref.rebuild();
    Alns alns;

    std::string bases, qual;
    for (int i = 0; i < 1003; ++i) {
        bases += "ACGTN"[(i * 7 + i / 5) % 5];
        qual += (char)('!' + (i * 13) % 94);
    }
    const std::string name(250, 'q');
    const std::string comment = "CO:Z:" + std::string(5000, 'k');
    Read r(name, bases, qual, comment.c_str());

    // 400 ops: 1M1I repeated 200 times, then 603M -> 1003 query bases.
    std::string cigar;
    for (int i = 0; i < 200; ++i) cigar += "1M1I";
    cigar += "603M";
    std::string md = "1";
    for (int i = 0; i < 300; ++i) md += "A1";
    mem_aln_t p = alns.make(0, 0, false, cigar, md.c_str(), 500, 33, 700, 650);
    std::string xa;
    for (int i = 0; i < 400; ++i) xa += ref.names[1] + ",+" + std::to_string(i + 1) + ",50M,2;";
    p.XA = const_cast<char *>(xa.c_str());
    p.HN = 400;
    mem_aln_t supp = alns.make(1, 3, true, "1000M3S", "1000", 0, 1, 900, 0);
    supp.flag = 0x800;
    mem_aln_t mate = alns.make(1, 6, false, cigar, md.c_str(), 1, 44, 1, 0);
    std::vector<mem_aln_t> l{p, supp};
    memset(bwa_rg_id, 'g', 255);
    bwa_rg_id[255] = '\0';
    g.opt->flag |= MEM_F_REF_HDR;

    const std::string expected =
        name + "\t65\tchrA\t1\t33\t" + cigar + "\t" + ref.names[1] + "\t7\t0\t" + bases + "\t" + qual
        + "\tNM:i:500\tMD:Z:" + md + "\tMC:Z:" + cigar + "\tMQ:i:44\tAS:i:700\tXS:i:650\tRG:Z:"
        + std::string(255, 'g') + "\tSA:Z:" + ref.names[1] + ",4,-,1000M3S,1,0;\tXA:Z:" + xa
        + "\tHN:i:400\t" + comment + "\tXR:Z:" + ref.annos[0] + "\n";
    l[0].flag = 0x40;

    // Appended to a buffer that already holds text: the estimate is added to
    // the length already there (str->l + need), not to zero.
    const std::string before = "@previous record\n";
    __kstring_t str = {0, 0, NULL};
    str.m = before.size() + 1;
    str.s = (char *)malloc(str.m);
    REQUIRE(str.s != NULL);
    memcpy(str.s, before.c_str(), before.size() + 1);
    str.l = before.size();
    mem_aln2sam(g.opt, &ref.bns, &str, &r.s, 2, l.data(), 0, &mate);
    CHECK(str.l < str.m);
    CHECK(std::string(str.s, str.l) == before + expected);
    CHECK(strlen(str.s) == str.l);
    free(str.s);
}

TEST_CASE("mem_io record: each estimated field alone is long enough to overrun a short estimate"
          * doctest::test_suite("unit/mem_io")) {
    // One field at a time is made 100 kB or more while every other field stays
    // a few bytes, so the rest of the estimate is far below the line: had the
    // estimate dropped or undersized that field's term, the line would not fit
    // the capacity it set (sam_line checks str.l < str.m). The -R ID is at
    // most 255 characters, inside the estimate's fixed slack, so it has no
    // case of its own.
    OptGuard g;
    RgIdGuard rg_guard;
    bwa_rg_id[0] = '\0';
    Alns alns;
    const int kLong = 100000;
    // "1M1D" k times then "1M": 2k + 1 operations over k + 1 query bases.
    auto many_ops = [](int k) {
        std::string c;
        for (int i = 0; i < k; ++i) c += "1M1D";
        return c + "1M";
    };
    const std::string ops = many_ops(25000);
    const std::string ops_bases(25001, 'A'), ops_qual(25001, 'I');

    // A forward single-end record on chrA:100 with the given pieces.
    auto line = [](const std::string &name, int flag, const std::string &rname, const std::string &cigar,
                   const std::string &rnext, const std::string &pnext, const std::string &tlen,
                   const std::string &seq, const std::string &qual, const std::string &md,
                   const std::string &tags) {
        return name + "\t" + std::to_string(flag) + "\t" + rname + "\t100\t60\t" + cigar + "\t" + rnext
               + "\t" + pnext + "\t" + tlen + "\t" + seq + "\t" + qual + "\tNM:i:0\tMD:Z:" + md + tags
               + "\n";
    };
    const std::string as_xs = "\tAS:i:10\tXS:i:0";

    SUBCASE("QNAME") {
        Ref ref;
        const std::string name((size_t)kLong, 'q');
        Read r(name, kBases, kQual);
        std::vector<mem_aln_t> l{alns.make(0, 99, false, "10M", "10", 0, 60, 10, 0)};
        CHECK(sam_line(g.opt, ref, r, l, 0, NULL)
              == line(name, 0, "chrA", "10M", "*", "0", "0", kBases, kQual, "10", as_xs));
    }
    SUBCASE("SEQ and QUAL") {
        Ref ref;
        const std::string bases((size_t)kLong, 'C'), qual((size_t)kLong, 'J');
        Read r("r1", bases, qual);
        const std::string cigar = std::to_string(kLong) + "M";
        std::vector<mem_aln_t> l{alns.make(0, 99, false, cigar, "1", 0, 60, 10, 0)};
        CHECK(sam_line(g.opt, ref, r, l, 0, NULL)
              == line("r1", 0, "chrA", cigar, "*", "0", "0", bases, qual, "1", as_xs));
    }
    SUBCASE("CIGAR") {
        Ref ref;
        Read r("r1", ops_bases, ops_qual);
        std::vector<mem_aln_t> l{alns.make(0, 99, false, ops, "1", 0, 60, 10, 0)};
        CHECK(sam_line(g.opt, ref, r, l, 0, NULL)
              == line("r1", 0, "chrA", ops, "*", "0", "0", ops_bases, ops_qual, "1", as_xs));
    }
    SUBCASE("MD") {
        Ref ref;
        Read r("r1", kBases, kQual);
        const std::string md((size_t)kLong, '7');
        std::vector<mem_aln_t> l{alns.make(0, 99, false, "10M", md.c_str(), 0, 60, 10, 0)};
        CHECK(sam_line(g.opt, ref, r, l, 0, NULL)
              == line("r1", 0, "chrA", "10M", "*", "0", "0", kBases, kQual, md, as_xs));
    }
    SUBCASE("RNAME") {
        Ref ref;
        ref.names[0] = std::string((size_t)kLong, 'a');
        ref.rebuild();
        Read r("r1", kBases, kQual);
        std::vector<mem_aln_t> l{alns.make(0, 99, false, "10M", "10", 0, 60, 10, 0)};
        CHECK(sam_line(g.opt, ref, r, l, 0, NULL)
              == line("r1", 0, ref.names[0], "10M", "*", "0", "0", kBases, kQual, "10", as_xs));
    }
    SUBCASE("RNEXT, the mate's contig") {
        Ref ref;
        ref.names[1] = std::string((size_t)kLong, 'b');
        ref.rebuild();
        Read r("r1", kBases, kQual);
        mem_aln_t mate = alns.make(1, 9, false, "10M", "10", 0, 7, 10, 0);
        std::vector<mem_aln_t> l{alns.make(0, 99, false, "10M", "10", 0, 60, 10, 0)};
        l[0].flag = 0x40;
        CHECK(sam_line(g.opt, ref, r, l, 0, &mate)
              == line("r1", 65, "chrA", "10M", ref.names[1], "10", "0", kBases, kQual, "10",
                      "\tMC:Z:10M\tMQ:i:7" + as_xs));
    }
    SUBCASE("MC, the mate's CIGAR") {
        Ref ref;
        Read r("r1", kBases, kQual);
        mem_aln_t mate = alns.make(1, 9, false, ops, "1", 0, 7, 10, 0);
        std::vector<mem_aln_t> l{alns.make(0, 99, false, "10M", "10", 0, 60, 10, 0)};
        l[0].flag = 0x40;
        CHECK(sam_line(g.opt, ref, r, l, 0, &mate)
              == line("r1", 65, "chrA", "10M", "chrB", "10", "0", kBases, kQual, "10",
                      "\tMC:Z:" + ops + "\tMQ:i:7" + as_xs));
    }
    SUBCASE("SA, the other hits") {
        Ref ref;
        Read r("r1", ops_bases, ops_qual);
        mem_aln_t p = alns.make(0, 99, false, "25001M", "1", 0, 60, 10, 0);
        mem_aln_t supp = alns.make(1, 0, false, ops, "1", 0, 3, 9, 0);
        supp.flag = 0x800;
        std::vector<mem_aln_t> l{p, supp};
        CHECK(sam_line(g.opt, ref, r, l, 0, NULL)
              == line("r1", 0, "chrA", "25001M", "*", "0", "0", ops_bases, ops_qual, "1",
                      as_xs + "\tSA:Z:chrB,1,+," + ops + ",3,0;"));
    }
    SUBCASE("XA") {
        Ref ref;
        Read r("r1", kBases, kQual);
        std::string xa;
        while (xa.size() < (size_t)kLong) xa += "chrB,+" + std::to_string(xa.size() + 1) + ",10M,0;";
        mem_aln_t p = alns.make(0, 99, false, "10M", "10", 0, 60, 10, 0);
        p.XA = const_cast<char *>(xa.c_str());
        std::vector<mem_aln_t> l{p};
        CHECK(sam_line(g.opt, ref, r, l, 0, NULL)
              == line("r1", 0, "chrA", "10M", "*", "0", "0", kBases, kQual, "10", as_xs + "\tXA:Z:" + xa));
    }
    SUBCASE("the FASTQ comment (-C)") {
        Ref ref;
        const std::string comment = "CO:Z:" + std::string((size_t)kLong, 'k');
        Read r("r1", kBases, kQual, comment.c_str());
        std::vector<mem_aln_t> l{alns.make(0, 99, false, "10M", "10", 0, 60, 10, 0)};
        CHECK(sam_line(g.opt, ref, r, l, 0, NULL)
              == line("r1", 0, "chrA", "10M", "*", "0", "0", kBases, kQual, "10", as_xs + "\t" + comment));
    }
    SUBCASE("XR, the contig's FASTA comment (-V)") {
        Ref ref;
        ref.annos[0] = std::string((size_t)kLong, 'd');
        ref.rebuild();
        g.opt->flag |= MEM_F_REF_HDR;
        Read r("r1", kBases, kQual);
        std::vector<mem_aln_t> l{alns.make(0, 99, false, "10M", "10", 0, 60, 10, 0)};
        CHECK(sam_line(g.opt, ref, r, l, 0, NULL)
              == line("r1", 0, "chrA", "10M", "*", "0", "0", kBases, kQual, "10",
                      as_xs + "\tXR:Z:" + ref.annos[0]));
    }
}
