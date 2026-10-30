// test/unit/test_mem_io_bam_writer.cpp
//
// The --bam output file's header: bam_writer_open builds it, and the file is
// read back with htslib, so the check is on what any BAM reader sees. The
// header lines and their order: the default @HD, one @SQ per contig (or the
// index sidecar's @HD/@SQ when it has them), the -H/-R lines, then @PG, as
// the SAM text path (bwa_print_sam_hdr2) writes them for these inputs. The
// records are read back too, as a control that the file is whole.
//
// The rest of the BAM writer is covered end to end: --compat=bwa-mem2's
// header by test/regression/compat_byte_identical.sh, the BGZF end-of-file
// block by bam_roundtrip.sh's quickcheck, output independent of the BGZF
// thread count by bam_threads_byte_identity.sh, and each record's fields by
// test_mem_io_sam_record.cpp's BAM half.

// htslib headers first: they share the KSTRING_H include guard with
// bwa-mem3's kstring.h, and the BAM writer is compiled against htslib's.
#include "htslib/kstring.h"
#include "htslib/sam.h"

#include "doctest/doctest.h"
#include "mem_record_fixture.h"
#include "temp_dir.h"

#include "bam_writer.h"
#include "bwa.h"
#include "bwamem.h"

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <string>

namespace {

using bwa_tests::OptGuard;

// Read r1: ACGTNACGTT / ABCDEFGHIJ (reverse complement AACGTNACGT).
struct Read : bwa_tests::RecordRead {
    Read() : RecordRead("r1", "ACGTNACGTT", "ABCDEFGHIJ") {}
};
typedef bwa_tests::RecordRef Ref;

const char *kFwd = "r1\t0\tchrA\t100\t60\t10M\t*\t0\t0\tACGTNACGTT\tABCDEFGHIJ\tNM:i:1\tMD:Z:4A5\tAS:i:5\tXS:i:0";
const char *kRev = "r1\t16\tchrB\t1\t7\t10M\t*\t0\t0\tAACGTNACGT\tJIHGFEDCBA\tNM:i:0\tMD:Z:10\tAS:i:10\tXS:i:0";
const char *kUnm = "r1\t4\t*\t0\t0\t*\t*\t0\t0\tACGTNACGTT\tABCDEFGHIJ\tAS:i:0\tXS:i:0";

const char *kPg = "@PG\tID:bwa-mem3\tPN:bwa-mem3\tVN:test";

// Writes the three records above through a bam_writer and returns the file's
// header text and records as htslib reads them back (each line ending '\n').
struct Written {
    std::string header, records;
};

Written write_and_read(const std::string &path, const char *idx_hdr, const char *hdr_line) {
    Ref ref;
    Read read;
    mem_opt_t *opt = mem_opt_init();
    bwa_tests::RecordAlns alns;
    const mem_aln_t fwd = alns.make(0, 99, false, "10M", "4A5", 1, 60, 5, 0);
    const mem_aln_t rev = alns.make(1, 0, true, "10M", "10", 0, 7, 10, 0);
    mem_aln_t unm = mem_reg2aln(opt, &ref.bns, NULL, read.s.l_seq, read.s.seq, NULL);

    Written out;
    bam_writer_t *w = bam_writer_open(path.c_str(), &ref.bns, idx_hdr, hdr_line, kPg, 6, NULL, 0);
    REQUIRE(w != NULL);
    const mem_aln_t *recs[] = {&fwd, &rev, &unm};
    for (const mem_aln_t *a : recs) {
        bam1_t *b = bam_writer_alloc();
        REQUIRE(b != NULL);
        REQUIRE(mem_aln_to_bam(b, opt, &ref.bns, &read.s, 1, a, 0, NULL) == 0);
        REQUIRE(bam_writer_write(w, b) >= 0);
        bam_writer_free(b);
    }
    CHECK(bam_writer_close(w) == 0);
    free(opt);

    htsFile *in = hts_open(path.c_str(), "r");
    REQUIRE(in != NULL);
    CHECK(hts_get_format(in)->format == bam);
    sam_hdr_t *h = sam_hdr_read(in);
    REQUIRE(h != NULL);
    out.header = sam_hdr_str(h);
    bam1_t *b = bam_init1();
    kstring_t ks = {0, 0, NULL};
    while (sam_read1(in, h, b) >= 0) {
        ks.l = 0;
        REQUIRE(sam_format1(h, b, &ks) >= 0);
        out.records += std::string(ks.s, ks.l) + "\n";
    }
    free(ks.s);
    bam_destroy1(b);
    sam_hdr_destroy(h);
    hts_close(in);
    return out;
}

const std::string kRecords = std::string(kFwd) + "\n" + kRev + "\n" + kUnm + "\n";

}  // namespace

TEST_CASE("mem_io bam_writer: header order, with and without an index sidecar"
          * doctest::test_suite("unit/mem_io")) {
    bwa_tests::TempDir dir;
    const char *user = "@CO\tuser line\n@RG\tID:grp1\tSM:s";

    SUBCASE("default @HD, generated @SQ, user lines, @PG") {
        const Written w = write_and_read(dir.file("out.bam"), NULL, user);
        CHECK(w.header
              == std::string("@HD\tVN:1.5\tSO:unsorted\tGO:query\n")
                     + "@SQ\tSN:chrA\tLN:1000\n@SQ\tSN:chrB\tLN:500\n"
                     + "@CO\tuser line\n@RG\tID:grp1\tSM:s\n" + kPg + "\n");
        CHECK(w.records == kRecords);
    }
    SUBCASE("the index sidecar's @HD and @SQ replace the generated ones") {
        const char *sidecar = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chrA\tLN:1000\tM5:aa\n"
                              "@SQ\tSN:chrB\tLN:500\tM5:bb";
        const Written w = write_and_read(dir.file("out.bam"), sidecar, user);
        CHECK(w.header
              == std::string("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chrA\tLN:1000\tM5:aa\n"
                             "@SQ\tSN:chrB\tLN:500\tM5:bb\n@CO\tuser line\n@RG\tID:grp1\tSM:s\n")
                     + kPg + "\n");
        CHECK(w.records == kRecords);
    }
}
