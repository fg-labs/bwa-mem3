// test/unit/test_mem_io_bam_writer.cpp
//
// The --bam output file: bam_writer_open builds the header, bam_writer_write
// appends records and bam_writer_close finishes the BGZF stream; mem_aln2sam
// under --bam (opt->bam_mode) hands each record to the read instead of
// writing SAM text. The file is read back with htslib, so the checks are on
// what any BAM reader sees:
//
//   - the header lines and their order, which must match the SAM text path's
//     (bwa_print_sam_hdr2): the default @HD (none under --compat=bwa-mem2),
//     one @SQ per contig (or the index sidecar's @HD/@SQ when it has them),
//     the -H/-R lines, then @PG;
//   - the records, printed back as SAM text, field by field as worked out by
//     hand (the same rules as test_mem_io_sam_record.cpp);
//   - the file ends with the 28-byte BGZF end-of-file block the SAM
//     specification (section 4.1.2) gives byte for byte, at every compression
//     level;
//   - the content does not depend on the compression level or on the number
//     of BGZF threads.

// htslib headers first: they share the KSTRING_H include guard with
// bwa-mem3's kstring.h, and the BAM writer is compiled against htslib's.
#include "htslib/kstring.h"
#include "htslib/sam.h"

#include "doctest/doctest.h"
#include "mem_fixture.h"

#include "bam_writer.h"
#include "bwa.h"
#include "bwamem.h"
#include "compat_target.h"

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <string>
#include <vector>

// See test_mem_io_sam_record.cpp: bwa-mem3's own kstring_t, re-declared so
// the call resolves to libbwa.a's definition.
struct __kstring_t {
    size_t l, m;
    char *s;
};
void mem_aln2sam(const mem_opt_t *opt, const bntseq_t *bns, __kstring_t *str, bseq1_t *s,
                 int n, const mem_aln_t *list, int which, const mem_aln_t *m_);

namespace {

struct Ref {
    bntann1_t anns[2];
    bntseq_t bns;
    Ref() {
        memset(anns, 0, sizeof(anns));
        memset(&bns, 0, sizeof(bns));
        anns[0].offset = 0;
        anns[0].len = 1000;
        anns[0].name = const_cast<char *>("chrA");
        anns[0].anno = const_cast<char *>("");
        anns[1].offset = 1000;
        anns[1].len = 500;
        anns[1].name = const_cast<char *>("chrB");
        anns[1].anno = const_cast<char *>("");
        bns.l_pac = 1500;
        bns.n_seqs = 2;
        bns.anns = anns;
    }
};

// Read r1: ACGTNACGTT / ABCDEFGHIJ (reverse complement AACGTNACGT).
struct Read {
    std::string name{"r1"}, qual{"ABCDEFGHIJ"};
    std::vector<char> codes{0, 1, 2, 3, 4, 0, 1, 2, 3, 3};
    bseq1_t s;
    Read() {
        memset(&s, 0, sizeof(s));
        s.l_seq = 10;
        s.name = const_cast<char *>(name.c_str());
        s.seq = codes.data();
        s.qual = const_cast<char *>(qual.c_str());
        s.meth_base_ot = -1;
    }
};

// A full-length (10M) alignment whose CIGAR buffer is followed by its MD.
struct Aln {
    uint32_t buf[4];
    mem_aln_t a;
    Aln(int rid, int64_t pos0, bool rev, const char *md, int nm, int mapq, int score) {
        memset(buf, 0, sizeof(buf));
        buf[0] = 10u << 4;  // 10M
        strncpy((char *)(buf + 1), md, sizeof(buf) - sizeof(buf[0]) - 1);
        memset(&a, 0, sizeof(a));
        a.pos = pos0;
        a.rid = rid;
        a.is_rev = rev ? 1 : 0;
        a.mapq = (uint32_t)mapq;
        a.NM = (uint32_t)nm;
        a.n_cigar = 1;
        a.cigar = buf;
        a.HN = -1;
        a.score = score;
        a.sub = 0;
        a.meth_hypothesis = -1;
    }
};

const char *kFwd = "r1\t0\tchrA\t100\t60\t10M\t*\t0\t0\tACGTNACGTT\tABCDEFGHIJ\tNM:i:1\tMD:Z:4A5\tAS:i:5\tXS:i:0";
const char *kRev = "r1\t16\tchrB\t1\t7\t10M\t*\t0\t0\tAACGTNACGT\tJIHGFEDCBA\tNM:i:0\tMD:Z:10\tAS:i:10\tXS:i:0";
const char *kUnm = "r1\t4\t*\t0\t0\t*\t*\t0\t0\tACGTNACGTT\tABCDEFGHIJ\tAS:i:0\tXS:i:0";

const char *kPg = "@PG\tID:bwa-mem3\tPN:bwa-mem3\tVN:test";

// Writes the three records above through a bam_writer and returns the file's
// header text and records as htslib reads them back (each line ending '\n').
struct Written {
    std::string header, records, tail;
};

Written write_and_read(const std::string &path, const char *idx_hdr, const char *hdr_line,
                       int level, const compat_target_t *compat, int threads) {
    Ref ref;
    Read read;
    mem_opt_t *opt = mem_opt_init();
    Aln fwd(0, 99, false, "4A5", 1, 60, 5), rev(1, 0, true, "10", 0, 7, 10);
    mem_aln_t unm = mem_reg2aln(opt, &ref.bns, NULL, read.s.l_seq, read.s.seq, NULL);

    Written out;
    bam_writer_t *w = bam_writer_open(path.c_str(), &ref.bns, idx_hdr, hdr_line, kPg, level, compat,
                                      threads);
    REQUIRE(w != NULL);
    const mem_aln_t *alns[] = {&fwd.a, &rev.a, &unm};
    for (const mem_aln_t *a : alns) {
        bam1_t *b = bam_writer_alloc();
        REQUIRE(b != NULL);
        CHECK(mem_aln_to_bam(b, opt, &ref.bns, &read.s, 1, a, 0, NULL) == 0);
        CHECK(bam_writer_write(w, b) >= 0);
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

    FILE *fp = fopen(path.c_str(), "rb");
    REQUIRE(fp != NULL);
    std::string bytes;
    char chunk[4096];
    size_t n;
    while ((n = fread(chunk, 1, sizeof(chunk), fp)) > 0) bytes.append(chunk, n);
    fclose(fp);
    out.tail = bytes.size() >= 28 ? bytes.substr(bytes.size() - 28) : bytes;
    return out;
}

// The BGZF end-of-file marker: an empty BGZF block (SAM spec section 4.1.2).
std::string bgzf_eof() {
    const unsigned char eof[28] = {0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0, 0x42, 0x43,
                                   0x02, 0, 0x1b, 0, 0x03, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    return std::string((const char *)eof, sizeof(eof));
}

const std::string kRecords = std::string(kFwd) + "\n" + kRev + "\n" + kUnm + "\n";

}  // namespace

TEST_CASE("mem_io bam_writer: header order and records read back as written"
          * doctest::test_suite("unit/mem_io")) {
    bwa_tests::TempDir dir;
    const char *user = "@CO\tuser line\n@RG\tID:grp1\tSM:s";

    SUBCASE("default @HD, generated @SQ, user lines, @PG; levels 0 and 6") {
        const std::string expected_header = std::string("@HD\tVN:1.5\tSO:unsorted\tGO:query\n")
                                            + "@SQ\tSN:chrA\tLN:1000\n@SQ\tSN:chrB\tLN:500\n"
                                            + "@CO\tuser line\n@RG\tID:grp1\tSM:s\n" + kPg + "\n";
        const int levels[] = {0, 1, 6, 9};
        for (int level : levels) {
            CAPTURE(level);
            const Written w = write_and_read(dir.file("out.bam"), NULL, user, level, NULL, 0);
            CHECK(w.header == expected_header);
            CHECK(w.records == kRecords);
            CHECK(w.tail == bgzf_eof());
        }
        // BGZF threads: same content.
        const Written t = write_and_read(dir.file("t.bam"), NULL, user, 6, NULL, 2);
        CHECK(t.header == expected_header);
        CHECK(t.records == kRecords);
        CHECK(t.tail == bgzf_eof());
    }
    SUBCASE("the index sidecar's @HD and @SQ replace the generated ones") {
        const char *sidecar = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chrA\tLN:1000\tM5:aa\n"
                              "@SQ\tSN:chrB\tLN:500\tM5:bb";
        const Written w = write_and_read(dir.file("out.bam"), sidecar, user, 6, NULL, 0);
        CHECK(w.header
              == std::string("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chrA\tLN:1000\tM5:aa\n"
                             "@SQ\tSN:chrB\tLN:500\tM5:bb\n@CO\tuser line\n@RG\tID:grp1\tSM:s\n")
                     + kPg + "\n");
        CHECK(w.records == kRecords);
    }
    SUBCASE("--compat=bwa-mem2 emits no default @HD") {
        const compat_target_t *mem2 = compat_target_from_name("bwa-mem2");
        REQUIRE(mem2 != NULL);
        const Written w = write_and_read(dir.file("out.bam"), NULL, NULL, 6, mem2, 0);
        CHECK(w.header == std::string("@SQ\tSN:chrA\tLN:1000\n@SQ\tSN:chrB\tLN:500\n") + kPg + "\n");
        CHECK(w.records == kRecords);
    }
}

TEST_CASE("mem_io bam_writer: under --bam mem_aln2sam hands the record to the read"
          * doctest::test_suite("unit/mem_io")) {
    Ref ref;
    Read read;
    mem_opt_t *opt = mem_opt_init();
    opt->bam_mode = 1;
    Aln fwd(0, 99, false, "4A5", 1, 60, 5), rev(1, 0, true, "10", 0, 7, 10);

    __kstring_t str = {0, 0, NULL};
    mem_aln2sam(opt, &ref.bns, &str, &read.s, 1, &fwd.a, 0, NULL);
    mem_aln2sam(opt, &ref.bns, &str, &read.s, 1, &rev.a, 0, NULL);
    CHECK(str.l == 0);  // no SAM text under --bam
    REQUIRE(read.s.n_bams == 2);

    const char *htext = "@SQ\tSN:chrA\tLN:1000\n@SQ\tSN:chrB\tLN:500\n";
    sam_hdr_t *h = sam_hdr_parse(strlen(htext), htext);
    REQUIRE(h != NULL);
    std::string got;
    kstring_t ks = {0, 0, NULL};
    for (int i = 0; i < read.s.n_bams; ++i) {
        ks.l = 0;
        REQUIRE(sam_format1(h, (const bam1_t *)read.s.bams[i], &ks) >= 0);
        got += std::string(ks.s, ks.l) + "\n";
        bam_writer_free((bam1_t *)read.s.bams[i]);
    }
    CHECK(got == std::string(kFwd) + "\n" + kRev + "\n");
    free(ks.s);
    free(read.s.bams);
    free(str.s);
    sam_hdr_destroy(h);
    free(opt);
}
