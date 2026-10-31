// Unit tests for the bwa read-input port (src/aln/bwaseqio.{h,cpp}): name, comment,
// barcode and quality-trimming rules of bwa_read_seq.

#include "doctest/doctest.h"

#include "../../src/aln/bwaseqio.h"

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <unistd.h>

namespace {

namespace ba = bwa_aln;

struct Rec {
    std::string name, comment, seq, qual;
};

// Writes a FASTQ to a temp file that is removed on destruction.
struct TmpFastq {
    std::string path;
    explicit TmpFastq(const std::vector<Rec> &recs) {
        const char *d = getenv("TMPDIR");
        std::string tmpl = std::string(d && *d ? d : "/tmp") + "/aln_seqio_XXXXXX";
        std::vector<char> buf(tmpl.begin(), tmpl.end());
        buf.push_back('\0');
        int fd = mkstemp(buf.data());
        REQUIRE(fd >= 0);
        close(fd);
        path = buf.data();
        FILE *fp = fopen(path.c_str(), "w");
        REQUIRE(fp != nullptr);
        for (const Rec &r : recs) {
            fprintf(fp, "@%s%s%s\n%s\n+\n%s\n", r.name.c_str(), r.comment.empty() ? "" : " ",
                    r.comment.c_str(), r.seq.c_str(), r.qual.c_str());
        }
        fclose(fp);
    }
    ~TmpFastq() { remove(path.c_str()); }
};

Rec rec(const std::string &name, const std::string &seq, const std::string &qual,
        const std::string &comment = "") {
    return Rec{name, comment, seq, qual};
}

std::vector<int> codes(const ba::ubyte_t *p, int n) {
    return std::vector<int>(p, p + n);
}

}  // namespace

TEST_CASE("aln_seqio: only a trailing /1 or /2 is trimmed") {
    std::vector<Rec> recs;
    for (const char *n : {"r/1", "r/2", "r/3", "r_1", "r/12"}) recs.push_back(rec(n, "ACGT", "IIII"));
    TmpFastq f(recs);
    ba::bwa_seqio_t *bs = ba::bwa_seq_open(f.path.c_str());
    int n = 0;
    ba::bwa_seq_t *s = ba::bwa_read_seq(bs, 0x40000, &n, 0, 0);
    REQUIRE(n == 5);
    CHECK(std::string(s[0].name) == "r");
    CHECK(std::string(s[1].name) == "r");
    CHECK(std::string(s[2].name) == "r/3");
    CHECK(std::string(s[3].name) == "r_1");
    CHECK(std::string(s[4].name) == "r/12");
    ba::bwa_free_read_seq(n, s);
    ba::bwa_seq_close(bs);
}

TEST_CASE("aln_seqio: -Y skips reads whose comment's first field after ':' is Y") {
    TmpFastq f({rec("a", "ACGT", "IIII", "1:Y:0:ACGT"), rec("b", "ACGT", "IIII", "1:N:0:ACGT"),
                rec("c", "ACGT", "IIII")});
    ba::bwa_seqio_t *bs = ba::bwa_seq_open(f.path.c_str());
    int n = 0;
    ba::bwa_seq_t *s = ba::bwa_read_seq(bs, 0x40000, &n, BWA_MODE_CFY, 0);
    REQUIRE(n == 2);
    CHECK(std::string(s[0].name) == "b");
    CHECK(std::string(s[1].name) == "c");
    ba::bwa_free_read_seq(n, s);
    ba::bwa_seq_close(bs);
}

TEST_CASE("aln_seqio: -I subtracts 31 from every quality before trimming") {
    TmpFastq f({rec("a", "ACGT", "hhhh")});
    ba::bwa_seqio_t *bs = ba::bwa_seq_open(f.path.c_str());
    int n = 0;
    ba::bwa_seq_t *s = ba::bwa_read_seq(bs, 0x40000, &n, BWA_MODE_IL13, 0);
    REQUIRE(n == 1);
    CHECK(std::string((const char *)s[0].qual) == "IIII");
    ba::bwa_free_read_seq(n, s);
    ba::bwa_seq_close(bs);
}

TEST_CASE("aln_seqio: -B 4 lower-cases barcode bases with qual-33 < 13 and drops len <= 4") {
    TmpFastq f({rec("a", "ACGTAAAA", "#I#IIIII"), rec("b", "ACGT", "IIII")});
    ba::bwa_seqio_t *bs = ba::bwa_seq_open(f.path.c_str());
    int n = 0;
    ba::bwa_seq_t *s = ba::bwa_read_seq(bs, 0x40000, &n, 4 << 24, 0);
    REQUIRE(n == 1);
    CHECK(std::string(s[0].bc) == "aCgT");
    CHECK((int)s[0].len == 4);
    CHECK((int)s[0].full_len == 4);
    CHECK(std::string((const char *)s[0].qual) == "IIII");
    // seq is stored reversed: the remaining "AAAA" is all zeros.
    CHECK(codes(s[0].seq, 4) == std::vector<int>({0, 0, 0, 0}));
    ba::bwa_free_read_seq(n, s);
    ba::bwa_seq_close(bs);
}

TEST_CASE("aln_seqio: -q never trims below 35 bp") {
    TmpFastq f({rec("a", std::string(50, 'A'), std::string(50, '#'))});
    ba::bwa_seqio_t *bs = ba::bwa_seq_open(f.path.c_str());
    int n = 0;
    ba::bwa_seq_t *s = ba::bwa_read_seq(bs, 0x40000, &n, 0, 20);
    REQUIRE(n == 1);
    CHECK((int)s[0].len == 35);
    CHECK((int)s[0].full_len == 50);
    CHECK(s[0].clip_len == 35);
    ba::bwa_free_read_seq(n, s);
    ba::bwa_seq_close(bs);
}

TEST_CASE("aln_seqio: seq is stored reversed, rseq reversed-and-complemented when COMPREAD") {
    TmpFastq f({rec("a", "AACG", "IIII")});
    ba::bwa_seqio_t *bs = ba::bwa_seq_open(f.path.c_str());
    int n = 0;
    ba::bwa_seq_t *s = ba::bwa_read_seq(bs, 0x40000, &n, BWA_MODE_COMPREAD, 0);
    REQUIRE(n == 1);
    CHECK(codes(s[0].seq, 4) == std::vector<int>({2, 1, 0, 0}));
    CHECK(codes(s[0].rseq, 4) == std::vector<int>({1, 2, 3, 3}));
    ba::bwa_free_read_seq(n, s);
    ba::bwa_seq_close(bs);
}

TEST_CASE("aln_seqio: lowercase maps like uppercase, IUPAC maps to 4") {
    TmpFastq f({rec("a", "acgtRYN", "IIIIIII")});
    ba::bwa_seqio_t *bs = ba::bwa_seq_open(f.path.c_str());
    int n = 0;
    ba::bwa_seq_t *s = ba::bwa_read_seq(bs, 0x40000, &n, 0, 0);
    REQUIRE(n == 1);
    CHECK(codes(s[0].seq, 7) == std::vector<int>({4, 4, 4, 3, 2, 1, 0}));
    ba::bwa_free_read_seq(n, s);
    ba::bwa_seq_close(bs);
}

TEST_CASE("aln_seqio: batches stop at n_needed kept reads") {
    std::vector<Rec> recs;
    for (int i = 0; i < 5; ++i) recs.push_back(rec("r" + std::to_string(i), "ACGT", "IIII"));
    TmpFastq f(recs);
    ba::bwa_seqio_t *bs = ba::bwa_seq_open(f.path.c_str());
    std::vector<int> sizes;
    for (;;) {
        int n = 0;
        ba::bwa_seq_t *s = ba::bwa_read_seq(bs, 2, &n, 0, 0);
        if (s == nullptr) {
            CHECK(n == 0);
            break;
        }
        sizes.push_back(n);
        ba::bwa_free_read_seq(n, s);
    }
    CHECK(sizes == std::vector<int>({2, 2, 1}));
    ba::bwa_seq_close(bs);
}

TEST_CASE("aln_seqio: empty file returns NULL with n == 0") {
    TmpFastq f({});
    ba::bwa_seqio_t *bs = ba::bwa_seq_open(f.path.c_str());
    int n = -1;
    CHECK(ba::bwa_read_seq(bs, 0x40000, &n, 0, 0) == nullptr);
    CHECK(n == 0);
    ba::bwa_seq_close(bs);
}
