// test/unit/test_mem_io_reader.cpp
//
// The two FASTA/FASTQ batch readers of `bwa-mem3 mem` at the level `mem`
// calls them: the default one (fast_reader -> fr_fastq -> bseq_read_fast) and
// --legacy-reader's (kseq over zlib's gzFile, bseq_read_orig), on plain, gzip
// (one member and several) and BGZF files written by the test; and
// bseq_classify, the -p smart-pairing split.
//
// The legacy reader has no other test. The default reader's parser is held to
// kseq on plain input, with tiny buffers, by test/fr_fastq_diff_test.c; its
// codec layer by test/fast_reader_selftest.c; its -C comment gate and
// read-number trimming by test/bseq_read_comment_copy_test.c; and a gzip
// stream cut mid-member by test/bseq_read_truncated_gzip_test.c. What this
// file adds is the batch layer above the parser, for both readers: the records
// each bseq_read_* returns, worked out by hand from the FASTA/FASTQ grammar
// kseq.h and src/fr_fastq.h document, through each codec:
//
//   - outside a record, the reader skips to the next '>' or '@'; text before
//     the first header and blank lines are skipped;
//   - the name runs to the first whitespace; the comment is the rest of that
//     line, absent when empty, with a trailing '\r' removed only when the line
//     is longer than it (kseq's rule), so a comment that is just '\r' is "\r";
//   - a record that ends after its sequence (EOF, or the next header, before a
//     '+' line) has no quality, as a FASTA record;
//   - the sequence is every line up to the next header or '+', newlines (and
//     CRLF) removed, bases kept as given (case, IUPAC);
//   - FASTQ quality lines are read until they reach the sequence length, so a
//     quality line may start with '@' or '+';
//   - bwa's trim_readno drops a trailing "/<digit>" from a name longer than two
//     characters;
//   - a batch stops once its bases reach the chunk size AND it holds an even
//     number of records; paired files are read in lockstep, interleaved (end 1,
//     end 2, ...), and stop at the shorter file, with a warning on stderr;
//   - a quality string of the wrong length, in either file of a pair, is a
//     fatal error, and so is a gzip or BGZF stream cut at a member or block
//     boundary: the default reader fails as it reads, the legacy one when
//     `mem` closes the input (err_gzclose), since zlib's gzread reports the
//     cut as end of file.
//
// The two readers differ in one documented way: the default reader keeps a
// record's comment only when asked to (copy_comment, set by -C), since `mem`
// drops comments otherwise; the legacy reader always keeps it.
#include "doctest/doctest.h"
#include "run_in_child.h"
#include "stderr_capture.h"
#include "temp_dir.h"

#include "../gz_test_util.h"

#include "bwa.h"
#include "fast_reader.h"
#include "fast_reader_bseq.h"
#include "kseq.h"
#include "read_arena.h"
#include "utils.h"

#include <fcntl.h>
#include <sys/wait.h>
#include <unistd.h>
#include <zlib.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

KSEQ_DECLARE(gzFile)

namespace {

// One parsed record, NULL fields shown as "<null>".
struct Rec {
    std::string name, comment, seq, qual;
    int l_seq, id;
};

std::string show(const Rec &r) {
    return r.name + "|" + r.comment + "|" + r.seq + "|" + r.qual + "|" + std::to_string(r.l_seq)
           + "|" + std::to_string(r.id);
}

std::string show(const std::vector<std::vector<Rec> > &batches) {
    std::string out;
    for (size_t b = 0; b < batches.size(); ++b) {
        out += "[batch " + std::to_string(b) + "]\n";
        for (const Rec &r : batches[b]) out += "  " + show(r) + "\n";
    }
    return out;
}

const char *kNull = "<null>";

Rec rec(const char *name, const char *comment, const char *seq, const char *qual, int id) {
    Rec r;
    r.name = name;
    r.comment = comment ? comment : kNull;
    r.seq = seq;
    r.qual = qual ? qual : kNull;
    r.l_seq = (int)strlen(seq);
    r.id = id;
    return r;
}

// Takes a batch apart into Recs and frees it (comments are per-read heap
// strings; name/seq/qual live in the batch's arena).
std::vector<Rec> drain(bseq1_t *seqs, int n, read_arena_t *arena) {
    std::vector<Rec> out;
    for (int i = 0; i < n; ++i) {
        Rec r;
        r.name = seqs[i].name;
        r.comment = seqs[i].comment ? seqs[i].comment : kNull;
        r.seq = seqs[i].seq;
        r.qual = seqs[i].qual ? seqs[i].qual : kNull;
        r.l_seq = seqs[i].l_seq;
        r.id = seqs[i].id;
        CHECK(strlen(seqs[i].seq) == (size_t)seqs[i].l_seq);
        if (seqs[i].qual) CHECK(strlen(seqs[i].qual) == (size_t)seqs[i].l_seq);
        CHECK(seqs[i].sam == NULL);
        CHECK(seqs[i].meth_base_ot == -1);
        free(seqs[i].comment);
        out.push_back(r);
    }
    free(seqs);
    if (arena) read_arena_destroy(arena);
    return out;
}

// Every batch of the --legacy-reader path (kseq over zlib's gzFile), until EOF.
// This reader always keeps the comment (mem frees it afterwards unless -C).
std::vector<std::vector<Rec> > read_legacy(const std::string &f1, const std::string &f2,
                                           int64_t chunk) {
    gzFile g1 = gzopen(f1.c_str(), "r");
    REQUIRE(g1 != NULL);
    gzFile g2 = f2.empty() ? NULL : gzopen(f2.c_str(), "r");
    if (!f2.empty()) REQUIRE(g2 != NULL);
    kseq_t *k1 = kseq_init(g1), *k2 = g2 ? kseq_init(g2) : NULL;
    std::vector<std::vector<Rec> > batches;
    for (int guard = 0; guard < 1000; ++guard) {
        int n = 0;
        int64_t size = 0;
        read_arena_t *arena = NULL;
        bseq1_t *seqs = bseq_read_orig(chunk, &n, k1, k2, &size, &arena);
        if (n == 0) {
            CHECK(seqs == NULL);
            CHECK(arena == NULL);
            break;
        }
        int64_t bases = 0;
        for (int i = 0; i < n; ++i) bases += seqs[i].l_seq;
        CHECK(size == bases);
        batches.push_back(drain(seqs, n, arena));
    }
    kseq_destroy(k1);
    if (k2) kseq_destroy(k2);
    gzclose(g1);
    if (g2) gzclose(g2);
    return batches;
}

// Each member a complete gzip stream (test/gz_test_util.h's gz_member);
// concatenated they are one file.
std::string gzip_bytes(const std::vector<std::string> &members) {
    std::string out;
    for (const std::string &m : members) {
        std::vector<unsigned char> buf(compressBound((uLong)m.size()) + 64);
        const size_t len = gz_member((const unsigned char *)m.data(), m.size(), buf.data(), buf.size());
        out.append((const char *)buf.data(), len);
    }
    return out;
}

// BGZF (SAM spec section 4.1): gzip members of at most 64 KiB, each with a
// "BC" extra subfield holding the member's size - 1, then an empty member as
// the end-of-file marker. Written here with zlib's raw deflate in blocks of
// `step` input bytes, small by default so the file has several and records
// straddle them (37 bytes; 30000 for a text over 4 KiB, which keeps the
// long-record case fast while its record still spans many blocks).
std::string bgzf_block(const std::string &data) {
    z_stream z;
    memset(&z, 0, sizeof(z));
    REQUIRE(deflateInit2(&z, 6, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY) == Z_OK);
    std::vector<unsigned char> body(compressBound((uLong)data.size()) + 16);
    z.next_in = (Bytef *)data.data();
    z.avail_in = (uInt)data.size();
    z.next_out = body.data();
    z.avail_out = (uInt)body.size();
    REQUIRE(deflate(&z, Z_FINISH) == Z_STREAM_END);
    const size_t clen = body.size() - z.avail_out;
    deflateEnd(&z);
    const size_t bsize = 18 + clen + 8;  // header, deflate data, CRC32 + ISIZE
    REQUIRE(bsize <= 65536);
    const unsigned char hdr[18] = {0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, 'B', 'C', 2, 0,
                                   (unsigned char)((bsize - 1) & 0xff),
                                   (unsigned char)((bsize - 1) >> 8)};
    const uLong crc = crc32(crc32(0L, Z_NULL, 0), (const Bytef *)data.data(), (uInt)data.size());
    std::string out((const char *)hdr, sizeof(hdr));
    out.append((const char *)body.data(), clen);
    for (int k = 0; k < 4; ++k) out += (char)((crc >> (8 * k)) & 0xff);
    for (int k = 0; k < 4; ++k) out += (char)((data.size() >> (8 * k)) & 0xff);
    return out;
}

std::string bgzf_bytes(const std::string &text, size_t step = 0) {
    if (step == 0) step = text.size() > 4096 ? 30000 : 37;
    std::string out;
    for (size_t i = 0; i < text.size(); i += step) out += bgzf_block(text.substr(i, step));
    return out + bgzf_block("");  // the end-of-file marker
}

// Every batch of the default reader (fast_reader_dopen -> fast_kseq_init ->
// bseq_read_fast), until EOF.
std::vector<std::vector<Rec> > read_fast(const std::string &f1, const std::string &f2, int64_t chunk,
                                         int copy_comment) {
    const char *err = NULL;
    fast_reader_t *r1 = fast_reader_dopen(open(f1.c_str(), O_RDONLY), &err);
    CAPTURE(std::string(err ? err : ""));
    REQUIRE(r1 != NULL);
    fast_reader_t *r2 = NULL;
    if (!f2.empty()) {
        r2 = fast_reader_dopen(open(f2.c_str(), O_RDONLY), &err);
        CAPTURE(std::string(err ? err : ""));
        REQUIRE(r2 != NULL);
    }
    void *k1 = fast_kseq_init(r1), *k2 = r2 ? fast_kseq_init(r2) : NULL;
    std::vector<std::vector<Rec> > batches;
    for (int guard = 0; guard < 1000; ++guard) {
        int n = 0;
        int64_t size = 0;
        read_arena_t *arena = NULL;
        bseq1_t *seqs = bseq_read_fast(chunk, &n, k1, k2, &size, &arena, copy_comment);
        if (n == 0) {
            CHECK(seqs == NULL);
            CHECK(arena == NULL);
            break;
        }
        int64_t bases = 0;
        for (int i = 0; i < n; ++i) bases += seqs[i].l_seq;
        CHECK(size == bases);
        batches.push_back(drain(seqs, n, arena));
    }
    fast_kseq_destroy(k1);
    fast_reader_close(r1);
    if (r2) {
        fast_kseq_destroy(k2);
        fast_reader_close(r2);
    }
    return batches;
}

// `batches` with every comment removed: what the default reader gives
// without copy_comment.
std::vector<std::vector<Rec> > without_comments(std::vector<std::vector<Rec> > batches) {
    for (size_t b = 0; b < batches.size(); ++b)
        for (size_t i = 0; i < batches[b].size(); ++i) batches[b][i].comment = kNull;
    return batches;
}

// Both readers must give `expected` for one file; the default reader with and
// without copy_comment.
void check_readers(const std::string &file, int64_t chunk, const std::vector<std::vector<Rec> > &expected) {
    CAPTURE(file);
    const std::string want = show(expected);
    CHECK(show(read_legacy(file, "", chunk)) == want);
    CHECK(show(read_fast(file, "", chunk, 1)) == want);
    CHECK(show(read_fast(file, "", chunk, 0)) == show(without_comments(expected)));
}

// check_readers on `text` written as plain, gzip (one member and split into
// three) and BGZF.
void check_all_codecs(const std::string &text, int64_t chunk,
                      const std::vector<std::vector<Rec> > &expected) {
    bwa_tests::TempDir dir;
    dir.write("in.fq", text);
    dir.write("in.fq.gz", gzip_bytes({text}));
    const size_t a = text.size() / 3, b = 2 * text.size() / 3;
    dir.write("multi.fq.gz", gzip_bytes({text.substr(0, a), text.substr(a, b - a), text.substr(b)}));
    dir.write("in.fq.bgz", bgzf_bytes(text));
    const char *files[] = {"in.fq", "in.fq.gz", "multi.fq.gz", "in.fq.bgz"};
    for (const char *f : files) check_readers(dir.file(f), chunk, expected);
}

const int64_t kWholeFile = 1 << 30;

// How many times `text` repeats `line`, or -1 if it holds anything else: a
// reader called again after a paired input ran short warns again.
int warning_lines(const std::string &text, const std::string &line) {
    int n = 0;
    for (size_t at = 0; at < text.size(); at += line.size(), ++n)
        if (text.compare(at, line.size(), line) != 0) return -1;
    return n;
}

}  // namespace

TEST_CASE("mem_io reader: FASTQ names and comments and read-number trimming"
          * doctest::test_suite("unit/mem_io")) {
    const std::string text =
        "@r1/1 BC:Z:AC\tXX:i:1\nACGTN\n+\nIIIII\n"  // trimmed; comment keeps its TAB
        "@r2/12\nacgtRYKM\n+r2/12\n!!##$$%%\n"     // "/12": not a /<digit> suffix
        "@r3/a\nAC\n+\nII\n"                       // "/a": not trimmed
        "@/1\nA\n+\nI\n"                           // two characters: not trimmed
        "@x/2 \nAC\n+\nII\n"                       // empty comment is no comment
        "@tab\tc1 c2\nGG\n+\n#$\n";                // TAB ends the name too
    const std::vector<std::vector<Rec> > with_comments{{
        rec("r1", "BC:Z:AC\tXX:i:1", "ACGTN", "IIIII", 0),
        rec("r2/12", NULL, "acgtRYKM", "!!##$$%%", 1),
        rec("r3/a", NULL, "AC", "II", 2),
        rec("/1", NULL, "A", "I", 3),
        rec("x", NULL, "AC", "II", 4),
        rec("tab", "c1 c2", "GG", "#$", 5),
    }};
    check_all_codecs(text, kWholeFile, with_comments);
}

TEST_CASE("mem_io reader: FASTA and multi-line records and CRLF and blank lines"
          * doctest::test_suite("unit/mem_io")) {
    SUBCASE("multi-line FASTA, no quality") {
        check_all_codecs(">c1 first desc\nACGT\nAC\n\n>c2\nGG\nT\n", kWholeFile,
                         {{rec("c1", "first desc", "ACGTAC", NULL, 0), rec("c2", NULL, "GGT", NULL, 1)}});
    }
    SUBCASE("multi-line FASTQ sequence and quality") {
        check_all_codecs("@q1\nACGT\nAC\n+\nIIII\nJJ\n@q2\nT\n+\n5\n", kWholeFile,
                         {{rec("q1", NULL, "ACGTAC", "IIIIJJ", 0), rec("q2", NULL, "T", "5", 1)}});
    }
    SUBCASE("quality lines that start with @ or +") {
        check_all_codecs("@a\nACGT\n+\n@III\n@b\nAC\n+\n+I\n@c\nG\n+\n@\n", kWholeFile,
                         {{rec("a", NULL, "ACGT", "@III", 0), rec("b", NULL, "AC", "+I", 1),
                           rec("c", NULL, "G", "@", 2)}});
    }
    SUBCASE("CRLF line endings") {
        check_all_codecs("@c1 cm\r\nACGT\r\n+\r\nIIII\r\n@c2\r\nGG\r\nTT\r\n+\r\nABCD\r\n",
                         kWholeFile,
                         {{rec("c1", "cm", "ACGT", "IIII", 0), rec("c2", NULL, "GGTT", "ABCD", 1)}});
        // A comment that is only the CR keeps it (kseq strips a CR only from a
        // line longer than one byte).
        check_all_codecs("@b \r\nA\r\n+\r\nI\r\n", kWholeFile, {{rec("b", "\r", "A", "I", 0)}});
    }
    SUBCASE("text before the first record and blank lines between records") {
        check_all_codecs("\n\njunk line\n@a\nAC\n+\nII\n\n\n@b\nGT\n+\nJJ", kWholeFile,
                         {{rec("a", NULL, "AC", "II", 0), rec("b", NULL, "GT", "JJ", 1)}});
    }
    SUBCASE("a record that ends after its sequence has no quality") {
        check_all_codecs("@a\nAC\n+\nII\n@b\nACGT\n@c\n", kWholeFile,
                         {{rec("a", NULL, "AC", "II", 0), rec("b", NULL, "ACGT", NULL, 1),
                           rec("c", NULL, "", NULL, 2)}});
    }
    SUBCASE("a long record") {
        std::string seq, qual;
        for (int i = 0; i < 300000; ++i) {
            seq += "ACGT"[(i * 7 + i / 13) % 4];
            qual += (char)('!' + i % 41);
        }
        check_all_codecs("@long\n" + seq + "\n+\n" + qual + "\n@s\nA\n+\nI\n", kWholeFile,
                         {{rec("long", NULL, seq.c_str(), qual.c_str(), 0), rec("s", NULL, "A", "I", 1)}});
    }
    SUBCASE("an empty file") { check_all_codecs("", kWholeFile, {}); }
}

TEST_CASE("mem_io reader: batches end at the chunk size on an even record count"
          * doctest::test_suite("unit/mem_io")) {
    // Batching sits above the codec, so plain input is enough.
    bwa_tests::TempDir dir;
    std::string text;
    std::vector<Rec> all;
    for (int i = 0; i < 7; ++i) {
        const std::string name = "r" + std::to_string(i);
        const std::string seq = std::string(10, "ACGT"[i % 4]);
        text += "@" + name + "\n" + seq + "\n+\n" + std::string(10, (char)('A' + i)) + "\n";
        all.push_back(rec(name.c_str(), NULL, seq.c_str(), std::string(10, (char)('A' + i)).c_str(), 0));
    }
    // chunk 25: after r0..r2 the batch holds 30 bases but an odd 3 records, so
    // r3 joins (40 bases, 4 records); the rest (r4..r6) is the next batch. ids
    // count from 0 in each batch.
    std::vector<std::vector<Rec> > b25(2);
    for (int i = 0; i < 7; ++i) {
        Rec r = all[(size_t)i];
        r.id = i < 4 ? i : i - 4;
        b25[i < 4 ? 0 : 1].push_back(r);
    }
    dir.write("in.fq", text);
    check_readers(dir.file("in.fq"), 25, b25);
    // chunk 20: 2 records reach it exactly -> batches of 2, 2, 2, 1.
    std::vector<std::vector<Rec> > b20(4);
    for (int i = 0; i < 7; ++i) {
        Rec r = all[(size_t)i];
        r.id = i % 2;
        b20[(size_t)(i / 2)].push_back(r);
    }
    check_readers(dir.file("in.fq"), 20, b20);
}

TEST_CASE("mem_io reader: two files are read in lockstep and interleaved"
          * doctest::test_suite("unit/mem_io")) {
    bwa_tests::TempDir dir;
    dir.write("r1.fq", "@p1/1 x\nAAAA\n+\nIIII\n@p2/1\nCCCC\n+\nJJJJ\n@p3/1\nGG\n+\nKK\n");
    dir.write("r2.fq", "@p1/2 y\nTTTT\n+\n5555\n@p2/2\nGGGG\n+\n6666\n@p3/2\nTT\n+\n77\n");
    dir.write("short.fq", "@p1/2\nTTTT\n+\n5555\n");
    const std::vector<std::vector<Rec> > pairs{{
        rec("p1", "x", "AAAA", "IIII", 0), rec("p1", "y", "TTTT", "5555", 1),
        rec("p2", NULL, "CCCC", "JJJJ", 2), rec("p2", NULL, "GGGG", "6666", 3),
        rec("p3", NULL, "GG", "KK", 4), rec("p3", NULL, "TT", "77", 5),
    }};
    CHECK(show(read_legacy(dir.file("r1.fq"), dir.file("r2.fq"), kWholeFile)) == show(pairs));
    CHECK(show(read_fast(dir.file("r1.fq"), dir.file("r2.fq"), kWholeFile, 1)) == show(pairs));

    // Chunk 5: the first pair holds 8 bases, so each batch is one pair.
    std::vector<std::vector<Rec> > one_pair_each(3);
    for (size_t i = 0; i < 6; ++i) {
        Rec r = pairs[0][i];
        r.id = (int)(i % 2);
        one_pair_each[i / 2].push_back(r);
    }
    CHECK(show(read_legacy(dir.file("r1.fq"), dir.file("r2.fq"), 5)) == show(one_pair_each));
    CHECK(show(read_fast(dir.file("r1.fq"), dir.file("r2.fq"), 5, 1)) == show(one_pair_each));

    // A shorter file ends the input at its last record, either way round, and
    // the reader says which file ran short: that warning is the only sign that
    // the longer file's remaining records were dropped.
    const std::vector<std::vector<Rec> > first_pair{
        {rec("p1", "x", "AAAA", "IIII", 0), rec("p1", NULL, "TTTT", "5555", 1)}};
    const std::vector<std::vector<Rec> > first_pair_swapped{
        {rec("p1", NULL, "TTTT", "5555", 0), rec("p1", "x", "AAAA", "IIII", 1)}};
    struct Case {
        const char *f1, *f2, *which;
    };
    const Case cases[] = {{"r1.fq", "short.fq", "2nd"}, {"short.fq", "r1.fq", "1st"}};
    for (const Case &c : cases) {
        CAPTURE(std::string(c.f1));
        const std::vector<std::vector<Rec> > &want = c.which[0] == '2' ? first_pair : first_pair_swapped;
        const std::string warning = std::string(" the ") + c.which + " file has fewer sequences.\n";
        {
            bwa_tests::StderrCapture err;
            const std::string got = show(read_legacy(dir.file(c.f1), dir.file(c.f2), kWholeFile));
            CHECK(warning_lines(err.text(), "[W::bseq_read_orig]" + warning) > 0);
            CHECK(got == show(want));
        }
        {
            bwa_tests::StderrCapture err;
            const std::string got = show(read_fast(dir.file(c.f1), dir.file(c.f2), kWholeFile, 1));
            CHECK(warning_lines(err.text(), "[W::bseq_read_fast]" + warning) > 0);
            CHECK(got == show(want));
        }
    }
}

namespace {

// Runs a reader over a file, or a pair of files, in a forked child (stderr
// silenced, no core file) and closes the input as `mem` does, and returns the
// child's wait status: a reader that rejects the input, or an input whose
// close reports an error (err_gzclose, for the legacy reader), exits with
// EXIT_FAILURE.
int read_in_child(bool legacy, const std::string &f1, const std::string &f2 = "") {
    return bwa_tests::run_in_child([&] {
        int n = 0;
        int64_t size = 0;
        read_arena_t *arena = NULL;
        if (legacy) {
            gzFile g1 = gzopen(f1.c_str(), "r"), g2 = f2.empty() ? NULL : gzopen(f2.c_str(), "r");
            if (g1 == NULL || (!f2.empty() && g2 == NULL)) _exit(4);
            kseq_t *k1 = kseq_init(g1), *k2 = g2 ? kseq_init(g2) : NULL;
            while (bseq_read_orig(1 << 20, &n, k1, k2, &size, &arena) != NULL && n > 0) arena = NULL;
            kseq_destroy(k1);
            err_gzclose(g1);
            if (k2) {
                kseq_destroy(k2);
                err_gzclose(g2);
            }
        } else {
            fast_reader_t *r1 = fast_reader_dopen(open(f1.c_str(), O_RDONLY), NULL);
            fast_reader_t *r2 = f2.empty() ? NULL : fast_reader_dopen(open(f2.c_str(), O_RDONLY), NULL);
            if (r1 == NULL || (!f2.empty() && r2 == NULL)) _exit(4);
            void *k1 = fast_kseq_init(r1), *k2 = r2 ? fast_kseq_init(r2) : NULL;
            while (bseq_read_fast(1 << 20, &n, k1, k2, &size, &arena, 1) != NULL && n > 0) arena = NULL;
        }
    });
}

}  // namespace

TEST_CASE("mem_io reader: a wrong-length quality string or a cut gzip or BGZF stream is fatal"
          * doctest::test_suite("unit/mem_io")) {
    bwa_tests::TempDir dir;
    const std::string good = "@a\nAC\n+\nII\n@b\nACGT\n+\nIIII\n";
    dir.write("good.fq", good);
    // A quality shorter than its sequence at the end of input, a '+' line with
    // nothing after it, and a quality longer than its sequence.
    dir.write("short_qual.fq", "@a\nACGT\n+\nII\n");
    dir.write("no_qual.fq", "@a\nAC\n+\nII\n@b\nACGT\n+");
    dir.write("long_qual.fq", "@a\nAC\n+\nIIII\n");
    // Streams cut right after a complete member or block, so what does decode
    // is whole records: only the cut itself can fail the read. The gzip file
    // is a member of 400 records and the 10-byte header of a second one; the
    // BGZF file is one block of the same records and the 18-byte header of the
    // next, with no end-of-file block.
    std::string many;
    for (int i = 0; i < 400; ++i) many += good;
    dir.write("whole.fq.gz", gzip_bytes({many, many}));
    dir.write("cut.fq.gz", gzip_bytes({many}) + gzip_bytes({many}).substr(0, 10));
    dir.write("whole.fq.bgz", bgzf_bytes(many + many, many.size()));
    dir.write("cut.fq.bgz", bgzf_block(many) + bgzf_block(many).substr(0, 18));
    const char *bad[] = {"short_qual.fq", "no_qual.fq", "long_qual.fq", "cut.fq.gz", "cut.fq.bgz"};
    const char *ok_files[] = {"good.fq", "whole.fq.gz", "whole.fq.bgz"};
    const bool legacy[] = {false, true};
    for (bool leg : legacy) {
        CAPTURE(leg);
        for (const char *f : bad) {
            CAPTURE(std::string(f));
            const int st = read_in_child(leg, dir.file(f));
            CHECK(WIFEXITED(st));
            CHECK(WEXITSTATUS(st) == EXIT_FAILURE);
        }
        // The same faults in the second file of a pair: at its first record,
        // and after the first file has ended (the reader drains the second).
        const char *bad_mates[] = {"short_qual.fq", "long_qual.fq"};
        for (const char *f : bad_mates) {
            CAPTURE(std::string(f));
            const int st = read_in_child(leg, dir.file("good.fq"), dir.file(f));
            CHECK(WIFEXITED(st));
            CHECK(WEXITSTATUS(st) == EXIT_FAILURE);
        }
        dir.write("one.fq", "@a\nAC\n+\nII\n");
        dir.write("good_then_bad.fq", "@a\nAC\n+\nII\n@b\nACGT\n+\nII\n");
        const int drained = read_in_child(leg, dir.file("one.fq"), dir.file("good_then_bad.fq"));
        CHECK(WIFEXITED(drained));
        CHECK(WEXITSTATUS(drained) == EXIT_FAILURE);
        // The controls: well-formed files read to the end and exit 0.
        for (const char *f : ok_files) {
            CAPTURE(std::string(f));
            const int ok = read_in_child(leg, dir.file(f));
            CHECK(WIFEXITED(ok));
            CHECK(WEXITSTATUS(ok) == 0);
        }
        const int pair_ok = read_in_child(leg, dir.file("good.fq"), dir.file("good.fq"));
        CHECK(WIFEXITED(pair_ok));
        CHECK(WEXITSTATUS(pair_ok) == 0);
    }
}

TEST_CASE("mem_io reader: bseq_classify pairs adjacent equal names (-p)"
          * doctest::test_suite("unit/mem_io")) {
    // bwa's smart pairing: a record whose name equals the previous one's forms
    // a pair with it; everything else is single-end. Both lists keep input
    // order and the records' ids.
    const char *names[] = {"a", "a", "b", "c", "c", "d", "e", "e", "e", "f"};
    const int n = 10;
    std::vector<bseq1_t> seqs((size_t)n);
    for (int i = 0; i < n; ++i) {
        memset(&seqs[(size_t)i], 0, sizeof(bseq1_t));
        seqs[(size_t)i].name = const_cast<char *>(names[i]);
        seqs[(size_t)i].id = i;
    }
    int m[2] = {-1, -1};
    bseq1_t *sep[2] = {NULL, NULL};
    bseq_classify(n, seqs.data(), m, sep);
    REQUIRE(m[0] == 4);
    REQUIRE(m[1] == 6);
    std::string singles, pairs;
    for (int i = 0; i < m[0]; ++i) singles += std::string(sep[0][i].name) + std::to_string(sep[0][i].id) + " ";
    for (int i = 0; i < m[1]; ++i) pairs += std::string(sep[1][i].name) + std::to_string(sep[1][i].id) + " ";
    // e e e: the first two pair, the third is single.
    CHECK(singles == "b2 d5 e8 f9 ");
    CHECK(pairs == "a0 a1 c3 c4 e6 e7 ");
    free(sep[0]);
    free(sep[1]);
}
