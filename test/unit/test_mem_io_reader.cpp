// test/unit/test_mem_io_reader.cpp
//
// The --legacy-reader FASTA/FASTQ batch reader (kseq over zlib's gzFile,
// bseq_read_orig) on plain, gzip (one member and several) and BGZF files
// written by the test, and bseq_classify, the -p smart-pairing split. The
// records are worked out by hand from the FASTA/FASTQ grammar kseq.h and
// src/fr_fastq.h document:
//
//   - a record starts at '>' or '@'; anything before the first header is
//     skipped, as are blank lines;
//   - the name runs to the first whitespace; the comment is the rest of that
//     line (a lone trailing '\r' removed), absent when empty;
//   - the sequence is every line up to the next header or '+', newlines (and
//     CRLF) removed, bases kept as given (case, IUPAC);
//   - FASTQ quality lines are read until they reach the sequence length, so a
//     quality line may start with '@' or '+';
//   - bwa's trim_readno drops a trailing "/<digit>" from a name longer than two
//     characters;
//   - a batch stops once its bases reach the chunk size AND it holds an even
//     number of records; paired files are read in lockstep, interleaved (end 1,
//     end 2, ...), and stop at the shorter file;
//   - a quality string of the wrong length is a fatal error.
//
// The default reader (fast_reader -> fr_fastq -> bseq_read_fast) cannot be
// linked into this binary (fast_reader.c needs zlib-ng's zng_* inflate, which
// test/Makefile does not link), so test/regression/mem_io_cli.sh holds it to
// the same grammar through `bwa-mem3 mem`, and to byte-identical SAM against
// this reader.

#include "doctest/doctest.h"
#include "mem_fixture.h"

#include "bwa.h"
#include "kseq.h"
#include "read_arena.h"

#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>
#include <zlib.h>

#include <csignal>
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
    bool operator==(const Rec &o) const {
        return name == o.name && comment == o.comment && seq == o.seq && qual == o.qual
               && l_seq == o.l_seq && id == o.id;
    }
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

void write_gzip(const std::string &path, const std::vector<std::string> &members) {
    // Each member is a complete gzip stream; concatenated they are one file.
    FILE *out = fopen(path.c_str(), "wb");
    REQUIRE(out != NULL);
    for (const std::string &m : members) {
        z_stream z;
        memset(&z, 0, sizeof(z));
        REQUIRE(deflateInit2(&z, 6, Z_DEFLATED, 15 + 16, 8, Z_DEFAULT_STRATEGY) == Z_OK);
        std::vector<unsigned char> buf(compressBound((uLong)m.size()) + 64);
        z.next_in = (Bytef *)m.data();
        z.avail_in = (uInt)m.size();
        z.next_out = buf.data();
        z.avail_out = (uInt)buf.size();
        REQUIRE(deflate(&z, Z_FINISH) == Z_STREAM_END);
        fwrite(buf.data(), 1, buf.size() - z.avail_out, out);
        deflateEnd(&z);
    }
    REQUIRE(fclose(out) == 0);
}

// BGZF (SAM spec section 4.1): gzip members of at most 64 KiB, each with a
// "BC" extra subfield holding the member's size - 1, then an empty member as
// the end-of-file marker. Written here with zlib's raw deflate, small blocks so
// the file has several.
void write_bgzf(const std::string &path, const std::string &text) {
    FILE *out = fopen(path.c_str(), "wb");
    REQUIRE(out != NULL);
    const size_t step = 37;
    for (size_t i = 0; i <= text.size(); i += step) {
        // The final pass (i == size, n == 0) writes the empty EOF block.
        const size_t n = i + step <= text.size() ? step : text.size() - i;
        z_stream z;
        memset(&z, 0, sizeof(z));
        REQUIRE(deflateInit2(&z, 6, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY) == Z_OK);
        std::vector<unsigned char> body(compressBound((uLong)n) + 16);
        z.next_in = (Bytef *)(text.data() + i);
        z.avail_in = (uInt)n;
        z.next_out = body.data();
        z.avail_out = (uInt)body.size();
        REQUIRE(deflate(&z, Z_FINISH) == Z_STREAM_END);
        const size_t clen = body.size() - z.avail_out;
        deflateEnd(&z);
        const size_t bsize = 18 + clen + 8;  // header, deflate data, CRC32 + ISIZE
        const unsigned char hdr[18] = {0x1f, 0x8b, 8, 4, 0, 0, 0, 0, 0, 0xff, 6, 0, 'B', 'C', 2, 0,
                                       (unsigned char)((bsize - 1) & 0xff),
                                       (unsigned char)((bsize - 1) >> 8)};
        const uLong crc = crc32(crc32(0L, Z_NULL, 0), (const Bytef *)(text.data() + i), (uInt)n);
        unsigned char tail[8];
        for (int k = 0; k < 4; ++k) {
            tail[k] = (unsigned char)((crc >> (8 * k)) & 0xff);
            tail[4 + k] = (unsigned char)((n >> (8 * k)) & 0xff);
        }
        fwrite(hdr, 1, sizeof(hdr), out);
        fwrite(body.data(), 1, clen, out);
        fwrite(tail, 1, sizeof(tail), out);
        if (n == 0) break;
    }
    REQUIRE(fclose(out) == 0);
}

// The reader must give `expected` for `text` written as plain, gzip (one member
// and split into three) and BGZF.
void check_all_codecs(const std::string &text, int64_t chunk,
                      const std::vector<std::vector<Rec> > &expected) {
    bwa_tests::TempDir dir;
    dir.write("in.fq", text);
    write_gzip(dir.file("in.fq.gz"), std::vector<std::string>(1, text));
    const size_t a = text.size() / 3, b = 2 * text.size() / 3;
    write_gzip(dir.file("multi.fq.gz"), {text.substr(0, a), text.substr(a, b - a), text.substr(b)});
    write_bgzf(dir.file("in.fq.bgz"), text);
    const char *files[] = {"in.fq", "in.fq.gz", "multi.fq.gz", "in.fq.bgz"};
    for (const char *f : files) {
        CAPTURE(f);
        CHECK(show(read_legacy(dir.file(f), "", chunk)) == show(expected));
    }
}

const int64_t kWholeFile = 1 << 30;

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

TEST_CASE("mem_io reader: FASTA and multi-line records and CRLF and stray lines"
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
    }
    SUBCASE("text before the first record and blank lines between records") {
        check_all_codecs("\n\njunk line\n@a\nAC\n+\nII\n\n\n@b\nGT\n+\nJJ", kWholeFile,
                         {{rec("a", NULL, "AC", "II", 0), rec("b", NULL, "GT", "JJ", 1)}});
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
    check_all_codecs(text, 25, b25);
    // chunk 20: 2 records reach it exactly -> batches of 2, 2, 2, 1.
    std::vector<std::vector<Rec> > b20(4);
    for (int i = 0; i < 7; ++i) {
        Rec r = all[(size_t)i];
        r.id = i % 2;
        b20[(size_t)(i / 2)].push_back(r);
    }
    check_all_codecs(text, 20, b20);
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

    // Chunk 5: the first pair holds 8 bases, so each batch is one pair.
    std::vector<std::vector<Rec> > one_pair_each(3);
    for (size_t i = 0; i < 6; ++i) {
        Rec r = pairs[0][i];
        r.id = (int)(i % 2);
        one_pair_each[i / 2].push_back(r);
    }
    CHECK(show(read_legacy(dir.file("r1.fq"), dir.file("r2.fq"), 5)) == show(one_pair_each));

    // A shorter second file ends the input at its last record, either way round.
    {
        bwa_tests::QuietStderr quiet;  // the reader warns about the shorter file
        const std::vector<std::vector<Rec> > first_pair{
            {rec("p1", "x", "AAAA", "IIII", 0), rec("p1", NULL, "TTTT", "5555", 1)}};
        CHECK(show(read_legacy(dir.file("r1.fq"), dir.file("short.fq"), kWholeFile))
              == show(first_pair));
        const std::vector<std::vector<Rec> > first_pair_swapped{
            {rec("p1", NULL, "TTTT", "5555", 0), rec("p1", "x", "AAAA", "IIII", 1)}};
        CHECK(show(read_legacy(dir.file("short.fq"), dir.file("r1.fq"), kWholeFile))
              == show(first_pair_swapped));
    }
}

namespace {

// Runs the reader over a file in a forked child (stderr silenced, no core
// file) and returns its wait status: a reader that rejects the input exits
// through err_fatal with EXIT_FAILURE.
int read_in_child(const std::string &path) {
    fflush(stdout);
    fflush(stderr);
    const pid_t pid = fork();
    if (pid == 0) {
        signal(SIGABRT, SIG_DFL);
        struct rlimit no_core = {0, 0};
        setrlimit(RLIMIT_CORE, &no_core);
        if (freopen("/dev/null", "w", stderr) == NULL) _exit(3);
        gzFile g = gzopen(path.c_str(), "r");
        if (g == NULL) _exit(4);
        kseq_t *k = kseq_init(g);
        int n = 0;
        int64_t size = 0;
        read_arena_t *arena = NULL;
        while (bseq_read_orig(1 << 20, &n, k, NULL, &size, &arena) != NULL && n > 0) arena = NULL;
        _exit(0);
    }
    REQUIRE(pid > 0);
    int status = 0;
    REQUIRE(waitpid(pid, &status, 0) == pid);
    return status;
}

}  // namespace

TEST_CASE("mem_io reader: a quality string of the wrong length is fatal"
          * doctest::test_suite("unit/mem_io")) {
    bwa_tests::TempDir dir;
    // A shorter quality at the end of input, and a '+' line with nothing after it.
    dir.write("short_qual.fq", "@a\nACGT\n+\nII\n");
    dir.write("no_qual.fq", "@a\nAC\n+\nII\n@b\nACGT\n+");
    dir.write("good.fq", "@a\nAC\n+\nII\n@b\nACGT\n+\nIIII\n");
    const char *bad[] = {"short_qual.fq", "no_qual.fq"};
    for (const char *f : bad) {
        CAPTURE(f);
        const int st = read_in_child(dir.file(f));
        CHECK(WIFEXITED(st));
        CHECK(WEXITSTATUS(st) == EXIT_FAILURE);
    }
    // The control: a well-formed file reads to the end and exits 0.
    const int ok = read_in_child(dir.file("good.fq"));
    CHECK(WIFEXITED(ok));
    CHECK(WEXITSTATUS(ok) == 0);
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

    // An empty batch yields two empty lists.
    bseq_classify(0, NULL, m, sep);
    CHECK(m[0] == 0);
    CHECK(m[1] == 0);
}
