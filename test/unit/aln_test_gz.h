// gzip FASTQ inputs whose read fails partway, for the aln driver's read-error tests
// (test_aln_se.cpp, test_aln_pe.cpp), and an oracle for where a serial reader stops on them.

#ifndef BWA_MEM3_TEST_ALN_TEST_GZ_H
#define BWA_MEM3_TEST_ALN_TEST_GZ_H

#include <cstdio>
#include <string>
#include <vector>

#include <zlib.h>

#include "doctest/doctest.h"

namespace aln_test {

// Writes `reads`, named `prefix`0, `prefix`1, ..., as one gzip member, then a second member whose
// deflate data is invalid (a block of the reserved type), so gzread fails once the first member's
// reads are read.
inline std::string write_fastq_gz_corrupt_tail(const std::string &path, const std::vector<std::string> &reads,
                                               const std::string &prefix)
{
    gzFile gz = gzopen(path.c_str(), "wb");
    REQUIRE(gz != nullptr);
    for (size_t i = 0; i < reads.size(); ++i) {
        std::string rec = "@" + prefix + std::to_string(i) + "\n" + reads[i] + "\n+\n" + std::string(reads[i].size(), 'I') + "\n";
        REQUIRE(gzwrite(gz, rec.data(), (unsigned)rec.size()) == (int)rec.size());
    }
    REQUIRE(gzclose(gz) == Z_OK);
    const unsigned char bad[] = {0x1f, 0x8b, 0x08, 0, 0, 0, 0, 0, 0, 0x03, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff};
    FILE *fp = fopen(path.c_str(), "ab");
    REQUIRE(fp != nullptr);
    REQUIRE(fwrite(bad, 1, sizeof bad, fp) == sizeof bad);
    REQUIRE(fclose(fp) == 0);
    return path;
}

// The number of reads a serial loop of batches of `batch` reads writes before the read that fails:
// whole batches of the complete four-line records that gzread returns before its error. Read with
// zlib directly, as kseq calls it (16384-byte reads, src/kseq.h's KSEQ_INIT2: zlib drops the bytes
// of the call that fails, so the buffer size sets where the reads end), not with bwa_read_seq, so
// that a reader that drops a batch it should keep cannot shrink the expected count with it.
inline int gz_reads_before_error(const std::string &path, int batch)
{
    gzFile gz = gzopen(path.c_str(), "rb");
    REQUIRE(gz != nullptr);
    std::vector<char> buf(16384);
    long n_lines = 0;
    int got;
    while ((got = gzread(gz, buf.data(), (unsigned)buf.size())) > 0)
        for (int i = 0; i < got; ++i) n_lines += buf[i] == '\n';
    REQUIRE(got < 0); // the corrupt member fails the read, not the end of the file
    gzclose(gz);
    const long n_reads = n_lines / 4;
    return (int)(n_reads / batch * batch);
}

} // namespace aln_test

#endif
