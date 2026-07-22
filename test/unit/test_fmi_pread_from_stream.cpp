// test/unit/test_fmi_pread_from_stream.cpp
//
// Byte-exactness and stream-postcondition tests for the parallel index read.
//
// FMI_search::load_index reads each big index array by handing the stream's
// current offset to a pool of pread workers, then re-anchoring the FILE* past
// the array (src/FMI_search.cpp, fmi_pread_from_stream). Two things have to
// hold for the loaded index to be correct, and neither is visible in the
// aligner's output until an alignment silently goes wrong:
//
//   1. The concatenated chunks reproduce the file's bytes exactly, at every
//      worker count -- a chunk-offset or chunk-length slip corrupts the array
//      across an interior boundary while both ends still look right.
//   2. The stream lands exactly past the array, because the trailing
//      sentinel_index is still read sequentially with err_fread_noeof().
//
// The reads run against a synthetic temp file so no index build is needed.

#include "doctest/doctest.h"
#include "../../src/FMI_search.h"

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <unistd.h>   /* getpid -- per-process temp path */
#include <vector>

namespace {

// Deterministic byte pattern with a long period, so a chunk read at the wrong
// offset -- or a chunk of the wrong length -- cannot coincidentally match.
uint8_t pattern_byte(size_t i)
{
    return (uint8_t)((i * 1103515245u + 12345u) >> 16);
}

// Index of the first byte of `buf` that departs from the pattern starting at
// file offset `file_off`, or `len` if every byte matches. Returning the offset
// rather than asserting per byte keeps the suite's assertion count sane while
// still naming the chunk boundary a failure landed on.
size_t first_pattern_mismatch(const uint8_t *buf, size_t len, size_t file_off)
{
    for (size_t i = 0; i < len; ++i) {
        if (buf[i] != pattern_byte(file_off + i)) return i;
    }
    return len;
}

// A temp file holding `nbytes` of pattern_byte(), removed on destruction.
class PatternFile {
  public:
    explicit PatternFile(size_t nbytes) : path_(make_path())
    {
        FILE *fp = fopen(path_.c_str(), "wb");
        REQUIRE(fp != NULL);
        std::vector<uint8_t> buf(nbytes);
        for (size_t i = 0; i < nbytes; ++i) buf[i] = pattern_byte(i);
        if (nbytes > 0) REQUIRE(fwrite(buf.data(), 1, nbytes, fp) == nbytes);
        REQUIRE(fclose(fp) == 0);
    }

    ~PatternFile() { remove(path_.c_str()); }

    const char *path() const { return path_.c_str(); }

  private:
    static std::string make_path()
    {
        const char *dir = getenv("TMPDIR");
        std::string p = std::string(dir != NULL ? dir : "/tmp");
        if (!p.empty() && p[p.size() - 1] != '/') p += '/';
        // Distinct per process so concurrent test runs don't collide.
        char suffix[64];
        snprintf(suffix, sizeof(suffix), "bwa3_pread_test_%ld", (long)getpid());
        return p + suffix;
    }

    std::string path_;
};

}  // namespace

TEST_CASE("fmi_pread_from_stream: chunked read reproduces the file bytes") {
    // Spans several FMI_PREAD_MIN_CHUNK floors so worker counts above 1 are
    // actually honored, plus a ragged tail so the remainder distribution
    // (base + 1 for the first `rem` chunks) is exercised rather than an even
    // split that would hide an off-by-one in the chunk arithmetic.
    const size_t nbytes = 5 * FMI_PREAD_MIN_CHUNK + 4097;
    PatternFile file(nbytes);

    const int worker_counts[] = {1, 2, 3, 5, 8, 64};
    for (size_t k = 0; k < sizeof(worker_counts) / sizeof(worker_counts[0]); ++k) {
        const int nthreads = worker_counts[k];
        CAPTURE(nthreads);

        FILE *fp = fopen(file.path(), "rb");
        REQUIRE(fp != NULL);

        std::vector<uint8_t> dst(nbytes, 0);
        fmi_pread_from_stream(fp, dst.data(), nbytes, nthreads);

        CHECK(first_pattern_mismatch(dst.data(), nbytes, 0) == nbytes);

        REQUIRE(fclose(fp) == 0);
    }
}

TEST_CASE("fmi_pread_from_stream: resumes from and restores the stream position") {
    // Mirrors load_index's actual shape: a sequential header read, then the
    // parallel array read, then a sequential trailer read. If the FILE* is not
    // re-anchored exactly past the array, the trailer comes back wrong -- this
    // is what keeps sentinel_index correct.
    const size_t header  = sizeof(int64_t);
    const size_t array   = 2 * FMI_PREAD_MIN_CHUNK + 123;
    const size_t trailer = sizeof(int64_t);
    PatternFile file(header + array + trailer);

    const int worker_counts[] = {1, 2, 4};
    for (size_t k = 0; k < sizeof(worker_counts) / sizeof(worker_counts[0]); ++k) {
        const int nthreads = worker_counts[k];
        CAPTURE(nthreads);

        FILE *fp = fopen(file.path(), "rb");
        REQUIRE(fp != NULL);

        uint8_t head[sizeof(int64_t)];
        REQUIRE(fread(head, 1, header, fp) == header);
        CHECK(first_pattern_mismatch(head, header, 0) == header);

        std::vector<uint8_t> dst(array, 0);
        fmi_pread_from_stream(fp, dst.data(), array, nthreads);
        CHECK(first_pattern_mismatch(dst.data(), array, header) == array);

        // The stream must sit exactly at header + array, both by report...
        CHECK((size_t)ftello(fp) == header + array);
        // ...and by what the next sequential read actually returns.
        uint8_t tail[sizeof(int64_t)];
        REQUIRE(fread(tail, 1, trailer, fp) == trailer);
        CHECK(first_pattern_mismatch(tail, trailer, header + array) == trailer);

        REQUIRE(fclose(fp) == 0);
    }
}

TEST_CASE("fmi_pread_from_stream: zero-length read is a no-op") {
    // cp_occ/sa_* are never empty in a real index, but parallel_pread's
    // nbytes == 0 short-circuit must leave the stream untouched rather than
    // divide by a zero worker count.
    PatternFile file(FMI_PREAD_MIN_CHUNK);

    FILE *fp = fopen(file.path(), "rb");
    REQUIRE(fp != NULL);

    uint8_t head[4];
    REQUIRE(fread(head, 1, sizeof(head), fp) == sizeof(head));

    fmi_pread_from_stream(fp, NULL, 0, 4);
    CHECK((size_t)ftello(fp) == sizeof(head));

    uint8_t next[4];
    REQUIRE(fread(next, 1, sizeof(next), fp) == sizeof(next));
    CHECK(first_pattern_mismatch(next, sizeof(next), sizeof(head)) == sizeof(next));

    REQUIRE(fclose(fp) == 0);
}
