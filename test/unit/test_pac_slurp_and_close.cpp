// test/unit/test_pac_slurp_and_close.cpp
//
// Byte-exactness and postcondition tests for pac_slurp_and_close (src/
// read_index_ele.cpp), the shared .pac slurp used by BOTH the seed-index loader
// (indexEle::bwa_idx_load_ele) and the --meth original-reference loader
// (meth_orig_ref_load_handles). The underlying fmi_pread_from_stream splitting
// is covered by test_fmi_pread_from_stream; this pins the wrapper the two call
// sites actually invoke:
//
//   1. the whole file is read into `dst` exactly, at both a serial worker count
//      and a multi-worker count large enough to force a real chunk split
//      (>= FMI_PREAD_MIN_CHUNK per worker) -- the property the aligner relies on
//      for byte-identical output regardless of -t / BWA3_LOAD_THREADS; and
//   2. the stream is closed and the caller's FILE* is NULLed, so the loaders'
//      subsequent `fp_pac == NULL` bookkeeping holds.
//
// Runs against a synthetic temp file so no index build is needed.

#include "doctest/doctest.h"
#include "../../src/read_index_ele.h"   /* pac_slurp_and_close */
#include "../../src/FMI_search.h"       /* FMI_PREAD_MIN_CHUNK */

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <unistd.h>   /* getpid -- per-process temp path */
#include <vector>

namespace {

// Deterministic long-period byte pattern: a chunk read at the wrong offset or
// of the wrong length cannot coincidentally match.
uint8_t pattern_byte(size_t i)
{
    return (uint8_t)((i * 1103515245u + 12345u) >> 16);
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
    const std::string &path() const { return path_; }

  private:
    static std::string make_path()
    {
        static unsigned seq = 0;   // unique per file within this process
        return "/tmp/bwa_mem3_pac_slurp_test." + std::to_string(getpid()) + "." +
               std::to_string(seq++);
    }
    std::string path_;
};

// Slurp `nbytes` from a fresh pattern file with `workers`, asserting the buffer
// reproduces the pattern exactly and the stream handle was NULLed.
void check_slurp(size_t nbytes, int workers)
{
    PatternFile file(nbytes);
    FILE *fp = fopen(file.path().c_str(), "rb");
    REQUIRE(fp != NULL);

    std::vector<uint8_t> dst(nbytes ? nbytes : 1, 0xAB);
    pac_slurp_and_close(&fp, dst.data(), (int64_t)nbytes, workers);

    CHECK(fp == NULL);   // stream closed and pointer cleared for the caller
    size_t mismatch = nbytes;
    for (size_t i = 0; i < nbytes; ++i) {
        if (dst[i] != pattern_byte(i)) { mismatch = i; break; }
    }
    CHECK(mismatch == nbytes);
}

}  // namespace

TEST_CASE("pac_slurp_and_close: serial read reproduces the file exactly"
          * doctest::test_suite("unit/pac_slurp_and_close")) {
    check_slurp(4096, /*workers=*/1);
}

TEST_CASE("pac_slurp_and_close: multi-worker read is byte-identical to serial"
          * doctest::test_suite("unit/pac_slurp_and_close")) {
    // Large enough that fmi_pread_worker_count actually hands out >1 chunk at
    // this worker count, exercising an interior chunk boundary.
    const size_t n = (size_t)FMI_PREAD_MIN_CHUNK * 4 + 12345;
    check_slurp(n, /*workers=*/4);
}

TEST_CASE("pac_slurp_and_close: sub-floor size falls back to one worker cleanly"
          * doctest::test_suite("unit/pac_slurp_and_close")) {
    // Below FMI_PREAD_MIN_CHUNK a >1 request collapses to a single worker; the
    // result must still be exact and the handle still NULLed.
    check_slurp(1000, /*workers=*/8);
}
