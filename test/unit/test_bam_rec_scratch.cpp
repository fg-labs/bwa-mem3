// test/unit/test_bam_rec_scratch.cpp
//
// Unit tests for BamRecScratch's grow-only buffers (src/bam_rec_scratch.h).
//
// The scratch is shared by both BAM emitters (bam_writer.cpp, meth_bam.cpp) and
// is reused across every emitted record, so its growth rule has to hold for the
// sequence of sizes a real run produces, not just one call:
//
//   * grow-only  — a shorter record must NOT shrink the buffer, or the next long
//                  record silently reallocates on every emission and the
//                  optimization this scratch exists for disappears.
//   * stable     — when the request already fits, the SAME pointer comes back;
//                  callers write through it immediately after.
//   * bounded    — the size is checked for overflow before malloc, so a wrapped
//                  `n * sizeof(T)` cannot under-allocate a buffer the caller
//                  then writes `n` elements into.
//   * non-destructive on failure — an allocation that cannot be satisfied must
//                  leave the existing buffer and capacity intact rather than
//                  freeing a buffer that was merely too small.

#include "doctest/doctest.h"
#include "../../src/bam_rec_scratch.h"

#include <cstdint>
#include <cstring>

using bwamem3::BamRecScratch;

TEST_CASE("BamRecScratch: buffers start empty") {
    BamRecScratch bs;
    // Nothing is allocated until first grown -- one unused instance per thread
    // in either mode is supposed to cost nothing.
    CHECK(bs.ensure_cigar(0) == nullptr);
    CHECK(bs.ensure_seq(0) == nullptr);
    CHECK(bs.ensure_qual(0) == nullptr);
}

TEST_CASE("BamRecScratch: grows on demand and is writable to the full request") {
    BamRecScratch bs;

    uint32_t *cig = bs.ensure_cigar(8);
    REQUIRE(cig != nullptr);
    for (int i = 0; i < 8; ++i) cig[i] = (uint32_t)(i + 1);
    for (int i = 0; i < 8; ++i) CHECK(cig[i] == (uint32_t)(i + 1));

    char *seq = bs.ensure_seq(16);
    REQUIRE(seq != nullptr);
    memset(seq, 'A', 15);
    seq[15] = '\0';
    CHECK(strlen(seq) == 15);

    char *qual = bs.ensure_qual(16);
    REQUIRE(qual != nullptr);
    memset(qual, 30, 16);
    CHECK(qual[15] == 30);
}

TEST_CASE("BamRecScratch: a request that already fits returns the same pointer") {
    BamRecScratch bs;
    uint32_t *first = bs.ensure_cigar(64);
    REQUIRE(first != nullptr);

    // Equal to capacity, and below it: both must reuse, not reallocate.
    CHECK(bs.ensure_cigar(64) == first);
    CHECK(bs.ensure_cigar(1) == first);
    CHECK(bs.ensure_cigar(0) == first);
}

TEST_CASE("BamRecScratch: growth is monotonic across a realistic record sequence") {
    BamRecScratch bs;

    // A long record, then a short one, then a long one again. The middle request
    // must not shrink the buffer -- if it does, the third call reallocates.
    char *big = bs.ensure_seq(4096);
    REQUIRE(big != nullptr);
    memset(big, 'C', 4096);

    char *small = bs.ensure_seq(32);
    CHECK(small == big);            // reused, not shrunk

    char *again = bs.ensure_seq(4096);
    CHECK(again == big);            // still the original allocation

    // Only a genuinely larger request may move it, and the data it returns must
    // be writable to the new length.
    char *bigger = bs.ensure_seq(8192);
    REQUIRE(bigger != nullptr);
    memset(bigger, 'G', 8192);
    CHECK(bigger[8191] == 'G');
}

TEST_CASE("BamRecScratch: an overflowing element count is refused, not under-allocated") {
    BamRecScratch bs;

    // n * sizeof(uint32_t) would wrap; the guard must refuse instead of
    // allocating a short buffer that the caller then writes n elements into.
    CHECK(bs.ensure_cigar(SIZE_MAX) == nullptr);
    CHECK(bs.ensure_cigar(SIZE_MAX / sizeof(uint32_t) + 1) == nullptr);

    // The refusal must not have disturbed a healthy existing buffer.
    uint32_t *ok = bs.ensure_cigar(4);
    REQUIRE(ok != nullptr);
    ok[3] = 0xDEADBEEFu;
    CHECK(bs.ensure_cigar(SIZE_MAX) == nullptr);
    CHECK(bs.ensure_cigar(4) == ok);        // same buffer, capacity kept
    CHECK(ok[3] == 0xDEADBEEFu);            // contents untouched
}

TEST_CASE("BamRecScratch: a failed allocation keeps the previous buffer usable") {
    BamRecScratch bs;
    char *seq = bs.ensure_seq(128);
    REQUIRE(seq != nullptr);
    memset(seq, 'T', 128);

    // sizeof(char) == 1, so SIZE_MAX passes the overflow guard and reaches
    // malloc, which cannot satisfy it. That exercises the ALLOCATION-FAILURE path
    // specifically (the overflow path is covered above): the request is refused,
    // and -- the point of allocating before freeing -- the existing buffer must
    // survive rather than be destroyed on the way to failing.
    CHECK(bs.ensure_seq(SIZE_MAX) == nullptr);

    // Same pointer, capacity retained, contents intact. Pre-fix the buffer had
    // already been freed here and cap zeroed, so this read was a use-after-free.
    CHECK(bs.ensure_seq(128) == seq);
    CHECK(seq[127] == 'T');
}
