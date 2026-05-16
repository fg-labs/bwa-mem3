// test/unit/test_mm_realloc_byte_count.cpp
//
// Regression test for the _mm_realloc copy-length bug: the historical body
// allocates `nsize * dsize` bytes but copies only `csize` bytes, so any caller
// passing dsize > 1 with csize meaning "elements" gets a truncated copy
// (csize/dsize elements survive, the rest is uninitialized scratch from
// _mm_malloc). Real callers in bwamem.cpp pass sizeof(int64_t) and
// sizeof(SMEM) — both > 1.

#include "doctest/doctest.h"
#include "../../src/bwamem.h"

#include <cstdint>
#include <cstring>

TEST_CASE("_mm_realloc: int64 buffer growth preserves every element") {
    // csize/nsize are element counts in the caller convention (see
    // bwamem.cpp:903 sa_coord realloc). dsize is sizeof(int64_t) = 8.
    const int64_t csize = 4;     // 4 int64 elements pre-grow
    const int64_t nsize = 8;     // 8 int64 elements post-grow
    const int16_t dsize = sizeof(int64_t);

    int64_t *buf = (int64_t *)_mm_malloc((size_t)csize * dsize, 64);
    REQUIRE(buf != nullptr);
    for (int64_t i = 0; i < csize; ++i) buf[i] = 0x1234567800000000LL + i;

    int64_t *grown = (int64_t *)_mm_realloc(buf, csize, nsize, dsize);
    REQUIRE(grown != nullptr);

    // Every pre-existing element must survive the copy. Pre-fix, this checked
    // only buf[0]: the body copied csize=4 BYTES of an 8-byte element, so
    // buf[0]'s low 4 bytes were correct and the high 4 bytes were scratch.
    for (int64_t i = 0; i < csize; ++i) {
        CHECK(grown[i] == 0x1234567800000000LL + i);
    }

    _mm_free(grown);
}

TEST_CASE("_mm_realloc: lazy-init from NULL with csize=0 is safe") {
    // Mirrors bwamem.cpp:1061 first-grow path: matchArray[tid] starts NULL
    // with wsize_mem[tid]=0, then grows to a positive count. _mm_realloc
    // must not segfault on memcpy(nptr, NULL, 0) and must return a usable
    // buffer.
    const int64_t csize = 0;
    const int64_t nsize = 8;
    const int16_t dsize = sizeof(int64_t);

    int64_t *grown = (int64_t *)_mm_realloc(nullptr, csize, nsize, dsize);
    REQUIRE(grown != nullptr);
    // Buffer must be writable across the full nsize range.
    for (int64_t i = 0; i < nsize; ++i) grown[i] = i;
    for (int64_t i = 0; i < nsize; ++i) CHECK(grown[i] == i);

    _mm_free(grown);
}

TEST_CASE("_mm_realloc: byte-buffer growth (dsize=1) preserves elements") {
    // dsize=1 hits an arm of the bug-free path (csize == csize*dsize) — sanity
    // check that the fix doesn't regress the seqBuf* callsites that pass
    // sizeof(uint8_t).
    const int64_t csize = 16;
    const int64_t nsize = 32;
    const int16_t dsize = sizeof(uint8_t);

    uint8_t *buf = (uint8_t *)_mm_malloc((size_t)csize * dsize, 64);
    REQUIRE(buf != nullptr);
    for (int64_t i = 0; i < csize; ++i) buf[i] = (uint8_t)(0xA0 + i);

    uint8_t *grown = (uint8_t *)_mm_realloc(buf, csize, nsize, dsize);
    REQUIRE(grown != nullptr);
    for (int64_t i = 0; i < csize; ++i) CHECK(grown[i] == (uint8_t)(0xA0 + i));

    _mm_free(grown);
}
