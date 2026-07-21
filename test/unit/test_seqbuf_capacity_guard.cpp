// test/unit/test_seqbuf_capacity_guard.cpp
//
// Regression test for silent int32 truncation of seqBuf* offsets.
//
// The extension buffers (seqBufLeftRef/Qer, seqBufRightRef/Qer) grow by
// doubling as a read block accumulates per-seed spans. The offsets into those
// buffers are handed to the SW kernels through SeqPair.idr/.idq, which are
// int32_t (src/bandedSWA.h). The accumulators feeding them
// (leftRefOffset/leftQerOffset, bwamem.cpp) are int64_t.
//
// Once a buffer doubles past INT32_MAX the narrowing assignment silently
// produces a negative index, and `seqBufLeftRef + sp.idr` then reads and
// writes outside the allocation. Long reads reach this: the ref buffer starts
// at BATCH_SIZE * SEEDS_PER_READ * MAX_SEQ_LEN_REF and a HiFi/ONT block needs
// far more span per seed than the 256-byte-per-seed model assumes.
//
// seqbuf_grow_capacity() is the single choke point: it refuses to report a
// capacity that cannot be addressed by an int32 offset.

#include "doctest/doctest.h"
#include "../../src/bwamem.h"
#include "../../src/bandedSWA.h"
#include "../../src/macro.h"

#include <cstdint>

TEST_CASE("seqbuf_grow_capacity: ordinary growth doubles") {
    CHECK(seqbuf_grow_capacity(1024) == 2048);
    CHECK(seqbuf_grow_capacity(131072000LL) == 262144000LL);  // initial ref buf
}

TEST_CASE("seqbuf_grow_capacity: refuses growth past the int32 offset range") {
    // Largest capacity that still leaves every in-range offset representable.
    const int64_t kMax = (int64_t)INT32_MAX;

    // Doubling from just under half the range is still fine.
    CHECK(seqbuf_grow_capacity(kMax / 4) == kMax / 4 * 2);

    // Doubling past INT32_MAX must be refused, not silently returned.
    CHECK(seqbuf_grow_capacity((int64_t)1 << 30) == SEQBUF_CAPACITY_OVERFLOW);
    CHECK(seqbuf_grow_capacity(kMax) == SEQBUF_CAPACITY_OVERFLOW);
}

TEST_CASE("seqbuf_grow_capacity: SeqPair offset fields are the constraint") {
    // The guard exists because these fields are 32-bit. If they are ever
    // widened, the guard threshold must be revisited (and this test updated)
    // rather than left silently over-restrictive.
    SeqPair sp;
    CHECK(sizeof(sp.idr) == 4);
    CHECK(sizeof(sp.idq) == 4);
}

TEST_CASE("seqbuf_grow_capacity: the overflow is reachable from real sizes") {
    // Documents reachability rather than asserting a policy: starting from the
    // shipped initial ref-buffer size, count doublings until the guard trips.
    // This is what a long-read block actually walks through.
    int64_t cap = (int64_t)BATCH_SIZE * SEEDS_PER_READ * MAX_SEQ_LEN_REF;
    int doublings = 0;
    while (cap != SEQBUF_CAPACITY_OVERFLOW) {
        cap = seqbuf_grow_capacity(cap);
        ++doublings;
        REQUIRE(doublings < 64);  // must terminate
    }
    // Reached in a handful of doublings -- not a theoretical concern.
    CHECK(doublings <= 6);
}
