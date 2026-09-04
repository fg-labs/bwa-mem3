// test/bam_rec_scratch_test.cpp -- pins the geometric-growth contract of
// BamRecScratch::grow().
//
// grow() doubles the capacity (floored at the requested n) instead of the old
// exact-fit cap = n. The observable contract is buffer stability: after a grow,
// a subsequent request up to the doubled capacity reuses the SAME buffer with
// no reallocation. The old exact-fit implementation would reallocate on the
// very next larger request, so this test fails against it and passes against
// the geometric version -- giving the change actual regression protection.

#include "bam_rec_scratch.h"
#include "grow_capacity.h"

#include <cstdio>
#include <cstdlib>

#define CHECK(cond) do {                                                    \
    if (!(cond)) {                                                          \
        std::fprintf(stderr, "CHECK failed at %s:%d: %s\n",                 \
                     __FILE__, __LINE__, #cond);                           \
        std::abort();                                                      \
    }                                                                      \
} while (0)

using bwamem3::BamRecScratch;

int main() {
    BamRecScratch s;

    // First allocation: cap 0 -> want max(0*2, 100) = 100.
    char *p1 = s.ensure_seq(100);
    CHECK(p1 != nullptr);

    // Grow: n=150 > cap=100 -> want = max(200, 150) = 200. New buffer.
    char *p2 = s.ensure_seq(150);
    CHECK(p2 != nullptr);
    CHECK(s.seq_cap == 200);

    // Geometric growth: cap is now 200, so a request up to 200 reuses the same
    // buffer. The old exact-fit impl (cap=150 here) would reallocate.
    CHECK(s.ensure_seq(200) == p2);
    CHECK(s.ensure_seq(1)   == p2);   // shrink request never reallocates

    // Beyond the doubled cap: n=500 > cap=200, want=max(400,500)=500. Grows.
    char *p3 = s.ensure_seq(500);
    CHECK(p3 != nullptr);
    // n=700 > cap=500, want=1000. Grows; then up to 1000 is stable.
    char *p4 = s.ensure_seq(700);
    CHECK(p4 != nullptr);
    CHECK(s.ensure_seq(1000) == p4);

    // The other two buffers share the same grow(): exercise a mid-range request
    // that forces a doubling grow, then confirm the doubled capacity and reuse.
    uint32_t *c1 = s.ensure_cigar(64);
    CHECK(c1 != nullptr);
    uint32_t *c2 = s.ensure_cigar(96);   // 96 > cap=64 -> want = max(128,96) = 128
    CHECK(c2 != nullptr);
    CHECK(s.cigar_cap == 128);
    CHECK(s.ensure_cigar(128) == c2);    // within doubled cap, stable

    char *q1 = s.ensure_qual(64);
    CHECK(q1 != nullptr);
    char *q2 = s.ensure_qual(96);        // 96 > cap=64 -> want = max(128,96) = 128
    CHECK(q2 != nullptr);
    CHECK(s.qual_cap == 128);
    CHECK(s.ensure_qual(128) == q2);     // within doubled cap, stable

    // grow_capacity() itself: no BamRecScratch buffer ever grows anywhere near
    // SIZE_MAX, so the overflow-clamp branch (elem_size > 1 && the doubled
    // value would overflow once scaled) is unreachable through BamRecScratch
    // alone. Call the pure function directly -- no malloc needed -- mirroring
    // the direct-call pattern the seqbuf_grow_capacity() guard test uses for
    // its own overflow branch.
    CHECK(grow_capacity(0,   100, 1) == 100);   // first alloc: floor at need
    CHECK(grow_capacity(100, 50,  1) == 200);   // shrink request still doubles

    const size_t huge_cur = SIZE_MAX / 4;
    // elem_size == 1: no scale-overflow check applies; ordinary doubling.
    CHECK(grow_capacity(huge_cur, 100, 1) == huge_cur * 2);
    // elem_size > 1: the doubled value is a valid size_t but would overflow
    // once scaled by elem_size -- clamp back to need.
    CHECK(grow_capacity(huge_cur, 100, 8) == 100);

    std::fprintf(stderr, "bam_rec_scratch_test OK\n");
    return 0;
}
