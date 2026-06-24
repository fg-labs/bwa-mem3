// Issue 173 / Task 2: mat-aware make_kswv freed-cell detection.
//
// The kswv factory learns a per-strand asymmetric (rank-1 freed-cell)
// substitution for bisulfite (--meth). A symmetric matrix (or nullptr)
// constructs the existing kernel unchanged; a matrix that frees exactly one
// ordered off-diagonal cell to a match (OT: ref-C/read-T, OB: ref-G/read-A)
// records that freed cell and is handled by the kernel; any other asymmetric
// matrix (>=2 freed cells, changed diagonal, non-match freed value) routes to
// the scalar fallback (needsScalar() == true).
//
// This test only exercises the ctor-side detection (needsScalar); the kernels
// still ignore the freed cell at this task, so non-meth output stays
// byte-identical.

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"
#include "kswv.h"
#include <cstring>

static void fill_sym(int8_t m[25], int a, int b) {
    for (int i=0;i<5;i++) for (int j=0;j<5;j++)
        m[i*5+j] = (i==4||j==4) ? 0 : (i==j ? (int8_t)a : (int8_t)-b);
}
TEST_CASE("kswv detects OT freed cell (ref-C x read-T)") {
    int8_t m[25]; fill_sym(m,1,4); m[1*5+3] = 1;          // free ref-C x read-T to match
    auto k = make_kswv(6,1,6,1, 1,-4, 1, 256,256, m);
    CHECK(k->needsScalar() == false);                      // rank-1 is handled
}
TEST_CASE("kswv detects OB freed cell (ref-G x read-A)") {
    int8_t m[25]; fill_sym(m,1,4); m[2*5+0] = 1;
    auto k = make_kswv(6,1,6,1, 1,-4, 1, 256,256, m);
    CHECK(k->needsScalar() == false);
}
TEST_CASE("kswv symmetric matrix is not a freed-cell case") {
    int8_t m[25]; fill_sym(m,1,4);
    auto k = make_kswv(6,1,6,1, 1,-4, 1, 256,256, m);
    CHECK(k->needsScalar() == false);                      // symmetric ⇒ normal kernel
}
TEST_CASE("kswv non-rank1 asymmetric falls back to scalar") {
    int8_t m[25]; fill_sym(m,1,4); m[1*5+3]=1; m[2*5+0]=1; // two freed cells
    auto k = make_kswv(6,1,6,1, 1,-4, 1, 256,256, m);
    CHECK(k->needsScalar() == true);
}
