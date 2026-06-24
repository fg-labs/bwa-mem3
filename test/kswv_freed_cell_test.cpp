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
#include "ksw.h"
#include "kswv.h"
#include <cstring>
#include <vector>
#include <cstdint>

static void fill_sym(int8_t m[25], int a, int b) {
    for (int i=0;i<5;i++) for (int j=0;j<5;j++)
        m[i*5+j] = (i==4||j==4) ? 0 : (i==j ? (int8_t)a : (int8_t)-b);
}

// Run getScores8 over a single ref/read pair through the batched NEON 8-bit
// kernel and return aln[0].score. `mat25 == nullptr` constructs the symmetric
// kernel (no freed cell); a meth matrix sets has_freed and applies the rank-1
// override. Ref/read are base-encoded (A=0,C=1,G=2,T=3). The pair is laid out
// so the optimal local alignment is the full diagonal — every column is a
// match except the deliberately planted converted base.
static int score8(const std::vector<uint8_t>& ref,
                  const std::vector<uint8_t>& read,
                  int a, int b, const int8_t* mat25) {
    const int32_t maxRefLen = 256, maxQerLen = 256;
    SeqPair sp = {};
    sp.idr = 0; sp.idq = 0; sp.id = 0;
    sp.len1 = (int)ref.size();
    sp.len2 = (int)read.size();
    sp.h0 = 0;                 // no KSW_XSUBO/XSTOP: plain best-score scan
    sp.seqid = 0; sp.regid = 0;

    std::vector<SeqPair> pairs(SIMD_WIDTH8 + 1);
    pairs[0] = sp;

    std::vector<kswr_t> aln(SIMD_WIDTH8 + 1, g_defr);
    // Mutable copies for the (non-const) buffer pointers getScores8 expects.
    std::vector<uint8_t> rbuf = ref, qbuf = read;

    kswv k(6, 1, 6, 1, (int8_t)a, (int8_t)-b, 1, maxRefLen, maxQerLen, mat25);
    // phase 0 fills aln[regid].score with the de-biased best SW score.
    k.getScores8(pairs.data(), rbuf.data(), qbuf.data(), aln.data(), 1, 1, 0);
    return aln[0].score;
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

// Task 4: the NEON kernel must actually APPLY the freed cell. Build a clean
// matching alignment with exactly one ref-C×read-T (OT) cell, score it through
// getScores8 under the symmetric vs OT matrix, and assert the OT score exceeds
// the symmetric score by exactly (a + b): the converted cell flips from a
// mismatch (-b) to a match (+a). OB is the ref-G×read-A analog.
TEST_CASE("kswv NEON 8-bit applies the OT freed cell (ref-C x read-T)") {
    const int a = 1, b = 4;
    // 20-base identical diagonal, then plant ref=C(1)/read=T(3) at one column.
    std::vector<uint8_t> ref  = {0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3};
    std::vector<uint8_t> read = ref;
    const int p = 10;                 // converted column
    ref[p]  = 1;                      // ref C
    read[p] = 3;                      // read T (bisulfite C->T on read)

    int8_t sym[25]; fill_sym(sym, a, b);
    int8_t ot[25];  fill_sym(ot, a, b); ot[1*5+3] = (int8_t)a;  // free C×T -> match

    int s_sym = score8(ref, read, a, b, sym);
    int s_ot  = score8(ref, read, a, b, ot);
    // exactly one freed cell on the optimal path: delta == a + b
    CHECK(s_ot - s_sym == a + b);
}

TEST_CASE("kswv NEON 8-bit applies the OB freed cell (ref-G x read-A)") {
    const int a = 1, b = 4;
    std::vector<uint8_t> ref  = {0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3};
    std::vector<uint8_t> read = ref;
    const int p = 10;
    ref[p]  = 2;                      // ref G
    read[p] = 0;                      // read A (bisulfite G->A on read, OB strand)

    int8_t sym[25]; fill_sym(sym, a, b);
    int8_t ob[25];  fill_sym(ob, a, b); ob[2*5+0] = (int8_t)a;  // free G×A -> match

    int s_sym = score8(ref, read, a, b, sym);
    int s_ob  = score8(ref, read, a, b, ob);
    CHECK(s_ob - s_sym == a + b);
}

// Two freed cells on the optimal path must each contribute (a + b). This also
// guards against an override that only fires once per row/strip.
TEST_CASE("kswv NEON 8-bit applies multiple OT freed cells additively") {
    const int a = 1, b = 4;
    std::vector<uint8_t> ref  = {0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3};
    std::vector<uint8_t> read = ref;
    ref[6]  = 1; read[6]  = 3;        // OT cell #1 (different rows)
    ref[14] = 1; read[14] = 3;        // OT cell #2

    int8_t sym[25]; fill_sym(sym, a, b);
    int8_t ot[25];  fill_sym(ot, a, b); ot[1*5+3] = (int8_t)a;

    int s_sym = score8(ref, read, a, b, sym);
    int s_ot  = score8(ref, read, a, b, ot);
    CHECK(s_ot - s_sym == 2 * (a + b));
}
