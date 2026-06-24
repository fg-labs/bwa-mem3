// Unit tests for mem_opt_fill_meth_mat under --meth-scoring.
//
// The bisulfite matrices mat_ot/mat_ob derive from the symmetric mat (diagonal
// +a, off-diagonal -b, N row/col 0) by freeing the conversion cell to +a:
//   OT frees mat[C][T] = mat[1*5+3]   (nt4 A=0,C=1,G=2,T=3,N=4)
//   OB frees mat[G][A] = mat[2*5+0]
// GENOMIC frees ONLY that cell (the mirror stays a -b mismatch -> variant-aware,
// rank-1 fast path). COLLAPSED also frees the MIRROR cell so C/T (and G/A) are
// interchangeable (bwameth-compatible, two freed cells).

#include "doctest/doctest.h"
#include "bwamem.h"

#include <cstdlib>

namespace {

// Build the meth matrices for a given scoring mode from a fresh default opt
// (a=1, b=4). Returns the heap opt; caller frees.
mem_opt_t *opt_for(int scoring)
{
    mem_opt_t *o = mem_opt_init();   // a=1, b=4, symmetric mat filled
    o->meth_scoring = scoring;
    mem_opt_fill_meth_mat(o);
    return o;
}

}  // namespace

TEST_CASE("genomic frees only the conversion cell; mirror stays a mismatch"
          * doctest::test_suite("unit/meth_scoring"))
{
    mem_opt_t *o = opt_for(MEM_METH_SCORING_GENOMIC);
    const int a = o->a, b = o->b;

    // OT: ref-C x read-T freed to +a; mirror ref-T x read-C stays -b.
    CHECK(o->mat_ot[1 * 5 + 3] == a);
    CHECK(o->mat_ot[3 * 5 + 1] == -b);
    // OB: ref-G x read-A freed; mirror ref-A x read-G stays -b.
    CHECK(o->mat_ob[2 * 5 + 0] == a);
    CHECK(o->mat_ob[0 * 5 + 2] == -b);
    // A real diagonal match and an unrelated off-diagonal are untouched.
    CHECK(o->mat_ot[0 * 5 + 0] == a);     // A/A match
    CHECK(o->mat_ot[0 * 5 + 1] == -b);    // ref-A x read-C real mismatch

    free(o);
}

TEST_CASE("collapsed frees both the conversion cell and its mirror"
          * doctest::test_suite("unit/meth_scoring"))
{
    mem_opt_t *o = opt_for(MEM_METH_SCORING_COLLAPSED);
    const int a = o->a;

    // OT: both ref-C x read-T and ref-T x read-C freed (C/T interchangeable).
    CHECK(o->mat_ot[1 * 5 + 3] == a);
    CHECK(o->mat_ot[3 * 5 + 1] == a);
    // OB: both ref-G x read-A and ref-A x read-G freed (G/A interchangeable).
    CHECK(o->mat_ob[2 * 5 + 0] == a);
    CHECK(o->mat_ob[0 * 5 + 2] == a);
    // Unrelated off-diagonal still a mismatch (e.g. ref-A x read-C).
    CHECK(o->mat_ot[0 * 5 + 1] == -o->b);

    free(o);
}

TEST_CASE("collapsed is the default scoring mode"
          * doctest::test_suite("unit/meth_scoring"))
{
    mem_opt_t *o = mem_opt_init();
    CHECK(o->meth_scoring == MEM_METH_SCORING_COLLAPSED);
    // The default-built matrices therefore have the mirror cell freed.
    CHECK(o->mat_ot[3 * 5 + 1] == o->a);
    free(o);
}
