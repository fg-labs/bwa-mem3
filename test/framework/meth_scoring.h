// test/framework/meth_scoring.h
//
// The production --meth OT / OB scoring matrices for test fixtures.
//
// Built through the same two calls the CLI makes -- bwa_fill_scmat, then
// mem_opt_fill_meth_mat -- so a test's oracle takes the matrix production
// uses. Restating the freed-cell layout in a test instead would stop testing
// mem_opt_fill_meth_mat and start testing a copy of it. Needs libbwa.a, which
// the unit and integration binaries link.

#ifndef BWA_TESTS_METH_SCORING_H
#define BWA_TESTS_METH_SCORING_H

#include "scoring.h"

namespace bwa_tests {

// The OT (ot = true: reference C / read T freed) or OB (reference G / read A
// freed) matrix for --meth-scoring `scoring` (MEM_METH_SCORING_*) at the given
// match / mismatch magnitudes.
ScoringMatrix meth_scoring_matrix(int scoring, bool ot, int match, int mismatch);

} // namespace bwa_tests

#endif
