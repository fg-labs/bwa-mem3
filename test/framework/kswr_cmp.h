// test/framework/kswr_cmp.h
//
// Comparators for kswr_t results produced by the scalar reference
// (ksw_align2) vs the batched SIMD kernel (kswv::getScores8). These
// embed the semantic exceptions the production code tolerates so that
// tests aren't tripped by known-unimportant discrepancies (e.g. scalar
// short-circuiting phase 2 when score < minsc).

#ifndef BWA_TESTS_KSWR_CMP_H
#define BWA_TESTS_KSWR_CMP_H

#include "ksw.h"
#include "ksw_runner.h"

namespace bwa_tests {

// Primary correctness bar: kernels must produce identical alignment scores.
bool kswr_score_eq(const kswr_t &a, const kswr_t &b);

// Coordinate equality (tb/qb). Phase-1 output.
// Two scalar-side exceptions match production behavior:
//   1. scalar.score == 0 — degenerate input, skip.
//   2. scalar.score < min_seed_len — scalar short-circuits phase 2
//      (ksw.cpp:374), leaving tb/qb at -1. Production's filter
//      `aln.score >= opt->min_seed_len` drops these pairs before tb/qb is
//      consulted, so the comparison is meaningless.
bool kswr_coords_eq(const kswr_t &scalar, const kswr_t &batched,
                    int min_seed_len = DEFAULT_MIN_SEED_LEN);

// Suboptimal-score equality (score2). KSW_XSUBO output. Scalar sets -1
// when absent; batched kernel may store 0 or -1; treat both as "no
// suboptimal." scalar.score == 0 skips the check.
bool kswr_score2_eq(const kswr_t &scalar, const kswr_t &batched);

} // namespace bwa_tests

#endif
