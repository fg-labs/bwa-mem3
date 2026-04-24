// test/framework/kswv_runner.h
//
// Batched SIMD kswv runner — two-pass (phase 0 + phase 1) execution
// matching the production mate-rescue path at bwamem_pair.cpp:660-676.

#ifndef BWA_TESTS_KSWV_RUNNER_H
#define BWA_TESTS_KSWV_RUNNER_H

#include <vector>

#include "ksw.h"
#include "ksw_runner.h"  // DEFAULT_* constants
#include "scoring.h"
#include "seqpair.h"

namespace bwa_tests {

// Run kswv::getScores8 over `pairs` through the full two-pass sequence
// production uses for mate rescue. Returns one kswr_t per pair (same order
// as input). xtra_flags = 0 is a sentinel meaning "derive from mat[0]" via
// default_xtra_flags(match_score); pass an explicit non-zero value to
// override.
std::vector<kswr_t> run_kswv_batch(const std::vector<TestPair> &pairs,
                                   const ScoringMatrix &mat,
                                   int gap_open   = DEFAULT_GAP_OPEN,
                                   int gap_extend = DEFAULT_GAP_EXTEND,
                                   int xtra_flags = 0);

} // namespace bwa_tests

#endif
