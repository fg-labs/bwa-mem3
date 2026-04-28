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

// Mirrors the ifdef logic at src/kswv.cpp:175 / 1182 / 1692 that gates
// kswv::getScores8 definitions. On SSE4.1-only x86 no SIMD branch compiles,
// so libbwa.a omits the symbol. Tests that would call run_kswv_batch must
// #if on BWA_TESTS_HAVE_KSWV. Test build and libbwa.a MUST be compiled
// with the same arch flags so this macro evaluates consistently.
#if defined(__ARM_NEON) || defined(__aarch64__) || defined(APPLE_SILICON) || \
    defined(__AVX2__) || defined(__AVX512BW__)
#define BWA_TESTS_HAVE_KSWV 1
#else
#define BWA_TESTS_HAVE_KSWV 0
#endif

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
