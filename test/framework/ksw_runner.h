// test/framework/ksw_runner.h
//
// Scalar ksw_align2 wrapper for tests that need a reference result to
// compare SIMD output against.

#ifndef BWA_TESTS_KSW_RUNNER_H
#define BWA_TESTS_KSW_RUNNER_H

#include "ksw.h"
#include "scoring.h"
#include "seqpair.h"

namespace bwa_tests {

// Default minimum-seed-length used by the kswv-family tests. Matches
// bwa-mem3's default mem_opt_t::min_seed_len.
constexpr int DEFAULT_MIN_SEED_LEN = 19;

// Default gap penalties for bwa-mem3: open=6, extend=1 (both strands).
constexpr int DEFAULT_GAP_OPEN   = 6;
constexpr int DEFAULT_GAP_EXTEND = 1;

// Default xtra flag mask for mate-rescue-style scoring: KSW_XSUBO |
// KSW_XSTART, optionally OR'd with KSW_XBYTE when minsc =
// min_seed_len * match_score < 250 (matching the production gate at
// src/bwamem_pair.cpp:209), and with the low 16 bits carrying minsc.
int default_xtra_flags(int match_score = 1,
                       int min_seed_len = DEFAULT_MIN_SEED_LEN);

// Run scalar ksw_align2 against a single TestPair. xtra_flags = 0 is a
// sentinel meaning "derive from mat[0]" via default_xtra_flags(match_score);
// pass an explicit non-zero value to override.
kswr_t run_scalar_ksw(const TestPair &p,
                      const ScoringMatrix &mat,
                      int gap_open   = DEFAULT_GAP_OPEN,
                      int gap_extend = DEFAULT_GAP_EXTEND,
                      int xtra_flags = 0);

} // namespace bwa_tests

#endif
