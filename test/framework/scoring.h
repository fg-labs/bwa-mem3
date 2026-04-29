// test/framework/scoring.h
//
// Scoring-matrix helpers for bwa-mem3 test fixtures.
//
// The bwa-mem3 production code fills a 5x5 substitution matrix (ACGTN) using
// `bwa_fill_scmat(a, b, ambig, mat[25])` (see test/integration/main_banded.cpp
// pre-migration, or the inlined callers throughout). Tests need the same
// layout so that `ksw_align2` and `kswv::getScores8` see identical scoring.
//
// This header exposes that construction as a value-returning helper (no
// caller-supplied buffer) so test code can stay terse.

#ifndef BWA_TESTS_SCORING_H
#define BWA_TESTS_SCORING_H

#include <array>
#include <cstdint>

namespace bwa_tests {

// 5x5 match/mismatch/ambig matrix indexed as mat[i*5 + j] for bases 0..4
// (A,C,G,T,N). Diagonal = +match, off-diagonal = -mismatch, N-rows/cols =
// -ambig. The sign convention matches bwa_fill_scmat and the internal
// w_ambig = -1 convention used by kswv.cpp.
using ScoringMatrix = std::array<int8_t, 25>;

// Build a matrix with the given match/mismatch/ambig weights. All three
// arguments are magnitudes — the returned matrix stores them with the
// appropriate signs.
ScoringMatrix build_scoring_matrix(int match, int mismatch, int ambig);

// bwa-mem3's default scoring (match=1, mismatch=4, ambig=1). Provided for
// tests that don't want to repeat the magic numbers.
ScoringMatrix default_scoring_matrix();

} // namespace bwa_tests

#endif
