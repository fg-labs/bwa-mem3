// test/framework/ksw_runner.cpp
#include "ksw_runner.h"

namespace bwa_tests {

int default_xtra_flags(int match_score, int min_seed_len) {
    // Mirror production's gate at src/bwamem_pair.cpp:209 — byte mode is
    // only enabled when minsc < 250. Keeps the test runner consistent with
    // the production flag selection for non-default match_score / seed_len.
    const int minsc = min_seed_len * match_score;
    const int xbyte = (minsc < 250) ? KSW_XBYTE : 0;
    return KSW_XSUBO | KSW_XSTART | xbyte | minsc;
}

kswr_t run_scalar_ksw(const TestPair &p,
                      const ScoringMatrix &mat,
                      int gap_open,
                      int gap_extend,
                      int xtra_flags) {
    if (xtra_flags == 0) {
        // Derive the KSW_XSUBO threshold (low 16 bits) from the matrix's
        // match score rather than the default=1, so tests that build a
        // non-unit-match ScoringMatrix keep scalar and batched flags
        // consistent.
        xtra_flags = default_xtra_flags(static_cast<int>(mat[0]));
    }
    // ksw_align2 declares non-const pointers; under KSW_XSTART it transiently
    // reverses the consumed query/target prefixes in place via revseq() and
    // restores them before return (see src/ksw.cpp:376-382). Net effect is
    // observationally const, but callers running concurrently against the
    // same buffers must not assume the data is undisturbed mid-call.
    return ksw_align2(static_cast<int>(p.qry.size()),
                      const_cast<uint8_t *>(p.qry.data()),
                      static_cast<int>(p.ref.size()),
                      const_cast<uint8_t *>(p.ref.data()),
                      5,
                      const_cast<int8_t *>(mat.data()),
                      gap_open, gap_extend, gap_open, gap_extend,
                      xtra_flags,
                      nullptr);
}

} // namespace bwa_tests
