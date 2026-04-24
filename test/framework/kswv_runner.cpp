// test/framework/kswv_runner.cpp
#include "kswv_runner.h"

#include <cstdio>
#include <cstdlib>

#include "kswv.h"
#include "seqpair_batch.h"

namespace bwa_tests {

std::vector<kswr_t> run_kswv_batch(const std::vector<TestPair> &pairs,
                                   const ScoringMatrix &mat,
                                   int gap_open,
                                   int gap_extend,
                                   int xtra_flags) {
    if (xtra_flags == 0) {
        // Derive the KSW_XSUBO threshold (low 16 bits) from the matrix's
        // match score so the batched runner agrees with run_scalar_ksw
        // when tests use non-unit match matrices.
        xtra_flags = default_xtra_flags(static_cast<int>(mat[0]));
    }

    BatchBuffers bb(pairs, xtra_flags);

    // Determine kswv sizing parameters.
    int32_t maxRefLen = 0, maxQerLen = 0;
    for (const auto &p : pairs) {
        if (static_cast<int32_t>(p.ref.size()) > maxRefLen) {
            maxRefLen = static_cast<int32_t>(p.ref.size());
        }
        if (static_cast<int32_t>(p.qry.size()) > maxQerLen) {
            maxQerLen = static_cast<int32_t>(p.qry.size());
        }
    }

    kswv *pwsw = new kswv(gap_open, gap_extend, gap_open, gap_extend,
                          mat[0],   // match weight (diagonal entry)
                          -mat[1],  // mismatch magnitude (off-diagonal sign)
                          1, maxRefLen, maxQerLen);

    // Phase 0: forward pass.
    pwsw->getScores8(bb.pairs(), bb.ref_buf(), bb.qer_buf(),
                     bb.aln(), bb.n(), 1, 0);

    if (const char *dbg = std::getenv("BWA_TESTS_DEBUG_PHASE0")) {
        (void)dbg;
        int nprint = (bb.n() < 5) ? bb.n() : 5;
        for (int i = 0; i < nprint; i++) {
            kswr_t r = bb.aln()[i];
            std::fprintf(stderr,
                         "[debug-phase0] pair %d: score=%d te=%d qe=%d tb=%d qb=%d\n",
                         i, r.score, r.te, r.qe, r.tb, r.qb);
        }
    }

    // Phase 0 → 1 transition: revseq + compact survivors.
    int pos = bb.prepare_phase1();

    // Snapshot phase-0 result to allow phase-1 diagnostics to compare.
    std::vector<kswr_t> pre_phase1(bb.aln(), bb.aln() + bb.n());

    // Phase 1: reverse pass fills tb/qb.
    pwsw->getScores8(bb.pairs(), bb.ref_buf(), bb.qer_buf(),
                     bb.aln(), pos, 1, 1);

    if (const char *dbg = std::getenv("BWA_TESTS_DEBUG_PHASE1")) {
        (void)dbg;
        int dumped = 0;
        for (int i = 0; i < bb.n() && dumped < 10; i++) {
            kswr_t r = bb.aln()[i];
            if (r.tb == -1 && pre_phase1[i].score > 0) {
                std::fprintf(stderr,
                             "[debug-phase1] pair %d: phase0{score=%d te=%d qe=%d} "
                             "post_phase1{score=%d te=%d qe=%d tb=%d qb=%d}\n",
                             i, pre_phase1[i].score, pre_phase1[i].te, pre_phase1[i].qe,
                             r.score, r.te, r.qe, r.tb, r.qb);
                dumped++;
            }
        }
    }

    delete pwsw;

    std::vector<kswr_t> out(bb.aln(), bb.aln() + bb.n());
    return out;
}

} // namespace bwa_tests
