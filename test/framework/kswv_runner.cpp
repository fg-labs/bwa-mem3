// test/framework/kswv_runner.cpp
#include "kswv_runner.h"

#if BWA_TESTS_HAVE_KSWV

#include <cstdio>
#include <cstdlib>
#include <memory>

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

    // Short-circuit for an empty pair set so we never construct kswv with
    // maxRefLen=maxQerLen=0 (the SIMD inner kernel's allocation/loop
    // bounds are not validated for n==0). Current callers don't hit this,
    // but defensive: future tests that build pair sets dynamically might.
    if (pairs.empty()) {
        return {};
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

    // C++11 throughout this tree, so spell out the unique_ptr ctor
    // rather than std::make_unique. RAII guarantees pwsw is destroyed
    // even if a later step (prepare_phase1, vector copy) throws.
    std::unique_ptr<kswv> pwsw(new kswv(
        gap_open, gap_extend, gap_open, gap_extend,
        mat[0],   // match weight (diagonal entry, +match)
        mat[1],   // mismatch weight (off-diagonal, -mismatch)
        1, maxRefLen, maxQerLen));

    // Phase 0: forward pass.
    pwsw->getScores8(bb.pairs(), bb.ref_buf(), bb.qer_buf(),
                     bb.aln(), bb.n(), 1, 0);

    if (std::getenv("BWA_TESTS_DEBUG_PHASE0")) {
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

    // Snapshot phase-0 result so phase-1 diagnostics can compare. Only
    // populated when BWA_TESTS_DEBUG_PHASE1 is set — copying ~10k kswr_t
    // per call adds up in the bulk-random unit test.
    const char *dbg1 = std::getenv("BWA_TESTS_DEBUG_PHASE1");
    std::vector<kswr_t> pre_phase1;
    if (dbg1) {
        pre_phase1.assign(bb.aln(), bb.aln() + bb.n());
    }

    // Phase 1: reverse pass fills tb/qb.
    pwsw->getScores8(bb.pairs(), bb.ref_buf(), bb.qer_buf(),
                     bb.aln(), pos, 1, 1);

    if (dbg1) {
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

    std::vector<kswr_t> out(bb.aln(), bb.aln() + bb.n());
    return out;
}

} // namespace bwa_tests

#endif // BWA_TESTS_HAVE_KSWV
