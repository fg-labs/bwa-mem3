// test/framework/seqpair_batch.h
//
// RAII buffer owner that flattens a batch of TestPair into the flat byte
// layout kswv::getScores8 reads. Handles SIMD_WIDTH8 tail padding and the
// between-phase revseq / pair-compaction dance that mate-rescue production
// (bwamem_pair.cpp:660-676) performs between phase 0 and phase 1.

#ifndef BWA_TESTS_SEQPAIR_BATCH_H
#define BWA_TESTS_SEQPAIR_BATCH_H

#include <cstdint>
#include <vector>

#include "bandedSWA.h"  // SeqPair
#include "ksw.h"        // kswr_t
#include "seqpair.h"

namespace bwa_tests {

// Tail padding is sized for the WIDEST 8-bit lane count of any tier (64, the
// AVX-512BW SIMD_WIDTH8), not this translation unit's compile-time
// SIMD_WIDTH8. The kernel rounds the batch up to a whole SIMD group and writes
// those dummy lanes, so a test that reaches a kernel wider than the tier the
// framework was compiled at -- which any test going through make_kswv() can
// do, since that dispatches on the runtime tier and honors BWAMEM3_FORCE_TIER
// -- would otherwise write past the end of these buffers.
//
// It backs getScores16 batches too: the widest SIMD_WIDTH16 is 32, so the
// 8-bit bound covers both and there is no reason to carry two constants.
constexpr int MAX_SIMD_WIDTH8 = 64;

class BatchBuffers {
public:
    // Pack `pairs` into internal seqBufRef/seqBufQer with tail padding for the
    // widest supported SIMD group. After construction,
    // pairs()/ref_buf()/qer_buf()/aln() are the correctly-shaped arguments for
    // kswv::getScores8.
    //
    // xtra_flags is written into every pair's h0 field — pass the exact
    // value production uses at bwamem_pair.cpp:1004 so that the minsc
    // threshold in the low 16 bits filters consistently.
    BatchBuffers(const std::vector<TestPair> &pairs, int xtra_flags);

    // Accessors for passing into kswv::getScores8.
    SeqPair  *pairs()    { return pairs_.data(); }
    uint8_t  *ref_buf()  { return seqBufRef_.data(); }
    uint8_t  *qer_buf()  { return seqBufQer_.data(); }
    kswr_t   *aln()      { return aln_.data(); }
    int       n() const  { return n_; }

    // Phase-0 → phase-1 transition. For each pair whose phase-0 result
    // passes the KSW_XSTART / KSW_XSUBO gate, reverse the consumed ref/qry
    // prefixes in place, set h0 = KSW_XSTOP | score, trim len2 to qe+1,
    // and compact the survivors to the front of pairs_. Returns the number
    // of surviving pairs, i.e. the `nPairs` argument for the phase-1 call.
    int prepare_phase1();

private:
    std::vector<SeqPair> pairs_;
    std::vector<uint8_t> seqBufRef_;
    std::vector<uint8_t> seqBufQer_;
    std::vector<kswr_t>  aln_;
    int                  n_;
};

} // namespace bwa_tests

#endif
