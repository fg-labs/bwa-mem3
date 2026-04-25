// Regression test for issue 38 / upstream PR 289.
//
// When every real pair in a SIMD-aligned batch has len1 == 0, the kswv DP
// loop never executes (nrow == 0), yet the post-loop rowMax store still
// writes SIMD_WIDTH* bytes to `rowMax + (i - 1) * SIMD_WIDTH*` with i == 0.
// That stomps the allocation prefix, producing `free(): invalid pointer`
// on glibc and asan-reportable heap-buffer-overflow (negative offset).
//
// This test exercises the 8-bit and 16-bit entry points with batches of
// all-len1==0 pairs. Without the kernel-side `if (i > 0)` guard, each
// kswv destructor / allocator touch after the overwrite is likely to
// abort. Exits 0 only when all five kernels short-circuit cleanly.
//
// Arch coverage is implicit in the build: NEON builds compile kswv_neon_u8
// and kswv_neon_16; AVX2 builds compile kswv256_u8 (getScores16 is a scalar
// fallback on AVX2 and won't trigger the crash); AVX-512BW builds compile
// kswv512_u8 and kswv512_16.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <vector>

#include "ksw.h"
#include "kswv.h"

static const int8_t SCORE_MATCH    = 1;
static const int8_t SCORE_MISMATCH = 4;
static const int    GAP_OPEN       = 6;
static const int    GAP_EXT        = 1;
static const int    MIN_SEED_LEN   = 19;

// Match production mate-rescue xtra: forces the kernel through the same
// KSW_XSTART/KSW_XSUBO paths that the crash reporter hit.
static const int XTRA_FLAGS = KSW_XSUBO | KSW_XSTART | KSW_XBYTE
                              | (MIN_SEED_LEN * SCORE_MATCH);

static void run_8(int n_pairs, int phase) {
    const int32_t maxRefLen = 128;
    const int32_t maxQerLen = 128;

    std::vector<SeqPair> pairs(n_pairs + SIMD_WIDTH8);
    std::vector<uint8_t> refBuf(64, 0);
    std::vector<uint8_t> qerBuf(64, 0);
    std::vector<kswr_t>  aln(n_pairs + SIMD_WIDTH8, g_defr);

    // Every real pair has len1 == 0 so maxLen1 across the batch is 0 and
    // the kernel runs with nrow == 0. len2 is irrelevant for the guard
    // but set to zero to keep the SoA packing trivial.
    for (int i = 0; i < n_pairs; i++) {
        SeqPair sp = {};
        sp.idr    = 0;
        sp.idq    = 0;
        sp.id     = i;
        sp.len1   = 0;
        sp.len2   = 0;
        sp.h0     = XTRA_FLAGS;
        sp.seqid  = i;
        sp.regid  = i;
        pairs[i] = sp;
    }

    kswv *pwsw = new kswv(GAP_OPEN, GAP_EXT, GAP_OPEN, GAP_EXT,
                          SCORE_MATCH, -SCORE_MISMATCH, 1, maxRefLen, maxQerLen);
    pwsw->getScores8(pairs.data(), refBuf.data(), qerBuf.data(),
                     aln.data(), n_pairs, 1, phase);
    delete pwsw;
}

static void run_16(int n_pairs, int phase) {
    const int32_t maxRefLen = 128;
    const int32_t maxQerLen = 128;

    std::vector<SeqPair> pairs(n_pairs + SIMD_WIDTH16);
    std::vector<uint8_t> refBuf(64, 0);
    std::vector<uint8_t> qerBuf(64, 0);
    std::vector<kswr_t>  aln(n_pairs + SIMD_WIDTH16, g_defr);

    for (int i = 0; i < n_pairs; i++) {
        SeqPair sp = {};
        sp.idr    = 0;
        sp.idq    = 0;
        sp.id     = i;
        sp.len1   = 0;
        sp.len2   = 0;
        sp.h0     = XTRA_FLAGS;
        sp.seqid  = i;
        sp.regid  = i;
        pairs[i] = sp;
    }

    kswv *pwsw = new kswv(GAP_OPEN, GAP_EXT, GAP_OPEN, GAP_EXT,
                          SCORE_MATCH, -SCORE_MISMATCH, 1, maxRefLen, maxQerLen);
    pwsw->getScores16(pairs.data(), refBuf.data(), qerBuf.data(),
                      aln.data(), n_pairs, 1, phase);
    delete pwsw;
}

// Exercise both phases. Production batched mate-rescue calls getScores8/16
// with phase == 0 followed by phase == 1 on the same compacted batch; the
// guarded post-loop store fires on both, so cover both explicitly.
int main() {
    for (int phase = 0; phase <= 1; phase++) {
        fprintf(stderr, "[nrow-zero-test] 8-bit numPairs=1 phase=%d ...\n", phase);
        run_8(1, phase);
        fprintf(stderr, "[nrow-zero-test] 8-bit numPairs=SIMD_WIDTH8+1 (=%d) phase=%d ...\n",
                SIMD_WIDTH8 + 1, phase);
        run_8(SIMD_WIDTH8 + 1, phase);

        fprintf(stderr, "[nrow-zero-test] 16-bit numPairs=1 phase=%d ...\n", phase);
        run_16(1, phase);
        fprintf(stderr, "[nrow-zero-test] 16-bit numPairs=SIMD_WIDTH16+1 (=%d) phase=%d ...\n",
                SIMD_WIDTH16 + 1, phase);
        run_16(SIMD_WIDTH16 + 1, phase);
    }

    fprintf(stderr, "kswv_nrow_zero_test: OK\n");
    return 0;
}
