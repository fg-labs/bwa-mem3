/*************************************************************************************
                           The MIT License

   BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),
   Copyright (C) 2019  Intel Corporation, Heng Li.

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.

Authors: Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>;
*****************************************************************************************/

#ifndef SCALAR_BANDEDSWA_HPP
#define SCALAR_BANDEDSWA_HPP

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include "macro.h"

/* SIMD compatibility layer for ARM/x86 */
#if defined(__ARM_NEON) || defined(__aarch64__) || defined(APPLE_SILICON)
    /* ARM/Apple Silicon - use sse2neon for SSE translation */
    #include "simd_compat.h"
#elif (__AVX512BW__ || __AVX2__)
    #include <immintrin.h>
#else
    #include <smmintrin.h>  // for SSE4.1
    #ifndef __mmask8
    #define __mmask8 uint8_t
    #endif
    #ifndef __mmask16
    #define __mmask16 uint16_t
    #endif
#endif

#define MAX_SEQ_LEN_REF 256
#define MAX_SEQ_LEN_QER 128
#define MAX_SEQ_LEN_EXT 256
#define MAX_NUM_PAIRS 10000000
#define MAX_NUM_PAIRS_ALLOC 20000

/* Whether any vector banded-SW variant (getScores8/getScores16 +
 * smithWatermanBatchWrapper8/16) is declared and defined. ARM gets the
 * SSE2/NEON path via sse2neon. x86 gets the SSE2 path only with SSSE3
 * (the SSE2/NEON kernels use _mm_shuffle_epi8/PSHUFB), the AVX2 path with
 * AVX2, or the AVX512 path with AVX512BW. The Makefile passes -mssse3 on
 * every x86 arch target, so this is true in every CI build today — but
 * any caller-side dispatch must consult this macro rather than hand-rolled
 * subset checks (e.g. `!__SSE2__`) so a hypothetical SSE2-only-no-SSSE3
 * build still picks scalarBandedSWAWrapper instead of link-failing. */
#if (defined(__ARM_NEON) || defined(__aarch64__) || defined(APPLE_SILICON)) || \
    (__AVX512BW__) || (__AVX2__) || ((__SSE2__) && (__SSSE3__))
#define HAVE_BSW_VECTOR_8_16 1
#else
#define HAVE_BSW_VECTOR_8_16 0
#endif

// used in BSW and SAM-SW
#define DEFAULT_AMBIG -1


// SIMD_WIDTH in bits
// ARM64/NEON (128-bit vectors)
#if defined(__ARM_NEON) || defined(__aarch64__) || defined(APPLE_SILICON)
#ifndef SIMD_WIDTH8
#define SIMD_WIDTH8 16    // 128-bit / 8-bit = 16 elements
#endif
#ifndef SIMD_WIDTH16
#define SIMD_WIDTH16 8    // 128-bit / 16-bit = 8 elements
#endif

// AVX2
#elif ((!__AVX512BW__) & (__AVX2__))
#define SIMD_WIDTH8 32
#define SIMD_WIDTH16 16

// AVX512
#elif __AVX512BW__
#define SIMD_WIDTH8 64
#define SIMD_WIDTH16 32

// SSE2
#elif ((!__AVX512BW__) & (!__AVX2__) & (__SSE2__))
#define SIMD_WIDTH8 16
#define SIMD_WIDTH16 8

// Scalar
#else
#define SIMD_WIDTH8 1
#define SIMD_WIDTH16 1
#endif

#define MAX_LINE_LEN 256
#define MAX_SEQ_LEN8 128
#define MAX_SEQ_LEN16 32768
#define MATRIX_MIN_CUTOFF 0
#define LOW_INIT_VALUE -128
#define SORT_BLOCK_SIZE 16384
#define min_(x, y) ((x)>(y)?(y):(x))
#define max_(x, y) ((x)>(y)?(x):(y))

typedef struct dnaSeqPair
{
    int32_t idr, idq, id;
    int32_t len1, len2;
    int32_t h0;
    int seqid, regid;
    int32_t score, tle, gtle, qle;
    int32_t gscore, max_off;
    // PR 26c.1: per-pair SW band upper bound derived from the ungapped
    // score. 0 means "use default opt->w"; any positive value is an
    // upper bound on useful band offset for this pair. Default-initialized
    // here so stack-local SeqPair instances that don't reach every
    // construction-site assignment (e.g. RIGHT-only paths that skip the
    // LEFT tight_band assignment, or seeding paths that don't run
    // ungapped_analyze) start at the safe sentinel rather than indeterminate.
    int32_t tight_band = 0;
    // Q3 instrumentation: would-be ungapped extension score (full diagonal
    // walk, mirrors the HIT-path walk semantics). Computed at LEFT queue
    // time for non-HIT pairs; carried through SW retry-collect; read at
    // commit. -1 means undefined (e.g., zero-length).
    int32_t ugp_walk_score = -1;
    // Set by the construction-time RIGHT ungapped attempt (when
    // a->score != -1 made fp_h0 known) so the post-left-SW pass doesn't
    // double-count UGP_R_ATTEMPT / UGP_R_TIGHT / UGP_R_HIT.
    uint8_t ugp_r_attempted = 0;
}SeqPair;


typedef struct dnaOutScore
{
    int32_t score, tle, gtle, qle;
    int32_t gscore, max_off;
} OutScore;

typedef struct {
    int32_t h, e;
} eh_t;


class BandedPairWiseSW {
    
public:
    uint64_t SW_cells;

    BandedPairWiseSW(const int o_del, const int e_del, const int o_ins,
                     const int e_ins, const int zdrop,
                     const int end_bonus, const int8_t *mat_,
                     const int8_t w_match, const int8_t w_mismatch, int numThreads);
    ~BandedPairWiseSW();

    // Owns dp_slab_ via raw pointer; copying or moving would alias the
    // allocation and double-free on destruction. Disable both.
    BandedPairWiseSW(const BandedPairWiseSW&)            = delete;
    BandedPairWiseSW& operator=(const BandedPairWiseSW&) = delete;
    BandedPairWiseSW(BandedPairWiseSW&&)                 = delete;
    BandedPairWiseSW& operator=(BandedPairWiseSW&&)      = delete;

    // Scalar code section
    int scalarBandedSWA(int qlen, const uint8_t *query, int tlen,
                        const uint8_t *target, int32_t w,
                        int h0, int *_qle, int *_tle,
                        int *_gtle, int *_gscore,
                        int *_max_off);

    void scalarBandedSWAWrapper(SeqPair *seqPairArray,
                                uint8_t *seqBufRef,
                                uint8_t *seqBufQer,
                                int numPairs,
                                int nthreads,
                                int32_t w);

#if (defined(__ARM_NEON) || defined(__aarch64__) || defined(APPLE_SILICON)) || ((!__AVX512BW__) && (!__AVX2__) && (__SSE2__) && (__SSSE3__))
    // On ARM: use SSE2 path via sse2neon translation
    // On x86 SSE2: native SSE2 implementation (requires SSSE3 for
    // _mm_shuffle_epi8 / PSHUFB used by SBT_PREPASS8_LUT). The declaration
    // guard here MUST stay in lockstep with the matching definition guard
    // in bandedSWA.cpp — the `#if ((!__AVX512BW__) && (!__AVX2__) &&
    // (__SSE2__) && (__SSSE3__))` block that opens the SSE2/NEON section
    // (search for "SSE2 code"); otherwise the SSE2-only build will link-fail.
    // 8 bit vector code section
    void getScores8(SeqPair *pairArray,
                    uint8_t *seqBufRef,
                    uint8_t *seqBufQer,
                    int32_t numPairs,
                    uint16_t numThreads,
                    int32_t w);

    void smithWatermanBatchWrapper8(SeqPair *pairArray,
                                   uint8_t *seqBufRef,
                                   uint8_t *seqBufQer,
                                   int32_t numPairs,
                                   uint16_t numThreads,
                                   int32_t w);

    void smithWaterman128_8(uint8_t seq1SoA[],
                            uint8_t seq2SoA[],
                            uint8_t nrow,
                            uint8_t ncol,
                            SeqPair *p,
                            uint8_t h0[],
                            uint16_t tid,
                            int32_t numPairs,
                            int zdrop,
                            int32_t w,
                            uint8_t qlen[],
                            uint8_t myband[]);
    // 16 bit vector code section
    void getScores16(SeqPair *pairArray,
                     uint8_t *seqBufRef,
                     uint8_t *seqBufQer,
                     int32_t numPairs,
                     uint16_t numThreads,
                     int32_t w);

    void smithWatermanBatchWrapper16(SeqPair *pairArray,
                                     uint8_t *seqBufRef,
                                     uint8_t *seqBufQer,
                                     int32_t numPairs,
                                     uint16_t numThreads,
                                     int32_t w);
    
    void smithWaterman128_16(uint16_t seq1SoA[],
                             uint16_t seq2SoA[],
                             uint16_t nrow,
                             uint16_t ncol,
                             SeqPair *p,
                             uint16_t h0[],
                             uint16_t tid,
                             int32_t numPairs,
                             int zdrop,
                             int32_t w,
                             uint16_t qlen[],
                             uint16_t myband[]);
    
#endif  // ARM/SSE2

#if !defined(__ARM_NEON) && !defined(__aarch64__) && !defined(APPLE_SILICON)
#if ((!__AVX512BW__) & (__AVX2__))
    // AVX256 is not updated for banding and separate ins/del in the inner loop.
    // 8 bit vector code section    
    void getScores8(SeqPair *pairArray,
                    uint8_t *seqBufRef,
                    uint8_t *seqBufQer,
                    int32_t numPairs,
                    uint16_t numThreads,
                    int32_t w);

    void smithWatermanBatchWrapper8(SeqPair *pairArray,
                                   uint8_t *seqBufRef,
                                   uint8_t *seqBufQer,
                                   int32_t numPairs,
                                   uint16_t numThreads,
                                   int32_t w);

    void smithWaterman256_8(uint8_t seq1SoA[],
                            uint8_t seq2SoA[],
                            uint8_t nrow,
                            uint8_t ncol,
                            SeqPair *p,
                            uint8_t h0[],
                            uint16_t tid,
                            int32_t numPairs,
                            int zdrop,
                            int32_t w,
                            uint8_t qlen[],
                            uint8_t myband[]);
    // 16 bit vector code section
    void getScores16(SeqPair *pairArray,
                     uint8_t *seqBufRef,
                     uint8_t *seqBufQer,
                     int32_t numPairs,
                     uint16_t numThreads,
                     int32_t w);

    void smithWatermanBatchWrapper16(SeqPair *pairArray,
                                     uint8_t *seqBufRef,
                                     uint8_t *seqBufQer,
                                     int32_t numPairs,
                                     uint16_t numThreads,
                                     int32_t w);
    
    void smithWaterman256_16(uint16_t seq1SoA[],
                             uint16_t seq2SoA[],
                             uint16_t nrow,
                             uint16_t ncol,
                             SeqPair *p,
                             uint16_t h0[],
                             uint16_t tid,
                             int32_t numPairs,
                             int zdrop,
                             int32_t w,
                             uint16_t qlen[],
                             uint16_t myband[]);
    
#endif  //avx2
#endif  // !ARM guard for AVX2

#if !defined(__ARM_NEON) && !defined(__aarch64__) && !defined(APPLE_SILICON)
#if __AVX512BW__
    // 8 bit vector code section    
    void getScores8(SeqPair *pairArray,
                    uint8_t *seqBufRef,
                    uint8_t *seqBufQer,
                    int32_t numPairs,
                    uint16_t numThreads,
                    int32_t w);

    void smithWatermanBatchWrapper8(SeqPair *pairArray,
                                   uint8_t *seqBufRef,
                                   uint8_t *seqBufQer,
                                   int32_t numPairs,
                                   uint16_t numThreads,
                                   int32_t w);

    void smithWaterman512_8(uint8_t seq1SoA[],
                            uint8_t seq2SoA[],
                            uint8_t nrow,
                            uint8_t ncol,
                            SeqPair *p,
                            uint8_t h0[],
                            uint16_t tid,
                            int32_t numPairs,
                            int zdrop,
                            int32_t w,
                            uint8_t qlen[],
                            uint8_t myband[]);

    // 16 bit vector code section
    void getScores16(SeqPair *pairArray,
                     uint8_t *seqBufRef,
                     uint8_t *seqBufQer,
                     int32_t numPairs,
                     uint16_t numThreads,
                     int32_t w);

    void smithWatermanBatchWrapper16(SeqPair *pairArray,
                                     uint8_t *seqBufRef,
                                     uint8_t *seqBufQer,
                                     int32_t numPairs,
                                     uint16_t numThreads,
                                     int32_t w);
    
    void smithWaterman512_16(uint16_t seq1SoA[],
                             uint16_t seq2SoA[],
                             uint16_t nrow,
                             uint16_t ncol,
                             SeqPair *p,
                             uint16_t h0[],
                             uint16_t tid,
                             int32_t numPairs,
                             int zdrop,
                             int32_t w,
                             uint16_t qlen[],
                             uint16_t myband[]);
#endif  // __AVX512BW__
#endif  // !ARM guard for AVX512

    int64_t getTicks();
    
private:
    int m;
    int end_bonus, zdrop;
    int o_del, o_ins, e_del, e_ins;
    const int8_t *mat;

    int8_t w_match;
    int8_t w_mismatch;
    int8_t w_open;
    int8_t w_extend;
    int8_t w_ambig;
    int8_t *F8_;
    int8_t *H8_, *H8__;

    int16_t *F16_;
    int16_t *H16_, *H16__;

    // Single 64-byte-aligned slab backing F8_/H8_/H8__/F16_/H16_/H16__
    // and the per-thread SBT pre-pass scratch (sbt8_/sbt16_). The eight
    // member pointers above/below are views into this slab; the destructor
    // frees only `dp_slab_`. sbt8_ / sbt16_ moved off the stack — the
    // 16-bit variant (MAX_SEQ_LEN16 * SIMD_WIDTH16 * 2 bytes) is multi-MB
    // and was a stack-overflow risk on small-stack threads.
    int8_t  *sbt8_;
    int16_t *sbt16_;
    void *dp_slab_;

    int64_t sort1Ticks;
    int64_t setupTicks;
    int64_t swTicks;
    int64_t sort2Ticks;
};


#define DP  4
#define DP1 5
#define DP2 6
#define DP3 7

#endif
