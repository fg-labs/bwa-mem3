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

Authors: Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>
*****************************************************************************************/

#include <stddef.h>
#include <string.h>
#include <unistd.h>
#include "kswv.h"
#include "limits.h"

/* ARM/NEON support */
#if defined(__ARM_NEON) || defined(__aarch64__) || defined(APPLE_SILICON)
#include "neon_utils.h"
#include "simd_compat.h"
#endif


// ------------------------------------------------------------------------------------
// MACROs for vector code
#if MAINY
uint64_t prof[10][112];
#else
extern uint64_t prof[10][112];
#endif

#define AMBIG_ 4  // ambiguous base
// for 16 bit
#define DUMMY1_ 4
#define DUMMY2_ 5
#define DUMMY3 26
#define AMBR16 15
#define AMBQ16 16
// for 8-bit
#define DUMMY8 8
#define DUMMY5 5
#define AMBRQ 0xFF
#define AMBR 4
#define AMBQ 8


// -----------------------------------------------------------------------------------
#if __AVX512BW__

#define MAIN_SAM_CODE8_OPT(s1, s2, h00, h11, e11, f11, f21, max512, sft512) \
    {                                                                   \
        __m512i sbt11, xor11, or11;                                     \
        xor11 = _mm512_xor_si512(s1, s2);                               \
        sbt11 = _mm512_shuffle_epi8(permSft512, xor11);                 \
        __mmask64 cmpq = _mm512_cmpeq_epu8_mask(s2, five512);           \
        sbt11 = _mm512_mask_blend_epi8(cmpq, sbt11, sft512);            \
        or11 =  _mm512_or_si512(s1, s2);                                \
        __mmask64 cmp = _mm512_movepi8_mask(or11);                      \
        __m512i m11 = _mm512_adds_epu8(h00, sbt11);                     \
        m11 = _mm512_mask_blend_epi8(cmp, m11, zero512);                \
        m11 = _mm512_subs_epu8(m11, sft512);                            \
        h11 = _mm512_max_epu8(m11, e11);                                \
        h11 = _mm512_max_epu8(h11, f11);                                \
        __mmask64 cmp0 = _mm512_cmpgt_epu8_mask(h11, imax512);          \
        imax512 = _mm512_max_epu8(imax512, h11);                        \
        iqe512 = _mm512_mask_blend_epi8(cmp0, iqe512, l512);            \
        __m512i gapE512 = _mm512_subs_epu8(h11, oe_ins512);             \
        e11 = _mm512_subs_epu8(e11, e_ins512);                          \
        e11 = _mm512_max_epu8(gapE512, e11);                            \
        __m512i gapD512 = _mm512_subs_epu8(h11, oe_del512);             \
        f21 = _mm512_subs_epu8(f11, e_del512);                          \
        f21 = _mm512_max_epu8(gapD512, f21);                            \
    }

#define MAIN_SAM_CODE16_OPT(s1, s2, h00, h11, e11, f11, f21, max512)    \
    {                                                                   \
        __m512i sbt11, xor11, or11;                                     \
        xor11 = _mm512_xor_si512(s1, s2);                               \
        sbt11 = _mm512_permutexvar_epi16(xor11, perm512);               \
        __m512i m11 = _mm512_add_epi16(h00, sbt11);                     \
        or11 =  _mm512_or_si512(s1, s2);                                \
        __mmask64 cmp = _mm512_movepi8_mask(or11);                      \
        m11 = _mm512_mask_blend_epi8(cmp, m11, zero512);                \
        h11 = _mm512_max_epi16(m11, e11);                               \
        h11 = _mm512_max_epi16(h11, f11);                               \
        h11 = _mm512_max_epi16(h11, zero512);                           \
        __mmask32 cmp0 = _mm512_cmpgt_epi16_mask(h11, imax512);         \
        imax512 = _mm512_max_epi16(imax512, h11);                       \
        iqe512 = _mm512_mask_blend_epi16(cmp0, iqe512, l512);           \
        __m512i gapE512 = _mm512_sub_epi16(h11, oe_ins512);             \
        e11 = _mm512_sub_epi16(e11, e_ins512);                          \
        e11 = _mm512_max_epi16(gapE512, e11);                           \
        __m512i gapD512 = _mm512_sub_epi16(h11, oe_del512);             \
        f21 = _mm512_sub_epi16(f11, e_del512);                          \
        f21 = _mm512_max_epi16(gapD512, f21);                           \
    }

#endif

// constructor
kswv::kswv(const int o_del, const int e_del, const int o_ins,
           const int e_ins, const int8_t w_match, const int8_t w_mismatch,
           int numThreads, int32_t maxRefLen = MAX_SEQ_LEN_REF_SAM,
           int32_t maxQerLen = MAX_SEQ_LEN_QER_SAM)
{

    this->m = 5;
    this->o_del = o_del;
    this->o_ins = o_ins;
    this->e_del = e_del;
    this->e_ins = e_ins;
    
    this->w_match    = w_match;
    this->w_mismatch = w_mismatch;
    this->w_open     = o_del;  // redundant, used in vector code.
    this->w_extend   = e_del;  // redundant, used in vector code.
    this->w_ambig    = DEFAULT_AMBIG;
    this->g_qmax = max_(w_match, w_mismatch);
    this->g_qmax = max_(this->g_qmax, w_ambig);

    this->maxRefLen = maxRefLen + 16;
    this->maxQerLen = maxQerLen + 16;
    
    this->swTicks = 0;
    setupTicks = 0;
    sort1Ticks = 0;
    swTicks = 0;
    sort2Ticks = 0;

    F16     = (int16_t *)_mm_malloc(this->maxQerLen * SIMD_WIDTH16 * numThreads * sizeof(int16_t), 64);
    H16_0   = (int16_t *)_mm_malloc(this->maxQerLen * SIMD_WIDTH16 * numThreads * sizeof(int16_t), 64);
    H16_1   = (int16_t *)_mm_malloc(this->maxQerLen * SIMD_WIDTH16 * numThreads * sizeof(int16_t), 64);
    H16_max = (int16_t *)_mm_malloc(this->maxQerLen * SIMD_WIDTH16 * numThreads * sizeof(int16_t), 64);
    rowMax16 = (int16_t *)_mm_malloc(this->maxRefLen * SIMD_WIDTH16 * numThreads * sizeof(int16_t), 64);

    F8 = (uint8_t*) F16;
    H8_0 = (uint8_t*) H16_0;
    H8_1 = (uint8_t*) H16_1;
    H8_max = (uint8_t*) H16_max;
    rowMax8 = (uint8_t*) rowMax16;
}

// destructor 
kswv::~kswv() {
    _mm_free(F16); _mm_free(H16_0); _mm_free(H16_max); _mm_free(H16_1);
    _mm_free(rowMax16);
}


/*******************************************************************************
 * ARM/NEON Implementation
 * Native NEON for 128-bit vectors (16 x 8-bit or 8 x 16-bit elements)
 * This provides optimized SIMD on Apple Silicon where sse2neon can't help
 * (AVX-512 has no sse2neon translation)
 ******************************************************************************/
#if defined(__ARM_NEON) || defined(__aarch64__) || defined(APPLE_SILICON)

void kswv::getScores8(SeqPair *pairArray,
                      uint8_t *seqBufRef,
                      uint8_t *seqBufQer,
                      kswr_t* aln,
                      int32_t numPairs,
                      uint16_t numThreads,
                      int phase)
{
    kswvBatchWrapper8(pairArray, seqBufRef, seqBufQer, aln,
                      numPairs, numThreads, phase);
}

void kswv::kswvBatchWrapper8(SeqPair *pairArray,
                             uint8_t *seqBufRef,
                             uint8_t *seqBufQer,
                             kswr_t* aln,
                             int32_t numPairs,
                             uint16_t numThreads,
                             int phase)
{
    uint8_t *seq1SoA = NULL;
    seq1SoA = (uint8_t *)_mm_malloc(this->maxRefLen * SIMD_WIDTH8 * numThreads * sizeof(uint8_t), 128);

    uint8_t *seq2SoA = NULL;
    seq2SoA = (uint8_t *)_mm_malloc(this->maxQerLen * SIMD_WIDTH8 * numThreads * sizeof(uint8_t), 128);

    assert(seq1SoA != NULL);
    assert(seq2SoA != NULL);

    int32_t ii;
    int32_t roundNumPairs = ((numPairs + SIMD_WIDTH8 - 1) / SIMD_WIDTH8) * SIMD_WIDTH8;

    for (ii = numPairs; ii < roundNumPairs; ii++)
    {
        pairArray[ii].regid = ii;
        pairArray[ii].id = ii;
        pairArray[ii].len1 = 0;
        pairArray[ii].len2 = 0;
    }

    {
        int32_t i;
        uint16_t tid = 0;   // forcing single thread
        uint8_t *mySeq1SoA = seq1SoA + tid * this->maxRefLen * SIMD_WIDTH8;
        uint8_t *mySeq2SoA = seq2SoA + tid * this->maxQerLen * SIMD_WIDTH8;
        uint8_t *seq1;
        uint8_t *seq2;

        int nstart = 0, nend = numPairs;

        for (i = nstart; i < nend; i += SIMD_WIDTH8)
        {
            int32_t j, k;
            int maxLen1 = 0;
            int maxLen2 = 0;

            for (j = 0; j < SIMD_WIDTH8; j++)
            {
                SeqPair sp = pairArray[i + j];
#if MAINY
                seq1 = seqBufRef + (int64_t)sp.id * this->maxRefLen;
#else
                seq1 = seqBufRef + sp.idr;
#endif
                for (k = 0; k < sp.len1; k++)
                {
                    mySeq1SoA[k * SIMD_WIDTH8 + j] = (seq1[k] == AMBIG_ ? AMBR : seq1[k]);
                }
                if (maxLen1 < sp.len1) maxLen1 = sp.len1;
            }
            for (j = 0; j < SIMD_WIDTH8; j++)
            {
                SeqPair sp = pairArray[i + j];
                for (k = sp.len1; k <= maxLen1; k++)
                {
                    mySeq1SoA[k * SIMD_WIDTH8 + j] = 0xFF;
                }
            }

            for (j = 0; j < SIMD_WIDTH8; j++)
            {
                SeqPair sp = pairArray[i + j];
#if MAINY
                seq2 = seqBufQer + (int64_t)sp.id * this->maxQerLen;
#else
                seq2 = seqBufQer + sp.idq;
#endif
                int quanta = (sp.len2 + 16 - 1) / 16;
                quanta *= 16;
                for (k = 0; k < sp.len2; k++)
                {
                    mySeq2SoA[k * SIMD_WIDTH8 + j] = (seq2[k] == AMBIG_ ? AMBQ : seq2[k]);
                }

                for (k = sp.len2; k < quanta; k++)
                {
                    mySeq2SoA[k * SIMD_WIDTH8 + j] = DUMMY5;
                }
                if (maxLen2 < quanta) maxLen2 = quanta;
            }

            for (j = 0; j < SIMD_WIDTH8; j++)
            {
                SeqPair sp = pairArray[i + j];
                int quanta = (sp.len2 + 16 - 1) / 16;
                quanta *= 16;
                for (k = quanta; k <= maxLen2; k++)
                {
                    mySeq2SoA[k * SIMD_WIDTH8 + j] = 0xFF;
                }
            }

            kswv_neon_u8(mySeq1SoA, mySeq2SoA,
                         maxLen1, maxLen2,
                         pairArray + i,
                         aln, i,
                         tid,
                         numPairs,
                         phase);
        }
    }

    _mm_free(seq1SoA);
    _mm_free(seq2SoA);

    return;
}

int kswv::kswv_neon_u8(uint8_t seq1SoA[],
                       uint8_t seq2SoA[],
                       int16_t nrow,
                       int16_t ncol,
                       SeqPair *p,
                       kswr_t *aln,
                       int po_ind,
                       uint16_t tid,
                       int32_t numPairs,
                       int phase)
{
    int m_b, n_b;
    uint8_t minsc[SIMD_WIDTH8] __attribute__((aligned(128))) = {0};
    uint8_t endsc[SIMD_WIDTH8] __attribute__((aligned(128))) = {0};
    uint64_t *b;

    uint8x16_t zero_vec = vdupq_n_u8(0);
    uint8x16_t one_vec = vdupq_n_u8(1);

    m_b = n_b = 0; b = 0;

    int8_t temp[SIMD_WIDTH8] __attribute((aligned(128))) = {0};

    uint8_t shift = 127, mdiff = 0;
    mdiff = max_(this->w_match, (int8_t)this->w_mismatch);
    mdiff = max_(mdiff, (int8_t)this->w_ambig);
    shift = min_(this->w_match, (int8_t)this->w_mismatch);
    shift = min_((int8_t)shift, this->w_ambig);

    shift = 256 - (uint8_t)shift;
    mdiff += shift;

    /* Build scoring table */
    temp[0] = this->w_match;
    temp[1] = temp[2] = temp[3] = this->w_mismatch;
    temp[4] = temp[5] = temp[6] = temp[7] = this->w_ambig;
    temp[8] = temp[9] = temp[10] = temp[11] = this->w_ambig;
    temp[12] = this->w_ambig;

    for (int i = 0; i < 16; i++)
        temp[i] += shift;

    uint8x16_t permSft = vld1q_u8((const uint8_t*)temp);
    uint8x16_t sft_vec = vdupq_n_u8(shift);

    uint16_t minsc_msk_a = 0x0000, endsc_msk_a = 0x0000;
    int val = 0;
    for (int i = 0; i < SIMD_WIDTH8; i++)
    {
        int xtra = p[i].h0;
        val = (xtra & KSW_XSUBO) ? xtra & 0xffff : 0x10000;
        if (val <= 255) {
            minsc[i] = val;
            minsc_msk_a |= (0x1 << i);
        }
        val = (xtra & KSW_XSTOP) ? xtra & 0xffff : 0x10000;
        if (val <= 255) {
            endsc[i] = val;
            endsc_msk_a |= (0x1 << i);
        }
    }

    uint8x16_t minsc_vec = vld1q_u8(minsc);
    uint8x16_t endsc_vec = vld1q_u8(endsc);

    uint8x16_t e_del_vec = vdupq_n_u8(this->e_del);
    uint8x16_t oe_del_vec = vdupq_n_u8(this->o_del + this->e_del);
    uint8x16_t e_ins_vec = vdupq_n_u8(this->e_ins);
    uint8x16_t oe_ins_vec = vdupq_n_u8(this->o_ins + this->e_ins);
    uint8x16_t five_vec = vdupq_n_u8(DUMMY5);
    uint8x16_t gmax_vec = zero_vec;
    // 8-bit SW processes 16 pairs per SIMD batch, but te must be tracked in
    // 16-bit precision (ref positions > 255 are possible). Use two
    // int16x8 vectors — _lo for pairs 0-7, _hi for pairs 8-15 — so every
    // pair has its own te slot. The earlier single int16x8 te_vec plus
    // hacky "duplicate on store" only tracked pairs 0-7 correctly; pairs
    // 8-15 got copies of the wrong pair's te.
    int16x8_t te_vec_lo = vdupq_n_s16(-1);
    int16x8_t te_vec_hi = vdupq_n_s16(-1);
    uint8x16_t cmax_vec = vdupq_n_u8(255);
    // "Frozen" mask: tracks which pairs have hit their KSW_XSTOP target
    // in a prior row. Once frozen, a pair's gmax/te/qe must not be updated
    // — otherwise later rows (which may contain untouched ref bytes when
    // the caller only partially reversed the buffer) can add spurious
    // matches, pushing phase-1's score above phase-0's and breaking the
    // write-back gate `aln[ind].score == score[l]` at line ~538.
    uint8x16_t frozen_vec = zero_vec;

    /* Per-lane mask of "this lane has a real KSW_XSTOP target set" (derived
     * from endsc_msk_a). endsc[i] defaults to 0 for lanes without a target,
     * which would otherwise make vcgeq_u8(gmax, endsc_vec) trivially true
     * and cause premature freezing. */
    uint8_t _has_endsc_bytes[SIMD_WIDTH8];
    for (int _i = 0; _i < SIMD_WIDTH8; _i++) {
        _has_endsc_bytes[_i] = (endsc_msk_a & (1 << _i)) ? 0xFF : 0x00;
    }
    uint8x16_t has_endsc_vec = vld1q_u8(_has_endsc_bytes);

    uint16_t exit0 = 0xFFFF;

    tid = 0;
    uint8_t *H0 = H8_0 + tid * SIMD_WIDTH8 * this->maxQerLen;
    uint8_t *H1 = H8_1 + tid * SIMD_WIDTH8 * this->maxQerLen;
    uint8_t *Hmax = H8_max + tid * SIMD_WIDTH8 * this->maxQerLen;
    uint8_t *F = F8 + tid * SIMD_WIDTH8 * this->maxQerLen;
    uint8_t *rowMax = rowMax8 + tid * SIMD_WIDTH8 * this->maxRefLen;

    /* Initialize arrays */
    for (int i = 0; i <= ncol; i++)
    {
        vst1q_u8(H0 + i * SIMD_WIDTH8, zero_vec);
        vst1q_u8(Hmax + i * SIMD_WIDTH8, zero_vec);
        vst1q_u8(F + i * SIMD_WIDTH8, zero_vec);
    }

    uint8x16_t max_vec = zero_vec, imax_vec, pimax_vec = zero_vec;
    uint16_t mask16 = 0x0000;
    uint16_t minsc_msk = 0x0000;

    uint8x16_t qe_vec = vdupq_n_u8(0);
    vst1q_u8(H0, zero_vec);
    vst1q_u8(H1, zero_vec);

    int i, limit = nrow;
    for (i = 0; i < nrow; i++)
    {
        uint8x16_t e11 = zero_vec;
        uint8x16_t h00, h11, h10, s1;
        int16x8_t i_vec = vdupq_n_s16(i);
        int j;

        s1 = vld1q_u8(seq1SoA + (i + 0) * SIMD_WIDTH8);
        h10 = zero_vec;
        imax_vec = zero_vec;
        uint8x16_t iqe_vec = vdupq_n_u8(0xFF);

        uint8x16_t l_vec = zero_vec;
        for (j = 0; j < ncol; j++)
        {
            uint8x16_t f11, s2, f21;
            h00 = vld1q_u8(H0 + j * SIMD_WIDTH8);
            s2 = vld1q_u8(seq2SoA + j * SIMD_WIDTH8);
            f11 = vld1q_u8(F + (j + 1) * SIMD_WIDTH8);

            /* Core Smith-Waterman computation using NEON */
            uint8x16_t xor_val = veorq_u8(s1, s2);
            uint8x16_t sbt = vqtbl1q_u8(permSft, xor_val);

            /* Check for ambiguous query base (DUMMY5) */
            uint8x16_t cmpq = vceqq_u8(s2, five_vec);
            sbt = vbslq_u8(cmpq, sft_vec, sbt);

            /* Check for boundary (high bit set) */
            uint8x16_t or_val = vorrq_u8(s1, s2);
            uint8x16_t high_bit = vshrq_n_u8(or_val, 7);
            uint8x16_t is_boundary = vceqq_u8(high_bit, one_vec);

            uint8x16_t m11 = vqaddq_u8(h00, sbt);
            m11 = vbslq_u8(is_boundary, zero_vec, m11);
            m11 = vqsubq_u8(m11, sft_vec);

            h11 = vmaxq_u8(m11, e11);
            h11 = vmaxq_u8(h11, f11);

            /* Update imax tracking */
            uint8x16_t cmp0 = vcgtq_u8(h11, imax_vec);
            imax_vec = vmaxq_u8(imax_vec, h11);
            iqe_vec = vbslq_u8(cmp0, l_vec, iqe_vec);

            /* Gap extension for E */
            uint8x16_t gapE = vqsubq_u8(h11, oe_ins_vec);
            e11 = vqsubq_u8(e11, e_ins_vec);
            e11 = vmaxq_u8(gapE, e11);

            /* Gap extension for F */
            uint8x16_t gapD = vqsubq_u8(h11, oe_del_vec);
            f21 = vqsubq_u8(f11, e_del_vec);
            f21 = vmaxq_u8(gapD, f21);

            vst1q_u8(H1 + (j + 1) * SIMD_WIDTH8, h11);
            vst1q_u8(F + (j + 1) * SIMD_WIDTH8, f21);
            l_vec = vaddq_u8(l_vec, one_vec);
        }

        /* Block I - row max tracking */
        if (i > 0)
        {
            uint8x16_t cmp_gt = vcgtq_u8(imax_vec, pimax_vec);
            uint16_t msk16 = neon_movemask_u8(cmp_gt);
            msk16 |= mask16;

            /* Apply masks */
            uint8x16_t msk_vec = vld1q_u8((uint8_t*)&msk16); // simplified
            pimax_vec = vbslq_u8(cmp_gt, zero_vec, pimax_vec);

            vst1q_u8(rowMax + (i - 1) * SIMD_WIDTH8, pimax_vec);
            mask16 = ~msk16;
        }
        pimax_vec = imax_vec;

        /* Check minsc threshold */
        uint8x16_t cmp_ge = vcgeq_u8(imax_vec, minsc_vec);
        minsc_msk = neon_movemask_u8(cmp_ge);
        minsc_msk &= minsc_msk_a;

        /* Block II: gmax, te */
        uint8x16_t cmp0 = vcgtq_u8(imax_vec, gmax_vec);
        uint16_t cmp0_msk = neon_movemask_u8(cmp0);
        cmp0_msk &= exit0;

        /* 16-bit-wide comparison mirrors cmp0 but with 0xFFFF/0x0000 per
         * 16-bit lane, suitable for updating the two halves of te. Must be
         * computed BEFORE the vmaxq_u8 update of gmax_vec (same as cmp0). */
        uint16x8_t cmp_lo_16 = vcgtq_u16(vmovl_u8(vget_low_u8(imax_vec)),
                                         vmovl_u8(vget_low_u8(gmax_vec)));
        uint16x8_t cmp_hi_16 = vcgtq_u16(vmovl_u8(vget_high_u8(imax_vec)),
                                         vmovl_u8(vget_high_u8(gmax_vec)));

        /* 16-bit frozen masks (0xFFFF per lane where frozen, else 0x0000).
         * Widen 0xFF/0x00 bytes: vtstq_u8-style check, but we can use
         * vmovl_u8 + vcgtq to turn non-zero into 0xFFFF. */
        uint16x8_t frozen_lo_16 = vcgtq_u16(
            vmovl_u8(vget_low_u8(frozen_vec)), vdupq_n_u16(0));
        uint16x8_t frozen_hi_16 = vcgtq_u16(
            vmovl_u8(vget_high_u8(frozen_vec)), vdupq_n_u16(0));

        /* Combine "imax > gmax" with "not frozen" for the final update masks */
        cmp_lo_16 = vbicq_u16(cmp_lo_16, frozen_lo_16);
        cmp_hi_16 = vbicq_u16(cmp_hi_16, frozen_hi_16);
        uint8x16_t cmp0_active = vbicq_u8(cmp0, frozen_vec);

        /* gmax: pick new max for active lanes, keep frozen lanes' stored
         * value. This is what prevents later rows from inflating a lane's
         * score beyond its phase-0 target once it's already "done". */
        uint8x16_t new_gmax = vmaxq_u8(gmax_vec, imax_vec);
        gmax_vec = vbslq_u8(frozen_vec, gmax_vec, new_gmax);

        /* Update te/qe only for active lanes (imax > gmax AND not frozen) */
        te_vec_lo = vbslq_s16(cmp_lo_16, i_vec, te_vec_lo);
        te_vec_hi = vbslq_s16(cmp_hi_16, i_vec, te_vec_hi);
        qe_vec = vbslq_u8(cmp0_active, iqe_vec, qe_vec);

        /* Check end score threshold */
        uint8x16_t cmp_end = vcgeq_u8(gmax_vec, endsc_vec);
        uint16_t cmp_end_msk = neon_movemask_u8(cmp_end);
        cmp_end_msk &= endsc_msk_a;

        /* Freeze any lane that just hit its KSW_XSTOP target. From the
         * next row onwards, frozen lanes will skip the gmax/te/qe updates
         * above. has_endsc_vec ensures we only freeze lanes that actually
         * set a target (endsc_vec defaults to 0 otherwise, which would
         * match cmp_end trivially). */
        uint8x16_t just_hit = vandq_u8(cmp_end, has_endsc_vec);
        frozen_vec = vorrq_u8(frozen_vec, just_hit);

        /* Check for overflow */
        uint8x16_t left_vec = vqaddq_u8(gmax_vec, sft_vec);
        uint8x16_t cmp2 = vcgeq_u8(left_vec, cmax_vec);
        uint16_t cmp2_msk = neon_movemask_u8(cmp2);

        exit0 = (~(cmp_end_msk | cmp2_msk)) & exit0;
        if (exit0 == 0)
        {
            limit = i++;
            break;
        }

        uint8_t *S = H1; H1 = H0; H0 = S;
    }

    /* Store final row max. Guard on i > 0: when every pair in the batch
     * has len1 == 0, nrow == 0, the DP loop never runs, and `i - 1`
     * underflows into the allocation prefix. Upstream PR 289 / issue 38. */
    if (i > 0) {
        vst1q_u8(rowMax + (i - 1) * SIMD_WIDTH8, pimax_vec);
    }

    /* Extract results */
    uint8_t score[SIMD_WIDTH8] __attribute((aligned(128)));
    int16_t te1[SIMD_WIDTH8] __attribute((aligned(128)));
    uint8_t qe[SIMD_WIDTH8] __attribute((aligned(128)));

    vst1q_u8(score, gmax_vec);
    vst1q_s16(te1, te_vec_lo);       // pairs 0-7
    vst1q_s16(te1 + 8, te_vec_hi);   // pairs 8-15
    vst1q_u8(qe, qe_vec);

    int live = 0;
    for (int l = 0; l < SIMD_WIDTH8 && (po_ind + l) < numPairs; l++) {
        int ind = po_ind + l;
        int16_t *te = te1;
#if !MAINY
        ind = p[l].regid;
        if (phase) {
            if (aln[ind].score == score[l]) {
                aln[ind].tb = aln[ind].te - te[l];
                aln[ind].qb = aln[ind].qe - qe[l];
            }
        } else {
            aln[ind].score = score[l] + shift < 255 ? score[l] : 255;
            aln[ind].te = te[l];
            aln[ind].qe = qe[l];
            if (aln[ind].score != 255) {
                qe[l] = 1;
                live++;
            }
            else qe[l] = 0;
        }
#else
        aln[ind].score = score[l] + shift < 255 ? score[l] : 255;
        aln[ind].te = te[l];
        aln[ind].qe = qe[l];
        if (aln[ind].score != 255) {
            qe[l] = 1;
            live++;
        }
        else qe[l] = 0;
#endif
    }

#if !MAINY
    if (phase) return 1;
#endif

    if (live == 0) return 1;

    /* Score2 and te2 computation.
     *
     * Per-lane scalar emulation of ksw_u8's b[] build + score2 scan
     * (src/ksw.cpp:201). scalar ksw_u8 collapses consecutive rows >=
     * minsc into a single b[] entry anchored at the max-score row; a
     * plateau straddling the primary region [te - val, te + val] is
     * excluded entirely when its anchor lies inside. A dense SIMD scan
     * treats every row independently and pulls out-of-region rows of a
     * boundary-straddling plateau in as false suboptimals, inflating `sub`.
     * Emulate scalar exactly, per lane, on the precomputed
     * rowMax. Cost is O(SIMD_WIDTH8 * limit) scalar ops per batch,
     * negligible next to the SIMD DP that built rowMax. */
    int qmax = this->g_qmax;
    int16_t low[SIMD_WIDTH8]  __attribute__((aligned(64)));
    int16_t high[SIMD_WIDTH8] __attribute__((aligned(64)));
    for (int j = 0; j < SIMD_WIDTH8; j++) {
        int val = (score[j] + qmax - 1) / qmax;
        low[j]  = te1[j] - val;
        high[j] = te1[j] + val;
    }

    /* rowMax entries produced: i == limit on normal completion, and
     * i == limit + 1 after early exit (limit = i++ stores old i). Using i
     * directly gives a single expression that covers both exit paths and
     * includes the terminating row stored at rowMax[limit]. */
    const int processed_rows = i;

    for (int l = 0; l < SIMD_WIDTH8 && (po_ind + l) < numPairs; l++) {
        int ind = p[l].regid;
        /* Match scalar ksw_u8: when KSW_XSUBO is absent, minsc = 0x10000
         * so the b[] list never starts and score2/te2 stay at -1.
         * minsc_msk_a tracks lanes with KSW_XSUBO set; lanes not in the
         * mask have minsc[l] zero-initialized and must skip the scan. */
        if (!qe[l] || !(minsc_msk_a & (1u << l))) {
            aln[ind].score2 = -1;
            aln[ind].te2    = -1;
            continue;
        }

        int len1_l  = (int)p[l].len1;
        int low_l   = (int)low[l];
        int high_l  = (int)high[l];
        int minsc_l = (int)minsc[l];
        int score2  = -1;
        int te2     = -1;
        /* b_pos = -2 sentinel: (-2 + 1) != any real row, forces an
         * "append" on the first qualifying row. */
        int b_score = -1;
        int b_pos   = -2;

        int nrows = processed_rows < len1_l ? processed_rows : len1_l;
        for (int i2 = 0; i2 < nrows; i2++) {
            int imax = (int)rowMax[i2 * SIMD_WIDTH8 + l];
            if (imax < minsc_l) continue;

            if (b_pos + 1 != i2) {
                /* APPEND: flush the outgoing b[] entry to score2 first. */
                if (b_pos >= 0 &&
                    (b_pos < low_l || b_pos > high_l) &&
                    b_score > score2) {
                    score2 = b_score;
                    te2    = b_pos;
                }
                b_score = imax;
                b_pos   = i2;
            } else if (b_score < imax) {
                /* UPDATE: strict greater extends the run's anchor. Equal
                 * imax leaves b_pos in place (matches scalar). */
                b_score = imax;
                b_pos   = i2;
            }
        }

        /* Flush trailing b[] entry. */
        if (b_pos >= 0 &&
            (b_pos < low_l || b_pos > high_l) &&
            b_score > score2) {
            score2 = b_score;
            te2    = b_pos;
        }

        aln[ind].score2 = score2;
        aln[ind].te2    = te2;
    }

    return 1;
}

/* 16-bit NEON implementation */
void kswv::getScores16(SeqPair *pairArray,
                       uint8_t *seqBufRef,
                       uint8_t *seqBufQer,
                       kswr_t* aln,
                       int32_t numPairs,
                       uint16_t numThreads,
                       int phase)
{
    kswvBatchWrapper16(pairArray, seqBufRef, seqBufQer, aln,
                       numPairs, numThreads, phase);
}

void kswv::kswvBatchWrapper16(SeqPair *pairArray,
                              uint8_t *seqBufRef,
                              uint8_t *seqBufQer,
                              kswr_t* aln,
                              int32_t numPairs,
                              uint16_t numThreads,
                              int phase)
{
    int16_t *seq1SoA = NULL;
    seq1SoA = (int16_t *)_mm_malloc(this->maxRefLen * SIMD_WIDTH16 * numThreads * sizeof(int16_t), 128);

    int16_t *seq2SoA = NULL;
    seq2SoA = (int16_t *)_mm_malloc(this->maxQerLen * SIMD_WIDTH16 * numThreads * sizeof(int16_t), 128);

    assert(seq1SoA != NULL);
    assert(seq2SoA != NULL);

    int32_t ii;
    int32_t roundNumPairs = ((numPairs + SIMD_WIDTH16 - 1) / SIMD_WIDTH16) * SIMD_WIDTH16;

    for (ii = numPairs; ii < roundNumPairs; ii++)
    {
        pairArray[ii].regid = ii;
        pairArray[ii].id = ii;
        pairArray[ii].len1 = 0;
        pairArray[ii].len2 = 0;
    }

    {
        int32_t i;
        uint16_t tid = 0;
        int16_t *mySeq1SoA = seq1SoA + tid * this->maxRefLen * SIMD_WIDTH16;
        int16_t *mySeq2SoA = seq2SoA + tid * this->maxQerLen * SIMD_WIDTH16;
        uint8_t *seq1;
        uint8_t *seq2;

        int nstart = 0, nend = numPairs;

        for (i = nstart; i < nend; i += SIMD_WIDTH16)
        {
            int32_t j, k;
            int maxLen1 = 0;
            int maxLen2 = 0;

            for (j = 0; j < SIMD_WIDTH16; j++)
            {
                SeqPair sp = pairArray[i + j];
#if MAINY
                seq1 = seqBufRef + (int64_t)sp.id * this->maxRefLen;
#else
                seq1 = seqBufRef + sp.idr;
#endif
                for (k = 0; k < sp.len1; k++)
                {
                    mySeq1SoA[k * SIMD_WIDTH16 + j] = (seq1[k] == AMBIG_ ? AMBR16 : seq1[k]);
                }
                if (maxLen1 < sp.len1) maxLen1 = sp.len1;
            }
            for (j = 0; j < SIMD_WIDTH16; j++)
            {
                SeqPair sp = pairArray[i + j];
                for (k = sp.len1; k <= maxLen1; k++)
                {
                    mySeq1SoA[k * SIMD_WIDTH16 + j] = 0xFFFF;
                }
            }

            for (j = 0; j < SIMD_WIDTH16; j++)
            {
                SeqPair sp = pairArray[i + j];
#if MAINY
                seq2 = seqBufQer + (int64_t)sp.id * this->maxQerLen;
#else
                seq2 = seqBufQer + sp.idq;
#endif
                int quanta = (sp.len2 + 8 - 1) / 8;
                quanta *= 8;
                for (k = 0; k < sp.len2; k++)
                {
                    mySeq2SoA[k * SIMD_WIDTH16 + j] = (seq2[k] == AMBIG_ ? AMBQ16 : seq2[k]);
                }

                for (k = sp.len2; k < quanta; k++)
                {
                    mySeq2SoA[k * SIMD_WIDTH16 + j] = DUMMY3;
                }
                if (maxLen2 < quanta) maxLen2 = quanta;
            }

            for (j = 0; j < SIMD_WIDTH16; j++)
            {
                SeqPair sp = pairArray[i + j];
                int quanta = (sp.len2 + 8 - 1) / 8;
                quanta *= 8;
                for (k = quanta; k <= maxLen2; k++)
                {
                    mySeq2SoA[k * SIMD_WIDTH16 + j] = 0xFFFF;
                }
            }

            kswv_neon_16(mySeq1SoA, mySeq2SoA,
                         maxLen1, maxLen2,
                         pairArray + i,
                         aln, i,
                         tid,
                         numPairs,
                         phase);
        }
    }

    _mm_free(seq1SoA);
    _mm_free(seq2SoA);
    return;
}

int kswv::kswv_neon_16(int16_t seq1SoA[],
                       int16_t seq2SoA[],
                       int16_t nrow,
                       int16_t ncol,
                       SeqPair *p,
                       kswr_t *aln,
                       int po_ind,
                       uint16_t tid,
                       int32_t numPairs,
                       int phase)
{
    /* NEON 16-bit mate-rescue kernel.
     *
     * Structural port of kswv_neon_u8: real rowMax writes each row,
     * per-lane freeze once a lane hits its KSW_XSTOP target (gated by
     * has_endsc_vec so lanes without a target aren't frozen prematurely),
     * per-row frozen_bits collapsed from frozen_vec for the batched early
     * exit, and a scalar b[]-emulation score2 scan over rowMax matching
     * scalar ksw_i16 exactly. The
     * pre-rewrite kernel never wrote rowMax, never computed score2
     * (hardcoded to -1), and global-broke the DP as soon as any lane's
     * gmax >= endsc_vec (trivially true for mate-rescue pairs whose
     * endsc_vec defaults to 0 since KSW_XSTOP is unset) — which made
     * `score`/`te`/`qe` almost garbage too. Replaced with proper SIMD
     * semantics matching the 8-bit NEON kernel.
     *
     * 16-bit specifics vs 8-bit:
     *   - 8 lanes (SIMD_WIDTH16), single te_vec (no _lo/_hi split)
     *   - int16 arithmetic; no shift offset, clamp intermediate H to 0
     *     via vmaxq_s16 (scalar ksw_i16 does the same)
     *   - No u8 overflow check; int16 has plenty of headroom for SW scores
     */
    int16_t minsc[SIMD_WIDTH16] __attribute((aligned(128))) = {0};
    int16_t endsc[SIMD_WIDTH16] __attribute((aligned(128))) = {0};

    int16x8_t zero_vec = vdupq_n_s16(0);
    int16x8_t one_vec  = vdupq_n_s16(1);

    /* Scoring table: 32-byte int8 table indexed by (s1 ^ s2) in [0..31].
     *
     * IMPORTANT: The table must be byte-wide (not int16-wide) AND must
     * cover all 32 possible xor values. The pre-rewrite kernel used an
     * 8-entry int16 table (16 bytes) and did `vqtbl1q_s8(perm, xor_val)`
     * + reinterpret to int16 — treating each int16 xor value's low+high
     * bytes as two byte indices, producing sbt[i] = (w_match << 8) |
     * perm[xor_low], inflating match scores by 256. That bug was latent
     * while smoke-1M never exercised the 16-bit path (l_ms*opt->a < 250
     * routes to 8-bit), but `-A 2` makes l_ms*2 = 300 >= 250 and
     * surfaces it.
     *
     * The 16-bit kernel's SoA uses wider ambig / tail-padding constants
     * than 8-bit (AMBR16=15, AMBQ16=16, DUMMY3=26 vs AMBR=4, AMBQ=8,
     * DUMMY5=5), so xor values can reach 31. A 32-byte table built via
     * vqtbl2_s8 handles all of them:
     *   [0]     -> w_match       (base == base; also AMBR^AMBR, AMBQ^AMBQ,
     *                             but those pairings don't occur)
     *   [1..3]  -> w_mismatch    (base != base)
     *   [21,24..27] -> 0         (anything ^ DUMMY3 = query tail padding,
     *                             must contribute 0 like scalar would)
     *   other   -> w_ambig       (any xor involving AMBR16 or AMBQ16)
     * Boundary (high bit of s1|s2) is handled separately by zeroing m11. */
    int8_t temp8[32] __attribute((aligned(16)));
    for (int i = 0; i < 32; i++) temp8[i] = this->w_ambig;
    temp8[0] = this->w_match;
    temp8[1] = temp8[2] = temp8[3] = this->w_mismatch;
    /* Query-tail DUMMY3 (=26): s1 ∈ {0..3, 15} xor 26 = {21, 24..27}.
     * Those indices must contribute 0 to match scalar's tail semantics. */
    temp8[21] = 0;
    temp8[24] = temp8[25] = temp8[26] = temp8[27] = 0;
    int8x16x2_t perm_vec;
    perm_vec.val[0] = vld1q_s8(temp8);
    perm_vec.val[1] = vld1q_s8(temp8 + 16);

    /* Per-lane minsc / endsc + msk bitmasks. Same convention as the 8-bit
     * kernel: minsc_msk_a bit set iff lane has a real KSW_XSUBO target that
     * fits in s16; same for endsc_msk_a / KSW_XSTOP. */
    uint8_t minsc_msk_a = 0x00, endsc_msk_a = 0x00;
    for (int i = 0; i < SIMD_WIDTH16; i++) {
        int xtra = p[i].h0;
        int val = (xtra & KSW_XSUBO) ? xtra & 0xffff : 0x10000;
        if (val <= SHRT_MAX) { minsc[i] = (int16_t)val; minsc_msk_a |= (0x1 << i); }
        val = (xtra & KSW_XSTOP) ? xtra & 0xffff : 0x10000;
        if (val <= SHRT_MAX) { endsc[i] = (int16_t)val; endsc_msk_a |= (0x1 << i); }
    }

    int16x8_t minsc_vec = vld1q_s16(minsc);
    int16x8_t endsc_vec = vld1q_s16(endsc);

    int16x8_t e_del_vec  = vdupq_n_s16(this->e_del);
    int16x8_t oe_del_vec = vdupq_n_s16(this->o_del + this->e_del);
    int16x8_t e_ins_vec  = vdupq_n_s16(this->e_ins);
    int16x8_t oe_ins_vec = vdupq_n_s16(this->o_ins + this->e_ins);

    int16x8_t gmax_vec = zero_vec;
    int16x8_t te_vec   = vdupq_n_s16(-1);
    int16x8_t qe_vec   = zero_vec;

    /* Per-lane freeze once a pair hits its KSW_XSTOP target. Gated by
     * has_endsc_vec so lanes with endsc[i] == 0 (no XSTOP) never freeze. */
    int16x8_t frozen_vec = zero_vec;
    int16_t _has_endsc_arr[SIMD_WIDTH16] __attribute((aligned(128)));
    for (int i = 0; i < SIMD_WIDTH16; i++)
        _has_endsc_arr[i] = (endsc_msk_a & (1 << i)) ? (int16_t)0xFFFF : (int16_t)0x0000;
    int16x8_t has_endsc_vec = vld1q_s16(_has_endsc_arr);

    tid = 0;
    int16_t *H0     = H16_0    + tid * SIMD_WIDTH16 * this->maxQerLen;
    int16_t *H1     = H16_1    + tid * SIMD_WIDTH16 * this->maxQerLen;
    int16_t *Hmax   = H16_max  + tid * SIMD_WIDTH16 * this->maxQerLen;
    int16_t *F      = F16      + tid * SIMD_WIDTH16 * this->maxQerLen;
    int16_t *rowMax = rowMax16 + tid * SIMD_WIDTH16 * this->maxRefLen;

    for (int i = 0; i <= ncol; i++) {
        vst1q_s16(H0   + i * SIMD_WIDTH16, zero_vec);
        vst1q_s16(Hmax + i * SIMD_WIDTH16, zero_vec);
        vst1q_s16(F    + i * SIMD_WIDTH16, zero_vec);
    }
    vst1q_s16(H0, zero_vec);
    vst1q_s16(H1, zero_vec);

    int16x8_t imax_vec, pimax_vec = zero_vec;

    int i, limit = nrow;
    for (i = 0; i < nrow; i++) {
        int16x8_t e11 = zero_vec;
        int16x8_t s1  = vld1q_s16(seq1SoA + i * SIMD_WIDTH16);
        imax_vec = zero_vec;
        int16x8_t iqe_vec = vdupq_n_s16(-1);
        int16x8_t l_vec   = zero_vec;
        int16x8_t i_vec   = vdupq_n_s16((int16_t)i);

        for (int j = 0; j < ncol; j++) {
            int16x8_t h00 = vld1q_s16(H0 + j * SIMD_WIDTH16);
            int16x8_t s2  = vld1q_s16(seq2SoA + j * SIMD_WIDTH16);
            int16x8_t f11 = vld1q_s16(F + (j + 1) * SIMD_WIDTH16);

            int16x8_t xor_val = veorq_s16(s1, s2);
            /* 8-lane int16 table lookup. Narrow xor (int16x8, values in
             * [0..31]) to 8 u8 indices (low byte of each int16; high byte
             * is always 0 for these xor values). Look up in the 32-byte
             * int8 score table via vqtbl2, sign-extend int8x8 back to
             * int16x8. */
            uint8x8_t  idx8 = vmovn_u16(vreinterpretq_u16_s16(xor_val));
            int8x8_t   sbt8 = vqtbl2_s8(perm_vec, idx8);
            int16x8_t  sbt  = vmovl_s8(sbt8);

            /* Boundary: high bit set in (s1 | s2) indicates padding. */
            int16x8_t or_val = vorrq_s16(s1, s2);
            uint16x8_t high_bit = vshrq_n_u16(vreinterpretq_u16_s16(or_val), 15);
            uint16x8_t is_boundary = vceqq_u16(high_bit, vdupq_n_u16(1));

            int16x8_t m11 = vaddq_s16(h00, sbt);
            m11 = vbslq_s16(is_boundary, zero_vec, m11);

            int16x8_t h11 = vmaxq_s16(m11, e11);
            h11 = vmaxq_s16(h11, f11);
            h11 = vmaxq_s16(h11, zero_vec);

            uint16x8_t cmp0 = vcgtq_s16(h11, imax_vec);
            imax_vec = vmaxq_s16(imax_vec, h11);
            iqe_vec  = vbslq_s16(cmp0, l_vec, iqe_vec);

            int16x8_t gapE = vsubq_s16(h11, oe_ins_vec);
            e11 = vsubq_s16(e11, e_ins_vec);
            e11 = vmaxq_s16(gapE, e11);

            int16x8_t gapD = vsubq_s16(h11, oe_del_vec);
            int16x8_t f21  = vsubq_s16(f11, e_del_vec);
            f21 = vmaxq_s16(gapD, f21);

            vst1q_s16(H1 + (j + 1) * SIMD_WIDTH16, h11);
            vst1q_s16(F  + (j + 1) * SIMD_WIDTH16, f21);
            l_vec = vaddq_s16(l_vec, one_vec);
        }

        /* Block I: write prior-row's pimax to rowMax. No intermediate
         * masking — frozen lanes' imax after the freeze row is still
         * produced by the DP but the score2 scan filters by per-lane
         * primary region / minsc / qe, which matches scalar. */
        if (i > 0) {
            vst1q_s16(rowMax + (i - 1) * SIMD_WIDTH16, pimax_vec);
        }
        pimax_vec = imax_vec;

        /* Block II: gmax / te / qe update with per-lane freeze mask. */
        uint16x8_t cmp0 = vcgtq_s16(imax_vec, gmax_vec);
        uint16x8_t frozen_u16 = vreinterpretq_u16_s16(frozen_vec);
        uint16x8_t cmp0_active = vbicq_u16(cmp0, frozen_u16);

        int16x8_t new_gmax = vmaxq_s16(gmax_vec, imax_vec);
        gmax_vec = vbslq_s16(frozen_u16, gmax_vec, new_gmax);

        te_vec = vbslq_s16(cmp0_active, i_vec, te_vec);
        qe_vec = vbslq_s16(cmp0_active, iqe_vec, qe_vec);

        /* Freeze newly endsc-qualifying lanes (has_endsc gate prevents
         * trivial freeze for lanes without a KSW_XSTOP target). */
        uint16x8_t cmp_end   = vcgeq_s16(gmax_vec, endsc_vec);
        uint16x8_t just_hit  = vandq_u16(cmp_end, vreinterpretq_u16_s16(has_endsc_vec));
        frozen_vec = vorrq_s16(frozen_vec, vreinterpretq_s16_u16(just_hit));

        /* Collapse frozen_vec (lanes are 0x0000 or 0xFFFF) to an 8-bit
         * mask: bit l set iff lane l has hit its KSW_XSTOP target. Uses
         * the existing _mm_movemask_epi16 helper (extracts each u16
         * MSB), avoiding the scalar store+loop round-trip. */
        uint8_t frozen_bits = (uint8_t)_mm_movemask_epi16(vreinterpretq_m128i_s16(frozen_vec));

        /* Early exit only when every lane that *could* freeze has frozen
         * — i.e. all KSW_XSTOP-carrying lanes are done. Non-XSTOP lanes
         * (endsc_msk_a bit clear) never contribute to frozen_bits, so the
         * batched break matches scalar ksw_i16 (no batched global exit). */
        if (endsc_msk_a != 0 && (frozen_bits & endsc_msk_a) == endsc_msk_a) {
            limit = i++;
            break;
        }

        int16_t *S = H1; H1 = H0; H0 = S;
    }

    /* Store final row's pimax. Guard on i > 0 for all-padding batches
     * (nrow == 0 → i stays at 0 → underflow). See issue 38 / PR 289. */
    if (i > 0) {
        vst1q_s16(rowMax + (i - 1) * SIMD_WIDTH16, pimax_vec);
    }

    /* Extract primary results. */
    int16_t score[SIMD_WIDTH16] __attribute((aligned(128)));
    int16_t te1[SIMD_WIDTH16]   __attribute((aligned(128)));
    int16_t qe[SIMD_WIDTH16]    __attribute((aligned(128)));
    vst1q_s16(score, gmax_vec);
    vst1q_s16(te1,   te_vec);
    vst1q_s16(qe,    qe_vec);

    int live = 0;
    for (int l = 0; l < SIMD_WIDTH16 && (po_ind + l) < numPairs; l++) {
        int ind = po_ind + l;
#if !MAINY
        ind = p[l].regid;
        if (phase) {
            if (aln[ind].score == score[l]) {
                aln[ind].tb = aln[ind].te - te1[l];
                aln[ind].qb = aln[ind].qe - qe[l];
            }
        } else {
            aln[ind].score = score[l];
            aln[ind].te    = te1[l];
            aln[ind].qe    = qe[l];
            if (score[l] > 0) { qe[l] = 1; live++; } else qe[l] = 0;
        }
#else
        aln[ind].score = score[l];
        aln[ind].te    = te1[l];
        aln[ind].qe    = qe[l];
        if (score[l] > 0) { qe[l] = 1; live++; } else qe[l] = 0;
#endif
    }

#if !MAINY
    if (phase) return 1;
#endif
    if (live == 0) return 1;

    /* Score2 / te2 via per-lane scalar b[]-emulation over the now-
     * populated rowMax. Identical semantics to the 8-bit kernel's
     * score2 scan; differs only in int16_t types and SIMD_WIDTH16 stride. */
    int qmax = this->g_qmax;
    int16_t low[SIMD_WIDTH16]  __attribute((aligned(128)));
    int16_t high[SIMD_WIDTH16] __attribute((aligned(128)));
    for (int j = 0; j < SIMD_WIDTH16; j++) {
        int val = (score[j] + qmax - 1) / qmax;
        low[j]  = te1[j] - val;
        high[j] = te1[j] + val;
    }

    /* rowMax entries produced: i == nrow on normal completion, and
     * i == limit + 1 after early exit via KSW_XSTOP (limit = i++ stores
     * the old i). Using i directly covers both paths and includes the
     * terminating row stored at rowMax[i-1]. */
    const int processed_rows = i;

    for (int l = 0; l < SIMD_WIDTH16 && (po_ind + l) < numPairs; l++) {
        int ind = p[l].regid;
        /* Match scalar ksw_i16: when KSW_XSUBO is absent, minsc = 0x10000
         * so the b[] list never starts and score2/te2 stay at -1.
         * minsc_msk_a tracks lanes with KSW_XSUBO set; lanes not in the
         * mask have minsc[l] zero-initialized and must skip the scan.
         *
         * The !qe[l] liveness check is equivalent to kswv512_16's
         * te[l] < 0: the extraction block above overwrites qe[l] with
         * 1 (score > 0) or 0 as a liveness flag before this scan runs,
         * so !qe[l] cleanly marks lanes with no primary alignment. Do
         * not rewrite to te[l] < 0 without also reverting that
         * liveness-flag overwrite. */
        if (!qe[l] || !(minsc_msk_a & (1u << l))) {
            aln[ind].score2 = -1;
            aln[ind].te2    = -1;
            continue;
        }

        int len1_l  = (int)p[l].len1;
        int low_l   = (int)low[l];
        int high_l  = (int)high[l];
        int minsc_l = (int)minsc[l];
        int score2  = -1;
        int te2     = -1;
        int b_score = -1;
        int b_pos   = -2;

        int nrows = processed_rows < len1_l ? processed_rows : len1_l;
        for (int i2 = 0; i2 < nrows; i2++) {
            int imax = (int)rowMax[i2 * SIMD_WIDTH16 + l];
            if (imax < minsc_l) continue;

            if (b_pos + 1 != i2) {
                if (b_pos >= 0 &&
                    (b_pos < low_l || b_pos > high_l) &&
                    b_score > score2) {
                    score2 = b_score;
                    te2    = b_pos;
                }
                b_score = imax;
                b_pos   = i2;
            } else if (b_score < imax) {
                b_score = imax;
                b_pos   = i2;
            }
        }

        if (b_pos >= 0 &&
            (b_pos < low_l || b_pos > high_l) &&
            b_score > score2) {
            score2 = b_score;
            te2    = b_pos;
        }

        aln[ind].score2 = score2;
        aln[ind].te2    = te2;
    }

    return 1;
}

#elif ((!__AVX512BW__) & (__AVX2__))
/* AVX2 kswv kernel — 256-bit vectors, 32 u8 lanes per batch.
 *
 * Direct port of the corrected NEON kernel (kswv_neon_u8 above) with all
 * four bug fixes pre-applied (te split, freeze mask, per-lane score2
 * exclusion, minsc filter in score2 scan).
 *
 * Intrinsic translation notes:
 *   - No unsigned 8-bit compare in AVX2 → sign-flip trick:
 *       cmpgt_u8(a, b) = cmpgt_i8(a ^ 0x80, b ^ 0x80)
 *   - No 8-bit shift → isolate high bit via sign compare against zero.
 *   - 16-byte lookup table → broadcast to 256 bits then shuffle
 *     (shuffle_epi8 is 128-bit-lane; broadcast makes both lanes
 *     identical so any 0..15 index lands on the correct byte).
 *   - te uses SIMD_WIDTH8=32 lanes worth of int16 → two __m256i
 *     (te_vec_lo = pairs 0..15, te_vec_hi = pairs 16..31).
 */

#include <immintrin.h>

static inline __m256i avx2_cmpgt_u8(__m256i a, __m256i b)
{
    const __m256i sign = _mm256_set1_epi8((char)0x80);
    return _mm256_cmpgt_epi8(_mm256_xor_si256(a, sign),
                             _mm256_xor_si256(b, sign));
}

static inline __m256i avx2_cmpge_u8(__m256i a, __m256i b)
{
    /* a >= b  <=>  b <= a  <=>  max(a,b) == a */
    return _mm256_cmpeq_epi8(_mm256_max_epu8(a, b), a);
}

/* NEON-style select: pick src where mask byte is 0xFF, else dst.
 * _mm256_blendv_epi8(dst, src, mask) uses the SIGN bit of each mask
 * byte: 0xFF (sign=1) picks src, 0x00 picks dst. Matches NEON's
 * vbslq_u8(mask, src, dst). */
static inline __m256i avx2_blendv_u8(__m256i mask, __m256i src, __m256i dst)
{
    return _mm256_blendv_epi8(dst, src, mask);
}

/* Widen the low/high 128 bits of a u8x32 into u16x16, zero-extended. */
static inline __m256i avx2_widen_u8_lo(__m256i v)
{
    return _mm256_cvtepu8_epi16(_mm256_extracti128_si256(v, 0));
}
static inline __m256i avx2_widen_u8_hi(__m256i v)
{
    return _mm256_cvtepu8_epi16(_mm256_extracti128_si256(v, 1));
}

/* Signed 16-bit compares are provided directly by AVX2 (_mm256_cmpgt_epi16,
 * _mm256_cmpeq_epi16). No cmplt / cmpge; synthesize. */
static inline __m256i avx2_cmplt_s16(__m256i a, __m256i b)
{
    return _mm256_cmpgt_epi16(b, a);
}
static inline __m256i avx2_cmpgt_s16(__m256i a, __m256i b)
{
    return _mm256_cmpgt_epi16(a, b);
}

int kswv::kswv256_u8(uint8_t seq1SoA[],
                     uint8_t seq2SoA[],
                     int16_t nrow,
                     int16_t ncol,
                     SeqPair *p,
                     kswr_t *aln,
                     int po_ind,
                     uint16_t tid,
                     int32_t numPairs,
                     int phase)
{
    uint8_t minsc[SIMD_WIDTH8] __attribute__((aligned(64))) = {0};
    uint8_t endsc[SIMD_WIDTH8] __attribute__((aligned(64))) = {0};

    const __m256i zero_vec = _mm256_setzero_si256();
    const __m256i one_vec  = _mm256_set1_epi8(1);

    /* Score lookup table (16 bytes). Indexed by (s1 ^ s2). */
    int8_t temp[SIMD_WIDTH8] __attribute__((aligned(64))) = {0};

    uint8_t shift;
    shift = min_(this->w_match, (int8_t)this->w_mismatch);
    shift = min_((int8_t)shift, this->w_ambig);
    shift = 256 - (uint8_t)shift;

    temp[0] = this->w_match;
    temp[1] = temp[2] = temp[3] = this->w_mismatch;
    temp[4] = temp[5] = temp[6] = temp[7] = this->w_ambig;
    temp[8] = temp[9] = temp[10] = temp[11] = this->w_ambig;
    temp[12] = this->w_ambig;
    for (int i = 0; i < 16; i++) temp[i] += shift;

    /* Broadcast 16-byte table across both 128-bit lanes of a 256-bit
     * vector so _mm256_shuffle_epi8 (128-bit-lane) still hits the right
     * entry for any 0..15 index. */
    __m128i permSft_128 = _mm_loadu_si128((const __m128i*)temp);
    __m256i permSft     = _mm256_broadcastsi128_si256(permSft_128);
    __m256i sft_vec     = _mm256_set1_epi8((char)shift);

    uint32_t minsc_msk_a = 0, endsc_msk_a = 0;
    for (int i = 0; i < SIMD_WIDTH8; i++) {
        int xtra = p[i].h0;
        int val  = (xtra & KSW_XSUBO) ? (xtra & 0xffff) : 0x10000;
        if (val <= 255) { minsc[i] = (uint8_t)val; minsc_msk_a |= (1u << i); }
        val = (xtra & KSW_XSTOP) ? (xtra & 0xffff) : 0x10000;
        if (val <= 255) { endsc[i] = (uint8_t)val; endsc_msk_a |= (1u << i); }
    }

    const __m256i minsc_vec  = _mm256_loadu_si256((const __m256i*)minsc);
    const __m256i endsc_vec  = _mm256_loadu_si256((const __m256i*)endsc);
    const __m256i e_del_vec  = _mm256_set1_epi8((char)this->e_del);
    const __m256i oe_del_vec = _mm256_set1_epi8((char)(this->o_del + this->e_del));
    const __m256i e_ins_vec  = _mm256_set1_epi8((char)this->e_ins);
    const __m256i oe_ins_vec = _mm256_set1_epi8((char)(this->o_ins + this->e_ins));
    const __m256i five_vec   = _mm256_set1_epi8((char)DUMMY5);
    const __m256i cmax_vec   = _mm256_set1_epi8((char)255);

    __m256i gmax_vec   = zero_vec;
    /* Fix 1: split te into two int16 vectors covering 32 pairs. */
    __m256i te_vec_lo  = _mm256_set1_epi16(-1);  /* pairs 0..15 */
    __m256i te_vec_hi  = _mm256_set1_epi16(-1);  /* pairs 16..31 */
    /* Fix 2: per-lane freeze once pair hits KSW_XSTOP. */
    __m256i frozen_vec = zero_vec;

    /* Fix 2 helper: only freeze lanes that actually set an endsc target.
     * endsc_vec defaults to 0 for lanes without a target (compares
     * trivially true with vcgeq_u8) — has_endsc_vec masks that out. */
    uint8_t _has_endsc_bytes[SIMD_WIDTH8];
    for (int i = 0; i < SIMD_WIDTH8; i++)
        _has_endsc_bytes[i] = (endsc_msk_a & (1u << i)) ? 0xFF : 0x00;
    __m256i has_endsc_vec = _mm256_loadu_si256((const __m256i*)_has_endsc_bytes);

    uint32_t exit0 = 0xFFFFFFFFu;

    tid = 0;
    uint8_t *H0     = H8_0    + tid * SIMD_WIDTH8 * this->maxQerLen;
    uint8_t *H1     = H8_1    + tid * SIMD_WIDTH8 * this->maxQerLen;
    uint8_t *Hmax   = H8_max  + tid * SIMD_WIDTH8 * this->maxQerLen;
    uint8_t *F      = F8      + tid * SIMD_WIDTH8 * this->maxQerLen;
    uint8_t *rowMax = rowMax8 + tid * SIMD_WIDTH8 * this->maxRefLen;

    for (int i = 0; i <= ncol; i++) {
        _mm256_storeu_si256((__m256i*)(H0   + i * SIMD_WIDTH8), zero_vec);
        _mm256_storeu_si256((__m256i*)(Hmax + i * SIMD_WIDTH8), zero_vec);
        _mm256_storeu_si256((__m256i*)(F    + i * SIMD_WIDTH8), zero_vec);
    }
    _mm256_storeu_si256((__m256i*)H0, zero_vec);
    _mm256_storeu_si256((__m256i*)H1, zero_vec);

    __m256i pimax_vec = zero_vec;
    uint32_t mask32   = 0;

    __m256i imax_vec;
    __m256i qe_vec = zero_vec;

    int i, limit = nrow;
    for (i = 0; i < nrow; i++) {
        __m256i e11 = zero_vec;
        __m256i s1  = _mm256_loadu_si256((const __m256i*)(seq1SoA + i * SIMD_WIDTH8));
        imax_vec    = zero_vec;
        __m256i iqe_vec = _mm256_set1_epi8((char)0xFF);
        __m256i l_vec = zero_vec;
        __m256i i_vec_s16 = _mm256_set1_epi16((int16_t)i);

        for (int j = 0; j < ncol; j++) {
            __m256i h00 = _mm256_loadu_si256((const __m256i*)(H0 + j * SIMD_WIDTH8));
            __m256i s2  = _mm256_loadu_si256((const __m256i*)(seq2SoA + j * SIMD_WIDTH8));
            __m256i f11 = _mm256_loadu_si256((const __m256i*)(F + (j + 1) * SIMD_WIDTH8));

            __m256i xor_val = _mm256_xor_si256(s1, s2);
            __m256i sbt     = _mm256_shuffle_epi8(permSft, xor_val);

            __m256i cmpq = _mm256_cmpeq_epi8(s2, five_vec);
            sbt = avx2_blendv_u8(cmpq, sft_vec, sbt);

            /* High bit of (s1 | s2) indicates boundary (padding 0xFF).
             * Sign compare against zero gives 0xFF wherever high bit is
             * set in the u8 byte. */
            __m256i or_val      = _mm256_or_si256(s1, s2);
            __m256i is_boundary = _mm256_cmpgt_epi8(zero_vec, or_val);

            __m256i m11 = _mm256_adds_epu8(h00, sbt);
            m11 = avx2_blendv_u8(is_boundary, zero_vec, m11);
            m11 = _mm256_subs_epu8(m11, sft_vec);

            __m256i h11 = _mm256_max_epu8(m11, e11);
            h11 = _mm256_max_epu8(h11, f11);

            __m256i cmp0 = avx2_cmpgt_u8(h11, imax_vec);
            imax_vec = _mm256_max_epu8(imax_vec, h11);
            iqe_vec  = avx2_blendv_u8(cmp0, l_vec, iqe_vec);

            __m256i gapE = _mm256_subs_epu8(h11, oe_ins_vec);
            e11 = _mm256_subs_epu8(e11, e_ins_vec);
            e11 = _mm256_max_epu8(gapE, e11);

            __m256i gapD = _mm256_subs_epu8(h11, oe_del_vec);
            __m256i f21  = _mm256_subs_epu8(f11, e_del_vec);
            f21 = _mm256_max_epu8(gapD, f21);

            _mm256_storeu_si256((__m256i*)(H1 + (j + 1) * SIMD_WIDTH8), h11);
            _mm256_storeu_si256((__m256i*)(F  + (j + 1) * SIMD_WIDTH8), f21);
            l_vec = _mm256_add_epi8(l_vec, one_vec);
        }

        /* Block I - rowMax tracking. Same shape as NEON. */
        if (i > 0) {
            __m256i cmp_gt = avx2_cmpgt_u8(imax_vec, pimax_vec);
            uint32_t msk32 = (uint32_t)_mm256_movemask_epi8(cmp_gt);
            msk32 |= mask32;
            pimax_vec = avx2_blendv_u8(cmp_gt, zero_vec, pimax_vec);
            _mm256_storeu_si256((__m256i*)(rowMax + (i - 1) * SIMD_WIDTH8), pimax_vec);
            mask32 = ~msk32;
        }
        pimax_vec = imax_vec;

        /* Block II: gmax, te with freeze mask (fix 2) */
        __m256i cmp0 = avx2_cmpgt_u8(imax_vec, gmax_vec);
        uint32_t cmp0_msk = (uint32_t)_mm256_movemask_epi8(cmp0);
        cmp0_msk &= exit0;

        /* 16-bit mirrors of cmp0 for te update (fix 1 split) */
        __m256i cmp_lo_16 = avx2_cmpgt_s16(
            avx2_widen_u8_lo(imax_vec), avx2_widen_u8_lo(gmax_vec));
        __m256i cmp_hi_16 = avx2_cmpgt_s16(
            avx2_widen_u8_hi(imax_vec), avx2_widen_u8_hi(gmax_vec));

        /* 16-bit mirrors of frozen_vec to mask te updates. */
        __m256i frozen_lo_16 = avx2_cmpgt_s16(
            avx2_widen_u8_lo(frozen_vec), _mm256_setzero_si256());
        __m256i frozen_hi_16 = avx2_cmpgt_s16(
            avx2_widen_u8_hi(frozen_vec), _mm256_setzero_si256());

        cmp_lo_16 = _mm256_andnot_si256(frozen_lo_16, cmp_lo_16);
        cmp_hi_16 = _mm256_andnot_si256(frozen_hi_16, cmp_hi_16);
        __m256i cmp0_active = _mm256_andnot_si256(frozen_vec, cmp0);

        __m256i new_gmax = _mm256_max_epu8(gmax_vec, imax_vec);
        gmax_vec = avx2_blendv_u8(frozen_vec, gmax_vec, new_gmax);

        te_vec_lo = avx2_blendv_u8(cmp_lo_16, i_vec_s16, te_vec_lo);
        te_vec_hi = avx2_blendv_u8(cmp_hi_16, i_vec_s16, te_vec_hi);
        qe_vec    = avx2_blendv_u8(cmp0_active, iqe_vec, qe_vec);

        /* End-score check + freeze update. */
        __m256i cmp_end = avx2_cmpge_u8(gmax_vec, endsc_vec);
        uint32_t cmp_end_msk = (uint32_t)_mm256_movemask_epi8(cmp_end);
        cmp_end_msk &= endsc_msk_a;

        __m256i just_hit = _mm256_and_si256(cmp_end, has_endsc_vec);
        frozen_vec = _mm256_or_si256(frozen_vec, just_hit);

        /* Overflow check. */
        __m256i left_vec = _mm256_adds_epu8(gmax_vec, sft_vec);
        __m256i cmp2     = avx2_cmpge_u8(left_vec, cmax_vec);
        uint32_t cmp2_msk = (uint32_t)_mm256_movemask_epi8(cmp2);

        exit0 = (~(cmp_end_msk | cmp2_msk)) & exit0;
        if (exit0 == 0) { limit = i++; break; }

        uint8_t *S = H1; H1 = H0; H0 = S;
    }

    /* Guard on i > 0: when every pair in the batch has len1 == 0, nrow
     * is 0, the DP loop never runs, and `i - 1` underflows into the
     * rowMax allocation prefix. See issue 38 / upstream PR 289. */
    if (i > 0) {
        _mm256_storeu_si256((__m256i*)(rowMax + (i - 1) * SIMD_WIDTH8), pimax_vec);
    }

    /* Extract results. */
    uint8_t score[SIMD_WIDTH8] __attribute__((aligned(64)));
    int16_t te1[SIMD_WIDTH8]    __attribute__((aligned(64)));
    uint8_t qe[SIMD_WIDTH8]     __attribute__((aligned(64)));

    _mm256_storeu_si256((__m256i*)score, gmax_vec);
    _mm256_storeu_si256((__m256i*)te1, te_vec_lo);        /* pairs 0..15 */
    _mm256_storeu_si256((__m256i*)(te1 + 16), te_vec_hi); /* pairs 16..31 */
    _mm256_storeu_si256((__m256i*)qe, qe_vec);

    int live = 0;
    for (int l = 0; l < SIMD_WIDTH8 && (po_ind + l) < numPairs; l++) {
        int ind = p[l].regid;
        int16_t *te = te1;
        if (phase) {
            if (aln[ind].score == score[l]) {
                aln[ind].tb = aln[ind].te - te[l];
                aln[ind].qb = aln[ind].qe - qe[l];
            }
        } else {
            aln[ind].score = score[l] + shift < 255 ? score[l] : 255;
            aln[ind].te = te[l];
            aln[ind].qe = qe[l];
            if (aln[ind].score != 255) { qe[l] = 1; live++; }
            else qe[l] = 0;
        }
    }

    if (phase) return 1;
    if (live == 0) return 1;

    /* Score2 and te2 computation.
     *
     * The SIMD dense-rowMax approach (max2 over every row outside the
     * primary region) diverges from scalar ksw_u8 on plateaus: scalar
     * collapses consecutive rows >= minsc into a single b[] entry
     * anchored at the max-score row, so a plateau straddling the primary
     * region boundary (anchor inside, later rows outside) is excluded
     * entirely. Dense rowMax would otherwise pull partial-plateau rows
     * in and inflate `sub`. Emulate scalar exactly, per lane, on the
     * precomputed rowMax. Cost is O(SIMD_WIDTH8 * limit) scalar ops per
     * batch, negligible next to the SIMD DP. */
    int qmax = this->g_qmax;
    int16_t low[SIMD_WIDTH8]  __attribute__((aligned(64)));
    int16_t high[SIMD_WIDTH8] __attribute__((aligned(64)));
    for (int j = 0; j < SIMD_WIDTH8; j++) {
        int val = (score[j] + qmax - 1) / qmax;
        low[j]  = te1[j] - val;
        high[j] = te1[j] + val;
    }

    /* rowMax entries produced: i == limit on normal completion, and
     * i == limit + 1 after early exit (limit = i++ stores old i). Using i
     * directly gives a single expression that covers both exit paths and
     * includes the terminating row stored at rowMax[limit]. */
    const int processed_rows = i;

    for (int l = 0; l < SIMD_WIDTH8 && (po_ind + l) < numPairs; l++) {
        int ind = p[l].regid;
        /* Match scalar ksw_u8: when KSW_XSUBO is absent, minsc = 0x10000
         * so the b[] list never starts and score2/te2 stay at -1.
         * minsc_msk_a tracks lanes with KSW_XSUBO set; lanes not in the
         * mask have minsc[l] zero-initialized and must skip the scan. */
        if (!qe[l] || !(minsc_msk_a & (1u << l))) {
            aln[ind].score2 = -1;
            aln[ind].te2    = -1;
            continue;
        }

        int len1_l  = (int)p[l].len1;
        int low_l   = (int)low[l];
        int high_l  = (int)high[l];
        int minsc_l = (int)minsc[l];
        int score2  = -1;
        int te2     = -1;
        /* b_pos = -2 sentinel: (-2 + 1) != any real row index, forces
         * an "append" on first qualifying row. */
        int b_score = -1;
        int b_pos   = -2;

        int nrows = processed_rows < len1_l ? processed_rows : len1_l;
        for (int i2 = 0; i2 < nrows; i2++) {
            int imax = (int)rowMax[i2 * SIMD_WIDTH8 + l];
            if (imax < minsc_l) continue;

            if (b_pos + 1 != i2) {
                /* APPEND: flush the outgoing b[] entry to score2 first. */
                if (b_pos >= 0 &&
                    (b_pos < low_l || b_pos > high_l) &&
                    b_score > score2) {
                    score2 = b_score;
                    te2    = b_pos;
                }
                b_score = imax;
                b_pos   = i2;
            } else if (b_score < imax) {
                /* UPDATE: strict greater extends the run's anchor. Equal
                 * imax leaves b_pos in place (matches scalar). */
                b_score = imax;
                b_pos   = i2;
            }
        }

        /* Flush trailing b[] entry. */
        if (b_pos >= 0 &&
            (b_pos < low_l || b_pos > high_l) &&
            b_score > score2) {
            score2 = b_score;
            te2    = b_pos;
        }

        aln[ind].score2 = score2;
        aln[ind].te2    = te2;
    }

    return 1;
}

void kswv::kswvBatchWrapper8_avx2(SeqPair *pairArray,
                                  uint8_t *seqBufRef,
                                  uint8_t *seqBufQer,
                                  kswr_t* aln,
                                  int32_t numPairs,
                                  uint16_t numThreads,
                                  int phase)
{
    uint8_t *seq1SoA = (uint8_t*)_mm_malloc(
        (size_t)this->maxRefLen * SIMD_WIDTH8 * numThreads * sizeof(uint8_t), 128);
    uint8_t *seq2SoA = (uint8_t*)_mm_malloc(
        (size_t)this->maxQerLen * SIMD_WIDTH8 * numThreads * sizeof(uint8_t), 128);
    assert(seq1SoA && seq2SoA);

    int32_t roundNumPairs = ((numPairs + SIMD_WIDTH8 - 1) / SIMD_WIDTH8) * SIMD_WIDTH8;
    for (int32_t ii = numPairs; ii < roundNumPairs; ii++) {
        pairArray[ii].regid = ii;
        pairArray[ii].id    = ii;
        pairArray[ii].len1  = 0;
        pairArray[ii].len2  = 0;
    }

    uint16_t tid = 0;
    uint8_t *mySeq1SoA = seq1SoA + tid * this->maxRefLen * SIMD_WIDTH8;
    uint8_t *mySeq2SoA = seq2SoA + tid * this->maxQerLen * SIMD_WIDTH8;

    for (int32_t i = 0; i < numPairs; i += SIMD_WIDTH8) {
        int maxLen1 = 0, maxLen2 = 0;

        for (int j = 0; j < SIMD_WIDTH8; j++) {
            SeqPair sp = pairArray[i + j];
            uint8_t *seq1 = seqBufRef + sp.idr;
            for (int k = 0; k < sp.len1; k++)
                mySeq1SoA[k * SIMD_WIDTH8 + j] = (seq1[k] == AMBIG_ ? AMBR : seq1[k]);
            if (maxLen1 < sp.len1) maxLen1 = sp.len1;
        }
        for (int j = 0; j < SIMD_WIDTH8; j++) {
            SeqPair sp = pairArray[i + j];
            for (int k = sp.len1; k <= maxLen1; k++)
                mySeq1SoA[k * SIMD_WIDTH8 + j] = 0xFF;
        }

        for (int j = 0; j < SIMD_WIDTH8; j++) {
            SeqPair sp = pairArray[i + j];
            uint8_t *seq2 = seqBufQer + sp.idq;
            int quanta = ((sp.len2 + 16 - 1) / 16) * 16;
            for (int k = 0; k < sp.len2; k++)
                mySeq2SoA[k * SIMD_WIDTH8 + j] = (seq2[k] == AMBIG_ ? AMBQ : seq2[k]);
            for (int k = sp.len2; k < quanta; k++)
                mySeq2SoA[k * SIMD_WIDTH8 + j] = DUMMY5;
            if (maxLen2 < quanta) maxLen2 = quanta;
        }
        for (int j = 0; j < SIMD_WIDTH8; j++) {
            SeqPair sp = pairArray[i + j];
            int quanta = ((sp.len2 + 16 - 1) / 16) * 16;
            for (int k = quanta; k <= maxLen2; k++)
                mySeq2SoA[k * SIMD_WIDTH8 + j] = 0xFF;
        }

        kswv256_u8(mySeq1SoA, mySeq2SoA,
                   (int16_t)maxLen1, (int16_t)maxLen2,
                   pairArray + i, aln, i, tid,
                   numPairs, phase);
    }

    _mm_free(seq1SoA);
    _mm_free(seq2SoA);
}

void kswv::getScores8(SeqPair *pairArray,
                      uint8_t *seqBufRef,
                      uint8_t *seqBufQer,
                      kswr_t *aln,
                      int32_t numPairs,
                      uint16_t numThreads,
                      int phase)
{
    kswvBatchWrapper8_avx2(pairArray, seqBufRef, seqBufQer, aln,
                           numPairs, numThreads, phase);
}

/* getScores16 remains a scalar fallback for now. The 16-bit NEON kernel
 * also has a full SIMD port — a parallel AVX2 port is a mechanical
 * follow-up after 8-bit lands. Most production mate-rescue traffic takes
 * the 8-bit path (KSW_XBYTE set when l_ms * w_match < 250). */
void kswv::getScores16(SeqPair *pairArray,
                       uint8_t *seqBufRef,
                       uint8_t *seqBufQer,
                       kswr_t *aln,
                       int32_t numPairs,
                       uint16_t /*numThreads*/,
                       int phase)
{
    if (phase != 0) return;

    int8_t mat[25];
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            mat[i*5 + j] = (i == j) ? this->w_match : this->w_mismatch;
        }
        mat[i*5 + 4] = this->w_ambig;
        mat[4*5 + i] = this->w_ambig;
    }
    mat[24] = this->w_ambig;

    for (int i = 0; i < numPairs; i++) {
        SeqPair *p = pairArray + i;
        kswr_t  *myaln = aln + p->regid;
        uint8_t *target = seqBufRef + p->idr;
        uint8_t *query  = seqBufQer + p->idq;
        kswr_t ks = ksw_align2(p->len2, query, p->len1, target, 5, mat,
                               this->o_del, this->e_del,
                               this->o_ins, this->e_ins,
                               p->h0, nullptr);
        myaln->score  = ks.score;
        myaln->qe     = ks.qe;
        myaln->te     = ks.te;
        myaln->qb     = ks.qb;
        myaln->tb     = ks.tb;
        myaln->score2 = ks.score2;
        myaln->te2    = ks.te2;
    }
}

#elif __AVX512BW__
/* AVX-512 Implementation for x86 */
void kswv::getScores8(SeqPair *pairArray,
                      uint8_t *seqBufRef,
                      uint8_t *seqBufQer,
                      kswr_t* aln,                    
                      int32_t numPairs,
                      uint16_t numThreads,
                      int phase)
{
    kswvBatchWrapper8(pairArray, seqBufRef, seqBufQer, aln,
                      numPairs, numThreads, phase);
}

#define PFD_ 2
void kswv::kswvBatchWrapper8(SeqPair *pairArray,
                             uint8_t *seqBufRef,
                             uint8_t *seqBufQer,
                             kswr_t* aln,
                             int32_t numPairs,
                             uint16_t numThreads,
                             int phase)
{
    int64_t st1, st2, st3, st4, st5;
#if RDT
    st1 = __rdtsc();
#endif
    uint8_t *seq1SoA = NULL;
    seq1SoA = (uint8_t *)_mm_malloc(this->maxRefLen * SIMD_WIDTH8 * numThreads * sizeof(uint8_t), 64);
    
    uint8_t *seq2SoA = NULL;
    seq2SoA = (uint8_t *)_mm_malloc(this->maxQerLen * SIMD_WIDTH8 * numThreads * sizeof(uint8_t), 64);
    
    assert(seq1SoA != NULL);
    assert(seq2SoA != NULL);

    int32_t ii;
    int32_t roundNumPairs = ((numPairs + SIMD_WIDTH8 - 1) / SIMD_WIDTH8 ) * SIMD_WIDTH8;
    // assert(roundNumPairs < BATCH_SIZE * SEEDS_PER_READ);
    for(ii = numPairs; ii < roundNumPairs; ii++)
    {
        pairArray[ii].regid = ii;
        pairArray[ii].id = ii;
        pairArray[ii].len1 = 0;
        pairArray[ii].len2 = 0;
    }

#if RDT
    st2 = __rdtsc();
#endif
    
#if SORT_PAIRS     // disbaled in bwa-mem2 (only used in separate benchmark sw code)
    {
        // Sort the sequences according to decreasing order of lengths
        SeqPair *tempArray = (SeqPair *)_mm_malloc(SORT_BLOCK_SIZE * numThreads *
                                                   sizeof(SeqPair), 64);
        int16_t *hist = (int16_t *)_mm_malloc((this->maxQerLen + 32) * numThreads *
                                              sizeof(int16_t), 64);
        
        #pragma omp parallel num_threads(numThreads)
        {
            int32_t tid = omp_get_thread_num();
            SeqPair *myTempArray = tempArray + tid * SORT_BLOCK_SIZE;
            int16_t *myHist = hist + tid * (this->maxQerLen + 32);
            
            #pragma omp for
            for(ii = 0; ii < roundNumPairs; ii+=SORT_BLOCK_SIZE)
            {
                int32_t first, last;
                first = ii;
                last  = ii + SORT_BLOCK_SIZE;
                if(last > roundNumPairs) last = roundNumPairs;
                sortPairsLen(pairArray + first, last - first, myTempArray, myHist);
            }
        }
        _mm_free(hist);
    }
#endif

#if RDT
    st3 = __rdtsc();
#endif
    
    //#pragma omp parallel num_threads(numThreads)
    {
        int32_t i;
        // uint16_t tid = omp_get_thread_num();
        uint16_t tid = 0;   // forcing single thread
        uint8_t *mySeq1SoA = seq1SoA + tid * this->maxRefLen * SIMD_WIDTH8;
        uint8_t *mySeq2SoA = seq2SoA + tid * this->maxQerLen * SIMD_WIDTH8;
        uint8_t *seq1;
        uint8_t *seq2;
                
        int nstart = 0, nend = numPairs;
        
        // #pragma omp for schedule(dynamic, 128)
        for(i = nstart; i < nend; i+=SIMD_WIDTH8)
        {
            int32_t j, k;
            int maxLen1 = 0;
            int maxLen2 = 0;

            for(j = 0; j < SIMD_WIDTH8; j++)
            {
                SeqPair sp = pairArray[i + j];
#if MAINY               
                seq1 = seqBufRef + (int64_t)sp.id * this->maxRefLen;
#else
                seq1 = seqBufRef + sp.idr;
#endif
                for(k = 0; k < sp.len1; k++)
                {
                    mySeq1SoA[k * SIMD_WIDTH8 + j] = (seq1[k] == AMBIG_ ? AMBR:seq1[k]);
                }
                if(maxLen1 < sp.len1) maxLen1 = sp.len1;
            }
            for(j = 0; j < SIMD_WIDTH8; j++)
            {
                SeqPair sp = pairArray[i + j];
                for(k = sp.len1; k <= maxLen1; k++) //removed "="
                {
                    mySeq1SoA[k * SIMD_WIDTH8 + j] = 0xFF;
                }
            }
            
            for(j = 0; j < SIMD_WIDTH8; j++)
            {               
                SeqPair sp = pairArray[i + j];
#if MAINY
                seq2 = seqBufQer + (int64_t)sp.id * this->maxQerLen;
#else
                seq2 = seqBufQer + sp.idq;
#endif
                // assert(sp.len2 < this->maxQerLen);
                int quanta = (sp.len2 + 16 - 1) / 16;  // based on SSE-8 bit lane
                quanta *= 16;                          // for matching the output of bwa-mem
                for(k = 0; k < sp.len2; k++)
                {
                    mySeq2SoA[k * SIMD_WIDTH8 + j] = (seq2[k]==AMBIG_? AMBQ:seq2[k]);
                }

                for(k = sp.len2; k < quanta; k++) {
                    mySeq2SoA[k * SIMD_WIDTH8 + j] = DUMMY5;  // SSE quanta
                }
                if(maxLen2 < (quanta)) maxLen2 = quanta;
            }
            
            for(j = 0; j < SIMD_WIDTH8; j++)
            {
                SeqPair sp = pairArray[i + j];
                int quanta = (sp.len2 + 16 - 1) / 16;  // based on SSE2-8 bit lane
                quanta *= 16;
                for(k = quanta; k <= maxLen2; k++)
                {
                    mySeq2SoA[k * SIMD_WIDTH8 + j] = 0xFF;
                }
            }

            kswv512_u8(mySeq1SoA, mySeq2SoA,
                       maxLen1, maxLen2,
                       pairArray + i,
                       aln, i,
                       tid,
                       numPairs,
                       phase);
        }
    }

#if RDT 
    st4 = __rdtsc();
#endif
    
#if SORT_PAIRS     // disbaled in bwa-mem2 (only used in separate benchmark sw code)
    {
        // Sort the sequences according to increasing order of id
        #pragma omp parallel num_threads(numThreads)
        {
            int32_t tid = omp_get_thread_num();
            SeqPair *myTempArray = tempArray + tid * SORT_BLOCK_SIZE;
            
            #pragma omp for
            for(ii = 0; ii < roundNumPairs; ii+=SORT_BLOCK_SIZE)
            {
                int32_t first, last;
                first = ii;
                last  = ii + SORT_BLOCK_SIZE;
                if(last > roundNumPairs) last = roundNumPairs;
                sortPairsId(pairArray + first, first, last - first, myTempArray);
            }
        }
        _mm_free(tempArray);
    }
#endif

#if RDT
    st5 = __rdtsc();
    setupTicks = st2 - st1;
    sort1Ticks = st3 - st2;
    swTicks = st4 - st3;
    sort2Ticks = st5 - st4;
#endif
    
    // free mem
    _mm_free(seq1SoA);
    _mm_free(seq2SoA);
    
    return;
}

int kswv::kswv512_u8(uint8_t seq1SoA[],
                     uint8_t seq2SoA[],
                     int16_t nrow,
                     int16_t ncol,
                     SeqPair *p,
                     kswr_t *aln,
                     int po_ind,
                     uint16_t tid,
                     int32_t numPairs,
                     int phase)
{
    int m_b, n_b;
    uint8_t minsc[SIMD_WIDTH8] __attribute__((aligned(64))) = {0};
    uint8_t endsc[SIMD_WIDTH8] __attribute__((aligned(64))) = {0};
    uint64_t *b;

    __m512i zero512 = _mm512_setzero_si512();
    __m512i one512  = _mm512_set1_epi8(1);
    
    m_b = n_b = 0; b = 0;

    int8_t temp[SIMD_WIDTH8] __attribute((aligned(64))) = {0};

    uint8_t shift = 127, mdiff = 0, qmax_;
    mdiff = max_(this->w_match, (int8_t) this->w_mismatch);
    mdiff = max_(mdiff, (int8_t) this->w_ambig);
    shift = min_(this->w_match, (int8_t) this->w_mismatch);
    shift = min_((int8_t) shift, this->w_ambig);

    qmax_ = mdiff;
    shift = 256 - (uint8_t) shift;
    mdiff += shift;
    
    temp[0] = this->w_match;                                   // states: 1. matches
    temp[1] = temp[2] = temp[3] =  this->w_mismatch;           // 2. mis-matches
    temp[4] = temp[5] = temp[6] = temp[7] =  this->w_ambig;    // 3. beyond boundary
    temp[8] = temp[9] = temp[10] = temp[11] = this->w_ambig;   // 4. 0 - sse2 region
    temp[12] = this->w_ambig;                                  // 5. ambig

    for (int i=0; i<16; i++) // for shuffle_epi8
        temp[i] += shift;

    int pos = 0;
    for (int i=16; i<SIMD_WIDTH8; i++) {
        temp[i] = temp[pos++];
        if (pos % 16 == 0) pos = 0;
    }
    
    __m512i permSft512 = _mm512_load_si512(temp);
    __m512i sft512 = _mm512_set1_epi8(shift);
    __m512i cmax512 = _mm512_set1_epi8(255);
    
    // __m512i minsc512, endsc512;
    __mmask64 minsc_msk_a = 0x0000, endsc_msk_a = 0x0000;
    int val = 0;
    for (int i=0; i<SIMD_WIDTH8; i++)
    {
        int xtra = p[i].h0;
        val = (xtra & KSW_XSUBO)? xtra & 0xffff : 0x10000;
        if (val <= 255) {
            minsc[i] = val;
            minsc_msk_a |= (0x1L << i);
        }
        // msc_mask;
        val = (xtra & KSW_XSTOP)? xtra & 0xffff : 0x10000;
        if (val <= 255) {
            endsc[i] = val;
            endsc_msk_a |= (0x1L << i);
        }
    }

    __m512i minsc512 = _mm512_load_si512((__m512i*) minsc);
    __m512i endsc512 = _mm512_load_si512((__m512i*) endsc);
       
    __m512i mismatch512 = _mm512_set1_epi8(this->w_mismatch + shift);
    __m512i e_del512    = _mm512_set1_epi8(this->e_del);
    __m512i oe_del512   = _mm512_set1_epi8(this->o_del + this->e_del);
    __m512i e_ins512    = _mm512_set1_epi8(this->e_ins);
    __m512i oe_ins512   = _mm512_set1_epi8(this->o_ins + this->e_ins);
    __m512i five512     = _mm512_set1_epi8(DUMMY5); // ambig mapping element
    __m512i gmax512     = zero512; // exit1 = zero512;
    __m512i te512       = _mm512_set1_epi16(-1);  // changed to -1
    __m512i te512_      = _mm512_set1_epi16(-1);  // changed to -1
    
    __mmask64 exit0 = 0xFFFFFFFFFFFFFFFF;

    tid = 0;  // no threading for now !!
    uint8_t *H0     = H8_0 + tid * SIMD_WIDTH8 * this->maxQerLen;
    uint8_t *H1     = H8_1 + tid * SIMD_WIDTH8 * this->maxQerLen;
    uint8_t *Hmax   = H8_max + tid * SIMD_WIDTH8 * this->maxQerLen;
    uint8_t *F      = F8 + tid * SIMD_WIDTH8 * this->maxQerLen;
    uint8_t *rowMax = rowMax8 + tid * SIMD_WIDTH8 * this->maxRefLen;
    
    
    for (int i=0; i <=ncol; i++)
    {
        _mm512_store_si512((__m512*) (H0 + i * SIMD_WIDTH8), zero512);
        _mm512_store_si512((__m512*) (Hmax + i * SIMD_WIDTH8), zero512);
        _mm512_store_si512((__m512*) (F + i * SIMD_WIDTH8), zero512);
    }

    __m512i max512 = zero512, imax512, pimax512 = zero512;
    __mmask64 mask512 = 0x0000;
    __mmask64 minsc_msk = 0x0000;

    __m512i qe512 = _mm512_set1_epi8(0);
    _mm512_store_si512((__m512i *)(H0), zero512);
    _mm512_store_si512((__m512i *)(H1), zero512);

    int i, limit = nrow;
    for (i=0; i < nrow; i++)
    {
        __m512i e11 = zero512;
        __m512i h00, h11, h10, s1;
        __m512i i512 = _mm512_set1_epi16(i);
        int j ;
        
        s1 = _mm512_load_si512((__m512i *)(seq1SoA + (i + 0) * SIMD_WIDTH8));
        h10 = zero512;
        imax512 = zero512;
        __m512i iqe512 = _mm512_set1_epi8(-1);

        __m512i l512 = zero512;
        for (j=0; j<ncol; j++)
        {
            __m512i f11, s2, f21;
            h00 = _mm512_load_si512((__m512i *)(H0 + j * SIMD_WIDTH8));  // check for col "0"
            s2  = _mm512_load_si512((__m512i *)(seq2SoA + (j) * SIMD_WIDTH8));
            f11 = _mm512_load_si512((__m512i *)(F + (j+1) * SIMD_WIDTH8));

            MAIN_SAM_CODE8_OPT(s1, s2, h00, h11, e11, f11, f21, max512, sft512);

            _mm512_store_si512((__m512i *)(H1 + (j + 1) * SIMD_WIDTH8), h11);  // check for col "0"
            _mm512_store_si512((__m512i *)(F + (j + 1)* SIMD_WIDTH8), f21);
            l512 = _mm512_add_epi8(l512, one512);
        }

        // Block I
        if (i > 0)
        {
            __mmask64 msk64 = _mm512_cmpgt_epu8_mask(imax512, pimax512);
            msk64 |= mask512;
            pimax512 = _mm512_mask_blend_epi8(msk64, pimax512, zero512);
            pimax512 = _mm512_mask_blend_epi8(minsc_msk, zero512, pimax512);
            pimax512 = _mm512_mask_blend_epi8(exit0, zero512, pimax512);
            
            _mm512_store_si512((__m512i *) (rowMax + (i-1)*SIMD_WIDTH8), pimax512);
            mask512 = ~msk64;
        }
        pimax512 = imax512;
        minsc_msk = _mm512_cmpge_epu8_mask(imax512, minsc512);
        minsc_msk &= minsc_msk_a;

        // Block II: gmax, te
        __mmask64 cmp0 = _mm512_cmpgt_epu8_mask(imax512, gmax512);
        cmp0 &= exit0;
        gmax512 = _mm512_mask_blend_epi8(cmp0, gmax512, imax512);
        te512 = _mm512_mask_blend_epi16((__mmask32)cmp0, te512,  i512);
        te512_ = _mm512_mask_blend_epi16((__mmask32) (cmp0 >> SIMD_WIDTH16), te512_,  i512);
        qe512 = _mm512_mask_blend_epi8(cmp0, qe512,  iqe512);
        
        cmp0 = _mm512_cmpge_epu8_mask(gmax512, endsc512);
        cmp0 &= endsc_msk_a;
        
        __m512i left512 = _mm512_adds_epu8(gmax512, sft512);
        __mmask64 cmp2 = _mm512_cmpge_epu8_mask(left512, cmax512);

        exit0 = (~(cmp0 | cmp2)) & exit0;
        if (exit0 == 0)
        {
            limit = i++;
            break;
        }       

        uint8_t *S = H1; H1 = H0; H0 = S;
        i512 = _mm512_add_epi16(i512, one512);
    } // for nrow

    pimax512 = _mm512_mask_blend_epi8(mask512, pimax512, zero512);
    pimax512 = _mm512_mask_blend_epi8(minsc_msk, zero512, pimax512);
    pimax512 = _mm512_mask_blend_epi8(exit0, zero512, pimax512);
    /* Guard on i > 0: all-len1==0 batches give nrow == 0 and would
     * otherwise underflow into the rowMax allocation. issue 38 / PR 289. */
    if (i > 0) {
        _mm512_store_si512((__m512i *) (rowMax + (i-1) * SIMD_WIDTH8), pimax512);
    }

    /******************* DP loop over *****************************/
    /**************** Partial output setting **********************/
    uint8_t score[SIMD_WIDTH8] __attribute((aligned(64)));
    int16_t te1[SIMD_WIDTH8] __attribute((aligned(64)));
    uint8_t qe[SIMD_WIDTH8] __attribute((aligned(64)));
    /* low[]/high[] are declared inside the score2 block below. */

    _mm512_store_si512((__m512i *) score, gmax512);
    _mm512_store_si512((__m512i *) te1, te512);
    _mm512_store_si512((__m512i *) (te1 + SIMD_WIDTH16), te512_);
    _mm512_store_si512((__m512i *) qe, qe512);

    int live = 0;
    for (int l=0; l<SIMD_WIDTH8 && (po_ind + l) < numPairs; l++) {
        int ind = po_ind + l;
        int16_t *te;
        if (i < SIMD_WIDTH16) te = te1;
        else te = te1;
#if !MAINY
        ind = p[l].regid;    // index of corr. aln
        if (phase) {
            if (aln[ind].score == score[l]) {
                aln[ind].tb = aln[ind].te - te[l];
                aln[ind].qb = aln[ind].qe - qe[l];
            }
        } else {
            aln[ind].score = score[l] + shift < 255? score[l] : 255;
            aln[ind].te = te[l];
            aln[ind].qe = qe[l];
            if (aln[ind].score != 255) {
                qe[l] = 1;
                live ++;                
            }
            else qe[l] = 0;
        }
#else
        aln[ind].score = score[l] + shift < 255? score[l] : 255;
        aln[ind].te = te[l];
        aln[ind].qe = qe[l];
        if (aln[ind].score != 255) {
            qe[l] = 1;
            live ++;                
        }
        else qe[l] = 0;     
#endif
    }
    
#if !MAINY
    if (phase) return 1;
#endif

    if (live == 0) return 1;

    /* Score2 and te2 computation.
     *
     * Per-lane scalar emulation of ksw_u8's b[] build + score2 scan
     * (src/ksw.cpp:201). scalar ksw_u8 collapses consecutive rows >=
     * minsc into a single b[] entry anchored at the max-score row; a
     * plateau straddling the primary region [te - val, te + val] is
     * excluded entirely when its anchor lies inside. A dense SIMD scan
     * treats every row independently and pulls out-of-region rows of a
     * boundary-straddling plateau in as false suboptimals, inflating `sub`.
     * Emulate scalar exactly, per lane, on the precomputed
     * rowMax. Cost is O(SIMD_WIDTH8 * limit) scalar ops per batch,
     * negligible next to the SIMD DP that built rowMax. */
    int qmax = this->g_qmax;
    int16_t low[SIMD_WIDTH8]  __attribute__((aligned(64)));
    int16_t high[SIMD_WIDTH8] __attribute__((aligned(64)));
    for (int j = 0; j < SIMD_WIDTH8; j++) {
        int val = (score[j] + qmax - 1) / qmax;
        low[j]  = te1[j] - val;
        high[j] = te1[j] + val;
    }

    /* rowMax entries produced: i == limit on normal completion, and
     * i == limit + 1 after early exit (limit = i++ stores old i). Using i
     * directly gives a single expression that covers both exit paths and
     * includes the terminating row stored at rowMax[limit]. */
    const int processed_rows = i;

    for (int l = 0; l < SIMD_WIDTH8 && (po_ind + l) < numPairs; l++) {
        int ind = p[l].regid;
        /* Match scalar ksw_u8: when KSW_XSUBO is absent, minsc = 0x10000
         * so the b[] list never starts and score2/te2 stay at -1.
         * minsc_msk_a tracks lanes with KSW_XSUBO set; lanes not in the
         * mask have minsc[l] zero-initialized and must skip the scan. */
        if (!qe[l] || !(minsc_msk_a & ((__mmask64)1 << l))) {
            aln[ind].score2 = -1;
            aln[ind].te2    = -1;
            continue;
        }

        int len1_l  = (int)p[l].len1;
        int low_l   = (int)low[l];
        int high_l  = (int)high[l];
        int minsc_l = (int)minsc[l];
        int score2  = -1;
        int te2     = -1;
        /* b_pos = -2 sentinel: (-2 + 1) != any real row, forces an
         * "append" on the first qualifying row. */
        int b_score = -1;
        int b_pos   = -2;

        int nrows = processed_rows < len1_l ? processed_rows : len1_l;
        for (int i2 = 0; i2 < nrows; i2++) {
            int imax = (int)rowMax[i2 * SIMD_WIDTH8 + l];
            if (imax < minsc_l) continue;

            if (b_pos + 1 != i2) {
                /* APPEND: flush the outgoing b[] entry to score2 first. */
                if (b_pos >= 0 &&
                    (b_pos < low_l || b_pos > high_l) &&
                    b_score > score2) {
                    score2 = b_score;
                    te2    = b_pos;
                }
                b_score = imax;
                b_pos   = i2;
            } else if (b_score < imax) {
                /* UPDATE: strict greater extends the run's anchor. Equal
                 * imax leaves b_pos in place (matches scalar). */
                b_score = imax;
                b_pos   = i2;
            }
        }

        /* Flush trailing b[] entry. */
        if (b_pos >= 0 &&
            (b_pos < low_l || b_pos > high_l) &&
            b_score > score2) {
            score2 = b_score;
            te2    = b_pos;
        }

        aln[ind].score2 = score2;
        aln[ind].te2    = te2;
    }

    return 1;   
}

/*********************************** Vectorized Code 16 bit *****************************/
/// 16 bit lanes
void kswv::getScores16(SeqPair *pairArray,
                       uint8_t *seqBufRef,
                       uint8_t *seqBufQer,
                       kswr_t* aln,
                       int32_t numPairs,
                       uint16_t numThreads,
                       int phase)
{
    kswvBatchWrapper16(pairArray, seqBufRef, seqBufQer, aln,
                       numPairs, numThreads, phase);
}

void kswv::kswvBatchWrapper16(SeqPair *pairArray,
                              uint8_t *seqBufRef,
                              uint8_t *seqBufQer,
                              kswr_t* aln,
                              int32_t numPairs,                           
                              uint16_t numThreads,
                              int phase)
{
    int64_t st1, st2, st3, st4, st5;
#if RDT
    st1 = __rdtsc();
#endif
    
    int16_t *seq1SoA = NULL;
    seq1SoA = (int16_t *)_mm_malloc(this->maxRefLen * SIMD_WIDTH16 * numThreads * sizeof(int16_t), 64);
    
    int16_t *seq2SoA = NULL;
    seq2SoA = (int16_t *)_mm_malloc(this->maxQerLen * SIMD_WIDTH16 * numThreads * sizeof(int16_t), 64);

    assert(seq1SoA != NULL);
    assert(seq2SoA != NULL);    
    
    int32_t ii;
    int32_t roundNumPairs = ((numPairs + SIMD_WIDTH16 - 1) / SIMD_WIDTH16 ) * SIMD_WIDTH16;
    // assert(roundNumPairs < BATCH_SIZE * SEEDS_PER_READ);
    for(ii = numPairs; ii < roundNumPairs; ii++)
    {
        pairArray[ii].regid = ii;
        pairArray[ii].id = ii;
        pairArray[ii].len1 = 0;
        pairArray[ii].len2 = 0;
    }

#if RDT
    st2 = __rdtsc();
#endif
    
#if SORT_PAIRS     // disbaled in bwa-mem2 (only used in separate benchmark sw code)
    {
        // Sort the sequences according to decreasing order of lengths
        SeqPair *tempArray = (SeqPair *)_mm_malloc(SORT_BLOCK_SIZE * numThreads *
                                                   sizeof(SeqPair), 64);
        int16_t *hist = (int16_t *)_mm_malloc((this->maxQerLen + 32) * numThreads *
                                              sizeof(int16_t), 64);
        
        #pragma omp parallel num_threads(numThreads)
        {
            int32_t tid = omp_get_thread_num();
            SeqPair *myTempArray = tempArray + tid * SORT_BLOCK_SIZE;
            int16_t *myHist = hist + tid * (this->maxQerLen + 32);
            
            #pragma omp for
            for(ii = 0; ii < roundNumPairs; ii+=SORT_BLOCK_SIZE)
            {
                int32_t first, last;
                first = ii;
                last  = ii + SORT_BLOCK_SIZE;
                if(last > roundNumPairs) last = roundNumPairs;
                sortPairsLen(pairArray + first, last - first, myTempArray, myHist);
            }
        }
        _mm_free(hist);
    }
#endif

#if RDT
    st3 = __rdtsc();
#endif
    
// #pragma omp parallel num_threads(numThreads)
    {
        int32_t i;
        // uint16_t tid = omp_get_thread_num();
        uint16_t tid = 0;  // no threading here.
        int16_t *mySeq1SoA = NULL;
        mySeq1SoA = seq1SoA + tid * this->maxRefLen * SIMD_WIDTH16;
        assert(mySeq1SoA != NULL);
            
        int16_t *mySeq2SoA = NULL;
        mySeq2SoA = seq2SoA + tid * this->maxQerLen * SIMD_WIDTH16;
        assert(mySeq1SoA != NULL);
        
        uint8_t *seq1;
        uint8_t *seq2;
        
        int nstart = 0, nend = numPairs;


//#pragma omp for schedule(dynamic, 128)
        for(i = nstart; i < nend; i+=SIMD_WIDTH16)
        {
            int32_t j, k;
            int maxLen1 = 0;
            int maxLen2 = 0;

            for(j = 0; j < SIMD_WIDTH16; j++)
            {
                SeqPair sp = pairArray[i + j];
#if MAINY
                seq1 = seqBufRef + (int64_t)sp.id * this->maxRefLen;
#else
                seq1 = seqBufRef + sp.idr;
#endif
                assert(sp.len1 < this->maxRefLen);
                for(k = 0; k < sp.len1; k++)
                {
                    mySeq1SoA[k * SIMD_WIDTH16 + j] = (seq1[k] == AMBIG_ ? AMBR16:seq1[k]);
                }
                if(maxLen1 < sp.len1) maxLen1 = sp.len1;
            }
            for(j = 0; j < SIMD_WIDTH16; j++)
            {
                SeqPair sp = pairArray[i + j];
                for(k = sp.len1; k <= maxLen1; k++) //removed "="
                {
                    // mySeq1SoA[k * SIMD_WIDTH16 + j] = DUMMY1_;
                    mySeq1SoA[k * SIMD_WIDTH16 + j] = 0xFFFF;
                }
            }

            for(j = 0; j < SIMD_WIDTH16; j++)
            {
                SeqPair sp = pairArray[i + j];
#if MAINY
                seq2 = seqBufQer + (int64_t)sp.id * this->maxQerLen;
#else
                seq2 = seqBufQer + sp.idq;
#endif
                assert(sp.len2 <= this->maxQerLen);


                // int quanta = 8 - sp.len2 % 8;  // based on SSE-16 bit lane
                int quanta = (sp.len2 + 8 - 1)/8;  // based on SSE-16 bit lane
                quanta *= 8;
                assert(sp.len2 < this->maxQerLen);
                for(k = 0; k < sp.len2; k++)
                {
                    mySeq2SoA[k * SIMD_WIDTH16 + j] = (seq2[k]==AMBIG_? AMBQ16:seq2[k]);
                }
                
                assert(quanta < this->maxQerLen); 
                for(k = sp.len2; k < quanta; k++) {
                    mySeq2SoA[k * SIMD_WIDTH16 + j] = DUMMY3;
                }
                if(maxLen2 < (quanta)) maxLen2 = quanta;
            }
            
            for(j = 0; j < SIMD_WIDTH16; j++)
            {
                SeqPair sp = pairArray[i + j];
                int quanta = (sp.len2 + 8 - 1)/8;  // based on SSE2-16 bit lane
                quanta *= 8;
                for(k = quanta; k <= maxLen2; k++)
                {
                    // mySeq2SoA[k * SIMD_WIDTH16 + j] = DUMMY2_;
                    mySeq2SoA[k * SIMD_WIDTH16 + j] = 0xFFFF;
                }
            }

            kswv512_16(mySeq1SoA, mySeq2SoA,
                       maxLen1, maxLen2,
                       pairArray + i,
                       aln, i,
                       tid,
                       numPairs,
                       phase);
        }
    }

#if RDT
    st4 = __rdtsc();
#endif
    
#if SORT_PAIRS     // disbaled in bwa-mem2 (only used in separate benchmark sw code)
    {
        // Sort the sequences according to increasing order of id
        #pragma omp parallel num_threads(numThreads)
        {
            int32_t tid = omp_get_thread_num();
            SeqPair *myTempArray = tempArray + tid * SORT_BLOCK_SIZE;
            
            #pragma omp for
            for(ii = 0; ii < roundNumPairs; ii+=SORT_BLOCK_SIZE)
            {
                int32_t first, last;
                first = ii;
                last  = ii + SORT_BLOCK_SIZE;
                if(last > roundNumPairs) last = roundNumPairs;
                sortPairsId(pairArray + first, first, last - first, myTempArray);
            }
        }
        _mm_free(tempArray);
    }
#endif

#if RDT
    st5 = __rdtsc();
    setupTicks = st2 - st1;
    sort1Ticks = st3 - st2;
    swTicks = st4 - st3;
    sort2Ticks = st5 - st4;
#endif
    // free mem
    _mm_free(seq1SoA);
    _mm_free(seq2SoA);
    return; 
}

int kswv::kswv512_16(int16_t seq1SoA[],
                     int16_t seq2SoA[],
                     int16_t nrow,
                     int16_t ncol,
                     SeqPair *p,
                     kswr_t *aln,
                     int po_ind,
                     uint16_t tid,
                     int32_t numPairs,
                     int phase)
{
    int m_b, n_b;
    int16_t minsc[SIMD_WIDTH16] = {0}, endsc[SIMD_WIDTH16] = {0};
    uint64_t *b;
    int limit = nrow;
    
    __m512i zero512 = _mm512_setzero_si512();
    __m512i one512  = _mm512_set1_epi16(1);
    __m512i minus1  = _mm512_set1_epi16(-1);

    int16_t temp[SIMD_WIDTH16] __attribute((aligned(64))) = {0};
    int16_t temp1[SIMD_WIDTH16] __attribute((aligned(64)));
    int16_t temp2[SIMD_WIDTH16] __attribute((aligned(64)));

    // query profile, we use xor operations on the strings
    // temp[0] = this->w_match;
    // temp[1] = temp[2] = temp[3] =  this->w_mismatch;
    // temp[4] = temp[5] = temp[6] = temp[7] =  this->w_mismatch;
    // temp[16] = temp[17] = temp[18] = temp[19] = this->w_ambig;
    // temp[31] =  this->w_ambig;
    // temp[12] = temp[13] = temp[14] = temp[15] =  this->w_ambig;
    // temp[10] = temp[20] = temp[30] = this->w_mismatch;
    
    temp[0] = this->w_match;    // matching
    temp[1]  = temp[2]  = temp[3]  =  this->w_mismatch;  // mis-matching    
    temp[12] = temp[13] = temp[14] = temp[15] =  this->w_ambig;
    temp[16] = temp[17] = temp[18] = temp[19] = this->w_ambig;
    temp[31] = this->w_ambig;
       
    __m512i perm512 = _mm512_load_si512(temp);
        
    m_b = n_b = 0; b = 0;

    __mmask32 minsc_msk_a = 0x0000, endsc_msk_a = 0x0000;
    int val = 0;
    for (int i=0; i<SIMD_WIDTH16; i++) {
        int xtra = p[i].h0;
        val = (xtra & KSW_XSUBO)? xtra & 0xffff : 0x10000;
        if (val <= SHRT_MAX) {
            minsc[i] = val;
            //if (val < SHRT_MIN)
            //  minsc[i] = -1;
            minsc_msk_a |= (0x1 << i);
        }
        // msc_mask;
        val = (xtra & KSW_XSTOP)? xtra & 0xffff : 0x10000;
        if (val <= SHRT_MAX) {
            endsc[i] = val;
            //if (val < SHRT_MIN)
            //  endsc[i] = -1;
            endsc_msk_a |= (0x1 << i);
        }
    }

    __m512i minsc512 = _mm512_load_si512((__m512i*) minsc);
    __m512i endsc512 = _mm512_load_si512((__m512i*) endsc);
    
    __m512i e_del512    = _mm512_set1_epi16(this->e_del);
    __m512i oe_del512   = _mm512_set1_epi16(this->o_del + this->e_del);
    __m512i e_ins512    = _mm512_set1_epi16(this->e_ins);
    __m512i oe_ins512   = _mm512_set1_epi16(this->o_ins + this->e_ins);
    __m512i gmax512     = zero512; // exit1 = zero512;
    // __m512i te512       = zero512;  // change to -1
    __m512i te512       = _mm512_set1_epi16(-1);
    __mmask32 exit0 = 0xFFFFFFFF;

    
    tid = 0;  // no threading here.
    int16_t *H0     = H16_0 + tid * SIMD_WIDTH16 * this->maxQerLen;
    int16_t *H1     = H16_1 + tid * SIMD_WIDTH16 * this->maxQerLen;
    int16_t *Hmax   = H16_max + tid * SIMD_WIDTH16 * this->maxQerLen;
    int16_t *F      = F16 + tid * SIMD_WIDTH16 * this->maxQerLen;
    int16_t *rowMax = rowMax16 + tid * SIMD_WIDTH16 * this->maxRefLen;
    
    _mm_prefetch((const char*) (F + SIMD_WIDTH16), _MM_HINT_NTA);
    _mm_prefetch((const char*) seq2SoA, _MM_HINT_NTA);
    _mm_prefetch((const char*) seq1SoA, _MM_HINT_NTA);
    _mm_prefetch((const char*) (H1 + SIMD_WIDTH16), _MM_HINT_NTA);
    _mm_prefetch((const char*) (F + SIMD_WIDTH16), _MM_HINT_NTA);

    for (int i=ncol; i >= 0; i--) {
        _mm512_store_si512((__m512*) (H0 + i * SIMD_WIDTH16), zero512);
        _mm512_store_si512((__m512*) (Hmax + i * SIMD_WIDTH16), zero512);
        _mm512_store_si512((__m512*) (F + i * SIMD_WIDTH16), zero512);
    }

    __m512i max512 = zero512, imax512, pimax512 = zero512;
    __mmask32 mask512 = 0x0000;
    __mmask32 minsc_msk = 0x0000;

    __m512i qe512 = _mm512_set1_epi16(0);
    _mm512_store_si512((__m512i *)(H0), zero512);
    _mm512_store_si512((__m512i *)(H1), zero512);
    __m512i i512 = zero512;
    
    int i;
    for (i=0; i < nrow; i++)
    {
        __m512i e11 = zero512;
        __m512i h00, h11, h10, s1;
        int j;
        
        s1 = _mm512_load_si512((__m512i *)(seq1SoA + (i + 0) * SIMD_WIDTH16));
        h10 = zero512;
        imax512 = zero512;
        __m512i iqe512 = _mm512_set1_epi16(-1);

        __m512i l512 = zero512;
        for (j=0; j<ncol; j++)
        {
            __m512i f11, s2, f21;
            h00 = _mm512_load_si512((__m512i *)(H0 + j * SIMD_WIDTH16));
            s2  = _mm512_load_si512((__m512i *)(seq2SoA + (j) * SIMD_WIDTH16));
            f11 = _mm512_load_si512((__m512i *)(F + (j+1) * SIMD_WIDTH16));
            
            MAIN_SAM_CODE16_OPT(s1, s2, h00, h11, e11, f11, f21, max512);

            _mm512_store_si512((__m512i *)(H1 + (j+1) * SIMD_WIDTH16), h11);
            _mm512_store_si512((__m512i *)(F + (j+1) * SIMD_WIDTH16), f21);
            l512 = _mm512_add_epi16(l512, one512);
            // prof[DP2][0] += 22;
            
        }   /* Inner DP loop */
        
        // Block I
        if (i > 0) {
            __mmask32 msk32 = _mm512_cmpgt_epi16_mask(imax512, pimax512);
            msk32 |= mask512;
            // msk32 &= ~minsc_msk;
            pimax512 = _mm512_mask_blend_epi16(msk32, pimax512, minus1);            
            pimax512 = _mm512_mask_blend_epi16(minsc_msk, minus1, pimax512);
            //pimax512 = _mm512_mask_mov_epi16(minus1, minsc_msk, pimax512);

            // new
            // pimax512 = _mm512_mask_blend_epi16(exit0, zero512, pimax512);
            pimax512 = _mm512_mask_blend_epi16(exit0, minus1, pimax512);
            
            _mm512_store_si512((__m512i *) (rowMax + (i-1)*SIMD_WIDTH16), pimax512);
            mask512 = ~msk32;
        }
        pimax512 = imax512;
        minsc_msk = _mm512_cmpge_epi16_mask(imax512, minsc512);
        minsc_msk &= minsc_msk_a;
        
        // Block II: gmax, te
        __mmask32 cmp0 = _mm512_cmpgt_epi16_mask(imax512, gmax512);
        cmp0 &= exit0;
        // gmax512 = _mm512_max_epi16(gmax512, imax512);
        gmax512 = _mm512_mask_blend_epi16(cmp0, gmax512, imax512);
        // te512 = _mm512_mask_mov_epi16(te512, cmp0, i512);
        te512 = _mm512_mask_blend_epi16(cmp0, te512,  i512);
        qe512 = _mm512_mask_blend_epi16(cmp0, qe512,  iqe512);
        
        cmp0 = _mm512_cmpge_epi16_mask(gmax512, endsc512);
        cmp0 &= endsc_msk_a;
        
        exit0 = (~cmp0) & exit0;
        if (exit0 == 0) {
            limit = i++;
            break;
        }
        
        int16_t *S = H1; H1 = H0; H0 = S;
        i512 = _mm512_add_epi16(i512, one512);
    } // for nrow
    
    pimax512 = _mm512_mask_blend_epi16(mask512, pimax512, minus1);
    pimax512 = _mm512_mask_blend_epi16(minsc_msk, minus1, pimax512);
    pimax512 = _mm512_mask_blend_epi16(exit0, minus1, pimax512);
    /* Guard on i > 0: all-len1==0 batches give nrow == 0 and would
     * otherwise underflow into the rowMax allocation. issue 38 / PR 289. */
    if (i > 0) {
        _mm512_store_si512((__m512i *) (rowMax + (i-1) * SIMD_WIDTH16), pimax512);
    }
    // __m512i max512_ = max512;

    /******************* DP loop over *****************************/
    /*************** Partial output setting ***************/
    int16_t score[SIMD_WIDTH16] __attribute((aligned(64)));
    int16_t te[SIMD_WIDTH16] __attribute((aligned(64)));
    int16_t qe[SIMD_WIDTH16] __attribute((aligned(64)));
    int16_t low[SIMD_WIDTH16] __attribute((aligned(64)));
    int16_t high[SIMD_WIDTH16] __attribute((aligned(64)));  
    _mm512_store_si512((__m512i *) score, gmax512); 
    _mm512_store_si512((__m512i *) te, te512);
    _mm512_store_si512((__m512i *) qe, qe512);

    for (int l=0; l<SIMD_WIDTH16 && (po_ind + l) < numPairs; l++) {
        int ind = po_ind + l;
#if !MAINY
        ind = p[l].regid;    // index of corr. aln
        if (phase) {
            if (aln[ind].score == score[l]) {
                aln[ind].tb = aln[ind].te - te[l];
                aln[ind].qb = aln[ind].qe - qe[l];
            }
        } else {
            aln[ind].score = score[l];
            aln[ind].te = te[l];
            aln[ind].qe = qe[l];
        }
#else
        aln[ind].score = score[l];
        aln[ind].te = te[l];
        aln[ind].qe = qe[l];
#endif
    }
    
#if !MAINY
    if (phase) return 1;
#endif

    /*************** Score2 and te2 — per-lane scalar emulation *******************
     *
     * Same b[]-consolidation fix as the 8-bit kernels (see kswv256_u8 for
     * full rationale). Additionally fixes three pre-existing bugs the
     * pre-fix 16-bit kernel carried but 8-bit did not:
     *
     *   (a) aggregate maxl/minh bounds — the same fix-3 issue PR #21
     *       addressed on the AVX-512BW 8-bit kernel; this 16-bit one
     *       never got it. A lane with a tight primary region would see
     *       rows inside its own [low, high] count as score2 candidates
     *       because the outer loop used the batch-wide max(low) /
     *       min(high) instead of per-lane bounds.
     *   (b) no minsc filter — scalar ksw_i16 records b[] entries only
     *       when imax >= minsc; pre-fix 16-bit had no such gate,
     *       leaking sub-minsc plateau scores into max2.
     *   (c) no qe mask — lanes without a valid primary (qe == 0) could
     *       contribute spurious rmax values since their rowMax was
     *       never zeroed.
     *
     * All three fall out of the per-lane scalar loop naturally. */
    int qmax = this->g_qmax;
    for (int j = 0; j < SIMD_WIDTH16; j++) {
        int val = (score[j] + qmax - 1) / qmax;
        low[j]  = te[j] - val;
        high[j] = te[j] + val;
    }

    /* rowMax entries produced: i == limit on normal completion, and
     * i == limit + 1 after early exit (limit = i++ stores old i). Using i
     * directly gives a single expression that covers both exit paths and
     * includes the terminating row stored at rowMax[limit]. */
    const int processed_rows = i;

    for (int l = 0; l < SIMD_WIDTH16 && (po_ind + l) < numPairs; l++) {
        int ind = p[l].regid;
        /* Match scalar ksw_i16: when KSW_XSUBO is absent, minsc = 0x10000
         * (> SHRT_MAX) so b[] never starts and score2/te2 stay at -1.
         * minsc_msk_a tracks lanes with KSW_XSUBO set; lanes not in the
         * mask have minsc[l] zero-initialized and must skip the scan.
         *
         * Inactive-lane check is te[l] < 0, not !qe[l]: unlike the 8-bit
         * kernels, 16-bit qe[l] is the real query-end coordinate and
         * qe == 0 is a valid endpoint (best cell at query column 0).
         * te[l] is initialized to -1 and only overwritten when a new
         * gmax is recorded, so te[l] < 0 cleanly marks lanes with no
         * primary alignment. */
        if (te[l] < 0 || !(minsc_msk_a & ((__mmask32)1 << l))) {
            aln[ind].score2 = -1;
            aln[ind].te2    = -1;
            continue;
        }

        int len1_l  = (int)p[l].len1;
        int low_l   = (int)low[l];
        int high_l  = (int)high[l];
        int minsc_l = (int)minsc[l];
        int score2  = -1;
        int te2     = -1;
        /* b_pos = -2 sentinel: (-2 + 1) != any real row, forces "append"
         * on the first qualifying row. */
        int b_score = -1;
        int b_pos   = -2;

        int nrows = processed_rows < len1_l ? processed_rows : len1_l;
        for (int i2 = 0; i2 < nrows; i2++) {
            int imax = (int)rowMax[i2 * SIMD_WIDTH16 + l];
            if (imax < minsc_l) continue;

            if (b_pos + 1 != i2) {
                /* APPEND: flush outgoing b[] entry to score2 first. */
                if (b_pos >= 0 &&
                    (b_pos < low_l || b_pos > high_l) &&
                    b_score > score2) {
                    score2 = b_score;
                    te2    = b_pos;
                }
                b_score = imax;
                b_pos   = i2;
            } else if (b_score < imax) {
                /* UPDATE: strict greater extends the run's anchor. */
                b_score = imax;
                b_pos   = i2;
            }
        }

        /* Flush trailing b[] entry. */
        if (b_pos >= 0 &&
            (b_pos < low_l || b_pos > high_l) &&
            b_score > score2) {
            score2 = b_score;
            te2    = b_pos;
        }

        aln[ind].score2 = score2;
        aln[ind].te2    = te2;
    }

    return 1;
}

#endif // AVX512BW



/**************************************Scalar code***************************************/
/* This is the original SW code from bwa-mem. We are keeping both, 8-bit and 16-bit 
   implementations, here for benchmarking purpose. 
   The interface call to the code is very simple and similar to the one we used above.
   By default the code is disabled.
 */

#if ORIGINAL_SW_SAM
void kswv::bwa_fill_scmat(int8_t mat[25]) {
    int a = this->w_match;
    int b = this->w_mismatch;
    int ambig = this->w_ambig;
    
    int i, j, k;
    for (i = k = 0; i < 4; ++i) {
        for (j = 0; j < 4; ++j)
            mat[k++] = i == j? a : b;
        mat[k++] = ambig; // ambiguous base
    }
    for (j = 0; j < 5; ++j) mat[k++] = ambig;
}

/**
 * Initialize the query data structure
 *
 * @param size   Number of bytes used to store a score; valid valures are 1 or 2
 * @param qlen   Length of the query sequence
 * @param query  Query sequence
 * @param m      Size of the alphabet
 * @param mat    Scoring matrix in a one-dimension array
 *
 * @return       Query data structure
 */
kswq_t* kswv::ksw_qinit(int size, int qlen, uint8_t *query, int m, const int8_t *mat)
{
    kswq_t *q;
    int slen, a, tmp, p;

    size = size > 1? 2 : 1;
    p = 8 * (3 - size); // # values per __m128i
    slen = (qlen + p - 1) / p; // segmented length
    q = (kswq_t*)malloc(sizeof(kswq_t) + 256 + 16 * slen * (m + 4)); // a single block of memory
    q->qp = (__m128i*)(((size_t)q + sizeof(kswq_t) + 15) >> 4 << 4); // align memory
    q->H0 = q->qp + slen * m;
    q->H1 = q->H0 + slen;
    q->E  = q->H1 + slen;
    q->Hmax = q->E + slen;
    q->slen = slen; q->qlen = qlen; q->size = size;
    // compute shift
    tmp = m * m;
    for (a = 0, q->shift = 127, q->mdiff = 0; a < tmp; ++a) { // find the minimum and maximum score
        if (mat[a] < (int8_t)q->shift) q->shift = mat[a];
        if (mat[a] > (int8_t)q->mdiff) q->mdiff = mat[a];
    }
    
    q->max = q->mdiff;
    q->shift = 256 - q->shift; // NB: q->shift is uint8_t
    q->mdiff += q->shift; // this is the difference between the min and max scores

    // An example: p=8, qlen=19, slen=3 and segmentation:
    //  {{0,3,6,9,12,15,18,-1},{1,4,7,10,13,16,-1,-1},{2,5,8,11,14,17,-1,-1}}
    if (size == 1) {
        int8_t *t = (int8_t*)q->qp;
        for (a = 0; a < m; ++a) {
            int i, k, nlen = slen * p;
            const int8_t *ma = mat + a * m;
            for (i = 0; i < slen; ++i)
                for (k = i; k < nlen; k += slen) // p iterations
                    *t++ = (k >= qlen? 0 : ma[query[k]]) + q->shift;
        }
    } else {
        int16_t *t = (int16_t*)q->qp;
        for (a = 0; a < m; ++a) {
            int i, k, nlen = slen * p;
            const int8_t *ma = mat + a * m;
            for (i = 0; i < slen; ++i)
                for (k = i; k < nlen; k += slen) // p iterations
                    *t++ = (k >= qlen? 0 : ma[query[k]]);
        }
    }
    return q;
}

kswr_t kswv::kswvScalar_u8(kswq_t *q, int tlen, const uint8_t *target,
                          int _o_del, int _e_del, int _o_ins, int _e_ins,
                          int xtra) // the first gap costs -(_o+_e)
{
    int slen, i, m_b, n_b, te = -1, gmax = 0, minsc, endsc;
    uint64_t *b;
    __m128i zero, oe_del, e_del, oe_ins, e_ins, shift, *H0, *H1, *E, *Hmax;
    kswr_t r;

#define __max_16(ret, xx) do { \
        (xx) = _mm_max_epu8((xx), _mm_srli_si128((xx), 8)); \
        (xx) = _mm_max_epu8((xx), _mm_srli_si128((xx), 4)); \
        (xx) = _mm_max_epu8((xx), _mm_srli_si128((xx), 2)); \
        (xx) = _mm_max_epu8((xx), _mm_srli_si128((xx), 1)); \
        (ret) = _mm_extract_epi16((xx), 0) & 0x00ff;        \
    } while (0)
    
    // initialization
    r = g_defr;
    minsc = (xtra & KSW_XSUBO)? xtra & 0xffff : 0x10000;
    endsc = (xtra & KSW_XSTOP)? xtra & 0xffff : 0x10000;
    m_b = n_b = 0; b = 0;   zero = _mm_set1_epi32(0);
    oe_del = _mm_set1_epi8(_o_del + _e_del);
    e_del = _mm_set1_epi8(_e_del);
    oe_ins = _mm_set1_epi8(_o_ins + _e_ins);
    e_ins = _mm_set1_epi8(_e_ins);
    shift = _mm_set1_epi8(q->shift);
    H0 = q->H0; H1 = q->H1; E = q->E; Hmax = q->Hmax;
    slen = q->slen;
    for (i = 0; i < slen; ++i) {
        _mm_store_si128(E + i, zero);
        _mm_store_si128(H0 + i, zero);
        _mm_store_si128(Hmax + i, zero);
    }
    // the core loop
    for (i = 0; i < tlen; ++i) {
        int j, k, cmp, imax;
        __m128i e, h, t, f = zero, max = zero, *S = q->qp + target[i] * slen; // s is the 1st score vector
        h = _mm_load_si128(H0 + slen - 1); // h={2,5,8,11,14,17,-1,-1} in the above example
        h = _mm_slli_si128(h, 1); // h=H(i-1,-1); << instead of >> because x64 is little-endian
        for (j = 0; LIKELY(j < slen); ++j) {
            /* SW cells are computed in the following order:
             *   H(i,j)   = max{H(i-1,j-1)+S(i,j), E(i,j), F(i,j)}
             *   E(i+1,j) = max{H(i,j)-q, E(i,j)-r}
             *   F(i,j+1) = max{H(i,j)-q, F(i,j)-r}
             */
            // compute H'(i,j); note that at the beginning, h=H'(i-1,j-1)
            h = _mm_adds_epu8(h, _mm_load_si128(S + j));
            h = _mm_subs_epu8(h, shift); // h=H'(i-1,j-1)+S(i,j)
            e = _mm_load_si128(E + j); // e=E'(i,j)
            h = _mm_max_epu8(h, e);
            h = _mm_max_epu8(h, f); // h=H'(i,j)
            max = _mm_max_epu8(max, h); // set max
            _mm_store_si128(H1 + j, h); // save to H'(i,j)
            // now compute E'(i+1,j)
            e = _mm_subs_epu8(e, e_del); // e=E'(i,j) - e_del
            t = _mm_subs_epu8(h, oe_del); // h=H'(i,j) - o_del - e_del
            e = _mm_max_epu8(e, t); // e=E'(i+1,j)
            _mm_store_si128(E + j, e); // save to E'(i+1,j)
            // now compute F'(i,j+1)
            f = _mm_subs_epu8(f, e_ins);
            t = _mm_subs_epu8(h, oe_ins); // h=H'(i,j) - o_ins - e_ins
            f = _mm_max_epu8(f, t);
            // get H'(i-1,j) and prepare for the next j
            h = _mm_load_si128(H0 + j); // h=H'(i-1,j)
        }
        // NB: we do not need to set E(i,j) as we disallow adjecent insertion and then deletion
        for (k = 0; LIKELY(k < 16); ++k) { // this block mimics SWPS3; NB: H(i,j) updated in the lazy-F loop cannot exceed max
            f = _mm_slli_si128(f, 1);
            for (j = 0; LIKELY(j < slen); ++j) {
                h = _mm_load_si128(H1 + j);
                h = _mm_max_epu8(h, f); // h=H'(i,j)
                _mm_store_si128(H1 + j, h);
                h = _mm_subs_epu8(h, oe_ins);
                f = _mm_subs_epu8(f, e_ins);
                cmp = _mm_movemask_epi8(_mm_cmpeq_epi8(_mm_subs_epu8(f, h), zero));
                if (UNLIKELY(cmp == 0xffff)) goto end_loop16;
            }
        }
end_loop16:
        //int k;for (k=0;k<16;++k)printf("%d ", ((uint8_t*)&max)[k]);printf("\n");
        __max_16(imax, max); // imax is the maximum number in max
        if (imax >= minsc) { // write the b array; this condition adds branching unfornately
            if (n_b == 0 || (int32_t)b[n_b-1] + 1 != i) { // then append
                if (n_b == m_b) {
                    m_b = m_b? m_b<<1 : 8;
                    b = (uint64_t*) realloc (b, 8 * m_b);
                }
                b[n_b++] = (uint64_t)imax<<32 | i;
            } else if ((int)(b[n_b-1]>>32) < imax) b[n_b-1] = (uint64_t)imax<<32 | i; // modify the last
        }
        if (imax > gmax) {
            gmax = imax; te = i; // te is the end position on the target
            for (j = 0; LIKELY(j < slen); ++j) // keep the H1 vector
                _mm_store_si128(Hmax + j, _mm_load_si128(H1 + j));
            if (gmax + q->shift >= 255 || gmax >= endsc) break;
        }
        S = H1; H1 = H0; H0 = S; // swap H0 and H1
    }
    r.score = gmax + q->shift < 255? gmax : 255;
    r.te = te;
    if (r.score != 255) { // get a->qe, the end of query match; find the 2nd best score
        int max = -1, tmp, low, high, qlen = slen * 16;
        uint8_t *t = (uint8_t*) Hmax;
        for (i = 0; i < qlen; ++i, ++t)
            if ((int)*t > max) max = *t, r.qe = i / 16 + i % 16 * slen;
            else if ((int)*t == max && (tmp = i / 16 + i % 16 * slen) < r.qe) r.qe = tmp; 
        if (b) {
            assert(q->max != 0);
            i = (r.score + q->max - 1) / q->max;
            low = te - i; high = te + i;
            for (i = 0; i < n_b; ++i) {
                int e = (int32_t)b[i];
                if ((e < low || e > high) && (int)(b[i]>>32) > r.score2)
                    r.score2 = b[i]>>32, r.te2 = e;
            }
        }
    }
    
#if MAXI
    fprintf(stderr, "score: %d, te: %d, qe: %d, score2: %d, te2: %d\n",
            r.score, r.te, r.qe, r.score2, r.te2);
#endif
    
    free(b);
    return r;
}

kswr_t kswv::kswvScalar_i16(kswq_t *q, int tlen, const uint8_t *target,
                            int _o_del, int _e_del, int _o_ins, int _e_ins,
                            int xtra) // the first gap costs -(_o+_e)
{
    int slen, i, m_b, n_b, te = -1, gmax = 0, minsc, endsc;
    uint64_t *b;
    __m128i zero, oe_del, e_del, oe_ins, e_ins, *H0, *H1, *E, *Hmax;
    kswr_t r;
#define SIMD16 8

#define __max_8(ret, xx) do { \
        (xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 8)); \
        (xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 4)); \
        (xx) = _mm_max_epi16((xx), _mm_srli_si128((xx), 2)); \
        (ret) = _mm_extract_epi16((xx), 0); \
    } while (0)

    // initialization
    r = g_defr;
    minsc = (xtra&KSW_XSUBO)? xtra&0xffff : 0x10000;
    endsc = (xtra&KSW_XSTOP)? xtra&0xffff : 0x10000;
    m_b = n_b = 0; b = 0;
    zero = _mm_set1_epi32(0);
    oe_del = _mm_set1_epi16(_o_del + _e_del);
    e_del = _mm_set1_epi16(_e_del);
    oe_ins = _mm_set1_epi16(_o_ins + _e_ins);
    e_ins = _mm_set1_epi16(_e_ins);
    H0 = q->H0; H1 = q->H1; E = q->E; Hmax = q->Hmax;
    slen = q->slen;
    for (i = 0; i < slen; ++i) {
        _mm_store_si128(E + i, zero);
        _mm_store_si128(H0 + i, zero);
        _mm_store_si128(Hmax + i, zero);
    }
    // the core loop
    for (i = 0; i < tlen; ++i) {
        int j, k, imax;
        __m128i e, t, h, f = zero, max = zero, *S = q->qp + target[i] * slen; // s is the 1st score vector
        h = _mm_load_si128(H0 + slen - 1); // h={2,5,8,11,14,17,-1,-1} in the above example
        h = _mm_slli_si128(h, 2);
        for (j = 0; LIKELY(j < slen); ++j) {

            h = _mm_adds_epi16(h, *S++);
            e = _mm_load_si128(E + j);
            h = _mm_max_epi16(h, e);
            h = _mm_max_epi16(h, f);
            max = _mm_max_epi16(max, h);
            _mm_store_si128(H1 + j, h);
            e = _mm_subs_epu16(e, e_del);
            t = _mm_subs_epu16(h, oe_del);
            e = _mm_max_epi16(e, t);
            _mm_store_si128(E + j, e);
            f = _mm_subs_epu16(f, e_ins);
            t = _mm_subs_epu16(h, oe_ins);
            f = _mm_max_epi16(f, t);
            h = _mm_load_si128(H0 + j);
        }
        
        for (k = 0; LIKELY(k < 16); ++k) {
            f = _mm_slli_si128(f, 2);
            for (j = 0; LIKELY(j < slen); ++j) {
                h = _mm_load_si128(H1 + j);
                h = _mm_max_epi16(h, f);
                _mm_store_si128(H1 + j, h);
                h = _mm_subs_epu16(h, oe_ins);
                f = _mm_subs_epu16(f, e_ins);
                if(UNLIKELY(!_mm_movemask_epi8(_mm_cmpgt_epi16(f, h)))) goto end_loop8;
            }
        }
        
end_loop8:
        __max_8(imax, max);
        if (imax >= minsc) {
            if (n_b == 0 || (int32_t)b[n_b-1] + 1 != i) {
                if (n_b == m_b) {
                    m_b = m_b? m_b<<1 : 8;
                    b = (uint64_t*)realloc(b, 8 * m_b);
                }
                b[n_b++] = (uint64_t)imax<<32 | i;
            } else if ((int)(b[n_b-1]>>32) < imax) b[n_b-1] = (uint64_t)imax<<32 | i; // modify the last
        }
        if (imax > gmax) {
            gmax = imax; te = i;
            for (j = 0; LIKELY(j < slen); ++j)
                _mm_store_si128(Hmax + j, _mm_load_si128(H1 + j));
            if (gmax >= endsc) break;
        }
        S = H1; H1 = H0; H0 = S;
    }
    
    r.score = gmax; r.te = te;
    {
        int max = -1, tmp, low, high, qlen = slen * 8;
        uint16_t *t = (uint16_t*)Hmax;
        for (i = 0, r.qe = -1; i < qlen; ++i, ++t)
            if ((int)*t > max) max = *t, r.qe = i / 8 + i % 8 * slen;
            else if ((int)*t == max && (tmp = i / 8 + i % 8 * slen) < r.qe) r.qe = tmp; 
        if (b) {
            assert(q->max != 0);
            i = (r.score + q->max - 1) / q->max;
            low = te - i; high = te + i;
            for (i = 0; i < n_b; ++i) {
                int e = (int32_t)b[i];
                if ((e < low || e > high) && (int)(b[i]>>32) > r.score2)
                    r.score2 = b[i]>>32, r.te2 = e;
            }
        }
    }

    free(b);
    return r;
}

// -------------------------------------------------------------
// kswc scalar, wrapper function, the interface.
//-------------------------------------------------------------
void kswv::kswvScalarWrapper(SeqPair *seqPairArray,
                             uint8_t *seqBufRef,
                             uint8_t *seqBufQer,
                             kswr_t* aln,
                             int numPairs,
                             int nthreads,
                             bool sw, int tid) {
    
    int8_t mat[25];
    bwa_fill_scmat(mat);

    int st = 0, ed = numPairs;    
    for (int i = st; i < ed; i++)
    {
        SeqPair *p = seqPairArray + i;
        kswr_t *myaln = aln + p->regid;
            
        uint8_t *target = seqBufRef + p->idr;
        uint8_t *query = seqBufQer + p->idq;

        int tlen = p->len1;
        int qlen = p->len2;
        int xtra = p->h0;

        kswq_t *q;
        kswr_t ks;
            
        int vw = (xtra & KSW_XBYTE)? 1 : 2;
        if (sw == 0) {
            assert(vw == 1);
            q = ksw_qinit((xtra & KSW_XBYTE)? 1 : 2, qlen, query, this->m, mat);
            ks = kswvScalar_u8(q, tlen, target,
                               o_del, e_del,
                               o_ins, e_ins,
                               xtra, tid, i);
            free(q);
        } else {
            assert(vw == 2);
            q = ksw_qinit(2, qlen, query, this->m, mat);
            ks = kswvScalar_i16(q, tlen, target,
                                o_del, e_del,
                                o_ins, e_ins,
                                xtra);
            free(q);
        }

        myaln->score = ks.score;
        myaln->tb = ks.tb;
        myaln->te = ks.te;
        myaln->qb = ks.qb;
        myaln->qe = ks.qe;
        myaln->score2 = ks.score2;
        myaln->te2 = ks.te2;
    }

}
#endif

/************************** standalone benchmark code *************************************/
/** This is a standalone code for the benchmarking purpose 
    By default it is disabled.
**/

#if MAINY
#define DEFAULT_MATCH 1
#define DEFAULT_MISMATCH -4
#define DEFAULT_OPEN 6
#define DEFAULT_EXTEND 1
#define DEFAULT_AMBIG -1

// #define MAX_NUM_PAIRS 1000
// #define MATRIX_MIN_CUTOFF -100000000
// #define LOW_INIT_VALUE (INT32_MIN/2)
// #define AMBIG 52

int spot = 42419;
double freq = 2.3*1e9;
int32_t w_match, w_mismatch, w_open, w_extend, w_ambig;
uint64_t SW_cells;
char *pairFileName;
FILE *pairFile;
int8_t h0 = 0;
double clock_freq;
//uint64_t prof[10][112], data, SW_cells2;


void parseCmdLine(int argc, char *argv[])
{
    int i;
    w_match = DEFAULT_MATCH;
    w_mismatch = DEFAULT_MISMATCH;
    w_open = DEFAULT_OPEN;
    w_extend = DEFAULT_EXTEND;
    w_ambig = DEFAULT_AMBIG;
    
    int pairFlag = 0;
    for(i = 1; i < argc; i+=2)
    {
        if(strcmp(argv[i], "-match") == 0)
        {
            w_match = atoi(argv[i + 1]);
        }
        if(strcmp(argv[i], "-mismatch") == 0) //penalty, +ve number
        {
            w_mismatch = atoi(argv[i + 1]);
        }
        if(strcmp(argv[i], "-ambig") == 0)
        {
            w_ambig = atoi(argv[i + 1]);
        }

        if(strcmp(argv[i], "-gapo") == 0)
        {
            w_open = atoi(argv[i + 1]);
        }
        if(strcmp(argv[i], "-gape") == 0)
        {
            w_extend = atoi(argv[i + 1]);
        }
        if(strcmp(argv[i], "-pairs") == 0)
        {
            pairFileName = argv[i + 1];
            pairFlag = 1;
        }
        if(strcmp(argv[i], "-h0") == 0)
        {
            h0 = atoi(argv[i + 1]);
        }
    }
    if(pairFlag == 0)
    {
        printf("ERROR! pairFileName not specified.\n");
        exit(EXIT_FAILURE);
    }
}

int loadPairs(SeqPair *seqPairArray, uint8_t *seqBufRef, uint8_t* seqBufQer, FILE *pairFile)
{
    static int32_t cnt = 0;
    int32_t numPairs = 0;
    while(numPairs < MAX_NUM_PAIRS_ALLOC)
    {
        int32_t xtra = 0;
        char temp[10];
        fgets(temp, 10, pairFile);
        sscanf(temp, "%x", &xtra);
        //printf("xtra: %d, %x, %s\n", xtra, xtra, temp);

        //if(!fgets((char *)(seqBuf + numPairs * 2 * MAX_SEQ_LEN), MAX_SEQ_LEN, pairFile))
        if(!fgets((char *)(seqBufRef + numPairs * this->maxRefLen), this->maxRefLen, pairFile))
        {
            break;
        }
        //if(!fgets((char *)(seqBuf + (numPairs * 2 + 1) * MAX_SEQ_LEN), MAX_SEQ_LEN, pairFile))
        if(!fgets((char *)(seqBufQer + numPairs * this->maxQerLen), this->maxQerLen, pairFile)) 
        {
            printf("ERROR! Odd number of sequences in %s\n", pairFileName);
            break;
        }

        SeqPair sp;
        sp.id = numPairs;
        // sp.seq1 = seqBuf + numPairs * 2 * MAX_SEQ_LEN;
        // sp.seq2 = seqBuf + (numPairs * 2 + 1) * MAX_SEQ_LEN;
        sp.len1 = strnlen((char *)(seqBufRef + numPairs * this->maxRefLen), this->maxRefLen) - 1;
        sp.len2 = strnlen((char *)(seqBufQer + numPairs * this->maxQerLen), this->maxQerLen) - 1;
        sp.h0 = xtra;
        // sp.score = 0;

        //uint8_t *seq1 = seqBuf + numPairs * 2 * MAX_SEQ_LEN;
        //uint8_t *seq2 = seqBuf + (numPairs * 2 + 1) * MAX_SEQ_LEN;
        uint8_t *seq1 = seqBufRef + numPairs * this->maxRefLen;
        uint8_t *seq2 = seqBufQer + numPairs * this->maxQerLen;

        for (int l=0; l<sp.len1; l++)
            seq1[l] -= 48;
        for (int l=0; l<sp.len2; l++)
            seq2[l] -= 48;

        seqPairArray[numPairs] = sp;
        numPairs++;
        // SW_cells += (sp.len1 * sp.len2);
    }
    // fclose(pairFile);
    return numPairs;
}

// profiling stats
uint64_t find_stats(uint64_t *val, int nt, double &min, double &max, double &avg)
{
    min = 1e10;
    max = 0;
    avg = 0;
    for (int i=0; i<nt; i++) {
        avg += val[i];
        if (max < val[i]) max = val[i];
        if (min > val[i]) min = val[i];     
    }
    avg /= nt;
}

int main(int argc, char *argv[])
{

#ifdef VTUNE_ANALYSIS
    printf("Vtune analysis enabled....\n");
    __itt_pause();
#endif
    FILE *fsam = fopen("fsam.txt", "w");    
    kswv *mate;

    parseCmdLine(argc, argv);

    SeqPair *seqPairArray = (SeqPair *)_mm_malloc((MAX_NUM_PAIRS + SIMD_WIDTH8) * sizeof(SeqPair), 64);
    uint8_t *seqBufRef = NULL, *seqBufQer = NULL;
    seqBufRef = (uint8_t *)_mm_malloc((this->maxRefLen * MAX_NUM_PAIRS + MAX_LINE_LEN)
                                      * sizeof(int8_t), 64);
    seqBufQer = (uint8_t *)_mm_malloc((this->maxQerLen * MAX_NUM_PAIRS + MAX_LINE_LEN)
                                      * sizeof(int8_t), 64);

    kswr_t *aln = NULL;
    aln = (kswr_t *) _mm_malloc ((MAX_NUM_PAIRS + SIMD_WIDTH8) * sizeof(kswr_t), 64);
    
    if (seqBufRef == NULL || seqBufQer == NULL || aln == NULL)
    {
        printf("Memory not allocated\nExiting...\n");
        exit(EXIT_FAILURE);
    } else {
        printf("Memory allocated: %0.2lf MB\n",
               ((int64_t)(this->maxRefLen + this->maxQerLen ) * MAX_NUM_PAIRS + MAX_LINE_LEN)/1e6);
    }
    uint64_t tim = __rdtsc(), readTim = 0;
    
    int32_t numThreads = 1;
#pragma omp parallel
    {
        int32_t tid = omp_get_thread_num();
        int32_t nt = omp_get_num_threads();
        if(tid == (nt - 1))
        {
            numThreads = nt;
        }
    }
    numThreads =1 ;
    //printf("Done reading input file!!, numPairs: %d, nt: %d\n",
    //     numPairs, numThreads);

#if SORT_PAIRS     // disbaled in bwa-mem2 (only used in separate benchmark sw code)
    printf("\tSorting is enabled !!\n");
#endif

#if MAXI
    printf("\tResults printing (on console) is enabled.....\n");
#endif

    tim = __rdtsc();
    sleep(1);
    freq = __rdtsc() - tim;

    int numPairs = 0, totNumPairs = 0;
    FILE *pairFile = fopen(pairFileName, "r");
    if(pairFile == NULL)
    {
        fprintf(stderr, "Could not open file: %s\n", pairFileName);
        exit(EXIT_FAILURE);
    }

    // BandedPairWiseSW *pwsw = new BandedPairWiseSW(w_match, w_mismatch, w_open,
    //                                            w_extend, w_ambig, end_bonus);
    kswv *pwsw = new kswv(w_open, w_extend, w_open, w_extend, w_match, w_mismatch, numThreads);
    

    int64_t myTicks = 0;
    printf("Processor freq: %0.2lf MHz\n", freq/1e6);
    printf("Executing int8 code....\n");
#if VEC
    printf("Executing AVX512 vectorized code!!\n");
    
    while (1) {
        uint64_t tim = __rdtsc();

        //printf("Loading current batch of pairs..\n");
        numPairs = loadPairs(seqPairArray, seqBufRef, seqBufQer, pairFile);
        totNumPairs += numPairs;
#if STAT
        printf("Loading done, numPairs: %d\n", numPairs);      
        if (totNumPairs > spot) spot -= totNumPairs - numPairs;
        else continue;
#endif
        readTim += __rdtsc() - tim;
        if (numPairs == 0) break;

        tim = __rdtsc();
        int phase = 0;
        pwsw->getScores8(seqPairArray, seqBufRef, seqBufQer, aln, numPairs, numThreads, phase);
        // pwsw->getScores16(seqPairArray, seqBufRef, seqBufQer, aln, numPairs, numThreads, phase);
        myTicks += __rdtsc() - tim;

#if STAT
        for (int l=0; l<10; l++)
        {
            // SeqPair r = seqPairArray[l];
            kswr_t r = aln[l];
            fprintf(stderr, "%d %d %d %d %d %d %d\n", r.score, r.tb, r.te, r.qb, r.qe, r.score2, r.te2);
        }       
        break;
#endif
        
#if MAXI
        // printf("Execution complete!!, writing output!\n");
        //for (int l=0; l<numPairs; l++)
        //{
        //  // SeqPair r = seqPairArray[l];
        //  kswr_t r = aln[l];
        //  fprintf(stderr, "%d %d %d %d %d %d %d\n", r.score, r.tb, r.te, r.qb, r.qe, r.score2, r.te2);
        //}
        // printf("Vector code: Writing output completed!!!\n\n");
        // printf("Vector code -- Wrote output to the file\n");
#endif
    }
    
    // int64_t myTicks = 0;//pwsw->getTicks();
    printf("Read time  = %0.2lf\n", readTim/freq);
    printf("Overall SW cycles = %ld, %0.2lf\n", myTicks, myTicks*1.0/freq);
    printf("SW cells(T)  = %ld\n", SW_cells);
    printf("SW cells(||)  = %ld\n", SW_cells2);
    printf("SW GCUPS  = %lf\n", SW_cells * freq / myTicks);
    
    {
        printf("More stats:\n");
        // double freq = 2.3*1e9;
        double min, max, avg;
        find_stats(prof[1], numThreads, min, max, avg);
        printf("Time in pre-processing: %0.2lf (%0.2lf, %0.2lf)\n",
               avg*1.0/freq, min*1.0/freq, max*1.0/freq);
        find_stats(prof[0], numThreads, min, max, avg);
        printf("Time spent in smithWaterman(): %0.2lf (%0.2lf, %0.2lf)\n",
               avg*1.0/freq, min*1.0/freq, max*1.0/freq);
        printf("\nTotal cells computed: %ld, %0.2lf\n",
               prof[2][0], prof[2][0]*1.0/(totNumPairs));
        printf("\tTotal useful cells computed: %ld, %0.2lf\n",
               prof[3][0], prof[3][0]*1.0/(totNumPairs));
        // printf("Total bp read from memory: %ld\n", data);
        printf("Computations after exit: %ld\n", prof[4][0]);
        printf("Cumulative seq1/2 len: %ld, %ld\n", prof[4][1], prof[4][2]);
    
#if 1
        printf("\nDebugging info:\n");
        printf("Time taken for DP loop: %0.2lf\n", prof[DP][0]*1.0/freq);
        printf("Time taken for DP loop upper part: %0.2lf\n", prof[DP3][0]*1.0/freq);   
        printf("Time taken for DP inner loop: %0.2lf\n", prof[DP1][0]*1.0/freq);
        printf("Time taken for DP loop lower part: %ld\n", prof[DP2][0]);
#endif
    }

    
#else
    printf("Executing scalar code!!\n");
    while (1) {
        uint64_t tim = __rdtsc();
        numPairs = loadPairs(seqPairArray, seqBufRef, seqBufQer, pairFile);
        totNumPairs += numPairs;
#if STAT
        printf("Loading done, numPairs: %d\n", numPairs);
        if (totNumPairs > spot) spot -= totNumPairs - numPairs;
        else continue;
#endif
        readTim += __rdtsc() - tim;
        if (numPairs == 0) break;

        tim = __rdtsc();
        pwsw->kswvScalarWrapper(seqPairArray,
                               seqBufRef,
                               seqBufQer,
                               aln,
                               numPairs,
                               numThreads);
        myTicks += __rdtsc() - tim;
#if STAT
        for (int l=0; l<10; l++)
        {
            // SeqPair r = seqPairArray[l];
            kswr_t r = aln[l];
            fprintf(stderr, "%d %d %d %d %d %d %d\n", r.score, r.tb, r.te, r.qb, r.qe, r.score2, r.te2);
        }       
        break;
#endif
        
    } // while
    
    printf("Read time  = %0.2lf\n", readTim/freq);
    printf("Overall SW cycles = %ld, %0.2lf\n", myTicks, myTicks*1.0/freq);
    printf("Time taken for DP loop lower part: %ld \n", prof[DP2][0]);

#if 1
        printf("\nDebugging info:\n");
        printf("Time taken for DP loop: %0.2lf\n", prof[DP][0]*1.0/freq);
        printf("Time taken for DP loop upper part: %0.2lf\n", prof[DP3][0]*1.0/freq);   
        printf("Time taken for DP inner loop: %0.2lf\n", prof[DP1][0]*1.0/freq);
        printf("Time taken for DP loop lower part: %ld\n", prof[DP2][0]);
#endif

#endif


#ifdef VTUNE_ANALYSIS
    printf("Vtune analysis enabled....\n");
#endif
        
    // free memory
    _mm_free(seqPairArray);
    _mm_free(seqBufRef);
    _mm_free(seqBufQer);
    _mm_free(aln);
    
    fclose(pairFile);
    fclose(fsam);
    return 1;
}
#endif  // MAINY
