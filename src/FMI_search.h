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

Authors: Sanchit Misra <sanchit.misra@intel.com>; Vasimuddin Md <vasimuddin.md@intel.com>;
*****************************************************************************************/

#ifndef _FMI_SEARCH_H
#define _FMI_SEARCH_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>

/* SIMD compatibility for ARM/x86 */
#if defined(__ARM_NEON) || defined(__aarch64__) || defined(APPLE_SILICON)
    #include "simd_compat.h"
#else
    #include <immintrin.h>
#endif
#include <fstream>

#include "read_index_ele.h"
#include "bwa.h"

#define DUMMY_CHAR 6

#define assert_not_null(x, size, cur_alloc) \
        if (x == NULL) { fprintf(stderr, "Allocation of %0.2lf GB for " #x " failed.\nCurrent Allocation = %0.2lf GB\n", size * 1.0 /(1024*1024*1024), cur_alloc * 1.0 /(1024*1024*1024)); exit(EXIT_FAILURE); }

#define CP_BLOCK_SIZE 64
#define CP_FILENAME_SUFFIX ".bwt.2bit.64"
#define CP_MASK 63
#define CP_SHIFT 6

typedef struct checkpoint_occ_scalar
{
    int64_t cp_count[4];
    uint64_t one_hot_bwt_str[4];
}CP_OCC;

#if defined(__clang__) || defined(__GNUC__)
static inline int _mm_countbits_64(unsigned long x) {
    return __builtin_popcountl(x);
}
#endif

#define \
GET_OCC(pp, c, occ_id_pp, y_pp, occ_pp, one_hot_bwt_str_c_pp, match_mask_pp) \
                int64_t occ_id_pp = pp >> CP_SHIFT; \
                int64_t y_pp = pp & CP_MASK; \
                int64_t occ_pp = cp_occ[occ_id_pp].cp_count[c]; \
                uint64_t one_hot_bwt_str_c_pp = cp_occ[occ_id_pp].one_hot_bwt_str[c]; \
                uint64_t match_mask_pp = one_hot_bwt_str_c_pp & one_hot_mask_array[y_pp]; \
                occ_pp += _mm_countbits_64(match_mask_pp);

typedef struct smem_struct
{
#ifdef DEBUG
    uint64_t info; // for debug
#endif
    uint32_t rid;
    uint32_t m, n;
    int64_t k, l, s;
}SMEM;

#define SAL_PFD 16

#ifndef SMEM_LOCKSTEP_N
#define SMEM_LOCKSTEP_N 16
#endif

/* Lockstep depth for the third-pass (bwtSeedStrategy) re-seeding, tuned
 * separately from the phase-2 SMEM depth above. The third-pass lockstep is
 * gated to arm64 at its call site (bwamem.cpp) because it only wins on non-SMT
 * cores; N=8 is the measured whole-aligner optimum on Graviton4 (-1.7% vs the
 * scalar third pass, -16.7% on the seeding stage). Set to 1 to fall back to the
 * scalar third-pass path even on arm. */
#ifndef BWTSEED_LOCKSTEP_N
#define BWTSEED_LOCKSTEP_N 8
#endif

class FMI_search: public indexEle
{
    public:
    FMI_search(const char *fname);
    ~FMI_search();
    //int64_t beCalls;

    /* Read-only size accessors. Used by load_index_from_shm for section
     * validation, and by the round-trip test for byte-equality checks. */
    int64_t cp_occ_size_bytes()      const;
    int64_t sa_sample_count()        const;
    int64_t sa_ms_byte_size_bytes()  const { return sa_sample_count() * (int64_t)sizeof(int8_t); }
    int64_t sa_ls_word_size_bytes() const { return sa_sample_count() * (int64_t)sizeof(uint32_t); }

    /* Read-only data accessors. Used by the round-trip test to byte-compare
     * the packed segment payload against the in-memory loader's buffers. */
    const void     *cp_occ_data()     const { return cp_occ; }
    const int8_t   *sa_ms_byte_data() const { return sa_ms_byte; }
    const uint32_t *sa_ls_word_data() const { return sa_ls_word; }
    const int64_t  *count_data()      const { return count; }

    /* Return the shm segment base if load_index attached from shm, else NULL.
     * fastmap reuses this so it doesn't have to re-attach for the ref string. */
    uint8_t *shm_attached_base()     const { return shm_base; }

    /* emit_unpacked_ref defaults false: skip writing the unpacked `<prefix>.0123`.
     * `mem` pac-fetches the original reference from `.pac`, so `.0123` is never
     * read; pass true only to emit it for an external consumer (e.g. bwa-mem2). */
    int build_index(bool emit_unpacked_ref = false);
    /* load_pac=false skips loading the 2-bit packed reference (BNS only). D3
     * --meth uses this for the SEED index: seeding needs the FM-index + bns
     * (for the seed->original remap) but never the seed pac — extension/scoring
     * runs against the ORIGINAL pac (meth_orig_pac). Saves ~1.6 GB on hg38. */
    void load_index(bool load_pac = true);

    /* Attach to a packed bwa-mem3 index segment from bwa_shm_attach. Sets
     * scalars and the cp_occ / sa_ms_byte / sa_ls_word pointers; the
     * destructor munmaps `base` and leaves the aliased buffers untouched. */
    void load_index_from_shm(uint8_t *base, size_t len);

    /* matchArray sizing contract (applies to all four SMEM-emitting methods
     * below). The previous internal `max_smem` capacity guard was removed;
     * the caller MUST pre-size matchArray to hold at least:
     *
     *   - getSMEMsOnePosOneThread, getSMEMsOnePosOneThread_lockstep,
     *     getSMEMsAllPosOneThread:   numReads * max_readlength SMEMs,
     *     where max_readlength is the function parameter passed in.
     *
     *   - bwtSeedStrategyAllPosOneThread:   numReads * max_seq_length SMEMs,
     *     where max_seq_length = max_i(seq_[i].l_seq), computed by the
     *     caller (this method takes no max_readlength parameter).
     *
     * Writing past the pre-sized capacity is undefined behavior. The actual
     * number of SMEMs written is reported via *__numTotalSmem (or the
     * int64_t return value on bwtSeedStrategyAllPosOneThread). */
    /* matchArray must hold at least numReads * max_readlength SMEMs (caller-sized). */
    void getSMEMsOnePosOneThread(uint8_t *enc_qdb,
                                 int32_t *query_pos_array,
                                 int32_t *min_intv_array,
                                 int32_t *rid_array,
                                 int32_t numReads,
                                 int32_t batch_size,
                                 const bseq1_t *seq_,
                                 int32_t *query_cum_len_ar,
                                 int32_t  max_readlength,
                                 int32_t minSeedLen,
                                 SMEM *matchArray,
                                 int64_t *__numTotalSmem);

    /* Lockstep-batched variant of getSMEMsOnePosOneThread: advances
     * SMEM_LOCKSTEP_N reads' SMEM walks in slot-interleaved order to expose
     * N independent backwardExt dependency chains to the CPU's out-of-order
     * engine. Per-read algorithm is byte-identical to the scalar path; only
     * the cross-read interleaving is new. Output in matchArray is written
     * in the same (read, smem) order as the scalar path via per-slot match
     * buffers flushed by input-index cursor.
     * matchArray must hold at least numReads * max_readlength SMEMs (caller-sized). */
    void getSMEMsOnePosOneThread_lockstep(uint8_t *enc_qdb,
                                          int32_t *query_pos_array,
                                          int32_t *min_intv_array,
                                          int32_t *rid_array,
                                          int32_t numReads,
                                          int32_t batch_size,
                                          const bseq1_t *seq_,
                                          int32_t *query_cum_len_ar,
                                          int32_t  max_readlength,
                                          int32_t minSeedLen,
                                          SMEM *matchArray,
                                          int64_t *__numTotalSmem);

    /* matchArray must hold at least numReads * max_readlength SMEMs (caller-sized). */
    void getSMEMsAllPosOneThread(uint8_t *enc_qdb,
                                 int32_t *min_intv_array,
                                 int32_t *rid_array,
                                 int32_t numReads,
                                 int32_t batch_size,
                                 const bseq1_t *seq_,
                                 int32_t *query_cum_len_ar,
                                 int32_t max_readlength,
                                 int32_t minSeedLen,
                                 SMEM *matchArray,
                                 int64_t *__numTotalSmem);


    /* matchArray must hold at least numReads * (longest l_seq in seq_) SMEMs;
     * returns the actual count written. */
    int64_t bwtSeedStrategyAllPosOneThread(uint8_t *enc_qdb,
                                           int32_t *max_intv_array,
                                           int32_t numReads,
                                           const bseq1_t *seq_,
                                           int32_t *query_cum_len_ar,
                                           int32_t minSeedLen,
                                           SMEM *matchArray);

    /* Lockstep batched variant of bwtSeedStrategyAllPosOneThread. Processes
     * BWTSEED_LOCKSTEP_N reads' forward-extension walks in slot-interleaved
     * order so the CP_OCC cache-line misses from independent reads issue
     * concurrently. Same SMEM emission order as the scalar (reads in input
     * order; within a read, outer-x ascending). matchArray must hold at
     * least numReads * max_readlength SMEMs. */
    int64_t bwtSeedStrategyAllPosOneThread_lockstep(uint8_t *enc_qdb,
                                                    int32_t *max_intv_array,
                                                    int32_t numReads,
                                                    const bseq1_t *seq_,
                                                    int32_t *query_cum_len_ar,
                                                    int32_t minSeedLen,
                                                    SMEM *matchArray,
                                                    int32_t max_readlength);
        
    void sortSMEMs(SMEM *matchArray,
                   int64_t numTotalSmem[],
                   int32_t numReads,
                   int32_t readlength,
                   int nthreads);
    int64_t get_sa_entry(int64_t pos);
    void get_sa_entries(int64_t *posArray,
                        int64_t *coordArray,
                        uint32_t count,
                        int32_t nthreads);
    void get_sa_entries(SMEM *smemArray,
                        int64_t *coordArray,
                        int32_t *coordCountArray,
                        uint32_t count,
                        int32_t max_occ);
    int64_t get_sa_entry_compressed(int64_t pos, int tid=0);
    void get_sa_entries(SMEM *smemArray,
                        int64_t *coordArray,
                        int32_t *coordCountArray,
                        uint32_t count,
                        int32_t max_occ,
                        int tid);
    int64_t call_one_step(int64_t pos, int64_t &sa_entry, int64_t &offset);
    void get_sa_entries_prefetch(SMEM *smemArray, int64_t *coordArray,
                                 int64_t *coordCountArray, int64_t count,
                                 const int32_t max_occ, int tid, int64_t &id_);
    
    int64_t reference_seq_len;
    int64_t sentinel_index;
private:
        char file_name[PATH_MAX];
        int64_t count[5];
        uint32_t *sa_ls_word;
        int8_t *sa_ms_byte;
        CP_OCC *cp_occ;

        uint64_t *one_hot_mask_array;

        /* If non-NULL, cp_occ / sa_ms_byte / sa_ls_word point into a shared
         * memory mapping owned by the shm segment rather than being
         * _mm_malloc'd. shm_len is the byte length of the mapping so the
         * destructor can munmap it; the destructor also skips _mm_free on
         * shm-backed buffers. */
        uint8_t *shm_base;
        size_t   shm_len;

        /* Defined inline so all hot callers (getSMEMs* and ls_advance_*)
         * fully absorb the body and pay no struct-by-value pass / return
         * cost. backwardExt is called ~10^9 times per 5M-pair WGS run, so
         * the SysV-ABI hidden-pointer dance for the 24-byte SMEM struct
         * dominates self-time on gcc 12+ (issue #87): a single output store
         * `vmovdqu %ymm0, (%r8)` writing the return SMEM hits 42% of the
         * function's CPU samples on c7a (Zen 4) with gcc-14, with the
         * matching argument load close behind. always_inline removes the
         * call boundary and recovers the full gcc-11 baseline (and beats
         * it by ~4% wall-clock on c7a wgs-5M shm-warmed). */
        __attribute__((always_inline)) inline
        SMEM backwardExt(SMEM smem, uint8_t a) const
        {
            /* sp/ep and their checkpoint block + one-hot mask do not depend on
             * the base b, so hoist them out of the per-base occ computation.
             * Only k[a] and s[a] are consumed for the result, but all four s[b]
             * feed the l cumulation below, so the full 4-lane occ is computed. */
            const int64_t sp = (int64_t)smem.k;
            const int64_t ep = (int64_t)smem.k + (int64_t)smem.s;
            const CP_OCC &blk_sp = cp_occ[sp >> CP_SHIFT];
            const CP_OCC &blk_ep = cp_occ[ep >> CP_SHIFT];
            const uint64_t mask_sp = one_hot_mask_array[sp & CP_MASK];
            const uint64_t mask_ep = one_hot_mask_array[ep & CP_MASK];

            int64_t occ_sp[4], s[4], l[4];

#if defined(__ARM_NEON) || defined(__aarch64__) || defined(APPLE_SILICON)
            /* arm64 has no scalar popcount instruction, so the scalar GET_OCC
             * path pays a GPR<->SIMD round-trip for every __builtin_popcountl
             * (eight per call). Compute all four bases in-lane instead: AND the
             * one-hot BWT words with the position mask, popcount each 64-bit
             * lane (cnt + pairwise-add reduction), and add the checkpoint
             * counts. Bit-identical to the scalar popcount, no cross-domain
             * moves, one load per checkpoint block. */
            const uint64x2_t msp = vdupq_n_u64(mask_sp);
            const uint64x2_t mep = vdupq_n_u64(mask_ep);
            #define BWA3_PC64(v) vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(vcntq_u8(vreinterpretq_u8_u64(v)))))
            const uint64x2_t psp01 = BWA3_PC64(vandq_u64(vld1q_u64(&blk_sp.one_hot_bwt_str[0]), msp));
            const uint64x2_t psp23 = BWA3_PC64(vandq_u64(vld1q_u64(&blk_sp.one_hot_bwt_str[2]), msp));
            const uint64x2_t pep01 = BWA3_PC64(vandq_u64(vld1q_u64(&blk_ep.one_hot_bwt_str[0]), mep));
            const uint64x2_t pep23 = BWA3_PC64(vandq_u64(vld1q_u64(&blk_ep.one_hot_bwt_str[2]), mep));
            #undef BWA3_PC64
            const uint64x2_t occ_sp01 = vaddq_u64(vld1q_u64((const uint64_t*)&blk_sp.cp_count[0]), psp01);
            const uint64x2_t occ_sp23 = vaddq_u64(vld1q_u64((const uint64_t*)&blk_sp.cp_count[2]), psp23);
            const uint64x2_t occ_ep01 = vaddq_u64(vld1q_u64((const uint64_t*)&blk_ep.cp_count[0]), pep01);
            const uint64x2_t occ_ep23 = vaddq_u64(vld1q_u64((const uint64_t*)&blk_ep.cp_count[2]), pep23);
            vst1q_u64((uint64_t*)&occ_sp[0], occ_sp01);
            vst1q_u64((uint64_t*)&occ_sp[2], occ_sp23);
            vst1q_u64((uint64_t*)&s[0], vsubq_u64(occ_ep01, occ_sp01));
            vst1q_u64((uint64_t*)&s[2], vsubq_u64(occ_ep23, occ_sp23));
#else
            for (uint8_t b = 0; b < 4; b++) {
                int64_t occ_s = blk_sp.cp_count[b] + _mm_countbits_64(blk_sp.one_hot_bwt_str[b] & mask_sp);
                int64_t occ_e = blk_ep.cp_count[b] + _mm_countbits_64(blk_ep.one_hot_bwt_str[b] & mask_ep);
                occ_sp[b] = occ_s;
                s[b]      = occ_e - occ_s;
            }
#endif

            int64_t sentinel_offset = 0;
            if ((smem.k <= sentinel_index) && ((smem.k + smem.s) > sentinel_index))
                sentinel_offset = 1;
            l[3] = smem.l + sentinel_offset;
            l[2] = l[3] + s[3];
            l[1] = l[2] + s[2];
            l[0] = l[1] + s[1];

            smem.k = count[a] + occ_sp[a];
            smem.l = l[a];
            smem.s = s[a];
            return smem;
        }

    // ----- Lockstep SMEM batching internals -----
    // Defined in FMI_search.cpp. Phases and per-slot state are scoped to
    // the lockstep driver; not used anywhere else.
    struct BatchSlot;

    void ls_init_slot(BatchSlot *s, int32_t input_idx,
                      const int32_t *query_pos_array,
                      const int32_t *min_intv_array,
                      const int32_t *rid_array,
                      const bseq1_t *seq_,
                      const int32_t *query_cum_len_ar,
                      const uint8_t *enc_qdb);
    void ls_prefetch_cp_occ(const BatchSlot *s);
    void ls_prefetch_cp_occ_t1(const BatchSlot *s);
    void ls_advance_forward_step(BatchSlot *s, const uint8_t *enc_qdb);
    void ls_prepare_backward(BatchSlot *s);
    void ls_advance_backward_step(BatchSlot *s,
                                  const uint8_t *enc_qdb,
                                  int32_t minSeedLen);

    // ----- Lockstep bwtSeed batching internals -----
    struct BwtSeedSlot;

    void bsd_init_slot(BwtSeedSlot *s, int32_t input_idx,
                       const int32_t *max_intv_array,
                       const bseq1_t *seq_,
                       const int32_t *query_cum_len_ar,
                       const uint8_t *enc_qdb);
    void bsd_prefetch_cp_occ(const BwtSeedSlot *s);
    void bsd_prefetch_cp_occ_t1(const BwtSeedSlot *s);
    void bsd_advance_step(BwtSeedSlot *s,
                          const uint8_t *enc_qdb,
                          int32_t minSeedLen);
};

#endif
