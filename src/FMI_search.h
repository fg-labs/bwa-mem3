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
#define SMEM_LOCKSTEP_N 8
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

    int build_index();
    void load_index();

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
    void getSMEMs(uint8_t *enc_qdb,
                  int32_t numReads,
                  int32_t batch_size,
                  int32_t readlength,
                  int32_t minSeedLengh,
                  int32_t numthreads,
                  SMEM *matchArray,
                  int64_t *numTotalSmem);

    /* matchArray must hold at least numReads * max_readlength SMEMs (caller-sized). */
    void getSMEMsOnePosOneThread(uint8_t *enc_qdb,
                                 int16_t *query_pos_array,
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
                                          int16_t *query_pos_array,
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

        SMEM backwardExt(SMEM smem, uint8_t a);

    // ----- Lockstep SMEM batching internals -----
    // Defined in FMI_search.cpp. Phases and per-slot state are scoped to
    // the lockstep driver; not used anywhere else.
    struct BatchSlot;

    void ls_init_slot(BatchSlot *s, int32_t input_idx,
                      const int16_t *query_pos_array,
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
};

#endif
