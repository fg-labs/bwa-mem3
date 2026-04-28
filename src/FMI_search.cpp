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

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <climits>
#include <sys/mman.h>     /* munmap */
#include "bwa_madvise.h"
#include "bwa_shm.h"
#include "FMI_search.h"
#include "memcpy_bwamem.h"
#include "profiling.h"
#include "libsais_build.h"

#include "safestringlib.h"

FMI_search::FMI_search(const char *fname)
{
    fprintf(stderr, "* Entering FMI_search\n");
    strcpy_s(file_name, PATH_MAX, fname);
    reference_seq_len = 0;
    sentinel_index = 0;
    sa_ls_word = NULL;
    sa_ms_byte = NULL;
    cp_occ = NULL;
    shm_base = NULL;
    shm_len = 0;

    /* one_hot_mask_array is constant across the FMI's lifetime and is
     * identical on disk and shm paths; initialize it once at construction. */
    int64_t one_hot_bytes = 64 * (int64_t)sizeof(uint64_t);
    one_hot_mask_array = (uint64_t *)_mm_malloc(one_hot_bytes, 64);
    assert_not_null(one_hot_mask_array, one_hot_bytes, one_hot_bytes);
    one_hot_mask_array[0] = 0;
    uint64_t base = 0x8000000000000000L;
    one_hot_mask_array[1] = base;
    for (int64_t i = 2; i < 64; ++i) {
        one_hot_mask_array[i] = (one_hot_mask_array[i - 1] >> 1) | base;
    }
}

FMI_search::~FMI_search()
{
    /* When attached from shm, cp_occ/sa_*_byte/sa_*_word alias mmap'd pages
     * owned by the shm segment, so we don't _mm_free them. We do munmap the
     * mapping itself — leaving it would leak VA in long-lived processes that
     * construct and destroy FMI_search repeatedly. */
    if (shm_base != NULL) {
        munmap(shm_base, shm_len);
        shm_base = NULL;
        shm_len  = 0;
    } else {
        if (sa_ms_byte) _mm_free(sa_ms_byte);
        if (sa_ls_word) _mm_free(sa_ls_word);
        if (cp_occ)     _mm_free(cp_occ);
    }
    if (one_hot_mask_array) _mm_free(one_hot_mask_array);
}

int64_t FMI_search::cp_occ_size_bytes() const {
    return ((reference_seq_len >> CP_SHIFT) + 1) * (int64_t)sizeof(CP_OCC);
}

int64_t FMI_search::sa_sample_count() const {
    return (reference_seq_len >> SA_COMPX) + 1;
}

void FMI_search::load_index_from_shm(uint8_t *base, size_t len)
{
    if (base == NULL) {
        fprintf(stderr, "ERROR! load_index_from_shm called with NULL base\n");
        exit(EXIT_FAILURE);
    }
    uint64_t off = 0, sz = 0;

    if (bwa_shm_section_find(base, BWA_SHM_SEC_FMI_SCALARS, &off, &sz) != 0
        || sz != BWA_SHM_FMI_SCALARS_BYTES) {
        fprintf(stderr, "ERROR! shm segment missing or malformed FMI_SCALARS\n");
        exit(EXIT_FAILURE);
    }
    memcpy(&reference_seq_len, base + off,                                    sizeof(int64_t));
    memcpy(count,              base + off + sizeof(int64_t),                  sizeof(int64_t) * 5);
    memcpy(&sentinel_index,    base + off + sizeof(int64_t) * 6,              sizeof(int64_t));

    /* Validate scalars before we use reference_seq_len in cp_occ_size_bytes()
     * and the SA size accessors. Bounds match the disk path's asserts in
     * load_index() and the writer-side checks in bwa_shm_compute(). A corrupt
     * segment would otherwise drive negative or overflowing section sizes. */
    if (reference_seq_len <= 0 || reference_seq_len > 0x7fffffffffLL) {
        fprintf(stderr,
            "ERROR! shm FMI_SCALARS: reference_seq_len=%lld out of bounds\n",
            (long long)reference_seq_len);
        exit(EXIT_FAILURE);
    }
    /* count[] is +1-adjusted by bwa_shm_compute; range [1, ref_seq_len+1]. */
    for (int i = 0; i < 5; ++i) {
        if (count[i] < 0 || count[i] > reference_seq_len + 1) {
            fprintf(stderr,
                "ERROR! shm FMI_SCALARS: count[%d]=%lld out of bounds (ref_seq_len=%lld)\n",
                i, (long long)count[i], (long long)reference_seq_len);
            exit(EXIT_FAILURE);
        }
    }
    if (sentinel_index < 0 || sentinel_index >= reference_seq_len) {
        fprintf(stderr,
            "ERROR! shm FMI_SCALARS: sentinel_index=%lld out of bounds (ref_seq_len=%lld)\n",
            (long long)sentinel_index, (long long)reference_seq_len);
        exit(EXIT_FAILURE);
    }

    if (bwa_shm_section_find(base, BWA_SHM_SEC_FMI_CP_OCC, &off, &sz) != 0
        || (int64_t)sz != cp_occ_size_bytes()) {
        fprintf(stderr, "ERROR! shm segment missing or sized FMI_CP_OCC\n");
        exit(EXIT_FAILURE);
    }
    cp_occ = (CP_OCC *)(base + off);

    if (bwa_shm_section_find(base, BWA_SHM_SEC_FMI_SA_MS, &off, &sz) != 0
        || (int64_t)sz != sa_ms_byte_size_bytes()) {
        fprintf(stderr, "ERROR! shm segment missing or sized FMI_SA_MS\n");
        exit(EXIT_FAILURE);
    }
    sa_ms_byte = (int8_t *)(base + off);

    if (bwa_shm_section_find(base, BWA_SHM_SEC_FMI_SA_LS, &off, &sz) != 0
        || (int64_t)sz != sa_ls_word_size_bytes()) {
        fprintf(stderr, "ERROR! shm segment missing or sized FMI_SA_LS\n");
        exit(EXIT_FAILURE);
    }
    sa_ls_word = (uint32_t *)(base + off);

    shm_base = base;
    shm_len  = len;

    fprintf(stderr, "* FMI attached from shm: ref_seq_len=%ld sentinel_index=%ld\n",
            (long)reference_seq_len, (long)sentinel_index);
}

int FMI_search::build_index() {

    char *prefix = file_name;

    // Read single-strand length from .ann to compute doubled pac_len.
    char ann_path[PATH_MAX];
    strcpy_s(ann_path, PATH_MAX, prefix);
    strcat_s(ann_path, PATH_MAX, ".ann");
    FILE* fann = fopen(ann_path, "r");
    if (fann == NULL) {
        fprintf(stderr, "ERROR: cannot open '%s'\n", ann_path);
        return 1;
    }
    int64_t l_pac = 0;
    int n_seqs = 0, seed = 0;
    if (fscanf(fann, "%" SCNd64 " %d %d", &l_pac, &n_seqs, &seed) != 3) {
        fprintf(stderr, "ERROR: malformed '%s'\n", ann_path);
        fclose(fann);
        return 1;
    }
    fclose(fann);
    // Defensive: a non-positive l_pac (corrupt or zero-length .ann) would
    // pass a bad pac_len into libsais and only fail much later. Catch it
    // up front with an actionable message, mirroring bntseq_restore_core.
    if (l_pac <= 0 || n_seqs < 0) {
        fprintf(stderr, "ERROR: malformed '%s' (l_pac=%" PRId64 ", n_seqs=%d)\n",
                ann_path, l_pac, n_seqs);
        return 1;
    }
    int64_t pac_len = 2 * l_pac;

    auto parse_ll = [](const char* s, const char* name,
                       long long max_val = LLONG_MAX) -> long long {
        char* end = nullptr;
        errno = 0;
        long long v = strtoll(s, &end, 10);
        if (errno || end == s || *end != '\0' || v <= 0 || v > max_val) {
            fprintf(stderr, "ERROR: invalid %s='%s' (expected positive integer)\n",
                    name, s);
            exit(1);
        }
        return v;
    };
    LibsaisBuildOpts opts;
    if (const char* th = getenv("BWA_INDEX_THREADS"))
        opts.num_threads      = (int)parse_ll(th, "BWA_INDEX_THREADS", INT_MAX);
    if (const char* mm = getenv("BWA_INDEX_MAX_MEMORY"))
        opts.max_memory_bytes = parse_ll(mm, "BWA_INDEX_MAX_MEMORY");
    if (const char* td = getenv("BWA_INDEX_TMPDIR"))
        opts.tmpdir           = td;
    return libsais_build_fm_index(prefix, pac_len, opts);
}

void FMI_search::load_index()
{
    /* Try the staged shm segment first. On hit, both the FMI internals
     * (cp_occ / sa_*) and the BNS+PAC are attached as views into the
     * mapping; the segment lifetime belongs to the shm-loading process. */
    {
        size_t   shm_attach_len = 0;
        uint8_t *shm_base_local = bwa_shm_attach(file_name, &shm_attach_len);
        if (shm_base_local != NULL) {
            load_index_from_shm(shm_base_local, shm_attach_len);
            bwa_idx_load_ele_from_shm(shm_base_local, shm_attach_len);
            fprintf(stderr, "* FMI+BNS+PAC attached from shm; "
                    "skipping disk load.\n");
            return;
        }
    }

    // Running total of index bytes allocated so far in this function.
    // Passed to assert_not_null so a failed allocation reports both the
    // attempted size and how much we'd already committed before failing.
    int64_t index_alloc = 0;

    char *ref_file_name = file_name;
    //beCalls = 0;
    char cp_file_name[PATH_MAX];
    strcpy_s(cp_file_name, PATH_MAX, ref_file_name);
    strcat_s(cp_file_name, PATH_MAX, CP_FILENAME_SUFFIX);

    // Read the BWT and FM index of the reference sequence
    FILE *cpstream = NULL;
    cpstream = fopen(cp_file_name,"rb");
    if (cpstream == NULL)
    {
        fprintf(stderr, "ERROR! Unable to open the file: %s\n", cp_file_name);
        exit(EXIT_FAILURE);
    }
    else
    {
        fprintf(stderr, "* Index file found. Loading index from %s\n", cp_file_name);
    }

    err_fread_noeof(&reference_seq_len, sizeof(int64_t), 1, cpstream);
    assert(reference_seq_len > 0);
    assert(reference_seq_len <= 0x7fffffffffL);

    fprintf(stderr, "* Reference seq len for bi-index = %ld\n", reference_seq_len);

    // create checkpointed occ
    int64_t cp_occ_size = (reference_seq_len >> CP_SHIFT) + 1;
    cp_occ = NULL;

    err_fread_noeof(&count[0], sizeof(int64_t), 5, cpstream);
    int64_t cp_occ_bytes = cp_occ_size * sizeof(CP_OCC);
    cp_occ = (CP_OCC *)_mm_malloc(cp_occ_bytes, 64);
    index_alloc += cp_occ_bytes;
    assert_not_null(cp_occ, cp_occ_bytes, index_alloc);
    bwamem_madv_hugepage(cp_occ, cp_occ_bytes);

    err_fread_noeof(cp_occ, sizeof(CP_OCC), cp_occ_size, cpstream);
    int64_t ii = 0;
    for(ii = 0; ii < 5; ii++)// update read count structure
    {
        count[ii] = count[ii] + 1;
    }

    #if SA_COMPRESSION

    int64_t reference_seq_len_ = (reference_seq_len >> SA_COMPX) + 1;
    int64_t sa_ms_bytes = reference_seq_len_ * sizeof(int8_t);
    int64_t sa_ls_bytes = reference_seq_len_ * sizeof(uint32_t);
    sa_ms_byte = (int8_t *)_mm_malloc(sa_ms_bytes, 64);
    index_alloc += sa_ms_bytes;
    assert_not_null(sa_ms_byte, sa_ms_bytes, index_alloc);
    sa_ls_word = (uint32_t *)_mm_malloc(sa_ls_bytes, 64);
    index_alloc += sa_ls_bytes;
    assert_not_null(sa_ls_word, sa_ls_bytes, index_alloc);
    bwamem_madv_hugepage(sa_ms_byte, sa_ms_bytes);
    bwamem_madv_hugepage(sa_ls_word, sa_ls_bytes);
    err_fread_noeof(sa_ms_byte, sizeof(int8_t), reference_seq_len_, cpstream);
    err_fread_noeof(sa_ls_word, sizeof(uint32_t), reference_seq_len_, cpstream);

    #else

    int64_t sa_ms_bytes = reference_seq_len * sizeof(int8_t);
    int64_t sa_ls_bytes = reference_seq_len * sizeof(uint32_t);
    sa_ms_byte = (int8_t *)_mm_malloc(sa_ms_bytes, 64);
    index_alloc += sa_ms_bytes;
    assert_not_null(sa_ms_byte, sa_ms_bytes, index_alloc);
    sa_ls_word = (uint32_t *)_mm_malloc(sa_ls_bytes, 64);
    index_alloc += sa_ls_bytes;
    assert_not_null(sa_ls_word, sa_ls_bytes, index_alloc);
    bwamem_madv_hugepage(sa_ms_byte, sa_ms_bytes);
    bwamem_madv_hugepage(sa_ls_word, sa_ls_bytes);
    err_fread_noeof(sa_ms_byte, sizeof(int8_t), reference_seq_len, cpstream);
    err_fread_noeof(sa_ls_word, sizeof(uint32_t), reference_seq_len, cpstream);

    #endif

    sentinel_index = -1;
    #if SA_COMPRESSION
    err_fread_noeof(&sentinel_index, sizeof(int64_t), 1, cpstream);
    fprintf(stderr, "* sentinel-index: %ld\n", sentinel_index);
    #endif
    fclose(cpstream);

    int64_t x;
    #if !SA_COMPRESSION
    for(x = 0; x < reference_seq_len; x++)
    {
        // fprintf(stderr, "x: %ld\n", x);
        #if SA_COMPRESSION
        if(get_sa_entry_compressed(x) == 0) {
            sentinel_index = x;
            break;
        }
        #else
        if(get_sa_entry(x) == 0) {
            sentinel_index = x;
            break;
        }
        #endif
    }
    fprintf(stderr, "\nsentinel_index: %ld\n", x);    
    #endif

    fprintf(stderr, "* Count:\n");
    for(x = 0; x < 5; x++)
    {
        fprintf(stderr, "%ld,\t%lu\n", x, (unsigned long)count[x]);
    }
    fprintf(stderr, "\n");  

    fprintf(stderr, "* Reading other elements of the index from files %s\n",
            ref_file_name);
    bwa_idx_load_ele(ref_file_name, BWA_IDX_ALL);

    fprintf(stderr, "* Done reading Index!!\n");
}

void FMI_search::getSMEMsOnePosOneThread(uint8_t *enc_qdb,
                                         int16_t *query_pos_array,
                                         int32_t *min_intv_array,
                                         int32_t *rid_array,
                                         int32_t numReads,
                                         int32_t batch_size,
                                         const bseq1_t *seq_,
                                         int32_t *query_cum_len_ar,
                                         int32_t max_readlength,
                                         int32_t minSeedLen,
                                         SMEM *matchArray,
                                         int64_t *__numTotalSmem)
{
    int64_t numTotalSmem = *__numTotalSmem;
    // Heap-allocate to avoid stack overflow on long reads (PacBio HiFi /
    // ONT 1 Mbp+); mirrors the lockstep path's max_readlength sizing.
    int64_t prevArray_bytes = (int64_t)max_readlength * sizeof(SMEM);
    SMEM *prevArray = (SMEM *)_mm_malloc((size_t)prevArray_bytes, 64);
    assert_not_null(prevArray, prevArray_bytes, prevArray_bytes);

    uint32_t i;
    // Perform SMEM for original reads
    for(i = 0; i < numReads; i++)
    {
        int x = query_pos_array[i];
        int32_t rid = rid_array[i];
        int next_x = x + 1;

        int readlength = seq_[rid].l_seq;
        int offset = query_cum_len_ar[rid];
        // uint8_t a = enc_qdb[rid * readlength + x];
        uint8_t a = enc_qdb[offset + x];

        if(a < 4)
        {
            SMEM smem;
            smem.rid = rid;
            smem.m = x;
            smem.n = x;
            smem.k = count[a];
            smem.l = count[3 - a];
            smem.s = count[a+1] - count[a];
            int numPrev = 0;
            
            int j;
            for(j = x + 1; j < readlength; j++)
            {
                // a = enc_qdb[rid * readlength + j];
                a = enc_qdb[offset + j];
                next_x = j + 1;
                if(a < 4)
                {
                    SMEM smem_ = smem;

                    // Forward extension is backward extension with the BWT of reverse complement
                    smem_.k = smem.l;
                    smem_.l = smem.k;
                    SMEM newSmem_ = backwardExt(smem_, 3 - a);
                    //SMEM newSmem_ = forwardExt(smem_, 3 - a);
                    SMEM newSmem = newSmem_;
                    newSmem.k = newSmem_.l;
                    newSmem.l = newSmem_.k;
                    newSmem.n = j;

                    int32_t s_neq_mask = newSmem.s != smem.s;

                    prevArray[numPrev] = smem;
                    numPrev += s_neq_mask;
                    if(newSmem.s < min_intv_array[i])
                    {
                        next_x = j;
                        break;
                    }
                    smem = newSmem;
#ifdef ENABLE_PREFETCH
                    _mm_prefetch((const char *)(&cp_occ[(smem.k) >> CP_SHIFT]), _MM_HINT_T0);
                    _mm_prefetch((const char *)(&cp_occ[(smem.l) >> CP_SHIFT]), _MM_HINT_T0);
#endif
                }
                else
                {
                    break;
                }
            }
            if(smem.s >= min_intv_array[i])
            {

                prevArray[numPrev] = smem;
                numPrev++;
            }

            SMEM *prev;
            prev = prevArray;

            int p;
            for(p = 0; p < (numPrev/2); p++)
            {
                SMEM temp = prev[p];
                prev[p] = prev[numPrev - p - 1];
                prev[numPrev - p - 1] = temp;
            }

            // Backward search
            int cur_j = readlength;
            for(j = x - 1; j >= 0; j--)
            {
                int numCurr = 0;
                int curr_s = -1;
                // a = enc_qdb[rid * readlength + j];
                a = enc_qdb[offset + j];

                if(a > 3)
                {
                    break;
                }
                for(p = 0; p < numPrev; p++)
                {
                    SMEM smem = prev[p];
                    SMEM newSmem = backwardExt(smem, a);
                    newSmem.m = j;

                    if((newSmem.s < min_intv_array[i]) && ((smem.n - smem.m + 1) >= minSeedLen))
                    {
                        cur_j = j;

                        matchArray[numTotalSmem++] = smem;
                        break;
                    }
                    if((newSmem.s >= min_intv_array[i]) && (newSmem.s != curr_s))
                    {
                        curr_s = newSmem.s;
                        prev[numCurr++] = newSmem;
#ifdef ENABLE_PREFETCH
                        _mm_prefetch((const char *)(&cp_occ[(newSmem.k) >> CP_SHIFT]), _MM_HINT_T0);
                        _mm_prefetch((const char *)(&cp_occ[(newSmem.k + newSmem.s) >> CP_SHIFT]), _MM_HINT_T0);
#endif
                        break;
                    }
                }
                p++;
                for(; p < numPrev; p++)
                {
                    SMEM smem = prev[p];

                    SMEM newSmem = backwardExt(smem, a);
                    newSmem.m = j;


                    if((newSmem.s >= min_intv_array[i]) && (newSmem.s != curr_s))
                    {
                        curr_s = newSmem.s;
                        prev[numCurr++] = newSmem;
#ifdef ENABLE_PREFETCH
                        _mm_prefetch((const char *)(&cp_occ[(newSmem.k) >> CP_SHIFT]), _MM_HINT_T0);
                        _mm_prefetch((const char *)(&cp_occ[(newSmem.k + newSmem.s) >> CP_SHIFT]), _MM_HINT_T0);
#endif
                    }
                }
                numPrev = numCurr;
                if(numCurr == 0)
                {
                    break;
                }
            }
            if(numPrev != 0)
            {
                SMEM smem = prev[0];
                if(((smem.n - smem.m + 1) >= minSeedLen))
                {

                    matchArray[numTotalSmem++] = smem;
                }
                numPrev = 0;
            }
        }
        query_pos_array[i] = next_x;
    }
    (*__numTotalSmem) = numTotalSmem;
    _mm_free(prevArray);
}

// ===== Lockstep SMEM batching =====
// Per-slot state for the lockstep SMEM walk. One instance per in-flight read
// in the batch. Every field mirrors a per-read local in the scalar
// getSMEMsOnePosOneThread body, so parity with the scalar path follows from
// composing the existing primitives (backwardExt, count[], cp_occ) in the
// same sequence the scalar would.
//
// The prev[] and match_buf[] buffers are heap-allocated by the driver
// (getSMEMsOnePosOneThread_lockstep), sized from the batch's
// max_readlength. No compile-time cap — long reads (PacBio HiFi, ONT)
// fit cleanly.

enum LockstepPhase : uint8_t {
    PH_FWD      = 0,  // forward extension inner loop active
    PH_BWD_INIT = 1,  // between-phases housekeeping pending
    PH_BWD      = 2,  // backward search outer loop active
    PH_DONE     = 3   // slot finished; match_buf ready for flush
};

// LISA trick #4: hybrid SoA. The per-slot "hot" state stays compact
// (~80 bytes) so cross-slot access (e.g. T1 prefetch lookahead against
// slots[(s+N/2)%N]) hits a tight L1-resident array instead of jumping
// 32KB strides across embedded buffers. The bulk arrays (prev[],
// match_buf[]; ~32 KB each) live in separate per-slot allocations,
// referenced by pointer. Accessor syntax (s->prev[i], s->match_buf[i])
// is unchanged from the embedded-array version.
struct FMI_search::BatchSlot {
    // Input identity — copied at init, never mutated thereafter.
    int32_t input_idx;           // index into the caller's input arrays
    int32_t rid;                 // rid_array[input_idx]
    int16_t start_pos;           // query_pos_array[input_idx] (saved for bwd init)
    int32_t min_intv;            // min_intv_array[input_idx]
    int32_t readlength;          // seq_[rid].l_seq
    int32_t offset;              // query_cum_len_ar[rid]

    // Output for query_pos_array[input_idx] write-back at flush time.
    int16_t next_x;

    // Walk state (mirrors scalar locals).
    SMEM smem;                   // current SA interval
    int32_t j;                   // current query position in the active phase's inner loop
    int32_t cur_j;               // backward phase bookkeeping (scalar's cur_j)
    LockstepPhase phase;

    // Per-slot bulk-buffer state.
    int32_t numPrev;
    int32_t match_count;
    bool    ready;

    // Pointers into per-slot bulk arrays (allocated by the driver, sized
    // by max_readlength).
    SMEM    *prev;
    SMEM    *match_buf;
};

void FMI_search::ls_prefetch_cp_occ(const BatchSlot *s)
{
#ifdef ENABLE_PREFETCH
    _mm_prefetch((const char *)(&cp_occ[(s->smem.k) >> CP_SHIFT]), _MM_HINT_T0);
    _mm_prefetch((const char *)(&cp_occ[(s->smem.l) >> CP_SHIFT]), _MM_HINT_T0);
    _mm_prefetch((const char *)(&cp_occ[(s->smem.k + s->smem.s) >> CP_SHIFT]), _MM_HINT_T0);
    _mm_prefetch((const char *)(&cp_occ[(s->smem.l + s->smem.s) >> CP_SHIFT]), _MM_HINT_T0);
#else
    (void)s;
#endif
}

/* T1 (L2) variant — used for cross-slot N/2-step lookahead. The T0 prefetch
 * above lands at "this slot's next step" granularity; T1 here lands at
 * "this slot's step ~N/2 from now". For N=8 that's 4 stepping-passes ahead,
 * giving DRAM-latency-class lookahead when the cp_occ working set spills
 * out of L2/L3. */
void FMI_search::ls_prefetch_cp_occ_t1(const BatchSlot *s)
{
#ifdef ENABLE_PREFETCH
    _mm_prefetch((const char *)(&cp_occ[(s->smem.k) >> CP_SHIFT]), _MM_HINT_T1);
    _mm_prefetch((const char *)(&cp_occ[(s->smem.l) >> CP_SHIFT]), _MM_HINT_T1);
    _mm_prefetch((const char *)(&cp_occ[(s->smem.k + s->smem.s) >> CP_SHIFT]), _MM_HINT_T1);
    _mm_prefetch((const char *)(&cp_occ[(s->smem.l + s->smem.s) >> CP_SHIFT]), _MM_HINT_T1);
#else
    (void)s;
#endif
}

// Populate a slot from the caller's input arrays at the given input_idx.
// After this call either:
//   phase == PH_FWD  — slot is ready to step the forward-extension inner loop
//   phase == PH_DONE — the first base is non-ACGT; scalar skips this read, so
//                      we match that (zero matches emitted, ready to flush).
void FMI_search::ls_init_slot(BatchSlot *s,
                              int32_t input_idx,
                              const int16_t *query_pos_array,
                              const int32_t *min_intv_array,
                              const int32_t *rid_array,
                              const bseq1_t *seq_,
                              const int32_t *query_cum_len_ar,
                              const uint8_t *enc_qdb)
{
    s->input_idx   = input_idx;
    s->rid         = rid_array[input_idx];
    s->start_pos   = query_pos_array[input_idx];
    s->min_intv    = min_intv_array[input_idx];
    s->readlength  = seq_[s->rid].l_seq;
    s->offset      = query_cum_len_ar[s->rid];
    s->next_x      = s->start_pos + 1;
    s->numPrev     = 0;
    s->match_count = 0;
    s->ready       = false;

    int32_t x = s->start_pos;
    uint8_t a = enc_qdb[s->offset + x];
    if (a < 4) {
        s->smem.rid = s->rid;
        s->smem.m   = x;
        s->smem.n   = x;
        s->smem.k   = count[a];
        s->smem.l   = count[3 - a];
        s->smem.s   = count[a + 1] - count[a];
        s->j        = x + 1;
        s->phase    = PH_FWD;
    } else {
        // Scalar path skips the whole read when the first base is N.
        // Match that by going straight to DONE with zero matches.
        s->phase = PH_DONE;
        s->ready = true;
    }
}
// Advance one slot through one step of forward extension.
// Mirrors the inner j-loop body of the scalar getSMEMsOnePosOneThread
// (src/FMI_search.cpp — the for(j = x+1; j < readlength; j++) loop).
// On the step that exits forward-ext (end-of-read, non-ACGT at j, or
// s < min_intv), transitions phase to PH_BWD_INIT.
void FMI_search::ls_advance_forward_step(BatchSlot *s, const uint8_t *enc_qdb)
{
    if (s->j >= s->readlength) {
        // Ran off the end of the read; keep the still-valid smem if any.
        if (s->smem.s >= s->min_intv) {
            s->prev[s->numPrev++] = s->smem;
        }
        s->phase = PH_BWD_INIT;
        return;
    }

    uint8_t a = enc_qdb[s->offset + s->j];
    s->next_x = s->j + 1;
    if (a >= 4) {
        // Non-ACGT base terminates forward ext.
        if (s->smem.s >= s->min_intv) {
            s->prev[s->numPrev++] = s->smem;
        }
        s->phase = PH_BWD_INIT;
        return;
    }

    SMEM smem_ = s->smem;
    smem_.k = s->smem.l;
    smem_.l = s->smem.k;
    SMEM newSmem_ = backwardExt(smem_, 3 - a);
    SMEM newSmem  = newSmem_;
    newSmem.k = newSmem_.l;
    newSmem.l = newSmem_.k;
    newSmem.n = s->j;

    int32_t s_neq_mask = (newSmem.s != s->smem.s);
    s->prev[s->numPrev] = s->smem;
    s->numPrev += s_neq_mask;

    if (newSmem.s < s->min_intv) {
        s->next_x = s->j;
        s->phase = PH_BWD_INIT;
        return;
    }

    s->smem = newSmem;
    s->j++;
#ifdef ENABLE_PREFETCH
    _mm_prefetch((const char *)(&cp_occ[(s->smem.k) >> CP_SHIFT]), _MM_HINT_T0);
    _mm_prefetch((const char *)(&cp_occ[(s->smem.l) >> CP_SHIFT]), _MM_HINT_T0);
#endif
}

// Between-phases housekeeping: reverse prev[] and set up backward-phase
// cursors. After this call the slot is either PH_BWD (ready for
// ls_advance_backward_step) or PH_DONE (nothing to do — backward outer
// loop would terminate on the first step, so we short-circuit to the
// final prev[0] emit and done).
void FMI_search::ls_prepare_backward(BatchSlot *s)
{
    // Reverse prev[] in place (matches scalar behavior before backward loop).
    for (int p = 0; p < (s->numPrev / 2); p++) {
        SMEM tmp = s->prev[p];
        s->prev[p] = s->prev[s->numPrev - p - 1];
        s->prev[s->numPrev - p - 1] = tmp;
    }
    s->cur_j = s->readlength;
    s->j = s->start_pos - 1;  // first position for the backward outer loop

    if (s->numPrev == 0) {
        // Nothing to emit; scalar's final `if (numPrev != 0)` block is a no-op.
        s->phase = PH_DONE;
        s->ready = true;
        return;
    }
    if (s->j < 0) {
        // Backward outer loop cannot execute (start_pos == 0). Transition to
        // PH_BWD so the first ls_advance_backward_step call sees j < 0, jumps
        // to DONE, and runs the final prev[0] emit with the correct minSeedLen.
        s->phase = PH_BWD;
        return;
    }
    s->phase = PH_BWD;
}

// Advance one slot through ONE outer-j iteration of backward search.
// Mirrors the body of the scalar `for (j = x-1; j >= 0; j--)` outer loop
// in getSMEMsOnePosOneThread. On terminating conditions (numCurr == 0,
// j < 0 after decrement, or non-ACGT at j), runs the scalar's final
// prev[0] emit and transitions to PH_DONE.
void FMI_search::ls_advance_backward_step(BatchSlot *s,
                                          const uint8_t *enc_qdb,
                                          int32_t minSeedLen)
{
    if (s->j < 0) goto DONE;

    {
        int32_t numCurr = 0;
        int32_t curr_s = -1;
        uint8_t a = enc_qdb[s->offset + s->j];
        if (a > 3) goto DONE;

        int p;
        for (p = 0; p < s->numPrev; p++) {
            SMEM smem = s->prev[p];
            SMEM newSmem = backwardExt(smem, a);
            newSmem.m = s->j;
            if ((newSmem.s < s->min_intv) && ((smem.n - smem.m + 1) >= minSeedLen)) {
                s->cur_j = s->j;
                s->match_buf[s->match_count++] = smem;
                break;
            }
            if ((newSmem.s >= s->min_intv) && (newSmem.s != curr_s)) {
                curr_s = newSmem.s;
                s->prev[numCurr++] = newSmem;
#ifdef ENABLE_PREFETCH
                _mm_prefetch((const char *)(&cp_occ[(newSmem.k) >> CP_SHIFT]), _MM_HINT_T0);
                _mm_prefetch((const char *)(&cp_occ[(newSmem.k + newSmem.s) >> CP_SHIFT]), _MM_HINT_T0);
#endif
                break;
            }
        }
        p++;
        for (; p < s->numPrev; p++) {
            SMEM smem = s->prev[p];
            SMEM newSmem = backwardExt(smem, a);
            newSmem.m = s->j;
            if ((newSmem.s >= s->min_intv) && (newSmem.s != curr_s)) {
                curr_s = newSmem.s;
                s->prev[numCurr++] = newSmem;
#ifdef ENABLE_PREFETCH
                _mm_prefetch((const char *)(&cp_occ[(newSmem.k) >> CP_SHIFT]), _MM_HINT_T0);
                _mm_prefetch((const char *)(&cp_occ[(newSmem.k + newSmem.s) >> CP_SHIFT]), _MM_HINT_T0);
#endif
            }
        }
        s->numPrev = numCurr;
        s->j--;                // advance for the next outer step
        if (numCurr == 0) goto DONE;
        return;                // stay in PH_BWD; next call continues
    }

DONE:
    // Scalar end-of-function final prev[0] emit (lines ~650-659 of scalar).
    if (s->numPrev != 0) {
        SMEM smem = s->prev[0];
        if ((smem.n - smem.m + 1) >= minSeedLen) {
            s->match_buf[s->match_count++] = smem;
        }
        s->numPrev = 0;
    }
    s->phase = PH_DONE;
    s->ready = true;
}
// ===== End lockstep SMEM batching =====

void FMI_search::getSMEMsAllPosOneThread(uint8_t *enc_qdb,
                                         int32_t *min_intv_array,
                                         int32_t *rid_array,
                                         int32_t numReads,
                                         int32_t batch_size,
                                         const bseq1_t *seq_,
                                         int32_t *query_cum_len_ar,
                                         int32_t max_readlength,
                                         int32_t minSeedLen,
                                         SMEM *matchArray,
                                         int64_t *__numTotalSmem)
{
    int16_t *query_pos_array = (int16_t *)_mm_malloc(numReads * sizeof(int16_t), 64);

    int32_t i;
    for(i = 0; i < numReads; i++)
        query_pos_array[i] = 0;

    int32_t numActive = numReads;
    (*__numTotalSmem) = 0;

    do
    {
        int32_t head = 0;
        int32_t tail = 0;
        for(head = 0; head < numActive; head++)
        {
            int readlength = seq_[rid_array[head]].l_seq;
            if(query_pos_array[head] < readlength)
            {
                rid_array[tail] = rid_array[head];
                query_pos_array[tail] = query_pos_array[head];
                min_intv_array[tail] = min_intv_array[head];
                tail++;
            }
        }
#if SMEM_LOCKSTEP_N > 1
        getSMEMsOnePosOneThread_lockstep(enc_qdb,
                                         query_pos_array,
                                         min_intv_array,
                                         rid_array,
                                         tail,
                                         batch_size,
                                         seq_,
                                         query_cum_len_ar,
                                         max_readlength,
                                         minSeedLen,
                                         matchArray,
                                         __numTotalSmem);
#else
        getSMEMsOnePosOneThread(enc_qdb,
                                query_pos_array,
                                min_intv_array,
                                rid_array,
                                tail,
                                batch_size,
                                seq_,
                                query_cum_len_ar,
                                max_readlength,
                                minSeedLen,
                                matchArray,
                                __numTotalSmem);
#endif
        numActive = tail;
    } while(numActive > 0);

    _mm_free(query_pos_array);
}

void FMI_search::getSMEMsOnePosOneThread_lockstep(uint8_t *enc_qdb,
                                                   int16_t *query_pos_array,
                                                   int32_t *min_intv_array,
                                                   int32_t *rid_array,
                                                   int32_t numReads,
                                                   int32_t batch_size,
                                                   const bseq1_t *seq_,
                                                   int32_t *query_cum_len_ar,
                                                   int32_t max_readlength,
                                                   int32_t minSeedLen,
                                                   SMEM *matchArray,
                                                   int64_t *__numTotalSmem)
{
    (void)batch_size;

    if (numReads == 0) return;

    const int32_t N = SMEM_LOCKSTEP_N;
    // LISA trick #4: hybrid SoA layout. `slots[]` holds only the small hot
    // state (~80 B per slot, full array fits in 1-2 cache lines for N=8).
    // Bulk per-slot buffers (prev/match_buf) live separately and are reused
    // across calls via thread_local caches sized from the batch's
    // max_readlength so reads of any length fit (see issue #44 / PR #55).
    //
    // The outer driver getSMEMsAllPosOneThread runs this in a do/while
    // loop, so allocating per-call (~2*N*max_readlength*sizeof(SMEM) ≈
    // 384 KB at N=8, max_readlength=1500) imposed measurable allocator
    // pressure. Cache per-thread (so OMP workers don't share) and grow
    // monotonically — max_readlength is bounded by the driver batch and
    // increases rarely in practice.
    //
    // TODO(memory): the cache only grows. For mixed-length workloads (e.g.
    // a single ONT-class read with max_readlength≈1e6 mid-stream — sized
    // to ~640 MB per thread, ~10 GB across 16 OMP workers — followed by
    // short reads), the high-water mark is held until thread/process exit.
    // Acceptable for the smoke1M Illumina PE150 benchmark (max_readlength≈
    // 150 → ~38 KB/thread). Revisit if a streaming aligner or long-running
    // service pipeline appears: gate the realloc on a configured upper
    // bound (e.g. MAX_SMEM_PER_SLOT) or shrink when cached_per_slot greatly
    // exceeds the current batch. Also flagged by leak-sanitizer/Valgrind
    // because the buffers are released only at process exit.
    BatchSlot slots[SMEM_LOCKSTEP_N] = {};
    static thread_local SMEM   *cached_prev  = NULL;
    static thread_local SMEM   *cached_match = NULL;
    static thread_local size_t  cached_per_slot = 0;
    const size_t per_slot_smems = (size_t)max_readlength;
    if (per_slot_smems > cached_per_slot) {
        if (cached_prev  != NULL) _mm_free(cached_prev);
        if (cached_match != NULL) _mm_free(cached_match);
        const size_t total_slot_bytes = (size_t)N * per_slot_smems * sizeof(SMEM);
        cached_prev  = (SMEM *)_mm_malloc(total_slot_bytes, 64);
        assert_not_null(cached_prev, total_slot_bytes, total_slot_bytes);
        cached_match = (SMEM *)_mm_malloc(total_slot_bytes, 64);
        assert_not_null(cached_match, total_slot_bytes, (size_t)2 * total_slot_bytes);
        cached_per_slot = per_slot_smems;
    }
    for (int32_t s = 0; s < N; s++) {
        slots[s].prev      = cached_prev  + (size_t)s * cached_per_slot;
        slots[s].match_buf = cached_match + (size_t)s * cached_per_slot;
    }

    // Seed the first min(N, numReads) slots.
    int32_t initial = (numReads < N) ? numReads : N;
    for (int32_t s = 0; s < initial; s++) {
        ls_init_slot(&slots[s], s, query_pos_array, min_intv_array,
                     rid_array, seq_, query_cum_len_ar, enc_qdb);
        if (slots[s].phase == PH_FWD) ls_prefetch_cp_occ(&slots[s]);
    }
    // Any unused slots (numReads < N) are marked DONE with invalid input_idx
    // so the stepping pass skips them and the flush pass ignores them.
    for (int32_t s = initial; s < N; s++) {
        slots[s].phase = PH_DONE;
        slots[s].ready = false;
        slots[s].input_idx = -1;
    }

    int32_t next_input   = initial;
    int32_t flush_cursor = 0;
    int64_t numTotalSmem = *__numTotalSmem;

    while (flush_cursor < numReads) {
        // --- Stepping pass: advance each non-DONE slot by one phase-step. ---
        for (int32_t s = 0; s < N; s++) {
            switch (slots[s].phase) {
                case PH_FWD:
                    ls_advance_forward_step(&slots[s], enc_qdb);
                    // If forward just transitioned to PH_BWD_INIT, run the
                    // between-phases housekeeping immediately. prepare_backward
                    // may transition straight to PH_DONE (e.g. empty prev[]).
                    if (slots[s].phase == PH_BWD_INIT) {
                        ls_prepare_backward(&slots[s]);
                    }
                    break;
                case PH_BWD_INIT:
                    ls_prepare_backward(&slots[s]);
                    break;
                case PH_BWD:
                    ls_advance_backward_step(&slots[s], enc_qdb, minSeedLen);
                    break;
                case PH_DONE:
                    break;
            }
            // LISA trick #5: cross-slot T1 (L2) prefetch with N/2-step lookahead.
            // The same-slot T0 prefetch issued inside ls_advance_*_step covers
            // the next single-step access; this T1 prefetch on slot[s+N/2] keeps
            // a copy in L2 for that slot's access ~N/2 stepping-passes from now,
            // hiding DRAM-class latency when cp_occ entries spill out of L3.
            const int32_t s_la = (s + (N / 2)) % N;
            if (slots[s_la].phase == PH_FWD || slots[s_la].phase == PH_BWD) {
                ls_prefetch_cp_occ_t1(&slots[s_la]);
            }
        }

        // --- Flush pass: in-order emit + slot recycle. ---
        bool progress = true;
        while (progress) {
            progress = false;
            for (int32_t s = 0; s < N; s++) {
                if (slots[s].phase == PH_DONE &&
                    slots[s].ready &&
                    slots[s].input_idx == flush_cursor) {
                    for (int32_t m = 0; m < slots[s].match_count; m++) {
                        matchArray[numTotalSmem++] = slots[s].match_buf[m];
                    }
                    query_pos_array[slots[s].input_idx] = slots[s].next_x;
                    flush_cursor++;

                    if (next_input < numReads) {
                        ls_init_slot(&slots[s], next_input,
                                     query_pos_array, min_intv_array,
                                     rid_array, seq_, query_cum_len_ar, enc_qdb);
                        if (slots[s].phase == PH_FWD) ls_prefetch_cp_occ(&slots[s]);
                        next_input++;
                    } else {
                        slots[s].input_idx = -1;  // retired
                        slots[s].ready = false;
                    }
                    progress = true;
                }
            }
        }
    }

    *__numTotalSmem = numTotalSmem;
    /* prev/match buffers are owned by thread_local caches above; intentionally
     * not freed here so the next call (same thread) reuses them without a
     * round-trip through the allocator. The caches leak at thread/process exit. */
}

int64_t FMI_search::bwtSeedStrategyAllPosOneThread(uint8_t *enc_qdb,
                                                   int32_t *max_intv_array,
                                                   int32_t numReads,
                                                   const bseq1_t *seq_,
                                                   int32_t *query_cum_len_ar,
                                                   int32_t minSeedLen,
                                                   SMEM *matchArray)
{
    int32_t i;

    int64_t numTotalSeed = 0;

    for(i = 0; i < numReads; i++)
    {
        int readlength = seq_[i].l_seq;
        int16_t x = 0;
        while(x < readlength)
        {
            int next_x = x + 1;

            // Forward search
            SMEM smem;
            smem.rid = i;
            smem.m = x;
            smem.n = x;
            
            int offset = query_cum_len_ar[i];
            uint8_t a = enc_qdb[offset + x];
            // uint8_t a = enc_qdb[i * readlength + x];

            if(a < 4)
            {
                smem.k = count[a];
                smem.l = count[3 - a];
                smem.s = count[a+1] - count[a];


                int j;
                for(j = x + 1; j < readlength; j++)
                {
                    next_x = j + 1;
                    // a = enc_qdb[i * readlength + j];
                    a = enc_qdb[offset + j];
                    if(a < 4)
                    {
                        SMEM smem_ = smem;

                        // Forward extension is backward extension with the BWT of reverse complement
                        smem_.k = smem.l;
                        smem_.l = smem.k;
                        SMEM newSmem_ = backwardExt(smem_, 3 - a);
                        //SMEM smem = backwardExt(smem, 3 - a);
                        //smem.n = j;
                        SMEM newSmem = newSmem_;
                        newSmem.k = newSmem_.l;
                        newSmem.l = newSmem_.k;
                        newSmem.n = j;
                        smem = newSmem;
#ifdef ENABLE_PREFETCH
                        _mm_prefetch((const char *)(&cp_occ[(smem.k) >> CP_SHIFT]), _MM_HINT_T0);
                        _mm_prefetch((const char *)(&cp_occ[(smem.l) >> CP_SHIFT]), _MM_HINT_T0);
#endif


                        if((smem.s < max_intv_array[i]) && ((smem.n - smem.m + 1) >= minSeedLen))
                        {

                            if(smem.s > 0)
                            {
                                matchArray[numTotalSeed++] = smem;
                            }
                            break;
                        }
                    }
                    else
                    {

                        break;
                    }
                }

            }
            x = next_x;
        }
    }
    return numTotalSeed;
}


void FMI_search::getSMEMs(uint8_t *enc_qdb,
        int32_t numReads,
        int32_t batch_size,
        int32_t readlength,
        int32_t minSeedLen,
        int32_t nthreads,
        SMEM *matchArray,
        int64_t *numTotalSmem)
{
    const size_t smem_arr_bytes = (size_t)nthreads * (size_t)readlength * sizeof(SMEM);
    SMEM *prevArray = (SMEM *)_mm_malloc(smem_arr_bytes, 64);
    assert_not_null(prevArray, smem_arr_bytes, smem_arr_bytes);
    SMEM *currArray = (SMEM *)_mm_malloc(smem_arr_bytes, 64);
    assert_not_null(currArray, smem_arr_bytes, (size_t)2 * smem_arr_bytes);


// #pragma omp parallel num_threads(nthreads)
    {
        int tid = 0; //omp_get_thread_num();   // removed omp
        numTotalSmem[tid] = 0;
        SMEM *myPrevArray = prevArray + tid * readlength;
        SMEM *myCurrArray = currArray + tid * readlength;

        int32_t perThreadQuota = (numReads + (nthreads - 1)) / nthreads;
        int32_t first = tid * perThreadQuota;
        int32_t last  = (tid + 1) * perThreadQuota;
        if(last > numReads) last = numReads;
        SMEM *myMatchArray = matchArray + first * readlength;

        uint32_t i;
        // Perform SMEM for original reads
        for(i = first; i < last; i++)
        {
            int x = readlength - 1;
            int numPrev = 0;
            int numSmem = 0;

            while (x >= 0)
            {
                // Forward search
                SMEM smem;
                smem.rid = i;
                smem.m = x;
                smem.n = x;
                uint8_t a = enc_qdb[i * readlength + x];

                if(a > 3)
                {
                    x--;
                    continue;
                }
                smem.k = count[a];
                smem.l = count[3 - a];
                smem.s = count[a+1] - count[a];

                int j;
                for(j = x + 1; j < readlength; j++)
                {
                    a = enc_qdb[i * readlength + j];
                    if(a < 4)
                    {
                        SMEM smem_ = smem;

                        // Forward extension is backward extension with the BWT of reverse complement
                        smem_.k = smem.l;
                        smem_.l = smem.k;
                        SMEM newSmem_ = backwardExt(smem_, 3 - a);
                        SMEM newSmem = newSmem_;
                        newSmem.k = newSmem_.l;
                        newSmem.l = newSmem_.k;
                        newSmem.n = j;

                        if(newSmem.s != smem.s)
                        {
                            myPrevArray[numPrev] = smem;
                            numPrev++;
                        }
                        smem = newSmem;
                        if(newSmem.s == 0)
                        {
                            break;
                        }
                    }
                    else
                    {
                        myPrevArray[numPrev] = smem;
                        numPrev++;
                        break;
                    }
                }
                if(smem.s != 0)
                {
                    myPrevArray[numPrev++] = smem;
                }

                SMEM *curr, *prev;
                prev = myPrevArray;
                curr = myCurrArray;

                int p;
                for(p = 0; p < (numPrev/2); p++)
                {
                    SMEM temp = prev[p];
                    prev[p] = prev[numPrev - p - 1];
                    prev[numPrev - p - 1] = temp;
                }

                int next_x = x - 1;

                // Backward search
                int cur_j = readlength;
                for(j = x - 1; j >= 0; j--)
                {
                    int numCurr = 0;
                    int curr_s = -1;
                    a = enc_qdb[i * readlength + j];
                    //printf("a = %d\n", a);
                    if(a > 3)
                    {
                        next_x = j - 1;
                        break;
                    }
                    for(p = 0; p < numPrev; p++)
                    {
                        SMEM smem = prev[p];
                        SMEM newSmem = backwardExt(smem, a);
                        newSmem.m = j;

                        if(newSmem.s == 0)
                        {
                            if((numCurr == 0) && (j < cur_j))
                            {
                                cur_j = j;
                                if((smem.n - smem.m + 1) >= minSeedLen)
                                    myMatchArray[numTotalSmem[tid] + numSmem++] = smem;
                            }
                        }
                        if((newSmem.s != 0) && (newSmem.s != curr_s))
                        {
                            curr_s = newSmem.s;
                            curr[numCurr++] = newSmem;
                        }
                    }
                    SMEM *temp = prev;
                    prev = curr;
                    curr = temp;
                    numPrev = numCurr;
                    if(numCurr == 0)
                    {
                        next_x = j;
                        break;
                    }
                    else
                    {
                        next_x = j - 1;
                    }
                }
                if(numPrev != 0)
                {
                    SMEM smem = prev[0];
                    if((smem.n - smem.m + 1) >= minSeedLen)
                        myMatchArray[numTotalSmem[tid] + numSmem++] = smem;
                    numPrev = 0;
                }
                x = next_x;
            }
            numTotalSmem[tid] += numSmem;
        }
    }

    _mm_free(prevArray);
    _mm_free(currArray);
}


int compare_smem(const void *a, const void *b)
{
    SMEM *pa = (SMEM *)a;
    SMEM *pb = (SMEM *)b;

    if(pa->rid < pb->rid)
        return -1;
    if(pa->rid > pb->rid)
        return 1;

    if(pa->m < pb->m)
        return -1;
    if(pa->m > pb->m)
        return 1;
    if(pa->n > pb->n)
        return -1;
    if(pa->n < pb->n)
        return 1;
    return 0;
}

void FMI_search::sortSMEMs(SMEM *matchArray,
        int64_t numTotalSmem[],
        int32_t numReads,
        int32_t readlength,
        int nthreads)
{
    int tid;
    int32_t perThreadQuota = (numReads + (nthreads - 1)) / nthreads;
    for(tid = 0; tid < nthreads; tid++)
    {
        int32_t first = tid * perThreadQuota;
        SMEM *myMatchArray = matchArray + first * readlength;
        qsort(myMatchArray, numTotalSmem[tid], sizeof(SMEM), compare_smem);
    }
}


SMEM FMI_search::backwardExt(SMEM smem, uint8_t a)
{
    //beCalls++;
    uint8_t b;

    int64_t k[4], l[4], s[4];
    for(b = 0; b < 4; b++)
    {
        int64_t sp = (int64_t)(smem.k);
        int64_t ep = (int64_t)(smem.k) + (int64_t)(smem.s);
        GET_OCC(sp, b, occ_id_sp, y_sp, occ_sp, one_hot_bwt_str_c_sp, match_mask_sp);
        GET_OCC(ep, b, occ_id_ep, y_ep, occ_ep, one_hot_bwt_str_c_ep, match_mask_ep);
        k[b] = count[b] + occ_sp;
        s[b] = occ_ep - occ_sp;
    }

    int64_t sentinel_offset = 0;
    if((smem.k <= sentinel_index) && ((smem.k + smem.s) > sentinel_index)) sentinel_offset = 1;
    l[3] = smem.l + sentinel_offset;
    l[2] = l[3] + s[3];
    l[1] = l[2] + s[2];
    l[0] = l[1] + s[1];

    smem.k = k[a];
    smem.l = l[a];
    smem.s = s[a];
    return smem;
}

int64_t FMI_search::get_sa_entry(int64_t pos)
{
    int64_t sa_entry = sa_ms_byte[pos];
    sa_entry = sa_entry << 32;
    sa_entry = sa_entry + sa_ls_word[pos];
    return sa_entry;
}

void FMI_search::get_sa_entries(int64_t *posArray, int64_t *coordArray, uint32_t count, int32_t nthreads)
{
    uint32_t i;
// #pragma omp parallel for num_threads(nthreads)
    for(i = 0; i < count; i++)
    {
        /* Prefetch the SAL slot SAL_PFD iterations ahead. Both
         * sa_ms_byte and sa_ls_word are large random-access arrays
         * (separate cache lines), so issue two prefetches per slot. */
        if (i + SAL_PFD < count) {
            int64_t pf_pos = posArray[i + SAL_PFD];
            _mm_prefetch((const char *)(sa_ms_byte + pf_pos), _MM_HINT_T0);
            _mm_prefetch((const char *)(sa_ls_word + pf_pos), _MM_HINT_T0);
        }
        int64_t pos = posArray[i];
        int64_t sa_entry = sa_ms_byte[pos];
        sa_entry = sa_entry << 32;
        sa_entry = sa_entry + sa_ls_word[pos];
        coordArray[i] = sa_entry;
    }
}

void FMI_search::get_sa_entries(SMEM *smemArray, int64_t *coordArray, int32_t *coordCountArray, uint32_t count, int32_t max_occ)
{
    uint32_t i;
    int32_t totalCoordCount = 0;
    for(i = 0; i < count; i++)
    {
        int32_t c = 0;
        SMEM smem = smemArray[i];
        int64_t hi = smem.k + smem.s;
        int64_t step = (smem.s > max_occ) ? smem.s / max_occ : 1;
        int64_t j;
        for(j = smem.k; (j < hi) && (c < max_occ); j+=step, c++)
        {
            int64_t pos = j;
            /* Prefetch SAL_PFD iterations ahead. Both arrays sit on
             * separate cache lines, so issue both prefetches per slot. */
            int64_t pf_pos = pos + SAL_PFD * step;
            if (pf_pos < hi) {
                _mm_prefetch((const char *)(sa_ms_byte + pf_pos), _MM_HINT_T0);
                _mm_prefetch((const char *)(sa_ls_word + pf_pos), _MM_HINT_T0);
            }
            int64_t sa_entry = sa_ms_byte[pos];
            sa_entry = sa_entry << 32;
            sa_entry = sa_entry + sa_ls_word[pos];
            coordArray[totalCoordCount + c] = sa_entry;
        }
        coordCountArray[i] = c;
        totalCoordCount += c;
    }
}

// sa_compression
int64_t FMI_search::get_sa_entry_compressed(int64_t pos, int tid)
{
    if ((pos & SA_COMPX_MASK) == 0) {
        
        #if  SA_COMPRESSION
        int64_t sa_entry = sa_ms_byte[pos >> SA_COMPX];
        #else
        int64_t sa_entry = sa_ms_byte[pos];     // simulation
        #endif
        
        sa_entry = sa_entry << 32;
        
        #if  SA_COMPRESSION
        sa_entry = sa_entry + sa_ls_word[pos >> SA_COMPX];
        #else
        sa_entry = sa_entry + sa_ls_word[pos];   // simulation
        #endif
        
        return sa_entry;        
    }
    else {
        // tprof[MEM_CHAIN][tid] ++;
        int64_t offset = 0; 
        int64_t sp = pos;
        while(true)
        {
            int64_t occ_id_pp_ = sp >> CP_SHIFT;
            int64_t y_pp_ = CP_BLOCK_SIZE - (sp & CP_MASK) - 1; 
            uint64_t *one_hot_bwt_str = cp_occ[occ_id_pp_].one_hot_bwt_str;
            uint8_t b;

            if((one_hot_bwt_str[0] >> y_pp_) & 1)
                b = 0;
            else if((one_hot_bwt_str[1] >> y_pp_) & 1)
                b = 1;
            else if((one_hot_bwt_str[2] >> y_pp_) & 1)
                b = 2;
            else if((one_hot_bwt_str[3] >> y_pp_) & 1)
                b = 3;
            else
                b = 4;

            if (b == 4) {
                return offset;
            }

            GET_OCC(sp, b, occ_id_sp, y_sp, occ_sp, one_hot_bwt_str_c_sp, match_mask_sp);

            sp = count[b] + occ_sp;
            
            offset ++;
            // tprof[ALIGN1][tid] ++;
            if ((sp & SA_COMPX_MASK) == 0) break;
        }
        // assert((reference_seq_len >> SA_COMPX) - 1 >= (sp >> SA_COMPX));
        #if  SA_COMPRESSION
        int64_t sa_entry = sa_ms_byte[sp >> SA_COMPX];
        #else
        int64_t sa_entry = sa_ms_byte[sp];      // simultion
        #endif
        
        sa_entry = sa_entry << 32;

        #if  SA_COMPRESSION
        sa_entry = sa_entry + sa_ls_word[sp >> SA_COMPX];
        #else
        sa_entry = sa_entry + sa_ls_word[sp];      // simulation
        #endif
        
        sa_entry += offset;
        return sa_entry;
    }
}

void FMI_search::get_sa_entries(SMEM *smemArray, int64_t *coordArray, int32_t *coordCountArray, uint32_t count, int32_t max_occ, int tid)
{
    
    uint32_t i;
    int32_t totalCoordCount = 0;
    for(i = 0; i < count; i++)
    {
        int32_t c = 0;
        SMEM smem = smemArray[i];
        int64_t hi = smem.k + smem.s;
        int64_t step = (smem.s > max_occ) ? smem.s / max_occ : 1;
        int64_t j;
        for(j = smem.k; (j < hi) && (c < max_occ); j+=step, c++)
        {
            int64_t pos = j;
            int64_t sa_entry = get_sa_entry_compressed(pos, tid);
            coordArray[totalCoordCount + c] = sa_entry;
        }
        // coordCountArray[i] = c;
        *coordCountArray += c;
        totalCoordCount += c;
    }
}

// SA_COPMRESSION w/ PREFETCH
int64_t FMI_search::call_one_step(int64_t pos, int64_t &sa_entry, int64_t &offset)
{
    if ((pos & SA_COMPX_MASK) == 0) {        
        sa_entry = sa_ms_byte[pos >> SA_COMPX];        
        sa_entry = sa_entry << 32;        
        sa_entry = sa_entry + sa_ls_word[pos >> SA_COMPX];        
        // return sa_entry;
        return 1;
    }
    else {
        // int64_t offset = 0; 
        int64_t sp = pos;

        int64_t occ_id_pp_ = sp >> CP_SHIFT;
        int64_t y_pp_ = CP_BLOCK_SIZE - (sp & CP_MASK) - 1; 
        uint64_t *one_hot_bwt_str = cp_occ[occ_id_pp_].one_hot_bwt_str;
        uint8_t b;

        if((one_hot_bwt_str[0] >> y_pp_) & 1)
            b = 0;
        else if((one_hot_bwt_str[1] >> y_pp_) & 1)
            b = 1;
        else if((one_hot_bwt_str[2] >> y_pp_) & 1)
            b = 2;
        else if((one_hot_bwt_str[3] >> y_pp_) & 1)
            b = 3;
        else
            b = 4;
        if (b == 4) {
            sa_entry = 0;
            return 1;
        }
        
        GET_OCC(sp, b, occ_id_sp, y_sp, occ_sp, one_hot_bwt_str_c_sp, match_mask_sp);
        
        sp = count[b] + occ_sp;
        
        offset ++;
        if ((sp & SA_COMPX_MASK) == 0) {
    
            sa_entry = sa_ms_byte[sp >> SA_COMPX];        
            sa_entry = sa_entry << 32;
            sa_entry = sa_entry + sa_ls_word[sp >> SA_COMPX];
            
            sa_entry += offset;
            // return sa_entry;
            return 1;
        }
        else {
            sa_entry = sp;
            return 0;
        }
    } // else
}

void FMI_search::get_sa_entries_prefetch(SMEM *smemArray, int64_t *coordArray,
                                         int64_t *coordCountArray, int64_t count,
                                         const int32_t max_occ, int tid, int64_t &id_)
{
    
    // uint32_t i;
    int32_t totalCoordCount = 0;
    int32_t mem_lim = 0, id = 0;
    
    for(int i = 0; i < count; i++)
    {
        int32_t c = 0;
        SMEM smem = smemArray[i];
        mem_lim += smem.s;
    }

    int64_t *pos_ar = (int64_t *) _mm_malloc( mem_lim * sizeof(int64_t), 64);
    int64_t *map_ar = (int64_t *) _mm_malloc( mem_lim * sizeof(int64_t), 64);

    for(int i = 0; i < count; i++)
    {
        int32_t c = 0;
        SMEM smem = smemArray[i];
        int64_t hi = smem.k + smem.s;
        int64_t step = (smem.s > max_occ) ? smem.s / max_occ : 1;
        int64_t j;
        for(j = smem.k; (j < hi) && (c < max_occ); j+=step, c++)
        {
            int64_t pos = j;
             pos_ar[id]  = pos;
             map_ar[id++] = totalCoordCount + c;
            // int64_t sa_entry = get_sa_entry_compressed(pos, tid);
            // coordArray[totalCoordCount + c] = sa_entry;
        }
        //coordCountArray[i] = c;
        *coordCountArray += c;
        totalCoordCount += c;
    }
    
    id_ += id;
    
    const int32_t sa_batch_size = 20;
    int64_t working_set[sa_batch_size], map_pos[sa_batch_size];;
    int64_t offset[sa_batch_size] = {-1};
    
    int i = 0, j = 0;    
    while(i<id && j<sa_batch_size)
    {
        int64_t pos =  pos_ar[i];
        working_set[j] = pos;
        map_pos[j] = map_ar[i];
        offset[j] = 0;
        
        if (pos & SA_COMPX_MASK == 0) {
            _mm_prefetch(&sa_ms_byte[pos >> SA_COMPX], _MM_HINT_T0);
            _mm_prefetch(&sa_ls_word[pos >> SA_COMPX], _MM_HINT_T0);
        }
        else {
            int64_t occ_id_pp_ = pos >> CP_SHIFT;
            _mm_prefetch(&cp_occ[occ_id_pp_], _MM_HINT_T0);
        }
        i++;
        j++;
    }
        
    int lim = j, all_quit = 0;
    while (all_quit < id)
    {
        
        for (int k=0; k<lim; k++)
        {
            int64_t sp = 0, pos = 0;
            bool quit;
            if (offset[k] >= 0) {
                quit = call_one_step(working_set[k], sp, offset[k]);
            }
            else
                continue;
            
            if (quit) {
                coordArray[map_pos[k]] = sp;
                all_quit ++;
                
                if (i < id)
                {
                    pos = pos_ar[i];
                    working_set[k] = pos;
                    map_pos[k] = map_ar[i++];
                    offset[k] = 0;
                    
                    if (pos & SA_COMPX_MASK == 0) {
                        _mm_prefetch(&sa_ms_byte[pos >> SA_COMPX], _MM_HINT_T0);
                        _mm_prefetch(&sa_ls_word[pos >> SA_COMPX], _MM_HINT_T0);
                    }
                    else {
                        int64_t occ_id_pp_ = pos >> CP_SHIFT;
                        _mm_prefetch(&cp_occ[occ_id_pp_], _MM_HINT_T0);
                    }
                }
                else
                    offset[k] = -1;
            }
            else {
                working_set[k] = sp;
                if (sp & SA_COMPX_MASK == 0) {
                    _mm_prefetch(&sa_ms_byte[sp >> SA_COMPX], _MM_HINT_T0);
                    _mm_prefetch(&sa_ls_word[sp >> SA_COMPX], _MM_HINT_T0);
                }
                else {
                    int64_t occ_id_pp_ = sp >> CP_SHIFT;
                    _mm_prefetch(&cp_occ[occ_id_pp_], _MM_HINT_T0);
                }                
            }
        }
    }
    
    _mm_free(pos_ar);
    _mm_free(map_ar);
}
