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
#include "io_utils.h"       /* parallel index-array read family (fmi_pread_*, index_load_threads) */
#include "lockstep_width.h"  /* SMEM_LOCKSTEP_N (+ _MAX), runtime width + probe */
#include "fmi_seed_api.h"

#define DUMMY_CHAR 6

#define assert_not_null(x, size, cur_alloc) \
        if (x == NULL) { fprintf(stderr, "Allocation of %0.2lf GB for " #x " failed.\nCurrent Allocation = %0.2lf GB\n", size * 1.0 /(1024*1024*1024), cur_alloc * 1.0 /(1024*1024*1024)); exit(EXIT_FAILURE); }

#define \
GET_OCC(pp, c, occ_id_pp, y_pp, occ_pp, one_hot_bwt_str_c_pp, match_mask_pp) \
                int64_t occ_id_pp = pp >> CP_SHIFT; \
                int64_t y_pp = pp & CP_MASK; \
                int64_t occ_pp = cp_occ[occ_id_pp].cp_count[c]; \
                uint64_t one_hot_bwt_str_c_pp = cp_occ[occ_id_pp].one_hot_bwt_str[c]; \
                uint64_t match_mask_pp = one_hot_bwt_str_c_pp & one_hot_mask_array[y_pp]; \
                occ_pp += _mm_countbits_64(match_mask_pp);

/* (sort key, source index) pair used only inside smem_sort_rid_mn; the full
 * definition is private to FMI_search.cpp. */
typedef struct smem_key_idx SmemKeyIdx;

/* Reusable, single-threaded scratch for the SMEM rid counting sort. The caller
 * owns one per worker thread (see mem_cache in bwamem.h) and passes it in on
 * every batch (audit SEED-15), so the buffers are allocated once and grown on
 * demand instead of allocated and freed per call:
 *   - cnt:   rid count/offset array, sized to the batch's rid range
 *            (smem_sort_rid_mn and FMI_search::sortSMEMs);
 *   - tmp:   SMEM staging buffer, sized to the SMEM count -- the rid scatter
 *            target in sortSMEMs, the gather target in smem_sort_rid_mn;
 *   - pairs: (key, idx) buffer, sized to the SMEM count (smem_sort_rid_mn only).
 * Buffers are never shrunk, so per thread the peak is 48 B (tmp) + 16 B (pairs)
 * per SMEM of the largest batch, plus 8 B per rid of its rid range.
 * Zero-initialize before first use (all NULL, zero caps); the owner frees
 * cnt/tmp/pairs with _mm_free at teardown (all NULL-safe). Never share one
 * instance across threads. */
typedef struct smem_sort_scratch
{
    int64_t *cnt;     /* rid count/offset array; >= observed rid range          */
    int64_t  cntCap;  /* allocated capacity of cnt, in int64_t entries           */
    SMEM    *tmp;     /* SMEM staging (scatter or gather) buffer; >= SMEM count  */
    int64_t  tmpCap;  /* allocated capacity of tmp, in SMEM entries              */
    SmemKeyIdx *pairs;    /* smem_sort_rid_mn's (key, idx) buffer; >= SMEM count  */
    int64_t     pairsCap; /* allocated capacity of pairs, in entries             */
}SmemSortScratch;

/* Sort n SMEMs into the order mem_collect_smem needs: ascending rid, and within
 * one rid by (m << 32 | n) with exactly the permutation that FMI_search::sortSMEMs
 * (stable counting sort by rid) followed by a per-rid ks_introsort(mem_intv1)
 * produces -- ties on (rid, m, n) included. The rid scatter and the per-rid
 * sort run on 16-byte (key, index) pairs with the same introsort and the same
 * comparisons, so the 48-byte SMEMs are not shifted by the sort; each is moved
 * twice (gathered into scratch, then copied back). scratch is single-threaded,
 * caller-owned, grown on demand. */
void smem_sort_rid_mn(SMEM *a, int64_t n, SmemSortScratch &scratch);

/* The bwtSeed pass's first walk on a read (from x = 0), recorded by the SMEM
 * pass's first forward walk on it. Both start at read base 0 and step every
 * base with the same forward extension, so the bwtSeed walk is a prefix of the
 * SMEM walk (min_intv 1), which goes on while s >= 1: the SMEM walk notes the
 * step at which the bwtSeed walk stops (s < max_intv and n + 1 >= min_len), or
 * where it would end when the SMEM walk stops first, and the bwtSeed pass
 * emits that seed and starts the read at x instead of walking it again.
 * mem_collect_smem fills `rec` (one entry per read, indexed by rid) through
 * getSMEMsAllPosOneThread, which first sets x = -1 for every rid it is given,
 * and reads it back in bwtSeedStrategyAllPosOneThread_lockstep, whose read i
 * must be rid i. A direct getSMEMsOnePosOneThread_lockstep caller must set
 * those x to -1 itself. */
typedef struct bsd_first_walk
{
    int64_t k, l, s;  /* the walk's last interval, emitted as its seed when s > 0;
                       * l is stale (frozen at the unique-match step) when the
                       * SMEM walk ended in its text compare (ls_s1_compare) */
    int32_t n;        /* the seed's last base (m = 0)                          */
    int32_t x;        /* where the bwtSeed pass goes on; < 0: not recorded     */
}BsdFirstWalk;

typedef struct bsd_first_walk_hint
{
    BsdFirstWalk *rec;     /* one per read, indexed by rid                     */
    int32_t       max_intv;/* the bwtSeed pass's max_intv (>= 1, every read)    */
    int32_t       min_len; /* the bwtSeed pass's minSeedLen                     */
}BsdFirstWalkHint;

#define SAL_PFD 16

/* SMEM_LOCKSTEP_N (phase-2 SMEM lockstep width, the compile-time floor/default
 * and the enable guard) and SMEM_LOCKSTEP_N_MAX now live in lockstep_width.h,
 * included above, alongside the runtime width g_smem_lockstep_n and the startup
 * MLP probe that resolves it. */

/* Third-pass (bwtSeedStrategy) re-seeding lockstep width + on/off policy now live
 * in lockstep_width.h (included above): BWTSEED_LOCKSTEP_N (compile default + the
 * `#if BWTSEED_LOCKSTEP_N > 1` enable guard), the runtime width
 * g_bwtseed_lockstep_n, and the on/off switch g_bwtseed_lockstep. */

/* Smallest slice of an index array a parallel-load worker is given. The load is
 * memory-bandwidth bound, so below this the per-thread create/join overhead is
 * a larger share of the work than the bandwidth it unlocks. */
#define FMI_PREAD_MIN_CHUNK (8UL << 20)

/* The parallel index-array read family (fmi_pread_worker_count,
 * fmi_pread_request_size, fmi_pread_from_stream) and the load-thread policy
 * (index_load_threads) are declared in io_utils.h (included above), the shared
 * low-level IO header, so the index-element loader can reach them without
 * depending on this derived FM-index header. Their definitions live in
 * FMI_search.cpp. FMI_PREAD_MIN_CHUNK stays here with fmi_pread_worker_count's
 * definition, which is the only user of it. */

/* Detect the trailing sa_compx tag on a loaded FM-index (.bwt.2bit.64) file,
 * falling back to `default_compx` when no valid tag is present (a legacy
 * index built before this field existed).
 *
 * `fd` is the open file descriptor to pread() from; `file_size` is its size
 * on disk (e.g. from fstat); `ref_seq_len` is the reference_seq_len already
 * read from the file's header, which the implied layout size depends on.
 *
 * Aborts the process (EXIT_FAILURE) on a tail-read error or short/EOF read:
 * silently degrading to `default_compx` there would mis-size the SA-sample
 * arrays of a genuinely non-default index with no diagnostic. `default_compx`
 * is returned only for a fully-read tag that fails the layout test (a legacy,
 * pre-tag index) or a file too small to carry a tail.
 *
 * Shared by FMI_search::load_index (the disk loader, FMI_search.cpp) and
 * bwa_shm_layout_load (the shm sizing path, bwa_shm.cpp) -- both must agree
 * on the detected rate, since bwa_shm_compute sizes the shm segment the
 * loader then attaches to.
 *
 * This tail-detection heuristic is provably unambiguous and MUST NOT be
 * replaced by a dedicated magic-number/version field in the on-disk format:
 * see the implementation comment in FMI_search.cpp for the proof (the size
 * delta a false positive would require, 5*Delta == 8, has no integer
 * solution since 5 does not divide 8). Adding a version field would change
 * the on-disk format and break this branch's byte-identity guarantee with
 * stock indexes for no benefit. */
int64_t detect_sa_compx(int fd, int64_t file_size, int64_t ref_seq_len, int64_t default_compx);

/* Configure the cross-read SMEM-interval (k,s) SA-dedup: --ks-dedup CLI value
 * (off|on|auto) > env BWA3_KS_DEDUP > default 'auto'; fatal on a bad value.
 * Mirrors mem_dedup_configure. Call once from fastmap.cpp after getopt (before
 * the first alignment batch); an explicitly-empty CLI value is rejected by the
 * getopt handler, an explicitly-empty env value is fatal here. Alignment output
 * is byte-identical in every mode -- the flag only trades SA-resolve LF-walks
 * for a per-chunk (k,s) dedup, and 'auto' latches it on/off by measured net. */
void ks_dedup_configure(const char *mode_arg);

/* Read the process-global cross-read (k,s) dedup position counters (the same
 * ones BWA3_KS_DEDUP_STATS prints at exit): *total = coordinate positions staged
 * across every resolved chunk, *distinct = those actually LF-walked (the reps).
 * distinct < total iff cross-read duplicates were collapsed. The counters only
 * accumulate while BWA3_KS_DEDUP_STATS is set (they are gated off the hot
 * staging loop otherwise), so a reader must enable that env before resolving.
 * Exposed so a test can assert the dedup staging path actually executed --
 * the coordinate oracle alone cannot, since the plain path returns the same
 * coordinates. */
void ks_dedup_position_counts(uint64_t *total, uint64_t *distinct);

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
     * read; pass true only to emit it for an external consumer (e.g. bwa-mem2).
     * sa_compx defaults to 3 (the historical rate); `bwa-mem3 index -u INT`
     * threads its validated [0,6] value down to here via bwa_idx_build. */
    int build_index(bool emit_unpacked_ref = false, int sa_compx = 3);
    /* load_pac=false skips loading the 2-bit packed reference (BNS only). D3
     * --meth uses this for the SEED index: seeding needs the FM-index + bns
     * (for the seed->original remap) but never the seed pac — extension/scoring
     * runs against the ORIGINAL pac (meth_orig_pac). Saves ~1.6 GB on hg38. */
    /* n_threads: worker count for the big-array disk reads (capped internally,
     * overridable via BWA3_LOAD_THREADS). Default 1 preserves prior behavior.
     * build_seed_jump=false skips the seed jump table (seed_jump stays NULL,
     * seed_jump_depth 0, as under BWA3_BWTSEED_JUMP=0): for a loader that never
     * runs the lockstep seeding drivers, such as the fmi_seed_api facade. The
     * table is built on both the disk and the shm-attach path. */
    void load_index(bool load_pac = true, int n_threads = 1, bool build_seed_jump = true);

    /* Attach to a packed bwa-mem3 index segment from bwa_shm_attach. Sets
     * scalars and the cp_occ / sa_ms_byte / sa_ls_word pointers; the
     * destructor munmaps `base` and leaves the aliased buffers untouched.
     * Builds no seed jump table: load_index builds it after the attach. */
    void load_index_from_shm(uint8_t *base, size_t len);

    /* Densify the SA sample table on the fly: given an existing FM-index's
     * checkpoint-occ block, C[] array, and SA samples taken at stride
     * `src_compx` (shift), synthesize the samples for a DENSER stride
     * `dst_compx` (dst_compx < src_compx). For each destination-sampled BWT row
     * whose value is not already present in the source table, the value is
     * recovered by the same LF-walk the resolver uses (get_sa_entry_compressed
     * evaluated against the SOURCE stride), so every synthesized sample is
     * byte-identical to what an index physically built at `dst_compx` would
     * hold. Used by `bwa-mem3 shm -u INT` to stage a denser SA table without
     * rebuilding the on-disk index. All pointer arguments are BORROWED: they
     * are cleared before return so this object's destructor frees nothing.
     * dst_ms/dst_ls must be sized for (ref_seq_len >> dst_compx) + 1 samples;
     * src_ms/src_ls for (ref_seq_len >> src_compx) + 1. Aborts on invalid
     * strides (must satisfy 0 <= dst_compx < src_compx <= 6).
     * n_threads: the per-destination-sample LF-walks are independent and write
     * disjoint slots, so the fill loop parallelizes with no locking; <=1 runs
     * serially. get_sa_entry_compressed only reads shared index state. */
    void densify_sa_into(const CP_OCC *cp_occ_src, const int64_t count_src[5],
                         int64_t ref_seq_len, int64_t sentinel_idx,
                         const int8_t *src_ms, const uint32_t *src_ls,
                         int64_t src_compx,
                         int8_t *dst_ms, uint32_t *dst_ls, int64_t dst_compx,
                         int n_threads = 1);

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
     * the cross-read interleaving is new. A slot is refilled as soon as its
     * read finishes; finished reads are emitted in input order through a
     * per-thread parked list (smem_parked.h), so matchArray is written in the
     * same (read, smem) order as the scalar path.
     * matchArray must hold at least numReads * max_readlength SMEMs (caller-sized).
     * With bsd_hint, a walk from query position 0 at min_intv 1 records the
     * bwtSeed pass's first walk on its read (BsdFirstWalk). The caller must
     * first set bsd_hint->rec[rid].x = -1 for every rid it passes, as
     * getSMEMsAllPosOneThread does; a read whose walk records nothing keeps
     * that mark. */
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
                                          int64_t *__numTotalSmem,
                                          BsdFirstWalkHint *bsd_hint = nullptr);

    /* matchArray must hold at least numReads * max_readlength SMEMs (caller-sized).
     * bsd_hint (first round only, lockstep build) as in
     * getSMEMsOnePosOneThread_lockstep; this call sets x = -1 in every rid's
     * record itself. */
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
                                 int64_t *__numTotalSmem,
                                 BsdFirstWalkHint *bsd_hint = nullptr);


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
     * least numReads * max_readlength SMEMs. With bsd_hint (filled by the
     * SMEM pass, BsdFirstWalk), a read whose first walk was recorded emits
     * that walk's seed and starts at its x; a hint whose min_len or max_intv
     * differs from this call's is not used. rec is indexed by read index,
     * which must equal the rid the SMEM pass recorded it under, and any
     * record with x > 0 is trusted as this batch's. A seed taken from a
     * record may carry a stale l (frozen at the SMEM walk's unique-match
     * step), so it matches the scalar driver's seed in (rid, m, n, k, s). */
    int64_t bwtSeedStrategyAllPosOneThread_lockstep(uint8_t *enc_qdb,
                                                    int32_t *max_intv_array,
                                                    int32_t numReads,
                                                    const bseq1_t *seq_,
                                                    int32_t *query_cum_len_ar,
                                                    int32_t minSeedLen,
                                                    SMEM *matchArray,
                                                    int32_t max_readlength,
                                                    const BsdFirstWalkHint *bsd_hint = nullptr);

    /* No production caller: mem_collect_smem uses smem_sort_rid_mn. Kept for the
     * integration tests (smem2_test, bwt_seed_strategy_test, fmi_test); the
     * smem_sort_rid_mn unit test also pins it to a stable sort by rid. */
    void sortSMEMs(SMEM *matchArray,
                   int64_t numTotalSmem[],
                   int32_t numReads,
                   int32_t readlength,
                   int nthreads,
                   SmemSortScratch &scratch);
    // Uncompressed SA accessors -- index sa_ms_byte/sa_ls_word by raw BWT row, so
    // they read out of bounds against the SA_COMPRESSION-sized arrays. Fenced to
    // match their definitions in FMI_search.cpp; the compressed equivalents below
    // are the shipped paths. Every in-tree caller is likewise `#if !SA_COMPRESSION`.
#if !SA_COMPRESSION
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
#endif // !SA_COMPRESSION
    int64_t get_sa_entry_compressed(int64_t pos, int tid=0);
    void get_sa_entries(SMEM *smemArray,
                        int64_t *coordArray,
                        int32_t *coordCountArray,
                        uint32_t count,
                        int32_t max_occ,
                        int tid);
    /* One step of a lane's LF walk in the pipelined resolver. `offset` is the
     * lane's in/out walk length: 0 when a row is first staged, then passed back
     * unchanged with each continuation. Returns 1 with sa_entry = the resolved
     * coordinate (sample(pos) + offset, or the sentinel value below). Returns 0
     * with sa_entry = the next row to pass back in as `pos`; that row may itself
     * be sampled (the step does not read the sample, so the caller can prefetch
     * it), and the next call resolves it.
     * drop_sentinel_offset selects the coordinate reported when the LF walk
     * reaches the sentinel row: 0 = the accumulated walk offset (correct, bwa's
     * behaviour), 1 = 0 (bwa-mem2 v2.2.1's behaviour, selected by
     * --compat=bwa-mem2; see compat_target_t::sa_sentinel_drop_offset). */
    int64_t call_one_step(int64_t pos, int64_t &sa_entry, int64_t &offset,
                          int drop_sentinel_offset = 0);
    /* trust_s1_tag: smemArray comes from bwa-mem3's seeding path
     * (mem_collect_smem), so an s == 1 SMEM carrying the lockstep driver's
     * text-position tag in info is written from it without an LF walk
     * (FMI_search.cpp, ls_s1_sa_tag); used only at drop_sentinel_offset 0.
     * External callers (the fmi_seed_sa_prefetch facade) leave it false. */
    void get_sa_entries_prefetch(SMEM *smemArray, int64_t *coordArray,
                                 int64_t *coordCountArray, int64_t count,
                                 const int32_t max_occ, int tid, int64_t &id_,
                                 int drop_sentinel_offset = 0,
                                 bool trust_s1_tag = false);
    /* One read's SMEMs a[0..n), sorted by (m, n) as mem_collect_smem sorts
     * them: an untagged s == 1 SMEM that lies inside a tagged one gives it the
     * same tag (FMI_search.cpp, ls_s1_sa_tag), so get_sa_entries_prefetch
     * writes its coordinate without an LF walk too. */
    static void s1_tag_contained(SMEM *a, int64_t n);
    /* s1_tag_contained over each read's run of a[0..n), an array sorted by rid
     * and then (m, n), as mem_collect_smem sorts it (and calls this). */
    static inline void s1_tag_contained_reads(SMEM *a, int64_t n)
    {
        for (int64_t b0 = 0, b1; b0 < n; b0 = b1) {
            for (b1 = b0 + 1; b1 < n && a[b1].rid == a[b0].rid; b1++) {}
            if (b1 - b0 > 1)
                s1_tag_contained(&a[b0], b1 - b0);
        }
    }
    /* True when the lockstep SMEM driver compares unique matches with the
     * reference text and so can tag SMEMs (the index's pac is loaded and
     * matches the BWT text length). When false no SMEM carries a tag. */
    bool s1_text_available() const;

    int64_t reference_seq_len;
    int64_t sentinel_index;
    int64_t sa_compx;        // SA sample shift; sample rate 1/(1<<sa_compx). Loaded per-index.
    int64_t sa_compx_mask;   // (1<<sa_compx)-1
private:
        char file_name[PATH_MAX];
        int64_t count[5];
        uint32_t *sa_ls_word;
        int8_t *sa_ms_byte;
        CP_OCC *cp_occ;

        /* If non-NULL, cp_occ / sa_ms_byte / sa_ls_word point into a shared
         * memory mapping owned by the shm segment rather than being
         * _mm_malloc'd. shm_len is the byte length of the mapping so the
         * destructor can munmap it; the destructor also skips _mm_free on
         * shm-backed buffers. */
        uint8_t *shm_base;
        size_t   shm_len;

        /* Seed jump table: the bidirectional interval (k, l, s) of every
         * string of length 1..seed_jump_depth, level L at offset
         * seed_jump_level_off(L), indexed within a level by the string's
         * 2-bit code, first base most significant. Built once per index load by
         * build_seed_jump_table after cp_occ, count[] and sentinel_index are in
         * (load_index, after the disk load or the shm attach); NULL (seed_jump_depth 0) when
         * BWA3_BWTSEED_JUMP=0, the loader passed build_seed_jump=false, or
         * the table could not be allocated.
         * Level 1 is bsd_init_slot's 1-base interval and every child is the
         * walk's own forward step (backwardExt on the k<->l-swapped parent with
         * 3 - c, swapped back), so an entry is exactly the interval a bwtSeed
         * walk reaches by stepping, s == 0 strings included. A walk cannot
         * leave before length minSeedLen except at an N or the read end, so
         * bwtSeedStrategyAllPosOneThread_lockstep jumps from its first base
         * straight to length min(depth, minSeedLen - 1) when the read has
         * that many ACGT bases there. The SMEM lockstep driver's forward walks
         * step the same way, so a walk's first visit takes its steps to depths
         * 2..min(depth, ACGT bases left) from the table, each with the forward
         * step's own s < min_intv exit (ls_advance_forward_table). Owned
         * (_mm_malloc), also on the shm path. */
        struct SeedJumpEntry { int64_t k, l, s; };
        SeedJumpEntry *seed_jump;
        int32_t        seed_jump_depth;
        static int64_t seed_jump_level_off(int L) { return (((int64_t)1 << (2 * L)) - 4) / 3; }
        static inline void seed_jump_prefetch(const SeedJumpEntry *e);
        void build_seed_jump_table(int n_threads);
        void seed_jump_children(int L, int64_t first, int64_t n);
        static void *seed_jump_worker(void *arg);

        /* Defined inline so all hot callers (getSMEMs* and ls_advance_*)
         * fully absorb the body and pay no struct-by-value pass / return
         * cost. backwardExt is called ~10^10 times per 5M-pair WGS run, so
         * the SysV-ABI hidden-pointer dance for the 24-byte SMEM struct
         * dominates self-time on gcc 12+ (issue #87): a single output store
         * `vmovdqu %ymm0, (%r8)` writing the return SMEM hits 42% of the
         * function's CPU samples on c7a (Zen 4) with gcc-14, with the
         * matching argument load close behind. always_inline removes the
         * call boundary and recovers the full gcc-11 baseline (and beats
         * it by ~4% wall-clock on c7a with a 5M-read WGS slice, shm-warmed). */
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

        /* Population count of (*word & *mask) for one 64-bit BWT word.
         * arm64 has no scalar popcount: __builtin_popcountll lowers to
         * fmov d,x; cnt.8b; addv.8b; fmov x,d -- two GPR<->SIMD moves on the
         * k' dependency chain. Loading the word and the mask straight into d
         * registers (both are memory operands here) leaves one move, the
         * result. Bit-identical. Other targets keep the scalar builtin. */
        static inline int64_t occ_popcount64(const uint64_t *word, const uint64_t *mask)
        {
#if defined(__ARM_NEON) || defined(__aarch64__) || defined(APPLE_SILICON)
            const uint8x8_t bits = vreinterpret_u8_u64(vand_u64(vld1_u64(word), vld1_u64(mask)));
            return (int64_t) vaddv_u8(vcnt_u8(bits));
#else
            return (int64_t) __builtin_popcountll(*word & *mask);
#endif
        }

        /* K-only backward extension for the pure-backward SMEM phase
         * (ls_advance_backward_step). Computes ONLY smem.k and smem.s, leaving
         * smem.l untouched: in the backward phase the reverse-complement interval
         * start `l` is dead — it is consumed only to compute the next step's `l`,
         * and every downstream reader of an emitted seed's `l` is redundant by the
         * (rid,m,n) => (k,l,s) invariant. Those readers were the smem_dedup
         * comparator (smem_dedup_inplace) and the scalar/lockstep field comparator
         * in the parity harness (test/smem_lockstep_parity_test.cpp); both have
         * been updated to drop `l`, so an emitted seed's identity is (rid,m,n,k,s).
         * NOTE this means the l this phase emits is intentionally stale (frozen at
         * its forward-phase value); so is the l of every seed the forward phase's
         * unique-match text compare finishes (ls_s1_compare in FMI_search.cpp
         * never steps l). It must never be resurrected as an observable field
         * without also restoring the l-chain here and in that compare. Dropping the l-chain
         * removes the sentinel test, the l cumulation, and (on x86) six of eight
         * popcounts — same two cache lines read, far fewer instructions.
         * Byte-identical in observable (SAM) output.
         *
         * s==1 fast path: when the interval is a single suffix (common once the
         * match is unique), occ(a, sp+1) - occ(a, sp) is exactly the one BWT bit
         * at sp, so s' is a single bit test on the sp block — no ep block, no
         * second popcount. one_hot_mask_array is MSB-first (mask[i] = top i bits),
         * so BWT position `pos` lives at bit (CP_MASK - pos). On the s'->0 exit the
         * new k is left stale (unset): this is safe ONLY because every caller
         * discards a seed with s < min_intv while emitting the OLD smem, and every
         * caller uses min_intv >= 1 (enforced by an xassert in ls_init_slot), so an
         * s'==0 result is always discarded and its k is never read. Under that
         * precondition the remaining popcount is skipped too. */
        __attribute__((always_inline)) inline
        SMEM backwardExt_konly(SMEM smem, uint8_t a) const
        {
            const int64_t sp = (int64_t)smem.k;
            if (smem.s == 1) {
                const CP_OCC &blk = cp_occ[sp >> CP_SHIFT];
                const int pos = (int)(sp & CP_MASK);
                const uint64_t oh = blk.one_hot_bwt_str[a];
                if ((oh >> (CP_MASK - pos)) & 1ULL) {  /* BWT[sp] == a  => s' = 1 */
                    smem.k = count[a] + blk.cp_count[a] +
                             occ_popcount64(&blk.one_hot_bwt_str[a], &one_hot_mask_array[pos]);
                    smem.s = 1;
                } else {                               /* s' = 0; new k is dead */
                    smem.s = 0;
                }
                return smem;
            }
            const int64_t ep = sp + (int64_t)smem.s;
            const CP_OCC &blk_sp = cp_occ[sp >> CP_SHIFT];
            const CP_OCC &blk_ep = cp_occ[ep >> CP_SHIFT];
            const int64_t occ_s = blk_sp.cp_count[a] +
                occ_popcount64(&blk_sp.one_hot_bwt_str[a], &one_hot_mask_array[(sp & CP_MASK)]);
            const int64_t occ_e = blk_ep.cp_count[a] +
                occ_popcount64(&blk_ep.one_hot_bwt_str[a], &one_hot_mask_array[(ep & CP_MASK)]);
            smem.k = count[a] + occ_s;
            smem.s = occ_e - occ_s;
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
                      const uint8_t *enc_qdb,
                      BsdFirstWalkHint *bsd_hint);
    void ls_prefetch_cp_occ(const BatchSlot *s);
    inline void ls_bsd_note(BatchSlot *s, int64_t k, int64_t l, int64_t sz);
    void ls_bsd_resolve(BatchSlot *s, const uint8_t *enc_qdb, int32_t p);
    void ls_advance_forward_step(BatchSlot *s, const uint8_t *enc_qdb,
                                 bool s1_text);
    inline void ls_prefetch_sa_row(int64_t p) const;
    void ls_s1_lf_step(BatchSlot *s, const uint8_t *pac, int64_t l_pac);
    void ls_s1_compare(BatchSlot *s, const uint8_t *enc_qdb,
                       const uint8_t *pac, int64_t l_pac);
    void ls_advance_forward_table(BatchSlot *s, const uint8_t *enc_qdb,
                                  bool s1_text);
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
                       const uint8_t *enc_qdb,
                       const SeedJumpEntry *jlev, int32_t jd,
                       int32_t minSeedLen,
                       const BsdFirstWalkHint *bsd_hint);
    void bsd_prefetch_cp_occ(const BwtSeedSlot *s);
    void bsd_plan_jump(BwtSeedSlot *s, const uint8_t *enc_qdb,
                       const SeedJumpEntry *jlev, int32_t jd);
    void bsd_advance_step(BwtSeedSlot *s,
                          const uint8_t *enc_qdb,
                          int32_t minSeedLen,
                          const SeedJumpEntry *jlev, int32_t jd);
};

#endif
