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
#include "bwa_hugepages.h"  /* bwamem_hugetlb_buf_t (--huge-pages occurrence table) */

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

/* Structural index checks shared by disk load, staging and re-sa. They compare
 * header lengths and file sizes, not sequence content or provenance: same-length
 * mixes can pass. A tagged file must also carry its actual SA-rate shift.
 * `fd` names the already-open .bwt.2bit.64. False reports a diagnostic.
 * from_shm validates only segment-internal lengths, not freshness against disk. */
bool fmi_bwt_size_ok(const char *func, const char *cp_path, int64_t file_size, int64_t ref_seq_len,
                     int64_t sa_compx, int fd);
bool fmi_bwt_matches_ann(const char *func, const char *prefix, int64_t ref_seq_len, int64_t l_pac,
                         bool from_shm = false);

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
    /* The high part of each SA sample (bits 32.. of a 2L-space position) is stored packed,
     * (1 << sa_hi_log2) bits per sample, a whole number of samples per 64-bit word (never
     * straddling one). On disk it is still one byte per sample; GRCh38-sized references need
     * one bit (positions < 2^33), which cuts the array 8-fold. */
    int64_t sa_hi_size_bytes()       const { return sa_hi_bytes_for(reference_seq_len, sa_compx); }
    int64_t sa_ls_word_size_bytes() const { return sa_sample_count() * (int64_t)sizeof(uint32_t); }

    /* The same sizes from the scalars alone, for a check made before any
     * FMI_search exists (bwa_shm's probe of a 1 GB stage). */
    static int64_t cp_occ_bytes_for(int64_t ref_seq_len);
    static int64_t sa_count_for(int64_t ref_seq_len, int64_t sa_compx);

    /* The in-memory layout of the occurrence checkpoints, chosen before
     * load_index() by the caller (main_mem: --compact-occ or BWA3_COMPACT_OCC;
     * the constructor always starts in the default layout):
     *   - default: the file's own 64-byte CP_OCC records, read as they are;
     *   - compact: 32-byte CP_OCC_PLANES records converted at load (half the
     *     memory: about 3 GB less on hg38), see fmi_seed_api.h.
     * Alignment records are byte-identical in both for the same SA sampling
     * (only @PG CL: differs, recording --compact-occ). Under --compat=bwa-mem2 a
     * run attached to a denser staged SA table can report other coordinates
     * than one reading the disk table: see docs/src/cli/mem.md, --compact-occ.
     * A staged shm segment holds one layout (`shm --compact-occ` stages the
     * compact one). By default
     * load_index() attaches only a segment of this layout and loads from disk
     * otherwise, so a caller that never calls prefer_staged_layout() (the seed
     * facade, the tests) always gets the layout it set.
     * set_compact_occ() after an index is loaded is refused (exit). */
    void set_compact_occ(bool on);
    bool compact_occ() const { return cp_compact; }
    /* main_mem only: load_index() attaches an index staged in shared memory
     * whatever its checkpoint layout, which then replaces set_compact_occ()'s
     * choice (that choice applies to a disk load). `requested`: the user chose
     * the layout explicitly (--compact-occ or BWA3_COMPACT_OCC), so a staged
     * index of the other layout is reported with a [W::load_index] line; with
     * no request, attaching a compact one is noted with an [M::load_index] line.
     * Before load_index(). */
    void prefer_staged_layout(bool requested);
    /* main_mem's (and the tests') reading of BWA3_COMPACT_OCC: on for a value
     * starting with '1'; a non-empty value starting with neither '0' nor '1' is
     * reported once and ignored. *set (if given): the value starts with '0' or
     * '1', so it chose a layout. The library never reads it on its own. */
    static bool compact_occ_env(bool *set = NULL);

    /* The compact layout's checkpoint conversion, for `shm --compact-occ`
     * (bwa_shm.cpp), which stages what a compact load_index() holds.
     * cp_compact_form: the span size and count form a compact load picks for
     * the index whose RAW (on-disk, unadjusted) cumulative counts are
     * raw_count. cp_compact_convert: converts the checkpoint records of the
     * .bwt.2bit.64 open at fd (`path`, for messages) into dst ((ref_seq_len >>
     * CP_SHIFT) + 1 records) with up to nthreads workers, and when form.wide
     * reads the span bases into a new int64_t[*n_spans][4] at *base (free()
     * it; NULL otherwise). raw_count: the file's C[] array, before
     * load_index's +1, as cp_load_from_file takes it. The conversion and checks
     * are load_index()'s own: false (nothing allocated) where load_index()
     * would refuse the file, a record that is not one-hot or out of range, or
     * whose counts are not the running totals, or a base-less row other than
     * sentinel_index's. */
    struct CpCompactForm { int span_shift; bool wide; };
    static CpCompactForm cp_compact_form(const int64_t raw_count[5]);
    static bool cp_compact_convert(int fd, const char *path, int64_t ref_seq_len,
                                   const int64_t raw_count[5], int64_t sentinel_index,
                                   CpCompactForm form,
                                   CP_OCC_PLANES *dst, int64_t (**base)[4], int64_t *n_spans,
                                   int nthreads);
    /* The span sizes (log2 rows) the compact layout takes: at least one block,
     * and at most 2^32 rows, so a span-relative count fits 32 bits. */
    static bool cp_span_shift_ok(int64_t span_shift) { return span_shift >= CP_SHIFT && span_shift <= 32; }
    /* Tests only: the span size (cp_span_shift_ok; 0 = the default, 32) and
     * whether to rebase the counts even on a narrow reference, applied by the
     * next compact load or cp_compact_form. Process-wide. */
    static void set_cp_layout_for_test(int span_shift, bool force_wide);

    /* log2 of the bits per packed SA-sample high part for a BWT of ref_seq_len rows
     * (positions run up to ref_seq_len - 1): 1, 2, 4 or 8 bits -> 0..3. */
    static int sa_hi_log2_for(int64_t ref_seq_len) {
        uint64_t hi = (uint64_t)(ref_seq_len > 0 ? ref_seq_len - 1 : 0) >> 32;
        int bits = 0;
        while (hi) { bits++; hi >>= 1; }
        return bits <= 1 ? 0 : bits <= 2 ? 1 : bits <= 4 ? 2 : 3;
    }
    static int64_t sa_hi_words(int64_t n, int log2) { return ((n << log2) + 63) >> 6; }
    static uint64_t sa_hi_mask_for(int log2) { return (1ull << (1u << log2)) - 1; }
    /* Bytes of the packed high parts (the FMI_SA_MS section) of a BWT of ref_seq_len rows
     * sampled at 1 / 2^sa_compx: the one size the loader, the stager and the attach check use. */
    static int64_t sa_hi_bytes_for(int64_t ref_seq_len, int64_t sa_compx) {
        return sa_hi_words(sa_count_for(ref_seq_len, sa_compx), sa_hi_log2_for(ref_seq_len))
               * (int64_t)sizeof(uint64_t);
    }
    /* Pack n one-byte high parts (src[0..n)) into dst, starting at sample index `first`, which
     * must be a multiple of 64 >> log2 (asserted) so the call owns whole words. Each word the
     * call touches is assigned, not OR-ed, with unused high fields of a partial last word zero,
     * so dst need not be zeroed. Returns false if a byte does not fit in (1 << log2) bits (a
     * corrupt index). Parallel over words when nthreads > 1. */
    static bool sa_hi_pack(const int8_t *src, int64_t first, int64_t n, int log2, uint64_t *dst,
                           int nthreads);
    /* Samples each worker of sa_hi_pread_pack reads and packs at a time (its private
     * buffer). Default 2^20, so a piece is packed from cache. The setter requires a positive
     * multiple of 64 (whole words at every width); tests lower it to cross many piece
     * boundaries on a small index. */
    static int64_t sa_hi_read_chunk();
    static void set_sa_hi_read_chunk(int64_t n);
    /* Read n one-byte high parts at offset `off` of the file open as `fd` (named `path`, for
     * the workers' private descriptors) and pack them into dst, assigning all of its words.
     * Up to nthreads workers (the index load's worker count for n bytes) each read a
     * whole-word range a piece at a time and pack it at once. Returns 0, SA_HI_READ_OOM, or
     * SA_HI_READ_RANGE if a byte does not fit (a corrupt index); an I/O error or a short
     * file exits, as the rest of the index load does. Does not move fd's file offset. */
    enum { SA_HI_READ_OOM = -1, SA_HI_READ_RANGE = -2 };
    static int sa_hi_pread_pack(int fd, const char *path, off_t off, int64_t n, int log2,
                                uint64_t *dst, int nthreads);

    /* Checks the packed FMI_SCALARS at `scalars` (BWA_SHM_FMI_SCALARS_BYTES
     * bytes) the way load_index_from_shm does. Returns true when they pass,
     * with sizes[0..2] set to the FMI_CP_OCC, FMI_SA_MS and FMI_SA_LS section
     * sizes the load then requires; otherwise false, with what is wrong
     * written to `why` (if not NULL). */
    static bool shm_scalars_ok(const uint8_t *scalars, int64_t sizes[3],
                               char *why, size_t why_len);

    /* Read-only data accessors. Used by the round-trip test to byte-compare
     * the packed segment payload against the in-memory loader's buffers.
     * cp_occ_data(): CP_OCC[] by default, CP_OCC_PLANES[] when compact. */
    const void     *cp_occ_data()     const { return cp_compact ? (const void *)cp_planes : (const void *)cp_occ; }
    /* Compact layout only: int64_t[4] per 2^cp_span_shift rows when the counts
     * are span-relative (cp_wide), else NULL. */
    const int64_t  *cp_count_base_data() const { return cp_compact && cp_wide ? &cp_base[0][0] : NULL; }
    const uint64_t *sa_hi_data()      const { return sa_hi; }

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
     * table is built on both the disk and the shm-attach path.
     * huge_pages: BWAMEM_HUGETLB_ON (main_mem's --huge-pages) maps the
     * occurrence table on 1 GB pages where it can, falling back to the default
     * allocation when it cannot; the SA samples always use the default
     * allocation. Disk path only (an shm attach prints that it has no effect).
     * Production callers pass OFF or ON; the other values are unit-test seams.
     * See bwa_hugepages.h. */
    void load_index(bool load_pac = true, int n_threads = 1, bool build_seed_jump = true,
                    bwamem_hugetlb_mode_t huge_pages = BWAMEM_HUGETLB_OFF);

    /* Attach to a packed bwa-mem3 index segment from bwa_shm_attach: `len` is
     * the packed index size and `map_len` the length of the mapping. Sets
     * scalars and the cp_occ / sa_hi / sa_ls_word pointers; the
     * destructor unmaps `map_len` bytes at `base` and leaves the aliased
     * buffers untouched. Builds no seed jump table: load_index builds it
     * after the attach. */
    void load_index_from_shm(uint8_t *base, size_t len, size_t map_len);

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
     * The high parts go to dst_ms (one byte per sample, the disk format) and/or
     * dst_hi (packed as sa_hi, sa_hi_words(count, sa_hi_log2_for(ref_seq_len))
     * words, written directly with no byte temporary); either may be NULL.
     * dst_ms/dst_ls must be sized for (ref_seq_len >> dst_compx) + 1 samples;
     * src_ms/src_ls for (ref_seq_len >> src_compx) + 1. Aborts on invalid
     * strides (must satisfy 0 <= dst_compx < src_compx <= 6). Returns false,
     * with a message on stderr, if memory for the packed source runs out or a
     * source or (with dst_hi) densified high part does not fit its width (a
     * corrupt index), so a caller holding a segment can clean up.
     * n_threads: the per-destination-sample LF-walks are independent and write
     * disjoint slots, so the fill loop parallelizes with no locking; <=1 runs
     * serially. get_sa_entry_compressed only reads shared index state.
     * cp_occ_src is the file's own checkpoint layout: this object is set to the
     * default layout for the walk, whatever BWA3_COMPACT_OCC says. */
    bool densify_sa_into(const CP_OCC *cp_occ_src, const int64_t count_src[5],
                         int64_t ref_seq_len, int64_t sentinel_idx,
                         const int8_t *src_ms, const uint32_t *src_ls,
                         int64_t src_compx,
                         int8_t *dst_ms, uint64_t *dst_hi, uint32_t *dst_ls,
                         int64_t dst_compx, int n_threads = 1);

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
    // Uncompressed SA accessors -- index sa_hi_at()/sa_ls_word by raw BWT row, so
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
        uint64_t *sa_hi;      /* packed SA-sample high parts; see sa_hi_size_bytes() */
        int sa_hi_log2;
        uint64_t sa_hi_mask;
        void set_sa_hi_geometry() {
            sa_hi_log2 = sa_hi_log2_for(reference_seq_len);
            sa_hi_mask = sa_hi_mask_for(sa_hi_log2);
        }
        /* The high part of SA sample i, and the address of the word holding it (for prefetch). */
        inline int64_t sa_hi_at(int64_t i) const {
            const uint64_t b = (uint64_t)i << sa_hi_log2;
            return (int64_t)((sa_hi[b >> 6] >> (b & 63)) & sa_hi_mask);
        }
        inline const char *sa_hi_addr(int64_t i) const {
            return (const char *)&sa_hi[((uint64_t)i << sa_hi_log2) >> 6];
        }
        /* The occurrence checkpoints, in one of two layouts (set_compact_occ):
         * cp_occ holds the file's CP_OCC records (default; NULL when compact),
         * cp_planes the CP_OCC_PLANES form (compact; NULL otherwise). */
        bool cp_compact;
        CP_OCC *cp_occ;
        CP_OCC_PLANES *cp_planes;
        /* Compact layout only. cp_wide: some base occurs 2^32 times or more, so
         * a 32-bit count cannot hold it. Then each record's counts are relative
         * to the full counts at the first block of its 2^cp_span_shift-row span,
         * cp_base[span] (cp_n_spans entries; cp_span_shift is 32 outside the
         * tests). cp_base is read (and allocated) only when cp_wide; NULL
         * otherwise. */
        bool cp_wide;
        int cp_span_shift;
        int64_t cp_n_spans;
        int64_t (*cp_base)[4];
        /* Compact load only: the one row below reference_seq_len whose record
         * sets no plane bit, checked against sentinel_index once it is read. */
        int64_t cp_empty_row;
        /* Compact layout: the rows of the sentinel's block after the sentinel,
         * the only rows whose A count the planes over-count (the sentinel row
         * reads as A), tested as (uint64_t)(row - cp_sent_after) < cp_sent_len.
         * Set by set_cp_sentinel() wherever sentinel_index is. */
        int64_t cp_sent_after;
        uint64_t cp_sent_len;
        void set_cp_sentinel() {
            cp_sent_after = sentinel_index + 1;
            cp_sent_len = (uint64_t)(((sentinel_index >> CP_SHIFT) + 1) * CP_BLOCK_SIZE - 1 - sentinel_index);
        }

        /* The compact layout's file-to-memory checkpoint transform (see
         * CP_OCC_PLANES), used by load_index and cp_compact_convert. cp_is_wide:
         * from the RAW (unadjusted) cumulative counts, whether some base occurs
         * 2^32 times or more. */
        static bool cp_is_wide(const int64_t raw_count[5]);
        /* n file records at file offset rec0_off of `path` (open as shared_fd) ->
         * dst, read and converted by up to nthreads workers, each on its own range
         * and descriptor. base[span] holds the full counts of the span's first
         * record (span = block >> (span_shift - CP_SHIFT)) and is read only when
         * wide. False if a record is not one-hot: a row with two bases, or a row
         * below ref_seq_len with none other than exactly one (the sentinel, whose
         * row is stored in *empty_row), or a count that does not fit; or if the
         * counts are not the running totals: the first record's are not zero,
         * another's are not the previous record's plus its one-hot popcounts, or
         * the last record's plus its popcounts are not the base totals that
         * raw_count (the file's C[] array, before load_index's +1) gives. */
        static bool cp_load_from_file(int shared_fd, const char *path, off_t rec0_off, int64_t n,
                                      int64_t ref_seq_len, const int64_t raw_count[5], CP_OCC_PLANES *dst,
                                      const int64_t (*base)[4], int span_shift, bool wide,
                                      int nthreads, int64_t *empty_row);
        /* Read the span bases (the file records at blocks span << (span_shift - CP_SHIFT))
         * from the .bwt.2bit.64 at `fd` into a new int64_t[n_spans][4]. */
        static int64_t (*cp_read_span_bases(int fd, const char *path, int64_t ref_seq_len,
                                            int span_shift, int64_t *n_spans))[4];
        /* The conversion load_index() and cp_compact_convert share: the span
         * bases when form.wide (*base, else NULL) and the n records at rec0_off,
         * converted into dst (raw_count and *empty_row as cp_load_from_file
         * takes and sets them). On false
         * *base is freed and NULL. cp_sentinel_row_ok is the check that follows
         * once sentinel_index is known. */
        static bool cp_convert_records(int fd, const char *path, off_t rec0_off, int64_t n,
                                       int64_t ref_seq_len, const int64_t raw_count[5],
                                       CpCompactForm form, CP_OCC_PLANES *dst,
                                       int64_t (**base)[4], int64_t *n_spans, int nthreads,
                                       int64_t *empty_row);
        static bool cp_sentinel_row_ok(int64_t empty_row, int64_t sentinel_index) {
            return empty_row == sentinel_index;
        }

        /* If non-NULL, cp_occ / sa_hi / sa_ls_word point into a shared
         * memory mapping owned by the shm segment rather than being
         * _mm_malloc'd. shm_len is the byte length of the mapping so the
         * destructor can munmap it; the destructor also skips _mm_free on
         * shm-backed buffers. */
        uint8_t *shm_base;
        size_t   shm_len;

        /* --huge-pages: when mapped, cp_occ is cp_occ_huge.base and the
         * destructor unmaps it instead of calling _mm_free. Unmapped (base
         * NULL) on every other path. The suffix-array samples are always on the
         * default allocation. */
        bwamem_hugetlb_buf_t cp_occ_huge;

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
        /* prefer_staged_layout(): attach a staged index of either layout, and
         * whether the layout set_compact_occ() chose was the user's request. */
        bool           cp_prefer_staged;
        bool           cp_layout_requested;
        static int64_t seed_jump_level_off(int L) { return (((int64_t)1 << (2 * L)) - 4) / 3; }
        static inline void seed_jump_prefetch(const SeedJumpEntry *e);
        void build_seed_jump_table(int n_threads);
        void seed_jump_children(int L, int64_t first, int64_t n);
        template <bool PL> void seed_jump_children_t(int L, int64_t first, int64_t n);
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
         * it by ~4% wall-clock on c7a with a 5M-read WGS slice, shm-warmed).
         * The default (file) layout's step; backwardExt_planes is the compact
         * layout's, and backwardExt_t<PL> / backwardExt pick one. */
        __attribute__((always_inline)) inline
        SMEM backwardExt_file(SMEM smem, uint8_t a) const
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
        SMEM backwardExt_konly_file(SMEM smem, uint8_t a) const
        {
            return backwardExt_konly_cp(cp_occ, count[a], smem, a);
        }

        /* backwardExt_konly on a caller-held cp_occ pointer and count[a], so the
         * general loops of ls_advance_backward_step keep both in registers across
         * their prev[] stores. Same arithmetic and same contract: k is left stale
         * on the s'->0 exit, so the caller must discard any s' < min_intv. */
        __attribute__((always_inline)) static inline
        SMEM backwardExt_konly_cp(const CP_OCC *cpo, int64_t count_a, SMEM smem, uint8_t a)
        {
            const int64_t sp = (int64_t)smem.k;
            if (smem.s == 1) {
                const CP_OCC &blk = cpo[sp >> CP_SHIFT];
                const int pos = (int)(sp & CP_MASK);
                const uint64_t oh = blk.one_hot_bwt_str[a];
                if ((oh >> (CP_MASK - pos)) & 1ULL) {  /* BWT[sp] == a  => s' = 1 */
                    smem.k = count_a + blk.cp_count[a] +
                             occ_popcount64(&blk.one_hot_bwt_str[a], &one_hot_mask_array[pos]);
                    smem.s = 1;
                } else {                               /* s' = 0; new k is dead */
                    smem.s = 0;
                }
                return smem;
            }
            const int64_t ep = sp + (int64_t)smem.s;
            const CP_OCC &blk_sp = cpo[sp >> CP_SHIFT];
            const CP_OCC &blk_ep = cpo[ep >> CP_SHIFT];
            const int64_t occ_s = blk_sp.cp_count[a] +
                occ_popcount64(&blk_sp.one_hot_bwt_str[a], &one_hot_mask_array[(sp & CP_MASK)]);
            const int64_t occ_e = blk_ep.cp_count[a] +
                occ_popcount64(&blk_ep.one_hot_bwt_str[a], &one_hot_mask_array[(ep & CP_MASK)]);
            smem.k = count_a + occ_s;
            smem.s = occ_e - occ_s;
            return smem;
        }


        /* ---- The compact layout's primitives (CP_OCC_PLANES, cp_planes) ---- */
        /* The base code (0..3) at `row` from its block's planes. The sentinel
         * row reads as A (0): callers test row == sentinel_index first. */
        static inline int cp_base_at(const CP_OCC_PLANES &b, int64_t row)
        {
            const int y = CP_MASK - (int)(row & CP_MASK);
            return (int)(((b.hi >> y) & 1) << 1 | ((b.lo >> y) & 1));
        }

        /* cp_occ_t / backwardExt_konly_ab's A for a base C, G or T (1..3) given at run
         * time in `a`. Never base A (0): this form skips the sentinel-row correction
         * that only A counts need, so an A count through it would be off by one past
         * the sentinel. cp_step_dispatch is its only user and sends a == 0 to the
         * A = 0 copy. Build with -DBWA_MEM3_DEBUG_CP_CGT_RUNTIME to check a at run
         * time (off by default; release code is unchanged). */
        static constexpr int CP_CGT_RUNTIME = -1;

        /* Everything an occurrence count reads besides the block, copied into
         * locals by the hot loops (cp_view()): read through `this`, clang reloads
         * each field after every store the loop makes, since it cannot tell the
         * stores do not alias the object. sent_after/sent_len: the rows of the
         * sentinel's block after the sentinel, the only rows whose A count the
         * planes over-count (the sentinel reads as A), as one unsigned range
         * test: (uint64_t)(pp - sent_after) < sent_len. */
        struct CpView {
            const CP_OCC_PLANES *cp;
            const int64_t (*base)[4];
            int64_t sent_after;
            uint64_t sent_len;
            int span_shift;
            bool wide;
        };
        inline CpView cp_view() const
        {
            CpView v;
            v.cp = cp_planes;
            v.base = cp_base;
            v.sent_after = cp_sent_after;
            v.sent_len = cp_sent_len;
            v.span_shift = cp_span_shift;
            v.wide = cp_wide;
            return v;
        }

        /* The plane masks of base c: all ones where c has that bit clear, so
         * (hi ^ hx) & (lo ^ lx) marks the rows holding c. */
        static inline uint64_t cp_hx(int c) { return (uint64_t)((c >> 1) & 1) - 1; }
        static inline uint64_t cp_lx(int c) { return (uint64_t)(c & 1) - 1; }

        /* Rows of `b`'s block before `pp` holding the base whose plane masks are
         * hx/lx. On arm64 the planes and the row mask are loaded straight into d
         * registers and only the count crosses back: arm64 has no scalar
         * popcount, so a word built in general registers would cross twice. */
        static inline int64_t cp_rows(const CP_OCC_PLANES &b, int64_t pp, uint64_t hx, uint64_t lx)
        {
#if defined(__ARM_NEON) || defined(__aarch64__) || defined(APPLE_SILICON)
            const uint64x1_t x = vand_u64(
                vand_u64(veor_u64(vld1_u64(&b.hi), vcreate_u64(hx)),
                         veor_u64(vld1_u64(&b.lo), vcreate_u64(lx))),
                vld1_u64(&one_hot_mask_array[pp & CP_MASK]));
            return (int64_t) vaddv_u8(vcnt_u8(vreinterpret_u8_u64(x)));
#else
            return (int64_t) __builtin_popcountll((b.hi ^ hx) & (b.lo ^ lx) & one_hot_mask_array[pp & CP_MASK]);
#endif
        }

#if defined(__ARM_NEON) || defined(__aarch64__) || defined(APPLE_SILICON)
        /* All four occurrence counts of bwt[0, row), `b` being row's block, in one vector, for
         * the narrow compact layout (!cp_wide: 32 bits hold every count) and before the sentinel
         * correction (cp_occ4_neon_sent): lanes [A, C] are ~hi & [~lo, lo] and [G, T] are
         * hi & [~lo, lo] under the row mask; cnt counts their bytes and three pairwise adds fold
         * each lane's eight bytes onto the checkpoint's 32-bit counts. A scalar popcount would
         * cross to the vector unit and back per base. */
        __attribute__((always_inline)) static inline uint32x4_t cp_occ4_neon(const CP_OCC_PLANES &b, int64_t row)
        {
            const uint64x2_t flip = vcombine_u64(vcreate_u64(~0ULL), vcreate_u64(0));
            const uint64x2_t pl = vld1q_u64(&b.hi);   /* [hi, lo] */
            const uint64x2_t m  = vld1q_dup_u64(&one_hot_mask_array[row & CP_MASK]);
            const uint64x2_t hi = vdupq_laneq_u64(pl, 0);
            const uint64x2_t lo = veorq_u64(vdupq_laneq_u64(pl, 1), flip);
            const uint64x2_t ac = vandq_u64(vbicq_u64(m, hi), lo);
            const uint64x2_t gt = vandq_u64(vandq_u64(m, hi), lo);
            const uint8x16_t bytes = vpaddq_u8(vcntq_u8(vreinterpretq_u8_u64(ac)),
                                               vcntq_u8(vreinterpretq_u8_u64(gt)));
            return vpadalq_u16(vld1q_u32(b.cp_count), vpaddlq_u8(bytes));
        }
        /* cp_occ4_neon's sentinel correction: the sentinel reads as A, so lane A is one less
         * for a row after it in its block ((uint64_t)(row - sent_after) < sent_len). */
        __attribute__((always_inline)) static inline uint32x4_t cp_occ4_neon_sent(uint32x4_t o, int64_t row,
                                                                                int64_t sent_after,
                                                                                uint64_t sent_len)
        {
            const uint32x4_t lane_a = vsetq_lane_u32(1, vdupq_n_u32(0), 0);
            if (__builtin_expect((uint64_t)(row - sent_after) < sent_len, 0)) o = vsubq_u32(o, lane_a);
            return o;
        }
#endif

        /* The occurrences of c in bwt[0, pp), `b` being pp's block. */
        static inline int64_t cp_occ_v(const CpView &v, const CP_OCC_PLANES &b, int64_t pp, int c,
                                       uint64_t hx, uint64_t lx)
        {
            int64_t occ = (int64_t)b.cp_count[c] + cp_rows(b, pp, hx, lx);
            if (__builtin_expect(v.wide, 0)) occ += v.base[pp >> v.span_shift][c];
            /* The range test comes before anything about c: c is the genome's own
             * base, so a branch on it would mispredict. */
            if (__builtin_expect((uint64_t)(pp - v.sent_after) < v.sent_len, 0)) occ -= (int64_t)(c == 0);
            return occ;
        }

        /* cp_occ_v for a base A and a count layout fixed at compile time, for
         * loops that extend many intervals by one base (the backward SMEM
         * phase): A's rows are two and/bic ops on the planes with no plane
         * masks to apply, there is no wide test, and only an A count tests the
         * sentinel range (the sentinel reads as A, so no other base's count
         * sees it). A = CP_CGT_RUNTIME takes base C, G or T from `a` at run time
         * (never A: see CP_CGT_RUNTIME), so it skips the sentinel test too. */
        template <int A, bool WIDE>
        static inline int64_t cp_occ_t(const CpView &v, const CP_OCC_PLANES &b, int64_t pp, int a)
        {
            if (A == CP_CGT_RUNTIME) {
#ifdef BWA_MEM3_DEBUG_CP_CGT_RUNTIME
                xassert(a >= 1 && a <= 3, "cp_occ_t<CP_CGT_RUNTIME> takes base C, G or T only");
#endif
                int64_t occ = (int64_t)b.cp_count[a] + cp_rows(b, pp, cp_hx(a), cp_lx(a));
                if (WIDE) occ += v.base[pp >> v.span_shift][a];
                return occ;
            }
            constexpr int C = A < 0 ? 0 : A;
#if defined(__ARM_NEON) || defined(__aarch64__) || defined(APPLE_SILICON)
            const uint64x1_t h = vld1_u64(&b.hi), l = vld1_u64(&b.lo);
            const uint64x1_t m = vld1_u64(&one_hot_mask_array[pp & CP_MASK]);
            const uint64x1_t mh = (C >> 1) ? vand_u64(m, h) : vbic_u64(m, h);
            const uint64x1_t x  = (C & 1) ? vand_u64(mh, l) : vbic_u64(mh, l);
            int64_t occ = (int64_t)b.cp_count[C] + (int64_t)vaddv_u8(vcnt_u8(vreinterpret_u8_u64(x)));
#else
            const uint64_t m = one_hot_mask_array[pp & CP_MASK];
            const uint64_t mh = (C >> 1) ? (m & b.hi) : (m & ~b.hi);
            const uint64_t x  = (C & 1) ? (mh & b.lo) : (mh & ~b.lo);
            int64_t occ = (int64_t)b.cp_count[C] + (int64_t)__builtin_popcountll(x);
#endif
            if (WIDE) occ += v.base[pp >> v.span_shift][C];
            if (A == 0 && __builtin_expect((uint64_t)(pp - v.sent_after) < v.sent_len, 0)) occ -= 1;
            return occ;
        }

        /* Return f.run<A, WIDE>(a) for base a (0..3) and count layout wide: the one
         * (wide, a) -> <A, WIDE> table the compact layout's per-base steps are picked
         * through (backwardExt_konly_planes, ls_advance_backward_step). F is a functor with
         * a member template `template <int A, bool WIDE> R run(int a)`. x86 takes a copy
         * per base, where each one's plane masks, count and checkpoint column are
         * constants. arm64 takes two, for A (the only base that tests the sentinel) and
         * for the other three (A = CP_CGT_RUNTIME, the base `a` at run time): the base is
         * the read's, so a four-way pick mispredicts on most rounds, and arm64 keeps the
         * run-time base's values in registers. */
        template <class R, class F>
        __attribute__((always_inline)) static inline R cp_step_dispatch(bool wide, int a, const F &f)
        {
#if defined(__aarch64__)
            /* a == 0 must take the A copy: CP_CGT_RUNTIME skips A's sentinel correction. */
            if (wide) return a == 0 ? f.template run<0, true>(0) : f.template run<CP_CGT_RUNTIME, true>(a);
            return a == 0 ? f.template run<0, false>(0) : f.template run<CP_CGT_RUNTIME, false>(a);
#else
            switch ((wide ? 4 : 0) | a) {
            case 0: return f.template run<0, false>(0);
            case 1: return f.template run<1, false>(1);
            case 2: return f.template run<2, false>(2);
            case 3: return f.template run<3, false>(3);
            case 4: return f.template run<0, true>(0);
            case 5: return f.template run<1, true>(1);
            case 6: return f.template run<2, true>(2);
            default: return f.template run<3, true>(3);
            }
#endif
        }

        /* backwardExt_konly on the compact layout for a base and count layout fixed at
         * compile time (or base C, G or T in `a` when A is CP_CGT_RUNTIME): the step the
         * seeding drivers run (ls_backward_one, ls_backward_general) and
         * backwardExt_konly_planes. k is left stale on the s'->0 exit, so the caller must
         * discard any s' < min_intv. The one-row case tests the row's base on its planes;
         * the sentinel row reads as A there but has no base, so it extends nothing. */
        template <int A, bool WIDE>
        __attribute__((always_inline)) static inline
        SMEM backwardExt_konly_ab(const CpView &v, int64_t count_a, SMEM smem, int a)
        {
            const int64_t sp = (int64_t)smem.k;
            if (smem.s == 1) {
                const CP_OCC_PLANES &blk = v.cp[sp >> CP_SHIFT];
                /* the sentinel row reads as A: only an A extension tests it */
                if (cp_base_at(blk, sp) == (A == CP_CGT_RUNTIME ? a : A) &&
                    (A != 0 || sp != v.sent_after - 1)) {
                    smem.k = count_a + cp_occ_t<A, WIDE>(v, blk, sp, a);
                    smem.s = 1;
                } else {
                    smem.s = 0;
                }
                return smem;
            }
            const int64_t ep = sp + (int64_t)smem.s;
            const int64_t occ_s = cp_occ_t<A, WIDE>(v, v.cp[sp >> CP_SHIFT], sp, a);
            const int64_t occ_e = cp_occ_t<A, WIDE>(v, v.cp[ep >> CP_SHIFT], ep, a);
            smem.k = count_a + occ_s;
            smem.s = occ_e - occ_s;
            return smem;
        }

        /* backwardExt_konly_ab with the base known only at run time: the plane masks
         * (cp_hx/cp_lx) select a's rows and the count layout is tested per count (it never
         * changes, so the branch predicts), instead of a dispatch on the query base, which
         * is unpredictable. For ls_backward_one, which extends a single interval, so a
         * dispatch there is not amortised over a loop. Same results as
         * backwardExt_konly_ab<a, wide>, stale k included. */
        __attribute__((always_inline)) static inline
        SMEM backwardExt_konly_rt(const CpView &v, int64_t count_a, SMEM smem, int a)
        {
            const int64_t sp = (int64_t)smem.k;
            const uint64_t hx = cp_hx(a), lx = cp_lx(a);
            if (smem.s == 1) {
                const CP_OCC_PLANES &blk = v.cp[sp >> CP_SHIFT];
                /* the sentinel row reads as A: only an A extension tests it */
                if ((cp_base_at(blk, sp) == a) & ((a != 0) | (sp != v.sent_after - 1))) {
                    smem.k = count_a + cp_occ_v(v, blk, sp, a, hx, lx);
                    smem.s = 1;
                } else {
                    smem.s = 0;
                }
                return smem;
            }
            const int64_t ep = sp + (int64_t)smem.s;
            const int64_t occ_s = cp_occ_v(v, v.cp[sp >> CP_SHIFT], sp, a, hx, lx);
            const int64_t occ_e = cp_occ_v(v, v.cp[ep >> CP_SHIFT], ep, a, hx, lx);
            smem.k = count_a + occ_s;
            smem.s = occ_e - occ_s;
            return smem;
        }

        inline int64_t occ_at(const CP_OCC_PLANES &b, int64_t pp, int c) const
        {
            const CpView v = cp_view();
            return cp_occ_v(v, b, pp, c, cp_hx(c), cp_lx(c));
        }

        /* The four occurrence counts of bwt[0, row), `b` being row's block, from
         * three popcounts: T is where both planes are set, G and C the rest of
         * each plane, A what remains of the rows. */
        inline void occ4_at(const CP_OCC_PLANES &b, int64_t row, int64_t occ[4]) const
        {
            const int64_t y = row & CP_MASK;
            const uint64_t m = one_hot_mask_array[y];
            const int64_t t = popcount64(b.hi & b.lo & m);
            const int64_t g = popcount64(b.hi & m) - t;
            const int64_t c = popcount64(b.lo & m) - t;
            int64_t a = y - t - g - c;
            if (__builtin_expect((uint64_t)(row - cp_sent_after) < cp_sent_len, 0)) a -= 1;
            occ[0] = (int64_t)b.cp_count[0] + a;
            occ[1] = (int64_t)b.cp_count[1] + c;
            occ[2] = (int64_t)b.cp_count[2] + g;
            occ[3] = (int64_t)b.cp_count[3] + t;
            if (__builtin_expect(cp_wide, 0)) {
                const int64_t *base = cp_base[row >> cp_span_shift];
                for (int k = 0; k < 4; k++) occ[k] += base[k];
            }
        }

        static inline int64_t popcount64(uint64_t x)
        {
#if defined(__ARM_NEON) || defined(__aarch64__) || defined(APPLE_SILICON)
            return (int64_t) vaddv_u8(vcnt_u8(vcreate_u8(x)));
#else
            return (int64_t) __builtin_popcountll(x);
#endif
        }

        /* The compact layout's full backward step (backwardExt_file's contract). */
        __attribute__((always_inline)) inline
        SMEM backwardExt_planes(SMEM smem, uint8_t a) const
        {
            /* sp/ep and their checkpoint blocks do not depend on the base b, so
             * they are loaded once. Only k[a] and s[a] are consumed for the
             * result, but all four s[b] feed the l cumulation below. */
            const int64_t sp = (int64_t)smem.k;
            const int64_t ep = (int64_t)smem.k + (int64_t)smem.s;
            const CP_OCC_PLANES &blk_sp = cp_planes[sp >> CP_SHIFT];
            const CP_OCC_PLANES &blk_ep = cp_planes[ep >> CP_SHIFT];

            int64_t occ_sp[4], s[4], l[4];

#if defined(__ARM_NEON) || defined(__aarch64__) || defined(APPLE_SILICON)
            if (__builtin_expect(!cp_wide, 1)) {
                /* All four bases of each block in one vector (cp_occ4_neon). */
                const uint32x4_t o_sp = cp_occ4_neon_sent(cp_occ4_neon(blk_sp, sp), sp, cp_sent_after, cp_sent_len);
                const uint32x4_t o_ep = cp_occ4_neon_sent(cp_occ4_neon(blk_ep, ep), ep, cp_sent_after, cp_sent_len);
                /* 32-bit lanes hold every count here (!cp_wide); widen on store. */
                const uint32x4_t s32 = vsubq_u32(o_ep, o_sp);
                vst1q_s64(&occ_sp[0], vreinterpretq_s64_u64(vmovl_u32(vget_low_u32(o_sp))));
                vst1q_s64(&occ_sp[2], vreinterpretq_s64_u64(vmovl_high_u32(o_sp)));
                vst1q_s64(&s[0], vreinterpretq_s64_u64(vmovl_u32(vget_low_u32(s32))));
                vst1q_s64(&s[2], vreinterpretq_s64_u64(vmovl_high_u32(s32)));
            } else
#endif
            {
                int64_t occ_ep[4];
                occ4_at(blk_sp, sp, occ_sp);
                occ4_at(blk_ep, ep, occ_ep);
                for (int b = 0; b < 4; b++) s[b] = occ_ep[b] - occ_sp[b];
            }

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

        /* Forward extension of `smem` by query base `a`: the interval of the
         * pattern with `a` appended on the right (n unchanged). The same result
         * as swapping k and l, calling backwardExt with the complement 3 - a and
         * swapping back, but only what the swapped-back result reads is
         * computed: base c = 3 - a's count at sp = smem.l and its size, and the
         * summed sizes of the bases above c (the l recurrence, which becomes k).
         * backwardExt builds all four k, s and l and picks one lane through the
         * stack. A one-row interval needs one block and one count: the row's own
         * base decides s, and whether it is above c (or the sentinel) moves k. */
        __attribute__((always_inline)) inline
        SMEM forwardExt_planes(SMEM smem, uint8_t a) const
        {
            const int c = 3 - (int)a;
            if (__builtin_expect(cp_wide, 0)) {
                SMEM r = smem;
                r.k = smem.l; r.l = smem.k;
                SMEM e = backwardExt_planes(r, (uint8_t)c);
                SMEM out = e;
                out.k = e.l; out.l = e.k;
                return out;
            }
            const int64_t sp = smem.l;
            const CpView v = cp_view();
            if (smem.s == 1) {
                const CP_OCC_PLANES &b = cp_planes[sp >> CP_SHIFT];
                const int sym = cp_base_at(b, sp);
                const bool sent = sp == sentinel_index;
                const int64_t kc = cp_occ_v(v, b, sp, c, cp_hx(c), cp_lx(c));
                smem.k = smem.k + (int64_t)(sent || sym > c);
                smem.l = count[c] + kc;
                smem.s = (int64_t)(!sent && sym == c);
                return smem;
            }
            const int64_t ep = sp + smem.s;
            const CP_OCC_PLANES &bsp = cp_planes[sp >> CP_SHIFT];
            const CP_OCC_PLANES &bep = cp_planes[ep >> CP_SHIFT];
            const int64_t so = (int64_t)((sp <= sentinel_index) & (ep > sentinel_index));
            int64_t kc, sc, sgt;
#if defined(__ARM_NEON) || defined(__aarch64__) || defined(APPLE_SILICON)
            {
                /* the sentinel reads as A: only lane A is corrected, and only
                 * kc/sc read it (A is never above c) */
                const uint32x4_t o_sp = cp_occ4_neon_sent(cp_occ4_neon(bsp, sp), sp, v.sent_after, v.sent_len);
                const uint32x4_t o_ep = cp_occ4_neon_sent(cp_occ4_neon(bep, ep), ep, v.sent_after, v.sent_len);
                static const uint32_t lanes[4] = {0, 1, 2, 3};
                const uint32x4_t lane = vld1q_u32(lanes);
                const uint32x4_t cv = vdupq_n_u32((uint32_t)c);
                const uint32x4_t eq = vceqq_u32(lane, cv), gt = vcgtq_u32(lane, cv);
                const uint32x4_t d = vsubq_u32(o_ep, o_sp);
                kc  = (int64_t)vaddvq_u32(vandq_u32(o_sp, eq));
                sc  = (int64_t)vaddvq_u32(vandq_u32(d, eq));
                sgt = (int64_t)vaddlvq_u32(vandq_u32(d, gt));
            }
#else
            {
                /* Base c's rows are (hi ^ hx) & (lo ^ lx); the rows above c are
                 * (hi & mh) | (lo & ml) | (hi & lo & mhl). Masks from a table: c
                 * is a random base, so a branch on it would mispredict. */
                static const uint64_t SEL[4][5] = {
                    {~0ULL, ~0ULL, ~0ULL, ~0ULL, 0},    /* A: above = C, G, T = hi | lo */
                    {~0ULL, 0, ~0ULL, 0, 0},            /* C: above = G, T = hi */
                    {0, ~0ULL, 0, 0, ~0ULL},            /* G: above = T = hi & lo */
                    {0, 0, 0, 0, 0},                    /* T: none */
                };
                const uint64_t *sel = SEL[c & 3];
                const uint64_t msp = one_hot_mask_array[sp & CP_MASK], mep = one_hot_mask_array[ep & CP_MASK];
                kc = cp_occ_v(v, bsp, sp, c, sel[0], sel[1]);
                sc = cp_occ_v(v, bep, ep, c, sel[0], sel[1]) - kc;
                /* the checkpoint counts of the bases above c, as a suffix sum picked
                 * by c rather than a loop whose trip count depends on c */
                const uint32_t d1 = bep.cp_count[1] - bsp.cp_count[1];
                const uint32_t d2 = bep.cp_count[2] - bsp.cp_count[2];
                const uint32_t d3 = bep.cp_count[3] - bsp.cp_count[3];
                const uint32_t above[4] = { d1 + d2 + d3, d2 + d3, d3, 0 };
                const uint32_t above_cnt = above[c & 3];
                const uint64_t ae = ((bep.hi & sel[2]) | (bep.lo & sel[3]) | (bep.hi & bep.lo & sel[4])) & mep;
                const uint64_t as = ((bsp.hi & sel[2]) | (bsp.lo & sel[3]) | (bsp.hi & bsp.lo & sel[4])) & msp;
                /* wrapping u32 is exact: the rows above c in [sp, ep) are at most s */
                sgt = (int64_t)(uint32_t)(above_cnt + (uint32_t)__builtin_popcountll(ae) - (uint32_t)__builtin_popcountll(as));
            }
#endif
            smem.k = smem.k + so + sgt;
            smem.l = count[c] + kc;
            smem.s = sc;
            return smem;
        }
        /* backwardExt_konly_file's contract, on the compact layout: the drivers' own step
         * (backwardExt_konly_ab), so the tests that call backwardExt_konly check the code the
         * aligner runs. */
        struct KonlyStep {
            const CpView &v;
            int64_t count_a;
            SMEM smem;
            template <int A, bool WIDE> SMEM run(int a) const { return backwardExt_konly_ab<A, WIDE>(v, count_a, smem, a); }
        };
        SMEM backwardExt_konly_planes(SMEM smem, uint8_t a) const
        {
            const CpView v = cp_view();
            return cp_step_dispatch<SMEM>(cp_wide, a, KonlyStep{v, count[a], smem});
        }


        /* ---- One layout fixed at compile time (PL: the compact layout) ----
         * The hot drivers are templates on PL, instantiated once per layout and
         * picked once per call of a public entry point, so no lookup tests the
         * layout. For PL = false each of these is exactly the default layout's
         * own code. */
        template <bool PL> inline const char *cp_line(int64_t row) const
        {
            return PL ? (const char *)&cp_planes[row >> CP_SHIFT] : (const char *)&cp_occ[row >> CP_SHIFT];
        }
        template <bool PL> __attribute__((always_inline)) inline
        SMEM backwardExt_t(SMEM smem, uint8_t a) const
        {
            return PL ? backwardExt_planes(smem, a) : backwardExt_file(smem, a);
        }
        /* Forward extension of `smem` by query base `a`: backward extension with
         * the BWT of the reverse complement (swap k <-> l, backwardExt by 3 - a,
         * swap back). `l` is LIVE here: after the swap it is the k that drives
         * the next forward step, so the full backwardExt's l-cumulation is
         * required and backwardExt_konly does not apply. The compact layout
         * computes only what the swapped-back result reads (forwardExt_planes). */
        template <bool PL> __attribute__((always_inline)) inline
        SMEM forwardExt_t(SMEM smem, uint8_t a) const
        {
            if (PL) return forwardExt_planes(smem, a);
            SMEM smem_ = smem;
            smem_.k = smem.l;
            smem_.l = smem.k;
            SMEM newSmem_ = backwardExt_file(smem_, 3 - a);
            SMEM newSmem = newSmem_;
            newSmem.k = newSmem_.l;
            newSmem.l = newSmem_.k;
            return newSmem;
        }

        /* The same three steps on this object's layout, tested per call: for the
         * tests (the seed jump table build dispatches once, seed_jump_children_t),
         * never a hot loop. */
        SMEM backwardExt(SMEM smem, uint8_t a) const
        {
            return cp_compact ? backwardExt_t<true>(smem, a) : backwardExt_t<false>(smem, a);
        }
        SMEM forwardExt(SMEM smem, uint8_t a) const
        {
            return cp_compact ? forwardExt_t<true>(smem, a) : forwardExt_t<false>(smem, a);
        }
        SMEM backwardExt_konly(SMEM smem, uint8_t a) const
        {
            return cp_compact ? backwardExt_konly_planes(smem, a) : backwardExt_konly_file(smem, a);
        }

        /* The SA resolver's LF walks (FMI_search.cpp): one template body per walk, its LF
         * step the only layout-specific part, and one out-of-line entry per layout. */
        inline int64_t sa_sample_value(int64_t sp) const;
        inline bool sa_lf_step_planes(int64_t &sp) const;
        template <bool PL> inline bool sa_lf_step_walk(int64_t &sp) const;
        template <bool PL> inline bool sa_lf_step_one(int64_t &sp) const;
        template <bool PL> inline int64_t get_sa_entry_walk(int64_t pos);
        template <bool PL> inline int64_t call_one_step_walk(int64_t pos, int64_t &sa_entry, int64_t &offset,
                                                             int drop_sentinel_offset);
        int64_t get_sa_entry_compressed_file(int64_t pos, int tid);
        int64_t get_sa_entry_compressed_planes(int64_t pos, int tid);
        template <bool PL> inline int64_t get_sa_entry_compressed_t(int64_t pos, int tid)
        {
            return PL ? get_sa_entry_compressed_planes(pos, tid) : get_sa_entry_compressed_file(pos, tid);
        }
        int64_t call_one_step_file(int64_t pos, int64_t &sa_entry, int64_t &offset, int drop_sentinel_offset);
        int64_t call_one_step_planes(int64_t pos, int64_t &sa_entry, int64_t &offset, int drop_sentinel_offset);
        template <bool PL> inline int64_t call_one_step_t(int64_t pos, int64_t &sa_entry, int64_t &offset,
                                                          int drop_sentinel_offset)
        {
            return PL ? call_one_step_planes(pos, sa_entry, offset, drop_sentinel_offset)
                      : call_one_step_file(pos, sa_entry, offset, drop_sentinel_offset);
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
    template <bool PL> void ls_prefetch_cp_occ(const BatchSlot *s);
    inline void ls_bsd_note(BatchSlot *s, int64_t k, int64_t l, int64_t sz);
    void ls_bsd_resolve(BatchSlot *s, const uint8_t *enc_qdb, int32_t p);
    template <bool PL> void ls_advance_forward_step(BatchSlot *s, const uint8_t *enc_qdb,
                                 bool s1_text);
    template <bool PL> inline void ls_prefetch_sa_row(int64_t p) const;
    template <bool PL> void ls_s1_lf_step(BatchSlot *s, const uint8_t *pac, int64_t l_pac);
    void ls_s1_compare(BatchSlot *s, const uint8_t *enc_qdb,
                       const uint8_t *pac, int64_t l_pac);
    template <bool PL> void ls_advance_forward_table(BatchSlot *s, const uint8_t *enc_qdb,
                                  bool s1_text);
    template <bool PL> void ls_prepare_backward(BatchSlot *s, int32_t minSeedLen);
    template <int A, bool WIDE>
    int32_t ls_backward_general(BatchSlot *s, SMEM *prev, int32_t numPrev, int32_t min_intv,
                                int32_t j, int32_t minSeedLen, int64_t count_a, int a);
    inline int32_t ls_backward_one(BatchSlot *s, int a, int32_t minSeedLen);
    /* ls_backward_general as a cp_step_dispatch functor. */
    struct LsBackwardGeneral {
        FMI_search *f; BatchSlot *s; SMEM *prev; int32_t numPrev, min_intv, j, minSeedLen; int64_t count_a;
        template <int A, bool WIDE> int32_t run(int a) const
        {
            return f->ls_backward_general<A, WIDE>(s, prev, numPrev, min_intv, j, minSeedLen, count_a, a);
        }
    };
    template <bool PL> void ls_advance_backward_step(BatchSlot *s,
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
    template <bool PL> void bsd_prefetch_cp_occ(const BwtSeedSlot *s);
    void bsd_plan_jump(BwtSeedSlot *s, const uint8_t *enc_qdb,
                       const SeedJumpEntry *jlev, int32_t jd);
    template <bool PL> bool bsd_start_walk(BwtSeedSlot *s, const uint8_t *enc_qdb,
                                           const SeedJumpEntry *jlev, int32_t jd);
    template <bool PL> void bsd_advance_step(BwtSeedSlot *s,
                          const uint8_t *enc_qdb,
                          int32_t minSeedLen,
                          const SeedJumpEntry *jlev, int32_t jd);
    // ----- The public drivers, one instantiation per checkpoint layout -----
    // (PL: the compact layout). The public entry points of the same names pick
    // one by cp_compact, once per call.
    template <bool PL>
    void getSMEMsOnePosOneThread_t(uint8_t *enc_qdb, int32_t *query_pos_array, int32_t *min_intv_array,
                                   int32_t *rid_array, int32_t numReads, int32_t batch_size,
                                   const bseq1_t *seq_, int32_t *query_cum_len_ar, int32_t max_readlength,
                                   int32_t minSeedLen, SMEM *matchArray, int64_t *__numTotalSmem);
    template <bool PL>
    void getSMEMsOnePosOneThread_lockstep_t(uint8_t *enc_qdb, int32_t *query_pos_array,
                                            int32_t *min_intv_array, int32_t *rid_array, int32_t numReads,
                                            int32_t batch_size, const bseq1_t *seq_,
                                            int32_t *query_cum_len_ar, int32_t max_readlength,
                                            int32_t minSeedLen, SMEM *matchArray, int64_t *__numTotalSmem,
                                            BsdFirstWalkHint *bsd_hint);
    template <bool PL>
    int64_t bwtSeedStrategyAllPosOneThread_t(uint8_t *enc_qdb, int32_t *max_intv_array, int32_t numReads,
                                             const bseq1_t *seq_, int32_t *query_cum_len_ar,
                                             int32_t minSeedLen, SMEM *matchArray);
    template <bool PL>
    int64_t bwtSeedStrategyAllPosOneThread_lockstep_t(uint8_t *enc_qdb, int32_t *max_intv_array,
                                                      int32_t numReads, const bseq1_t *seq_,
                                                      int32_t *query_cum_len_ar, int32_t minSeedLen,
                                                      SMEM *matchArray, int32_t max_readlength,
                                                      const BsdFirstWalkHint *bsd_hint);
    template <bool PL>
    void get_sa_entries_t(SMEM *smemArray, int64_t *coordArray, int32_t *coordCountArray,
                          uint32_t count, int32_t max_occ, int tid);
    template <bool PL>
    void get_sa_entries_prefetch_t(SMEM *smemArray, int64_t *coordArray, int64_t *coordCountArray,
                                   int64_t count, const int32_t max_occ, int tid, int64_t &id_,
                                   int drop_sentinel_offset, bool trust_s1_tag);
};

#endif
