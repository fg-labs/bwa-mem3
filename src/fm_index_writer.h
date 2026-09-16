#ifndef BWA_FM_INDEX_WRITER_H
#define BWA_FM_INDEX_WRITER_H

#include <cstdint>

#include "fmi_seed_api.h"   // CP_OCC, CP_SHIFT

// Fixed .bwt.2bit.64 header: ref_seq_len (int64) + count[5] (5×int64).
static const int64_t FMI_BWT2BIT_HDR_BYTES = 6 * (int64_t)sizeof(int64_t);

// The byte layout of a .bwt.2bit.64 file for a given reference length and SA
// sample-rate shift. Single source of truth for the offset arithmetic that was
// previously hand-duplicated across the writer, the resampler, the shm packer,
// and the disk loader — every one of those must agree bit-for-bit or reads
// mis-align, so the formula lives in one place. Section order on disk:
//   [header] [cp_occ] [sa_ms_byte] [sa_ls_word] [sentinel_index] [sa_compx tag]
struct Bwt2bitLayout {
    int64_t cp_occ_count;      // (ref_seq_len >> CP_SHIFT) + 1
    int64_t sa_sample_count;   // (ref_seq_len >> sa_compx) + 1
    int64_t off_cp_occ;        // == FMI_BWT2BIT_HDR_BYTES
    int64_t off_ms_byte;
    int64_t off_ls_word;
    int64_t off_sentinel;
    int64_t off_sa_compx_tag;  // sentinel + int64
    int64_t total_bytes;       // through the trailing sa_compx tag
};

// Compute the section offsets for `ref_seq_len` at SA sample-rate shift
// `sa_compx`. Pure arithmetic, no I/O. `sa_compx` must be in [0, CP_SHIFT].
Bwt2bitLayout fmi_bwt2bit_layout(int64_t ref_seq_len, int sa_compx);

// Write a bwa-mem2 .bwt.2bit.64 FM index directly from the libsais GSA
// output, without materialising the full BWT byte array or the dense
// SA-sample array in RAM. The writer walks the SA twice:
//
//   Pass 1 (parallel): per-stripe BWT-base histograms + sentinel location.
//   Pass 2 (parallel): per-stripe emit of cp_occ blocks + every-8th SA
//                      sample through small per-thread ring buffers,
//                      pwritten to disjoint file ranges.
//
// Peak residency inside the writer is O(KB) per thread regardless of
// input size: no bwt[], no samples[].
//
// `buf`            — libsais-alphabet input buffer of length pac_len+1
//                    ({0=$, 1=A, 2=C, 3=G, 4=T}). The writer never
//                    writes to it; it's the source of every BWT byte
//                    via BWT[i] = alphabet_shift(buf[(SA[i]-1) mod N+1]).
// `sa`             — pointer to the GSA output. Either int32_t* (when
//                    sa_is_64bit is false) or int64_t*, length pac_len+1.
// `sa_is_64bit`    — true when `sa` is int64_t*, false for int32_t*.
//                    libsais picks based on whether N+1 fits in INT32_MAX.
// `pac_len`        — doubled-text length N; the file header stores N+1.
// `count`          — 5-entry prefix-sum histogram over the doubled text
//                    (count[0]=0, count[c] = sum(freq[<c]) for c in 1..4).
//                    Typically from compute_counts() on the same text.
// `out_sentinel_index` — set by pass 1 to the single BWT row where
//                    SA[i]=0 (i.e. the full text comes right after $).
// `sa_compx`       — SA sample-rate shift: one SA row in (1<<sa_compx) is
//                    stored. Persisted as a trailing int64_t appended after
//                    sentinel_index so the loader can recover it; existing
//                    header/cp_occ/SA-sample/sentinel offsets are unchanged
//                    by its presence.
// `num_threads`    — 1 keeps the serial path; > 1 partitions work into
//                    CP_BLOCK_SIZE-aligned stripes. OpenMP-based; reuses
//                    the existing libomp thread pool (shared with libsais)
//                    so no new arenas / stacks are allocated.
//
// Fatal-errors via err_fatal on any failure; returns normally on success.
//
// Threading contract: must be invoked from a serial context (no enclosing
// OpenMP parallel region). The implementation issues
// `#pragma omp parallel num_threads(T)` which would otherwise nest a fresh
// team of size T inside the caller's region, regardless of
// OMP_NESTED / OMP_MAX_ACTIVE_LEVELS. The current call site
// (libsais_build_fm_index, after `libsais_gsa_omp` has joined) satisfies
// this; future callers must too.
void write_fm_index_streaming(const char* out_path,
                              const uint8_t* buf,
                              const void* sa,
                              bool sa_is_64bit,
                              int64_t pac_len,
                              const int64_t count[5],
                              int64_t* out_sentinel_index,
                              int sa_compx,
                              int num_threads = 1);

// Rewrite an existing .bwt.2bit.64 with a resampled SA sample table, reusing
// the already-materialised cp_occ block and header of the source index. Unlike
// write_fm_index_streaming (which recomputes everything from the raw suffix
// array during `index`), this only re-lays-out the SA-sample layer at a new
// sample rate -- the `re-sa` command persists on disk what `shm -u` synthesises
// per-stage. The caller supplies:
//
//   `ref_seq_len`      — text length N+1, written verbatim to the header.
//   `count`            — the RAW on-disk cumulative C[] (as read back, i.e.
//                        NOT the +1-adjusted in-memory form); written verbatim.
//   `cp_occ`           — the source checkpoint-occ block, `cp_occ_size` entries,
//                        copied through unchanged (it is independent of the SA
//                        sample rate).
//   `sa_ms_byte`/`sa_ls_word` — the resampled SA samples, `sa_sample_cnt` =
//                        (ref_seq_len >> sa_compx) + 1 entries, already computed
//                        for the target `sa_compx` (denser: FMI_search::
//                        densify_sa_into; coarser: decimation of the disk table).
//   `sentinel_index`   — unchanged from the source; written verbatim.
//   `sa_compx`         — the new SA sample-rate shift, persisted as the trailing
//                        int64_t tag so the loader auto-detects the new rate.
//
// Writes to `<out_path>.tmp` and rename(2)s over `out_path` on success, so a
// concurrent reader sees either the old or new index, never a partial. Sibling
// index files (.pac/.ann/.amb/.0123) are untouched and must not be moved.
// Fatal-errors via err_fatal on any failure; returns normally on success.
void rewrite_fm_index_resampled_sa(const char* out_path,
                                   int64_t ref_seq_len,
                                   const int64_t count[5],
                                   const CP_OCC* cp_occ,
                                   int64_t cp_occ_size,
                                   const int8_t* sa_ms_byte,
                                   const uint32_t* sa_ls_word,
                                   int64_t sa_sample_cnt,
                                   int64_t sentinel_index,
                                   int sa_compx);

#endif
