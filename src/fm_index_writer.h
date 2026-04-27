#ifndef BWA_FM_INDEX_WRITER_H
#define BWA_FM_INDEX_WRITER_H

#include <cstdint>

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
                              int num_threads = 1);

#endif
