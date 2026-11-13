#ifndef BWA_LIBSAIS_BUILD_H
#define BWA_LIBSAIS_BUILD_H

#include "fm_index_writer.h"
#include "packed_text.h"

#include <chrono>
#include <cstdint>
#include <functional>
#include <memory>
#include <string>

struct LibsaisBuildOpts {
    int64_t     max_memory_bytes = 0;
    int         num_threads      = 0;
    std::string tmpdir;
    /* Emit the unpacked `<prefix>.0123` reference. Defaults false: `mem`
     * pac-fetches the original reference from `<prefix>.pac` on demand and
     * never reads `.0123` (~8x the `.pac`; ~6.4 GB on hg38, ~13 GB for a
     * doubled --meth seed). The FM build itself never consumes the `.0123`.
     * Set true only to emit it for an external consumer (e.g. bwa-mem2). */
    bool        emit_unpacked_ref = false;
    /* True when --max-memory came from the command line rather than from
     * host detection. Only affects diagnostics: a "retry on a bigger host"
     * hint is noise when the binding constraint is the user's own flag. */
    bool        max_memory_user_specified = false;
    /* SA sample-rate shift passed to the streaming FM-index writer: one SA
     * row in (1<<sa_compx) is stored. Defaults to the historical rate (3,
     * i.e. 1/8) so a caller that doesn't set it builds byte-identical
     * indexes to before this field existed; a CLI flag to override it is
     * wired up by the caller in bwtindex.cpp, not here. */
    int         sa_compx           = 3;
};

// Forward text `fwd` (l_pac bases) followed by its reverse complement, 2-bit
// packed into (2*l_pac+3)/4 bytes: exactly the packed bytes write_doubled_pac
// writes. Parallelised over output bytes with `num_threads` OpenMP threads.
std::unique_ptr<uint8_t[]> pack_doubled(const PackedText& fwd, int64_t l_pac, int num_threads);

// Write `pack_doubled` of the forward .pac at `fwd_pac_path` to `out_path` in
// the on-disk .pac format (packed bytes, pad byte when 2*l_pac % 4 == 0, count byte).
void write_doubled_pac(const char* fwd_pac_path, int64_t l_pac,
                       const char* out_path, int num_threads);

// True when a doubled text of length `N` forces libsais's 64-bit SA path,
// which doubles the per-suffix cost.
bool libsais_sa_is_64bit(int64_t N);

// Estimated peak bytes for a libsais build over a doubled text of length `N`
// (= 2 * l_pac). 8 B/base on the int32 SA path, 12 B/base once `N` forces the
// int64 SA; both cover measured peaks with margin for libsais aux arrays and
// OMP/allocator overhead. Exposed so callers can preflight a build they have
// not started yet.
int64_t libsais_estimate_peak_bytes(int64_t N);

// Report a memory-budget shortfall on stderr in one standard form, naming what
// is needed, the budget, and the remedies. `what` labels the build being
// refused (e.g. "libsais"); `ref_bases` is the reference length in bp.
void libsais_report_budget_shortfall(const char* what, int64_t need_bytes,
                                     int64_t budget_bytes, bool budget_user_specified,
                                     int64_t ref_bases);

// Build the bwa-mem3 FM index via libsais's generalized-suffix-array
// construction. Precondition: `<prefix>.pac` and `<prefix>.ann` already
// exist, with the .pac encoding the forward-only bases emitted by
// bns_fasta2bntseq (l_pac bases, 2-bit, alphabet A=0 C=1 G=2 T=3; N was
// replaced by a pseudo-random base).
//
// Emits `<prefix>.bwt.2bit.64` (and, only when opts.emit_unpacked_ref is set,
// `<prefix>.0123`), byte-identical to the historical sais-lite-based build.
int libsais_build_fm_index(const char* prefix, int64_t pac_len,
                           const LibsaisBuildOpts& opts);

// Shared core of every build: phases 1-3 (unpack the doubled text `doubled`, libsais GSA,
// streaming FM-index write) over N = doubled.length() bases, T threads, `count` = the
// compute_counts() histogram of the text. The caller does the budget check and thread
// resolution beforehand. `open_sink` is called once, after the GSA is built, and returns the
// sink the index is written through; deferring it means a failure in phase 1 or 2 has
// created nothing (e.g. no partial .tmp file). `t_start` is when the caller's build began; the
// phase-3 "total" is measured from it. Returns 0; other failures are fatal.
int libsais_build_fm_core(const PackedText& doubled, const LibsaisBuildOpts& opts,
                          const int64_t count[5], int T,
                          std::chrono::steady_clock::time_point t_start,
                          const std::function<FmSink&()>& open_sink);

// In-memory build of a forward-only 2-bit text (l_pac bases): pack_doubled + core + MemSink.
// Uses no temp dir and no files. Returns 0 on success, 3 on a budget shortfall.
int libsais_build_fm_memory(const uint8_t* fwd_pac, int64_t l_pac, const LibsaisBuildOpts& opts,
                            FmArrays* out);

#endif
