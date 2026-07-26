#ifndef BWA_LIBSAIS_BUILD_H
#define BWA_LIBSAIS_BUILD_H

#include <cstdint>
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
};

// True when a doubled text of length `N` forces libsais's 64-bit SA path,
// which doubles the per-suffix cost.
bool libsais_sa_is_64bit(int64_t N);

// Estimated peak bytes for a libsais build over a doubled text of length `N`
// (= 2 * l_pac). 6 B/base on the int32 SA path, 12 B/base once `N` forces the
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

#endif
