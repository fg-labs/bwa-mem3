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
};

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
