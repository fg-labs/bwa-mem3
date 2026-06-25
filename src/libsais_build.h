#ifndef BWA_LIBSAIS_BUILD_H
#define BWA_LIBSAIS_BUILD_H

#include <cstdint>
#include <string>

struct LibsaisBuildOpts {
    int64_t     max_memory_bytes = 0;
    int         num_threads      = 0;
    std::string tmpdir;
    /* Emit the unpacked `<prefix>.0123` reference. The D3 --meth SEED index sets
     * this false: `mem --meth` extends against the ORIGINAL `.0123`, never the
     * seed's, so emitting the seed `.0123` (~13 GB on hg38) would only write a
     * file nothing reads. The FM build itself never consumes the `.0123`. */
    bool        emit_unpacked_ref = true;
};

// Build the bwa-mem3 FM index via libsais's generalized-suffix-array
// construction. Precondition: `<prefix>.pac` and `<prefix>.ann` already
// exist, with the .pac encoding the forward-only bases emitted by
// bns_fasta2bntseq (l_pac bases, 2-bit, alphabet A=0 C=1 G=2 T=3; N was
// replaced by a pseudo-random base).
//
// Emits `<prefix>.bwt.2bit.64` (and `<prefix>.0123` unless opts.emit_unpacked_ref
// is false), byte-identical to the historical sais-lite-based build.
int libsais_build_fm_index(const char* prefix, int64_t pac_len,
                           const LibsaisBuildOpts& opts);

#endif
