#ifndef BWA_LIBSAIS_BUILD_H
#define BWA_LIBSAIS_BUILD_H

#include <cstdint>
#include <string>

struct LibsaisBuildOpts {
    int64_t     max_memory_bytes = 0;
    int         num_threads      = 0;
    std::string tmpdir;
};

// Build the bwa-mem3 FM index via libsais's generalized-suffix-array
// construction. Precondition: `<prefix>.pac` and `<prefix>.ann` already
// exist, with the .pac encoding the forward-only bases emitted by
// bns_fasta2bntseq (l_pac bases, 2-bit, alphabet A=0 C=1 G=2 T=3; N was
// replaced by a pseudo-random base).
//
// Emits `<prefix>.0123` and `<prefix>.bwt.2bit.64`, byte-identical to the
// historical sais-lite-based build.
int libsais_build_fm_index(const char* prefix, int64_t pac_len,
                           const LibsaisBuildOpts& opts);

#endif
