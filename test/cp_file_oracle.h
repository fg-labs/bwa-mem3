/* Test oracle: the occurrence checkpoints exactly as <prefix>.bwt.2bit.64 stores
 * them (CP_OCC: 64-bit counts and one one-hot word per base), read straight from
 * the file, independent of the in-memory form the aligner converts them to.
 * occ(b, p) = cp_count[b] + popcount(one_hot_bwt_str[b] & mask[p & CP_MASK]);
 * a row with no plane set is the sentinel. */
#ifndef CP_FILE_ORACLE_H
#define CP_FILE_ORACLE_H

#include "fmi_seed_api.h"
#include "fm_index_writer.h"   /* FMI_BWT2BIT_HDR_BYTES */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

static inline std::vector<CP_OCC> cp_file_records(const char *prefix, int64_t ref_seq_len)
{
    const std::string path = std::string(prefix) + CP_FILENAME_SUFFIX;
    FILE *f = fopen(path.c_str(), "rb");
    if (f == NULL) { fprintf(stderr, "cannot open %s\n", path.c_str()); exit(2); }
    std::vector<CP_OCC> v((size_t) ((ref_seq_len >> CP_SHIFT) + 1));
    if (fseek(f, (long) FMI_BWT2BIT_HDR_BYTES, SEEK_SET) != 0 ||
        fread(v.data(), sizeof(CP_OCC), v.size(), f) != v.size()) {
        fprintf(stderr, "cannot read the checkpoints of %s\n", path.c_str());
        exit(2);
    }
    fclose(f);
    return v;
}

/* The SA samples exactly as <prefix>.bwt.2bit.64 stores them, after the checkpoints: one high
 * byte, then one low word, per sample at rate 1/2^compx. The aligner packs the high bytes to a
 * few bits in memory (sa_hi_at); an oracle built from these does not go through that packing. */
static inline void sa_file_samples(const char *prefix, int64_t ref_seq_len, int64_t compx,
                                   std::vector<int8_t> &ms, std::vector<uint32_t> &ls)
{
    const std::string path = std::string(prefix) + CP_FILENAME_SUFFIX;
    FILE *f = fopen(path.c_str(), "rb");
    if (f == NULL) { fprintf(stderr, "cannot open %s\n", path.c_str()); exit(2); }
    const size_t n = (size_t) ((ref_seq_len >> compx) + 1);
    const long at = (long) FMI_BWT2BIT_HDR_BYTES + (long) (((ref_seq_len >> CP_SHIFT) + 1) * (int64_t) sizeof(CP_OCC));
    ms.resize(n);
    ls.resize(n);
    if (fseek(f, at, SEEK_SET) != 0 || fread(ms.data(), sizeof(int8_t), n, f) != n ||
        fread(ls.data(), sizeof(uint32_t), n, f) != n) {
        fprintf(stderr, "cannot read the SA samples of %s\n", path.c_str());
        exit(2);
    }
    fclose(f);
}

static inline int64_t cp_file_occ(const std::vector<CP_OCC> &cp, int b, int64_t p)
{
    const CP_OCC &blk = cp[(size_t) (p >> CP_SHIFT)];
    return blk.cp_count[b] + (int64_t) __builtin_popcountll(blk.one_hot_bwt_str[b] & one_hot_mask_array[p & CP_MASK]);
}

/* The base at row p, or 4 for the sentinel row (no plane set). */
static inline int cp_file_base(const std::vector<CP_OCC> &cp, int64_t p)
{
    const CP_OCC &blk = cp[(size_t) (p >> CP_SHIFT)];
    const int y = (int) (CP_MASK - (p & CP_MASK));
    for (int c = 0; c < 4; c++)
        if ((blk.one_hot_bwt_str[c] >> y) & 1ULL) return c;
    return 4;
}

#endif
