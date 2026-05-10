// Shared test fixture for tests that need a loaded meth_orig_ref_t. Builds
// an in-memory doubled-c2t bns + pac (matching what bwa-mem3 index --meth
// produces) from a forward-strand sequence, then constructs a cmap and
// calls meth_orig_ref_load. Same code path the production runtime uses.

#ifndef BWAMEM3_TEST_METH_ORIG_REF_FIXTURE_H
#define BWAMEM3_TEST_METH_ORIG_REF_FIXTURE_H

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "bntseq.h"
#include "meth_bam.h"
#include "meth_orig_ref.h"

namespace meth_test {

inline uint8_t base_to_2bit(char c) {
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default:  return 0;  // N or unknown -> A (ambs side-table marks N)
    }
}

inline void pac_pack_set(uint8_t *pac, int64_t i, uint8_t v) {
    pac[i >> 2] |= (uint8_t)(v << ((~i & 3) << 1));
}

/* RAII fixture that owns:
 *   - a doubled-c2t bns (one chrom, contigs r-X then f-X)
 *   - a packed c2t pac (G->A on r-X, C->T on f-X)
 *   - the cmap built from the doubled bns
 *   - a meth_orig_ref_t loaded against bns + pac + cmap
 *
 * Construct from a forward-strand sequence (one chrom). Test code uses
 * .orig (the loaded handle) and .real_tid (always 0) for slice calls. */
struct OrigRefFixture {
    bntseq_t                  bns{};
    std::vector<bntann1_t>    anns;
    std::vector<bntamb1_t>    ambs;
    std::vector<std::string>  name_storage;
    std::vector<uint8_t>      pac;
    meth_chrom_map_t         *cmap = nullptr;
    meth_orig_ref_t          *orig = nullptr;
    int                       real_tid = 0;

    explicit OrigRefFixture(const std::string &forward_seq) {
        const int len = (int)forward_seq.size();
        bns.l_pac = (int64_t)len * 2;
        bns.n_seqs = 2;
        anns.resize(2);
        name_storage = {"rX", "fX"};
        for (int i = 0; i < 2; ++i) {
            anns[i] = bntann1_t{};
            anns[i].offset = (int64_t)i * len;
            anns[i].len    = len;
            anns[i].name   = const_cast<char *>(name_storage[i].c_str());
            anns[i].anno   = const_cast<char *>("");
        }
        bns.anns = anns.data();

        /* ambs: any 'N' in the input marks a 1-bp interval on BOTH halves
         * of the doubled pac (matching what bwa-mem3 index --meth produces).
         * bns_iter_ambi binary-searches bns->ambs assuming non-decreasing
         * offset+len, so emit all r-X intervals (offsets 0..len-1) before
         * any f-X interval (offsets len..2*len-1) to preserve that order. */
        for (int half = 0; half < 2; ++half) {
            for (int i = 0; i < len; ++i) {
                if (forward_seq[i] != 'N') continue;
                bntamb1_t a;
                a.offset = anns[half].offset + i;
                a.len = 1;
                a.amb = 'N';
                ambs.push_back(a);
            }
        }
        bns.ambs = ambs.data();
        bns.n_holes = (int)ambs.size();

        const int pac_bytes = ((int)bns.l_pac + 3) / 4;
        pac.assign(pac_bytes, 0);
        /* r-X: G->A applied to the forward strand. */
        for (int i = 0; i < len; ++i) {
            char b = forward_seq[i];
            char projected = (b == 'G') ? 'A' : b;
            pac_pack_set(pac.data(), anns[0].offset + i, base_to_2bit(projected));
        }
        /* f-X: C->T applied to the forward strand. */
        for (int i = 0; i < len; ++i) {
            char b = forward_seq[i];
            char projected = (b == 'C') ? 'T' : b;
            pac_pack_set(pac.data(), anns[1].offset + i, base_to_2bit(projected));
        }

        cmap = meth_chrom_map_build_from_bns(&bns);
        if (cmap == nullptr) std::abort();
        orig = meth_orig_ref_load(&bns, pac.data(), cmap);
        if (orig == nullptr) std::abort();
    }

    ~OrigRefFixture() {
        if (orig) meth_orig_ref_free(orig);
        if (cmap) meth_chrom_map_free(cmap);
    }

    OrigRefFixture(const OrigRefFixture &) = delete;
    OrigRefFixture &operator=(const OrigRefFixture &) = delete;
};

}  // namespace meth_test

#endif
