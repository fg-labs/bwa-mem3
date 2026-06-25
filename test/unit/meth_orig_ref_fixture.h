// Shared test fixture for tests that need an ORIGINAL (un-converted) bns + pac
// — the same handles the D3 (PR-5) runtime loads via
// meth_orig_ref_load_handles (bns_restore + slurped 2-bit pac). Builds them
// in memory from a forward-strand sequence (one contig), storing the bases
// un-converted (no c2t fold) so meth_build_xm can decode forward-genome
// bases inline. Same code path the production runtime exercises.

#ifndef BWAMEM3_TEST_METH_ORIG_REF_FIXTURE_H
#define BWAMEM3_TEST_METH_ORIG_REF_FIXTURE_H

#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>

#include "bntseq.h"

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

/* RAII fixture that owns an original (un-converted) single-contig bns + pac:
 *   - bns: one contig "X" of length len at offset 0
 *   - pac: the forward sequence packed 2-bit, un-converted
 *   - ambs: one 1-bp 'N' interval per 'N' in the input
 *
 * Construct from a forward-strand sequence (one contig). Test code uses
 * .orig_bns / .orig_pac (the loaded handles) and .real_tid (always 0). */
struct OrigRefFixture {
    bntseq_t                  bns{};
    bntann1_t                 ann{};
    std::vector<bntamb1_t>    ambs;
    std::string               name_storage = "X";
    std::vector<uint8_t>      pac;
    bntseq_t                 *orig_bns = nullptr;
    uint8_t                  *orig_pac = nullptr;
    int                       real_tid = 0;

    explicit OrigRefFixture(const std::string &forward_seq) {
        const int len = (int)forward_seq.size();
        bns.l_pac  = (int64_t)len;
        bns.n_seqs = 1;

        ann.offset = 0;
        ann.len    = len;
        ann.n_ambs = 0;
        ann.name   = const_cast<char *>(name_storage.c_str());
        ann.anno   = const_cast<char *>("");
        bns.anns   = &ann;

        /* ambs: any 'N' in the input marks a 1-bp interval. bns_iter_ambi
         * binary-searches bns->ambs assuming non-decreasing offset+len, which
         * is satisfied by emitting them in input order. */
        for (int i = 0; i < len; ++i) {
            if (forward_seq[i] != 'N') continue;
            bntamb1_t a{};
            a.offset = i;
            a.len    = 1;
            a.amb    = 'N';
            ambs.push_back(a);
            ++ann.n_ambs;
        }
        bns.ambs    = ambs.data();
        bns.n_holes = (int)ambs.size();

        const int pac_bytes = (len + 3) / 4;
        pac.assign(pac_bytes, 0);
        /* Un-converted: store the forward sequence as-is (no c2t / g2a fold). */
        for (int i = 0; i < len; ++i) {
            pac_pack_set(pac.data(), i, base_to_2bit(forward_seq[i]));
        }

        orig_bns = &bns;
        orig_pac = pac.data();
    }

    ~OrigRefFixture() = default;

    OrigRefFixture(const OrigRefFixture &) = delete;
    OrigRefFixture &operator=(const OrigRefFixture &) = delete;
};

}  // namespace meth_test

#endif
