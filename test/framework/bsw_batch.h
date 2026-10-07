// test/framework/bsw_batch.h
//
// Score a batch of extension jobs (bwa_tests::ExtPair) through one banded-SW
// kernel copy's getScores8 / getScores16, results in input order. Shared by the
// bandedSWA kernel tests that compare kernel copies, groupings and the scalar
// oracle.

#ifndef BWA_TESTS_BSW_BATCH_H
#define BWA_TESTS_BSW_BATCH_H

#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

#include "bandedSWA.h"
#include "ext_ladder.h"

#if defined(__aarch64__)
/* Defined by src/bandedSWA.rowalt.o (bandedSWA.cpp compiled with
 * -DKERNEL_VARIANT=_rowalt -DBSW8_ROW_LEAN_INVERT=1), which the aarch64 unit
 * binary links: the kernel built with the other BSW8_ROW_LEAN setting. */
extern "C" IBandedPairWiseSW *make_bsw_kernel_rowalt(int, int, int, int, int, int,
                                                     const int8_t *, int8_t, int8_t, int);
#endif

namespace bwa_tests {

/// Score `pairs` as one batch through getScores8 (`width` 8) or getScores16
/// (`width` 16) with band `w`; returns each pair's six result fields in input
/// order. The SeqPair array has 64 slots of padding-lane slack (the widest tier's
/// lane count; the kernel writes the pads) and both sequence buffers carry 4 KiB
/// of zero tail. Throws if the kernel hands back a pair id out of range or twice.
inline std::vector<ExtResult> score_bsw_batch(IBandedPairWiseSW &bsw, const std::vector<const ExtPair *> &pairs,
                                              int width, int w)
{
    const int n = (int)pairs.size();
    std::vector<uint8_t> ref, qer;
    std::vector<SeqPair> sp(n + 64);
    for (int q = 0; q < n; q++) {
        SeqPair p;
        memset(&p, 0, sizeof p);
        p.id = p.seqid = p.regid = q;
        p.idr = (int)ref.size();
        p.idq = (int)qer.size();
        p.len1 = (int)pairs[q]->target.size();
        p.len2 = (int)pairs[q]->query.size();
        p.h0 = pairs[q]->h0;
        p.score = p.tle = p.gtle = p.qle = p.gscore = p.max_off = -1;
        ref.insert(ref.end(), pairs[q]->target.begin(), pairs[q]->target.end());
        qer.insert(qer.end(), pairs[q]->query.begin(), pairs[q]->query.end());
        sp[q] = p;
    }
    ref.resize(ref.size() + 4096, 0);
    qer.resize(qer.size() + 4096, 0);
    if (width == 8) bsw.getScores8(sp.data(), ref.data(), qer.data(), n, 1, w);
    else            bsw.getScores16(sp.data(), ref.data(), qer.data(), n, 1, w);
    std::vector<ExtResult> out(n);
    std::vector<char> seen(n, 0);
    for (int q = 0; q < n; q++) {
        const SeqPair &p = sp[q];
        if (p.id < 0 || p.id >= n || seen[p.id])   // every pair's result read exactly once
            throw std::logic_error("getScores" + std::to_string(width) + " returned pair id " +
                                   std::to_string(p.id) + " at slot " + std::to_string(q) + " of " +
                                   std::to_string(n));
        seen[p.id] = 1;
        ExtResult &o = out[p.id];
        o.score = p.score; o.qle = p.qle; o.tle = p.tle; o.gtle = p.gtle;
        o.gscore = p.gscore; o.max_off = p.max_off;
    }
    return out;
}

} // namespace bwa_tests

#endif
