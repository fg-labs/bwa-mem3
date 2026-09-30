// test/framework/ext_ladder.h
//
// Test-side model of bwamem.cpp's seed-extension retry ladder over a scalar
// extension kernel, for differential tests of the ladder's rung count. Two
// independent kernels are
// wrapped: bwa-mem3's scalarBandedSWA and upstream bwa's ksw_extend2 (src/ksw.cpp,
// unchanged from bwa), so a claim of "what the full ladder would commit" is
// checked against both.
//
// The ladder here is the exact one (ACCEPT_PAIR's non-cert branch with
// band_start <= 0): rung i runs at w0 << i and is accepted when it is the last
// rung, or the score equals the previous rung's, or max_off < (w>>1)+(w>>2).
// Upstream bwa and bwa-mem2 use two rungs (MAX_BAND_TRY 2).

#ifndef BWA_TESTS_EXT_LADDER_H
#define BWA_TESTS_EXT_LADDER_H

#include <cstdint>
#include <cstring>
#include <vector>

#include "bandedSWA.h"
#include "ksw.h"
#include "scoring.h"

namespace bwa_tests {

struct ExtScoring {
    int a, b, o_del, e_del, o_ins, e_ins, pen_clip, zdrop;
    int8_t mat[25];
    // bwa_fill_scmat(a, b): match a, mismatch -b, N vs anything -1.
    void fill_mat() { std::memcpy(mat, build_scoring_matrix(a, b, 1).data(), sizeof mat); }
    int o_min() const { return o_del < o_ins ? o_del : o_ins; }
    int e_min() const { return e_del < e_ins ? e_del : e_ins; }
};

struct ExtResult {
    int score, qle, tle, gtle, gscore, max_off;
    bool operator==(const ExtResult &o) const {
        return score == o.score && qle == o.qle && tle == o.tle && gtle == o.gtle &&
               gscore == o.gscore && max_off == o.max_off;
    }
};

// One extension pair as the ladder sees it: query (len2 = qlen), target (len1 =
// tlen >= qlen), seed score h0.
struct ExtPair {
    std::vector<uint8_t> query, target;
    int h0;
};

// bwa-mem3's scalar kernel. A fresh BandedPairWiseSW per call keeps the test
// independent of the object's scratch state.
inline ExtResult run_mem3_scalar(const ExtScoring &sc, const ExtPair &p, int w) {
    BandedPairWiseSW bsw(sc.o_del, sc.e_del, sc.o_ins, sc.e_ins, sc.zdrop, sc.pen_clip,
                         sc.mat, (int8_t)sc.a, (int8_t)sc.b, 1);
    ExtResult r;
    r.score = bsw.scalarBandedSWA((int)p.query.size(), p.query.data(), (int)p.target.size(),
                                  p.target.data(), w, p.h0, &r.qle, &r.tle, &r.gtle, &r.gscore,
                                  &r.max_off);
    return r;
}

// Upstream bwa's kernel (src/ksw.cpp ksw_extend2, verbatim from bwa).
inline ExtResult run_upstream_ksw(const ExtScoring &sc, const ExtPair &p, int w) {
    ExtResult r;
    r.score = ksw_extend2((int)p.query.size(), p.query.data(), (int)p.target.size(),
                          p.target.data(), 5, sc.mat, sc.o_del, sc.e_del, sc.o_ins, sc.e_ins, w,
                          sc.pen_clip, sc.zdrop, p.h0, &r.qle, &r.tle, &r.gtle, &r.gscore,
                          &r.max_off);
    return r;
}

typedef ExtResult (*ExtKernel)(const ExtScoring &, const ExtPair &, int);

struct LadderOutcome {
    ExtResult res;   // the accepted rung's result
    int rung;        // index of the accepted rung
    int w;           // its band width; the region records a->w = max(opt->w, w)
};

// Run the exact ladder [w0, 2*w0, ...] with nband rungs. prev0 is the a->score
// the ladder starts from: -1 for a left extension (the region's initial score),
// the left result (== h0) for a right extension.
inline LadderOutcome run_ladder(ExtKernel k, const ExtScoring &sc, const ExtPair &p, int w0,
                                int nband, int prev0) {
    LadderOutcome out;
    int prev = prev0;
    for (int i = 0; i < nband; ++i) {
        const int w = w0 << i;
        ExtResult r = k(sc, p, w);
        out.res = r; out.rung = i; out.w = w;
        if (i + 1 == nband || r.score == prev || r.max_off < ((w >> 1) + (w >> 2))) break;
        prev = r.score;
    }
    return out;
}

} // namespace bwa_tests

#endif
