// Regression test: an 8-bit banded-SW lane with a NEGATIVE seed score h0 must
// score the same whatever other pairs share its SIMD group.
//
// Under --meth a seed is rescored in original space, so the extension seed
// score h0 can be negative. The 8-bit kernels clamp the byte seed to 0 but keep
// the raw h0 in the wide per-lane best_abs side channel. Their per-row wide
// update ran best_abs = max(best_abs, byte max) over every lane of the group
// whenever ANY lane's max advanced that row, on the assumption that best_abs is
// never below the byte max. A negative h0 breaks that assumption: a lane whose
// cells never score above 0 was lifted from h0 to 0 exactly when a neighbouring
// lane advanced, so its score depended on which pairs were batched with it.
// Because the extension-DP dedup (default `auto`, latched ON/OFF from measured
// wall time) changes that grouping, --meth output differed from run to run.
//
// Checks, for every tier the host builds natively (-march=native):
//   1. a "dead" negative-h0 pair (all cells <= 0) scores its h0 alone, next to
//      an advancing pair in one group, and among many pairs, and that equals
//      the scalar oracle;
//   2. random mixed batches (negative and positive h0) give every pair the same
//      result as that pair scored alone and as the scalar oracle.
// Exits non-zero on any mismatch.

#include <cstdio>
#include <cstdint>
#include <random>
#include <vector>

#include "bandedSWA.h"

namespace {

struct Pair { int h0; std::vector<uint8_t> ref, qer; };
struct Out { int score, tle, qle, gscore, gtle, max_off; };

bool same(const Out &a, const Out &b)
{
    return a.score == b.score && a.tle == b.tle && a.qle == b.qle && a.gscore == b.gscore
        && a.gtle == b.gtle && a.max_off == b.max_off;
}

// Score `pairs` in one batch (vector or scalar); results in input order.
std::vector<Out> score(BandedPairWiseSW &bsw, const std::vector<Pair> &pairs, bool scalar, int w)
{
    std::vector<uint8_t> ref, qer;
    const int n = (int)pairs.size();
    std::vector<SeqPair> sp(n + 256);   // padding-lane slack, as production allocates
    for (int i = 0; i < n; i++) {
        SeqPair p = {};
        p.id = i; p.seqid = i; p.regid = 0;
        p.idr = (int)ref.size(); p.idq = (int)qer.size();
        p.len1 = (int)pairs[i].ref.size(); p.len2 = (int)pairs[i].qer.size();
        p.h0 = pairs[i].h0;
        p.score = p.tle = p.gtle = p.qle = p.gscore = p.max_off = -1;
        ref.insert(ref.end(), pairs[i].ref.begin(), pairs[i].ref.end());
        qer.insert(qer.end(), pairs[i].qer.begin(), pairs[i].qer.end());
        sp[i] = p;
    }
    ref.resize(ref.size() + 256);
    qer.resize(qer.size() + 256);
    if (scalar) bsw.scalarBandedSWAWrapper(sp.data(), ref.data(), qer.data(), n, 1, w);
    else        bsw.getScores8(sp.data(), ref.data(), qer.data(), n, 1, w);
    std::vector<Out> out(n);
    for (int i = 0; i < n; i++) {
        const SeqPair &p = sp[i];
        out[p.id] = Out{p.score, p.tle, p.qle, p.gscore, p.gtle, p.max_off};
    }
    return out;
}

Pair make_pair(std::mt19937 &rng, int h0, int len1, int len2, bool copy)
{
    Pair p;
    p.h0 = h0;
    p.ref.resize(len1);
    p.qer.resize(len2);
    for (auto &b : p.ref) b = (uint8_t)(rng() & 3);
    for (int j = 0; j < len2; j++) p.qer[j] = copy && j < len1 ? p.ref[j] : (uint8_t)(rng() & 3);
    return p;
}

}  // namespace

int main()
{
    const int8_t a = 1, b = 4, ambig = -1;
    const int o = 6, e = 1, zdrop = 100, end_bonus = 5, w = 100;
    int8_t mat[25];
    { int k = 0;
      for (int i = 0; i < 4; ++i) { for (int j = 0; j < 4; ++j) mat[k++] = (i == j) ? a : -b; mat[k++] = ambig; }
      for (int j = 0; j < 5; ++j) mat[k++] = ambig; }
    BandedPairWiseSW bsw(o, e, o, e, zdrop, end_bonus, mat, a, b, 1);

    int fails = 0;
    std::mt19937 rng(20260926);

    // 1. A dead negative-h0 pair: reference all A, query all C, so every cell
    //    scores <= 0 and its maximum never advances past the clamped seed.
    Pair dead;
    dead.h0 = -5;
    dead.ref.assign(30, 0);
    dead.qer.assign(20, 1);
    const Out want = score(bsw, {dead}, true, w)[0];
    if (want.score != dead.h0) {
        fprintf(stderr, "scalar oracle: dead pair scored %d, expected its h0 %d\n", want.score, dead.h0);
        fails++;
    }
    // Alone, next to one advancing pair, and inside many-pair batches that put
    // it in a group whose other lanes advance on every row.
    const Pair live = make_pair(rng, 20, 40, 40, true);
    std::vector<std::vector<Pair>> batches = {{dead}, {live, dead}, {dead, live}};
    for (int n : {31, 63, 64, 65, 130}) {
        std::vector<Pair> v;
        for (int i = 0; i < n; i++) v.push_back(i == n / 2 ? dead : make_pair(rng, 10 + (int)(rng() % 30), 40, 30, true));
        batches.push_back(v);
    }
    for (size_t t = 0; t < batches.size(); t++) {
        const std::vector<Out> got = score(bsw, batches[t], false, w);
        for (size_t i = 0; i < batches[t].size(); i++) {
            if (batches[t][i].h0 != dead.h0) continue;
            if (!same(got[i], want)) {
                fprintf(stderr, "batch %zu (n=%zu): dead pair got score=%d tle=%d qle=%d gscore=%d gtle=%d max_off=%d, "
                        "scalar score=%d tle=%d qle=%d gscore=%d gtle=%d max_off=%d\n", t, batches[t].size(),
                        got[i].score, got[i].tle, got[i].qle, got[i].gscore, got[i].gtle, got[i].max_off,
                        want.score, want.tle, want.qle, want.gscore, want.gtle, want.max_off);
                fails++;
            }
        }
    }

    // 2. Random mixed batches: every pair must score as it does alone and as the
    //    scalar oracle does.
    for (int round = 0; round < 20; round++) {
        const int n = 1 + (int)(rng() % 200);
        std::vector<Pair> v;
        for (int i = 0; i < n; i++) {
            const int len1 = 5 + (int)(rng() % 60), len2 = 5 + (int)(rng() % 40);
            const int h0 = (rng() % 3 == 0) ? -(int)(1 + rng() % 30) : (int)(1 + rng() % 40);
            v.push_back(make_pair(rng, h0, len1, len2, rng() % 2 == 0));
        }
        const std::vector<Out> got = score(bsw, v, false, w);
        const std::vector<Out> oracle = score(bsw, v, true, w);
        for (int i = 0; i < n; i++) {
            const Out alone = score(bsw, {v[i]}, false, w)[0];
            if (!same(got[i], oracle[i]) || !same(got[i], alone)) {
                fprintf(stderr, "round %d pair %d (h0=%d len1=%zu len2=%zu): batch score=%d gscore=%d tle=%d qle=%d, "
                        "alone score=%d gscore=%d tle=%d qle=%d, scalar score=%d gscore=%d tle=%d qle=%d\n",
                        round, i, v[i].h0, v[i].ref.size(), v[i].qer.size(),
                        got[i].score, got[i].gscore, got[i].tle, got[i].qle, alone.score, alone.gscore, alone.tle, alone.qle,
                        oracle[i].score, oracle[i].gscore, oracle[i].tle, oracle[i].qle);
                if (++fails > 20) return 1;
            }
        }
    }

    if (fails) { fprintf(stderr, "bandedswa_negative_h0_test: %d mismatches\n", fails); return 1; }
    printf("bandedswa_negative_h0_test: OK\n");
    return 0;
}
