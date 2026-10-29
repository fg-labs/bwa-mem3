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
//      result as that pair scored alone and as the scalar oracle;
//   3. with a match score of 1, 2 and 3, a pair routed as production routes it
//      (bsw8_envelope_ok picks getScores8 or getScores16) matches the scalar
//      oracle for h0 from -(a+1) to 0. The 8-bit kernels clamp the seed byte at
//      0, which matches scalar only while h0 + a <= 0 on the seed's diagonal, so
//      the gate must keep a negative h0 off the 8-bit kernel;
//   4. the production classifier (sortPairsLenExt) sends a negative-h0 pair
//      whose h0 + min(len1,len2)*a is negative to the scalar tier, not 16-bit.
// Exits non-zero on any mismatch.

#include <cstdio>
#include <cstdint>
#include <random>
#include <vector>

#include "bandedSWA.h"
#include "bwamem.h"
#include "framework/ext_result.h"

namespace {

struct Pair { int h0; std::vector<uint8_t> ref, qer; };
using bwa_tests::ExtResult;

enum Kernel { SCALAR, K8, K16 };

// Score `pairs` in one batch with `kernel`; results in input order.
std::vector<ExtResult> score_with(BandedPairWiseSW &bsw, const std::vector<Pair> &pairs, Kernel kernel, int w)
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
    if (kernel == SCALAR)  bsw.scalarBandedSWAWrapper(sp.data(), ref.data(), qer.data(), n, 1, w);
    else if (kernel == K8) bsw.getScores8(sp.data(), ref.data(), qer.data(), n, 1, w);
    else                   bsw.getScores16(sp.data(), ref.data(), qer.data(), n, 1, w);
    std::vector<ExtResult> out(n);
    for (int i = 0; i < n; i++) {
        const SeqPair &p = sp[i];
        ExtResult &o = out[p.id];
        o.score = p.score; o.tle = p.tle; o.qle = p.qle;
        o.gscore = p.gscore; o.gtle = p.gtle; o.max_off = p.max_off;
    }
    return out;
}

std::vector<ExtResult> score(BandedPairWiseSW &bsw, const std::vector<Pair> &pairs, bool scalar, int w)
{
    return score_with(bsw, pairs, scalar ? SCALAR : K8, w);
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
    const ExtResult want = score(bsw, {dead}, true, w)[0];
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
        const std::vector<ExtResult> got = score(bsw, batches[t], false, w);
        for (size_t i = 0; i < batches[t].size(); i++) {
            if (batches[t][i].h0 != dead.h0) continue;
            if (!(got[i] == want)) {
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
        const std::vector<ExtResult> got = score(bsw, v, false, w);
        const std::vector<ExtResult> oracle = score(bsw, v, true, w);
        for (int i = 0; i < n; i++) {
            const ExtResult alone = score(bsw, {v[i]}, false, w)[0];
            if (!(got[i] == oracle[i]) || !(got[i] == alone)) {
                fprintf(stderr, "round %d pair %d (h0=%d len1=%zu len2=%zu): batch score=%d gscore=%d tle=%d qle=%d, "
                        "alone score=%d gscore=%d tle=%d qle=%d, scalar score=%d gscore=%d tle=%d qle=%d\n",
                        round, i, v[i].h0, v[i].ref.size(), v[i].qer.size(),
                        got[i].score, got[i].gscore, got[i].tle, got[i].qle, alone.score, alone.gscore, alone.tle, alone.qle,
                        oracle[i].score, oracle[i].gscore, oracle[i].tle, oracle[i].qle);
                if (++fails > 20) return 1;
            }
        }
    }

    // 3. Routed through the production gate at a = 1, 2, 3. Each query copies the
    //    start of its reference, so the seed's diagonal scores a per base and an
    //    h0 in (-a, 0) can start a positive extension.
    for (int8_t ra = 1; ra <= 3; ra++) {
        int8_t rmat[25];
        { int k = 0;
          for (int i = 0; i < 4; ++i) { for (int j = 0; j < 4; ++j) rmat[k++] = (i == j) ? ra : -b; rmat[k++] = ambig; }
          for (int j = 0; j < 5; ++j) rmat[k++] = ambig; }
        BandedPairWiseSW rbsw(o, e, o, e, zdrop, end_bonus, rmat, ra, b, 1);
        for (int h0 = -(ra + 1); h0 <= 0; h0++) {
            std::vector<Pair> v8, v16;
            for (int i = 0; i < 100; i++) {
                const int len2 = 20 + (int)(rng() % 20), len1 = len2 + (int)(rng() % 20);
                Pair p = make_pair(rng, h0, len1, len2, true);
                (bsw8_envelope_ok(len1, len2, w, ra, zdrop, h0) ? v8 : v16).push_back(p);
            }
            // The gate rejects every negative seed by design, but must admit the
            // h0 = 0 control, or the check below never reaches getScores8.
            if (h0 == 0 && v8.empty()) {
                fprintf(stderr, "routed a=%d h0=0: no pair reached getScores8\n", ra);
                fails++;
            }
            for (const std::vector<Pair> *v : {&v8, &v16}) {
                if (v->empty()) continue;
                const Kernel k = v == &v8 ? K8 : K16;
                const std::vector<ExtResult> got = score_with(rbsw, *v, k, w);
                const std::vector<ExtResult> oracle = score_with(rbsw, *v, SCALAR, w);
                for (size_t i = 0; i < v->size(); i++) {
                    if (got[i] == oracle[i]) continue;
                    fprintf(stderr, "routed a=%d h0=%d pair %zu (len1=%zu len2=%zu) on %s: score=%d qle=%d tle=%d gscore=%d, "
                            "scalar score=%d qle=%d tle=%d gscore=%d\n", ra, h0, i, (*v)[i].ref.size(), (*v)[i].qer.size(),
                            k == K8 ? "getScores8" : "getScores16", got[i].score, got[i].qle, got[i].tle, got[i].gscore,
                            oracle[i].score, oracle[i].qle, oracle[i].tle, oracle[i].gscore);
                    if (++fails > 20) return 1;
                }
            }
        }
    }

    // 4. The production classifier (sortPairsLenExt) at a = 1, len1 = len2 = 20:
    //    h0 = 0 takes the 8-bit tier; h0 = -1 (minval 19) leaves it for 16-bit;
    //    h0 = -30 (minval -10) fails the 16-bit range check (minval < 0), so goes to the
    //    scalar tier.
    {
        const int ra = 1, len = 20;
        const int h0s[3] = {0, -1, -30};
        std::vector<SeqPair> sp(3), tmp(3);
        for (int i = 0; i < 3; i++) {
            SeqPair p = {};
            p.id = i; p.seqid = i;
            p.len1 = p.len2 = len;
            p.h0 = h0s[i];
            sp[i] = p;
        }
        int n8 = -1, n16 = -1, n1 = -1;
        sortPairsLenExt(sp.data(), 3, tmp.data(), n8, n16, n1, ra, w, zdrop);
        if (n8 != 1 || n16 != 1 || n1 != 1) {
            fprintf(stderr, "classifier: tier counts 8-bit=%d 16-bit=%d scalar=%d, expected 1/1/1\n", n8, n16, n1);
            fails++;
        }
        // Tiers are laid out 8-bit, then 16-bit, then scalar.
        for (int i = 0; i < 3; i++) {
            if (sp[i].h0 == h0s[i]) continue;
            fprintf(stderr, "classifier: slot %d holds h0=%d, expected %d\n", i, sp[i].h0, h0s[i]);
            fails++;
        }
    }

    if (fails) { fprintf(stderr, "bandedswa_negative_h0_test: %d mismatches\n", fails); return 1; }
    printf("bandedswa_negative_h0_test: OK\n");
    return 0;
}
