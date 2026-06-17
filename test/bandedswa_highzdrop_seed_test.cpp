// Regression test for the unsigned [0,255] h0-prefix seed (issue #146 item:
// "make the 8-bit h0-prefix seed unsigned").
//
// The 8-bit banded-SW DP body is unsigned [0,255]. The h0-prefix column/row
// seed (smithWaterman*_8 wrapper setup) is now also unsigned-saturating
// (subs_epu8); it previously used signed int8 ops, which required the seeded
// byte min(h0, REBASE_KEEP) <= 127 and so capped zdrop at 126
// (BSW8_MAX_ZDROP_STEP). With the unsigned seed the gate admits zdrop up to
// ~252, and seed prefix bytes range over (127, 254].
//
// This test exercises exactly that new regime: zdrop = 200 with h0 in
// (127, 200], so every admitted pair seeds prefix bytes > 127 — the range the
// old signed seed could not represent. getScores8 must remain byte-identical to
// the scalar oracle for every pair the routing gate (bsw8_envelope_ok) would
// admit (re-baseline provably inert), so we mirror that gate here, compact the
// admitted pairs into a tight batch, and score only that batch -- out-of-envelope
// pairs are documented as requiring 16-bit routing and never reach the 8-bit
// kernel. Deterministic (fixed RNG seed); exits non-zero on any score-field
// mismatch.

#include <cstdio>
#include <cstdint>
#include <vector>
#include <random>

#include "bandedSWA.h"

namespace {

struct Out { int score, tle, gtle, qle, gscore, max_off; };

// Mirror of bwamem.cpp bsw8_envelope_ok for a = maxStep = 1 (this test's
// scoring). Only admitted pairs are required to match scalar.
bool envelope_ok(int len1, int len2, int w, int zdrop, int h0, int maxStep) {
    int shorter = len1 < len2 ? len1 : len2;
    return len1 < MAX_SEQ_LEN8 && len2 < MAX_SEQ_LEN8 && len1 >= len2 &&
           w <= 127 && zdrop + maxStep <= 253 && h0 <= zdrop + 1 &&
           (h0 + shorter * 1) < 255 - maxStep;
}

} // namespace

int main() {
    const long n      = 50000;
    const int  maxlen = 110;
    const int  zdrop  = 200;        // > 126: the new regime
    const int  h0min  = 128, h0max = 200;  // seed bytes strictly > 127
    const int  a = 1, b = 4, ambig = -1, o = 6, e = 1, end_bonus = 5, w = 100;
    const int  maxStep = a > 1 ? a : 1;
    const int  STRIDE = 1280;

    int8_t mat[25];
    { int k = 0;
      for (int i = 0; i < 4; ++i) { for (int j = 0; j < 4; ++j) mat[k++] = (i == j) ? a : -b; mat[k++] = ambig; }
      for (int j = 0; j < 5; ++j) mat[k++] = ambig; }

    BandedPairWiseSW bsw(o, e, o, e, zdrop, end_bonus, mat, a, b, 1);

    std::mt19937_64 rng(12345);  // fixed seed -> deterministic
    std::vector<uint8_t> ref((size_t)STRIDE * n, 0), qer((size_t)STRIDE * n, 0);
    std::vector<SeqPair> pairs(n);
    std::vector<Out> oracle(n);
    std::uniform_int_distribution<int> lenD(20, maxlen);
    std::uniform_int_distribution<int> hD(h0min, h0max);
    std::uniform_int_distribution<int> unit(3, 12);

    for (long c = 0; c < n; ++c) {
        int aa = lenD(rng), bb = lenD(rng);
        int len1 = aa > bb ? aa : bb, len2 = aa > bb ? bb : aa;  // len1 >= len2
        int h0 = hD(rng);
        uint8_t *s1 = &ref[(size_t)c * STRIDE];
        uint8_t *s2 = &qer[(size_t)c * STRIDE];
        int u = unit(rng); uint8_t ubuf[16];
        for (int i = 0; i < u; i++) ubuf[i] = (uint8_t)(rng() % 4);
        for (int i = 0; i < len1; i++) s1[i] = ubuf[i % u];   // tandem repeat
        int ti = 0;
        for (int i = 0; i < len2; i++) {                       // query: derived w/ noise
            uint8_t base = (ti < len1) ? s1[ti] : (uint8_t)(rng() % 4);
            if ((int)(rng() % 100) < 5) base = (uint8_t)(rng() % 4);
            s2[i] = base;
            if ((rng() % 40) == 0) { if (rng() & 1) ti += 2; } else ti++;
        }
        SeqPair &p = pairs[c];
        p.id = c; p.len1 = len1; p.len2 = len2; p.h0 = h0;
        p.idr = (int)((size_t)c * STRIDE); p.idq = (int)((size_t)c * STRIDE);
        p.seqid = c; p.regid = c;
        p.score = p.tle = p.gtle = p.qle = p.gscore = p.max_off = -1;
        Out &O = oracle[c];
        O.score = bsw.scalarBandedSWA(len2, s2, len1, s1, w, h0,
                                      &O.qle, &O.tle, &O.gtle, &O.gscore, &O.max_off);
    }

    // Routing contract: only envelope-admitted pairs may traverse the 8-bit
    // kernel — out-of-envelope pairs are documented as requiring 16-bit routing,
    // so feeding them to getScores8 exercises an unspecified path. Mirror the
    // gate here and compact the admitted pairs into a tight batch before
    // scoring. Each compacted SeqPair keeps its original idr/idq, so it still
    // references the correct slice of the shared ref/qer buffers.
    std::vector<long> admittedIdx;
    admittedIdx.reserve(n);
    for (long c = 0; c < n; ++c) {
        const SeqPair &p = pairs[c];
        if (envelope_ok(p.len1, p.len2, w, zdrop, p.h0, maxStep)) admittedIdx.push_back(c);
    }
    const long admitted = (long)admittedIdx.size();

    // getScores8 pads the batch to a whole number of SIMD lanes and WRITES the
    // trailing padding lanes pairArray[numPairs .. roundup(numPairs, SIMD_WIDTH8)).
    // Size the compact batch to that rounded capacity so the kernel's padding
    // writes stay in-bounds (AVX2 SIMD_WIDTH8=32, AVX-512 =64; the scalar/NEON
    // tiers round to 1/16), then bound the kernel call and comparison to the
    // first `admitted` lanes.
    const long round_admitted = ((admitted + SIMD_WIDTH8 - 1) / SIMD_WIDTH8) * SIMD_WIDTH8;
    std::vector<SeqPair> apairs(round_admitted);
    for (long k = 0; k < admitted; ++k) apairs[k] = pairs[admittedIdx[k]];

    bsw.getScores8(apairs.data(), ref.data(), qer.data(), (int32_t)admitted, 1, w);

    long seed_hi = 0, diffs = 0;
    for (long k = 0; k < admitted; ++k) {
        const SeqPair &q = apairs[k];
        long c = admittedIdx[k];
        int h0p = q.h0 < zdrop + 1 ? q.h0 : zdrop + 1;   // seeded byte min(h0,REBASE_KEEP)
        if (h0p > 127) seed_hi++;
        const Out &O = oracle[c];
        if (O.score != q.score || O.tle != q.tle || O.gtle != q.gtle ||
            O.qle != q.qle || O.gscore != q.gscore || O.max_off != q.max_off) {
            if (diffs < 10)
                fprintf(stderr,
                    "[highzdrop-seed] DIFF c=%ld len1=%d len2=%d h0=%d | scalar %d/%d/%d/%d/%d/%d | 8 %d/%d/%d/%d/%d/%d\n",
                    c, q.len1, q.len2, q.h0,
                    O.score, O.tle, O.gtle, O.qle, O.gscore, O.max_off,
                    q.score, q.tle, q.gtle, q.qle, q.gscore, q.max_off);
            diffs++;
        }
    }

    fprintf(stderr, "[highzdrop-seed] zdrop=%d h0=[%d,%d]: admitted=%ld (seed byte >127: %ld) diffs=%ld\n",
            zdrop, h0min, h0max, admitted, seed_hi, diffs);

    if (admitted == 0 || seed_hi == 0) {
        fprintf(stderr, "bandedswa_highzdrop_seed_test: FAIL — fixture admitted no >127-seed pairs\n");
        return 2;
    }
    if (diffs != 0) {
        fprintf(stderr, "bandedswa_highzdrop_seed_test: FAIL — %ld getScores8/scalar mismatches\n", diffs);
        return 1;
    }
    fprintf(stderr, "bandedswa_highzdrop_seed_test: OK\n");
    return 0;
}
