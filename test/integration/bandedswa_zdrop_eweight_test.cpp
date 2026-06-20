// Regression test for z-drop gap-extension weighting in the 8-bit banded-SW
// vectorized epilogue (finding: the SIMD z-drop term was |Δrow-Δcol| with an
// implicit gap-extend penalty of 1, while the scalar reference weights the
// drift by e_del / e_ins). At the default `-E 1` the two are identical, which
// is why the existing byte-identity tests (all default-params) never caught
// it. This probe exercises NON-DEFAULT gap-extend so the weighting matters.
//
// It mirrors the high-zdrop seed test's fixture (tandem-repeat targets + noisy
// query producing real indels/drift) but uses:
//   - a SMALL zdrop, so z-drop early-exit actually FIRES for many pairs, and
//   - ASYMMETRIC, non-unit gap-extend (e_del != e_ins, both > 1), so the
//     deletion vs insertion branch of the weighted term is exercised.
// getScores8 must remain byte-identical to the scalar oracle for every pair the
// routing gate (bsw8_envelope_ok, mirrored here) would admit. Deterministic
// (fixed RNG seed); exits non-zero on any score-field mismatch.

#include <cstdio>
#include <cstdint>
#include <vector>
#include <random>

#include "bandedSWA.h"

namespace {

struct Out { int score, tle, gtle, qle, gscore, max_off; };

// Mirror of bwamem.cpp bsw8_envelope_ok for maxStep = 1 (this test's scoring,
// a = 1). The gate is independent of the gap penalties, so it is unchanged from
// the high-zdrop seed test. Only admitted pairs are required to match scalar.
bool envelope_ok(int len1, int len2, int w, int zdrop, int h0, int maxStep) {
    int shorter = len1 < len2 ? len1 : len2;
    return len1 < MAX_SEQ_LEN8 && len2 < MAX_SEQ_LEN8 && len1 >= len2 &&
           w <= 127 && zdrop + maxStep <= 127 && h0 <= zdrop + 1 &&
           (h0 + shorter * maxStep) < 255 - maxStep;   // mirrors bwamem.cpp: h0 + shorter*score_a
}

// Run one (o_del, e_del, o_ins, e_ins, zdrop) configuration: build the fixture,
// score with the scalar oracle and getScores8, diff the output fields over
// admitted pairs (gtle only when scalar gscore > 0). Returns the number of
// mismatches; sets *admitted_out to the count of envelope-admitted pairs and
// *gpos_out to how many of those had gscore > 0 (the regime where gtle is
// compared — used by main() to fail if that coverage is ever vacuous).
long run_config(const char *name, int o_del, int e_del, int o_ins, int e_ins,
                int zdrop, long *admitted_out, long *gpos_out) {
    // Multiple of 64 = max SIMD_WIDTH8 (avx512bw). getScores8 processes pairs in
    // SIMD_WIDTH8-wide batches; a non-multiple leaves a partial final batch whose
    // padding lanes read/write past the array (heap corruption on avx2/avx512).
    const long n      = 49920;
    const int  maxlen = 110;
    const int  a = 1, b = 4, ambig = -1, end_bonus = 5, w = 100;
    const int  maxStep = a > 1 ? a : 1;
    const int  STRIDE = 1280;
    const int  h0min = 1, h0max = zdrop + 1;  // admissible seed range

    int8_t mat[25];
    { int k = 0;
      for (int i = 0; i < 4; ++i) { for (int j = 0; j < 4; ++j) mat[k++] = (i == j) ? a : -b; mat[k++] = ambig; }
      for (int j = 0; j < 5; ++j) mat[k++] = ambig; }

    BandedPairWiseSW bsw(o_del, e_del, o_ins, e_ins, zdrop, end_bonus, mat, a, b, 1);

    std::mt19937_64 rng(0xC0FFEE ^ (uint64_t)zdrop ^ ((uint64_t)e_del << 8) ^ ((uint64_t)e_ins << 16));
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
        for (int i = 0; i < len1; i++) s1[i] = ubuf[i % u];     // tandem repeat
        int ti = 0;
        for (int i = 0; i < len2; i++) {                        // query: derived w/ noise
            uint8_t base = (ti < len1) ? s1[ti] : (uint8_t)(rng() % 4);
            if ((int)(rng() % 100) < 12) base = (uint8_t)(rng() % 4);  // mismatch noise
            s2[i] = base;
            if ((rng() % 20) == 0) { if (rng() & 1) ti += 2; } else ti++;  // indel noise
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

    bsw.getScores8(pairs.data(), ref.data(), qer.data(), (int32_t)n, 1, w);

    long admitted = 0, diffs = 0, gpos = 0;
    for (long c = 0; c < n; ++c) {
        const SeqPair &q = pairs[c];
        if (!envelope_ok(q.len1, q.len2, w, zdrop, q.h0, maxStep)) continue;
        admitted++;
        const Out &O = oracle[c];
        if (O.gscore > 0) gpos++;   // pairs in the regime where gtle IS compared
        // gtle (gscore row = max_ie+1) CONTRACT: it is byte-identical to scalar
        // whenever gscore > 0, and may differ only in the gscore == 0 query-end
        // tail. Root cause: the 8-bit vector row loop ends at the e-dependent band
        // bound mlenw = min(qlen+myband, tlen), while scalar runs to its dynamic
        // m==0 break and keeps updating max_ie on the trailing end==qlen rows
        // (all h1==0, so gscore stays 0). bwa-mem only consumes gtle under
        // `gscore > 0` (every gtle read in bwamem.cpp is in the else-branch of an
        // `if (gscore <= 0 || ...) {...} else { ...gtle...}` guard), so the
        // gscore==0 divergence never reaches a SAM record. We ENFORCE the safe
        // half of that contract here: gtle must match for gscore > 0. If a future
        // change widens the divergence into the gscore > 0 regime (the consumed
        // one), this fails. The gscore == 0 tail is intentionally not compared.
        const bool gtle_must_match = (O.gscore > 0);
        if (O.score != q.score || O.tle != q.tle ||
            O.qle != q.qle || O.gscore != q.gscore || O.max_off != q.max_off ||
            (gtle_must_match && O.gtle != q.gtle)) {
            if (diffs < 10)
                fprintf(stderr,
                    "[zdrop-eweight:%s] DIFF c=%ld len1=%d len2=%d h0=%d | scalar %d/%d/%d/%d/%d/%d | 8 %d/%d/%d/%d/%d/%d\n",
                    name, c, q.len1, q.len2, q.h0,
                    O.score, O.tle, O.gtle, O.qle, O.gscore, O.max_off,
                    q.score, q.tle, q.gtle, q.qle, q.gscore, q.max_off);
            diffs++;
        }
    }
    fprintf(stderr, "[zdrop-eweight:%s] o_del=%d e_del=%d o_ins=%d e_ins=%d zdrop=%d: admitted=%ld (gscore>0: %ld) diffs=%ld\n",
            name, o_del, e_del, o_ins, e_ins, zdrop, admitted, gpos, diffs);
    *admitted_out = admitted;
    *gpos_out = gpos;
    return diffs;
}

} // namespace

int main() {
    long adm = 0, gpos = 0, total_diffs = 0, min_admitted = -1, total_gpos = 0;

    // Config A: asymmetric, non-unit gap-extend, small zdrop (z-drop fires
    // often, deletion/insertion weighting both exercised). This is the config
    // that diverges pre-fix.
    total_diffs += run_config("asym", /*o_del*/6, /*e_del*/3, /*o_ins*/6, /*e_ins*/1, /*zdrop*/25, &adm, &gpos);
    if (min_admitted < 0 || adm < min_admitted) min_admitted = adm;
    total_gpos += gpos;

    // Config B: symmetric non-unit gap-extend.
    total_diffs += run_config("sym2", 6, 2, 6, 2, 30, &adm, &gpos);
    if (adm < min_admitted) min_admitted = adm;
    total_gpos += gpos;

    // Config C: default gap-extend (e = 1). Must match both before and after
    // the fix — a guard that the weighting change is a no-op at default params.
    total_diffs += run_config("e1", 6, 1, 6, 1, 25, &adm, &gpos);
    if (adm < min_admitted) min_admitted = adm;
    total_gpos += gpos;

    if (min_admitted <= 0) {
        fprintf(stderr, "bandedswa_zdrop_eweight_test: FAIL — a config admitted no pairs\n");
        return 2;
    }
    // The gtle invariant (compared iff gscore > 0) is the point of this probe;
    // fail loudly if the fixture never exercises a gscore > 0 pair, else the
    // gtle check is silently vacuous (cf. seed_hi guard in the high-zdrop test).
    if (total_gpos <= 0) {
        fprintf(stderr, "bandedswa_zdrop_eweight_test: FAIL — no gscore>0 pairs; gtle invariant untested\n");
        return 2;
    }
    if (total_diffs != 0) {
        fprintf(stderr, "bandedswa_zdrop_eweight_test: FAIL — %ld getScores8/scalar mismatches\n", total_diffs);
        return 1;
    }
    fprintf(stderr, "bandedswa_zdrop_eweight_test: OK\n");
    return 0;
}
