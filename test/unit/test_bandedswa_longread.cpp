// test/unit/test_bandedswa_longread.cpp
//
// Byte-identity test for the banded Smith-Waterman batched SIMD kernels
// (BandedPairWiseSW::getScores8 / getScores16) versus the scalar oracle
// (scalarBandedSWA), across SHORT and LONG reads.
//
// Locks in the recovered long-read 8-bit path (reads >=128bp): every one of the
// six outputs (score, tle, gtle, qle, gscore, max_off) must match the scalar
// reference on the host SIMD tier. Pairs are generated with target >= query
// (len1 >= len2) -- the saturation-safe domain the 8-bit routing envelope
// targets. Runs on every CI matrix row that defines a vector kernel.

#include <cstdint>
#include <random>
#include <vector>

#include "doctest/doctest.h"
#include "bandedSWA.h"

#if HAVE_BSW_VECTOR_8_16

namespace {

struct Out { int score, tle, gtle, qle, gscore, max_off; };

// The batched kernels (getScores8/16) round numPairs up to a multiple of the
// SIMD width (roundNumPairs in smithWatermanBatchWrapper{8,16}) and initialize
// + prefetch-read padding lanes PAST numPairs (the production pipeline
// over-allocates seqPairArray by MAX_LINE_LEN to cover this). A direct caller
// that allocates exactly numPairs SeqPairs therefore gets an out-of-bounds
// write/read on any tier whose SIMD width does not divide numPairs (e.g.
// AVX2=32, AVX-512=64). Pad the SeqPair array: round up to SIMD_WIDTH8 -- the
// 8-bit lane count of the tier this TU compiles for, which equals the width the
// linked kernel rounds to because the test build is given the same arch macros
// as libbwa.a (see test/Makefile ARCH_FLAGS_FROM_PARENT) -- plus MAX_LINE_LEN
// slack for the prefetch look-ahead, mirroring the production over-allocation.
// Padding pairs are value-initialized (len1=len2=0), inert dummy lanes.
inline int padPairs(int n) {
    return ((n + SIMD_WIDTH8 - 1) / SIMD_WIDTH8) * SIMD_WIDTH8 + MAX_LINE_LEN;
}

// Build the bench-default 5x5 nucleotide scoring matrix (match=a, mismatch=-b,
// ambiguous=-1).
void build_mat(int8_t mat[25], int a, int b, int ambig) {
    int k = 0;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) mat[k++] = (i == j) ? (int8_t)a : (int8_t)(-b);
        mat[k++] = (int8_t)ambig;
    }
    for (int j = 0; j < 5; ++j) mat[k++] = (int8_t)ambig;
}

// Drive width-bit batched kernel (8 or 16) vs the scalar oracle over n random
// target>=query pairs up to maxlen, and tally per-field mismatches. Returns the
// total number of pairs with any mismatch.
int run_parity(int width, int n, int maxlen, unsigned long seed, int zdrop,
               int &bs, int &bt, int &bg, int &bq, int &bgs, int &bm) {
    const int a = 1, b = 4, ambig = -1, o = 6, e = 1, end_bonus = 5, w = 100;
    const int STRIDE = maxlen + MAX_LINE_LEN;   // per-pair seq slot; scales with read
                                                // length, + slack for prefetch look-ahead
                                                // (16-bit handles up to MAX_SEQ_LEN16=32768)
    int8_t mat[25];
    build_mat(mat, a, b, ambig);

    BandedPairWiseSW bsw(o, e, o, e, zdrop, end_bonus, mat, a, b, 1);

    std::vector<uint8_t> seqBufRef((size_t)STRIDE * n, 0);
    std::vector<uint8_t> seqBufQer((size_t)STRIDE * n, 0);
    std::vector<SeqPair> pairs(padPairs(n));   // padded for SIMD-width round-up + prefetch
    std::vector<Out> oracle(n);

    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<int> lenD(5, maxlen);
    for (int c = 0; c < n; ++c) {
        // target (len1) >= query (len2): the 8-bit routing envelope's domain.
        int aa = lenD(rng), bb = lenD(rng);
        int len1 = aa > bb ? aa : bb, len2 = aa > bb ? bb : aa;
        int h0 = (rng() & 1) ? 19 : (int)(rng() % 100 + 1);
        uint8_t *s1 = &seqBufRef[(size_t)c * STRIDE];
        uint8_t *s2 = &seqBufQer[(size_t)c * STRIDE];
        for (int i = 0; i < len1; ++i) s1[i] = (uint8_t)(rng() % 4);
        // mutated copy of the target prefix -> realistic (not random) alignments
        for (int i = 0, ti = 0; i < len2; ++i, ++ti) {
            uint8_t base = (ti < len1) ? s1[ti] : (uint8_t)(rng() % 4);
            if ((int)(rng() % 100) < 5) base = (uint8_t)(rng() % 4);   // 5% mismatch
            s2[i] = base;
        }
        SeqPair &sp = pairs[c];
        sp.id = c; sp.len1 = len1; sp.len2 = len2; sp.h0 = h0;
        sp.idr = (int)((size_t)c * STRIDE); sp.idq = (int)((size_t)c * STRIDE);
        sp.seqid = c; sp.regid = c;
        sp.score = sp.tle = sp.gtle = sp.qle = sp.gscore = sp.max_off = -1;
        Out &O = oracle[c];
        O.score = bsw.scalarBandedSWA(len2, s2, len1, s1, w, h0,
                                      &O.qle, &O.tle, &O.gtle, &O.gscore, &O.max_off);
    }

    if (width == 8)
        bsw.getScores8(pairs.data(), seqBufRef.data(), seqBufQer.data(), (int32_t)n, 1, w);
    else
        bsw.getScores16(pairs.data(), seqBufRef.data(), seqBufQer.data(), (int32_t)n, 1, w);

    int bad = 0;
    bs = bt = bg = bq = bgs = bm = 0;
    for (int c = 0; c < n; ++c) {
        const Out &O = oracle[c]; const SeqPair &p = pairs[c];
        bool sd = O.score != p.score, td = O.tle != p.tle, gd = O.gtle != p.gtle,
             qd = O.qle != p.qle, gsd = O.gscore != p.gscore, md = O.max_off != p.max_off;
        if (sd) bs++; if (td) bt++; if (gd) bg++; if (qd) bq++; if (gsd) bgs++; if (md) bm++;
        if (sd || td || gd || qd || gsd || md) bad++;
    }
    return bad;
}

void check_width(int width, int n, int maxlen, unsigned long seed, int zdrop = 100) {
    int bs, bt, bg, bq, bgs, bm;
    int bad = run_parity(width, n, maxlen, seed, zdrop, bs, bt, bg, bq, bgs, bm);
    MESSAGE("bandedSWA getScores" << width << " vs scalar: maxlen=" << maxlen
            << " zdrop=" << zdrop
            << " n=" << n << " ANY=" << bad
            << " (score=" << bs << " tle=" << bt << " gtle=" << bg
            << " qle=" << bq << " gscore=" << bgs << " max_off=" << bm << ")");
    CHECK(bs == 0);
    CHECK(bt == 0);
    CHECK(bg == 0);
    CHECK(bq == 0);
    CHECK(bgs == 0);
    CHECK(bm == 0);
}

// --- Repeat-rich, large-h0 adversarial parity (re-baseline regression) -------
//
// run_parity above uses a clean high-identity prefix copy: a single dominant
// diagonal, so even when the score exceeds 254 and re-baselining fires, the
// saturating-subtract only zeroes genuinely-dead cells and the 8-bit result
// stays exact. The 8-bit re-baseline truncation bug needed TANDEM-REPEAT
// structure -- strong shifted (off-diagonal) suboptimal alignments whose
// still-positive cells get wrongly zeroed (then misread as the h00==0
// local-restart sentinel). This generator reproduces that class: a short
// repeated unit + occasional indel + wide h0.
//
// Contract under test (what the bwamem.cpp routing gate relies on):
//   * Within the 8-bit envelope (re-baseline provably inert) getScores8 is
//     byte-identical to scalar -- INCLUDING scores in [128,254], which only the
//     unsigned [0,255] recurrence can represent (a revert to the old signed
//     cap-127 recurrence would fail here).
//   * getScores16 is byte-identical to scalar over the WHOLE set (it has no
//     re-baseline); the out-of-envelope repeat pairs route to it.

// Mirror of bwamem.cpp:bsw8_envelope_ok for the fixed scoring used here
// (a=1 => max_step=1, zdrop=100, w=100). Keep in sync with bwamem.cpp.
bool envelope_ok_8(int len1, int len2, int h0) {
    const int a = 1, zdrop = 100, w = 100, max_step = 1;
    int shorter = len1 < len2 ? len1 : len2;
    int max_score = h0 + shorter * a;
    return len1 < MAX_SEQ_LEN8 && len2 < MAX_SEQ_LEN8 && len1 >= len2 &&
           w <= 127 && zdrop + max_step <= 127 &&
           h0 <= zdrop + 1 && max_score < 255 - max_step;
}

struct RepeatStats {
    int bad8_in_env;   // 8-bit mismatches among envelope-admitted pairs (must be 0)
    int bad16_all;     // 16-bit mismatches over all pairs (must be 0)
    int n_in_env;      // # envelope-admitted pairs
    int n_in_env_hi;   // # admitted pairs scoring in [128,254] (exercises unsigned range)
    int n_out_env;     // # pairs the gate would route to 16-bit
};

// Generate n tandem-repeat target>=query pairs (wide h0), run getScores8 and
// getScores16 vs the scalar oracle, and tally the contract above.
RepeatStats run_repeat_parity(int n, int maxlen, int maxh0, unsigned long seed) {
    const int a = 1, b = 4, ambig = -1, o = 6, e = 1, zdrop = 100, end_bonus = 5, w = 100;
    const int STRIDE = maxlen + MAX_LINE_LEN;   // per-pair seq slot (cf. run_parity):
                                                // read length + prefetch look-ahead slack
    int8_t mat[25];
    build_mat(mat, a, b, ambig);
    BandedPairWiseSW bsw(o, e, o, e, zdrop, end_bonus, mat, a, b, 1);

    std::vector<uint8_t> ref((size_t)STRIDE * n, 0), qer((size_t)STRIDE * n, 0);
    std::vector<SeqPair> p8(padPairs(n)), p16(padPairs(n));   // padded (see padPairs)
    std::vector<Out> oracle(n);

    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<int> lenD(50, maxlen);
    std::uniform_int_distribution<int> hD(1, maxh0);
    std::uniform_int_distribution<int> unitD(3, 12);   // tandem-repeat unit length
    for (int c = 0; c < n; ++c) {
        int aa = lenD(rng), bb = lenD(rng);
        int len1 = aa > bb ? aa : bb, len2 = aa > bb ? bb : aa;   // target >= query
        int h0 = hD(rng);
        uint8_t *s1 = &ref[(size_t)c * STRIDE];
        uint8_t *s2 = &qer[(size_t)c * STRIDE];
        int u = unitD(rng);
        uint8_t ubuf[16];
        for (int i = 0; i < u; ++i) ubuf[i] = (uint8_t)(rng() % 4);
        for (int i = 0; i < len1; ++i) s1[i] = ubuf[i % u];        // tandem repeat target
        int ti = 0;
        for (int i = 0; i < len2; ++i) {
            uint8_t base = (ti < len1) ? s1[ti] : (uint8_t)(rng() % 4);
            if ((int)(rng() % 100) < 5) base = (uint8_t)(rng() % 4);   // 5% mismatch
            s2[i] = base;
            if ((rng() % 40) == 0) { if (rng() & 1) ti += 2; } else ++ti;  // occasional indel
        }
        for (SeqPair *P : {&p8[c], &p16[c]}) {
            P->id = c; P->len1 = len1; P->len2 = len2; P->h0 = h0;
            P->idr = (int)((size_t)c * STRIDE); P->idq = (int)((size_t)c * STRIDE);
            P->seqid = c; P->regid = c;
            P->score = P->tle = P->gtle = P->qle = P->gscore = P->max_off = -1;
        }
        Out &O = oracle[c];
        O.score = bsw.scalarBandedSWA(len2, s2, len1, s1, w, h0,
                                      &O.qle, &O.tle, &O.gtle, &O.gscore, &O.max_off);
    }

    bsw.getScores8(p8.data(), ref.data(), qer.data(), (int32_t)n, 1, w);
    bsw.getScores16(p16.data(), ref.data(), qer.data(), (int32_t)n, 1, w);

    RepeatStats st{0, 0, 0, 0, 0};
    for (int c = 0; c < n; ++c) {
        const Out &O = oracle[c];
        const SeqPair &q8 = p8[c], &q16 = p16[c];
        auto diff = [&](const SeqPair &p) {
            return O.score != p.score || O.tle != p.tle || O.gtle != p.gtle ||
                   O.qle != p.qle || O.gscore != p.gscore || O.max_off != p.max_off;
        };
        if (diff(q16)) st.bad16_all++;
        if (envelope_ok_8(q8.len1, q8.len2, q8.h0)) {
            st.n_in_env++;
            if (O.score >= 128 && O.score <= 254) st.n_in_env_hi++;
            if (diff(q8)) st.bad8_in_env++;
        } else {
            st.n_out_env++;
        }
    }
    return st;
}

} // namespace

TEST_CASE("bandedSWA getScores8 byte-identical to scalar (short + long reads)"
          * doctest::test_suite("unit/bandedswa")) {
    SUBCASE("short reads (maxlen 120)") { check_width(8, 3000, 120, 12345); }
    SUBCASE("long reads (maxlen 1000)") { check_width(8, 1500, 1000, 12345); }
}

TEST_CASE("bandedSWA getScores16 byte-identical to scalar (short + long reads)"
          * doctest::test_suite("unit/bandedswa")) {
    SUBCASE("short reads (maxlen 120)") { check_width(16, 3000, 120, 999); }
    SUBCASE("long reads (maxlen 1000)") { check_width(16, 1500, 1000, 999); }
}

TEST_CASE("bandedSWA getScores16 byte-identical to scalar past MAX_SEQ_LEN8 + high zdrop"
          * doctest::test_suite("unit/bandedswa")) {
    // The 8-bit path caps at MAX_SEQ_LEN8 (1088); reads beyond that use the 16-bit
    // kernel exclusively, a regime that previously had no byte-identity coverage.
    SUBCASE("16-bit long reads (maxlen 2000)") { check_width(16, 600, 2000, 4242); }
    SUBCASE("16-bit long reads (maxlen 4000)") { check_width(16, 300, 4000, 4243); }
    // zdrop beyond the 8-bit gate's <=126 cap: exercises the 16-bit path's
    // handling of a large z-drop horizon (getScores16 has no re-baseline, so it
    // takes zdrop directly; this also covers what the 8-bit gate routes here).
    SUBCASE("16-bit high zdrop (maxlen 1000, zdrop 150)") { check_width(16, 1500, 1000, 4244, 150); }
}

TEST_CASE("bandedSWA 8-bit byte-identical to scalar on repeat-rich, large-h0 reads"
          " within the routing envelope (re-baseline regression)"
          * doctest::test_suite("unit/bandedswa")) {
    RepeatStats st = run_repeat_parity(10000, 450, 300, 1);
    MESSAGE("repeat-rich parity: in_env=" << st.n_in_env
            << " (hi-score[128,254]=" << st.n_in_env_hi << ")"
            << " out_env=" << st.n_out_env
            << " bad8_in_env=" << st.bad8_in_env
            << " bad16_all=" << st.bad16_all);
    // 8-bit is byte-exact for every envelope-admitted pair (incl. [128,254]).
    CHECK(st.bad8_in_env == 0);
    // 16-bit is byte-exact for the whole set -- the path out-of-envelope pairs
    // route to.
    CHECK(st.bad16_all == 0);
    // Non-vacuity: the case must actually straddle the envelope boundary and
    // exercise the unsigned [128,254] range, else a generation change could
    // silently make the parity checks trivial.
    CHECK(st.n_in_env > 0);
    CHECK(st.n_in_env_hi > 0);
    CHECK(st.n_out_env > 0);
}

// --- Asymmetric (--meth OT/OB) query-end gscore/gtle capture --------------------
//
// Regression for the 8-bit "exit due to zero score by a row" break dropping the
// current row's query-end (gscore/gtle) capture. Under --meth's asymmetric
// substitution matrix (OT frees ref-C/read-T), a real captured extension pair
// reaches the query end at a late target row where the whole band has decayed to
// 0 (gscore==0). scalar records that row's gscore before its m==0 break; the 8-bit
// kernel's zero-row break used to exit BEFORE folding the capture into the wide
// gbest_abs/ierow, so getScores8 reported gtle from an early row (1) instead of the
// scalar row (67). Benign for symmetric scoring (gscore==0 gtle is unused), but
// consumed under --meth, where it produced spurious soft-clips / placement drift.
//
// The captured pair (from a 1M-pair meth slice) scores with the bisulfite params
// a=1, b=2, gap open/extend 6/1, end_bonus (pen_clip) 10, band 100, seed h0 27.
TEST_CASE("bandedSWA 8-bit query-end gscore/gtle byte-identical to scalar under an"
          " asymmetric --meth (OT) matrix (zero-score-row break capture)"
          * doctest::test_suite("unit/bandedswa")) {
    const int a = 1, b = 2, o = 6, e = 1, end_bonus = 10, w = 100, zdrop = 100, h0 = 27;
    // OT matrix: symmetric a/-b plus freed ref-C/read-T (and the collapsed mirror
    // ref-T/read-C), ambiguous = -1 — exactly mem_opt_fill_meth_mat's OT output.
    int8_t mat[25];
    build_mat(mat, a, b, -1);
    mat[1 * 5 + 3] = (int8_t)a;   // ref C / read T -> match
    mat[3 * 5 + 1] = (int8_t)a;   // ref T / read C -> match (collapsed)

    const uint8_t ref[] = {2,2,0,0,0,0,3,0,0,3,3,0,1,3,3,3,3,0,1,0,0,2,3,0,0,0,3,0,1,0,3,1,2,0,1,2,2,0,3,0,3,0,0,0,2,2,3,1,1,1,3,3,3,2,3,0,3,1,1,3,2,3,3,1,2,2,3,1,0,3,3,0,3,0,0,3,1,0,3,3,0,2,0,0,2,3,0,0,3,0,1,2,3,3,3,0,1,0,0,2,3,0,3};
    const uint8_t qer[] = {3,3,3,3,0,0,0,3,2,3,3,3,0,3,0,3,3,2,0,0,3,3,2,0,3,0,3,0,3,2,2,3,3,3,3,3,0,3,3,3,3,3,3,3,3,3,3,0,3,3,3,0,0,2};
    const int len1 = (int)(sizeof(ref) / sizeof(ref[0]));
    const int len2 = (int)(sizeof(qer) / sizeof(qer[0]));

    BandedPairWiseSW bsw(o, e, o, e, zdrop, end_bonus, mat, a, b, 1);

    const int STRIDE = len1 + MAX_LINE_LEN;
    std::vector<uint8_t> seqRef(STRIDE, 0), seqQer(STRIDE, 0);
    for (int i = 0; i < len1; ++i) seqRef[i] = ref[i];
    for (int i = 0; i < len2; ++i) seqQer[i] = qer[i];

    Out oracle;
    oracle.score = bsw.scalarBandedSWA(len2, seqQer.data(), len1, seqRef.data(), w, h0,
                                       &oracle.qle, &oracle.tle, &oracle.gtle,
                                       &oracle.gscore, &oracle.max_off);

    std::vector<SeqPair> pairs(padPairs(1));
    pairs[0].id = 0; pairs[0].len1 = len1; pairs[0].len2 = len2; pairs[0].h0 = h0;
    pairs[0].idr = 0; pairs[0].idq = 0; pairs[0].seqid = 0; pairs[0].regid = 0;
    pairs[0].score = pairs[0].tle = pairs[0].gtle = pairs[0].qle =
        pairs[0].gscore = pairs[0].max_off = -1;
    bsw.getScores8(pairs.data(), seqRef.data(), seqQer.data(), 1, 1, w);

    MESSAGE("meth-OT query-end: scalar gtle=" << oracle.gtle << " gscore=" << oracle.gscore
            << " | get8 gtle=" << pairs[0].gtle << " gscore=" << pairs[0].gscore);
    CHECK(pairs[0].score   == oracle.score);
    CHECK(pairs[0].qle     == oracle.qle);
    CHECK(pairs[0].tle     == oracle.tle);
    CHECK(pairs[0].gtle    == oracle.gtle);    // was 1 vs scalar 67 before the fix
    CHECK(pairs[0].gscore  == oracle.gscore);
    CHECK(pairs[0].max_off == oracle.max_off);
    // Non-vacuity: this input genuinely exercises the query-end tail (gtle deep in
    // the target, well past the local end) so the check is not trivially satisfied.
    CHECK(oracle.gtle > oracle.tle + 10);
}

// Regression test for the --meth per-hypothesis sub-slice overshoot bug.
//
// The batched kernels round numPairs up to their SIMD width and initialize the
// padding lanes pairArray[numPairs .. roundUp) in place (len=0 dummy pairs).
// A caller may legitimately invoke the kernel on a SUB-SLICE of a larger array
// -- the --meth dispatch (bsw_run_tier) runs three contiguous slices
// (OT | OB | SYM) of one pair array back-to-back. If the kernel scribbles
// padding past numPairs it zeroes the START of the next slice; that slice's
// kernel call then scores those real pairs as empty (score==h0, gscore==-1),
// producing a spurious soft-clip / confident mismap on --meth reads. So the
// kernels must not modify pairArray outside [0, numPairs).
TEST_CASE("bandedSWA getScores16/getScores8 leave pairArray[numPairs..] untouched"
          " (sub-slice overshoot safety)"
          * doctest::test_suite("unit/bandedswa")) {
    const int a = 1, b = 2, o = 6, e = 1, end_bonus = 10, w = 100, zdrop = 100;
    int8_t mat[25];
    build_mat(mat, a, b, -1);
    BandedPairWiseSW bsw(o, e, o, e, zdrop, end_bonus, mat, a, b, 1);
    // This test emulates a sub-slice caller (the --meth OT/OB kernels), for which
    // the overshoot guard is enabled; it is off by default for whole-array callers.
    bsw.set_guard_overshoot(true);

    // n deliberately NOT a multiple of any SIMD width so roundUp(n) > n and the
    // kernel has padding lanes to initialize past numPairs.
    const int n = 5;
    const int GUARD = SIMD_WIDTH8;          // widest tier's overshoot
    const int STRIDE = 64;
    std::vector<uint8_t> ref((size_t)(n + GUARD) * STRIDE + MAX_LINE_LEN, 0);
    std::vector<uint8_t> qer((size_t)(n + GUARD) * STRIDE + MAX_LINE_LEN, 0);
    std::vector<SeqPair> pairs((size_t)(n + GUARD) + MAX_LINE_LEN, SeqPair());

    // n real all-match pairs (ref A-run vs read A-run).
    for (int c = 0; c < n; ++c) {
        uint8_t *rs = &ref[(size_t)c * STRIDE];
        uint8_t *qs = &qer[(size_t)c * STRIDE];
        for (int i = 0; i < 20; ++i) rs[i] = 0;
        for (int i = 0; i < 15; ++i) qs[i] = 0;
        pairs[c].id = c; pairs[c].len1 = 20; pairs[c].len2 = 15; pairs[c].h0 = 30;
        pairs[c].idr = (int64_t)c * STRIDE; pairs[c].idq = (int64_t)c * STRIDE;
        pairs[c].seqid = c; pairs[c].regid = 0;
        pairs[c].score = pairs[c].tle = pairs[c].gtle = pairs[c].qle =
            pairs[c].gscore = pairs[c].max_off = -1;
    }
    // GUARD sentinel pairs occupying pairArray[n .. n+GUARD): stand-ins for the
    // next slice's real pairs. Distinctive values so any kernel write is caught.
    auto set_sentinel = [](SeqPair &p, int c) {
        p.id = 1000 + c; p.len1 = 37; p.len2 = 21; p.h0 = 129;
        p.idr = 700000 + c; p.idq = 800000 + c;
        p.seqid = 500 + c; p.regid = 3;
        p.score = 111; p.qle = 112; p.tle = 113;
        p.gscore = 114; p.gtle = 115; p.max_off = 116;
    };
    std::vector<SeqPair> expected(GUARD);
    for (int c = 0; c < GUARD; ++c) { set_sentinel(pairs[n + c], c); expected[c] = pairs[n + c]; }

    auto sentinels_intact = [&](const char *which) {
        for (int c = 0; c < GUARD; ++c) {
            const SeqPair &g = pairs[n + c], &x = expected[c];
            INFO("call=" << which << " sentinel=" << c);
            CHECK(g.id == x.id); CHECK(g.len1 == x.len1); CHECK(g.len2 == x.len2);
            CHECK(g.h0 == x.h0); CHECK(g.idr == x.idr); CHECK(g.idq == x.idq);
            CHECK(g.seqid == x.seqid); CHECK(g.regid == x.regid);
            CHECK(g.score == x.score); CHECK(g.gscore == x.gscore);
            CHECK(g.qle == x.qle); CHECK(g.tle == x.tle);
            CHECK(g.gtle == x.gtle); CHECK(g.max_off == x.max_off);
        }
    };

    bsw.getScores16(pairs.data(), ref.data(), qer.data(), n, 1, w);
    sentinels_intact("getScores16");

    for (int c = 0; c < GUARD; ++c) set_sentinel(pairs[n + c], c);  // reset
    for (int c = 0; c < n; ++c) {
        pairs[c].score = pairs[c].tle = pairs[c].gtle = pairs[c].qle =
            pairs[c].gscore = pairs[c].max_off = -1;
    }
    bsw.getScores8(pairs.data(), ref.data(), qer.data(), n, 1, w);
    sentinels_intact("getScores8");
}

#endif // HAVE_BSW_VECTOR_8_16
