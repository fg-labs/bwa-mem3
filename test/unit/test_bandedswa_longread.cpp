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
int run_parity(int width, int n, int maxlen, unsigned long seed,
               int &bs, int &bt, int &bg, int &bq, int &bgs, int &bm) {
    const int a = 1, b = 4, ambig = -1, o = 6, e = 1, zdrop = 100, end_bonus = 5, w = 100;
    const int STRIDE = 1280;   // > MAX_SEQ_LEN8 (1088); generous per-pair slot
    int8_t mat[25];
    build_mat(mat, a, b, ambig);

    BandedPairWiseSW bsw(o, e, o, e, zdrop, end_bonus, mat, a, b, 1);

    std::vector<uint8_t> seqBufRef((size_t)STRIDE * n, 0);
    std::vector<uint8_t> seqBufQer((size_t)STRIDE * n, 0);
    std::vector<SeqPair> pairs(n);
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

void check_width(int width, int n, int maxlen, unsigned long seed) {
    int bs, bt, bg, bq, bgs, bm;
    int bad = run_parity(width, n, maxlen, seed, bs, bt, bg, bq, bgs, bm);
    MESSAGE("bandedSWA getScores" << width << " vs scalar: maxlen=" << maxlen
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

#endif // HAVE_BSW_VECTOR_8_16
