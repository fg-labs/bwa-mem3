// test/unit/test_meth_asym_directed.cpp
//
// D3 B5 — directed γ (asymmetric bisulfite) scoring tests.
//
// test_bandedswa_asym.cpp already proves the SIMD kernel HONORS an arbitrary
// asymmetric matrix (getScores8/16 == scalar oracle). These tests assert the
// directed SEMANTICS the D3 design depends on, exercising the production 8-bit
// kernel (getScores8) directly:
//
//   1. The OT matrix CHANGES the selected score by exactly the freed penalty —
//      a C->T conversion costs nothing under OT but (a+b) under symmetric.
//   2. Hypothesis selection FLIPS by conversion direction — a C->T read scores
//      higher under OT than OB, a G->A read higher under OB than OT.
//   3. #136-style segdup substrate: a competitor differing from the true window
//      only at a C/T (freed) position stays TIED with it (paralogs are not
//      falsely separated post-conversion), while a competitor differing at a
//      non-freed position separates — AND the reverse cell mat[T][C] (read-C at
//      ref-T) is NOT freed (negative control against freeing too much).
//
// Assertions are on score DELTAS between matrices/windows on the same pair, which
// cancel any h0 / end-bonus baseline and isolate the matrix effect.
//
// A,C,G,T = 0,1,2,3; target-major mat[ref*5 + read].

#include <cstdint>
#include <random>
#include <vector>

#include "doctest/doctest.h"
#include "bandedSWA.h"

#if HAVE_BSW_VECTOR_8_16

namespace {

constexpr int A_ = 0, C_ = 1, G_ = 2, T_ = 3;
constexpr int MATCH = 1, MISMATCH = 4, AMBIG = -1;
constexpr int FREED = MATCH + MISMATCH;  // score gained when a -b cell becomes +a

void build_sym_mat(int8_t mat[25]) {
    int k = 0;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) mat[k++] = (i == j) ? (int8_t)MATCH : (int8_t)(-MISMATCH);
        mat[k++] = (int8_t)AMBIG;
    }
    for (int j = 0; j < 5; ++j) mat[k++] = (int8_t)AMBIG;
}
// OT: free ref-C x read-T (an unmethylated C read as T).
void build_ot_mat(int8_t mat[25]) { build_sym_mat(mat); mat[C_ * 5 + T_] = (int8_t)MATCH; }
// OB: free ref-G x read-A (bottom-strand C->T appears as G->A on the top strand).
void build_ob_mat(int8_t mat[25]) { build_sym_mat(mat); mat[G_ * 5 + A_] = (int8_t)MATCH; }

// Score ONE (ref, read) pair through the production 8-bit kernel under `mat`.
int score_one(const int8_t mat[25], const std::vector<uint8_t> &ref,
              const std::vector<uint8_t> &read) {
    const int o = 6, e = 1, end_bonus = 5, w = 100, zdrop = 100;
    BandedPairWiseSW bsw(o, e, o, e, zdrop, end_bonus, const_cast<int8_t*>(mat),
                         MATCH, MISMATCH, 1);
    const int len1 = (int)ref.size(), len2 = (int)read.size();
    const int stride = (len1 > len2 ? len1 : len2) + 256;
    std::vector<uint8_t> refbuf(stride, 0), qerbuf(stride, 0);
    for (int i = 0; i < len1; ++i) refbuf[i] = ref[i];
    for (int i = 0; i < len2; ++i) qerbuf[i] = read[i];
    std::vector<SeqPair> pairs(64);   // SIMD-width round-up (getScores8 writes padding lanes)
    SeqPair &sp = pairs[0];
    // h0 is the seed baseline carried into the DP; the 8-bit kernel's
    // diagonal-offset encoding needs a positive h0 (h0=0 degenerates). The score
    // deltas these tests assert on are h0-independent, so any positive value works.
    sp.id = 0; sp.len1 = len1; sp.len2 = len2; sp.h0 = 40;
    sp.idr = 0; sp.idq = 0; sp.seqid = 0; sp.regid = 0;
    sp.score = sp.tle = sp.gtle = sp.qle = sp.gscore = sp.max_off = -1;
    bsw.getScores8(pairs.data(), refbuf.data(), qerbuf.data(), 1, 1, w);
    return sp.score;
}

// A deterministic ref window with no homopolymer artifacts, then forced bases.
std::vector<uint8_t> make_ref(int len, unsigned long seed) {
    std::vector<uint8_t> v(len);
    std::mt19937_64 rng(seed);
    for (int i = 0; i < len; ++i) v[i] = (uint8_t)(rng() % 4);
    return v;
}

} // namespace

TEST_CASE("γ directed: OT matrix changes the selected score by exactly the freed penalty"
          * doctest::test_suite("unit/meth-asym-directed")) {
    const int L = 60, conv[3] = {15, 30, 45};
    auto ref = make_ref(L, 11);
    for (int p : conv) ref[p] = C_;          // force ref-C at the conversion sites
    auto read = ref;                          // identical...
    for (int p : conv) read[p] = T_;          // ...except 3 C->T conversions

    int8_t sym[25], ot[25]; build_sym_mat(sym); build_ot_mat(ot);
    int s_sym = score_one(sym, ref, read);
    int s_ot  = score_one(ot,  ref, read);

    // OT frees all 3 conversions; symmetric pays (a+b) for each.
    CHECK(s_ot - s_sym == 3 * FREED);
    CHECK(s_ot > s_sym);
}

TEST_CASE("γ directed: hypothesis selection flips with conversion direction"
          * doctest::test_suite("unit/meth-asym-directed")) {
    const int L = 60, sites[3] = {15, 30, 45};
    int8_t ot[25], ob[25]; build_ot_mat(ot); build_ob_mat(ob);

    SUBCASE("C->T read favors OT over OB") {
        auto ref = make_ref(L, 21);
        for (int p : sites) ref[p] = C_;
        auto read = ref; for (int p : sites) read[p] = T_;
        CHECK(score_one(ot, ref, read) - score_one(ob, ref, read) == 3 * FREED);
    }
    SUBCASE("G->A read favors OB over OT") {
        auto ref = make_ref(L, 22);
        for (int p : sites) ref[p] = G_;
        auto read = ref; for (int p : sites) read[p] = A_;
        CHECK(score_one(ob, ref, read) - score_one(ot, ref, read) == 3 * FREED);
    }
}

TEST_CASE("γ directed: #136 segdup substrate — C/T-only competitor ties, others separate"
          * doctest::test_suite("unit/meth-asym-directed")) {
    const int L = 60, conv[3] = {15, 30, 45};
    auto trueref = make_ref(L, 31);
    for (int p : conv) trueref[p] = C_;        // true locus has C at the conversion sites
    // Force a known non-C site for the negative control / over-freeing checks.
    const int site_T = 40; trueref[site_T] = T_;
    const int site_A = 50; trueref[site_A] = A_;

    auto read = trueref;
    for (int p : conv) read[p] = T_;           // the bisulfite-converted read (C->T)

    int8_t ot[25]; build_ot_mat(ot);
    int s_true = score_one(ot, trueref, read);

    SUBCASE("paralog differing ONLY at the C/T sites stays tied (not falsely separated)") {
        auto comp = trueref; for (int p : conv) comp[p] = T_;   // paralog has T where true has C
        // read has T at conv sites: matches comp exactly (+a), matches true via the
        // freed ref-C x read-T cell (+a). Equal score => genuine MAPQ ambiguity kept.
        CHECK(score_one(ot, comp, read) == s_true);
    }
    SUBCASE("competitor differing at a non-freed base separates") {
        auto comp = trueref; comp[site_A] = G_;   // read has A at site_A; comp has G => -b
        CHECK(s_true - score_one(ot, comp, read) == FREED);
    }
    SUBCASE("over-freeing guard: reverse cell mat[T][C] (read-C at ref-T) is NOT freed") {
        auto read2 = trueref;            // start from a perfect-match read
        read2[site_T] = C_;              // read-C against ref-T: must be penalized under OT
        CHECK(score_one(ot, trueref, trueref) - score_one(ot, trueref, read2) == FREED);
    }
}

#endif // HAVE_BSW_VECTOR_8_16
