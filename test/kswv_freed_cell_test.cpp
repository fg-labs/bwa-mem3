// Issue 173 / Task 2: mat-aware make_kswv freed-cell detection.
//
// The kswv factory learns a per-strand asymmetric (freed-cell)
// substitution for bisulfite (--meth). A symmetric matrix (or nullptr)
// constructs the existing kernel unchanged; a matrix that frees exactly one
// ordered off-diagonal cell (OT: ref-C/read-T, OB: ref-G/read-A) to ANY value
// records that freed cell and its value and is handled by the kernel — a match
// (+a) under GENOMIC, 0 under NEUTRAL; an exact mirrored pair (collapsed
// --meth) is likewise handled; any other asymmetric matrix (a non-mirror
// multi-cell free, changed diagonal) routes to the scalar fallback
// (needsScalar() == true).
//
// Two layers are tested: (1) ctor-side detection (needsScalar) on every tier,
// and (2) the kernel actually APPLYING the freed cell via getScores8. Layer (2)
// runs only on tiers with the freed-cell kernel override (NEON/AVX2/AVX-512BW);
// on SSE-only tiers make_kswv reports needsScalar()==true and layer (1) verifies
// that fallback contract instead (see tier_supports_freed_cell below).

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"
#include "ksw.h"
#include "kswv.h"
#include "simd_dispatch.h"
#include <cstring>
#include <vector>
#include <cstdint>

static void fill_sym(int8_t m[25], int a, int b) {
    for (int i=0;i<5;i++) for (int j=0;j<5;j++)
        m[i*5+j] = (i==4||j==4) ? 0 : (i==j ? (int8_t)a : (int8_t)-b);
}

// True when the host's native kswv tier implements the freed-cell kernel
// override (NEON / AVX2 / AVX-512BW). On SSE41/SSE42/AVX there is no override:
// make_kswv reports needsScalar() == true for a freed-cell matrix (so the caller
// routes to scalar ksw_align2), and getScores8/16 are unreachable exit() stubs.
// The freed-cell DETECTION cases below therefore assert needsScalar() == !FREED
// for single-cell (any freed value) and mirror matrices, and the score-delta
// cases (which call getScores8) run only when FREED — on a scalar-only tier the
// fallback contract is what the detection cases verify, and there is no kernel
// score to compare.
//
// This gate is derived from the INDEPENDENT SIMD dispatch tier
// (bwamem3_simd_tier, which is the same tier make_kswv dispatches on and which
// honors BWAMEM3_FORCE_TIER), NOT from needsScalar() on the matrix under test.
// A needsScalar()-based probe would be self-referential: a constructor
// regression that wrongly returned needsScalar()==true on a freed-cell tier
// would also flip this gate to false, so the line-79-style assertions
// (needsScalar() == !tier_supports_freed_cell()) would still pass and the
// score-delta cases would silently become no-ops. Keying off the tier oracle
// instead leaves needsScalar() strictly as the assertion target, so such a
// regression FAILS the suite instead of hiding in it.
static bool tier_supports_freed_cell() {
    bwamem3_simd_init();
    switch (bwamem3_simd_tier()) {
        case BWAMEM3_TIER_NEON:
        case BWAMEM3_TIER_AVX2:
        case BWAMEM3_TIER_AVX512BW:
            return true;
        default:
            return false;
    }
}

// Run getScores8 over a single ref/read pair through the batched NEON 8-bit
// kernel and return aln[0].score. `mat25 == nullptr` constructs the symmetric
// kernel (no freed cell); a meth matrix sets has_freed and blends the freed
// cell to fr_val. Ref/read are base-encoded (A=0,C=1,G=2,T=3). The pair is laid out
// so the optimal local alignment is the full diagonal — every column is a
// match except the deliberately planted converted base.
static int score8(const std::vector<uint8_t>& ref,
                  const std::vector<uint8_t>& read,
                  int a, int b, const int8_t* mat25) {
    const int32_t maxRefLen = 256, maxQerLen = 256;
    SeqPair sp = {};
    sp.idr = 0; sp.idq = 0; sp.id = 0;
    sp.len1 = (int)ref.size();
    sp.len2 = (int)read.size();
    sp.h0 = 0;                 // no KSW_XSUBO/XSTOP: plain best-score scan
    sp.seqid = 0; sp.regid = 0;

    std::vector<SeqPair> pairs(SIMD_WIDTH8 + 1);
    pairs[0] = sp;

    std::vector<kswr_t> aln(SIMD_WIDTH8 + 1, g_defr);
    // Mutable copies for the (non-const) buffer pointers getScores8 expects.
    std::vector<uint8_t> rbuf = ref, qbuf = read;

    kswv k(6, 1, 6, 1, (int8_t)a, (int8_t)-b, 1, maxRefLen, maxQerLen, mat25);
    // phase 0 fills aln[regid].score with the de-biased best SW score.
    k.getScores8(pairs.data(), rbuf.data(), qbuf.data(), aln.data(), 1, 1, 0);
    return aln[0].score;
}
TEST_CASE("kswv detects OT freed cell (ref-C x read-T)") {
    int8_t m[25]; fill_sym(m,1,4); m[1*5+3] = 1;          // free ref-C x read-T to match
    auto k = make_kswv(6,1,6,1, 1,-4, 1, 256,256, m);
    // rank-1 is handled on a freed-cell tier; scalar fallback on SSE-only tiers.
    CHECK(k->needsScalar() == !tier_supports_freed_cell());
}
TEST_CASE("kswv detects OB freed cell (ref-G x read-A)") {
    int8_t m[25]; fill_sym(m,1,4); m[2*5+0] = 1;
    auto k = make_kswv(6,1,6,1, 1,-4, 1, 256,256, m);
    CHECK(k->needsScalar() == !tier_supports_freed_cell());
}
TEST_CASE("kswv symmetric matrix is not a freed-cell case") {
    int8_t m[25]; fill_sym(m,1,4);
    auto k = make_kswv(6,1,6,1, 1,-4, 1, 256,256, m);
    CHECK(k->needsScalar() == false);                      // symmetric ⇒ normal kernel
}
TEST_CASE("kswv non-mirror multi-cell asymmetric falls back to scalar") {
    int8_t m[25]; fill_sym(m,1,4); m[1*5+3]=1; m[2*5+0]=1; // two freed, not mirrors
    auto k = make_kswv(6,1,6,1, 1,-4, 1, 256,256, m);
    CHECK(k->needsScalar() == true);
}

// Task 4: the NEON kernel must actually APPLY the freed cell. Build a clean
// matching alignment with exactly one ref-C×read-T (OT) cell, score it through
// getScores8 under the symmetric vs OT matrix, and assert the OT score exceeds
// the symmetric score by exactly (a + b): the converted cell flips from a
// mismatch (-b) to a match (+a). OB is the ref-G×read-A analog.
TEST_CASE("kswv NEON 8-bit applies the OT freed cell (ref-C x read-T)") {
    if (!tier_supports_freed_cell()) return;   // scalar-only tier: no kernel score
    const int a = 1, b = 4;
    // 20-base identical diagonal, then plant ref=C(1)/read=T(3) at one column.
    std::vector<uint8_t> ref  = {0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3};
    std::vector<uint8_t> read = ref;
    const int p = 10;                 // converted column
    ref[p]  = 1;                      // ref C
    read[p] = 3;                      // read T (bisulfite C->T on read)

    int8_t sym[25]; fill_sym(sym, a, b);
    int8_t ot[25];  fill_sym(ot, a, b); ot[1*5+3] = (int8_t)a;  // free C×T -> match

    int s_sym = score8(ref, read, a, b, sym);
    int s_ot  = score8(ref, read, a, b, ot);
    // exactly one freed cell on the optimal path: delta == a + b
    CHECK(s_ot - s_sym == a + b);
}

TEST_CASE("kswv NEON 8-bit applies the OB freed cell (ref-G x read-A)") {
    if (!tier_supports_freed_cell()) return;   // scalar-only tier: no kernel score
    const int a = 1, b = 4;
    std::vector<uint8_t> ref  = {0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3};
    std::vector<uint8_t> read = ref;
    const int p = 10;
    ref[p]  = 2;                      // ref G
    read[p] = 0;                      // read A (bisulfite G->A on read, OB strand)

    int8_t sym[25]; fill_sym(sym, a, b);
    int8_t ob[25];  fill_sym(ob, a, b); ob[2*5+0] = (int8_t)a;  // free G×A -> match

    int s_sym = score8(ref, read, a, b, sym);
    int s_ob  = score8(ref, read, a, b, ob);
    CHECK(s_ob - s_sym == a + b);
}

// Two freed cells on the optimal path must each contribute (a + b). This also
// guards against an override that only fires once per row/strip.
TEST_CASE("kswv NEON 8-bit applies multiple OT freed cells additively") {
    if (!tier_supports_freed_cell()) return;   // scalar-only tier: no kernel score
    const int a = 1, b = 4;
    std::vector<uint8_t> ref  = {0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3};
    std::vector<uint8_t> read = ref;
    ref[6]  = 1; read[6]  = 3;        // OT cell #1 (different rows)
    ref[14] = 1; read[14] = 3;        // OT cell #2

    int8_t sym[25]; fill_sym(sym, a, b);
    int8_t ot[25];  fill_sym(ot, a, b); ot[1*5+3] = (int8_t)a;

    int s_sym = score8(ref, read, a, b, sym);
    int s_ot  = score8(ref, read, a, b, ot);
    CHECK(s_ot - s_sym == 2 * (a + b));
}

// ---------------------------------------------------------------------------
// Issue 173 follow-up: COLLAPSED bisulfite scoring (the --meth default) frees
// the conversion cell AND its mirror, so C/T (and G/A) are interchangeable —
// two freed cells forming a symmetric (i,j)/(j,i) pair. The batched kernel must
// handle this (not route to scalar / not assert), freeing BOTH cells. A
// NON-mirror two-cell matrix (e.g. OT+OB combined) is still scalar.
// ---------------------------------------------------------------------------
TEST_CASE("kswv detects collapsed OT mirror pair (C<->T) — handled, not scalar") {
    int8_t m[25]; fill_sym(m,1,4); m[1*5+3]=1; m[3*5+1]=1;  // free (C,T) and (T,C)
    auto k = make_kswv(6,1,6,1, 1,-4, 1, 256,256, m);
    // mirror pair handled on a freed-cell tier; scalar fallback on SSE-only.
    CHECK(k->needsScalar() == !tier_supports_freed_cell());
}
TEST_CASE("kswv detects collapsed OB mirror pair (G<->A) — handled, not scalar") {
    int8_t m[25]; fill_sym(m,1,4); m[2*5+0]=1; m[0*5+2]=1;  // free (G,A) and (A,G)
    auto k = make_kswv(6,1,6,1, 1,-4, 1, 256,256, m);
    CHECK(k->needsScalar() == !tier_supports_freed_cell());
}
// The decisive behavior: collapsed frees the MIRROR cell that genomic leaves as
// a mismatch. Plant ONE mirror cell on the optimal diagonal (clean arithmetic,
// same as the single-cell tests above). Under genomic OT the mirror (T,C) stays
// a mismatch (delta 0); under collapsed OT it is freed (delta a+b).
TEST_CASE("kswv NEON 8-bit: collapsed OT frees the mirror cell (ref-T x read-C)") {
    if (!tier_supports_freed_cell()) return;   // scalar-only tier: no kernel score
    const int a = 1, b = 4;
    std::vector<uint8_t> ref  = {0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3};
    std::vector<uint8_t> read = ref;
    const int p = 10;
    ref[p]  = 3;                      // ref T
    read[p] = 1;                      // read C  → the (T,C) mirror of the OT cell

    int8_t sym[25]; fill_sym(sym, a, b);
    int8_t gen[25]; fill_sym(gen, a, b); gen[1*5+3]=(int8_t)a;                       // genomic OT: (C,T) only
    int8_t col[25]; fill_sym(col, a, b); col[1*5+3]=(int8_t)a; col[3*5+1]=(int8_t)a; // collapsed OT: + mirror (T,C)

    int s_sym = score8(ref, read, a, b, sym);
    int s_gen = score8(ref, read, a, b, gen);
    int s_col = score8(ref, read, a, b, col);
    CHECK(s_gen - s_sym == 0);          // genomic does NOT free the mirror
    CHECK(s_col - s_sym == a + b);      // collapsed frees the mirror cell
}
TEST_CASE("kswv NEON 8-bit: collapsed OB frees the mirror cell (ref-A x read-G)") {
    if (!tier_supports_freed_cell()) return;   // scalar-only tier: no kernel score
    const int a = 1, b = 4;
    std::vector<uint8_t> ref  = {0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3};
    std::vector<uint8_t> read = ref;
    const int p = 10;
    ref[p]  = 0;                      // ref A
    read[p] = 2;                      // read G  → the (A,G) mirror of the OB cell

    int8_t sym[25]; fill_sym(sym, a, b);
    int8_t gen[25]; fill_sym(gen, a, b); gen[2*5+0]=(int8_t)a;                       // genomic OB: (G,A) only
    int8_t col[25]; fill_sym(col, a, b); col[2*5+0]=(int8_t)a; col[0*5+2]=(int8_t)a; // collapsed OB: + mirror (A,G)

    int s_sym = score8(ref, read, a, b, sym);
    int s_gen = score8(ref, read, a, b, gen);
    int s_col = score8(ref, read, a, b, col);
    CHECK(s_gen - s_sym == 0);
    CHECK(s_col - s_sym == a + b);
}
// Additive guard: with BOTH the conversion cell and its mirror on the optimal
// path, collapsed must score exactly (a+b) above genomic — i.e. it frees the
// one extra (mirror) cell, and the kernel's second blend fires independently of
// the first. (Robust to the local-alignment baseline: compares col vs gen.)
TEST_CASE("kswv NEON 8-bit: collapsed OT beats genomic by exactly the mirror cell") {
    if (!tier_supports_freed_cell()) return;   // scalar-only tier: no kernel score
    const int a = 1, b = 4;
    std::vector<uint8_t> ref  = {0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3};
    std::vector<uint8_t> read = ref;
    ref[6]  = 1; read[6]  = 3;        // (C,T) conversion cell — freed by both
    ref[14] = 3; read[14] = 1;        // (T,C) mirror cell — freed only by collapsed

    int8_t gen[25]; fill_sym(gen, a, b); gen[1*5+3]=(int8_t)a;
    int8_t col[25]; fill_sym(col, a, b); col[1*5+3]=(int8_t)a; col[3*5+1]=(int8_t)a;

    int s_gen = score8(ref, read, a, b, gen);
    int s_col = score8(ref, read, a, b, col);
    CHECK(s_col - s_gen == a + b);
}
TEST_CASE("kswv non-mirror two-cell matrix (OT+OB) still routes to scalar") {
    int8_t m[25]; fill_sym(m,1,4); m[1*5+3]=1; m[2*5+0]=1;  // (C,T)+(G,A): not a mirror pair
    auto k = make_kswv(6,1,6,1, 1,-4, 1, 256,256, m);
    CHECK(k->needsScalar() == true);
}

// ---------------------------------------------------------------------------
// NEUTRAL scoring (the --meth=taps default): the conversion cell is freed to 0
// — tolerated but NOT rewarded — while the mirror stays a real mismatch. The
// kernel blends the freed cell to the matrix value (fr_val), not to a match, so
// the freed value is what these cases pin down. A freed value other than the
// match score used to route to scalar; it must now run batched on a freed-cell
// tier, exactly like GENOMIC.
// ---------------------------------------------------------------------------
TEST_CASE("kswv detects neutral OT freed cell (ref-C x read-T scored 0)") {
    int8_t m[25]; fill_sym(m,1,4); m[1*5+3] = 0;           // free (C,T) to 0, not to +a
    auto k = make_kswv(6,1,6,1, 1,-4, 1, 256,256, m);
    CHECK(k->needsScalar() == !tier_supports_freed_cell());
}
TEST_CASE("kswv detects neutral OB freed cell (ref-G x read-A scored 0)") {
    int8_t m[25]; fill_sym(m,1,4); m[2*5+0] = 0;
    auto k = make_kswv(6,1,6,1, 1,-4, 1, 256,256, m);
    CHECK(k->needsScalar() == !tier_supports_freed_cell());
}
// The 8-bit kernel blends fr_val into a biased u8 domain that only represents
// [w_mismatch, w_match]. A single freed cell outside that range would wrap the
// biased byte, so the ctor must route it to the scalar fallback instead.
TEST_CASE("kswv freed value outside [w_mismatch, w_match] routes to scalar") {
    int8_t lo[25]; fill_sym(lo,1,4); lo[1*5+3] = -5;   // below the mismatch score
    CHECK(make_kswv(6,1,6,1, 1,-4, 1, 256,256, lo)->needsScalar() == true);
    int8_t hi[25]; fill_sym(hi,1,4); hi[1*5+3] = 2;    // above the match score
    CHECK(make_kswv(6,1,6,1, 1,-4, 1, 256,256, hi)->needsScalar() == true);
}
// The decisive behavior: neutral scores the conversion cell 0, so it recovers
// exactly the mismatch penalty b relative to the symmetric matrix — and stays a
// full match score a BELOW genomic, which rewards the same cell with +a.
TEST_CASE("kswv NEON 8-bit: neutral OT scores the conversion cell 0") {
    if (!tier_supports_freed_cell()) return;   // scalar-only tier: no kernel score
    const int a = 1, b = 4;
    std::vector<uint8_t> ref  = {0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3};
    std::vector<uint8_t> read = ref;
    const int p = 10;
    ref[p]  = 1;                      // ref C
    read[p] = 3;                      // read T (the C->T conversion)

    int8_t sym[25]; fill_sym(sym, a, b);
    int8_t gen[25]; fill_sym(gen, a, b); gen[1*5+3] = (int8_t)a;   // genomic: freed to +a
    int8_t neu[25]; fill_sym(neu, a, b); neu[1*5+3] = 0;           // neutral: freed to 0

    int s_sym = score8(ref, read, a, b, sym);
    int s_gen = score8(ref, read, a, b, gen);
    int s_neu = score8(ref, read, a, b, neu);
    CHECK(s_neu - s_sym == b);        // mismatch penalty removed, no reward
    CHECK(s_gen - s_neu == a);        // genomic additionally rewards the match
}
TEST_CASE("kswv NEON 8-bit: neutral OB scores the conversion cell 0") {
    if (!tier_supports_freed_cell()) return;   // scalar-only tier: no kernel score
    const int a = 1, b = 4;
    std::vector<uint8_t> ref  = {0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3};
    std::vector<uint8_t> read = ref;
    const int p = 10;
    ref[p]  = 2;                      // ref G
    read[p] = 0;                      // read A (the G->A conversion)

    int8_t sym[25]; fill_sym(sym, a, b);
    int8_t neu[25]; fill_sym(neu, a, b); neu[2*5+0] = 0;

    CHECK(score8(ref, read, a, b, neu) - score8(ref, read, a, b, sym) == b);
}
// Neutral leaves the MIRROR cell a real mismatch (that is what keeps NM/MD
// variant-aware): freeing (C,T) to 0 must not touch (T,C).
TEST_CASE("kswv NEON 8-bit: neutral leaves the mirror cell a mismatch") {
    if (!tier_supports_freed_cell()) return;   // scalar-only tier: no kernel score
    const int a = 1, b = 4;
    std::vector<uint8_t> ref  = {0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3};
    std::vector<uint8_t> read = ref;
    const int p = 10;
    ref[p]  = 3;                      // ref T
    read[p] = 1;                      // read C → the (T,C) mirror, a real variant

    int8_t sym[25]; fill_sym(sym, a, b);
    int8_t neu[25]; fill_sym(neu, a, b); neu[1*5+3] = 0;

    CHECK(score8(ref, read, a, b, neu) == score8(ref, read, a, b, sym));
}
