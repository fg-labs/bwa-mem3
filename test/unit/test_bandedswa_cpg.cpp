// test/unit/test_bandedswa_cpg.cpp
//
// The --meth-context cpg arm of the SIMD banded extension kernels (getScores8 /
// getScores16, every tier, with and without 8-bit lane compaction). Under cpg
// the OT/OB kernel objects are built on the 7-code matrices mat7_ot / mat7_ob
// (meth_model.h) and score CpG-labelled targets: codes 0..3, plus C5 (a C whose
// next frame base is G) and G6 (a G whose previous frame base is C). Rows and
// columns 0..3 of mat7 are the plain matrix; the one freed cell is C5 x read T
// (OT) or G6 x read A (OB), worth 0 under NEUTRAL and +a under GENOMIC. The
// object learns that cell through set_meth_cpg(fr_ref_base, fr_read, fr_val).
//
// The oracle is the same object's scalarBandedSWA, which scores the labelled
// target through mat7 at stride 7 (itself checked against ksw_extend2 at m = 7
// in test_meth_model.cpp). Every SeqPair result field (score, qle, tle, gtle,
// gscore, max_off) must equal it, for OT and OB, NEUTRAL and GENOMIC, 8- and
// 16-bit, over forward and reverse-strand windows of a CpG-rich contig labelled
// by meth_label_cpg, with C->T (OT) or G->A (OB) conversions at CpG and at
// non-CpG positions. The scoring sets include the high-h0 / small-zdrop regime
// of bandedswa_high_h0_zdrop_test.cpp. 8-bit batches hold only pairs inside
// the 8-bit routing envelope (bsw8_envelope_ok in bwamem.cpp), as production
// does; the 16-bit batches take every pair.
//
// The kernel is made through make_banded_pair_wise_sw, so it runs the tier the
// host (or BWAMEM3_FORCE_TIER) selects. On aarch64 the test also runs the
// second kernel copy (src/bandedSWA.rowalt.o, BSW8_ROW_LEAN / BSW8_SKEW flipped).

#include <cstdint>
#include <cstdlib>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "doctest/doctest.h"
#include "bandedSWA.h"
#include "bntseq.h"
#include "bwa.h"
#include "bwamem.h"
#include "ext_result.h"
#include "meth_model.h"
#include "meth_orig_ref_fixture.h"
#include "simd_dispatch.h"

#if HAVE_BSW_VECTOR_8_16

#if defined(__aarch64__)
/* Defined by src/bandedSWA.rowalt.o: the kernel copy with BSW8_ROW_LEAN / BSW8_SKEW flipped. */
extern "C" IBandedPairWiseSW *make_bsw_kernel_m_rowalt(int, int, int, int, int, int,
                                                       const int8_t *, int8_t, int8_t, int, int);
#endif

static_assert(BSW_CPG_C5 == METH_C5 && BSW_CPG_G6 == METH_G6,
              "the kernels' CpG codes must be meth_model.h's labels");

namespace {

using bwa_tests::ExtResult;

struct Params {
    int a, b, o, e, zdrop, w, end_bonus;
    bool high_h0;     // draw h0 above zdrop + 1 (the high-h0 / small-zdrop regime)
    const char *label;
};

const Params kParams[] = {
    {1, 4,  6, 1, 100, 100, 5, false, "defaults"},
    {1, 4,  6, 1,   5, 100, 5, true,  "zdrop=5 high h0"},
    {1, 4,  6, 1,  20, 100, 5, true,  "zdrop=20 high h0"},
    {1, 4,  6, 1,  10,  20, 5, true,  "zdrop=10 narrow band"},
    {2, 8, 12, 2,  20, 100, 5, true,  "-A2 -B8 -O12 -E2 zdrop=20"},
    // gap open + extend byte 128: the AVX-512BW 8-bit kernel's unsigned score domain
    {1, 4, 127, 1, 100, 100, 5, false, "-O127 -E1 (gap byte 128)"},
};

const int kQlens[] = {30, 100, 150, 250};

// A CpG-rich contig: at each step a "CG" with probability 0.15, otherwise one
// random base, so well over 10% of positions are labelled C5 or G6.
std::string cpg_rich_contig(std::mt19937_64 &rng, int len)
{
    static const char *B = "ACGT";
    std::string s;
    while ((int)s.size() < len) {
        if (rng() % 100 < 15) s += "CG";
        else                  s += B[rng() % 4];
    }
    s.resize((size_t)len);
    return s;
}

// The 8-bit routing envelope (bsw8_envelope_ok in bwamem.cpp).
bool envelope8_ok(int len1, int len2, int w, int a, int zdrop, int h0)
{
    const int64_t max_step = a > 1 ? a : 1;
    const int64_t shorter = len1 < len2 ? len1 : len2;
    return len1 < MAX_SEQ_LEN8 && len2 < MAX_SEQ_LEN8 && w <= 124 &&
           zdrop + max_step <= 253 && h0 + shorter * max_step < 255 - max_step;
}

struct Corpus {
    std::vector<uint8_t> ref, qer;
    std::vector<SeqPair> pairs;
    int n = 0;
    long labelled = 0, bases = 0;   // labelled (C5/G6) target positions, all target positions
};

// n pairs: a labelled window (forward or reverse strand) as the target and a
// converted copy of its prefix as the query. hyp 1 = OT (C->T), 0 = OB (G->A).
Corpus make_corpus(const meth_test::OrigRefFixture &f, const Params &P, int width, int hyp,
                   int n, uint64_t seed)
{
    const int STRIDE = 512;
    std::mt19937_64 rng(seed);
    Corpus c;
    c.ref.assign((size_t)STRIDE * n, 0);
    c.qer.assign((size_t)STRIDE * n, 0);
    c.pairs.assign((size_t)n + 64, SeqPair{});
    const int64_t L = f.bns.l_pac;
    while (c.n < n) {
        const int qlen0 = kQlens[rng() % 4];
        const int tlen = qlen0 + (int)(rng() % 41);
        const bool rev = rng() & 1;
        const int64_t rb = (rev ? L : 0) + (int64_t)(rng() % (uint64_t)(L - tlen));
        std::vector<uint8_t> t((size_t)tlen + 1);
        int64_t got = 0;
        bns_get_seq_into(L, f.orig_pac, rb, rb + tlen, t.data(), &got);
        REQUIRE(got == tlen);
        meth_label_cpg(t.data(), tlen, rb, 0, f.orig_bns, f.orig_pac);

        // The query: the target's plain bases, with conversions at about half of
        // the convertible positions (CpG or not), ~3% substitutions, rare
        // indels and query Ns.
        std::vector<uint8_t> q;
        for (int i = 0; i < tlen && (int)q.size() < qlen0; ++i) {
            uint8_t b = meth_code_base(t[(size_t)i]);
            const unsigned u = (unsigned)(rng() % 1000);
            if (hyp == 1 && b == METH_C && u < 500) b = METH_T;
            else if (hyp == 0 && b == METH_G && u < 500) b = METH_A;
            else if (u < 530) b = (uint8_t)(rng() % 4);
            else if (u < 535) b = METH_N;
            else if (u < 542) continue;                       // deletion from the read
            else if (u < 549) q.push_back((uint8_t)(rng() % 4));   // insertion
            q.push_back(b);
        }
        while ((int)q.size() < qlen0) q.push_back((uint8_t)(rng() % 4));
        q.resize((size_t)qlen0);
        const int qlen = qlen0;

        int h0;
        if (P.high_h0) h0 = P.zdrop + 2 + (int)(rng() % 120);
        else           h0 = 1 + (int)(rng() % 60);
        if (width == 8) {
            const int shorter = tlen < qlen ? tlen : qlen;
            const int ms = P.a > 1 ? P.a : 1;
            const int hmax = 255 - ms - shorter * ms - 1;   // largest h0 inside the envelope
            if (hmax < 0) continue;                          // no 8-bit pair of this length
            if (h0 > hmax) {
                const int lo = P.high_h0 && hmax >= P.zdrop + 2 ? P.zdrop + 2 : 0;
                h0 = lo + (int)(rng() % (uint64_t)(hmax - lo + 1));
            }
            REQUIRE(envelope8_ok(tlen, qlen, P.w, P.a, P.zdrop, h0));
        }

        const int k = c.n;
        uint8_t *s1 = &c.ref[(size_t)k * STRIDE];
        uint8_t *s2 = &c.qer[(size_t)k * STRIDE];
        for (int i = 0; i < tlen; ++i) {
            s1[i] = t[(size_t)i];
            c.labelled += t[(size_t)i] >= METH_C5;
        }
        c.bases += tlen;
        for (int i = 0; i < qlen; ++i) s2[i] = q[(size_t)i];
        SeqPair &sp = c.pairs[(size_t)k];
        sp.id = k; sp.seqid = k; sp.regid = k;
        sp.len1 = tlen; sp.len2 = qlen; sp.h0 = h0;
        sp.idr = (int)((size_t)k * STRIDE); sp.idq = (int)((size_t)k * STRIDE);
        sp.score = sp.tle = sp.gtle = sp.qle = sp.gscore = sp.max_off = -1;
        ++c.n;
    }
    return c;
}

mem_opt_t *opt_for(int scoring, int a, int b)
{
    mem_opt_t *o = mem_opt_init();
    o->a = a; o->b = b;
    bwa_fill_scmat(a, b, o->mat);
    o->meth_scoring = scoring;
    mem_opt_fill_meth_mat(o);
    return o;
}

using MakeM = IBandedPairWiseSW *(*)(int, int, int, int, int, int, const int8_t *, int8_t, int8_t, int, int);

struct Outcome {
    int diffs = 0;          // pairs whose SIMD record differs from the m = 7 oracle
    int freed_matters = 0;  // pairs whose oracle score differs from the unlabelled target's
};

// One corpus through one kernel copy: the SIMD records against the scalar m = 7
// oracle on the same object.
Outcome run_one(const meth_test::OrigRefFixture &f, const Params &P, int scoring, int hyp,
                int width, int compact, MakeM make, uint64_t seed)
{
    mem_opt_t *o = opt_for(scoring, P.a, P.b);
    const int8_t *mat7 = hyp ? o->mat7_ot : o->mat7_ob;
    const int fr_ref = hyp ? METH_C : METH_G, fr_read = hyp ? METH_T : METH_A;
    const int fr_val = mat7[(fr_ref + 4) * METH_M + fr_read];
    REQUIRE(fr_val == (scoring == MEM_METH_SCORING_NEUTRAL ? 0 : P.a));

    std::unique_ptr<IBandedPairWiseSW> k;
    if (make)
        k.reset(make(P.o, P.e, P.o, P.e, P.zdrop, P.end_bonus, mat7, (int8_t)P.a, (int8_t)P.b, 1, METH_M));
    else
        k = make_banded_pair_wise_sw(P.o, P.e, P.o, P.e, P.zdrop, P.end_bonus, mat7,
                                     (int8_t)P.a, (int8_t)P.b, 1, METH_M);
    k->set_meth_cpg(fr_ref, fr_read, fr_val);
    k->set_lane_compaction(compact, BSW_COMPACT_MIN_W_DEFAULT);

    Corpus c = make_corpus(f, P, width, hyp, 400, seed);
    CHECK((double)c.labelled >= 0.10 * (double)c.bases);

    std::vector<ExtResult> oracle((size_t)c.n);
    Outcome out;
    for (int i = 0; i < c.n; ++i) {
        const SeqPair &sp = c.pairs[(size_t)i];
        const uint8_t *t = &c.ref[(size_t)sp.idr], *q = &c.qer[(size_t)sp.idq];
        ExtResult &R = oracle[(size_t)i];
        R.score = k->scalarBandedSWA(sp.len2, q, sp.len1, t, P.w, sp.h0,
                                     &R.qle, &R.tle, &R.gtle, &R.gscore, &R.max_off);
        std::vector<uint8_t> tp(t, t + sp.len1);
        meth_unlabel(tp.data(), sp.len1);
        int u[5];
        if (k->scalarBandedSWA(sp.len2, q, sp.len1, tp.data(), P.w, sp.h0,
                                   &u[0], &u[1], &u[2], &u[3], &u[4]) != R.score)
            ++out.freed_matters;
    }
    if (width == 8) k->getScores8(c.pairs.data(), c.ref.data(), c.qer.data(), c.n, 1, P.w);
    else            k->getScores16(c.pairs.data(), c.ref.data(), c.qer.data(), c.n, 1, P.w);
    // On the tier with two 8-bit score domains (AVX-512BW), the batch ran in the one its gap
    // bytes select, so both domains' CpG arms are covered across kParams.
    const BswScoreDomainCounts dc = k->bsw8_score_domain_counts();
    if (width == 8 && dc.signed_batches + dc.unsigned_batches > 0)
        CHECK((bsw8_gap_bytes_fit_int8(P.o, P.e, P.o, P.e) ? dc.signed_batches : dc.unsigned_batches) > 0);
    for (int i = 0; i < c.n; ++i) {
        const SeqPair &p = c.pairs[(size_t)i];
        const ExtResult &R = oracle[(size_t)p.id];
        if (p.score != R.score || p.qle != R.qle || p.tle != R.tle || p.gtle != R.gtle ||
            p.gscore != R.gscore || p.max_off != R.max_off) {
            if (out.diffs < 3)
                MESSAGE("diff pair " << p.id << " len1 " << p.len1 << " len2 " << p.len2
                        << " h0 " << p.h0 << ": simd " << p.score << "/" << p.qle << "/" << p.tle
                        << "/" << p.gtle << "/" << p.gscore << "/" << p.max_off << " oracle "
                        << R.score << "/" << R.qle << "/" << R.tle << "/" << R.gtle << "/"
                        << R.gscore << "/" << R.max_off);
            ++out.diffs;
        }
    }
    free(o);
    return out;
}

// The grid for one scoring set: every kernel copy, OT and OB, NEUTRAL and GENOMIC, 8-bit with
// and without compaction, and 16-bit. One case per scoring set keeps each case short.
void check_params(const Params &P, uint64_t seed)
{
    bwamem3_simd_init();
    MESSAGE("bandedSWA cpg arm (" << std::string(P.label) << "): active SIMD tier "
            << std::string(bwamem3_simd_tier_name(bwamem3_simd_tier())));
    std::mt19937_64 rng(20261028);
    const std::string contig = cpg_rich_contig(rng, 200000);
    meth_test::OrigRefFixture f(contig);

    struct Copy { MakeM make; const char *name; };
    std::vector<Copy> copies = {{nullptr, "dispatched tier"}};
#if defined(__aarch64__)
    copies.push_back({make_bsw_kernel_m_rowalt, "neon (other BSW8_ROW_LEAN / BSW8_SKEW setting)"});
#endif
    const int scorings[] = {MEM_METH_SCORING_NEUTRAL, MEM_METH_SCORING_GENOMIC};
    long total_freed = 0;
    for (const Copy &cp : copies)
    for (int sc : scorings)
    for (int hyp = 0; hyp < 2; ++hyp)
    for (int width : {8, 16})
    for (int compact : {0, 1}) {
        if (width == 16 && compact) continue;   // compaction is an 8-bit driver
        // compaction off, or on at 8 groups per superblock (the tiers without a
        // compaction driver ignore it) at the production band gate
        const int groups = !compact ? 0 : BSW_COMPACT_GROUPS_NEON;
        CAPTURE(cp.name); CAPTURE(P.label); CAPTURE(sc); CAPTURE(hyp); CAPTURE(width);
        CAPTURE(compact);
        const Outcome r = run_one(f, P, sc, hyp, width, groups, cp.make, seed++);
        CHECK(r.diffs == 0);
        total_freed += r.freed_matters;
    }
    MESSAGE("bandedSWA cpg arm (" << std::string(P.label) << "): " << total_freed << " pairs scored by the CpG cell");
    // The freed cell is exercised: on many pairs the labelled target scores
    // differently from the same target with its labels removed.
    CHECK(total_freed > 500);
}

}  // namespace

TEST_CASE("bandedSWA cpg arm == scalar m=7 oracle: defaults" * doctest::test_suite("unit/bandedswa_cpg"))
{
    check_params(kParams[0], 1);
}

TEST_CASE("bandedSWA cpg arm == scalar m=7 oracle: zdrop=5 high h0" * doctest::test_suite("unit/bandedswa_cpg"))
{
    check_params(kParams[1], 101);
}

TEST_CASE("bandedSWA cpg arm == scalar m=7 oracle: zdrop=20 high h0" * doctest::test_suite("unit/bandedswa_cpg"))
{
    check_params(kParams[2], 201);
}

TEST_CASE("bandedSWA cpg arm == scalar m=7 oracle: zdrop=10 narrow band"
          * doctest::test_suite("unit/bandedswa_cpg"))
{
    check_params(kParams[3], 301);
}

TEST_CASE("bandedSWA cpg arm == scalar m=7 oracle: -A2 -B8 -O12 -E2 zdrop=20"
          * doctest::test_suite("unit/bandedswa_cpg"))
{
    check_params(kParams[4], 401);
}

TEST_CASE("bandedSWA cpg arm == scalar m=7 oracle: -O127 -E1 (gap byte 128)"
          * doctest::test_suite("unit/bandedswa_cpg"))
{
    check_params(kParams[5], 501);
}

#endif // HAVE_BSW_VECTOR_8_16
