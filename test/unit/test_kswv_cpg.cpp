// test/unit/test_kswv_cpg.cpp
//
// The --meth-context cpg gate of the batched mate-rescue kernels (kswv
// getScores8 / getScores16, at the runtime tier: NEON, AVX2, AVX-512BW).
//
// Under cpg the rescue window is labelled by meth_label_cpg in its fetch frame:
// reference codes 0..3, N (4), plus C5 (a C whose next frame base is G) and G6
// (a G whose previous frame base is C). The kswv object is built from the
// 7-code matrix (mat7_ot / mat7_ob, meth_model.h) with make_kswv(..., mat7,
// meth_cpg = true). Rows and columns 0..4 of mat7 are the plain matrix; the one
// freed cell is C5 x read T (OT) or G6 x read A (OB), worth 0 under NEUTRAL and
// +a under GENOMIC. A conversion against a plain C (OT) or G (OB) is an
// ordinary mismatch.
//
// The oracle is scalar ksw_align2 at m = 7 with the same matrix on the same
// labelled window. The batch runs both production passes (phase 0, then the
// in-place reversal and phase 1 on the survivors, as mem_sam_pe_batch_run does)
// and every result field must equal the oracle's: score, te, qe, score2 (with
// te2 where both report one) and tb / qb, with the tolerances the existing kswv
// scalar-parity tests apply (kswr_cmp.h) to tb / qb: none below min_seed_len.
// score2 is compared exactly (see fields_eq).
//
// Grid: OT and OB, NEUTRAL and GENOMIC, 8- and 16-bit, several scorings
// (including unequal gap costs, which take the non-FScan bodies), and the
// BWA3_RESCUE_* toggles that select the other kernel bodies. Windows come from
// the forward and reverse strands of a CpG-rich contig with internal repeats
// (so score2 is exercised), with reference Ns in some windows; queries carry
// conversions at CpG and non-CpG positions, substitutions, indels and Ns.
// On SSE4.1 / SSE4.2 / AVX the kernels have no freed-cell override, so the
// object must report needsScalar() (the caller's scalar fallback) instead.

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "doctest/doctest.h"
#include "bntseq.h"
#include "bwa.h"
#include "bwamem.h"
#include "ksw.h"
#include "kswv.h"
#include "kswr_cmp.h"
#include "kswv_runner.h"
#include "meth_model.h"
#include "meth_orig_ref_fixture.h"
#include "scoped_env.h"
#include "seqpair_batch.h"
#include "simd_dispatch.h"

#if BWA_TESTS_HAVE_KSWV

namespace {

using bwa_tests::BatchBuffers;
using bwa_tests::ScopedEnv;
using bwa_tests::TestPair;

struct Params {
    int a, b, o_del, e_del, o_ins, e_ins;
    const char *label;
};

const Params kParams[] = {
    {1, 4,  6, 1,  6, 1, "defaults"},
    {2, 8, 12, 2, 12, 2, "-A2 -B8 -O12 -E2"},
    {1, 4,  5, 2,  6, 1, "equal open+extend sums, different costs (FScan)"},
    {1, 4,  6, 1,  4, 1, "unequal gap sums (original bodies)"},
};

// Kernel-body selectors read on every kernel call (rescue_env_on).
struct Toggle { const char *name; const char *value; };
const Toggle kToggles[] = {
    {nullptr, nullptr},                 // the default bodies
    {"BWA3_RESCUE_FSCAN", "0"},
    {"BWA3_RESCUE_USQADD", "0"},
    {"BWA3_RESCUE_ROWPAIR", "0"},
    {"BWA3_RESCUE_LAZYQE", "0"},
};

const int kMinSeedLen = 19;

bool tier_has_freed_kernel()
{
    bwamem3_simd_init();
    const int t = bwamem3_simd_tier();
    return t == BWAMEM3_TIER_NEON || t == BWAMEM3_TIER_AVX2 || t == BWAMEM3_TIER_AVX512BW;
}

// A CpG-rich contig with internal repeats: at each step a "CG" with
// probability 0.15, otherwise one random base; every 1500 bases the previous
// 300 are repeated with ~3% substitutions, so a rescue window that holds both
// copies has a real second-best alignment.
std::string cpg_rich_contig(std::mt19937_64 &rng, int len)
{
    static const char *B = "ACGT";
    std::string s;
    while ((int)s.size() < len) {
        if (s.size() >= 1500 && s.size() % 1500 < 2) {
            const std::string blk = s.substr(s.size() - 300, 300);
            for (char c : blk) s += (rng() % 100 < 3) ? B[rng() % 4] : c;
        }
        if (rng() % 100 < 15) s += "CG";
        else                  s += B[rng() % 4];
    }
    s.resize((size_t)len);
    return s;
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

// The rescue xtra word: production's (bwamem_pair.cpp), KSW_XBYTE for 8-bit.
int xtra_for(int a, bool u8)
{
    return KSW_XSUBO | KSW_XSTART | (u8 ? KSW_XBYTE : 0) | (kMinSeedLen * a);
}

// The largest 8-bit query length the router admits (matesw_use_u8, with
// shift = b for these scorings).
int max_u8_qlen(int a, int b)
{
    int q = 0;
    while ((q + 1) * a + b <= 254 && (q + 1) * a < 250) ++q;
    return q;
}

struct Corpus {
    std::vector<TestPair> pairs;
    long labelled = 0, bases = 0, ref_n = 0;
};

// n rescue jobs: a labelled window (forward or reverse strand) as the target
// and, as the query, a converted copy of a segment of it (or, for one job in
// ten, an unrelated random read). hyp 1 = OT (C->T), 0 = OB (G->A): most CpG
// sites and ~1.5% of the other Cs (OT) or Gs (OB) are converted, plus ~3%
// substitutions, rare indels and Ns.
Corpus make_corpus(const meth_test::OrigRefFixture &f, const Params &P, int width, int hyp,
                   int n, uint64_t seed)
{
    std::mt19937_64 rng(seed);
    Corpus c;
    const int64_t L = f.bns.l_pac;
    const int qmax8 = max_u8_qlen(P.a, P.b);
    const int qlens8[] = {30, 100, 150, qmax8};
    const int qlens16[] = {30, 100, 150, 250, 300};
    while ((int)c.pairs.size() < n) {
        int qlen = width == 8 ? qlens8[rng() % 4] : qlens16[rng() % 5];
        if (width == 8 && qlen > qmax8) qlen = qmax8;
        const int tlen = qlen + 50 + (int)(rng() % 600);
        const bool rev = rng() & 1;
        const int64_t rb = (rev ? L : 0) + (int64_t)(rng() % (uint64_t)(L - tlen));
        std::vector<uint8_t> t((size_t)tlen);
        int64_t got = 0;
        bns_get_seq_into(L, f.orig_pac, rb, rb + tlen, t.data(), &got);
        REQUIRE(got == tlen);
        meth_label_cpg(t.data(), tlen, rb, 0, f.orig_bns, f.orig_pac);
        // Reference Ns (code 4) in one window in eight: --meth references hold
        // none, but the kernels must keep scoring them as w_ambig.
        if (rng() % 8 == 0)
            for (int i = 0; i < tlen; ++i)
                if (rng() % 100 < 2) t[(size_t)i] = METH_N;

        std::vector<uint8_t> q;
        if (rng() % 10 == 0) {
            for (int i = 0; i < qlen; ++i) q.push_back((uint8_t)(rng() % 4));
        } else {
            const int start = (int)(rng() % (uint64_t)(tlen - qlen));
            for (int i = start; i < tlen && (int)q.size() < qlen; ++i) {
                const uint8_t code = t[(size_t)i];
                uint8_t b = meth_code_base(code);
                if (b == METH_N) b = (uint8_t)(rng() % 4);
                // Conversions: most labelled (CpG) sites, a few plain ones.
                const unsigned u = (unsigned)(rng() % 1000);
                const unsigned conv = code >= METH_C5 ? 600 : 15;
                if (hyp == 1 && b == METH_C && u < conv) b = METH_T;
                else if (hyp == 0 && b == METH_G && u < conv) b = METH_A;
                else if (u < 630 && u >= 600) b = (uint8_t)(rng() % 4);
                else if (u >= 630 && u < 635) b = METH_N;
                else if (u >= 635 && u < 641) continue;                 // deletion from the read
                else if (u >= 641 && u < 647) q.push_back((uint8_t)(rng() % 4));   // insertion
                q.push_back(b);
            }
            while ((int)q.size() < qlen) q.push_back((uint8_t)(rng() % 4));
            q.resize((size_t)qlen);
        }
        for (uint8_t x : t) { c.labelled += x >= METH_C5; c.ref_n += x == METH_N; }
        c.bases += tlen;
        TestPair tp;
        tp.ref = t;
        tp.qry = q;
        tp.tag = std::string(rev ? "rev" : "fwd") + " rb=" + std::to_string((long long)rb);
        c.pairs.push_back(tp);
    }
    return c;
}

// Both production passes through one kswv object: phase 0, then the survivors'
// in-place reversal and phase 1 (BatchBuffers::prepare_phase1).
std::vector<kswr_t> run_batch(Ikswv &k, const std::vector<TestPair> &pairs, int xtra, int width)
{
    BatchBuffers bb(pairs, xtra);
    if (width == 8) k.getScores8(bb.pairs(), bb.ref_buf(), bb.qer_buf(), bb.aln(), bb.n(), 1, 0);
    else            k.getScores16(bb.pairs(), bb.ref_buf(), bb.qer_buf(), bb.aln(), bb.n(), 1, 0);
    const int survivors = bb.prepare_phase1();
    if (survivors > 0) {
        if (width == 8) k.getScores8(bb.pairs(), bb.ref_buf(), bb.qer_buf(), bb.aln(), survivors, 1, 1);
        else            k.getScores16(bb.pairs(), bb.ref_buf(), bb.qer_buf(), bb.aln(), survivors, 1, 1);
    }
    return std::vector<kswr_t>(bb.aln(), bb.aln() + pairs.size());
}

// score2 is compared exactly, not through kswr_score2_eq (which takes a batched 0 as the scalar's
// -1): the rescue copies score2 into csub, and a 0 there raises the MAPQ cap
// raw_mapq(score - csub) by one step against the -1 that ksw_align2 reports.
bool fields_eq(const kswr_t &scalar, const kswr_t &simd)
{
    using namespace bwa_tests;
    if (!kswr_score_eq(scalar, simd) || !kswr_ends_eq(scalar, simd) ||
        !kswr_coords_eq(scalar, simd, kMinSeedLen))
        return false;
    if (scalar.score > 0 && scalar.score2 != simd.score2) return false;
    if (scalar.score > 0 && scalar.score2 > 0 && scalar.te2 != simd.te2)
        return false;
    return true;
}

struct Outcome {
    int diffs = 0;          // jobs (summed over the kernel-body toggles) whose kernel result
                            // differs from the m = 7 oracle
    int freed_matters = 0;  // jobs whose oracle score differs from the unlabelled window's
    int with_score2 = 0;    // jobs whose oracle reports a second-best score
};

// One corpus, scored once by the scalar oracle, then through the kernel under every kernel-body
// toggle: the toggles select only the SIMD bodies, so the oracle is shared.
Outcome run_one(const meth_test::OrigRefFixture &f, const Params &P, int scoring, int hyp,
                int width, uint64_t seed)
{
    mem_opt_t *o = opt_for(scoring, P.a, P.b);
    const int8_t *mat7 = hyp ? o->mat7_ot : o->mat7_ob;
    const int fr_ref = hyp ? METH_C5 : METH_G6, fr_read = hyp ? METH_T : METH_A;
    REQUIRE(mat7[fr_ref * METH_M + fr_read] == (scoring == MEM_METH_SCORING_NEUTRAL ? 0 : P.a));

    Corpus c = make_corpus(f, P, width, hyp, 150, seed);
    CHECK((double)c.labelled >= 0.10 * (double)c.bases);

    int32_t maxr = 0, maxq = 0;
    for (const TestPair &tp : c.pairs) {
        if ((int32_t)tp.ref.size() > maxr) maxr = (int32_t)tp.ref.size();
        if ((int32_t)tp.qry.size() > maxq) maxq = (int32_t)tp.qry.size();
    }
    std::unique_ptr<Ikswv> k = make_kswv(P.o_del, P.e_del, P.o_ins, P.e_ins,
                                         (int8_t)P.a, (int8_t)-P.b, 1, maxr, maxq, mat7, true);
    Outcome out;
    if (!tier_has_freed_kernel()) {
        // No freed-cell override: the caller must take the scalar ksw_align2.
        CHECK(k->needsScalar());
        free(o);
        return out;
    }
    REQUIRE_FALSE(k->needsScalar());

    const int xtra = xtra_for(P.a, width == 8);
    std::vector<kswr_t> oracle;
    for (const TestPair &tp : c.pairs) {
        std::vector<uint8_t> t = tp.ref, q = tp.qry;   // ksw_align2 reverses in place
        kswr_t r = ksw_align2((int)q.size(), q.data(), (int)t.size(), t.data(), METH_M, mat7,
                              P.o_del, P.e_del, P.o_ins, P.e_ins, xtra, 0);
        oracle.push_back(r);
        out.with_score2 += r.score > 0 && r.score2 > 0;
        std::vector<uint8_t> tu = tp.ref, qu = tp.qry;
        meth_unlabel(tu.data(), (int64_t)tu.size());
        kswr_t u = ksw_align2((int)qu.size(), qu.data(), (int)tu.size(), tu.data(), METH_M, mat7,
                              P.o_del, P.e_del, P.o_ins, P.e_ins, xtra, 0);
        out.freed_matters += u.score != r.score;
    }

    for (const Toggle &tg : kToggles) {
        std::unique_ptr<ScopedEnv> env;
        if (tg.name) env.reset(new ScopedEnv(tg.name, tg.value));
        const char *body = tg.name ? tg.name : "default bodies";
        const std::vector<kswr_t> simd = run_batch(*k, c.pairs, xtra, width);
        REQUIRE(simd.size() == oracle.size());
        for (size_t i = 0; i < simd.size(); ++i) {
            if (fields_eq(oracle[i], simd[i])) continue;
            if (out.diffs < 3) {
                const kswr_t &a = oracle[i], &b = simd[i];
                MESSAGE("diff job " << i << " (" << body << ", " << c.pairs[i].tag << ", tlen "
                        << c.pairs[i].ref.size() << ", qlen " << c.pairs[i].qry.size() << "): oracle "
                        << a.score << "/" << a.te << "/" << a.qe << "/" << a.score2 << "/" << a.te2
                        << "/" << a.tb << "/" << a.qb << " simd " << b.score << "/" << b.te << "/"
                        << b.qe << "/" << b.score2 << "/" << b.te2 << "/" << b.tb << "/" << b.qb);
            }
            ++out.diffs;
        }
    }
    free(o);
    return out;
}

// The grid for one scoring set: OT and OB, NEUTRAL and GENOMIC, 8- and 16-bit, each corpus
// through every kernel-body toggle. One case per scoring set keeps each case short.
void check_params(const Params &P, uint64_t seed)
{
    bwamem3_simd_init();
    MESSAGE("kswv cpg gate (" << std::string(P.label) << "): active SIMD tier "
            << std::string(bwamem3_simd_tier_name(bwamem3_simd_tier())));
    std::mt19937_64 rng(20261028);
    const std::string contig = cpg_rich_contig(rng, 200000);
    meth_test::OrigRefFixture f(contig);

    const int scorings[] = {MEM_METH_SCORING_NEUTRAL, MEM_METH_SCORING_GENOMIC};
    long total_freed = 0, total_score2 = 0;
    for (int sc : scorings)
    for (int hyp = 0; hyp < 2; ++hyp)
    for (int width : {8, 16}) {
        CAPTURE(P.label); CAPTURE(sc); CAPTURE(hyp); CAPTURE(width);
        const Outcome r = run_one(f, P, sc, hyp, width, seed++);
        CHECK(r.diffs == 0);
        total_freed += r.freed_matters;
        total_score2 += r.with_score2;
    }
    MESSAGE("kswv cpg gate (" << std::string(P.label) << "): " << total_freed
            << " jobs scored by the CpG cell, " << total_score2 << " with a second-best score");
    if (tier_has_freed_kernel()) {
        // Non-vacuity: the CpG cell decides many scores, and score2 is exercised.
        CHECK(total_freed > 500);
        CHECK(total_score2 > 250);
    }
}

}  // namespace

TEST_CASE("kswv cpg gate == ksw_align2 m=7: defaults" * doctest::test_suite("unit/kswv_cpg"))
{
    check_params(kParams[0], 1);
}

TEST_CASE("kswv cpg gate == ksw_align2 m=7: -A2 -B8 -O12 -E2" * doctest::test_suite("unit/kswv_cpg"))
{
    check_params(kParams[1], 101);
}

TEST_CASE("kswv cpg gate == ksw_align2 m=7: equal open+extend sums, different costs (FScan)"
          * doctest::test_suite("unit/kswv_cpg"))
{
    check_params(kParams[2], 201);
}

TEST_CASE("kswv cpg gate == ksw_align2 m=7: unequal gap sums (original bodies)"
          * doctest::test_suite("unit/kswv_cpg"))
{
    check_params(kParams[3], 301);
}

TEST_CASE("kswv cpg: matrices the gate cannot express take the scalar fallback"
          * doctest::test_suite("unit/kswv_cpg"))
{
    mem_opt_t *o = opt_for(MEM_METH_SCORING_GENOMIC, 1, 4);
    const bool simd = tier_has_freed_kernel();
    for (int hyp = 0; hyp < 2; ++hyp) {
        const int8_t *mat7 = hyp ? o->mat7_ot : o->mat7_ob;
        CAPTURE(hyp);
        CHECK(make_kswv(6, 1, 6, 1, 1, -4, 1, 256, 256, mat7, true)->needsScalar() == !simd);

        int8_t m[49];
        const int frow = hyp ? METH_C5 : METH_G6;   // the freed (labelled) row
        // A freed cell on a plain base row: the ACGT block is not symmetric.
        memcpy(m, mat7, sizeof m);
        m[(hyp ? METH_C : METH_G) * METH_M + (hyp ? METH_T : METH_A)] = 1;
        CHECK(make_kswv(6, 1, 6, 1, 1, -4, 1, 256, 256, m, true)->needsScalar());
        // An N row / N column entry that is not w_ambig.
        memcpy(m, mat7, sizeof m);
        m[METH_N * METH_M + METH_A] = 0;
        CHECK(make_kswv(6, 1, 6, 1, 1, -4, 1, 256, 256, m, true)->needsScalar());
        memcpy(m, mat7, sizeof m);
        m[METH_C * METH_M + METH_N] = 0;
        CHECK(make_kswv(6, 1, 6, 1, 1, -4, 1, 256, 256, m, true)->needsScalar());
        // A labelled row whose read-N column differs from its base row.
        memcpy(m, mat7, sizeof m);
        m[frow * METH_M + METH_N] = 0;
        CHECK(make_kswv(6, 1, 6, 1, 1, -4, 1, 256, 256, m, true)->needsScalar());
        // Two freed cells on the CpG rows.
        memcpy(m, mat7, sizeof m);
        m[METH_C5 * METH_M + METH_T] = 1;
        m[METH_G6 * METH_M + METH_A] = 1;
        CHECK(make_kswv(6, 1, 6, 1, 1, -4, 1, 256, 256, m, true)->needsScalar());
        // The freed row differs from its base row in a second cell as well.
        memcpy(m, mat7, sizeof m);
        m[frow * METH_M + (hyp ? METH_A : METH_C)] = 0;
        CHECK(make_kswv(6, 1, 6, 1, 1, -4, 1, 256, 256, m, true)->needsScalar());
        // A freed value outside [w_mismatch, w_match].
        memcpy(m, mat7, sizeof m);
        m[(hyp ? METH_C5 : METH_G6) * METH_M + (hyp ? METH_T : METH_A)] = 2;
        CHECK(make_kswv(6, 1, 6, 1, 1, -4, 1, 256, 256, m, true)->needsScalar());
        // No freed cell at all.
        memcpy(m, mat7, sizeof m);
        m[(hyp ? METH_C5 : METH_G6) * METH_M + (hyp ? METH_T : METH_A)] = -4;
        CHECK(make_kswv(6, 1, 6, 1, 1, -4, 1, 256, 256, m, true)->needsScalar());
    }
    free(o);
}

#endif // BWA_TESTS_HAVE_KSWV
