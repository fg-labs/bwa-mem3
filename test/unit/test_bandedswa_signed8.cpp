// test/unit/test_bandedswa_signed8.cpp
//
// The x86 8-bit banded extension kernel (AVX-512BW; see the score-domain note above its
// MAIN_CODE8_CORE in bandedSWA.cpp) carries its H/E/F bytes as score ^ 0x80 when every gap
// penalty byte fits an int8, and as the plain unsigned score otherwise. The two are separate
// instantiations of one kernel template, chosen per batch by bsw8_gap_bytes_fit_int8(). This
// suite pins the chooser at its boundaries, runs synthetic extension batches through both
// instantiations on every x86 8-bit tier the host has (lane compaction off and on) against the
// scalar reference, and checks through the object's score-domain counters that each scoring
// really ran the instantiation the chooser names. Every scoring keeps each gap penalty byte
// equal to the penalty, so the scalar reference scores the same problem the kernel does, and
// every pair stays inside the 8-bit routing envelope (h0 + min(len) * a < 254, tlen >= qlen).

#include <cstdint>
#include <cstdlib>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "doctest/doctest.h"
#include "bandedSWA.h"
#include "bsw_batch.h"
#include "ext_ladder.h"

TEST_CASE("bsw8: the signed-domain chooser follows each gap penalty byte at the int8 boundary"
          * doctest::test_suite("unit/bandedswa-signed8")) {
    // production scoring and the largest penalties that still fit: signed domain
    CHECK(bsw8_gap_bytes_fit_int8(6, 1, 6, 1));
    CHECK(bsw8_gap_bytes_fit_int8(126, 1, 126, 1));      // o + e == 127 on both sides
    CHECK(bsw8_gap_bytes_fit_int8(0, 127, 0, 127));      // e == 127 == o + e
    // one byte past 127 on any single clause, deletion or insertion side: unsigned domain.
    // Each case is rejected by exactly one of the four clauses, so losing any clause fails here.
    CHECK_FALSE(bsw8_gap_bytes_fit_int8(127, 1, 6, 1));      // o_del + e_del == 128, e_del fits
    CHECK_FALSE(bsw8_gap_bytes_fit_int8(6, 1, 127, 1));      // o_ins + e_ins == 128, e_ins fits
    CHECK_FALSE(bsw8_gap_bytes_fit_int8(250, 130, 6, 1));    // e_del byte 130; the sum 380 wraps to 124
    CHECK_FALSE(bsw8_gap_bytes_fit_int8(6, 1, 250, 130));    // e_ins byte 130; the sum wraps to 124
    // several clauses at once
    CHECK_FALSE(bsw8_gap_bytes_fit_int8(6, 128, 6, 1));      // e_del == 128 and o + e == 134
    CHECK_FALSE(bsw8_gap_bytes_fit_int8(200, 50, 200, 50));  // both sums 250
    // the decision is on the byte of each value: o + e == 260 is the byte 4 and fits
    CHECK(bsw8_gap_bytes_fit_int8(250, 10, 250, 10));
}

#if (defined(__x86_64__) || defined(__i386__)) && HAVE_BSW_VECTOR_8_16

namespace {

using bwa_tests::ExtPair;
using bwa_tests::ExtResult;
using bwa_tests::ExtScoring;
using bwa_tests::score_bsw_batch;
using Factory = IBandedPairWiseSW *(*)(int, int, int, int, int, int, const int8_t *, int8_t, int8_t, int);

struct Tier { const char *name; Factory make; };

std::vector<Tier> tiers_8bit()
{
    std::vector<Tier> ts;
    __builtin_cpu_init();
    if (__builtin_cpu_supports("avx2")) ts.push_back({"avx2", make_bsw_kernel_avx2});
    if (__builtin_cpu_supports("avx512bw")) ts.push_back({"avx512bw", make_bsw_kernel_avx512bw});
    return ts;
}

/// One scoring regime: the gap penalties (each byte equal to the penalty), the matrix shape
/// (0 symmetric, 1 one freed cell = the rank-1 prepass, 2 two freed cells = the generic
/// prepass), the seed and length ranges, and whether the chooser must take the signed domain.
struct Row {
    const char *what;
    int o_del, e_del, o_ins, e_ins;
    int matrix;
    int h0_lo, h0_hi, len_lo, len_hi;
    bool fits;
};

const int kA = 1, kB = 4, kZdrop = 100, kEndBonus = 5, kW = 100, kPairs = 128;

ExtScoring scoring(const Row &r)
{
    ExtScoring sc;
    sc.a = kA; sc.b = kB; sc.o_del = r.o_del; sc.e_del = r.e_del; sc.o_ins = r.o_ins; sc.e_ins = r.e_ins;
    sc.zdrop = kZdrop; sc.pen_clip = kEndBonus;
    sc.fill_mat();                                  // match a, mismatch -b, N -1
    if (r.matrix >= 1) sc.mat[1 * 5 + 3] = (int8_t) sc.a;   // ref C x read T freed (rank-1)
    if (r.matrix >= 2) sc.mat[2 * 5 + 0] = (int8_t) sc.a;   // and ref G x read A (generic)
    return sc;
}

/// kPairs synthetic pairs: the read copies a ref prefix with sparse mismatches, N and short
/// indels (an indel in the first fifth of a long read leaves enough matches after it for a
/// gap of cost 128 to pay, so the unsigned instantiation's gap paths are not idle). h0 is
/// drawn from the row's range, with every eighth pair at the range's low end so a seed of 0
/// is covered where the range allows it. Every pair satisfies the 8-bit routing envelope.
std::vector<ExtPair> make_pairs(const Row &r, unsigned long seed)
{
    std::vector<ExtPair> ps(kPairs);
    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<int> lenD(r.len_lo, r.len_hi), h0D(r.h0_lo, r.h0_hi);
    for (int c = 0; c < kPairs; ++c) {
        ExtPair &p = ps[c];
        const int la = lenD(rng), lb = lenD(rng);
        const int len1 = la > lb ? la : lb, len2 = la > lb ? lb : la;   // target >= query
        p.h0 = (c % 8 == 0) ? r.h0_lo : h0D(rng);
        REQUIRE(p.h0 + len2 * kA < 254);   // bsw8_envelope_ok: no row max can reach the byte ceiling
        p.target.resize(len1); p.query.resize(len2);
        for (int i = 0; i < len1; ++i) p.target[i] = (uint8_t) (rng() % 4);
        int i1 = 0;   // the ref position the read copies
        for (int i = 0; i < len2; ++i) {
            const int roll = (int) (rng() % 1000);
            const bool early = i < len2 / 5;
            if (roll < 10 && early && i1 + 1 < len1) i1++;        // 1%: a ref base skipped (deletion)
            uint8_t base = p.target[i1 < len1 ? i1 : len1 - 1];
            if (roll >= 10 && roll < 40) base = (uint8_t) (rng() % 4); // 3%: mismatch
            else if (roll >= 40 && roll < 60) base = 4;                // 2%: N
            if (roll < 990 || !early) i1++;                            // else 1%: the ref base repeats (insertion)
            p.query[i] = base;
        }
    }
    return ps;
}

std::vector<ExtResult> scalar_reference(IBandedPairWiseSW &bsw, const std::vector<ExtPair> &ps)
{
    std::vector<ExtResult> out(ps.size());
    for (size_t c = 0; c < ps.size(); ++c) {
        ExtResult &o = out[c];
        o.score = bsw.scalarBandedSWA((int) ps[c].query.size(), ps[c].query.data(),
                                      (int) ps[c].target.size(), ps[c].target.data(), kW, ps[c].h0,
                                      &o.qle, &o.tle, &o.gtle, &o.gscore, &o.max_off);
    }
    return out;
}

bool env_set(const char *name)
{
    const char *v = std::getenv(name);
    return v != nullptr && v[0] != '\0' && std::string(v) != "0";
}

}  // namespace

// Rows 1-3 take the signed instantiation, 4-7 the unsigned one. Seeds of 0 and of 140 (a byte
// above 127) are covered where the lengths allow; the rank-1 and generic matrices run in both
// domains; row 4 is the long-read regime in which gaps of cost 128 engage. The scalar
// reference is computed once per row and shared by every (tier, compaction) leg.
TEST_CASE("bsw8: both score-domain instantiations match the scalar reference on every x86 8-bit tier"
          * doctest::test_suite("unit/bandedswa-signed8")) {
    const std::vector<Tier> ts = tiers_8bit();
    bool has512 = false;
    for (const Tier &t : ts) has512 = has512 || std::string(t.name) == "avx512bw";
    const bool require512 = env_set("BWA_REQUIRE_SIGNED8");
    if (!has512 && require512)
        FAIL("bandedSWA signed8: BWA_REQUIRE_SIGNED8 is set but this host has no AVX-512BW; the "
             "signed-domain instantiation cannot run");
    if (ts.empty()) {
        MESSAGE("bandedSWA signed8: this host has no x86 8-bit SIMD tier (needs avx2 or avx512bw); skipped");
        return;
    }
    if (!has512)
        MESSAGE("bandedSWA signed8: no AVX-512BW on this host; the signed-domain instantiation is not exercised");

    const Row rows[] = {
        {"production -O 6 -E 1, symmetric",            6,   1,   6,   1, 0, 40,  40,  40, 150, true},
        {"-O 126 -E 1 (byte 127), rank-1 matrix",      126, 1,   126, 1, 1, 0,   60,  40, 150, true},
        {"-O 6 -E 1, generic matrix, seeds up to 140", 6,   1,   6,   1, 2, 0,   140, 40, 100, true},
        {"-O 127 -E 1 (byte 128), long reads, gaps",   127, 1,   127, 1, 0, 0,   10,  160, 240, false},
        {"-O 127,6 -E 1 (deletion side), seeds 0..40", 127, 1,   6,   1, 0, 0,   40,  40, 150, false},
        {"-E 1,128 (insertion extend), generic matrix", 6,  1,   6,   128, 2, 0,  140, 40, 100, false},
        {"-O 200 -E 50 (byte 250) both sides, rank-1", 200, 50,  200, 50, 1, 0,  60,  40, 150, false},
    };
    int rowi = 0;
    for (const Row &r : rows) {
        rowi++;
        CHECK(bsw8_gap_bytes_fit_int8(r.o_del, r.e_del, r.o_ins, r.e_ins) == r.fits);
        const ExtScoring sc = scoring(r);
        const unsigned long seed = 0x5EEDULL * rowi + 1000003ULL * r.o_del + 7919ULL * r.e_del +
                                   104729ULL * r.o_ins + 1299709ULL * r.e_ins;
        const std::vector<ExtPair> ps = make_pairs(r, seed);
        std::vector<const ExtPair *> pp;
        for (const ExtPair &p : ps) pp.push_back(&p);

        // the scalar reference once per row, through the first tier's object (scalar code)
        std::unique_ptr<IBandedPairWiseSW> ref_obj(ts[0].make(sc.o_del, sc.e_del, sc.o_ins, sc.e_ins, sc.zdrop,
                                                              sc.pen_clip, sc.mat, (int8_t) sc.a, (int8_t) sc.b, 1));
        const std::vector<ExtResult> ref = scalar_reference(*ref_obj, ps);

        for (const Tier &t : ts) {
            for (int compact = 0; compact < 2; ++compact) {
                std::unique_ptr<IBandedPairWiseSW> k(t.make(sc.o_del, sc.e_del, sc.o_ins, sc.e_ins, sc.zdrop,
                                                            sc.pen_clip, sc.mat, (int8_t) sc.a, (int8_t) sc.b, 1));
                if (!compact) k->set_lane_compaction(0, BSW_COMPACT_MIN_W_DEFAULT);   // else the tier's default
                const std::vector<ExtResult> got = score_bsw_batch(*k, pp, 8, kW);
                int diff = 0;
                for (int c = 0; c < kPairs; ++c) diff += !(got[c] == ref[c]);
                CHECK_MESSAGE(diff == 0, "tier " << t.name << ", " << r.what << ", compaction "
                              << (compact ? "default" : "off") << ": " << diff << " of " << kPairs
                              << " pairs differ from the scalar reference");
                // the compaction-on leg must have compacted: superblocks run and lanes moved
                if (compact) {
                    const BswCompactCounts cc = k->lane_compaction_counts();
                    const bool compacted = cc.superblocks > 0 && cc.moves > 0;
                    CHECK_MESSAGE(compacted, "tier " << t.name << ", " << r.what
                                  << ": lane compaction did not run (superblocks " << cc.superblocks
                                  << ", moves " << cc.moves << ")");
                }
                // the instantiation that ran, from the object's per-batch domain counters
                const BswScoreDomainCounts dc = k->bsw8_score_domain_counts();
                if (std::string(t.name) == "avx512bw") {
                    CHECK_MESSAGE((r.fits ? dc.signed_batches : dc.unsigned_batches) > 0,
                                  r.what << ": the " << (r.fits ? "signed" : "unsigned")
                                  << " instantiation did not run (compaction " << (compact ? "default" : "off") << ")");
                    CHECK_MESSAGE((r.fits ? dc.unsigned_batches : dc.signed_batches) == 0,
                                  r.what << ": the " << (r.fits ? "unsigned" : "signed")
                                  << " instantiation ran (compaction " << (compact ? "default" : "off") << ")");
                }
            }
        }
    }
}

#else

TEST_CASE("bsw8: both score-domain instantiations match the scalar reference on every x86 8-bit tier"
          * doctest::test_suite("unit/bandedswa-signed8")) {
    const char *v = std::getenv("BWA_REQUIRE_SIGNED8");
    if (v != nullptr && v[0] != '\0' && std::string(v) != "0")
        FAIL("bandedSWA signed8: BWA_REQUIRE_SIGNED8 is set but this is not an x86 build with the vector kernels");
    MESSAGE("bandedSWA signed8: not an x86 build; the x86 8-bit kernels are not exercised here");
}

#endif
