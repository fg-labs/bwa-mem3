// Shared synthetic batches and six-field oracles for 16-bit lane compaction tests.
#ifndef BWA_TESTS_BSW_COMPACT16_TEST_H
#define BWA_TESTS_BSW_COMPACT16_TEST_H

#include <algorithm>
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

#if HAVE_BSW_VECTOR_8_16

namespace {

using bwa_tests::ExtPair;
using bwa_tests::ExtResult;
using bwa_tests::ExtScoring;

ExtScoring scoring(int a, int b, int o_del, int e_del, int o_ins, int e_ins, int zdrop, int end_bonus)
{
    ExtScoring sc;
    sc.a = a; sc.b = b; sc.o_del = o_del; sc.e_del = e_del; sc.o_ins = o_ins; sc.e_ins = e_ins;
    sc.zdrop = zdrop; sc.pen_clip = end_bonus;
    sc.fill_mat();
    return sc;
}

/// `sc` with an asymmetric --meth-style matrix: kind 1 frees reference C against read T (the
/// rank-1 path); kind 2 also G against A (the generic matrix path). Kind 0 is symmetric.
ExtScoring with_matrix(ExtScoring sc, int kind)
{
    if (kind >= 1) sc.mat[1 * 5 + 3] = (int8_t)sc.a;
    if (kind >= 2) sc.mat[2 * 5 + 0] = (int8_t)sc.a;
    return sc;
}

using Factory = IBandedPairWiseSW *(*)(int, int, int, int, int, int, const int8_t *, int8_t, int8_t, int);

/// A kernel tier, with its getScores16 lane count and its default groups per superblock (0: it
/// does not compact the 16-bit kernel).
struct Tier {
    std::string name;
    Factory make;
    int lanes;
    int default_groups;
};

#if defined(__aarch64__)
IBandedPairWiseSW *make_neon(int o_del, int e_del, int o_ins, int e_ins, int zdrop, int end_bonus,
                             const int8_t *mat, int8_t a, int8_t b, int nthreads)
{
    return new BandedPairWiseSW(o_del, e_del, o_ins, e_ins, zdrop, end_bonus, mat, a, b, nthreads);
}
#endif

/// The tiers with 16-bit compaction that this host runs. As in test_bandedswa_compact.cpp, an
/// empty list records a WARN, and BWA_REQUIRE_COMPACTION_TIER (1, or a tier name) makes a
/// missing tier a failure.
std::vector<Tier> compacting_tiers()
{
    std::vector<Tier> ts;
#if defined(__aarch64__)
    ts.push_back({"neon", make_neon, 8, BSW16_COMPACT_GROUPS_NEON});
#elif defined(__x86_64__) || defined(__i386__)
    __builtin_cpu_init();
    if (__builtin_cpu_supports("avx2")) ts.push_back({"avx2", make_bsw_kernel_avx2, 16, BSW16_COMPACT_GROUPS_AVX2});
    if (__builtin_cpu_supports("avx512bw"))
        ts.push_back({"avx512bw", make_bsw_kernel_avx512bw, 32, BSW16_COMPACT_GROUPS_AVX512});
#endif
    const char *req = std::getenv("BWA_REQUIRE_COMPACTION_TIER");
    if (req != nullptr && req[0] != '\0' && std::string(req) != "0") {
        const std::string want(req);
        const bool named = want != "1";
        const bool found = named ? std::any_of(ts.begin(), ts.end(), [&](const Tier &t) { return t.name == want; })
                                 : !ts.empty();
        if (!found)
            FAIL("lane compaction 16-bit: BWA_REQUIRE_COMPACTION_TIER=" << want << " but this host has "
                 << (named ? "no such compacting SIMD tier" : "no compacting SIMD tier"));
    } else if (ts.empty()) {
        WARN_MESSAGE(false, "lane compaction 16-bit: this host has no compacting SIMD tier (needs avx2 or "
                            "avx512bw on x86); the equivalence tests are skipped here");
    }
    return ts;
}

std::unique_ptr<IBandedPairWiseSW> make(const Tier &t, const ExtScoring &sc, int groups, int rows)
{
    std::unique_ptr<IBandedPairWiseSW> k(
        t.make(sc.o_del, sc.e_del, sc.o_ins, sc.e_ins, sc.zdrop, sc.pen_clip, sc.mat, (int8_t)sc.a, (int8_t)sc.b, 1));
    k->set_lane_compaction16(groups, rows);
    return k;
}

std::string show(const ExtResult &o)
{
    return std::to_string(o.score) + "/" + std::to_string(o.tle) + "/" + std::to_string(o.qle) + "/" +
           std::to_string(o.gscore) + "/" + std::to_string(o.gtle) + "/" + std::to_string(o.max_off);
}

std::vector<uint8_t> random_bases(std::mt19937_64 &rng, int n)
{
    std::vector<uint8_t> v(n);
    for (auto &c : v) c = (uint8_t)(rng() % 4);
    return v;
}

/// A pair that matches perfectly for `len` bases; the target runs `extra` bases further.
ExtPair perfect(std::mt19937_64 &rng, int len, int extra, int h0)
{
    ExtPair p;
    p.target = random_bases(rng, len + extra);
    p.query.assign(p.target.begin(), p.target.begin() + len);
    p.h0 = h0;
    return p;
}

/// A pair with no useful alignment.
ExtPair divergent(std::mt19937_64 &rng, int len1, int len2, int h0)
{
    ExtPair p;
    p.target = random_bases(rng, len1);
    p.query = random_bases(rng, len2);
    p.h0 = h0;
    return p;
}

/// A pair aligned `d` columns right of the main diagonal (the first `d` query bases inserted).
ExtPair shifted(std::mt19937_64 &rng, int d, int len, int h0)
{
    ExtPair p;
    p.target = random_bases(rng, len + 40);
    p.query = random_bases(rng, d);
    p.query.insert(p.query.end(), p.target.begin(), p.target.begin() + len);
    p.h0 = h0;
    return p;
}

/// Staggered finishes: `groups` * `lanes` pairs, mostly matches of lengths spread over [8, 400]
/// with a few mismatches and 16-bit-range seeds, a quarter divergent or z-dropping.
std::vector<ExtPair> staggered_batch(std::mt19937_64 &rng, int lanes, int groups)
{
    std::vector<ExtPair> v;
    for (int q = 0; q < lanes * groups; q++) {
        const int kind = (int)(rng() % 8);
        if (kind == 0) {
            const int len2 = 10 + (int)(rng() % 120);
            v.push_back(divergent(rng, len2 + (int)(rng() % 100), len2, 10 + (int)(rng() % 60)));
        } else if (kind == 1) {   // a match that z-drops part-way
            ExtPair p = perfect(rng, 40 + (int)(rng() % 200), (int)(rng() % 40), 20 + (int)(rng() % 200));
            for (size_t c = p.query.size() / 2; c < p.query.size(); c++) p.query[c] = (uint8_t)(rng() % 4);
            v.push_back(p);
        } else {
            const int len = rng() % 3 == 0 ? 200 + (int)(rng() % 200) : 8 + (int)(rng() % 120);
            ExtPair p = perfect(rng, len, (int)(rng() % 60), 5 + (int)(rng() % 300));
            for (int m = 0; m < len / 50; m++) p.query[rng() % p.query.size()] = (uint8_t)(rng() % 4);
            v.push_back(p);
        }
    }
    return v;
}

/// Score `batch` with compaction at `groups` / `rows` and with it off (the whole batch, then
/// each pair alone) and against the scalar recurrence; every field must agree. Returns the
/// number of disagreements and adds the compacted run's counters to `moves` and `retired`.
int check_batch(const Tier &t, const ExtScoring &sc, const std::vector<ExtPair> &batch, int w, int groups,
                int rows, const std::string &what, uint64_t &moves, uint64_t &retired)
{
    std::vector<const ExtPair *> ptrs;
    for (const ExtPair &p : batch) ptrs.push_back(&p);
    const std::unique_ptr<IBandedPairWiseSW> cmp = make(t, sc, groups, rows);
    const std::unique_ptr<IBandedPairWiseSW> plain = make(t, sc, 0, 0);
    const std::vector<ExtResult> got = bwa_tests::score_bsw_batch(*cmp, ptrs, 16, w);
    const std::vector<ExtResult> ref = bwa_tests::score_bsw_batch(*plain, ptrs, 16, w);
    const BswCompactCounts c = cmp->lane_compaction16_counts();
    moves += c.moves;
    retired += c.retired;
    int bad = 0;
    for (size_t q = 0; q < batch.size(); q++) {
        const ExtResult alone = bwa_tests::score_bsw_batch(*plain, {ptrs[q]}, 16, w)[0];
        const ExtResult scalar = bwa_tests::run_mem3_scalar(sc, batch[q], w);
        if (got[q] == ref[q] && got[q] == alone && got[q] == scalar) continue;
        if (bad++ < 4)
            MESSAGE(t.name << " " << what << " K=" << groups << " R=" << rows << " w=" << w << " pair " << q << ": compacted "
                         << show(got[q]) << ", plain " << show(ref[q]) << ", alone " << show(alone)
                         << ", scalar " << show(scalar)
                         << " (len1=" << batch[q].target.size() << " len2=" << batch[q].query.size()
                         << " h0=" << batch[q].h0 << "; score/tle/qle/gscore/gtle/max_off)");
    }
    return bad;
}

} // namespace


#endif   // HAVE_BSW_VECTOR_8_16
#endif   // BWA_TESTS_BSW_COMPACT16_TEST_H
