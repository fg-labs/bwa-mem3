// test/unit/test_chain_flt_survivors.cpp
//
// On x86-64, mem_chain_flt's pairwise pass reads each survivor through a
// contiguous record (bounds, weight, the float product w * drop_ratio, is_alt),
// and the shadow promotion and extend-cap count run on the pass's first[]/kept[]
// scratch with the write-back fused into the compaction; other targets compile
// the original form. Both must meet the same contract, so these tests run on
// every target and check whichever form it compiles. The contract is the naive filter:
// the O(n^2) overlap scan over the weight-sorted chains, shadow promotion of
// each survivor's first shadower, the extend cap, and compaction of the kept
// chains in weight order with every other field unchanged. This pins
// mem_chain_flt against an independent transcription of that filter over
// generated reads and a spread of filter options.
//
// The oracle sorts with the same klib introsort and comparator as the filter,
// so chains of equal weight come out in the same order and tied reads are
// covered. Some chains keep their seeds on the heap (m > SEEDS_PER_CHAIN), so
// the filter frees the dropped ones; build with ASAN=1 to check that path (the
// rebuild-variants CI job runs this suite under ASan, where LeakSanitizer
// fails it on a missed free).

#include "doctest/doctest.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <random>
#include <vector>

#include "bwamem.h"
#include "ksort.h"
#include "macro.h"

// mem_chain_flt is non-static in bwamem.cpp but not declared in bwamem.h.
int mem_chain_flt(const mem_opt_t *opt, int n_chn_, mem_chain_t *a_, int tid);

namespace {

#define oracle_flt_lt(a, b) ((a).w > (b).w)  // the filter's mem_flt comparator
KSORT_INIT(oracle_flt, mem_chain_t, oracle_flt_lt)

struct Paths {
    long dominated = 0;   // chains the overlap scan drops
    long phase2 = 0;      // first[] set by a survivor's later overlaps (after its first)
};

int beg(const mem_chain_t &c) { return c.seeds->qbeg; }
int end(const mem_chain_t &c) { return c.seeds[c.n - 1].qbeg + c.seeds[c.n - 1].len; }

// The naive filter, past the weight filter, for the default target: the
// surviving chains, in order, with first/kept set and every other field as given.
std::vector<mem_chain_t> oracle(const mem_opt_t *opt, const std::vector<mem_chain_t> &in, Paths &paths)
{
    std::vector<mem_chain_t> f;
    for (auto c : in) {
        c.first = -1; c.kept = 0;
        c.w = mem_chain_weight(&c);
        if (opt->min_chain_weight <= 0 || (int)c.w >= opt->min_chain_weight) f.push_back(c);
    }
    if (f.empty()) return f;
    if (f.size() == 1) { f[0].kept = 3; return f; }
    ks_introsort(oracle_flt, f.size(), f.data());
    const int n = (int)f.size();

    const float mask_level = opt->mask_level, drop_ratio = opt->drop_ratio;
    const int max_chain_gap = opt->max_chain_gap, min_seed_gap = opt->min_seed_len << 1;
    std::vector<int> first(n, -1), kept(n, 0), surv;
    kept[0] = 3; surv.push_back(0);
    for (int i = 1; i < n; ++i) {
        int large_ovlp = 0;
        const int cbi = beg(f[i]), cei = end(f[i]), li = cei - cbi, wi = f[i].w;
        const int ai_alt = f[i].is_alt;
        size_t k;
        for (k = 0; k < surv.size(); ++k) {
            const int j = surv[k];
            const int b_max = std::max(beg(f[j]), cbi), e_min = std::min(end(f[j]), cei);
            if (e_min > b_max && (!f[j].is_alt || ai_alt)) {
                const int lj = end(f[j]) - beg(f[j]);
                const int min_l = std::min(li, lj);
                if (e_min - b_max >= min_l * mask_level && min_l < max_chain_gap) {
                    if (large_ovlp && first[j] < 0) ++paths.phase2;
                    large_ovlp = 1;
                    if (first[j] < 0) first[j] = i;
                    if (wi < f[j].w * drop_ratio && (int)f[j].w - wi >= min_seed_gap) break;
                }
            }
        }
        if (k == surv.size()) { surv.push_back(i); kept[i] = large_ovlp ? 2 : 3; }
        else ++paths.dominated;
    }

    // Promote every survivor's first shadower, count the extend cap over the
    // kept=1/2 chains, demote those past it and keep the rest in order.
    for (int s : surv) if (first[s] >= 0) kept[first[s]] = 1;
    int i, k;
    for (i = k = 0; i < n; ++i) {
        if (kept[i] == 0 || kept[i] == 3) continue;
        if (++k >= opt->max_chain_extend) break;
    }
    const int cap_cutoff = i;
    std::vector<mem_chain_t> out;
    for (i = 0; i < n; ++i) {
        int kp = kept[i];
        if (i >= cap_cutoff && kp < 3) kp = 0;
        if (kp == 0) continue;
        f[i].first = first[i]; f[i].kept = kp;
        out.push_back(f[i]);
    }
    return out;
}

// The oracle's chains share their seed arrays with the filter's input, so the
// seed pointers alone cannot show a retained seed rewritten in place: copy the
// seeds out before mem_chain_flt runs and compare their contents after.
std::vector<std::vector<mem_seed_t>> snapshot_seeds(const std::vector<mem_chain_t> &chains)
{
    std::vector<std::vector<mem_seed_t>> out;
    for (const auto &c : chains) out.emplace_back(c.seeds, c.seeds + c.n);
    return out;
}

void check_same_chain(const mem_chain_t &got, const mem_chain_t &want, const std::vector<mem_seed_t> &want_seeds)
{
    CHECK(got.seeds == want.seeds);
    CHECK(got.n == want.n);
    REQUIRE(got.n == (int)want_seeds.size());
    for (int j = 0; j < got.n; ++j) {
        CAPTURE(j);
        CHECK(got.seeds[j].rbeg == want_seeds[j].rbeg);
        CHECK(got.seeds[j].qbeg == want_seeds[j].qbeg);
        CHECK(got.seeds[j].len == want_seeds[j].len);
        CHECK(got.seeds[j].score == want_seeds[j].score);
        CHECK(got.seeds[j].aln == want_seeds[j].aln);
        CHECK(got.seeds[j].n_hits == want_seeds[j].n_hits);
    }
    CHECK(got.m == want.m);
    CHECK(got.first == want.first);
    CHECK(got.rid == want.rid);
    CHECK(got.seqid == want.seqid);
    CHECK(got.cseed == want.cseed);
    CHECK((int)got.w == (int)want.w);
    CHECK((int)got.kept == (int)want.kept);
    CHECK((int)got.is_alt == (int)want.is_alt);
    CHECK((int)got.meth_hypothesis == (int)want.meth_hypothesis);
    CHECK(got.frac_rep == want.frac_rep);
    CHECK(got.pos == want.pos);
}

// A read's chains: 1-3 seeds each; `heap` chains own a calloc'd seed array
// (m > SEEDS_PER_CHAIN), the rest point into `pool`. Seed lengths come from a
// narrow range, so equal weights are common.
std::vector<mem_chain_t> make_read(std::mt19937 &rng, int n, std::vector<mem_seed_t> &pool)
{
    std::uniform_int_distribution<int> qb(0, 300), len(5, 60), ns(1, 3), pct(0, 99), rid(0, 3);
    std::vector<mem_chain_t> a;
    pool.assign((size_t)n * 3, mem_seed_t());
    for (int c = 0; c < n; ++c) {
        const int nseed = ns(rng);
        const bool heap = pct(rng) < 30;
        mem_seed_t *s = heap ? (mem_seed_t *)calloc(SEEDS_PER_CHAIN + 3, sizeof(mem_seed_t))
                             : &pool[(size_t)c * 3];
        REQUIRE(s != nullptr);
        int q = qb(rng);
        for (int j = 0; j < nseed; ++j) {
            s[j].qbeg = q;
            s[j].len = len(rng);
            s[j].rbeg = 100000 + q;
            s[j].score = s[j].len;
            q += s[j].len / 2 + 1;
        }
        mem_chain_t ch{};
        ch.n = nseed;
        ch.m = heap ? SEEDS_PER_CHAIN + 3 : 1;
        ch.seeds = s;
        ch.rid = rid(rng);
        ch.cseed = pct(rng);
        ch.is_alt = pct(rng) < 20;
        ch.meth_hypothesis = pct(rng) < 10 ? 1 : -1;
        ch.frac_rep = (float)pct(rng) / 100.0f;
        ch.pos = 5000 * c + pct(rng);
        a.push_back(ch);
    }
    return a;
}

void free_survivors(mem_chain_t *a, int k)
{
    for (int t = 0; t < k; ++t)
        if (a[t].m > SEEDS_PER_CHAIN) free(a[t].seeds);
}

}  // namespace

TEST_CASE("mem_chain_flt matches the naive overlap filter, shadow promotion and extend cap"
          * doctest::test_suite("unit/chain_flt"))
{
    mem_opt_t *opt = mem_opt_init();
    // Each option axis cycles at its own period, so all 576 combinations occur.
    const float drop_ratios[] = {0.5f, 0.3f, 0.77f, 0.9f};
    const float mask_levels[] = {0.5f, 0.2f, 0.8f};
    const int extends[] = {1, 2, 3, 1 << 30};
    const int seed_lens[] = {19, 3, 40};
    const int gaps[] = {10000, 60};
    const int min_ws[] = {0, 30};

    std::mt19937 rng(20261022);
    std::uniform_int_distribution<int> nchn(2, 40);
    std::vector<mem_seed_t> pool;
    Paths paths;
    long promoted = 0, overlapping = 0, ties = 0;
    for (int round = 0; round < 1728; ++round) {
        opt->drop_ratio = drop_ratios[round % 4];
        opt->mask_level = mask_levels[(round / 4) % 3];
        opt->max_chain_extend = extends[(round / 12) % 4];
        opt->min_seed_len = seed_lens[(round / 48) % 3];
        opt->max_chain_gap = gaps[(round / 144) % 2];
        opt->min_chain_weight = min_ws[(round / 288) % 2];
        CAPTURE(round);

        std::vector<mem_chain_t> a = make_read(rng, nchn(rng), pool);
        const std::vector<mem_chain_t> want = oracle(opt, a, paths);
        const auto want_seeds = snapshot_seeds(want);
        const int k = mem_chain_flt(opt, (int)a.size(), a.data(), /*tid=*/0);
        REQUIRE(k == (int)want.size());
        for (int t = 0; t < k; ++t) {
            CAPTURE(t);
            check_same_chain(a[t], want[t], want_seeds[t]);
            if (want[t].kept == 1) ++promoted;
            if (want[t].kept == 2) ++overlapping;
            if (t > 0 && want[t].w == want[t - 1].w) ++ties;
        }
        free_survivors(a.data(), k);
    }
    // The generated reads reach every path the filter's result depends on.
    CHECK(paths.dominated > 0);
    CHECK(paths.phase2 > 0);
    CHECK(promoted > 0);
    CHECK(overlapping > 0);
    CHECK(ties > 0);
    free(opt);
}

TEST_CASE("mem_chain_flt domination uses the float product w * drop_ratio"
          * doctest::test_suite("unit/chain_flt"))
{
    // Weights 10, 3, 2 on query spans [0,10), [8,11), [10,12) at drop_ratio
    // 0.3f. Chain 1 overlaps chain 0, and the domination test is 3 < 10 * 0.3f:
    // the float product is exactly 3 (10 * 0.3f lies half an ulp above 3 and
    // rounds to even), so chain 1 survives. Chain 2 then shadows chain 1, and
    // promotion gives kept {3, 1, 1}, first {1, 2, -1}. Computed in double the
    // product is above 3, chain 1 is dominated, chain 2 keeps 3 and chain 1 is
    // only resurrected by promotion with first = -1: a different result.
    mem_opt_t *opt = mem_opt_init();
    opt->drop_ratio = 0.3f;
    opt->min_seed_len = 3;   // min_seed_gap 6, so the weight gap 10 - 3 qualifies
    REQUIRE(10 * opt->drop_ratio == 3.0f);
    REQUIRE(10 * (double)opt->drop_ratio > 3.0);
    mem_seed_t s[3];
    const int qbeg[3] = {0, 8, 10}, len[3] = {10, 3, 2};
    std::vector<mem_chain_t> a(3);
    for (int c = 0; c < 3; ++c) {
        s[c].qbeg = qbeg[c]; s[c].len = len[c];
        s[c].rbeg = 1000 + qbeg[c]; s[c].score = len[c];
        a[c] = mem_chain_t{};
        a[c].n = 1; a[c].m = 1; a[c].seeds = &s[c]; a[c].meth_hypothesis = -1;
    }
    Paths paths;
    const std::vector<mem_chain_t> want = oracle(opt, a, paths);
    const auto want_seeds = snapshot_seeds(want);
    CHECK(paths.dominated == 0);
    const int k = mem_chain_flt(opt, 3, a.data(), /*tid=*/0);
    REQUIRE(k == 3);
    REQUIRE(want.size() == 3);
    const int kept[3] = {3, 1, 1}, first[3] = {1, 2, -1};
    for (int t = 0; t < 3; ++t) {
        CAPTURE(t);
        CHECK(a[t].seeds == &s[t]);
        CHECK((int)a[t].kept == kept[t]);
        CHECK(a[t].first == first[t]);
        check_same_chain(a[t], want[t], want_seeds[t]);
    }
    free(opt);
}
