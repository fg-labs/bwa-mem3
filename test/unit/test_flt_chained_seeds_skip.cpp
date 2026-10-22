// test/unit/test_flt_chained_seeds_skip.cpp — the seed-SW screen's per-read skip.
//
// mem_flt_chained_seeds rescores each chain's seeds with a local SW, but only for
// reads long enough that the screen's minimum seed length min_l is at most
// MEM_SEEDSW_COEF * l_query. mem_kernel1_core asks mem_flt_chained_seeds_noop once
// per read and does not call mem_flt_chained_seeds when it answers 1. That is
// byte-identical only if:
//
//   1. with --min-ext-len off, mem_flt_chained_seeds_noop is exactly the length
//      test mem_flt_chained_seeds applies per chain (pinned against an independent
//      copy of the expression, including the strict `>` at an exact float tie),
//   2. with --min-ext-len on it answers 0 for every length, because the short-seed
//      filter runs on every chain before the length test, and
//   3. the per-read memo returns what a fresh call would for every read, across
//      changes of read length, and
//   4. on a read it skips, mem_flt_chained_seeds changes nothing, and in
//      particular never touches the reference (bns/pac) or decodes the --meth
//      original read. The fixtures pass NULL bns/pac: a skipped chain that
//      reached the --meth block (which reads bns->l_pac) or the seed SW would
//      crash instead of passing.
//
// Fixtures are built in memory, no data files.

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "doctest/doctest.h"
#include "bwamem.h"

namespace {

// mem_seed_t has a user constructor that sets only its members, so its padding is
// indeterminate: compare it field by field, never with memcmp.
bool seed_eq(const mem_seed_t &a, const mem_seed_t &b)
{
    return a.rbeg == b.rbeg && a.qbeg == b.qbeg && a.len == b.len && a.score == b.score
           && a.aln == b.aln && a.n_hits == b.n_hits;
}

// Independent copy of the screen's length test as it stood in mem_flt_chained_seeds
// (MEM_HSP_COEF 1.1f, MEM_MINSC_COEF 5.5f, MEM_SEEDSW_COEF 0.05f; same float/double
// promotions). Kept separate from the code under test so a drift in either fails.
int reference_skips(int min_chain_weight, int l_query)
{
    double min_l = min_chain_weight ? 1.1f * min_chain_weight : 5.5f * log(l_query);
    return min_l > 0.05f * l_query;
}

struct OptGuard {
    mem_opt_t *opt;
    OptGuard() : opt(mem_opt_init()) {}
    ~OptGuard() { free(opt); }
};

// One read with `n_chains` chains of `seeds_per_chain` seeds each, all stamped
// seqid 0. Seed fields are distinct so a modified seed is detectable.
struct ReadFixture {
    std::string seq;
    std::string orig;
    bseq1_t read;
    std::vector<std::vector<mem_seed_t>> seeds;
    std::vector<mem_chain_t> chains;
    std::vector<std::vector<mem_seed_t>> seeds0;  // copies taken at construction
    std::vector<mem_chain_t> chains0;

    ReadFixture(int l_query, int n_chains, int seeds_per_chain, bool meth)
        : seq(static_cast<size_t>(l_query), '\0'), orig(static_cast<size_t>(l_query), 'T')
    {
        for (int i = 0; i < l_query; ++i) seq[static_cast<size_t>(i)] = static_cast<char>(i & 3);
        memset(&read, 0, sizeof(read));
        read.l_seq = l_query;
        read.seq = &seq[0];
        read.meth_orig_seq = meth ? &orig[0] : nullptr;
        seeds.resize(static_cast<size_t>(n_chains));
        chains.resize(static_cast<size_t>(n_chains));
        for (int c = 0; c < n_chains; ++c) {
            auto &s = seeds[static_cast<size_t>(c)];
            s.resize(static_cast<size_t>(seeds_per_chain));
            for (int j = 0; j < seeds_per_chain; ++j) {
                s[static_cast<size_t>(j)].rbeg = 1000 * c + 10 * j;
                s[static_cast<size_t>(j)].qbeg = j;
                s[static_cast<size_t>(j)].len = 5 + j;
                s[static_cast<size_t>(j)].score = 100 + j;
            }
            mem_chain_t &ch = chains[static_cast<size_t>(c)];
            memset(&ch, 0, sizeof(ch));
            ch.n = ch.m = seeds_per_chain;
            ch.seeds = s.empty() ? nullptr : s.data();
            ch.seqid = 0;
            ch.meth_hypothesis = meth ? 1 : -1;
            ch.w = 40 + c;
            ch.pos = 5000 * c;
        }
        seeds0 = seeds;
        chains0 = chains;
    }

    // Every chain record is byte-for-byte what it was at construction (chains are
    // memset, so their padding compares equal too), and every seed field is unchanged.
    bool untouched() const
    {
        if (memcmp(chains.data(), chains0.data(), chains.size() * sizeof(mem_chain_t)) != 0)
            return false;
        for (size_t c = 0; c < seeds.size(); ++c)
            for (size_t j = 0; j < seeds[c].size(); ++j)
                if (!seed_eq(seeds[c][j], seeds0[c][j])) return false;
        return true;
    }
};

}  // namespace

TEST_CASE("--min-ext-len off: mem_flt_chained_seeds_noop matches the per-chain length test"
          * doctest::test_suite("unit/flt_chained_seeds"))
{
    OptGuard g;
    for (int mcw : {0, 1, 5, 20, 40, 100, 1000}) {
        g.opt->min_chain_weight = mcw;
        for (int l = 0; l <= 25000; ++l) {
            if (mem_flt_chained_seeds_noop(g.opt, l) != reference_skips(mcw, l)) {
                FAIL_CHECK("mismatch at min_chain_weight=" << mcw << " l_query=" << l);
                break;
            }
        }
    }
}

TEST_CASE("default options: short reads skip the screen, long reads do not"
          * doctest::test_suite("unit/flt_chained_seeds"))
{
    OptGuard g;
    REQUIRE(g.opt->min_chain_weight == 0);
    // 5.5 * ln(l) > 0.05 * l holds for 2 <= l <= 724; ln(0) = -inf and ln(1) = 0 fail it.
    CHECK(mem_flt_chained_seeds_noop(g.opt, 0) == 0);
    CHECK(mem_flt_chained_seeds_noop(g.opt, 1) == 0);
    CHECK(mem_flt_chained_seeds_noop(g.opt, 2) == 1);
    CHECK(mem_flt_chained_seeds_noop(g.opt, 100) == 1);
    CHECK(mem_flt_chained_seeds_noop(g.opt, 150) == 1);
    CHECK(mem_flt_chained_seeds_noop(g.opt, 250) == 1);
    CHECK(mem_flt_chained_seeds_noop(g.opt, 724) == 1);
    CHECK(mem_flt_chained_seeds_noop(g.opt, 725) == 0);
    CHECK(mem_flt_chained_seeds_noop(g.opt, 10000) == 0);
}

TEST_CASE("min_chain_weight branch: strict > at an exact float tie"
          * doctest::test_suite("unit/flt_chained_seeds"))
{
    OptGuard g;
    g.opt->min_chain_weight = 5;  // 1.1f * 5 == 0.05f * 110 == 5.5 exactly
    CHECK(mem_flt_chained_seeds_noop(g.opt, 109) == 1);
    CHECK(mem_flt_chained_seeds_noop(g.opt, 110) == 0);
    CHECK(mem_flt_chained_seeds_noop(g.opt, 111) == 0);
}

TEST_CASE("--min-ext-len on: no read length is a no-op"
          * doctest::test_suite("unit/flt_chained_seeds"))
{
    OptGuard g;
    for (int mel : {1, 30}) {
        g.opt->min_ext_len = mel;
        for (int l = 0; l <= 25000; ++l) {
            if (mem_flt_chained_seeds_noop(g.opt, l) != 0) {
                FAIL_CHECK("no-op at min_ext_len=" << mel << " l_query=" << l);
                break;
            }
        }
    }
    g.opt->min_ext_len = -1;  // <= 0 is off, like 0
    CHECK(mem_flt_chained_seeds_noop(g.opt, 150) == 1);
}

TEST_CASE("a skipped read's chains are left untouched without reading the reference"
          * doctest::test_suite("unit/flt_chained_seeds"))
{
    OptGuard g;
    for (int l : {2, 101, 150, 724}) {
        REQUIRE(mem_flt_chained_seeds_noop(g.opt, l) == 1);
        ReadFixture f(l, 4, 3, false);
        mem_flt_chained_seeds(g.opt, nullptr, nullptr, &f.read, static_cast<int>(f.chains.size()),
                              f.chains.data());
        CHECK_MESSAGE(f.untouched(), "l_query=" << l);
    }
    // The min_chain_weight branch, up to the last skipped length before the tie.
    for (int mcw : {5, 40}) {
        g.opt->min_chain_weight = mcw;
        for (int l : {50, 109}) {
            REQUIRE(mem_flt_chained_seeds_noop(g.opt, l) == 1);
            ReadFixture f(l, 2, 3, false);
            mem_flt_chained_seeds(g.opt, nullptr, nullptr, &f.read, static_cast<int>(f.chains.size()),
                                  f.chains.data());
            CHECK_MESSAGE(f.untouched(), "min_chain_weight=" << mcw << " l_query=" << l);
        }
    }
}

TEST_CASE("--meth: a skipped read is not decoded and its chains are untouched"
          * doctest::test_suite("unit/flt_chained_seeds"))
{
    // The --meth block reads seeds[0].rbeg against bns->l_pac; with bns == NULL it
    // would crash, so passing proves the length test runs before it.
    OptGuard g;
    g.opt->meth_mode = 1;
    for (int l : {150, 724}) {
        ReadFixture f(l, 3, 2, true);
        mem_flt_chained_seeds(g.opt, nullptr, nullptr, &f.read, static_cast<int>(f.chains.size()),
                              f.chains.data());
        CHECK_MESSAGE(f.untouched(), "l_query=" << l);
    }
}

TEST_CASE("--min-ext-len still drops short seeds on a read the screen skips"
          * doctest::test_suite("unit/flt_chained_seeds"))
{
    // The caller keeps calling mem_flt_chained_seeds when min_ext_len > 0, because the
    // short-seed filter runs before the length test.
    OptGuard g;
    g.opt->min_ext_len = 6;
    ReadFixture f(150, 2, 3, false);  // seed lengths 5, 6, 7
    mem_flt_chained_seeds(g.opt, nullptr, nullptr, &f.read, static_cast<int>(f.chains.size()),
                          f.chains.data());
    for (size_t c = 0; c < f.chains.size(); ++c) {
        REQUIRE(f.chains[c].n == 2);
        // Survivors are the original seeds 1 and 2, copied whole.
        for (int k = 0; k < 2; ++k)
            CHECK(seed_eq(f.chains[c].seeds[k], f.seeds0[c][static_cast<size_t>(k + 1)]));
    }
}

TEST_CASE("the per-read memo answers like a fresh call across length changes"
          * doctest::test_suite("unit/flt_chained_seeds"))
{
    OptGuard g;
    const int lens[] = {150, 150, 2000, 2000, 150, 724, 725, 724, 1, 2, 0, 150};
    for (int mel : {0, 30}) {
        g.opt->min_ext_len = mel;
        mem_flt_noop_memo_t m = MEM_FLT_NOOP_MEMO_INIT;
        for (int l : lens)
            CHECK_MESSAGE(mem_flt_chained_seeds_noop_memo(&m, g.opt, l)
                              == mem_flt_chained_seeds_noop(g.opt, l),
                          "min_ext_len=" << mel << " l_query=" << l);
    }
}
