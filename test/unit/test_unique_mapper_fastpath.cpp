// test/unit/test_unique_mapper_fastpath.cpp
//
// Pin the single-item fast paths added to mem_chain_flt and
// mem_mark_primary_se. Both are algebraically equivalent to the general path
// for the single-element case, but that equivalence rests on a manual trace
// (mem_mark_primary_se_core's inner loop starting at i=1, the
// secondary_all/secondary normalization sequencing, mem_chain_flt's one-range
// [0,1) split collapsing to `kept = 3`). These tests lock the exact
// post-conditions so a future refactor of either function preserves the
// single-element semantics.

#include "doctest/doctest.h"

#include <climits>
#include <cstdint>
#include <cstdlib>

#include "bwamem.h"
#include "utils.h"  // hash_64

// mem_chain_flt is a non-static (external) symbol in bwamem.cpp but is not
// declared in bwamem.h; forward-declare it here with the matching C++
// signature so the linker resolves it against libbwa.a.
int mem_chain_flt(const mem_opt_t *opt, int n_chn_, mem_chain_t *a_, int tid);

namespace {

// Build a single chain carrying one seed. With opt->min_chain_weight == 0 the
// chain always survives the weight filter, so mem_chain_flt reaches its
// n_chn_ == 1 fast path.
mem_chain_t make_single_chain(mem_seed_t *seed_slot, int is_alt) {
    seed_slot->qbeg = 0;
    seed_slot->rbeg = 1000;
    seed_slot->len = 100;
    seed_slot->score = 100;

    mem_chain_t chain{};
    chain.seqid = 0;
    chain.rid = 0;
    chain.n = 1;
    chain.m = 1;  // <= SEEDS_PER_CHAIN, so the drop path never frees seeds
    chain.seeds = seed_slot;
    chain.first = 7;  // sentinel: the filter loop must reset this to -1
    chain.kept = 0;
    chain.is_alt = is_alt ? 1 : 0;
    return chain;
}

}  // namespace

TEST_CASE("mem_chain_flt: single-chain fast path keeps the sole chain") {
    mem_opt_t *opt = mem_opt_init();
    REQUIRE(opt->min_chain_weight == 0);

    for (int is_alt = 0; is_alt <= 1; ++is_alt) {
        mem_seed_t seed{};  // abc() ctor: n_hits = 1, rest zeroed
        mem_chain_t chain = make_single_chain(&seed, is_alt);

        int ret = mem_chain_flt(opt, 1, &chain, /*tid=*/0);

        int kept = chain.kept;  // copy out of the bitfield: doctest binds a ref
        CHECK(ret == 1);            // sole surviving chain is returned
        CHECK(kept == 3);           // marked kept unconditionally
        CHECK(chain.first == -1);   // sentinel cleared by the weight-filter loop
    }

    free(opt);
}

// bwa-mem2 parity for the degenerate "weight filter dropped everything" case.
// The array is never compacted, so slot 0 still holds chain 0, and the legacy
// seqid-range scan built the range {0,1} over that uncompacted slot -- returning
// 1 with the chain marked kept = 3 rather than reporting zero survivors. That is
// bwa-mem2's observable behaviour and this pins it.
//
// bwa-mem2 free()s slot 0's seeds in the filter and then returns the chain
// pointing at the freed block, so its extension reads whatever the allocator
// left behind. bwa-mem3 defers that one free, giving the same chain with its own
// seeds intact. Reachable via -W, -x ont2d and -x pacbio/pbref, the only modes
// that set a nonzero min_chain_weight. Build with ASAN=1 to catch a regression
// here as a use-after-free rather than a silent data change.
TEST_CASE("mem_chain_flt: dropping every chain resurrects slot 0 with live seeds") {
    mem_opt_t *opt = mem_opt_init();
    opt->min_chain_weight = 1000;  // far above the weight of the chains below

    // m > SEEDS_PER_CHAIN puts the seeds on the heap rather than in the caller's
    // seedBuf arena, so the weight filter's free(c->seeds) branch is the one taken.
    const int m = SEEDS_PER_CHAIN + 1;
    mem_chain_t chains[2]{};
    for (int i = 0; i < 2; ++i) {
        mem_seed_t *seeds = (mem_seed_t *)calloc(m, sizeof(mem_seed_t));
        REQUIRE(seeds != NULL);
        seeds[0].qbeg = 0;
        seeds[0].rbeg = 1000 + 10000 * i;
        seeds[0].len = 100;   // weight 100, well under min_chain_weight
        seeds[0].score = 100;

        chains[i].seqid = 0;
        chains[i].rid = 0;
        chains[i].n = 1;
        chains[i].m = m;
        chains[i].seeds = seeds;
        chains[i].first = 7;
        chains[i].kept = 0;
    }
    const mem_seed_t *slot0_seeds = chains[0].seeds;

    CHECK(mem_chain_flt(opt, 2, chains, /*tid=*/0) == 1);

    int kept = chains[0].kept;  // copy out of the bitfield: doctest binds a ref
    CHECK(kept == 3);                          // resurrected slot is marked kept
    CHECK(chains[0].first == -1);              // sentinel cleared by the filter loop
    CHECK(chains[0].seeds == slot0_seeds);     // slot 0's seeds survived the filter
    CHECK(chains[0].seeds[0].len == 100);      // ...and are readable (ASAN gate)

    // Slot 0 comes back live, so its seeds are the caller's to release -- exactly
    // what mem_kernel2_core's chain teardown does. Chain 1 was dropped and freed.
    free(chains[0].seeds);
    free(opt);
}

TEST_CASE("mem_mark_primary_se: single-region fast path normalizes the sole hit") {
    mem_opt_t *opt = mem_opt_init();
    const int64_t id = 12345;

    for (int is_alt = 0; is_alt <= 1; ++is_alt) {
        mem_alnreg_t a{};
        a.is_alt = is_alt ? 1 : 0;
        // Sentinels: the fast path must overwrite sub/alt_sc/secondary/
        // secondary_all/hash and must leave sub_n untouched.
        a.sub = 999;
        a.alt_sc = 999;
        a.secondary = 999;
        a.secondary_all = 999;
        a.hash = 0;
        a.sub_n = 42;

        int ret = mem_mark_primary_se(opt, 1, &a, id);

        CHECK(ret == (is_alt ? 0 : 1));         // n_pri: primary hits only
        CHECK(a.sub == 0);
        CHECK(a.alt_sc == 0);
        CHECK(a.secondary == -1);
        CHECK(a.secondary_all == -1);
        CHECK(a.hash == hash_64((uint64_t)id));  // id, not id + rank
        CHECK(a.sub_n == 42);                    // untouched by the fast path
    }

    free(opt);
}
