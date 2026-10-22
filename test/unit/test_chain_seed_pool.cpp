// Unit tests for the chaining seed pool (src/chain_seed_pool.h). While a read
// is chained, each chain's seeds are a list linked by index in a per-thread
// pool (ChainSeedPool::start / append), and chain_pool_drain copies a finished
// chain's seeds, contiguous and in merge order, into the per-thread seed window,
// or into a heap array of its own (m == SEEDS_PER_CHAIN + 1) once the window is
// full. The checks compare against an independent model that keeps each chain's
// seeds in a std::vector of its own.

#include "doctest/doctest.h"
#include "chain_seed_pool.h"

#include <algorithm>
#include <random>
#include <vector>

namespace {

mem_seed_t make_seed(int64_t id)
{
    mem_seed_t s;
    s.rbeg = 1000 * id;
    s.qbeg = static_cast<int32_t>(id);
    s.len = static_cast<int32_t>(19 + id % 40);
    s.score = static_cast<int32_t>(id % 97);
    s.aln = static_cast<int>(id % 3);
    s.n_hits = static_cast<int32_t>(1 + id % 5);
    return s;
}

// mem_seed_t has a user constructor that sets only its members, so its padding
// is indeterminate: compare field by field.
bool seed_eq(const mem_seed_t &a, const mem_seed_t &b)
{
    return a.rbeg == b.rbeg && a.qbeg == b.qbeg && a.len == b.len && a.score == b.score
           && a.aln == b.aln && a.n_hits == b.n_hits;
}

// One read: n_chains chains built by n_seeds random start/append calls through
// the pool and through the owned-array model, then drained in a random order into
// a window of window_size slots starting at window_used. Checks every drained
// chain against the model and the window/heap accounting.
void check_read(std::mt19937 &rng, ChainSeedPool &pool, int n_chains, int n_seeds,
                int64_t window_size, int64_t window_used)
{
    pool.clear();
    std::vector<mem_chain_t> chains(static_cast<size_t>(n_chains));
    std::vector<std::vector<mem_seed_t>> model(static_cast<size_t>(n_chains));
    int64_t id = 0;
    for (int c = 0; c < n_chains; ++c) {
        const mem_seed_t s = make_seed(id++);
        chains[static_cast<size_t>(c)].pos = s.rbeg;
        chains[static_cast<size_t>(c)].rid = c;
        pool.start(&chains[static_cast<size_t>(c)], &s);
        model[static_cast<size_t>(c)].push_back(s);
    }
    std::uniform_int_distribution<int> pick(0, n_chains - 1);
    for (int k = 0; k < n_seeds; ++k) {
        const size_t c = static_cast<size_t>(pick(rng));
        const mem_seed_t s = make_seed(id++);
        pool.append(&chains[c], &s);
        model[c].push_back(s);
    }

    std::vector<int> order(static_cast<size_t>(n_chains));
    for (int c = 0; c < n_chains; ++c) order[static_cast<size_t>(c)] = c;
    std::shuffle(order.begin(), order.end(), rng);

    // Every window slot starts as a sentinel, so a drain that writes past its
    // chain's n slots is seen at the slot just after them.
    mem_seed_t sentinel = make_seed(-1);
    sentinel.rbeg = -12345;
    std::vector<mem_seed_t> window(static_cast<size_t>(window_size), sentinel);
    int64_t count = window_used;
    for (int c : order) {
        const auto &want = model[static_cast<size_t>(c)];
        const int64_t before = count;
        const bool fits = before + static_cast<int64_t>(want.size()) <= window_size;
        const mem_chain_t out = chain_pool_drain(chains[static_cast<size_t>(c)], pool,
                                                 window.data(), &count, window_size);
        CAPTURE(c);
        CAPTURE(fits);
        REQUIRE(out.n == static_cast<int32_t>(want.size()));
        CHECK(out.pos == chains[static_cast<size_t>(c)].pos);
        CHECK(out.rid == c);
        for (size_t k = 0; k < want.size(); ++k) CHECK(seed_eq(out.seeds[k], want[k]));
        if (fits) {
            CHECK(out.m == SEEDS_PER_CHAIN);
            CHECK(out.seeds == window.data() + before);
            CHECK(count == before + out.n);
            if (count < window_size) CHECK(seed_eq(window[static_cast<size_t>(count)], sentinel));
        } else {
            CHECK(out.m == SEEDS_PER_CHAIN + 1);
            CHECK(count == before);
            free(out.seeds);
        }
    }
}

}  // namespace

TEST_CASE("chain seed pool: drained chains hold their seeds in merge order in the window"
          * doctest::test_suite("unit/chain"))
{
    std::mt19937 rng(7);
    ChainSeedPool pool;
    for (int rep = 0; rep < 200; ++rep) {
        const int n_chains = 1 + static_cast<int>(rng() % 40);
        const int n_seeds = static_cast<int>(rng() % 60);
        CAPTURE(rep);
        // The window always has room: every chain lands in it.
        check_read(rng, pool, n_chains, n_seeds, n_chains + n_seeds + 8, 3);
    }
}

TEST_CASE("chain seed pool: a full window sends a chain to its own heap array"
          * doctest::test_suite("unit/chain"))
{
    std::mt19937 rng(11);
    ChainSeedPool pool;
    for (int rep = 0; rep < 200; ++rep) {
        const int n_chains = 1 + static_cast<int>(rng() % 20);
        const int n_seeds = static_cast<int>(rng() % 40);
        const int64_t size = static_cast<int64_t>(rng() % 24);
        const int64_t used = size > 0 ? static_cast<int64_t>(rng() % static_cast<uint32_t>(size + 1)) : 0;
        CAPTURE(rep);
        CAPTURE(size);
        CAPTURE(used);
        // Some chains fit and some do not, in any order; a smaller chain may still
        // fit after a larger one did not.
        check_read(rng, pool, n_chains, n_seeds, size, used);
    }
}

TEST_CASE("chain seed pool: clear() starts the next read from an empty pool"
          * doctest::test_suite("unit/chain"))
{
    ChainSeedPool pool;
    mem_chain_t a{}, b{};
    const mem_seed_t s0 = make_seed(0), s1 = make_seed(1), s2 = make_seed(2);
    pool.start(&a, &s0);
    pool.append(&a, &s1);
    CHECK(pool.seed.size() == 2);
    pool.clear();
    CHECK(pool.seed.empty());
    CHECK(pool.next.empty());
    pool.start(&b, &s2);
    CHECK(b.n == 1);
    CHECK(b.link.head == 0);
    CHECK(b.link.tail == 0);
    mem_seed_t window[1];
    int64_t count = 0;
    const mem_chain_t out = chain_pool_drain(b, pool, window, &count, 1);
    CHECK(out.m == SEEDS_PER_CHAIN);
    CHECK(seed_eq(out.seeds[0], s2));
}

TEST_CASE("chain seed pool: clear() releases a pool grown past KEEP_SEEDS"
          * doctest::test_suite("unit/chain"))
{
    ChainSeedPool pool;
    mem_chain_t c{};
    const mem_seed_t s = make_seed(3);
    pool.start(&c, &s);
    pool.clear();
    const size_t small_cap = pool.seed.capacity();
    CHECK(small_cap > 0);   // an ordinary read keeps its capacity
    pool.seed.reserve(ChainSeedPool::KEEP_SEEDS + 1);
    pool.next.reserve(ChainSeedPool::KEEP_SEEDS + 1);
    pool.clear();
    CHECK(pool.seed.capacity() == 0);
    CHECK(pool.next.capacity() == 0);
    pool.start(&c, &s);   // and is usable again
    CHECK(pool.seed.size() == 1);
    CHECK(c.link.head == 0);
}
