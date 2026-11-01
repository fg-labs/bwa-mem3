// test/unit/test_pair64_sort.cpp -- pdqsort_128 reproduces ks_introsort_128.
//
// mem_pair (src/bwamem_pair.cpp) sorts two pair64_t arrays whose (x, y) keys
// are pairwise distinct by construction: `v` packs the region index and read
// number into y, `u` packs the (k, i) candidate pair into y. On distinct keys
// the sorted order under pair64_lt is unique, so any correct sort produces the
// same array, and pdqsort_128 can stand in for the klib introsort. This pins
// that claim against the introsort itself, on arrays shaped like mem_pair's
// (many equal x -- same position, different regions -- and distinct y), and
// asserts the distinctness the argument rests on. (Duplicate pair64_t keys are
// byte-identical, so a permutation difference among them could not be observed
// anyway; the uniqueness argument, not the algorithm, carries the parity.) A
// second case pins that for duplicate-heavy input too, which `aln`'s pairing
// (src/aln/bwape.cpp) sorts with pdqsort_128: its rows repeat when a read has
// the same position and tag twice, e.g. an unplaced (uint64_t)-1 position.
//
// Fixtures are generated in-memory; no test data files are read.
#include <cstdint>
#include <cstring>
#include <vector>
#include "doctest/doctest.h"
#include "utils.h"

namespace {

struct Rng {
    uint64_t s;
    explicit Rng(uint64_t seed) : s(seed) {}
    uint64_t next() { s = s * 6364136223846793005ULL + 1442695040888963407ULL; return s >> 11; }
    int in(int lo, int hi) { return lo + static_cast<int>(next() % static_cast<uint64_t>(hi - lo + 1)); }
};

// n keys in mem_pair's `v` shape: x = position drawn from `x_pool` slots (so
// equal x is common), y = score << 32 | i << 2 | strand << 1 | r with the index
// i distinct per element, so (x, y) is always distinct.
std::vector<pair64_t> mem_pair_like(Rng &rng, int n, int x_pool) {
    std::vector<pair64_t> v;
    for (int i = 0; i < n; ++i) {
        pair64_t k;
        k.x = static_cast<uint64_t>(rng.in(0, x_pool - 1)) * 1000;
        k.y = static_cast<uint64_t>(rng.in(20, 300)) << 32 |
              static_cast<uint64_t>(i) << 2 | static_cast<uint64_t>(rng.in(0, 1)) << 1 |
              static_cast<uint64_t>(rng.in(0, 1));
        v.push_back(k);
    }
    return v;
}

// n keys in mem_pair's `u` shape: y = k << 32 | i over distinct (k, i) pairs,
// x = q << 32 | 32-bit hash with q drawn from few values so x ties are common.
std::vector<pair64_t> mem_pair_u_like(Rng &rng, int n, int q_pool) {
    std::vector<pair64_t> v;
    for (int t = 0; t < n; ++t) {
        pair64_t u;
        // Distinct (k, i): walk a grid so no pair repeats.
        const uint64_t k = static_cast<uint64_t>(t / 7), i = static_cast<uint64_t>(t % 7);
        u.y = k << 32 | i;
        u.x = static_cast<uint64_t>(rng.in(0, q_pool - 1)) << 32 | (rng.next() & 0xffffffffULL);
        v.push_back(u);
    }
    return v;
}

bool same(const std::vector<pair64_t> &a, const std::vector<pair64_t> &b) {
    return a.size() == b.size() &&
           (a.empty() || memcmp(a.data(), b.data(), a.size() * sizeof(pair64_t)) == 0);
}

}  // namespace

TEST_CASE("pdqsort_128 equals ks_introsort_128 on distinct (x, y) keys"
          * doctest::test_suite("unit/pair64_sort")) {
    // Each variant is its own SUBCASE so doctest reports and --subcase-filters
    // them independently; doctest re-enters the body per subcase, so each holds
    // its own Rng with a distinct seed.
    SUBCASE("v-shape keys: heavy x ties (realistic case), y always distinct") {
        Rng rng(0x5eed1234ULL);
        for (int t = 0; t < 300; ++t) {
            const int n = 1 + t;
            std::vector<pair64_t> in = mem_pair_like(rng, n, 1 + t / 10);
            std::vector<pair64_t> p(in), k(in);
            pdqsort_128(p.size(), p.data());
            ks_introsort_128(k.size(), k.data());
            CHECK(same(p, k));
            // Sorted, and every adjacent pair strictly ordered (distinct keys).
            for (size_t i = 1; i < k.size(); ++i)
                CHECK((k[i - 1].x < k[i].x || (k[i - 1].x == k[i].x && k[i - 1].y < k[i].y)));
        }
    }
    // The `u` shape (candidate pairs): equal pairing scores are common, the
    // low-32 hash rarely collides, and (k, i) is distinct.
    SUBCASE("u-shape candidate pairs: common equal scores, rare hash collisions") {
        Rng rng(0x5eed5678ULL);
        for (int t = 0; t < 200; ++t) {
            const int n = 1 + t;
            std::vector<pair64_t> in = mem_pair_u_like(rng, n, 1 + t / 20);
            std::vector<pair64_t> p(in), k(in);
            pdqsort_128(p.size(), p.data());
            ks_introsort_128(k.size(), k.data());
            CHECK(same(p, k));
            for (size_t i = 1; i < k.size(); ++i)
                CHECK((k[i - 1].x < k[i].x || (k[i - 1].x == k[i].x && k[i - 1].y < k[i].y)));
        }
    }
}

// Duplicates are allowed: any input sorts to the one array ks_introsort_128
// gives, because rows equal under (x, y) are the same bytes. Small x/y pools make
// duplicates common; a pool of 1 makes every row equal; the extreme words include
// (uint64_t)-1, `aln`'s unplaced position.
TEST_CASE("pdqsort_128 equals ks_introsort_128 on duplicate-heavy keys"
          * doctest::test_suite("unit/pair64_sort")) {
    Rng rng(0xd0b1e5ULL);
    auto check = [](std::vector<pair64_t> in) {
        std::vector<pair64_t> p(in), k(in);
        p.reserve(1);  // non-null data() for n == 0
        k.reserve(1);
        pdqsort_128(p.size(), p.data());
        ks_introsort_128(k.size(), k.data());
        return same(p, k);
    };
    auto pooled = [&rng](int n, int x_pool, int y_pool) {
        std::vector<pair64_t> v(static_cast<size_t>(n));
        for (auto &r : v) {
            r.x = static_cast<uint64_t>(rng.in(0, x_pool - 1)) * 977;
            r.y = static_cast<uint64_t>(rng.in(0, y_pool - 1));
        }
        return v;
    };
    for (int n = 0; n <= 300; ++n) {
        CHECK(check(pooled(n, 3, 4)));
        CHECK(check(pooled(n, 1 + n / 8, 2)));
        CHECK(check(pooled(n, 1, 1 + n)));
    }
    for (int n = 301; n <= 5000; n += 53) {
        CHECK(check(pooled(n, 5, 8)));
        CHECK(check(pooled(n, 1000, 3)));
    }
    for (int n : {0, 1, 2, 3, 23, 24, 25, 100, 5000}) CHECK(check(pooled(n, 1, 1)));
    const uint64_t vals[] = {0, 1, 2, 0x7fffffffffffffffULL, 0xfffffffffffffffeULL, 0xffffffffffffffffULL};
    for (int n = 0; n <= 200; ++n) {
        std::vector<pair64_t> v(static_cast<size_t>(n));
        for (auto &r : v) { r.x = vals[rng.in(0, 5)]; r.y = vals[rng.in(0, 5)]; }
        CHECK(check(v));
    }
}

// mem_pair reaches the pairing sort with zero candidates when a thread's first
// call (or a released scratch buffer) leaves the buffer null and n == 0. The
// call site guards `if (v.n > 1)` so PDQSORT_INIT never forms `a + n` on a null
// pointer; this pins the empty/single-element sort as a no-op that leaves the
// (possibly one) element in place. Buffers are reserved so data() is non-null.
TEST_CASE("pdqsort_128 is a no-op on empty and single-element inputs"
          * doctest::test_suite("unit/pair64_sort")) {
    std::vector<pair64_t> empty;
    empty.reserve(1);  // non-null data() for the n == 0 no-op
    pdqsort_128(empty.size(), empty.data());
    CHECK(empty.empty());

    std::vector<pair64_t> one(1);
    one[0].x = 42;
    one[0].y = 7;
    pdqsort_128(one.size(), one.data());
    CHECK(one[0].x == 42u);
    CHECK(one[0].y == 7u);
}

TEST_CASE("pair64_lt_bf agrees with pair64_lt"
          * doctest::test_suite("unit/pair64_sort")) {
    // pdqsort_128's branch-free comparator must be the same order as the
    // short-circuit one the strict-order check uses. Tiny ranges plus extremes
    // take every equal-x branch.
    const uint64_t vals[] = {0, 1, 2, 0x7fffffffffffffffULL, 0xffffffffffffffffULL};
    int n = 0;
    for (const uint64_t ax : vals) for (const uint64_t ay : vals)
        for (const uint64_t bx : vals) for (const uint64_t by : vals) {
            pair64_t a, b;
            a.x = ax; a.y = ay; b.x = bx; b.y = by;
            CHECK(!!pair64_lt(a, b) == !!pair64_lt_bf(a, b));
            ++n;
        }
    CHECK(n == 625);
}
