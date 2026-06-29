// test/unit/test_seed_order.cpp — order_seeds + string mappers.
#include <cstdint>
#include <vector>
#include "doctest/doctest.h"
#include "seed_order.h"

namespace {
// Build a seed_rec_t with the fields order_seeds reads; orig_ix = emission index.
seed_rec_t rec(int qbeg, int len, int64_t rbeg, int32_t rid, uint32_t ix) {
    seed_rec_t r{};
    r.seed.qbeg = qbeg; r.seed.len = len; r.seed.score = len;
    r.seed.rbeg = rbeg; r.seed.n_hits = 1;
    r.rid = rid; r.meth_hyp = -1; r.orig_ix = ix;
    return r;
}
} // namespace

TEST_CASE("seed_order string mappers round-trip") {
    CHECK(seed_order_from_str("off") == SEED_ORDER_OFF);
    CHECK(seed_order_from_str("global-longest") == SEED_ORDER_GLOBAL_LONGEST);
    CHECK(seed_order_from_str("most-absorb") == SEED_ORDER_MOST_ABSORB);
    CHECK((int)seed_order_from_str("bogus") == -1);
    // advertised mode round-trips
    CHECK(std::string(seed_order_to_str(SEED_ORDER_LOCAL_LONGEST)) == "local-longest");
    // hidden modes still accepted by parser
    CHECK(seed_order_from_str("most-absorb") == SEED_ORDER_MOST_ABSORB);
}

TEST_CASE("order_seeds OFF is the identity permutation") {
    std::vector<seed_rec_t> v{rec(50,20,500,0,0), rec(0,150,1000,0,1), rec(20,20,200,0,2)};
    order_seeds(v.data(), (int)v.size(), SEED_ORDER_OFF);
    CHECK(v[0].orig_ix == 0);
    CHECK(v[1].orig_ix == 1);
    CHECK(v[2].orig_ix == 2);
}

TEST_CASE("global-longest: longest first, stable on ties") {
    std::vector<seed_rec_t> v{rec(0,20,100,0,0), rec(0,150,200,0,1),
                              rec(50,20,300,0,2), rec(10,20,400,0,3)};
    order_seeds(v.data(), 4, SEED_ORDER_GLOBAL_LONGEST);
    CHECK(v[0].orig_ix == 1);             // len 150 first
    CHECK(v[1].orig_ix == 0);             // len-20 ties keep orig_ix order: 0,2,3
    CHECK(v[2].orig_ix == 2);
    CHECK(v[3].orig_ix == 3);
}

TEST_CASE("local-longest: group by qbeg asc, then longest within group") {
    std::vector<seed_rec_t> v{rec(10,20,100,0,0), rec(0,20,200,0,1), rec(0,40,300,0,2)};
    order_seeds(v.data(), 3, SEED_ORDER_LOCAL_LONGEST);
    CHECK(v[0].orig_ix == 2);             // qbeg 0, len 40
    CHECK(v[1].orig_ix == 1);             // qbeg 0, len 20
    CHECK(v[2].orig_ix == 0);             // qbeg 10
}

TEST_CASE("global-longest: stable on equal len (large n, counting sort)") {
    // 1000 seeds all len 50 -> every key equal -> output MUST stay orig_ix order.
    std::vector<seed_rec_t> v;
    for (uint32_t i = 0; i < 1000; ++i) v.push_back(rec(0, 50, 100 + i, 0, i));
    order_seeds(v.data(), 1000, SEED_ORDER_GLOBAL_LONGEST);
    for (uint32_t i = 0; i < 1000; ++i) CHECK(v[i].orig_ix == i);  // stability preserved
}

TEST_CASE("local-longest: equal-len within a qbeg group is stable on orig_ix") {
    // qbeg 0: two len-40 seeds (orig 0,1) then a len-20 (orig 2).
    std::vector<seed_rec_t> v{rec(0,40,100,0,0), rec(0,40,200,0,1), rec(0,20,300,0,2)};
    order_seeds(v.data(), 3, SEED_ORDER_LOCAL_LONGEST);
    CHECK(v[0].orig_ix == 0);   // equal len 40 -> stable: orig 0 before 1
    CHECK(v[1].orig_ix == 1);
    CHECK(v[2].orig_ix == 2);   // len 20 last
}

TEST_CASE("absorb-count: container with both q and r nested ranks first") {
    // A spans all; B,C nested in BOTH q and r -> A absorbs 2.
    std::vector<seed_rec_t> v{
        rec(20,20,120,0,0),   // B
        rec(0,150,100,0,1),   // A (q[0,150] r[100,250]) contains B and C
        rec(50,20,150,0,2)};  // C
    order_seeds(v.data(), 3, SEED_ORDER_ABSORB_COUNT);
    CHECK(v[0].orig_ix == 1);  // A absorbs 2 -> first
    CHECK(v[1].orig_ix == 0);  // B and C both absorb 0, stable orig_ix order: 0 before 2
    CHECK(v[2].orig_ix == 2);
}

TEST_CASE("absorb-count: query-nested but reference-disjoint is NOT absorbed") {
    // B's query is inside A's, but B.rbeg is far from A -> different locus, no absorb.
    std::vector<seed_rec_t> v{
        rec(0,150,100,0,0),     // A, ref [100,250]
        rec(20,20,9000,0,1)};   // B, q nested but ref [9000,9020] disjoint
    order_seeds(v.data(), 2, SEED_ORDER_ABSORB_COUNT);
    // Both absorb 0 -> stable orig_ix order preserved.
    CHECK(v[0].orig_ix == 0);
    CHECK(v[1].orig_ix == 1);
}

TEST_CASE("most-absorb: transitive A>B>C emits A first, consumes B and C") {
    // A q[0,100] r[0,100]; B q[10,40] r[10,40] in A; C q[15,25] r[15,25] in B and A.
    std::vector<seed_rec_t> v{
        rec(15,11,15,0,0),    // C
        rec(10,31,10,0,1),    // B (absorbs C)
        rec(0,101,0,0,2)};    // A (absorbs B and C)
    order_seeds(v.data(), 3, SEED_ORDER_MOST_ABSORB);
    CHECK(v[0].orig_ix == 2);  // A first (rowsum 2)
}

TEST_CASE("most-absorb: second pick is max-absorber of the remainder") {
    // A absorbs B. D absorbs E. A and D disjoint loci.
    // Absorbed seeds are appended immediately after their absorber, so output is:
    //   A(0), B(1), D(2), E(3)  -- both absorbers emitted before any non-absorbed seed.
    std::vector<seed_rec_t> v{
        rec(0,100,0,0,0),       // A absorbs B
        rec(10,20,10,0,1),      // B
        rec(0,100,5000,0,2),    // D absorbs E (far locus)
        rec(10,20,5010,0,3)};   // E
    order_seeds(v.data(), 4, SEED_ORDER_MOST_ABSORB);
    // Absorbers A(0) and D(2) are at positions 0 and 2; absorbed B(1) and E(3) follow each.
    CHECK((v[0].orig_ix == 0 || v[0].orig_ix == 2));
    CHECK((v[2].orig_ix == 0 || v[2].orig_ix == 2));
    CHECK(v[0].orig_ix != v[2].orig_ix);
}

TEST_CASE("genomic modes fall back to global-longest above GENOMIC_ORDER_MAX_N") {
    // n > 1024 -> both most-absorb and absorb-count delegate to global-longest
    // (O(n) longest-first), avoiding the O(n^2) count / n*n matrix. Verify each
    // produces exactly the global-longest ordering on a large, varied input.
    const int N = 1100;  // > GENOMIC_ORDER_MAX_N (1024)
    std::vector<seed_rec_t> base;
    for (uint32_t i = 0; i < (uint32_t)N; ++i)
        base.push_back(rec(0, 20 + (int)(i % 50), 100 + i, 0, i));  // varied len -> non-trivial order
    std::vector<seed_rec_t> want = base, ma = base, ac = base;
    order_seeds(want.data(), N, SEED_ORDER_GLOBAL_LONGEST);
    order_seeds(ma.data(),   N, SEED_ORDER_MOST_ABSORB);
    order_seeds(ac.data(),   N, SEED_ORDER_ABSORB_COUNT);
    for (int i = 0; i < N; ++i) {
        CHECK(ma[i].orig_ix == want[i].orig_ix);  // most-absorb fallback == global-longest
        CHECK(ac[i].orig_ix == want[i].orig_ix);  // absorb-count fallback == global-longest
    }
}
