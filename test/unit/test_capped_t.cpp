// T1: capped_t helpers (capped_push / capped_overlaps / chain_qspan).
#include "doctest/doctest.h"
#include "bwamem.h"
#include <cstring>

namespace {
capped_chain_t cc(int w, int rid = 0) {
    capped_chain_t d; memset(&d, 0, sizeof d); d.w = w; d.est = w; d.rid = rid; return d;
}
capped_chain_t ccq(int qb, int qe) {
    capped_chain_t d; memset(&d, 0, sizeof d); d.qb = qb; d.qe = qe; return d;
}
void push(capped_t *c, int w, int rid = 0) { capped_chain_t d = cc(w, rid); capped_push(c, &d); }
}  // namespace

TEST_CASE("capped_push: keeps top-CAPPED_KEEP by weight, desc, stable on ties") {
    capped_t c; capped_clear(&c);
    push(&c, 10);
    push(&c, 50, /*rid=*/1);   // first 50
    push(&c, 90);
    push(&c, 50, /*rid=*/2);   // second 50
    push(&c, 30);
    push(&c, 70);
    REQUIRE(c.n == CAPPED_KEEP);            // 4
    CHECK(c.d[0].w == 90);
    CHECK(c.d[1].w == 70);
    CHECK(c.d[2].w == 50);                  // the two 50s survive over 30/10
    CHECK(c.d[3].w == 50);
    CHECK(c.d[2].rid == 1);                 // stable: first-inserted 50 ranks ahead
    CHECK(c.d[3].rid == 2);
}

TEST_CASE("capped_push: lighter-than-lightest when full is dropped") {
    capped_t c; capped_clear(&c);
    for (int w : {40, 30, 20, 10}) push(&c, w);
    push(&c, 5);                             // below the lightest (10) -> ignored
    CHECK(c.n == 4);
    CHECK(c.d[3].w == 10);
    push(&c, 25);                            // beats 10 -> replaces, re-sorts
    CHECK(c.d[0].w == 40); CHECK(c.d[1].w == 30); CHECK(c.d[2].w == 25); CHECK(c.d[3].w == 20);
}

TEST_CASE("capped_overlaps: mask_level boundary is inclusive") {
    capped_chain_t d = ccq(0, 80);           // span 80
    CHECK(capped_overlaps(&d, 40, 120, 0.5f) == 1);   // overlap [40,80)=40 -> exactly 40, true
    CHECK(capped_overlaps(&d, 41, 120, 0.5f) == 0);   // overlap 39 -> false
    CHECK(capped_overlaps(&d, 200, 280, 0.5f) == 0);  // disjoint
    CHECK(capped_overlaps(&d, 0, 20, 0.5f) == 1);     // smaller region span 20 -> need >=10, overlap 20
}

TEST_CASE("chain_qspan: min qbeg .. max qbeg+len over unsorted seeds") {
    mem_seed_t s[3];
    memset(s, 0, sizeof s);
    s[0].qbeg = 30; s[0].len = 20;           // [30,50)
    s[1].qbeg = 5;  s[1].len = 10;           // [5,15)   <- min
    s[2].qbeg = 60; s[2].len = 25;           // [60,85)  <- max end
    mem_chain_t c; memset(&c, 0, sizeof c);
    c.n = 3; c.seeds = s;
    int qb = -1, qe = -1; chain_qspan(&c, &qb, &qe);
    CHECK(qb == 5);
    CHECK(qe == 85);
}
