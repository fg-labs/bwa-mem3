// test/unit/test_extend_csub.cpp — mem_seed_capped_csub (--extend-csub), T3 per-region.
//
// When --max-extend-chains / --extend-tie-frac / --extend-dom-gap prunes chains
// before banded-SW, the pruned competitor is never scored, so it never reaches
// `csub`/`sub_n` in mem_approx_mapq_se and MAPQ comes out overstated.
// mem_seed_capped_csub seeds them back in from the dropped-chain records (capped_t),
// PER REGION: a dropped chain demotes a region only where it OVERLAPS it on the
// query (the mask_level rule mem_mark_primary_se_core applies to real regions), a
// non-ALT competitor never demotes an ALT region, the estimate is CLAMPED strictly
// below the region's own score (else `sub >= score -> mapq 0`), and sub_n increments
// on an estimated near-tie by the same band as `sub`.
//
// Every expectation is derived from the contract, not recorded from a run. Fixtures
// are built in memory; no test data files are read.

#include <cstring>
#include <vector>

#include "doctest/doctest.h"
#include "bwamem.h"

namespace {

// A region spanning [qb,qe) with the given score, pre-existing csub, secondary flag,
// and ALT flag. secondary < 0 marks it primary (what mem_mark_primary_se leaves).
mem_alnreg_t reg(int qb, int qe, int score, int csub = 0, int secondary = -1, int is_alt = 0) {
    mem_alnreg_t a;
    memset(&a, 0, sizeof(a));
    a.qb = qb; a.qe = qe;
    a.score = a.truesc = score;
    a.csub = csub;
    a.secondary = secondary;
    a.is_alt = is_alt;
    a.sub_n = 0;
    a.meth_hypothesis = a.meth_strand_hyp = -1;
    return a;
}

// A dropped-chain record spanning [qb,qe) with estimated score `est` and ALT flag.
capped_chain_t cap(int qb, int qe, int est, int is_alt = 0) {
    capped_chain_t d; memset(&d, 0, sizeof d);
    d.qb = qb; d.qe = qe; d.w = est; d.est = est; d.is_alt = (int8_t) is_alt;
    return d;
}

struct RegV {
    std::vector<mem_alnreg_t> storage;
    mem_alnreg_v v;
    RegV(std::vector<mem_alnreg_t> regs, std::vector<capped_chain_t> caps) : storage(std::move(regs)) {
        memset(&v, 0, sizeof(v));
        v.a = storage.data();
        v.n = v.m = storage.size();
        capped_clear(&v.capped);
        for (auto &c : caps) capped_push(&v.capped, &c);
    }
};

// Defaults: a=1 b=4 o_del=6 e_del=1 o_ins=6 e_ins=1 mask_level=0.5, so the sub_n
// near-tie band tmp = max(1+4, 6+1, 6+1) = 7.
struct Opt {
    mem_opt_t *o;
    Opt() { o = mem_opt_init(); }
    ~Opt() { free(o); }
};

}  // namespace

TEST_CASE("csub: no-op when nothing was dropped (off-state byte-identity)") {
    Opt opt;
    RegV r({reg(0, 100, 80, /*csub=*/7)}, /*caps=*/{});
    mem_seed_capped_csub(&r.v, opt.o);
    CHECK(r.v.a[0].csub == 7);
    CHECK(r.v.a[0].sub_n == 0);
}

TEST_CASE("csub: an overlapping dropped chain seeds csub, clamped below score") {
    Opt opt;
    // region [0,100) score 80; dropped record spans [0,100) est 30 -> csub=30 (< score).
    RegV r({reg(0, 100, 80)}, {cap(0, 100, 30)});
    mem_seed_capped_csub(&r.v, opt.o);
    CHECK(r.v.a[0].csub == 30);
    CHECK(r.v.a[0].sub_n == 0);            // 80-30 = 50 > 7
}

TEST_CASE("csub: est at/above score is clamped to score-1 and counts as a near-tie") {
    Opt opt;
    RegV r({reg(0, 100, 80, /*csub=*/3)}, {cap(0, 100, 10000)});
    mem_seed_capped_csub(&r.v, opt.o);
    CHECK(r.v.a[0].csub == 79);            // clamped
    CHECK(r.v.a[0].sub_n == 1);            // 80-79 = 1 <= 7
}

TEST_CASE("csub: per-region — a dropped chain seeds only the region it overlaps") {
    Opt opt;
    // split read: L=[0,80) score 80, R=[80,150) score 70 (both primary).
    // dropped record spans [0,75) est 60 -> overlaps L (75>=0.5*80) but not R.
    RegV r({reg(0, 80, 80), reg(80, 150, 70)}, {cap(0, 75, 60)});
    mem_seed_capped_csub(&r.v, opt.o);
    CHECK(r.v.a[0].csub == 60);            // L seeded
    CHECK(r.v.a[1].csub == 0);             // R untouched
    CHECK(r.v.a[1].sub_n == 0);
}

TEST_CASE("csub: sub_n increments iff est >= score - tmp (tmp = 7)") {
    Opt opt;
    {   RegV r({reg(0, 100, 80)}, {cap(0, 100, 73)}); // 80-73 = 7 <= 7 -> +1
        mem_seed_capped_csub(&r.v, opt.o);
        CHECK(r.v.a[0].sub_n == 1); CHECK(r.v.a[0].csub == 73); }
    {   RegV r({reg(0, 100, 80)}, {cap(0, 100, 72)}); // 80-72 = 8 > 7 -> +0
        mem_seed_capped_csub(&r.v, opt.o);
        CHECK(r.v.a[0].sub_n == 0); CHECK(r.v.a[0].csub == 72); }
}

TEST_CASE("csub: a non-ALT competitor never demotes an ALT region (ALT protection)") {
    Opt opt;
    // rule mirrors mem_mark_primary_se_core: seed iff (competitor.is_alt || !region.is_alt).
    {   RegV r({reg(0, 100, 80, /*csub=*/0, /*secondary=*/-1, /*is_alt=*/1)},   // ALT region
               {cap(0, 100, 40, /*is_alt=*/0)});                                // non-ALT competitor
        mem_seed_capped_csub(&r.v, opt.o);
        CHECK(r.v.a[0].csub == 0);         // protected: non-ALT competitor cannot demote an ALT region
        CHECK(r.v.a[0].sub_n == 0); }
    {   RegV r({reg(0, 100, 80, /*csub=*/0, /*secondary=*/-1, /*is_alt=*/1)},   // ALT region
               {cap(0, 100, 40, /*is_alt=*/1)});                                // ALT competitor
        mem_seed_capped_csub(&r.v, opt.o);
        CHECK(r.v.a[0].csub == 40); }      // both ALT -> seeded
    {   RegV r({reg(0, 100, 80, /*csub=*/0, /*secondary=*/-1, /*is_alt=*/0)},   // non-ALT region
               {cap(0, 100, 40, /*is_alt=*/1)});                                // ALT competitor
        mem_seed_capped_csub(&r.v, opt.o);
        CHECK(r.v.a[0].csub == 40); }      // a non-ALT region can be demoted by anyone
}

TEST_CASE("csub: secondary regions are never seeded") {
    Opt opt;
    RegV r({reg(0, 100, 80, /*csub=*/3, /*secondary=*/0)}, {cap(0, 100, 40)});
    mem_seed_capped_csub(&r.v, opt.o);
    CHECK(r.v.a[0].csub == 3);             // untouched (secondary >= 0)
    CHECK(r.v.a[0].sub_n == 0);
}

TEST_CASE("csub: only raises csub, never lowers a larger existing value") {
    Opt opt;
    RegV r({reg(0, 100, 80, /*csub=*/50)}, {cap(0, 100, 30)});
    mem_seed_capped_csub(&r.v, opt.o);
    CHECK(r.v.a[0].csub == 50);            // 30 < 50, unchanged
}
