// test/unit/test_mem_extend_patch.cpp
//
// Region patching (mem_patch_reg, src/bwamem.cpp, reached through
// mem_sort_dedup_patch): when one read yields two alignment regions a and b
// that are collinear on one strand of one contig (a before b on both the query
// and the reference), bwa-mem joins them into one region if a global alignment
// of query[a.qb, b.qe) to reference[a.rb, b.re) scores at least 0.90 of what
// the two regions predict for that span. The gates, as bwa-mem defines them:
//
//   - w = |(a.re - b.rb) - (a.qe - b.qb)|, the diagonal shift between them,
//     and r = |(a.re - b.rb)/(b.re - a.rb) - (a.qe - b.qb)/(b.qe - a.qb)|, the
//     shift relative to the span;
//   - regions apart on the query or the reference: w <= 2*opt->w and r < 0.05;
//     overlapping on both: w <= 4*opt->w and r < 0.10 (the limits are the
//     floats 0.05f and 0.10f, a little above 0.05 and 0.10, so an r of
//     exactly 0.05 passes);
//   - the global alignment (bwa_gen_cigar2) is given a band of
//     min(w + a.w + b.w, 4*opt->w), which it may narrow or floor;
//   - predicted score: the two scores scaled by span / summed lengths, on the
//     query (q_s) and on the reference (r_s), each rounded; the patch needs
//     score / max(q_s, r_s) >= 0.90 (the float 0.90f, a little under 0.9, so
//     no integer score can tell >= from >);
//   - mem_sort_dedup_patch then gives the later region a's start, the patched
//     score as score and truesc, the band as w, the larger seedcov, sub and
//     csub of the two, and drops a.
//
// Every expected value below is worked out by hand in its comment, on a random
// reference packed in the test (no index needed: patching reads only l_pac and
// the packed bases).

#include "doctest/doctest.h"

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <vector>

#include "bntseq.h"
#include "bwamem.h"

// Non-static in bwamem.cpp but not declared in bwamem.h.
int mem_patch_reg(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, uint8_t *query,
                  const mem_alnreg_t *a, const mem_alnreg_t *b, int *_w, const int8_t *mat);
int mem_sort_dedup_patch(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, uint8_t *query, int n,
                         mem_alnreg_t *a, const int8_t *mat);

namespace {

struct Fixture {
    std::vector<uint8_t> fwd, pac;
    bntseq_t bns;
    mem_opt_t *opt;
    Fixture() : opt(mem_opt_init()) {
        std::mt19937 rng(27182u);
        fwd.resize(2000);
        for (auto &c : fwd) c = (uint8_t)(rng() % 4);
        pac.assign(fwd.size() / 4 + 1, 0);
        for (size_t k = 0; k < fwd.size(); ++k) _set_pac(pac, k, fwd[k]);
        memset(&bns, 0, sizeof bns);
        bns.l_pac = (int64_t)fwd.size();
    }
    ~Fixture() { free(opt); }
    Fixture(const Fixture &) = delete;
    Fixture &operator=(const Fixture &) = delete;
    std::vector<uint8_t> ref(size_t b, size_t e) const {
        return std::vector<uint8_t>(fwd.begin() + (long)b, fwd.begin() + (long)e);
    }
};

mem_alnreg_t region(int qb, int qe, int64_t rb, int64_t re, int score, int w) {
    mem_alnreg_t r;
    memset(&r, 0, sizeof r);
    r.qb = qb, r.qe = qe, r.rb = rb, r.re = re, r.score = r.truesc = score, r.w = w;
    r.secondary = -1;
    return r;
}

// Substitute query bases so they differ from the reference there.
void substitute(std::vector<uint8_t> &q, std::initializer_list<int> at) {
    for (int k : at) q[(size_t)k] = (uint8_t)((q[(size_t)k] + 1) % 4);
}

}  // namespace

TEST_CASE("mem_patch_reg: regions that are not collinear are not joined"
          * doctest::test_suite("unit/alnreg_sort_dedup")) {
    Fixture f;
    int w = -1;
    // Both pairs pass the band, relative-band and score tests (w = 0, r = 0, a
    // perfect global alignment), so only the collinearity guard refuses them.
    // b starting on the query where a does: not collinear. a = [10,60)x
    // [500,550), b = [10,100)x[500,590); both overlap gates give w = r = 0.
    std::vector<uint8_t> q = f.ref(490, 590);
    const mem_alnreg_t a0 = region(10, 60, 500, 550, 50, 10), b0 = region(10, 100, 500, 590, 90, 10);
    CHECK(mem_patch_reg(f.opt, &f.bns, f.pac.data(), q.data(), &a0, &b0, &w, f.opt->mat) == 0);
    // b ending on the query and the reference where a does: not collinear.
    const mem_alnreg_t a1 = region(0, 90, 490, 580, 90, 10), b1 = region(10, 90, 500, 580, 80, 10);
    CHECK(mem_patch_reg(f.opt, &f.bns, f.pac.data(), q.data(), &a1, &b1, &w, f.opt->mat) == 0);
    CHECK(w == -1);
}

TEST_CASE("mem_patch_reg: the band and relative-band gates, at their edges"
          * doctest::test_suite("unit/alnreg_sort_dedup")) {
    Fixture f;
    // The read is reference 500-699 then reference 700+d to 899+d: a d-base
    // deletion. a = [0,200)x[500,700), b = [200,400)x[700+d,900+d), scores 200
    // and 200, band 10 each. Apart on the reference, so w = |(700 - (700+d)) -
    // (200 - 200)| = d must be <= 2*opt->w and r = d / (400 + d) under 0.05.
    // Global score 400 - (6 + d); q_s = round(400/400 * 400) = 400 and
    // r_s = round((400+d)/400 * 400) = 400 + d.
    auto patch = [&](int d, int opt_w, int *w) {
        std::vector<uint8_t> q = f.ref(500, 700), tail = f.ref(700 + (size_t)d, 900 + (size_t)d);
        q.insert(q.end(), tail.begin(), tail.end());
        f.opt->w = opt_w;
        const mem_alnreg_t a = region(0, 200, 500, 700, 200, 10), b = region(200, 400, 700 + d, 900 + d, 200, 10);
        *w = -1;
        return mem_patch_reg(f.opt, &f.bns, f.pac.data(), q.data(), &a, &b, w, f.opt->mat);
    };
    int w = -1;
    // d = 4: r = 4/404; score 400 - 10 = 390, 390 / 404 = 0.965. Band
    // min(4 + 10 + 10, 4*100) = 24.
    CHECK(patch(4, 100, &w) == 390);
    CHECK(w == 24);
    // At -w 1 the band gate allows w <= 2: d = 2 passes (392 / 402 = 0.975,
    // band min(2 + 20, 4) = 4), d = 3 does not.
    CHECK(patch(2, 1, &w) == 392);
    CHECK(w == 4);
    CHECK(patch(3, 1, &w) == 0);
    CHECK(w == -1);

    // The relative band alone, with region scores low enough that the score
    // test would pass: a = [0,30)x[500,530) and b = [30,60)x[530+d,560+d),
    // scores 25 and 25. q_s = round(60/60 * 50) = 50, r_s = round((60+d)/60 * 50).
    auto patch_short = [&](int d) {
        std::vector<uint8_t> q = f.ref(500, 530), tail = f.ref(530 + (size_t)d, 560 + (size_t)d);
        q.insert(q.end(), tail.begin(), tail.end());
        f.opt->w = 100;
        const mem_alnreg_t a = region(0, 30, 500, 530, 25, 10), b = region(30, 60, 530 + d, 560 + d, 25, 10);
        int ww = -1;
        return mem_patch_reg(f.opt, &f.bns, f.pac.data(), q.data(), &a, &b, &ww, f.opt->mat);
    };
    // d = 4: r = 4/64 = 0.0625 >= 0.05: refused, although 50 / 53 = 0.94.
    CHECK(patch_short(4) == 0);

    // r of exactly 0.05: a = [0,30)x[500,530) scoring 25 and b =
    // [30,57)x[533,560) scoring 22, a 3-base deletion over a 60-base reference
    // span: r = 3/60 = 0.05, under the float limit 0.05f. Score 57 - 9 = 48,
    // q_s = 47, r_s = round(60/57 * 47) = 49; 48 / 49 = 0.98: joined.
    {
        std::vector<uint8_t> q = f.ref(500, 530), tail = f.ref(533, 560);
        q.insert(q.end(), tail.begin(), tail.end());
        f.opt->w = 100;
        const mem_alnreg_t a = region(0, 30, 500, 530, 25, 10), b = region(30, 57, 533, 560, 22, 10);
        int ww = -1;
        CHECK(mem_patch_reg(f.opt, &f.bns, f.pac.data(), q.data(), &a, &b, &ww, f.opt->mat) == 48);
    }

    // An insertion makes the query span the longer one, so q_s decides: the
    // read is reference 500-699, d inserted bases, then reference 700-899.
    // a = [0,200)x[500,700), b = [200+d,400+d)x[700,900), scores 200 each.
    // w = d, r = d / (400 + d). Score 400 - (6 + d); q_s = 400 + d, r_s = 400.
    // d = 17: 377 / 417 = 0.904, joined. d = 19: 375 / 419 = 0.895, refused,
    // although 375 / r_s = 0.94 would pass.
    auto patch_ins = [&](int d) {
        std::vector<uint8_t> q = f.ref(500, 700), tail = f.ref(700, 900);
        for (int k = 0; k < d; ++k) q.push_back((uint8_t)((q.back() + 1 + k) % 4));
        q.insert(q.end(), tail.begin(), tail.end());
        f.opt->w = 100;
        const mem_alnreg_t a = region(0, 200, 500, 700, 200, 10), b = region(200 + d, 400 + d, 700, 900, 200, 10);
        int ww = -1;
        return mem_patch_reg(f.opt, &f.bns, f.pac.data(), q.data(), &a, &b, &ww, f.opt->mat);
    };
    CHECK(patch_ins(17) == 377);
    CHECK(patch_ins(19) == 0);
}

TEST_CASE("mem_patch_reg: regions overlapping on both axes use the wider gates"
          * doctest::test_suite("unit/alnreg_sort_dedup")) {
    Fixture f;
    // The read is reference 500-799 then 820-1119 (a 20-base deletion), with
    // regions that overlap on both axes: a = [0,300)x[500,800) and
    // b = [280,600)x[800,1120), each scoring 300 with band 0. a.qe 300 > b.qb
    // 280 and a.re 800 is not below b.rb 800, so the overlap gates apply:
    // w = |(800 - 800) - (300 - 280)| = 20 <= 4*opt->w needs -w 5 or more;
    // r = |0/620 - 20/600| = 0.033 < 0.10. Band min(20 + 0 + 0, 4*opt->w).
    // Global score 600 - (6 + 20) = 574; q_s = round(600/620 * 600) = 581,
    // r_s = round(620/620 * 600) = 600; 574 / 600 = 0.957.
    std::vector<uint8_t> q = f.ref(500, 800), tail = f.ref(820, 1120);
    q.insert(q.end(), tail.begin(), tail.end());
    const mem_alnreg_t a = region(0, 300, 500, 800, 300, 0), b = region(280, 600, 800, 1120, 300, 0);
    int w = -1;
    f.opt->w = 5;
    CHECK(mem_patch_reg(f.opt, &f.bns, f.pac.data(), q.data(), &a, &b, &w, f.opt->mat) == 574);
    CHECK(w == 20);
    f.opt->w = 4;
    w = -1;
    CHECK(mem_patch_reg(f.opt, &f.bns, f.pac.data(), q.data(), &a, &b, &w, f.opt->mat) == 0);
}

TEST_CASE("mem_sort_dedup_patch: a patched pair becomes one region spanning both"
          * doctest::test_suite("unit/alnreg_sort_dedup")) {
    Fixture f;
    // The read is reference 500-599 with read bases 50 and 51 changed; the
    // regions either side are a = [0,50)x[500,550) scoring 50 and
    // b = [52,100)x[552,600) scoring 48, each with band 10. w = r = 0, band
    // min(0 + 10 + 10, 400) = 20. Global score 98 - 2*4 = 90; q_s = r_s =
    // round(100 / (48 + 50) * 98) = 100; 90 / 100 = 0.90 reaches the join
    // threshold. b takes a's start, score and truesc 90, band 20; a is
    // dropped. The larger seedcov, sub and csub of the two are kept.
    std::vector<uint8_t> q = f.ref(500, 600);
    substitute(q, {50, 51});
    std::vector<mem_alnreg_t> regs = {region(0, 50, 500, 550, 50, 10), region(52, 100, 552, 600, 48, 10)};
    regs[0].seedcov = 50, regs[1].seedcov = 48;
    regs[0].sub = 7, regs[1].sub = 3;
    regs[0].csub = 2, regs[1].csub = 9;
    const int n = mem_sort_dedup_patch(f.opt, &f.bns, f.pac.data(), q.data(), 2, regs.data(), f.opt->mat);
    REQUIRE(n == 1);
    const mem_alnreg_t &m = regs[0];
    CHECK(m.qb == 0);
    CHECK(m.qe == 100);
    CHECK(m.rb == 500);
    CHECK(m.re == 600);
    CHECK(m.score == 90);
    CHECK(m.truesc == 90);
    CHECK(m.w == 20);
    CHECK(m.seedcov == 50);
    CHECK(m.sub == 7);
    CHECK(m.csub == 9);

    // With three changed bases (global 97 - 3*4 = 85, q_s = 100; 0.85 < 0.90)
    // the pair stays two regions.
    std::vector<uint8_t> q3 = f.ref(500, 600);
    substitute(q3, {50, 51, 52});
    std::vector<mem_alnreg_t> regs3 = {region(53, 100, 553, 600, 47, 10), region(0, 50, 500, 550, 50, 10)};
    CHECK(mem_sort_dedup_patch(f.opt, &f.bns, f.pac.data(), q3.data(), 2, regs3.data(), f.opt->mat) == 2);
}
