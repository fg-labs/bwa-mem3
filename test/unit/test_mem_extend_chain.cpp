// test/unit/test_mem_extend_chain.cpp
//
// Chain weighting and the seed SW screen (src/bwamem.cpp), steps between
// seeding and extension, checked against expectations derived here from bwa's
// documented model. (The chain filter between them, mem_chain_flt, is held to
// an independent transcription in test_chain_flt_survivors.cpp.)
//
//   - mem_chain_weight: a chain's weight is the number of query bases its seeds
//     cover or the number of reference bases they cover, whichever is smaller,
//     overlapping seeds counted once (bwa-mem's "seeded bases", -W);
//   - mem_flt_chained_seeds / mem_seed_sw: on a read long enough for the
//     screen (min_l at most 0.05 * read length, where min_l is 1.1 * -W, or
//     5.5 * ln(read length) without -W), each seed under 200 bases is rescored
//     by a local SW (ksw_align2) of the read and the reference within 50 bases
//     around it, the reference window clamped to the contig and strand of the
//     seed's middle; it is dropped when that score is under round(a * min_l),
//     and a kept seed's score becomes the SW score (a seed of 200 bases or
//     more, or whose read or reference window is, is kept with score
//     a * length instead).
//
// Seeds are given in query and reference order, as chaining produces them.

#include "doctest/doctest.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <vector>

#include "bwamem.h"
#include "ksw.h"
#include "temp_index.h"

namespace {

struct OptGuard {
    mem_opt_t *opt;
    OptGuard() : opt(mem_opt_init()) {}
    ~OptGuard() { free(opt); }
};

mem_seed_t seed(int64_t rbeg, int32_t qbeg, int32_t len) {
    mem_seed_t s;
    s.rbeg = rbeg;
    s.qbeg = qbeg;
    s.len = len;
    s.score = len;
    return s;
}

// A chain over caller-owned seeds; m = 1 (SEEDS_PER_CHAIN) marks the seeds as
// not the chain's to free, so the filter never frees them. pos tags the chain.
mem_chain_t chain(std::vector<mem_seed_t> &seeds, int64_t pos, int is_alt = 0) {
    mem_chain_t c;
    memset(&c, 0, sizeof c);
    c.n = (int32_t)seeds.size();
    c.m = 1;
    c.seeds = seeds.data();
    c.pos = pos;
    c.is_alt = (uint32_t)is_alt;
    c.meth_hypothesis = -1;
    return c;
}

}  // namespace

TEST_CASE("mem_chain_weight: the smaller of query and reference coverage, overlaps counted once"
          * doctest::test_suite("unit/chain")) {
    struct Case {
        const char *what;
        std::vector<mem_seed_t> seeds;
        int want;
    };
    const Case cases[] = {
        {"one seed", {seed(1000, 0, 30)}, 30},
        {"two seeds overlapping on both axes: 30 + 30 - 10", {seed(1000, 0, 30), seed(1020, 20, 30)}, 50},
        {"a deletion between seeds: both axes 60", {seed(1000, 0, 30), seed(1040, 30, 30)}, 60},
        {"overlap on the query only: query 50, reference 60", {seed(1000, 0, 30), seed(1040, 20, 30)}, 50},
        {"overlap on the reference only: query 60, reference 50", {seed(1000, 0, 30), seed(1020, 40, 30)}, 50},
        {"a seed inside another: 40", {seed(1000, 0, 40), seed(1010, 10, 20)}, 40},
    };
    for (const Case &c : cases) {
        CAPTURE(c.what);
        std::vector<mem_seed_t> seeds = c.seeds;
        const mem_chain_t ch = chain(seeds, 0);
        CHECK(mem_chain_weight(&ch) == c.want);
    }
}

TEST_CASE("mem_flt_chained_seeds: the seed SW screen keeps a seed scoring exactly the threshold"
          * doctest::test_suite("unit/flt_chained_seeds")) {
    // Two contigs indexed in a temporary directory: ctgA (3000 bases) then ctgB
    // (1500), so ctgB starts at forward coordinate 3000 and l_pac is 4500. The
    // 1000-base read is random except for:
    //   [0, 400)    ctgA 400-799                      seed L: qbeg 5, len 40
    //   [600, 615)  ctgA 2500-2514                    seed S: a 15-base seed alone
    //   [680, 720)  ctgB 1460-1499, then [720, 750) the reverse complement of
    //               ctgB 1470-1499 (the reverse strand just past l_pac)
    //                                                 seed E: next to l_pac
    //   [750, 800)  ctgA 2950-2999 (ctgA's tail)
    //   [800, 845)  ctgB 5-49                         seed B: near ctgB's start
    //   [900, 950)  reverse complement of ctgA 1500-1549   seed R: reverse strand
    std::mt19937 rng(4242u);
    std::string a(3000, 'A'), b(1500, 'A'), read(1000, 'A');
    for (auto &c : a) c = "ACGT"[rng() % 4];
    for (auto &c : b) c = "ACGT"[rng() % 4];
    for (auto &c : read) c = "ACGT"[rng() % 4];
    bwa_tests::TempIndex idx({{"ctgA", a}, {"ctgB", b}});
    const int64_t L = 4500;
    REQUIRE(idx.bns()->l_pac == L);
    read.replace(0, 400, a.substr(400, 400));
    read.replace(600, 15, a.substr(2500, 15));
    read.replace(680, 40, b.substr(1460, 40));
    read.replace(720, 30, bwa_tests::revcomp(b.substr(1470, 30)));
    read.replace(750, 50, a.substr(2950, 50));
    read.replace(800, 45, b.substr(5, 45));
    read.replace(900, 50, bwa_tests::revcomp(a.substr(1500, 50)));
    std::vector<uint8_t> enc = bwa_tests::nt4(read);
    bseq1_t rec;
    memset(&rec, 0, sizeof rec);
    rec.l_seq = (int)read.size();
    rec.seq = (char *)enc.data();

    const std::string fwd = a + b;
    // Doubled coordinate x: forward base below L, complement of 2L-1-x above.
    auto base = [&](int64_t x) -> char {
        if (x < L) return fwd[(size_t)x];
        const char c = fwd[(size_t)(2 * L - 1 - x)];
        return c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : 'A';
    };
    // The window mem_seed_sw scores: 50 bases around the seed on both sides,
    // clamped to the read, and on the reference to the contig and strand of the
    // seed's middle base. Score: a local SW with gaps opened from H.
    auto window_score = [&](const mem_seed_t &sd) {
        const int qb = std::max(0, sd.qbeg - 50), qe = std::min((int)read.size(), sd.qbeg + sd.len + 50);
        const int64_t mid = sd.rbeg + sd.len / 2;
        const bool rev = mid >= L;
        const int64_t fmid = rev ? 2 * L - 1 - mid : mid;
        const int64_t cb = fmid < 3000 ? 0 : 3000, ce = fmid < 3000 ? 3000 : L;  // contig, forward
        const int64_t lo = rev ? 2 * L - ce : cb, hi = rev ? 2 * L - cb : ce;   // on the seed's strand
        const int64_t rb = std::max(lo, sd.rbeg - 50), re = std::min(hi, sd.rbeg + sd.len + 50);
        std::string t;
        for (int64_t x = rb; x < re; ++x) t += base(x);
        const std::vector<uint8_t> qv = bwa_tests::nt4(read.substr((size_t)qb, (size_t)(qe - qb)));
        const std::vector<uint8_t> tv = bwa_tests::nt4(t);
        int best = 0;
        const size_t W = qv.size() + 1;
        std::vector<int> H(W * (tv.size() + 1), 0), E(H), F(H);
        for (size_t i = 1; i <= tv.size(); ++i)
            for (size_t j = 1; j <= qv.size(); ++j) {
                const int sc = tv[i - 1] == qv[j - 1] ? 1 : -4;
                E[i * W + j] = std::max(H[(i - 1) * W + j] - 7, E[(i - 1) * W + j] - 1);
                F[i * W + j] = std::max(H[i * W + j - 1] - 7, F[i * W + j - 1] - 1);
                H[i * W + j] = std::max(std::max(0, H[(i - 1) * W + j - 1] + sc), std::max(E[i * W + j], F[i * W + j]));
                best = std::max(best, H[i * W + j]);
            }
        return best;
    };
    const mem_seed_t sL = seed(405, 5, 40), sS = seed(2500, 600, 15), sB = seed(3005, 800, 45),
                     sR = seed(2 * L - 1550, 900, 50);
    const int scL = window_score(sL), scS = window_score(sS), scB = window_score(sB), scR = window_score(sR);
    CHECK(scL >= 95);                    // read[0, 95) lies inside ctgA 355-494
    CHECK((scS >= 15 && scS < 30));      // the 15 bases plus at most chance matches
    // Seed B's window stops at ctgB's first base: read [750, 800) matches
    // ctgA's tail, so an unclamped window would score about 80.
    CHECK((scB >= 45 && scB < 60));
    CHECK((scR >= 50 && scR < 65));
    OptGuard g;

    // -W values whose threshold round(1.1 * W) brackets seed S's score; both
    // leave the screen on (1.1 * W <= 0.05 * 1000).
    int w_keep = 0;
    for (int w = 1; w <= 45; ++w)
        if ((int)(1.1f * w + .499) <= scS) w_keep = w;
    REQUIRE(w_keep > 0);
    REQUIRE((int)(1.1f * (w_keep + 1) + .499) > scS);
    CHECK((int)(1.1f * w_keep + .499) == scS);  // exactly the threshold
    for (int w : {w_keep, w_keep + 1}) {
        CAPTURE(w);
        g.opt->min_chain_weight = w;
        std::vector<mem_seed_t> seeds = {sL, sS, sB, sR};
        mem_chain_t c = chain(seeds, 0);
        mem_flt_chained_seeds(g.opt, idx.bns(), idx.pac(), &rec, 1, &c);
        // The survivors keep their order and carry their SW score, as the
        // window model computes it.
        std::vector<int> qbeg, score;
        for (int k = 0; k < c.n; ++k) qbeg.push_back(c.seeds[k].qbeg), score.push_back(c.seeds[k].score);
        if (w == w_keep) {
            CHECK(qbeg == std::vector<int>{5, 600, 800, 900});
            CHECK(score == std::vector<int>{scL, scS, scB, scR});
        } else {
            CHECK(qbeg == std::vector<int>{5, 800, 900});
            CHECK(score == std::vector<int>{scL, scB, scR});
        }
    }
    // Without -W the threshold is round(5.5 * ln 1000) = 38: seed S goes.
    g.opt->min_chain_weight = 0;
    {
        std::vector<mem_seed_t> seeds = {sL, sS, sB, sR};
        mem_chain_t c = chain(seeds, 0);
        mem_flt_chained_seeds(g.opt, idx.bns(), idx.pac(), &rec, 1, &c);
        REQUIRE(c.n == 3);
        CHECK(c.seeds[1].qbeg == 800);
    }
    // Seed E (ctgB 1480-1499) has its 50-base flank past l_pac, where the
    // reverse strand continues with the read's next 30 bases: the window stops
    // at l_pac, so it scores the 40 forward bases (about 40, not about 70).
    {
        const mem_seed_t sE = seed(4480, 700, 20);
        const int scE = window_score(sE);
        CHECK((scE >= 40 && scE < 55));
        std::vector<mem_seed_t> seeds = {sE};
        mem_chain_t c = chain(seeds, 0);
        mem_flt_chained_seeds(g.opt, idx.bns(), idx.pac(), &rec, 1, &c);
        REQUIRE(c.n == 1);
        CHECK(c.seeds[0].score == scE);
    }
    // A 150-base seed with its full flanks has a 250-base window, at or over
    // 200: it is not rescored, and is kept with score a * 150 (rescored, its
    // window would score 250: read [0, 400) lies inside ctgA 400-799).
    {
        mem_seed_t mid = seed(500, 100, 150);
        mid.score = 7;
        std::vector<mem_seed_t> seeds = {mid};
        mem_chain_t c = chain(seeds, 0);
        mem_flt_chained_seeds(g.opt, idx.bns(), idx.pac(), &rec, 1, &c);
        REQUIRE(c.n == 1);
        CHECK(c.seeds[0].score == 150);
    }
    // A seed of 200 bases is not rescored: it is kept with score a * 200
    // (its window would score 250: read [0, 250) lies inside ctgA 400-649).
    {
        mem_seed_t longseed = seed(400, 0, 200);
        longseed.score = 1;  // so that a kept, unchanged score shows
        std::vector<mem_seed_t> seeds = {longseed, sS};
        mem_chain_t c = chain(seeds, 0);
        mem_flt_chained_seeds(g.opt, idx.bns(), idx.pac(), &rec, 1, &c);
        REQUIRE(c.n == 1);
        CHECK(c.seeds[0].qbeg == 0);
        CHECK(c.seeds[0].score == 200);
    }
}
