// test/unit/test_mem_extend_chain.cpp
//
// Chain weighting and filtering (src/bwamem.cpp), the steps between seeding and
// extension that decide which candidate locations of a read are extended at
// all, checked against expectations derived here from bwa's documented model:
//
//   - mem_chain_weight: a chain's weight is the number of query bases its seeds
//     cover or the number of reference bases they cover, whichever is smaller,
//     overlapping seeds counted once (bwa-mem's "seeded bases", -W);
//   - mem_chain_flt: chains lighter than -W are discarded; the rest, heaviest
//     first, are compared with the chains already kept. Two chains overlap
//     significantly when their query spans share at least mask_level (-X) of
//     the shorter span, and the shorter is under max_chain_gap (-G) long; a
//     kept ALT chain never shadows a primary one. A chain lighter than -D times
//     a significantly overlapping kept chain, and lighter by at least twice the
//     minimum seed length, is dropped; a chain kept with a significant overlap
//     is marked kept=2, one without kept=3. The first chain found overlapping
//     each kept chain is kept as well, with kept=1, so its score can inform the
//     mapping quality. Counting kept=1/2 chains heaviest first (kept=3 chains
//     are not counted), the -N-th one and every kept=1/2 chain after it are
//     dropped. Survivors stay heaviest first. What leaves the filter, and what
//     is compared, is which chains survive and in what order (the order they
//     are extended in); the kept marks and the stored weight decide that inside
//     the filter, and nothing after it reads them, so they are not compared;
//   - mem_flt_chained_seeds / mem_seed_sw: on a read long enough for the
//     screen (min_l at most 0.05 * read length, where min_l is 1.1 * -W, or
//     5.5 * ln(read length) without -W), each seed under 200 bases is rescored
//     by a local SW (ksw_align2) of the read and the reference within 50 bases
//     around it, the reference window clamped to the contig and strand of the
//     seed's middle; it is dropped when that score is under round(a * min_l),
//     and a kept seed's score becomes the SW score (a seed of 200 bases or
//     more, or whose window is, is kept with score a * length).
//
// The filter's order on equal weights is that of an unstable sort, so the
// generated chains here have distinct weights. Seeds are given in query and
// reference order, as chaining produces them.
//
// Random inputs use std::mt19937 raw output (rng() % n).

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
#include "mem_fixture.h"

// Non-static in bwamem.cpp but not declared in bwamem.h.
int mem_chain_flt(const mem_opt_t *opt, int n_chn_, mem_chain_t *a_, int tid);

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

// Union length of half-open intervals.
int covered(std::vector<std::pair<int64_t, int64_t> > iv) {
    std::sort(iv.begin(), iv.end());
    int64_t total = 0, end = INT64_MIN;
    for (const auto &x : iv) {
        if (x.second <= end) continue;
        total += x.second - std::max(x.first, end);
        end = x.second;
    }
    return (int)total;
}

int weight_oracle(const std::vector<mem_seed_t> &seeds) {
    std::vector<std::pair<int64_t, int64_t> > q, r;
    for (const mem_seed_t &s : seeds) {
        q.push_back(std::make_pair((int64_t)s.qbeg, (int64_t)s.qbeg + s.len));
        r.push_back(std::make_pair(s.rbeg, s.rbeg + s.len));
    }
    return std::min(covered(q), covered(r));
}

struct FltParams {
    int min_chain_weight, max_chain_gap, min_seed_len, max_chain_extend;
    float mask_level, drop_ratio;
};

struct Kept {
    int64_t pos;
    int kept, w;
};

// The filter as the header describes it, over chains given by their seeds.
std::vector<Kept> filter_oracle(const std::vector<std::vector<mem_seed_t> > &seeds, const std::vector<int> &alt,
                                const FltParams &f) {
    struct C {
        int idx, w, beg, end, alt;
    };
    std::vector<C> cs;
    for (size_t i = 0; i < seeds.size(); ++i) {
        const int w = weight_oracle(seeds[i]);
        if (f.min_chain_weight > 0 && w < f.min_chain_weight) continue;
        cs.push_back({(int)i, w, seeds[i].front().qbeg, seeds[i].back().qbeg + seeds[i].back().len, alt[i]});
    }
    std::sort(cs.begin(), cs.end(), [](const C &a, const C &b) { return a.w > b.w; });
    const int n = (int)cs.size();
    std::vector<int> kept(n, 0), first(n, -1), surv;
    if (n > 0) {
        kept[0] = 3;
        surv.push_back(0);
    }
    for (int i = 1; i < n; ++i) {
        bool large = false, dropped = false;
        for (int j : surv) {
            const int b = std::max(cs[j].beg, cs[i].beg), e = std::min(cs[j].end, cs[i].end);
            if (e <= b || (cs[j].alt && !cs[i].alt)) continue;
            const int min_l = std::min(cs[i].end - cs[i].beg, cs[j].end - cs[j].beg);
            if (e - b >= min_l * f.mask_level && min_l < f.max_chain_gap) {
                large = true;
                if (first[j] < 0) first[j] = i;
                if (cs[i].w < cs[j].w * f.drop_ratio && cs[j].w - cs[i].w >= 2 * f.min_seed_len) {
                    dropped = true;
                    break;
                }
            }
        }
        if (!dropped) {
            surv.push_back(i);
            kept[i] = large ? 2 : 3;
        }
    }
    for (int j : surv)
        if (first[j] >= 0) kept[first[j]] = 1;
    int cutoff = n;
    for (int i = 0, k = 0; i < n; ++i) {
        if (kept[i] == 0 || kept[i] == 3) continue;
        if (++k >= f.max_chain_extend) {
            cutoff = i;
            break;
        }
    }
    std::vector<Kept> out;
    for (int i = 0; i < n; ++i) {
        if (i >= cutoff && kept[i] < 3) kept[i] = 0;
        if (kept[i] != 0) out.push_back({(int64_t)cs[i].idx, kept[i], cs[i].w});
    }
    return out;
}

std::vector<Kept> run_filter(std::vector<std::vector<mem_seed_t> > &seeds, const std::vector<int> &alt,
                             const FltParams &f) {
    OptGuard g;
    g.opt->min_chain_weight = f.min_chain_weight;
    g.opt->max_chain_gap = f.max_chain_gap;
    g.opt->min_seed_len = f.min_seed_len;
    g.opt->max_chain_extend = f.max_chain_extend;
    g.opt->mask_level = f.mask_level;
    g.opt->drop_ratio = f.drop_ratio;
    std::vector<mem_chain_t> cs;
    for (size_t i = 0; i < seeds.size(); ++i) cs.push_back(chain(seeds[i], (int64_t)i, alt[i]));
    const int n = mem_chain_flt(g.opt, (int)cs.size(), cs.data(), 0);
    std::vector<Kept> out;
    for (int i = 0; i < n; ++i) out.push_back({cs[(size_t)i].pos, (int)cs[(size_t)i].kept, (int)cs[(size_t)i].w});
    return out;
}

std::string show(const std::vector<Kept> &v) {
    std::string s;
    for (const Kept &k : v)
        s += "(" + std::to_string(k.pos) + " kept " + std::to_string(k.kept) + " w " + std::to_string(k.w) + ")";
    return s;
}

// The filter's output contract: the same chains survive, in the same order.
// kept and w are listed by show() for a failing case but not compared (see the
// header).
bool same(const std::vector<Kept> &a, const std::vector<Kept> &b) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i)
        if (a[i].pos != b[i].pos) return false;
    return true;
}

FltParams defaults() { return {0, 10000, 19, 1 << 30, 0.50f, 0.50f}; }

// One seed per chain spanning query [qb, qe) on reference 1000 * id + qb.
std::vector<mem_seed_t> span(int id, int qb, int qe) { return {seed(1000 * (int64_t)id + qb, qb, qe - qb)}; }

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
        CHECK(weight_oracle(c.seeds) == c.want);
    }
    // Random collinear chains: seeds in query and reference order, with gaps,
    // overlaps and indels between them.
    std::mt19937 rng(1234u);
    for (int n = 0; n < 300; ++n) {
        std::vector<mem_seed_t> seeds;
        int32_t q = (int32_t)(rng() % 20);
        int64_t r = 5000 + (int64_t)(rng() % 20);
        for (int k = 1 + (int)(rng() % 6); k > 0; --k) {
            const int32_t len = 5 + (int32_t)(rng() % 40);
            seeds.push_back(seed(r, q, len));
            q += (int32_t)(rng() % (uint32_t)(len + 10));
            r += (int64_t)(rng() % (uint32_t)(len + 10));
        }
        CAPTURE(n);
        const mem_chain_t ch = chain(seeds, 0);
        CHECK(mem_chain_weight(&ch) == weight_oracle(seeds));
    }
}

TEST_CASE("mem_chain_flt: hand-worked drops, shadows and caps"
          * doctest::test_suite("unit/chain")) {
    // A spans query [0, 100) (weight 100). Each case adds lighter chains over
    // the same read; the expected survivors are listed heaviest first.
    SUBCASE("-D: a chain under half of an overlapping heavier one is dropped, the first is kept=1") {
        // B [0, 49) w 49 < 50 and 100 - 49 >= 38: dropped, but B is the first
        // chain to overlap A, so it comes back with kept=1. C [0, 30) w 30: dropped.
        std::vector<std::vector<mem_seed_t> > s = {span(0, 0, 100), span(1, 0, 49), span(2, 0, 30)};
        const std::vector<Kept> got = run_filter(s, {0, 0, 0}, defaults());
        const std::vector<Kept> want = {{0, 3, 100}, {1, 1, 49}};
        INFO(show(got));
        CHECK(same(got, want));
        CHECK(same(filter_oracle(s, {0, 0, 0}, defaults()), want));
    }
    SUBCASE("-D: at exactly half the weight a chain survives (the test is strict)") {
        // B0 [0, 90) w 90 survives and is A's first shadow. B [0, 50) w 50 is not
        // < 100 * 0.5, so it survives (kept=2); a non-strict test would drop it,
        // and it would not come back, A's first shadow being B0. C [0, 49) is
        // dropped by A (49 < 50 and 51 >= 38) and is nobody's first shadow.
        std::vector<std::vector<mem_seed_t> > s = {span(0, 0, 100), span(1, 0, 90), span(2, 0, 50),
                                                   span(3, 0, 49)};
        const std::vector<Kept> got = run_filter(s, {0, 0, 0, 0}, defaults());
        const std::vector<Kept> want = {{0, 3, 100}, {1, 1, 90}, {2, 1, 50}};
        INFO(show(got));
        CHECK(same(got, want));
        CHECK(same(filter_oracle(s, {0, 0, 0, 0}, defaults()), want));
    }
    SUBCASE("-D 0.3 keeps the lighter chains; each is the first shadow of the one before") {
        // -D is a float: 100 * 0.3f is 30.0000012, so a weight of 30 would
        // still be dropped; 31 is not.
        FltParams f = defaults();
        f.drop_ratio = 0.3f;
        std::vector<std::vector<mem_seed_t> > s = {span(0, 0, 100), span(1, 0, 49), span(2, 0, 31)};
        const std::vector<Kept> got = run_filter(s, {0, 0, 0}, f);
        const std::vector<Kept> want = {{0, 3, 100}, {1, 1, 49}, {2, 1, 31}};
        INFO(show(got));
        CHECK(same(got, want));
    }
    SUBCASE("a weight gap under twice the seed length keeps a light chain") {
        // A' w 60, B w 25: 25 < 30 but 60 - 25 = 35 < 2 * 19 = 38.
        std::vector<std::vector<mem_seed_t> > s = {span(0, 0, 60), span(1, 0, 25)};
        const std::vector<Kept> got = run_filter(s, {0, 0}, defaults());
        CHECK(same(got, {{0, 3, 60}, {1, 1, 25}}));
        // With -k 17 the gap is 35 >= 34: dropped, then kept=1 as the first shadow.
        FltParams f = defaults();
        f.min_seed_len = 17;
        CHECK(same(run_filter(s, {0, 0}, f), {{0, 3, 60}, {1, 1, 25}}));
        // A third chain shows the difference: with -k 17 it is dropped for good.
        std::vector<std::vector<mem_seed_t> > s3 = {span(0, 0, 60), span(1, 0, 25), span(2, 0, 24)};
        CHECK(same(run_filter(s3, {0, 0, 0}, defaults()), {{0, 3, 60}, {1, 1, 25}, {2, 1, 24}}));
        CHECK(same(run_filter(s3, {0, 0, 0}, f), {{0, 3, 60}, {1, 1, 25}}));
    }
    SUBCASE("-X: an overlap under mask_level of the shorter span is not significant") {
        // B [70, 110) overlaps A [0, 100) by 30 of its 40 bases: 30 >= 40*0.5,
        // significant, so B is dropped and comes back as A's first shadow
        // (kept=1); with -X 0.8 (needs 32) it is not, and B is kept=3. -N 1 makes
        // the difference visible: it cuts the first kept=1/2 chain but not a
        // kept=3 one.
        std::vector<std::vector<mem_seed_t> > s = {span(0, 0, 100), span(1, 70, 110)};
        FltParams f = defaults();
        f.max_chain_extend = 1;
        CHECK(same(run_filter(s, {0, 0}, f), {{0, 3, 100}}));
        f.mask_level = 0.8f;
        CHECK(same(run_filter(s, {0, 0}, f), {{0, 3, 100}, {1, 3, 40}}));
    }
    SUBCASE("-G: a shorter span of max_chain_gap or more is never a significant overlap") {
        std::vector<std::vector<mem_seed_t> > s = {span(0, 0, 100), span(1, 0, 45), span(2, 0, 44)};
        FltParams f = defaults();
        f.max_chain_gap = 45;  // B (45 long) is exempt, C (44) is not
        CHECK(same(run_filter(s, {0, 0, 0}, f), {{0, 3, 100}, {1, 3, 45}, {2, 1, 44}}));
    }
    SUBCASE("-W discards light chains first; when none is left, nothing is extended") {
        FltParams f = defaults();
        f.min_chain_weight = 45;
        std::vector<std::vector<mem_seed_t> > s = {span(0, 0, 100), span(1, 200, 245), span(2, 300, 344)};
        CHECK(same(run_filter(s, {0, 0, 0}, f), {{0, 3, 100}, {1, 3, 45}}));
        f.min_chain_weight = 101;
        CHECK(run_filter(s, {0, 0, 0}, f).empty());
    }
    SUBCASE("a kept ALT chain does not shadow a primary chain") {
        // A primary B under an ALT A is not overlapped as far as the filter is
        // concerned (kept=3); an ALT B under a primary A is dropped and comes
        // back as A's first shadow (kept=1). Under -N 1 only the kept=3 one stays.
        std::vector<std::vector<mem_seed_t> > s = {span(0, 0, 100), span(1, 0, 40)};
        FltParams f = defaults();
        f.max_chain_extend = 1;
        CHECK(same(run_filter(s, {1, 0}, f), {{0, 3, 100}, {1, 3, 40}}));
        CHECK(same(run_filter(s, {0, 1}, f), {{0, 3, 100}}));
        CHECK(same(filter_oracle(s, {0, 1}, f), {{0, 3, 100}}));
    }
    SUBCASE("-N: the N-th kept=1/2 chain and later ones are dropped; kept=3 chains are not counted") {
        // B and C overlap A and survive (-D 0.1); D is alone on the read (kept=3).
        FltParams f = defaults();
        f.drop_ratio = 0.1f;
        std::vector<std::vector<mem_seed_t> > s = {span(0, 0, 100), span(1, 0, 90), span(2, 0, 80),
                                                   span(3, 200, 270)};
        CHECK(same(run_filter(s, {0, 0, 0, 0}, f), {{0, 3, 100}, {1, 1, 90}, {2, 1, 80}, {3, 3, 70}}));
        f.max_chain_extend = 2;
        CHECK(same(run_filter(s, {0, 0, 0, 0}, f), {{0, 3, 100}, {1, 1, 90}, {3, 3, 70}}));
        f.max_chain_extend = 1;
        CHECK(same(run_filter(s, {0, 0, 0, 0}, f), {{0, 3, 100}, {3, 3, 70}}));
    }
}

TEST_CASE("mem_chain_flt: random chain sets equal the documented filter"
          * doctest::test_suite("unit/chain")) {
    std::mt19937 rng(777u);
    int checked = 0, dropped_some = 0;
    for (int n = 0; n < 400; ++n) {
        const int n_chn = 2 + (int)(rng() % 14);
        std::vector<std::vector<mem_seed_t> > seeds;
        std::vector<int> alt, weights;
        for (int c = 0; c < n_chn; ++c) {
            // One to three seeds in query/reference order on a 200 bp read.
            std::vector<mem_seed_t> s;
            int32_t q = (int32_t)(rng() % 150);
            int64_t r = 100000 * (int64_t)(c + 1) + q;
            for (int k = 1 + (int)(rng() % 3); k > 0 && q < 200; --k) {
                const int32_t len = std::min<int32_t>(19 + (int32_t)(rng() % 60), 200 - q);
                s.push_back(seed(r, q, len));
                q += len + (int32_t)(rng() % 10);
                r += len + (int64_t)(rng() % 10);
            }
            const int w = weight_oracle(s);
            if (std::find(weights.begin(), weights.end(), w) != weights.end()) continue;  // distinct weights
            weights.push_back(w);
            seeds.push_back(s);
            alt.push_back(rng() % 8 == 0 ? 1 : 0);
        }
        FltParams f = defaults();
        f.drop_ratio = (float)(rng() % 10) / 10.0f;
        f.mask_level = (float)(rng() % 10) / 10.0f;
        f.min_seed_len = 10 + (int)(rng() % 20);
        f.min_chain_weight = rng() % 3 == 0 ? (int)(rng() % 80) : 0;
        f.max_chain_gap = rng() % 4 == 0 ? 40 + (int)(rng() % 100) : 10000;
        f.max_chain_extend = rng() % 3 == 0 ? 1 + (int)(rng() % 4) : 1 << 30;
        CAPTURE(n);
        const std::vector<Kept> want = filter_oracle(seeds, alt, f);
        const std::vector<Kept> got = run_filter(seeds, alt, f);
        INFO("got " << show(got) << "; want " << show(want));
        CHECK(same(got, want));
        if (want.size() < seeds.size()) ++dropped_some;
        ++checked;
    }
    CHECK(checked == 400);
    CHECK(dropped_some >= 100);
}

TEST_CASE("mem_flt_chained_seeds: the seed SW screen keeps a seed scoring exactly the threshold"
          * doctest::test_suite("unit/flt_chained_seeds")) {
    // Two contigs indexed in a temporary directory: ctgA (3000 bases) then ctgB
    // (1500), so ctgB starts at forward coordinate 3000 and l_pac is 4500. The
    // 1000-base read is random except for:
    //   [0, 400)    ctgA 400-799                      seed L: qbeg 5, len 40
    //   [600, 615)  ctgA 2500-2514                    seed S: a 15-base seed alone
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
    CHECK((scB >= 45 && scB < 60));      // the window stops at ctgB's first base
    CHECK((scR >= 50 && scR < 65));
    OptGuard g;

    // -W values whose threshold round(1.1 * W) brackets seed S's score; both
    // leave the screen on (1.1 * W <= 0.05 * 1000).
    int w_keep = 0;
    for (int w = 1; w <= 45; ++w)
        if ((int)(1.1f * w + .499) <= scS) w_keep = w;
    REQUIRE(w_keep > 0);
    REQUIRE((int)(1.1f * (w_keep + 1) + .499) > scS);
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
    // On a 150-base read the screen is off (5.5 * ln 150 = 27.6 > 0.05 * 150):
    // nothing is rescored or dropped.
    {
        bseq1_t short_rec = rec;
        short_rec.l_seq = 150;
        std::vector<mem_seed_t> seeds = {sL, seed(2500, 120, 15)};
        mem_chain_t c = chain(seeds, 0);
        mem_flt_chained_seeds(g.opt, idx.bns(), idx.pac(), &short_rec, 1, &c);
        REQUIRE(c.n == 2);
        CHECK(c.seeds[0].score == 40);
        CHECK(c.seeds[1].score == 15);
    }
}
