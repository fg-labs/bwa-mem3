// test/unit/test_mem_pair_pairing.cpp -- pair selection (mem_pair) and the
// paired-end MAPQ and proper-pair decision (mem_pair_resolve_batch_post).
//
// mem_pair picks, among every combination of the two ends' hits, the pair with
// the best pair score, and reports the runner-up score and how many pairs come
// close to it. mem_pair_resolve_batch_post turns that into the records' FLAG 0x2
// and MAPQ. Neither had a test: a pair-window bound made exclusive, a changed
// constant in the pair score, a tie broken differently, or a slip in the MAPQ
// recombination would change output on every paired-end run and pass the suite.
//
// mem_pair is checked against a brute-force oracle written from bwa's
// documented model (bwa 0.7.x bwamem_pair.c, mem_pair): every two hits of
// different ends, ordered by forward-strand position, whose orientation (strand
// of the upstream hit, strand of the downstream hit) is fitted and whose
// distance lies inside that orientation's [low, high], scores
// s1 + s2 + .721 * log(2 * erfc(|z| / sqrt 2)) * a + .499, truncated and floored
// at 0, with z = (distance - mean) / std; ties go to the larger
// hash_64((k << 32 | i) ^ id << 8), k < i the two hits' ranks in that position
// order. hash_64 is bwa's 64-bit integer hash (Thomas Wang's mix, bwa utils.h);
// the oracle carries its own copy, written out below, so that a change to the
// library's hash, which would move every tie-break in the output, is caught
// rather than followed. The oracle enumerates all pairs instead of sweeping, so
// a sweep that stops early or skips a hit disagrees with it.
//
// The MAPQ cases are worked by hand from the same model (bwa's mem_sam_pe):
// q_pe = 6.02 * (o - max(sub, s1 + s2 - U)) / a, less 4.343 * ln(n_sub + 1),
// capped at 60 and scaled by 1 - (frac_rep1 + frac_rep2) / 2; each end's MAPQ is
// raised to min(q_pe, its single-end MAPQ + 40) and capped at
// 6.02 * (score - csub) / a; a pair scoring no better than the two ends apart
// (o <= s1 + s2 - U) is not proper; under --meth the margin o - subo is also
// capped at o + 4a minus both ends' second-best scores. When no pair is
// resolved, 0x2 follows the unpaired rule (mem_proper_pair_extra_flag), checked
// at its window edges. Mate rescue is off in these cases
// (MEM_F_NO_RESCUE, i.e. -S), which is what lets them run without an index:
// the resolve then reads only the regions, the insert-size model and the
// contig offsets.

#include "doctest/doctest.h"
#include "bwamem.h"
#include "bntseq.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <random>
#include <vector>

int mem_pair(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac,
             const mem_pestat_t pes[4], bseq1_t s[2], mem_alnreg_v a[2], int id,
             int *sub, int *n_sub, int z[2], int n_pri[2]);

namespace {

struct OptGuard {
    mem_opt_t *o;
    OptGuard() : o(mem_opt_init()) {}
    ~OptGuard() { free(o); }
};

// A hand-built bntseq_t: `n` contigs of `len` bases laid end to end. Pairing
// reads only l_pac and each contig's offset.
struct Bns {
    bntseq_t b;
    std::vector<bntann1_t> anns;
    Bns(int n, int64_t len) : anns((size_t)n) {
        memset(&b, 0, sizeof(b));
        for (int i = 0; i < n; ++i) {
            memset(&anns[(size_t)i], 0, sizeof(bntann1_t));
            anns[(size_t)i].offset = (int64_t)i * len;
            anns[(size_t)i].len = (int32_t)len;
        }
        b.l_pac = (int64_t)n * len;
        b.n_seqs = n;
        b.anns = anns.data();
    }
    // A forward region at contig-relative [pos, pos + len), query [qb, qb + len).
    mem_alnreg_t fwd(int rid, int64_t pos, int score, int len = 100, int qb = 0) const {
        mem_alnreg_t a;
        memset(&a, 0, sizeof(a));
        a.rid = rid;
        a.rb = anns[(size_t)rid].offset + pos;
        a.re = a.rb + len;
        a.qb = qb;
        a.qe = qb + len;
        a.score = a.truesc = score;
        a.seedcov = len;
        a.secondary = a.secondary_all = -1;
        return a;
    }
    // A reverse region whose forward-strand span is contig-relative [pos, pos + len).
    mem_alnreg_t rev(int rid, int64_t pos, int score, int len = 100, int qb = 0) const {
        mem_alnreg_t a = fwd(rid, pos, score, len, qb);
        const int64_t fwd_end = anns[(size_t)rid].offset + pos + len;
        a.rb = 2 * b.l_pac - fwd_end;
        a.re = a.rb + len;
        return a;
    }
};

// The -I 300,30 model: FR fitted with mean 300, std 30 and window
// [(int)(300 - 120 + .499), (int)(300 + 120 + .499)] = [180, 420]; the other
// three orientations failed.
void pes_i300(mem_pestat_t pes[4]) {
    memset(pes, 0, 4 * sizeof(mem_pestat_t));
    for (int d = 0; d < 4; ++d) pes[d].failed = 1;
    pes[1].failed = 0;
    pes[1].avg = 300;
    pes[1].std = 30;
    pes[1].low = 180;
    pes[1].high = 420;
}

// ---- mem_pair oracle ---------------------------------------------------------

// bwa's hash_64 (Thomas Wang's 64-bit integer mix), the tie-break hash of the
// pair model, written out here rather than taken from the library.
uint64_t wang_hash64(uint64_t key) {
    key += ~(key << 32);
    key ^= (key >> 22);
    key += ~(key << 13);
    key ^= (key >> 8);
    key += (key << 3);
    key ^= (key >> 15);
    key += ~(key << 27);
    key ^= (key >> 31);
    return key;
}

struct Hit {
    uint64_t x;    // rid << 32 | contig-relative forward position (reverse: last base)
    uint64_t y;    // score << 32 | index << 2 | strand << 1 | end
};

struct PairResult {
    int ret, sub, n_sub;
    int z[2];
};

PairResult oracle_pair(const mem_opt_t *opt, const Bns &B, const mem_pestat_t pes[4],
                       const std::vector<mem_alnreg_t> regs[2], int id) {
    const int64_t l_pac = B.b.l_pac;
    std::vector<Hit> h;
    for (int r = 0; r < 2; ++r) {
        for (size_t i = 0; i < regs[r].size(); ++i) {
            const mem_alnreg_t &e = regs[r][i];
            const bool rev = e.rb >= l_pac;
            const int64_t f = rev ? 2 * l_pac - 1 - e.rb : e.rb;
            Hit t;
            t.x = (uint64_t)e.rid << 32 | (uint64_t)(f - B.anns[(size_t)e.rid].offset);
            t.y = (uint64_t)e.score << 32 | (uint64_t)i << 2 | (uint64_t)rev << 1 | (uint64_t)r;
            h.push_back(t);
        }
    }
    std::sort(h.begin(), h.end(), [](const Hit &p, const Hit &q) {
        return p.x < q.x || (p.x == q.x && p.y < q.y);
    });
    struct Cand {
        uint64_t key;    // q << 32 | tie-break hash
        size_t k, i;
    };
    std::vector<Cand> u;
    for (size_t i = 0; i < h.size(); ++i) {
        for (size_t k = 0; k < i; ++k) {
            if ((h[k].y & 1) == (h[i].y & 1)) continue;  // same end
            const int dir = (int)((h[k].y >> 1 & 1) << 1 | (h[i].y >> 1 & 1));
            if (pes[dir].failed) continue;
            const int64_t dist = (int64_t)(h[i].x - h[k].x);
            if (dist < pes[dir].low || dist > pes[dir].high) continue;
            const double ns = pes[dir].std > 0 ? (dist - pes[dir].avg) / pes[dir].std : 0.;
            int q = (int)((double)(h[i].y >> 32) + (double)(h[k].y >> 32)
                          + .721 * log(2. * erfc(fabs(ns) * M_SQRT1_2)) * opt->a + .499);
            if (q < 0) q = 0;
            const uint64_t y = (uint64_t)k << 32 | i;
            Cand c;
            c.key = (uint64_t)q << 32 | (wang_hash64(y ^ (uint64_t)(int64_t)(id << 8)) & 0xffffffffu);
            c.k = k;
            c.i = i;
            u.push_back(c);
        }
    }
    PairResult res;
    res.z[0] = res.z[1] = -7;
    if (u.empty()) {
        res.ret = res.sub = res.n_sub = 0;
        return res;
    }
    // Equal keys (a 32-bit hash collision) are ordered by (k, i), as the model's
    // (key, k << 32 | i) pairs are, so the order is total.
    std::sort(u.begin(), u.end(), [](const Cand &p, const Cand &q) {
        return p.key < q.key || (p.key == q.key && (p.k < q.k || (p.k == q.k && p.i < q.i)));
    });
    const Cand &best = u.back();
    res.ret = (int)(best.key >> 32);
    res.sub = u.size() > 1 ? (int)(u[u.size() - 2].key >> 32) : 0;
    const int tmp = std::max(opt->a + opt->b,
                             std::max(opt->o_del + opt->e_del, opt->o_ins + opt->e_ins));
    res.n_sub = 0;
    for (size_t j = 0; j + 1 < u.size(); ++j)
        if (res.sub - (int)(u[j].key >> 32) <= tmp) ++res.n_sub;
    res.z[h[best.i].y & 1] = (int)(h[best.i].y >> 2 & 0x3fffffff);
    res.z[h[best.k].y & 1] = (int)(h[best.k].y >> 2 & 0x3fffffff);
    return res;
}

PairResult run_pair(const mem_opt_t *opt, const Bns &B, const mem_pestat_t pes[4],
                    std::vector<mem_alnreg_t> regs[2], int id) {
    mem_alnreg_v a[2];
    for (int r = 0; r < 2; ++r) {
        a[r].n = a[r].m = regs[r].size();
        a[r].a = regs[r].empty() ? NULL : regs[r].data();
        a[r].capped_w = 0;
    }
    bseq1_t s[2];
    memset(s, 0, sizeof(s));
    int n_pri[2] = {(int)regs[0].size(), (int)regs[1].size()};
    PairResult res;
    res.z[0] = res.z[1] = -7;
    res.ret = mem_pair(opt, &B.b, NULL, pes, s, a, id, &res.sub, &res.n_sub, res.z, n_pri);
    return res;
}

void check_pair(const PairResult &got, const PairResult &want) {
    CHECK(got.ret == want.ret);
    CHECK(got.sub == want.sub);
    CHECK(got.n_sub == want.n_sub);
    CHECK(got.z[0] == want.z[0]);
    CHECK(got.z[1] == want.z[1]);
}

// ---- the resolve half ---------------------------------------------------------

struct Resolved {
    int paired, extra_flag;
    int z[2], q_se[2], n_pri[2];
};

// Runs mem_pair_resolve_batch_post with mate rescue off. `regs` is resolved in
// place (mem_mark_primary_se reorders it), so z indexes the reordered vectors.
Resolved resolve(mem_opt_t *opt, const Bns &B, const mem_pestat_t pes[4],
                 std::vector<mem_alnreg_t> regs[2], uint64_t id) {
    opt->flag |= MEM_F_NO_RESCUE;
    mem_alnreg_v a[2];
    for (int r = 0; r < 2; ++r) {
        a[r].n = a[r].m = regs[r].size();
        a[r].a = regs[r].data();
        a[r].capped_w = 0;
    }
    bseq1_t s[2];
    memset(s, 0, sizeof(s));
    mem_cache *mmc = (mem_cache *)calloc(1, sizeof(mem_cache));
    REQUIRE(mmc != NULL);
    int32_t gcnt = 0;
    Resolved out;
    memset(&out, 0, sizeof(out));
    out.z[0] = out.z[1] = -7;
    out.q_se[0] = out.q_se[1] = -7;
    mem_pair_resolve_batch_post(opt, &B.b, NULL, pes, id, s, a, NULL, mmc, gcnt, 0, out.n_pri,
                                out.z, out.q_se, &out.extra_flag, &out.paired);
    free(mmc);
    return out;
}

}  // namespace

TEST_CASE("mem_pair: random hit sets agree with an all-pairs oracle"
          * doctest::test_suite("unit/pair")) {
    OptGuard g;
    Bns B(3, 20000);
    std::mt19937_64 rng(2610);
    int paired = 0, with_sub = 0;
    for (int iter = 0; iter < 400; ++iter) {
        CAPTURE(iter);
        mem_pestat_t pes[4];
        for (int d = 0; d < 4; ++d) {
            pes[d].failed = (int)(rng() % 3 == 0);
            pes[d].avg = 150 + (double)(rng() % 300);
            pes[d].std = (double)(rng() % 60);          // 0 included: the std == 0 guard
            pes[d].low = 1 + (int)(rng() % 150);
            pes[d].high = pes[d].low + (int)(rng() % 500);
        }
        std::vector<mem_alnreg_t> regs[2];
        for (int r = 0; r < 2; ++r) {
            const int n = 1 + (int)(rng() % 6);
            for (int i = 0; i < n; ++i) {
                const int rid = (int)(rng() % 2);
                const int64_t pos = 1000 + (int64_t)(rng() % 900);
                const int score = 30 + (int)(rng() % 90);
                regs[r].push_back(rng() % 2 ? B.rev(rid, pos, score) : B.fwd(rid, pos, score));
            }
        }
        const int id = (int)(rng() % 100000);
        const PairResult want = oracle_pair(g.o, B, pes, regs, id);
        const PairResult got = run_pair(g.o, B, pes, regs, id);
        check_pair(got, want);
        paired += want.ret > 0;
        with_sub += want.sub > 0;
    }
    // Both outcomes must be common, or the comparison proves little.
    CHECK(paired > 100);
    CHECK(with_sub > 50);
}

TEST_CASE("mem_pair: the pair window is inclusive at both ends"
          * doctest::test_suite("unit/pair")) {
    // Read 1 forward at [1000, 1100); read 2 reverse with its last base at
    // 1000 + d, so the FR distance is d. The window is [180, 420].
    OptGuard g;
    Bns B(1, 100000);
    mem_pestat_t pes[4];
    pes_i300(pes);
    const struct {
        int64_t d;
        bool pairs;
    } cases[] = {{179, false}, {180, true}, {420, true}, {421, false}};
    for (const auto &c : cases) {
        CAPTURE(c.d);
        std::vector<mem_alnreg_t> regs[2];
        regs[0].push_back(B.fwd(0, 1000, 100));
        regs[1].push_back(B.rev(0, 1000 + c.d - 99, 100));
        const PairResult got = run_pair(g.o, B, pes, regs, 1);
        CHECK((got.ret > 0) == c.pairs);
        if (c.pairs) {
            CHECK(got.z[0] == 0);
            CHECK(got.z[1] == 0);
        } else {
            CHECK(got.z[0] == -7);   // z is left alone when nothing pairs
        }
    }
}

TEST_CASE("mem_proper_pair_extra_flag: the unpaired rule's window is inclusive"
          * doctest::test_suite("unit/pair")) {
    // When no pair is resolved (an end with two primary hits, or no pair in the
    // window), 0x2 is set iff the two ends' emitted regions lie at a distance
    // inside their orientation's window: low <= d <= high (bwa's mem_sam_pe).
    // Read 1 forward at [1000, 1100), read 2 reverse with its last base at
    // 1000 + d: FR at distance d, window [180, 420].
    OptGuard g;
    Bns B(1, 100000);
    mem_pestat_t pes[4];
    pes_i300(pes);
    const struct {
        int64_t d;
        int flag;
    } cases[] = {{179, 0}, {180, 2}, {300, 2}, {420, 2}, {421, 0}};
    for (const auto &c : cases) {
        CAPTURE(c.d);
        mem_alnreg_t r0 = B.fwd(0, 1000, 100), r1 = B.rev(0, 1000 + c.d - 99, 100);
        mem_alnreg_v a[2];
        a[0].n = a[0].m = 1;
        a[0].a = &r0;
        a[1].n = a[1].m = 1;
        a[1].a = &r1;
        const int which[2] = {0, 0};
        CHECK(mem_proper_pair_extra_flag(g.o, B.b.l_pac, a, which, pes) == c.flag);
        // A failed orientation is never proper.
        pes[1].failed = 1;
        CHECK(mem_proper_pair_extra_flag(g.o, B.b.l_pac, a, which, pes) == 0);
        pes[1].failed = 0;
    }
}

TEST_CASE("mem_pair: the pair score's insert-size term, by hand"
          * doctest::test_suite("unit/pair")) {
    // Two 100-point ends, mean 300, std 30 (or 0). q = (int)(200 + .721 *
    // ln(2 erfc(|z| / sqrt 2)) + .499):
    //   z = 0:   .721 * ln 2 = .49976; 200.99876 -> 200
    //   z = 1:   .721 * ln(0.634697) = -.32774; 200.17126 -> 200
    //   z = 2:   .721 * ln(0.0455003) = -2.22790; 198.27110 -> 198
    //   z = 4:   .721 * ln(1.26685e-4) = -6.47004; 194.02896 -> 194
    //   std = 0: z is taken as 0 (no division by zero) -> 200
    OptGuard g;
    Bns B(1, 100000);
    const struct {
        int64_t d;
        double std_;
        int q;
    } cases[] = {{300, 30, 200}, {330, 30, 200}, {360, 30, 198}, {420, 30, 194}, {300, 0, 200}};
    for (const auto &c : cases) {
        CAPTURE(c.d);
        CAPTURE(c.std_);
        mem_pestat_t pes[4];
        pes_i300(pes);
        pes[1].std = c.std_;
        std::vector<mem_alnreg_t> regs[2];
        regs[0].push_back(B.fwd(0, 1000, 100));
        regs[1].push_back(B.rev(0, 1000 + c.d - 99, 100));
        const PairResult got = run_pair(g.o, B, pes, regs, 1);
        CHECK(got.ret == c.q);
        CHECK(got.sub == 0);
        CHECK(got.n_sub == 0);
    }
}

TEST_CASE("mem_pair: runner-up score and the count of near runners-up"
          * doctest::test_suite("unit/pair")) {
    // One mate (100) and three hits of the other end (100, 96 and 90 or 88), all
    // pairing at distances 299, 294 and 289 under std 1000, where the insert term
    // is .721 * ln(2 erfc(< .0011)) + .499 = .996 and so adds nothing: pair scores
    // 200, 196 and 190 (or 188). sub = 196; n_sub counts the pairs other than the
    // best within max(a + b, o + e) = 7 of sub: 196 and 190 (6 <= 7) give 2, 196
    // and 188 (8 > 7) give 1.
    OptGuard g;
    Bns B(1, 100000);
    mem_pestat_t pes[4];
    pes_i300(pes);
    pes[1].std = 1000;
    pes[1].low = 1;
    pes[1].high = 4300;
    for (int third = 90; third >= 88; third -= 2) {
        CAPTURE(third);
        std::vector<mem_alnreg_t> regs[2];
        regs[0].push_back(B.fwd(0, 1000, 100));
        regs[0].push_back(B.fwd(0, 1005, 96));
        regs[0].push_back(B.fwd(0, 1010, third));
        regs[1].push_back(B.rev(0, 1200, 100));
        const PairResult got = run_pair(g.o, B, pes, regs, 3);
        CHECK(got.ret == 200);
        CHECK(got.sub == 196);
        CHECK(got.n_sub == (third == 90 ? 2 : 1));
        CHECK(got.z[0] == 0);
        CHECK(got.z[1] == 0);
    }
}

TEST_CASE("mem_pair_resolve_batch_post: a unique FR pair is proper with MAPQ 60"
          * doctest::test_suite("unit/pair")) {
    // Ends at [1000, 1100) forward and [1200, 1300) reverse: distance 299,
    // z = -1/30, o = (int)(200 + .4804 + .499) = 200. Unpaired score
    // 200 - U = 183; q_pe = (int)(6.02 * 17 + .499) = 102, capped at 60. Each end
    // alone: sub = 19 (min_seed_len * a), MAPQ (int)(6.02 * 81 * (3/ln 100)^2 +
    // .499) = 207 -> 60.
    OptGuard g;
    Bns B(3, 100000);
    mem_pestat_t pes[4];
    pes_i300(pes);
    std::vector<mem_alnreg_t> regs[2];
    regs[0].push_back(B.fwd(0, 1000, 100));
    regs[1].push_back(B.rev(0, 1200, 100));
    const Resolved r = resolve(g.o, B, pes, regs, 11);
    CHECK(r.paired == 1);
    CHECK(r.extra_flag == 3);
    CHECK(r.n_pri[0] == 1);
    CHECK(r.n_pri[1] == 1);
    CHECK(r.z[0] == 0);
    CHECK(r.z[1] == 0);
    CHECK(r.q_se[0] == 60);
    CHECK(r.q_se[1] == 60);
}

TEST_CASE("mem_pair_resolve_batch_post: a repeat end placed by its mate gets min(q_pe, 40)"
          * doctest::test_suite("unit/pair")) {
    // End 0 has three equal 100-point hits, on contigs 0, 1 and 2; its mate is
    // unique on contig 1, 299 bp downstream of the contig-1 copy. Only that copy
    // pairs: o = 200, no runner-up, so subo = 200 - U and q_pe = (int)(6.02 * U +
    // .499): 102 -> 60 at U 17, (int)(24.08 + .499) = 24 at U 4, (int)(12.04 +
    // .499) = 12 at U 2. End 0 alone has MAPQ 0 (an equal second hit), so its
    // MAPQ is min(q_pe, 0 + 40): 40, 24, 12. The mate keeps 60.
    const struct {
        int U, q0;
    } cases[] = {{17, 40}, {4, 24}, {2, 12}};
    int switched = 0, kept = 0;
    for (const auto &c : cases) {
        for (uint64_t id = 0; id < 12; ++id) {
            CAPTURE(c.U);
            CAPTURE(id);
            OptGuard g;
            g.o->pen_unpaired = c.U;
            Bns B(3, 100000);
            mem_pestat_t pes[4];
            pes_i300(pes);
            std::vector<mem_alnreg_t> regs[2];
            regs[0].push_back(B.fwd(0, 1000, 100));
            regs[0].push_back(B.fwd(1, 1000, 100));
            regs[0].push_back(B.fwd(2, 1000, 100));
            regs[1].push_back(B.rev(1, 1200, 100));
            const Resolved r = resolve(g.o, B, pes, regs, id);
            REQUIRE(r.paired == 1);
            CHECK(r.extra_flag == 3);
            REQUIRE(r.z[0] >= 0);
            REQUIRE(r.z[0] < 3);
            CHECK(regs[0][(size_t)r.z[0]].rid == 1);
            CHECK(r.q_se[0] == c.q0);
            CHECK(r.q_se[1] == 60);
            // The chosen copy heads the secondary_all tree the XA tag is built
            // from: it has none, and the other two point at it.
            CHECK(regs[0][(size_t)r.z[0]].secondary_all == -1);
            for (int j = 0; j < 3; ++j)
                if (j != r.z[0]) CHECK(regs[0][(size_t)j].secondary_all == r.z[0]);
            // Whether the copy was already the single-end primary depends on the
            // read id's hash; count both so the switch is known to have run.
            if (regs[0][0].rid == 1) ++kept;
            else ++switched;
        }
    }
    CHECK(switched > 0);
    CHECK(kept > 0);
}

TEST_CASE("mem_pair_resolve_batch_post: near runners-up lower q_pe"
          * doctest::test_suite("unit/pair")) {
    // End 0: A at [1000, 1100) scoring 100 and B at [1010, 1110) scoring 96 over
    // the same read bases (B is A's secondary: A.sub = 96, A.sub_n = 1 since
    // 100 - 96 <= 7). The mate at [1200, 1300) reverse pairs with both: A at 299
    // (o = 200), B at 289 (z = -11/30: .721 * ln(2 erfc(.25927)) = .2567,
    // (int)(196.7557) = 196). subo = max(196, 183) = 196; q_pe = (int)(6.02 * 4 +
    // .499) = 24, less (int)(4.343 * ln 2 + .499) = 3 for n_sub = 1: 21.
    // A alone: (int)(6.02 * 4 * (3/ln 100)^2 + .499) = 10, less 3 for sub_n: 7,
    // raised to min(21, 47) = 21. The mate stays at 60.
    OptGuard g;
    Bns B(1, 100000);
    mem_pestat_t pes[4];
    pes_i300(pes);
    std::vector<mem_alnreg_t> regs[2];
    regs[0].push_back(B.fwd(0, 1010, 96));
    regs[0].push_back(B.fwd(0, 1000, 100));
    regs[1].push_back(B.rev(0, 1200, 100));
    const Resolved r = resolve(g.o, B, pes, regs, 5);
    REQUIRE(r.paired == 1);
    CHECK(r.extra_flag == 3);
    CHECK(r.z[0] == 0);
    CHECK(regs[0][0].score == 100);
    CHECK(r.q_se[0] == 21);
    CHECK(r.q_se[1] == 60);
}

TEST_CASE("mem_pair_resolve_batch_post: repetitive seeds and a tandem hit scale MAPQ"
          * doctest::test_suite("unit/pair")) {
    OptGuard g;
    Bns B(1, 100000);
    mem_pestat_t pes[4];
    pes_i300(pes);
    SUBCASE("frac_rep 0.5 on end 0") {
        // q_pe = (int)(60 * (1 - .5 * (.5 + 0)) + .499) = 45. End 0 alone:
        // (int)(60 * .5 + .499) = 30, raised to min(45, 70) = 45. Mate: 60.
        std::vector<mem_alnreg_t> regs[2];
        regs[0].push_back(B.fwd(0, 1000, 100));
        regs[0][0].frac_rep = 0.5f;
        regs[1].push_back(B.rev(0, 1200, 100));
        const Resolved r = resolve(g.o, B, pes, regs, 7);
        REQUIRE(r.paired == 1);
        CHECK(r.q_se[0] == 45);
        CHECK(r.q_se[1] == 60);
    }
    SUBCASE("csub 95 on end 0") {
        // End 0 alone: sub = max(19, 95) = 95, (int)(6.02 * 5 * .424376 + .499) =
        // 13, raised to min(60, 53) = 53, then capped at the tandem score
        // (int)(6.02 * (100 - 95) + .499) = 30. Mate: 60.
        std::vector<mem_alnreg_t> regs[2];
        regs[0].push_back(B.fwd(0, 1000, 100));
        regs[0][0].csub = 95;
        regs[1].push_back(B.rev(0, 1200, 100));
        const Resolved r = resolve(g.o, B, pes, regs, 7);
        REQUIRE(r.paired == 1);
        CHECK(r.q_se[0] == 30);
        CHECK(r.q_se[1] == 60);
    }
}

TEST_CASE("mem_pair_resolve_batch_post: a pair no better than its ends apart is not proper"
          * doctest::test_suite("unit/pair")) {
    // The mate's last base is 420 bp past end 0's first (z = 4): o = 194 (see
    // the pair-score case). At U 17 the ends apart score 183 < 194, so the pair
    // is proper and q_pe = (int)(6.02 * 11 + .499) = 66 -> 60; at U 7, 193 < 194,
    // still proper (q_pe (int)(6.02 + .499) = 6, below each end's own 60). At U
    // 6 the ends apart score 194, equal to o, and at U 3 197 > 194: the pair is
    // preferred only when strictly better, so no 0x2, and both ends keep their
    // own MAPQ (60) at index 0.
    const struct {
        int U, flag;
    } cases[] = {{17, 3}, {7, 3}, {6, 1}, {3, 1}};
    for (const auto &c : cases) {
        CAPTURE(c.U);
        OptGuard g;
        g.o->pen_unpaired = c.U;
        Bns B(1, 100000);
        mem_pestat_t pes[4];
        pes_i300(pes);
        std::vector<mem_alnreg_t> regs[2];
        regs[0].push_back(B.fwd(0, 1000, 100));
        regs[1].push_back(B.rev(0, 1321, 100));
        const Resolved r = resolve(g.o, B, pes, regs, 9);
        REQUIRE(r.paired == 1);
        CHECK(r.extra_flag == c.flag);
        CHECK(r.z[0] == 0);
        CHECK(r.z[1] == 0);
        CHECK(r.q_se[0] == 60);
        CHECK(r.q_se[1] == 60);
    }
}

TEST_CASE("mem_pair_resolve_batch_post: no pair, a split end and -P leave the ends unpaired"
          * doctest::test_suite("unit/pair")) {
    Bns B(3, 100000);
    mem_pestat_t pes[4];
    pes_i300(pes);
    SUBCASE("mate on another contig") {
        OptGuard g;
        std::vector<mem_alnreg_t> regs[2];
        regs[0].push_back(B.fwd(0, 1000, 100));
        regs[1].push_back(B.rev(1, 1200, 100));
        const Resolved r = resolve(g.o, B, pes, regs, 1);
        CHECK(r.paired == 0);
        CHECK(r.extra_flag == 1);
    }
    SUBCASE("an end with a second, non-overlapping primary hit (a split read)") {
        // End 0: [1000, 1060) for read bases 0-59 (60) and contig 2 for bases
        // 60-99 (40 >= -T 30). The pair with the mate exists (o = 160), but an end
        // with two primary hits is left to the single-end path.
        OptGuard g;
        std::vector<mem_alnreg_t> regs[2];
        regs[0].push_back(B.fwd(0, 1000, 60, 60, 0));
        regs[0].push_back(B.fwd(2, 5000, 40, 40, 60));
        regs[1].push_back(B.rev(0, 1200, 100));
        const Resolved r = resolve(g.o, B, pes, regs, 1);
        CHECK(r.n_pri[0] == 2);
        CHECK(r.paired == 0);
        CHECK(r.extra_flag == 1);
    }
    SUBCASE("-P (MEM_F_NOPAIRING) skips pairing after marking primaries") {
        OptGuard g;
        g.o->flag |= MEM_F_NOPAIRING;
        std::vector<mem_alnreg_t> regs[2];
        regs[0].push_back(B.fwd(0, 1000, 100));
        regs[1].push_back(B.rev(0, 1200, 100));
        const Resolved r = resolve(g.o, B, pes, regs, 1);
        CHECK(r.paired == 0);
        CHECK(r.extra_flag == 1);
        CHECK(r.n_pri[0] == 1);
        CHECK(r.n_pri[1] == 1);
    }
}

TEST_CASE("mem_pair_resolve_batch_post: --meth caps q_pe by both ends' second hits"
          * doctest::test_suite("unit/pair")) {
    // Each end has a 99-point second hit over the same read bases on contig 2,
    // too far apart to pair; only the contig-0 hits pair (o = 200, subo = 183,
    // q_pe 60). End MAPQ alone: (int)(6.02 * 1 * .424376 + .499) = 3, less 3 for
    // sub_n = 1: 0. Without --meth each end is raised to min(60, 40) = 40. Under
    // --meth the pair margin is capped at o + 4a - (99 + 99) = 6:
    // q_pe = (int)(36.12 + .499) = 36, and each end gets min(36, 40) = 36.
    for (int meth = 0; meth <= 1; ++meth) {
        CAPTURE(meth);
        OptGuard g;
        g.o->meth_mode = meth;
        Bns B(3, 100000);
        mem_pestat_t pes[4];
        pes_i300(pes);
        std::vector<mem_alnreg_t> regs[2];
        regs[0].push_back(B.fwd(0, 1000, 100));
        regs[0].push_back(B.fwd(2, 1000, 99));
        regs[1].push_back(B.rev(0, 1200, 100));
        regs[1].push_back(B.rev(2, 50000, 99));
        const Resolved r = resolve(g.o, B, pes, regs, 4);
        REQUIRE(r.paired == 1);
        CHECK(r.extra_flag == 3);
        CHECK(r.q_se[0] == (meth ? 36 : 40));
        CHECK(r.q_se[1] == (meth ? 36 : 40));
    }
}
