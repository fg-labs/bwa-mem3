// test/unit/test_ungapped_fastpath.cpp
//
// Differential proof-check of the ungapped extension fast path (src/ungapped_ext.h).
// When ungapped_analyze returns FP_STATUS_HIT the aligner skips the banded-SW retry
// ladder and commits the region's score / qb / rb / truesc / a->w from the diagonal
// walk. The HIT envelope comment in ungapped_ext.h argues those equal what the
// full-width ladder would have committed. This file checks that argument against
// two independent kernels -- bwa-mem3's scalarBandedSWA and upstream bwa's
// ksw_extend2 -- over a parameter grid (-A/-B/-O/-E/-L/-d/-w) with random and
// adversarial pairs. It pins the two guards (E1: -w below 2, E2: z-drop) and the
// strict tie-break as load-bearing -- each has a pair on which the walk's fields
// differ from the kernels', so a HIT there would have changed output -- and checks
// that a dying diagonal (E3) needs no guard.
//
// Adversarial shapes:
//   * tie at the query end: a mismatch followed by exactly b/a matches, so the
//     diagonal returns to its record without exceeding it (the kernels keep the
//     EARLIER tied row; a `>=` walk took the later one and moved qb/rb under -L 0);
//   * mismatches first: the diagonal dips to h0 - b*X (E3, the zero-row break);
//   * a mismatch cluster after a record: drop b*X below the running max (guard E2,
//     z-drop under a small -d).
#include <cstdint>
#include <random>
#include <vector>

#include "doctest/doctest.h"
#include "ext_ladder.h"
#include "ungapped_ext.h"

using bwa_tests::ExtPair;
using bwa_tests::ExtResult;
using bwa_tests::ExtScoring;
using bwa_tests::LadderOutcome;

namespace {

struct FpOut { int status, score, qle, gscore, gtle, tb; };

FpOut analyze(const ExtScoring &sc, const ExtPair &p, int w) {
    FpOut o;
    o.status = ungapped_analyze(p.query.data(), p.target.data(), (int)p.query.size(), p.h0,
                                sc.a, sc.b, sc.o_min(), sc.e_min(), sc.x_threshold(w), w,
                                sc.zdrop, &o.score, &o.qle, &o.gscore, &o.gtle, &o.tb);
    return o;
}

// The diagonal walk's would-be HIT fields, with no guards: what a HIT would commit
// if the guard under test were absent.
FpOut walk_only(const ExtScoring &sc, const ExtPair &p) {
    const int N = (int)p.query.size();
    uint64_t mis[FP_MIS_NWORDS] = {0};
    for (int j = 0; j < N; ++j)
        if (p.query[j] != p.target[j]) mis[j >> 6] |= 1ULL << (j & 63);
    FpOut o; int cur;
    ungapped_walk_mis(mis, N, p.h0, sc.a, sc.b, &o.score, &o.qle, &cur);
    o.status = FP_STATUS_HIT; o.gscore = cur; o.gtle = N; o.tb = 0;
    return o;
}

uint8_t other_base(std::mt19937_64 &rng, uint8_t b) {
    return (uint8_t)((b + 1 + (int)(rng() % 3)) % 4);
}

// query: N random bases; target: the query on the diagonal with mismatches at
// `mis_pos`, then a random tail so tlen >= qlen.
ExtPair make_pair(std::mt19937_64 &rng, int N, int tail, const std::vector<int> &mis_pos, int h0) {
    ExtPair p;
    p.h0 = h0;
    p.query.resize(N);
    for (int j = 0; j < N; ++j) p.query[j] = (uint8_t)(rng() % 4);
    p.target = p.query;
    for (size_t k = 0; k < mis_pos.size(); ++k) {
        const int j = mis_pos[k];
        if (j >= 0 && j < N) p.target[j] = other_base(rng, p.query[j]);
    }
    for (int t = 0; t < tail; ++t) p.target.push_back((uint8_t)(rng() % 4));
    return p;
}

std::vector<int> random_positions(std::mt19937_64 &rng, int N, int k) {
    std::vector<int> v;
    for (int i = 0; i < k; ++i) v.push_back((int)(rng() % (uint64_t)N));
    return v;
}

// One scoring grid shared by the cases below. Mismatch penalties, gap costs, clip
// penalties, z-drops and band widths chosen so that x_threshold spans 0..6 and every
// guard has parameter cells on both sides of its boundary.
struct GridCell { int a, b, o_del, e_del, o_ins, e_ins, pen_clip, zdrop, w; };

std::vector<GridCell> grid() {
    static const int A[] = {1, 2};
    static const int B[] = {0, 1, 2, 4, 9};
    // {o_del, e_del, o_ins, e_ins}: symmetric, and two asymmetric schemes (the bound
    // reads o_min/e_min across the two gap directions).
    static const int OE[][4] = {{1, 1, 1, 1}, {6, 1, 6, 1}, {16, 1, 16, 1}, {20, 1, 20, 1},
                                {25, 2, 25, 2}, {30, 5, 30, 5}, {6, 1, 10, 2}, {20, 1, 6, 1}};
    static const int L[] = {0, 5};
    static const int Z[] = {0, 3, 4, 100};
    static const int W[] = {2, 100};
    std::vector<GridCell> g;
    for (size_t ia = 0; ia < 2; ++ia) for (size_t ib = 0; ib < 5; ++ib)
    for (size_t io = 0; io < 8; ++io) for (size_t il = 0; il < 2; ++il)
    for (size_t iz = 0; iz < 4; ++iz) for (size_t iw = 0; iw < 2; ++iw) {
        GridCell c; c.a = A[ia]; c.b = B[ib];
        c.o_del = OE[io][0]; c.e_del = OE[io][1]; c.o_ins = OE[io][2]; c.e_ins = OE[io][3];
        c.pen_clip = L[il]; c.zdrop = Z[iz]; c.w = W[iw];
        g.push_back(c);
    }
    return g;
}

ExtScoring scoring_of(const GridCell &c) {
    ExtScoring sc;
    sc.a = c.a; sc.b = c.b; sc.o_del = c.o_del; sc.e_del = c.e_del; sc.o_ins = c.o_ins; sc.e_ins = c.e_ins;
    sc.pen_clip = c.pen_clip; sc.zdrop = c.zdrop; sc.fill_mat();
    return sc;
}

// All pairs a grid cell is exercised with: random mismatch counts around
// x_threshold plus the three adversarial shapes, at several seed scores.
std::vector<ExtPair> pairs_for(std::mt19937_64 &rng, const ExtScoring &sc, int w) {
    static const int H0[] = {1, 5, 19, 200};
    const int X = sc.x_threshold(w) < 0 ? 0 : sc.x_threshold(w);
    std::vector<ExtPair> v;
    for (size_t ih = 0; ih < 4; ++ih) {
        const int h0 = H0[ih];
        // one short pair per seed score, and one pair up to FP_N_MAX (bitmap words
        // 2..7, the run walk across word boundaries) per cell
        for (int rep = 0; rep < (ih == 0 ? 2 : 1); ++rep) {
            const int N = rep == 0 ? 1 + (int)(rng() % 128) : 129 + (int)(rng() % (FP_N_MAX - 128));
            const int tail = (int)(rng() % 60);
            // random: k in [0, X+1]
            v.push_back(make_pair(rng, N, tail, random_positions(rng, N, (int)(rng() % (uint64_t)(X + 2))), h0));
            // tie at the end: one mismatch, then exactly b/a matches (b % a == 0)
            if (sc.b % sc.a == 0 && N > sc.b / sc.a) {
                std::vector<int> pos(1, N - 1 - sc.b / sc.a);
                v.push_back(make_pair(rng, N, tail, pos, h0));
            }
            // mismatches first (E3): X consecutive mismatches at the start
            { std::vector<int> pos; for (int k = 0; k < X && k < N; ++k) pos.push_back(k);
              v.push_back(make_pair(rng, N, tail, pos, h0)); }
            // cluster after a record (E2): X consecutive mismatches mid-way
            { std::vector<int> pos; const int s = N / 2;
              for (int k = 0; k < X && s + k < N; ++k) pos.push_back(s + k);
              v.push_back(make_pair(rng, N, tail, pos, h0)); }
        }
    }
    return v;
}

} // namespace

TEST_CASE("unit/ungapped_fastpath: a HIT commits exactly the full ladder's fields, on both kernels") {
    std::mt19937_64 rng(0x5EEDF00Dull);
    long n_hit = 0, n_tie_hit = 0, n_pairs = 0, n_x0 = 0, n_dead_hit = 0;
    const std::vector<GridCell> g = grid();
    for (size_t ic = 0; ic < g.size(); ++ic) {
        const ExtScoring sc = scoring_of(g[ic]);
        const int w = g[ic].w;
        const std::vector<ExtPair> ps = pairs_for(rng, sc, w);
        for (size_t ip = 0; ip < ps.size(); ++ip) {
            const ExtPair &p = ps[ip];
            ++n_pairs;
            const FpOut fp = analyze(sc, p, w);
            if (fp.status != FP_STATUS_HIT) continue;
            ++n_hit;
            if (sc.x_threshold(w) == 0) ++n_x0;
            // Both kernels, the exact 2-rung ladder, left (prev0 = -1) and right
            // (prev0 = h0) starting scores. Every raw field must match and the ladder
            // must accept at its first rung, so the region records a->w = opt->w.
            const bwa_tests::ExtKernel kernels[2] = {bwa_tests::run_mem3_scalar, bwa_tests::run_upstream_ksw};
            const int prev0s[2] = {-1, p.h0};
            for (int ik = 0; ik < 2; ++ik) for (int iv = 0; iv < 2; ++iv) {
                const LadderOutcome L = bwa_tests::run_ladder(kernels[ik], sc, p, w, 2, prev0s[iv]);
                CAPTURE(ic); CAPTURE(ip); CAPTURE(ik); CAPTURE(iv);
                CAPTURE(sc.a); CAPTURE(sc.b); CAPTURE(sc.o_del); CAPTURE(sc.e_del);
                CAPTURE(sc.pen_clip); CAPTURE(sc.zdrop); CAPTURE(w); CAPTURE(p.h0);
                CAPTURE((int)p.query.size());
                REQUIRE(L.rung == 0);
                REQUIRE(L.res.max_off == 0);
                REQUIRE(L.res.score == fp.score);
                REQUIRE(L.res.qle == fp.qle);
                REQUIRE(L.res.tle == fp.qle);
                if (fp.gscore > 0) {
                    REQUIRE(L.res.gscore == fp.gscore);
                    REQUIRE(L.res.gtle == fp.gtle);
                } else {
                    // the diagonal died (HIT envelope E3): the kernel broke on the
                    // same all-zero row and captured no positive query-end score, so
                    // both sides take branch A and gtle is not read.
                    ++n_dead_hit;
                    REQUIRE(L.res.gscore <= 0);
                }
                // and hence the committed region fields
                const bwa_tests::LeftCommit ck = bwa_tests::left_commit(L.res.score, L.res.qle, L.res.tle,
                                                                        L.res.gscore, L.res.gtle, sc.pen_clip, 1000);
                const bwa_tests::LeftCommit cf = bwa_tests::left_commit(fp.score, fp.qle, fp.qle,
                                                                        fp.gscore, fp.gtle, sc.pen_clip, 1000);
                REQUIRE(ck == cf);
                // a tie at the end is the shape the old >= walk got wrong: count the
                // HITs whose final diagonal value equals the record without a strict
                // record at the end (kernel qle < N with gscore == score).
                if (ik == 0 && iv == 0 && L.res.gscore == L.res.score && L.res.qle < (int)p.query.size()) ++n_tie_hit;
            }
        }
    }
    MESSAGE("pairs " << n_pairs << ", HITs " << n_hit << ", tie-at-end HITs " << n_tie_hit
            << ", HITs at x_threshold 0 " << n_x0 << ", HITs whose diagonal died " << n_dead_hit);
    REQUIRE(n_hit > 10000);
    REQUIRE(n_tie_hit > 100);
    REQUIRE(n_dead_hit > 100);
}

TEST_CASE("unit/ungapped_fastpath: strict tie-break -- the L0 pattern commits the kernels' qle") {
    // Left extension reversed prefix [X][m m m m] at default scoring, h0 = 145
    // (test/ungapped_hit_parity_test.sh's -L 0 read). The diagonal returns to h0
    // at step 5 without exceeding it: the kernels report qle = tle = 0, gscore =
    // h0. Under -L 0 that is branch A, so qle decides qb: a `>=` walk gave 5.
    ExtScoring sc; sc.a = 1; sc.b = 4; sc.o_del = sc.o_ins = 6; sc.e_del = sc.e_ins = 1;
    sc.pen_clip = 0; sc.zdrop = 100; sc.fill_mat();
    std::mt19937_64 rng(7);
    ExtPair p = make_pair(rng, 5, 20, std::vector<int>(1, 0), 145);
    const FpOut fp = analyze(sc, p, 100);
    REQUIRE(fp.status == FP_STATUS_HIT);
    const ExtResult k3 = bwa_tests::run_mem3_scalar(sc, p, 100);
    const ExtResult ku = bwa_tests::run_upstream_ksw(sc, p, 100);
    CHECK(k3 == ku);
    CHECK(k3.qle == 0);
    CHECK(fp.qle == k3.qle);
    CHECK(fp.score == k3.score);
    CHECK(fp.gscore == k3.gscore);
    CHECK(fp.gtle == k3.gtle);
    // and the -O 20 shape [X][m m m m][X][X]: tie at step 5, then two mismatches;
    // gscore = h0 - 8 <= h0 - 5 is branch A under the default -L 5, so qle decides
    // qb even at the default clip penalty.
    ExtScoring sc2 = sc; sc2.o_del = sc2.o_ins = 20; sc2.pen_clip = 5;
    int pos2[3] = {0, 5, 6};
    ExtPair p2 = make_pair(rng, 7, 20, std::vector<int>(pos2, pos2 + 3), 143);
    const FpOut fp2 = analyze(sc2, p2, 100);
    REQUIRE(fp2.status == FP_STATUS_HIT);
    const ExtResult k2 = bwa_tests::run_mem3_scalar(sc2, p2, 100);
    CHECK(k2 == bwa_tests::run_upstream_ksw(sc2, p2, 100));
    CHECK(k2.qle == 0);
    CHECK(fp2.qle == k2.qle);
    CHECK(fp2.gscore == k2.gscore);
    CHECK(fp2.gscore <= fp2.score - sc2.pen_clip);
}

TEST_CASE("unit/ungapped_fastpath: guard E2 (z-drop) rejects and is load-bearing; a dying diagonal (E3) needs none") {
    std::mt19937_64 rng(0xE2E3ull);
    long n_e2 = 0, n_e2_diverge = 0, n_dead = 0;
    const std::vector<GridCell> g = grid();
    for (size_t ic = 0; ic < g.size(); ++ic) {
        const ExtScoring sc = scoring_of(g[ic]);
        const int w = g[ic].w;
        const int X = sc.x_threshold(w);
        if (X < 1) continue;
        for (int rep = 0; rep < 4; ++rep) {
            const int N = 8 + (int)(rng() % 120);
            // E2: zdrop in (0, b*X) with a mismatch cluster right after the record.
            if (sc.zdrop > 0 && sc.b * X > sc.zdrop) {
                std::vector<int> pos; const int s = N / 2;
                for (int k = 0; k < X; ++k) pos.push_back(s + k);
                ExtPair p = make_pair(rng, N, 30, pos, 200);
                int actual = 0;
                for (int j = 0; j < N; ++j) actual += (p.query[j] != p.target[j]);
                if (sc.b * actual <= sc.zdrop) continue;   // cluster truncated by N: no drop
                const FpOut fp = analyze(sc, p, w);
                CAPTURE(ic); CAPTURE(rep);
                REQUIRE(fp.status == FP_STATUS_FALLBACK);
                ++n_e2;
                const FpOut wo = walk_only(sc, p);
                const ExtResult k = bwa_tests::run_mem3_scalar(sc, p, w);
                CHECK(k == bwa_tests::run_upstream_ksw(sc, p, w));
                if (!(k.score == wo.score && k.qle == wo.qle && k.gscore == wo.gscore)) ++n_e2_diverge;
            }
            // E3: h0 <= b*X with X mismatches first, so the diagonal reaches 0. Still
            // a HIT: the kernels break on the same row with no positive query-end
            // score, so every committed field agrees (branch A on both sides).
            if (sc.b > 0 && (sc.zdrop <= 0 || sc.b * X <= sc.zdrop)) {
                std::vector<int> pos; for (int k = 0; k < X; ++k) pos.push_back(k);
                const int h0 = 1 + (int)(rng() % (uint64_t)(sc.b * X));
                ExtPair p = make_pair(rng, N, 30, pos, h0);
                int actual = 0;
                for (int j = 0; j < N; ++j) actual += (p.query[j] != p.target[j]);
                if (h0 > sc.b * actual) continue;   // run truncated by N: the diagonal survives
                const FpOut fp = analyze(sc, p, w);
                CAPTURE(ic); CAPTURE(rep); CAPTURE(h0);
                REQUIRE(fp.status == FP_STATUS_HIT);
                REQUIRE(fp.gscore == 0);
                ++n_dead;
                const ExtResult k = bwa_tests::run_mem3_scalar(sc, p, w);
                CHECK(k == bwa_tests::run_upstream_ksw(sc, p, w));
                REQUIRE(k.score == fp.score);
                REQUIRE(k.qle == fp.qle);
                REQUIRE(k.tle == fp.qle);
                REQUIRE(k.gscore <= 0);
                const bwa_tests::LeftCommit ck = bwa_tests::left_commit(k.score, k.qle, k.tle, k.gscore, k.gtle, sc.pen_clip, 1000);
                const bwa_tests::LeftCommit cf = bwa_tests::left_commit(fp.score, fp.qle, fp.qle, fp.gscore, fp.gtle, sc.pen_clip, 1000);
                REQUIRE(ck == cf);
            }
        }
    }
    MESSAGE("E2 rejects " << n_e2 << " (kernel differs from the walk on " << n_e2_diverge
            << "); dying-diagonal HITs " << n_dead);
    REQUIRE(n_e2 > 100);
    REQUIRE(n_e2_diverge > 0);
    REQUIRE(n_dead > 100);
}

TEST_CASE("unit/ungapped_fastpath: guard E1 -- at -w 1 the ladder records a->w = 2, so the fast path is off") {
    ExtScoring sc; sc.a = 1; sc.b = 4; sc.o_del = sc.o_ins = 6; sc.e_del = sc.e_ins = 1;
    sc.pen_clip = 5; sc.zdrop = 100; sc.fill_mat();
    CHECK(sc.x_threshold(1) == -1);
    CHECK(sc.x_threshold(2) == 1);
    std::mt19937_64 rng(11);
    ExtPair p = make_pair(rng, 40, 10, std::vector<int>(), 30);   // perfect diagonal
    const LadderOutcome L1 = bwa_tests::run_ladder(bwa_tests::run_mem3_scalar, sc, p, 1, 2, -1);
    CHECK(L1.rung == 1);   // (1>>1)+(1>>2) == 0: max_off 0 is not < 0, so rung 0 is not accepted
    CHECK(L1.w == 2);
    const LadderOutcome L2 = bwa_tests::run_ladder(bwa_tests::run_mem3_scalar, sc, p, 2, 2, -1);
    CHECK(L2.rung == 0);
    CHECK(analyze(sc, p, 1).status == FP_STATUS_FALLBACK);
    CHECK(analyze(sc, p, 2).status == FP_STATUS_HIT);
}
