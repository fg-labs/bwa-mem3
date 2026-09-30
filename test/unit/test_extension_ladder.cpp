// test/unit/test_extension_ladder.cpp
//
// The seed-extension retry ladder's rung count (bwamem.cpp MAX_BAND_TRY). Upstream
// bwa and bwa-mem2 run w = opt->w, then 2*opt->w, and keep the 2w result whatever
// the per-rung stop test says there. bwa-mem3 carried a 4-rung ladder
// [w, 2w, 4w, 8w] from #58 (a leftover of an initial-band experiment reverted in
// the same change), and a pair that fails the stop test at 2w then ran on: a
// different alignment when a gap lies beyond offset 2w, and always a wider a->w
// (the contained-seed purge, the mem_reg2aln CIGAR band, the mem_patch_reg band).
//
// Two facts about a model of the exact ladder (ext_ladder.h run_ladder, with the
// rung count passed in) are checked here, differentially over both scalar kernels
// (bwa-mem3's scalarBandedSWA and upstream ksw_extend2). The production rung
// count itself is pinned end to end by test/ladder_rungs_test.sh.
//
//  1. When a third rung can run. Failing the stop test at 2w needs a
//     record-setting cell at diagonal offset >= T = (2w>>1)+(2w>>2), which costs a
//     gap of >= o_min + e_min*T and pays back at most a per aligned column, so it
//     needs a*L > o_min + e_min*T query bases: at defaults L >= 157. Below that
//     bound no pair, random or adversarial, fails the 2w test.
//
//  2. Past the bound a third rung changes what is committed: a right extension
//     with insertions at offsets 80, 150 and 250 (the read of
//     test/ladder_rungs_test.sh) stops at 2w = 200 on the 2-rung ladder and at
//     4w = 400 on the 4-rung one, with a higher score and a longer alignment.
#include <cstdint>
#include <random>
#include <vector>

#include "doctest/doctest.h"
#include "ext_ladder.h"

using bwa_tests::ExtPair;
using bwa_tests::ExtScoring;
using bwa_tests::LadderOutcome;

namespace {

ExtScoring scoring(int a, int b, int o, int e, int pen_clip, int zdrop) {
    ExtScoring sc;
    sc.a = a; sc.b = b; sc.o_del = sc.o_ins = o; sc.e_del = sc.e_ins = e;
    sc.pen_clip = pen_clip; sc.zdrop = zdrop; sc.fill_mat();
    return sc;
}

void append_random(std::mt19937_64 &rng, std::vector<uint8_t> &v, int n) {
    for (int i = 0; i < n; ++i) v.push_back((uint8_t)(rng() % 4));
}

// Copy of `seg` with a deterministic substitution every `every` bases: kills exact
// k-mers of length >= every so a read built from it seeds only on its exact part.
std::vector<uint8_t> mutated(const std::vector<uint8_t> &seg, int every) {
    std::vector<uint8_t> m = seg;
    for (size_t p = (size_t)every - 1; p < m.size(); p += (size_t)every) m[p] = (uint8_t)((m[p] + 1) % 4);
    return m;
}

// The query-side length bound for a third rung: the largest L with a*L <= o_min + e_min*T.
int third_rung_length_bound(const ExtScoring &sc, int w) {
    const int w2 = 2 * w;
    const int T = (w2 >> 1) + (w2 >> 2);
    return (sc.o_min() + sc.e_min() * T) / sc.a;
}

// Random pair shapes that push the alignment off the diagonal: a single insertion
// or deletion of `gap` bases after `pre` diagonal bases, then matches to the end.
ExtPair gapped_pair(std::mt19937_64 &rng, int qlen, int pre, int gap, bool insertion, int h0) {
    ExtPair p; p.h0 = h0;
    std::vector<uint8_t> common;
    append_random(rng, common, qlen + gap + 8);
    if (insertion) {
        // query: pre common bases, gap random bases, then the rest of common
        // (pre <= qlen - gap keeps the tail range forward and the query qlen long)
        if (pre > qlen - gap) pre = qlen - gap;
        p.query.assign(common.begin(), common.begin() + pre);
        append_random(rng, p.query, gap);
        p.query.insert(p.query.end(), common.begin() + pre, common.begin() + (qlen - gap));
        p.target.assign(common.begin(), common.begin() + qlen);       // no gap
    } else {
        // target: pre common bases, gap random bases, then the rest of common
        p.query.assign(common.begin(), common.begin() + qlen);
        p.target.assign(common.begin(), common.begin() + pre);
        append_random(rng, p.target, gap);
        p.target.insert(p.target.end(), common.begin() + pre, common.end());
    }
    while ((int)p.target.size() < (int)p.query.size()) p.target.push_back((uint8_t)(rng() % 4));
    return p;
}

} // namespace

TEST_CASE("unit/extension_ladder: no pair within the length bound fails the 2w stop test") {
    std::mt19937_64 rng(0x1ADDE12ull);
    static const int A[] = {1, 2}, B[] = {1, 4, 9};
    static const int OE[][2] = {{1, 1}, {6, 1}, {20, 1}, {25, 2}};
    static const int W[] = {2, 10, 20, 100};
    long n_pairs = 0, n_second_rung = 0;
    for (size_t ia = 0; ia < 2; ++ia) for (size_t ib = 0; ib < 3; ++ib)
    for (size_t io = 0; io < 4; ++io) for (size_t iw = 0; iw < 4; ++iw) {
        const ExtScoring sc = scoring(A[ia], B[ib], OE[io][0], OE[io][1], 5, 100);
        const int w = W[iw];
        const int L = third_rung_length_bound(sc, w);
        const int T = ((2 * w) >> 1) + ((2 * w) >> 2);
        for (int rep = 0; rep < 18; ++rep) {
            // qlen at or below the bound (a few far below). Gap shapes: at or just
            // beyond T (the offset a third rung would need); in [3w/4, w] right after
            // the start, so rung 0 fails its max_off test and the pair reaches the 2w
            // rung, where the same record sits below T; or random up to 2w.
            const int qlen = (rep < 12) ? L - (int)(rng() % 4) : 1 + (int)(rng() % (uint64_t)L);
            if (qlen < 2) continue;
            const int T0 = (w >> 1) + (w >> 2);
            int gap, pre;
            if (rep % 3 == 0) { gap = T + (int)(rng() % 3); pre = (int)(rng() % (uint64_t)qlen); }
            else if (rep % 3 == 1) { gap = T0 + (int)(rng() % (uint64_t)(w - T0 + 1)); pre = (int)(rng() % 3); }
            else { gap = 1 + (int)(rng() % (uint64_t)(2 * w)); pre = (int)(rng() % (uint64_t)qlen); }
            const bool ins = (rng() & 1) != 0;
            const int h0 = 1 + (int)(rng() % 300);
            ExtPair p = gapped_pair(rng, qlen, pre, gap, ins && gap < qlen, h0);
            REQUIRE((int)p.query.size() == qlen);
            ++n_pairs;
            const bwa_tests::ExtKernel kernels[2] = {bwa_tests::run_mem3_scalar, bwa_tests::run_upstream_ksw};
            for (int ik = 0; ik < 2; ++ik) {
                // The 4-rung ladder must stop by rung 1 (w or 2w): the stop test at 2w
                // never fails inside the bound.
                const LadderOutcome L4 = bwa_tests::run_ladder(kernels[ik], sc, p, w, 4, -1);
                CAPTURE(ia); CAPTURE(ib); CAPTURE(io); CAPTURE(iw); CAPTURE(rep); CAPTURE(qlen); CAPTURE(gap);
                REQUIRE(L4.rung <= 1);
                if (L4.rung == 1) {
                    ++n_second_rung;
                    // and it stopped because of the test, not because rung 1 was last
                    const bwa_tests::ExtResult r2 = kernels[ik](sc, p, 2 * w);
                    REQUIRE(r2.max_off < T);   // the max_off clause alone suffices inside the bound
                }
            }
        }
    }
    MESSAGE("pairs " << n_pairs << ", of which reached the 2w rung " << n_second_rung);
    REQUIRE(n_pairs > 1000);
    REQUIRE(n_second_rung > 100);
}

TEST_CASE("unit/extension_ladder: past the bound a third rung changes the committed alignment") {
    // The right extension of test/ladder_rungs_test.sh's read at default scoring:
    // h0 = 210 (the 210 bp exact seed), query = 80 inserted bases, Q2a (140 bp, a
    // substitution every 18), 70 inserted, Q2b (130 bp), 100 inserted, Q3 (170 bp);
    // target = Q2a Q2b Q3 and a 450 bp tail, so the target spans qlen + 200 like
    // the reference window stage_seed_extension passes (cal_max_gap allows 200). Each segment nets +19 over its insertion, at
    // diagonal offsets 80, 150 and 250.
    std::mt19937_64 rng(0x900ull);
    std::vector<uint8_t> q2a, q2b, q3, tail;
    append_random(rng, q2a, 140); append_random(rng, q2b, 130); append_random(rng, q3, 170);
    append_random(rng, tail, 450);
    ExtPair p; p.h0 = 210;
    append_random(rng, p.query, 80);
    { std::vector<uint8_t> m = mutated(q2a, 18); p.query.insert(p.query.end(), m.begin(), m.end()); }
    append_random(rng, p.query, 70);
    { std::vector<uint8_t> m = mutated(q2b, 18); p.query.insert(p.query.end(), m.begin(), m.end()); }
    append_random(rng, p.query, 100);
    { std::vector<uint8_t> m = mutated(q3, 18); p.query.insert(p.query.end(), m.begin(), m.end()); }
    p.target = q2a;
    p.target.insert(p.target.end(), q2b.begin(), q2b.end());
    p.target.insert(p.target.end(), q3.begin(), q3.end());
    p.target.insert(p.target.end(), tail.begin(), tail.end());
    REQUIRE(p.query.size() == 690u);

    const ExtScoring sc = scoring(1, 4, 6, 1, 5, 100);
    REQUIRE((int)p.query.size() > third_rung_length_bound(sc, 100));
    const bwa_tests::ExtKernel kernels[2] = {bwa_tests::run_mem3_scalar, bwa_tests::run_upstream_ksw};
    for (int ik = 0; ik < 2; ++ik) {
        CAPTURE(ik);
        // right extension: the ladder starts from the left score, prev0 = h0
        const LadderOutcome L2 = bwa_tests::run_ladder(kernels[ik], sc, p, 100, 2, p.h0);
        const LadderOutcome L4 = bwa_tests::run_ladder(kernels[ik], sc, p, 100, 4, p.h0);
        // upstream: rung 0 finds Q2a (offset 80 >= 75, score changed), rung 1 finds
        // Q2b (offset 150 >= 150, score changed) and is the last rung.
        CHECK(L2.rung == 1);
        CHECK(L2.w == 200);
        CHECK(L2.res.score == 210 + 19 + 19);
        CHECK(L2.res.max_off == 150);
        // Q2b's last substitution (position 126) is followed by 4 matches, which only
        // tie the record; the kernels keep the earlier row, so the alignment ends 5
        // bases short of Q2b's end (bwa's 125M275S).
        CHECK(L2.res.qle == 80 + 140 + 70 + 130 - 5);
        // the former 4-rung ladder: rung 2 (w = 400) reaches Q3 at offset 250 and
        // is accepted there (max_off 250 < 300).
        CHECK(L4.rung == 2);
        CHECK(L4.w == 400);
        CHECK(L4.res.score == 210 + 19 + 19 + 19);
        CHECK(L4.res.qle == 690);
        CHECK(L4.res.max_off == 250);
        // so score, qe/re and a->w all differ from upstream's
        CHECK(L4.res.score != L2.res.score);
        CHECK(L4.res.qle != L2.res.qle);
        CHECK(L4.w != L2.w);
    }
    // and both kernels agree with each other rung for rung
    for (int i = 0; i < 4; ++i) {
        CAPTURE(i);
        CHECK(bwa_tests::run_mem3_scalar(sc, p, 100 << i) == bwa_tests::run_upstream_ksw(sc, p, 100 << i));
    }
}
