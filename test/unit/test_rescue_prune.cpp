// test/unit/test_rescue_prune.cpp
//
// Exact mate-rescue pruning (src/rescue_prune.h) must never change a consumed
// rescue result. Checked on generated jobs (random, mutated copies of the
// window, planted near-threshold matches with a weaker second copy, tandem
// repeats, poly-A, tiny lengths, windows past the NEON filter's capacity, and
// N bases in both sequences):
//
//   1. Against the independent scalar Smith-Waterman, ksw_align2, at bwa's
//      default scoring and several rescue thresholds (min_seed_len * a):
//        - B1 (proven failure): the full window scores below the threshold;
//        - B2 (hull [hb, he]): whenever the full window reaches the threshold,
//          the full window and the hull alone agree on score, qe, te, tb, qb
//          and score2 (te/tb shifted by hb).
//   2. The guard exits the lemma depends on: any N, a threshold below 5, and a
//      query or window beyond the scratch capacity all return FULL.
//   3. On aarch64, the NEON filter (rescue_prune_window at threshold 19)
//      against the int32 scalar filter (rescue_prune_window_scalar) and, where
//      it fits, the int16 reference lean(): identical (kind, hb, he).
//
// Mates stay at <= 250 bases so every job is in the 8-bit kernel's domain,
// the only one the filter runs in.

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <memory>
#include <random>
#include <vector>

#include "doctest/doctest.h"
#include "ksw.h"
#include "rescue_prune.h"
#include "scoring.h"

namespace {

constexpr int kNeonMinsc = 19;       // min_seed_len * a at the defaults; the NEON filter's threshold
constexpr int kDefaultMaxHits = 400; // BWA3_RESCUE_PRUNE_MAX_HITS default
constexpr int kNoGate = 1 << 30;
constexpr int kGapOpen = 6, kGapExtend = 1;  // the only gaps the pruning lemma is derived for
constexpr uint8_t kN = 4;

struct Job {
    std::vector<uint8_t> ref, q;
    const char *tag;
};

std::vector<uint8_t> random_bases(std::mt19937 &rng, int n)
{
    std::vector<uint8_t> v((size_t)n);
    for (auto &b : v) b = (uint8_t)(rng() & 3);
    return v;
}

// A copy of ref[at, at + len) with `subs` substitutions and `indels` 1-base
// insertions or deletions, so the mate shares 5-mer runs with the window.
std::vector<uint8_t> mutated_copy(std::mt19937 &rng, const std::vector<uint8_t> &ref, int at, int len,
                                  int subs, int indels)
{
    std::vector<uint8_t> q(ref.begin() + at, ref.begin() + at + len);
    for (int k = 0; k < subs && !q.empty(); k++) {
        auto &b = q[rng() % q.size()];
        b = (uint8_t)((b + 1 + rng() % 3) & 3);
    }
    for (int k = 0; k < indels && q.size() > 2; k++) {
        const size_t p = 1 + rng() % (q.size() - 2);
        if (rng() & 1) q.erase(q.begin() + (long)p);
        else q.insert(q.begin() + (long)p, (uint8_t)(rng() & 3));
    }
    return q;
}

// A random mate carrying one exact copy of a window segment (score ~seg) and a
// weaker, mismatched copy of part of it elsewhere in the window, so the best
// score sits just above the threshold and score2 near it: the hull's tail
// slack is what keeps that second copy, and with it score2.
Job planted(std::mt19937 &rng, int len1, int len2, int seg)
{
    auto ref = random_bases(rng, len1);
    auto q = random_bases(rng, len2);
    const int at = (int)(rng() % (unsigned)(len1 - seg)), qat = (int)(rng() % (unsigned)(len2 - seg));
    std::copy(ref.begin() + at, ref.begin() + at + seg, q.begin() + qat);
    const int seg2 = seg - 2 - (int)(rng() % 4);
    const int at2 = (int)(rng() % (unsigned)(len1 - seg2));
    std::copy(q.begin() + qat, q.begin() + qat + seg2, ref.begin() + at2);
    ref[(size_t)(at2 + seg2 / 2)] = (uint8_t)((ref[(size_t)(at2 + seg2 / 2)] + 1) & 3);
    return {std::move(ref), std::move(q), "planted near-threshold"};
}

std::vector<Job> build_jobs(std::mt19937 &rng)
{
    std::vector<Job> jobs;
    for (int k = 0; k < 80; k++) {  // mutated copies: B2 hulls of every width
        const int len1 = 150 + (int)(rng() % 1100), len2 = 40 + (int)(rng() % 210);
        auto ref = random_bases(rng, len1);
        const int len = std::min(len2, len1 - 1), at = (int)(rng() % (unsigned)(len1 - len));
        auto q = mutated_copy(rng, ref, at, len, (int)(rng() % 12), (int)(rng() % 4));
        jobs.push_back({std::move(ref), std::move(q), "mutated copy"});
    }
    for (int k = 0; k < 120; k++)  // near the thresholds tested below, score2 included
        jobs.push_back(planted(rng, 200 + (int)(rng() % 1000), 60 + (int)(rng() % 190), 12 + (int)(rng() % 26)));
    for (int k = 0; k < 60; k++) {  // unrelated: mostly B1
        const int len1 = 100 + (int)(rng() % 1200), len2 = 30 + (int)(rng() % 220);
        jobs.push_back({random_bases(rng, len1), random_bases(rng, len2), "random"});
    }
    for (int k = 0; k < 24; k++) {  // tandem repeats: repeated codes, multi-occurrence layers
        const int unit = 2 + (int)(rng() % 6), len1 = 200 + (int)(rng() % 600), len2 = 60 + (int)(rng() % 120);
        const auto motif = random_bases(rng, unit);
        std::vector<uint8_t> ref((size_t)len1), q((size_t)len2);
        for (int i = 0; i < len1; i++) ref[(size_t)i] = motif[(size_t)(i % unit)];
        for (int j = 0; j < len2; j++) q[(size_t)j] = motif[(size_t)((j + k) % unit)];
        for (int m = 0; m < 3; m++) ref[rng() % ref.size()] = (uint8_t)(rng() & 3);
        jobs.push_back({std::move(ref), std::move(q), "tandem repeat"});
    }
    for (int k = 0; k < 6; k++) {  // poly-A: the hit gate, and the int16 fallback without it
        jobs.push_back({std::vector<uint8_t>((size_t)(300 + k * 150), 0),
                        std::vector<uint8_t>((size_t)(60 + k * 30), 0), "poly-A"});
    }
    for (int k = 0; k < 30; k++) {  // tiny lengths around the 5-mer and threshold edges
        const int len1 = 5 + (int)(rng() % 40), len2 = 5 + (int)(rng() % 30);
        auto ref = random_bases(rng, len1);
        auto q = (rng() & 1) ? mutated_copy(rng, ref, 0, std::min(len1, len2), 0, 0) : random_bases(rng, len2);
        jobs.push_back({std::move(ref), std::move(q), "tiny"});
    }
    for (int k = 0; k < 3; k++)  // past the NEON filter's capacity (~4000 rows): the scalar filter decides
        jobs.push_back(planted(rng, 4200 + k * 400, 150, 30));
    return jobs;
}

// Jobs with an N in the window or the mate, placed in the vector-scanned
// prefix and in the scalar tail of each: the filter must refuse all of them.
std::vector<Job> build_n_jobs(std::mt19937 &rng)
{
    std::vector<Job> jobs;
    for (int k = 0; k < 16; k++) {
        Job jb = planted(rng, 300 + (int)(rng() % 500), 70 + (int)(rng() % 150), 30);
        std::vector<uint8_t> &seq = (k & 1) ? jb.q : jb.ref;
        const size_t n = seq.size();
        seq[(k & 2) ? (size_t)(rng() % 16) : n - 1 - (size_t)(rng() % (n % 16 + 1))] = kN;
        jb.tag = "N base";
        jobs.push_back(std::move(jb));
    }
    return jobs;
}

// ksw_align2 at the rescue settings for threshold `minsc`. It reverses its
// target in place for the start pass, so both sequences are copied per call.
kswr_t scalar_sw(const std::vector<uint8_t> &q, const uint8_t *ref, int len1, int minsc,
                 const bwa_tests::ScoringMatrix &mat)
{
    std::vector<uint8_t> qq(q), rr(ref, ref + len1);
    const int xtra = KSW_XSUBO | KSW_XSTART | KSW_XBYTE | minsc;
    return ksw_align2((int)qq.size(), qq.data(), len1, rr.data(), 5, mat.data(),
                      kGapOpen, kGapExtend, kGapOpen, kGapExtend, xtra, nullptr);
}

// Exact 5-mer hits between window and mate (the int16 reference lean() is exact only below 32000).
long count_hits(const Job &jb)
{
    std::vector<int> qcnt(1024, 0);
    const int len1 = (int)jb.ref.size(), len2 = (int)jb.q.size();
    for (int j = 4; j < len2; j++) {
        int c = 0;
        for (int t = j - 4; t <= j; t++) c = (c << 2) | jb.q[(size_t)t];
        qcnt[(size_t)c]++;
    }
    long hits = 0;
    for (int i = 4; i < len1; i++) {
        int c = 0;
        for (int t = i - 4; t <= i; t++) c = (c << 2) | jb.ref[(size_t)t];
        hits += qcnt[(size_t)c];
    }
    return hits;
}

struct Tally { int b1 = 0, b2 = 0, full = 0; };

// One job at one threshold and hit gate, against the full-window oracle.
void check_against_oracle(const Job &jb, int minsc, int max_hits, const bwa_tests::ScoringMatrix &mat,
                          Tally &t)
{
    const int len1 = (int)jb.ref.size(), len2 = (int)jb.q.size();
    int hb = -2, he = -2;
    const int kind = rescue_prune_window(jb.ref.data(), len1, jb.q.data(), len2, minsc, max_hits, &hb, &he);
    CAPTURE(jb.tag); CAPTURE(len1); CAPTURE(len2); CAPTURE(minsc); CAPTURE(max_hits); CAPTURE(kind);
    if (kind == RESCUE_PRUNE_FULL) { ++t.full; return; }
    const kswr_t full = scalar_sw(jb.q, jb.ref.data(), len1, minsc, mat);
    CAPTURE(full.score); CAPTURE(full.te); CAPTURE(full.score2);
    if (kind == RESCUE_PRUNE_B1) {
        ++t.b1;
        CHECK(full.score < minsc);
        return;
    }
    REQUIRE(kind == RESCUE_PRUNE_B2);
    ++t.b2;
    CAPTURE(hb); CAPTURE(he);
    REQUIRE(hb >= 0);
    REQUIRE(hb <= he);
    REQUIRE(he < len1);
    // Below the threshold the rescue gate fails and nothing else is consumed; a
    // sub-window can only score lower, so there is nothing to compare there.
    if (full.score < minsc) return;
    const kswr_t sub = scalar_sw(jb.q, jb.ref.data() + hb, he - hb + 1, minsc, mat);
    CAPTURE(sub.score); CAPTURE(sub.te); CAPTURE(sub.score2);
    CHECK(sub.score == full.score);
    CHECK(sub.qe == full.qe);
    CHECK(sub.te + hb == full.te);
    CHECK(sub.qb == full.qb);
    CHECK(sub.tb + hb == full.tb);
    CHECK(sub.score2 == full.score2);
}

} // namespace

TEST_CASE("rescue prune: B1 and B2 decisions reproduce every consumed ksw_align2 field"
          * doctest::test_suite("unit/rescue")) {
    std::mt19937 rng(20260928);
    const auto jobs = build_jobs(rng);
    const auto mat = bwa_tests::default_scoring_matrix();
    Tally t19;
    for (const int max_hits : {kDefaultMaxHits, kNoGate})
        for (const Job &jb : jobs) check_against_oracle(jb, kNeonMinsc, max_hits, mat, t19);
    MESSAGE("threshold 19: B1=" << t19.b1 << " B2=" << t19.b2 << " FULL=" << t19.full);
    CHECK(t19.b1 > 0);
    CHECK(t19.b2 > 0);
    CHECK(t19.full > 0);

    // Every other threshold (-k with a = 1) runs the scalar filter, whose hull
    // slack is ub - minsc - 2; pin it away from 19 too.
    for (const int minsc : {5, 10, 32}) {
        Tally t;
        for (const Job &jb : jobs) check_against_oracle(jb, minsc, kNoGate, mat, t);
        CAPTURE(minsc);
        MESSAGE("threshold " << minsc << ": B1=" << t.b1 << " B2=" << t.b2 << " FULL=" << t.full);
        CHECK(t.b2 > 0);
    }
}

TEST_CASE("rescue prune: N, a threshold below 5 and oversized inputs keep the full window"
          * doctest::test_suite("unit/rescue")) {
    std::mt19937 rng(31337);
    int hb = -2, he = -2;
    for (const Job &jb : build_n_jobs(rng)) {
        CAPTURE(jb.ref.size()); CAPTURE(jb.q.size());
        CHECK(rescue_prune_window(jb.ref.data(), (int)jb.ref.size(), jb.q.data(), (int)jb.q.size(),
                                  kNeonMinsc, kNoGate, &hb, &he) == RESCUE_PRUNE_FULL);
        CHECK(rescue_prune_window(jb.ref.data(), (int)jb.ref.size(), jb.q.data(), (int)jb.q.size(),
                                  10, kNoGate, &hb, &he) == RESCUE_PRUNE_FULL);
    }
    const Job small = planted(rng, 400, 120, 30);
    CHECK(rescue_prune_window(small.ref.data(), 400, small.q.data(), 120, 4, kNoGate, &hb, &he)
          == RESCUE_PRUNE_FULL);
    const auto long_q = random_bases(rng, rescue_prune_scratch::QCAP + 1);
    CHECK(rescue_prune_window(small.ref.data(), 400, long_q.data(), (int)long_q.size(), kNeonMinsc, kNoGate,
                              &hb, &he) == RESCUE_PRUNE_FULL);
    const auto long_ref = random_bases(rng, 30001);
    CHECK(rescue_prune_window(long_ref.data(), (int)long_ref.size(), small.q.data(), 120, kNeonMinsc, kNoGate,
                              &hb, &he) == RESCUE_PRUNE_FULL);
}

TEST_CASE("rescue prune: the NEON filter decides exactly as the scalar filters"
          * doctest::test_suite("unit/rescue")) {
#if !defined(__aarch64__)
    MESSAGE("skipped: the NEON rescue-prune filter is aarch64-only");
    return;
#else
    std::mt19937 rng(8675309);
    auto jobs = build_jobs(rng);
    const auto n_jobs = build_n_jobs(rng);
    jobs.insert(jobs.end(), n_jobs.begin(), n_jobs.end());
    std::unique_ptr<rescue_prune_scratch> scalar(new rescue_prune_scratch());
    std::unique_ptr<rescue_prune_neon::Scratch> ref16(new rescue_prune_neon::Scratch());
    int compared_lean = 0;
    for (const int max_hits : {kDefaultMaxHits, kNoGate}) {
        for (size_t i = 0; i < jobs.size(); i++) {
            const Job &jb = jobs[i];
            const int len1 = (int)jb.ref.size(), len2 = (int)jb.q.size();
            CAPTURE(jb.tag); CAPTURE(i); CAPTURE(len1); CAPTURE(len2); CAPTURE(max_hits);
            int hb = -2, he = -2, shb = -2, she = -2;
            const int kind = rescue_prune_window(jb.ref.data(), len1, jb.q.data(), len2, kNeonMinsc, max_hits,
                                                 &hb, &he);
            const int skind = rescue_prune_window_scalar(jb.ref.data(), len1, jb.q.data(), len2, kNeonMinsc,
                                                         max_hits, *scalar, &shb, &she);
            CHECK(kind == skind);
            if (kind == RESCUE_PRUNE_B2 && skind == RESCUE_PRUNE_B2) {
                CHECK(hb == shb);
                CHECK(he == she);
            }
            // lean() is the int16 reference the NEON rewrite was derived from; it has no hit gate,
            // a fixed capacity and int16 sums, so compare it only without the gate, where it
            // fits, and below 32000 hits.
            if (max_hits == kNoGate && len1 + kswv_query_quantum8(len2) + 1 <= 4096 && len1 >= 5 && len2 >= 5
                && count_hits(jb) <= 32000) {
                const rescue_prune_neon::Job lj{len1, len2, 0, 0, -1, -1, jb.ref.data(), jb.q.data()};
                int lhb = -2, lhe = -2;
                const rescue_prune_neon::Kind lk = rescue_prune_neon::lean(lj, *ref16, lhb, lhe);
                ++compared_lean;
                CHECK((int)lk == skind);
                if (lk == rescue_prune_neon::B2 && skind == RESCUE_PRUNE_B2) {
                    CHECK(lhb == shb);
                    CHECK(lhe == she);
                }
            }
        }
    }
    CHECK(compared_lean > 0);
#endif
}
