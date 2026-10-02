// test/unit/test_rescue_prune.cpp
//
// Exact mate-rescue pruning (src/rescue_prune.h) must never change a consumed
// rescue result. Checked on generated jobs (random, mutated copies of the
// window, planted near-threshold matches with a weaker second copy, tandem
// repeats, poly-A, tiny lengths, windows past the SIMD filter's capacity, and
// N bases in both sequences):
//
//   1. Against the independent scalar Smith-Waterman, ksw_align2, at bwa's
//      default scoring and several rescue thresholds (min_seed_len * a):
//        - B1 (proven failure): the full window scores below the threshold;
//        - B2 (hull [hb, he]): whenever the full window reaches the threshold,
//          the full window and the hull alone agree on score, qe, te, tb, qb
//          and score2 (te/tb shifted by hb).
//      The same at other scorings the lemma admits (rescue_prune_params):
//      -B 6, -O 8 -E 2, split gap costs and -A 2 -B 8 -O 12 -E 2, and the
//      refused ones (-B 3, a cheap insertion) keep every window in full.
//      And under --meth (rescue_prune_params::set_meth), against ksw_align2
//      with the genomic, neutral and collapsed OT / OB matrices, on mates
//      with converted bases.
//   2. The guard exits the lemma depends on: any N, a threshold below 5, and a
//      query or window beyond the scratch capacity all return FULL.
//   3. Where a SIMD filter is compiled in (NEON on aarch64, the SSE4.1 / SSSE3
//      port on x86 AVX2 builds), that filter (rescue_prune_window, at 19 and
//      at thresholds 5 / 10 / 32 in rotation, and at scorings with other bound
//      weights: c = 2, a = 2, -x intractg's tail) against the int32 scalar
//      filter (rescue_prune_window_scalar) and, where it fits, the int16
//      reference lean(): identical (kind, hb, he).
//   4. There too, the SIMD filter's repeat memo and its component list: a
//      repeated job returns the first call's decision and view, the same job
//      at another threshold is not answered from the memo, and band
//      planning's components taken from the filter's list equal a rescan of
//      the view, at every cap.
//
// Most mates stay at <= 250 bases, in the 8-bit kernel's domain; a long-mate
// category (251 to 600 bases) and the a = 2 scoring put jobs on the 16-bit
// kernels, which the filter runs on too, and the oracle then runs ksw_align2's
// 16-bit path (scalar_sw).

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "doctest/doctest.h"
#include "bwamem.h"
#include "ksw.h"
#include "matesw_u8.h"
#include "meth_scoring.h"
#include "rescue_band.h"
#include "rescue_prune.h"
#include "rescue_prune_test.h"
#include "scoring.h"
#include "simd_dispatch.h"

namespace {

constexpr int kSimdMinsc = 19;       // min_seed_len * a at the defaults
constexpr int kDefaultMaxHits = 400; // BWA3_RESCUE_PRUNE_MAX_HITS default
constexpr int kNoGate = 1 << 30;
constexpr int kGapOpen = 6, kGapExtend = 1;  // the default gap costs
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
    for (int k = 0; k < 3; k++)  // past the SIMD filter's capacity (~4000 rows): the scalar filter decides
        jobs.push_back(planted(rng, 4200 + k * 400, 150, 30));
    for (int k = 0; k < 12; k++) {  // long mates, past the 8-bit kernels: the 16-bit path at every a
        const int len2 = 251 + (int)(rng() % 350), len1 = len2 + 50 + (int)(rng() % 500);
        if (k & 1) {
            jobs.push_back(planted(rng, len1, len2, 20 + (int)(rng() % 30)));
            jobs.back().tag = "long mate, planted";
        } else {
            auto ref = random_bases(rng, len1);
            const int at = (int)(rng() % (unsigned)(len1 - len2));
            auto q = mutated_copy(rng, ref, at, len2, (int)(rng() % 20), (int)(rng() % 4));
            jobs.push_back({std::move(ref), std::move(q), "long mate, mutated copy"});
        }
    }
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

// Whether production runs this mate on the 8-bit kernels (matesw_use_u8, with the kernels'
// bias, the mismatch penalty b, which is the matrices' minimum here). Longer mates take the
// 16-bit kernels, which pruning covers too.
bool u8_job(int len2, const rescue_prune_params &p)
{
    return matesw_use_u8(len2, p.a, std::max(1, p.b)) != 0;
}

// ksw_align2 at the rescue settings for threshold p.minsc and p's gap costs
// (the matrix carries a and b), on its 8-bit path for an 8-bit job and its 16-bit
// path otherwise, as production routes the job (u8_job). It reverses its target in
// place for the start pass, so both sequences are copied per call.
kswr_t scalar_sw(const std::vector<uint8_t> &q, const uint8_t *ref, int len1, const rescue_prune_params &p,
                 const bwa_tests::ScoringMatrix &mat)
{
    std::vector<uint8_t> qq(q), rr(ref, ref + len1);
    const int xtra = KSW_XSUBO | KSW_XSTART | (u8_job((int)q.size(), p) ? KSW_XBYTE : 0) | p.minsc;
    return ksw_align2((int)qq.size(), qq.data(), len1, rr.data(), 5, mat.data(),
                      p.o_del, p.e_del, p.o_ins, p.e_ins, xtra, nullptr);
}

bool has_n(const std::vector<uint8_t> &seq)
{
    return std::any_of(seq.begin(), seq.end(), [](uint8_t b) { return b >= kN; });
}

// Exact 5-mer hits between window and mate (the int16 reference lean() is exact only below 32000).
// A 5-mer with an N matches nothing and is skipped, which also keeps every code inside the
// 1024-entry table.
long count_hits(const Job &jb)
{
    std::vector<int> qcnt(1024, 0);
    const int len1 = (int)jb.ref.size(), len2 = (int)jb.q.size();
    auto code = [](const std::vector<uint8_t> &seq, int end) {
        int c = 0;
        for (int t = end - 4; t <= end; t++) {
            if (seq[(size_t)t] >= kN) return -1;
            c = (c << 2) | seq[(size_t)t];
        }
        return c;
    };
    for (int j = 4; j < len2; j++) {
        const int c = code(jb.q, j);
        if (c >= 0) qcnt[(size_t)c]++;
    }
    long hits = 0;
    for (int i = 4; i < len1; i++) {
        const int c = code(jb.ref, i);
        if (c >= 0) hits += qcnt[(size_t)c];
    }
    return hits;
}

struct Tally { int b1 = 0, b2 = 0, full = 0; };

// ksw_align2 over one job's full window, and over its last B2 hull, at one
// scoring and threshold: each computed at most once however many hit gates the
// job is checked under (the gate changes the decision, never the oracle).
struct Oracle {
    const Job &jb;
    const rescue_prune_params p;
    const int minsc;
    const bwa_tests::ScoringMatrix &mat;
    bool have_full = false;
    kswr_t full{};
    int sub_hb = -1, sub_he = -1;
    kswr_t sub{};

    Oracle(const Job &j, int m, const bwa_tests::ScoringMatrix &sm)
        : jb(j), p(rescue_prune_params::defaults(m)), minsc(m), mat(sm) {}
    Oracle(const Job &j, const rescue_prune_params &sp, const bwa_tests::ScoringMatrix &sm)
        : jb(j), p(sp), minsc(sp.minsc), mat(sm) {}
    const kswr_t &full_window()
    {
        if (!have_full) {
            full = scalar_sw(jb.q, jb.ref.data(), (int)jb.ref.size(), p, mat);
            have_full = true;
        }
        return full;
    }
    const kswr_t &hull(int hb, int he)
    {
        if (hb != sub_hb || he != sub_he) {
            sub = scalar_sw(jb.q, jb.ref.data() + hb, he - hb + 1, p, mat);
            sub_hb = hb;
            sub_he = he;
        }
        return sub;
    }
};

// One job at one hit gate, against the full-window oracle at o.minsc.
void check_against_oracle(int max_hits, Oracle &o, Tally &t)
{
    const Job &jb = o.jb;
    const int minsc = o.minsc;
    const int len1 = (int)jb.ref.size(), len2 = (int)jb.q.size();
    int hb = -2, he = -2;
    const int kind = rescue_prune_window(jb.ref.data(), len1, jb.q.data(), len2, o.p, max_hits, &hb, &he);
    CAPTURE(o.p.a); CAPTURE(o.p.b); CAPTURE(o.p.o_del); CAPTURE(o.p.e_del); CAPTURE(o.p.o_ins); CAPTURE(o.p.e_ins);
    CAPTURE(jb.tag); CAPTURE(len1); CAPTURE(len2); CAPTURE(minsc); CAPTURE(max_hits); CAPTURE(kind);
    if (kind == RESCUE_PRUNE_FULL) { ++t.full; return; }
    const kswr_t &full = o.full_window();
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
    const kswr_t &sub = o.hull(hb, he);
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
          * doctest::test_suite("unit/pair")) {
    std::mt19937 rng(20260928);
    const auto jobs = build_jobs(rng);
    const auto mat = bwa_tests::default_scoring_matrix();
    Tally t19;
    for (const Job &jb : jobs) {
        Oracle o(jb, kSimdMinsc, mat);
        for (const int max_hits : {kDefaultMaxHits, kNoGate}) check_against_oracle(max_hits, o, t19);
    }
    MESSAGE("threshold 19: B1=" << t19.b1 << " B2=" << t19.b2 << " FULL=" << t19.full);
    CHECK(t19.b1 > 0);
    CHECK(t19.b2 > 0);
    CHECK(t19.full > 0);

    // Every other threshold (-k with a = 1) runs the scalar filter, whose hull
    // slack is ub - minsc - 2; pin it away from 19 too. Each job is checked at
    // one of them in rotation: every job category (each at least three jobs,
    // built contiguously) still reaches every threshold, at a third of the
    // oracle calls, which keeps this case inside the unit-test time budget.
    const int other_minsc[3] = {5, 10, 32};
    Tally t[3];
    for (size_t i = 0; i < jobs.size(); i++) {
        Oracle o(jobs[i], other_minsc[i % 3], mat);
        check_against_oracle(kNoGate, o, t[i % 3]);
    }
    for (int k = 0; k < 3; k++) {
        CAPTURE(other_minsc[k]);
        MESSAGE("threshold " << other_minsc[k] << ": B1=" << t[k].b1 << " B2=" << t[k].b2
                             << " FULL=" << t[k].full);
        CHECK(t[k].b2 > 0);
    }
}

// Other scorings: the lemma with K-mers, the charge c and the tail scaled to the deletion costs
// (rescue_prune_params). Each job is checked at one scoring in rotation, at the default -k 19
// (minsc = 19 a), which keeps the case inside the unit-test time budget. Mates that leave the
// 8-bit kernels (matesw_use_u8 false: at a = 2, mates over about 120 bases) are pruned on the
// 16-bit kernels, and the oracle takes ksw_align2's 16-bit path for them (scalar_sw); the case
// requires some of them. The scorings the lemma refuses must keep every window in full.
TEST_CASE("rescue prune: decisions at other scorings reproduce every consumed ksw_align2 field"
          * doctest::test_suite("unit/pair")) {
    std::mt19937 rng(4242);
    const auto jobs = build_jobs(rng);
    struct Sc { int a, b, o_del, e_del, o_ins, e_ins; const char *name; };
    const Sc admitted[4] = {{1, 6, 6, 1, 6, 1, "-B 6"}, {1, 4, 8, 2, 8, 2, "-O 8 -E 2"},
                            {1, 4, 6, 1, 7, 2, "-O 6,7 -E 1,2"}, {2, 8, 12, 2, 12, 2, "-A 2 -B 8 -O 12 -E 2"}};
    const Sc refused[2] = {{1, 3, 6, 1, 6, 1, "-B 3"}, {1, 4, 6, 1, 2, 1, "-O 6,2 -E 1"}};
    Tally t[4], t16;
    for (size_t i = 0; i < jobs.size(); i++) {
        const Sc &s = admitted[i % 4];
        const std::string name(s.name);
        CAPTURE(name);
        const rescue_prune_params p = rescue_prune_params::from(s.a, s.b, s.o_del, s.e_del, s.o_ins, s.e_ins,
                                                                kSimdMinsc * s.a);
        REQUIRE(p.valid);
        const auto mat = bwa_tests::build_scoring_matrix(s.a, s.b, 1);
        Oracle o(jobs[i], p, mat);
        check_against_oracle(kNoGate, o, u8_job((int)jobs[i].q.size(), p) ? t[i % 4] : t16);
    }
    MESSAGE("16-bit jobs: B1=" << t16.b1 << " B2=" << t16.b2 << " FULL=" << t16.full);
    CHECK(t16.b1 > 0);
    CHECK(t16.b2 > 0);
    for (int k = 0; k < 4; k++) {
        const std::string name(admitted[k].name);
        CAPTURE(name);
        MESSAGE(std::string(admitted[k].name) << ": B1=" << t[k].b1 << " B2=" << t[k].b2 << " FULL=" << t[k].full);
        CHECK(t[k].b1 > 0);
        CHECK(t[k].b2 > 0);
    }
    for (const Sc &s : refused) {
        const std::string name(s.name);
        CAPTURE(name);
        const rescue_prune_params p = rescue_prune_params::from(s.a, s.b, s.o_del, s.e_del, s.o_ins, s.e_ins,
                                                                kSimdMinsc * s.a);
        CHECK(!p.valid);
        int hb = -2, he = -2;
        for (size_t i = 0; i < jobs.size(); i += 7)
            CHECK(rescue_prune_window(jobs[i].ref.data(), (int)jobs[i].ref.size(), jobs[i].q.data(),
                                      (int)jobs[i].q.size(), p, kNoGate, &hb, &he) == RESCUE_PRUNE_FULL);
    }
}

// K > 5: a scoring whose mismatch and gap costs admit longer K-mers filters with them, up to the
// caller's k_max. Three scorings in rotation: -B 6 (K = 7), -x intractg (its costs admit K = 10, the
// cap 8 takes K = 8) and -B 8 under a cap of 6 (K = 6). Each job runs through the ksw_align2 oracle,
// and through the wrapper (the SIMD K-mer filter where the build has one, which must decide some
// jobs at every K) against a scalar filter whose tables are cleared in full every 8 jobs, which
// checks the touched-code reset of the tables the wrapper keeps across queries.
TEST_CASE("rescue prune: K-mer decisions up to K = 8 reproduce ksw_align2 and the scalar filter"
          * doctest::test_suite("unit/pair")) {
    std::mt19937 rng(9001);
    const auto jobs = build_jobs(rng);
    struct Sc { int a, b, o_del, e_del, o_ins, e_ins, kmax, K; const char *name; };
    const Sc sc[3] = {{1, 6, 6, 1, 6, 1, 8, 7, "-B 6"}, {1, 9, 16, 1, 16, 1, 8, 8, "-x intractg"},
                      {1, 8, 6, 1, 6, 1, 6, 6, "-B 8 at k_max 6"}};
    Tally t[3];
    std::unique_ptr<rescue_prune_scratch> fresh(new rescue_prune_scratch());
    int n_simd[3] = {0, 0, 0};
    for (size_t i = 0; i < jobs.size(); i++) {
        const Sc &s = sc[i % 3];
        const std::string name(s.name);
        CAPTURE(name);
        const rescue_prune_params p = rescue_prune_params::from(s.a, s.b, s.o_del, s.e_del, s.o_ins, s.e_ins,
                                                                kSimdMinsc * s.a, 5, s.kmax);
        REQUIRE(p.valid);
        REQUIRE(p.K == s.K);
        const auto mat = bwa_tests::build_scoring_matrix(s.a, s.b, 1);
        Oracle o(jobs[i], p, mat);
        check_against_oracle(kNoGate, o, t[i % 3]);
        const Job &jb = jobs[i];
        const int len1 = (int)jb.ref.size(), len2 = (int)jb.q.size();
        if (i % 8 == 0) fresh->clear_tables();
        int hb = -2, he = -2, shb = -2, she = -2;
        rescue_prune_view v;
        const int kind = rescue_prune_window(jb.ref.data(), len1, jb.q.data(), len2, p, kDefaultMaxHits, &hb, &he, &v);
        const int skind = rescue_prune_window_scalar(jb.ref.data(), len1, jb.q.data(), len2, p, kDefaultMaxHits,
                                                     *fresh, &shb, &she);
        CAPTURE(jb.tag); CAPTURE(len1); CAPTURE(len2);
        CHECK(kind == skind);
        if (kind == RESCUE_PRUNE_B2 && skind == RESCUE_PRUNE_B2) {
            CHECK(hb == shb);
            CHECK(he == she);
        }
        n_simd[i % 3] += kind == RESCUE_PRUNE_B2 && v.bnd16 != nullptr;
    }
    for (int k = 0; k < 3; k++) {
        const std::string name(sc[k].name);
        CAPTURE(name);
        MESSAGE(name << ": B1=" << t[k].b1 << " B2=" << t[k].b2 << " FULL=" << t[k].full);
        CHECK(t[k].b1 > 0);
        CHECK(t[k].b2 > 0);
#if RESCUE_PRUNE_HAVE_SIMD
        CHECK(n_simd[k] > 0);   // the SIMD filter's K-mer instantiation decided some, not only the scalar fallback
#endif
    }
}

// rescue_prune_params::from takes the largest valid K up to its cap: where the costs admit a K that
// fails the charge (c > 0) or the threshold (minsc > (K - 1) a), a smaller K that passes is taken
// rather than refusing the scoring, so raising the cap never turns pruning off.
TEST_CASE("rescue prune: the K cap takes the largest valid K, never refusing a scoring a smaller K admits"
          * doctest::test_suite("unit/pair")) {
    // -B 6 -O 6,5 -E 1: the costs admit K = 7, whose insertion charge 5 + 1 - 6 is 0; K = 6 passes.
    rescue_prune_params p = rescue_prune_params::from(1, 6, 6, 1, 5, 1, kSimdMinsc, 5, 8);
    CHECK(p.valid);
    CHECK(p.K == 6);
    CHECK(p.c == 1);
    // -B 6 at minsc 6 (-k 6): K = 7 needs minsc > 6; K = 6 passes, as K = 5 did before the cap rose.
    p = rescue_prune_params::from(1, 6, 6, 1, 6, 1, 6, 5, 8);
    CHECK(p.valid);
    CHECK(p.K == 6);
    // The same scoring at minsc 19 takes K = 7, and a cap of 5 keeps K = 5.
    CHECK(rescue_prune_params::from(1, 6, 6, 1, 6, 1, kSimdMinsc, 5, 8).K == 7);
    CHECK(rescue_prune_params::from(1, 6, 6, 1, 6, 1, kSimdMinsc, 5, 5).K == 5);
    // Every K refused (minsc 4 at -B 6 fails from K = 7 down to k_min 5): not valid.
    CHECK(!rescue_prune_params::from(1, 6, 6, 1, 6, 1, 4, 5, 8).valid);
}

// --meth: the filter matches C -> T (OT) or G -> A (OB) converted copies of the window and the mate,
// which over-counts hits under every meth matrix, so its decisions must reproduce ksw_align2 with the
// matrix itself (mat[ref * 5 + read], built by mem_opt_fill_meth_mat through
// bwa_tests::meth_scoring_matrix: the conversion cell freed to +a, or to 0 under neutral, and under
// collapsed its mirror too), at -B 4 so the lemma admits collapsed. Each job is checked at one (matrix, hypothesis) pair in rotation, with half of its
// mate's convertible bases converted.
TEST_CASE("rescue prune: --meth decisions reproduce ksw_align2 with the meth matrix"
          * doctest::test_suite("unit/pair")) {
    std::mt19937 rng(5150);
    const auto jobs = build_jobs(rng);
    const char *names[3] = {"genomic", "neutral", "collapsed"};
    std::vector<bwa_tests::ScoringMatrix> mats;
    const int scorings[3] = {MEM_METH_SCORING_GENOMIC, MEM_METH_SCORING_NEUTRAL, MEM_METH_SCORING_COLLAPSED};
    for (int kind = 0; kind < 3; kind++)
        for (int hyp = 0; hyp < 2; hyp++)   // hyp 1 = OT: ref C / read T; 0 = OB: ref G / read A
            mats.push_back(bwa_tests::meth_scoring_matrix(scorings[kind], hyp == 1, 1, 4));
    Tally t[3];
    for (size_t i = 0; i < jobs.size(); i++) {
        const int kind = (int)(i % 3), hyp = (int)((i / 3) % 2);
        const std::string name(names[kind]);
        CAPTURE(name); CAPTURE(hyp);
        Job jb = jobs[i];
        rescue_convert_mate(jb.q, hyp, 2, rng);
        rescue_prune_params p = rescue_prune_params::defaults(kSimdMinsc);
        p.set_meth(hyp);
        Oracle o(jb, p, mats[(size_t)(kind * 2 + hyp)]);
        check_against_oracle(kNoGate, o, t[kind]);
    }
    for (int k = 0; k < 3; k++) {
        const std::string name(names[k]);
        CAPTURE(name);
        MESSAGE("--meth " << name << ": B1=" << t[k].b1 << " B2=" << t[k].b2 << " FULL=" << t[k].full);
        CHECK(t[k].b1 > 0);
        CHECK(t[k].b2 > 0);
    }
}

// The --meth relation (set_meth_rel): a read T relates to a reference T or C (OT), a read A to a
// reference A or G (OB), every other pair only to itself -- the exact match cells of the genomic
// matrix and a superset of the neutral one's. Mates mostly unconverted (a TAPS-like 5 % conversion),
// where converted copies would relate nearly everything.
TEST_CASE("rescue prune: --meth relation decisions reproduce ksw_align2 with the genomic and neutral matrices"
          * doctest::test_suite("unit/pair")) {
    std::mt19937 rng(7117);
    const auto jobs = build_jobs(rng);
    const char *names[2] = {"genomic", "neutral"};
    std::vector<bwa_tests::ScoringMatrix> mats;
    const int scorings[2] = {MEM_METH_SCORING_GENOMIC, MEM_METH_SCORING_NEUTRAL};
    for (int kind = 0; kind < 2; kind++)
        for (int hyp = 0; hyp < 2; hyp++)   // hyp 1 = OT: ref C / read T; 0 = OB: ref G / read A
            mats.push_back(bwa_tests::meth_scoring_matrix(scorings[kind], hyp == 1, 1, 4));
    Tally t[2];
    for (size_t i = 0; i < jobs.size(); i++) {
        const int kind = (int)(i % 2), hyp = (int)((i / 2) % 2);
        const std::string name(names[kind]);
        CAPTURE(name); CAPTURE(hyp);
        Job jb = jobs[i];
        rescue_convert_mate(jb.q, hyp, 20, rng);
        rescue_prune_params p = rescue_prune_params::defaults(kSimdMinsc);
        p.set_meth_rel(hyp);
        Oracle o(jb, p, mats[(size_t)(kind * 2 + hyp)]);
        check_against_oracle(kNoGate, o, t[kind]);
    }
    for (int k = 0; k < 2; k++) {
        const std::string name(names[k]);
        CAPTURE(name);
        MESSAGE("--meth relation " << name << ": B1=" << t[k].b1 << " B2=" << t[k].b2 << " FULL=" << t[k].full);
        CHECK(t[k].b1 > 0);
        CHECK(t[k].b2 > 0);
    }
}

// The relation is exact only for a matrix that frees the single conversion cell: rescue_meth_rel
// checks the run's own matrices (rescue_meth_rel_matrix_ok) rather than --meth-scoring. The
// production matrices (mem_opt_fill_meth_mat, through meth_scoring_matrix): genomic and neutral pass
// under both hypotheses, at -B 4 and -B 6; collapsed, which frees the mirror cell too, fails even at
// -B 4, where the lemma admits it; a matrix checked against the other hypothesis fails; and so does
// any other cell raised above -b.
TEST_CASE("rescue prune: the relation admits the genomic and neutral matrices, not collapsed"
          * doctest::test_suite("unit/pair")) {
    for (const int b : {4, 6})
        for (const int hyp : {0, 1}) {
            CAPTURE(b); CAPTURE(hyp);
            const auto gen = bwa_tests::meth_scoring_matrix(MEM_METH_SCORING_GENOMIC, hyp == 1, 1, b);
            const auto neu = bwa_tests::meth_scoring_matrix(MEM_METH_SCORING_NEUTRAL, hyp == 1, 1, b);
            const auto col = bwa_tests::meth_scoring_matrix(MEM_METH_SCORING_COLLAPSED, hyp == 1, 1, b);
            CHECK(rescue_meth_rel_matrix_ok(gen.data(), hyp, 1, b));
            CHECK(rescue_meth_rel_matrix_ok(neu.data(), hyp, 1, b));
            CHECK(!rescue_meth_rel_matrix_ok(col.data(), hyp, 1, b));
            CHECK(!rescue_meth_rel_matrix_ok(gen.data(), !hyp, 1, b));
            for (int r = 0; r < 4; r++)
                for (int q = 0; q < 4; q++) {
                    if (r == q || gen[(size_t)(r * 5 + q)] != -b) continue;
                    auto m = gen;
                    m[(size_t)(r * 5 + q)] = (int8_t)(-b + 1);
                    CAPTURE(r); CAPTURE(q);
                    CHECK(!rescue_meth_rel_matrix_ok(m.data(), hyp, 1, b));
                }
        }
}

TEST_CASE("rescue prune: N, a threshold below 5 and oversized inputs keep the full window"
          * doctest::test_suite("unit/pair")) {
    std::mt19937 rng(31337);
    int hb = -2, he = -2;
    for (const Job &jb : build_n_jobs(rng)) {
        CAPTURE(jb.ref.size()); CAPTURE(jb.q.size());
        CHECK(rescue_prune_window(jb.ref.data(), (int)jb.ref.size(), jb.q.data(), (int)jb.q.size(),
                                  kSimdMinsc, kNoGate, &hb, &he) == RESCUE_PRUNE_FULL);
        CHECK(rescue_prune_window(jb.ref.data(), (int)jb.ref.size(), jb.q.data(), (int)jb.q.size(),
                                  10, kNoGate, &hb, &he) == RESCUE_PRUNE_FULL);
    }
    const Job small = planted(rng, 400, 120, 30);
    CHECK(rescue_prune_window(small.ref.data(), 400, small.q.data(), 120, 4, kNoGate, &hb, &he)
          == RESCUE_PRUNE_FULL);
    const auto long_q = random_bases(rng, rescue_prune_scratch::QCAP + 1);
    CHECK(rescue_prune_window(small.ref.data(), 400, long_q.data(), (int)long_q.size(), kSimdMinsc, kNoGate,
                              &hb, &he) == RESCUE_PRUNE_FULL);
    const auto long_ref = random_bases(rng, 30001);
    CHECK(rescue_prune_window(long_ref.data(), (int)long_ref.size(), small.q.data(), 120, kSimdMinsc, kNoGate,
                              &hb, &he) == RESCUE_PRUNE_FULL);
}

TEST_CASE("rescue prune: the SIMD filter decides exactly as the scalar filters"
          * doctest::test_suite("unit/pair")) {
#if !RESCUE_PRUNE_HAVE_SIMD
    MESSAGE("skipped: no SIMD rescue-prune filter in this build (aarch64, or x86 with AVX2)");
    return;
#else
    std::mt19937 rng(8675309);
    auto jobs = build_jobs(rng);
    const auto n_jobs = build_n_jobs(rng);
    jobs.insert(jobs.end(), n_jobs.begin(), n_jobs.end());
    std::unique_ptr<rescue_prune_scratch> scalar(new rescue_prune_scratch());
    std::unique_ptr<rescue_prune_neon::Scratch> ref16(new rescue_prune_neon::Scratch());
    int compared_lean = 0;
    // The SIMD filters take the threshold at run time: every other job at 19, the rest at 5, 10 and 32
    // in rotation (the tail slack ub - minsc - 2 and the mw threshold move with it).
    const int other_minsc[3] = {5, 10, 32};
    // And the bound weights, every fourth job: -O 8 -E 2 (c = 2, so a single-hit diagonal weighs
    // a - c < 0), -A 2 -B 8 -O 12 -E 2 (a = 2) and -x intractg (the tail offset 12), at -k 19.
    const rescue_prune_params weighted[3] = {rescue_prune_params::from(1, 4, 8, 2, 8, 2, kSimdMinsc),
                                             rescue_prune_params::from(2, 8, 12, 2, 12, 2, 2 * kSimdMinsc),
                                             rescue_prune_params::from(1, 9, 16, 1, 16, 1, kSimdMinsc)};
    for (const rescue_prune_params &w : weighted) REQUIRE(w.simd_ok());
    int n_weighted = 0;
    for (const int max_hits : {kDefaultMaxHits, kNoGate}) {
        for (size_t i = 0; i < jobs.size(); i++) {
            const Job &jb = jobs[i];
            const int len1 = (int)jb.ref.size(), len2 = (int)jb.q.size();
            const int minsc = i % 2 ? kSimdMinsc : other_minsc[(i / 2) % 3];
            const bool wtd = i % 4 == 3;
            const rescue_prune_params p = wtd ? weighted[(i / 4) % 3] : rescue_prune_params::defaults(minsc);
            n_weighted += wtd;
            CAPTURE(jb.tag); CAPTURE(i); CAPTURE(len1); CAPTURE(len2); CAPTURE(max_hits); CAPTURE(p.minsc);
            CAPTURE(p.a); CAPTURE(p.c);
            int hb = -2, he = -2, shb = -2, she = -2;
            const int kind = rescue_prune_window(jb.ref.data(), len1, jb.q.data(), len2, p, max_hits, &hb, &he);
            const int skind = rescue_prune_window_scalar(jb.ref.data(), len1, jb.q.data(), len2, p, max_hits,
                                                         *scalar, &shb, &she);
            CHECK(kind == skind);
            if (kind == RESCUE_PRUNE_B2 && skind == RESCUE_PRUNE_B2) {
                CHECK(hb == shb);
                CHECK(he == she);
            }
            // lean() is the int16 reference the SIMD filters were derived from; it has no hit gate,
            // no N guard, a fixed capacity and int16 sums, so compare it only without the gate, on
            // N-free jobs (the filters refuse the rest, checked above), where it fits, and below
            // 32000 hits.
            if (max_hits == kNoGate && !wtd && !has_n(jb.ref) && !has_n(jb.q)
                && len1 + kswv_query_quantum8(len2) + 1 <= 4096 && len1 >= 5 && len2 >= 5
                && count_hits(jb) <= 32000) {
                const rescue_prune_neon::Job lj{len1, len2, 0, 0, -1, -1, jb.ref.data(), jb.q.data()};
                int lhb = -2, lhe = -2;
                const rescue_prune_neon::Kind lk = rescue_prune_neon::lean(lj, *ref16, lhb, lhe, minsc);
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
    CHECK(n_weighted > 50);
#endif
}

// The SIMD filter under the --meth relation (the Rel instantiation of lean_neon_core / lean_x86_core)
// against the scalar filter: mates mostly unconverted (a TAPS-like 5 % conversion), each job under its
// hypothesis, then the other one and its own again on the same bytes (the query table and the repeat
// memo are keyed on relx, so the middle call must not be answered from the first), at the default
// gate and with it open; and relx-rich mates whose entries exceed REL_ECAP, which both filters must
// send to the full window. A decision counts as the SIMD filter's when its view carries bnd16.
TEST_CASE("rescue prune: the SIMD filter decides exactly as the scalar filter under the --meth relation"
          * doctest::test_suite("unit/pair")) {
#if !RESCUE_PRUNE_HAVE_SIMD
    MESSAGE("skipped: no SIMD rescue-prune filter in this build (aarch64, or x86 with AVX2)");
    return;
#else
    std::mt19937 rng(424242);
    auto jobs = build_jobs(rng);
    std::unique_ptr<rescue_prune_scratch> scalar(new rescue_prune_scratch());
    int simd_b2 = 0, over = 0;
    auto check = [&](const Job &jb, int hyp, int max_hits) {
        const int len1 = (int)jb.ref.size(), len2 = (int)jb.q.size();
        rescue_prune_params p = rescue_prune_params::defaults(kSimdMinsc);
        p.set_meth_rel(hyp);
        CAPTURE(jb.tag); CAPTURE(len1); CAPTURE(len2); CAPTURE(hyp); CAPTURE(max_hits);
        int hb = -2, he = -2, shb = -2, she = -2;
        rescue_prune_view v;
        const int kind = rescue_prune_window(jb.ref.data(), len1, jb.q.data(), len2, p, max_hits, &hb, &he, &v);
        const int skind = len1 < 5 || len2 < 5 ? RESCUE_PRUNE_FULL
                                               : rescue_prune_window_scalar(jb.ref.data(), len1, jb.q.data(), len2,
                                                                            p, max_hits, *scalar, &shb, &she);
        CHECK(kind == skind);
        if (kind == RESCUE_PRUNE_B2 && skind == RESCUE_PRUNE_B2) {
            CHECK(hb == shb);
            CHECK(he == she);
            simd_b2 += v.bnd16 != nullptr;
        }
        return skind;
    };
    for (size_t i = 0; i < jobs.size(); i++) {
        Job jb = jobs[i];
        const int hyp = (int)(i & 1);
        rescue_convert_mate(jb.q, hyp, 20, rng);
        const int mh = i % 3 ? kDefaultMaxHits : kNoGate;
        check(jb, hyp, mh);
        check(jb, !hyp, mh);
        check(jb, hyp, mh);
    }
    // relx-rich mates: every relx position doubles the entries of the 5-mers covering it, so a long
    // poly-relx stretch passes REL_ECAP.
    for (const int hyp : {0, 1}) {
        Job jb = planted(rng, 600, 1000, 30);
        for (size_t j = 0; j < jb.q.size(); j++)
            if (rng() % 10) jb.q[j] = (uint8_t)(hyp ? 3 : 0);   // ~90 %: about 25k entries
        over += check(jb, hyp, kNoGate) == RESCUE_PRUNE_FULL;
    }
    CHECK(simd_b2 > 0);
    CHECK(over == 2);
#endif
}

// The SIMD filter's two stateful shortcuts, against the scalar filter:
//   - the hit gate, decided from the per-code occurrence counts before any hit is accumulated:
//     at max_hits = hits - 1 the job must go to the full window, at max_hits = hits it must be
//     decided exactly as the scalar filter decides it;
//   - the per-query cache (query table, presence map, occurrence chains), rebuilt only when the
//     query bytes change: one query run against several windows in a row, then a different
//     query of the same length, must be decided exactly as the scalar filter decides each.
TEST_CASE("rescue prune: the SIMD filter's hit gate and query cache match the scalar filter"
          * doctest::test_suite("unit/pair")) {
#if !RESCUE_PRUNE_HAVE_SIMD
    MESSAGE("skipped: no SIMD rescue-prune filter in this build (aarch64, or x86 with AVX2)");
    return;
#else
    std::mt19937 rng(24601);
    const auto jobs = build_jobs(rng);
    std::unique_ptr<rescue_prune_scratch> scalar(new rescue_prune_scratch());
    auto same_as_scalar = [&](const std::vector<uint8_t> &ref, const std::vector<uint8_t> &q, int max_hits) {
        int hb = -2, he = -2, shb = -2, she = -2;
        const int kind = rescue_prune_window(ref.data(), (int)ref.size(), q.data(), (int)q.size(), kSimdMinsc,
                                             max_hits, &hb, &he);
        const int skind = rescue_prune_window_scalar(ref.data(), (int)ref.size(), q.data(), (int)q.size(),
                                                     kSimdMinsc, max_hits, *scalar, &shb, &she);
        CHECK(kind == skind);
        if (kind == RESCUE_PRUNE_B2 && skind == RESCUE_PRUNE_B2) {
            CHECK(hb == shb);
            CHECK(he == she);
        }
        return kind;
    };
    int gated = 0, at_limit = 0;
    for (size_t i = 0; i < jobs.size(); i++) {
        const Job &jb = jobs[i];
        if (has_n(jb.ref) || has_n(jb.q)) continue;
        const long hits = count_hits(jb);
        if (hits == 0 || hits > 32000) continue;
        CAPTURE(jb.tag); CAPTURE(i); CAPTURE(hits);
        int hb = -2, he = -2;
        CHECK(rescue_prune_window(jb.ref.data(), (int)jb.ref.size(), jb.q.data(), (int)jb.q.size(), kSimdMinsc,
                                  (int)hits - 1, &hb, &he) == RESCUE_PRUNE_FULL);
        ++gated;
        at_limit += same_as_scalar(jb.ref, jb.q, (int)hits) != RESCUE_PRUNE_FULL;
    }
    MESSAGE("hit gate: " << gated << " jobs at hits - 1, " << at_limit << " decided (not FULL) at hits");
    CHECK(gated > 50);
    CHECK(at_limit > 20);
    // Query cache: each job's query against the next few windows, in a row.
    int cached = 0;
    for (size_t i = 0; i + 4 < jobs.size(); i += 5) {
        const auto &q = jobs[i].q;
        for (size_t k = i; k < i + 4; k++) { same_as_scalar(jobs[k].ref, q, kNoGate); ++cached; }
        // Same length, one base changed: the cache must notice and rebuild.
        auto q2 = q;
        q2[q2.size() / 2] = (uint8_t)((q2[q2.size() / 2] + 1) & 3);
        same_as_scalar(jobs[i].ref, q2, kNoGate);
        same_as_scalar(jobs[i].ref, q, kNoGate);
    }
    CHECK(cached > 20);
#endif
}

// The two shortcuts band planning takes from the SIMD filter (rescue_prune_neon.h / rescue_prune_x86.h):
//   - an exact repeat of the previous job returns the previous decision without recomputing it, so
//     it must return the same (kind, hb, he) and a view with the same components, and the same job
//     at another threshold must be decided afresh, as the scalar filter decides it;
//   - the filter lists its components at the call's threshold, and rescue_band_components reuses the list for the
//     whole view instead of rescanning it. With the list hidden (ncomp = -1) it rescans; both must
//     give the same components and the same success at every cap, including a window with more
//     components than the filter stores (COMP_CAP) and an output vector that is not empty.
TEST_CASE("rescue prune: the SIMD filter's repeat memo and component list match a rescan"
          * doctest::test_suite("unit/pair")) {
#if !RESCUE_PRUNE_HAVE_SIMD
    MESSAGE("skipped: no SIMD rescue-prune filter in this build (aarch64, or x86 with AVX2)");
    return;
#else
    std::mt19937 rng(314159);
    auto jobs = build_jobs(rng);
    std::unique_ptr<rescue_prune_scratch> scalar(new rescue_prune_scratch());
    {   // one window with more components than COMP_CAP: a 20-base mate segment planted every 26 rows
        // over nearly the SIMD filter's capacity (about 150 components), from its own generator so
        // the job set above can grow without changing this window
        std::mt19937 crng(271828);
        const int len1 = 3950, len2 = 100;
        auto ref = random_bases(crng, len1);
        auto q = random_bases(crng, len2);
        for (int at = 0; at + 20 <= len1; at += 26) std::copy(q.begin() + 40, q.begin() + 60, ref.begin() + at);
        jobs.push_back({std::move(ref), std::move(q), "many components"});
    }
    const int comp_cap = rescue_prune_simd_scratch_t::COMP_CAP;
    auto same = [](const std::vector<rb_comp> &x, const std::vector<rb_comp> &y) {
        if (x.size() != y.size()) return false;
        for (size_t c = 0; c < x.size(); c++)
            if (x[c].ub != y[c].ub || x[c].i0 != y[c].i0 || x[c].dlo != y[c].dlo || x[c].dhi != y[c].dhi
                || x[c].dmaxhit != y[c].dmaxhit)
                return false;
        return true;
    };
    int n_b2 = 0, n_over_cap = 0;
    std::vector<rb_comp> listed, rescanned;
    for (size_t i = 0; i < jobs.size(); i++) {
        const Job &jb = jobs[i];
        const int len1 = (int)jb.ref.size(), len2 = (int)jb.q.size();
        CAPTURE(jb.tag); CAPTURE(i); CAPTURE(len1); CAPTURE(len2);
        int hb = -2, he = -2, rhb = -2, rhe = -2;
        rescue_prune_view v, rv;
        const int kind = rescue_prune_window(jb.ref.data(), len1, jb.q.data(), len2, kSimdMinsc, kNoGate,
                                             &hb, &he, &v);
        // The repeat: same bytes from a different buffer, so only the contents can match.
        const std::vector<uint8_t> ref2(jb.ref), q2(jb.q);
        const int rkind = rescue_prune_window(ref2.data(), len1, q2.data(), len2, kSimdMinsc, kNoGate,
                                              &rhb, &rhe, &rv);
        CHECK(rkind == kind);
        CHECK(rhb == hb);
        CHECK(rhe == he);
        // The view reports the repeat (a caller reuses the job's result on it) and only the repeat;
        // a job the NEON filter handed to the scalar one (past its int16 range) has no memo.
        if (v.bnd16) CHECK(rv.repeat);
        CHECK(!v.repeat);
        if (i % 8 == 0) {   // the same bytes at another threshold: not a repeat
            const uint64_t hits0 = rescue_prune_memo_hits();
            int ohb = -2, ohe = -2, shb = -2, she = -2;
            const int okind = rescue_prune_window(ref2.data(), len1, q2.data(), len2, 25, kNoGate, &ohb, &ohe);
            CHECK(rescue_prune_memo_hits() == hits0);
            const int skind = rescue_prune_window_scalar(ref2.data(), len1, q2.data(), len2, 25, kNoGate,
                                                         *scalar, &shb, &she);
            CHECK(okind == skind);
            if (okind == RESCUE_PRUNE_B2 && skind == RESCUE_PRUNE_B2) { CHECK(ohb == shb); CHECK(ohe == she); }
            // restore the view at kSimdMinsc for the component checks below
            rescue_prune_window(ref2.data(), len1, q2.data(), len2, kSimdMinsc, kNoGate, &rhb, &rhe, &rv);
        }
        if (kind != RESCUE_PRUNE_B2 || !v.bnd16) continue;
        ++n_b2;
        REQUIRE(v.ncomp >= 0);
        CHECK(rv.ncomp == v.ncomp);
        n_over_cap += v.ncomp > comp_cap;
        rescue_prune_view hidden = v;
        hidden.ncomp = -1;
        for (const int cap : {1 << 20, 0, 1, 2, comp_cap - 1, comp_cap, comp_cap + 1, v.ncomp - 1, v.ncomp,
                              v.ncomp + 1}) {
            for (const int prefill : {0, 3}) {
                CAPTURE(cap); CAPTURE(prefill);
                listed.assign((size_t)prefill, rb_comp{0, 0, 0, 0, 0});
                rescanned = listed;
                const bool ok_l = rescue_band_components(rv, kSimdMinsc, 0, rv.nd, listed, cap);
                const bool ok_r = rescue_band_components(hidden, kSimdMinsc, 0, hidden.nd, rescanned, cap);
                CHECK(ok_l == ok_r);
                CHECK(same(listed, rescanned));
            }
        }
    }
    MESSAGE("component list: " << n_b2 << " B2 jobs, " << n_over_cap << " with more than COMP_CAP components");
    CHECK(n_b2 > 50);
    CHECK(n_over_cap > 0);
#endif
}

TEST_CASE("rescue band: the pass-0 cost gate's default follows the kswv tier"
          * doctest::test_suite("unit/pair")) {
    // Where kswv sweeps 64 lanes (the AVX-512BW tier) the 32-lane band kernel cannot undercut the
    // hull, so the default bands no pass-0 parent there; every other tier keeps the 85 % margin
    // (rescue_band.cpp, rb_cost_pct). BWA3_RESCUE_BAND_COST overrides either.
    CHECK(rescue_band_cost_pct_default(BWAMEM3_TIER_AVX512BW) == 0);
    for (const int tier : {BWAMEM3_TIER_NONE, BWAMEM3_TIER_SSE41, BWAMEM3_TIER_SSE42, BWAMEM3_TIER_AVX,
                           BWAMEM3_TIER_AVX2, BWAMEM3_TIER_NEON}) {
        CAPTURE(tier);
        CHECK(rescue_band_cost_pct_default(tier) == 85);
    }
}

// The run-level cost gate (rescue_prune_cost_ok): every outcome gives the same output, so no identity
// check can see a wrong one; pinned here per architecture. aarch64 prunes everywhere but --meth with
// chemistry other than EM-seq (TAPS) on converted copies, which prunes under the relation; x86 never
// under --meth, relation or not, only where the SIMD filter takes the
// weights (K = 5, a <= 16), and not at the AVX-512BW tier from seed length 25 (minsc >= 25 a), so
// -A 2 at the default -k 19 (minsc 38) still prunes there and -A 2 -k 25 (minsc 50) does not.
TEST_CASE("rescue prune: the cost gate per architecture, --meth chemistry and kswv tier"
          * doctest::test_suite("unit/pair")) {
    const rescue_prune_params dflt = rescue_prune_params::defaults(19);
    const rescue_prune_params k25 = rescue_prune_params::defaults(25);
    const rescue_prune_params a17 = rescue_prune_params::from(17, 68, 102, 17, 102, 17, 19 * 17);
    REQUIRE(dflt.valid);
    REQUIRE(k25.valid);
    const rescue_prune_params a2 = rescue_prune_params::from(2, 8, 12, 2, 12, 2, 19 * 2);
    const rescue_prune_params a2k25 = rescue_prune_params::from(2, 8, 12, 2, 12, 2, 25 * 2);
    REQUIRE(a17.valid);
    REQUIRE(!a17.simd_ok());
    REQUIRE(a2.simd_ok());
    REQUIRE(a2k25.simd_ok());
    // meth_fits: the run's --meth filter fits its reads (EM-seq on converted copies, or the relation).
    for (const bool avx512 : {false, true}) {
        CAPTURE(avx512);
#if defined(__aarch64__)
        CHECK(rescue_prune_cost_ok(dflt, false, false, avx512));
        CHECK(rescue_prune_cost_ok(dflt, true, true, avx512));    // --meth: EM-seq, or TAPS under the relation
        CHECK(!rescue_prune_cost_ok(dflt, true, false, avx512));  // --meth=taps on converted copies
        CHECK(rescue_prune_cost_ok(k25, false, false, avx512));
        CHECK(rescue_prune_cost_ok(a17, false, false, avx512));   // the scalar filter pays here
        CHECK(rescue_prune_cost_ok(a2k25, false, false, avx512));
#else
        CHECK(rescue_prune_cost_ok(dflt, false, false, avx512));
        CHECK(!rescue_prune_cost_ok(dflt, true, true, avx512));   // no --meth pruning on x86, fitting or not
        CHECK(!rescue_prune_cost_ok(dflt, true, false, avx512));
        CHECK(rescue_prune_cost_ok(k25, false, false, avx512) == !avx512);
        CHECK(rescue_prune_cost_ok(rescue_prune_params::defaults(24), false, false, avx512));
        CHECK(!rescue_prune_cost_ok(a17, false, false, avx512));  // the scalar filter would decide
        CHECK(rescue_prune_cost_ok(a2, false, false, avx512));    // -A 2 -k 19: the seed length gates
        CHECK(rescue_prune_cost_ok(a2k25, false, false, avx512) == !avx512);
#endif
    }
}

// The default hit gate (rescue_prune_max_hits_default): like the cost gate, every value gives the
// same output, so it is pinned here. aarch64 takes 1000 only where banding can turn the pruned
// windows into savings (banding on, and not a --meth run that leaves its pruned windows unbanded);
// everything else, and all of x86, takes 400.
TEST_CASE("rescue prune: the default hit gate per architecture, banding and --meth"
          * doctest::test_suite("unit/pair")) {
    for (const bool banding : {false, true})
        for (const bool meth : {false, true}) {
            CAPTURE(banding);
            CAPTURE(meth);
#if defined(__aarch64__)
            CHECK(rescue_prune_max_hits_default(banding, meth) == (banding && !meth ? 1000 : kDefaultMaxHits));
#else
            CHECK(rescue_prune_max_hits_default(banding, meth) == kDefaultMaxHits);
#endif
        }
}

// kswv's 8-bit kernels load each gap type's open plus extend as a byte, so a sum past 255 wraps there
// and neither pruning nor the band kernels may take such a scoring (kswv8_scoring_ok).
TEST_CASE("rescue prune: scorings past kswv's 8-bit gap byte are refused"
          * doctest::test_suite("unit/pair")) {
    CHECK(kswv8_scoring_ok(1, 4, 249, 6, 249, 6));
    CHECK(rescue_prune_params::from(1, 4, 249, 6, 249, 6, 19).valid);
    for (const int s : {0, 1}) {
        CAPTURE(s);
        const int od = s ? 6 : 250, ed = s ? 1 : 6, oi = s ? 250 : 6, ei = s ? 6 : 1;
        CHECK(!kswv8_scoring_ok(1, 4, od, ed, oi, ei));
        CHECK(!rescue_prune_params::from(1, 4, od, ed, oi, ei, 19).valid);
        rb_scoring r;
        r.o_del = od; r.e_del = ed; r.o_ins = oi; r.e_ins = ei;
        CHECK(!r.valid());
    }
}
