// test/unit/test_ksw_global2_oracle.cpp
//
// ksw_global2 (src/ksw.cpp) checked against an independent oracle, on every
// SIMD tier this host can run. `mem` reaches ksw_global2 through
// bwa_gen_cigar2/bwa_gen_cigar3 for every record it emits: it is what turns an
// alignment region into the CIGAR, NM and MD of the SAM line.
//
// The existing ksw_global2 tests (test_ksw_global2_wave.cpp,
// test_ksw_global2_scalar.cpp) are differential: each fast path is held
// byte-identical to the in-tree scalar reference, so a fault the reference
// shares would pass them. This file adds the other half: a textbook
// full-matrix Gotoh DP written here from the scoring model ksw.h and ksw.cpp
// document, not from their code. The model:
//
//   - global: all of query q[0,qlen) against all of target t[0,tlen);
//   - a substitution scores mat[t*5+q] (bwa_fill_scmat: a on a match, -b on a
//     mismatch, -1 against N);
//   - a run of k deleted target bases costs o_del + k*e_del, a run of k
//     inserted query bases o_ins + k*e_ins;
//   - a deletion never directly follows an insertion or the reverse (ksw.cpp:
//     E and F open only from M, so "10M3I3D10M" is not a global() CIGAR).
//
// With a band w >= max(qlen, tlen) the band excludes no cell, so the returned
// score must equal the oracle's optimum. The CIGAR is one of possibly several
// optimal paths, so it is checked by rescoring it under the same model: it must
// consume exactly qlen and tlen, use only M/I/D, never put I next to D, and
// score exactly the returned value. The hand-worked cases at the top pin the
// CIGAR and the direction of each gap's price, on inputs with a single optimal
// path. The int16 wavefront takes every small-score pair, so the random cases
// reach it; a high-gap case rejects int16 and checks the int32 wavefront
// directly against the oracle (score, full rescored CIGAR, and the wave
// counters). Both are also held to the scalar reference in
// test_ksw_global2_wave.cpp.
//
// Tiers: x86 calls the per-tier kernels by name (ksw_global2_avx2 and, where
// the host has AVX-512BW, ksw_global2_avx512bw; libbwa.a carries every tier's
// objects on every x86 build), so one unit binary covers each tier the runner
// has; arm64 calls the one NEON build. A tier the host lacks is reported with
// MESSAGE and skipped, never failed. The sub-AVX2 tiers are left out: `mem`
// refuses them at run time. On a tier with wavefront kernels the random cases
// must also have entered them, so the comparison cannot pass by running only
// the scalar path.
//
// Random inputs use std::mt19937 raw output (rng() % n), which the standard
// fixes; std::uniform_int_distribution differs between libstdc++ and libc++.

#include "doctest/doctest.h"
#include "scoring.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <random>
#include <string>
#include <vector>

#include "ksw.h"

#if defined(__x86_64__) || defined(__i386__)
// The per-tier builds of the kernel (src/kernel_dispatch.h renames each tier's
// symbols with its suffix); ksw.h declares only the unsuffixed dispatcher.
extern "C" {
int ksw_global2_avx2(int, const uint8_t *, int, const uint8_t *, int, const int8_t *, int, int, int, int, int,
                     int *, uint32_t **);
int ksw_global2_avx512bw(int, const uint8_t *, int, const uint8_t *, int, const int8_t *, int, int, int, int,
                         int, int *, uint32_t **);
unsigned long ksw_g2_wave_exec_count_avx2(void);
unsigned long ksw_g2_wave_exec_count_avx512bw(void);
unsigned long ksw_g2_wave16_exec_count_avx2(void);
unsigned long ksw_g2_wave16_exec_count_avx512bw(void);
int ksw_g2_wave_wmin_avx2(void);
int ksw_g2_wave_wmin_avx512bw(void);
int ksw_g2_wave16_wmin_avx2(void);
int ksw_g2_wave16_wmin_avx512bw(void);
}
#endif

namespace {

extern "C" {
using G2Fn = int (*)(int, const uint8_t *, int, const uint8_t *, int, const int8_t *, int, int, int, int,
                     int, int *, uint32_t **);
using CountFn = unsigned long (*)(void);
using WminFn = int (*)(void);
}

struct Tier {
    const char *name;
    G2Fn g2;
    CountFn wave_count;    // int16 and int32 wavefront entries
    CountFn wave16_count;  // int16 entries only
    WminFn wave_wmin;
    WminFn wave16_wmin;
};

std::vector<Tier> runnable_tiers() {
    std::vector<Tier> ts;
#if defined(__x86_64__) || defined(__i386__)
    __builtin_cpu_init();
    if (__builtin_cpu_supports("avx2"))
        ts.push_back({"avx2", ksw_global2_avx2, ksw_g2_wave_exec_count_avx2, ksw_g2_wave16_exec_count_avx2,
                      ksw_g2_wave_wmin_avx2, ksw_g2_wave16_wmin_avx2});
    else
        MESSAGE("ksw_global2 oracle: no AVX2 on this host; the avx2 tier is skipped");
    if (__builtin_cpu_supports("avx512bw"))
        ts.push_back({"avx512bw", ksw_global2_avx512bw, ksw_g2_wave_exec_count_avx512bw,
                      ksw_g2_wave16_exec_count_avx512bw, ksw_g2_wave_wmin_avx512bw, ksw_g2_wave16_wmin_avx512bw});
    else
        MESSAGE("ksw_global2 oracle: no AVX-512BW on this host; the avx512bw tier is skipped");
#else
    ts.push_back({"neon", ksw_global2, ksw_g2_wave_exec_count, ksw_g2_wave16_exec_count, ksw_g2_wave_wmin,
                  ksw_g2_wave16_wmin});
#endif
    return ts;
}

struct Scoring {
    int a, b, o_del, e_del, o_ins, e_ins;
    int8_t mat[25];
};

// bwa_fill_scmat's matrix, written out: a on the diagonal, -b off it, -1 for N.
Scoring scoring(int a, int b, int o_del, int e_del, int o_ins, int e_ins) {
    Scoring s{a, b, o_del, e_del, o_ins, e_ins, {}};
    const bwa_tests::ScoringMatrix m = bwa_tests::build_scoring_matrix(a, b, 1);
    std::copy(m.begin(), m.end(), s.mat);
    return s;
}

// Full-matrix three-state Gotoh: M ends in a substitution, E in a deletion
// (target base against nothing), F in an insertion. E and F open only from M.
int oracle_score(const std::vector<uint8_t> &q, const std::vector<uint8_t> &t, const Scoring &s) {
    const int NEG = -(1 << 28);
    const size_t Q = q.size(), T = t.size();
    std::vector<int> M((T + 1) * (Q + 1), NEG), E(M), F(M);
    auto at = [Q](size_t i, size_t j) { return i * (Q + 1) + j; };
    M[at(0, 0)] = 0;
    for (size_t j = 1; j <= Q; ++j) F[at(0, j)] = -(s.o_ins + (int)j * s.e_ins);
    for (size_t i = 1; i <= T; ++i) E[at(i, 0)] = -(s.o_del + (int)i * s.e_del);
    for (size_t i = 1; i <= T; ++i) {
        for (size_t j = 1; j <= Q; ++j) {
            const int h = std::max(M[at(i - 1, j - 1)], std::max(E[at(i - 1, j - 1)], F[at(i - 1, j - 1)]));
            M[at(i, j)] = h + s.mat[t[i - 1] * 5 + q[j - 1]];
            E[at(i, j)] = std::max(M[at(i - 1, j)] - s.o_del - s.e_del, E[at(i - 1, j)] - s.e_del);
            F[at(i, j)] = std::max(M[at(i, j - 1)] - s.o_ins - s.e_ins, F[at(i, j - 1)] - s.e_ins);
        }
    }
    return std::max(M[at(T, Q)], std::max(E[at(T, Q)], F[at(T, Q)]));
}

struct G2 {
    int score = 0;
    std::vector<uint32_t> cigar;
};

G2 run(G2Fn fn, const std::vector<uint8_t> &q, const std::vector<uint8_t> &t, const Scoring &s, int w) {
    int n_cigar = 0;
    uint32_t *cigar = nullptr;
    G2 r;
    r.score = fn((int)q.size(), q.data(), (int)t.size(), t.data(), 5, s.mat, s.o_del, s.e_del, s.o_ins,
                 s.e_ins, w, &n_cigar, &cigar);
    REQUIRE((n_cigar > 0) == (cigar != nullptr));
    if (n_cigar > 0) r.cigar.assign(cigar, cigar + n_cigar);
    free(cigar);
    return r;
}

std::string cigar_string(const std::vector<uint32_t> &cigar) {
    std::string out;
    for (uint32_t c : cigar) {
        const uint32_t op = c & 0xf;
        out += std::to_string(c >> 4) + (op < 9 ? "MIDNSHP=X"[op] : '?');
    }
    return out;
}

// Score a CIGAR under the model; false (and a FAIL) if it is not a valid
// global() path for q and t.
bool rescore(const std::vector<uint32_t> &cigar, const std::vector<uint8_t> &q, const std::vector<uint8_t> &t,
             const Scoring &s, int *score) {
    size_t i = 0, j = 0;  // target, query
    int sc = 0, prev = -1;
    for (uint32_t c : cigar) {
        const int op = (int)(c & 0xf), len = (int)(c >> 4);
        if (len <= 0 || (op != 0 && op != 1 && op != 2)) return false;
        if ((op == 1 && prev == 2) || (op == 2 && prev == 1)) return false;
        if (op == 0) {
            for (int k = 0; k < len; ++k, ++i, ++j) {
                if (i >= t.size() || j >= q.size()) return false;
                sc += s.mat[t[i] * 5 + q[j]];
            }
        } else if (op == 1) {
            j += (size_t)len;
            sc -= s.o_ins + len * s.e_ins;
        } else {
            i += (size_t)len;
            sc -= s.o_del + len * s.e_del;
        }
        prev = op;
    }
    *score = sc;
    return i == t.size() && j == q.size();
}

std::vector<uint8_t> random_seq(std::mt19937 &rng, size_t n, bool with_n) {
    std::vector<uint8_t> v(n);
    for (size_t k = 0; k < n; ++k) v[k] = (uint8_t)(with_n && rng() % 50 == 0 ? 4 : rng() % 4);
    return v;
}

// t derived from q by substitutions and short indels, so most pairs align
// with gaps rather than as noise.
std::vector<uint8_t> mutate(std::mt19937 &rng, const std::vector<uint8_t> &q) {
    std::vector<uint8_t> t;
    for (size_t k = 0; k < q.size(); ++k) {
        const uint32_t r = rng() % 100;
        if (r < 5) {
            t.push_back((uint8_t)((q[k] + 1 + rng() % 3) % 4));
        } else if (r < 8) {
            continue;  // deleted from t: an insertion in q
        } else if (r < 11) {
            const uint32_t n = 1 + rng() % 4;
            for (uint32_t m = 0; m < n; ++m) t.push_back((uint8_t)(rng() % 4));
            t.push_back(q[k]);
        } else {
            t.push_back(q[k]);
        }
    }
    if (t.empty()) t.push_back(0);
    return t;
}

}  // namespace

TEST_CASE("ksw_global2 oracle: hand-worked alignments on every runnable tier"
          * doctest::test_suite("unit/ksw_global2")) {
    // q is 40 bases with no base equal to its neighbour at the indel sites, so
    // each case below has exactly one optimal path. Each case runs on every
    // tier (the tier loop sits inside the case, not around it, so doctest
    // enters every case once per run).
    const std::vector<uint8_t> q = {0, 1, 2, 3, 1, 0, 3, 2, 0, 2, 1, 3, 0, 1, 2, 0, 3, 1, 2, 0,
                                    1, 3, 2, 0, 1, 2, 3, 0, 2, 1, 0, 3, 1, 0, 2, 3, 0, 1, 3, 2};
    const Scoring def = scoring(1, 4, 6, 1, 6, 1);
    const std::vector<Tier> tiers = runnable_tiers();
    if (tiers.empty()) return;  // no tier to check; runnable_tiers() reported each one skipped

    SUBCASE("three target bases the query lacks: one deletion, o_del + 3*e_del") {
        // q[19] = 0 (A) and q[20] = 1 (C); inserting G G G between them in
        // the target leaves no shifted placement of the gap.
        std::vector<uint8_t> t(q.begin(), q.begin() + 20);
        t.insert(t.end(), {2, 2, 2});
        t.insert(t.end(), q.begin() + 20, q.end());
        for (const Tier &tier : tiers) {
            CAPTURE(tier.name);
            const G2 r = run(tier.g2, q, t, def, 43);
            CHECK(r.score == 40 - (6 + 3 * 1));
            CHECK(cigar_string(r.cigar) == "20M3D20M");
            // Under -O 6,5 -E 1,2 the deletion keeps its own price, 6 + 3*1
            // (the insertion price would be 5 + 3*2 = 11, giving 29).
            const G2 ra = run(tier.g2, q, t, scoring(1, 4, 6, 1, 5, 2), 43);
            CHECK(ra.score == 40 - (6 + 3 * 1));
            CHECK(cigar_string(ra.cigar) == "20M3D20M");
        }
    }
    SUBCASE("two query bases the target lacks: one insertion, o_ins + 2*e_ins") {
        const Scoring asym = scoring(1, 4, 6, 1, 5, 2);
        std::vector<uint8_t> qi(q.begin(), q.begin() + 20);
        qi.insert(qi.end(), {2, 2});
        qi.insert(qi.end(), q.begin() + 20, q.end());
        for (const Tier &tier : tiers) {
            CAPTURE(tier.name);
            const G2 r = run(tier.g2, qi, q, asym, 42);
            CHECK(r.score == 40 - (5 + 2 * 2));
            CHECK(cigar_string(r.cigar) == "20M2I20M");
        }
    }
}

TEST_CASE("ksw_global2 oracle: random pairs score the Gotoh optimum on every runnable tier"
          * doctest::test_suite("unit/ksw_global2")) {
    const Scoring scorings[] = {
        scoring(1, 4, 6, 1, 6, 1),   // bwa-mem3 defaults
        scoring(1, 4, 5, 2, 7, 1),   // -O 5,7 -E 2,1: deletions and insertions priced apart
        scoring(2, 8, 12, 2, 12, 2), // -A 2, which scales -B -O -E with it
    };
    const std::vector<Tier> tiers = runnable_tiers();
    if (tiers.empty()) return;  // no tier to check; runnable_tiers() reported each one skipped
    for (const Tier &tier : tiers) {
        CAPTURE(tier.name);
        const unsigned long wave_before = tier.wave_count();
        std::mt19937 rng(20261024u);
        for (const Scoring &s : scorings) {
            for (int n = 0; n < 40; ++n) {
                const size_t qlen = 1 + rng() % 120;
                const std::vector<uint8_t> q = random_seq(rng, qlen, n % 3 == 0);
                const std::vector<uint8_t> t = n % 10 == 9 ? random_seq(rng, 1 + rng() % 120, false) : mutate(rng, q);
                const int w = (int)std::max(q.size(), t.size());
                CAPTURE(n);
                CAPTURE(q.size());
                CAPTURE(t.size());
                const G2 r = run(tier.g2, q, t, s, w);
                CHECK(r.score == oracle_score(q, t, s));
                int sc = 0;
                const bool valid = rescore(r.cigar, q, t, s, &sc);
                CAPTURE(cigar_string(r.cigar));
                CHECK(valid);
                CHECK(sc == r.score);
            }
        }
        // A tier with a wavefront kernel must have run it on some of these
        // (every w >= its crossover goes there); otherwise this compared only
        // the scalar path.
        if (tier.wave_wmin() > 0 || tier.wave16_wmin() > 0) CHECK(tier.wave_count() > wave_before);
    }
}

TEST_CASE("ksw_global2 oracle: high gap opens reject int16, run the int32 wavefront and score the Gotoh optimum"
          * doctest::test_suite("unit/ksw_global2")) {
    // A gap open of 3000 puts o + e past the int16 kernel's sentinel headroom,
    // so ksw_g2_wave16_safe rejects the pair at any length and, at w >= the
    // int32 crossover, the int32 wavefront runs. The other cases here use small
    // scores, which the int16 kernel always takes, so without this case the
    // independent oracle never checks the int32 kernel. Each pair must leave the
    // int16 entry count alone, raise the aggregate wavefront count by one, score
    // the oracle's optimum and return a CIGAR that rescores to that score. The
    // hand-made pairs pin one deletion and one insertion under the first gap
    // set. The second keeps deletions cheap under a still-rejected insertion
    // open, so its random pairs with extra target bases align through the
    // int32 kernel's cheap-deletion path; some must (a D op in a random pair's
    // CIGAR), or a change to the generated pairs left that path unchecked.
    // (With an insertion open of 3000 no optimum takes gap beyond what the
    // length difference forces, so the floor counts D ops, not extra gap.)
    const Scoring big = scoring(1, 4, 3000, 1, 3000, 1);
    const Scoring cheap_del = scoring(1, 4, 5, 2, 3000, 3);
    // ksw_g2_wave16_safe (src/ksw_global2_wave.h, which only ksw.cpp may
    // include) rejects int16 when KSW_VMINF16 - (o_max + e_max) <= -32768 +
    // SLACK, with KSW_VMINF16 = -30000 and SLACK = 64 there: o + e >= 2704.
    // The wave16 counter below shows the rejection directly. These are copies
    // of those constants, so a change there is caught by that counter, not here;
    // the copies only document why this scoring set is rejected.
    const int kVminf16 = -30000, kSlack = 64;
    const int int16_reject_oe = kVminf16 + 32768 - kSlack;   // 2704
    const std::vector<Tier> tiers = runnable_tiers();
    if (tiers.empty()) return;  // no tier to check; runnable_tiers() reported each one skipped
    std::mt19937 rng(6870013u);
    struct P { std::vector<uint8_t> q, t; };
    std::vector<P> pairs;
    {
        const std::vector<uint8_t> q = random_seq(rng, 64, false);
        // Three target bases the query lacks, all c, between q[29] and q[30]:
        // c differs from both, so the gap has one optimal placement.
        uint8_t c = 0;
        while (c == q[29] || c == q[30]) ++c;
        std::vector<uint8_t> t(q.begin(), q.begin() + 30);
        t.insert(t.end(), {c, c, c});
        t.insert(t.end(), q.begin() + 30, q.end());
        pairs.push_back({q, t});
        pairs.push_back({t, q});                             // the same, as an insertion
    }
    for (int n = 0; n < 12; ++n) {
        const std::vector<uint8_t> q = random_seq(rng, 60 + rng() % 60, n % 3 == 0);
        pairs.push_back({q, mutate(rng, q)});
    }
    int n_cheap_del_d = 0, n_tiers = 0;
    for (const Tier &tier : tiers) {
        CAPTURE(tier.name);
        if (tier.wave_wmin() == 0) {
            MESSAGE("ksw_global2 oracle: tier " << tier.name << " has no int32 wavefront; case skipped");
            continue;
        }
        ++n_tiers;
        for (const Scoring *s : {&big, &cheap_del}) {
            REQUIRE(s->o_ins + s->e_ins >= int16_reject_oe);   // the int16 rejection this case depends on
            for (size_t k = 0; k < pairs.size(); ++k) {
                const P &p = pairs[k];
                const int w = std::max({(int)p.q.size(), (int)p.t.size(), tier.wave_wmin()});
                CAPTURE(k);
                CAPTURE(s->o_del);
                CAPTURE(p.q.size());
                CAPTURE(p.t.size());
                const unsigned long agg0 = tier.wave_count(), w16_0 = tier.wave16_count();
                const G2 r = run(tier.g2, p.q, p.t, *s, w);
                CHECK(tier.wave16_count() == w16_0);    // int16 rejected
                CHECK(tier.wave_count() == agg0 + 1);   // the int32 wavefront ran
                CHECK(r.score == oracle_score(p.q, p.t, *s));
                int sc = 0;
                const bool valid = rescore(r.cigar, p.q, p.t, *s, &sc);
                CAPTURE(cigar_string(r.cigar));
                CHECK(valid);
                CHECK(sc == r.score);
                if (s == &cheap_del && k >= 2)
                    n_cheap_del_d += std::any_of(r.cigar.begin(), r.cigar.end(),
                                                 [](uint32_t c) { return (c & 0xf) == 2; });
                if (s == &big && k == 0) CHECK(cigar_string(r.cigar) == "30M3D34M");
                if (s == &big && k == 1) CHECK(cigar_string(r.cigar) == "30M3I34M");
            }
        }
    }
    // Each gap set's floor: under big, the pinned 30M3D34M / 30M3I34M above;
    // under cheap_del, a random pair (not a hand-made one) whose CIGAR has a D.
    if (n_tiers > 0) CHECK(n_cheap_del_d > 0);
}
