// test/unit/test_mem_extend_bsw_oracle.cpp
//
// The banded seed-extension kernels (bandedSWA.cpp: getScores8, getScores16 and
// the scalar fallback scalarBandedSWA) checked against an oracle written here, on
// every SIMD tier this host can run. `mem` scores every seed extension with these
// kernels: the 8-bit kernel for pairs inside its routing envelope, the 16-bit
// kernel for the rest, the scalar kernel for pairs too long for either. Each pair
// returns six fields (score, qle, tle, gscore, gtle, max_off) that decide where
// the alignment ends, whether it is clipped and whether the band is retried, so
// every field is compared, for every pair.
//
// The existing kernel tests (test_bandedswa_*.cpp) are differential: each vector
// kernel is held to scalarBandedSWA, so a fault shared with the scalar kernel, or
// a change made to all of them together, would pass them. This file adds an
// oracle written from the model the kernels document: bwa's ksw_extend2
// (ksw.h's "Extend alignment" comment and the lane-independence note at the top
// of bandedSWA.cpp). The model is an extension, not a local alignment, and
// carries bwa's heuristics, which the oracle implements as documented because
// they change the result:
//
//   - the alignment is anchored at the seed: the cell before the first query
//     and target base scores h0, a cell whose diagonal predecessor scores 0 is
//     dead (no local restart), and the seed row and column are the gap scores
//     max(0, h0 - o - k*e);
//   - a gap opens from the diagonal score M, not from the cell maximum H, so an
//     insertion never directly follows a deletion or the reverse;
//   - the band w is first clamped to the longest insertion and deletion that
//     could pay for themselves, (qlen*max(mat) + end_bonus - o)/e + 1;
//   - row i computes the columns [beg, end) of its band, where beg and end
//     shrink to the nonzero cells of the row before (end may run two columns
//     past the last nonzero one), and a column left of a row's range keeps the
//     value it last held (the seed row where never computed);
//   - the row maximum takes its last column on ties, the overall maximum only a
//     strictly larger row maximum; gscore takes the last row on ties;
//   - the extension stops at the first row whose maximum is 0, or, with
//     zdrop > 0, whose maximum has dropped more than zdrop below the best score
//     after discounting the gap length (e_del or e_ins per base of diagonal
//     drift).
//
// Expected values come from the hand-worked cases at the top (which also check
// the oracle), and from the oracle on random and adversarial pairs: mutated,
// homopolymer, all-mismatch, long insertions and deletions, N bases, seeds at
// the 8-bit kernel's score ceiling, short queries against long targets with a
// large seed (where the band clamp decides gscore), and pairs placed at the exact
// z-drop and band width where their result changes. Each pair is scored in
// batches of several sizes and orders, so a result that depends on its lane
// group, lane-compaction move or batch position fails.
//
// The 8-bit kernel is given only pairs inside the routing envelope documented at
// bsw8_envelope_ok (bwamem.cpp): len1 >= len2, both < MAX_SEQ_LEN8, w <= 124,
// zdrop + max(a, 1) <= 253 and h0 + min(len1, len2)*a < 255 - max(a, 1). The
// 16-bit kernel is given every pair.
//
// Tiers: x86 builds each tier's kernel by name (make_bsw_kernel_avx2 and, where
// the host has AVX-512BW, make_bsw_kernel_avx512bw; libbwa.a carries every
// tier's objects on every x86 build), as constructed for `mem`, lane compaction
// at the tier's default included. arm64 builds the one NEON kernel. A tier the
// host lacks is reported with MESSAGE and skipped. The sub-AVX2 tiers are left
// out: `mem` refuses them at run time.
//
// Random inputs use std::mt19937 raw output (rng() % n), which the standard
// fixes; std::uniform_int_distribution differs between libstdc++ and libc++.

#include "doctest/doctest.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "bandedSWA.h"
#include "bsw_batch.h"
#include "ext_ladder.h"
#include "ext_result.h"

#if HAVE_BSW_VECTOR_8_16

namespace {

using bwa_tests::ExtPair;
using bwa_tests::ExtResult;

struct Scoring {
    int a, b, o_del, e_del, o_ins, e_ins, zdrop, end_bonus;
    int8_t mat[25];
};

// bwa_fill_scmat's matrix written out (a on the diagonal, -b off it, -1 for N);
// `ot` also scores reference C against read T as a match, the --meth OT matrix.
Scoring scoring(int a, int b, int o_del, int e_del, int o_ins, int e_ins, int zdrop, int end_bonus,
                bool ot = false) {
    Scoring s{a, b, o_del, e_del, o_ins, e_ins, zdrop, end_bonus, {}};
    for (int i = 0; i < 5; ++i)
        for (int j = 0; j < 5; ++j) s.mat[i * 5 + j] = (int8_t)(i == 4 || j == 4 ? -1 : i == j ? a : -b);
    if (ot) s.mat[1 * 5 + 3] = (int8_t)a;
    return s;
}

std::string show(const ExtResult &r) {
    return "score " + std::to_string(r.score) + " qle " + std::to_string(r.qle) + " tle " +
           std::to_string(r.tle) + " gscore " + std::to_string(r.gscore) + " gtle " + std::to_string(r.gtle) +
           " max_off " + std::to_string(r.max_off);
}

// The ksw_extend2 model (see the header). Rows are target bases, columns query
// bases. diag[j] holds the H value cell (i, j) reads diagonally, H(i-1, j-1),
// as last stored; vgap[j] holds E(i, j), the deletion score carried down
// column j.
//
// With `drops` non-null, every row that does not raise the best score appends
// the z-drop statistic it computed (also when zdrop is 0 and nothing is cut).
ExtResult oracle(const std::vector<uint8_t> &q, const std::vector<uint8_t> &t, const Scoring &s, int w,
                 int h0, std::vector<int> *drops = nullptr) {
    const int Q = (int)q.size(), T = (int)t.size();
    int maxsc = 0;
    for (int k = 0; k < 25; ++k) maxsc = std::max(maxsc, (int)s.mat[k]);
    const int max_ins = std::max(1, (int)((double)(Q * maxsc + s.end_bonus - s.o_ins) / s.e_ins + 1.0));
    const int max_del = std::max(1, (int)((double)(Q * maxsc + s.end_bonus - s.o_del) / s.e_del + 1.0));
    w = std::min(w, std::min(max_ins, max_del));

    std::vector<int> diag(Q + 1, 0), vgap(Q + 1, 0);
    diag[0] = h0;
    for (int j = 1; j <= Q; ++j) diag[j] = std::max(0, h0 - s.o_ins - j * s.e_ins);

    int best = h0, best_i = -1, best_j = -1, gscore = -1, g_i = -1, max_off = 0;
    int beg = 0, end = Q;
    for (int i = 0; i < T; ++i) {
        beg = std::max(beg, i - w);
        end = std::min(end, std::min(i + w + 1, Q));
        int left = beg == 0 ? std::max(0, h0 - (s.o_del + s.e_del * (i + 1))) : 0;  // H(i, beg-1)
        int f = 0, rmax = 0, rarg = -1;
        for (int j = beg; j < end; ++j) {
            const int hd = diag[j];
            diag[j] = left;
            const int m = hd != 0 ? hd + s.mat[t[(size_t)i] * 5 + q[(size_t)j]] : 0;
            const int h = std::max(m, std::max(vgap[j], f));
            left = h;
            if (h >= rmax) rmax = h, rarg = j;
            vgap[j] = std::max(vgap[j] - s.e_del, std::max(0, m - s.o_del - s.e_del));
            f = std::max(f - s.e_ins, std::max(0, m - s.o_ins - s.e_ins));
        }
        const int jlast = std::max(beg, end);  // where the column loop stopped
        diag[end] = left;
        vgap[end] = 0;
        if (jlast == Q && left >= gscore) gscore = left, g_i = i;
        if (rmax == 0) break;
        if (rmax > best) {
            best = rmax, best_i = i, best_j = rarg;
            max_off = std::max(max_off, std::abs(rarg - i));
        } else if (s.zdrop > 0 || drops != nullptr) {
            const int di = i - best_i, dj = rarg - best_j;
            const int drop = di > dj ? best - rmax - (di - dj) * s.e_del : best - rmax - (dj - di) * s.e_ins;
            if (drops != nullptr) drops->push_back(drop);
            if (s.zdrop > 0 && drop > s.zdrop) break;
        }
        int nb = beg;
        while (nb < end && diag[nb] == 0 && vgap[nb] == 0) ++nb;
        int ne = end;
        while (ne >= nb && diag[ne] == 0 && vgap[ne] == 0) --ne;
        beg = nb;
        end = std::min(ne + 2, Q);
    }
    ExtResult r;
    r.score = best;
    r.qle = best_j + 1;
    r.tle = best_i + 1;
    r.gtle = g_i + 1;
    r.gscore = gscore;
    r.max_off = max_off;
    return r;
}

// The 8-bit routing envelope, as bsw8_envelope_ok documents it.
bool in_bsw8_envelope(const ExtPair &p, const Scoring &s, int w) {
    const int len1 = (int)p.target.size(), len2 = (int)p.query.size();
    const int max_step = std::max(s.a, 1);
    return len1 < MAX_SEQ_LEN8 && len2 < MAX_SEQ_LEN8 && len1 >= len2 && w <= 124 &&
           s.zdrop + max_step <= 253 && p.h0 + std::min(len1, len2) * s.a < 255 - max_step;
}

struct Kernel {
    const char *name;
    IBandedPairWiseSW *(*make)(int, int, int, int, int, int, const int8_t *, int8_t, int8_t, int);
};

#if defined(__aarch64__)
IBandedPairWiseSW *make_native(int o_del, int e_del, int o_ins, int e_ins, int zdrop, int end_bonus,
                               const int8_t *mat, int8_t a, int8_t b, int nthreads) {
    return new BandedPairWiseSW(o_del, e_del, o_ins, e_ins, zdrop, end_bonus, mat, a, b, nthreads);
}
#endif

std::vector<Kernel> runnable_kernels() {
    std::vector<Kernel> ks;
#if defined(__x86_64__) || defined(__i386__)
    __builtin_cpu_init();
    if (__builtin_cpu_supports("avx2")) ks.push_back({"avx2", make_bsw_kernel_avx2});
    if (__builtin_cpu_supports("avx512bw")) ks.push_back({"avx512bw", make_bsw_kernel_avx512bw});
    else MESSAGE("banded SW oracle: no AVX-512BW on this host; the avx512bw tier is skipped");
#elif defined(__aarch64__)
    ks.push_back({"neon", make_native});
#endif
    return ks;
}

std::unique_ptr<IBandedPairWiseSW> build(const Kernel &k, const Scoring &s) {
    return std::unique_ptr<IBandedPairWiseSW>(
        k.make(s.o_del, s.e_del, s.o_ins, s.e_ins, s.zdrop, s.end_bonus, s.mat, (int8_t)s.a, (int8_t)s.b, 1));
}

ExtResult scalar_result(IBandedPairWiseSW &bsw, const ExtPair &p, int w) {
    ExtResult r;
    r.score = bsw.scalarBandedSWA((int)p.query.size(), p.query.data(), (int)p.target.size(), p.target.data(), w,
                                  p.h0, &r.qle, &r.tle, &r.gtle, &r.gscore, &r.max_off);
    return r;
}

// Score `pairs` through every kernel width that may take them (8-bit only inside
// its envelope) in three batchings: as given, rotated by a third of the batch,
// and one pair per batch. Every result must equal want[k]. Returns the number of
// (pair, width, batching) results compared.
int check_batches(IBandedPairWiseSW &bsw, const char *tier, const std::vector<ExtPair> &pairs,
                  const std::vector<ExtResult> &want, const Scoring &s, int w) {
    int compared = 0;
    for (int width : {8, 16}) {
        std::vector<int> idx;
        for (int k = 0; k < (int)pairs.size(); ++k)
            if (width == 16 || in_bsw8_envelope(pairs[(size_t)k], s, w)) idx.push_back(k);
        if (idx.empty()) continue;
        const int n = (int)idx.size();
        for (int rot : {0, n / 3 + 1}) {
            std::vector<const ExtPair *> batch;
            std::vector<int> order;
            for (int r = 0; r < n; ++r) {
                const int k = idx[(size_t)((r + rot) % n)];
                batch.push_back(&pairs[(size_t)k]);
                order.push_back(k);
            }
            const std::vector<ExtResult> got = bwa_tests::score_bsw_batch(bsw, batch, width, w);
            for (int r = 0; r < n; ++r) {
                const int k = order[(size_t)r];
                CAPTURE(tier);
                CAPTURE(width);
                CAPTURE(w);
                CAPTURE(rot);
                CAPTURE(k);
                CAPTURE(pairs[(size_t)k].query.size());
                CAPTURE(pairs[(size_t)k].target.size());
                CAPTURE(pairs[(size_t)k].h0);
                INFO("got " << show(got[(size_t)r]) << "; want " << show(want[(size_t)k]));
                CHECK(got[(size_t)r] == want[(size_t)k]);
                ++compared;
            }
        }
        for (int k : idx) {
            const ExtResult got = bwa_tests::score_bsw_batch(bsw, {&pairs[(size_t)k]}, width, w)[0];
            CAPTURE(tier);
            CAPTURE(width);
            CAPTURE(w);
            CAPTURE(k);
            INFO("alone: got " << show(got) << "; want " << show(want[(size_t)k]));
            CHECK(got == want[(size_t)k]);
            ++compared;
        }
    }
    return compared;
}

std::vector<uint8_t> random_seq(std::mt19937 &rng, size_t n) {
    std::vector<uint8_t> v(n);
    for (size_t k = 0; k < n; ++k) v[k] = (uint8_t)(rng() % 4);
    return v;
}

// A query and a target that differ by the edits of one family.
enum Family {
    kMutated,       // substitutions and short indels
    kHomopolymer,   // runs of one base with a few others
    kAllMismatch,   // every aligned base differs
    kLongDeletion,  // the target carries 8-40 extra bases mid-way
    kLongInsertion, // the query carries 8-30 extra bases mid-way
    kWithN,         // N bases in both
    kUnrelated,     // independent sequences
    kNumFamilies
};

ExtPair make_pair(std::mt19937 &rng, Family fam, int max_len) {
    ExtPair p;
    const size_t len = 8 + rng() % (uint32_t)(max_len - 8);
    std::vector<uint8_t> q = random_seq(rng, len);
    std::vector<uint8_t> t;
    switch (fam) {
        case kMutated:
            for (size_t k = 0; k < q.size(); ++k) {
                const uint32_t r = rng() % 100;
                if (r < 6) t.push_back((uint8_t)((q[k] + 1 + rng() % 3) % 4));
                else if (r < 9) continue;
                else if (r < 12) {
                    for (uint32_t m = 1 + rng() % 3; m > 0; --m) t.push_back((uint8_t)(rng() % 4));
                    t.push_back(q[k]);
                } else t.push_back(q[k]);
            }
            break;
        case kHomopolymer: {
            const uint8_t base = (uint8_t)(rng() % 4);
            for (auto &c : q) c = rng() % 10 == 0 ? (uint8_t)((base + 1) % 4) : base;
            t = q;
            for (int m = (int)(rng() % 4); m > 0; --m) t.insert(t.begin() + (long)(rng() % t.size()), base);
            break;
        }
        case kAllMismatch:
            for (uint8_t c : q) t.push_back((uint8_t)((c + 1 + rng() % 3) % 4));
            break;
        case kLongDeletion: {
            t = q;
            const size_t at = 1 + rng() % (q.size() - 1);
            const std::vector<uint8_t> extra = random_seq(rng, 8 + rng() % 33);
            t.insert(t.begin() + (long)at, extra.begin(), extra.end());
            break;
        }
        case kLongInsertion: {
            t = q;
            const size_t at = 1 + rng() % (q.size() - 1);
            const std::vector<uint8_t> extra = random_seq(rng, 8 + rng() % 23);
            q.insert(q.begin() + (long)at, extra.begin(), extra.end());
            break;
        }
        case kWithN:
            for (auto &c : q)
                if (rng() % 15 == 0) c = 4;
            t = q;
            for (auto &c : t)
                if (rng() % 15 == 0) c = rng() % 3 == 0 ? (uint8_t)4 : (uint8_t)((c + 1) % 4);
            break;
        default:
            t = random_seq(rng, 8 + rng() % (uint32_t)(max_len - 8));
            break;
    }
    if (t.empty()) t.push_back(0);
    // A target tail beyond the alignment, as `mem` gives the extension a
    // reference window longer than the read segment.
    const std::vector<uint8_t> tail = random_seq(rng, rng() % 24);
    t.insert(t.end(), tail.begin(), tail.end());
    p.query = q;
    p.target = t;
    p.h0 = 1 + (int)(rng() % 60);
    return p;
}

// Pairs whose seed puts them at the top of the 8-bit kernel's score range:
// h0 + len2*a lands on 255 - max(a,1) - 1 - k for k in {0, 1, 2}.
ExtPair ceiling_pair(std::mt19937 &rng, const Scoring &s) {
    ExtPair p;
    const int max_step = std::max(s.a, 1);
    const int len2 = 20 + (int)(rng() % 60);
    const int len1 = len2 + (int)(rng() % 30);
    p.target = random_seq(rng, (size_t)len1);
    p.query.assign(p.target.begin(), p.target.begin() + len2);
    for (int m = (int)(rng() % 3); m > 0; --m) {
        const size_t at = rng() % (size_t)len2;
        p.query[at] = (uint8_t)((p.query[at] + 1) % 4);
    }
    p.h0 = 255 - max_step - 1 - (int)(rng() % 3) - len2 * s.a;
    return p;
}

// A short query matching the target only after `off` target bases, with a seed
// large enough to pay for the deletion: whether the band (clamped through
// max_del) reaches the matching diagonal decides gscore and gtle.
ExtPair clamp_pair(std::mt19937 &rng) {
    ExtPair p;
    const int qlen = 3 + (int)(rng() % 8);
    const int off = 2 + (int)(rng() % 20);
    p.target = random_seq(rng, (size_t)(off + qlen + rng() % 10));
    p.query.assign(p.target.begin() + off, p.target.begin() + off + qlen);
    p.h0 = 30 + (int)(rng() % 60);
    return p;
}

std::vector<uint8_t> digits(const char *s) {
    std::vector<uint8_t> v;
    for (; *s; ++s) v.push_back((uint8_t)(*s - '0'));
    return v;
}

ExtResult res(int score, int qle, int tle, int gscore, int gtle, int max_off) {
    ExtResult r;
    r.score = score;
    r.qle = qle;
    r.tle = tle;
    r.gscore = gscore;
    r.gtle = gtle;
    r.max_off = max_off;
    return r;
}

// One hand-worked case: the oracle, the scalar kernel and each vector width
// that may take the pair must return `want`.
void check_hand(const std::vector<Kernel> &kernels, const Scoring &s, const ExtPair &p, int w,
                const ExtResult &want) {
    INFO("want " << show(want));
    const ExtResult o = oracle(p.query, p.target, s, w, p.h0);
    INFO("oracle " << show(o));
    CHECK(o == want);
    for (const Kernel &k : kernels) {
        CAPTURE(k.name);
        std::unique_ptr<IBandedPairWiseSW> bsw = build(k, s);
        const ExtResult sc = scalar_result(*bsw, p, w);
        INFO("scalar " << show(sc));
        CHECK(sc == want);
        for (int width : {8, 16}) {
            if (width == 8 && !in_bsw8_envelope(p, s, w)) continue;
            CAPTURE(width);
            const ExtResult got = bwa_tests::score_bsw_batch(*bsw, {&p}, width, w)[0];
            INFO("kernel " << show(got));
            CHECK(got == want);
        }
    }
}

}  // namespace

TEST_CASE("banded SW oracle: hand-worked extensions on every runnable tier"
          * doctest::test_suite("unit/bandedswa")) {
    const std::vector<Kernel> kernels = runnable_kernels();
    REQUIRE(!kernels.empty());
    const Scoring def = scoring(1, 4, 6, 1, 6, 1, 100, 5);
    // 40 bases, no base equal to its neighbour, so a gap has one placement.
    const std::vector<uint8_t> q = digits("0123102302130120312013201230210310231032");
    REQUIRE(q.size() == 40);

    SUBCASE("exact match: h0 + 40, every field at the query end") {
        ExtPair p{q, q, 20};
        check_hand(kernels, def, p, 100, res(60, 40, 40, 60, 40, 0));
    }
    SUBCASE("a target longer than the query: gtle at the query's last base") {
        ExtPair p{q, q, 20};
        p.target.insert(p.target.end(), {0, 0, 1, 1, 2, 2, 3, 3});
        check_hand(kernels, def, p, 100, res(60, 40, 40, 60, 40, 0));
    }
    SUBCASE("a mismatch two bases from the end: the maximum stops before it, gscore runs through") {
        // 20 + 37 = 57 at (37,37); 57 - 4 = 53 at the mismatch; 53 + 2 = 55 at the end.
        ExtPair p{q, q, 20};
        p.target[37] = (uint8_t)((p.target[37] + 2) % 4);
        check_hand(kernels, def, p, 100, res(57, 37, 37, 55, 40, 0));
    }
    SUBCASE("a mismatch at base 10: the full extension wins") {
        // 20 + 39 - 4 = 55, above the 30 reached before the mismatch.
        ExtPair p{q, q, 20};
        p.target[10] = (uint8_t)((p.target[10] + 2) % 4);
        check_hand(kernels, def, p, 100, res(55, 40, 40, 55, 40, 0));
    }
    SUBCASE("three extra target bases: one deletion, o_del + 3*e_del, max_off 3") {
        // 20 + 40 - (6 + 3) = 51; the query ends on target base 43.
        ExtPair p{q, q, 20};
        p.target.insert(p.target.begin() + 20, {2, 2, 2});
        check_hand(kernels, def, p, 100, res(51, 40, 43, 51, 43, 3));
    }
    SUBCASE("two extra query bases: one insertion priced -O 5 -E 2") {
        // 20 + 40 - (5 + 2*2) = 51; the query ends 2 columns right of the diagonal.
        const Scoring ins = scoring(1, 4, 6, 1, 5, 2, 100, 5);
        ExtPair p{q, q, 20};
        p.query.insert(p.query.begin() + 20, {1, 1});
        check_hand(kernels, ins, p, 100, res(51, 42, 40, 51, 40, 2));
    }
    SUBCASE("a band narrower than the deletion keeps the extension on its first diagonal") {
        // With w = 2 the 3-base deletion is out of reach; the best cell is the
        // last one before it, 20 + 20 = 40 at (19,19).
        ExtPair p{q, q, 20};
        p.target.insert(p.target.begin() + 20, {2, 2, 2});
        const ExtResult o = oracle(p.query, p.target, def, 2, p.h0);
        CHECK(o.score == 40);
        CHECK(o.qle == 20);
        CHECK(o.tle == 20);
        check_hand(kernels, def, p, 2, o);
        check_hand(kernels, def, p, 3, res(51, 40, 43, 51, 43, 3));
    }
    SUBCASE("-A 2 scales the match and the seed") {
        const Scoring a2 = scoring(2, 8, 12, 2, 12, 2, 200, 10);
        ExtPair p{q, q, 40};
        p.target[30] = (uint8_t)((p.target[30] + 1) % 4);
        // 40 + 39*2 - 8 = 110 at the end; 40 + 30*2 = 100 before the mismatch.
        check_hand(kernels, a2, p, 100, res(110, 40, 40, 110, 40, 0));
    }
    SUBCASE("a seed too small for the first mismatch ends the extension at once") {
        // h0 = 3: the first base mismatches (3 - 4 < 0) and no gap is affordable,
        // so row 0's maximum is the seed-row deletion/insertion score 0: stop.
        ExtPair p{q, q, 3};
        p.target[0] = (uint8_t)((p.target[0] + 1) % 4);
        p.query[0] = (uint8_t)((p.query[0] + 2) % 4);
        // Row 0 still reaches the last query column (with H = 0), so gscore is 0
        // at gtle 1.
        check_hand(kernels, def, p, 100, res(3, 0, 0, 0, 1, 0));
    }
}

TEST_CASE("banded SW oracle: z-drop stops at the documented threshold on every runnable tier"
          * doctest::test_suite("unit/bandedswa")) {
    // 30 matching bases, two mismatches at query/target 30 and 31, then 30
    // matching bases again. Along the diagonal the row maximum is 20 + 30 = 50
    // at row 29, then 46 and 42: a drop of 8 at row 31, after which the
    // alignment recovers to 50 - 8 + 30 = 72. A z-drop of 8 lets it through
    // (the test is "drop > zdrop"); 7 stops it at row 31 with the best cell
    // (29,29). The pair is chosen (seed 20261024) so that no off-diagonal cell
    // of rows 30 and 31 beats the diagonal; the oracle confirms both results.
    const std::vector<Kernel> kernels = runnable_kernels();
    REQUIRE(!kernels.empty());
    std::mt19937 rng(20261024u);
    ExtPair p;
    for (int attempt = 0;; ++attempt) {
        REQUIRE(attempt < 1000);
        p.query = random_seq(rng, 62);
        p.target = p.query;
        p.target[30] = (uint8_t)((p.target[30] + 1) % 4);
        p.target[31] = (uint8_t)((p.target[31] + 1) % 4);
        p.h0 = 20;
        const Scoring z8 = scoring(1, 4, 6, 1, 6, 1, 8, 5);
        const Scoring z7 = scoring(1, 4, 6, 1, 6, 1, 7, 5);
        if (oracle(p.query, p.target, z8, 100, 20) == res(72, 62, 62, 72, 62, 0) &&
            oracle(p.query, p.target, z7, 100, 20).score == 50)
            break;
    }
    // When stopped, gscore comes from the long insertion to the query end in
    // rows 29-31; its value depends on the bases, so the oracle supplies it.
    const Scoring z8 = scoring(1, 4, 6, 1, 6, 1, 8, 5);
    const Scoring z7 = scoring(1, 4, 6, 1, 6, 1, 7, 5);
    const ExtResult stopped = oracle(p.query, p.target, z7, 100, 20);
    CHECK(stopped.score == 50);
    CHECK(stopped.qle == 30);
    CHECK(stopped.tle == 30);
    check_hand(kernels, z8, p, 100, res(72, 62, 62, 72, 62, 0));
    check_hand(kernels, z7, p, 100, stopped);
    // zdrop = 0 disables the test: the alignment runs to the end.
    check_hand(kernels, scoring(1, 4, 6, 1, 6, 1, 0, 5), p, 100, res(72, 62, 62, 72, 62, 0));
}

namespace {

// For each scoring, 50 pairs (every family six times, four seeds at the 8-bit
// ceiling, four short queries at the band clamp), scored at bands 1, 33, 124
// and 200 (the wide-band retry; 16-bit only) through every tier's kernels and
// the scalar kernel, against the oracle.
void check_random(const Scoring *scorings, int n_scorings, uint32_t seed) {
    const std::vector<Kernel> kernels = runnable_kernels();
    REQUIRE(!kernels.empty());
    const int widths[] = {1, 33, 124, 200};
    std::mt19937 rng(seed);
    int compared = 0, compared8 = 0;
    for (int si = 0; si < n_scorings; ++si) {
        const Scoring &s = scorings[si];
        std::vector<ExtPair> pairs;
        for (int n = 0; n < 6 * kNumFamilies; ++n) pairs.push_back(make_pair(rng, (Family)(n % kNumFamilies), 90));
        for (int n = 0; n < 4; ++n) pairs.push_back(ceiling_pair(rng, s));
        for (int n = 0; n < 4; ++n) pairs.push_back(clamp_pair(rng));
        std::vector<std::unique_ptr<IBandedPairWiseSW> > bsw;
        for (const Kernel &k : kernels) bsw.push_back(build(k, s));
        for (int w : widths) {
            std::vector<ExtResult> want;
            for (const ExtPair &p : pairs) want.push_back(oracle(p.query, p.target, s, w, p.h0));
            for (const ExtPair &p : pairs)
                if (in_bsw8_envelope(p, s, w)) ++compared8;
            for (size_t ki = 0; ki < kernels.size(); ++ki) {
                compared += check_batches(*bsw[ki], kernels[ki].name, pairs, want, s, w);
                for (size_t i = 0; i < pairs.size(); ++i) {
                    CAPTURE(kernels[ki].name);
                    CAPTURE(w);
                    CAPTURE(i);
                    CHECK(scalar_result(*bsw[ki], pairs[i], w) == want[i]);
                }
            }
        }
    }
    // Both kernel widths saw many pairs, so the comparison cannot pass vacuously.
    CHECK(compared8 > 100 * n_scorings);
    CHECK(compared > 500 * n_scorings);
}

}  // namespace

TEST_CASE("banded SW oracle: random and adversarial pairs equal the ksw_extend2 model in every lane"
          * doctest::test_suite("unit/bandedswa")) {
    const Scoring scorings[] = {
        scoring(1, 4, 6, 1, 6, 1, 100, 5),  // bwa-mem3 defaults
        scoring(1, 4, 5, 2, 7, 1, 30, 5),   // -O 5,7 -E 2,1 -d 30: gaps priced apart
    };
    check_random(scorings, 2, 424242u);
}

TEST_CASE("banded SW oracle: random pairs under -A 2 and the --meth OT matrix equal the model in every lane"
          * doctest::test_suite("unit/bandedswa")) {
    const Scoring scorings[] = {
        scoring(2, 8, 12, 2, 12, 2, 200, 10),     // -A 2, which scales -B -O -E -d -L with it
        scoring(1, 2, 6, 1, 6, 1, 100, 5, true),  // --meth OT matrix (-B 2), asymmetric
    };
    check_random(scorings, 2, 515151u);
}

TEST_CASE("banded SW oracle: large batches through lane compaction equal the model"
          * doctest::test_suite("unit/bandedswa")) {
    // 150 pairs: more than two 64-lane groups, so the shipped lane compaction
    // (bands of 31 or more, batches of more than one group) moves lanes between
    // groups as their pairs finish at different rows.
    const std::vector<Kernel> kernels = runnable_kernels();
    REQUIRE(!kernels.empty());
    const Scoring s = scoring(1, 4, 6, 1, 6, 1, 100, 5);
    std::mt19937 rng(77u);
    std::vector<ExtPair> pairs;
    for (int n = 0; n < 150; ++n) {
        ExtPair p = make_pair(rng, (Family)(n % kNumFamilies), 60);
        if (p.target.size() < p.query.size()) p.target.resize(p.query.size(), 3);
        if (p.h0 + (int)p.query.size() >= 253) p.h0 = 1 + n % 20;
        pairs.push_back(p);
    }
    const int w = 100;
    std::vector<ExtResult> want;
    for (const ExtPair &p : pairs) want.push_back(oracle(p.query, p.target, s, w, p.h0));
    for (const Kernel &k : kernels) {
        CAPTURE(k.name);
        std::unique_ptr<IBandedPairWiseSW> bsw = build(k, s);
        for (int width : {8, 16}) {
            CAPTURE(width);
            std::vector<const ExtPair *> batch;
            for (const ExtPair &p : pairs) batch.push_back(&p);
            const std::vector<ExtResult> got = bwa_tests::score_bsw_batch(*bsw, batch, width, w);
            int bad = 0;
            for (size_t i = 0; i < pairs.size(); ++i) {
                REQUIRE(in_bsw8_envelope(pairs[i], s, w));
                if (!(got[i] == want[i]) && bad++ < 3)
                    MESSAGE(k.name << " getScores" << width << " pair " << i << ": " << show(got[i]) << "; want "
                                   << show(want[i]));
            }
            CHECK(bad == 0);
        }
    }
}

TEST_CASE("banded SW oracle: results at the exact z-drop and band width where they change"
          * doctest::test_suite("unit/bandedswa")) {
    // The rows of an extension do not depend on the z-drop until it stops one,
    // so one run with the statistic recorded gives every row's drop: a z-drop of
    // Z stops at the first row whose drop exceeds Z. For each pair the largest
    // drop D is the edge: at zdrop = D - 1 the extension stops early, at D it
    // does not. The band edge is the result's own max_off M: at w = M - 1 the
    // best cell found at w = 124 lies outside the band. The kernels are checked
    // on both sides of each edge, so a comparison off by one (">=" for ">", a
    // band one cell narrower or wider) fails here.
    const std::vector<Kernel> kernels = runnable_kernels();
    REQUIRE(!kernels.empty());
    std::mt19937 rng(9001u);
    int zedges = 0, wedges = 0;
    for (int n = 0; n < 24; ++n) {
        const Family fam = n % 3 == 0 ? kMutated : n % 3 == 1 ? kLongDeletion : kLongInsertion;
        ExtPair p = make_pair(rng, fam, 70);
        if (p.target.size() < p.query.size()) p.target.resize(p.query.size(), 2);
        p.h0 = 10 + (int)(rng() % 40);
        CAPTURE(n);
        const Scoring off = scoring(1, 4, 6, 1, 6, 1, 0, 5);
        std::vector<int> drops;
        const ExtResult free_run = oracle(p.query, p.target, off, 100, p.h0, &drops);
        const int dmax = drops.empty() ? 0 : *std::max_element(drops.begin(), drops.end());
        if (dmax >= 2 && dmax <= 250) {
            const Scoring below = scoring(1, 4, 6, 1, 6, 1, dmax - 1, 5);
            const Scoring at = scoring(1, 4, 6, 1, 6, 1, dmax, 5);
            const ExtResult rb = oracle(p.query, p.target, below, 100, p.h0);
            const ExtResult ra = oracle(p.query, p.target, at, 100, p.h0);
            CHECK(ra == free_run);
            if (!(rb == ra)) ++zedges;
            check_hand(kernels, below, p, 100, rb);
            check_hand(kernels, at, p, 100, ra);
        }
        const ExtResult wide = oracle(p.query, p.target, off, 124, p.h0);
        if (wide.max_off >= 1) {
            const ExtResult narrow = oracle(p.query, p.target, off, wide.max_off - 1, p.h0);
            if (!(narrow == wide)) ++wedges;
            for (int w = wide.max_off - 1; w <= wide.max_off + 1; ++w) {
                if (w < 1) continue;
                CAPTURE(w);
                check_hand(kernels, off, p, w, oracle(p.query, p.target, off, w, p.h0));
            }
        }
    }
    // The edges were really crossed: on enough pairs the result changes there.
    CHECK(zedges >= 6);
    CHECK(wedges >= 6);
}

#endif  // HAVE_BSW_VECTOR_8_16
