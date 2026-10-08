// test/unit/test_bandedswa_row_lean.cpp
//
// Exactness gate for BSW8_ROW_LEAN, the lean per-row bookkeeping in the 128-bit
// 8-bit banded extension kernel (smithWaterman128_8, bandedSWA.cpp): counted
// band-trim loop, one-pass band narrowing, and the narrow per-row epilogue with
// the z-drop test behind a gate. It is a pure speed change, so getScores8 must
// return every result field (score, tle, gtle, qle, gscore, max_off) unchanged
// for every pair of every batch.
//
// On aarch64 the unit binary links two copies of bandedSWA.cpp: libbwa.a's, and
// src/bandedSWA.rowalt.o, built with BSW8_ROW_LEAN set the other way
// (BSW8_ROW_LEAN_INVERT) and its symbols renamed with a _rowalt suffix
// (KERNEL_VARIANT). By default the libbwa.a kernel is the lean one on Linux arm64
// and the original on macOS; the first case checks that the two really differ.
// The original code is the reference: it is the kernel that shipped before
// BSW8_ROW_LEAN, an independent implementation of every per-row step. x86 always
// runs the original code, so there the file holds a single case that says so.
//
// Each batch is scored in several groupings (production's length sort, a random
// shuffle, reversed, and a random split into sub-batches), and for each grouping
// the lean kernel must equal the original kernel on that same grouping. A change
// to how one lane's result depends on its neighbours would show up as a
// difference here. Pairs scored alone are also checked against the scalar
// oracle. That a pair's result does not depend on its group at all (grouped ==
// alone == scalar) is test_bandedswa_lane_indep.cpp's property; this file
// compares the two kernel copies.
//
// Generated batches (fixed seeds): default and random scoring (match, mismatch,
// gap open / extend, end bonus), --meth OT / OB matrices under all three
// --meth-scoring modes with negative seed scores, N bases, bands 1 to 124,
// default, small, disabled and above-255 z-drop, gap extends large enough to
// wrap the z-drop drift, short and long targets, and query lengths up to the
// 8-bit routing envelope.

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "doctest/doctest.h"
#include "bandedSWA.h"
#include "bwamem.h"
#include "meth_scoring.h"
#include "bsw_batch.h"
#include "ext_ladder.h"
#include "scoring.h"

#if defined(__aarch64__)
#if !HAVE_BSW_VECTOR_8_16
#error "aarch64 builds have the 8-bit vector kernels; this test would otherwise compile to nothing"
#endif

/* Defined by src/bandedSWA.rowalt.o (make_bsw_kernel_rowalt is declared in
 * bsw_batch.h). */
extern "C" int bsw8_row_lean_enabled_rowalt(void);
/* Defined by libbwa.a's bandedSWA.o: the BSW8_ROW_LEAN setting it was built with. */
extern "C" int bsw8_row_lean_enabled(void);

namespace {

using bwa_tests::ExtPair;     // query, target, h0
using bwa_tests::ExtResult;   // the six result fields, operator==
using bwa_tests::ExtScoring;  // a, b, gap open / extend, pen_clip (end bonus), zdrop, mat

/// bwa-mem3's default extension scoring (-A1 -B4 -O6 -E1 -L5 -d100).
ExtScoring default_scoring()
{
    ExtScoring sc;
    sc.a = 1; sc.b = 4; sc.o_del = 6; sc.e_del = 1; sc.o_ins = 6; sc.e_ins = 1;
    sc.pen_clip = 5; sc.zdrop = 100;
    sc.fill_mat();
    return sc;
}

class Rng {
public:
    explicit Rng(uint64_t seed) : g_(seed) {}
    int below(int n) { return n <= 0 ? 0 : (int)(g_() % (uint64_t)n); }
    int range(int lo, int hi) { return lo + below(hi - lo + 1); }
    double unif() { return (g_() >> 11) * (1.0 / 9007199254740992.0); }
    std::mt19937_64 &engine() { return g_; }

private:
    std::mt19937_64 g_;
};

/// One batch of the given pairs through getScores8, results in input order.
std::vector<ExtResult> score8(IBandedPairWiseSW &bsw, const std::vector<const ExtPair *> &pairs, int w)
{
    return bwa_tests::score_bsw_batch(bsw, pairs, 8, w);
}

/// A query derived from `src` with substitutions, indels and N bases.
std::vector<uint8_t> mutate(Rng &r, const std::vector<uint8_t> &src, int len, double sub,
                            double indel, double npct)
{
    std::vector<uint8_t> dst;
    size_t k = 0;
    while ((int)dst.size() < len) {
        uint8_t c = k < src.size() ? src[k] : (uint8_t)r.below(4);
        const double u = r.unif();
        if (u < indel / 2) { k++; continue; }                            // deletion from the query
        if (u < indel) { dst.push_back((uint8_t)r.below(4)); continue; }  // insertion
        if (r.unif() < sub) c = (uint8_t)((c + 1 + r.below(3)) & 3);
        if (r.unif() < npct) c = 4;
        dst.push_back(c);
        k++;
    }
    return dst;
}

/// Shapes of generated pairs.
enum class Shape {
    Mixed,        // near-identical, divergent and random pairs, a few long targets
    ShortQuery,   // short queries on long targets: rows run far past the query end
    LongQuery,    // queries near the 8-bit envelope through deletion-rich targets: rows > 255
};

/// Pairs inside the 8-bit routing envelope (h0 + min(len1, len2) * max step
/// below 255 - max step).
std::vector<ExtPair> make_pairs(Rng &r, const ExtScoring &sc, int w, bool meth, double npct,
                                Shape shape)
{
    const int n = r.unif() < 0.3 ? r.range(1, 20) : r.range(16, 64);
    const int maxstep = sc.a > 1 ? sc.a : 1;
    const int lenmax = r.unif() < 0.15 ? 400 : 150;
    std::vector<ExtPair> pairs(n);
    for (ExtPair &p : pairs) {
        if (shape == Shape::LongQuery) {
            // A query near the envelope (match score 1, small seed) and a target that
            // is the query with 10-24 single-base insertions (deletions from the
            // query), so the best and query-end rows pass 255 and the int16 row side
            // channels carry values no byte could.
            const int len2 = r.range(232, 252 - maxstep - 4);
            p.h0 = r.range(1, 4);
            p.query.resize(len2);
            for (auto &c : p.query) c = (uint8_t)r.below(4);
            const int ndel = r.range(10, 24);
            for (int k = 0; k < len2; k++) {
                p.target.push_back(r.unif() < 0.01 ? (uint8_t)((p.query[k] + 1) & 3) : p.query[k]);
                if (r.below(len2) < ndel) p.target.push_back((uint8_t)r.below(4));
            }
            continue;
        }
        int len2 = r.range(1, lenmax);
        int len1 = len2 + r.below(w + 30);
        const double u = r.unif();
        if (u < 0.03) {
            len1 = r.range(len2, MAX_SEQ_LEN8 - 1);   // long target: rows past the int8 range
        } else if (shape == Shape::ShortQuery ? u < 0.97 : u < 0.10) {
            // Short query against a long target: the band keeps trimming while rows run
            // far past the query end, the case a band-trim off-by-one shows up in.
            len2 = r.range(1, 24);
            len1 = r.unif() < 0.7 ? len2 + r.below(w + 120) : r.range(len2, MAX_SEQ_LEN8 - 1);
        }
        if (len1 >= MAX_SEQ_LEN8) len1 = MAX_SEQ_LEN8 - 1;
        int hmax = 255 - maxstep - 1 - std::min(len1, len2) * sc.a;
        while (hmax < 0) {
            len2 = std::max(1, len2 / 2);
            len1 = std::min(len1, len2 + w + 30);
            hmax = 255 - maxstep - 1 - std::min(len1, len2) * sc.a;
        }
        const double hu = r.unif();
        p.h0 = hu < 0.2 ? r.below(std::min(hmax, 12) + 1)
             : hu < 0.4 ? hmax - r.below(std::min(hmax, 12) + 1)   // high seeds
                        : r.below(hmax + 1);
        if (meth && r.unif() < 0.1) p.h0 = -r.range(0, 6);   // --meth rescored seeds
        p.target.resize(len1);
        for (auto &c : p.target) c = r.unif() < npct ? 4 : (uint8_t)r.below(4);
        const double kind = r.unif();
        if (kind < 0.15) {
            p.query.resize(len2);
            for (auto &c : p.query) c = r.unif() < npct ? 4 : (uint8_t)r.below(4);
        } else {
            p.query = mutate(r, p.target, len2, kind < 0.6 ? 0.02 : 0.12, kind < 0.8 ? 0.01 : 0.05, npct);
        }
        if (meth)   // bisulfite-like C->T in the read, so the freed cell is hit
            for (auto &c : p.query)
                if (c == 1 && r.unif() < 0.5) c = 3;
    }
    return pairs;
}

std::unique_ptr<IBandedPairWiseSW> make_kernel(const ExtScoring &sc, bool rowalt)
{
    if (rowalt)
        return std::unique_ptr<IBandedPairWiseSW>(make_bsw_kernel_rowalt(
            sc.o_del, sc.e_del, sc.o_ins, sc.e_ins, sc.zdrop, sc.pen_clip, sc.mat,
            (int8_t)sc.a, (int8_t)sc.b, 1));
    return std::unique_ptr<IBandedPairWiseSW>(new BandedPairWiseSW(
        sc.o_del, sc.e_del, sc.o_ins, sc.e_ins, sc.zdrop, sc.pen_clip, sc.mat,
        (int8_t)sc.a, (int8_t)sc.b, 1));
}

/// The lean and the original kernel for one scoring (whichever copy is which).
struct KernelPair {
    std::unique_ptr<IBandedPairWiseSW> lean, orig;
    explicit KernelPair(const ExtScoring &sc)
    {
        const bool lib_is_lean = bsw8_row_lean_enabled() != 0;
        lean = make_kernel(sc, !lib_is_lean);
        orig = make_kernel(sc, lib_is_lean);
    }
};

std::string show(const ExtResult &x)
{
    return std::to_string(x.score) + "/" + std::to_string(x.tle) + "/" + std::to_string(x.gtle) +
           "/" + std::to_string(x.qle) + "/" + std::to_string(x.gscore) + "/" +
           std::to_string(x.max_off);
}

struct Totals {
    long batches = 0, groupings = 0, pairs = 0, diffs = 0, zdrop_effect = 0, gscore_set = 0;
    long tle_over_255 = 0, gtle_over_255 = 0;
};

/// Score `order` (one grouping of a batch) with both kernels; count differences.
void compare_grouping(KernelPair &k, const std::vector<const ExtPair *> &order, int w, Totals &t)
{
    const std::vector<ExtResult> a = score8(*k.lean, order, w), b = score8(*k.orig, order, w);
    for (size_t q = 0; q < order.size(); q++) {
        if (!(a[q] == b[q])) {
            if (t.diffs < 5)
                MESSAGE("len1=" << order[q]->target.size() << " len2=" << order[q]->query.size()
                        << " h0=" << order[q]->h0 << " w=" << w << " | lean " << show(a[q])
                        << " | original " << show(b[q]) << " (score/tle/gtle/qle/gscore/max_off)");
            t.diffs++;
        }
        if (a[q].gscore >= 0) t.gscore_set++;
        if (a[q].tle > 255) t.tle_over_255++;
        if (a[q].gtle > 255) t.gtle_over_255++;
    }
    t.pairs += (long)order.size();
    t.groupings++;
}

enum class Mode { Default, RandomScoring, Meth, NBases, SmallZdrop, Combined, ShortQuery,
                  TinyBand, LongQuery, HugeExtend };

/// `nbatch` generated batches of one mode, each scored under four groupings:
/// production's length sort, a random shuffle, reversed, and a random split
/// into sub-batches.
Totals run_mode(Mode mode, uint64_t seed, int nbatch)
{
    Rng r(seed);
    Totals t;
    static const int kBands[] = {20, 30, 50, 100};
    const bool short_q = mode == Mode::ShortQuery, tiny = mode == Mode::TinyBand;
    const bool long_q = mode == Mode::LongQuery, huge_e = mode == Mode::HugeExtend;
    const bool small_zd = mode == Mode::SmallZdrop || mode == Mode::Combined || huge_e;
    for (int bi = 0; bi < nbatch; bi++) {
        ExtScoring sc = default_scoring();
        const bool rscore = mode == Mode::RandomScoring || mode == Mode::Combined || short_q ||
                            tiny || long_q || huge_e;
        if (rscore) {
            sc.a = long_q ? 1 : r.range(1, 3);
            sc.b = r.range(1, 6);
            // Cheap gap opens are over-weighted: long-lived F (gap) cells are what a
            // band-trim off-by-one leaks into the next row.
            sc.o_del = short_q || tiny || long_q || r.unif() < 0.4 ? r.range(0, 2) : r.range(0, 10);
            sc.e_del = long_q ? 1 : r.range(1, 3);
            sc.o_ins = short_q || tiny || r.unif() < 0.4 ? r.range(0, 2) : r.range(0, 10);
            sc.e_ins = r.range(1, 3);
            sc.pen_clip = r.range(0, 10);
            const double u = r.unif();
            sc.zdrop = u < 0.1 ? 0 : u < 0.5 ? r.range(1, 40) : u < 0.9 ? r.range(1, 253 - sc.a)
                                                                  : r.range(254, 1000);
            sc.fill_mat();
        }
        if (huge_e) {
            // Gap extends so large that the z-drop step's int32 |drift| * e wraps
            // (INT32_MAX, 2^30, 2^32 / 3 rounded up, 2^30 + 1, 2^29, 2^30 - 1), negative
            // ones (the API does not reject them), mixed with small ones.
            static const int kE[] = {2147483647, 1 << 30, 0x55555556, 0x40000001, 1 << 29,
                                     0x3FFFFFFF, -1, -3};
            sc.e_del = r.unif() < 0.7 ? kE[r.below(8)] : r.range(1, 5);
            sc.e_ins = r.unif() < 0.7 ? kE[r.below(8)] : r.range(1, 5);
            sc.a = r.range(1, 5);
            sc.fill_mat();
        }
        // A wrapped dif is a small negative number, so the original kills a lane at
        // a drop just below zdrop; keep zdrop small there so such drops are common.
        if (small_zd) sc.zdrop = huge_e ? r.range(1, 4) : r.range(1, 30);
        bool meth = false;
        if ((mode == Mode::Meth || mode == Mode::Combined) && r.unif() < 0.8) {
            static const int kMeth[] = {MEM_METH_SCORING_GENOMIC, MEM_METH_SCORING_NEUTRAL,
                                        MEM_METH_SCORING_COLLAPSED};
            memcpy(sc.mat,
                   bwa_tests::meth_scoring_matrix(kMeth[r.below(3)], r.below(2) == 1, sc.a, sc.b).data(),
                   sizeof sc.mat);
            meth = true;
        }
        const double npct = mode == Mode::NBases ? 0.10 : mode == Mode::Combined ? 0.05 : 0.0;
        const int w = tiny ? r.range(1, 3) : rscore ? r.range(1, 124) : kBands[r.below(4)];
        const Shape shape = short_q ? Shape::ShortQuery : long_q ? Shape::LongQuery : Shape::Mixed;

        const std::vector<ExtPair> pairs = make_pairs(r, sc, w, meth, npct, shape);
        KernelPair k(sc);
        std::vector<const ExtPair *> order;
        for (const ExtPair &p : pairs) order.push_back(&p);
        // 1. production order: stable sort by max(len1, len2), as bsw_score_batch does
        std::stable_sort(order.begin(), order.end(), [](const ExtPair *x, const ExtPair *y) {
            return std::max(x->target.size(), x->query.size()) <
                   std::max(y->target.size(), y->query.size());
        });
        compare_grouping(k, order, w, t);
        // 2. a random shuffle, 3. reversed
        std::shuffle(order.begin(), order.end(), r.engine());
        compare_grouping(k, order, w, t);
        std::reverse(order.begin(), order.end());
        compare_grouping(k, order, w, t);
        // 4. a random split into sub-batches (ragged groups, different neighbours)
        for (size_t s = 0; s < order.size();) {
            const size_t len = std::min(order.size() - s, (size_t)r.range(1, 40));
            compare_grouping(k, std::vector<const ExtPair *>(order.begin() + s, order.begin() + s + len),
                             w, t);
            s += len;
        }
        // The z-drop modes must reach the lean kernel's z-drop path: count pairs
        // whose lean result changes when z-drop is disabled (zdrop = 0).
        if (small_zd) {
            ExtScoring off = sc;
            off.zdrop = 0;
            KernelPair koff(off);
            const std::vector<ExtResult> zd = score8(*k.lean, order, w), nz = score8(*koff.lean, order, w);
            for (size_t q = 0; q < order.size(); q++)
                if (!(zd[q] == nz[q])) t.zdrop_effect++;
        }
        t.batches++;
    }
    return t;
}

void check_mode(Mode mode, uint64_t seed, int nbatch)
{
    const Totals t = run_mode(mode, seed, nbatch);
    MESSAGE("batches=" << t.batches << " groupings=" << t.groupings << " pair-scorings=" << t.pairs
            << " gscore-set=" << t.gscore_set << " tle>255=" << t.tle_over_255 << " gtle>255="
            << t.gtle_over_255 << " zdrop-effect=" << t.zdrop_effect << " diffs=" << t.diffs);
    CHECK(t.batches == nbatch);
    CHECK(t.groupings >= 4 * t.batches);   // every batch scored under all four groupings
    CHECK(t.gscore_set > 0);               // query-end capture exercised
    if (mode == Mode::LongQuery) {         // rows beyond a byte reach the outputs
        CHECK(t.tle_over_255 > 0);
        CHECK(t.gtle_over_255 > 0);
    }
    if (mode == Mode::SmallZdrop || mode == Mode::Combined || mode == Mode::HugeExtend)
        CHECK(t.zdrop_effect > 0);
    CHECK(t.diffs == 0);
}

} // namespace

TEST_CASE("BSW8_ROW_LEAN: the two linked kernel copies differ in BSW8_ROW_LEAN"
          * doctest::test_suite("unit/bandedswa")) {
    // Otherwise the cases below would compare a kernel with itself.
    CHECK(bsw8_row_lean_enabled() != bsw8_row_lean_enabled_rowalt());
}

TEST_CASE("BSW8_ROW_LEAN: both kernels equal the scalar oracle for pairs scored alone"
          * doctest::test_suite("unit/bandedswa")) {
    // Independent oracle: a pair alone in its group (the other 15 lanes are
    // padding), with both the lean and the original kernel, must equal
    // scalarBandedSWA on all six fields.
    Rng r(0x5EB8A17ull);
    long checked = 0, diffs = 0, zdrop_low = 0;
    for (int bi = 0; bi < 24; bi++) {
        ExtScoring sc = default_scoring();
        sc.a = r.range(1, 3);
        sc.b = r.range(1, 6);
        sc.o_del = r.range(0, 10);
        sc.e_del = r.range(1, 3);
        sc.o_ins = r.range(0, 10);
        sc.e_ins = r.range(1, 3);
        sc.pen_clip = r.range(0, 10);
        sc.zdrop = r.unif() < 0.5 ? r.range(1, 30) : 100;
        if (sc.zdrop <= 30) zdrop_low++;
        sc.fill_mat();
        const int w = r.range(1, 124);
        const std::vector<ExtPair> pairs = make_pairs(r, sc, w, false, 0.02, Shape::Mixed);
        KernelPair k(sc);
        for (size_t q = 0; q < pairs.size() && q < 8; q++) {
            const ExtPair &p = pairs[q];
            const ExtResult o = bwa_tests::run_mem3_scalar(sc, p, w);
            const std::vector<const ExtPair *> one{&p};
            const ExtResult a = score8(*k.lean, one, w)[0], b = score8(*k.orig, one, w)[0];
            if (!(a == o) || !(b == o)) {
                if (diffs < 5)
                    MESSAGE("len1=" << p.target.size() << " len2=" << p.query.size() << " h0=" << p.h0
                            << " w=" << w << " | lean " << show(a) << " | original " << show(b)
                            << " | scalar " << show(o) << " (score/tle/gtle/qle/gscore/max_off)");
                diffs++;
            }
            checked++;
        }
    }
    CHECK(checked > 100);
    CHECK(zdrop_low > 0);
    CHECK(diffs == 0);
}

TEST_CASE("BSW8_ROW_LEAN: lean == original on every field, default scoring"
          * doctest::test_suite("unit/bandedswa")) {
    check_mode(Mode::Default, 0x5EB8A11ull, 40);
}
TEST_CASE("BSW8_ROW_LEAN: lean == original on every field, random scoring and bands"
          * doctest::test_suite("unit/bandedswa")) {
    check_mode(Mode::RandomScoring, 0x5EB8A12ull, 40);
}
TEST_CASE("BSW8_ROW_LEAN: lean == original on every field, --meth matrices"
          * doctest::test_suite("unit/bandedswa")) {
    check_mode(Mode::Meth, 0x5EB8A13ull, 40);
}
TEST_CASE("BSW8_ROW_LEAN: lean == original on every field, 10% N bases"
          * doctest::test_suite("unit/bandedswa")) {
    check_mode(Mode::NBases, 0x5EB8A14ull, 40);
}
TEST_CASE("BSW8_ROW_LEAN: lean == original on every field, small z-drop"
          * doctest::test_suite("unit/bandedswa")) {
    check_mode(Mode::SmallZdrop, 0x5EB8A15ull, 28);
}
TEST_CASE("BSW8_ROW_LEAN: lean == original on every field, combined"
          * doctest::test_suite("unit/bandedswa")) {
    check_mode(Mode::Combined, 0x5EB8A16ull, 40);
}
TEST_CASE("BSW8_ROW_LEAN: lean == original on every field, short queries on long targets"
          * doctest::test_suite("unit/bandedswa")) {
    // Cheap gap opens, high seeds and rows running far past the query end: the
    // band keeps trimming on every row, which is where a band-trim off-by-one shows.
    check_mode(Mode::ShortQuery, 0x5EB8A18ull, 150);
}
TEST_CASE("BSW8_ROW_LEAN: lean == original on every field, bands of 1 to 3"
          * doctest::test_suite("unit/bandedswa")) {
    check_mode(Mode::TinyBand, 0x5EB8A1Aull, 150);
}
TEST_CASE("BSW8_ROW_LEAN: lean == original on every field, best and query-end rows past 255"
          * doctest::test_suite("unit/bandedswa")) {
    check_mode(Mode::LongQuery, 0x5EB8A1Bull, 12);
}
/// A directed pair: scoring, band, seed and sequences (bases as digits 0-4).
struct DirectedCase {
    int a, b, o_del, o_ins, e_del, e_ins, pen_clip, zdrop, w, h0;
    const char *ref, *qry;
};

/// Score `c` alone and as a full 16-lane group of copies with both kernels;
/// every field must agree.
void check_directed(const DirectedCase &c)
{
    ExtScoring sc = default_scoring();
    sc.a = c.a; sc.b = c.b; sc.o_del = c.o_del; sc.o_ins = c.o_ins;
    sc.e_del = c.e_del; sc.e_ins = c.e_ins; sc.pen_clip = c.pen_clip; sc.zdrop = c.zdrop;
    sc.fill_mat();
    ExtPair p;
    p.h0 = c.h0;
    for (const char *q = c.ref; *q; q++) p.target.push_back((uint8_t)(*q - '0'));
    for (const char *q = c.qry; *q; q++) p.query.push_back((uint8_t)(*q - '0'));
    KernelPair k(sc);
    const std::vector<const ExtPair *> one{&p}, group(16, &p);
    for (const auto *batch : {&one, &group}) {
        const std::vector<ExtResult> x = score8(*k.lean, *batch, c.w), y = score8(*k.orig, *batch, c.w);
        for (size_t q = 0; q < batch->size(); q++) {
            INFO("lean " << show(x[q]) << " | original " << show(y[q]) << " (score/tle/gtle/qle/gscore/max_off)");
            CHECK(x[q] == y[q]);
        }
    }
}

TEST_CASE("BSW8_ROW_LEAN: lean == original on directed band-edge pairs"
          * doctest::test_suite("unit/bandedswa")) {
    // Pairs that a band-narrowing (nend one column too far) and a band-trim (F
    // left untrimmed) off-by-one each change, found by mutation testing.
    static const DirectedCase kCases[] = {
        {1, 4, 6, 6, 1, 1, 5, 100, 30, 133,
         "11130021110233310201330200223200002022101023333300201210102023200202110031213213"
         "13301321331330023102011233201002310",
         "21203031020323121002021302103021102212101332121132100212010022100313120201301130"
         "0102220033102312013111031200333"},
        {1, 6, 1, 10, 1, 3, 8, 14, 22, 223,
         "20122003110133410121302031032420312211302133003220", "003321320021433100302001030333"},
    };
    for (const DirectedCase &c : kCases) check_directed(c);
}

// The two cases below compare against the original's ungated epilogue, which
// runs the z-drop step on every row. On Apple silicon the original copy is the
// Apple-gated epilogue instead, whose need_z gate assumes dif >= 0 just as an
// ungated lean gate would, so there it is not a reference for these inputs.
#if !defined(__APPLE__)
TEST_CASE("BSW8_ROW_LEAN: lean == original on pairs whose z-drop drift wraps int32"
          * doctest::test_suite("unit/bandedswa")) {
    // Pairs where |drift| * e wraps to a negative dif in the z-drop step, so the
    // original kills the lane at a drop that does not exceed zdrop: a gate on
    // drop > zdrop alone would keep it alive and report a longer alignment. Both
    // gap extends huge, then only the insertion one, then only the deletion one.
    static const DirectedCase kCases[] = {
        {1, 0, 4, 2, 1073741824, 1073741825, 9, 1, 90, 166,
         "03301310222100301322122030211200123010030111203131210032120231212210231000203022"
         "10012323133133301001320230133101332312032032102001023313220230012311223022131212"
         "1300330011221",
         "121103333213222032232313010002222303313223102302000"},
        {3, 1, 4, 2, 1073741825, 2147483647, 11, 4, 25, 167,
         "12112121310033233023002321211031221102103002031202013", "2132301323"},
        {5, 1, 1, 7, 1, 1073741825, 9, 4, 52, 11,
         "20300200223333123312300133310101301002303303203232321330300223133221302131210311"
         "320320120112032231100110130202102211",
         "0003301032230323130113313103330111"},
        {5, 6, 3, 0, 1073741825, 4, 6, 4, 102, 100,
         "20333211321030232101322320213210031003001012200333030223033232000013311323200223"
         "133231303310130321303123030202130033322321221033211023110113100202021102",
         "0231201103330113002100301"},
    };
    for (const DirectedCase &c : kCases) check_directed(c);
}
TEST_CASE("BSW8_ROW_LEAN: lean == original on every field, random huge or negative gap extends"
          * doctest::test_suite("unit/bandedswa")) {
    // The z-drop step forms |drift| * e in int32, as the original does; for these
    // extends the lean gate steps aside (zgate_exact). Random batches rarely hit a
    // wrap that changes a result, so the directed pairs above are what pin it;
    // this case covers everything else such extends reach.
    check_mode(Mode::HugeExtend, 0x5EB8A19ull, 60);
}
#endif

#else

TEST_CASE("BSW8_ROW_LEAN: not built on this architecture (x86 runs only the original code)"
          * doctest::test_suite("unit/bandedswa")) {
    MESSAGE("SKIP: BSW8_ROW_LEAN is aarch64-only; nothing to compare");
}

#endif
