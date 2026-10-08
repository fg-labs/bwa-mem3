// test/unit/test_bandedswa_lane_indep.cpp
//
// Lane independence of the banded extension kernels (getScores8 / getScores16,
// bandedSWA.cpp). A pair's result must not depend on which other pairs share its
// SIMD lane group: scored inside any batch, every pair must return what it
// returns scored alone, and both must equal the scalar reference
// (scalarBandedSWA, bwa's ksw_extend2), on all six result fields. See the
// lane-independence note at the top of bandedSWA.cpp for the mechanism.
//
// Kernel copies checked: on x86 every tier the CPU supports (sse41, sse42, avx,
// avx2, avx512bw; the ones it lacks are reported); on aarch64 NEON with both
// BSW8_ROW_LEAN settings (libbwa.a and src/bandedSWA.rowalt.o; the second with
// lane compaction off, since the compaction driver runs the same lean row in
// both copies). The tiers with same-row lane compaction in getScores8 (NEON,
// avx512bw) are checked three ways: as shipped (the tier's default setting),
// with compaction off, and compacting every band in superblocks of three
// groups.
//
//   1. The regression guard: directed pairs that were coupled to their group in
//      the old kernels (found by fuzzing them; among them the 16-bit corner
//      len1=8 len2=28 h0=29 w=20), two pairs that pin the band trim, and one
//      that pins the phantom column's zero gap.
//      Each is scored alone and next to a full group of wide partners at the
//      first, a middle and the last slot. The old kernels fail here, and so does
//      the fix with either of its two parts removed (the seed-row write-back, or
//      the zero gap at the phantom column).
//   2. A broader backstop: generated batches (random and asymmetric scoring,
//      negative seed scores, N bases, small z-drop, bands 0 to 124) scored in the
//      generated order, reversed, shuffled and split into random sub-batches, and
//      a sample next to wide partners. Coupled pairs are rare, so these batches
//      are not what pins the fix; they keep the property honest on many shapes.

#include <algorithm>
#include <cstdint>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "doctest/doctest.h"
#include "bandedSWA.h"
#include "bsw_batch.h"
#include "ext_ladder.h"

#if HAVE_BSW_VECTOR_8_16

namespace {

using bwa_tests::ExtPair;     // query, target, h0
using bwa_tests::ExtResult;   // the six result fields, operator==
using bwa_tests::ExtScoring;  // a, b, gap open / extend, pen_clip (end bonus), zdrop, mat

/// The scoring with its matrix filled; `asym` additionally scores reference C
/// against read T as a match (an OT-like --meth matrix).
ExtScoring scoring(int a, int b, int o_del, int e_del, int o_ins, int e_ins, int zdrop, int end_bonus, bool asym)
{
    ExtScoring sc;
    sc.a = a; sc.b = b; sc.o_del = o_del; sc.e_del = e_del; sc.o_ins = o_ins; sc.e_ins = e_ins;
    sc.zdrop = zdrop; sc.pen_clip = end_bonus;
    sc.fill_mat();
    if (asym) sc.mat[1 * 5 + 3] = (int8_t)a;
    return sc;
}

std::string show(const ExtResult &o)
{
    return std::to_string(o.score) + "/" + std::to_string(o.tle) + "/" + std::to_string(o.qle) + "/" +
           std::to_string(o.gscore) + "/" + std::to_string(o.gtle) + "/" + std::to_string(o.max_off);
}

/// A kernel copy and the name it is reported under. The kernel keeps a pointer
/// to the scoring's matrix, so the ExtScoring must outlive it.
struct Kernel {
    std::string name;
    std::unique_ptr<IBandedPairWiseSW> (*make)(const ExtScoring &);
};

template <IBandedPairWiseSW *(*Factory)(int, int, int, int, int, int, const int8_t *, int8_t, int8_t, int)>
std::unique_ptr<IBandedPairWiseSW> make_with(const ExtScoring &sc)
{
    return std::unique_ptr<IBandedPairWiseSW>(Factory(sc.o_del, sc.e_del, sc.o_ins, sc.e_ins, sc.zdrop,
                                                      sc.pen_clip, sc.mat, (int8_t)sc.a, (int8_t)sc.b, 1));
}

/// make_with with getScores8's lane compaction set to `groups` / `min_w`.
template <IBandedPairWiseSW *(*Factory)(int, int, int, int, int, int, const int8_t *, int8_t, int8_t, int),
          int groups, int min_w>
std::unique_ptr<IBandedPairWiseSW> make_compact(const ExtScoring &sc)
{
    std::unique_ptr<IBandedPairWiseSW> k = make_with<Factory>(sc);
    k->set_lane_compaction(groups, min_w);
    return k;
}

#if defined(__aarch64__)
IBandedPairWiseSW *make_lib_raw(int o_del, int e_del, int o_ins, int e_ins, int zdrop, int end_bonus,
                                const int8_t *mat, int8_t a, int8_t b, int nthreads)
{
    return new BandedPairWiseSW(o_del, e_del, o_ins, e_ins, zdrop, end_bonus, mat, a, b, nthreads);
}
#endif

/// Every kernel copy this host can run; the x86 tiers it cannot run are reported.
std::vector<Kernel> host_kernels()
{
    std::vector<Kernel> ks;
#if defined(__aarch64__)
    ks.push_back({"neon", make_with<make_lib_raw>});
    ks.push_back({"neon (lane compaction off)", make_compact<make_lib_raw, 0, 0>});
    ks.push_back({"neon (lane compaction, 3 groups, every band)", make_compact<make_lib_raw, 3, 0>});
    ks.push_back({"neon (other BSW8_ROW_LEAN setting)", make_compact<make_bsw_kernel_rowalt, 0, 0>});
#elif defined(__x86_64__) || defined(__i386__)
    __builtin_cpu_init();
    std::string skipped;
    auto add = [&](bool supported, const char *name, std::unique_ptr<IBandedPairWiseSW> (*make)(const ExtScoring &)) {
        if (supported) ks.push_back({name, make});
        else skipped += std::string(" ") + name;
    };
    add(__builtin_cpu_supports("sse4.1"), "sse41", make_with<make_bsw_kernel_sse41>);
    add(__builtin_cpu_supports("sse4.2"), "sse42", make_with<make_bsw_kernel_sse42>);
    add(__builtin_cpu_supports("avx"), "avx", make_with<make_bsw_kernel_avx>);
    add(__builtin_cpu_supports("avx2"), "avx2", make_with<make_bsw_kernel_avx2>);
    add(__builtin_cpu_supports("avx512bw"), "avx512bw", make_with<make_bsw_kernel_avx512bw>);
    add(__builtin_cpu_supports("avx512bw"), "avx512bw (lane compaction off)",
        make_compact<make_bsw_kernel_avx512bw, 0, 0>);
    add(__builtin_cpu_supports("avx512bw"), "avx512bw (lane compaction, 3 groups, every band)",
        make_compact<make_bsw_kernel_avx512bw, 3, 0>);
    if (!skipped.empty()) MESSAGE("lane independence: this CPU cannot run tier(s)" << skipped << "; not checked here");
#endif
    return ks;
}

/// The largest query a perfect-match partner may have on the 8-bit kernel (the
/// routing envelope: h0 + len2 * a stays below the byte ceiling).
int partner_len2_cap8(const ExtScoring &sc, int h0) { return (250 - h0 - 2 * sc.a) / sc.a; }

/// A long perfect-match pair whose band stays wide for the whole run, so a group
/// that holds it has a column range well past the band of a divergent pair. The
/// target is longer than every pair's, and so is the query where the 8-bit
/// envelope allows. A higher seed on the 16-bit kernel keeps its band from
/// collapsing before the pair's does.
ExtPair wide_partner(std::mt19937_64 &rng, const ExtScoring &sc, const ExtPair &p, int width)
{
    ExtPair x;
    x.h0 = width == 8 ? 5 : 50;
    int len2 = (int)std::max(p.query.size(), p.target.size()) + 24;
    if (width == 8) len2 = std::min(len2, partner_len2_cap8(sc, x.h0));
    const int len1 = std::max<int>(len2, (int)p.target.size() + 24);
    x.target.resize(len1);
    for (auto &c : x.target) c = (uint8_t)(rng() % 4);
    x.query.assign(x.target.begin(), x.target.begin() + len2);
    return x;
}

/// Score `p` alone and next to wide partners filling a 64-lane (8-bit) / 32-lane
/// (16-bit) batch, the widest group of any tier, at the first, a middle and the
/// last slot. Each must equal `want`; returns the number of disagreements and
/// reports the first few.
int check_against_partners(IBandedPairWiseSW &bsw, const std::string &kname, const ExtScoring &sc,
                           const ExtPair &p, int width, int w, const ExtResult &want, std::mt19937_64 &rng)
{
    int bad = 0;
    auto expect = [&](const ExtResult &got, const std::string &how) {
        if (got == want) return;
        if (bad++ < 3)
            MESSAGE(kname << " getScores" << width << " " << how << ": " << show(got) << " vs expected "
                          << show(want) << " (len1=" << p.target.size() << " len2=" << p.query.size()
                          << " h0=" << p.h0 << " w=" << w << "; score/tle/qle/gscore/gtle/max_off)");
    };
    expect(bwa_tests::score_bsw_batch(bsw, {&p}, width, w)[0], "alone");
    const int lanes = width == 8 ? 64 : 32;
    std::vector<ExtPair> partners;
    for (int q = 0; q < lanes - 1; q++) partners.push_back(wide_partner(rng, sc, p, width));
    for (int slot : {0, lanes / 2 - 3, lanes - 1}) {
        std::vector<const ExtPair *> batch;
        for (const ExtPair &x : partners) batch.push_back(&x);
        batch.insert(batch.begin() + slot, &p);
        expect(bwa_tests::score_bsw_batch(bsw, batch, width, w)[slot],
               "next to wide partners, slot " + std::to_string(slot));
    }
    return bad;
}

std::vector<uint8_t> digits(const char *s)
{
    std::vector<uint8_t> v;
    for (; *s; s++) v.push_back((uint8_t)(*s - '0'));
    return v;
}

/// A directed pair with its scalar result.
struct Directed {
    int a, b, o_del, e_del, o_ins, e_ins, zdrop, end_bonus;
    bool asym;
    int w, width, h0;
    const char *ref, *qry;
    int score, tle, qle, gscore, gtle, max_off;   // scalarBandedSWA's result
};

// Bases as digits (0-3 = ACGT, 4 = N).
const Directed kDirected[] = {
    // 8-bit pairs that were coupled to their group
    {1, 4, 6, 1, 8, 1, 50, 20, true, 5, 8, 22,
     "013022211210002112",
     "23122313233300030",
     22, 0, 0, 0, 16, 0},
    {1, 4, 6, 1, 6, 1, 100, 5, false, 4, 8, 49,
     "03121020332203343331233433010310220322031131",
     "0332400202422133431324332200212213231300",
     51, 2, 2, 0, 36, 0},
    {1, 4, 6, 1, 6, 1, 100, 5, false, 16, 8, 69,
     "233232322100130014121243200231132303202122230113300031301231202121",
     "11313111313113002043314300013113232011112103113330003130",
     69, 0, 0, 0, 41, 0},
    {1, 4, 6, 1, 6, 1, 100, 5, false, 16, 8, 85,
     "11102232113210012021231022340101212222231323112324130201303021042122110230301041333002032213331210202302132201023",
     "31102231033241002013212022130000002023222333100312130203302210142123103",
     87, 7, 7, 5, 55, 0},
    {1, 8, 1, 2, 11, 2, 200, 5, true, 2, 8, 116,
     "1001100020301123213212123313020313033223342240102211213313042201120023023000132211404201042040023011100300011411132123310333002303220343311303",
     "3321121330023143300200003214123313033031",
     117, 1, 1, 13, 42, 0},
    // 8-bit pairs that pin the band trim (found against a trim variant that
    // changed the pair's result even alone)
    {2, 9, 10, 3, 11, 1, 200, 0, false, 20, 8, 51,
     "21130133311133100132224132202102231321102310102203",
     "03100332311333202013",
     51, 0, 0, 23, 9, 0},
    {1, 5, 9, 3, 11, 1, 150, 5, false, 50, 8, 69,
     "1022111222310220213131023002123311302",
     "131330112421014103313133433212341231",
     70, 1, 1, 18, 26, 0},
    // An 8-bit pair whose result depends on the phantom column's gap output being
    // stored as 0 (here even its score moves next to wide partners)
    {2, 9, 5, 1, 9, 1, 200, 100, true, 2, 8, 21,
     "2033000121330202100134311131204313033311103003210231113203313403303200202110231033240034233024330212414202023332003231130041112042313040302332123034331032023302300340",
     "133303012333202023401324313314244113043331101003030131330003334034020420214023103342400333",
     21, 0, 0, -1, 0, 0},
    // 16-bit pairs that were coupled to their group; the first is the corner
    // len1=8 len2=28 h0=29 w=20 (query longer than target)
    {1, 4, 6, 1, 6, 1, 100, 5, false, 20, 16, 29,
     "30003001",
     "0002300142323301230012003002",
     29, 0, 0, 2, 8, 0},
    {1, 4, 6, 1, 6, 1, 100, 5, false, 20, 16, 29,
     "11124202",
     "0112202213232130131330123400",
     29, 0, 0, 0, 8, 0},
    {2, 7, 12, 2, 12, 1, 200, 5, false, 2, 16, 95,
     "3222310013022012121110032111131123031320221302322232110213203202332312301330113321000131131221021332122321201223131301222200230013230103201200200001103021033303032033020010320011232233003301321133203232113313121022030300133231331222120101103002233121123202223310012111223230033321320313121303030113311102213310013213022001312111013333133301013301001333310023320222111311110110303121",
     "3021331111220020130011032010131023031333020",
     97, 1, 1, 12, 42, 0},
};

/// Generated pairs for one batch: three length regimes, 0-70 % substitutions,
/// up to 5 % indels and 8 % N, an occasional negative seed score. 8-bit pairs
/// stay inside the routing envelope (len1 >= len2, h0 + len2 * a below the byte
/// ceiling).
std::vector<ExtPair> make_batch(std::mt19937_64 &rng, const ExtScoring &sc, bool asym, int width, int n)
{
    std::vector<ExtPair> v(n);
    const int mode = (int)(rng() % 3);
    for (ExtPair &p : v) {
        int len1, len2;
        if (width == 8) {
            for (;;) {
                if (mode == 0) { len2 = 1 + (int)(rng() % 150); len1 = len2 + (int)(rng() % 120); }
                else if (mode == 1) { len2 = 1 + (int)(rng() % (200 / sc.a)); len1 = len2 + (int)(rng() % 400); }
                else { len2 = 1 + (int)(rng() % 40); len1 = len2 + (int)(rng() % 10); }
                const int hmax = 250 - 2 * sc.a - len2 * sc.a;
                if (hmax < 1) continue;
                p.h0 = 1 + (int)(rng() % hmax);
                break;
            }
        } else {
            if (mode == 0) { len1 = 1 + (int)(rng() % 150); len2 = 1 + (int)(rng() % 150); }
            else if (mode == 1) { len1 = 1 + (int)(rng() % 40); len2 = len1 + (int)(rng() % 60); }
            else { len1 = 1 + (int)(rng() % 400); len2 = 1 + (int)(rng() % 250); }
            p.h0 = 1 + (int)(rng() % 400);
        }
        if (rng() % 10 == 0) p.h0 = -(int)(rng() % 30);   // --meth rescored seeds
        const unsigned mm = (unsigned)(rng() % 70), npct = (unsigned)(rng() % 9), indel = (unsigned)(rng() % 6);
        p.target.resize(len1);
        for (auto &c : p.target) c = (uint8_t)(rng() % 100 < npct ? 4 : rng() % 4);
        p.query.resize(len2);
        int t = 0;
        for (int q = 0; q < len2; q++) {
            if (rng() % 100 < indel) {
                if (rng() & 1) t++;
                else { p.query[q] = (uint8_t)(rng() % 4); continue; }
            }
            p.query[q] = t < len1 ? p.target[t] : (uint8_t)(rng() % 4);
            t++;
            if (rng() % 100 < mm) p.query[q] = (uint8_t)(rng() % 4);
            if (rng() % 100 < npct) p.query[q] = 4;
            if (asym && p.query[q] == 1 && (rng() & 1)) p.query[q] = 3;   // C->T conversions
        }
    }
    return v;
}

/// Generated batches for one kernel width: every grouping of every batch must
/// give each pair its scalar result (or, for a negative seed on the 8-bit
/// kernel, its result scored alone; see below).
void check_generated(int width, uint64_t seed, int nbatch)
{
    static const int bands[] = {0, 1, 2, 5, 10, 20, 50, 100, 124};
    static const int zdrops[] = {1, 5, 20, 50, 100, 150, 200};
    static const int bonuses[] = {0, 5, 20, 100};
    const std::vector<Kernel> kernels = host_kernels();
    REQUIRE(!kernels.empty());
    for (const Kernel &k : kernels) {
        std::mt19937_64 rng(seed);
        long pairs = 0, bad = 0;
        for (int bi = 0; bi < nbatch; bi++) {
            // Every other batch uses bwa's default scoring, where most coupled pairs
            // were found; the rest draw random scoring.
            const bool dflt = bi % 2 == 0;
            const bool asym = !dflt && rng() % 3 == 0;
            const int a = dflt ? 1 : 1 + (int)(rng() % 3), b = dflt ? 4 : 1 + (int)(rng() % 9);
            const int o_del = dflt ? 6 : (int)(rng() % 17), e_del = dflt ? 1 : 1 + (int)(rng() % 4);
            const int o_ins = dflt ? 6 : (int)(rng() % 17), e_ins = dflt ? 1 : 1 + (int)(rng() % 4);
            const int zdrop = dflt ? 100 : zdrops[rng() % 7], end_bonus = dflt ? 5 : bonuses[rng() % 4];
            const ExtScoring sc = scoring(a, b, o_del, e_del, o_ins, e_ins, zdrop, end_bonus, asym);
            const int w = bands[rng() % 9];
            const std::vector<ExtPair> batch = make_batch(rng, sc, asym, width, 1 + (int)(rng() % 120));
            const std::unique_ptr<IBandedPairWiseSW> bsw = k.make(sc);
            // The expected result is the scalar one. For a negative seed the 8-bit
            // kernels seed the byte DP with max(h0, 0) where the scalar code can
            // start from h0 < 0 (PR #528), so there the pair scored alone is the
            // reference: the grouping property still holds exactly.
            std::vector<ExtResult> want;
            for (const ExtPair &p : batch)
                want.push_back(width == 8 && p.h0 < 0 ? bwa_tests::score_bsw_batch(*bsw, {&p}, width, w)[0]
                                                      : bwa_tests::run_mem3_scalar(sc, p, w));
            std::vector<int> order(batch.size());
            for (size_t q = 0; q < order.size(); q++) order[q] = (int)q;
            // generated order, reversed, two shuffles, then random sub-batches
            for (int g = 0; g < 5; g++) {
                if (g == 1) std::reverse(order.begin(), order.end());
                if (g == 2 || g == 3) std::shuffle(order.begin(), order.end(), rng);
                size_t s = 0;
                while (s < order.size()) {
                    const size_t len = g == 4 ? 1 + rng() % 40 : order.size();
                    std::vector<const ExtPair *> sub;
                    for (size_t q = s; q < std::min(order.size(), s + len); q++) sub.push_back(&batch[order[q]]);
                    const std::vector<ExtResult> got = bwa_tests::score_bsw_batch(*bsw, sub, width, w);
                    for (size_t q = 0; q < sub.size(); q++) {
                        const int id = order[s + q];
                        pairs++;
                        if (got[q] == want[id]) continue;
                        if (bad++ < 5)
                            MESSAGE(k.name << " getScores" << width << " batch " << bi << " grouping " << g
                                           << ": " << show(got[q]) << " vs expected " << show(want[id])
                                           << " (len1=" << batch[id].target.size() << " len2="
                                           << batch[id].query.size() << " h0=" << batch[id].h0 << " w=" << w
                                           << ")");
                    }
                    s += sub.size();
                }
            }
            // a sample of the batch next to wide partners
            for (size_t q = 0; q < batch.size() && q < 2; q++)
                bad += check_against_partners(*bsw, k.name, sc, batch[q], width, w, want[q], rng);
        }
        MESSAGE(k.name << " getScores" << width << ": " << pairs << " pair-scorings, " << bad << " disagreements");
        CHECK(pairs > (width == 8 ? 3000 : 2000));   // the generator did not collapse
        CHECK(bad == 0);
    }
}

} // namespace

TEST_CASE("bandedSWA lane independence: formerly coupled pairs score their scalar result in any group"
          * doctest::test_suite("unit/bandedswa")) {
    const std::vector<Kernel> kernels = host_kernels();
    REQUIRE(!kernels.empty());
    for (const Kernel &k : kernels) {
        std::mt19937_64 rng(0x1A9E5EEDull);
        int bad = 0;
        for (const Directed &d : kDirected) {
            const ExtScoring sc = scoring(d.a, d.b, d.o_del, d.e_del, d.o_ins, d.e_ins, d.zdrop, d.end_bonus, d.asym);
            ExtPair p;
            p.target = digits(d.ref);
            p.query = digits(d.qry);
            p.h0 = d.h0;
            ExtResult pinned;
            pinned.score = d.score; pinned.tle = d.tle; pinned.qle = d.qle;
            pinned.gscore = d.gscore; pinned.gtle = d.gtle; pinned.max_off = d.max_off;
            const ExtResult scalar = bwa_tests::run_mem3_scalar(sc, p, d.w);
            if (!(scalar == pinned)) {
                MESSAGE("scalar result " << show(scalar) << " differs from the pinned " << show(pinned));
                bad++;
            }
            const std::unique_ptr<IBandedPairWiseSW> bsw = k.make(sc);
            bad += check_against_partners(*bsw, k.name, sc, p, d.width, d.w, scalar, rng);
        }
        INFO("kernel " << k.name);
        CHECK(bad == 0);
    }
}

TEST_CASE("bandedSWA lane independence: generated 8-bit batches in any grouping"
          * doctest::test_suite("unit/bandedswa")) {
    check_generated(8, 0x1A9E0008ull, 16);
}

TEST_CASE("bandedSWA lane independence: generated 16-bit batches in any grouping"
          * doctest::test_suite("unit/bandedswa")) {
    check_generated(16, 0x1A9E0016ull, 10);
}

#endif // HAVE_BSW_VECTOR_8_16
