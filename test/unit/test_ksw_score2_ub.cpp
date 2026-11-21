// test/unit/test_ksw_score2_ub.cpp
//
// kswr_t::score2_ub: an upper bound on every column maximum at or above minsc that ends outside
// the primary zone [te - Z, te + Z], Z = ceil(score / max). score2 is not one: its b[] build keeps
// one entry per run of adjacent columns, so a column outside the zone can be dropped in favour of a
// higher neighbour inside it (the neighbour then excluded as in-zone). The --meth-mapq
// conversion-blind rescue gate bounds the mirror-free second-best from this value, so it must cover
// those columns.
//
// A scalar full-window DP gives the true column maxima; every rescue kernel (scalar ksw_align2 at
// 8 and 16 bits, kswv getScores8 and getScores16) must report score2_ub >= the true out-of-zone
// maximum, on generated windows in which a second copy of the query ends near the zone edges.
// The generated set must contain windows where that maximum is above score2, so the test fails
// if score2_ub is only score2.

#include <algorithm>
#include <cstdint>
#include <random>
#include <string>
#include <vector>

#include "doctest/doctest.h"

#include "ksw_runner.h"
#include "kswv_runner.h"
#include "scoring.h"
#include "seqpair.h"

namespace {

constexpr int kMinsc = 14;
constexpr int kOpen = 6, kExt = 1;

struct Truth {
    int score = 0, te = -1;
    int score2 = -1;      // ksw's merged score2
    int out_max = -1;     // the best column maximum >= minsc outside the zone, unmerged
};

// Smith-Waterman over the whole window (affine gaps, the first gap base costs open + extend, as
// in ksw), the column maximum per target position, and ksw's te (first column reaching the best).
// The kernels pad the query to a multiple of their lane count with cells that score 0 against any
// base, so an alignment reaching the query end carries its score into the next columns; `width`
// (16 for 8-bit scores, 8 for 16-bit) pads the query the same way.
Truth truth_of(const bwa_tests::TestPair &p, const bwa_tests::ScoringMatrix &mat, int width)
{
    const int qlen0 = (int)p.qry.size(), tlen = (int)p.ref.size();
    const int qlen = (qlen0 + width - 1) / width * width;
    std::vector<int> H(qlen + 1, 0), Hn(qlen + 1, 0), E(qlen + 1, 0), col(tlen, 0);
    Truth t;
    int amax = 0;
    for (int8_t v : mat) amax = std::max(amax, (int)v);
    for (int i = 0; i < tlen; ++i) {
        int f = 0, m = 0;
        Hn[0] = 0;
        for (int j = 1; j <= qlen; ++j) {
            int h = H[j - 1] + (j > qlen0 ? 0 : mat[p.ref[i] * 5 + p.qry[j - 1]]);
            h = std::max({ h, E[j], f, 0 });
            Hn[j] = h;
            m = std::max(m, h);
            E[j] = std::max({ E[j] - kExt, h - kOpen - kExt, 0 });
            f = std::max({ f - kExt, h - kOpen - kExt, 0 });
        }
        std::swap(H, Hn);
        col[i] = m;
        if (m > t.score) t.score = m, t.te = i;
    }
    const int z = (t.score + amax - 1) / amax, low = t.te - z, high = t.te + z;
    int bs = -1, bp = -2;
    auto flush = [&] {
        if (bp >= 0 && (bp < low || bp > high) && bs > t.score2) t.score2 = bs;
    };
    for (int i = 0; i < tlen; ++i) {
        if (col[i] < kMinsc) continue;
        if (i < low || i > high) t.out_max = std::max(t.out_max, col[i]);
        if (bp + 1 != i) flush(), bs = col[i], bp = i;
        else if (bs < col[i]) bs = col[i], bp = i;
    }
    flush();
    return t;
}

std::vector<uint8_t> mutate(std::mt19937 &rng, std::vector<uint8_t> s, int pct)
{
    for (auto &b : s)
        if ((int)(rng() % 100) < pct) b = (uint8_t)((b + 1 + rng() % 3) & 3);
    return s;
}

// A window holding a copy of the query and a second, partial copy that ends close to one edge of
// the best's zone, on either side, so that a run of adjacent columns can cross the edge.
std::vector<bwa_tests::TestPair> windows(uint64_t seed, int n)
{
    std::mt19937 rng((unsigned)seed);
    auto ri = [&](int lo, int hi) { return lo + (int)(rng() % (unsigned)(hi - lo + 1)); };
    std::vector<bwa_tests::TestPair> out;
    for (int k = 0; k < n; ++k) {
        bwa_tests::TestPair p;
        const int qlen = ri(30, 80), tlen = ri(3 * qlen, 4 * qlen);
        p.qry.resize(qlen);
        p.ref.resize(tlen);
        for (auto &b : p.qry) b = (uint8_t)(rng() & 3);
        for (auto &b : p.ref) b = (uint8_t)(rng() & 3);
        const std::vector<uint8_t> c1 = mutate(rng, p.qry, ri(0, 6));
        const int p1 = ri(qlen, tlen - 2 * qlen);
        std::copy(c1.begin(), c1.end(), p.ref.begin() + p1);
        // the second copy: a prefix or suffix of the query, its end placed about one zone away
        const int len2 = ri(kMinsc + 4, qlen - 2);
        const bool prefix = rng() & 1;
        std::vector<uint8_t> c2(p.qry.begin() + (prefix ? 0 : qlen - len2),
                                p.qry.begin() + (prefix ? len2 : qlen));
        c2 = mutate(rng, c2, ri(0, 10));
        const int te = p1 + qlen - 1, side = (rng() & 1) ? 1 : -1;
        const int end2 = te + side * (qlen + ri(-4, 4));
        const int st2 = end2 - len2 + 1;
        if (st2 < 0 || end2 >= tlen) { --k; continue; }
        if (st2 <= p1 + qlen - 1 && end2 >= p1) { --k; continue; } // keep the copies apart
        std::copy(c2.begin(), c2.end(), p.ref.begin() + st2);
        out.push_back(std::move(p));
    }
    return out;
}

// score2_ub >= the out-of-zone maximum in every window; and the kernel's score, te and score2 are
// the reference's, so the reference models this kernel's columns.
void check(const std::string &name, const std::vector<bwa_tests::TestPair> &pairs,
           const std::vector<Truth> &truth, const std::vector<kswr_t> &got)
{
    REQUIRE(got.size() == pairs.size());
    long bad = 0, differ = 0;
    for (size_t i = 0; i < pairs.size(); ++i) {
        if (got[i].score == 255) continue; // saturated
        if (got[i].score != truth[i].score || got[i].te != truth[i].te || got[i].score2 != truth[i].score2) {
            if (++differ <= 3)
                MESSAGE(name << " pair " << i << ": score/te/score2 " << got[i].score << "/" << got[i].te << "/"
                             << got[i].score2 << ", reference " << truth[i].score << "/" << truth[i].te << "/"
                             << truth[i].score2);
            continue;
        }
        if (truth[i].out_max > got[i].score2_ub) {
            if (++bad <= 5)
                MESSAGE(name << " pair " << i << ": out-of-zone max " << truth[i].out_max
                              << " > score2_ub " << got[i].score2_ub << " (score2 " << got[i].score2
                              << ", score " << got[i].score << ", te " << got[i].te << ")");
        }
    }
    CHECK_MESSAGE(differ == 0, name << ": " << differ << " windows differ from the reference DP");
    CHECK_MESSAGE(bad == 0, name << ": score2_ub below the true out-of-zone maximum in " << bad
                                 << " windows");
}

} // namespace

TEST_CASE("score2_ub bounds every out-of-zone column maximum, merged or not")
{
    const auto mat = bwa_tests::build_scoring_matrix(1, 4, 1);
    const auto pairs = windows(20261123, 4000);
    std::vector<Truth> t8, t16; // the 8-bit kernels' (16-lane query pad) and the 16-bit ones'
    long masked = 0;
    for (const auto &p : pairs) {
        t8.push_back(truth_of(p, mat, 16));
        t16.push_back(truth_of(p, mat, 8));
        masked += t8.back().out_max > t8.back().score2;
    }
    // the windows must exercise the merge: some out-of-zone column above score2, which score2 as
    // the bound would miss
    MESSAGE(masked << " of " << pairs.size() << " windows have an out-of-zone column above score2");
    REQUIRE(masked > 0);

    const int xtra = KSW_XSUBO | KSW_XSTART | kMinsc;
    std::vector<kswr_t> u8, i16;
    for (const auto &p : pairs) {
        u8.push_back(bwa_tests::run_scalar_ksw(p, mat, kOpen, kExt, xtra | KSW_XBYTE));
        i16.push_back(bwa_tests::run_scalar_ksw(p, mat, kOpen, kExt, xtra));
    }
    check("ksw_align2 u8", pairs, t8, u8);
    check("ksw_align2 i16", pairs, t16, i16);
#if BWA_TESTS_HAVE_KSWV
    check("kswv u8", pairs, t8, bwa_tests::run_kswv_batch(pairs, mat, kOpen, kExt, xtra | KSW_XBYTE, false));
    check("kswv i16", pairs, t16, bwa_tests::run_kswv_batch(pairs, mat, kOpen, kExt, xtra, true));
#endif
}
