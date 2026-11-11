// Unit tests for the --meth refutations (src/aln/aln_meth_prove.h): the text masks against the index
// text, the banded edit distance against full dynamic programming, and the row and window
// refutations against exhaustive edit distances on the index text.

#include "doctest/doctest.h"

#include "aln_meth_prove_fixture.h"

#include <vector>

namespace {

namespace ba = bwa_aln;
using aln_prove_test::full_ed;
using aln_prove_test::prove_index;
using aln_prove_test::Rng;

// A read-like string from T[y, y + n) with e random unit edits (mismatch, insertion or deletion).
std::vector<uint8_t> mutate(const std::vector<uint8_t> &text, int64_t y, int n, int e, Rng &rng)
{
    std::vector<uint8_t> s(text.begin() + y, text.begin() + y + n);
    for (int k = 0; k < e && s.size() > 2; ++k) {
        const int at = rng.below((int)s.size());
        switch (rng.below(3)) {
        case 0: s[at] = (uint8_t)((s[at] + 1 + rng.below(3)) & 3); break;
        case 1: s.insert(s.begin() + at, (uint8_t)rng.below(4)); break;
        default: s.erase(s.begin() + at); break;
        }
    }
    return s;
}

} // namespace

TEST_SUITE_BEGIN("unit/aln_meth_prove");

TEST_CASE("aln_meth_prove: every text mask holds the base the converted index holds there")
{
    auto &x = prove_index();
    const int64_t n = (int64_t)x.bwt->seq_len;
    std::vector<uint8_t> m(n);
    ba::aln_meth_text_masks(x.ref, 0, n, m.data());
    int64_t wide = 0;
    for (int64_t y = 0; y < n; ++y) {
        REQUIRE((m[y] >> x.text[y] & 1) == 1);
        wide += __builtin_popcount(m[y]) > 1;
    }
    // the masks are exact outside the converted bases, the soft-masked stretch and the N run
    CHECK(wide < n / 2);
    // any window agrees with the whole-text masks (contig, strand and junction boundaries)
    Rng rng(7);
    for (int t = 0; t < 2000; ++t) {
        const int64_t y0 = rng.below((int)n), len = rng.below(300);
        const int64_t y1 = std::min(n, y0 + len);
        std::vector<uint8_t> w(y1 - y0);
        ba::aln_meth_text_masks(x.ref, y0, y1, w.data());
        for (int64_t y = y0; y < y1; ++y) REQUIRE(w[y - y0] == m[y]);
    }
}

TEST_CASE("aln_meth_prove: the banded edit distance equals full dynamic programming within its band")
{
    Rng rng(11);
    std::vector<int16_t> dp;
    for (int t = 0; t < 4000; ++t) {
        const int n = 1 + rng.below(30), m = rng.below(40), kmax = rng.below(8);
        std::vector<uint8_t> q(n), tc(m), tm(m);
        for (auto &c : q) c = (uint8_t)(rng.below(20) == 0 ? 4 : rng.below(4));
        for (int b = 0; b < m; ++b) { tc[b] = (uint8_t)rng.below(4); tm[b] = (uint8_t)(1 << tc[b]); }
        if (rng.below(2) && n <= m) // make it close: t derived from q
            for (int b = 0; b < n; ++b) if (q[b] < 4 && rng.below(4)) { tc[b] = q[b]; tm[b] = (uint8_t)(1 << q[b]); }
        for (bool free_start : {false, true}) {
            const int want = full_ed(q, tc, !free_start, true);
            // a band of every diagonal: exact up to kmax
            const int got = ba::aln_meth_band_ed(q.data(), n, tm.data(), m, -n, m, free_start, kmax, &dp);
            CHECK(got == std::min(want, kmax + 1));
            // any band: never above the true distance when that is within kmax
            const int lo = -rng.below(n + 1), hi = rng.below(m + 1);
            const int banded = ba::aln_meth_band_ed(q.data(), n, tm.data(), m, lo, hi, free_start, kmax, &dp);
            CHECK(banded >= std::min(want, kmax + 1));
        }
    }
}

TEST_CASE("aln_meth_prove: a window refutation never refutes an alignment that starts in its band")
{
    auto &x = prove_index();
    const int64_t n = (int64_t)x.bwt->seq_len;
    Rng rng(17);
    ba::aln_meth_verify_scratch_t s;
    int found = 0, refuted = 0;
    for (int t = 0; t < 3000; ++t) {
        const int len = 40 + rng.below(80), kmax = 2 + rng.below(5);
        const int64_t y = rng.below((int)(n - 200));
        std::vector<uint8_t> seq = mutate(x.text, y, len, rng.below(kmax + 3), rng);
        // the window's start: the true start, shifted by up to kmax either way (as a piece's diagonal is)
        const int64_t y0 = y - kmax + rng.below(2 * kmax + 1);
        const int got = ba::aln_meth_refute_window(x.ref, n, seq.data(), (int)seq.size(), y0, kmax, &s);
        // exhaustive: an alignment of seq starting in [y0, y0 + 2 kmax]
        bool any = false;
        for (int64_t st = std::max<int64_t>(0, y0); st <= y0 + 2 * kmax && !any; ++st) {
            const int64_t end = std::min<int64_t>(n, st + (int64_t)seq.size() + kmax);
            std::vector<uint8_t> tt(x.text.begin() + st, x.text.begin() + end);
            any = full_ed(seq, tt, true, true) <= kmax;
        }
        if (any) CHECK(got == 0);
        (got ? refuted : found) += 1;
    }
    CHECK(found > 100);
    CHECK(refuted > 100);
}

TEST_SUITE_END();
