// Unit tests for the --meth refutations (src/aln/aln_meth_prove.h): the text masks against the index
// text, the banded edit distance against full dynamic programming, and the row and window
// refutations against exhaustive edit distances on the index text.

#include "doctest/doctest.h"

#include "aln_meth_prove_fixture.h"
#include "../../src/aln/bwtgap.h"

#include <functional>
#include <vector>

namespace {

namespace ba = bwa_aln;
using aln_prove_test::full_ed;
using aln_prove_test::prove_index;
using aln_prove_test::prove_index_plain;
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

// Semi-global unit edit distance against the permissive masks (a read N matches nothing).
// This is a lower bound on the indexed-text distance, not an oracle for actual indexed hits.
static int mask_ed_anywhere(const std::vector<uint8_t> &seq, const std::vector<uint8_t> &masks, int kmax)
{
    const size_t n = seq.size(), m = masks.size();
    std::vector<int> prev(m + 1, 0), cur(m + 1);
    for (size_t a = 1; a <= n; ++a) {
        cur[0] = (int)a;
        const int q = seq[a - 1];
        int rmin = cur[0];
        for (size_t b = 1; b <= m; ++b) {
            cur[b] = std::min({prev[b - 1] + (q < 4 && (masks[b - 1] >> q & 1) ? 0 : 1), prev[b] + 1, cur[b - 1] + 1});
            rmin = std::min(rmin, cur[b]);
        }
        if (rmin > kmax) return kmax + 1;
        std::swap(prev, cur);
    }
    return std::min(*std::min_element(prev.begin(), prev.end()), kmax + 1);
}

// The pieces aln_meth_pigeonhole cuts a read of len bases into for kmax, and their allowances.
static void pieces(int len, int kmax, std::vector<int> *pa, std::vector<int> *pe)
{
    const int n1 = (kmax + 1) / 2, p = n1 + (kmax + 1) % 2;
    pa->clear(); pe->clear();
    for (int j = 0; j <= p; ++j) pa->push_back((int)((int64_t)j * len / p));
    for (int j = 0; j < p; ++j) pe->push_back(j < n1 ? 1 : 0);
}

// A segment of exactly want bases aligned to the text from *y on with d edits (positions in the
// segment chosen at random, the first and last base favoured); *y moves past the text it used.
static std::vector<uint8_t> edited_segment(const std::vector<uint8_t> &text, int64_t *y, int want, int d, Rng &rng)
{
    std::vector<int> type(d);
    int net = 0; // inserted - deleted read bases
    for (int &t : type) { t = rng.below(3); net += t == 1 ? 1 : t == 2 ? -1 : 0; }
    const int take = want - net;
    std::vector<uint8_t> s(text.begin() + *y, text.begin() + *y + take);
    *y += take;
    for (int t : type) {
        const int r = rng.below(6);
        const int at = r == 0 ? 0 : r == 1 ? (int)s.size() - 1 : rng.below((int)s.size());
        if (t == 0) s[at] = (uint8_t)((s[at] + 1 + rng.below(3)) & 3);
        else if (t == 1) s.insert(s.begin() + at + (rng.below(2) && at + 1 <= (int)s.size() ? 1 : 0), (uint8_t)rng.below(4));
        else s.erase(s.begin() + at);
    }
    REQUIRE((int)s.size() == want);
    return s;
}

TEST_CASE("aln_meth_prove: a planted hit is found for every split of its edits over the pieces")
{
    auto &x = prove_index();
    const int64_t n = (int64_t)x.bwt->seq_len;
    Rng rng(19);
    ba::aln_meth_ph_scratch_t s;
    int planted = 0;
    for (int kmax = 3; kmax <= 7; ++kmax) {
        const int len = 40 + 16 * kmax;
        std::vector<int> pa, pe;
        pieces(len, kmax, &pa, &pe);
        const int p = (int)pe.size();
        // every composition of kmax edits into p pieces (d_0 + ... + d_{p-1} = kmax)
        std::vector<int> d(p, 0);
        std::function<void(int, int)> each = [&](int j, int left) {
            if (j == p - 1) {
                d[j] = left;
                for (int rep = 0; rep < 3; ++rep) {
                    int64_t y = 50 + rng.below((int)(n - len - 200));
                    std::vector<uint8_t> seq;
                    for (int q = 0; q < p; ++q) {
                        std::vector<uint8_t> seg = edited_segment(x.text, &y, pa[q + 1] - pa[q], d[q], rng);
                        seq.insert(seq.end(), seg.begin(), seg.end());
                    }
                    REQUIRE((int)seq.size() == len);
                    long steps = 0, rows = 0;
                    CHECK(ba::aln_meth_pigeonhole(x.bwt, x.ref, seq.data(), len, kmax, nullptr, &s, &steps, &rows) == 0);
                    CHECK(steps > 0);
                    ++planted;
                }
                return;
            }
            for (int v = 0; v <= left; ++v) { d[j] = v; each(j + 1, left - v); }
        };
        each(0, kmax);
    }
    CHECK(planted > 300);
}

TEST_CASE("aln_meth_prove: mask-only alignments do not prevent an indexed-text no-hit proof")
{
    auto &x = prove_index();
    // Regression: the permissive masks admit five edits, but the actual indexed bases need seven.
    const std::string bases = "TGTTTTTTTGAAAGGTATTGGGGTTATTCTAATAATTGTTGTATTTTTATAATTTAAATATGATATGAGTTTGTTAGATGTTTTGTTTATNTTTTGGGGTAAATTAAACAAATGTGTGTTG";
    std::vector<uint8_t> seq;
    for (char c : bases) seq.push_back(c == 'A' ? 0 : c == 'C' ? 1 : c == 'G' ? 2 : c == 'T' ? 3 : 4);
    std::vector<uint8_t> masks(x.text.size());
    ba::aln_meth_text_masks(x.ref, 0, (int64_t)x.text.size(), masks.data());
    CHECK(mask_ed_anywhere(seq, masks, 5) == 5);
    REQUIRE(full_ed(seq, x.text, false, true) == 7);
    ba::aln_meth_ph_scratch_t s;
    long steps = 0, rows = 0;
    CHECK(ba::aln_meth_pigeonhole(x.bwt, x.ref, seq.data(), (int)seq.size(), 5, nullptr, &s, &steps, &rows) == 1);
    // The same read has a real alignment at seven edits and must not be refuted at that limit.
    CHECK(ba::aln_meth_pigeonhole(x.bwt, x.ref, seq.data(), (int)seq.size(), 7, nullptr, &s, &steps, &rows) == 0);
}

TEST_CASE("aln_meth_prove: the pigeonhole proof never refutes an indexed hit and proves mask no-hits")
{
    auto &x = prove_index();
    const int64_t n = (int64_t)x.bwt->seq_len;
    std::vector<uint8_t> masks(n);
    ba::aln_meth_text_masks(x.ref, 0, n, masks.data());
    Rng rng(23);
    ba::aln_meth_ph_scratch_t s;
    int proven = 0, found = 0;
    for (int t = 0; t < 160; ++t) {
        const int kmax = 3 + rng.below(4), len = 56 + rng.below(70);
        std::vector<uint8_t> seq;
        if (t % 4 == 0) { // random bases: no hit
            for (int b = 0; b < len; ++b) seq.push_back((uint8_t)rng.below(4));
        } else { // from the text with kmax - 1 .. kmax + 2 edits: near misses either side
            const int64_t y = rng.below((int)(n - len - 20));
            seq = mutate(x.text, y, len, kmax - 1 + rng.below(4), rng);
        }
        if (t % 9 == 0) seq[rng.below((int)seq.size())] = 4;
        const int masked_ed = mask_ed_anywhere(seq, masks, kmax);
        const int indexed_ed = full_ed(seq, x.text, false, true);
        long steps = 0, rows = 0;
        const int got = ba::aln_meth_pigeonhole(x.bwt, x.ref, seq.data(), (int)seq.size(), kmax, nullptr, &s, &steps, &rows);
        INFO("t=" << t << " kmax=" << kmax << " indexed_ed=" << indexed_ed << " masked_ed=" << masked_ed);
        CHECK(masked_ed <= indexed_ed);
        // Safety is against the actual index, independently recovered through the suffix array.
        if (indexed_ed <= kmax) CHECK(got == 0);
        if (got) CHECK(indexed_ed > kmax);
        // A mask no-hit is sufficient but not necessary: piece enumeration uses actual bases.
        // No resource cap binds on this small index.
        if (masked_ed > kmax) CHECK(got == 1);
        (got ? proven : found) += 1;
    }
    CHECK(proven > 40);
    CHECK(found > 30);
}

TEST_CASE("aln_meth_prove: on a plain index the masks are the indexed base, or any base in an N run")
{
    auto &x = prove_index_plain();
    const int64_t n = (int64_t)x.bwt->seq_len;
    REQUIRE(n == 2 * (int64_t)x.ref.bns->l_pac);
    std::vector<uint8_t> m(n);
    ba::aln_meth_text_masks(x.ref, 0, n, m.data());
    int64_t any = 0;
    for (int64_t y = 0; y < n; ++y) {
        REQUIRE((m[y] >> x.text[y] & 1) == 1);
        REQUIRE((m[y] == 0xf || __builtin_popcount(m[y]) == 1));
        any += m[y] == 0xf;
    }
    CHECK(any == 2 * 30); // the N run, on both strands
    Rng rng(11);
    for (int t = 0; t < 2000; ++t) {
        const int64_t y0 = rng.below((int)n), len = rng.below(300);
        const int64_t y1 = std::min(n, y0 + len);
        std::vector<uint8_t> w(y1 - y0);
        ba::aln_meth_text_masks(x.ref, y0, y1, w.data());
        for (int64_t y = y0; y < y1; ++y) REQUIRE(w[y - y0] == m[y]);
    }
}

TEST_CASE("aln_meth_prove: plain-index ambiguous masks are not indexed occurrences")
{
    auto &x = prove_index_plain();
    // Derived from text [5480, 5552), changing only bases in the original N run [5500, 5530).
    const std::string bases = "GGACTGATCGTTTGACCTGTGTTTCCGCTGAGAGGAATCTGAGATATTTATCAATAACAGCAGATGCCACTT";
    std::vector<uint8_t> seq;
    for (char c : bases) seq.push_back(c == 'A' ? 0 : c == 'C' ? 1 : c == 'G' ? 2 : 3);
    std::vector<uint8_t> masks(x.text.size());
    ba::aln_meth_text_masks(x.ref, 0, (int64_t)x.text.size(), masks.data());
    const int masked_ed = mask_ed_anywhere(seq, masks, 3);
    const int indexed_ed = full_ed(seq, x.text, false, true);
    long steps = 0, rows = 0;
    ba::aln_meth_ph_scratch_t s;
    const int got = ba::aln_meth_pigeonhole(x.bwt, x.ref, seq.data(), (int)seq.size(), 3, nullptr, &s, &steps, &rows);
    CHECK(masked_ed == 0);
    REQUIRE(indexed_ed == 21);
    CHECK(got == 1);
    // An actual indexed substring crossing that same run must never be refuted.
    seq.assign(x.text.begin() + 5480, x.text.begin() + 5552);
    REQUIRE(full_ed(seq, x.text, false, true) == 0);
    CHECK(ba::aln_meth_pigeonhole(x.bwt, x.ref, seq.data(), (int)seq.size(), 3, nullptr, &s, &steps, &rows) == 0);
}

TEST_CASE("aln_meth_prove: on a plain index the proof never refutes an indexed hit and proves mask no-hits")
{
    auto &x = prove_index_plain();
    const int64_t n = (int64_t)x.bwt->seq_len;
    std::vector<uint8_t> masks(n);
    ba::aln_meth_text_masks(x.ref, 0, n, masks.data());
    Rng rng(29);
    ba::aln_meth_ph_scratch_t s;
    int proven = 0, found = 0;
    for (int t = 0; t < 160; ++t) {
        const int kmax = 3 + rng.below(4), len = 40 + rng.below(80);
        std::vector<uint8_t> seq;
        if (t % 4 == 0) {
            for (int b = 0; b < len; ++b) seq.push_back((uint8_t)rng.below(4));
        } else {
            const int64_t y = rng.below((int)(n - len - 20));
            seq = mutate(x.text, y, len, kmax - 1 + rng.below(4), rng);
        }
        if (t % 9 == 0) seq[rng.below((int)seq.size())] = 4;
        const int masked_ed = mask_ed_anywhere(seq, masks, kmax);
        const int indexed_ed = full_ed(seq, x.text, false, true);
        long steps = 0, rows = 0;
        const int got = ba::aln_meth_pigeonhole(x.bwt, x.ref, seq.data(), (int)seq.size(), kmax, nullptr, &s, &steps, &rows);
        INFO("t=" << t << " kmax=" << kmax << " indexed_ed=" << indexed_ed << " masked_ed=" << masked_ed);
        CHECK(masked_ed <= indexed_ed);
        // Ambiguous text masks are permissive; safety must be checked against actual indexed bases.
        if (indexed_ed <= kmax) CHECK(got == 0);
        if (got) CHECK(indexed_ed > kmax);
        if (masked_ed > kmax) CHECK(got == 1);
        (got ? proven : found) += 1;
    }
    CHECK(proven > 40);
    CHECK(found > 30);
}

TEST_CASE("aln_meth_prove: the proof uses the search's widths for its first piece and agrees without them")
{
    auto &x = prove_index();
    const int64_t n = (int64_t)x.bwt->seq_len;
    Rng rng(29);
    ba::aln_meth_ph_scratch_t s;
    for (int t = 0; t < 200; ++t) {
        const int kmax = 3 + rng.below(4), len = 60 + rng.below(60);
        const int64_t y = rng.below((int)(n - len - 20));
        std::vector<uint8_t> seq = mutate(x.text, y, len, kmax + rng.below(3) - 1, rng);
        const int L = (int)seq.size();
        // the search's widths: bwt_cal_width on the uncomplemented bases
        std::vector<uint8_t> v(L);
        for (int b = 0; b < L; ++b) v[b] = seq[b] > 3 ? 4 : (uint8_t)(3 - seq[b]);
        std::vector<ba::bwt_width_t> w(L + 1);
        ba::bwt_cal_width(x.bwt, L, v.data(), w.data());
        long s1 = 0, r1 = 0, s2 = 0, r2 = 0;
        const int a = ba::aln_meth_pigeonhole(x.bwt, x.ref, seq.data(), L, kmax, w.data(), &s, &s1, &r1);
        const int b = ba::aln_meth_pigeonhole(x.bwt, x.ref, seq.data(), L, kmax, nullptr, &s, &s2, &r2);
        CHECK(a == b);
        CHECK(r1 == r2);
        CHECK(s1 <= s2);
    }
}

TEST_SUITE_END();
