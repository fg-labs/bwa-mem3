// test/unit/test_mem_extend_cigar.cpp
//
// CIGAR, MD and NM generation (bwa_gen_cigar2, src/bwa.cpp), the step that turns
// an alignment region [rb, re) x query into the CIGAR, NM:i and MD:Z of every
// SAM record `mem` writes. Checked against expectations derived here, on a
// packed reference built in the test:
//
//   - the CIGAR is a global alignment of the whole query to the whole region:
//     it consumes exactly l_query query bases and re - rb reference bases, uses
//     only M, I and D, never puts I next to D (ksw_global2's model), and
//     rescored under the scoring it equals the returned score. That score is
//     the optimum of a full-matrix three-state DP written here wherever the
//     band bwa documents for this step holds that optimum: half of the longest
//     gap the read could pay for plus the length difference, at most the
//     caller's band and at least |tlen - qlen| + 3 (bwa's bwa_gen_cigar2);
//   - MD follows the SAM specification over that CIGAR: the run of matching
//     bases, then each mismatching reference base, and "^" plus the deleted
//     reference bases for a deletion; bwa leaves a deletion at either end of
//     the CIGAR out of MD and NM. Given the CIGAR, the query and the
//     reference, MD is fully determined, so it is rebuilt here independently
//     and compared as a string;
//   - NM is the number of mismatching aligned columns (a literal difference;
//     N against a base counts) plus every inserted and every interior deleted
//     base;
//   - on the reverse strand (rb >= l_pac) the region is the reverse complement
//     of the forward reference, and the CIGAR and MD describe the reverse
//     complement of the query against the forward reference, as SAM reports a
//     reverse-strand read;
//   - indels are placed at the leftmost position on the forward reference on
//     both strands (bwa_gen_cigar4: reversing both sequences on the reverse
//     strand "is to ensure indels to be placed at the leftmost position"),
//     checked on homopolymer and dinucleotide-repeat indels whose placement is
//     otherwise free.
//
// Random inputs use std::mt19937 raw output (rng() % n).

#include "doctest/doctest.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <vector>

#include "bwa.h"

namespace {

struct Scoring {
    int a, b, o_del, e_del, o_ins, e_ins;
    int8_t mat[25];
};

Scoring scoring(int a, int b, int o_del, int e_del, int o_ins, int e_ins) {
    Scoring s{a, b, o_del, e_del, o_ins, e_ins, {}};
    for (int i = 0; i < 5; ++i)
        for (int j = 0; j < 5; ++j) s.mat[i * 5 + j] = (int8_t)(i == 4 || j == 4 ? -1 : i == j ? a : -b);
    return s;
}

// A forward reference and its 2-bit packing (bntseq.h: base k sits in bits
// (~k & 3) * 2 of byte k >> 2).
struct Packed {
    std::vector<uint8_t> fwd, pac;
    int64_t l_pac() const { return (int64_t)fwd.size(); }
};

Packed pack(const std::vector<uint8_t> &fwd) {
    Packed p;
    p.fwd = fwd;
    p.pac.assign(fwd.size() / 4 + 1, 0);
    for (size_t k = 0; k < fwd.size(); ++k) p.pac[k >> 2] |= (uint8_t)(fwd[k] << ((~k & 3) << 1));
    return p;
}

// Doubled-coordinate region [rb, re) as `mem` sees it: forward bases below
// l_pac, the reverse complement above.
std::vector<uint8_t> region(const Packed &p, int64_t rb, int64_t re) {
    std::vector<uint8_t> out;
    const int64_t L = p.l_pac();
    for (int64_t x = rb; x < re; ++x)
        out.push_back(x < L ? p.fwd[(size_t)x] : (uint8_t)(3 - p.fwd[(size_t)(2 * L - 1 - x)]));
    return out;
}

std::vector<uint8_t> revcomp(const std::vector<uint8_t> &s) {
    std::vector<uint8_t> r(s.rbegin(), s.rend());
    for (auto &c : r) c = c < 4 ? (uint8_t)(3 - c) : c;
    return r;
}

std::string cigar_string(const uint32_t *cigar, int n) {
    std::string out;
    for (int k = 0; k < n; ++k) out += std::to_string(cigar[k] >> 4) + "MIDNSHP=X"[cigar[k] & 0xf];
    return out;
}

// The band bwa_gen_cigar2 documents for a query of qlen bases against a region
// of rlen under band cap w_ (bwa 0.7.x bwa.c): half of the longest insertion or
// deletion that half the read's matches could pay for, plus |rlen - qlen|, at
// most w_ and at least |rlen - qlen| + 3. A narrower band computed by the code
// is allowed only where it gives the same alignment, so where this band holds
// the unbanded optimum the result must be that optimum.
int documented_band(const Scoring &s, int qlen, int rlen, int w_) {
    const int half = (qlen + 1) >> 1;
    const int max_ins = (int)((double)(half * s.a - s.o_ins) / s.e_ins + 1.);
    const int max_del = (int)((double)(half * s.a - s.o_del) / s.e_del + 1.);
    const int max_gap = std::max(1, std::max(max_ins, max_del));
    const int dlen = std::abs(rlen - qlen);
    return std::max(std::min((max_gap + dlen + 1) >> 1, w_), dlen + 3);
}

// Full-matrix global three-state DP; gaps open only from M (ksw_global2).
// With band >= 0, only cells within `band` of the main diagonal are used.
int global_optimum(const std::vector<uint8_t> &q, const std::vector<uint8_t> &t, const Scoring &s, int band = -1) {
    const int NEG = -(1 << 28);
    const size_t Q = q.size(), T = t.size();
    std::vector<int> M((T + 1) * (Q + 1), NEG), E(M), F(M);
    const size_t W = Q + 1;
    M[0] = 0;
    for (size_t j = 1; j <= Q; ++j) F[j] = -(s.o_ins + (int)j * s.e_ins);
    for (size_t i = 1; i <= T; ++i) E[i * W] = -(s.o_del + (int)i * s.e_del);
    for (size_t i = 1; i <= T; ++i)
        for (size_t j = 1; j <= Q; ++j) {
            if (band >= 0 && (i > j ? i - j : j - i) > (size_t)band) continue;
            const size_t d = (i - 1) * W + j - 1;
            const int h = std::max(M[d], std::max(E[d], F[d]));
            M[i * W + j] = h + s.mat[t[i - 1] * 5 + q[j - 1]];
            E[i * W + j] = std::max(M[(i - 1) * W + j] - s.o_del - s.e_del, E[(i - 1) * W + j] - s.e_del);
            F[i * W + j] = std::max(M[i * W + j - 1] - s.o_ins - s.e_ins, F[i * W + j - 1] - s.e_ins);
        }
    return std::max(M[T * W + Q], std::max(E[T * W + Q], F[T * W + Q]));
}

// What the SAM specification and bwa's conventions make of a CIGAR over
// (query, ref), both in forward-reference orientation.
struct Derived {
    bool valid = false;
    int score = 0, nm = 0;
    std::string md;
};

Derived derive(const uint32_t *cigar, int n, const std::vector<uint8_t> &q, const std::vector<uint8_t> &r,
               const Scoring &s) {
    Derived d;
    size_t x = 0, y = 0;  // query, reference
    int run = 0, prev = -1;
    for (int k = 0; k < n; ++k) {
        const int op = (int)(cigar[k] & 0xf), len = (int)(cigar[k] >> 4);
        if (len <= 0 || op > 2) return d;
        if ((op == 1 && prev == 2) || (op == 2 && prev == 1)) return d;
        if (op == 0) {
            for (int i = 0; i < len; ++i, ++x, ++y) {
                if (x >= q.size() || y >= r.size()) return d;
                d.score += s.mat[r[y] * 5 + q[x]];
                if (q[x] != r[y]) {
                    d.md += std::to_string(run) + "ACGTN"[r[y]];
                    run = 0;
                    ++d.nm;
                } else {
                    ++run;
                }
            }
        } else if (op == 1) {
            x += (size_t)len;
            d.nm += len;
            d.score -= s.o_ins + len * s.e_ins;
        } else {
            if (k > 0 && k < n - 1) {
                d.md += std::to_string(run) + "^";
                for (int i = 0; i < len; ++i) d.md += "ACGTN"[r[y + (size_t)i]];
                run = 0;
                d.nm += len;
            }
            y += (size_t)len;
            d.score -= s.o_del + len * s.e_del;
        }
        prev = op;
    }
    d.md += std::to_string(run);
    d.valid = x == q.size() && y == r.size();
    return d;
}

struct Generated {
    int score = 0, n_cigar = 0, nm = 0;
    std::vector<uint32_t> cigar;
    std::string md;
};

Generated gen(const Packed &p, const Scoring &s, int w, std::vector<uint8_t> query, int64_t rb, int64_t re) {
    Generated g;
    const std::vector<uint8_t> before = query;
    uint32_t *cigar = bwa_gen_cigar2(s.mat, s.o_del, s.e_del, s.o_ins, s.e_ins, w, p.l_pac(), p.pac.data(),
                                     (int)query.size(), query.data(), rb, re, &g.score, &g.n_cigar, &g.nm);
    REQUIRE(cigar != nullptr);
    REQUIRE(g.n_cigar > 0);
    g.cigar.assign(cigar, cigar + g.n_cigar);
    g.md = (const char *)(cigar + g.n_cigar);  // MD is appended after the ops
    free(cigar);
    CHECK(query == before);  // the query is restored after the reverse-strand flip
    return g;
}

std::vector<uint8_t> random_seq(std::mt19937 &rng, size_t n) {
    std::vector<uint8_t> v(n);
    for (size_t k = 0; k < n; ++k) v[k] = (uint8_t)(rng() % 4);
    return v;
}

std::vector<uint8_t> mutate(std::mt19937 &rng, const std::vector<uint8_t> &r) {
    std::vector<uint8_t> q;
    for (size_t k = 0; k < r.size(); ++k) {
        const uint32_t x = rng() % 100;
        if (x < 4) q.push_back((uint8_t)((r[k] + 1 + rng() % 3) % 4));
        else if (x < 5) q.push_back(4);
        else if (x < 7) continue;
        else if (x < 9) {
            for (uint32_t m = 1 + rng() % 3; m > 0; --m) q.push_back((uint8_t)(rng() % 4));
            q.push_back(r[k]);
        } else q.push_back(r[k]);
    }
    if (q.empty()) q.push_back(r[0]);
    return q;
}

std::vector<uint8_t> digits(const char *str) {
    std::vector<uint8_t> v;
    for (; *str; ++str) v.push_back((uint8_t)(*str - '0'));
    return v;
}

}  // namespace

TEST_CASE("bwa_gen_cigar2: hand-worked CIGAR, MD and NM on both strands"
          * doctest::test_suite("unit/cigar")) {
    const Scoring def = scoring(1, 4, 6, 1, 6, 1);
    // 60 bases, no base equal to its neighbour; the region is [10, 50).
    const std::vector<uint8_t> fwd = digits("012310230213012031201320123021031023103201231023021301203120");
    REQUIRE(fwd.size() == 60);
    const Packed p = pack(fwd);
    const std::vector<uint8_t> ref(fwd.begin() + 10, fwd.begin() + 50);

    SUBCASE("exact: 40M, MD 40, NM 0, score 40") {
        const Generated g = gen(p, def, 100, ref, 10, 50);
        CHECK(cigar_string(g.cigar.data(), g.n_cigar) == "40M");
        CHECK(g.md == "40");
        CHECK(g.nm == 0);
        CHECK(g.score == 40);
    }
    SUBCASE("a substitution at query base 15: MD names the reference base") {
        std::vector<uint8_t> q = ref;
        q[15] = (uint8_t)((q[15] + 1) % 4);
        const Generated g = gen(p, def, 100, q, 10, 50);
        CHECK(cigar_string(g.cigar.data(), g.n_cigar) == "40M");
        CHECK(g.md == std::string("15") + "ACGT"[ref[15]] + "24");
        CHECK(g.nm == 1);
        CHECK(g.score == 39 - 4);
    }
    SUBCASE("two adjacent substitutions: a 0 between the two bases") {
        std::vector<uint8_t> q = ref;
        q[20] = (uint8_t)((q[20] + 1) % 4);
        q[21] = (uint8_t)((q[21] + 1) % 4);
        const Generated g = gen(p, def, 100, q, 10, 50);
        CHECK(g.md == std::string("20") + "ACGT"[ref[20]] + "0" + "ACGT"[ref[21]] + "18");
        CHECK(g.nm == 2);
        CHECK(g.score == 38 - 8);
    }
    SUBCASE("an N in the query is a mismatch for MD and NM and scores -1") {
        std::vector<uint8_t> q = ref;
        q[5] = 4;
        const Generated g = gen(p, def, 100, q, 10, 50);
        CHECK(g.md == std::string("5") + "ACGT"[ref[5]] + "34");
        CHECK(g.nm == 1);
        CHECK(g.score == 39 - 1);
    }
    SUBCASE("two reference bases missing from the query: 20M2D18M, MD ^ and the bases") {
        std::vector<uint8_t> q(ref.begin(), ref.begin() + 20);
        q.insert(q.end(), ref.begin() + 22, ref.end());
        const Generated g = gen(p, def, 100, q, 10, 50);
        CHECK(cigar_string(g.cigar.data(), g.n_cigar) == "20M2D18M");
        CHECK(g.md == std::string("20^") + "ACGT"[ref[20]] + "ACGT"[ref[21]] + "18");
        CHECK(g.nm == 2);
        CHECK(g.score == 38 - (6 + 2));
    }
    SUBCASE("three extra query bases: 20M3I20M, MD skips them, NM counts them") {
        std::vector<uint8_t> q(ref.begin(), ref.begin() + 20);
        q.insert(q.end(), {2, 2, 2});
        q.insert(q.end(), ref.begin() + 20, ref.end());
        const Generated g = gen(p, def, 100, q, 10, 50);
        CHECK(cigar_string(g.cigar.data(), g.n_cigar) == "20M3I20M");
        CHECK(g.md == "40");
        CHECK(g.nm == 3);
        CHECK(g.score == 40 - (6 + 3));
    }
    SUBCASE("reverse strand: the CIGAR and MD read forward, the score is the same") {
        // Doubled coordinates [2*60 - 50, 2*60 - 10) = [70, 110) are the reverse
        // complement of forward [10, 50). The query aligned there is the reverse
        // complement of a forward read with one substitution at forward 15 and a
        // deletion of forward bases 30-31.
        std::vector<uint8_t> fq(ref.begin(), ref.begin() + 30);
        fq.insert(fq.end(), ref.begin() + 32, ref.end());
        fq[15] = (uint8_t)((fq[15] + 2) % 4);
        const Generated g = gen(p, def, 100, revcomp(fq), 70, 110);
        CHECK(cigar_string(g.cigar.data(), g.n_cigar) == "30M2D8M");
        CHECK(g.md == std::string("15") + "ACGT"[ref[15]] + "14^" + "ACGT"[ref[30]] + "ACGT"[ref[31]] + "8");
        CHECK(g.nm == 3);
        CHECK(g.score == 37 - 4 - (6 + 2));
    }
}

TEST_CASE("bwa_gen_cigar2: random regions give an optimal CIGAR and the MD and NM the SAM spec derives"
          * doctest::test_suite("unit/cigar")) {
    const Scoring scorings[] = {
        scoring(1, 4, 6, 1, 6, 1),
        scoring(1, 4, 5, 2, 7, 1),
        scoring(2, 8, 12, 2, 12, 2),
    };
    std::mt19937 rng(5150u);
    const Packed p = pack(random_seq(rng, 1200));
    const int64_t L = p.l_pac();
    int checked = 0, gapped = 0, optimal = 0;
    for (const Scoring &s : scorings) {
        for (int n = 0; n < 50; ++n) {
            const int64_t len = 30 + (int64_t)(rng() % 120);
            const bool rev = n % 2 == 1;
            const int64_t fb = (int64_t)(rng() % (uint32_t)(L - len));
            // The region and the query in the doubled coordinates mem passes.
            const int64_t rb = rev ? 2 * L - (fb + len) : fb, re = rb + len;
            const std::vector<uint8_t> fref(p.fwd.begin() + fb, p.fwd.begin() + fb + len);
            const std::vector<uint8_t> fq = mutate(rng, fref);  // the read, forward orientation
            const std::vector<uint8_t> query = rev ? revcomp(fq) : fq;
            REQUIRE(region(p, rb, re) == (rev ? revcomp(fref) : fref));
            CAPTURE(n);
            CAPTURE(rev);
            CAPTURE(len);
            const Generated g = gen(p, s, 100, query, rb, re);
            const std::string cs = cigar_string(g.cigar.data(), g.n_cigar);
            CAPTURE(cs);
            const Derived d = derive(g.cigar.data(), g.n_cigar, fq, fref, s);
            CHECK(d.valid);
            CHECK(g.score == d.score);
            // Where the documented band holds the unbanded optimum, the score
            // must be that optimum.
            const int opt = global_optimum(fq, fref, s);
            if (global_optimum(fq, fref, s, documented_band(s, (int)fq.size(), (int)fref.size(), 100)) == opt) {
                CHECK(g.score == opt);
                ++optimal;
            }
            CHECK(g.md == d.md);
            CHECK(g.nm == d.nm);
            if (g.n_cigar > 1) ++gapped;
            ++checked;
        }
    }
    CHECK(checked == 150);
    CHECK(optimal >= 130);
    CHECK(gapped >= 60);  // the indel paths were exercised, not just ungapped regions
}

TEST_CASE("bwa_gen_cigar2: short reads in long regions keep the band floor of |tlen - qlen| + 3"
          * doctest::test_suite("unit/cigar")) {
    // A 10-24 base read placed in a region 3-12 bases longer (one or two
    // deletions from the read's point of view, sometimes also a 1-2 base
    // insertion and substitutions). For reads this short, the band bwa_gen_cigar4
    // derives from the gap penalties is below |tlen - qlen| + 3, so that floor
    // is the band, and the optimal path may stray up to two diagonals past the
    // end diagonal (an insertion before a deletion). Where the floor band holds
    // the unbanded optimum the score must be that optimum; the test also
    // requires cases where a band one or two narrower would not hold it, so the
    // floor itself is exercised.
    const Scoring def = scoring(1, 4, 6, 1, 6, 1);
    std::mt19937 rng(1919u);
    const Packed p = pack(random_seq(rng, 2000));
    const int64_t L = p.l_pac();
    int optimal = 0, floor_needed = 0;
    for (int n = 0; n < 600; ++n) {
        const int qlen = 10 + (int)(rng() % 15);
        std::vector<uint8_t> fq = random_seq(rng, (size_t)qlen);
        std::vector<uint8_t> fref;
        // The region: the read with extra bases at one or two places.
        const int cut1 = 1 + (int)(rng() % (uint32_t)(qlen - 1)), extra1 = 2 + (int)(rng() % 8);
        const int cut2 = rng() % 2 ? 1 + (int)(rng() % (uint32_t)(qlen - 1)) : -1, extra2 = 1 + (int)(rng() % 4);
        for (int k = 0; k < qlen; ++k) {
            if (k == cut1)
                for (int m = 0; m < extra1; ++m) fref.push_back((uint8_t)(rng() % 4));
            if (k == cut2)
                for (int m = 0; m < extra2; ++m) fref.push_back((uint8_t)(rng() % 4));
            fref.push_back(fq[(size_t)k]);
        }
        // Sometimes a short insertion into the read and a substitution.
        if (rng() % 2) {
            const size_t at = 1 + rng() % (fq.size() - 1);
            for (int m = 1 + (int)(rng() % 2); m > 0; --m) fq.insert(fq.begin() + (long)at, (uint8_t)(rng() % 4));
        }
        if (rng() % 2) fq[rng() % fq.size()] = (uint8_t)(rng() % 4);
        if (fref.size() <= fq.size()) continue;
        // Write the region into the packed reference at a fresh place.
        const int64_t fb = 20 + (int64_t)(rng() % (uint32_t)(L - 100));
        Packed q = p;
        for (size_t k = 0; k < fref.size(); ++k) q.fwd[(size_t)fb + k] = fref[k];
        q = pack(q.fwd);
        const int64_t rb = fb, re = fb + (int64_t)fref.size();
        CAPTURE(n);
        const Generated g = gen(q, def, 100, fq, rb, re);
        const Derived d = derive(g.cigar.data(), g.n_cigar, fq, fref, def);
        CHECK(d.valid);
        CHECK(g.score == d.score);
        CHECK(g.md == d.md);
        CHECK(g.nm == d.nm);
        const int dlen = (int)fref.size() - (int)fq.size();
        const int opt = global_optimum(fq, fref, def);
        if (global_optimum(fq, fref, def, dlen + 3) == opt) {
            CHECK(g.score == opt);
            ++optimal;
            if (global_optimum(fq, fref, def, dlen + 1) < opt) ++floor_needed;
        }
    }
    CHECK(optimal >= 400);
    CHECK(floor_needed >= 5);
}

TEST_CASE("bwa_gen_cigar2: an insertion and a deletion far apart on the diagonal are aligned in the documented band"
          * doctest::test_suite("unit/cigar")) {
    // The read is the region with a k-base insertion and, further on, a k-base
    // deletion (k = 8-20), so read and region have the same length and the
    // floor band (3) cannot hold the alignment: the band derived from the gap
    // penalties (half the payable gap; 23 for a 100-base read at the defaults)
    // must. Both strands; the score must be the unbanded optimum, which the
    // test also requires to need more than the floor band.
    const Scoring def = scoring(1, 4, 6, 1, 6, 1);
    std::mt19937 rng(2323u);
    const Packed p = pack(random_seq(rng, 3000));
    const int64_t L = p.l_pac();
    int beyond_floor = 0;
    for (int n = 0; n < 40; ++n) {
        const int len = 80 + (int)(rng() % 70), k = 8 + (int)(rng() % 13);
        const int64_t fb = (int64_t)(rng() % (uint32_t)(L - len));
        const std::vector<uint8_t> fref(p.fwd.begin() + fb, p.fwd.begin() + fb + len);
        // insertion after base i0, deletion of [d0, d0 + k), d0 >= i0 + 10
        const int i0 = 10 + (int)(rng() % (uint32_t)(len / 3)), d0 = i0 + 10 + (int)(rng() % 15);
        if (d0 + k + 10 > len) continue;
        std::vector<uint8_t> fq(fref.begin(), fref.begin() + i0);
        const std::vector<uint8_t> ins = random_seq(rng, (size_t)k);
        fq.insert(fq.end(), ins.begin(), ins.end());
        fq.insert(fq.end(), fref.begin() + i0, fref.begin() + d0);
        fq.insert(fq.end(), fref.begin() + d0 + k, fref.end());
        REQUIRE(fq.size() == fref.size());
        const bool rev = n % 2 == 1;
        const int64_t rb = rev ? 2 * L - (fb + len) : fb, re = rb + len;
        CAPTURE(n);
        CAPTURE(k);
        CAPTURE(rev);
        const Generated g = gen(p, def, 100, rev ? revcomp(fq) : fq, rb, re);
        const Derived d = derive(g.cigar.data(), g.n_cigar, fq, fref, def);
        CHECK(d.valid);
        CHECK(g.score == d.score);
        CHECK(g.md == d.md);
        CHECK(g.nm == d.nm);
        const int opt = global_optimum(fq, fref, def);
        const int band = documented_band(def, len, len, 100);
        REQUIRE(global_optimum(fq, fref, def, band) == opt);
        CHECK(g.score == opt);
        if (global_optimum(fq, fref, def, 3) < opt) ++beyond_floor;
    }
    CHECK(beyond_floor >= 25);
}

TEST_CASE("bwa_gen_cigar2: an indel in a repeat is placed leftmost on the forward reference, both strands"
          * doctest::test_suite("unit/cigar")) {
    const Scoring def = scoring(1, 4, 6, 1, 6, 1);
    // Forward reference: unique flank, a run of 8 A, unique flank, (AC)x5,
    // unique flank. Each read drops or adds one A, or one AC, inside a repeat,
    // so the indel can sit anywhere in it at equal score; leftmost means right
    // after the flank before the repeat.
    const std::vector<uint8_t> left = digits("1230132102312"), mid = digits("2310213"), right = digits("3123021321032");
    std::vector<uint8_t> fwd = left;
    const size_t homo = fwd.size();
    fwd.insert(fwd.end(), 8, 0);
    fwd.insert(fwd.end(), mid.begin(), mid.end());
    const size_t dinuc = fwd.size();
    for (int k = 0; k < 5; ++k) fwd.insert(fwd.end(), {0, 1});
    fwd.insert(fwd.end(), right.begin(), right.end());
    const Packed p = pack(fwd);
    const int64_t L = p.l_pac(), len = (int64_t)fwd.size();
    struct Edit {
        const char *what;
        std::vector<uint8_t> fq;
        std::string cigar;
    };
    std::vector<Edit> edits;
    {
        std::vector<uint8_t> q = fwd;
        q.erase(q.begin() + (long)homo + 5);  // any A of the run
        edits.push_back({"one A deleted from the run", q,
                         std::to_string(homo) + "M1D" + std::to_string(len - homo - 1) + "M"});
        q = fwd;
        q.insert(q.begin() + (long)homo + 3, 0);
        edits.push_back({"one A inserted into the run", q,
                         std::to_string(homo) + "M1I" + std::to_string(len - homo) + "M"});
        q = fwd;
        q.erase(q.begin() + (long)dinuc + 4, q.begin() + (long)dinuc + 6);  // one AC
        edits.push_back({"one AC deleted from the repeat", q,
                         std::to_string(dinuc) + "M2D" + std::to_string(len - dinuc - 2) + "M"});
    }
    for (const Edit &e : edits) {
        CAPTURE(e.what);
        for (bool rev : {false, true}) {
            CAPTURE(rev);
            const int64_t rb = rev ? 2 * L - len : 0, re = rb + len;
            const Generated g = gen(p, def, 100, rev ? revcomp(e.fq) : e.fq, rb, re);
            CHECK(cigar_string(g.cigar.data(), g.n_cigar) == e.cigar);
            const Derived d = derive(g.cigar.data(), g.n_cigar, e.fq, fwd, def);
            CHECK(g.md == d.md);
            CHECK(g.nm == d.nm);
        }
    }
}
