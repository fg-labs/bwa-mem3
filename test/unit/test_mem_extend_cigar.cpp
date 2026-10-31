// test/unit/test_mem_extend_cigar.cpp
//
// CIGAR generation (bwa_gen_cigar2, src/bwa.cpp), the step that turns an
// alignment region [rb, re) x query into the CIGAR of every SAM record `mem`
// writes, on the two cases nothing else pins (MD, NM and the score against a
// naive walk of the CIGAR are covered by test_gen_cigar_no_md.cpp, and the
// CIGAR against ksw_global2 at its band by test_global_cert_band.cpp):
//
//   - a read with an insertion and, further on, a deletion of 8-20 bases has
//     the same length as its region, so the band floor of |tlen - qlen| + 3
//     cannot hold the alignment and the band derived from the gap penalties
//     (half the gap the read could pay for, at most the caller's band; bwa's
//     bwa_gen_cigar2) must: the score is the unbanded optimum of a full-matrix
//     three-state DP written here;
//   - indels are placed at the leftmost position on the forward reference on
//     both strands (bwa_gen_cigar4: reversing both sequences on the reverse
//     strand "is to ensure indels to be placed at the leftmost position"),
//     checked on homopolymer and dinucleotide-repeat indels whose placement is
//     otherwise free. On the reverse strand (rb >= l_pac) the region is the
//     reverse complement of the forward reference, and the CIGAR describes the
//     reverse complement of the query against the forward reference.
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
#include "scoring.h"

namespace {

struct Scoring {
    int a, b, o_del, e_del, o_ins, e_ins;
    int8_t mat[25];
};

Scoring scoring(int a, int b, int o_del, int e_del, int o_ins, int e_ins) {
    Scoring s{a, b, o_del, e_del, o_ins, e_ins, {}};
    const bwa_tests::ScoringMatrix m = bwa_tests::build_scoring_matrix(a, b, 1);
    std::copy(m.begin(), m.end(), s.mat);
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
    for (size_t k = 0; k < fwd.size(); ++k) _set_pac(p.pac, k, fwd[k]);
    return p;
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
// of rlen under band cap w_ (bwa 0.7.x bwa.c): half of the sum of the longest
// insertion or deletion that half the read's matches could pay for and
// |rlen - qlen|, at most w_ and at least |rlen - qlen| + 3. A narrower band computed by the code
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

std::vector<uint8_t> digits(const char *str) {
    std::vector<uint8_t> v;
    for (; *str; ++str) v.push_back((uint8_t)(*str - '0'));
    return v;
}

}  // namespace

TEST_CASE("bwa_gen_cigar2: an insertion and a deletion far apart on the diagonal are aligned in the documented band"
          * doctest::test_suite("unit/cigar")) {
    // The read is the region with a k-base insertion and, further on, a k-base
    // deletion (k = 8-20), so read and region have the same length and the
    // floor band (3) cannot hold the alignment: the band derived from the gap
    // penalties (half the payable gap; 23 for a 100-base read at the defaults)
    // must. Both strands; the score must be the unbanded optimum, which the
    // test also requires to need more than the floor band. A draw whose
    // documented band would not hold the optimum (a read too short for k) is
    // skipped and counted.
    const Scoring def = scoring(1, 4, 6, 1, 6, 1);
    std::mt19937 rng(2323u);
    const Packed p = pack(random_seq(rng, 3000));
    const int64_t L = p.l_pac();
    int beyond_floor = 0, in_band = 0;
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
        const int opt = global_optimum(fq, fref, def);
        const int band = documented_band(def, len, len, 100);
        if (global_optimum(fq, fref, def, band) != opt) continue;
        ++in_band;
        CHECK(g.score == opt);
        if (global_optimum(fq, fref, def, 3) < opt) ++beyond_floor;
    }
    CHECK(in_band >= 25);
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
        }
    }
}
