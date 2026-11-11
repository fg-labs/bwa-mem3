// A converted (`index --meth`) test index for the --meth refutation tests (aln_meth_prove.h), and
// the index text T as the index itself holds it, read back through its suffix array: T[SA[r]] is the
// first base of row r. The tests compare the refutations with exhaustive computations on that
// text, so they do not rely on how the text masks are derived from the original reference.
//
// Contigs: random bases; random bases around an N run; a soft-masked (lower-case) stretch, which the
// converted FASTA does not convert; and a 120-base unit repeated with a few substitutions per copy,
// so that intervals hold several rows.
#ifndef BWA_MEM3_TEST_ALN_METH_PROVE_FIXTURE_H
#define BWA_MEM3_TEST_ALN_METH_PROVE_FIXTURE_H

#include "doctest/doctest.h"

#include "aln_test_ref.h"
#include "../framework/temp_dir.h"
#include "../../src/aln/aln_bwt.h"
#include "../../src/aln/aln_meth.h"
#include "../../src/aln/aln_meth_prove.h"
#include "../../src/FMI_search.h"
#include "../../src/main.h"

#include <getopt.h>
#include <algorithm>
#include <cctype>
#include <memory>
#include <string>
#include <vector>

namespace aln_prove_test {

namespace ba = bwa_aln;

struct ProveIndex {
    bwa_tests::TempDir tmp{"aln_meth_prove"};
    std::string fa;
    std::unique_ptr<FMI_search> fmi;
    ba::bwt_t *bwt = nullptr;
    ba::aln_meth_ref_t ref;
    std::vector<uint8_t> text; // T as indexed (0-3), length bwt->seq_len
    std::string repeat_unit;

    ProveIndex()
    {
        fa = tmp.path() + "/ref.fa";
        uint32_t st = 20261111;
        std::string c1 = aln_test::lcg_dna(st, 4000);
        // Sequence the stateful draws explicitly so every compiler builds the same index.
        std::string c2_right = aln_test::lcg_dna(st, 1500);
        std::string c2_left = aln_test::lcg_dna(st, 1500);
        std::string c2 = c2_left + std::string(30, 'N') + c2_right;
        std::string c3 = aln_test::lcg_dna(st, 1000);
        for (int i = 200; i < 700; ++i) c3[i] = (char)std::tolower(c3[i]);
        repeat_unit = aln_test::lcg_dna(st, 120);
        std::string c4 = aln_test::lcg_dna(st, 50);
        for (int k = 0; k < 20; ++k) {
            std::string u = repeat_unit;
            for (int s = 0; s < 3; ++s) { // three substitutions per copy (some copies share them)
                const int at = (int)((k * 37 + s * 41) % 120);
                u[at] = u[at] == 'A' ? 'G' : 'A';
            }
            c4 += u + aln_test::lcg_dna(st, 30);
        }
        REQUIRE(aln_test::write_test_fasta(fa, {{"c1", c1}, {"c2", c2}, {"c3", c3}, {"c4", c4}}));
        std::vector<std::string> a = {"index", "--meth", fa};
        std::vector<char *> argv;
        for (std::string &s : a) argv.push_back(&s[0]);
        argv.push_back(nullptr);
        optind = 1;
#if defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || defined(__NetBSD__)
        optreset = 1;
#endif
        REQUIRE(bwa_index((int)a.size(), argv.data()) == 0);
        fmi.reset(new FMI_search((fa + ".meth").c_str()));
        fmi->load_index(false, 1);
        bwt = ba::bwt_from_fmi(fmi.get());
        REQUIRE(ba::aln_meth_ref_load(&ref, fa.c_str(), fmi->idx->bns, 1) == 0);
        // T from the index: row r (1-based in bwa's numbering; row 0 is `$`) starts with base c
        // for L2[c] < r <= L2[c + 1].
        text.assign(bwt->seq_len, 0xff);
        for (ba::bwtint_t r = 1; r <= bwt->seq_len; ++r) {
            int c = 0;
            while (c < 3 && r > bwt->L2[c + 1]) ++c;
            const ba::bwtint_t y = ba::bwt_sa(bwt, r);
            REQUIRE(y < bwt->seq_len);
            text[y] = (uint8_t)c;
        }
        for (uint8_t b : text) REQUIRE(b < 4);
    }
    ~ProveIndex()
    {
        ba::aln_meth_ref_free(&ref);
        ba::bwt_destroy(bwt);
    }
    // The SA interval of the exact string s (codes 0-3) in T, by backward search; false when empty.
    bool interval(const std::vector<uint8_t> &s, ba::bwtint_t *k, ba::bwtint_t *l) const
    {
        ba::bwtint_t kk = 0, ll = bwt->seq_len;
        for (int i = (int)s.size() - 1; i >= 0; --i) {
            const uint8_t c = s[i];
            if (c > 3) return false;
            kk = bwt->L2[c] + ba::bwt_occ(bwt, kk - 1, c) + 1;
            ll = bwt->L2[c] + ba::bwt_occ(bwt, ll, c);
            if (kk > ll) return false;
        }
        *k = kk; *l = ll;
        return true;
    }
};

inline ProveIndex &prove_index()
{
    static ProveIndex idx;
    return idx;
}

// Unit edit distance of q against t (both codes; q's 4 matches nothing) by full dynamic programming:
// anchored at both ends of q; t's start anchored when anchor_start, its end free when free_end.
inline int full_ed(const std::vector<uint8_t> &q, const std::vector<uint8_t> &t, bool anchor_start, bool free_end)
{
    const size_t n = q.size(), m = t.size();
    std::vector<int> prev(m + 1), cur(m + 1);
    for (size_t b = 0; b <= m; ++b) prev[b] = anchor_start ? (int)b : 0;
    for (size_t a = 1; a <= n; ++a) {
        cur[0] = (int)a;
        for (size_t b = 1; b <= m; ++b)
            cur[b] = std::min({prev[b - 1] + (q[a - 1] < 4 && q[a - 1] == t[b - 1] ? 0 : 1), prev[b] + 1, cur[b - 1] + 1});
        std::swap(prev, cur);
    }
    return free_end ? *std::min_element(prev.begin(), prev.end()) : prev[m];
}

// A small deterministic generator for the tests' random choices.
struct Rng {
    uint64_t s;
    explicit Rng(uint64_t seed) : s(seed * 0x9e3779b97f4a7c15ULL + 1) {}
    uint32_t next() { s ^= s << 13; s ^= s >> 7; s ^= s << 17; return (uint32_t)(s >> 11); }
    int below(int n) { return (int)(next() % (uint32_t)n); }
};

} // namespace aln_prove_test

#endif
