// A converted (`index --meth`) test index for the --meth refutation tests (aln_meth_prove.h), or a
// plain one (`index`, aln_meth_ref_t::plain: plain aln's no-hit proof), and
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

    bool plain;

    explicit ProveIndex(bool plain_index = false) : plain(plain_index)
    {
        fa = tmp.unique_file(".fa"); // a unique basename: aln attaches a staged index by name
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
        std::vector<std::string> a = plain ? std::vector<std::string>{"index", fa}
                                           : std::vector<std::string>{"index", "--meth", fa};
        std::vector<char *> argv;
        for (std::string &s : a) argv.push_back(&s[0]);
        argv.push_back(nullptr);
        optind = 1;
#if defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || defined(__NetBSD__)
        optreset = 1;
#endif
        REQUIRE(bwa_index((int)a.size(), argv.data()) == 0);
        fmi.reset(new FMI_search((plain ? fa : fa + ".meth").c_str()));
        fmi->load_index(plain, 1);
        bwt = ba::bwt_from_fmi(fmi.get());
        if (plain) { // the index's own bns and pac, as plain aln sets them (plain_prove_scope_t)
            ref.conv_bns = fmi->idx->bns; ref.bns = fmi->idx->bns; ref.pac = fmi->idx->pac; ref.plain = 1;
        } else {
            REQUIRE(ba::aln_meth_ref_load(&ref, fa.c_str(), fmi->idx->bns, 1) == 0);
        }
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
        // Everything is loaded (the last file read is aln_meth_ref_load): the shared fixture keeps
        // no files on disk for the rest of the run.
        tmp.remove();
    }
    ~ProveIndex()
    {
        if (!plain) ba::aln_meth_ref_free(&ref);
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
inline ProveIndex &prove_index_plain()
{
    static ProveIndex idx(true);
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

// Myers' bit-parallel approximate search (one 128-bit word per text position) for a query of at
// most 128 bases: the fewest edits of q against any substring of a text of n positions, where
// eq(peq, b) gives the query positions that match text position b (peq[c]: those holding base c;
// a query code past 3 matches nothing). Lets a test afford a whole-text oracle.
template <class Eq> inline int myers_search(const std::vector<uint8_t> &q, size_t n, Eq eq_at)
{
    using word = unsigned __int128;
    const size_t m = q.size();
    REQUIRE(m <= 128);
    if (m == 0) return 0;
    word peq[4] = {0, 0, 0, 0};
    for (size_t i = 0; i < m; ++i)
        if (q[i] < 4) peq[q[i]] |= (word)1 << i;
    const word high = (word)1 << (m - 1);
    word pv = ~(word)0, mv = 0;
    int score = (int)m, best = score;
    for (size_t b = 0; b < n; ++b) {
        const word eq = eq_at(peq, b);
        const word xv = eq | mv, xh = (((eq & pv) + pv) ^ pv) | eq;
        word ph = mv | ~(xh | pv), mh = pv & xh;
        if (ph & high) ++score;
        if (mh & high) --score;
        ph <<= 1; mh <<= 1; // the first row is free (no carry in): a match may start anywhere
        pv = mh | ~(xv | ph); mv = ph & xv;
        best = std::min(best, score);
    }
    return best;
}

// full_ed(q, t, false, true), bit-parallel.
inline int search_ed(const std::vector<uint8_t> &q, const std::vector<uint8_t> &t)
{
    return myers_search(q, t.size(), [&](const unsigned __int128 *peq, size_t b) {
        return t[b] < 4 ? peq[t[b]] : (unsigned __int128)0;
    });
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
