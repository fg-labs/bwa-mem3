// test/unit/test_gen_cigar_no_md.cpp — the MD-free CIGAR regen used for XA hits.
//
// mem_gen_alt formats each XA hit from mem_reg2aln's CIGAR and NM and never
// reads the MD string mem_reg2aln builds after the CIGAR, so it asks for none
// (want_md 0). That is only a performance change if everything the XA entry
// reads is what the MD-building form gives. These tests run both forms on the
// same inputs and require the same score, CIGAR words and NM, from
// bwa_gen_cigar4 (both strands, the ungapped path and the global DP, leading
// and trailing deletions, query Ns, and the --meth NM rule that skips columns
// the matrix does not penalise) and from mem_reg2aln (the deletion squeeze and
// the clipping moves, which carry no MD bytes in the MD-free form while the MD
// form still carries its MD through them, with and without the --meth regen).
// Both forms share one NM rule, so hand-built alignments pin the NM values
// themselves, and mem_gen_alt's XA string is checked against one written out
// by hand. Both forms' score and NM, and the MD form's MD, are also checked
// against a naive scalar walk of the returned CIGAR (gen_cigar_oracle), so an
// error the two forms share (the vector M-op walks, or the reverse-strand
// flips) cannot pass by agreeing with itself.

#include "doctest/doctest.h"

#include "bwa.h"
#include "bwamem.h"
#include "meth_orig_ref_fixture.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <vector>

namespace {

constexpr int kRefLen = 3000;

// A deterministic, low-repeat reference: a 32-bit LCG over ACGT.
std::string make_reference(int len) {
    std::string ref(len, 'A');
    uint32_t x = 2024;
    for (int i = 0; i < len; ++i) {
        x = x * 1103515245u + 12345u;
        ref[i] = "ACGT"[(x >> 16) & 3];
    }
    return ref;
}

uint8_t nt4(char c) { return nst_nt4_table[(unsigned char)c]; }

// The default symmetric matrix (-A 1 -B 4, N scores -1), and a --meth-style
// asymmetric copy whose ref-C/query-T cell is free (`freed_score` >= 0), so
// that column is a literal difference that NM under nm_from_mat 1 skips.
struct Matrices {
    int8_t sym[25];
    int8_t meth[25];
    explicit Matrices(int freed_score) {
        bwa_fill_scmat(1, 4, sym);
        memcpy(meth, sym, sizeof(meth));
        meth[1 * 5 + 3] = (int8_t)freed_score;
    }
};

struct GenCigarOut {
    uint32_t *cigar = nullptr;
    int score = 0, n_cigar = 0, nm = -1;
    std::vector<uint8_t> query_after;  // the query buffer after the call
};

GenCigarOut run_gen_cigar(const meth_test::OrigRefFixture &fx, const int8_t *mat,
                          const std::vector<uint8_t> &query, int64_t rb, int64_t re, int band,
                          int nm_from_mat, int want_md) {
    GenCigarOut out;
    out.query_after = query;
    out.cigar = bwa_gen_cigar4(mat, 5, NULL, 6, 1, 6, 1, band, fx.bns.l_pac, fx.orig_pac,
                               (int)query.size(), out.query_after.data(), rb, re, &out.score,
                               &out.n_cigar, &out.nm, nm_from_mat, want_md);
    return out;
}

std::vector<uint32_t> words(const uint32_t *cigar, int n) {
    return cigar ? std::vector<uint32_t>(cigar, cigar + n) : std::vector<uint32_t>();
}

// The forward-strand bases ref[at, at + len) with random edits, in the read's
// own orientation (reverse-complemented when `rev`), and the doubled-strand
// window [rb, re) to align it to: the edited span plus `pad_l`/`pad_r` extra
// reference bases, which the global alignment closes with deletions at the
// ends.
struct Case {
    std::vector<uint8_t> query;
    int64_t rb = 0, re = 0;
    int n_edits = 0;
};

Case make_case(const std::string &ref, std::mt19937 &rng, bool rev, int pad_l, int pad_r,
               int len = 0) {
    if (len <= 0) len = 40 + (int)(rng() % 160);
    const int at = 50 + (int)(rng() % (kRefLen - len - 100));
    std::string fwd = ref.substr(at, len);
    std::string edited;
    Case c;
    for (int i = 0; i < len; ++i) {
        const uint32_t r = rng() % 100;
        if (i > 3 && i < len - 4 && r < 3) {  // substitution, often a C->T
            edited += (fwd[i] == 'C' && (rng() & 1)) ? 'T' : "ACGT"[(nt4(fwd[i]) + 1 + rng() % 3) & 3];
            ++c.n_edits;
        } else if (i > 3 && i < len - 4 && r < 4) {  // query N
            edited += 'N';
            ++c.n_edits;
        } else if (i > 3 && i < len - 4 && r < 5) {  // deletion of this base
            ++c.n_edits;
        } else if (i > 3 && i < len - 4 && r < 6) {  // insertion before this base
            edited += "ACGT"[rng() & 3];
            edited += fwd[i];
            ++c.n_edits;
        } else {
            edited += fwd[i];
        }
    }
    if (rev) {
        std::string rc(edited.rbegin(), edited.rend());
        for (char &ch : rc) ch = ch == 'A' ? 'T' : ch == 'C' ? 'G' : ch == 'G' ? 'C' : ch == 'T' ? 'A' : 'N';
        edited = rc;
    }
    for (char ch : edited) c.query.push_back(nt4(ch));
    const int64_t fb = at - pad_l, fe = at + len + pad_r;
    c.rb = rev ? 2 * (int64_t)kRefLen - fe : fb;
    c.re = rev ? 2 * (int64_t)kRefLen - fb : fe;
    return c;
}

// Op counts over a CIGAR, to check the random cases reach every NM rule.
struct OpSeen {
    int ungapped = 0, internal_del = 0, end_del = 0, ins = 0, freed = 0;
};

void tally(const std::vector<uint32_t> &cig, OpSeen &seen) {
    const int n = (int)cig.size();
    if (n == 1 && (cig[0] & 0xf) == 0) ++seen.ungapped;
    for (int k = 0; k < n; ++k) {
        const int op = cig[k] & 0xf;
        if (op == 1) ++seen.ins;
        if (op == 2) ++((k > 0 && k < n - 1) ? seen.internal_del : seen.end_del);
    }
}

// The reference bases an MD string spans (match runs, mismatched bases and
// deleted bases), or -1 if it is not well formed.
int md_ref_span(const char *md) {
    int span = 0;
    for (const char *p = md; *p;) {
        if (*p >= '0' && *p <= '9') {
            span += (int)strtol(p, (char **)&p, 10);
        } else if (*p == '^') {
            for (++p; *p >= 'A' && *p <= 'Z'; ++p) ++span;
        } else if (*p >= 'A' && *p <= 'Z') {
            ++span, ++p;
        } else {
            return -1;
        }
    }
    return span;
}

// The reference bases MD lists for a CIGAR: its M ops and its internal D ops.
// MD leaves out a deletion at either end of the aligned ops, and mem_reg2aln
// squeezes out only one of them, so the other can remain next to a clip.
int cigar_ref_span(const uint32_t *cigar, int n) {
    int first = 0, last = n - 1, span = 0;
    while (first < n && (cigar[first] & 0xf) == 3) ++first;
    while (last >= 0 && (cigar[last] & 0xf) == 3) --last;
    for (int k = first; k <= last; ++k) {
        const int op = cigar[k] & 0xf;
        if (op == 0 || (op == 2 && k > first && k < last)) span += (int)(cigar[k] >> 4);
    }
    return span;
}

// What bwa_gen_cigar4 must return for a CIGAR, computed the slow way. The
// window [rb, re) and the query are put in the frame the walk runs in: on the
// reverse strand (rb >= l_pac) the window is the complement of forward bases
// [2 l_pac - re, 2 l_pac - rb) read in forward order, and the query is reversed;
// MD then names a reference base by its forward-strand letter. The score is the
// global alignment score of the CIGAR's path (mat over M columns, o + e * len
// per gap). NM counts every literal mismatch (or, under nm_from_mat, only those
// mat penalises), every inserted base and the deleted bases of a D that is not
// the first or the last op; MD lists the M ops' literal mismatches and those
// internal deletions.
struct Oracle {
    int score = 0, nm = 0;
    int freed = 0;  // literal mismatches NM skipped under nm_from_mat
    std::string md;
};

Oracle gen_cigar_oracle(const std::string &ref, const int8_t *mat, int o_del, int e_del,
                        int o_ins, int e_ins, const std::vector<uint8_t> &query, int64_t rb,
                        int64_t re, const uint32_t *cigar, int n_cigar, int nm_from_mat) {
    const int64_t l_pac = (int64_t)ref.size();
    const bool rev = rb >= l_pac;
    const int lq = (int)query.size(), rlen = (int)(re - rb);
    std::vector<uint8_t> q(query), r(rlen);
    if (rev) std::reverse(q.begin(), q.end());
    for (int i = 0; i < rlen; ++i)
        r[i] = rev ? (uint8_t)(3 - nt4(ref[2 * l_pac - re + i])) : nt4(ref[rb + i]);
    const char *int2base = rev ? "TGCAN" : "ACGTN";
    Oracle o;
    int x = 0, y = 0, u = 0;
    for (int k = 0; k < n_cigar; ++k) {
        const int op = cigar[k] & 0xf, len = (int)(cigar[k] >> 4);
        if (op == 0) {
            for (int i = 0; i < len; ++i, ++x, ++y) {
                o.score += mat[r[y] * 5 + q[x]];
                if (q[x] != r[y]) {
                    o.md += std::to_string(u);
                    o.md += int2base[r[y]];
                    u = 0;
                    if (!nm_from_mat || mat[r[y] * 5 + q[x]] < 0) ++o.nm;
                    else ++o.freed;
                } else {
                    ++u;
                }
            }
        } else if (op == 1) {
            o.score -= o_ins + e_ins * len;
            o.nm += len;
            x += len;
        } else if (op == 2) {
            o.score -= o_del + e_del * len;
            if (k > 0 && k < n_cigar - 1) {
                o.md += std::to_string(u) + "^";
                for (int i = 0; i < len; ++i) o.md += int2base[r[y + i]];
                u = 0;
                o.nm += len;
            }
            y += len;
        }
    }
    o.md += std::to_string(u);
    CHECK(x == lq);
    CHECK(y == rlen);
    return o;
}

}  // namespace

TEST_CASE("bwa_gen_cigar4: score, NM and MD match a naive walk of the CIGAR" *
          doctest::test_suite("unit/gen_cigar")) {
#if defined(__AVX2__)
    MESSAGE(std::string("M-op walks: AVX2 32-column blocks"));
#else
    MESSAGE(std::string("M-op walks: scalar loops (no AVX2 on this target)"));
#endif
    const std::string ref = make_reference(kRefLen);
    meth_test::OrigRefFixture fx(ref);
    std::mt19937 rng(23);
    int n_rev = 0, n_gapped = 0, n_n = 0, n_short = 0, n_freed_zero = 0;
    // The freed --meth cell scores 1, then 0: NM under nm_from_mat skips a
    // column whose cell is not negative, so 0 is the edge of that test.
    for (int freed_score : {1, 0}) {
        const Matrices m(freed_score);
        // Every query length from 1 to 200, so the 32-column blocks, the
        // remainder loads at either end and the scalar fallback of short ops
        // all run.
        for (int len = 1; len <= 200; ++len) {
            for (int rep = 0; rep < 6; ++rep) {
                const bool rev = rep & 1;
                const int pad_l = (rep == 2) ? 1 + (int)(rng() % 3) : 0;
                const int pad_r = (rep == 3) ? 1 + (int)(rng() % 3) : 0;
                Case c = make_case(ref, rng, rev, pad_l, pad_r, len);
                // substitutions and Ns anywhere, including the first and last bases
                for (uint8_t &b : c.query) {
                    const uint32_t p = rng() % 100;
                    if (p < 4) b = (uint8_t)((b + 1 + rng() % 3) & 3);
                    else if (p < 6) b = 4;
                }
                static const int kBands[] = {0, 3, 10, 50, 100};
                const int band = kBands[rng() % 5];
                for (int combo = 0; combo < 3; ++combo) {
                    const int nm_from_mat = combo == 2;
                    const int8_t *mat = combo == 0 ? m.sym : m.meth;
                    CAPTURE(freed_score);
                    CAPTURE(len);
                    CAPTURE(rep);
                    CAPTURE(combo);
                    CAPTURE(band);
                    GenCigarOut md = run_gen_cigar(fx, mat, c.query, c.rb, c.re, band, nm_from_mat, 1);
                    REQUIRE(md.cigar != nullptr);
                    const Oracle o = gen_cigar_oracle(ref, mat, 6, 1, 6, 1, c.query, c.rb, c.re,
                                                      md.cigar, md.n_cigar, nm_from_mat);
                    CHECK(md.score == o.score);
                    CHECK(md.nm == o.nm);
                    CHECK(std::string((const char *)(md.cigar + md.n_cigar)) == o.md);
                    CHECK(md.query_after == c.query);
                    // the MD-free form counts NM with its own walks (a lane
                    // mask under AVX2), so it is held to the oracle as well
                    GenCigarOut no = run_gen_cigar(fx, mat, c.query, c.rb, c.re, band, nm_from_mat, 0);
                    REQUIRE(no.cigar != nullptr);
                    CHECK(no.score == o.score);
                    CHECK(no.nm == o.nm);
                    CHECK(words(no.cigar, no.n_cigar) == words(md.cigar, md.n_cigar));
                    CHECK(no.query_after == c.query);
                    n_rev += rev;
                    n_gapped += md.n_cigar > 1;
                    n_short += c.query.size() < 32;
                    if (freed_score == 0) n_freed_zero += o.freed > 0;
                    for (uint8_t b : c.query) n_n += b == 4;
                    free(md.cigar);
                    free(no.cigar);
                }
            }
        }
    }
    // the cases reached both strands, gapped CIGARs, query Ns, short queries and
    // freed columns whose cell scores exactly 0
    CHECK(n_rev > 1000);
    CHECK(n_gapped > 500);
    CHECK(n_n > 1000);
    CHECK(n_short > 300);
    CHECK(n_freed_zero > 20);
}

// Test-only hook in bwa.cpp: fills this thread's bwa_gen_cigar4 reference scratch
// (at least min_size bytes) with `byte`, so the next call's headroom past rlen holds it.
extern "C" void bwa_gen_cigar_fill_rseq_scratch(size_t min_size, int byte);

namespace {

// A query that aligns to a short reference window with long insertions: the
// edited forward bases ref[at, at + L), with `tail_l` and `tail_r` random bases
// before and after, in the read's own orientation, and the doubled-strand
// window [rb, re) of the L bases. The window is shorter than a 32-column block
// and the query longer, so an M op's forward remainder load reads the reference
// headroom past rlen.
Case make_short_window_case(const std::string &ref, std::mt19937 &rng, bool rev, int L,
                            int tail_l, int tail_r) {
    const int at = 50 + (int)(rng() % (kRefLen - 200));
    std::string edited;
    for (int i = 0; i < tail_l; ++i) edited += "ACGT"[rng() & 3];
    for (int i = 0; i < L; ++i) {
        const uint32_t p = rng() % 100;
        edited += p < 8 ? "ACGT"[(nt4(ref[at + i]) + 1 + rng() % 3) & 3] : p < 10 ? 'N' : ref[at + i];
    }
    for (int i = 0; i < tail_r; ++i) edited += "ACGT"[rng() & 3];
    if (rev) {
        std::string rc(edited.rbegin(), edited.rend());
        for (char &ch : rc) ch = ch == 'A' ? 'T' : ch == 'C' ? 'G' : ch == 'G' ? 'C' : ch == 'T' ? 'A' : 'N';
        edited = rc;
    }
    Case c;
    for (char ch : edited) c.query.push_back(nt4(ch));
    c.rb = rev ? 2 * (int64_t)kRefLen - (at + L) : at;
    c.re = rev ? 2 * (int64_t)kRefLen - at : at + L;
    return c;
}

// The M ops of a CIGAR (in the walk's frame) whose vector remainder is the
// forward load and reads past rlen: cig_blocks' condition, written out.
int forward_loads_past_rlen(const uint32_t *cigar, int n_cigar, int l_query, int rlen) {
    int x = 0, y = 0, n = 0;
    for (int k = 0; k < n_cigar; ++k) {
        const int op = cigar[k] & 0xf, len = (int)(cigar[k] >> 4);
        if (op == 0) {
            const int rem = len % 32, i = len - rem;
            if (rem && !(x + len >= 32 && y + len >= 32) && l_query - x - i >= 32 && y + i + 32 > rlen) ++n;
            x += len, y += len;
        } else if (op == 1) {
            x += len;
        } else if (op == 2) {
            y += len;
        }
    }
    return n;
}

}  // namespace

TEST_CASE("bwa_gen_cigar4: bytes with the high bits set in the reference headroom change nothing" *
          doctest::test_suite("unit/gen_cigar")) {
    // The vector walks' forward remainder load reads up to 31 reference bytes
    // past rlen, which no call writes, and relies on the lane masks alone to
    // keep them out of the score, NM and MD (cig_lut_get's r << 2 carries a high
    // byte's top bits into the next lane, which is masked as well). Fill that
    // headroom with such bytes before every call and hold both forms to the
    // naive walk.
    const std::string ref = make_reference(kRefLen);
    meth_test::OrigRefFixture fx(ref);
    std::mt19937 rng(31);
    int n_past = 0;
    for (int fill : {0xff, 0x80, 0x7f}) {
        const Matrices m(0);
        for (int iter = 0; iter < 600; ++iter) {
            const bool rev = iter & 1;
            const int L = 1 + (int)(rng() % 31);
            const int tail = 32 + (int)(rng() % 17);
            const int split = (iter / 2) % 3;  // the insertions after, before or around the window
            const int tail_l = split == 0 ? 0 : split == 1 ? tail : tail / 2;
            const Case c = make_short_window_case(ref, rng, rev, L, tail_l, tail - tail_l);
            const size_t need = (size_t)(c.re - c.rb) + 64;
            for (int combo = 0; combo < 3; ++combo) {
                const int nm_from_mat = combo == 2;
                const int8_t *mat = combo == 0 ? m.sym : m.meth;
                CAPTURE(fill);
                CAPTURE(iter);
                CAPTURE(combo);
                bwa_gen_cigar_fill_rseq_scratch(need, fill);
                GenCigarOut md = run_gen_cigar(fx, mat, c.query, c.rb, c.re, 100, nm_from_mat, 1);
                bwa_gen_cigar_fill_rseq_scratch(need, fill);
                GenCigarOut no = run_gen_cigar(fx, mat, c.query, c.rb, c.re, 100, nm_from_mat, 0);
                REQUIRE(md.cigar != nullptr);
                REQUIRE(no.cigar != nullptr);
                const Oracle o = gen_cigar_oracle(ref, mat, 6, 1, 6, 1, c.query, c.rb, c.re,
                                                  md.cigar, md.n_cigar, nm_from_mat);
                CHECK(md.score == o.score);
                CHECK(md.nm == o.nm);
                CHECK(std::string((const char *)(md.cigar + md.n_cigar)) == o.md);
                CHECK(no.score == o.score);
                CHECK(no.nm == o.nm);
                CHECK(words(no.cigar, no.n_cigar) == words(md.cigar, md.n_cigar));
                CHECK(md.query_after == c.query);
                CHECK(no.query_after == c.query);
                n_past += forward_loads_past_rlen(md.cigar, md.n_cigar, (int)c.query.size(),
                                                  (int)(c.re - c.rb));
                free(md.cigar);
                free(no.cigar);
            }
        }
    }
    bwa_gen_cigar_fill_rseq_scratch(0, 0);
    // the forward remainder load read past rlen in many of the calls
    CHECK(n_past > 1000);
}

TEST_CASE("bwa_gen_cigar4: want_md 0 gives the same score, CIGAR and NM as the MD form" *
          doctest::test_suite("unit/gen_cigar")) {
    const std::string ref = make_reference(kRefLen);
    meth_test::OrigRefFixture fx(ref);
    std::mt19937 rng(7);
    OpSeen seen;
    for (int freed_score : {1, 0}) {
        const Matrices m(freed_score);
        for (int iter = 0; iter < 1000; ++iter) {
            const bool rev = iter & 1;
            const int pad_l = (iter % 7 == 0) ? 1 + (int)(rng() % 3) : 0;
            const int pad_r = (iter % 11 == 0) ? 1 + (int)(rng() % 3) : 0;
            const Case c = make_case(ref, rng, rev, pad_l, pad_r);
            static const int kBands[] = {0, 3, 10, 50, 100};
            const int band = kBands[rng() % 5];
            // the literal NM under both matrices, and the --meth NM under the asymmetric one
            for (int combo = 0; combo < 3; ++combo) {
                const int nm_from_mat = combo == 2;
                const int8_t *mat = combo == 0 ? m.sym : m.meth;
                CAPTURE(iter);
                CAPTURE(freed_score);
                CAPTURE(combo);
                CAPTURE(band);
                GenCigarOut md = run_gen_cigar(fx, mat, c.query, c.rb, c.re, band, nm_from_mat, 1);
                GenCigarOut no = run_gen_cigar(fx, mat, c.query, c.rb, c.re, band, nm_from_mat, 0);
                REQUIRE(md.cigar != nullptr);
                REQUIRE(no.cigar != nullptr);
                CHECK(no.score == md.score);
                CHECK(no.n_cigar == md.n_cigar);
                CHECK(no.nm == md.nm);
                CHECK(words(no.cigar, no.n_cigar) == words(md.cigar, md.n_cigar));
                CHECK(no.query_after == c.query);  // the reverse-strand flip is undone
                CHECK(md.query_after == c.query);
                tally(words(md.cigar, md.n_cigar), seen);
                if (nm_from_mat) {
                    // NM differs from the literal count iff a freed C->T column was skipped
                    GenCigarOut lit = run_gen_cigar(fx, mat, c.query, c.rb, c.re, band, 0, 0);
                    if (lit.nm != no.nm) ++seen.freed;
                    free(lit.cigar);
                }
                free(md.cigar);
                free(no.cigar);
            }
        }
    }
    // The random cases reached every branch of the NM walk.
    CHECK(seen.ungapped > 100);
    CHECK(seen.internal_del > 100);
    CHECK(seen.end_del > 20);
    CHECK(seen.ins > 100);
    CHECK(seen.freed > 100);
}

TEST_CASE("bwa_gen_cigar4: want_md 0 NM values on hand-built alignments" *
          doctest::test_suite("unit/gen_cigar")) {
    const std::string ref = make_reference(kRefLen);
    meth_test::OrigRefFixture fx(ref);
    const Matrices m(1);
    const int at = 1000, len = 60;
    auto query_of = [&](const std::string &s) {
        std::vector<uint8_t> q;
        for (char ch : s) q.push_back(nt4(ch));
        return q;
    };
    std::string base = ref.substr(at, len);

    SUBCASE("two mismatches, ungapped") {
        std::string s = base;
        s[10] = s[10] == 'A' ? 'G' : 'A';
        s[40] = s[40] == 'A' ? 'G' : 'A';
        GenCigarOut no = run_gen_cigar(fx, m.sym, query_of(s), at, at + len, 0, 0, 0);
        REQUIRE(no.cigar != nullptr);
        CHECK(no.n_cigar == 1);
        CHECK(no.cigar[0] == (uint32_t)(len << 4));
        CHECK(no.nm == 2);
        CHECK(no.score == len - 2 - 2 * 4);
        free(no.cigar);
    }
    SUBCASE("query N counts as an edit") {
        std::string s = base;
        s[20] = 'N';
        GenCigarOut no = run_gen_cigar(fx, m.sym, query_of(s), at, at + len, 0, 0, 0);
        REQUIRE(no.cigar != nullptr);
        CHECK(no.nm == 1);
        free(no.cigar);
    }
    SUBCASE("a 2-base internal deletion counts 2") {
        std::string s = base.substr(0, 30) + base.substr(32);
        // the first M is 30 up to however far the deletion can slide in a repeat
        int lo = 30, hi = 30;
        while (lo > 0 && base[lo - 1] == base[lo + 1]) --lo;
        while (hi + 2 < len && base[hi] == base[hi + 2]) ++hi;
        GenCigarOut no = run_gen_cigar(fx, m.sym, query_of(s), at, at + len, 20, 0, 0);
        REQUIRE(no.cigar != nullptr);
        REQUIRE(no.n_cigar == 3);
        const int m1 = (int)(no.cigar[0] >> 4);
        CHECK(m1 >= lo);
        CHECK(m1 <= hi);
        CHECK(words(no.cigar, 3) == std::vector<uint32_t>{(uint32_t)m1 << 4, 2u << 4 | 2, (uint32_t)(len - 2 - m1) << 4});
        CHECK(no.nm == 2);
        free(no.cigar);
    }
    SUBCASE("a 3-base insertion counts 3") {
        std::string s = base.substr(0, 30) + "ACG" + base.substr(30);
        // the first M is 30 up to however far the insertion can slide in a repeat
        int lo = 30, hi = 30;
        while (lo > 0 && s[lo - 1] == s[lo + 2]) --lo;
        while (hi + 3 < (int)s.size() && s[hi] == s[hi + 3]) ++hi;
        GenCigarOut no = run_gen_cigar(fx, m.sym, query_of(s), at, at + len, 20, 0, 0);
        REQUIRE(no.cigar != nullptr);
        REQUIRE(no.n_cigar == 3);
        const int m1 = (int)(no.cigar[0] >> 4);
        CHECK(m1 >= lo);
        CHECK(m1 <= hi);
        CHECK(words(no.cigar, 3) == std::vector<uint32_t>{(uint32_t)m1 << 4, 3u << 4 | 1, (uint32_t)(len - m1) << 4});
        CHECK(no.nm == 3);
        free(no.cigar);
    }
    SUBCASE("a leading deletion is not counted") {
        // the window starts 2 bases before the read
        GenCigarOut no = run_gen_cigar(fx, m.sym, query_of(base), at - 2, at + len, 20, 0, 0);
        REQUIRE(no.cigar != nullptr);
        CHECK(no.n_cigar == 2);
        CHECK(no.cigar[0] == (2u << 4 | 2));
        CHECK(no.nm == 0);
        free(no.cigar);
    }
    SUBCASE("a freed C->T is an edit only for the literal NM") {
        size_t c = base.find('C', 5);
        REQUIRE(c != std::string::npos);
        REQUIRE(c < (size_t)len - 5);
        std::string s = base;
        s[c] = 'T';
        GenCigarOut lit = run_gen_cigar(fx, m.meth, query_of(s), at, at + len, 0, 0, 0);
        GenCigarOut mat = run_gen_cigar(fx, m.meth, query_of(s), at, at + len, 0, 1, 0);
        REQUIRE(lit.cigar != nullptr);
        REQUIRE(mat.cigar != nullptr);
        CHECK(lit.nm == 1);
        CHECK(mat.nm == 0);
        free(lit.cigar);
        free(mat.cigar);
    }
}

// Both mem_reg2aln forms of `r` must agree on every field mem_gen_alt reads,
// and the MD form must still carry its MD through the squeeze and clip moves.
// Returns the NM.
int check_reg2aln_forms(const mem_opt_t *opt, const meth_test::OrigRefFixture &fx,
                        const std::vector<char> &seq, const mem_alnreg_t &r,
                        const char *meth_orig_query, int *n_squeezed_gapped_clipped) {
    mem_aln_t md = mem_reg2aln(opt, fx.orig_bns, fx.orig_pac, (int)seq.size(), seq.data(), &r,
                               meth_orig_query, 1);
    mem_aln_t no = mem_reg2aln(opt, fx.orig_bns, fx.orig_pac, (int)seq.size(), seq.data(), &r,
                               meth_orig_query, 0);
    REQUIRE(md.cigar != nullptr);
    REQUIRE(no.cigar != nullptr);
    CHECK(no.rid == md.rid);
    CHECK(no.pos == md.pos);
    CHECK((int)no.is_rev == (int)md.is_rev);
    CHECK((int)no.flag == (int)md.flag);
    CHECK((int)no.mapq == (int)md.mapq);
    CHECK((int)no.NM == (int)md.NM);
    CHECK(no.score == md.score);
    CHECK(no.sub == md.sub);
    CHECK(no.n_cigar == md.n_cigar);
    CHECK(words(no.cigar, no.n_cigar) == words(md.cigar, md.n_cigar));
    CHECK(md_ref_span((const char *)(md.cigar + md.n_cigar)) == cigar_ref_span(md.cigar, md.n_cigar));
    if (r.rb < fx.bns.l_pac && md.pos > r.rb) ++n_squeezed_gapped_clipped[0];  // a leading D squeezed into pos
    for (int k = 0; k < md.n_cigar; ++k)
        if ((md.cigar[k] & 0xf) == 1 || (md.cigar[k] & 0xf) == 2) { ++n_squeezed_gapped_clipped[1]; break; }
    if ((md.cigar[0] & 0xf) == 3 || (md.cigar[md.n_cigar - 1] & 0xf) == 3) ++n_squeezed_gapped_clipped[2];
    const int nm = md.NM;
    free(md.cigar);
    free(no.cigar);
    return nm;
}

TEST_CASE("mem_reg2aln: want_md 0 matches the MD form in every field the XA entry reads" *
          doctest::test_suite("unit/gen_cigar")) {
    const std::string ref = make_reference(kRefLen);
    meth_test::OrigRefFixture fx(ref);
    mem_opt_t *opt = mem_opt_init();
    // --meth: the native regen (the read's original bases, the per-hypothesis matrix and the
    // matrix NM rule) that mem_gen_alt's XA hits take, one arm per scoring mode
    mem_opt_t *meth[3];
    for (int m = 0; m < 3; ++m) {
        meth[m] = mem_opt_init();
        meth[m]->meth_mode = 1;
        meth[m]->meth_scoring = m;
        mem_opt_fill_meth_mat(meth[m]);
    }
    std::mt19937 rng(11);
    int seen[3] = {0, 0, 0};  // squeezed, gapped, clipped
    int meth_nm_differs = 0;
    for (int iter = 0; iter < 900; ++iter) {
        CAPTURE(iter);
        const bool rev = iter & 1;
        const int pad_l = (iter % 5 == 0) ? 1 + (int)(rng() % 2) : 0;
        const int pad_r = (iter % 9 == 0) ? 1 + (int)(rng() % 2) : 0;
        const Case c = make_case(ref, rng, rev, pad_l, pad_r);
        // the aligned core flanked by unaligned bases, which mem_reg2aln soft-clips
        const int lclip = (iter % 3 == 0) ? (int)(rng() % 20) : 0;
        const int rclip = (iter % 4 == 0) ? (int)(rng() % 20) : 0;
        std::vector<char> seq;
        for (int i = 0; i < lclip; ++i) seq.push_back((char)(rng() & 3));
        for (uint8_t b : c.query) seq.push_back((char)b);
        for (int i = 0; i < rclip; ++i) seq.push_back((char)(rng() & 3));
        const int qlen = (int)c.query.size();
        mem_alnreg_t r{};
        r.rb = c.rb;
        r.re = c.re;
        r.qb = lclip;
        r.qe = lclip + qlen;
        r.rid = 0;
        // a deficit across both sides of the zero-band threshold
        r.score = r.truesc = qlen - (int)(rng() % 3) * 5 * c.n_edits;
        r.secondary = (iter % 6 == 0) ? 0 : -1;
        r.w = opt->w;
        r.seedcov = r.seedlen0 = qlen;
        r.meth_hypothesis = -1;
        const int nm_plain = check_reg2aln_forms(opt, fx, seq, r, NULL, seen);
        // the same region under --meth, as an OT read (hypothesis 1; C->T is the freed cell)
        CAPTURE(iter % 3);
        r.meth_hypothesis = 1;
        const int nm_meth = check_reg2aln_forms(meth[iter % 3], fx, seq, r, seq.data(), seen);
        if (nm_meth != nm_plain) ++meth_nm_differs;
        r.meth_hypothesis = 0;
        check_reg2aln_forms(meth[iter % 3], fx, seq, r, seq.data(), seen);
    }
    free(opt);
    for (mem_opt_t *m : meth) free(m);
    CHECK(seen[0] > 50);
    CHECK(seen[1] > 300);
    CHECK(seen[2] > 300);
    // the --meth arm took the native regen: its NM skipped freed C->T columns
    CHECK(meth_nm_differs > 100);
}

// mem_gen_alt's XA string for a 150-base read whose segment is copied three
// more times into the reference: with two mismatches, with a 2-base deletion and
// one mismatch, and reverse complemented with one mismatch. Under --meth (`meth`
// non-NULL) the read also carries two bisulfite C->T conversions, which the XA
// hits' NM must not count, on either strand.
void check_xa_string(mem_opt_t *meth) {
    std::string ref = make_reference(kRefLen);
    const int seg_at = 100, seg_len = 150;
    const std::string seg = ref.substr(seg_at, seg_len);
    auto flip = [](char b) { return b == 'A' ? 'C' : 'A'; };
    std::string c1 = seg;
    c1[20] = flip(c1[20]);
    c1[90] = flip(c1[90]);
    std::string c2 = seg.substr(0, 60) + "TT" + seg.substr(60);
    if (c2[59] == 'T' || c2[62] == 'T') c2 = seg.substr(0, 60) + "GG" + seg.substr(60);
    c2[120] = flip(c2[120]);
    std::string c3(seg.rbegin(), seg.rend());
    for (char &b : c3) b = b == 'A' ? 'T' : b == 'C' ? 'G' : b == 'G' ? 'C' : 'A';
    c3[40] = flip(c3[40]);
    const int at1 = 1000, at2 = 1600, at3 = 2200;
    ref.replace(at1, c1.size(), c1);
    ref.replace(at2, c2.size(), c2);
    ref.replace(at3, c3.size(), c3);
    meth_test::OrigRefFixture fx(ref);
    mem_opt_t *plain = mem_opt_init();
    mem_opt_t *opt = meth ? meth : plain;
    std::string read = seg;
    if (meth) {  // convert the first two Cs past the edits (read positions 125 on)
        int n = 0;
        for (int i = 125; i < seg_len && n < 2; ++i)
            if (read[i] == 'C') read[i] = 'T', ++n;
        REQUIRE(n == 2);
    }
    std::vector<char> seq;
    for (char b : read) seq.push_back((char)nt4(b));

    std::vector<mem_alnreg_t> regs(4);
    auto region = [&](mem_alnreg_t &r, int64_t rb, int64_t re, int score, int secondary_all) {
        r = mem_alnreg_t{};
        r.rb = rb;
        r.re = re;
        r.qb = 0;
        r.qe = seg_len;
        r.rid = 0;
        r.score = r.truesc = score;
        r.secondary = secondary_all;
        r.secondary_all = secondary_all;
        r.w = opt->w;
        r.seedcov = r.seedlen0 = seg_len;
        // the genome-strand hypothesis of a C->T read: OT on a forward hit, OB on a
        // reverse one, which mem_reg2aln flips back to the read's C->T matrix
        r.meth_hypothesis = meth ? (rb < fx.bns.l_pac ? 1 : 0) : -1;
    };
    region(regs[0], seg_at, seg_at + seg_len, 150, -1);
    region(regs[1], at1, at1 + seg_len, 150 - 2 * 5, 0);
    region(regs[2], at2, at2 + seg_len + 2, 150 - 5 - 8, 0);
    region(regs[3], 2 * (int64_t)kRefLen - (at3 + seg_len), 2 * (int64_t)kRefLen - at3, 150 - 5, 0);
    mem_alnreg_v a{regs.size(), regs.size(), regs.data()};
    int *hn = nullptr;
    char **xa = mem_gen_alt(opt, fx.orig_bns, fx.orig_pac, &a, seg_len, seq.data(), &hn,
                            meth ? seq.data() : NULL);
    REQUIRE(xa != nullptr);
    REQUIRE(xa[0] != nullptr);
    // the deletion's first M may slide left or right in a repeat
    int lo = 60, hi = 60;
    while (lo > 0 && c2[lo - 1] == c2[lo + 1]) --lo;
    while (hi + 2 < (int)c2.size() && c2[hi] == c2[hi + 2]) ++hi;
    bool matched = false;
    for (int m1 = lo; m1 <= hi && !matched; ++m1) {
        const std::string want = "X,+" + std::to_string(at1 + 1) + ",150M,2;" + "X,+" +
                                 std::to_string(at2 + 1) + "," + std::to_string(m1) + "M2D" +
                                 std::to_string(seg_len - m1) + "M,3;" + "X,-" +
                                 std::to_string(at3 + 1) + ",150M,1;";
        matched = std::string(xa[0]) == want;
    }
    const std::string got = xa[0];
    CAPTURE(got);
    CHECK(matched);
    CHECK(hn[0] == 3);
    for (size_t k = 0; k < regs.size(); ++k) free(xa[k]);
    free(xa);
    free(hn);
    free(plain);
}

TEST_CASE("mem_gen_alt: the XA string carries each hit's CIGAR and NM" *
          doctest::test_suite("unit/gen_cigar")) {
    SUBCASE("without --meth") { check_xa_string(nullptr); }
    SUBCASE("--meth: bisulfite conversions are not edits") {
        for (int scoring : {MEM_METH_SCORING_COLLAPSED, MEM_METH_SCORING_GENOMIC,
                            MEM_METH_SCORING_NEUTRAL}) {
            CAPTURE(scoring);
            mem_opt_t *meth = mem_opt_init();
            meth->meth_mode = 1;
            meth->meth_scoring = scoring;
            mem_opt_fill_meth_mat(meth);
            check_xa_string(meth);
            free(meth);
        }
    }
}
