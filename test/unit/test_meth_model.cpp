// Unit tests for the methylation model (src/meth_model.{h,cpp}): the 7-code
// reference alphabet, the 7x7 per-hypothesis matrices, and in-place CpG
// labelling of a fetched reference window.
//
// Labelling is defined in the window's fetch frame: a C whose next frame base
// is G becomes METH_C5, a G whose previous frame base is C becomes METH_G6.
// Reverse-strand windows (doubled coordinate rb >= l_pac) arrive already
// reverse-complemented from bns_get_seq_into, so their frame is the reverse
// complement of the forward contig. The oracle below labels the whole contig
// (or its reverse complement) as a string and slices the window out of it, so
// every window-edge neighbour is taken in the same frame by construction.

#include "doctest/doctest.h"
#include "bwamem.h"
#include "bwa.h"
#include "bntseq.h"
#include "meth_model.h"
#include "meth_orig_ref_fixture.h"

#include <cstdlib>
#include <string>
#include <vector>

namespace {

using meth_test::OrigRefFixture;

// Overwrite one 2-bit base in the packed pac (pac_pack_set only ORs).
void pac_put(std::vector<uint8_t> &pac, int64_t i, uint8_t v)
{
    const int sh = (int)((~i & 3) << 1);
    pac[i >> 2] = (uint8_t)((pac[i >> 2] & ~(3u << sh)) | (v << sh));
}

// Fetch window [rb, rb + len) of the doubled coordinate space, exactly as the
// aligner does, and label it.
std::vector<uint8_t> fetch_label(const bntseq_t *bns, const uint8_t *pac,
                                 int64_t rb, int64_t len, int rid)
{
    std::vector<uint8_t> w((size_t)len + 1, 0xff);
    int64_t got = 0;
    bns_get_seq_into(bns->l_pac, pac, rb, rb + len, w.data(), &got);
    REQUIRE(got == len);
    w.resize((size_t)len);
    meth_label_cpg(w.data(), len, rb, rid, bns, pac);
    return w;
}

char comp(char c)
{
    switch (c) {
        case 'A': return 'T'; case 'C': return 'G';
        case 'G': return 'C'; case 'T': return 'A';
        default:  return 'N';
    }
}

// Expected labels of one contig in a given frame (frame = forward contig
// string, or its reverse complement). N never forms a CpG.
std::vector<int> oracle_frame(const std::string &frame)
{
    std::vector<int> out(frame.size());
    for (size_t i = 0; i < frame.size(); ++i) {
        const char c = frame[i];
        int code = c == 'A' ? METH_A : c == 'C' ? METH_C : c == 'G' ? METH_G
                 : c == 'T' ? METH_T : -1;   // -1: N, code unchecked
        if (c == 'C' && i + 1 < frame.size() && frame[i + 1] == 'G') code = METH_C5;
        if (c == 'G' && i > 0 && frame[i - 1] == 'C') code = METH_G6;
        out[i] = code;
    }
    return out;
}

std::string revcomp(const std::string &s)
{
    std::string r(s.rbegin(), s.rend());
    for (char &c : r) c = comp(c);
    return r;
}

// Check a labelled window against the oracle. Forward window: frame index k
// is forward base rb + k. Reverse window: frame index k is reverse-complement
// index (rb - l_pac) + k. Positions that are N in the contig are skipped
// (their code is the random pac base); only their non-labelling matters, and
// that is covered by the neighbours' expected values plus the "no 5/6 on an
// N" check.
void check_window(const std::string &contig, const std::vector<uint8_t> &w,
                  int64_t rb, int64_t l_pac)
{
    const bool rev = rb >= l_pac;
    const std::string frame = rev ? revcomp(contig) : contig;
    const std::vector<int> exp = oracle_frame(frame);
    const int64_t start = rev ? rb - l_pac : rb;
    for (size_t k = 0; k < w.size(); ++k) {
        const int e = exp[(size_t)start + k];
        CAPTURE(rb); CAPTURE(k);
        if (e < 0) CHECK(w[k] < METH_C5);
        else       CHECK((int)w[k] == e);
    }
}

// A multi-contig in-memory bns + pac (contigs adjacent in pac, no gaps).
struct MultiRef {
    bntseq_t                bns{};
    std::vector<bntann1_t>  anns;
    std::vector<bntamb1_t>  ambs;
    std::vector<uint8_t>    pac;
    std::string             all;

    explicit MultiRef(const std::vector<std::string> &contigs)
    {
        for (const auto &c : contigs) {
            bntann1_t a{};
            a.offset = (int64_t)all.size();
            a.len    = (int32_t)c.size();
            a.name   = const_cast<char *>("c");
            a.anno   = const_cast<char *>("");
            anns.push_back(a);
            all += c;
        }
        bns.l_pac  = (int64_t)all.size();
        bns.n_seqs = (int32_t)anns.size();
        bns.anns   = anns.data();
        pac.assign((all.size() + 3) / 4, 0);
        for (size_t i = 0; i < all.size(); ++i) {
            if (all[i] == 'N') {
                bntamb1_t a{};
                a.offset = (int64_t)i; a.len = 1; a.amb = 'N';
                ambs.push_back(a);
            }
            meth_test::pac_pack_set(pac.data(), (int64_t)i, meth_test::base_to_2bit(all[i]));
        }
        bns.ambs    = ambs.data();
        bns.n_holes = (int32_t)ambs.size();
    }
};

mem_opt_t *opt_for(int scoring, int a, int b)
{
    mem_opt_t *o = mem_opt_init();
    o->a = a; o->b = b;
    bwa_fill_scmat(a, b, o->mat);
    o->meth_scoring = scoring;
    mem_opt_fill_meth_mat(o);
    return o;
}

}  // namespace

TEST_CASE("meth_code_base maps labelled codes to their base"
          * doctest::test_suite("unit/meth_model"))
{
    CHECK(meth_code_base(METH_A) == METH_A);
    CHECK(meth_code_base(METH_C) == METH_C);
    CHECK(meth_code_base(METH_G) == METH_G);
    CHECK(meth_code_base(METH_T) == METH_T);
    CHECK(meth_code_base(METH_N) == METH_N);
    CHECK(meth_code_base(METH_C5) == METH_C);
    CHECK(meth_code_base(METH_G6) == METH_G);
}

TEST_CASE("mat7 rows/cols 0..4 equal mat" * doctest::test_suite("unit/meth_model"))
{
    const int scorings[] = { MEM_METH_SCORING_COLLAPSED, MEM_METH_SCORING_GENOMIC,
                             MEM_METH_SCORING_NEUTRAL };
    const int ab[][2] = { {1, 4}, {2, 4}, {1, 1}, {3, 7}, {1, 0} };
    for (int s : scorings) {
        for (const auto &p : ab) {
            mem_opt_t *o = opt_for(s, p[0], p[1]);
            CAPTURE(s); CAPTURE(p[0]); CAPTURE(p[1]);
            for (int r = 0; r < METH_M; ++r) {
                for (int c = 0; c < METH_M; ++c) {
                    const int8_t base = o->mat[meth_code_base(r) * 5 + meth_code_base(c)];
                    const bool freed_ot = (r == METH_C5 && c == METH_T);
                    const bool freed_ob = (r == METH_G6 && c == METH_A);
                    if (!freed_ot) CHECK(o->mat7_ot[r * 7 + c] == base);
                    if (!freed_ob) CHECK(o->mat7_ob[r * 7 + c] == base);
                }
            }
            free(o);
        }
    }
}

TEST_CASE("mat7 freed cells" * doctest::test_suite("unit/meth_model"))
{
    mem_opt_t *o = opt_for(MEM_METH_SCORING_NEUTRAL, 1, 4);
    CHECK(o->mat7_ot[5 * 7 + 3] == 0);
    CHECK(o->mat7_ot[1 * 7 + 3] == -4);
    CHECK(o->mat7_ob[6 * 7 + 0] == 0);
    CHECK(o->mat7_ob[2 * 7 + 0] == -4);
    // The other hypothesis' cell is not freed.
    CHECK(o->mat7_ot[6 * 7 + 0] == -4);
    CHECK(o->mat7_ob[5 * 7 + 3] == -4);
    free(o);

    o = opt_for(MEM_METH_SCORING_GENOMIC, 1, 4);
    CHECK(o->mat7_ot[5 * 7 + 3] == 1);
    CHECK(o->mat7_ob[6 * 7 + 0] == 1);
    CHECK(o->mat7_ot[1 * 7 + 3] == -4);
    free(o);

    o = opt_for(MEM_METH_SCORING_GENOMIC, 2, 6);
    CHECK(o->mat7_ot[5 * 7 + 3] == 2);
    CHECK(o->mat7_ot[1 * 7 + 3] == -6);
    free(o);
}

TEST_CASE("meth_label_cpg label interior" * doctest::test_suite("unit/meth_model"))
{
    OrigRefFixture f("ACGTTCGA");
    std::vector<uint8_t> w = fetch_label(f.orig_bns, f.orig_pac, 0, 8, 0);
    CHECK(w == std::vector<uint8_t>({0, 5, 6, 3, 3, 5, 6, 0}));
    meth_unlabel(w.data(), (int64_t)w.size());
    CHECK(w == std::vector<uint8_t>({0, 1, 2, 3, 3, 1, 2, 0}));
}

TEST_CASE("meth_label_cpg window-edge lookahead/lookbehind"
          * doctest::test_suite("unit/meth_model"))
{
    //                  0123456789
    OrigRefFixture f("TTCGTTCGTT");
    // Window [3, 7): starts on the G of CG@2, ends on the C of CG@6.
    std::vector<uint8_t> w = fetch_label(f.orig_bns, f.orig_pac, 3, 4, 0);
    CHECK(w == std::vector<uint8_t>({METH_G6, METH_T, METH_T, METH_C5}));
    // Same as labelling one base wider on each side, then trimming.
    std::vector<uint8_t> wide = fetch_label(f.orig_bns, f.orig_pac, 2, 6, 0);
    CHECK(std::vector<uint8_t>(wide.begin() + 1, wide.end() - 1) == w);
    // One-base windows on each half of a CpG.
    CHECK(fetch_label(f.orig_bns, f.orig_pac, 2, 1, 0) == std::vector<uint8_t>({METH_C5}));
    CHECK(fetch_label(f.orig_bns, f.orig_pac, 3, 1, 0) == std::vector<uint8_t>({METH_G6}));
}

TEST_CASE("meth_label_cpg contig boundary" * doctest::test_suite("unit/meth_model"))
{
    MultiRef r({"TTAC", "GTTA"});
    // Last window of contig 0, first window of contig 1: no CpG across them.
    CHECK(fetch_label(&r.bns, r.pac.data(), 0, 4, 0)
          == std::vector<uint8_t>({3, 3, 0, 1}));
    CHECK(fetch_label(&r.bns, r.pac.data(), 4, 4, 1)
          == std::vector<uint8_t>({2, 3, 3, 0}));
    // Reverse strand: contig 1 reversed is TAAC, contig 0 reversed is GTAA;
    // in frame contig 1 precedes contig 0 and its C still meets no G.
    const int64_t L = r.bns.l_pac;
    CHECK(fetch_label(&r.bns, r.pac.data(), L, 4, 1)
          == std::vector<uint8_t>({3, 0, 0, 1}));
    CHECK(fetch_label(&r.bns, r.pac.data(), L + 4, 4, 0)
          == std::vector<uint8_t>({2, 3, 0, 0}));
}

TEST_CASE("meth_label_cpg ambs masking" * doctest::test_suite("unit/meth_model"))
{
    {
        // An N stored in pac as C, followed by G: no label on either.
        OrigRefFixture f("TTNGTT");
        pac_put(f.pac, 2, METH_C);
        std::vector<uint8_t> w = fetch_label(f.orig_bns, f.orig_pac, 0, 6, 0);
        CHECK(w == std::vector<uint8_t>({3, 3, 1, 2, 3, 3}));
        // The same N as a window-edge neighbour (lookbehind of a window
        // starting on the G).
        CHECK(fetch_label(f.orig_bns, f.orig_pac, 3, 2, 0)
              == std::vector<uint8_t>({2, 3}));
    }
    {
        // A C then an N stored as G, also as a lookahead neighbour.
        OrigRefFixture f("TTCNTT");
        pac_put(f.pac, 3, METH_G);
        CHECK(fetch_label(f.orig_bns, f.orig_pac, 0, 6, 0)
              == std::vector<uint8_t>({3, 3, 1, 2, 3, 3}));
        CHECK(fetch_label(f.orig_bns, f.orig_pac, 0, 3, 0)
              == std::vector<uint8_t>({3, 3, 1}));
    }
    {
        // A C followed by a G, both inside an N run.
        OrigRefFixture f("TNNNNT");
        pac_put(f.pac, 2, METH_C);
        pac_put(f.pac, 3, METH_G);
        CHECK(fetch_label(f.orig_bns, f.orig_pac, 0, 6, 0)
              == std::vector<uint8_t>({3, 0, 1, 2, 0, 3}));
        // And on the reverse strand (revcomp of pac TACGAT is ATCGTA).
        CHECK(fetch_label(f.orig_bns, f.orig_pac, 6, 6, 0)
              == std::vector<uint8_t>({0, 3, 1, 2, 3, 0}));
    }
}

TEST_CASE("meth_label_cpg reverse-strand frame" * doctest::test_suite("unit/meth_model"))
{
    {
        OrigRefFixture f("ACGT");
        CHECK(fetch_label(f.orig_bns, f.orig_pac, 4, 4, 0)
              == std::vector<uint8_t>({0, 5, 6, 3}));
    }
    {
        // Forward CG at the contig's first two bases is the last two bases
        // of the reverse window (revcomp TTAACG).
        OrigRefFixture f("CGTTAA");
        const int64_t L = f.bns.l_pac;
        CHECK(fetch_label(f.orig_bns, f.orig_pac, L, 6, 0)
              == std::vector<uint8_t>({3, 3, 0, 0, 5, 6}));
        CHECK(fetch_label(f.orig_bns, f.orig_pac, L + 4, 2, 0)
              == std::vector<uint8_t>({5, 6}));
    }
}

TEST_CASE("meth_label_cpg reverse window with a CpG at each edge"
          * doctest::test_suite("unit/meth_model"))
{
    // revcomp("AACGTTACGTT") = "AACGTAACGTT"; reverse window over revcomp
    // indices [3, 8) = "GTAAC": its first G follows frame C (revcomp[2]) and
    // its last C precedes frame G (revcomp[8]). The window covers forward
    // [3, 8); the lookbehind is the complement of forward base 8 (G), the one
    // AFTER the window's forward end, and the lookahead the complement of
    // forward base 2 (C), the one BEFORE its forward start.
    const std::string contig = "AACGTTACGTT";
    OrigRefFixture f(contig);
    const int64_t L = f.bns.l_pac;
    REQUIRE(revcomp(contig) == "AACGTAACGTT");
    std::vector<uint8_t> w = fetch_label(f.orig_bns, f.orig_pac, L + 3, 5, 0);
    CHECK(w == std::vector<uint8_t>({METH_G6, METH_T, METH_A, METH_A, METH_C5}));
    std::vector<uint8_t> wide = fetch_label(f.orig_bns, f.orig_pac, L + 2, 7, 0);
    CHECK(std::vector<uint8_t>(wide.begin() + 1, wide.end() - 1) == w);
    check_window(contig, w, L + 3, L);
}

TEST_CASE("meth_label_cpg matches the frame oracle on random windows"
          * doctest::test_suite("unit/meth_model"))
{
    srand48(11);
    const char alpha[] = "ACGTCGCGN";   // CpG-rich, some N
    for (int rep = 0; rep < 40; ++rep) {
        std::string contig;
        const int n = 8 + (int)(lrand48() % 40);
        for (int i = 0; i < n; ++i) contig += alpha[lrand48() % 9];
        OrigRefFixture f(contig);
        // Store each N as a random base, as the indexer does.
        for (int i = 0; i < n; ++i)
            if (contig[i] == 'N') pac_put(f.pac, i, (uint8_t)(lrand48() & 3));
        const int64_t L = f.bns.l_pac;
        for (int64_t beg = 0; beg < 2 * L; ++beg) {
            const int64_t half_end = beg < L ? L : 2 * L;
            for (int64_t end = beg + 1; end <= half_end; ++end) {
                std::vector<uint8_t> w = fetch_label(f.orig_bns, f.orig_pac, beg, end - beg, 0);
                check_window(contig, w, beg, L);
            }
        }
    }
}
