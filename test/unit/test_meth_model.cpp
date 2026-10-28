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
#include "matesw_u8.h"
#include "ksw.h"

#include <cstdlib>
#include <random>
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
// aligner does.
std::vector<uint8_t> fetch_raw(const bntseq_t *bns, const uint8_t *pac,
                               int64_t rb, int64_t len)
{
    std::vector<uint8_t> w((size_t)len + 1, 0xff);
    int64_t got = 0;
    bns_get_seq_into(bns->l_pac, pac, rb, rb + len, w.data(), &got);
    REQUIRE(got == len);
    w.resize((size_t)len);
    return w;
}

// Fetch the window and label it.
std::vector<uint8_t> fetch_label(const bntseq_t *bns, const uint8_t *pac,
                                 int64_t rb, int64_t len, int rid)
{
    std::vector<uint8_t> w = fetch_raw(bns, pac, rb, len);
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

TEST_CASE("meth_label_cpg matches the per-contig oracle on windows straddling contigs"
          * doctest::test_suite("unit/meth_model"))
{
    // Windows anywhere in the doubled coordinate space, possibly spanning a
    // contig boundary; rid is the contig of the window's midpoint, as the
    // aligner picks it. Positions outside rid stay unlabelled (raw codes);
    // positions inside match the oracle run on contig rid alone, in frame.
    srand48(23);
    const char alpha[] = "ACGTCGCGN";
    for (int rep = 0; rep < 30; ++rep) {
        std::vector<std::string> contigs;
        const int nc = 2 + (int)(lrand48() % 3);
        for (int j = 0; j < nc; ++j) {
            std::string c;
            const int n = 1 + (int)(lrand48() % 12);
            for (int i = 0; i < n; ++i) c += alpha[lrand48() % 9];
            contigs.push_back(c);
        }
        MultiRef r(contigs);
        for (size_t i = 0; i < r.all.size(); ++i)
            if (r.all[i] == 'N') pac_put(r.pac, (int64_t)i, (uint8_t)(lrand48() & 3));
        const int64_t L = r.bns.l_pac;
        for (int64_t beg = 0; beg < 2 * L; ++beg) {
            const int64_t half_end = beg < L ? L : 2 * L;
            for (int64_t end = beg + 1; end <= half_end; ++end) {
                const int64_t len = end - beg;
                int is_rev = 0;
                const int64_t mid_f = bns_depos(&r.bns, beg + len / 2, &is_rev);
                const int rid = bns_pos2rid(&r.bns, mid_f);
                const std::vector<uint8_t> raw = fetch_raw(&r.bns, r.pac.data(), beg, len);
                std::vector<uint8_t> w = raw;
                meth_label_cpg(w.data(), len, beg, rid, &r.bns, r.pac.data());

                const std::string &cg = contigs[(size_t)rid];
                const int64_t off = r.anns[(size_t)rid].offset, clen = (int64_t)cg.size();
                const std::vector<int> exp = oracle_frame(is_rev ? revcomp(cg) : cg);
                for (int64_t k = 0; k < len; ++k) {
                    const int64_t fc = is_rev ? (2 * L - 1 - (beg + k)) : beg + k;
                    CAPTURE(rep); CAPTURE(beg); CAPTURE(len); CAPTURE(k);
                    if (fc < off || fc >= off + clen) { CHECK(w[(size_t)k] == raw[(size_t)k]); continue; }
                    const int64_t t = fc - off;
                    const int e = exp[(size_t)(is_rev ? clen - 1 - t : t)];
                    if (e < 0) CHECK(w[(size_t)k] == raw[(size_t)k]);
                    else       CHECK((int)w[(size_t)k] == e);
                }
            }
        }
    }
}

TEST_CASE("meth_label_cpg leaves a window outside contig rid unchanged"
          * doctest::test_suite("unit/meth_model"))
{
    MultiRef r({"ACGCGT", "TCGCGA"});
    const int64_t L = r.bns.l_pac;
    // Forward window inside contig 0, labelled as contig 1; and the reverse
    // window of contig 1 labelled as contig 0.
    CHECK(fetch_label(&r.bns, r.pac.data(), 0, 6, 1) == fetch_raw(&r.bns, r.pac.data(), 0, 6));
    CHECK(fetch_label(&r.bns, r.pac.data(), L, 6, 0) == fetch_raw(&r.bns, r.pac.data(), L, 6));
    // Sanity: with the right rid both are labelled.
    CHECK(fetch_label(&r.bns, r.pac.data(), 0, 6, 0) != fetch_raw(&r.bns, r.pac.data(), 0, 6));
    CHECK(fetch_label(&r.bns, r.pac.data(), L, 6, 1) != fetch_raw(&r.bns, r.pac.data(), L, 6));
}

TEST_CASE("meth_label_cpg masks a merged multi-base ambs run exactly"
          * doctest::test_suite("unit/meth_model"))
{
    //                                   0123456789012345
    const std::string stored          = "TACGCGCGCGCGCGTA";   // pac bases
    const std::string contig_with_n   = "TACNNNNNNNNNCGTA";   // run [3, 12)
    MultiRef r({stored});
    // Replace the (empty) per-base ambs with one merged 9-bp run.
    r.ambs.clear();
    bntamb1_t a{};
    a.offset = 3; a.len = 9; a.amb = 'N';
    r.ambs.push_back(a);
    r.bns.ambs = r.ambs.data();
    r.bns.n_holes = 1;
    const int64_t L = r.bns.l_pac;

    // The full forward window: the C before the run and the CpG inside it
    // are unlabelled; the CpG after it is labelled.
    CHECK(fetch_label(&r.bns, r.pac.data(), 0, 16, 0)
          == std::vector<uint8_t>({3, 0, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 5, 6, 3, 0}));
    // A window strictly inside the run (the run extends past it both sides).
    CHECK(fetch_label(&r.bns, r.pac.data(), 5, 5, 0) == fetch_raw(&r.bns, r.pac.data(), 5, 5));
    CHECK(fetch_label(&r.bns, r.pac.data(), L + 16 - 10, 5, 0)
          == fetch_raw(&r.bns, r.pac.data(), L + 16 - 10, 5));
    // Every window on both strands against the oracle, with N positions raw.
    for (int64_t beg = 0; beg < 2 * L; ++beg) {
        const int64_t half_end = beg < L ? L : 2 * L;
        for (int64_t end = beg + 1; end <= half_end; ++end) {
            const std::vector<uint8_t> raw = fetch_raw(&r.bns, r.pac.data(), beg, end - beg);
            const std::vector<uint8_t> w = fetch_label(&r.bns, r.pac.data(), beg, end - beg, 0);
            check_window(contig_with_n, w, beg, L);
            const bool rev = beg >= L;
            for (int64_t k = 0; k < end - beg; ++k) {
                const int64_t fc = rev ? 2 * L - 1 - (beg + k) : beg + k;
                CAPTURE(beg); CAPTURE(k);
                if (fc >= 3 && fc < 12) CHECK(w[(size_t)k] == raw[(size_t)k]);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// --meth-context cpg on the scalar consumers: bwa_gen_cigar4 at m = 7, the
// scalar banded extension at m = 7, and the model helpers.
// ---------------------------------------------------------------------------

namespace {

std::vector<uint8_t> codes_of(const std::string &s)
{
    std::vector<uint8_t> v(s.size());
    for (size_t i = 0; i < s.size(); ++i) v[i] = meth_test::base_to_2bit(s[i]);
    return v;
}

struct Cigar4Out {
    int score = 0, nm = -1;
    std::string cigar, md;
};

// bwa_gen_cigar4 over doubled window [rb, re) with `query` (codes, in the
// window's fetch frame, as mem_reg2aln passes it), NM from the matrix and MD.
Cigar4Out gen_cigar4(const bntseq_t *bns, const uint8_t *pac, const int8_t *mat, int m,
                     std::vector<uint8_t> query, int64_t rb, int64_t re, int w)
{
    Cigar4Out o;
    int n_cigar = 0;
    uint32_t *c = bwa_gen_cigar4(mat, m, bns, 6, 1, 6, 1, w, bns->l_pac, pac,
                                 (int)query.size(), query.data(), rb, re, &o.score,
                                 &n_cigar, &o.nm, 1, 1);
    REQUIRE(c != NULL);
    for (int k = 0; k < n_cigar; ++k)
        o.cigar += std::to_string(c[k] >> 4) + "MIDSH"[c[k] & 0xf];
    o.md = (const char *)(c + n_cigar);
    free(c);
    return o;
}

}  // namespace

TEST_CASE("bwa_gen_cigar4 m=7 frees a conversion only at a CpG"
          * doctest::test_suite("unit/meth_model"))
{
    //                   0123456789012345
    const std::string g = "TTACGATTCATTGGAT";   // CpG C at 3, plain C at 8
    OrigRefFixture f(g);
    mem_opt_t *o = opt_for(MEM_METH_SCORING_NEUTRAL, 1, 4);
    std::string r = g;
    r[3] = 'T'; r[8] = 'T';
    // cpg: the CpG conversion is free (0), the plain one a mismatch (-4).
    Cigar4Out c7 = gen_cigar4(f.orig_bns, f.orig_pac, o->mat7_ot, 7, codes_of(r), 0, 16, 0);
    CHECK(c7.cigar == "16M");
    CHECK(c7.score == 14 * 1 + 0 - 4);
    CHECK(c7.nm == 1);
    CHECK(c7.md == "3C4C7");          // MD stays literal and prints the plain base
    // all (the 5-code matrix): both free.
    Cigar4Out c5 = gen_cigar4(f.orig_bns, f.orig_pac, o->mat_ot, 5, codes_of(r), 0, 16, 0);
    CHECK(c5.score == 14);
    CHECK(c5.nm == 0);
    CHECK(c5.md == "3C4C7");
    // The MD-free walk counts NM the same way.
    {
        std::vector<uint8_t> q = codes_of(r);
        int sc = 0, n_cigar = 0, nm = -1;
        uint32_t *c = bwa_gen_cigar4(o->mat7_ot, 7, f.orig_bns, 6, 1, 6, 1, 0, f.bns.l_pac,
                                     f.orig_pac, 16, q.data(), 0, 16, &sc, &n_cigar, &nm, 1, 0);
        REQUIRE(c != NULL);
        CHECK(nm == 1);
        CHECK(sc == c7.score);
        free(c);
    }
    free(o);
}

TEST_CASE("bwa_gen_cigar4 m=7 through the global DP (a deletion)"
          * doctest::test_suite("unit/meth_model"))
{
    //                   0         1         2         3
    //                   0123456789012345678901234567890123456789
    const std::string g = "GATTACGATTCAGTAGGATCCATAGTTAGCATAGGATTAC";
    OrigRefFixture f(g);
    mem_opt_t *o = opt_for(MEM_METH_SCORING_NEUTRAL, 1, 4);
    // Read: CpG C at 5 and plain C at 10 converted, base 22 (T) deleted.
    std::string r = g;
    r[5] = 'T'; r[10] = 'T';
    r.erase(22, 1);
    Cigar4Out c7 = gen_cigar4(f.orig_bns, f.orig_pac, o->mat7_ot, 7, codes_of(r), 0, 40, 10);
    CHECK(c7.cigar == "22M1D17M");
    CHECK(c7.score == 37 - 4 - 7);    // 37 matches, the plain C->T, a 1-base deletion
    CHECK(c7.nm == 2);
    CHECK(c7.md == "5C4C11^T17");
    Cigar4Out c5 = gen_cigar4(f.orig_bns, f.orig_pac, o->mat_ot, 5, codes_of(r), 0, 40, 10);
    CHECK(c5.cigar == c7.cigar);
    CHECK(c5.score == 37 - 7);
    CHECK(c5.nm == 1);
    // The MD-free walk counts the gapped NM the same way.
    {
        std::vector<uint8_t> q = codes_of(r);
        int sc = 0, n_cigar = 0, nm = -1;
        uint32_t *c = bwa_gen_cigar4(o->mat7_ot, 7, f.orig_bns, 6, 1, 6, 1, 10, f.bns.l_pac,
                                     f.orig_pac, (int)q.size(), q.data(), 0, 40, &sc, &n_cigar, &nm, 1, 0);
        REQUIRE(c != NULL);
        CHECK(n_cigar == 3);
        CHECK(nm == 2);
        CHECK(sc == c7.score);
        free(c);
    }
    free(o);
}

TEST_CASE("bwa_gen_cigar4 m=7: a reverse-strand window gets the forward twin's freed cells"
          * doctest::test_suite("unit/meth_model"))
{
    // Review focus: an OT read aligned on the reverse strand (hypothesis flipped
    // to OB) frees exactly the cells its forward-strand alignment frees.
    const std::string g = "TTACGATTCATTGGATCGAA";   // CpG Cs at 3 and 16, plain C at 8
    OrigRefFixture f(g);
    const int64_t L = f.bns.l_pac;
    mem_opt_t *o = opt_for(MEM_METH_SCORING_NEUTRAL, 1, 4);
    std::string r = g;
    r[3] = 'T'; r[8] = 'T'; r[16] = 'T';
    Cigar4Out fwd = gen_cigar4(f.orig_bns, f.orig_pac, o->mat7_ot, 7, codes_of(r), 0, L, 0);
    Cigar4Out rev = gen_cigar4(f.orig_bns, f.orig_pac, o->mat7_ob, 7, codes_of(revcomp(r)), L, 2 * L, 0);
    CHECK(fwd.score == 17 - 4);
    CHECK(fwd.nm == 1);
    CHECK(rev.score == fwd.score);
    CHECK(rev.nm == fwd.nm);
    CHECK(rev.cigar == fwd.cigar);
    CHECK(fwd.md == "3C4C7C3");
    // The reverse record's MD is in forward coordinates too: the same string.
    CHECK(rev.md == fwd.md);
    free(o);
}

TEST_CASE("bwa_gen_cigar4 m=7: a CpG straddling the window's last base is freed"
          * doctest::test_suite("unit/meth_model"))
{
    const std::string g = "AATTGGATCGTT";   // CpG at 8-9
    OrigRefFixture f(g);
    mem_opt_t *o = opt_for(MEM_METH_SCORING_NEUTRAL, 1, 4);
    // Window [0, 9) ends on the CpG's C; its G is outside the alignment.
    std::string r = g.substr(0, 9);
    r[8] = 'T';
    Cigar4Out c = gen_cigar4(f.orig_bns, f.orig_pac, o->mat7_ot, 7, codes_of(r), 0, 9, 0);
    CHECK(c.score == 8);
    CHECK(c.nm == 0);
    // And the reverse window that starts on the CpG's G (frame C before it).
    const int64_t L = f.bns.l_pac;
    // revcomp(g) = AACGATCCAATT: frame [3, 12) starts on the CpG's G, whose C
    // (frame 2) is outside the window.
    std::string rr = revcomp(g).substr(3, 9);
    REQUIRE(rr[0] == 'G');
    rr[0] = 'A';
    Cigar4Out cr = gen_cigar4(f.orig_bns, f.orig_pac, o->mat7_ob, 7, codes_of(rr), L + 3, L + 12, 0);
    CHECK(cr.score == 8);
    CHECK(cr.nm == 0);
    free(o);
}

TEST_CASE("scalarBandedSWA at m=7 equals ksw_extend2 at m=7 on labelled targets"
          * doctest::test_suite("unit/meth_model"))
{
    std::mt19937 rng(20261028);
    const char *bases = "ACGT";
    const int scorings[] = { MEM_METH_SCORING_GENOMIC, MEM_METH_SCORING_NEUTRAL };
    for (int sc : scorings) {
        mem_opt_t *o = opt_for(sc, 1, 4);
        for (int hyp = 0; hyp < 2; ++hyp) {
            const int8_t *mat = hyp ? o->mat7_ot : o->mat7_ob;
            auto k = make_banded_pair_wise_sw(o->o_del, o->e_del, o->o_ins, o->e_ins,
                                              o->zdrop, o->pen_clip5, mat, o->a, o->b, 1, 7);
            int n_diff_unlabelled = 0;
            for (int it = 0; it < 300; ++it) {
                const int tlen = 20 + (int)(rng() % 60), qlen = 10 + (int)(rng() % 50);
                std::string t;
                for (int i = 0; i < tlen; ++i) t += bases[rng() % 4];
                OrigRefFixture f(t);
                std::vector<uint8_t> raw = fetch_raw(f.orig_bns, f.orig_pac, 0, tlen);
                std::vector<uint8_t> lab = raw;
                meth_label_cpg(lab.data(), tlen, 0, 0, f.orig_bns, f.orig_pac);
                // Query: the target's prefix with conversions and noise.
                std::vector<uint8_t> q;
                for (int i = 0; i < qlen; ++i) {
                    uint8_t b = i < tlen ? raw[i] : (uint8_t)(rng() % 4);
                    const unsigned u = rng() % 100;
                    if (hyp && b == METH_C && u < 60) b = METH_T;
                    else if (!hyp && b == METH_G && u < 60) b = METH_A;
                    else if (u < 6) b = (uint8_t)(rng() % 5);
                    q.push_back(b);
                }
                const int h0 = 5 + (int)(rng() % 30), w = 5 + (int)(rng() % 40);
                int r1[5], r2[5];
                const int s1 = k->scalarBandedSWA(qlen, q.data(), tlen, lab.data(), w, h0,
                                                  &r1[0], &r1[1], &r1[2], &r1[3], &r1[4]);
                const int s2 = ksw_extend2(qlen, q.data(), tlen, lab.data(), 7, mat,
                                           o->o_del, o->e_del, o->o_ins, o->e_ins, w,
                                           o->pen_clip5, o->zdrop, h0,
                                           &r2[0], &r2[1], &r2[2], &r2[3], &r2[4]);
                CAPTURE(sc); CAPTURE(hyp); CAPTURE(it);
                CHECK(s1 == s2);
                for (int j = 0; j < 5; ++j) CHECK(r1[j] == r2[j]);
                int u[5];
                if (k->scalarBandedSWA(qlen, q.data(), tlen, raw.data(), w, h0,
                                       &u[0], &u[1], &u[2], &u[3], &u[4]) != s1)
                    ++n_diff_unlabelled;
            }
            // The labels matter: plain codes score the CpG conversions as mismatches.
            CHECK(n_diff_unlabelled > 0);
        }
        free(o);
    }
}

TEST_CASE("kernel matrix caches read a 7-code matrix at stride 7"
          * doctest::test_suite("unit/meth_model"))
{
    // Rows/cols 0..3 of mat7 are the plain symmetric matrix (the freed cells live
    // in rows 5/6), so at stride 7 the SIMD-side detectors see a symmetric matrix.
    mem_opt_t *o = opt_for(MEM_METH_SCORING_GENOMIC, 1, 4);
    const int8_t a = (int8_t)o->a, mm = (int8_t)-o->b;
    for (const int8_t *mat : { o->mat7_ot, o->mat7_ob }) {
        CHECK_FALSE(bsw_generic_matrix(mat, a, mm, 7));
        CHECK_FALSE(bsw_freed_cell(mat, a, mm, false, 7).single);
        CHECK_FALSE(bsw_freed_pair(mat, a, mm, 7).supported);
    }
    // ...while the 5-code OT matrix is the single freed cell, as before.
    CHECK(bsw_freed_cell(o->mat_ot, a, mm, false).rank1);
    // Constructing a kernel object on a 49-entry matrix is safe on every tier.
    auto k = make_banded_pair_wise_sw(o->o_del, o->e_del, o->o_ins, o->e_ins, o->zdrop,
                                      o->pen_clip5, o->mat7_ot, o->a, o->b, 1, 7);
    CHECK(k != nullptr);
    free(o);
}

TEST_CASE("meth_model_mat and the u8 rescue admission under cpg"
          * doctest::test_suite("unit/meth_model"))
{
    mem_opt_t *o = opt_for(MEM_METH_SCORING_NEUTRAL, 1, 4);
    o->meth_mode = 1;
    int m = 0;
    o->meth_context = METH_CTX_ALL;
    CHECK(meth_model_mat(o, 1, &m) == o->mat_ot);
    CHECK(m == 5);
    CHECK(meth_model_mat(o, -1, &m) == o->mat);
    CHECK(m == 5);
    const int shift_all = matesw_u8_shift_for(o, o->mat), qmax_all = matesw_u8_qmax_for(o, o->mat);
    o->meth_context = METH_CTX_CPG;
    CHECK(meth_model_mat(o, 1, &m) == o->mat7_ot);
    CHECK(m == 7);
    CHECK(meth_model_mat(o, 0, &m) == o->mat7_ob);
    CHECK(m == 7);
    CHECK(meth_model_mat(o, -1, &m) == o->mat);   // no hypothesis: plain codes
    CHECK(m == 5);
    CHECK(meth_mat_m(o, o->mat7_ob) == 7);
    CHECK(meth_mat_m(o, o->mat_ob) == 5);
    // mat7 holds the same values as mat/mat_ot/mat_ob, so the 8-bit admission is unchanged.
    CHECK(matesw_u8_shift(o->mat7_ot, 7) == matesw_u8_shift(o->mat_ot));
    CHECK(matesw_u8_qmax(o->mat7_ot, 7) == matesw_u8_qmax(o->mat_ot));
    CHECK(matesw_u8_shift_for(o, o->mat) == shift_all);
    CHECK(matesw_u8_qmax_for(o, o->mat) == qmax_all);
    o->meth_mode = 0;
    CHECK(meth_model_mat(o, 1, &m) == o->mat);   // cpg is a --meth option
    CHECK(m == 5);
    free(o);
}
