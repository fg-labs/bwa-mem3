// test/unit/test_xa_regen_reuse.cpp — mem_gen_alt's reuse of duplicate XA hits.
//
// 1. bns_pac_fwd_equal compares two forward-strand windows of the 2-bit pac; it
//    keys mem_gen_alt's reuse of an earlier hit's global alignment, so a false
//    "equal" would copy the wrong CIGAR/NM. It is checked against a plain
//    string compare, including windows at and near the end of the pac and
//    windows that span a contig boundary, on a pac that ends right before an
//    unreadable page so that any read past its last byte faults.
// 2. mem_reg2aln with `reuse` takes the given CIGAR ops and NM verbatim and still
//    squeezes end deletions, shifts the position and adds the clips; with `made`
//    it hands out the ops before that step. A made/reuse round trip equals a
//    plain mem_reg2aln.
// 3. mem_gen_alt reuses the CIGAR ops and NM of an earlier hit with an equal key.
//    Its XA string must be byte-identical to one written out from a fresh
//    mem_reg2aln of every hit. Besides exact, reverse, one-base-off and gapped
//    duplicates, the hits include pairs that agree on every key component but
//    one and regenerate differently (strand, query window, true score, band,
//    --meth hypothesis), so a key that drops a component fails; and clipped hits
//    whose global alignment ends in a deletion, reused and not.

#include "doctest/doctest.h"

#include "bwa.h"
#include "bwamem.h"
#include "bntseq.h"
#include "meth_opts.h"
#include "meth_orig_ref_fixture.h"

#include <sys/mman.h>
#include <unistd.h>

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <vector>

#ifndef MAP_ANONYMOUS
#define MAP_ANONYMOUS MAP_ANON
#endif

using meth_test::flip;
using meth_test::nt4;
using meth_test::random_bases;
using meth_test::revcomp;

namespace {

// The 2-bit pac of `s` in exactly as many bytes as the bases need, the last of
// them immediately followed by a PROT_NONE page.
struct GuardedPac {
    uint8_t *map = nullptr, *pac = nullptr;
    size_t map_len = 0;
    explicit GuardedPac(const std::string &s) {
        const size_t n = (s.size() + 3) / 4, page = (size_t)sysconf(_SC_PAGESIZE);
        const size_t body = (n + page - 1) / page * page + page;  // at least one readable page before
        map_len = body + page;
        void *m = mmap(nullptr, map_len, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
        REQUIRE(m != MAP_FAILED);
        map = (uint8_t *)m;
        REQUIRE(mprotect(map + body, page, PROT_NONE) == 0);
        pac = map + body - n;
        meth_test::pack_into(s, pac);
    }
    ~GuardedPac() { munmap(map, map_len); }
    GuardedPac(const GuardedPac &) = delete;
    GuardedPac &operator=(const GuardedPac &) = delete;
};

void check_pac_fwd_equal(const int len1, const int len2) {
    // two "contigs" back to back: the pac does not know about the boundary
    const int L = len1 + len2;
    std::string ref = random_bases(L, 7 + L);
    std::mt19937 rng(11);
    // plant equal windows (some spanning the boundary at len1, some ending at the
    // last base) so that "equal" is exercised as often as "different"
    struct Plant { int a, b, n; };
    std::vector<Plant> plants = {{10, 300, 150},      {len1 - 60, 100, 130}, {37, L - 123, 123},
                                 {L - 97, 503, 97},   {200, 201, 64},        {5, len1 + 1, 33},
                                 {650, L - 200, 200}, {300, L - 70, 70},     {400, L - 33, 33}};
    for (const Plant &p : plants) ref.replace(p.b, p.n, ref.substr(p.a, p.n));
    const GuardedPac pac(ref);
    auto oracle = [&](int64_t a, int64_t b, int64_t n) { return ref.compare(a, n, ref, b, n) == 0; };
    int n_eq = 0, n_ne = 0, n_slow = 0, n_fast_edge = 0;
    auto check = [&](int64_t a, int64_t b, int64_t n) {
        CAPTURE(a); CAPTURE(b); CAPTURE(n);
        const bool want = oracle(a, b, n);
        CHECK((bns_pac_fwd_equal(L, pac.pac, a, b, n) != 0) == want);
        CHECK((bns_pac_fwd_equal(L, pac.pac, b, a, n) != 0) == want);
        (want ? n_eq : n_ne)++;
        const int64_t hi = a > b ? a : b;
        if (hi + n + 64 >= L) ++n_slow;
        else if (hi + n + 64 >= L - 4) ++n_fast_edge;
    };
    for (const Plant &p : plants) {
        check(p.a, p.b, p.n);
        // every prefix and suffix length crosses the 32-base word and the tail mask
        for (int n = 0; n <= p.n; n += 1 + (n > 40) * 7) {
            check(p.a, p.b, n);
            check(p.a + p.n - n, p.b + p.n - n, n);
        }
        // one differing base at the first, last and every 13th position, alone
        for (int k = 0; k < p.n; k += (k == 0 ? 1 : 13)) {
            std::string saved = ref;
            ref[p.b + k] = flip(ref[p.b + k]);
            const GuardedPac pac2(ref);
            CHECK(bns_pac_fwd_equal(L, pac2.pac, p.a, p.b, p.n) == 0);
            ref = saved;
        }
        std::string saved = ref;
        ref[p.b + p.n - 1] = flip(ref[p.b + p.n - 1]);
        const GuardedPac pac3(ref);
        CHECK(bns_pac_fwd_equal(L, pac3.pac, p.a, p.b, p.n) == 0);
        ref = saved;
    }
    // every window that ends within the last 100 bases, and the windows on either
    // side of the fast path's bound (hi + len + 64 < l_pac), at every phase
    for (int64_t n = 0; n <= 140; ++n)
        for (int64_t end = L - 100; end <= L; ++end) {
            if (end - n < 0) continue;
            check(end - n, (end - n) * 7 % (L - n + 1), n);
        }
    for (int64_t n = 0; n <= 100; ++n)
        for (int64_t d = -4; d <= 4; ++d) {
            const int64_t b = L - 64 - n + d;
            if (b < 0 || b + n > L) continue;
            for (int64_t ph = 0; ph < 4; ++ph) check(ph + 3, b, n);
        }
    // random windows at every offset phase, many ending exactly at the last base
    for (int it = 0; it < 4000; ++it) {
        const int64_t n = rng() % 140;
        int64_t a = rng() % (L - n + 1), b = rng() % (L - n + 1);
        if (it % 5 == 0) b = L - n;  // ends at the last base of the pac
        if (it % 7 == 0) a = b ^ (rng() & 3);  // same byte, other phase
        if (a + n > L) a = L - n;
        check(a, b, n);
    }
    CHECK(n_eq > 100);
    CHECK(n_ne > 1000);
    CHECK(n_slow > 100);
    CHECK(n_fast_edge > 100);
}

}  // namespace

TEST_CASE("bns_pac_fwd_equal agrees with a string compare" * doctest::test_suite("unit/xa_regen")) {
    SUBCASE("full last byte") { check_pac_fwd_equal(700, 520); }
    SUBCASE("partial last byte") { check_pac_fwd_equal(701, 522); }
}

namespace {

struct Reads {
    meth_test::OrigRefFixture fx;
    std::vector<char> seq;  // nt4
    Reads(const std::string &ref, const std::string &read) : fx(ref) {
        for (char b : read) seq.push_back((char)nt4(b));
    }
};

// A region over forward bases [fb, fe) of a kRef-long reference, on the strand `rev`.
mem_alnreg_t region(int64_t kRef, int64_t fb, int64_t fe, bool rev, int qb, int qe, int score, int truesc, int w,
                    int secondary, int hyp) {
    mem_alnreg_t r{};
    r.rb = rev ? 2 * kRef - fe : fb;
    r.re = rev ? 2 * kRef - fb : fe;
    r.qb = qb;
    r.qe = qe;
    r.rid = 0;
    r.score = score;
    r.truesc = truesc;
    r.secondary = r.secondary_all = secondary;
    r.w = w;
    r.seedcov = r.seedlen0 = qe - qb;
    r.meth_hypothesis = hyp;
    return r;
}

std::string cigar_str(const mem_aln_t &t) {
    std::string s;
    for (int k = 0; k < t.n_cigar; ++k) s += std::to_string(t.cigar[k] >> 4) + "MIDSHN"[t.cigar[k] & 0xf];
    return s;
}

// "X,<strand><pos>,<CIGAR>,<NM>;" from a fresh want_md 0 regeneration: the XA
// entry format, written independently of mem_gen_alt.
std::string fresh_entry(const mem_opt_t *opt, const Reads &rd, const mem_alnreg_t &r, const char *meth_orig) {
    mem_aln_t t = mem_reg2aln(opt, rd.fx.orig_bns, rd.fx.orig_pac, (int)rd.seq.size(), rd.seq.data(), &r, meth_orig, 0);
    std::string s = std::string(rd.fx.orig_bns->anns[t.rid].name) + "," + (t.is_rev ? "-" : "+") +
                    std::to_string(t.pos + 1) + "," + cigar_str(t) + "," + std::to_string(t.NM) + ";";
    free(t.cigar);
    return s;
}

// Every field of a mem_aln_t that mem_reg2aln sets, the CIGAR ops included.
void check_same_aln(const mem_aln_t &x, const mem_aln_t &y) {
    CHECK(x.rid == y.rid);
    CHECK(x.pos == y.pos);
    CHECK(x.is_rev == y.is_rev);
    CHECK(x.flag == y.flag);
    CHECK(x.mapq == y.mapq);
    CHECK(x.NM == y.NM);
    CHECK(x.score == y.score);
    CHECK(x.sub == y.sub);
    CHECK(x.is_alt == y.is_alt);
    CHECK(x.alt_sc == y.alt_sc);
    CHECK(x.HN == y.HN);
    CHECK(x.meth_hypothesis == y.meth_hypothesis);
    CHECK(cigar_str(x) == cigar_str(y));
}

// mem_gen_alt over `regs` (regs[0] the primary, the rest its hits) must give the
// XA string written out from a fresh regeneration of every hit, reusing exactly
// `want_reused` of them; each `distinct` pair of hits must regenerate differently
// (that is what makes the pair discriminate a key component).
void check_gen_alt(const mem_opt_t *opt, const Reads &rd, const std::vector<mem_alnreg_t> &regs, const char *meth_orig,
                   long want_reused, const std::vector<std::pair<int, int>> &distinct) {
    mem_alnreg_v a{regs.size(), regs.size(), const_cast<mem_alnreg_t *>(regs.data())};
    std::vector<std::string> fresh(regs.size());
    std::string want;
    for (size_t k = 1; k < regs.size(); ++k) want += fresh[k] = fresh_entry(opt, rd, regs[k], meth_orig);
    for (const auto &d : distinct) {
        CAPTURE(d.first); CAPTURE(d.second);
        CHECK(fresh[d.first] != fresh[d.second]);
    }
    int *hn = nullptr;
    const long reused0 = mem_gen_alt_n_reused;
    char **xa = mem_gen_alt(opt, rd.fx.orig_bns, rd.fx.orig_pac, &a, (int)rd.seq.size(), rd.seq.data(), &hn, meth_orig);
    CHECK(mem_gen_alt_n_reused - reused0 == want_reused);
    REQUIRE(xa != nullptr);
    REQUIRE(xa[0] != nullptr);
    CHECK(std::string(xa[0]) == want);
    CHECK(hn[0] == (int)regs.size() - 1);
    for (size_t k = 0; k < regs.size(); ++k) free(xa[k]);
    free(xa);
    free(hn);
}

// A reference with copies of a 150-base segment and a read that is the segment
// with four substitutions (so gapped copies take the global DP).
struct Dups {
    static constexpr int kRef = 6000, seg_len = 150, seg_at = 200;
    std::string ref, seg, read;
    Dups() {
        ref = random_bases(kRef, 2026);
        seg = ref.substr(seg_at, seg_len);
        read = seg;
        for (int k : {20, 55, 100, 130}) read[k] = flip(read[k]);
    }
};

void check_duplicates(mem_opt_t *meth) {
    Dups d;
    const int kRef = Dups::kRef, seg_len = Dups::seg_len;
    const std::string &seg = d.seg;
    // gapped copy: two bases inserted in the reference at 70
    const std::string gapped = seg.substr(0, 70) + "GT" + seg.substr(70);
    std::string onediff = seg;
    onediff[90] = flip(onediff[90]);
    // equal length, but the read takes an insertion at 40 and a deletion at 110:
    // a zero band gives an ungapped CIGAR, a wide one the gapped alignment
    const std::string shifted = seg.substr(0, 40) + seg.substr(41, 70) + "A" + seg.substr(111);
    REQUIRE(shifted.size() == (size_t)seg_len);
    struct Copy { int at; std::string bases; bool rev; };
    const std::vector<Copy> copies = {
        {800, seg, false},          {1300, seg, false},           // exact forward duplicates
        {1800, seg, true},          {2300, seg, true},            // exact reverse duplicates
        {2800, onediff, false},                                   // equal key, one base differs
        {3300, gapped, false},      {3800, gapped, false},        // gapped duplicates
        {4300, gapped, true},       {4800, gapped, true},
        {5200, shifted, false},     {5500, shifted, false},       // equal-length gapped duplicates
    };
    for (const Copy &c : copies) d.ref.replace(c.at, c.bases.size(), c.rev ? revcomp(c.bases) : c.bases);
    Reads rd(d.ref, d.read);
    mem_opt_t *plain = mem_opt_init();
    mem_opt_t *opt = meth ? meth : plain;
    opt->max_XA_hits = opt->max_XA_hits_alt = 50;
    const char *meth_orig = meth ? rd.seq.data() : NULL;

    const int sc = 150 - 4 * 5;  // the hits' score, kept in the XA drop window
    std::vector<mem_alnreg_t> regs;
    auto add = [&](int64_t fb, int64_t fe, bool rev, int truesc, int qb = 0, int qe = Dups::seg_len, int w = -1) {
        const int64_t rb = rev ? 2 * (int64_t)kRef - fe : fb;
        regs.push_back(region(kRef, fb, fe, rev, qb, qe, sc, truesc, w < 0 ? opt->w : w, regs.empty() ? -1 : 0,
                              meth ? (rb < (int64_t)kRef ? 1 : 0) : -1));
        return (int)regs.size() - 1;
    };
    add(Dups::seg_at, Dups::seg_at + seg_len, false, sc);  // the primary
    std::vector<int> at;
    for (const Copy &c : copies) at.push_back(add(c.at, c.at + (int)c.bases.size(), c.rev, c.bases == gapped ? sc - 8 : sc));
    add(1300, 1300 + seg_len, false, sc - 1);  // same bases, other truesc: its own key
    // strand: the forward strand of the locus whose reverse strand is a duplicate
    const int fwd1800 = add(1800, 1800 + seg_len, false, sc);
    // query window: the same reference bases under a window clipped at either end;
    // the global alignment then starts (qb) or ends (qe) in a 10-base deletion
    const int qb10 = add(800, 800 + seg_len, false, sc, 10, seg_len);
    add(1300, 1300 + seg_len, false, sc, 10, seg_len);  // reused, clipped and squeezed
    const int qe140 = add(800, 800 + seg_len, false, sc, 0, seg_len - 10);
    add(1300, 1300 + seg_len, false, sc, 0, seg_len - 10);  // reused
    add(1800, 1800 + seg_len, true, sc, 10, seg_len);  // clipped on the reverse strand
    add(2300, 2300 + seg_len, true, sc, 10, seg_len);  // reused
    // true score: a near-perfect one infers a zero band (ungapped) on `shifted`, a
    // lower one a band wide enough for its two gaps
    const int ts_hi = add(5200, 5200 + seg_len, false, 140);
    const int ts_lo = add(5200, 5200 + seg_len, false, 120);
    add(5500, 5500 + seg_len, false, 140);  // reused
    // band: with a true score so low that the inferred band exceeds -w, the region's
    // own band w caps it: 0 keeps it ungapped, 100 lets the gaps in
    const int w0 = add(5200, 5200 + seg_len, false, 40, 0, seg_len, 0);
    const int w100 = add(5200, 5200 + seg_len, false, 40, 0, seg_len, 100);
    add(5500, 5500 + seg_len, false, 40, 0, seg_len, 100);  // reused
    std::vector<std::pair<int, int>> distinct = {
        {at[2], fwd1800},  // strand
        {at[0], qb10},     // qb
        {at[0], qe140},    // qe
        {ts_hi, ts_lo},    // truesc
        {w0, w100},        // w
    };
    // reused: the copies at 1300, 2300 (reverse), 3800, 4800 (gapped) and 5500, the two
    // clipped windows at 1300, the clipped reverse at 2300, and 5500 under truesc 140
    // and under w 100
    check_gen_alt(opt, rd, regs, meth_orig, 10, distinct);

    // the gapped and clipped hits did reach the DP, the deletion squeeze and the clips
    {
        bool saw_gap = false, saw_rev = false, saw_clip = false;
        for (size_t k = 1; k < regs.size(); ++k) {
            const std::string e = fresh_entry(opt, rd, regs[k], meth_orig);
            saw_gap |= e.find('I') != std::string::npos || e.find('D') != std::string::npos;
            saw_rev |= e.find(",-") != std::string::npos;
            saw_clip |= e.find('S') != std::string::npos;
        }
        CHECK(saw_gap);
        CHECK(saw_rev);
        CHECK(saw_clip);
    }
    free(plain);
}

// --meth: the hypothesis picks the matrix, so a read with C->T conversions
// regenerates with a different NM under hypotheses 0 and 1 over the same bases.
void check_meth_hypothesis(mem_opt_t *opt) {
    Dups d;
    const int kRef = Dups::kRef, seg_len = Dups::seg_len;
    for (int at : {1000, 2000}) d.ref.replace(at, seg_len, d.seg);
    std::string read = d.seg;
    int n_conv = 0;
    for (int k = 0; k < seg_len && n_conv < 8; ++k)
        if (read[k] == 'C') read[k] = 'T', ++n_conv;
    REQUIRE(n_conv == 8);
    Reads rd(d.ref, read);
    opt->max_XA_hits = opt->max_XA_hits_alt = 50;
    std::vector<mem_alnreg_t> regs;
    auto add = [&](int64_t fb, int hyp) {
        regs.push_back(region(kRef, fb, fb + seg_len, false, 0, seg_len, 140, 150, opt->w, regs.empty() ? -1 : 0, hyp));
        return (int)regs.size() - 1;
    };
    add(Dups::seg_at, 0);
    const int h0 = add(1000, 0), h1 = add(1000, 1);
    add(2000, 0);
    add(2000, 1);  // both reused
    check_gen_alt(opt, rd, regs, rd.seq.data(), 2, {{h0, h1}});
}

// --meth-context cpg: the CpG labels of a window read the reference base on
// either side of it, which the key does not hold. Two copies of a segment that
// ends in C, one followed by G (a CpG) and one by A, regenerate a read with a T
// there differently: a conversion at the CpG, a C/T mismatch at the other.
void check_cpg_flank(mem_opt_t *opt) {
    Dups d;
    const int kRef = Dups::kRef, seg_len = Dups::seg_len;
    std::string seg = d.seg;
    seg[seg_len - 1] = 'C';
    d.ref.replace(Dups::seg_at, seg_len, seg);
    d.ref[Dups::seg_at + seg_len] = 'G';
    d.ref.replace(1000, seg_len, seg);
    d.ref[1000 + seg_len] = 'G';
    d.ref.replace(2000, seg_len, seg);
    d.ref[2000 + seg_len] = 'A';
    // The same on the reverse strand: forward copies of revcomp(seg), whose window
    // (the reverse-complement frame) ends in that C and continues into the forward
    // base before the copy, complemented: C there makes the CpG, T does not.
    const std::string rc = revcomp(seg);
    d.ref.replace(3000, seg_len, rc);
    d.ref[3000 - 1] = 'C';
    d.ref.replace(4000, seg_len, rc);
    d.ref[4000 - 1] = 'T';
    std::string read = seg;
    read[seg_len - 1] = 'T';
    Reads rd(d.ref, read);
    opt->max_XA_hits = opt->max_XA_hits_alt = 50;
    std::vector<mem_alnreg_t> regs;
    auto add = [&](int64_t fb, bool rev) {
        // the hypothesis is the genome strand's: OT (1) on the forward half, OB (0) on the reverse
        regs.push_back(region(kRef, fb, fb + seg_len, rev, 0, seg_len, 140, 150, opt->w, regs.empty() ? -1 : 0,
                              rev ? 0 : 1));
        return (int)regs.size() - 1;
    };
    add(Dups::seg_at, false);
    const int cpg = add(1000, false), non_cpg = add(2000, false);
    const int rev_cpg = add(3000, true), rev_non_cpg = add(4000, true);
    for (const auto &pr : {std::make_pair(cpg, non_cpg), std::make_pair(rev_cpg, rev_non_cpg)}) {
        mem_core_key_t k0, k1;
        REQUIRE(mem_core_key(rd.fx.orig_bns, &regs[pr.first], &k0));
        REQUIRE(mem_core_key(rd.fx.orig_bns, &regs[pr.second], &k1));
        REQUIRE(mem_core_key_equal(rd.fx.orig_bns, rd.fx.orig_pac, &k0, &k1));  // equal key, yet ...
    }
    // ... distinct: no reuse, on either strand
    check_gen_alt(opt, rd, regs, rd.seq.data(), 0, {{cpg, non_cpg}, {rev_cpg, rev_non_cpg}});
}

}  // namespace

TEST_CASE("mem_gen_alt: a reused duplicate hit gives the XA entry of a fresh regeneration" *
          doctest::test_suite("unit/xa_regen")) {
    SUBCASE("without --meth") { check_duplicates(nullptr); }
    SUBCASE("--meth") {
        mem_opt_t *meth = mem_opt_init();
        meth->meth_mode = 1;
        meth->meth_scoring = MEM_METH_SCORING_GENOMIC;
        mem_opt_fill_meth_mat(meth);
        check_duplicates(meth);
        check_meth_hypothesis(meth);
        free(meth);
    }
    SUBCASE("--meth-context cpg") {
        mem_opt_t *meth = mem_opt_init();
        meth->meth_mode = 1;
        meth->meth_scoring = MEM_METH_SCORING_GENOMIC;
        meth->meth_context = METH_CTX_CPG;
        mem_opt_fill_meth_mat(meth);
        check_cpg_flank(meth);
        free(meth);
    }
}

TEST_CASE("mem_reg2aln: `made` hands out the global alignment and `reuse` takes it verbatim" *
          doctest::test_suite("unit/xa_regen")) {
    Dups d;
    const int kRef = Dups::kRef, seg_len = Dups::seg_len;
    for (int at : {800, 1300}) d.ref.replace(at, seg_len, d.seg);
    for (int at : {1800, 2300}) d.ref.replace(at, seg_len, revcomp(d.seg));
    Reads rd(d.ref, d.read);
    mem_opt_t *opt = mem_opt_init();
    const bntseq_t *bns = rd.fx.orig_bns;
    const uint8_t *pac = rd.fx.orig_pac;
    const int l = (int)rd.seq.size();
    const char *q = rd.seq.data();

    SUBCASE("a made/reuse round trip equals a plain regeneration") {
        // clipped at either end and on either strand: the global alignment starts or
        // ends in a 10-base deletion that the squeeze removes
        struct Win { int qb, qe; bool rev; };
        for (const Win &win : {Win{10, seg_len, false}, Win{0, seg_len - 10, false}, Win{10, seg_len, true},
                               Win{0, seg_len - 10, true}, Win{0, seg_len, false}}) {
            CAPTURE(win.qb); CAPTURE(win.qe); CAPTURE(win.rev);
            const int at = win.rev ? 1800 : 800;  // a copy of the read on its strand, and a duplicate of it
            const mem_alnreg_t miss = region(kRef, at, at + seg_len, win.rev, win.qb, win.qe, 130, 130, opt->w, 0, -1);
            const mem_alnreg_t dup = region(kRef, at + 500, at + 500 + seg_len, win.rev, win.qb, win.qe, 130, 130, opt->w, 0, -1);
            mem_core_key_t k0, k1;
            REQUIRE(mem_core_key(bns, &miss, &k0));
            REQUIRE(mem_core_key(bns, &dup, &k1));
            CHECK(mem_core_key_equal(bns, pac, &k0, &k1));
            std::vector<uint32_t> ops = {0xdeadbeef};  // `made` appends after what is there
            int nm = -1;
            mem_aln_t made = mem_reg2aln(opt, bns, pac, l, q, &miss, NULL, 0, NULL, &ops, &nm);
            mem_aln_t plain = mem_reg2aln(opt, bns, pac, l, q, &miss, NULL, 0);
            check_same_aln(made, plain);
            REQUIRE(ops.size() >= 2);
            CHECK(ops[0] == 0xdeadbeef);
            CHECK(nm == (int)plain.NM);
            // the ops are those before the squeeze: a clipped window's start in a deletion
            const bool end_d = (ops[1] & 0xf) == 2 || (ops.back() & 0xf) == 2;
            CHECK(end_d == (win.qb != 0 || win.qe != seg_len));
            const mem_cigar_core_t core = {ops.data() + 1, (int)ops.size() - 1, nm};
            mem_aln_t reused = mem_reg2aln(opt, bns, pac, l, q, &dup, NULL, 0, &core);
            mem_aln_t fresh = mem_reg2aln(opt, bns, pac, l, q, &dup, NULL, 0);
            check_same_aln(reused, fresh);
            CHECK(reused.pos == plain.pos + 500);
            free(made.cigar); free(plain.cigar); free(reused.cigar); free(fresh.cigar);
        }
    }
    SUBCASE("a reused core is taken verbatim, then squeezed and clipped") {
        // a core no DP would produce for these bases, so only `reuse` can explain it
        const uint32_t lead_d[] = {3u << 4 | 2, 140u << 4 | 0};
        const uint32_t trail_d[] = {140u << 4 | 0, 2u << 4 | 2};
        const mem_alnreg_t fwd = region(kRef, 1300, 1300 + seg_len, false, 10, seg_len, 130, 130, opt->w, 0, -1);
        const mem_alnreg_t rev = region(kRef, 1300, 1300 + seg_len, true, 10, seg_len, 130, 130, opt->w, 0, -1);
        const mem_cigar_core_t c_lead = {lead_d, 2, 77}, c_trail = {trail_d, 2, 66};
        mem_aln_t t = mem_reg2aln(opt, bns, pac, l, q, &fwd, NULL, 0, &c_lead);
        CHECK((int)t.NM == 77);
        CHECK((int)t.is_rev == 0);
        CHECK(t.pos == 1300 + 3);  // the leading deletion moves the start
        CHECK(cigar_str(t) == "10S140M");
        free(t.cigar);
        t = mem_reg2aln(opt, bns, pac, l, q, &fwd, NULL, 0, &c_trail);
        CHECK((int)t.NM == 66);
        CHECK(t.pos == 1300);
        CHECK(cigar_str(t) == "10S140M");
        free(t.cigar);
        t = mem_reg2aln(opt, bns, pac, l, q, &rev, NULL, 0, &c_lead);
        CHECK((int)t.NM == 77);
        CHECK((int)t.is_rev == 1);
        CHECK(t.pos == 1300 + 3);
        CHECK(cigar_str(t) == "140M10S");  // a reverse hit's query start is clipped at its 3' end
        free(t.cigar);
    }
    free(opt);
}

TEST_CASE("mem_core_key: no key for a window off one strand" * doctest::test_suite("unit/xa_regen")) {
    // mem_core_key reads only bns->l_pac; its rejects are what keep
    // bns_pac_fwd_equal inside [0, l_pac)
    const int64_t kRef = 6000;
    bntseq_t bns{};
    bns.l_pac = kRef;
    mem_core_key_t k;
    mem_alnreg_t r = region(kRef, 100, 250, false, 0, 150, 130, 130, 100, -1, -1);
    REQUIRE(mem_core_key(&bns, &r, &k));  // control: a forward window
    CHECK(k.fb == 100);
    CHECK(k.rlen == 150);
    r = region(kRef, 100, 250, true, 0, 150, 130, 130, 100, -1, -1);
    REQUIRE(mem_core_key(&bns, &r, &k));  // control: its reverse strand, on forward coordinates
    CHECK(k.rev == 1);
    CHECK(k.fb == 100);
    r.rb = kRef - 50; r.re = kRef + 100;  // bridges the forward/reverse boundary
    CHECK(!mem_core_key(&bns, &r, &k));
    r.rb = 2 * kRef - 100; r.re = 2 * kRef + 1;  // runs past the reverse strand
    CHECK(!mem_core_key(&bns, &r, &k));
    r.rb = 300; r.re = 300;  // empty
    CHECK(!mem_core_key(&bns, &r, &k));
    r.rb = 300; r.re = 200;  // inverted
    CHECK(!mem_core_key(&bns, &r, &k));
    r.rb = -1; r.re = 100;  // unmapped
    CHECK(!mem_core_key(&bns, &r, &k));
}
