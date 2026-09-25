// test/unit/test_reg2sam_anchored.cpp — mem_reg2sam_anchored must emit exactly
// what mem_reg2sam emits.
//
// On the no-pairing path, mem_sam_pe_batch_post converts each end's anchor
// region with mem_reg2aln for the mate fields, then hands that conversion to
// mem_reg2sam_anchored so the region is not converted (and globally aligned)
// a second time. That is only a performance change if the reused conversion is
// byte-for-byte what a fresh mem_reg2aln would produce, so these tests emit
// reads with and without each region anchored and require identical SAM:
// a chimeric read (forward primary plus a clipped reverse-strand
// supplementary), a gapped read whose CIGAR and MD come from the global
// alignment the change exists to skip, and an ALT hit. They also pin the
// ownership contract (the caller's anchor, which the pair path still reads as
// the other end's mate, is left untouched) and the guard: an anchor that is
// not the conversion of its region, such as the mate's or a tied region's,
// must be a fatal error instead of being emitted.

#include "doctest/doctest.h"

#include "bwamem.h"
#include "meth_orig_ref_fixture.h"

#include <algorithm>
#include <csignal>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

namespace {

// A deterministic, low-repeat reference: a 32-bit LCG over ACGT.
std::string make_reference(int len) {
    std::string ref(len, 'A');
    uint32_t x = 12345;
    for (int i = 0; i < len; ++i) {
        x = x * 1103515245u + 12345u;
        ref[i] = "ACGT"[(x >> 16) & 3];
    }
    return ref;
}

int nt4(char c) { return meth_test::base_to_2bit(c); }

constexpr int kRefLen = 2000;

// A read built from reference pieces; `seq` is nt4, the form mem_reg2aln and
// mem_aln2sam expect.
struct Read {
    std::string       name{"read1"};
    std::vector<char> seq;
    std::string       qual;
    bseq1_t           s{};

    // Forward-strand bases ref[at, at + len).
    void add_forward(const std::string &ref, int at, int len) {
        for (int i = 0; i < len; ++i) seq.push_back((char)nt4(ref[at + i]));
    }
    // The reverse complement of ref[at, at + len).
    void add_reverse(const std::string &ref, int at, int len) {
        for (int i = 0; i < len; ++i) seq.push_back((char)(3 - nt4(ref[at + len - 1 - i])));
    }
    void finish() {
        qual.assign(seq.size(), 'I');
        s.l_seq = (int)seq.size();
        s.name  = const_cast<char *>(name.c_str());
        s.seq   = seq.data();
        s.qual  = const_cast<char *>(qual.c_str());
    }
};

// A region aligning read [qb, qe) to doubled-strand reference [rb, re) with
// alignment score `score`.
mem_alnreg_t make_region(const mem_opt_t *opt, int64_t rb, int64_t re, int qb, int qe,
                         int score) {
    mem_alnreg_t r{};
    r.rb = rb;
    r.re = re;
    r.qb = qb;
    r.qe = qe;
    r.rid = 0;
    r.score = r.truesc = score;
    r.sub = r.csub = 0;
    r.secondary = -1;
    r.w = opt->w;
    r.seedcov = r.seedlen0 = qe - qb;
    r.meth_hypothesis = -1;
    return r;
}

// Doubled-strand start of the reverse-strand hit on ref[at, at + len).
int64_t reverse_rb(int at, int len) { return 2 * (int64_t)kRefLen - (at + len); }

// SAM text from mem_reg2sam_anchored (a NULL anchor means mem_reg2sam).
std::string emit(const mem_opt_t *opt, const meth_test::OrigRefFixture &fx, Read &read,
                 mem_alnreg_v *a, int extra_flag, const mem_aln_t *mate, int anchor_k,
                 const mem_aln_t *anchor) {
    if (anchor == NULL)
        mem_reg2sam(opt, fx.orig_bns, fx.orig_pac, &read.s, a, extra_flag, mate);
    else
        mem_reg2sam_anchored(opt, fx.orig_bns, fx.orig_pac, &read.s, a, extra_flag, mate,
                             anchor_k, anchor);
    REQUIRE(read.s.sam != NULL);
    std::string out(read.s.sam);
    free(read.s.sam);
    read.s.sam = NULL;
    return out;
}

mem_aln_t convert(const mem_opt_t *opt, const meth_test::OrigRefFixture &fx, const Read &read,
                  const mem_alnreg_t *region) {
    return mem_reg2aln(opt, fx.orig_bns, fx.orig_pac, read.s.l_seq, read.s.seq, region, NULL);
}

// The bytes of a mem_reg2aln CIGAR buffer: n_cigar ops, then the MD string.
std::string cigar_bytes(const mem_aln_t &aln) {
    const char *md = (const char *)(aln.cigar + aln.n_cigar);
    return std::string((const char *)aln.cigar, (size_t)aln.n_cigar * 4 + strlen(md) + 1);
}

/* Anchors each region of `a` in turn and requires the anchored SAM to equal
 * mem_reg2sam's, with the caller's anchor (record and CIGAR+MD buffer)
 * unchanged by the call. */
void check_each_anchor(const mem_opt_t *opt, const meth_test::OrigRefFixture &fx, Read &read,
                       mem_alnreg_v *a, int extra_flag, const mem_aln_t *mate) {
    const std::string expected = emit(opt, fx, read, a, extra_flag, mate, -1, NULL);
    for (int anchor_k = 0; anchor_k < (int)a->n; ++anchor_k) {
        CAPTURE(anchor_k);
        mem_aln_t anchor = convert(opt, fx, read, &a->a[anchor_k]);
        REQUIRE(anchor.cigar != NULL);
        const mem_aln_t before = anchor;
        const std::string before_cigar = cigar_bytes(anchor);

        CHECK(emit(opt, fx, read, a, extra_flag, mate, anchor_k, &anchor) == expected);

        CHECK(anchor.cigar == before.cigar);
        CHECK(anchor.n_cigar == before.n_cigar);
        CHECK(anchor.flag == before.flag);
        CHECK((int)anchor.mapq == (int)before.mapq);
        CHECK(anchor.pos == before.pos);
        CHECK(cigar_bytes(anchor) == before_cigar);
        free(anchor.cigar);
    }
}

// A mapped mate for the records' RNEXT/PNEXT/MC fields.
mem_aln_t make_mate(uint32_t *cigar) {
    mem_aln_t m{};
    m.pos = 900;
    m.rid = 0;
    m.is_rev = 1;
    m.mapq = 60;
    m.n_cigar = 1;
    m.cigar = cigar;
    m.HN = -1;
    m.meth_hypothesis = -1;
    return m;
}

/* Run `body` in a forked child with stderr and core dumps suppressed and
 * report how it terminated (as test_xassert_macro.cpp does, which explains
 * why the child restores SIGABRT's default disposition). A child that runs
 * to completion exits 0. */
int run_in_child(void (*body)()) {
    const pid_t pid = fork();
    if (pid == 0) {
        signal(SIGABRT, SIG_DFL);
        struct rlimit no_core = {0, 0};
        setrlimit(RLIMIT_CORE, &no_core);
        (void)freopen("/dev/null", "w", stderr);
        body();
        _exit(0);
    }
    REQUIRE(pid > 0);
    int status = 0;
    REQUIRE(waitpid(pid, &status, 0) == pid);
    return status;
}

/* The guard fails through err_fatal, which exits with EXIT_FAILURE; a signal
 * (a crash or a failed assert) is not the guard firing. */
bool exited_via_err_fatal(int status) {
    return WIFEXITED(status) && WEXITSTATUS(status) == EXIT_FAILURE;
}

/* Two equal-scoring forward hits for one read at different positions on the
 * same contig: the tie the guard's rid/score/is_alt checks alone could not
 * tell apart. Anchors region `anchor_k` with the conversion of region
 * `converted_k`. */
void emit_tied_regions(int anchor_k, int converted_k) {
    const std::string ref = make_reference(kRefLen);
    meth_test::OrigRefFixture fx(ref);
    mem_opt_t *opt = mem_opt_init();
    opt->flag |= MEM_F_ALL;  // keep the tied hit as a record rather than an XA entry
    Read read;
    read.add_forward(ref, 300, 60);
    read.finish();
    mem_alnreg_t regs[2] = {
        make_region(opt, 300, 360, 0, 60, 60),
        make_region(opt, 1300, 1360, 0, 60, 60),
    };
    mem_alnreg_v a{};
    a.n = a.m = 2;
    a.a = regs;
    mem_aln_t anchor = convert(opt, fx, read, &a.a[converted_k]);
    mem_reg2sam_anchored(opt, fx.orig_bns, fx.orig_pac, &read.s, &a, 0, NULL, anchor_k, &anchor);
    free(read.s.sam);
    free(anchor.cigar);
    free(opt);
}

void anchor_own_region() { emit_tied_regions(0, 0); }
void anchor_tied_region() { emit_tied_regions(0, 1); }

/* The pair path's mistake to guard against: end 0's region anchored with
 * end 1's conversion (the mate), both perfect forward hits of equal score. */
void anchor_the_mate() {
    const std::string ref = make_reference(kRefLen);
    meth_test::OrigRefFixture fx(ref);
    mem_opt_t *opt = mem_opt_init();
    Read end0, end1;
    end0.add_forward(ref, 300, 60);
    end0.finish();
    end1.add_forward(ref, 700, 60);
    end1.finish();
    mem_alnreg_t reg0 = make_region(opt, 300, 360, 0, 60, 60);
    mem_alnreg_t reg1 = make_region(opt, 700, 760, 0, 60, 60);
    mem_alnreg_v a0{};
    a0.n = a0.m = 1;
    a0.a = &reg0;
    mem_aln_t mate = convert(opt, fx, end1, &reg1);
    mem_reg2sam_anchored(opt, fx.orig_bns, fx.orig_pac, &end0.s, &a0, 0x41, &mate, 0, &mate);
    free(end0.s.sam);
    free(mate.cigar);
    free(opt);
}

}  // namespace

TEST_CASE("reg2sam_anchored: a chimeric read emits the same SAM with either region anchored"
          * doctest::test_suite("unit/pair")) {
    const std::string ref = make_reference(kRefLen);
    meth_test::OrigRefFixture fx(ref);
    mem_opt_t *opt = mem_opt_init();
    // [0, 70) is ref[300, 370) forward; [70, 110) is the reverse complement of
    // ref[1200, 1240). Each region soft-clips the other's bases, and the
    // supplementary is hard-clipped in the output.
    Read read;
    read.add_forward(ref, 300, 70);
    read.add_reverse(ref, 1200, 40);
    read.finish();
    mem_alnreg_t regs[2] = {
        make_region(opt, 300, 370, 0, 70, 70),
        make_region(opt, reverse_rb(1200, 40), reverse_rb(1200, 40) + 40, 70, 110, 40),
    };
    mem_alnreg_v a{};
    a.n = a.m = 2;
    a.a = regs;

    // A primary plus its supplementary: two records.
    const std::string sam = emit(opt, fx, read, &a, 0, NULL, -1, NULL);
    REQUIRE(std::count(sam.begin(), sam.end(), '\n') == 2);
    CHECK(sam.find("\tSA:Z:") != std::string::npos);

    uint32_t mate_cigar[2] = {(60u << 4) | 0, 0};
    mem_aln_t mate = make_mate(mate_cigar);
    check_each_anchor(opt, fx, read, &a, 0, NULL);
    check_each_anchor(opt, fx, read, &a, 0x41, &mate);

    free(opt);
}

TEST_CASE("reg2sam_anchored: a gapped region reuses the global alignment's CIGAR and MD"
          * doctest::test_suite("unit/pair")) {
    const std::string ref = make_reference(kRefLen);
    meth_test::OrigRefFixture fx(ref);
    mem_opt_t *opt = mem_opt_init();
    // ref[500, 540), a 2 bp deletion (ref[540, 542)), then ref[542, 582) with
    // one substitution: mem_reg2aln must run the banded global alignment, and
    // the CIGAR has three ops with a '^' deletion in MD.
    Read read;
    read.add_forward(ref, 500, 40);
    read.add_forward(ref, 542, 40);
    read.seq[60] = (char)((read.seq[60] + 1) & 3);
    read.finish();
    const int score = 79 * opt->a - opt->b - (opt->o_del + 2 * opt->e_del);
    mem_alnreg_t reg = make_region(opt, 500, 582, 0, 80, score);
    reg.sub = score - 10;  // a competing hit: MAPQ below 60
    mem_alnreg_v a{};
    a.n = a.m = 1;
    a.a = &reg;

    mem_aln_t probe = convert(opt, fx, read, &reg);
    REQUIRE(probe.n_cigar == 3);
    CHECK((probe.cigar[1] & 0xf) == 2);  // D
    CHECK(strchr((const char *)(probe.cigar + probe.n_cigar), '^') != NULL);
    CHECK((int)probe.mapq < 60);
    free(probe.cigar);

    uint32_t mate_cigar[2] = {(60u << 4) | 0, 0};
    mem_aln_t mate = make_mate(mate_cigar);
    check_each_anchor(opt, fx, read, &a, 0x81, &mate);

    free(opt);
}

TEST_CASE("reg2sam_anchored: an ALT hit emits the same SAM when anchored"
          * doctest::test_suite("unit/pair")) {
    const std::string ref = make_reference(kRefLen);
    meth_test::OrigRefFixture fx(ref);
    mem_opt_t *opt = mem_opt_init();
    Read read;
    read.add_forward(ref, 800, 60);
    read.finish();
    mem_alnreg_t reg = make_region(opt, 800, 860, 0, 60, 60);
    reg.is_alt = 1;
    reg.alt_sc = 55;
    mem_alnreg_v a{};
    a.n = a.m = 1;
    a.a = &reg;
    check_each_anchor(opt, fx, read, &a, 0, NULL);
    free(opt);
}

TEST_CASE("reg2sam_anchored: an anchor on a region that is not emitted changes nothing"
          * doctest::test_suite("unit/pair")) {
    const std::string ref = make_reference(kRefLen);
    meth_test::OrigRefFixture fx(ref);
    mem_opt_t *opt = mem_opt_init();
    Read read;
    read.add_forward(ref, 300, 70);
    read.add_reverse(ref, 1200, 40);
    read.finish();
    // Region 1 scores below -T, so the emit loop skips it.
    mem_alnreg_t regs[2] = {
        make_region(opt, 300, 370, 0, 70, 70),
        make_region(opt, reverse_rb(1200, 40), reverse_rb(1200, 40) + 40, 70, 110, opt->T - 1),
    };
    mem_alnreg_v a{};
    a.n = a.m = 2;
    a.a = regs;

    const std::string expected = emit(opt, fx, read, &a, 0, NULL, -1, NULL);
    REQUIRE(std::count(expected.begin(), expected.end(), '\n') == 1);
    mem_aln_t anchor = convert(opt, fx, read, &a.a[1]);
    CHECK(emit(opt, fx, read, &a, 0, NULL, 1, &anchor) == expected);
    free(anchor.cigar);
    free(opt);
}

TEST_CASE("reg2sam_anchored: a NULL anchor falls back to converting the region"
          * doctest::test_suite("unit/pair")) {
    const std::string ref = make_reference(kRefLen);
    meth_test::OrigRefFixture fx(ref);
    mem_opt_t *opt = mem_opt_init();
    Read read;
    read.add_forward(ref, 300, 70);
    read.finish();
    mem_alnreg_t reg = make_region(opt, 300, 370, 0, 70, 70);
    mem_alnreg_v a{};
    a.n = a.m = 1;
    a.a = &reg;

    const std::string expected = emit(opt, fx, read, &a, 0, NULL, -1, NULL);
    mem_reg2sam_anchored(opt, fx.orig_bns, fx.orig_pac, &read.s, &a, 0, NULL, 0, NULL);
    REQUIRE(read.s.sam != NULL);
    CHECK(std::string(read.s.sam) == expected);
    free(read.s.sam);
    read.s.sam = NULL;
    free(opt);
}

TEST_CASE("reg2sam_anchored: an anchor that is not its region's conversion is fatal"
          * doctest::test_suite("unit/pair")) {
    // Control: the harness itself runs to completion on a correct anchor.
    const int own = run_in_child(anchor_own_region);
    REQUIRE(WIFEXITED(own));
    CHECK(WEXITSTATUS(own) == 0);

    // Same contig, same score, same strand: only the position differs.
    CHECK(exited_via_err_fatal(run_in_child(anchor_tied_region)));
    // The pair path's own hazard: an end anchored with its mate's conversion.
    CHECK(exited_via_err_fatal(run_in_child(anchor_the_mate)));
}
