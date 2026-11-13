// test/unit/test_alt_lift.cpp — the ALT -> primary liftover behind --alt-project.
//
// alt_lift_load parses the .alt sidecar (bwakit SAM lines) against a bns,
// alt_lift_interval maps an interval of an ALT contig onto the primary assembly,
// and mem_alt_project rewrites a read's alignment regions with it. These cases
// pin the interval arithmetic (forward and reverse-complemented ALTs, deletions,
// insertions, clips, the band an interval needs), every loader skip rule, and
// the region rewrite: which regions are lifted, which stay on the ALT (including
// the band bwa_gen_cigar4 will give a lifted region, and the widest floor alt_bw
// can carry), the band floor a lifted region carries, which regions are one
// locus and which of them survive (a lifted hit replaces a primary hit only when
// it outscores it by more than -A times the primary hit's read bases it leaves
// out; one that cannot is tested only against the hits that survive; of two
// tied lifted hits, one that can replace beats one that cannot, and otherwise
// the lower ALT contig index wins), and the by-score output order. The expected
// values are worked out by hand from the CIGARs, not read back from the code.
// test/regression/alt_project.sh drives the option end to end, and
// test_gen_cigar_no_md.cpp pins bwa_gen_cigar_band's values.

#include "doctest/doctest.h"
#include "temp_index.h"
#include "temp_path.h"
#include "alt_lift.h"
#include "bntseq.h"
#include "bwamem.h"

#include <algorithm>
#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <random>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <unistd.h>
#include <vector>
#include <zlib.h>
#include <sys/stat.h>

namespace {
struct Ctg {
    const char *name;
    int len;
    int is_alt;
};

// A minimal bns: names, lengths, offsets and is_alt are all these paths read.
struct FakeBns {
    bntseq_t bns;
    std::vector<bntann1_t> anns;
    explicit FakeBns(const std::vector<Ctg> &ctgs) : anns(ctgs.size()) {
        memset(&bns, 0, sizeof(bns));
        int64_t off = 0;
        for (size_t i = 0; i < ctgs.size(); ++i) {
            memset(&anns[i], 0, sizeof(bntann1_t));
            anns[i].name = const_cast<char *>(ctgs[i].name);
            anns[i].len = ctgs[i].len;
            anns[i].is_alt = ctgs[i].is_alt;
            anns[i].offset = off;
            off += ctgs[i].len;
        }
        bns.n_seqs = (int32_t) ctgs.size();
        bns.anns = anns.data();
        bns.l_pac = off;
    }
};

struct TmpFile {
    std::string path;
    explicit TmpFile(const char *body) {
        char tmpl[4096];
        int fd = bwa_test_mkstemp(tmpl, sizeof(tmpl), "alt_lift_test");
        REQUIRE(fd >= 0);
        path = tmpl;
        FILE *fp = fdopen(fd, "w");
        REQUIRE(fp != nullptr);
        fputs(body, fp);
        fclose(fp);
    }
    ~TmpFile() { unlink(path.c_str()); }
};

struct LiftGuard {
    alt_lift_t *l;
    explicit LiftGuard(alt_lift_t *x) : l(x) {}
    ~LiftGuard() { alt_lift_destroy(l); }
};

// Compare lift data, never allocation addresses or struct padding.
void check_lift_tables(const alt_lift_t *actual, const alt_lift_t *expected) {
    REQUIRE(actual != nullptr);
    REQUIRE(expected != nullptr);
    REQUIRE(actual->n_seqs == expected->n_seqs);
    CHECK(actual->n_lift == expected->n_lift);
    for (int rid = 0; rid < actual->n_seqs; ++rid) {
        CAPTURE(rid);
        const auto &a = actual->a[rid], &b = expected->a[rid];
        REQUIRE(a.trid == b.trid);
        if (a.trid < 0) continue;
        CHECK(a.rev == b.rev);
        CHECK(a.tpos == b.tpos);
        CHECK(a.tlen == b.tlen);
        CHECK(a.qlen == b.qlen);
        REQUIRE(a.n_cigar == b.n_cigar);
        for (int i = 0; i < a.n_cigar; ++i) CHECK(a.cigar[i] == b.cigar[i]);
        for (int i = 0; i <= a.n_cigar; ++i) {
            CHECK(a.qoff[i] == b.qoff[i]);
            CHECK(a.roff[i] == b.roff[i]);
        }
    }
}

void write_gzip(const std::string &path, const std::string &text) {
    gzFile fp = gzopen(path.c_str(), "wb");
    REQUIRE(fp != nullptr);
    CHECK(gzwrite(fp, text.data(), (unsigned) text.size()) == (int) text.size());
    REQUIRE(gzclose(fp) == Z_OK);
}

struct OptGuard {
    mem_opt_t *o;
    OptGuard() : o(mem_opt_init()) {}
    ~OptGuard() { free(o); }
};

// chr1 (primary, 100 kb); alt_fwd and alt_rev: 1500 bases matching chr1:3000-4499
// (0-based), a 3-bp deletion, 497 matching; alt_rev is its reverse complement
// (FLAG 16). alt_clip: 100S200M at chr1:20001 (ALT-only first 100 bases).
// alt_bigdel: 100M500D100M at chr1:30001 (an ALT deletion wider than the band).
// alt_ins: 100M80I100M at chr1:50001 (an ALT insertion shorter than the band).
// alt_wdel: 100M100D100M at chr1:60001 (an ALT deletion of exactly the default
// band, -w 100). alt_bal: 100M5I30M5D100M at chr1:70001 (indels that cancel).
// alt_b35 and alt_b40: 100M{k}I30M{k}D100M at chr1:80001 and 81001, indels that
// cancel 30 bases apart with an excursion of 35 and 40. alt_maxdel and
// alt_hugedel: 100M1022D100M at chr1:90001 and 100M1100D100M at chr1:93001, ALT
// deletions of exactly MEM_ALT_BW_MAX and wider than it. alt_fwd2: a copy of
// alt_fwd, last in the index.
const std::vector<Ctg> kCtgs = {
    {"chr1", 100000, 0}, {"alt_fwd", 1997, 1}, {"alt_rev", 1997, 1},
    {"alt_clip", 300, 1}, {"alt_bigdel", 200, 1}, {"alt_ins", 280, 1},
    {"alt_wdel", 200, 1}, {"alt_bal", 235, 1}, {"alt_b35", 265, 1}, {"alt_b40", 270, 1},
    {"alt_maxdel", 200, 1}, {"alt_hugedel", 200, 1}, {"alt_fwd2", 1997, 1},
};
const char *kAlt =
    "@HD\tVN:1.6\n"
    "alt_fwd\t0\tchr1\t3001\t60\t1500M3D497M\t*\t0\t0\t*\t*\n"
    "alt_rev\t16\tchr1\t3001\t60\t1500M3D497M\t*\t0\t0\t*\t*\n"
    "alt_clip\t0\tchr1\t20001\t60\t100S200M\t*\t0\t0\t*\t*\n"
    "alt_bigdel\t0\tchr1\t30001\t60\t100M500D100M\t*\t0\t0\t*\t*\n"
    "alt_ins\t0\tchr1\t50001\t60\t100M80I100M\t*\t0\t0\t*\t*\n"
    "alt_wdel\t0\tchr1\t60001\t60\t100M100D100M\t*\t0\t0\t*\t*\n"
    "alt_bal\t0\tchr1\t70001\t60\t100M5I30M5D100M\t*\t0\t0\t*\t*\n"
    "alt_b35\t0\tchr1\t80001\t60\t100M35I30M35D100M\t*\t0\t0\t*\t*\n"
    "alt_b40\t0\tchr1\t81001\t60\t100M40I30M40D100M\t*\t0\t0\t*\t*\n"
    "alt_maxdel\t0\tchr1\t90001\t60\t100M1022D100M\t*\t0\t0\t*\t*\n"
    "alt_hugedel\t0\tchr1\t93001\t60\t100M1100D100M\t*\t0\t0\t*\t*\n"
    "alt_fwd2\t0\tchr1\t3001\t60\t1500M3D497M\t*\t0\t0\t*\t*\n";

enum { CHR1 = 0, ALT_FWD = 1, ALT_REV = 2, ALT_CLIP = 3, ALT_BIGDEL = 4, ALT_INS = 5, ALT_WDEL = 6,
       ALT_BAL = 7, ALT_B35 = 8, ALT_B40 = 9, ALT_MAXDEL = 10, ALT_HUGEDEL = 11, ALT_FWD2 = 12 };

// One alignment region: forward or reverse strand of [cs, ce) on contig rid, in
// bwa's doubled coordinate space. `tag` (stored in seedcov) identifies it after
// the rewrite.
mem_alnreg_t reg(const FakeBns &fb, int rid, int64_t cs, int64_t ce, int rev, int score, int tag,
                 int qb = 0, int qe = 150) {
    mem_alnreg_t r;
    memset(&r, 0, sizeof(r));
    const int64_t off = fb.bns.anns[rid].offset, l2 = fb.bns.l_pac << 1;
    r.rid = rid;
    r.rb = rev ? l2 - (off + ce) : off + cs;
    r.re = rev ? l2 - (off + cs) : off + ce;
    r.qb = qb; r.qe = qe;
    r.score = r.truesc = score;
    r.secondary = -1;
    r.is_alt = fb.bns.anns[rid].is_alt;
    r.seedcov = tag;
    return r;
}
// Resolve already-projected candidates without rescue, to isolate pairing/MAPQ.
struct PairOutcome { int paired, flag, q[2], score[2]; int64_t rb[2]; };
PairOutcome project_and_pair(mem_opt_t *opt, const FakeBns &fb,
                             const std::vector<mem_alnreg_t> input[2]) {
    std::vector<mem_alnreg_t> regions[2] = {input[0], input[1]};
    opt->flag |= MEM_F_NO_RESCUE;
    mem_alnreg_v a[2] = {};
    for (int i = 0; i < 2; ++i) {
        regions[i].resize(mem_alt_project(opt, &fb.bns, (int) regions[i].size(), regions[i].data()));
        a[i].a = regions[i].data();
        a[i].n = a[i].m = regions[i].size();
    }
    mem_pestat_t pes[4] = {};
    for (auto &p : pes) p.failed = 1;
    pes[1].failed = 0;
    pes[1].avg = 299; pes[1].std = 30; pes[1].low = 180; pes[1].high = 420;
    bseq1_t reads[2] = {};
    mem_cache *cache = (mem_cache *) calloc(1, sizeof(mem_cache));
    REQUIRE(cache != nullptr);
    int32_t gcnt = 0;
    int n_pri[2], selected[2] = {-1, -1};
    PairOutcome out = {};
    mem_pair_resolve_batch_post(opt, &fb.bns, nullptr, pes, 7, reads, a, nullptr, cache,
                               gcnt, 0, n_pri, selected, out.q, &out.flag, &out.paired);
    for (int i = 0; i < 2; ++i) {
        out.rb[i] = out.paired ? a[i].a[selected[i]].rb : -1;
        out.score[i] = out.paired ? a[i].a[selected[i]].score : -1;
    }
    free(cache);
    return out;
}
} // namespace

TEST_CASE("alt_lift_load keeps usable lines and skips every other kind" * doctest::test_suite("unit/alt_lift")) {
    FakeBns fb({{"chr1", 10000, 0}, {"a_ok", 100, 1}, {"a_dup", 100, 1}, {"a_unmapped", 100, 1},
                {"a_pos0", 100, 1}, {"a_star", 100, 1}, {"a_supp", 100, 1}, {"a_sec", 100, 1},
                {"a_alttarget", 100, 1}, {"a_badop", 100, 1}, {"a_trailing", 100, 1},
                {"a_len", 100, 1}, {"a_pastend", 100, 1}, {"a_eqx", 100, 1}, {"a_notarget", 100, 1},
                {"a_short", 100, 1}, {"a_end", 100, 1}, {"a_posmax", 100, 1}, {"a_plus", 100, 1},
                {"a_space", 100, 1}, {"a_posend", 100, 1}, {"a_posalt", 100, 1}});
    TmpFile f(
        "a_ok\t0\tchr1\t101\t60\t100M\t*\t0\t0\t*\t*\n"
        "a_dup\t0\tchr1\t201\t60\t100M\t*\t0\t0\t*\t*\n"
        "a_dup\t0\tchr1\t901\t60\t100M\t*\t0\t0\t*\t*\n"        // second line: ignored
        "a_unmapped\t4\tchr1\t1\t0\t100M\t*\t0\t0\t*\t*\n"
        "a_pos0\t0\tchr1\t0\t60\t100M\t*\t0\t0\t*\t*\n"
        "a_star\t0\tchr1\t1\t60\t*\t*\t0\t0\t*\t*\n"
        "a_supp\t2048\tchr1\t1\t60\t100M\t*\t0\t0\t*\t*\n"
        "a_sec\t256\tchr1\t1\t60\t100M\t*\t0\t0\t*\t*\n"
        "a_alttarget\t0\ta_ok\t1\t60\t100M\t*\t0\t0\t*\t*\n"
        "a_badop\t0\tchr1\t1\t60\t100Q\t*\t0\t0\t*\t*\n"
        "a_trailing\t0\tchr1\t1\t60\t100M5\t*\t0\t0\t*\t*\n"
        "a_len\t0\tchr1\t1\t60\t50M\t*\t0\t0\t*\t*\n"
        "a_pastend\t0\tchr1\t9951\t60\t100M\t*\t0\t0\t*\t*\n"
        "a_eqx\t0\tchr1\t501\t60\t40=1P10X50=\t*\t0\t0\t*\t*\n"
        "a_notarget\t0\tchrZ\t1\t60\t100M\t*\t0\t0\t*\t*\n"
        "a_short\t0\tchr1\t1\t60\n"
        "a_end\t0\tchr1\t9901\t60\t100M\t*\t0\t0\t*\t*\n"     // ends exactly at the target end
        "not_in_index\t0\tchr1\t1\t60\t10M\t*\t0\t0\t*\t*\n"
        // POS = LONG_MAX: POS - 1 + the CIGAR's 100 reference bases would overflow
        "a_posmax\t0\tchr1\t9223372036854775807\t60\t100M\t*\t0\t0\t*\t*\n"
        "a_plus\t0\tchr1\t1\t60\t10M+90M\t*\t0\t0\t*\t*\n"   // strtol alone takes the sign
        "a_space\t0\tchr1\t1\t60\t10M 90M\t*\t0\t0\t*\t*\n"   // and the space
        // POS one past the end with a CIGAR that consumes no reference: only the
        // POS check rejects it (the end check sees 0 reference bases fit in 0)
        "a_posend\t0\tchr1\t10001\t60\t100S\t*\t0\t0\t*\t*\n"
        // POS past the end is counted after the other skip reasons: an ALT
        // target, and a second line for a contig that already has one
        "a_posalt\t0\ta_ok\t500\t60\t100M\t*\t0\t0\t*\t*\n"
        "a_ok\t0\tchr1\t99999\t60\t100M\t*\t0\t0\t*\t*\n");
    LiftGuard g(alt_lift_load(f.path.c_str(), &fb.bns));
    REQUIRE(g.l != nullptr);
    CHECK(g.l->n_lift == 4);                     // a_ok, a_dup (first line), a_eqx, a_end
    CHECK(g.l->a[16].trid == 0);                 // a_end: one base short of a_pastend's overrun
    CHECK(g.l->a[16].tpos == 9900);
    CHECK(g.l->a[1].trid == 0);
    CHECK(g.l->a[1].tpos == 100);
    CHECK(g.l->a[1].tlen == 10000);
    CHECK(g.l->a[2].tpos == 200);                // first a_dup line wins
    CHECK(g.l->a[13].trid == 0);                 // =, X and P parse
    CHECK(g.l->a[13].tpos == 500);
    CHECK(g.l->a[13].n_cigar == 3);              // P dropped
    for (int rid : {0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 17, 18, 19, 20, 21}) {
        CAPTURE(rid);
        CHECK(g.l->a[rid].trid == -1);
    }
    // every skipped line, by the reason the startup warning gives
    CHECK(g.l->skip.no_ctg == 1);                // not_in_index
    CHECK(g.l->skip.no_target == 1);             // a_notarget
    CHECK(g.l->skip.alt_target == 2);            // a_alttarget, a_posalt
    CHECK(g.l->skip.unmapped == 3);              // a_unmapped, a_pos0, a_star
    CHECK(g.l->skip.not_primary == 2);           // a_supp, a_sec
    CHECK(g.l->skip.dup == 2);                   // a_dup's and a_ok's second lines
    CHECK(g.l->skip.bad == 9);                   // a_badop, a_trailing, a_len, a_pastend,
                                                 // a_short, a_posmax, a_plus, a_space,
                                                 // a_posend
    CHECK(alt_lift_load("/nonexistent/x.alt", &fb.bns) == nullptr);
}

TEST_CASE("alt_lift_load skips a line whose FLAG or POS is not a whole integer" * doctest::test_suite("unit/alt_lift")) {
    // Read as numbers with no end check, "abc" would be FLAG 0 (a forward-strand
    // primary alignment) and "16x" FLAG 16; both must be skipped as malformed.
    FakeBns fb({{"chr1", 10000, 0}, {"a_flagtext", 100, 1}, {"a_flagtrail", 100, 1},
                {"a_flagsign", 100, 1}, {"a_postext", 100, 1}, {"a_postrail", 100, 1},
                {"a_flagempty", 100, 1}, {"a_rev", 100, 1}});
    TmpFile f(
        "a_flagtext\tabc\tchr1\t101\t60\t100M\t*\t0\t0\t*\t*\n"
        "a_flagtrail\t16x\tchr1\t101\t60\t100M\t*\t0\t0\t*\t*\n"
        "a_flagsign\t-16\tchr1\t101\t60\t100M\t*\t0\t0\t*\t*\n"
        "a_postext\t0\tchr1\tx\t60\t100M\t*\t0\t0\t*\t*\n"
        "a_postrail\t0\tchr1\t101x\t60\t100M\t*\t0\t0\t*\t*\n"
        "a_flagempty\t\tchr1\t101\t60\t100M\t*\t0\t0\t*\t*\n"
        "a_rev\t16\tchr1\t301\t60\t100M\t*\t0\t0\t*\t*\n");
    LiftGuard g(alt_lift_load(f.path.c_str(), &fb.bns));
    REQUIRE(g.l != nullptr);
    CHECK(g.l->n_lift == 1);                     // only a_rev
    for (int rid : {1, 2, 3, 4, 5, 6}) {
        CAPTURE(rid);
        CHECK(g.l->a[rid].trid == -1);
    }
    CHECK(g.l->a[7].trid == 0);
    CHECK(g.l->a[7].rev == 1);
    CHECK(g.l->a[7].tpos == 300);
    CHECK(g.l->skip.bad == 6);
    CHECK(g.l->skip.unmapped + g.l->skip.no_ctg + g.l->skip.dup == 0);
}

TEST_CASE("alt_lift_interval: forward ALT, before and across a deletion" * doctest::test_suite("unit/alt_lift")) {
    FakeBns fb(kCtgs);
    TmpFile f(kAlt);
    LiftGuard g(alt_lift_load(f.path.c_str(), &fb.bns));
    REQUIRE(g.l != nullptr);
    alt_lift_iv_t iv;
    REQUIRE(alt_lift_interval(g.l, ALT_FWD, 100, 250, &iv) == 0);
    CHECK(iv.trid == CHR1);
    CHECK(iv.rev == 0);
    CHECK(iv.ps == 3100);
    CHECK(iv.pe == 3250);
    CHECK(iv.n_aln == 150);
    CHECK(iv.n_clip == 0);
    // [1450, 1600) spans the deletion: 50 bases before, 3 deleted ref bases, 100 after.
    REQUIRE(alt_lift_interval(g.l, ALT_FWD, 1450, 1600, &iv) == 0);
    CHECK(iv.ps == 4450);
    CHECK(iv.pe == 4603);
    CHECK(iv.n_aln == 150);
    // ending exactly at the deletion: the end stays before it
    REQUIRE(alt_lift_interval(g.l, ALT_FWD, 1400, 1500, &iv) == 0);
    CHECK(iv.ps == 4400);
    CHECK(iv.pe == 4500);
    // starting exactly after it: shifted by 3
    REQUIRE(alt_lift_interval(g.l, ALT_FWD, 1500, 1650, &iv) == 0);
    CHECK(iv.ps == 4503);
    CHECK(iv.pe == 4653);
    // empty, inverted and out-of-contig intervals, and a primary contig, do not lift
    CHECK(alt_lift_interval(g.l, ALT_FWD, 50, 50, &iv) == -1);
    CHECK(alt_lift_interval(g.l, ALT_FWD, 60, 50, &iv) == -1);
    CHECK(alt_lift_interval(g.l, ALT_FWD, 1900, 2000, &iv) == -1);
    CHECK(alt_lift_interval(g.l, CHR1, 0, 10, &iv) == -1);
}

TEST_CASE("alt_lift_interval: reverse-complemented ALT mirrors the interval" * doctest::test_suite("unit/alt_lift")) {
    FakeBns fb(kCtgs);
    TmpFile f(kAlt);
    LiftGuard g(alt_lift_load(f.path.c_str(), &fb.bns));
    REQUIRE(g.l != nullptr);
    alt_lift_iv_t iv;
    // contig [1847, 1997) is CIGAR-query [0, 150): chr1 [3000, 3150)
    REQUIRE(alt_lift_interval(g.l, ALT_REV, 1847, 1997, &iv) == 0);
    CHECK(iv.rev == 1);
    CHECK(iv.ps == 3000);
    CHECK(iv.pe == 3150);
    // contig [0, 150) is CIGAR-query [1847, 1997): chr1 [4850, 5000)
    REQUIRE(alt_lift_interval(g.l, ALT_REV, 0, 150, &iv) == 0);
    CHECK(iv.ps == 4850);
    CHECK(iv.pe == 5000);
}

TEST_CASE("alt_lift_interval: clips and insertions are reported, not hidden" * doctest::test_suite("unit/alt_lift")) {
    FakeBns fb({{"chr1", 10000, 0}, {"a_ins", 230, 1}});
    // 10S50M20I150M (230 query bases, 200 reference bases) at chr1:5001
    TmpFile f("a_ins\t0\tchr1\t5001\t60\t10S50M20I150M\t*\t0\t0\t*\t*\n");
    LiftGuard g(alt_lift_load(f.path.c_str(), &fb.bns));
    REQUIRE(g.l != nullptr);
    REQUIRE(g.l->n_lift == 1);
    alt_lift_iv_t iv;
    // wholly inside the clip: maps to the alignment start, all of it clipped
    REQUIRE(alt_lift_interval(g.l, 1, 0, 5, &iv) == 0);
    CHECK(iv.ps == 5000);
    CHECK(iv.pe == 5001);
    CHECK(iv.n_aln == 0);
    CHECK(iv.n_clip == 5);
    // across the clip edge: 5 clipped, 45 aligned
    REQUIRE(alt_lift_interval(g.l, 1, 5, 55, &iv) == 0);
    CHECK(iv.ps == 5000);
    CHECK(iv.pe == 5045);
    CHECK(iv.n_aln == 45);
    CHECK(iv.n_clip == 5);
    // across the insertion [60, 80): 20 query bases have no reference
    REQUIRE(alt_lift_interval(g.l, 1, 40, 100, &iv) == 0);
    CHECK(iv.ps == 5030);
    CHECK(iv.pe == 5070);
    CHECK(iv.n_aln == 40);
    CHECK(iv.n_clip == 0);
    CHECK(iv.max_diag == 20);
    CHECK(iv.edge_ins == 0);
    // wholly inside the insertion
    REQUIRE(alt_lift_interval(g.l, 1, 62, 78, &iv) == 0);
    CHECK(iv.ps == 5050);
    CHECK(iv.pe == 5051);
    CHECK(iv.n_aln == 0);
    CHECK(iv.edge_ins == 1);
    // starting in the insertion's first base, or 3 bases into it: no left anchor
    REQUIRE(alt_lift_interval(g.l, 1, 60, 120, &iv) == 0);
    CHECK(iv.edge_ins == 1);
    CHECK(iv.max_diag == 20);
    REQUIRE(alt_lift_interval(g.l, 1, 63, 120, &iv) == 0);
    CHECK(iv.edge_ins == 1);
    CHECK(iv.max_diag == 17);
    // ending in its last base, or 5 bases into it
    REQUIRE(alt_lift_interval(g.l, 1, 20, 80, &iv) == 0);
    CHECK(iv.edge_ins == 1);
    REQUIRE(alt_lift_interval(g.l, 1, 20, 65, &iv) == 0);
    CHECK(iv.edge_ins == 1);
    CHECK(iv.max_diag == 5);
    // ending just before it and starting just after it: anchored on both ends
    REQUIRE(alt_lift_interval(g.l, 1, 20, 60, &iv) == 0);
    CHECK(iv.edge_ins == 0);
    CHECK(iv.max_diag == 0);
    REQUIRE(alt_lift_interval(g.l, 1, 80, 120, &iv) == 0);
    CHECK(iv.edge_ins == 0);
    CHECK(iv.max_diag == 0);
}

TEST_CASE("alt_lift_interval: the band an interval needs is its widest diagonal excursion" * doctest::test_suite("unit/alt_lift")) {
    FakeBns fb(kCtgs);
    TmpFile f(kAlt);
    LiftGuard g(alt_lift_load(f.path.c_str(), &fb.bns));
    REQUIRE(g.l != nullptr);
    alt_lift_iv_t iv;
    // alt_bal (100M5I30M5D100M): across both indels the net span change is 0, but
    // the alignment runs 5 bases off the diagonal between them.
    REQUIRE(alt_lift_interval(g.l, ALT_BAL, 50, 200, &iv) == 0);
    CHECK(iv.ps == 70050);
    CHECK(iv.pe == 70200);
    CHECK(iv.max_diag == 5);
    // part of the insertion only (3 of its 5 bases in the interval)
    REQUIRE(alt_lift_interval(g.l, ALT_BAL, 102, 200, &iv) == 0);
    CHECK(iv.edge_ins == 1);
    CHECK(iv.max_diag == 3);
    // the deletion only: the excursion is the deletion
    REQUIRE(alt_lift_interval(g.l, ALT_BAL, 120, 200, &iv) == 0);
    CHECK(iv.max_diag == 5);
    CHECK(iv.pe - iv.ps == 85);
    // a reverse-complemented ALT walks the CIGAR in the primary's orientation:
    // alt_rev [0, 150) is CIGAR-query [1847, 1997), after the 3D, so no excursion;
    // its contig [400, 550) is CIGAR-query [1447, 1597), across it
    REQUIRE(alt_lift_interval(g.l, ALT_REV, 0, 150, &iv) == 0);
    CHECK(iv.max_diag == 0);
    REQUIRE(alt_lift_interval(g.l, ALT_REV, 400, 550, &iv) == 0);
    CHECK(iv.max_diag == 3);
    CHECK(iv.ps == 4447);
    CHECK(iv.pe == 4600);
}

TEST_CASE("mem_alt_project: a lifted hit and a primary hit at one locus reduce to one" * doctest::test_suite("unit/alt_lift")) {
    FakeBns fb(kCtgs);
    TmpFile f(kAlt);
    LiftGuard g(alt_lift_load(f.path.c_str(), &fb.bns));
    REQUIRE(g.l != nullptr);
    OptGuard o;
    o.o->alt_lift = g.l;

    SUBCASE("a tie keeps the primary-assembly hit, ALT region first") {
        mem_alnreg_t a[2] = {reg(fb, ALT_FWD, 100, 250, 0, 150, 2), reg(fb, CHR1, 3100, 3250, 0, 150, 1)};
        REQUIRE(mem_alt_project(o.o, &fb.bns, 2, a) == 1);
        CHECK(a[0].seedcov == 1);
        CHECK(a[0].rid == CHR1);
        CHECK(a[0].rb == 3100);
        CHECK(a[0].re == 3250);
        CHECK(int(a[0].is_alt) == 0);
        CHECK(int(a[0].alt_bw) == 0);
    }
    SUBCASE("a tie keeps the primary-assembly hit, primary region first") {
        mem_alnreg_t a[2] = {reg(fb, CHR1, 3100, 3250, 0, 150, 1), reg(fb, ALT_FWD, 100, 250, 0, 150, 2)};
        REQUIRE(mem_alt_project(o.o, &fb.bns, 2, a) == 1);
        CHECK(a[0].seedcov == 1);
        CHECK(int(a[0].alt_bw) == 0);
    }
    SUBCASE("a higher-scoring ALT hit wins and is rewritten onto the primary") {
        mem_alnreg_t a[2] = {reg(fb, ALT_FWD, 100, 250, 0, 155, 2), reg(fb, CHR1, 3100, 3250, 0, 150, 1)};
        REQUIRE(mem_alt_project(o.o, &fb.bns, 2, a) == 1);
        CHECK(a[0].seedcov == 2);
        CHECK(a[0].rid == CHR1);
        CHECK(a[0].rb == 3100);
        CHECK(a[0].re == 3250);
        CHECK(a[0].qb == 0);
        CHECK(a[0].qe == 150);
        CHECK(a[0].score == 155);
        CHECK(int(a[0].is_alt) == 0);
        CHECK(int(a[0].alt_bw) == 1);            // lifted, no ALT indel in [100, 250): floor 0
    }
    SUBCASE("opposite strands are different loci: both kept, by score") {
        mem_alnreg_t a[2] = {reg(fb, ALT_FWD, 100, 250, 1, 160, 2), reg(fb, CHR1, 3100, 3250, 0, 150, 1)};
        REQUIRE(mem_alt_project(o.o, &fb.bns, 2, a) == 2);
        const int64_t l2 = fb.bns.l_pac << 1;
        CHECK(a[0].seedcov == 2);
        CHECK(a[0].rb == l2 - 3250);
        CHECK(a[0].re == l2 - 3100);
        CHECK(a[1].seedcov == 1);
    }
    SUBCASE("a forward hit on a FLAG-16 ALT lifts to the reverse strand") {
        mem_alnreg_t a[1] = {reg(fb, ALT_REV, 1847, 1997, 0, 150, 2)};
        REQUIRE(mem_alt_project(o.o, &fb.bns, 1, a) == 1);
        const int64_t l2 = fb.bns.l_pac << 1;
        CHECK(a[0].rid == CHR1);
        CHECK(a[0].rb == l2 - 3150);
        CHECK(a[0].re == l2 - 3000);
        CHECK(int(a[0].is_alt) == 0);
    }
}

TEST_CASE("mem_alt_project: hits the lift cannot describe stay on the ALT" * doctest::test_suite("unit/alt_lift")) {
    FakeBns fb(kCtgs);
    TmpFile f(kAlt);
    LiftGuard g(alt_lift_load(f.path.c_str(), &fb.bns));
    REQUIRE(g.l != nullptr);
    OptGuard o;
    o.o->alt_lift = g.l;
    const int64_t clip_off = fb.bns.anns[ALT_CLIP].offset, del_off = fb.bns.anns[ALT_BIGDEL].offset;

    SUBCASE("inside the ALT-only clip") {
        mem_alnreg_t a[1] = {reg(fb, ALT_CLIP, 0, 100, 0, 100, 2, 0, 100)};
        REQUIRE(mem_alt_project(o.o, &fb.bns, 1, a) == 1);
        CHECK(a[0].rid == ALT_CLIP);
        CHECK(a[0].rb == clip_off);
        CHECK(int(a[0].is_alt) == 1);
    }
    SUBCASE("partly in the clip") {
        mem_alnreg_t a[1] = {reg(fb, ALT_CLIP, 50, 200, 0, 150, 2)};
        REQUIRE(mem_alt_project(o.o, &fb.bns, 1, a) == 1);
        CHECK(a[0].rid == ALT_CLIP);
        CHECK(int(a[0].is_alt) == 1);
    }
    SUBCASE("clear of the clip: lifted") {
        mem_alnreg_t a[1] = {reg(fb, ALT_CLIP, 120, 270, 0, 150, 2)};
        REQUIRE(mem_alt_project(o.o, &fb.bns, 1, a) == 1);
        CHECK(a[0].rid == CHR1);
        CHECK(a[0].rb == 20020);
        CHECK(a[0].re == 20170);
    }
    SUBCASE("across an ALT deletion wider than the band") {
        mem_alnreg_t a[1] = {reg(fb, ALT_BIGDEL, 25, 175, 0, 150, 2)};
        REQUIRE(mem_alt_project(o.o, &fb.bns, 1, a) == 1);
        CHECK(a[0].rid == ALT_BIGDEL);
        CHECK(a[0].rb == del_off + 25);
        CHECK(int(a[0].is_alt) == 1);
    }
    SUBCASE("mostly in an ALT insertion") {
        // alt_ins [60, 210): 40 + 30 aligned bases, 80 inserted; the excursion (80)
        // is within the band, so only the at-least-half-aligned rule keeps it
        mem_alnreg_t a[1] = {reg(fb, ALT_INS, 60, 210, 0, 150, 2)};
        REQUIRE(mem_alt_project(o.o, &fb.bns, 1, a) == 1);
        CHECK(a[0].rid == ALT_INS);
        CHECK(int(a[0].is_alt) == 1);
        CHECK(int(a[0].alt_bw) == 0);
    }
    SUBCASE("starting or ending inside an ALT insertion") {
        // alt_bal's 5I is ALT [100, 105): [102, 252) starts in it, [0, 103) ends in it
        mem_alnreg_t a[2] = {reg(fb, ALT_BAL, 102, 232, 0, 130, 2, 0, 130),
                             reg(fb, ALT_BAL, 0, 103, 1, 103, 3, 0, 103)};
        REQUIRE(mem_alt_project(o.o, &fb.bns, 2, a) == 2);
        CHECK(a[0].rid == ALT_BAL);
        CHECK(a[1].rid == ALT_BAL);
        CHECK(int(a[0].is_alt) == 1);
        CHECK(int(a[1].is_alt) == 1);
    }
}

TEST_CASE("mem_alt_project: band floor and the exact band-width edge" * doctest::test_suite("unit/alt_lift")) {
    FakeBns fb(kCtgs);
    TmpFile f(kAlt);
    LiftGuard g(alt_lift_load(f.path.c_str(), &fb.bns));
    REQUIRE(g.l != nullptr);
    OptGuard o;
    o.o->alt_lift = g.l;
    REQUIRE(o.o->w == 100);

    SUBCASE("indels that cancel inside the region: lifted with their excursion as the floor") {
        mem_alnreg_t a[1] = {reg(fb, ALT_BAL, 50, 200, 0, 150, 2)};
        a[0].dedup_re_rank = 5;                  // its `re` moves, so the rank must reset
        REQUIRE(mem_alt_project(o.o, &fb.bns, 1, a) == 1);
        CHECK(a[0].rid == CHR1);
        CHECK(a[0].rb == 70050);
        CHECK(a[0].re == 70200);
        CHECK(int(a[0].alt_bw) == 6);            // 1 + the 5-base excursion
        CHECK(a[0].dedup_re_rank == 0);
    }
    SUBCASE("a deletion of exactly -w lifts") {
        mem_alnreg_t a[1] = {reg(fb, ALT_WDEL, 25, 175, 0, 150, 2)};
        REQUIRE(mem_alt_project(o.o, &fb.bns, 1, a) == 1);
        CHECK(a[0].rid == CHR1);
        CHECK(a[0].rb == 60025);
        CHECK(a[0].re == 60275);
        CHECK(int(a[0].alt_bw) == 101);
    }
    SUBCASE("one base past -w stays") {
        o.o->w = 99;
        mem_alnreg_t a[1] = {reg(fb, ALT_WDEL, 25, 175, 0, 150, 2)};
        REQUIRE(mem_alt_project(o.o, &fb.bns, 1, a) == 1);
        CHECK(a[0].rid == ALT_WDEL);
        CHECK(int(a[0].is_alt) == 1);
    }
    // bwa_gen_cigar4 caps the band of a 150-base region at 35 when the lengths
    // match (see "bwa_gen_cigar_band" in test_gen_cigar_no_md.cpp), whatever -w
    // is: a 35-base
    // excursion fits it, a 40-base one does not.
    SUBCASE("balanced indels 35 bases off the diagonal: lifted at the read-length band") {
        // alt_b35 [50, 200): 50M, 35I, 30M, 35D, 35M; chr1 [80050, 80200)
        mem_alnreg_t a[1] = {reg(fb, ALT_B35, 50, 200, 0, 150, 2)};
        REQUIRE(mem_alt_project(o.o, &fb.bns, 1, a) == 1);
        CHECK(a[0].rid == CHR1);
        CHECK(a[0].rb == 80050);
        CHECK(a[0].re == 80200);
        CHECK(int(a[0].alt_bw) == 36);
    }
    SUBCASE("balanced indels 40 bases off the diagonal: past the read-length band, stays") {
        mem_alnreg_t a[1] = {reg(fb, ALT_B40, 50, 200, 0, 150, 2)};
        REQUIRE(mem_alt_project(o.o, &fb.bns, 1, a) == 1);
        CHECK(a[0].rid == ALT_B40);
        CHECK(int(a[0].is_alt) == 1);
        CHECK(int(a[0].alt_bw) == 0);
    }
    SUBCASE("a deletion of exactly MEM_ALT_BW_MAX lifts with -w above it") {
        o.o->w = 2000;
        mem_alnreg_t a[1] = {reg(fb, ALT_MAXDEL, 25, 175, 0, 150, 2)};
        REQUIRE(mem_alt_project(o.o, &fb.bns, 1, a) == 1);
        CHECK(a[0].rid == CHR1);
        CHECK(a[0].rb == 90025);
        CHECK(a[0].re == 90025 + 150 + 1022);
        CHECK(int(a[0].alt_bw) == MEM_ALT_BW_MAX + 1);
    }
    SUBCASE("a deletion wider than MEM_ALT_BW_MAX stays, even with -w above it") {
        // the read-length band (1103: the length difference plus 3) and -w both
        // allow it; only the width of alt_bw does not
        o.o->w = 2000;
        mem_alnreg_t a[1] = {reg(fb, ALT_HUGEDEL, 25, 175, 0, 150, 2)};
        REQUIRE(mem_alt_project(o.o, &fb.bns, 1, a) == 1);
        CHECK(a[0].rid == ALT_HUGEDEL);
        CHECK(int(a[0].is_alt) == 1);
        CHECK(int(a[0].alt_bw) == 0);
    }
}

TEST_CASE("mem_alt_project: what counts as one locus" * doctest::test_suite("unit/alt_lift")) {
    FakeBns fb(kCtgs);
    TmpFile f(kAlt);
    LiftGuard g(alt_lift_load(f.path.c_str(), &fb.bns));
    REQUIRE(g.l != nullptr);
    OptGuard o;
    o.o->alt_lift = g.l;

    SUBCASE("one repeat unit away in a tandem repeat: two loci, both kept") {
        // lifted [3100, 3250) and a primary hit 50 bases on (100 of 150 overlap,
        // below dedup's 0.95): ordinary dedup keeps both, and so does the lift
        mem_alnreg_t a[2] = {reg(fb, ALT_FWD, 100, 250, 0, 150, 2), reg(fb, CHR1, 3150, 3300, 0, 145, 1)};
        REQUIRE(mem_alt_project(o.o, &fb.bns, 2, a) == 2);
        CHECK(a[0].seedcov == 2);
        CHECK(a[1].seedcov == 1);
        CHECK(a[1].rb == 3150);
    }
    SUBCASE("a few bases off the lifted diagonal: one locus") {
        mem_alnreg_t a[2] = {reg(fb, ALT_FWD, 100, 250, 0, 150, 2), reg(fb, CHR1, 3104, 3254, 0, 145, 1)};
        REQUIRE(mem_alt_project(o.o, &fb.bns, 2, a) == 1);
        CHECK(a[0].seedcov == 2);
    }
    SUBCASE("a primary piece inside the lifted hit: one locus") {
        mem_alnreg_t a[2] = {reg(fb, ALT_FWD, 100, 250, 0, 150, 2),
                             reg(fb, CHR1, 3100, 3171, 0, 71, 1, 0, 71)};
        REQUIRE(mem_alt_project(o.o, &fb.bns, 2, a) == 1);
        CHECK(a[0].seedcov == 2);
    }
    SUBCASE("only a kept region drops another, whatever the input order") {
        // X lifted [3100, 3250) 150; Y primary [3104, 3254) 140 and Z primary
        // [3096, 3246) 150 each share X's locus. Z wins its tie with X, so X is
        // dropped; Y lost only to X, so Y is kept. Visited in input order, X would
        // first drop Y and then itself lose to Z, leaving only Z.
        mem_alnreg_t a[3] = {reg(fb, ALT_FWD, 100, 250, 0, 150, 2), reg(fb, CHR1, 3104, 3254, 0, 140, 3),
                             reg(fb, CHR1, 3096, 3246, 0, 150, 1)};
        REQUIRE(mem_alt_project(o.o, &fb.bns, 3, a) == 2);
        CHECK(a[0].seedcov == 1);
        CHECK(a[1].seedcov == 3);
    }
    SUBCASE("read intervals that share no base: two loci, both kept") {
        // the same reference window, but read bases [0, 150) and [150, 300)
        mem_alnreg_t a[2] = {reg(fb, ALT_FWD, 100, 250, 0, 150, 2),
                             reg(fb, CHR1, 3100, 3250, 0, 140, 1, 150, 300)};
        REQUIRE(mem_alt_project(o.o, &fb.bns, 2, a) == 2);
        CHECK(a[0].seedcov == 2);
        CHECK(a[0].rid == CHR1);
        CHECK(a[1].seedcov == 1);
    }
    SUBCASE("read intervals sharing exactly half: two loci, both kept") {
        // [0, 150) and [75, 225) share 75 bases, half of 150, not more than half
        mem_alnreg_t a[2] = {reg(fb, ALT_FWD, 100, 250, 0, 150, 2),
                             reg(fb, CHR1, 3100, 3250, 0, 140, 1, 75, 225)};
        REQUIRE(mem_alt_project(o.o, &fb.bns, 2, a) == 2);
        CHECK(a[0].seedcov == 2);
        CHECK(a[1].seedcov == 1);
    }
    SUBCASE("read intervals sharing one base more than half: one locus") {
        // [0, 150) and [74, 224) share 76 bases: one locus. The lifted hit leaves
        // out 74 of the primary hit's bases, and 150 is not more than 140 + 74, so
        // the primary hit is kept.
        mem_alnreg_t a[2] = {reg(fb, ALT_FWD, 100, 250, 0, 150, 2),
                             reg(fb, CHR1, 3100, 3250, 0, 140, 1, 74, 224)};
        REQUIRE(mem_alt_project(o.o, &fb.bns, 2, a) == 1);
        CHECK(a[0].seedcov == 1);
        CHECK(a[0].qb == 74);
    }
    SUBCASE("two lifted hits at one primary locus: the better one is kept") {
        // alt_fwd [100, 250) and a reverse-strand hit on alt_rev (FLAG 16) at its
        // contig [1747, 1897), CIGAR-query [100, 250): both lift to forward chr1
        // [3100, 3250)
        for (int order = 0; order < 2; ++order) {
            CAPTURE(order);
            mem_alnreg_t x = reg(fb, ALT_FWD, 100, 250, 0, 148, 2), y = reg(fb, ALT_REV, 1747, 1897, 1, 150, 3);
            mem_alnreg_t a[2] = {order ? y : x, order ? x : y};
            REQUIRE(mem_alt_project(o.o, &fb.bns, 2, a) == 1);
            CHECK(a[0].seedcov == 3);
            CHECK(a[0].rid == CHR1);
            CHECK(a[0].rb == 3100);
            CHECK(a[0].re == 3250);
            CHECK(int(a[0].alt_bw) == 1);
        }
    }
    SUBCASE("two lifted hits tied on score, rb and qb: the lower ALT contig index wins, in either order") {
        // both lift to chr1 [3100, 3250) with read [0, 150) and score 150; alt_fwd
        // (rid 1) comes before alt_rev (rid 2) in the index. The survivor's other
        // fields (here sub, which MAPQ reads) come with it.
        for (int order = 0; order < 2; ++order) {
            CAPTURE(order);
            mem_alnreg_t x = reg(fb, ALT_FWD, 100, 250, 0, 150, 2), y = reg(fb, ALT_REV, 1747, 1897, 1, 150, 3);
            x.sub = 10;
            y.sub = 140;
            mem_alnreg_t a[2] = {order ? y : x, order ? x : y};
            REQUIRE(mem_alt_project(o.o, &fb.bns, 2, a) == 1);
            CHECK(a[0].seedcov == 2);
            CHECK(a[0].sub == 10);
            CHECK(a[0].rid == CHR1);
            CHECK(a[0].rb == 3100);
            CHECK(int(a[0].alt_bw) == 1);
        }
    }
    SUBCASE("the tie goes by ALT contig index even where the pre-lift rb disagrees") {
        // a reverse-strand hit on alt_rev (rid 2) and a forward hit on alt_fwd2
        // (rid 12), both lifting to chr1 [3100, 3250) at 150: alt_rev comes first
        // in the index, though its reverse-strand rb before the lift is the larger
        for (int order = 0; order < 2; ++order) {
            CAPTURE(order);
            mem_alnreg_t x = reg(fb, ALT_REV, 1747, 1897, 1, 150, 3), y = reg(fb, ALT_FWD2, 100, 250, 0, 150, 4);
            REQUIRE(x.rb > y.rb);
            mem_alnreg_t a[2] = {order ? y : x, order ? x : y};
            REQUIRE(mem_alt_project(o.o, &fb.bns, 2, a) == 1);
            CHECK(a[0].seedcov == 3);
            CHECK(a[0].rb == 3100);
        }
    }
}

TEST_CASE("mem_alt_project: a lifted hit replaces a primary hit only if it outscores it by more than the read bases it leaves out" * doctest::test_suite("unit/alt_lift")) {
    FakeBns fb(kCtgs);
    TmpFile f(kAlt);
    LiftGuard g(alt_lift_load(f.path.c_str(), &fb.bns));
    REQUIRE(g.l != nullptr);
    OptGuard o;
    o.o->alt_lift = g.l;

    SUBCASE("a shorter, higher-scoring lifted hit leaves the full primary hit as it was") {
        // A read whose first 100 bases match the ALT's last 100 (score 100) and
        // whose 150 bases align to chr1 with mismatches (score 90). Replacing the
        // primary hit would soft-clip the 50 bases only it aligned, worth 50 at
        // -A 1, and 100 is not more than 90 + 50.
        for (int order = 0; order < 2; ++order) {
            CAPTURE(order);
            mem_alnreg_t l = reg(fb, ALT_FWD, 100, 200, 0, 100, 2, 0, 100), p = reg(fb, CHR1, 3100, 3250, 0, 90, 1);
            mem_alnreg_t a[2] = {order ? p : l, order ? l : p};
            REQUIRE(mem_alt_project(o.o, &fb.bns, 2, a) == 1);
            CHECK(a[0].seedcov == 1);
            CHECK(a[0].rid == CHR1);
            CHECK(a[0].qb == 0);
            CHECK(a[0].qe == 150);
            CHECK(a[0].rb == 3100);
            CHECK(a[0].re == 3250);
            CHECK(a[0].score == 90);
            CHECK(int(a[0].alt_bw) == 0);
        }
    }
    SUBCASE("the shorter lifted hit at the threshold: 140 keeps the primary hit, 141 replaces it") {
        // the same two hits: the 50 bases are worth 50, so 90 + 50 is a tie
        for (int sc = 140; sc <= 141; ++sc) {
            for (int order = 0; order < 2; ++order) {
                CAPTURE(sc);
                CAPTURE(order);
                mem_alnreg_t l = reg(fb, ALT_FWD, 100, 200, 0, sc, 2, 0, 100), p = reg(fb, CHR1, 3100, 3250, 0, 90, 1);
                mem_alnreg_t a[2] = {order ? p : l, order ? l : p};
                REQUIRE(mem_alt_project(o.o, &fb.bns, 2, a) == 1);
                CHECK(a[0].seedcov == (sc == 141 ? 2 : 1));
                CHECK(a[0].qe == (sc == 141 ? 100 : 150));
            }
        }
    }
    SUBCASE("at -A 2 the left-out bases are worth twice as much: 190 keeps the primary hit, 191 replaces it") {
        // the same two hits: the 50 bases are worth 100 at -A 2, so 90 + 100 is a
        // tie
        o.o->a = 2;
        bwa_fill_scmat(o.o->a, o.o->b, o.o->mat);
        for (int sc = 190; sc <= 191; ++sc) {
            for (int order = 0; order < 2; ++order) {
                CAPTURE(sc);
                CAPTURE(order);
                mem_alnreg_t l = reg(fb, ALT_FWD, 100, 200, 0, sc, 2, 0, 100), p = reg(fb, CHR1, 3100, 3250, 0, 90, 1);
                mem_alnreg_t a[2] = {order ? p : l, order ? l : p};
                REQUIRE(mem_alt_project(o.o, &fb.bns, 2, a) == 1);
                CHECK(a[0].seedcov == (sc == 191 ? 2 : 1));
                CHECK(a[0].qe == (sc == 191 ? 100 : 150));
            }
        }
    }
    SUBCASE("a lifted hit inside the primary hit on both ends is charged for both: 140 keeps it, 141 replaces it") {
        // lifted [5, 145), primary [0, 150): the lifted hit leaves out the primary
        // hit's bases [0, 5) and [145, 150), worth 10, so 130 + 10 is a tie
        for (int sc = 140; sc <= 141; ++sc) {
            for (int order = 0; order < 2; ++order) {
                CAPTURE(sc);
                CAPTURE(order);
                mem_alnreg_t l = reg(fb, ALT_FWD, 105, 245, 0, sc, 2, 5, 145), p = reg(fb, CHR1, 3100, 3250, 0, 130, 1);
                mem_alnreg_t a[2] = {order ? p : l, order ? l : p};
                REQUIRE(mem_alt_project(o.o, &fb.bns, 2, a) == 1);
                CHECK(a[0].seedcov == (sc == 141 ? 2 : 1));
                CHECK(a[0].qb == (sc == 141 ? 5 : 0));
                CHECK(a[0].qe == (sc == 141 ? 145 : 150));
            }
        }
    }
    SUBCASE("a much better lifted hit that misses a few clipped end bases replaces the primary hit") {
        // The read's first 4 bases differ from the ALT, so its ALT hit is
        // [4, 150) at 146; on chr1 it aligns [0, 150) at 120. Leaving out 4
        // bases costs 4: 146 > 124, so the lifted hit replaces it.
        for (int order = 0; order < 2; ++order) {
            CAPTURE(order);
            mem_alnreg_t l = reg(fb, ALT_FWD, 104, 250, 0, 146, 2, 4, 150), p = reg(fb, CHR1, 3100, 3250, 0, 120, 1);
            mem_alnreg_t a[2] = {order ? p : l, order ? l : p};
            REQUIRE(mem_alt_project(o.o, &fb.bns, 2, a) == 1);
            CHECK(a[0].seedcov == 2);
            CHECK(a[0].rid == CHR1);
            CHECK(a[0].rb == 3104);
            CHECK(a[0].re == 3250);
            CHECK(a[0].qb == 4);
            CHECK(a[0].score == 146);
        }
    }
    SUBCASE("a lifted hit that extends past the primary hit on one end is charged for the other") {
        // lifted [5, 150), primary [0, 145): one locus; the lifted hit leaves out
        // the primary hit's bases [0, 5), worth 5. 135 ties 130 + 5 and keeps the
        // primary hit; 136 replaces it.
        for (int sc = 135; sc <= 136; ++sc) {
            CAPTURE(sc);
            mem_alnreg_t a[2] = {reg(fb, ALT_FWD, 105, 250, 0, sc, 2, 5, 150),
                                 reg(fb, CHR1, 3100, 3245, 0, 130, 1, 0, 145)};
            REQUIRE(mem_alt_project(o.o, &fb.bns, 2, a) == 1);
            CHECK(a[0].seedcov == (sc == 136 ? 2 : 1));
            CHECK(a[0].qb == (sc == 136 ? 5 : 0));
        }
    }
    SUBCASE("a lifted hit covering a shorter primary hit replaces it") {
        // lifted [0, 150), primary [5, 150): the lifted hit covers it and scores more
        mem_alnreg_t a[2] = {reg(fb, CHR1, 3105, 3250, 0, 130, 1, 5, 150),
                             reg(fb, ALT_FWD, 100, 250, 0, 150, 2)};
        REQUIRE(mem_alt_project(o.o, &fb.bns, 2, a) == 1);
        CHECK(a[0].seedcov == 2);
        CHECK(a[0].qb == 0);
        CHECK(a[0].qe == 150);
        CHECK(a[0].rb == 3100);
    }
    SUBCASE("one primary hit it cannot replace keeps the lifted hit off the locus") {
        // P1 [0, 150) at 140 is covered and beaten; P2 [0, 160) at 145 is not
        // replaced (150 is not more than 145 + 10): the lifted hit replaces
        // neither, is dropped for P2, and both primary hits stay, as without the
        // flag
        mem_alnreg_t a[3] = {reg(fb, ALT_FWD, 100, 250, 0, 150, 2), reg(fb, CHR1, 3100, 3250, 0, 140, 1),
                             reg(fb, CHR1, 3100, 3258, 0, 145, 3, 0, 160)};
        REQUIRE(mem_alt_project(o.o, &fb.bns, 3, a) == 2);
        CHECK(a[0].seedcov == 3);
        CHECK(a[1].seedcov == 1);
    }
}

// Every input order of `in` must give `want` (tags, in output order).
static void check_all_orders(const mem_opt_t *opt, const bntseq_t *bns, std::vector<mem_alnreg_t> in,
                             const std::vector<int> &want) {
    std::vector<int> perm(in.size());
    for (size_t i = 0; i < perm.size(); ++i) perm[i] = (int) i;
    do {
        std::vector<mem_alnreg_t> a;
        for (int i : perm) a.push_back(in[i]);
        const int n = mem_alt_project(opt, bns, (int) a.size(), a.data());
        std::vector<int> got;
        for (int i = 0; i < n; ++i) got.push_back(a[i].seedcov);
        CAPTURE(perm);
        CHECK(got == want);
    } while (std::next_permutation(perm.begin(), perm.end()));
}

TEST_CASE("mem_alt_project: a lifted hit is dropped only for a hit that is kept" * doctest::test_suite("unit/alt_lift")) {
    FakeBns fb(kCtgs);
    TmpFile f(kAlt);
    LiftGuard g(alt_lift_load(f.path.c_str(), &fb.bns));
    REQUIRE(g.l != nullptr);
    OptGuard o;
    o.o->alt_lift = g.l;

    SUBCASE("a primary hit replaced by one lifted hit does not take down another lifted hit at a separate locus") {
        // L1 lifts to chr1 [3000, 3150) at 100, P is chr1 [3004, 3154) at 120 and
        // L2 lifts to chr1 [3008, 3158) at 130, all over read [0, 150). L1 and P
        // are one locus (146 > 0.95 * 150), so are P and L2, but L1 and L2 are not
        // (142). L1 cannot replace P; L2 can, and does. P is then gone, so L1 is
        // tested against nothing kept at its locus and stays: two placements, as
        // ordinary dedup would keep two hits 8 bases apart.
        check_all_orders(o.o, &fb.bns,
                         {reg(fb, ALT_FWD, 0, 150, 0, 100, 1), reg(fb, CHR1, 3004, 3154, 0, 120, 2),
                          reg(fb, ALT_FWD, 8, 158, 0, 130, 3)},
                         {3, 1});
    }
    SUBCASE("an unlifted hit is dropped only for a lifted hit that is kept") {
        // L1 lifts to chr1 [3000, 3150) at 150, L2 to chr1 [3004, 3154) at 145 and
        // P is chr1 [3008, 3158) at 130, all over read [0, 150). L1 and L2 are one
        // locus (146), so are L2 and P (146), but L1 and P are not (142). L2 could
        // replace P, but L1 is kept and drops L2, and no kept lifted hit is at P's
        // locus, so P stays.
        check_all_orders(o.o, &fb.bns,
                         {reg(fb, ALT_FWD, 0, 150, 0, 150, 1), reg(fb, ALT_FWD, 4, 154, 0, 145, 2),
                          reg(fb, CHR1, 3008, 3158, 0, 130, 3)},
                         {1, 3});
    }
    SUBCASE("a lifted hit that loses to a primary hit the other lifted hit replaced is kept") {
        // L1 lifts to chr1 [3100, 3240), read [0, 140), at 140; P is chr1 [3140,
        // 3240), read [40, 140), at 95; L2 lifts to chr1 [3140, 3280), read
        // [40, 180), at 90. L1 replaces P; L2 cannot (90 < 95), but P is gone and
        // L1 is not at L2's locus (100 of 140), so L2, the only hit for read bases
        // [140, 180), is kept.
        check_all_orders(o.o, &fb.bns,
                         {reg(fb, ALT_FWD, 100, 240, 0, 140, 1, 0, 140), reg(fb, CHR1, 3140, 3240, 0, 95, 2, 40, 140),
                          reg(fb, ALT_FWD, 140, 280, 0, 90, 3, 40, 180)},
                         {1, 3});
        // L3, a copy of L2 on alt_fwd2 at 89, cannot replace P either. L2 is kept
        // in step 2, so L3, at L2's locus, is dropped for it: a region kept
        // earlier in step 2 counts as kept.
        check_all_orders(o.o, &fb.bns,
                         {reg(fb, ALT_FWD, 100, 240, 0, 140, 1, 0, 140), reg(fb, CHR1, 3140, 3240, 0, 95, 2, 40, 140),
                          reg(fb, ALT_FWD, 140, 280, 0, 90, 3, 40, 180), reg(fb, ALT_FWD2, 140, 280, 0, 89, 4, 40, 180)},
                         {1, 3});
    }
    SUBCASE("a higher-scoring lifted hit that cannot replace the primary hit gives way to one that can") {
        // A 160-base read. L1 lifts to chr1 [3100, 3250), read [0, 150), at 150;
        // L2 (a reverse hit on the FLAG-16 alt_rev) lifts to chr1 [3100, 3260),
        // read [0, 160), at 145; P is chr1 [3100, 3260), read [0, 160), at 140.
        // L1 leaves out 10 of P's bases and 150 is not more than 140 + 10, so only
        // L2 can replace P. It does, and L1, at L2's locus, is dropped for it.
        check_all_orders(o.o, &fb.bns,
                         {reg(fb, ALT_FWD, 100, 250, 0, 150, 1, 0, 150), reg(fb, ALT_REV, 1737, 1897, 1, 145, 2, 0, 160),
                          reg(fb, CHR1, 3100, 3260, 0, 140, 3, 0, 160)},
                         {2});
        // without P, the two lifted hits are one locus and the better is kept
        check_all_orders(o.o, &fb.bns,
                         {reg(fb, ALT_FWD, 100, 250, 0, 150, 1, 0, 150), reg(fb, ALT_REV, 1737, 1897, 1, 145, 2, 0, 160)},
                         {1});
    }
    SUBCASE("of two tied lifted hits, one that can replace beats one that cannot, whatever the contig index") {
        // L1 (alt_fwd, rid 1) lifts to chr1 [3100, 3240), read [0, 140), at 140; L2
        // (alt_fwd2, rid 12) to chr1 [3100, 3250), read [0, 150), at 140; P is chr1
        // [3100, 3250) at 135. L1 and L2 tie on score, rb and qb. L1 leaves out 10
        // of P's bases and 140 is not more than 135 + 10, so only L2 can replace P:
        // L2 is kept and L1, at its locus, is dropped.
        check_all_orders(o.o, &fb.bns,
                         {reg(fb, ALT_FWD, 100, 240, 0, 140, 1, 0, 140), reg(fb, ALT_FWD2, 100, 250, 0, 140, 2),
                          reg(fb, CHR1, 3100, 3250, 0, 135, 3)},
                         {2});
        // without P both can replace, and the lower contig index wins
        check_all_orders(o.o, &fb.bns,
                         {reg(fb, ALT_FWD, 100, 240, 0, 140, 1, 0, 140), reg(fb, ALT_FWD2, 100, 250, 0, 140, 2)},
                         {1});
    }
    SUBCASE("a lifted hit that cannot replace is dropped for a kept primary hit even if it beats another") {
        // L lifts to chr1 [3100, 3250) at 150; P1 chr1 [3100, 3250) at 140 is
        // beaten and covered, P2 chr1 [3100, 3258), read [0, 160), at 145 is not
        // replaced. L replaces neither, so both primary hits stay.
        check_all_orders(o.o, &fb.bns,
                         {reg(fb, ALT_FWD, 100, 250, 0, 150, 1), reg(fb, CHR1, 3100, 3250, 0, 140, 2),
                          reg(fb, CHR1, 3100, 3258, 0, 145, 3, 0, 160)},
                         {3, 2});
    }
}

TEST_CASE("mem_alt_project: many regions, by-score order" * doctest::test_suite("unit/alt_lift")) {
    FakeBns fb(kCtgs);
    TmpFile f(kAlt);
    LiftGuard g(alt_lift_load(f.path.c_str(), &fb.bns));
    REQUIRE(g.l != nullptr);
    OptGuard o;
    o.o->alt_lift = g.l;
    // 69 non-overlapping primary hits in ascending score order plus one lifted
    // ALT hit: once a region lifts, the output is in by-score order (score
    // descending, then rb, then qb) whatever the input order. The lifted hit
    // ties the primary hit at chr1:48000 at 140 and sorts before it by rb.
    std::vector<mem_alnreg_t> a;
    for (int k = 0; k < 69; ++k)
        a.push_back(reg(fb, CHR1, 40000 + 200 * k, 40000 + 200 * k + 150, 0, 100 + k, 1));
    a.push_back(reg(fb, ALT_FWD, 100, 250, 0, 140, 2));
    const int n = mem_alt_project(o.o, &fb.bns, (int) a.size(), a.data());
    REQUIRE(n == 70);
    int n_lifted = 0, i_lifted = -1;
    for (int i = 0; i < n; ++i) {
        CAPTURE(i);
        CHECK(a[i].rid == CHR1);
        CHECK(int(a[i].is_alt) == 0);
        if (i) {
            const mem_alnreg_t &x = a[i - 1], &y = a[i];
            CHECK((x.score > y.score || (x.score == y.score && (x.rb < y.rb || (x.rb == y.rb && x.qb <= y.qb)))));
        }
        if (a[i].seedcov == 2) {
            ++n_lifted;
            i_lifted = i;
            CHECK(a[i].rb == 3100);
        }
    }
    CHECK(n_lifted == 1);
    // scores 168 down to 141 come first (28 hits), then the 140 tie
    REQUIRE(i_lifted == 28);
    CHECK(a[28].score == 140);
    CHECK(a[29].score == 140);
    CHECK(a[29].rb == 48000);
    CHECK(a[30].score == 139);
}

TEST_CASE("mem_alt_project: a read with no region to lift is left as it is" * doctest::test_suite("unit/alt_lift")) {
    FakeBns fb(kCtgs);
    TmpFile f(kAlt);
    LiftGuard g(alt_lift_load(f.path.c_str(), &fb.bns));
    REQUIRE(g.l != nullptr);
    OptGuard o;
    o.o->alt_lift = g.l;
    // primary hits out of score order and an ALT hit that stays (in the clip):
    // nothing lifts, so neither the order nor any field changes
    mem_alnreg_t a[3] = {reg(fb, CHR1, 3100, 3250, 0, 120, 1), reg(fb, CHR1, 3104, 3254, 0, 140, 2),
                         reg(fb, ALT_CLIP, 0, 100, 0, 100, 3, 0, 100)};
    mem_alnreg_t b[3];
    memcpy(b, a, sizeof(a));
    REQUIRE(mem_alt_project(o.o, &fb.bns, 3, a) == 3);
    CHECK(memcmp(a, b, sizeof(a)) == 0);
}

TEST_CASE("mem_alt_project: paired-confidence gap retains full candidate scores" * doctest::test_suite("unit/alt_lift")) {
    FakeBns fb({{"chr1", 20000, 0}, {"alt", 1000, 1}});
    TmpFile f("alt\t0\tchr1\t1001\t60\t1000M\t*\t0\t0\t*\t*\n");
    LiftGuard lift(alt_lift_load(f.path.c_str(), &fb.bns));
    OptGuard opt;
    opt.o->alt_lift = lift.l;
    opt.o->alt_project_min_gap = 2;
    opt.o->alt_project_pair_gap = 1;
    std::vector<mem_alnreg_t> regions[2] = {
        {reg(fb, 0, 1000, 1100, 0, 100, 100, 0, 100),
         reg(fb, 1, 500, 600, 0, 100, 100, 0, 100),
         reg(fb, 0, 1500, 1600, 0, 90, 100, 0, 100)},
        {reg(fb, 0, 1200, 1300, 1, 100, 100, 0, 100),
         reg(fb, 1, 700, 800, 1, 100, 100, 0, 100),
         reg(fb, 0, 1700, 1800, 1, 90, 100, 0, 100)}
    };
    SUBCASE("an equal paired alternative is not weakened by ordinary hits") {
        const auto out = project_and_pair(opt.o, fb, regions);
        REQUIRE(out.paired == 1);
        CHECK(out.q[0] == 0);
        CHECK(out.q[1] == 0);
    }
    SUBCASE("a three-point pair advantage minus gap two caps MAPQ at six") {
        regions[0][0] = reg(fb, 1, 0, 100, 0, 100, 100, 0, 100);
        regions[1][1].score = regions[1][1].truesc = 97;
        const auto out = project_and_pair(opt.o, fb, regions);
        REQUIRE(out.paired == 1);
        CHECK(out.rb[0] == 1000);
        CHECK(out.q[0] <= 6);
        CHECK(out.q[1] <= 6);
        CHECK(out.q[0] > 0);
    }
    SUBCASE("match-score scaling preserves the confidence margin") {
        regions[0][0] = reg(fb, 1, 0, 100, 0, 200, 100, 0, 100);
        regions[0][1].score = regions[0][1].truesc = 200;
        regions[0][2].score = regions[0][2].truesc = 180;
        regions[1][0].score = regions[1][0].truesc = 200;
        regions[1][1].score = regions[1][1].truesc = 194;
        regions[1][2].score = regions[1][2].truesc = 180;
        opt.o->a = 2;
        opt.o->pen_unpaired *= 2;
        const auto out = project_and_pair(opt.o, fb, regions);
        REQUIRE(out.paired == 1);
        CHECK(out.rb[0] == 1000);
        CHECK(out.q[0] > 0);
        CHECK(out.q[0] <= 6);
        CHECK(out.q[1] <= 6);
    }
    SUBCASE("a split-query primary bypasses paired confidence") {
        regions[0] = {reg(fb, 1, 0, 60, 0, 60, 60, 0, 60),
                      reg(fb, 0, 10000, 10040, 0, 40, 40, 60, 100)};
        regions[1].resize(1);
        const auto out = project_and_pair(opt.o, fb, regions);
        CHECK(out.paired == 0);
    }
    SUBCASE("a unique mate still resolves a single-end tie") {
        regions[1].resize(1);
        const auto out = project_and_pair(opt.o, fb, regions);
        REQUIRE(out.paired == 1);
        CHECK(out.rb[0] == 1000);
        CHECK(out.q[0] > 0);
    }
    SUBCASE("a pair with no selected lifted hit keeps ordinary confidence") {
        regions[0].resize(1);
        regions[1].resize(1);
        opt.o->alt_project_min_gap = 2147483647;
        const auto out = project_and_pair(opt.o, fb, regions);
        REQUIRE(out.paired == 1);
        CHECK(out.rb[0] == 1000);
        CHECK(out.q[0] == 60);
        CHECK(out.q[1] == 60);
    }
    SUBCASE("unpaired-preferred ends retain independent confidence") {
        regions[0] = {reg(fb, 1, 0, 100, 0, 100, 100, 0, 100)};
        regions[1] = {reg(fb, 0, 1317, 1417, 1, 100, 100, 0, 100)};
        opt.o->pen_unpaired = 0;
        const auto out = project_and_pair(opt.o, fb, regions);
        REQUIRE(out.paired == 1);
        CHECK((out.flag & 2) == 0);
        CHECK(out.q[0] == 60);
        CHECK(out.q[1] == 60);
    }
    SUBCASE("no-pairing and failed-pair paths do not read selected indices") {
        regions[0] = {reg(fb, 1, 0, 100, 0, 100, 100, 0, 100)};
        SUBCASE("no pairing requested") { opt.o->flag |= MEM_F_NOPAIRING; }
        SUBCASE("no compatible mate") { regions[1] = {reg(fb, 0, 10000, 10100, 1, 100, 100, 0, 100)}; }
        const auto out = project_and_pair(opt.o, fb, regions);
        CHECK(out.paired == 0);
    }
    SUBCASE("a wide gap cannot overflow or change the chosen placement") {
        opt.o->a = 127;
        regions[0][0] = reg(fb, 1, 0, 100, 0, 100, 100, 0, 100);
        regions[1].resize(1);
        opt.o->alt_project_min_gap = 2147483647;
        const auto out = project_and_pair(opt.o, fb, regions);
        REQUIRE(out.paired == 1);
        CHECK(out.rb[0] == 1000);
        CHECK(out.q[0] == 0);
        CHECK(out.q[1] == 0);
    }
}

TEST_CASE("mem_alt_project: paired-confidence gap only lowers confidence over varied candidates" * doctest::test_suite("unit/alt_lift")) {
    FakeBns fb({{"chr1", 20000, 0}, {"alt", 4000, 1}});
    TmpFile f("alt\t0\tchr1\t1001\t60\t4000M\t*\t0\t0\t*\t*\n");
    LiftGuard lift(alt_lift_load(f.path.c_str(), &fb.bns));
    OptGuard opt;
    opt.o->alt_lift = lift.l;
    opt.o->alt_project_pair_gap = 1;
    std::mt19937 random(775);
    int lowered = 0;
    for (int trial = 0; trial < 200; ++trial) {
        std::vector<mem_alnreg_t> regions[2];
        opt.o->a = 1 + (trial % 3);
        opt.o->pen_unpaired = (trial % 4 == 0 ? 0 : 17 * opt.o->a);
        for (int locus = 0; locus < 3; ++locus) {
            const int pos = 1000 + locus * 1000;
            for (int end = 0; end < 2; ++end) {
                const bool lifted = (random() % 2) != 0;
                const int start = pos + (end ? 200 : 0) + (int)(random() % 7);
                const int score = (80 + (int)(random() % 21)) * opt.o->a;
                regions[end].push_back(reg(fb, lifted ? 1 : 0,
                    lifted ? start - 1000 : start, lifted ? start - 900 : start + 100,
                    end, score, 100, 0, 100));
            }
        }
        opt.o->alt_project_min_gap = 0;
        const auto baseline = project_and_pair(opt.o, fb, regions);
        opt.o->alt_project_min_gap = 1 + trial % 5;
        const auto guarded = project_and_pair(opt.o, fb, regions);
        CHECK(guarded.paired == baseline.paired);
        CHECK(guarded.flag == baseline.flag);
        if (!baseline.paired) continue;
        for (int end = 0; end < 2; ++end) {
            CHECK(guarded.rb[end] == baseline.rb[end]);
            CHECK(guarded.score[end] == baseline.score[end]);
            CHECK(guarded.q[end] <= baseline.q[end]);
            lowered += guarded.q[end] < baseline.q[end];
        }
    }
    CHECK(lowered > 0);
}

TEST_CASE("mem_alt_project: score-gap guard refuses ambiguous cross-locus projections" * doctest::test_suite("unit/alt_lift")) {
    FakeBns fb(kCtgs);
    TmpFile tmp(kAlt);
    LiftGuard l(alt_lift_load(tmp.path.c_str(), &fb.bns));
    OptGuard o;
    o.o->alt_lift = l.l;
    o.o->alt_project_min_gap = 2;
    // The ALT projects to chr1:3000, while the same query also matches chr1:10000.
    // Removing the gate, changing <= to <, or ignoring -A must break this test.
    for (int match : {1, 2}) {
        o.o->a = match;
        for (int competitor : {150 - 2 * match - 1, 150 - 2 * match, 150, 151}) {
            for (int order = 0; order < 2; ++order) {
                mem_alnreg_t a[] = {reg(fb, ALT_FWD, 0, 150, 0, 150, 1),
                                    reg(fb, CHR1, 10000, 10150, 0, competitor, 2)};
                a[0].seedlen0 = 19;
                const auto original = a[0];
                if (order) std::swap(a[0], a[1]);
                REQUIRE(mem_alt_project(o.o, &fb.bns, 2, a) == 2);
                const auto &lift = a[a[0].seedcov == 1 ? 0 : 1];
                if (competitor >= 150 - 2 * match) CHECK(memcmp(&lift, &original, sizeof(lift)) == 0);
                CHECK(lift.rid == (competitor >= 150 - 2 * match ? ALT_FWD : CHR1));
                CHECK(lift.is_alt == (competitor >= 150 - 2 * match));
            }
        }
    }
}

TEST_CASE("mem_alt_project: gap decisions use projected competitors without suppressing split pieces" * doctest::test_suite("unit/alt_lift")) {
    FakeBns fb(kCtgs);
    TmpFile tmp(kAlt);
    LiftGuard l(alt_lift_load(tmp.path.c_str(), &fb.bns));
    OptGuard o;
    o.o->alt_lift = l.l;
    o.o->alt_project_min_gap = 2;
    SUBCASE("same-locus hits are not competitors") {
        mem_alnreg_t a[] = {reg(fb, ALT_FWD, 0, 150, 0, 150, 1), reg(fb, CHR1, 3000, 3150, 0, 149, 2)};
        REQUIRE(mem_alt_project(o.o, &fb.bns, 2, a) == 1);
        CHECK(a[0].seedcov == 1);
        CHECK(a[0].rid == CHR1);
    }
    SUBCASE("disjoint query pieces do not veto a lift") {
        mem_alnreg_t a[] = {reg(fb, ALT_FWD, 0, 100, 0, 100, 1, 0, 100),
                            reg(fb, CHR1, 10000, 10100, 0, 100, 2, 100, 200)};
        REQUIRE(mem_alt_project(o.o, &fb.bns, 2, a) == 2);
        CHECK(a[0].rid == CHR1);
        CHECK(a[1].rid == CHR1);
    }
    SUBCASE("equal projected competitors both refuse in every order") {
        for (int order = 0; order < 2; ++order) {
            mem_alnreg_t a[] = {reg(fb, ALT_FWD, 0, 150, 0, 150, 1),
                                reg(fb, ALT_FWD2, 500, 650, 0, 150, 2)};
            if (order) std::swap(a[0], a[1]);
            REQUIRE(mem_alt_project(o.o, &fb.bns, 2, a) == 2);
            CHECK((int) a[0].is_alt == 1);
            CHECK((int) a[1].is_alt == 1);
        }
    }
    SUBCASE("unprojectable ALT cannot veto a valid lift") {
        mem_alnreg_t a[] = {reg(fb, ALT_FWD, 0, 150, 0, 150, 1),
                            reg(fb, ALT_CLIP, 0, 150, 0, 150, 2)};
        REQUIRE(mem_alt_project(o.o, &fb.bns, 2, a) == 2);
        CHECK(a[0].rid == CHR1);
        CHECK(a[1].rid == ALT_CLIP);
    }
    SUBCASE("opposite-strand placements compete") {
        mem_alnreg_t a[] = {reg(fb, ALT_FWD, 0, 150, 0, 150, 1),
                            reg(fb, CHR1, 3000, 3150, 1, 149, 2)};
        REQUIRE(mem_alt_project(o.o, &fb.bns, 2, a) == 2);
        const auto &lift = a[a[0].seedcov == 1 ? 0 : 1];
        CHECK(lift.rid == ALT_FWD);
    }
    SUBCASE("half query overlap is not a competitor; one base more is") {
        for (int begin : {74, 75}) {
            mem_alnreg_t a[] = {reg(fb, ALT_FWD, 0, 150, 0, 150, 1),
                                reg(fb, CHR1, 10000, 10150, 0, 149, 2, begin, begin + 150)};
            REQUIRE(mem_alt_project(o.o, &fb.bns, 2, a) == 2);
            const auto &lift = a[a[0].seedcov == 1 ? 0 : 1];
            CHECK(lift.rid == (begin == 74 ? ALT_FWD : CHR1));
        }
    }
    SUBCASE("wide products do not wrap the margin") {
        o.o->a = 127;
        o.o->alt_project_min_gap = 2147483647;
        mem_alnreg_t a[] = {reg(fb, ALT_FWD, 0, 150, 0, 150, 1), reg(fb, CHR1, 10000, 10150, 0, 1, 2)};
        REQUIRE(mem_alt_project(o.o, &fb.bns, 2, a) == 2);
        const auto &lift = a[a[0].seedcov == 1 ? 0 : 1];
        CHECK(lift.rid == ALT_FWD);
    }
    SUBCASE("zero disables even an exact cross-locus tie") {
        o.o->alt_project_min_gap = 0;
        mem_alnreg_t a[] = {reg(fb, ALT_FWD, 0, 150, 0, 150, 1), reg(fb, CHR1, 10000, 10150, 0, 150, 2)};
        REQUIRE(mem_alt_project(o.o, &fb.bns, 2, a) == 2);
        CHECK(a[0].rid == CHR1);
        CHECK(a[1].rid == CHR1);
    }
}

TEST_CASE("mem_alt_project: provenance follows the winning source contig" * doctest::test_suite("unit/alt_lift")) {
    FakeBns fb(kCtgs);
    TmpFile tmp(kAlt);
    LiftGuard l(alt_lift_load(tmp.path.c_str(), &fb.bns));
    OptGuard o;
    o.o->alt_lift = l.l;
    for (int rid : {ALT_FWD, ALT_REV, ALT_FWD2}) {
        for (int strand : {0, 1}) {
            auto a = reg(fb, rid, 0, 150, strand, 150, 1);
            a.seedlen0 = 19;
            REQUIRE(mem_alt_project(o.o, &fb.bns, 1, &a) == 1);
            REQUIRE(a.rid == CHR1);
            CHECK(a.alt_source_rid == rid);
        }
    }
}

TEST_CASE("mem_alt_project: without --alt-project it is a no-op" * doctest::test_suite("unit/alt_lift")) {
    FakeBns fb(kCtgs);
    OptGuard o;
    REQUIRE(o.o->alt_lift == nullptr);
    // an ALT hit that would lift and a primary hit at the same locus: with the
    // option off neither is lifted or dropped, and no field changes
    mem_alnreg_t a[2] = {reg(fb, ALT_FWD, 100, 250, 0, 150, 1), reg(fb, CHR1, 3100, 3250, 0, 140, 2)};
    mem_alnreg_t b[2];
    memcpy(b, a, sizeof(a));
    REQUIRE(mem_alt_project(o.o, &fb.bns, 2, a) == 2);
    CHECK(memcmp(a, b, sizeof(a)) == 0);
}


// Omitting exclusion, using aligned blocks, or filtering only after first-wins
// selection must break these hand-derived tables.
TEST_CASE("exclude BED footprint boundaries and CIGAR" * doctest::test_suite("unit/alt_lift")) {
    struct Case { const char *cigar; int pos, qlen; bool excluded; };
    const Case cases[] = {
        {"100M", 1, 100, false}, {"100M", 100, 100, true},
        {"100M", 201, 100, false}, {"100M", 100, 100, true},
        {"50M100D50M", 1, 100, true}, {"50M100N50M", 1, 100, true},
        {"50=1X49=", 100, 100, true}, {"50M100I50M", 1, 200, false},
        {"100S100M", 1, 200, false}, {"100H100M", 1, 200, false},
        {"50M100P50M", 1, 100, false}, {"100S", 101, 100, false},
        {"100M", 901, 100, false},
    };
    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); ++i) {
        const auto &c = cases[i];
        CAPTURE(i);
        FakeBns fb({{"chr1", 1000, 0}, {"tile", c.qlen, 1}});
        // D/N cases intersect only the deleted/skipped portion.
        TmpFile bed(i == 4 || i == 5 ? "chr1\t75\t76\n" : "chr1\t100\t200\n");
        const int flag = i == 3 ? 16 : 0;
        const std::string suffix = "\tchr1\t" + std::to_string(c.pos) + "\t60\t" + c.cigar + "\t*\t0\t0\t*\t*\n";
        TmpFile source(("tile\t" + std::to_string(flag) + suffix).c_str());
        TmpFile masked(("tile\t" + std::to_string(flag | (c.excluded ? 4 : 0)) + suffix).c_str());
        alt_lift_exclude_stats_t stats = {};
        char error[1024] = {};
        LiftGuard actual(alt_lift_load_excluded(source.path.c_str(), &fb.bns, bed.path.c_str(), &stats, error, sizeof(error)));
        LiftGuard expected(alt_lift_load(masked.path.c_str(), &fb.bns));
        CAPTURE(error);
        check_lift_tables(actual.l, expected.l);
        CHECK(stats.input_intervals == 1);
        CHECK(stats.merged_intervals == 1);
        CHECK(stats.zero_width_intervals == 0);
        CHECK(stats.excluded == (c.excluded ? 1 : 0));
        CHECK(fb.bns.anns[1].is_alt == 1);
        if (!c.excluded) {
            REQUIRE(actual.l->n_lift == 1);
            CHECK(actual.l->a[1].tpos == c.pos - 1);
            CHECK(actual.l->a[1].rev == (flag == 16));
        }
    }
}

TEST_CASE("exclude BED duplicate fallback equals FLAG masking" * doctest::test_suite("unit/alt_lift")) {
    FakeBns fb({{"chr1", 1000, 0}, {"tile", 100, 1}});
    TmpFile bed("chr1\t100\t200\n");
    const char *source_cases[] = {
        "tile\t0\tchr1\t101\t60\t100M\ntile\t16\tchr1\t301\t60\t100M\n",
        "tile\t16\tchr1\t301\t60\t100M\ntile\t0\tchr1\t101\t60\t100M\n",
        "tile\t0\tchr1\t101\t60\t100M\ntile\t16\tchr1\t101\t60\t100M\n",
    };
    const char *masked_cases[] = {
        "tile\t4\tchr1\t101\t60\t100M\ntile\t16\tchr1\t301\t60\t100M\n",
        "tile\t16\tchr1\t301\t60\t100M\ntile\t4\tchr1\t101\t60\t100M\n",
        "tile\t4\tchr1\t101\t60\t100M\ntile\t20\tchr1\t101\t60\t100M\n",
    };
    for (int i = 0; i < 3; ++i) {
        CAPTURE(i);
        const std::string prefix =
            "tile\t4\tchr1\t101\t60\tbad\n"
            "tile\t0\t*\t0\t60\t*\n"
            "tile\t256\tchr1\t101\t60\t100M\n"
            "tile\t2048\tchr1\t101\t60\t100M\n"
            "tile\t0\tchr1\t101\t60\t50M\n"
            "tile\t0\tchr1\t101\t60\tbad\n";
        TmpFile source((prefix + source_cases[i]).c_str());
        TmpFile masked((prefix + masked_cases[i]).c_str());
        alt_lift_exclude_stats_t stats = {};
        char error[1024] = {};
        LiftGuard actual(alt_lift_load_excluded(source.path.c_str(), &fb.bns, bed.path.c_str(), &stats, error, sizeof(error)));
        LiftGuard expected(alt_lift_load(masked.path.c_str(), &fb.bns));
        check_lift_tables(actual.l, expected.l);
        CHECK(stats.excluded == (i == 0 ? 1 : i == 1 ? 0 : 2));
        CHECK(actual.l->skip.bad == 2);
        CHECK(actual.l->skip.unmapped == 2);
        CHECK(actual.l->skip.not_primary == 2);
        CHECK(actual.l->n_lift == (i == 2 ? 0 : 1));
        if (i != 2) CHECK(actual.l->a[1].tpos == 300);
        CHECK(fb.bns.anns[1].is_alt == 1);
    }
}

TEST_CASE("exclude BED reads every member of a multi-member gzip" * doctest::test_suite("unit/alt_lift")) {
    // bgzip output is a run of concatenated gzip members; build the same
    // structure with zlib (append mode starts a new member).
    FakeBns fb({{"chr1", 1000, 0}, {"tileA", 100, 1}, {"tileB", 100, 1}});
    TmpFile source("tileA\t0\tchr1\t101\t60\t100M\ntileB\t0\tchr1\t501\t60\t100M\n");
    TmpFile multi("");
    write_gzip(multi.path, "chr1\t150\t160\n");
    {
        gzFile fp = gzopen(multi.path.c_str(), "ab");
        REQUIRE(fp != nullptr);
        const std::string second = "chr1\t550\t560\n";
        CHECK(gzwrite(fp, second.data(), (unsigned) second.size()) == (int) second.size());
        REQUIRE(gzclose(fp) == Z_OK);
    }
    alt_lift_exclude_stats_t stats = {};
    char error[1024] = {'x'};
    LiftGuard actual(alt_lift_load_excluded(source.path.c_str(), &fb.bns, multi.path.c_str(), &stats, error, sizeof(error)));
    REQUIRE(actual.l != nullptr);
    CHECK(error[0] == 0);
    CHECK(stats.input_intervals == 2);
    CHECK(stats.merged_intervals == 2);
    CHECK(stats.excluded == 2);
    CHECK(actual.l->n_lift == 0);
}

TEST_CASE("exclude BED dictionary and compressed errors" * doctest::test_suite("unit/alt_lift")) {
    FakeBns fb({{"chr1", 1000, 0}, {"tile", 100, 1}, {"track_chr", 1000, 0}, {"browser_chr", 1000, 0}});
    TmpFile source("tile\t0\tchr1\t101\t60\t100M\n");
    const std::string text = "# comment\r\ntrack name=exclude\nbrowser position chr1\n\nchr1\t150\t200\textra\r\nchr1\t100\t150\nchr1\t125\t175\nchr1\t220\t220\ntrack_chr\t0\t1\nbrowser_chr\t0\t1";
    TmpFile plain(text.c_str()), compressed("");
    write_gzip(compressed.path, text);
    for (const auto &path : {plain.path, compressed.path}) {
        alt_lift_exclude_stats_t stats = {};
        char error[1024] = {'x'};
        LiftGuard actual(alt_lift_load_excluded(source.path.c_str(), &fb.bns, path.c_str(), &stats, error, sizeof(error)));
        REQUIRE(actual.l != nullptr);
        CHECK(actual.l->n_lift == 0);
        CHECK(error[0] == 0);
        CHECK(stats.input_intervals == 5);
        CHECK(stats.merged_intervals == 3);
        CHECK(stats.zero_width_intervals == 1);
        CHECK(stats.excluded == 1);
    }
    const char *bad_rows[] = {
        "chr1\t-1\t10", "chr1\t+1\t10", "chr1\t1x\t10", "chr1\t20\t10",
        "chr1\t1\t1001", "chr1\t1\t9223372036854775808", "chr1\t1", "chr1\t\t10",
        "chrZ\t1\t10", "tile\t1\t10", "chr1 1 10", "chr1\t1\t10x",
        "track_bad\t1\t10", "browser_bad\t1\t10",
    };
    for (const auto row : bad_rows) {
        CAPTURE(row);
        for (bool gzip : {false, true}) {
            const std::string content = std::string("# comment\n") + row + "\n";
            TmpFile bad(content.c_str());
            if (gzip) write_gzip(bad.path, content);
            alt_lift_exclude_stats_t stats = {};
            char error[1024] = {};
            LiftGuard actual(alt_lift_load_excluded(source.path.c_str(), &fb.bns, bad.path.c_str(), &stats, error, sizeof(error)));
            CHECK(actual.l == nullptr);
            CHECK(std::string(error).find(bad.path) != std::string::npos);
            CHECK(std::string(error).find("line 2") != std::string::npos);
        }
    }
    SUBCASE("truncated gzip after valid BED rows fails") {
        // Missing trailer can surface only after the last complete text line.
        struct stat info = {};
        REQUIRE(stat(compressed.path.c_str(), &info) == 0);
        REQUIRE(info.st_size > 8);
        REQUIRE(truncate(compressed.path.c_str(), info.st_size - 8) == 0);
        alt_lift_exclude_stats_t stats = {};
        char error[1024] = {};
        LiftGuard actual(alt_lift_load_excluded(source.path.c_str(), &fb.bns, compressed.path.c_str(), &stats, error, sizeof(error)));
        CHECK(actual.l == nullptr);
        CHECK(std::string(error).find(compressed.path) != std::string::npos);
    }
    SUBCASE("missing BED and missing sidecar diagnose correct path") {
        for (bool missing_bed : {false, true}) {
            const char *missing = "/nonexistent/alt_exclude_test";
            alt_lift_exclude_stats_t stats = {};
            char error[32];
            memset(error, 'x', sizeof(error));
            LiftGuard actual(alt_lift_load_excluded(missing_bed ? source.path.c_str() : missing, &fb.bns,
                                                  missing_bed ? missing : plain.path.c_str(), &stats, error, sizeof(error)));
            CHECK(actual.l == nullptr);
            CHECK(error[sizeof(error) - 1] == 0);
            CHECK(std::string(error).find(missing) != std::string::npos);
        }
    }
    SUBCASE("empty and no-effect BED remain inspectable by the loader") {
        for (const char *content : {"# empty\nchr1\t10\t10\n", "chr1\t500\t600\n"}) {
            TmpFile bed(content);
            alt_lift_exclude_stats_t stats = {};
            char error[1024] = {};
            LiftGuard actual(alt_lift_load_excluded(source.path.c_str(), &fb.bns, bed.path.c_str(), &stats, error, sizeof(error)));
            REQUIRE(actual.l != nullptr);
            CHECK(actual.l->n_lift == 1);
            CHECK(stats.excluded == 0);
        }
    }
}

namespace {
struct BnsGuard {
    bntseq_t *b;
    explicit BnsGuard(bntseq_t *x) : b(x) {}
    ~BnsGuard() { if (b) bns_destroy(b); }
    BnsGuard(const BnsGuard &) = delete;
    BnsGuard &operator=(const BnsGuard &) = delete;
};

void write_text(const std::string &path, const std::string &text) {
    FILE *fp = fopen(path.c_str(), "wb");
    REQUIRE(fp != nullptr);
    CHECK(fwrite(text.data(), 1, text.size(), fp) == text.size());
    REQUIRE(fclose(fp) == 0);
}

std::vector<std::vector<std::string>> sam_fields(const std::string &text) {
    std::vector<std::vector<std::string>> rows;
    size_t at = 0;
    while (at < text.size()) {
        size_t end = text.find('\n', at);
        if (end == std::string::npos) end = text.size();
        std::vector<std::string> row;
        size_t f = at;
        while (true) {
            const size_t tab = text.find('\t', f);
            if (tab == std::string::npos || tab > end) { row.push_back(text.substr(f, end - f)); break; }
            row.push_back(text.substr(f, tab - f));
            f = tab + 1;
        }
        rows.push_back(row);
        at = end + 1;
    }
    return rows;
}
} // namespace

// The real index loader marks ALT contigs from the sidecar's QNAMEs alone
// (bns_restore), so neither runtime exclusion nor FLAG-masking the excluded line
// may change which contigs are ALT. Generated index; one ALT name has the
// conventional alt_ prefix and one (an HLA-style name) has no prefix or suffix,
// so prefix-based filtering would be caught. Each tile is excluded in turn.
TEST_CASE("exclude BED preserves real ALT marking" * doctest::test_suite("unit/alt_lift")) {
    std::mt19937 rng(775);
    std::string chr1(2000, 'A');
    for (char &c : chr1) c = "ACGT"[rng() % 4];
    const std::string tile_a = chr1.substr(100, 200);
    const std::string tile_b = bwa_tests::revcomp(chr1.substr(1000, 200));
    const char *name_a = "alt_chr1_tile";
    const char *name_b = "HLA-A*01:01:01:01";
    // The index is built and loaded before any sidecar exists, so its own bns
    // marks nothing: every ALT bit below comes from a sidecar's names.
    bwa_tests::TempIndex idx({{"chr1", chr1}, {name_a, tile_a}, {name_b, tile_b}});
    for (int rid = 0; rid < idx.bns()->n_seqs; ++rid) CHECK(idx.bns()->anns[rid].is_alt == 0);
    const std::string prefix = idx.prefix();
    const std::string alt_path = prefix + ".alt";
    const std::string source = std::string(name_a) + "\t0\tchr1\t101\t60\t200M\t*\t0\t0\t*\t*\n" +
                               name_b + "\t16\tchr1\t1001\t60\t200M\t*\t0\t0\t*\t*\n";
    for (int excluded = 1; excluded <= 2; ++excluded) {
        CAPTURE(excluded);
        // One base inside the excluded tile's footprint (BED is 0-based).
        TmpFile bed(excluded == 1 ? "chr1\t250\t251\n" : "chr1\t1000\t1001\n");
        // The oracle: the excluded line's FLAG gains 0x4, nothing else changes.
        std::string masked = source;
        const std::string flag_from = excluded == 1 ? std::string(name_a) + "\t0\t" : std::string(name_b) + "\t16\t";
        const std::string flag_to = excluded == 1 ? std::string(name_a) + "\t4\t" : std::string(name_b) + "\t20\t";
        REQUIRE(masked.find(flag_from) != std::string::npos);
        masked.replace(masked.find(flag_from), flag_from.size(), flag_to);
        const auto src_rows = sam_fields(source), mask_rows = sam_fields(masked);
        REQUIRE(src_rows.size() == 2);
        REQUIRE(mask_rows.size() == 2);
        for (size_t i = 0; i < 2; ++i) {
            CAPTURE(i);
            CHECK(src_rows[i][0] == mask_rows[i][0]);  // same QNAME sequence
            CHECK((src_rows[i][1] != mask_rows[i][1]) == ((int) i + 1 == excluded));
            for (size_t f = 2; f < src_rows[i].size(); ++f) CHECK(src_rows[i][f] == mask_rows[i][f]);
        }

        write_text(alt_path, source);
        BnsGuard src_bns(bns_restore(prefix.c_str()));
        REQUIRE(src_bns.b != nullptr);
        REQUIRE(src_bns.b->n_seqs == 3);
        CHECK(std::string(src_bns.b->anns[0].name) == "chr1");
        CHECK(std::string(src_bns.b->anns[1].name) == name_a);
        CHECK(std::string(src_bns.b->anns[2].name) == name_b);
        const int want_alt[3] = {0, 1, 1};
        for (int rid = 0; rid < 3; ++rid) CHECK(src_bns.b->anns[rid].is_alt == want_alt[rid]);

        // Ordinary projection (no BED) lifts both and leaves the bits alone.
        LiftGuard off(alt_lift_load(alt_path.c_str(), src_bns.b));
        REQUIRE(off.l != nullptr);
        CHECK(off.l->n_lift == 2);
        CHECK(off.l->a[1].trid == 0);
        CHECK(off.l->a[2].trid == 0);
        alt_lift_exclude_stats_t stats = {};
        char error[1024] = {};
        LiftGuard runtime(alt_lift_load_excluded(alt_path.c_str(), src_bns.b, bed.path.c_str(), &stats, error, sizeof(error)));
        CAPTURE(error);
        REQUIRE(runtime.l != nullptr);
        CHECK(stats.excluded == 1);
        CHECK(runtime.l->n_lift == 1);
        CHECK(runtime.l->a[excluded].trid == -1);
        CHECK(runtime.l->a[3 - excluded].trid == 0);
        for (int rid = 0; rid < 3; ++rid) CHECK(src_bns.b->anns[rid].is_alt == want_alt[rid]);

        write_text(alt_path, masked);
        BnsGuard mask_bns(bns_restore(prefix.c_str()));
        REQUIRE(mask_bns.b != nullptr);
        REQUIRE(mask_bns.b->n_seqs == src_bns.b->n_seqs);
        for (int rid = 0; rid < 3; ++rid) {
            CAPTURE(rid);
            CHECK(std::string(mask_bns.b->anns[rid].name) == src_bns.b->anns[rid].name);
            CHECK(mask_bns.b->anns[rid].len == src_bns.b->anns[rid].len);
            CHECK(mask_bns.b->anns[rid].offset == src_bns.b->anns[rid].offset);
            CHECK(mask_bns.b->anns[rid].is_alt == src_bns.b->anns[rid].is_alt);
        }
        LiftGuard oracle(alt_lift_load(alt_path.c_str(), mask_bns.b));
        check_lift_tables(runtime.l, oracle.l);
        CHECK(oracle.l->skip.unmapped == 1);
        for (int rid = 0; rid < 3; ++rid) CHECK(mask_bns.b->anns[rid].is_alt == want_alt[rid]);
    }
}

// ---------------------------------------------------------------------------
// Table gate: restricted loading versus a FLAG-masked sidecar.
//
// The masked sidecar is the source with FLAG 0x4 OR-ed into every record whose
// primary footprint overlaps the BED, as an independent overlap oracle writes
// it. Loading source + BED must give the same lift table, field by field, as
// loading the masked sidecar without a BED, against a bns that each sidecar's
// own names marked. The gate also checks the inputs it relies on to equate
// per-record decisions with table entries: identical QNAME order and non-FLAG
// bytes, unique QNAMEs, and every source record usable when loaded unfiltered.
// ---------------------------------------------------------------------------
namespace {
struct AltGateInputs {
    std::string prefix, source, masked, bed;
};

struct AltGateCounts {
    long n_seqs = 0, alt_contigs = 0, rows = 0, excluded = 0, retained = 0;
    long input_intervals = 0, merged_intervals = 0;
};

struct SidecarRow {
    std::string qname;
    long flag;
    std::string rest;  // every byte after FLAG up to the newline
};

bool file_exists(const std::string &path) {
    struct stat st;
    return stat(path.c_str(), &st) == 0 && S_ISREG(st.st_mode);
}

std::vector<SidecarRow> read_sidecar(const std::string &path) {
    std::ifstream in(path, std::ios::binary);
    REQUIRE(in.good());
    std::vector<SidecarRow> rows;
    std::string line;
    while (std::getline(in, line)) {
        const size_t t1 = line.find('\t');
        const size_t t2 = t1 == std::string::npos ? t1 : line.find('\t', t1 + 1);
        REQUIRE(t2 != std::string::npos);
        REQUIRE(line[0] != '@');
        char *end = nullptr;
        const std::string flag = line.substr(t1 + 1, t2 - t1 - 1);
        const long f = strtol(flag.c_str(), &end, 10);
        REQUIRE((!flag.empty() && *end == '\0'));
        rows.push_back({line.substr(0, t1), f, line.substr(t2)});
    }
    return rows;
}

// bns_restore a private prefix whose .ann/.amb/.pac link to `prefix` and whose
// .alt links to `alt`, so ALT marking comes from that sidecar's names and the
// shared index directory is never written.
bntseq_t *restore_with_alt(const std::string &dir, const char *tag, const std::string &prefix,
                           const std::string &alt) {
    const std::string base = dir + "/" + tag;
    auto link = [](const std::string &target, const std::string &at) {
        char *abs = realpath(target.c_str(), nullptr);
        REQUIRE(abs != nullptr);
        const int rc = symlink(abs, at.c_str());
        free(abs);
        REQUIRE(rc == 0);
    };
    for (const char *ext : {".ann", ".amb", ".pac"}) link(prefix + ext, base + ext);
    link(alt, base + ".alt");
    return bns_restore(base.c_str());
}

// Compare contig dictionaries and ALT bits. Returns the number of differing contigs.
long count_bns_differences(const bntseq_t *a, const bntseq_t *b) {
    long diff = 0;
    for (int rid = 0; rid < a->n_seqs; ++rid) {
        const bntann1_t &x = a->anns[rid], &y = b->anns[rid];
        if (strcmp(x.name, y.name) != 0 || x.len != y.len || x.offset != y.offset || x.is_alt != y.is_alt) {
            if (diff == 0) MESSAGE("first differing contig: rid " << rid << " " << x.name);
            ++diff;
        }
    }
    return diff;
}

// The whole gate. `tamper`, if given, edits the oracle table after loading; the
// generated cases use it to prove a single wrong field is caught.
AltGateCounts run_alt_table_gate(const AltGateInputs &in,
                                 const std::function<void(alt_lift_t *)> &tamper = nullptr) {
    for (const char *ext : {".ann", ".amb", ".pac"}) {
        INFO("missing index file " << in.prefix << ext);
        REQUIRE(file_exists(in.prefix + ext));
    }
    for (const auto &path : {in.source, in.masked, in.bed}) {
        INFO("missing gate input " << path);
        REQUIRE(file_exists(path));
    }
    AltGateCounts counts;

    // Real ALT marking from each sidecar's names, against the same index.
    char dir[4096];
    REQUIRE(bwa_test_mkdtemp(dir, sizeof(dir), "alt_gate") == 0);
    struct DirGuard {
        const char *d;
        ~DirGuard() { bwa_test_remove_flat_dir(d); }
    } dir_guard{dir};
    BnsGuard src_bns(restore_with_alt(dir, "source", in.prefix, in.source));
    BnsGuard mask_bns(restore_with_alt(dir, "masked", in.prefix, in.masked));
    REQUIRE(src_bns.b != nullptr);
    REQUIRE(mask_bns.b != nullptr);
    REQUIRE(src_bns.b->n_seqs == mask_bns.b->n_seqs);
    CHECK(count_bns_differences(src_bns.b, mask_bns.b) == 0);
    // The deployed layout: the index's own sidecar, read by bns_restore(prefix).
    if (file_exists(in.prefix + ".alt")) {
        BnsGuard own_bns(bns_restore(in.prefix.c_str()));
        REQUIRE(own_bns.b != nullptr);
        REQUIRE(own_bns.b->n_seqs == src_bns.b->n_seqs);
        CHECK(count_bns_differences(own_bns.b, src_bns.b) == 0);
    }
    counts.n_seqs = src_bns.b->n_seqs;
    std::vector<int> alt_bits(counts.n_seqs);
    for (int rid = 0; rid < counts.n_seqs; ++rid) {
        alt_bits[rid] = src_bns.b->anns[rid].is_alt;
        counts.alt_contigs += alt_bits[rid] != 0;
    }

    // Row-level oracle shape: same QNAME order and non-FLAG bytes, FLAG either
    // unchanged or exactly FLAG | 4 on a record that lacked 0x4, unique ALT names.
    const auto src_rows = read_sidecar(in.source), mask_rows = read_sidecar(in.masked);
    REQUIRE(src_rows.size() == mask_rows.size());
    std::unordered_map<std::string, int> rid_of;
    rid_of.reserve(counts.n_seqs);
    for (int rid = 0; rid < counts.n_seqs; ++rid) rid_of.emplace(src_bns.b->anns[rid].name, rid);
    std::vector<int> row_rid(src_rows.size());
    std::vector<char> row_masked(src_rows.size());
    std::unordered_set<std::string> seen;
    seen.reserve(src_rows.size());
    long bad_shape = 0, bad_name = 0;
    for (size_t i = 0; i < src_rows.size(); ++i) {
        const auto &s = src_rows[i], &m = mask_rows[i];
        const bool shape_ok = s.qname == m.qname && s.rest == m.rest && !(s.flag & 4) &&
                              (m.flag == s.flag || m.flag == (s.flag | 4));
        if (!shape_ok) {
            if (bad_shape == 0) MESSAGE("first bad masked row " << i << " " << s.qname);
            ++bad_shape;
        }
        const auto it = rid_of.find(s.qname);
        const bool name_ok = it != rid_of.end() && alt_bits[it->second] && seen.insert(s.qname).second;
        if (!name_ok) {
            if (bad_name == 0) MESSAGE("first unknown, non-ALT or repeated QNAME at row " << i << " " << s.qname);
            ++bad_name;
        }
        row_rid[i] = it == rid_of.end() ? -1 : it->second;
        row_masked[i] = m.flag != s.flag;
    }
    CHECK(bad_shape == 0);
    REQUIRE(bad_name == 0);
    counts.rows = (long) src_rows.size();

    // Usability: unfiltered, every record is a table entry, so with unique names
    // each per-record decision is exactly one table entry's presence.
    LiftGuard plain(alt_lift_load(in.source.c_str(), src_bns.b));
    REQUIRE(plain.l != nullptr);
    const alt_lift_skips_t &ps = plain.l->skip;
    CHECK(ps.no_ctg + ps.no_target + ps.alt_target + ps.unmapped + ps.not_primary + ps.dup + ps.bad == 0);
    REQUIRE(plain.l->n_lift == counts.rows);

    alt_lift_exclude_stats_t stats = {};
    char error[1024] = {};
    LiftGuard runtime(alt_lift_load_excluded(in.source.c_str(), src_bns.b, in.bed.c_str(), &stats, error,
                                             sizeof(error)));
    INFO("restricted load error: " << error);
    REQUIRE(runtime.l != nullptr);  // unknown or out-of-bounds BED rows stop here
    LiftGuard oracle(alt_lift_load(in.masked.c_str(), mask_bns.b));
    REQUIRE(oracle.l != nullptr);
    if (tamper) tamper(oracle.l);

    // Per-record decisions against the oracle's FLAG bit.
    long bad_decision = 0;
    for (size_t i = 0; i < src_rows.size(); ++i) {
        const int rid = row_rid[i];
        const bool runtime_excluded = runtime.l->a[rid].trid < 0;
        const bool oracle_excluded = oracle.l->a[rid].trid < 0;
        if (runtime_excluded != (bool) row_masked[i] || oracle_excluded != (bool) row_masked[i]) {
            if (bad_decision == 0) MESSAGE("first differing decision at row " << i << " " << src_rows[i].qname);
            ++bad_decision;
        }
        if (row_masked[i]) ++counts.excluded; else ++counts.retained;
    }
    CHECK(bad_decision == 0);
    CHECK(stats.excluded == counts.excluded);
    CHECK(runtime.l->n_lift == counts.retained);
    CHECK(oracle.l->skip.unmapped == counts.excluded);

    // The full table, every contig, field by field.
    check_lift_tables(runtime.l, oracle.l);

    // Loading never changes ALT marking.
    for (const bntseq_t *b : {(const bntseq_t *) src_bns.b, (const bntseq_t *) mask_bns.b}) {
        long changed = 0;
        for (int rid = 0; rid < counts.n_seqs; ++rid) changed += b->anns[rid].is_alt != alt_bits[rid];
        CHECK(changed == 0);
    }
    counts.input_intervals = stats.input_intervals;
    counts.merged_intervals = stats.merged_intervals;
    return counts;
}

// A generated index with four ALT tiles: a forward tile and a tile whose
// deletion holds the only BED base in its footprint are excluded; a reverse tile
// and a tile with an insertion whose footprint ends where a BED interval starts
// are retained. Names include the alt_/altc_ prefixes and an HLA-style name.
struct AltGateFixture {
    std::string chr1, chr2;
    bwa_tests::TempIndex idx;
    TmpFile source, masked, bed;
    static std::string random_seq(std::mt19937 &rng, size_t n) {
        std::string s(n, 'A');
        for (char &c : s) c = "ACGT"[rng() % 4];
        return s;
    }
    static std::string sidecar(bool masked, int wrong_flag_row) {
        const int flags[4] = {0, 16, 0, 0};
        const bool excluded[4] = {true, false, true, false};
        const char *rest[4] = {
            "\tchr1\t101\t60\t200M\t*\t0\t0\t*\t*\n",
            "\tchr1\t1001\t60\t200M\t*\t0\t0\t*\t*\n",
            "\tchr2\t501\t60\t100M10D100M\t*\t0\t0\t*\t*\n",
            "\tchr1\t2001\t60\t75M5I75M\t*\t0\t0\t*\t*\n",
        };
        const char *names[4] = {"alt_chr1_t1", "alt_chr1_t2", "altc_chr2_t3", "HLA-A*02:01:01:01"};
        std::string out;
        for (int i = 0; i < 4; ++i) {
            const bool mask = masked && (excluded[i] != (i == wrong_flag_row));
            out += std::string(names[i]) + "\t" + std::to_string(flags[i] | (mask ? 4 : 0)) + rest[i];
        }
        return out;
    }
    static std::vector<bwa_tests::Contig> contigs(const std::string &c1, const std::string &c2) {
        return {{"chr1", c1},
                {"chr2", c2},
                {"alt_chr1_t1", c1.substr(100, 200)},
                {"alt_chr1_t2", bwa_tests::revcomp(c1.substr(1000, 200))},
                {"altc_chr2_t3", c2.substr(500, 100) + c2.substr(610, 100)},
                {"HLA-A*02:01:01:01", c1.substr(2000, 75) + "ACGTA" + c1.substr(2075, 75)}};
    }
    explicit AltGateFixture(int wrong_flag_row = -1, std::mt19937 rng = std::mt19937(7754))
        : chr1(random_seq(rng, 3000)), chr2(random_seq(rng, 2000)), idx(contigs(chr1, chr2)),
          source(sidecar(false, -1).c_str()), masked(sidecar(true, wrong_flag_row).c_str()),
          bed("chr1\t250\t260\nchr2\t605\t606\nchr1\t2150\t2200\n") {}
    AltGateInputs inputs() const { return {idx.prefix(), source.path, masked.path, bed.path}; }
};

const char *const kAltGateEnv[4] = {"BWAMEM3_ALT_GATE_PREFIX", "BWAMEM3_ALT_GATE_SOURCE",
                                    "BWAMEM3_ALT_GATE_MASKED", "BWAMEM3_ALT_GATE_BED"};

bool alt_gate_any_env_set() {
    for (const char *name : kAltGateEnv) {
        const char *v = getenv(name);
        if (v && *v) return true;
    }
    return false;
}

// Counts failed CHECK/REQUIRE assertions (not WARNs) over the whole run, so a
// test case can tell whether any assertion failed since a given point. doctest
// reports every assertion to listeners, including non-fatal CHECK failures that
// do not stop the case.
std::atomic<long> g_alt_gate_failed_asserts{0};

struct AltGateFailureCounter : doctest::IReporter {
    explicit AltGateFailureCounter(const doctest::ContextOptions &) {}
    void report_query(const doctest::QueryData &) override {}
    void test_run_start() override {}
    void test_run_end(const doctest::TestRunStats &) override {}
    void test_case_start(const doctest::TestCaseData &) override {}
    void test_case_reenter(const doctest::TestCaseData &) override {}
    void test_case_end(const doctest::CurrentTestCaseStats &) override {}
    void test_case_exception(const doctest::TestCaseException &) override {}
    void subcase_start(const doctest::SubcaseSignature &) override {}
    void subcase_end() override {}
    void log_assert(const doctest::AssertData &ad) override {
        if (ad.m_failed && !(ad.m_at & doctest::assertType::is_warn)) ++g_alt_gate_failed_asserts;
    }
    void log_message(const doctest::MessageData &) override {}
    void test_case_skipped(const doctest::TestCaseData &) override {}
};

// Writes the ALT_GATE_FROZEN_COMPLETE line to `out` only if no assertion has
// failed since `failed_before` was read; returns whether it wrote the line.
bool alt_gate_emit_marker(FILE *out, long failed_before, const AltGateCounts &c) {
    if (g_alt_gate_failed_asserts.load() != failed_before) return false;
    fprintf(out,
            "ALT_GATE_FROZEN_COMPLETE n_seqs=%ld alt_contigs=%ld rows=%ld excluded=%ld retained=%ld "
            "bed_intervals=%ld merged_intervals=%ld\n",
            c.n_seqs, c.alt_contigs, c.rows, c.excluded, c.retained, c.input_intervals, c.merged_intervals);
    fflush(out);
    return true;
}

// Reads back what alt_gate_emit_marker wrote to a temporary file.
std::string alt_gate_marker_text(long failed_before, const AltGateCounts &c, bool *emitted) {
    FILE *f = tmpfile();
    REQUIRE(f != nullptr);
    *emitted = alt_gate_emit_marker(f, failed_before, c);
    rewind(f);
    std::string text;
    char buf[256];
    while (fgets(buf, sizeof(buf), f)) text += buf;
    fclose(f);
    return text;
}
} // namespace

DOCTEST_REGISTER_LISTENER("alt_gate_failure_counter", 1, AltGateFailureCounter);

TEST_CASE("exclude BED table gate on a generated index" * doctest::test_suite("unit/alt_lift")) {
    AltGateFixture fx;
    const long failed_before = g_alt_gate_failed_asserts.load();
    const AltGateCounts c = run_alt_table_gate(fx.inputs());
    CHECK(g_alt_gate_failed_asserts.load() == failed_before);
    bool emitted = false;
    const std::string marker = alt_gate_marker_text(failed_before, c, &emitted);
    CHECK(emitted);
    CHECK(marker.rfind("ALT_GATE_FROZEN_COMPLETE n_seqs=6 alt_contigs=4 rows=4 excluded=2 retained=2 ", 0) == 0);
    CHECK(c.n_seqs == 6);
    CHECK(c.alt_contigs == 4);
    CHECK(c.rows == 4);
    CHECK(c.excluded == 2);
    CHECK(c.retained == 2);
    CHECK(c.input_intervals == 3);
    CHECK(c.merged_intervals == 3);
}

// The gate must fail on a single wrong masked FLAG, in either direction. Each
// direction fails exactly five assertions: the row decision, the excluded and
// retained counts, the comparator's n_lift and its fatal per-entry trid check.
TEST_CASE("exclude BED table gate catches a masked FLAG wrongly set" * doctest::test_suite("unit/alt_lift") *
          doctest::expected_failures(5)) {
    AltGateFixture fx(/*wrong_flag_row=*/1);
    run_alt_table_gate(fx.inputs());
}

TEST_CASE("exclude BED table gate catches a masked FLAG wrongly clear" * doctest::test_suite("unit/alt_lift") *
          doctest::expected_failures(5)) {
    AltGateFixture fx(/*wrong_flag_row=*/2);
    run_alt_table_gate(fx.inputs());
}

// One differing selected target or offset is exactly one failed assertion.
TEST_CASE("exclude BED table gate catches one differing target" * doctest::test_suite("unit/alt_lift") *
          doctest::expected_failures(1)) {
    AltGateFixture fx;
    // alt_chr1_t2 (rid 3) is retained on chr1; were it not, more checks would fail.
    run_alt_table_gate(fx.inputs(), [](alt_lift_t *l) { l->a[3].trid = 1; });
}

TEST_CASE("exclude BED table gate catches one differing offset" * doctest::test_suite("unit/alt_lift") *
          doctest::expected_failures(1)) {
    AltGateFixture fx;
    // The HLA-style tile (rid 5, 75M5I75M) is retained: reference bases before op 2.
    run_alt_table_gate(fx.inputs(), [](alt_lift_t *l) {
        if (l->a[5].n_cigar == 3) l->a[5].roff[2] += 1;
    });
}

TEST_CASE("exclude BED table gate catches one differing position" * doctest::test_suite("unit/alt_lift") *
          doctest::expected_failures(1)) {
    AltGateFixture fx;
    run_alt_table_gate(fx.inputs(), [](alt_lift_t *l) { l->a[5].tpos += 1; });
}

// A non-fatal (CHECK) gate failure must withhold the completion marker. The
// offset tamper fails exactly one CHECK and lets the gate run to the end; any
// extra failure here (counter missed it, or the marker was written) makes the
// case fail its expected count.
TEST_CASE("exclude BED frozen marker is withheld after a non-fatal gate failure" *
          doctest::test_suite("unit/alt_lift") * doctest::expected_failures(1)) {
    AltGateFixture fx;
    const long failed_before = g_alt_gate_failed_asserts.load();
    const AltGateCounts c = run_alt_table_gate(fx.inputs(), [](alt_lift_t *l) {
        if (l->a[5].n_cigar == 3) l->a[5].roff[2] += 1;
    });
    const long failed_after = g_alt_gate_failed_asserts.load();
    bool emitted = true;
    const std::string marker = alt_gate_marker_text(failed_before, c, &emitted);
    if (failed_after - failed_before != 1) FAIL("failure counter saw " << failed_after - failed_before);
    if (emitted || !marker.empty()) FAIL("marker written after a failed gate: " << marker);
}

// Release gate on a real index, configured only by environment: unset, the case
// is skipped; set, all four inputs are required. Only when no assertion in the
// case has failed does it print one ALT_GATE_FROZEN_COMPLETE line with the counts
// it checked, so a release run can tell a completed gate from a skipped or
// failed one.
TEST_CASE("exclude BED frozen table gate" * doctest::test_suite("unit/alt_lift_frozen") *
          doctest::skip(!alt_gate_any_env_set())) {
    const long failed_before = g_alt_gate_failed_asserts.load();
    std::string value[4];
    for (int i = 0; i < 4; ++i) {
        const char *v = getenv(kAltGateEnv[i]);
        INFO("partial configuration: " << std::string(kAltGateEnv[i]) << " is unset");
        REQUIRE((v != nullptr && *v != '\0'));
        value[i] = v;
    }
    const AltGateCounts c = run_alt_table_gate({value[0], value[1], value[2], value[3]});
    REQUIRE(c.rows > 0);
    REQUIRE(c.excluded > 0);
    REQUIRE(c.retained > 0);
    REQUIRE(c.input_intervals > 0);
    alt_gate_emit_marker(stdout, failed_before, c);
}
