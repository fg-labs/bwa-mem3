// test/unit/test_mem_io_reg2sam.cpp
//
// mem_reg2sam: from a read's marked alignment regions to its SAM records. This
// is the step that decides which regions become records (primary, -a
// secondaries, supplementary, unmapped), sets the supplementary FLAG (0x800,
// or 0x100 under -M), caps a supplementary's MAPQ at the primary's (unless -q
// or -5), and builds the XA list and HN count of the hits a primary shadows
// (mem_gen_alt: -h, -z, -u). The regions are built by hand on a small indexed
// reference (bwa_tests::TempIndex) whose repeats are planted, so every record
// below is worked out from the construction:
//
//   c1  400 bp random.
//   c2  400 bp random, with c1[100, 160) copied to c2[200, 260) and a
//       substitution at copy offset 30 (the base becomes the next in ACGT).
//   c3  300 bp random, with the reverse complement of c1[100, 160) at
//       c3[50, 110).
//
// The regions are given, not found, so only the bases inside each region's
// windows matter. Coordinates in comments are 0-based half-open unless written
// as SAM POS. MAPQ comes from mem_approx_mapq_se (bwa's formula, default -A 1
// -B 4, mapQ_coef_len 50): for an exact hit shorter than 50 bp it is
// 6.02 * (score - s) rounded, capped at 60, where s is the larger of csub and
// sub (19, the minimum seed length, when sub is 0), and 0 when s reaches the
// hit's own score; a hit of l >= 50 bp is scaled by (3 / ln l)^2, where 3 is
// bwa's integer-typed coefficient (int)ln 50.

#include "doctest/doctest.h"
#include "mem_fixture.h"

#include "bwamem.h"

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <random>
#include <string>
#include <vector>

namespace {

std::string random_bases(std::mt19937 &rng, int n) {
    std::string s((size_t)n, 'A');
    for (char &c : s) c = "ACGT"[rng() % 4];
    return s;
}

char next_base(char c) { return c == 'A' ? 'C' : c == 'C' ? 'G' : c == 'G' ? 'T' : 'A'; }

struct Fixture {
    std::vector<bwa_tests::Contig> contigs;
    std::string read1;  // c1[100, 160)
    std::string read2;  // c1[300, 345) + c3[200, 235): a split read
    std::string read3;  // c3[200, 235) + c1[300, 345): the same, other way round
    std::string read4;  // c3[200, 230) + c1[100, 160): a split read whose longer part has copies
    std::unique_ptr<bwa_tests::TempIndex> idx;
    int64_t l_pac = 0;

    Fixture() {
        std::mt19937 rng(20261024u);
        std::string c1 = random_bases(rng, 400), c2 = random_bases(rng, 400),
                    c3 = random_bases(rng, 300);
        std::string copy = c1.substr(100, 60);
        copy[30] = next_base(copy[30]);
        c2.replace(200, 60, copy);
        c3.replace(50, 60, bwa_tests::revcomp(c1.substr(100, 60)));
        contigs = {{"c1", c1}, {"c2", c2}, {"c3", c3}};
        read1 = c1.substr(100, 60);
        read2 = c1.substr(300, 45) + c3.substr(200, 35);
        read3 = c3.substr(200, 35) + c1.substr(300, 45);
        read4 = c3.substr(200, 30) + c1.substr(100, 60);
        idx.reset(new bwa_tests::TempIndex(contigs));
        l_pac = idx->bns()->l_pac;
    }
};

const Fixture &fixture() {
    static Fixture fx;
    return fx;
}

// One read handed to mem_reg2sam: nt4 codes, a quality string with a distinct
// character per position (so a trimming or reversal error moves a visible
// character), and the name.
struct Read {
    std::string name, bases, qual;
    std::vector<char> codes;
    bseq1_t s{};

    Read(const std::string &n, const std::string &b) : name(n), bases(b) {
        for (size_t i = 0; i < b.size(); ++i) qual += (char)('!' + (int)(i % 90));
        std::vector<uint8_t> v = bwa_tests::nt4(b);
        codes.assign(v.begin(), v.end());
        s.l_seq = (int)codes.size();
        s.name = const_cast<char *>(name.c_str());
        s.seq = codes.data();
        s.qual = const_cast<char *>(qual.c_str());
        s.meth_base_ot = -1;
    }
};

// A region of read [qb, qe) on contig `rid` at contig coordinate pos0, with
// the given strand and scores; marked primary (secondary -1) by default.
mem_alnreg_t region(const mem_opt_t *opt, int rid, int64_t pos0, bool rev, int qb, int qe,
                    int score, int sub) {
    const Fixture &fx = fixture();
    const int64_t off = fx.idx->bns()->anns[rid].offset, len = qe - qb;
    mem_alnreg_t r;
    memset(&r, 0, sizeof(r));
    if (!rev) {
        r.rb = off + pos0;
    } else {  // the reverse strand is [l_pac, 2 * l_pac), mirrored
        r.rb = 2 * fx.l_pac - (off + pos0 + len);
    }
    r.re = r.rb + len;
    r.qb = qb;
    r.qe = qe;
    r.rid = rid;
    r.score = r.truesc = score;
    r.sub = sub;
    r.csub = 0;
    r.secondary = r.secondary_all = -1;
    r.w = opt->w;
    r.seedcov = r.seedlen0 = (int)len;
    r.chain_n_hits = 1;
    r.meth_hypothesis = r.meth_strand_hyp = -1;
    return r;
}

std::string emit(const mem_opt_t *opt, Read &r, std::vector<mem_alnreg_t> regs,
                 int extra_flag = 0) {
    const Fixture &fx = fixture();
    mem_alnreg_v a;
    a.n = a.m = regs.size();
    a.a = regs.data();
    a.capped_w = 0;
    mem_reg2sam(opt, fx.idx->bns(), fx.idx->pac(), &r.s, &a, extra_flag, NULL);
    REQUIRE(r.s.sam != NULL);
    std::string out(r.s.sam);
    free(r.s.sam);
    r.s.sam = NULL;
    return out;
}

struct OptGuard {
    mem_opt_t *opt;
    OptGuard() : opt(mem_opt_init()) {}
    ~OptGuard() { free(opt); }
};

// read1's regions: the c1 hit (primary), the c2 copy (one mismatch: score
// 60 - 1 - 4 = 55) and the c3 reverse-complement copy (score 60), both
// shadowed by the primary (secondary = secondary_all = 0).
std::vector<mem_alnreg_t> read1_regions(const mem_opt_t *opt) {
    std::vector<mem_alnreg_t> regs{region(opt, 0, 100, false, 0, 60, 60, 60),
                                   region(opt, 1, 200, false, 0, 60, 55, 0),
                                   region(opt, 2, 50, true, 0, 60, 60, 0)};
    for (size_t k = 1; k < regs.size(); ++k) regs[k].secondary = regs[k].secondary_all = 0;
    return regs;
}

// The MD of the c2 copy: 30 matches, the reference base there, 29 matches.
std::string c2_md() {
    return "30" + std::string(1, fixture().contigs[1].seq[230]) + "29";
}

std::string reverse(std::string s) { return std::string(s.rbegin(), s.rend()); }

}  // namespace

TEST_CASE("mem_io reg2sam: shadowed hits become the primary's XA list and HN count"
          * doctest::test_suite("unit/mem_io")) {
    OptGuard g;
    Read r("r1", fixture().read1);
    // Primary: c1 POS 101, 60M. sub 60 >= score 60, so MAPQ 0; XS 60. XA lists
    // the shadowed hits in region order as contig,strand+POS,CIGAR,NM: the c2
    // copy at POS 201 with one mismatch, the c3 copy on the reverse strand at
    // c3[50, 110) -> POS 51.
    const std::string head = "r1\t0\tc1\t101\t0\t60M\t*\t0\t0\t" + r.bases + "\t" + r.qual
                             + "\tNM:i:0\tMD:Z:60\tAS:i:60\tXS:i:60";

    SUBCASE("default: one record") {
        CHECK(emit(g.opt, r, read1_regions(g.opt))
              == head + "\tXA:Z:c2,+201,60M,1;c3,-51,60M,0;\tHN:i:2\n");
    }
    SUBCASE("-u: each XA entry also carries its score and MAPQ (0 for a secondary)") {
        g.opt->flag |= MEM_F_XB;
        CHECK(emit(g.opt, r, read1_regions(g.opt))
              == head + "\tXA:Z:c2,+201,60M,1,55,0;c3,-51,60M,0,60,0;\tHN:i:2\n");
    }
    SUBCASE("-h: more hits than the cap drop XA, but HN still counts them") {
        g.opt->max_XA_hits = 1;
        CHECK(emit(g.opt, r, read1_regions(g.opt)) == head + "\tHN:i:2\n");
        g.opt->max_XA_hits = 2;
        CHECK(emit(g.opt, r, read1_regions(g.opt))
              == head + "\tXA:Z:c2,+201,60M,1;c3,-51,60M,0;\tHN:i:2\n");
    }
    SUBCASE("a clipped shadowed hit keeps its clip in its XA entry") {
        // The c2 copy aligned from read base 5 only: c2[205, 260), 55 bases
        // with the one mismatch, score 54 - 4 = 50 (>= 0.8 * 60, so listed).
        // Its entry is c2,+206 with the 5 unaligned read bases as 5S: 5S55M.
        std::vector<mem_alnreg_t> regs = read1_regions(g.opt);
        regs[1] = region(g.opt, 1, 205, false, 5, 60, 50, 0);
        regs[1].secondary = regs[1].secondary_all = 0;
        CHECK(emit(g.opt, r, regs) == head + "\tXA:Z:c2,+206,5S55M,1;c3,-51,60M,0;\tHN:i:2\n");
    }
    SUBCASE("-z: a hit below the fraction of the primary's score is not counted") {
        // 55 < 0.95 * 60 = 57: the c2 copy drops out of XA and HN.
        g.opt->XA_drop_ratio = 0.95f;
        CHECK(emit(g.opt, r, read1_regions(g.opt)) == head + "\tXA:Z:c3,-51,60M,0;\tHN:i:1\n");
    }
    SUBCASE("-a: shadowed hits are 0x100 records instead, with no XA or HN") {
        g.opt->flag |= MEM_F_ALL;
        // A secondary prints SEQ/QUAL as *, MAPQ 0 and no XS (sub is reset).
        CHECK(emit(g.opt, r, read1_regions(g.opt))
              == head + "\n"
                 + "r1\t256\tc2\t201\t0\t60M\t*\t0\t0\t*\t*\tNM:i:1\tMD:Z:" + c2_md()
                 + "\tAS:i:55\n"
                 + "r1\t272\tc3\t51\t0\t60M\t*\t0\t0\t*\t*\tNM:i:0\tMD:Z:60\tAS:i:60\n");
        // -a still drops a secondary below drop_ratio (0.5) of its parent.
        std::vector<mem_alnreg_t> regs = read1_regions(g.opt);
        g.opt->drop_ratio = 0.95f;  // 55 < 0.95 * 60
        CHECK(emit(g.opt, r, regs)
              == head + "\n" + "r1\t272\tc3\t51\t0\t60M\t*\t0\t0\t*\t*\tNM:i:0\tMD:Z:60\tAS:i:60\n");
    }
    SUBCASE("the caller's extra FLAG bits reach every record") {
        g.opt->flag |= MEM_F_ALL;
        const std::string out = emit(g.opt, r, read1_regions(g.opt), 0x200);
        CHECK(out.find("r1\t512\tc1\t101\t") == 0);
        CHECK(out.find("\nr1\t768\tc2\t201\t") != std::string::npos);
        CHECK(out.find("\nr1\t784\tc3\t51\t") != std::string::npos);
    }
}

TEST_CASE("mem_io reg2sam: a split read gives a primary and a supplementary record"
          * doctest::test_suite("unit/mem_io")) {
    OptGuard g;
    Read r("r2", fixture().read2);
    // Read bases 0-44 on c1[300, 345) (score 45, sub 40: MAPQ 6.02 * 5 = 30.1
    // -> 30) and bases 45-79 on c3[200, 235) (score 35, sub 0 -> 19: MAPQ
    // 6.02 * 16 = 96.3 -> 60, capped at the primary's 30).
    std::vector<mem_alnreg_t> regs{region(g.opt, 0, 300, false, 0, 45, 45, 40),
                                   region(g.opt, 2, 200, false, 45, 80, 35, 0)};
    const std::string seq_tail = r.bases.substr(45), qual_tail = r.qual.substr(45);
    const std::string primary = "r2\t0\tc1\t301\t30\t45M35S\t*\t0\t0\t" + r.bases + "\t" + r.qual
                                + "\tNM:i:0\tMD:Z:45\tAS:i:45\tXS:i:40\tSA:Z:c3,201,+,45S35M,";

    SUBCASE("default: 0x800, hard clips, MAPQ capped at the primary's") {
        CHECK(emit(g.opt, r, regs)
              == primary + "30,0;\tHN:i:0\n"
                 + "r2\t2048\tc3\t201\t30\t45H35M\t*\t0\t0\t" + seq_tail + "\t" + qual_tail
                 + "\tNM:i:0\tMD:Z:35\tAS:i:35\tXS:i:0\tSA:Z:c1,301,+,45M35S,30,0;\tHN:i:0\n");
    }
    SUBCASE("-q: the supplementary keeps its own MAPQ") {
        g.opt->flag |= MEM_F_KEEP_SUPP_MAPQ;
        CHECK(emit(g.opt, r, regs)
              == primary + "60,0;\tHN:i:0\n"
                 + "r2\t2048\tc3\t201\t60\t45H35M\t*\t0\t0\t" + seq_tail + "\t" + qual_tail
                 + "\tNM:i:0\tMD:Z:35\tAS:i:35\tXS:i:0\tSA:Z:c1,301,+,45M35S,30,0;\tHN:i:0\n");
    }
    SUBCASE("-M: the shorter hit is FLAG 256 instead") {
        g.opt->flag |= MEM_F_NO_MULTI;
        CHECK(emit(g.opt, r, regs)
              == primary + "30,0;\tHN:i:0\n"
                 + "r2\t256\tc3\t201\t30\t45H35M\t*\t0\t0\t" + seq_tail + "\t" + qual_tail
                 + "\tNM:i:0\tMD:Z:35\tAS:i:35\tXS:i:0\tSA:Z:c1,301,+,45M35S,30,0;\tHN:i:0\n");
    }
    SUBCASE("-Y: soft clips and the whole read on the supplementary") {
        g.opt->flag |= MEM_F_SOFTCLIP;
        CHECK(emit(g.opt, r, regs)
              == primary + "30,0;\tHN:i:0\n"
                 + "r2\t2048\tc3\t201\t30\t45S35M\t*\t0\t0\t" + r.bases + "\t" + r.qual
                 + "\tNM:i:0\tMD:Z:35\tAS:i:35\tXS:i:0\tSA:Z:c1,301,+,45M35S,30,0;\tHN:i:0\n");
    }
    SUBCASE("-T above the shorter hit's score: no supplementary, no SA") {
        g.opt->T = 40;
        CHECK(emit(g.opt, r, regs)
              == "r2\t0\tc1\t301\t30\t45M35S\t*\t0\t0\t" + r.bases + "\t" + r.qual
                 + "\tNM:i:0\tMD:Z:45\tAS:i:45\tXS:i:40\tHN:i:0\n");
    }
    SUBCASE("-T above every hit: one unmapped record") {
        g.opt->T = 50;
        CHECK(emit(g.opt, r, regs)
              == "r2\t4\t*\t0\t0\t*\t*\t0\t0\t" + r.bases + "\t" + r.qual + "\tAS:i:0\tXS:i:0\n");
    }
    SUBCASE("no regions at all: one unmapped record") {
        CHECK(emit(g.opt, r, std::vector<mem_alnreg_t>())
              == "r2\t4\t*\t0\t0\t*\t*\t0\t0\t" + r.bases + "\t" + r.qual + "\tAS:i:0\tXS:i:0\n");
    }
    SUBCASE("a reverse-strand supplementary") {
        // Read bases 45-79 as the reverse complement of c3[200, 235): the hit is
        // on the reverse strand, its CIGAR in reference orientation puts the 45
        // clipped read bases last (35M45H), SEQ is c3[200, 235) itself and QUAL
        // the reversed quality of read bases 45-79.
        const Fixture &fx = fixture();
        Read rr("r2", fx.read2.substr(0, 45) + bwa_tests::revcomp(fx.contigs[2].seq.substr(200, 35)));
        std::vector<mem_alnreg_t> rv{region(g.opt, 0, 300, false, 0, 45, 45, 40),
                                     region(g.opt, 2, 200, true, 45, 80, 35, 0)};
        CHECK(emit(g.opt, rr, rv)
              == "r2\t0\tc1\t301\t30\t45M35S\t*\t0\t0\t" + rr.bases + "\t" + rr.qual
                 + "\tNM:i:0\tMD:Z:45\tAS:i:45\tXS:i:40\tSA:Z:c3,201,-,35M45S,30,0;\tHN:i:0\n"
                 + "r2\t2064\tc3\t201\t30\t35M45H\t*\t0\t0\t" + fx.contigs[2].seq.substr(200, 35)
                 + "\t" + reverse(rr.qual.substr(45))
                 + "\tNM:i:0\tMD:Z:35\tAS:i:35\tXS:i:0\tSA:Z:c1,301,+,45M35S,30,0;\tHN:i:0\n");
    }
}

TEST_CASE("mem_io reg2sam: -5 makes the hit nearest the read's 5' end primary"
          * doctest::test_suite("unit/mem_io")) {
    // read3 is c3[200, 235) + c1[300, 345). Marking puts the higher-scoring c1
    // hit (read bases 35-79) first; mem_reorder_primary5 moves the hit with the
    // smallest query start (the c3 hit, read bases 0-34) to the front, and -5
    // also keeps supplementary MAPQ (it implies -q).
    OptGuard g;
    g.opt->flag |= MEM_F_PRIMARY5 | MEM_F_KEEP_SUPP_MAPQ;
    Read r("r3", fixture().read3);
    std::vector<mem_alnreg_t> regs{region(g.opt, 0, 300, false, 35, 80, 45, 40),
                                   region(g.opt, 2, 200, false, 0, 35, 35, 0)};
    mem_alnreg_v a;
    a.n = a.m = regs.size();
    a.a = regs.data();
    a.capped_w = 0;
    mem_reorder_primary5(g.opt->T, &a);
    REQUIRE(regs[0].rid == 2);
    REQUIRE(regs[1].rid == 0);
    // Primary: c3 POS 201, 35M45S, MAPQ 60 (sub 19: 6.02 * 16 -> 60). The c1
    // hit keeps its own MAPQ 30 under -q.
    CHECK(emit(g.opt, r, regs)
          == "r3\t0\tc3\t201\t60\t35M45S\t*\t0\t0\t" + r.bases + "\t" + r.qual
             + "\tNM:i:0\tMD:Z:35\tAS:i:35\tXS:i:0\tSA:Z:c1,301,+,35S45M,30,0;\tHN:i:0\n"
             + "r3\t2048\tc1\t301\t30\t35H45M\t*\t0\t0\t" + r.bases.substr(35) + "\t"
             + r.qual.substr(35)
             + "\tNM:i:0\tMD:Z:45\tAS:i:45\tXS:i:40\tSA:Z:c3,201,+,35M45S,60,0;\tHN:i:0\n");

    SUBCASE("the XA list and HN stay with the hit they shadow") {
        // read4 is c3[200, 230) + c1[100, 160). Regions in marking order: the
        // c1 hit (read bases 30-89, score 60, sub 55), the c3 hit (bases 0-29,
        // score 30) and the c2 copy of the c1 part (score 55, one mismatch),
        // shadowed by the c1 hit (secondary = secondary_all = 0). -5 moves the
        // c3 hit to the front; the c2 copy must still be listed under the c1
        // hit, now the supplementary record: XA:Z:c2,+201,30S60M,1; (read bases
        // 0-29 unaligned there) and HN 1, with HN 0 and no XA on the c3 primary.
        // MAPQ: c3 hit 6.02 * (30 - 19) = 66.2 -> 60. c1 hit: 60 bases >= 50,
        // so the length factor applies, with bwa's integer coefficient
        // (int)ln 50 = 3: 6.02 * (60 - 55) * (3 / ln 60)^2 = 30.1 * 0.537
        // = 16.2 -> 16, kept under -5 (which implies -q).
        const Fixture &fx = fixture();
        Read r4("r4", fx.read4);
        std::vector<mem_alnreg_t> regs{region(g.opt, 0, 100, false, 30, 90, 60, 55),
                                       region(g.opt, 2, 200, false, 0, 30, 30, 0),
                                       region(g.opt, 1, 200, false, 30, 90, 55, 0)};
        regs[2].secondary = regs[2].secondary_all = 0;
        mem_alnreg_v v;
        v.n = v.m = regs.size();
        v.a = regs.data();
        v.capped_w = 0;
        mem_reorder_primary5(g.opt->T, &v);
        REQUIRE(regs[0].rid == 2);
        CHECK(emit(g.opt, r4, regs)
              == "r4\t0\tc3\t201\t60\t30M60S\t*\t0\t0\t" + r4.bases + "\t" + r4.qual
                 + "\tNM:i:0\tMD:Z:30\tAS:i:30\tXS:i:0\tSA:Z:c1,101,+,30S60M,16,0;\tHN:i:0\n"
                 + "r4\t2048\tc1\t101\t16\t30H60M\t*\t0\t0\t" + r4.bases.substr(30) + "\t"
                 + r4.qual.substr(30)
                 + "\tNM:i:0\tMD:Z:60\tAS:i:60\tXS:i:55\tSA:Z:c3,201,+,30M60S,60,0;"
                   "\tXA:Z:c2,+201,30S60M,1;\tHN:i:1\n");
    }
    SUBCASE("a single primary is left alone") {
        std::vector<mem_alnreg_t> one{region(g.opt, 0, 300, false, 35, 80, 45, 40)};
        mem_alnreg_v b;
        b.n = b.m = 1;
        b.a = one.data();
        b.capped_w = 0;
        mem_reorder_primary5(g.opt->T, &b);
        CHECK(one[0].rid == 0);
        CHECK(one[0].qb == 35);
    }
}
