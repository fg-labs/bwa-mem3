// Unit tests for `aln --meth`: orientation, chemistry and tags, reward ranking,
// bounded candidate discovery, conservative caps, duplicate/thread invariance,
// and stock alignment against a --meth-built index.

#include "doctest/doctest.h"

#include "aln_test_ref.h"
#include "../framework/stderr_capture.h"
#include "../../src/aln/aln_main.h"
#include "../../src/aln/aln_meth.h"
#include "../../src/aln/aln_meth_search.h"
#include "../../src/aln/bwaseqio.h"
#include "../../src/aln/bwase.h"
#include "../../src/aln/aln_out.h"
#include "../../src/FMI_search.h"
#include "../../src/bwa.h"
#include "../../src/main.h"

#include <getopt.h>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <unistd.h>

namespace {

std::string conv(std::string s, char from, char to)
{
    for (char &c : s)
        if (c == from) c = to;
    return s;
}

// The reference: c1 and c2 random; c3 holds the -R counterexample's loci (see below).
struct MethRef {
    bwa_tests::TempDir tmp{"aln_meth_unit"}; // removed (with its files) at process exit
    std::string dir, fa, c1, c2, read_r, c3, affine_read;
    int x_pos = 0; // 1-based position of the counterexample's true locus in c3

    MethRef()
    {
        dir = tmp.path();
        fa = dir + "/ref.fa";
        uint32_t st = 1234567;
        c1 = aln_test::lcg_dna(st, 3000);
        c2 = aln_test::lcg_dna(st, 2000);
        c2 += c1.substr(2500, 50); // a copy of a c1 segment: a read from it has an XA hit on another contig
        // The -R counterexample (taps, neutral, read 1 forward, so OT: ref T under read C is a
        // strict mismatch the converted search frees). R is the read. X = R with one A/G
        // substitution: collapsed 1, strict 1. Y (31 copies) = R with two of R's C made T:
        // collapsed 0, strict 2. Z = X with two more of R's C made T: collapsed 1, strict 3.
        // Stock's -R on the collapsed counts would stop on Y's 31 equally best hits; on the
        // rescored counts the read is unique at X.
        std::string r;
        for (;;) {
            r = aln_test::lcg_dna(st, 40);
            int nc = 0;
            for (char c : r) nc += c == 'C';
            if (nc >= 6 && (r[20] == 'A' || r[20] == 'G')) break;
        }
        read_r = r;
        std::vector<int> cpos;
        for (int i = 0; i < 40; ++i)
            if (r[i] == 'C') cpos.push_back(i);
        std::string x = r;
        x[20] = r[20] == 'A' ? 'G' : 'A';
        std::string y = r;
        y[cpos[1]] = 'T'; y[cpos[3]] = 'T';
        std::string z = x;
        z[cpos[2]] = 'T'; z[cpos[4]] = 'T';
        c3 = aln_test::lcg_dna(st, 100);
        x_pos = (int)c3.size() + 1;
        c3 += x + aln_test::lcg_dna(st, 100) + z + aln_test::lcg_dna(st, 100);
        for (int k = 0; k < 31; ++k) c3 += y + aln_test::lcg_dna(st, 60);
        // First locus: two reverse-conversion mismatches, collapsed cost zero, strict cost 12.
        // Second locus: a 5-base read insertion, five user differences but affine cost 5.
        affine_read = aln_test::lcg_dna(st, 60);
        std::string expensive = affine_read;
        int changed = 0;
        for (int i = 10; i < 50 && changed < 2; ++i)
            if (expensive[i] == 'C') { expensive[i] = 'T'; ++changed; }
        REQUIRE(changed == 2);
        std::string cheap = affine_read;
        cheap.erase(28, 5);
        REQUIRE(aln_test::write_test_fasta(fa, {{"c1", c1}, {"c2", c2}, {"c3", c3},
                                              {"affine_expensive", expensive}, {"affine_cheap", cheap}}));
        std::vector<std::string> a = {"index", "--meth", fa};
        std::vector<char *> argv;
        for (std::string &s : a) argv.push_back(&s[0]);
        argv.push_back(nullptr);
        optind = 1;
#if defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || defined(__NetBSD__)
        optreset = 1;
#endif
        REQUIRE(bwa_index((int)a.size(), argv.data()) == 0);
    }
};

MethRef &ref()
{
    static MethRef r;
    return r;
}

std::string write_fq(const std::string &name, const std::vector<std::string> &reads)
{
    std::string path = ref().dir + "/" + name;
    FILE *fp = fopen(path.c_str(), "w");
    REQUIRE(fp != nullptr);
    for (size_t i = 0; i < reads.size(); ++i)
        fprintf(fp, "@r%zu\n%s\n+\n%s\n", i, reads[i].c_str(), std::string(reads[i].size(), 'I').c_str());
    REQUIRE(fclose(fp) == 0);
    return path;
}

struct Rec {
    std::string name, rname, cigar;
    int flag = 0, pos = 0, mapq = 0, pnext = 0, tlen = 0;
    std::string rnext;
    std::map<std::string, std::string> tag; // "XG:Z" -> "CT"
    std::vector<std::string> raw_tags;
};

// The alignment records of a SAM stream.
std::vector<Rec> parse_records(std::istream &in)
{
    std::vector<Rec> v;
    for (std::string l; std::getline(in, l);) {
        if (l.empty() || l[0] == '@') continue;
        std::vector<std::string> f;
        std::stringstream ss(l);
        for (std::string x; std::getline(ss, x, '\t');) f.push_back(x);
        Rec r;
        REQUIRE(f.size() >= 11);
        r.name = f[0]; r.flag = atoi(f[1].c_str()); r.rname = f[2]; r.pos = atoi(f[3].c_str());
        r.mapq = atoi(f[4].c_str()); r.cigar = f[5];
        r.rnext=f[6];r.pnext=atoi(f[7].c_str());r.tlen=atoi(f[8].c_str());
        for (size_t i = 11; i < f.size(); ++i) {
            r.raw_tags.push_back(f[i]);
            REQUIRE(f[i].size() >= 5);
            REQUIRE(f[i][2] == ':');
            REQUIRE(f[i][4] == ':');
            r.tag[f[i].substr(0, 4)] = f[i].substr(5);
        }
        v.push_back(r);
    }
    return v;
}

std::vector<Rec> run_aln(std::vector<std::string> args, const std::vector<std::string> &inputs,
                         const std::string &reference = "")
{
    const std::string out = ref().dir + "/out.sam";
    std::vector<std::string> a = {"aln", "-f", out};
    a.insert(a.end(), args.begin(), args.end());
    a.push_back(reference.empty()?ref().fa:reference);
    a.insert(a.end(), inputs.begin(), inputs.end());
    std::vector<char *> argv;
    for (std::string &s : a) argv.push_back(&s[0]);
    argv.push_back(nullptr);
    REQUIRE(main_aln((int)a.size(), argv.data()) == 0);
    std::ifstream in(out);
    return parse_records(in);
}

std::string slurp_out()
{
    std::ifstream in(ref().dir + "/out.sam");
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

bool has_tag(const Rec &r, const char *tag_type) // e.g. "XM:i"
{
    for (const std::string &t : r.raw_tags)
        if (t.compare(0, 4, tag_type) == 0) return true;
    return false;
}

const std::string &tag_value(const Rec &r, const char *key) {
    REQUIRE(r.tag.count(key) == 1);
    return r.tag.at(key);
}

std::string expected_md(const std::string &genome, int start, const std::string &aligned) {
    std::string md;
    int run = 0;
    for (size_t i = 0; i < aligned.size(); ++i) {
        const char base = genome[start + i];
        if (base == aligned[i]) ++run;
        else { md += std::to_string(run) + base; run = 0; }
    }
    return md + std::to_string(run);
}

std::string expected_xm(const std::string &genome, int start, const std::string &aligned, bool top, bool taps) {
    std::string xm(aligned.size(), '.');
    for (size_t i = 0; i < aligned.size(); ++i) {
        const int pos = start + (int)i;
        if (genome[pos] != (top ? 'C' : 'G')) continue;
        const char retained = top ? 'C' : 'G', converted = top ? 'T' : 'A';
        if (aligned[i] != retained && aligned[i] != converted) continue;
        const int step = top ? 1 : -1;
        const char partner = top ? 'G' : 'C';
        auto base_at=[&](int at) { return at<0 || at>=(int)genome.size()?'N':genome[at]; };
        const char first=base_at(pos+step),second=base_at(pos+2*step);
        const char context = first=='N'?'u':first==partner?'z':second=='N'?'u':second==partner?'x':'h';
        const bool methylated = taps ? aligned[i] == converted : aligned[i] == retained;
        xm[i] = methylated ? context - 'a' + 'A' : context;
    }
    return xm;
}

} // namespace

TEST_SUITE_BEGIN("unit/aln_meth");

TEST_CASE("aln_meth: read 1 on both genome strands: position, strand, XR/XG, NM and MD")
{
    const std::string x = ref().c1.substr(500, 50);
    const std::string fwd = conv(x, 'C', 'T');                       // OT, forward
    const std::string rev = conv(aln_test::revcomp(x), 'C', 'T');    // OB, reverse
    auto v = run_aln({"--meth"}, {write_fq("se.fq", {fwd, rev})});
    REQUIRE(v.size() == 2);
    CHECK(v[0].flag == 0);
    CHECK(v[0].rname == "c1");
    CHECK(v[0].pos == 501);
    CHECK(tag_value(v[0], "XR:Z") == "CT");
    CHECK(tag_value(v[0], "XG:Z") == "CT");
    CHECK(tag_value(v[0], "NM:i") == "0");
    CHECK(tag_value(v[0], "MD:Z") == expected_md(ref().c1, 500, fwd));
    CHECK(tag_value(v[0], "XT:A") == "U");
    CHECK(v[1].flag == 16);
    CHECK(v[1].pos == 501);
    CHECK(tag_value(v[1], "XR:Z") == "CT");
    CHECK(tag_value(v[1], "XG:Z") == "GA");
    CHECK(tag_value(v[1], "NM:i") == "0");
    const std::string aligned_rev = aln_test::revcomp(rev);
    CHECK(tag_value(v[1], "MD:Z") == expected_md(ref().c1, 500, aligned_rev));
    CHECK(tag_value(v[0], "XM:Z") == expected_xm(ref().c1, 500, fwd, true, false));
    CHECK(tag_value(v[1], "XM:Z") == expected_xm(ref().c1, 500, aligned_rev, false, false));
    const auto taps = run_aln({"--meth=taps"}, {write_fq("se-taps.fq", {fwd, rev})});
    REQUIRE(taps.size() == 2);
    for (size_t i = 0; i < taps.size(); ++i) {
        REQUIRE(taps[i].pos == 501);
        CHECK(taps[i].flag == (i ? 16 : 0));
        const std::string aligned = i ? aligned_rev : fwd;
        CHECK(tag_value(taps[i], "MD:Z") == expected_md(ref().c1, 500, aligned));
        CHECK(tag_value(taps[i], "XM:Z") == expected_xm(ref().c1, 500, aligned, i == 0, true));
    }
    for (const Rec &r : v) {
        CHECK_FALSE(has_tag(r, "XM:i"));
        CHECK_FALSE(has_tag(r, "XO:i"));
        CHECK_FALSE(has_tag(r, "XG:i"));
        CHECK(tag_value(r, "XM:Z").size() == 50);
    }
}

TEST_CASE("aln_meth: read 2 converts G->A and pairs on the fragment's genome strand")
{
    const std::string f = ref().c1.substr(1000, 200);
    // OT fragment: read 1 from the top strand, read 2 from its reverse complement
    const std::string r1_ot = conv(f.substr(0, 50), 'C', 'T');
    const std::string r2_ot = conv(aln_test::revcomp(f).substr(0, 50), 'G', 'A');
    // OB fragment: read 1 from the bottom strand
    const std::string r1_ob = conv(aln_test::revcomp(f).substr(0, 50), 'C', 'T');
    const std::string r2_ob = conv(f.substr(0, 50), 'G', 'A');
    auto v = run_aln({"--meth"}, {write_fq("p1.fq", {r1_ot, r1_ob}), write_fq("p2.fq", {r2_ot, r2_ob})});
    REQUIRE(v.size() == 4);
    // pair 0: OT
    CHECK((v[0].flag & 0x10) == 0);
    CHECK((v[1].flag & 0x10) != 0);
    CHECK((v[0].flag & 0x2) != 0);
    CHECK(v[0].tag["XG:Z"] == "CT");
    CHECK(v[1].tag["XG:Z"] == "CT");
    CHECK(v[0].tag["XR:Z"] == "CT");
    CHECK(v[1].tag["XR:Z"] == "GA");
    CHECK(v[0].pos == 1001);
    CHECK(v[1].pos == 1151);
    CHECK(v[0].flag == 99); CHECK(v[1].flag == 147);
    CHECK(v[0].mapq == 60); CHECK(v[1].mapq == 60);
    CHECK(v[0].rnext=="=");CHECK(v[1].rnext=="=");
    CHECK(v[0].pnext==1151);CHECK(v[1].pnext==1001);
    CHECK(v[0].tlen==200);CHECK(v[1].tlen==-200);
    CHECK(v[0].cigar == "50M"); CHECK(v[1].cigar == "50M");
    CHECK(tag_value(v[0],"MD:Z") == expected_md(ref().c1,1000,r1_ot));
    CHECK(tag_value(v[1],"MD:Z") == expected_md(ref().c1,1150,aln_test::revcomp(r2_ot)));
    // pair 1: OB
    CHECK((v[2].flag & 0x10) != 0);
    CHECK((v[3].flag & 0x10) == 0);
    CHECK(v[2].tag["XG:Z"] == "GA");
    CHECK(v[3].tag["XG:Z"] == "GA");
    CHECK(v[3].tag["XR:Z"] == "GA");
    CHECK(v[2].tag["NM:i"] == "0");
    CHECK(v[3].tag["NM:i"] == "0");
}

TEST_CASE("aln_meth: scoring modes: a read C over a reference T counts unless collapsed")
{
    const std::string x = ref().c2.substr(300, 40);
    std::string r = x;
    int t_at = -1, c_at = -1;
    for (int i = 5; i < 35; ++i) {
        if (t_at < 0 && x[i] == 'T') t_at = i;
        if (c_at < 0 && x[i] == 'C') c_at = i;
    }
    REQUIRE(t_at >= 0);
    REQUIRE(c_at >= 0);
    r[t_at] = 'C'; // the reverse direction: read C over reference T (OT)
    r[c_at] = 'T'; // a conversion
    const std::string fq = write_fq("modes.fq", {r});
    auto neutral = run_aln({"--meth=taps"}, {fq});
    auto genomic = run_aln({"--meth=taps", "--meth-scoring", "genomic"}, {fq});
    auto collapsed = run_aln({"--meth=taps", "--meth-scoring", "collapsed"}, {fq});
    REQUIRE(neutral.size() == 1);
    REQUIRE(genomic.size() == 1);
    REQUIRE(collapsed.size() == 1);
    CHECK(neutral[0].pos == 301);
    CHECK(tag_value(neutral[0], "NM:i") == "1");
    CHECK(tag_value(genomic[0], "NM:i") == "1");
    CHECK(tag_value(collapsed[0], "NM:i") == "0");
}

TEST_CASE("aln_meth: --meth-context cpg frees a conversion only at a CpG")
{
    const std::string &c2 = ref().c2;
    int cpg = -1, other = -1;
    for (int i = 820; i < 980; ++i) {
        if (c2[i] != 'C') continue;
        if (c2[i + 1] == 'G' && cpg < 0) cpg = i;
        if (c2[i + 1] != 'G' && other < 0 && i > cpg + 3) other = i;
    }
    REQUIRE(cpg >= 0);
    REQUIRE(other >= 0);
    const int st = (cpg < other ? cpg : other) - 10;
    std::string x = c2.substr(st, 50);
    std::string at_cpg = x, at_other = x;
    at_cpg[cpg - st] = 'T';
    at_other[other - st] = 'T';
    const std::string fq = write_fq("cpg.fq", {at_cpg, at_other});
    auto all = run_aln({"--meth=taps"}, {fq});
    auto ctx = run_aln({"--meth=5base", "--meth-context", "cpg"}, {fq});
    REQUIRE(all.size() == 2);
    REQUIRE(ctx.size() == 2);
    CHECK(all[0].tag["NM:i"] == "0");
    CHECK(all[1].tag["NM:i"] == "0");
    CHECK(ctx[0].tag["NM:i"] == "0");
    CHECK(ctx[1].tag["NM:i"] == "1");
    CHECK(ctx[0].pos == st + 1);
    CHECK(ctx[1].pos == st + 1);
}

TEST_CASE("aln_meth: cpg's partner base beyond the read end still makes a CpG")
{
    // a read ending on the C of a CpG: the G is outside the read, read from the reference (as mem)
    const std::string &c1 = ref().c1;
    int cpg = -1;
    for (int i = 900; i < 1400; ++i)
        if (c1[i] == 'C' && c1[i + 1] == 'G') { cpg = i; break; }
    REQUIRE(cpg >= 0);
    std::string r = c1.substr(cpg - 39, 40); // its last base is the CpG's C
    r[39] = 'T';
    auto v = run_aln({"--meth=taps", "--meth-context", "cpg", "-n", "2", "-l", "1024"}, {write_fq("cpgedge.fq", {r})});
    REQUIRE(v.size() == 1);
    CHECK(v[0].pos == cpg - 38);
    CHECK(v[0].tag["NM:i"] == "0");
}

TEST_CASE("aln_meth: XM call polarity follows the chemistry")
{
    const std::string &c2 = ref().c2;
    int cpg = -1;
    for (int i = 1200; i < 1400; ++i)
        if (c2[i] == 'C' && c2[i + 1] == 'G') { cpg = i; break; }
    REQUIRE(cpg >= 0);
    const int st = cpg - 20;
    const std::string x = c2.substr(st, 45);
    // em-seq: every non-CpG C converts; the CpG C is kept, so it is methylated (Z)
    std::string em = conv(x, 'C', 'T');
    em[cpg - st] = 'C';
    // TAPS: only the methylated CpG C converts (to T), so it is methylated (Z) too
    std::string taps = x;
    taps[cpg - st] = 'T';
    auto e = run_aln({"--meth"}, {write_fq("xm_em.fq", {em})});
    auto t = run_aln({"--meth=taps"}, {write_fq("xm_taps.fq", {taps})});
    REQUIRE(e.size() == 1);
    REQUIRE(t.size() == 1);
    REQUIRE(tag_value(e[0], "XM:Z").size() == x.size());
    REQUIRE(tag_value(t[0], "XM:Z").size() == x.size());
    CHECK(tag_value(e[0], "XM:Z")[cpg - st] == 'Z');
    CHECK(tag_value(t[0], "XM:Z")[cpg - st] == 'Z');
    // the same base read the other way round is unmethylated under each chemistry
    std::string em_u = conv(x, 'C', 'T');
    std::string taps_u = x;
    auto eu = run_aln({"--meth"}, {write_fq("xm_emu.fq", {em_u})});
    auto tu = run_aln({"--meth=taps"}, {write_fq("xm_tapsu.fq", {taps_u})});
    REQUIRE(eu.size() == 1);
    REQUIRE(tu.size() == 1);
    REQUIRE(tag_value(eu[0], "XM:Z").size() == x.size());
    REQUIRE(tag_value(tu[0], "XM:Z").size() == x.size());
    CHECK(tag_value(eu[0], "XM:Z")[cpg - st] == 'z');
    CHECK(tag_value(tu[0], "XM:Z")[cpg - st] == 'z');
}

TEST_CASE("aln_meth: --meth-tags selects the tags; none drops all three")
{
    const std::string fq = write_fq("tags.fq", {conv(ref().c1.substr(700, 40), 'C', 'T')});
    auto none = run_aln({"--meth", "--meth-tags", "none"}, {fq});
    REQUIRE(none.size() == 1);
    CHECK(none[0].tag.count("XR:Z") == 0);
    CHECK(none[0].tag.count("XG:Z") == 0);
    CHECK(none[0].tag.count("XM:Z") == 0);
    auto some = run_aln({"--meth", "--meth-tags", "XR,XG"}, {fq});
    REQUIRE(some.size() == 1);
    CHECK(some[0].tag.count("XR:Z") == 1);
    CHECK(some[0].tag.count("XG:Z") == 1);
    CHECK(some[0].tag.count("XM:Z") == 0);
}

TEST_CASE("aln_meth: converted repeats do not trigger a rescored -R early stop")
{
    const std::string fq = write_fq("rstop.fq", {ref().read_r});
    auto v = run_aln({"--meth=taps", "-n", "2", "-l", "1024"}, {fq});
    REQUIRE(v.size() == 1);
    CHECK(v[0].rname == "c3");
    CHECK(v[0].pos == ref().x_pos);
    CHECK(tag_value(v[0], "XT:A") == "U");
    CHECK(tag_value(v[0], "X0:i") == "1");
    CHECK(tag_value(v[0], "X1:i") == "31"); // actual next reward stratum among generated loci
    CHECK(v[0].mapq == 8); // 23 - round(4.343*ln(31))
    // under collapsed scoring the repeats are exact: a repeat, MAPQ 0
    auto c = run_aln({"--meth=taps", "--meth-scoring", "collapsed", "-n", "2", "-l", "1024"}, {fq});
    REQUIRE(c.size() == 1);
    CHECK(tag_value(c[0], "XT:A") == "R");
    CHECK(c[0].mapq == 0);
}

TEST_CASE("aln_meth: -N preserves best placements while retaining more worse costs")
{
    // reads from c1/c2 with conversions, substitutions and reverse-direction cells
    std::vector<std::string> reads;
    uint32_t st = 99;
    for (int i = 0; i < 150; ++i) {
        const std::string &c = i & 1 ? ref().c2 : ref().c1;
        st = st * 1103515245u + 12345u;
        const int pos = (int)((st >> 8) % (c.size() - 60));
        const int len = 25 + (int)((st >> 4) % 30);
        std::string r = c.substr(pos, len);
        for (char &b : r) {
            st = st * 1103515245u + 12345u;
            const unsigned u = (st >> 16) % 100;
            if (b == 'C' && u < 40) b = 'T';
            else if (b == 'T' && u < 3) b = 'C';
            else if (u >= 98) b = "ACGT"[(st >> 9) & 3];
        }
        if (i % 3 == 0) r = aln_test::revcomp(r);
        reads.push_back(r);
    }
    const std::string fq = write_fq("oracle.fq", reads);
    for (const char *mode : {"--meth", "--meth=taps"}) {
        for (const char *n : {"2", "0.04"}) {
            auto a = run_aln({mode, "-n", n, "-l", "1024", "-o", "2"}, {fq});
            auto b = run_aln({mode, "-n", n, "-l", "1024", "-o", "2", "-N"}, {fq});
            REQUIRE(a.size() == b.size());
            for (size_t i = 0; i < a.size(); ++i) {
                CAPTURE(mode); CAPTURE(n); CAPTURE(i);
                CHECK(a[i].flag == b[i].flag);
                CHECK(a[i].rname == b[i].rname);
                CHECK(a[i].pos == b[i].pos);
                CHECK(a[i].cigar == b[i].cigar);
                CHECK(a[i].mapq >= b[i].mapq);
                REQUIRE(a[i].tag.count("X0:i") == b[i].tag.count("X0:i"));
                if (a[i].tag.count("X0:i")) CHECK(tag_value(a[i], "X0:i") == tag_value(b[i], "X0:i"));
                REQUIRE(a[i].tag.count("XT:A") == b[i].tag.count("XT:A"));
                if (a[i].tag.count("XT:A")) CHECK(tag_value(a[i], "XT:A") == tag_value(b[i], "XT:A"));
                REQUIRE(a[i].tag.count("X1:i") == b[i].tag.count("X1:i"));
                if (a[i].tag.count("X1:i"))
                    CHECK(std::stoi(tag_value(a[i], "X1:i")) <= std::stoi(tag_value(b[i], "X1:i")));
            }
        }
    }
}

TEST_CASE("aln_meth: the records do not depend on the thread count")
{
    // 192 distinct reads: three of the 64-read chunks the search threads claim, so -t 3 has work
    // for every thread
    std::vector<std::string> reads;
    for (int i = 0; i < 192; ++i) reads.push_back(conv(ref().c1.substr(15 * i, 36), 'C', 'T'));
    const std::string fq = write_fq("threads.fq", reads);
    run_aln({"--meth", "-t", "1"}, {fq});
    const std::string one = slurp_out();
    run_aln({"--meth", "-t", "3"}, {fq});
    const std::string three = slurp_out();
    // the @PG lines name the thread count; compare the records
    auto body = [](const std::string &s) {
        std::string b;
        std::stringstream ss(s);
        for (std::string l; std::getline(ss, l);)
            if (l.compare(0, 3, "@PG") != 0) b += l + "\n";
        return b;
    };
    CHECK(body(one) == body(three));
}

TEST_CASE("aln_meth: an index built with --meth still aligns without --meth (stock records)")
{
    const std::string fq = write_fq("stock.fq", {ref().c1.substr(200, 40)});
    auto v = run_aln({}, {fq});
    REQUIRE(v.size() == 1);
    CHECK(v[0].pos == 201);
    CHECK(has_tag(v[0], "XM:i"));
    CHECK(v[0].tag.count("XR:Z") == 0);
}

namespace {

// A single-end --meth run through the library, so the test hooks of aln_opt_t can be set.
std::string run_se_api(const std::string &fq, int share_dup_search, int batch_size, const char *chem_taps,
                       const std::string &mate = "", int barcode_len = 0)
{
    namespace ba = bwa_aln;
    ba::aln_opt_t o;
    ba::aln_opt_init(&o);
    o.gap.mode |= barcode_len << 24; // -B, unchecked: the library entry points take any length
    o.share_dup_search = share_dup_search;
    o.batch_size = batch_size;
    o.gap.max_diff = 2;
    o.gap.fnr = -1.0;
    o.gap.seed_len = 1024;
    ba::aln_meth_opt_t mo;
    mo.m.meth_mode = 1;
    mo.m.meth_chem = chem_taps ? 1 : 0;
    mo.m.meth_scoring = chem_taps ? MEM_METH_SCORING_NEUTRAL : MEM_METH_SCORING_COLLAPSED;
    ba::aln_meth_opt_finish(&mo);
    const std::string conv = ref().fa + ".meth";
    FMI_search fmi(conv.c_str());
    fmi.load_index(false, 1);
    ba::bwt_t *bwt = ba::bwt_from_fmi(&fmi);
    ba::aln_meth_ref_t mr;
    REQUIRE(ba::aln_meth_ref_load(&mr, ref().fa.c_str(), fmi.idx->bns, 1) == 0);
    o.meth = &mo;
    o.meth_ref = &mr;
    const std::string out_path = ref().dir + "/api.sam";
    ba::aln_out_t *out = ba::aln_out_open(out_path.c_str(), 0, 0, "", 1);
    REQUIRE(out != nullptr);
    if (mate.empty()) REQUIRE(ba::aln_run_se_meth(&o, bwt, fq.c_str(), out) == 0);
    else REQUIRE(ba::aln_run_pe_meth(&o, bwt, fq.c_str(), mate.c_str(), out) == 0);
    REQUIRE(ba::aln_out_close(out) == 0);
    ba::aln_meth_ref_free(&mr);
    ba::bwt_destroy(bwt);
    std::ifstream in(out_path);
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

} // namespace

TEST_CASE("aln_meth: aln_run_pe_meth bounds the joined barcode to BWA_MAX_BCLEN")
{
    // main_aln limits -B to BWA_MAX_BCLEN / 2 with two read files, but aln_run_pe_meth is a
    // public entry point: two 40-base barcodes join to 80 bases, more than the buffer holds, and
    // the join keeps the first BWA_MAX_BCLEN.
    const int l_bc = 40;
    REQUIRE(2 * l_bc > BWA_MAX_BCLEN);
    std::vector<std::string> one, two;
    for (int i = 0; i < 5; ++i) {
        const int start = 100 + 300 * i;
        one.push_back(conv(ref().c1.substr(start, 90), 'C', 'T'));
        two.push_back(conv(aln_test::revcomp(ref().c1.substr(start + 150, 90)), 'C', 'T'));
    }
    std::stringstream sam(run_se_api(write_fq("bc-one.fq", one), 1, 0x40000, nullptr,
                                     write_fq("bc-two.fq", two), l_bc));
    int n_rec = 0;
    for (std::string l; std::getline(sam, l);) {
        if (l.empty() || l[0] == '@') continue;
        const int pair = n_rec / 2;
        ++n_rec;
        CAPTURE(pair);
        const size_t at = l.find("\tBC:Z:");
        REQUIRE(at != std::string::npos);
        const size_t st = at + 6, en = l.find('\t', st);
        const std::string bc = l.substr(st, en == std::string::npos ? std::string::npos : en - st);
        const std::string want = (one[pair].substr(0, l_bc) + two[pair].substr(0, l_bc)).substr(0, BWA_MAX_BCLEN);
        CHECK(bc == want); // both mates carry the same joined barcode
    }
    CHECK(n_rec == 10);
}

TEST_CASE("aln_meth: identical reads share a search and still draw their own main hit (§4.4)")
{
    // a batch of mostly duplicated reads, so the duplicate pass runs, with repeats (several best
    // hits) so the draw matters
    std::vector<std::string> reads;
    for (int k = 0; k < 40; ++k) {
        reads.push_back(ref().read_r);                                   // a read with 31 repeat loci
        reads.push_back(conv(ref().c1.substr(37 * (k % 20), 36), 'C', 'T'));
    }
    const std::string fq = write_fq("dups.fq", reads);
    for (const char *chem : {(const char *)nullptr, "taps"}) {
        const std::string shared = run_se_api(fq, 1, 0x40000, chem);
        const std::string alone = run_se_api(fq, 0, 0x40000, chem);
        CHECK(shared == alone);
        // and the result does not depend on the batch boundaries
        const std::string small = run_se_api(fq, 1, 7, chem);
        CHECK(small == alone);
    }
}

TEST_CASE("aln_meth: aln_opt_init clears the --meth pointers of a reused aln_opt_t")
{
    namespace ba = bwa_aln;
    ba::aln_opt_t o;
    ba::aln_meth_opt_t mo;
    ba::aln_meth_ref_t mr;
    o.meth = &mo;
    o.meth_ref = &mr;
    ba::aln_opt_init(&o);
    CHECK(o.meth == nullptr);
    CHECK(o.meth_ref == nullptr);
}

TEST_CASE("aln_meth: a duplicate's search work is counted once, by the read it copies")
{
    // the stats line of a run (-v 4): SA lookups, windows scored, flushes
    auto work = [](const std::string &fq) {
        const int verbose = bwa_verbose;
        bwa_verbose = 4;
        std::string err;
        {
            bwa_tests::StderrCapture cap;
            run_se_api(fq, 1, 0x40000, nullptr);
            err = cap.text();
        }
        bwa_verbose = verbose;
        const size_t at = err.find("[M::aln] --meth read");
        REQUIRE(at != std::string::npos);
        REQUIRE(err.find("[M::aln] --meth read", at + 1) == std::string::npos); // one batch
        int n = 0;
        long lookups = 0, wrong = 0, scored = 0, dropped = 0, flushes = 0;
        REQUIRE(sscanf(err.c_str() + at, "[M::aln] --meth read %*d: %d reads, %ld SA lookups (%ld wrong-copy), "
                       "%ld unique windows scored (%ld dropped as unrepresentable), %ld flushes",
                       &n, &lookups, &wrong, &scored, &dropped, &flushes) == 6);
        return std::vector<long>{lookups, wrong, scored, dropped, flushes};
    };
    const std::string a = ref().read_r, b = conv(ref().c1.substr(500, 36), 'C', 'T');
    std::vector<std::string> dups;
    for (int k = 0; k < 20; ++k) { dups.push_back(a); dups.push_back(b); }
    const std::vector<long> once = work(write_fq("work_once.fq", {a, b}));
    CHECK(once[0] > 0);
    CHECK(work(write_fq("work_dups.fq", dups)) == once);
}

TEST_CASE("aln_meth: XM is built on the record's own contig when XA names another")
{
    std::string x = ref().c1.substr(2500, 50);
    int nc = 0;
    for (char c : x) nc += c == 'C';
    REQUIRE(nc > 3);
    auto v = run_aln({"--meth"}, {write_fq("xa_xm.fq", {conv(x, 'C', 'T')})});
    REQUIRE(v.size() == 1);
    CHECK(v[0].tag.count("XA:Z") == 1);
    const auto cpg=run_aln({"--meth=taps","--meth-context","cpg"},
                           {write_fq("xa-cpg.fq",{x})});
    REQUIRE(cpg.size()==1);
    const std::string &original=cpg[0].rname=="c1"?ref().c1:ref().c2;
    CHECK(tag_value(cpg[0],"XM:Z")==expected_xm(original,cpg[0].pos-1,x,true,true));
    CHECK(cpg[0].tag.count("XA:Z")==1);
    const std::string &xm = tag_value(v[0], "XM:Z");
    int calls = 0;
    for (char c : xm) calls += c != '.';
    CHECK(calls == nc); // every C of the top strand gets a call (z/x/h), wherever the read is placed
}

TEST_CASE("aln_meth: jobs never allocate stock hit arrays, including immediate completion")
{
    namespace ba = bwa_aln;
    const std::string conv_prefix = ref().fa + ".meth";
    FMI_search fmi(conv_prefix.c_str());
    fmi.load_index(false, 1);
    ba::bwt_t *bwt = ba::bwt_from_fmi(&fmi);
    ba::gap_pool_t pool;
    ba::gap_init_pool(&pool, 2);
    ba::gap_opt_t *opt = ba::gap_init_opt();
    opt->fnr = -1; opt->max_diff = 0;
    ba::gap_job_t job = {};
    job.opt = *opt;
    job.stack = ba::gap_init_cstack(0, opt->max_gapo, opt->max_gape, opt, &pool);
    const unsigned char ambiguous[] = {4, 4, 4, 4};
    const unsigned char bases[] = {0, 1, 2, 3};
    job.len = 4;
    std::vector<ba::bwt_width_t> width(5);
    job.width = width.data();
    for (int i = 0; i < 100; ++i) {
        job.seq = i % 2 ? bases : ambiguous;
        const int done = ba::gap_job_start_meth(&job, bwt);
        CHECK(done == (i % 2 == 0));
        CHECK(job.aln == nullptr);
        CHECK(job.m_aln == 0);
        CHECK(job.n_aln == 0);
        // Also exercise cleanup of a started but abandoned job, not just normal completion.
        ba::gap_release_cstack(job.stack);
    }
    ba::gap_destroy_cstack(job.stack);
    CHECK(pool.n_alloc <= pool.max_free);
    ba::gap_destroy_pool(&pool);
    free(opt);
    ba::bwt_destroy(bwt);
}

TEST_CASE("aln_meth: slot reuse handles mapped unmapped and early-complete reads")
{
    // Every read is distinct, so duplicate sharing collapses nothing and each read takes a search
    // slot of its own; mapped, early-complete (all but ten bases N) and mapped reads alternate, so
    // reused slots go from mapped to unmapped and from unmapped to mapped.
    std::vector<std::string> reads;
    for (int k = 0; k < 90; ++k) {
        reads.push_back(conv(ref().c1.substr(500 + k, 50), 'C', 'T'));
        reads.push_back(ref().c1.substr(1500 + k, 10) + std::string(40, 'N'));
        reads.push_back(conv(aln_test::revcomp(ref().c2.substr(600 + k, 50)), 'C', 'T'));
    }
    const std::string fq = write_fq("slot-reuse.fq", reads);
    auto check = [&](const std::vector<Rec> &records) {
        REQUIRE(records.size() == reads.size());
        for (size_t i = 0; i < records.size(); ++i) {
            CAPTURE(i);
            const int k = (int)(i / 3);
            CHECK(records[i].name == "r" + std::to_string(i));
            if (i % 3 == 1) {
                CHECK((records[i].flag & 4) != 0);
                CHECK(records[i].rname == "*");
                CHECK(tag_value(records[i], "XR:Z") == "CT");
                CHECK(records[i].tag.count("XM:Z") == 0);
                CHECK(records[i].tag.count("XG:Z") == 0);
            } else {
                const bool fwd = i % 3 == 0;
                CHECK(records[i].flag == (fwd ? 0 : 16));
                CHECK(records[i].rname == (fwd ? "c1" : "c2"));
                CHECK(records[i].pos == (fwd ? 501 : 601) + k);
                CHECK(records[i].cigar == "50M");
            }
        }
    };
    SUBCASE("three threads") { check(run_aln({"--meth", "-t", "3"}, {fq})); }
    SUBCASE("duplicate sharing off") {
        std::stringstream sam(run_se_api(fq, 0, 0x40000, nullptr));
        check(parse_records(sam));
    }
}

TEST_CASE("aln_meth: stock cost selects a lower reward refined CIGAR")
{
    const std::string fq = write_fq("cost-primary.fq", {ref().affine_read});
    for (bool nonstop : {false, true}) {
        std::vector<std::string> options = {"--meth=taps", "-n", "6", "-l", "1024",
                                            "-M", "6", "-O", "1", "-E", "1"};
        if (nonstop) options.push_back("-N");
        const auto hits = run_aln(options, {fq});
        REQUIRE(hits.size() == 1);
        CHECK(hits[0].rname == "affine_cheap"); // cost 5 vs 12, despite rewards 45 vs 52
        CHECK(hits[0].cigar.find("5I") != std::string::npos);
        CHECK(tag_value(hits[0], "NM:i") == "5");
        CHECK(tag_value(hits[0], "MD:Z") == "55");
        CHECK(tag_value(hits[0], "X0:i") == "1");
        // Ordinary cost window ends at 11; -N also retains four worse-cost loci,
        // including the cost-12 mismatch placement and shifted gapped placements.
        CHECK(tag_value(hits[0], "X1:i") == (nonstop ? "4" : "0"));
        CHECK(hits[0].mapq == (nonstop ? 17 : 37));
    }
    for (const auto &penalties : {std::vector<std::string>{"0","1","1"},
                                   std::vector<std::string>{"6","0","0"}}) {
        const auto zero = run_aln({"--meth=taps", "-n", "6", "-l", "1024",
                                   "-M", penalties[0], "-O", penalties[1], "-E", penalties[2]}, {fq});
        REQUIRE(zero.size()==1);
        if (penalties[0]=="0") {
            CHECK(zero[0].rname=="affine_expensive");
            CHECK(zero[0].cigar=="60M");
            CHECK(tag_value(zero[0],"NM:i")=="2");
        } else {
            CHECK(zero[0].rname=="affine_cheap");
            CHECK(zero[0].cigar.find("5I")!=std::string::npos);
            CHECK(tag_value(zero[0],"NM:i")=="5");
        }
        CHECK(zero[0].mapq==37);
    }
    // With expensive gaps the same refined candidates instead prefer two mismatches.
    const auto costly_gap = run_aln({"--meth=taps", "-n", "6", "-l", "1024",
                                    "-M", "3", "-O", "11", "-E", "4"}, {fq});
    REQUIRE(costly_gap.size() == 1);
    CHECK(costly_gap[0].rname == "affine_expensive");
    CHECK(costly_gap[0].cigar == "60M");
    CHECK(tag_value(costly_gap[0], "NM:i") == "2");
    CHECK(costly_gap[0].mapq == 37);
}

TEST_CASE("aln_meth: stock MAPQ mismatch limit unique and tied placements")
{
    const std::string exact = conv(ref().c1.substr(500, 50), 'C', 'T');
    auto unique = run_aln({"--meth", "-n", "2", "-l", "1024"},
                          {write_fq("mapq-unique.fq", {exact})});
    REQUIRE(unique.size() == 1);
    CHECK(tag_value(unique[0], "X0:i") == "1");
    CHECK(tag_value(unique[0], "X1:i") == "0");
    CHECK(unique[0].mapq == 37);
    std::string changed = exact;
    changed[25] = (changed[25] == 'A' ? 'G' : 'A');
    const auto limit = run_aln({"--meth", "-n", "1", "-l", "1024"},
                              {write_fq("mapq-limit.fq", {changed})});
    REQUIRE(limit.size() == 1);
    CHECK(tag_value(limit[0], "NM:i") == "1");
    CHECK(limit[0].mapq == 25);
    const auto tied = run_aln({"--meth", "-n", "2", "-l", "1024"},
                             {write_fq("mapq-tied.fq", {conv(ref().c1.substr(2500,50),'C','T')})});
    REQUIRE(tied.size() == 1);
    CHECK(tag_value(tied[0], "X0:i") == "2");
    CHECK(tied[0].mapq == 0);
}

TEST_CASE("aln_meth: parallel finish preserves nonzero second-stratum MAPQ")
{
    const std::string fq=write_fq("mapq-workers.fq",std::vector<std::string>(160,ref().read_r));
    run_aln({"--meth=taps","-n","2","-l","1024","-t","1"},{fq});
    const std::string serial=slurp_out();
    const auto parallel=run_aln({"--meth=taps","-n","2","-l","1024","-t","4"},{fq});
    CHECK(slurp_out()==serial);
    REQUIRE(parallel.size()==160);
    for(const auto &record:parallel) {
        CHECK(record.mapq==8);
        CHECK(tag_value(record,"X1:i")=="31");
    }
}

TEST_CASE("aln_meth: trim merges existing terminal clips safely for primary and XA")
{
    namespace ba=bwa_aln;
    ba::aln_meth_ref_t mr;
    const std::string prefix=ref().fa+".meth";
    FMI_search fmi(prefix.c_str());fmi.load_index(false,1);
    REQUIRE(ba::aln_meth_ref_load(&mr,ref().fa.c_str(),fmi.idx->bns,1)==0);
    ba::aln_meth_opt_t mo;ba::aln_meth_opt_finish(&mo);
    ba::gap_opt_t *opt=ba::gap_init_opt();ba::bwase_initialize();
    for(int reverse : {0,1}) {
        ba::bwa_seqio_t *input=ba::bwa_seq_open(write_fq("owned-trim.fq",{std::string(17060,'A')}).c_str());
        int n=0;auto *read=ba::bwa_read_seq(input,1,&n,BWA_MODE_COMPREAD,0);
        REQUIRE(n==1);read->len=read->clip_len=60;
        ba::aln_meth_read_t result;result.main=0;result.c1=2;
        for(int strand : {reverse,1-reverse}) {
            ba::aln_meth_occ_t o;o.rid=0;o.hyp=1;o.strand=strand;
            o.pos=mr.bns->anns[0].offset+500;o.alignment.rb=o.pos;o.alignment.re=o.pos+57;
            o.alignment.cigar=strand ? std::vector<uint32_t>{3u<<4|4,57u<<4}
                                     : std::vector<uint32_t>{57u<<4,3u<<4|4};
            result.occ.push_back(std::move(o));
        }
        ba::aln_meth_finish_se(&mr,&mo,opt,3,1,read,&result);
        auto span=[](const ba::bwa_cigar_t *c,int count) {
            int total=0;for(int k=0;k<count;++k) if(__cigar_op(c[k])!=2) total+=__cigar_len(c[k]);
            return total;
        };
        CHECK(span(read->cigar,read->n_cigar)==17060);
        REQUIRE(read->n_multi==1);
        CHECK(span(read->multi[0].cigar,read->multi[0].n_cigar)==17060);
        auto clipped=[](const ba::bwa_cigar_t *c,int count) {
            int total=0;for(int k=0;k<count;++k) if(__cigar_op(c[k])==3) total+=__cigar_len(c[k]);
            return total;
        };
        CHECK(clipped(read->cigar,read->n_cigar)==17003);
        CHECK(clipped(read->multi[0].cigar,read->multi[0].n_cigar)==17003);
        CHECK((int)read->len==17060);
        ba::bwa_free_read_seq(n,read);ba::bwa_seq_close(input);
    }
    free(opt);ba::aln_meth_ref_free(&mr);
}

TEST_CASE("aln_meth: an XA alternative too long for the legacy edit field drops XA, not the primary")
{
    namespace ba=bwa_aln;
    ba::aln_meth_ref_t mr;
    const std::string prefix=ref().fa+".meth";
    FMI_search fmi(prefix.c_str());fmi.load_index(false,1);
    REQUIRE(ba::aln_meth_ref_load(&mr,ref().fa.c_str(),fmi.idx->bns,1)==0);
    ba::aln_meth_opt_t mo;ba::aln_meth_opt_finish(&mo);
    ba::gap_opt_t *opt=ba::gap_init_opt();ba::bwase_initialize();
    ba::bwa_seqio_t *input=ba::bwa_seq_open(write_fq("xa-edits.fq",{ref().c1.substr(100,400)}).c_str());
    int n=0;auto *read=ba::bwa_read_seq(input,1,&n,BWA_MODE_COMPREAD,0);
    REQUIRE(n==1);
    ba::aln_meth_read_t result;result.main=0;result.c1=1;result.c2=1;
    const int64_t c1_off=mr.bns->anns[0].offset;
    ba::aln_meth_occ_t primary;primary.rid=0;primary.hyp=1;primary.strand=0;
    primary.pos=c1_off+100;primary.alignment.rb=primary.pos;primary.alignment.re=primary.pos+400;
    primary.alignment.cigar={400u<<4};
    result.occ.push_back(std::move(primary));
    // An internal 256-base insertion: one gap open and 255 extensions fit a record, but the
    // alternative's edit distance (at least 256) does not fit XA's 8-bit field.
    ba::aln_meth_occ_t alt;alt.rid=0;alt.hyp=1;alt.strand=0;alt.cost=1;
    alt.pos=c1_off+1000;alt.alignment.rb=alt.pos;alt.alignment.re=alt.pos+144;
    alt.alignment.cigar={72u<<4,256u<<4|1,72u<<4};
    result.occ.push_back(std::move(alt));
    ba::aln_meth_finish_se(&mr,&mo,opt,3,1,read,&result);
    CHECK((int)read->type==BWA_TYPE_UNIQUE);
    CHECK((int64_t)read->pos==c1_off+100);
    REQUIRE(read->n_cigar==1);
    CHECK(__cigar_op(read->cigar[0])==0);
    CHECK(__cigar_len(read->cigar[0])==400);
    CHECK((int)read->nm==0);
    CHECK(read->n_multi==0);
    CHECK(read->multi==nullptr);
    ba::bwa_free_read_seq(n,read);ba::bwa_seq_close(input);
    free(opt);ba::aln_meth_ref_free(&mr);
}

TEST_CASE("aln_meth: retained cost window counts all worse strata with wide arithmetic")
{
    namespace ba = bwa_aln;
    ba::gap_opt_t *opt=ba::gap_init_opt();
    opt->s_mm=3;
    ba::aln_meth_group_t group; group.opt=opt;
    ba::aln_meth_read_t result;
    ba::aln_meth_slot_t slot; slot.grp=&group; slot.res=&result; slot.keep=1;
    auto populate=[&](int64_t best) {
        result=ba::aln_meth_read_t();
        for(int k=0;k<5;++k) {
            ba::aln_meth_occ_t o;
            o.rid=0; o.hyp=1; o.strand=0; o.pos=100+k;
            o.cost=best+k;
            o.alignment.reward=100+k; // reward deliberately points in the opposite direction
            result.occ.push_back(std::move(o));
        }
    };
    populate(0);ba::aln_meth_slot_finish(&slot);
    CHECK(result.c1==1); CHECK(result.c2==3); CHECK(result.occ.size()==4);
    CHECK(result.occ[result.main].pos==100);
    opt->mode|=BWA_MODE_NONSTOP;
    populate(0);ba::aln_meth_slot_finish(&slot);
    CHECK(result.c1==1); CHECK(result.c2==4); CHECK(result.occ.size()==5);
    opt->mode&=~BWA_MODE_NONSTOP;opt->s_mm=INT_MAX;
    populate(INT_MAX-4);ba::aln_meth_slot_finish(&slot);
    CHECK(result.c1==1); CHECK(result.c2==4); CHECK(result.occ.size()==5);
    free(opt);
}

TEST_CASE("aln_meth: a dropped unrepresentable window that could compete marks the read incomplete")
{
    namespace ba = bwa_aln;
    ba::gap_opt_t *opt=ba::gap_init_opt();
    opt->s_mm=3;
    ba::aln_meth_group_t group; group.opt=opt;
    ba::aln_meth_read_t result;
    ba::aln_meth_slot_t slot; slot.grp=&group; slot.res=&result; slot.keep=1;
    // Retained costs 10 and 11, so the retained window ends at 10 + 3. A dropped window inside it
    // could have been a best or retained worse locus; one beyond it could not.
    for(const auto &c : std::vector<std::pair<int64_t,int>>{{9,ba::ALN_METH_BIND_UNREP},{13,ba::ALN_METH_BIND_UNREP},{14,0}}) {
        CAPTURE(c.first);
        result=ba::aln_meth_read_t();
        for(int k=0;k<2;++k) {
            ba::aln_meth_occ_t o; o.rid=0; o.hyp=1; o.strand=0; o.pos=100+k; o.cost=10+k;
            result.occ.push_back(std::move(o));
        }
        result.n_unrepresentable=1; slot.unrep_cost=c.first;
        ba::aln_meth_slot_finish(&slot);
        CHECK(result.flags==c.second);
        CHECK(result.c1==1); CHECK(result.c2==1);
    }
    // No dropped window: complete, even under -N's unbounded window.
    opt->mode|=BWA_MODE_NONSTOP;
    result=ba::aln_meth_read_t(); slot.unrep_cost=INT64_MAX;
    ba::aln_meth_occ_t o; o.rid=0; o.hyp=1; o.strand=0; o.pos=100; result.occ.push_back(std::move(o));
    ba::aln_meth_slot_finish(&slot);
    CHECK(result.flags==0);
    free(opt);
}

TEST_CASE("aln_meth: window cache scores spans once and caps cannot claim unique strata")
{
    namespace ba = bwa_aln;
    const std::string fq = write_fq("cache.fq", {conv(ref().c1.substr(500,50),'C','T')});
    ba::bwa_seqio_t *input = ba::bwa_seq_open(fq.c_str());
    REQUIRE(input != nullptr);
    int n=0;ba::bwa_seq_t *reads=ba::bwa_read_seq(input,1,&n,BWA_MODE_COMPREAD,0);
    REQUIRE(n==1);REQUIRE(reads!=nullptr);
    const std::string prefix=ref().fa+".meth";
    FMI_search fmi(prefix.c_str());fmi.load_index(false,1);
    ba::bwt_t *bwt=ba::bwt_from_fmi(&fmi);
    ba::aln_meth_ref_t mr;
    REQUIRE(ba::aln_meth_ref_load(&mr,ref().fa.c_str(),fmi.idx->bns,1)==0);
    ba::aln_meth_opt_t mo;ba::aln_meth_opt_finish(&mo);
    ba::gap_opt_t *opt=ba::gap_init_opt();
    ba::aln_meth_group_t group;group.bwt=bwt;group.ref=&mr;group.mo=&mo;group.opt=opt;
    group.n_slots=1;group.slot.resize(1);
    ba::gap_slot_core_t core={};group.cores=&core;
    auto &slot=group.slot[0];slot.grp=&group;slot.job=&core.job;
    ba::aln_meth_read_t result;
    ba::aln_meth_slot_start(&slot,reads,&result,0,3,1);
    int strand=0; bwtint_t row=(bwtint_t)-1;
    // Find a row for the forward C->T copy at c1 offset 500; independent of search hit order.
    for(bwtint_t k=0;k<=bwt->seq_len;++k) {
        auto pos=ba::bwa_pos_from_sa(mr.conv_bns,ba::bwt_sa(bwt,k),50,&strand);
        if(pos!=(bwtint_t)-1 && strand==0 && bns_pos2rid(mr.conv_bns,pos)==1 &&
           pos==(bwtint_t)(mr.conv_bns->anns[1].offset+500)){row=k;break;}
    }
    REQUIRE(row!=(bwtint_t)-1);
    ba::gap_entry_t exact={},longer={};longer.n_del=1;
    ba::aln_meth_hit(&slot,&exact,row,row);
    ba::aln_meth_hit(&slot,&exact,row,row);
    ba::aln_meth_hit(&slot,&longer,row,row);
    REQUIRE(slot.pend.size()==2);
    core.state=ba::GAP_SLOT_PARKED;ba::aln_meth_flush_parked(&group);
    CHECK(result.n_rescored==2);CHECK(slot.windows.size()==2);
    // The two spans have the same normalized start; their minimum-cost alignment wins.
    REQUIRE(result.occ.size()==1);
    CHECK(result.occ[0].cost==0);
    CHECK(result.occ[0].n_mm==0);
    CHECK(result.occ[0].n_gapo==0);
    CHECK(result.occ[0].n_gape==0);
    CHECK(result.occ[0].n_conv==(int)std::count(ref().c1.begin()+500,ref().c1.begin()+550,'C'));
    ba::aln_meth_hit(&slot,&exact,row,row); // repeated interval across flushes
    CHECK(slot.pend.empty());
    ba::aln_meth_slot_finish(&slot);
    REQUIRE(result.occ.size()==1);CHECK(result.c1==1);
    // A terminal deletion is charged on the raw CIGAR even though display removes it.
    for (int penalty : {0, 11, INT_MAX}) {
        opt->s_gapo=penalty;
        ba::aln_meth_slot_start(&slot,reads,&result,3,3,1);
        ba::aln_meth_hit(&slot,&longer,row,row);
        core.state=ba::GAP_SLOT_PARKED;ba::aln_meth_flush_parked(&group);
        REQUIRE(result.occ.size()==1);
        CHECK(result.occ[0].n_mm==0);
        CHECK(result.occ[0].n_gapo==1);
        CHECK(result.occ[0].n_gape==0);
        CHECK(result.occ[0].n_del==1);
        CHECK(result.occ[0].cost==penalty);
        CHECK(((result.occ[0].alignment.raw_cigar.front()&15)==2 ||
               (result.occ[0].alignment.raw_cigar.back()&15)==2));
        CHECK(result.occ[0].alignment.cigar.size()==1);
        CHECK(result.occ[0].alignment.cigar[0]==(50u<<4));
    }
    opt->s_gapo=11;
    mo.p_best=0;mo.p_second=0;
    ba::aln_meth_slot_start(&slot,reads,&result,1,3,1);
    ba::aln_meth_hit(&slot,&exact,row,row);core.state=ba::GAP_SLOT_PARKED;
    ba::aln_meth_flush_parked(&group);ba::aln_meth_slot_finish(&slot);
    CHECK(result.flags!=0);CHECK(result.main==-1);CHECK(result.c1==0);
    mo.p_best=1;mo.p_second=0;
    ba::aln_meth_slot_start(&slot,reads,&result,2,3,1);
    ba::aln_meth_hit(&slot,&longer,row,row); // lower reward evaluated first
    ba::aln_meth_hit(&slot,&exact,row,row); // better window hidden behind the cap
    core.state=ba::GAP_SLOT_PARKED;
    ba::aln_meth_flush_parked(&group);ba::aln_meth_slot_finish(&slot);
    REQUIRE(result.main>=0);CHECK(result.flags!=0);CHECK(result.n_rescored==1);
    CHECK(result.occ[result.main].alignment.reward < 50);
    ba::bwase_initialize(); // standalone caller, before any finish workers
    ba::aln_meth_finish_se(&mr,&mo,opt,3,1,reads,&result);
    CHECK((int)reads[0].type==BWA_TYPE_REPEAT);CHECK((int)reads[0].mapQ==0);
    CHECK((int)reads[0].c1==0);CHECK((int)reads[0].c2==0);CHECK(reads[0].n_multi==0);
    kstring_t text={0,0,nullptr};
    ba::aln_meth_print_sam1(&text,&mr,&mo,reads,&result,1,nullptr,BWA_MODE_COMPREAD,30);
    std::string record(text.s,text.l);free(text.s);
    CHECK(record.find("X0:i:")==std::string::npos);
    CHECK(record.find("X1:i:")==std::string::npos);
    CHECK(record.find("XA:Z:")==std::string::npos);
    ba::aln_meth_ref_free(&mr);ba::bwt_destroy(bwt);free(opt);
    ba::bwa_free_read_seq(n,reads);CHECK(ba::bwa_seq_close(input)==0);
}

TEST_CASE("aln_meth: the CIGAR byte cap retains an incomplete placement and the next read still maps")
{
    namespace ba = bwa_aln;
    const std::string fq = write_fq("cigar-cap.fq", std::vector<std::string>(2, conv(ref().c1.substr(500, 50), 'C', 'T')));
    ba::bwa_seqio_t *input = ba::bwa_seq_open(fq.c_str());
    REQUIRE(input != nullptr);
    int n = 0;
    ba::bwa_seq_t *reads = ba::bwa_read_seq(input, 2, &n, BWA_MODE_COMPREAD, 0);
    REQUIRE(reads != nullptr);
    REQUIRE(n == 2);
    const std::string prefix = ref().fa + ".meth";
    FMI_search fmi(prefix.c_str());
    fmi.load_index(false, 1);
    ba::bwt_t *bwt = ba::bwt_from_fmi(&fmi);
    ba::aln_meth_ref_t mr;
    REQUIRE(ba::aln_meth_ref_load(&mr, ref().fa.c_str(), fmi.idx->bns, 1) == 0);
    ba::aln_meth_opt_t mo;
    ba::aln_meth_opt_finish(&mo);
    ba::gap_opt_t *opt = ba::gap_init_opt();
    ba::aln_meth_group_t group;
    group.bwt = bwt; group.ref = &mr; group.mo = &mo; group.opt = opt;
    group.n_slots = 1; group.slot.resize(1);
    ba::gap_slot_core_t core = {};
    group.cores = &core;
    auto &slot = group.slot[0];
    slot.grp = &group; slot.job = &core.job;
    ba::aln_meth_read_t results[2];
    ba::aln_meth_slot_start(&slot, reads, results, 0, 3, 1);
    bwtint_t row = (bwtint_t)-1;
    for (bwtint_t k = 0; k <= bwt->seq_len; ++k) {
        int strand = 0;
        const auto pos = ba::bwa_pos_from_sa(mr.conv_bns, ba::bwt_sa(bwt, k), 50, &strand);
        if (pos != (bwtint_t)-1 && strand == 0 && bns_pos2rid(mr.conv_bns, pos) == 1 &&
            pos == (bwtint_t)(mr.conv_bns->anns[1].offset + 500)) { row = k; break; }
    }
    REQUIRE(row < bwt->seq_len);
    // Seed the existing counter rather than allocating 64 MiB. An exact 50M candidate owns
    // one raw and one display operation (8 bytes); these are the adjacent byte boundaries.
    uint64_t remaining = 8;
    SUBCASE("one byte below the cap") { remaining = 9; }
    SUBCASE("exactly at the cap") { remaining = 8; }
    SUBCASE("one byte over the cap, without a retained placement") { remaining = 7; }
    const uint64_t cap = 64ULL << 20;
    slot.alignment_bytes = cap - remaining;
    ba::gap_entry_t exact = {}, longer = {};
    longer.n_del = 1;
    ba::aln_meth_hit(&slot, &exact, row, row);
    core.state = ba::GAP_SLOT_PARKED;
    ba::aln_meth_flush_parked(&group);
    const bool retained = remaining >= 8;
    if (retained) {
        CHECK(results[0].flags == 0);
        CHECK(slot.stop == 0);
        CHECK(slot.alignment_bytes == cap - remaining + 8);
        REQUIRE(results[0].occ.size() == 1);
        // Overflow on the first row of an interval: do not process its later rows or the
        // next pending interval once the cap binds.
        ba::aln_meth_hit(&slot, &longer, row, row + 1);
        ba::aln_meth_hit(&slot, &longer, row, row);
        core.state = ba::GAP_SLOT_PARKED;
        ba::aln_meth_flush_parked(&group);
    }
    CHECK(results[0].flags == ba::ALN_METH_BIND_CIGAR);
    CHECK(slot.stop == 1);
    CHECK(results[0].n_lookups == (retained ? 2 : 1));
    CHECK(slot.pend.empty());
    ba::aln_meth_slot_finish(&slot);
    CHECK((results[0].main >= 0) == retained);
    if (retained) CHECK(results[0].occ[results[0].main].pos == 500);
    // Reuse the same slot for a later read: the cap and stop state belong to one read only.
    ba::aln_meth_slot_start(&slot, reads + 1, results + 1, 1, 3, 1);
    CHECK(slot.stop == 0);
    CHECK(slot.alignment_bytes == 0);
    ba::aln_meth_hit(&slot, &exact, row, row);
    core.state = ba::GAP_SLOT_PARKED;
    ba::aln_meth_flush_parked(&group);
    ba::aln_meth_slot_finish(&slot);
    REQUIRE(results[1].main >= 0);
    CHECK(results[1].flags == 0);
    ba::bwase_initialize();
    ba::aln_meth_finish_se(&mr, &mo, opt, 3, 2, reads, results);
    CHECK((int)reads[0].mapQ == 0);
    CHECK((int)reads[0].c1 == 0);
    CHECK((int)reads[0].c2 == 0);
    CHECK(reads[0].n_multi == 0);
    CHECK((int)reads[1].type == BWA_TYPE_UNIQUE);
    CHECK((int)reads[1].mapQ > 0);
    kstring_t text = {0, 0, nullptr};
    ba::aln_meth_print_sam1(&text, &mr, &mo, reads, results, 1, nullptr, BWA_MODE_COMPREAD, 30);
    const std::string record(text.s, text.l);
    free(text.s);
    CHECK(record.find(retained ? "r0\t0\tc1\t501\t0\t50M\t" : "r0\t4\t*\t0\t0\t*\t") == 0);
    CHECK(record.find("X0:i:") == std::string::npos);
    CHECK(record.find("X1:i:") == std::string::npos);
    CHECK(record.find("XA:Z:") == std::string::npos);
    ba::aln_meth_ref_free(&mr);
    ba::bwt_destroy(bwt);
    free(opt);
    ba::bwa_free_read_seq(n, reads);
    CHECK(ba::bwa_seq_close(input) == 0);
}

TEST_CASE("aln_meth: paired-end records do not depend on threads and -N retains competitors")
{
    std::vector<std::string> r1, r2;
    uint32_t st = 5;
    for (int i = 0; i < 60; ++i) {
        st = st * 1103515245u + 12345u;
        const std::string &c = i & 1 ? ref().c2 : ref().c1;
        const int pos = (int)((st >> 8) % (c.size() - 260));
        const std::string f = c.substr(pos, 150 + (int)((st >> 4) % 100));
        if (i % 2 == 0) { // OT
            r1.push_back(conv(f.substr(0, 40), 'C', 'T'));
            r2.push_back(conv(aln_test::revcomp(f).substr(0, 40), 'G', 'A'));
        } else {          // OB
            r1.push_back(conv(aln_test::revcomp(f).substr(0, 40), 'C', 'T'));
            r2.push_back(conv(f.substr(0, 40), 'G', 'A'));
        }
    }
    // Repeat pairs also exercise compaction and independent ordinal draws in both streams.
    const std::vector<std::string> copy1(r1.begin(),r1.begin()+10),copy2(r2.begin(),r2.begin()+10);
    r1.insert(r1.end(),copy1.begin(),copy1.end());
    r2.insert(r2.end(),copy2.begin(),copy2.end());
    const std::string f1 = write_fq("pt1.fq", r1), f2 = write_fq("pt2.fq", r2);
    CHECK(run_se_api(f1,1,0x40000,nullptr,f2)==run_se_api(f1,0,0x40000,nullptr,f2));
    auto body = [](const std::string &s) {
        std::string b;
        std::stringstream ss(s);
        for (std::string l; std::getline(ss, l);)
            if (l.compare(0, 3, "@PG") != 0) b += l + "\n";
        return b;
    };
    run_aln({"--meth", "-t", "1"}, {f1, f2});
    const std::string one = body(slurp_out());
    auto v = run_aln({"--meth", "-t", "3"}, {f1, f2});
    CHECK(one == body(slurp_out()));
    int proper = 0;
    for (size_t i = 0; i + 1 < v.size(); i += 2) {
        CAPTURE(i);
        if (v[i].flag & 0x2) {
            ++proper;
            CHECK(v[i].tag["XG:Z"] == v[i + 1].tag["XG:Z"]); // both ends on the fragment's genome strand
        }
        CHECK(v[i].tag["XR:Z"] == "CT");
        CHECK(v[i + 1].tag["XR:Z"] == "GA");
    }
    CHECK(proper >= 50);
    auto a = run_aln({"--meth", "-n", "2", "-l", "1024"}, {f1, f2});
    auto b = run_aln({"--meth", "-n", "2", "-l", "1024", "-N"}, {f1, f2});
    REQUIRE(a.size() == b.size());
    for (size_t i = 0; i < a.size(); ++i) {
        CHECK(a[i].pos == b[i].pos);
        CHECK(a[i].flag == b[i].flag);
        CHECK(a[i].mapq >= b[i].mapq);
    }
}

TEST_CASE("aln_meth: pairing moves to stored second-locus CIGAR and counts atomically")
{
    namespace ba=bwa_aln;
    const std::string prefix=ref().fa+".meth";
    FMI_search fmi(prefix.c_str());fmi.load_index(false,1);
    ba::aln_meth_ref_t mr;REQUIRE(ba::aln_meth_ref_load(&mr,ref().fa.c_str(),fmi.idx->bns,1)==0);
    ba::aln_meth_opt_t mo;ba::aln_meth_opt_finish(&mo);
    auto *opt=ba::gap_init_opt();opt->mode|=BWA_MODE_NONSTOP;
    auto *pe=ba::bwa_init_pe_opt();pe->is_sw=0;
    ba::bwase_initialize();
    ba::bwa_seq_t *seqs[2];ba::bwa_seqio_t *input[2];int n[2];
    for(int j=0;j<2;++j) {
        input[j]=ba::bwa_seq_open(write_fq(j?"pair-move2.fq":"pair-move1.fq",{std::string(50,'A')}).c_str());
        seqs[j]=ba::bwa_read_seq(input[j],1,&n[j],BWA_MODE_COMPREAD,0);REQUIRE(n[j]==1);
    }
    ba::aln_meth_read_t result[2];ba::aln_meth_read_t *r[2]={result,result+1};
    const int64_t base=mr.bns->anns[0].offset;
    for(int j=0;j<2;++j) {
        ba::aln_meth_occ_t o;o.rid=0;o.hyp=1;o.strand=j;o.pos=base+(j?700:100);
        o.alignment.rb=o.pos;o.alignment.re=o.pos+50;
        o.alignment.cigar={50u<<4};result[j].occ.push_back(o);
        result[j].main=0;result[j].c1=1;
    }
    auto alternate=result[0].occ[0];alternate.pos=base+550;alternate.cost=3;
    alternate.n_mm=1;alternate.n_gapo=1;alternate.n_ins=1;
    alternate.alignment.rb=alternate.pos;alternate.alignment.re=alternate.pos+49;
    alternate.alignment.cigar={25u<<4,1u<<4|1,24u<<4};
    result[0].occ.push_back(alternate);result[0].c2=1;
    bool incomplete=false;
    SUBCASE("complete loci may move") {}
    SUBCASE("incomplete end is not promoted") {incomplete=true;result[0].flags=ba::ALN_METH_BIND_BEST;}
    ba::aln_meth_pe_stage1(opt,1,seqs,r);
    ba::isize_info_t ii={};ii.avg=ii.std=-1;
    uint64_t total[2]={},mapped[2]={};
    ba::aln_meth_pe_stage2(&mr,&mo,opt,pe,&ii,1,seqs,r,total,mapped);
    if(incomplete) {
        CHECK(result[0].main==0);CHECK((int)seqs[0]->pos==base+100);
        CHECK((int)seqs[0]->mapQ==0);CHECK((int)seqs[0]->c1==0);
        CHECK((seqs[0]->extra_flag&SAM_FPP)==0);CHECK(seqs[0]->n_multi==0);
    } else {
    CHECK(result[0].main==1);CHECK(result[0].hyp==1);
    CHECK((int)seqs[0]->pos==base+550);CHECK(seqs[0]->score==3);
    CHECK((int)seqs[0]->n_mm==1);CHECK((int)seqs[0]->n_gapo==1);
    REQUIRE(seqs[0]->n_cigar==3); // the stored 25M1I24M, operation by operation
    const int ops[3]={FROM_M,FROM_I,FROM_M},lens[3]={25,1,24};
    for(int k=0;k<3;++k) {
        CAPTURE(k);CHECK(__cigar_op(seqs[0]->cigar[k])==ops[k]);CHECK(__cigar_len(seqs[0]->cigar[k])==lens[k]);
    }
    CHECK((int)seqs[0]->mapQ==29);CHECK((int)seqs[0]->seQ==0);
    CHECK((seqs[0]->extra_flag&SAM_FPP)!=0);
    CHECK(std::string(seqs[0]->md).find('^')==std::string::npos);
    }
    for(int j=0;j<2;++j) {ba::bwa_free_read_seq(1,seqs[j]);ba::bwa_seq_close(input[j]);}
    free(pe);free(opt);ba::aln_meth_ref_free(&mr);
}

TEST_CASE("aln_meth: ordinary mate rescue owns its CIGAR and chemistry tags")
{
    std::vector<std::string> one,two;
    for(int i=0;i<45;++i) {
        const int start=100+50*i,span=180+i%21;
        one.push_back(ref().c1.substr(start,50));
        two.push_back(aln_test::revcomp(ref().c1.substr(start+span-50,50)));
    }
    // Search cannot propose the last mate; stock-shaped local rescue can recover it.
    for(int at : {15,25,35}) two.back()[at]=two.back()[at]=='A'?'C':'A';
    const auto v=run_aln({"--meth","-n","0","-l","1024"},
                        {write_fq("rescue-one.fq",one),write_fq("rescue-two.fq",two)});
    REQUIRE(v.size()==90);
    CHECK((v[89].flag&4)==0);CHECK((v[89].flag&2)!=0);
    CHECK(v[89].cigar=="50M");CHECK(v[89].pos==2433);
    CHECK(tag_value(v[89],"XT:A")=="M");
    CHECK(tag_value(v[89],"NM:i")=="2"); // one planted substitution is chemistry-free
    CHECK(tag_value(v[89],"XR:Z")=="GA");CHECK(tag_value(v[89],"XG:Z")=="CT");
    CHECK(tag_value(v[89],"MD:Z")==expected_md(ref().c1,2432,aln_test::revcomp(two.back())));
    CHECK(tag_value(v[89],"XM:Z")==expected_xm(ref().c1,2432,aln_test::revcomp(two.back()),true,false));
    const auto cpg=run_aln({"--meth=taps","--meth-context","cpg","-n","0","-l","1024"},
                          {write_fq("rescue-cpg-one.fq",one),write_fq("rescue-cpg-two.fq",two)});
    REQUIRE(cpg.size()==90);CHECK(tag_value(cpg[89],"XT:A")=="M");
    CHECK(cpg[89].pos==2433);CHECK(cpg[89].cigar=="50M");
    const std::string aligned=aln_test::revcomp(two.back());
    int mismatches=0;
    for(int k=0;k<50;++k) if(aligned[k]!=ref().c1[2432+k]) {
        const bool free=ref().c1[2432+k]=='C'&&aligned[k]=='T'&&ref().c1[2433+k]=='G';
        if(!free) ++mismatches;
    }
    CHECK(tag_value(cpg[89],"NM:i")==std::to_string(mismatches));
    CHECK(tag_value(cpg[89],"MD:Z")==expected_md(ref().c1,2432,aligned));
    CHECK(tag_value(cpg[89],"XM:Z")==expected_xm(ref().c1,2432,aligned,true,true));
}

TEST_CASE("aln_meth: -v 0 silences the paired-end insert-size and rescue reports, not the work")
{
    std::vector<std::string> one,two;
    for(int i=0;i<45;++i) {
        const int start=100+50*i,span=180+i%21;
        one.push_back(ref().c1.substr(start,50));
        two.push_back(aln_test::revcomp(ref().c1.substr(start+span-50,50)));
    }
    for(int at : {15,25,35}) two.back()[at]=two.back()[at]=='A'?'C':'A'; // only rescue places it
    const std::string f1=write_fq("quiet-one.fq",one),f2=write_fq("quiet-two.fq",two);
    // both paired drivers: --meth and plain (stock sampe) aln
    for(const bool meth : {true,false}) {
        CAPTURE(meth);
        const int verbose=bwa_verbose;
        std::vector<Rec> rec[2];std::string err[2];
        for(int quiet=0;quiet<2;++quiet) {
            bwa_tests::StderrCapture capture;
            std::vector<std::string> args={"-n","0","-l","1024","-v",quiet?"0":"3"};
            if(meth) args.insert(args.begin(),"--meth");
            rec[quiet]=run_aln(args,{f1,f2});
            err[quiet]=capture.text();
        }
        bwa_verbose=verbose;
        CHECK(err[0].find("[infer_isize] inferred external isize")!=std::string::npos);
        CHECK(err[0].find("[bwa_paired_sw]")!=std::string::npos);
        CHECK(err[1].find("[infer_isize]")==std::string::npos);
        CHECK(err[1].find("[bwa_paired_sw]")==std::string::npos);
        CHECK(err[1].find("[M::aln]")==std::string::npos);
        // inference and rescue still ran: the rescued mate is placed and the records match
        REQUIRE(rec[1].size()==90);REQUIRE(rec[0].size()==90);
        CHECK(tag_value(rec[1][89],"XT:A")=="M");CHECK((rec[1][89].flag&2)!=0);
        for(size_t i=0;i<rec[0].size();++i) {
            CAPTURE(i);
            CHECK(rec[1][i].flag==rec[0][i].flag);CHECK(rec[1][i].pos==rec[0][i].pos);
            CHECK(rec[1][i].cigar==rec[0][i].cigar);CHECK(rec[1][i].tlen==rec[0][i].tlen);
            CHECK(rec[1][i].raw_tags==rec[0][i].raw_tags);
        }
    }
}

TEST_CASE("aln_meth: insert inference counts only same-contig pairs, as --meth pairs them")
{
    // 30 same-contig pairs on c1 near 400 bases, and 30 pairs whose ends straddle the c1/c2
    // boundary with a narrow ~200-base coordinate distance. Stock counts both; --meth, which
    // never pairs across contigs, must estimate from the same-contig pairs alone.
    namespace ba=bwa_aln;
    FMI_search fmi((ref().fa+".meth").c_str());fmi.load_index(false,1);
    ba::aln_meth_ref_t mr;REQUIRE(ba::aln_meth_ref_load(&mr,ref().fa.c_str(),fmi.idx->bns,1)==0);
    const int64_t c2=mr.bns->anns[1].offset;
    REQUIRE(c2==(int64_t)ref().c1.size());
    std::vector<ba::bwa_seq_t> ends[2];ends[0].resize(60);ends[1].resize(60);
    for(int i=0;i<60;++i) for(int j=0;j<2;++j) {
        auto &p=ends[j][i];p.mapQ=37;p.len=p.full_len=50;
        if(i<30) p.pos=j?c2+50+i%7:c2-150-i%5; // straddles: c1 tail, c2 head
        else p.pos=j?100+(i-30)*60+350+(i%9):100+(i-30)*60;
    }
    ba::bwa_seq_t *all[2]={ends[0].data(),ends[1].data()};
    ba::bwa_seq_t *same[2]={ends[0].data()+30,ends[1].data()+30};
    ba::isize_info_t want,stock,meth;
    REQUIRE(ba::infer_isize(30,same,&want,1e-5,mr.bns->l_pac,nullptr,0)==0);
    REQUIRE(ba::infer_isize(60,all,&stock,1e-5,mr.bns->l_pac,nullptr,0)==0);
    REQUIRE(ba::infer_isize(60,all,&meth,1e-5,mr.bns->l_pac,mr.bns,0)==0);
    CHECK(stock.avg<want.avg-50); // not vacuous: the straddling pairs pull stock's estimate down
    CHECK(meth.avg==want.avg);CHECK(meth.std==want.std);
    CHECK(meth.low==want.low);CHECK(meth.high==want.high);CHECK(meth.high_bayesian==want.high_bayesian);
    ba::aln_meth_ref_free(&mr);
}

TEST_CASE("aln_meth: stock insert inference rejects mixed short overlaps via global max read length")
{
    namespace ba=bwa_aln;
    std::vector<ba::bwa_seq_t> ends[2];ends[0].resize(41);ends[1].resize(41);
    for(int i=0;i<41;++i) for(int j=0;j<2;++j) {
        auto &p=ends[j][i];p.mapQ=37;p.pos=100+i*200+(j?5:0);
        p.len=p.full_len=(i==40?150:20+i%16);
    }
    ba::bwa_seq_t *seqs[2]={ends[0].data(),ends[1].data()};
    ba::isize_info_t ii;
    CHECK(ba::infer_isize(41,seqs,&ii,1e-5,100000)==-1);
    CHECK(ii.avg<0);CHECK(ii.low==150);CHECK(ii.high<150);
    // Removing the unrelated long-read pair permits the same short overlaps.
    CHECK(ba::infer_isize(40,seqs,&ii,1e-5,100000)==0);
    CHECK(ii.avg>20);CHECK(ii.avg<40);
}

TEST_CASE("aln_meth: rescue windows belong to the anchor contig at both boundaries")
{
    namespace ba=bwa_aln;
    FMI_search fmi((ref().fa+".meth").c_str());fmi.load_index(false,1);
    ba::aln_meth_ref_t mr;REQUIRE(ba::aln_meth_ref_load(&mr,ref().fa.c_str(),fmi.idx->bns,1)==0);
    ba::aln_meth_opt_t mo;ba::aln_meth_opt_finish(&mo);
    auto *opt=ba::gap_init_opt();auto *pe=ba::bwa_init_pe_opt();ba::bwase_initialize();
    struct Case {int rid,strand,anchor,target_rid,target;bool rescued;};
    const Case cases[]={{0,0,2920,1,70,false},{1,1,110,1,10,true},
                        {1,0,50,1,180,true},{0,1,2900,0,2800,true},
                        {0,0,2950,1,70,false},{1,1,0,0,2900,false},
                        {0,1,150,0,0,true}}; // the unmapped mate's pos is 0 too
    for(const auto &c:cases) {
        CAPTURE(c.rid);CAPTURE(c.strand);CAPTURE(c.anchor);
        const std::string &genome=c.target_rid?ref().c2:ref().c1;
        std::string mate=genome.substr(c.target,50);
        if(!c.strand) mate=aln_test::revcomp(mate);
        ba::bwa_seq_t *seqs[2];ba::bwa_seqio_t *inputs[2];int n[2];
        for(int j=0;j<2;++j) {
            inputs[j]=ba::bwa_seq_open(write_fq(j?"boundary2.fq":"boundary1.fq",{j?mate:std::string(50,'A')}).c_str());
            seqs[j]=ba::bwa_read_seq(inputs[j],1,&n[j],BWA_MODE_COMPREAD,0);REQUIRE(n[j]==1);
        }
        ba::aln_meth_read_t result[2];ba::aln_meth_read_t *r[2]={result,result+1};
        ba::aln_meth_occ_t o;o.rid=c.rid;o.strand=c.strand;o.hyp=1-c.strand;
        o.pos=mr.bns->anns[c.rid].offset+c.anchor;o.alignment.rb=o.pos;o.alignment.re=o.pos+50;
        o.alignment.cigar={50u<<4};result[0].occ.push_back(o);result[0].main=0;result[0].c1=1;
        ba::aln_meth_pe_stage1(opt,1,seqs,r);
        ba::isize_info_t ii={};ii.avg=190;ii.std=10;ii.ap_prior=1e-5;ii.low=100;ii.high=250;ii.high_bayesian=250;
        uint64_t total[2]={},mapped[2]={};ba::aln_meth_pe_stage2(&mr,&mo,opt,pe,&ii,1,seqs,r,total,mapped);
        if(c.rescued) {
            CHECK((int)seqs[1]->type==BWA_TYPE_MATESW);
            CHECK(seqs[1]->pos==(bwtint_t)(mr.bns->anns[c.rid].offset+c.target));
            CHECK((seqs[1]->extra_flag&SAM_FPP)!=0);
            CHECK(result[1].hyp==o.hyp);
            if(seqs[1]->md) CHECK(std::string(seqs[1]->md)=="50");
            else CHECK(seqs[1]->md!=nullptr);
        } else {
            CHECK((int)seqs[1]->type==BWA_TYPE_NO_MATCH);
            CHECK((seqs[0]->extra_flag&SAM_FPP)==0);CHECK((seqs[1]->extra_flag&SAM_FPP)==0);
        }
        for(int j=0;j<2;++j){ba::bwa_free_read_seq(1,seqs[j]);ba::bwa_seq_close(inputs[j]);}
    }
    free(pe);free(opt);ba::aln_meth_ref_free(&mr);
}

TEST_CASE("aln_meth: eligible pairs aggregate ties across separate contigs")
{
    namespace ba=bwa_aln;
    FMI_search fmi((ref().fa+".meth").c_str());fmi.load_index(false,1);
    ba::aln_meth_ref_t mr;REQUIRE(ba::aln_meth_ref_load(&mr,ref().fa.c_str(),fmi.idx->bns,1)==0);
    ba::bwase_initialize();auto *pe=ba::bwa_init_pe_opt();
    for(int second_cost:{3,6,9}) {
        ba::bwa_seq_t reads[2]={};ba::bwa_seq_t *p[2]={reads,reads+1};
        ba::bwt_aln1_t hits[2][2]={};
        for(int j=0;j<2;++j){p[j]->len=p[j]->full_len=50;p[j]->type=BWA_TYPE_UNIQUE;
            p[j]->mapQ=p[j]->seQ=37;p[j]->aln=hits[j];}
        p[0]->pos=100;p[0]->strand=0;p[1]->pos=3100;p[1]->strand=1;
        hits[1][0].score=3;hits[1][1].score=second_cost;
        pair64_t rows[4]={{100,0},{200,3},{3100,4},{3200,7}};
        pair64_v v;v.a=rows;v.n=v.m=4;
        ba::isize_info_t ii={};ii.avg=ii.std=-1;
        ba::bwa_pe_pairing_same_contig(mr.bns,p,&v,pe,3,&ii);
        CHECK(p[0]->pos==100);CHECK(p[1]->pos==200);
        CHECK((p[0]->extra_flag&SAM_FPP)!=0);
        // One competing pair on another contig: tie0, score-gap30/2=15, poor competitor23.
        CHECK((int)p[1]->mapQ==(second_cost==3?0:second_cost==6?15:23));
    }
    free(pe);ba::aln_meth_ref_free(&mr);
}

TEST_CASE("aln_meth: valid same-contig pair wins over cheaper cross-contig placement")
{
    namespace ba=bwa_aln;
    FMI_search fmi((ref().fa+".meth").c_str());fmi.load_index(false,1);
    ba::aln_meth_ref_t mr;REQUIRE(ba::aln_meth_ref_load(&mr,ref().fa.c_str(),fmi.idx->bns,1)==0);
    ba::aln_meth_opt_t mo;ba::aln_meth_opt_finish(&mo);
    auto *opt=ba::gap_init_opt();auto *pe=ba::bwa_init_pe_opt();pe->is_sw=0;ba::bwase_initialize();
    for(int max_insert:{200,500}) {
        pe->max_isize=max_insert;
        ba::bwa_seq_t *seqs[2];ba::bwa_seqio_t *inputs[2];int n[2];
        for(int j=0;j<2;++j){inputs[j]=ba::bwa_seq_open(write_fq(j?"valid2.fq":"valid1.fq",{std::string(50,'A')}).c_str());
            seqs[j]=ba::bwa_read_seq(inputs[j],1,&n[j],BWA_MODE_COMPREAD,0);REQUIRE(n[j]==1);}
        ba::aln_meth_read_t result[2];ba::aln_meth_read_t *r[2]={result,result+1};
        for(int j=0;j<2;++j){ba::aln_meth_occ_t o;o.rid=j;o.hyp=1;o.strand=j;o.pos=j?3100:2800;
            o.alignment.rb=o.pos;o.alignment.re=o.pos+50;o.alignment.cigar={50u<<4};
            result[j].occ.push_back(o);result[j].main=0;result[j].c1=1;}
        auto valid=result[1].occ[0];valid.rid=0;valid.pos=2900;valid.cost=3;valid.n_mm=1;
        valid.alignment.rb=2900;valid.alignment.re=2950;result[1].occ.push_back(valid);result[1].c2=1;
        ba::aln_meth_pe_stage1(opt,1,seqs,r);ba::isize_info_t ii={};ii.avg=ii.std=-1;
        uint64_t total[2]={},mapped[2]={};ba::aln_meth_pe_stage2(&mr,&mo,opt,pe,&ii,1,seqs,r,total,mapped);
        CHECK(result[1].main==1);CHECK(seqs[1]->pos==2900);CHECK(seqs[1]->score==3);
        CHECK((seqs[0]->extra_flag&SAM_FPP)!=0);CHECK((seqs[1]->extra_flag&SAM_FPP)!=0);
        CHECK((int)seqs[1]->mapQ==29);CHECK(ba::pos_5(seqs[1])-ba::pos_5(seqs[0])==150);
        for(int j=0;j<2;++j){ba::bwa_free_read_seq(1,seqs[j]);ba::bwa_seq_close(inputs[j]);}
    }
    free(pe);free(opt);ba::aln_meth_ref_free(&mr);
}

TEST_CASE("aln_meth: CpG per-base tags and costs follow original reference on both strands")
{
    const int start=500;const std::string window=ref().c1.substr(start,50);
    std::vector<std::string> reads;
    for(int strand=0;strand<2;++strand) {
        std::string aligned=window;bool inside=false,outside=false;
        for(int k=3;k<47;++k) if(window[k]==(strand?'G':'C')) {
            const bool cpg=strand?ref().c1[start+k-1]=='C':ref().c1[start+k+1]=='G';
            if((cpg&&!inside)||(!cpg&&!outside)) {
                aligned[k]=strand?'A':'T';if(cpg) inside=true;else outside=true;
            }
        }
        REQUIRE(inside);REQUIRE(outside);
        reads.push_back(strand?aln_test::revcomp(aligned):aligned);
    }
    const std::string fq=write_fq("cpg-both.fq",reads);
    for(const char *mode:{"genomic","neutral"}) {
        const auto hits=run_aln({"--meth=taps","--meth-context","cpg","--meth-scoring",mode,"-n","20","-l","1024"},{fq});
        REQUIRE(hits.size()==2);
        for(int strand=0;strand<2;++strand) {
            const std::string aligned=strand?aln_test::revcomp(reads[strand]):reads[strand];
            int expected_nm=0;
            for(int k=0;k<50;++k) if(aligned[k]!=window[k]) {
                const int pos=start+k;
                const bool free=strand ? window[k]=='G'&&ref().c1[pos-1]=='C'
                                       : window[k]=='C'&&ref().c1[pos+1]=='G';
                if(!free) ++expected_nm;
            }
            CAPTURE(strand);CAPTURE(mode);
            CHECK(hits[strand].rname=="c1");CHECK(hits[strand].pos==501);
            CHECK(hits[strand].cigar=="50M");
            CHECK(tag_value(hits[strand],"NM:i")==std::to_string(expected_nm));
            CHECK(tag_value(hits[strand],"MD:Z")==expected_md(ref().c1,start,aligned));
            CHECK(tag_value(hits[strand],"XM:Z")==expected_xm(ref().c1,start,aligned,!strand,true));
            CHECK(tag_value(hits[strand],"XG:Z")== (strand?"GA":"CT"));
        }
    }
}

TEST_CASE("aln_meth: CpG indel shifts query tags without shifting reference context")
{
    const int start=700;
    std::string query=ref().c1.substr(start,60);
    for(int k:{12,42}) {
        int at=k;
        while(at<k+12&&!(ref().c1[start+at]=='C'&&ref().c1[start+at+1]=='G')) ++at;
        REQUIRE(at<k+12);query[at]='T';
    }
    query.insert(30,"AAA");
    const auto hits=run_aln({"--meth=taps","--meth-context","cpg","-n","5","-l","1024","-o","1"},
                            {write_fq("cpg-insertion.fq",{query})});
    REQUIRE(hits.size()==1);CHECK(hits[0].rname=="c1");CHECK(hits[0].pos==start+1);
    CHECK(hits[0].cigar.find("3I")!=std::string::npos);
    // Walk the emitted alignment independently; reference context always uses original coordinates.
    int q=0,t=hits[0].pos-1,run=0,nm=0,matches=0;
    std::string xm(query.size(),'.'),md;
    for(char op:hits[0].cigar) {
        if(op>='0'&&op<='9'){run=10*run+op-'0';continue;}
        if(op=='M') for(int k=0;k<run;++k) {
            const char base=ref().c1[t],read=query[q];
            if(base!=read) {
                md+=std::to_string(matches)+base;matches=0;
                if(!(base=='C'&&read=='T'&&ref().c1[t+1]=='G')) ++nm;
            } else ++matches;
            if(base=='C'&&(read=='C'||read=='T')) {
                const char context=ref().c1[t+1]=='G'?'z':ref().c1[t+2]=='G'?'x':'h';
                xm[q]=read=='T'?context-'a'+'A':context;
            }
            ++q;++t;
        }
        else if(op=='I'){q+=run;nm+=run;}
        else if(op=='S')q+=run;
        else {FAIL("unexpected deletion in planted insertion fixture");}
        run=0;
    }
    md+=std::to_string(matches);
    CHECK(q==(int)query.size());CHECK(tag_value(hits[0],"NM:i")==std::to_string(nm));
    CHECK(nm==3);CHECK(tag_value(hits[0],"MD:Z")==md);CHECK(tag_value(hits[0],"XM:Z")==xm);
}

TEST_CASE("aln_meth: CpG edge and original ambiguity tags never borrow neighbor contigs")
{
    const std::string fa=ref().dir+"/context-boundaries.fa";
    uint32_t seed=514;
    const std::string left=aln_test::lcg_dna(seed,59)+"C";
    const std::string right="G"+aln_test::lcg_dna(seed,59);
    std::string ambiguous=aln_test::lcg_dna(seed,80);ambiguous[29]='C';ambiguous[30]='N';ambiguous[31]='G';
    REQUIRE(aln_test::write_test_fasta(fa,{{"left",left},{"right",right},{"ambiguous",ambiguous}}));
    std::vector<std::string> cmd={"index","--meth",fa};std::vector<char*> argv;
    for(auto &x:cmd)argv.push_back(&x[0]);argv.push_back(nullptr);
    optind=1;
#if defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || defined(__NetBSD__)
    optreset=1;
#endif
    REQUIRE(bwa_index((int)cmd.size(),argv.data())==0);
    std::string l=left;l.back()='T';
    std::string r=right;r[0]='A';r=aln_test::revcomp(r);
    const auto edge=run_aln({"--meth=taps","--meth-context","cpg","-n","2","-l","1024"},
                            {write_fq("context-boundaries.fq",{l,r})},fa);
    REQUIRE(edge.size()==2);
    CHECK(edge[0].rname=="left");CHECK(edge[1].rname=="right");
    CHECK(tag_value(edge[0],"NM:i")=="1");CHECK(tag_value(edge[1],"NM:i")=="1");
    CHECK(tag_value(edge[0],"MD:Z")=="59C0");CHECK(tag_value(edge[1],"MD:Z")=="0G59");
    CHECK(tag_value(edge[0],"XM:Z").back()=='U');
    CHECK(tag_value(edge[1],"XM:Z").front()=='U');
    const auto n=run_aln({"--meth=taps","--meth-context","cpg","-n","4","-l","1024"},
                         {write_fq("context-ambiguity.fq",{ambiguous})},fa);
    REQUIRE(n.size()==1);CHECK(n[0].rname=="ambiguous");
    CHECK(tag_value(n[0],"NM:i")=="1");
    const auto &xm=tag_value(n[0],"XM:Z");REQUIRE(xm.size()==80);
    CHECK(xm[29]=='u');CHECK(xm[30]=='.');
    CHECK(xm==expected_xm(ambiguous,0,ambiguous,true,true));
    // Stock MD describes the packed original reference, whose N replacement is index-seeded.
    FMI_search fmi((fa+".meth").c_str());fmi.load_index(false,1);
    bwa_aln::aln_meth_ref_t original;
    REQUIRE(bwa_aln::aln_meth_ref_load(&original,fa.c_str(),fmi.idx->bns,1)==0);
    const int64_t pos=original.bns->anns[2].offset+30;
    std::string packed=ambiguous;
    packed[30]="ACGT"[(original.pac[pos>>2]>>((~pos&3)*2))&3];
    CHECK(tag_value(n[0],"MD:Z")==expected_md(packed,0,ambiguous));
    bwa_aln::aln_meth_ref_free(&original);
}

TEST_CASE("aln_meth: a refined window whose counts overflow a record is dropped, not fatal")
{
    // Read 1 C->T converts a poly-C read to poly-T, which hits a homopolymer T run in the converted
    // index with few differences. Refined against the original reference (where read C over
    // reference T is not free for TAPS/5-base), the best global alignment needs more than 255 gap
    // extensions: it cannot be stored in a record, so it must not compete, and must not abort.
    const std::string fa = ref().dir + "/homopolymer.fa";
    uint32_t st = 424242;
    const std::string left = aln_test::lcg_dna(st, 2000), right = aln_test::lcg_dna(st, 2000);
    REQUIRE(aln_test::write_test_fasta(fa, {{"poly", left + std::string(400, 'T') + right}}));
    std::vector<std::string> args = {"index", "--meth", fa};
    std::vector<char *> argv;
    for (auto &s : args) argv.push_back(&s[0]);
    argv.push_back(nullptr);
    optind = 1;
#if defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || defined(__NetBSD__)
    optreset = 1;
#endif
    REQUIRE(bwa_index((int)args.size(), argv.data()) == 0);
    const std::string polyc(150, 'C');
    const std::string adapter = "GAAGAGCACACGTCTGAACTCCAGTCACATGTAAGTATCTAGTGTGTCGTCTTTTTTTTT";
    for (const char *mode : {"--meth=taps", "--meth=5base"}) {
        CAPTURE(mode);
        // Single-end: no representable candidate, so the read is unmapped.
        auto se = run_aln({mode}, {write_fq("polyc.fq", {polyc})}, fa);
        REQUIRE(se.size() == 1);
        CHECK((se[0].flag & 4) != 0);
        // Paired: both mates are reported, neither aborts the run. One pair cannot infer an
        // insert size, so mate rescue does not run here (see the next case for rescue).
        auto pe = run_aln({mode}, {write_fq("polyc_1.fq", {adapter}), write_fq("polyc_2.fq", {std::string(150, 'G')})}, fa);
        REQUIRE(pe.size() == 2);
    }
    // EM-seq frees read C over reference T after conversion, so the same read still places.
    auto em = run_aln({"--meth=emseq", "--meth-scoring", "collapsed"}, {write_fq("polyc_em.fq", {polyc})}, fa);
    REQUIRE(em.size() == 1);
    CHECK((em[0].flag & 4) == 0);
}

TEST_CASE("aln_meth: mate rescue drops a candidate whose counts overflow a record, not fatal")
{
    // Mate rescue with a valid insert size and a unique anchor: the mate is 800 bases of c1 with
    // a read insertion. A 200-base insertion is rescued (the rescue path runs); a 260-base one
    // needs 259 gap extensions, more than a record holds, so the mate stays unmapped.
    namespace ba=bwa_aln;
    FMI_search fmi((ref().fa+".meth").c_str());fmi.load_index(false,1);
    ba::aln_meth_ref_t mr;REQUIRE(ba::aln_meth_ref_load(&mr,ref().fa.c_str(),fmi.idx->bns,1)==0);
    ba::aln_meth_opt_t mo;ba::aln_meth_opt_finish(&mo);
    auto *opt=ba::gap_init_opt();auto *pe=ba::bwa_init_pe_opt();ba::bwase_initialize();
    const int anchor=100,start=600,span=800;
    uint32_t st=97531;
    for(int insertion : {200,260}) {
        CAPTURE(insertion);
        const std::string seg=ref().c1.substr(start,span);
        const std::string mate=aln_test::revcomp(seg.substr(0,400)+aln_test::lcg_dna(st,insertion)+seg.substr(400));
        ba::bwa_seq_t *seqs[2];ba::bwa_seqio_t *inputs[2];int n[2];
        for(int j=0;j<2;++j) {
            inputs[j]=ba::bwa_seq_open(write_fq(j?"overflow-rescue2.fq":"overflow-rescue1.fq",{j?mate:std::string(50,'A')}).c_str());
            seqs[j]=ba::bwa_read_seq(inputs[j],1,&n[j],BWA_MODE_COMPREAD,0);REQUIRE(n[j]==1);
        }
        ba::aln_meth_read_t result[2];ba::aln_meth_read_t *r[2]={result,result+1};
        ba::aln_meth_occ_t o;o.rid=0;o.strand=0;o.hyp=1;
        o.pos=mr.bns->anns[0].offset+anchor;o.alignment.rb=o.pos;o.alignment.re=o.pos+50;
        o.alignment.cigar={50u<<4};result[0].occ.push_back(o);result[0].main=0;result[0].c1=1;
        ba::aln_meth_pe_stage1(opt,1,seqs,r);
        REQUIRE((int)seqs[0]->mapQ>=17);REQUIRE((int)seqs[1]->type==BWA_TYPE_NO_MATCH);
        ba::isize_info_t ii={};ii.avg=start+span-anchor;ii.std=10;ii.ap_prior=1e-5;
        ii.low=100;ii.high=2000;ii.high_bayesian=2000;
        uint64_t total[2]={},mapped[2]={};ba::aln_meth_pe_stage2(&mr,&mo,opt,pe,&ii,1,seqs,r,total,mapped);
        CHECK(total[1]==1); // the singleton was offered to mate rescue
        if(insertion==200) {
            CHECK(mapped[1]==1);
            REQUIRE((int)seqs[1]->type==BWA_TYPE_MATESW);
            CHECK(seqs[1]->pos==(bwtint_t)(mr.bns->anns[0].offset+start));
            CHECK((int)seqs[1]->strand==1);
            CHECK((int)seqs[1]->n_gapo==1);CHECK((int)seqs[1]->n_gape==insertion-1);
            REQUIRE(seqs[1]->n_cigar==3);
            CHECK(__cigar_op(seqs[1]->cigar[0])==FROM_M);CHECK(__cigar_op(seqs[1]->cigar[2])==FROM_M);
            CHECK(__cigar_op(seqs[1]->cigar[1])==FROM_I);CHECK(__cigar_len(seqs[1]->cigar[1])==insertion);
            CHECK(__cigar_len(seqs[1]->cigar[0])+__cigar_len(seqs[1]->cigar[2])==span);
            CHECK((seqs[1]->extra_flag&SAM_FPP)!=0);
        } else {
            CHECK(mapped[1]==0);
            CHECK((int)seqs[1]->type==BWA_TYPE_NO_MATCH);CHECK(seqs[1]->cigar==nullptr);
            CHECK((seqs[0]->extra_flag&SAM_FPP)==0);CHECK((seqs[1]->extra_flag&SAM_FPP)==0);
        }
        for(int j=0;j<2;++j){ba::bwa_free_read_seq(1,seqs[j]);ba::bwa_seq_close(inputs[j]);}
    }
    free(pe);free(opt);ba::aln_meth_ref_free(&mr);
}
