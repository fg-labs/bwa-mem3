// Unit tests for `aln --meth`: orientation, chemistry and tags, stock cost ranking,
// bounded candidate discovery, conservative caps, duplicate/thread invariance,
// and stock alignment against a --meth-built index.

#include "doctest/doctest.h"

#include "aln_test_ref.h"
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
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <unistd.h>
#include <sys/wait.h>
#include <limits>

namespace {

std::string conv(std::string s, char from, char to)
{
    for (char &c : s)
        if (c == from) c = to;
    return s;
}

// The reference: c1 and c2 random; c3 holds the -R counterexample's loci (see below).
struct MethRef {
    std::string dir, fa, c1, c2, read_r, c3, affine_read;
    std::string c4, read_u, read_v, span_read, span_reference;
    std::string c5,c6,c7,c8,seg_p,frag_m,seg_w,frag_n;
    std::string c9,seg_q,frag_r,c10,cap_read;int q_pos=0,r_pos=0;
    int p_pos=0,m_pos=0,w_pos=0,n_pos=0;
    int u_pos = 0, v_pos[2] = {0, 0};
    int x_pos = 0; // 1-based position of the counterexample's true locus in c3

    MethRef()
    {
        std::string tmpl_s = aln_test::tmp_template("aln_meth_unit");
        char *tmpl = &tmpl_s[0];
        REQUIRE(mkdtemp(tmpl) != nullptr);
        dir = tmpl;
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
        auto with_ts = [&](int n_t) { // a random 40-mer with at least n_t T in [5, 35)
            for (;;) {
                std::string w = aln_test::lcg_dna(st, 40);
                int nt = 0;
                for (int i = 5; i < 35; ++i) nt += w[i] == 'T';
                if (nt >= n_t) return w;
            }
        };
        auto t_to_c = [](std::string w, int nth, int n) { // the nth.. (nth+n-1)th T of [5, 35) made C
            for (int i = 5, k = 0; i < 35 && n > 0; ++i)
                if (w[i] == 'T' && k++ >= nth) w[i] = 'C', --n;
            return w;
        };
        read_u = with_ts(4);
        read_v = with_ts(2);
        c4 = aln_test::lcg_dna(st, 100);
        u_pos = (int)c4.size() + 1;
        c4 += read_u + aln_test::lcg_dna(st, 80) + t_to_c(read_u, 0, 1) + aln_test::lcg_dna(st, 80)
            + t_to_c(read_u, 1, 1) + aln_test::lcg_dna(st, 80) + t_to_c(read_u, 2, 2) + aln_test::lcg_dna(st, 80);
        v_pos[0] = (int)c4.size() + 1;
        c4 += read_v + aln_test::lcg_dna(st, 80);
        v_pos[1] = (int)c4.size() + 1;
        c4 += read_v + aln_test::lcg_dna(st, 80) + t_to_c(read_v, 0, 1) + aln_test::lcg_dna(st, 100);

        span_read=aln_test::lcg_dna(st,20)+"T"+std::string(19,'A');
        span_reference=span_read.substr(0,20)+"CT"+std::string(19,'A');
        seg_p = aln_test::lcg_dna(st, 36);
        frag_m = aln_test::lcg_dna(st, 36);
        seg_w = aln_test::lcg_dna(st, 36);
        frag_n = aln_test::lcg_dna(st, 36);
        c5 = aln_test::lcg_dna(st, 50);
        for (int k = 0; k < 600; ++k) c5 += seg_p + aln_test::lcg_dna(st, 64);
        c8 = aln_test::lcg_dna(st, 50);
        for (int k = 0; k < 600; ++k) {
            std::string w = seg_w;
            if (k >= 300) w[10 + k % 16] = w[10 + k % 16] == 'A' ? 'G' : 'A'; // one substitution
            c8 += w + aln_test::lcg_dna(st, 64);
        }
        c6 = aln_test::lcg_dna(st, 100);
        p_pos = (int)c6.size() + 1;
        c6 += seg_p + aln_test::lcg_dna(st, 150);
        m_pos = (int)c6.size() + 1;
        c6 += frag_m + aln_test::lcg_dna(st, 2000);
        w_pos = (int)c6.size() + 1;
        c6 += seg_w + aln_test::lcg_dna(st, 150);
        n_pos = (int)c6.size() + 1;
        c6 += frag_n + aln_test::lcg_dna(st, 100);
        c7 = aln_test::lcg_dna(st, 100) + frag_m + aln_test::lcg_dna(st, 2000) + frag_n + aln_test::lcg_dna(st, 100);
        seg_q = aln_test::lcg_dna(st, 36);
        frag_r = aln_test::lcg_dna(st, 36);
        c6 += aln_test::lcg_dna(st,2000);
        q_pos = (int)c6.size() + 1;
        c6 += seg_q + aln_test::lcg_dna(st, 150);
        r_pos = (int)c6.size() + 1;
        c6 += frag_r + aln_test::lcg_dna(st, 100);
        c7 += aln_test::lcg_dna(st,2000)+frag_r+aln_test::lcg_dna(st,100);
        c9 = aln_test::lcg_dna(st, 50);
        for (int k = 0; k < 4800; ++k) c9 += seg_q + aln_test::lcg_dna(st, 30);
        cap_read=aln_test::lcg_dna(st,50);cap_read[20]='T';cap_read[21]='G';
        std::string converted_cap=cap_read;converted_cap[20]='C';
        c10=aln_test::lcg_dna(st,250)+converted_cap+aln_test::lcg_dna(st,350)+cap_read+aln_test::lcg_dna(st,100);
        REQUIRE(aln_test::write_test_fasta(fa, {{"c1", c1}, {"c2", c2}, {"c3", c3},
                                              {"c4", c4}, {"c5",c5},{"c6",c6},{"c7",c7},{"c8",c8},{"c9",c9},{"cap_choice",c10}, {"span",span_reference}, {"affine_expensive", expensive}, {"affine_cheap", cheap}}));
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
    ~MethRef()
    {
        std::string cmd = "rm -rf '" + dir + "'";
        if (system(cmd.c_str()) != 0) std::fprintf(stderr, "warning: could not remove %s\n", dir.c_str());
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
                CHECK(a[i].pos == b[i].pos);
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
    std::vector<std::string> reads;
    for (int i = 0; i < 40; ++i) reads.push_back(conv(ref().c1.substr(37 * i, 36), 'C', 'T'));
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
                       const std::string &mate = "")
{
    namespace ba = bwa_aln;
    ba::aln_opt_t o;
    ba::aln_opt_init(&o);
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

TEST_CASE("aln_meth: identical reads share a search and still draw their own main hit")
{
    // a batch of mostly duplicated reads, so the duplicate pass runs, with repeats (several best
    // hits) so the draw matters
    std::vector<std::string> reads;
    for (int k = 0; k < 600; ++k) {
        reads.push_back(ref().read_r);                                   // a read with 31 repeat loci
        reads.push_back(ref().read_v); // the neutral conversion split must replay per read
        reads.push_back(conv(ref().c1.substr(37 * (k % 40), 36), 'C', 'T'));
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
    std::vector<std::string> reads;
    for (int i = 0; i < 90; ++i) {
        reads.push_back(conv(ref().c1.substr(500, 50), 'C', 'T'));
        reads.push_back(std::string(50, 'N'));
        reads.push_back(conv(aln_test::revcomp(ref().c2.substr(600, 50)), 'C', 'T'));
    }
    const auto records = run_aln({"--meth", "-t", "3"}, {write_fq("slot-reuse.fq", reads)});
    REQUIRE(records.size() == reads.size());
    for (size_t i = 0; i < records.size(); ++i) {
        CAPTURE(i);
        if (i % 3 == 1) {
            CHECK((records[i].flag & 4) != 0);
            CHECK(records[i].rname == "*");
            CHECK(tag_value(records[i], "XR:Z") == "CT");
            CHECK(records[i].tag.count("XM:Z") == 0);
            CHECK(records[i].tag.count("XG:Z") == 0);
        } else CHECK((records[i].flag & 4) == 0);
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
    for(int duplicate=0;duplicate<601;++duplicate) ba::aln_meth_hit(&slot,&exact,row,row);
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
    // Metadata, not the independently random-filled original/converted PACs, owns ambiguity.
    const int old_holes=mr.bns->n_holes;
    bntamb1_t *old_ambs=mr.bns->ambs;
    bntamb1_t ambiguous={}; ambiguous.offset=mr.bns->anns[0].offset+520; ambiguous.len=1;
    mr.bns->n_holes=1; mr.bns->ambs=&ambiguous;
    for (int context : {METH_CTX_ALL, METH_CTX_CPG}) {
        mo.m.meth_context=context; ba::aln_meth_opt_finish(&mo);
        ba::aln_meth_slot_start(&slot,reads,&result,0,3,1);
        ba::aln_meth_hit(&slot,&exact,row,row); core.state=ba::GAP_SLOT_PARKED;
        ba::aln_meth_flush_parked(&group);
        REQUIRE(result.occ.size()==1);
        int expected_mm=1;
        if (context==METH_CTX_CPG)
            for (int k=0;k<50;++k)
                if (k!=20 && ref().c1[500+k]=='C' &&
                    (ref().c1[501+k]!='G' || k==19)) ++expected_mm;
        CHECK(result.occ[0].n_mm==expected_mm);
        CHECK(result.occ[0].cost==expected_mm*opt->s_mm);
    }
    mo.m.meth_context=METH_CTX_ALL; ba::aln_meth_opt_finish(&mo);
    mr.bns->n_holes=old_holes; mr.bns->ambs=old_ambs;
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
    mr.bns->n_holes=1; mr.bns->ambs=&ambiguous;
    ba::aln_meth_finish_se(&mr,&mo,opt,3,1,reads,&result);
    CHECK((int)reads[0].nm==1);
    CHECK(std::string(reads[0].md).find('N')!=std::string::npos);
    CHECK((int)reads[0].type==BWA_TYPE_REPEAT);CHECK((int)reads[0].mapQ==0);
    CHECK((int)reads[0].c1==0);CHECK((int)reads[0].c2==0);CHECK(reads[0].n_multi==0);
    kstring_t text={0,0,nullptr};
    ba::aln_meth_print_sam1(&text,&mr,&mo,reads,&result,1,nullptr,BWA_MODE_COMPREAD,30);
    std::string record(text.s,text.l);free(text.s);
    CHECK(record.find("X0:i:")==std::string::npos);
    CHECK(record.find("X1:i:")==std::string::npos);
    CHECK(record.find("XA:Z:")==std::string::npos);
    mr.bns->n_holes=old_holes; mr.bns->ambs=old_ambs;
    ba::aln_meth_ref_free(&mr);ba::bwt_destroy(bwt);free(opt);
    ba::bwa_free_read_seq(n,reads);CHECK(ba::bwa_seq_close(input)==0);
}

TEST_CASE("aln_meth: ineligible searches batch pending hits instead of flushing score boundaries"
          * doctest::test_suite("unit/aln_meth"))
{
    namespace ba=bwa_aln;
    ba::gap_opt_t *opt=ba::gap_init_opt();
    ba::gap_job_t job={}; job.opt=*opt;
    ba::aln_meth_opt_t mo;
    ba::aln_meth_group_t group; group.mo=&mo;
    ba::aln_meth_read_t result;
    ba::aln_meth_slot_t slot; slot.grp=&group; slot.job=&job; slot.res=&result;
    ba::bwt_width_t seed={};
    for (int mode=0;mode<7;++mode) {
        CAPTURE(mode);
        job.opt=*opt; job.seed_width=nullptr; result.flags=0;
        if (mode==0) job.opt.mode|=BWA_MODE_NONSTOP;
        if (mode==1) job.seed_width=&seed;
        if (mode==2) job.opt.s_mm=0;
        if (mode==3) job.opt.s_gape=-1;
        if (mode==4) result.flags=ba::ALN_METH_BIND_M;
        if (mode==5) result.flags=ba::ALN_METH_BIND_BEST;
        if (mode==6) job.opt.s_gapo=job.opt.s_mm;
        slot.flushed_score=3;
        slot.pend.clear();
        CHECK_FALSE(ba::aln_meth_needs_flush(&slot,6));
        slot.pend.resize(31);
        CHECK_FALSE(ba::aln_meth_needs_flush(&slot,6));
        slot.pend.resize(32);
        CHECK(ba::aln_meth_needs_flush(&slot,6));
        CHECK(ba::aln_meth_needs_flush(&slot,3));
    }
    // No scored candidates yet: boundary refinement must still be able to find the first best.
    job.opt=*opt; job.seed_width=nullptr; result.flags=0;
    slot.flushed_score=3; slot.pend.resize(1);
    CHECK(result.occ.empty());
    CHECK(ba::aln_meth_needs_flush(&slot,6));
    CHECK_FALSE(ba::aln_meth_needs_flush(&slot,3));
    // Reusing the same slot for an eligible read must not inherit a disabled gate.
    job.opt.mode|=BWA_MODE_NONSTOP;
    CHECK_FALSE(ba::aln_meth_needs_flush(&slot,6));
    job.opt.mode&=~BWA_MODE_NONSTOP; slot.flushed_score=3;
    CHECK(ba::aln_meth_needs_flush(&slot,6));
    // An already exhausted lookup budget must be detected at the original boundary,
    // before a later expansion can hit the work cap.
    job.opt.mode|=BWA_MODE_NONSTOP;
    slot.pend.resize(3);
    slot.pend[0].k=10; slot.pend[0].l=11;
    slot.pend[1].k=20; slot.pend[1].l=22;
    slot.pend[2].k=30; slot.pend[2].l=33;
    mo.p_best=9; mo.p_second=0;
    CHECK_FALSE(ba::aln_meth_needs_flush(&slot,6));
    result.n_lookups=1; slot.flushed_score=3;
    CHECK(ba::aln_meth_needs_flush(&slot,6));
    CHECK_FALSE(ba::aln_meth_needs_flush(&slot,3));
    result.n_lookups=0; mo.p_best=0; slot.flushed_score=3;
    CHECK(ba::aln_meth_needs_flush(&slot,6));
    free(opt);
}

TEST_CASE("aln_meth: suppressed boundaries do not move lookup-cap detection into a stratum"
          * doctest::test_suite("unit/aln_meth"))
{
    namespace ba=bwa_aln;
    auto *opt=ba::gap_init_opt();
    ba::gap_job_t job={}; job.opt=*opt; job.opt.mode|=BWA_MODE_NONSTOP;
    ba::aln_meth_opt_t mo; mo.p_best=1; mo.p_second=0;
    ba::aln_meth_group_t group; group.mo=&mo;
    ba::aln_meth_read_t result;
    ba::aln_meth_slot_t slot; slot.grp=&group; slot.job=&job; slot.res=&result;
    CHECK_FALSE(ba::aln_meth_needs_flush(&slot,0));
    slot.pend.resize(1); slot.pend[0].k=10; slot.pend[0].l=11;
    // The old empty score-zero flush would already have marked this score examined.
    CHECK_FALSE(ba::aln_meth_needs_flush(&slot,0));
    CHECK(ba::aln_meth_needs_flush(&slot,3));
    // The bounded-buffer override remains independent of stratum boundaries.
    slot.pend.resize(32);
    CHECK(ba::aln_meth_needs_flush(&slot,0));
    free(opt);
}

TEST_CASE("aln_meth: mixed seed eligibility and caps survive real slot reuse and final flushing"
          * doctest::test_suite("unit/aln_meth"))
{
    namespace ba=bwa_aln;
    FMI_search fmi((ref().fa+".meth").c_str()); fmi.load_index(false,1);
    auto *bwt=ba::bwt_from_fmi(&fmi);
    ba::aln_meth_ref_t mr;
    REQUIRE(ba::aln_meth_ref_load(&mr,ref().fa.c_str(),fmi.idx->bns,1)==0);
    ba::aln_meth_opt_t mo; ba::aln_meth_opt_finish(&mo);
    auto *opt=ba::gap_init_opt(); opt->fnr=-1; opt->max_diff=2; opt->seed_len=30;
    std::vector<std::string> queries;
    // Both slot-group sizes (8 and 12) revisit each slot on the opposite side of seed_len.
    for (int i=0;i<40;++i) {
        const int len=(i/4)%2 ? 50 : 20;
        queries.push_back(conv(ref().c1.substr(400+i*3,len),'C','T'));
    }
    const auto fq=write_fq("flush-slot-reuse.fq",queries);
    for (bool nonstop : {false,true}) for (int budget : {1,4596})
        for (int max_entries : {20,2000000}) {
        CAPTURE(nonstop); CAPTURE(budget); CAPTURE(max_entries);
        if (nonstop) opt->mode|=BWA_MODE_NONSTOP;
        else opt->mode&=~BWA_MODE_NONSTOP;
        opt->max_entries=max_entries; mo.p_best=budget; mo.p_second=0;
        auto *input=ba::bwa_seq_open(fq.c_str()); int n=0;
        auto *reads=ba::bwa_read_seq(input,100,&n,BWA_MODE_COMPREAD,0);
        REQUIRE(n==40);
        std::vector<ba::aln_meth_read_t> batch(n);
        ba::aln_dup_state_t dup;
        ba::aln_meth_search_batch(bwt,&mr,&mo,n,reads,opt,1,1,0,1,1000,&dup,batch.data());
        // Fresh single-read groups provide an independent lifecycle reference for reused slots.
        for (int i=0;i<n;++i) {
            ba::aln_meth_read_t single; ba::aln_dup_state_t single_dup;
            ba::aln_meth_search_batch(bwt,&mr,&mo,1,reads+i,opt,1,1,i,1,1000,&single_dup,&single);
            const auto &actual=batch[i];
            CHECK(actual.flags==single.flags);
            CHECK(actual.n_lookups==single.n_lookups);
            CHECK(actual.n_rescored==single.n_rescored);
            CHECK(actual.c1==single.c1); CHECK(actual.c2==single.c2);
            CHECK(actual.main==single.main);
            REQUIRE(actual.occ.size()==single.occ.size());
            for (size_t k=0;k<actual.occ.size();++k) {
                const auto &a=actual.occ[k], &b=single.occ[k];
                CHECK(std::make_tuple(a.rid,a.hyp,a.strand,a.pos,a.cost)==
                      std::make_tuple(b.rid,b.hyp,b.strand,b.pos,b.cost));
                CHECK(a.alignment.raw_cigar==b.alignment.raw_cigar);
                CHECK(a.alignment.cigar==b.alignment.cigar);
            }
            if (budget==4596 && max_entries==2000000) {
                CHECK(actual.flags==0);
                CHECK_FALSE(actual.occ.empty()); // Pending exact hits must reach final refinement.
            }
        }
        ba::bwa_free_read_seq(n,reads); ba::bwa_seq_close(input);
    }
    free(opt); ba::aln_meth_ref_free(&mr); ba::bwt_destroy(bwt);
}

TEST_CASE("aln_meth: retained-cost stop waits for complete strata and handles decreasing best costs")
{
    namespace ba=bwa_aln;
    ba::gap_opt_t *opt=ba::gap_init_opt();
    ba::aln_meth_group_t group; group.opt=opt;
    ba::gap_job_t job={}; job.opt=*opt; job.opt.max_diff=4;
    ba::aln_meth_read_t result;
    ba::aln_meth_occ_t exact; exact.cost=0; result.occ.push_back(exact);
    ba::aln_meth_slot_t slot; slot.grp=&group; slot.job=&job; slot.res=&result;
    CHECK_FALSE(ba::aln_meth_can_stop(&slot,opt->s_mm));
    CHECK(ba::aln_meth_can_stop(&slot,opt->s_mm+1));
    job.e.info=(uint32_t)(opt->s_mm+1)<<21;
    ba::gap_entry_t held={};
    CHECK(ba::aln_meth_resume(&slot,&job,&held)==1);
    CHECK(slot.flushed_score==opt->s_mm+1);
    result.occ[0].cost=3;
    CHECK_FALSE(ba::aln_meth_can_stop(&slot,6));
    CHECK(ba::aln_meth_can_stop(&slot,7));
    // A better result discovered while flushing lowers the necessary stratum boundary.
    result.occ[0].cost=6; CHECK_FALSE(ba::aln_meth_can_stop(&slot,7));
    CHECK(ba::aln_meth_can_stop(&slot,10));
    result.occ[0].cost=3; CHECK(ba::aln_meth_can_stop(&slot,7));
    job.opt.max_diff=1; CHECK_FALSE(ba::aln_meth_can_stop(&slot,7)); job.opt.max_diff=4;
    job.opt.s_gapo=6; CHECK_FALSE(ba::aln_meth_can_stop(&slot,7)); job.opt.s_gapo=11;
    result.occ[0].cost=INT_MAX; CHECK_FALSE(ba::aln_meth_can_stop(&slot,INT_MAX));
    job.opt.s_mm=INT_MAX; job.opt.s_gapo=INT_MAX;
    result.occ[0].cost=INT_MAX; CHECK_FALSE(ba::aln_meth_can_stop(&slot,INT_MAX));
    job.opt.s_mm=3;job.opt.s_gapo=11;
    // Completing score 3 must not suppress a later flush at 6/7 with a nonzero best.
    slot.flushed_score=3; CHECK(ba::aln_meth_needs_flush(&slot,6));
    job.e.info=6u<<21; result.occ[0].cost=3;
    CHECK(ba::aln_meth_resume(&slot,&job,&held)==0);
    CHECK_FALSE(ba::aln_meth_needs_flush(&slot,6));
    CHECK(ba::aln_meth_needs_flush(&slot,7));
    job.e.info=7u<<21; CHECK(ba::aln_meth_resume(&slot,&job,&held)==1);
    result.occ[0].cost=0;
    for (int flag : {ba::ALN_METH_BIND_BEST,ba::ALN_METH_BIND_M}) {
        result.flags=flag; CHECK_FALSE(ba::aln_meth_can_stop(&slot,opt->s_mm+1));
    }
    result.flags=0;
    ba::bwt_width_t seed={}; job.seed_width=&seed;
    CHECK_FALSE(ba::aln_meth_can_stop(&slot,opt->s_mm+1)); job.seed_width=nullptr;
    job.opt.mode|=BWA_MODE_NONSTOP; CHECK_FALSE(ba::aln_meth_can_stop(&slot,opt->s_mm+1));
    job.opt.mode&=~BWA_MODE_NONSTOP;
    job.opt.max_diff=0; CHECK_FALSE(ba::aln_meth_can_stop(&slot,opt->s_mm+1)); job.opt.max_diff=4;
    job.opt.s_gapo=job.opt.s_mm; CHECK_FALSE(ba::aln_meth_can_stop(&slot,opt->s_mm+1));
    job.opt.s_gapo=11; job.opt.s_mm=0; CHECK_FALSE(ba::aln_meth_can_stop(&slot,1));
    job.opt.s_mm=3; job.opt.s_gape=-1; CHECK_FALSE(ba::aln_meth_can_stop(&slot,4));
    free(opt);
}

TEST_CASE("aln_meth: zero-cost stopping preserves an independently enumerated tiny candidate set")
{
    namespace ba=bwa_aln;
    const std::string fa=ref().dir+"/zero-bound.fa";
    const std::string query="ACGTCGATCGTA";
    std::vector<std::pair<std::string,std::string>> contigs={{"exact",query}};
    for (int i=0;i<(int)query.size();++i) for(char base : std::string("ACGTN")) {
        std::string target=query; target[i]=base;
        contigs.emplace_back("v"+std::to_string(contigs.size()),target);
    }
    REQUIRE(aln_test::write_test_fasta(fa,contigs));
    std::vector<std::string> args={"index","--meth",fa};std::vector<char*> argv;
    for(auto &s:args)argv.push_back(&s[0]);argv.push_back(nullptr);optind=1;
#if defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || defined(__NetBSD__)
    optreset=1;
#endif
    REQUIRE(bwa_index((int)args.size(),argv.data())==0);
    FMI_search fmi((fa+".meth").c_str());fmi.load_index(false,1);
    auto *bwt=ba::bwt_from_fmi(&fmi);ba::aln_meth_ref_t mr;
    REQUIRE(ba::aln_meth_ref_load(&mr,fa.c_str(),fmi.idx->bns,1)==0);
    auto *opt=ba::gap_init_opt();opt->fnr=-1;opt->max_diff=2;opt->seed_len=1024;
    using Key=std::tuple<int,int,int,int64_t,int64_t>;
    for(int scoring : {MEM_METH_SCORING_COLLAPSED,MEM_METH_SCORING_GENOMIC,MEM_METH_SCORING_NEUTRAL})
        for(int context : {METH_CTX_ALL,METH_CTX_CPG}) for(int read_no : {0,1}) for(int hyp : {0,1}) {
            ba::aln_meth_opt_t mo;mo.m.meth_scoring=scoring;mo.m.meth_context=context;
            mo.m.meth_chem=METH_CHEM_TAPS;ba::aln_meth_opt_finish(&mo);
            const int strand=hyp^read_no;
            const std::string input_read=strand?aln_test::revcomp(query):query;
            std::set<Key> expected;
            // Every contig has exactly one full-length window; any gap costs more than the cutoff.
            for(int rid=0;rid<(int)contigs.size();++rid) {
                int mismatches=0;
                const auto &target=contigs[rid].second;
                for(int k=0;k<(int)query.size();++k) {
                    const char r=target[k],q=query[k];
                    const char from=hyp?'C':'G',to=hyp?'T':'A';
                    const bool cpg=hyp?(r=='C' && k+1<(int)target.size() && target[k+1]=='G')
                                      :(r=='G' && k>0 && target[k-1]=='C');
                    const bool conversion=r==from && q==to && (context==METH_CTX_ALL || cpg);
                    const bool mirror=context==METH_CTX_ALL && scoring==MEM_METH_SCORING_COLLAPSED && r==to && q==from;
                    mismatches+=r=='N' || (r!=q && !conversion && !mirror);
                }
                // A source N may have a seed fill mismatch; at most one N still fits discovery.
                if(mismatches<=1) expected.emplace(rid,hyp,strand,mr.bns->anns[rid].offset,mismatches*opt->s_mm);
            }
            std::set<Key> actual;
            const auto fq=write_fq("zero-bound.fq",{input_read});auto *input=ba::bwa_seq_open(fq.c_str());
            int n=0;auto *reads=ba::bwa_read_seq(input,1,&n,BWA_MODE_COMPREAD,0);REQUIRE(n==1);
            ba::aln_meth_read_t result;ba::aln_dup_state_t dup;
            REQUIRE(ba::aln_meth_search_batch(bwt,&mr,&mo,1,reads,opt,2,read_no,0,1,1000,&dup,&result)==0);
            CHECK(result.flags==0);
            for(const auto &o:result.occ) actual.emplace(o.rid,o.hyp,o.strand,o.pos,o.cost);
            CHECK(actual==expected);
            ba::bwa_free_read_seq(n,reads);ba::bwa_seq_close(input);
        }
    free(opt);ba::aln_meth_ref_free(&mr);ba::bwt_destroy(bwt);
}

TEST_CASE("aln_meth: nonzero stopping preserves independent candidate strata with custom penalties")
{
    namespace ba=bwa_aln;
    const std::string fa=ref().dir+"/nonzero-bound.fa";
    const std::string query="ACGTCGATCGTA";
    std::vector<std::pair<std::string,std::string>> contigs={{"best","C"+query.substr(1)}};
    for (int i=1;i<(int)query.size();++i) for(char base : std::string("ACGTN")) {
        std::string target="C"+query.substr(1); target[i]=base;
        contigs.emplace_back("v"+std::to_string(contigs.size()),target);
    }
    REQUIRE(aln_test::write_test_fasta(fa,contigs));
    std::vector<std::string> args={"index","--meth",fa};std::vector<char*> argv;
    for(auto &s:args)argv.push_back(&s[0]);argv.push_back(nullptr);optind=1;
#if defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || defined(__NetBSD__)
    optreset=1;
#endif
    REQUIRE(bwa_index((int)args.size(),argv.data())==0);
    FMI_search fmi((fa+".meth").c_str());fmi.load_index(false,1);
    auto *bwt=ba::bwt_from_fmi(&fmi);ba::aln_meth_ref_t mr;
    REQUIRE(ba::aln_meth_ref_load(&mr,fa.c_str(),fmi.idx->bns,1)==0);
    auto *opt=ba::gap_init_opt();opt->fnr=-1;opt->max_diff=3;opt->seed_len=1024;
    using Key=std::tuple<int,int,int,int64_t,int64_t>;
    for (const auto &penalties : {std::vector<int>{3,11,4},std::vector<int>{2,7,0},std::vector<int>{4,13,2}})
    for(int scoring : {MEM_METH_SCORING_COLLAPSED,MEM_METH_SCORING_GENOMIC,MEM_METH_SCORING_NEUTRAL})
        for(int context : {METH_CTX_ALL,METH_CTX_CPG}) for(int read_no : {0,1}) for(int hyp : {0,1}) {
            opt->s_mm=penalties[0];opt->s_gapo=penalties[1];opt->s_gape=penalties[2];
            ba::aln_meth_opt_t mo;mo.m.meth_scoring=scoring;mo.m.meth_context=context;
            mo.m.meth_chem=METH_CHEM_TAPS;ba::aln_meth_opt_finish(&mo);
            const int strand=hyp^read_no;
            const std::string input_read=strand?aln_test::revcomp(query):query;
            std::set<Key> expected;
            // Every contig has exactly one full-length window; any gap costs more than the cutoff.
            for(int rid=0;rid<(int)contigs.size();++rid) {
                int mismatches=0;
                const auto &target=contigs[rid].second;
                for(int k=0;k<(int)query.size();++k) {
                    const char r=target[k],q=query[k];
                    const char from=hyp?'C':'G',to=hyp?'T':'A';
                    const bool cpg=hyp?(r=='C' && k+1<(int)target.size() && target[k+1]=='G')
                                      :(r=='G' && k>0 && target[k-1]=='C');
                    const bool conversion=r==from && q==to && (context==METH_CTX_ALL || cpg);
                    const bool mirror=context==METH_CTX_ALL && scoring==MEM_METH_SCORING_COLLAPSED && r==to && q==from;
                    mismatches+=r=='N' || (r!=q && !conversion && !mirror);
                }
                // A source N may have a seed fill mismatch; at most two differences still fit discovery.
                if(mismatches<=2) expected.emplace(rid,hyp,strand,mr.bns->anns[rid].offset,mismatches*opt->s_mm);
            }
            std::set<Key> actual;
            const auto fq=write_fq("nonzero-bound.fq",{input_read});auto *input=ba::bwa_seq_open(fq.c_str());
            int n=0;auto *reads=ba::bwa_read_seq(input,1,&n,BWA_MODE_COMPREAD,0);REQUIRE(n==1);
            ba::aln_meth_read_t result;ba::aln_dup_state_t dup;
            REQUIRE(ba::aln_meth_search_batch(bwt,&mr,&mo,1,reads,opt,2,read_no,0,1,1000,&dup,&result)==0);
            CHECK(result.flags==0);
            for(const auto &o:result.occ) actual.emplace(o.rid,o.hyp,o.strand,o.pos,o.cost);
            CHECK(actual==expected);
            const int complete_c1=result.c1,complete_c2=result.c2;
            ba::bwa_free_read_seq(n,reads);ba::bwa_seq_close(input);
            // -N disables stopping but must produce exactly the same retained <=best+M strata.
            opt->mode|=BWA_MODE_NONSTOP;
            input=ba::bwa_seq_open(fq.c_str());reads=ba::bwa_read_seq(input,1,&n,BWA_MODE_COMPREAD,0);
            ba::aln_meth_read_t nonstop;ba::aln_dup_state_t dup_nonstop;
            REQUIRE(ba::aln_meth_search_batch(bwt,&mr,&mo,1,reads,opt,2,read_no,0,1,1000,&dup_nonstop,&nonstop)==0);
            CHECK(nonstop.flags==0);actual.clear();int within_c2=0;
            for(const auto &o:nonstop.occ) if(o.cost<=2*opt->s_mm) {
                actual.emplace(o.rid,o.hyp,o.strand,o.pos,o.cost);
                within_c2+=o.cost>opt->s_mm;
            }
            // Neutral conversion-displaced best ties also contribute to X1.
            within_c2+=nonstop.n_conv_moved;
            CHECK(actual==expected);CHECK(nonstop.c1==complete_c1);CHECK(within_c2==complete_c2);
            opt->mode&=~BWA_MODE_NONSTOP;
            ba::bwa_free_read_seq(n,reads);ba::bwa_seq_close(input);
            if(scoring==MEM_METH_SCORING_GENOMIC && context==METH_CTX_ALL && read_no==1 && hyp==1) {
                mo.p_best=1;mo.p_second=0;
                std::set<Key> capped_loci;
                long first_lookups=-1,first_windows=-1;int first_flags=-1;
                for(bool no_stop : {false,true}) {
                    if(no_stop) opt->mode|=BWA_MODE_NONSTOP;
                    input=ba::bwa_seq_open(fq.c_str());reads=ba::bwa_read_seq(input,1,&n,BWA_MODE_COMPREAD,0);
                    ba::aln_meth_read_t capped;ba::aln_dup_state_t capped_dup;
                    REQUIRE(ba::aln_meth_search_batch(bwt,&mr,&mo,1,reads,opt,2,read_no,0,1,1000,&capped_dup,&capped)==0);
                    actual.clear();for(const auto &o:capped.occ)actual.emplace(o.rid,o.hyp,o.strand,o.pos,o.cost);
                    CHECK(capped.flags!=0);
                    if(!no_stop) {
                        capped_loci=actual;first_flags=capped.flags;
                        first_lookups=capped.n_lookups;first_windows=capped.n_rescored;
                    } else {
                        CHECK(actual==capped_loci);CHECK(capped.flags==first_flags);
                        CHECK(capped.n_lookups==first_lookups);CHECK(capped.n_rescored==first_windows);
                    }
                    ba::bwa_free_read_seq(n,reads);ba::bwa_seq_close(input);
                }
                opt->mode&=~BWA_MODE_NONSTOP;
            }
        }
    free(opt);ba::aln_meth_ref_free(&mr);ba::bwt_destroy(bwt);
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
    SUBCASE("complete ends beyond 500 move their stored CIGAR") {
        const int rid=4;
        for(int k=0;k<600;++k) {
            auto extra=result[0].occ[0];extra.rid=rid;extra.pos=mr.bns->anns[rid].offset+50+100*k;
            extra.cost=3;extra.n_mm=1;extra.alignment.rb=extra.pos;extra.alignment.re=extra.pos+50;
            result[0].occ.push_back(extra);
        }
        result[0].c2+=600;
    }
    SUBCASE("work-incomplete end is not promoted") {incomplete=true;result[0].flags=ba::ALN_METH_BIND_M;}
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
    REQUIRE(seqs[0]->n_cigar==3);
    CHECK(__cigar_op(seqs[0]->cigar[1])==1);
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
    CHECK(tag_value(v[89],"XM:Z").size()==50);
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
    bool conversion_guard=false,cap_bound=false;
    SUBCASE("ordinary anchors") {}
    SUBCASE("cap-bound no-hit ends use bounded rescue without proper-pair claims") {cap_bound=true;}
    SUBCASE("conversion-selected anchors retain their ceiling") {conversion_guard=true;}
    struct Case {int rid,strand,anchor,target_rid,target;bool rescued;};
    const Case cases[]={{0,0,2920,1,70,false},{1,1,110,1,10,true},
                        {1,0,50,1,180,true},{0,1,2900,0,2800,true},
                        {0,0,2950,1,70,false},{1,1,0,0,2900,false}};
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
        if(cap_bound) result[1].flags=ba::ALN_METH_BIND_BEST;
        if(conversion_guard) {result[0].n_conv_moved=1;result[0].conv_gap=3;result[0].c2=1;}
        ba::aln_meth_pe_stage1(opt,1,seqs,r);
        if(conversion_guard) CHECK((int)seqs[0]->mapQ==18);
        ba::isize_info_t ii={};ii.avg=190;ii.std=10;ii.ap_prior=1e-5;ii.low=100;ii.high=250;ii.high_bayesian=250;
        uint64_t total[2]={},mapped[2]={};ba::aln_meth_pe_stage2(&mr,&mo,opt,pe,&ii,1,seqs,r,total,mapped);
        if(c.rescued) {
            CHECK((int)seqs[1]->type==BWA_TYPE_MATESW);
            CHECK(seqs[1]->pos==(bwtint_t)(mr.bns->anns[c.rid].offset+c.target));
            CHECK(((seqs[1]->extra_flag&SAM_FPP)!=0)==!cap_bound);
            CHECK(result[1].hyp==o.hyp);
            if(cap_bound) {
                CHECK((int)seqs[1]->mapQ==0);CHECK((int)seqs[1]->c1==0);CHECK((int)seqs[1]->c2==0);
                CHECK((int)seqs[0]->mapQ==37);CHECK((seqs[0]->extra_flag&SAM_FPP)==0);
                REQUIRE(result[1].occ.size()==1);CHECK(result[1].occ[0].cost==0);
                CHECK(result[1].occ[0].alignment.raw_cigar==std::vector<uint32_t>{50u<<4});
            }
            if(conversion_guard) {CHECK((int)seqs[0]->mapQ==10);CHECK((int)seqs[1]->mapQ==10);}
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
        std::vector<pair64_t> rows={{100,0},{200,3},{3100,4},{3200,7}};
        for(int k=0;k<600;++k) rows.push_back({(uint64_t)(mr.bns->anns[4].offset+50+100*k),0});
        pair64_v v;v.a=rows.data();v.n=v.m=rows.size();
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
    bool large=false;
    SUBCASE("ordinary retained strata") {}
    SUBCASE("more than 500 retained distinct loci") {large=true;}
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
        for(int k=0;k<(large?600:0);++k) {
            auto extra=result[0].occ[0];extra.rid=4;extra.pos=mr.bns->anns[4].offset+50+100*k;
            extra.cost=3;extra.n_mm=1;extra.alignment.rb=extra.pos;extra.alignment.re=extra.pos+50;
            result[0].occ.push_back(extra);
        }
        result[0].c2+=large?600:0;
        ba::aln_meth_pe_stage1(opt,1,seqs,r);ba::isize_info_t ii={};ii.avg=ii.std=-1;
        uint64_t total[2]={},mapped[2]={};ba::aln_meth_pe_stage2(&mr,&mo,opt,pe,&ii,1,seqs,r,total,mapped);
        CHECK(result[1].main==1);CHECK(seqs[1]->pos==2900);CHECK(seqs[1]->score==3);
        CHECK((seqs[0]->extra_flag&SAM_FPP)!=0);CHECK((seqs[1]->extra_flag&SAM_FPP)!=0);
        // The large worse stratum makes the unchanged anchor MAPQ 0 under stock accounting.
        CHECK((int)seqs[1]->mapQ==(large?0:29));CHECK(ba::pos_5(seqs[1])-ba::pos_5(seqs[0])==150);
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
    // Meth MD reports original ambiguity, independently of the index's random PAC fill.
    CHECK(tag_value(n[0],"MD:Z")=="30N49");
}

TEST_CASE("aln_meth: TAPS neutral prefers the locus needing the fewest conversions; MAPQ is capped")
{
    const std::string fq = write_fq("conv.fq", {ref().read_u, ref().read_v});
    for (const char *chem : {"--meth=taps", "--meth=5base"}) {
        auto v = run_aln({chem, "-n", "2", "-l", "1024"}, {fq});
        REQUIRE(v.size() == 2);
        // read_u: the exact locus is unique; the three needing conversions move to the second
        // stratum; one more conversion than the fewest caps MAPQ at 6
        CHECK(v[0].rname == "c4");
        CHECK(v[0].pos == ref().u_pos);
        CHECK(v[0].tag["XT:A"] == "U");
        CHECK(v[0].tag["X0:i"] == "1");
        CHECK(v[0].tag["X1:i"] == "3");
        CHECK(v[0].mapq == 6);
        // read_v: two exact loci stay tied (a repeat, MAPQ 0), and the main hit is one of them
        CHECK(v[1].tag["XT:A"] == "R");
        CHECK(v[1].tag["X0:i"] == "2");
        CHECK(v[1].tag["X1:i"] == "1");
        CHECK(v[1].mapq == 0);
        CHECK((v[1].pos == ref().v_pos[0] || v[1].pos == ref().v_pos[1]));
    }
    const auto repeated=run_aln({"--meth=taps","-n","2","-l","1024","-t","4"},
        {write_fq("conv-primary.fq",std::vector<std::string>(64,ref().read_u))});
    REQUIRE(repeated.size()==64);
    for(const auto &record:repeated) {CHECK(record.pos==ref().u_pos);CHECK(record.mapq==6);}
    // em-seq (any scoring), and TAPS genomic, keep the plain stratum: four (three) tied loci
    for (std::vector<std::string> args : {std::vector<std::string>{"--meth", "--meth-scoring", "neutral"},
                                          std::vector<std::string>{"--meth=taps", "--meth-scoring", "genomic"}}) {
        args.insert(args.end(), {"-n", "2", "-l", "1024"});
        auto v = run_aln(args, {fq});
        REQUIRE(v.size() == 2);
        CHECK(v[0].tag["X0:i"] == "4");
        CHECK(v[0].mapq == 0);
        CHECK(v[1].tag["X0:i"] == "3");
    }
}

TEST_CASE("aln_meth: TAPS neutral's tiebreak in pairs: the end's main hit and MAPQ before pairing")
{
    // read 1 is read_u (forward, OT); read 2 a G->A read from the reverse strand 150 bp on
    const std::string &c4 = ref().c4;
    const int p2 = ref().u_pos - 1 + 150;
    REQUIRE(p2 + 40 <= (int)c4.size());
    const std::string r2 = conv(aln_test::revcomp(c4.substr(p2, 40)), 'G', 'A');
    auto v = run_aln({"--meth=taps", "-n", "2", "-l", "1024"},
                     {write_fq("cp1.fq", {ref().read_u}), write_fq("cp2.fq", {r2})});
    REQUIRE(v.size() == 2);
    CHECK(v[0].pos == ref().u_pos);
    CHECK(v[0].tag["X0:i"] == "1");
    CHECK(v[0].tag["X1:i"] == "3");
    CHECK(v[1].pos == p2 + 1);
    CHECK(v[0].flag == 99); CHECK(v[1].flag == 147);
    CHECK(v[0].mapq == 6); CHECK(v[1].mapq == 43);
    CHECK(v[0].cigar == "40M"); CHECK(v[1].cigar == "40M");
    CHECK(v[0].tlen == 190); CHECK(v[1].tlen == -190);
}


TEST_CASE("aln_meth: neutral lexicographic strata honor stock cost and incomplete state")
{
    namespace ba=bwa_aln;
    auto *opt=ba::gap_init_opt(); opt->s_mm=6;opt->s_gapo=6;opt->s_gape=0;
    ba::aln_meth_opt_t mo;mo.m.meth_chem=METH_CHEM_TAPS;
    mo.m.meth_scoring=MEM_METH_SCORING_NEUTRAL;ba::aln_meth_opt_finish(&mo);
    ba::aln_meth_group_t group;group.mo=&mo;group.opt=opt;
    for(int gap : {1,2,9,10}) for(bool incomplete : {false,true}) {
        ba::aln_meth_read_t result;
        for(int k=0;k<4;++k) {
            ba::aln_meth_occ_t o;o.rid=0;o.hyp=1;o.strand=0;o.pos=100+100*k;
            o.cost=k==3?12:6;o.n_conv=k==0?2:k==1?2+gap:1;
            o.n_mm=k==0?1:0;o.n_gapo=k==1?1:0;
            o.alignment.rb=o.pos;o.alignment.re=o.pos+40;
            o.alignment.cigar={40u<<4};o.alignment.raw_cigar=o.alignment.cigar;
            result.occ.push_back(o);
        }
        result.flags=incomplete?ba::ALN_METH_BIND_M:0;
        ba::aln_meth_slot_t slot;slot.grp=&group;slot.res=&result;slot.keep=1;
        ba::aln_meth_slot_finish(&slot);
        CHECK(result.occ[result.main].cost==6);
        if(incomplete) { CHECK(result.n_conv_moved==0);CHECK(result.conv_gap==0); }
        else {
            // The third cost-6 locus has fewer conversions than both other compositions.
            CHECK(result.occ[result.main].pos==300);
            CHECK(result.c1==1); CHECK(result.c2==3);
            CHECK(result.n_conv_moved==2); CHECK(result.conv_gap==1);
            for(int ordinal=0;ordinal<40;++ordinal) {
                ba::aln_meth_redraw(&result,ordinal);
                CHECK(result.occ[result.main].pos==300);
            }
        }
    }
    free(opt);
}

TEST_CASE("aln_meth: conversion-only unchanged pairs retain both guards while moved loci use pair evidence")
{
    namespace ba=bwa_aln;
    const std::string prefix=ref().fa+".meth";
    FMI_search fmi(prefix.c_str());fmi.load_index(false,1);
    ba::aln_meth_ref_t mr;REQUIRE(ba::aln_meth_ref_load(&mr,ref().fa.c_str(),fmi.idx->bns,1)==0);
    ba::aln_meth_opt_t mo;mo.m.meth_chem=METH_CHEM_TAPS;
    mo.m.meth_scoring=MEM_METH_SCORING_NEUTRAL;ba::aln_meth_opt_finish(&mo);
    auto *opt=ba::gap_init_opt();auto *pe=ba::bwa_init_pe_opt();pe->is_sw=0;
    ba::bwase_initialize();
    bool zero_initial=false,large=false;
    SUBCASE("nonzero single-end confidence") {}
    SUBCASE("unchanged zero single-end confidence still retains conversion cap") {zero_initial=true;}
    SUBCASE("more than 500 retained loci preserve conversion guards") {large=true;}
    for(bool moved : {false,true}) {
        ba::bwa_seq_t *seqs[2];ba::bwa_seqio_t *input[2];int n[2];
        ba::aln_meth_read_t result[2];ba::aln_meth_read_t *r[2]={result,result+1};
        for(int j=0;j<2;++j) {
            input[j]=ba::bwa_seq_open(write_fq(j?"guard2.fq":"guard1.fq",{std::string(40,'A')}).c_str());
            seqs[j]=ba::bwa_read_seq(input[j],1,&n[j],BWA_MODE_COMPREAD,0);REQUIRE(n[j]==1);
            ba::aln_meth_occ_t o;o.rid=0;o.hyp=1;o.strand=j;o.pos=j?500:350;
            o.alignment.rb=o.pos;o.alignment.re=o.pos+40;o.alignment.cigar={40u<<4};
            result[j].occ.push_back(o);result[j].main=0;result[j].c1=1;result[j].c2=1;
            result[j].n_conv_moved=1;result[j].conv_gap=j?2:1;
            if(j==0 && moved) {
                result[j].occ[0].pos=20;
                auto alternate=o;alternate.n_conv=1;alternate.n_mm=1;alternate.cost=3;
                alternate.alignment.cigar={20u<<4,1u<<4|1,19u<<4};alternate.n_ins=1;alternate.n_gapo=1;
                alternate.alignment.re=alternate.pos+39;result[j].occ.push_back(alternate);
            }
        }
        if(large) for(int k=0;k<600;++k) {
            auto extra=result[0].occ[0];extra.rid=4;extra.pos=mr.bns->anns[4].offset+50+100*k;
            extra.cost=3;extra.n_mm=1;extra.alignment.rb=extra.pos;extra.alignment.re=extra.pos+40;
            result[0].occ.push_back(extra);
        }
        if(large) result[0].c2+=600;
        if(zero_initial) result[0].c2=300;
        ba::aln_meth_pe_stage1(opt,1,seqs,r);
        CHECK((int)seqs[0]->mapQ==((large||zero_initial)?0:6));CHECK((int)seqs[1]->mapQ==12);
        ba::isize_info_t ii={};ii.avg=ii.std=-1;
        uint64_t total[2]={},mapped[2]={};
        ba::aln_meth_pe_stage2(&mr,&mo,opt,pe,&ii,1,seqs,r,total,mapped);
        CHECK((int)seqs[0]->pos==350);CHECK((int)seqs[1]->pos==500);
        CHECK((seqs[0]->extra_flag&SAM_FPP)!=0);CHECK((seqs[1]->extra_flag&SAM_FPP)!=0);
        if(moved) {
            CHECK(result[0].main==1);CHECK((int)seqs[0]->mapQ==12);CHECK((int)seqs[0]->seQ==0);
            REQUIRE(seqs[0]->n_cigar==3);CHECK(__cigar_op(seqs[0]->cigar[1])==1);
            CHECK((int)seqs[0]->n_mm==1);CHECK((int)seqs[0]->n_gapo==1);
        } else {
            CHECK(result[0].main==0);CHECK((int)seqs[0]->mapQ==6);
            CHECK((int)seqs[0]->seQ==((large||zero_initial)?0:6));
        }
        CHECK((int)seqs[1]->mapQ==12);
        for(int j=0;j<2;++j) {ba::bwa_free_read_seq(1,seqs[j]);ba::bwa_seq_close(input[j]);}
    }
    free(pe);free(opt);ba::aln_meth_ref_free(&mr);
}

TEST_CASE("aln_meth: neutral conversion advantage applies the original six-per-conversion ceiling")
{
    namespace ba=bwa_aln;
    const std::string prefix=ref().fa+".meth";
    FMI_search fmi(prefix.c_str());fmi.load_index(false,1);
    ba::aln_meth_ref_t mr;REQUIRE(ba::aln_meth_ref_load(&mr,ref().fa.c_str(),fmi.idx->bns,1)==0);
    ba::aln_meth_opt_t mo;mo.m.meth_chem=METH_CHEM_TAPS;
    mo.m.meth_scoring=MEM_METH_SCORING_NEUTRAL;ba::aln_meth_opt_finish(&mo);
    auto *opt=ba::gap_init_opt();opt->max_diff=5;opt->fnr=-1;
    ba::aln_meth_group_t group;group.mo=&mo;group.opt=opt;ba::bwase_initialize();
    for(int gap : {1,2,3,4,10}) for(bool incomplete : {false,true}) {
        auto *input=ba::bwa_seq_open(write_fq("cap.fq",{ref().read_u}).c_str());int n=0;
        auto *reads=ba::bwa_read_seq(input,1,&n,BWA_MODE_COMPREAD,0);REQUIRE(n==1);
        ba::aln_meth_read_t result;
        for(int k=0;k<2;++k) {
            ba::aln_meth_occ_t o;o.rid=3;o.hyp=1;o.strand=0;
            o.pos=mr.bns->anns[3].offset+ref().u_pos-1+120*k;o.n_conv=k?gap:0;
            o.alignment.rb=o.pos;o.alignment.re=o.pos+40;
            o.alignment.cigar={40u<<4};o.alignment.raw_cigar=o.alignment.cigar;
            result.occ.push_back(o);
        }
        result.flags=incomplete?ba::ALN_METH_BIND_BEST:0;
        ba::aln_meth_slot_t slot;slot.grp=&group;slot.res=&result;slot.keep=1;
        ba::aln_meth_slot_finish(&slot);
        ba::aln_meth_finish_se(&mr,&mo,opt,3,1,reads,&result);
        CHECK((int)reads->mapQ==(incomplete?0:std::min(23,6*gap)));
        CHECK((int)reads->seQ==(int)reads->mapQ);
        CHECK((int)reads->c1==(incomplete?0:1));CHECK((int)reads->c2==(incomplete?0:1));
        CHECK((int)reads->n_multi==(incomplete?0:1));
        ba::bwa_free_read_seq(1,reads);ba::bwa_seq_close(input);
    }
    free(opt);ba::aln_meth_ref_free(&mr);
}

TEST_CASE("aln_meth: equal-cost same-locus windows retain the fewest-conversion raw CIGAR atomically")
{
    namespace ba=bwa_aln;
    const std::string prefix=ref().fa+".meth";
    FMI_search fmi(prefix.c_str());fmi.load_index(false,1);auto *bwt=ba::bwt_from_fmi(&fmi);
    ba::aln_meth_ref_t mr;REQUIRE(ba::aln_meth_ref_load(&mr,ref().fa.c_str(),fmi.idx->bns,1)==0);
    int rid=-1;for(int k=0;k<mr.bns->n_seqs;++k) if(std::string(mr.bns->anns[k].name)=="span") rid=k;
    REQUIRE(rid>=0);
    for(bool neutral : {false,true}) for(bool reverse_order : {false,true}) {
        auto *input=ba::bwa_seq_open(write_fq("span.fq",{ref().span_read}).c_str());int n=0;
        auto *reads=ba::bwa_read_seq(input,1,&n,BWA_MODE_COMPREAD,0);REQUIRE(n==1);
        ba::aln_meth_opt_t mo;mo.m.meth_chem=METH_CHEM_TAPS;
        mo.m.meth_scoring=neutral?MEM_METH_SCORING_NEUTRAL:MEM_METH_SCORING_GENOMIC;
        ba::aln_meth_opt_finish(&mo);
        auto *opt=ba::gap_init_opt();opt->s_mm=opt->s_gapo=3;opt->s_gape=0;
        ba::aln_meth_group_t group;group.bwt=bwt;group.ref=&mr;group.mo=&mo;group.opt=opt;
        group.read_no=1;group.n_slots=1;group.slot.resize(1);
        ba::gap_slot_core_t core={};group.cores=&core;
        auto &slot=group.slot[0];slot.grp=&group;slot.job=&core.job;
        ba::aln_meth_read_t result;ba::aln_meth_slot_start(&slot,reads,&result,0,3,1);
        for(int iteration=0;iteration<2;++iteration) {
            const int delta=(iteration^(int)reverse_order);
            const int span=40+delta;int strand=0;bwtint_t row=(bwtint_t)-1;
            for(bwtint_t k=0;k<=bwt->seq_len;++k) {
                auto pos=ba::bwa_pos_from_sa(mr.conv_bns,ba::bwt_sa(bwt,k),span,&strand);
                if(pos==(bwtint_t)mr.conv_bns->anns[2*rid+1].offset && strand==0) {row=k;break;}
            }
            REQUIRE(row!=(bwtint_t)-1);
            ba::gap_entry_t entry={};entry.n_del=delta;
            ba::aln_meth_hit(&slot,&entry,row,row);core.state=ba::GAP_SLOT_PARKED;
            ba::aln_meth_flush_parked(&group);
        }
        REQUIRE(result.occ.size()==1);
        const auto &o=result.occ[0];
        CHECK(o.cost==3);
        if(neutral) {
            CHECK(o.n_conv==0);CHECK(o.n_mm==0);CHECK(o.n_gapo==1);CHECK(o.n_del==1);
            CHECK(o.alignment.raw_cigar==std::vector<uint32_t>{20u<<4,1u<<4|2,20u<<4});
            CHECK(o.alignment.cigar==o.alignment.raw_cigar);
            CHECK(o.alignment.re-o.alignment.rb==41);
        }
        ba::bwa_free_read_seq(1,reads);ba::bwa_seq_close(input);free(opt);
    }
    ba::aln_meth_ref_free(&mr);ba::bwt_destroy(bwt);
}

TEST_CASE("aln_meth: an end with more than 500 loci still pairs when its strata are complete")
{
    // read 1 has 601 equally good loci, read 2 two; only c6 holds a proper pair. Mate rescue has no
    // anchor (both ends MAPQ 0), so only pairing on all 601 loci finds it.
    const std::string r1 = conv(ref().seg_p, 'C', 'T');
    const std::string r2 = conv(aln_test::revcomp(ref().frag_m), 'G', 'A');
    auto v = run_aln({"--meth"}, {write_fq("big1.fq", {r1}), write_fq("big2.fq", {r2})});
    REQUIRE(v.size() == 2);
    CHECK((v[0].flag & 0x2) != 0);
    CHECK(v[0].rname == "c6");
    CHECK(v[0].pos == ref().p_pos);
    CHECK(v[0].tag["X0:i"] == "601");
    CHECK(v[1].rname == "c6");
    CHECK(v[1].pos == ref().m_pos);
}


TEST_CASE("aln_meth: resolved distinct loci obey max-occ and incomplete lookup eligibility")
{
    namespace ba=bwa_aln;
    FMI_search fmi((ref().fa+".meth").c_str());fmi.load_index(false,1);
    auto *bwt=ba::bwt_from_fmi(&fmi);ba::aln_meth_ref_t mr;
    REQUIRE(ba::aln_meth_ref_load(&mr,ref().fa.c_str(),fmi.idx->bns,1)==0);
    auto *opt=ba::gap_init_opt();opt->seed_len=1024;opt->fnr=-1;opt->max_diff=2;
    auto *pe=ba::bwa_init_pe_opt();pe->is_sw=0;ba::bwase_initialize();
    for(bool neutral : {false,true}) for(int budget : {500,4596}) for(int limit : {600,601}) {
        CAPTURE(neutral);CAPTURE(budget);CAPTURE(limit);
        ba::aln_meth_opt_t mo;mo.m.meth_chem=neutral?METH_CHEM_TAPS:METH_CHEM_EMSEQ;
        mo.m.meth_scoring=neutral?MEM_METH_SCORING_NEUTRAL:MEM_METH_SCORING_COLLAPSED;
        mo.p_best=budget;mo.p_second=0;ba::aln_meth_opt_finish(&mo);pe->max_occ=limit;
        ba::bwa_seq_t *seqs[2];ba::bwa_seqio_t *input[2];int n[2];
        ba::aln_meth_read_t result[2];ba::aln_meth_read_t *res[2]={result,result+1};
        for(int j=0;j<2;++j) {
            const auto read=j?conv(aln_test::revcomp(ref().frag_m),'G','A'):conv(ref().seg_p,'C','T');
            input[j]=ba::bwa_seq_open(write_fq(j?"resolved2.fq":"resolved1.fq",{read}).c_str());
            seqs[j]=ba::bwa_read_seq(input[j],1,&n[j],BWA_MODE_COMPREAD,0);REQUIRE(n[j]==1);
            ba::aln_dup_state_t dup;
            REQUIRE(ba::aln_meth_search_batch(bwt,&mr,&mo,1,seqs[j],opt,2,1-j,0,1,3,&dup,res[j])==0);
        }
        std::set<std::tuple<int,int,int,bwtint_t>> keys;
        for(const auto &o:result[0].occ) keys.emplace(o.rid,o.hyp,o.strand,o.pos);
        CHECK(keys.size()==result[0].occ.size());
        const bool complete=budget==4596;
        if(complete) {CHECK(result[0].flags==0);CHECK(result[0].occ.size()==601);CHECK(result[0].c1==601);}
        else {CHECK(result[0].flags!=0);CHECK(result[0].occ.size()<=500);}
        // Force both single-end choices away from the sole compatible pair.
        for(int j=0;j<2;++j) for(size_t k=0;k<result[j].occ.size();++k) {
            const auto &o=result[j].occ[k];
            if(std::string(mr.bns->anns[o.rid].name)==(j?"c7":"c5")) {result[j].main=(int)k;break;}
        }
        ba::aln_meth_pe_stage1(opt,1,seqs,res);
        ba::isize_info_t ii={};ii.avg=ii.std=-1;uint64_t total[2]={},mapped[2]={};
        ba::aln_meth_pe_stage2(&mr,&mo,opt,pe,&ii,1,seqs,res,total,mapped);
        if(complete && limit==601) {
            CHECK((seqs[0]->extra_flag&SAM_FPP)!=0);CHECK((seqs[1]->extra_flag&SAM_FPP)!=0);
            CHECK(seqs[0]->pos==(bwtint_t)(mr.bns->anns[5].offset+ref().p_pos-1));
            CHECK(seqs[1]->pos==(bwtint_t)(mr.bns->anns[5].offset+ref().m_pos-1));
            CHECK((int)seqs[0]->mapQ==9);CHECK((int)seqs[1]->mapQ==9);
            CHECK((int)seqs[0]->seQ==0);CHECK((int)seqs[1]->seQ==0);
            CHECK(result[0].occ[result[0].main].pos==seqs[0]->pos);
            CHECK(seqs[0]->n_cigar==1);CHECK(__cigar_len(seqs[0]->cigar[0])==36);
        } else {
            CHECK((seqs[0]->extra_flag&SAM_FPP)==0);CHECK((seqs[1]->extra_flag&SAM_FPP)==0);
            CHECK((int)seqs[0]->mapQ==0);CHECK((int)seqs[1]->mapQ==0);
            CHECK((int)seqs[0]->type==BWA_TYPE_REPEAT);
            if(!complete) {
                CHECK((int)seqs[0]->c1==0);CHECK((int)seqs[0]->c2==0);CHECK(seqs[0]->n_multi==0);
                kstring_t text={0,0,nullptr};
                ba::aln_meth_print_sam1(&text,&mr,&mo,seqs[0],result,1,seqs[1],BWA_MODE_COMPREAD,30);
                const std::string record(text.s,text.l);free(text.s);
                CHECK(record.find("X0:i:")==std::string::npos);CHECK(record.find("X1:i:")==std::string::npos);
            }
        }
        for(int j=0;j<2;++j){ba::bwa_free_read_seq(1,seqs[j]);ba::bwa_seq_close(input[j]);}
    }
    free(pe);free(opt);ba::aln_meth_ref_free(&mr);ba::bwt_destroy(bwt);
}

TEST_CASE("aln_meth: resolved worse-cost loci beyond 255 pair without obsolete survivor stops")
{
    const std::string r1=conv(ref().seg_w,'C','T');
    const std::string r2=conv(aln_test::revcomp(ref().frag_n),'G','A');
    const auto hits=run_aln({"--meth","-n","1","-o","0","-l","1024","-R","1000","--max-occ","601"},
        {write_fq("second601-1.fq",{r1}),write_fq("second601-2.fq",{r2})});
    REQUIRE(hits.size()==2);
    CHECK(hits[0].rname=="c6");CHECK(hits[0].pos==ref().w_pos);
    CHECK(hits[1].rname=="c6");CHECK(hits[1].pos==ref().n_pos);
    CHECK(hits[0].flag==99);CHECK(hits[1].flag==147);
    CHECK(tag_value(hits[0],"X0:i")=="301");CHECK(tag_value(hits[0],"X1:i")=="300");
}

TEST_CASE("aln_meth: a cap-bound end is rescued at its resolved mate's loci")
{
    // 40 ordinary OT pairs (196-256 bp fragments) for the insert-size estimate, then the pair whose
    // read 1 is cap-bound (4801 equally good loci; XT:A:R, no X0, not paired) and whose read 2 has
    // two loci, MAPQ 0 (no anchor for ordinary rescue). Only c6 holds both.
    std::vector<std::string> r1s, r2s;
    for (int k = 0; k < 40; ++k) {
        const std::string f = ref().c1.substr(100 + 50 * k, 196 + (k * 7) % 61); // 196-256 bp
        r1s.push_back(conv(f.substr(0, 36), 'C', 'T'));
        r2s.push_back(conv(aln_test::revcomp(f).substr(0, 36), 'G', 'A'));
    }
    r1s.push_back(conv(ref().seg_q, 'C', 'T'));
    r2s.push_back(conv(aln_test::revcomp(ref().frag_r), 'G', 'A'));
    auto v = run_aln({"--meth"}, {write_fq("cap1.fq", r1s), write_fq("cap2.fq", r2s)});
    REQUIRE(v.size() == 82);
    Rec a = v[80], b = v[81];
    CHECK(a.rname == "c6");
    CHECK(a.pos == ref().q_pos);
    CHECK(a.tag["XT:A"] == "M");
    CHECK(a.tag.count("X0:i") == 0);
    CHECK((a.flag & 0x2) == 0);
    CHECK(a.tag["XG:Z"] == "CT");
    CHECK(b.rname == "c6");
    CHECK(b.pos == ref().r_pos);
    CHECK((b.flag & 0x2) == 0);
    CHECK(a.mapq == 0); // omitted genomic competitors still prohibit confidence
    CHECK(b.mapq == 0);
}


TEST_CASE("aln_meth: cap rescue accounts all bounded anchors and preserves unique neutral confidence")
{
    namespace ba=bwa_aln;
    FMI_search fmi((ref().fa+".meth").c_str());fmi.load_index(false,1);
    ba::aln_meth_ref_t mr;REQUIRE(ba::aln_meth_ref_load(&mr,ref().fa.c_str(),fmi.idx->bns,1)==0);
    auto *opt=ba::gap_init_opt();opt->s_mm=6;opt->s_gapo=5;opt->s_gape=1;
    auto *pe=ba::bwa_init_pe_opt();ba::bwase_initialize();
    int bad_range=0;bool changed_identity=false;
    SUBCASE("valid bounded insert model") {}
    SUBCASE("same-position different-hypothesis owned anchor is a moved selection") {changed_identity=true;}
    SUBCASE("nonfinite insert model is fatal") {bad_range=1;}
    SUBCASE("oversized bounded window is fatal") {bad_range=2;}
    SUBCASE("anchor outside its contig is fatal") {bad_range=3;}
    for(int anchors:{1,2,255,256}) for(bool cpg:{false,true}) {
        CAPTURE(anchors);CAPTURE(cpg);
        ba::aln_meth_opt_t mo;mo.m.meth_chem=METH_CHEM_TAPS;
        mo.m.meth_scoring=MEM_METH_SCORING_NEUTRAL;mo.m.meth_context=cpg?METH_CTX_CPG:METH_CTX_ALL;
        ba::aln_meth_opt_finish(&mo);
        ba::bwa_seq_t *seqs[2];ba::bwa_seqio_t *inputs[2];int n[2];
        const std::string mate=aln_test::revcomp(ref().c1.substr(700,50));
        for(int j=0;j<2;++j) {
            inputs[j]=ba::bwa_seq_open(write_fq(j?"anchors2.fq":"anchors1.fq",{j?mate:ref().c1.substr(550,50)}).c_str());
            seqs[j]=ba::bwa_read_seq(inputs[j],1,&n[j],BWA_MODE_COMPREAD,0);REQUIRE(n[j]==1);
        }
        ba::aln_meth_read_t result[2];ba::aln_meth_read_t *r[2]={result,result+1};
        for(int k=0;k<anchors;++k) {
            ba::aln_meth_occ_t a;a.rid=0;a.hyp=1;a.strand=0;a.pos=550+k;
            a.alignment.rb=a.pos;a.alignment.re=a.pos+50;a.alignment.cigar={50u<<4};
            result[0].occ.push_back(a);
        }
        result[0].main=anchors-1;result[0].c1=anchors;
        if(anchors==1) {result[0].c2=1;result[0].n_conv_moved=1;result[0].conv_gap=3;}
        if(changed_identity && anchors==1) {
            auto alternate=result[0].occ[0];alternate.hyp=0;
            alternate.alignment.cigar={25u<<4,1u<<4|1,24u<<4};alternate.alignment.re=alternate.pos+49;
            result[0].occ.insert(result[0].occ.begin(),alternate);result[0].main=1;
        }
        result[1].flags=ba::ALN_METH_BIND_BEST;
        ba::aln_meth_pe_stage1(opt,1,seqs,r);
        ba::isize_info_t ii={};ii.avg=200;ii.std=20;ii.high=0;ii.ap_prior=1e-5;
        if(bad_range && anchors<=255) {
            if(bad_range==1) ii.std=std::numeric_limits<double>::infinity();
            if(bad_range==2) {
                int rid=-1;for(int k=0;k<mr.bns->n_seqs;++k) if(std::string(mr.bns->anns[k].name)=="c9") rid=k;
                REQUIRE(rid>=0);ii.std=1000000;
                for(auto &a:result[0].occ) {a.rid=rid;a.pos=mr.bns->anns[rid].offset+100;
                    a.alignment.rb=a.pos;a.alignment.re=a.pos+50;}
            }
            if(bad_range==3) result[0].occ[0].alignment.re=mr.bns->anns[0].len+1;
            const pid_t child=fork();REQUIRE(child>=0);
            if(child==0) {
                // A fatal exit must not run the inherited fixture destructor in the child.
                std::atexit([] { _exit(1); });
                uint64_t total[2]={},mapped[2]={};
                ba::aln_meth_pe_stage2(&mr,&mo,opt,pe,&ii,1,seqs,r,total,mapped);
                _exit(0);
            }
            int status=0;REQUIRE(waitpid(child,&status,0)==child);
            CHECK(WIFEXITED(status));CHECK(WEXITSTATUS(status)==1);
            for(int j=0;j<2;++j){ba::bwa_free_read_seq(1,seqs[j]);ba::bwa_seq_close(inputs[j]);}
            continue;
        }
        uint64_t total[2]={},mapped[2]={};ba::aln_meth_pe_stage2(&mr,&mo,opt,pe,&ii,1,seqs,r,total,mapped);
        if(anchors<=255) {
            CHECK((int)seqs[1]->type==BWA_TYPE_MATESW);CHECK(seqs[1]->pos==700);
            CHECK((int)seqs[1]->mapQ==0);CHECK((int)seqs[1]->seQ==0);
            CHECK((int)seqs[1]->c1==0);CHECK((int)seqs[1]->c2==0);CHECK(seqs[1]->n_multi==0);
            CHECK((seqs[0]->extra_flag&SAM_FPP)==0);CHECK((seqs[1]->extra_flag&SAM_FPP)==0);
            CHECK(result[0].main==0);CHECK(seqs[0]->pos==550);
            CHECK((int)seqs[0]->mapQ==(anchors==1 && !changed_identity?18:0));
            if(changed_identity && anchors==1) {
                CHECK(result[0].hyp==0);REQUIRE(seqs[0]->n_cigar==3);
                CHECK(__cigar_op(seqs[0]->cigar[1])==1);
            }
            REQUIRE(result[1].occ.size()==1);CHECK(result[1].occ[0].cost==0);
            CHECK(result[1].occ[0].n_conv==0);CHECK(result[1].occ[0].n_mm==0);
        } else {
            CHECK((int)seqs[1]->type==BWA_TYPE_NO_MATCH);CHECK(result[1].main==-1);
            CHECK(result[0].main==255);CHECK((int)seqs[0]->mapQ==0);
        }
        for(int j=0;j<2;++j){ba::bwa_free_read_seq(1,seqs[j]);ba::bwa_seq_close(inputs[j]);}
    }
    free(pe);free(opt);ba::aln_meth_ref_free(&mr);
}

TEST_CASE("aln_meth: cap rescue ranks raw chemistry cost then neutral conversions independent of anchor order")
{
    namespace ba=bwa_aln;
    FMI_search fmi((ref().fa+".meth").c_str());fmi.load_index(false,1);
    ba::aln_meth_ref_t mr;REQUIRE(ba::aln_meth_ref_load(&mr,ref().fa.c_str(),fmi.idx->bns,1)==0);
    int rid=-1;for(int k=0;k<mr.bns->n_seqs;++k) if(std::string(mr.bns->anns[k].name)=="cap_choice") rid=k;
    REQUIRE(rid>=0);const int64_t base=mr.bns->anns[rid].offset;
    auto *opt=ba::gap_init_opt();opt->s_mm=6;opt->s_gapo=5;opt->s_gape=1;
    auto *pe=ba::bwa_init_pe_opt();ba::bwase_initialize();
    for(bool neutral:{false,true}) for(bool order:{false,true}) for(bool cpg:{false,true}) {
        ba::aln_meth_opt_t mo;mo.m.meth_chem=METH_CHEM_TAPS;
        mo.m.meth_scoring=neutral?MEM_METH_SCORING_NEUTRAL:MEM_METH_SCORING_GENOMIC;
        mo.m.meth_context=cpg?METH_CTX_CPG:METH_CTX_ALL;ba::aln_meth_opt_finish(&mo);
        ba::bwa_seq_t *seqs[2];ba::bwa_seqio_t *input[2];int n[2];
        for(int j=0;j<2;++j) {
            const auto seq=j?aln_test::revcomp(ref().cap_read):ref().c10.substr(100,50);
            input[j]=ba::bwa_seq_open(write_fq(j?"cap-choice2.fq":"cap-choice1.fq",{seq}).c_str());
            seqs[j]=ba::bwa_read_seq(input[j],1,&n[j],BWA_MODE_COMPREAD,0);REQUIRE(n[j]==1);
        }
        ba::aln_meth_read_t result[2];ba::aln_meth_read_t *r[2]={result,result+1};
        for(int k=0;k<2;++k) {
            ba::aln_meth_occ_t a;a.rid=rid;a.hyp=1;a.strand=0;a.pos=base+100+400*(k^(int)order);
            a.alignment.rb=a.pos;a.alignment.re=a.pos+50;a.alignment.cigar={50u<<4};result[0].occ.push_back(a);
        }
        result[0].main=0;result[0].c1=2;result[1].flags=ba::ALN_METH_BIND_BEST;
        // Existing best-known cost zero still accepts neutral conversions at zero cost.
        ba::aln_meth_occ_t known;known.rid=rid;known.hyp=1;known.strand=1;known.pos=base;
        known.alignment.rb=base;known.alignment.re=base+50;known.alignment.cigar={50u<<4};
        result[1].occ.push_back(known);result[1].main=0;
        ba::aln_meth_pe_stage1(opt,1,seqs,r);
        ba::isize_info_t ii={};ii.avg=200;ii.std=5;ii.ap_prior=1e-5;
        uint64_t total[2]={},mapped[2]={};ba::aln_meth_pe_stage2(&mr,&mo,opt,pe,&ii,1,seqs,r,total,mapped);
        REQUIRE(result[1].occ.size()==1);const auto &o=result[1].occ[0];
        CHECK(o.pos==(bwtint_t)(base+(neutral?650:250)));
        CHECK(o.n_conv==(neutral?0:1));CHECK(o.n_mm==0);CHECK(o.cost==0);
        CHECK(o.alignment.raw_cigar==std::vector<uint32_t>{50u<<4});
        CHECK(seqs[0]->pos==(bwtint_t)(base+(neutral?500:100)));
        CHECK((int)seqs[1]->type==BWA_TYPE_MATESW);CHECK((int)seqs[1]->mapQ==0);
        CHECK((int)seqs[0]->mapQ==0);CHECK((seqs[0]->extra_flag&SAM_FPP)==0);
        for(int j=0;j<2;++j){ba::bwa_free_read_seq(1,seqs[j]);ba::bwa_seq_close(input[j]);}
    }
    free(pe);free(opt);ba::aln_meth_ref_free(&mr);
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
        // Paired (mate rescue runs too): both mates are reported, neither aborts the run.
        auto pe = run_aln({mode}, {write_fq("polyc_1.fq", {adapter}), write_fq("polyc_2.fq", {std::string(150, 'G')})}, fa);
        REQUIRE(pe.size() == 2);
    }
    // EM-seq frees read C over reference T after conversion, so the same read still places.
    auto em = run_aln({"--meth=emseq"}, {write_fq("polyc_em.fq", {polyc})}, fa);
    REQUIRE(em.size() == 1);
    CHECK((em[0].flag & 4) == 0);
}
