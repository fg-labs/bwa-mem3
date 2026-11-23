// Unit tests for `aln --meth`: orientation, chemistry and tags, stock cost ranking,
// bounded candidate discovery, conservative caps, duplicate/thread invariance,
// and stock alignment against a --meth-built index.

#include "doctest/doctest.h"

#include "aln_test_ref.h"
#include "../framework/stderr_capture.h"
#include "../framework/run_in_child.h"
#include "../framework/temp_dir.h"
#include "../../src/aln/aln_main.h"
#include "../../src/aln/aln_meth.h"
#include "../../src/aln/aln_meth_search.h"
#include "../../src/aln/aln_known_snvs.h"
#include "../../src/aln/aln_meth_realign.h"
#include "../../src/aln/bwaseqio.h"
#include "../../src/aln/bwase.h"
#include "../../src/aln/aln_out.h"
#include "../../src/FMI_search.h"
#include "../../src/bwa.h"
#include "../../src/main.h"

#include <getopt.h>
#include <algorithm>
#include <cmath>
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

// A finished read's stored CIGAR, as a vector.
std::vector<uint32_t> cigar_of(const bwa_aln::aln_meth_read_t &r, const bwa_aln::aln_meth_hit_t &h)
{
    const uint32_t *c = bwa_aln::aln_meth_hit_cigar(r, h);
    return std::vector<uint32_t>(c, c + h.n_cig);
}

// Stores hand-built reads' loci as hits, as finishing does. It skips finishing's drop of an
// incomplete read's worse-cost loci: the tests that use it pick each end's loci themselves.
void store_pair(bwa_aln::aln_meth_read_t *r[2])
{
    bwa_aln::aln_meth_store_loci(r[0]); bwa_aln::aln_meth_store_loci(r[1]);
}

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
    std::string c4, read_u, read_v, span_read, span_reference;
    std::string c5,c6,c7,c8,seg_p,frag_m,seg_w,frag_n;
    std::string c10,cap_read;
    int p_pos=0,m_pos=0,w_pos=0,n_pos=0;
    int u_pos = 0, v_pos[2] = {0, 0};
    int x_pos = 0; // 1-based position of the counterexample's true locus in c3

    MethRef()
    {
        dir = tmp.path();
        fa = tmp.unique_file();
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
        cap_read=aln_test::lcg_dna(st,50);cap_read[20]='T';cap_read[21]='G';
        std::string converted_cap=cap_read;converted_cap[20]='C';
        c10=aln_test::lcg_dna(st,250)+converted_cap+aln_test::lcg_dna(st,350)+cap_read+aln_test::lcg_dna(st,100);
        REQUIRE(aln_test::write_test_fasta(fa, {{"c1", c1}, {"c2", c2}, {"c3", c3},
                                              {"c4", c4}, {"c5",c5},{"c6",c6},{"c7",c7},{"c8",c8},{"cap_choice",c10}, {"span",span_reference}, {"affine_expensive", expensive}, {"affine_cheap", cheap}}));
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

// Writes `reads` as FASTQ records r0, r1, ... to dir/name and returns the path.
std::string write_fq_in(const std::string &dir, const std::string &name, const std::vector<std::string> &reads)
{
    std::string path = dir + "/" + name;
    FILE *fp = fopen(path.c_str(), "w");
    REQUIRE(fp != nullptr);
    for (size_t i = 0; i < reads.size(); ++i)
        fprintf(fp, "@r%zu\n%s\n+\n%s\n", i, reads[i].c_str(), std::string(reads[i].size(), 'I').c_str());
    REQUIRE(fclose(fp) == 0);
    return path;
}

std::string write_fq(const std::string &name, const std::vector<std::string> &reads)
{
    return write_fq_in(ref().dir, name, reads);
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

// Runs aln and returns its records. The output goes to out_dir/out.sam (default ref().dir).
std::vector<Rec> run_aln(std::vector<std::string> args, const std::vector<std::string> &inputs,
                         const std::string &reference = "", const std::string &out_dir = "")
{
    const std::string out = (out_dir.empty() ? ref().dir : out_dir) + "/out.sam";
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
    // every mode places the read at the same alignment, so the NM difference is the scoring's
    for (const auto *v : {&neutral, &genomic, &collapsed}) {
        CHECK((*v)[0].flag == 0);
        CHECK((*v)[0].rname == "c2");
        CHECK((*v)[0].pos == 301);
        CHECK((*v)[0].cigar == "40M");
    }
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
    // an exact forward copy of c1[200, 240), unique in the reference: the stock record is known
    CHECK(v[0].flag == 0);
    CHECK(v[0].rname == "c1");
    CHECK(v[0].pos == 201);
    CHECK(v[0].cigar == "40M");
    CHECK(tag_value(v[0], "XT:A") == "U");
    CHECK(tag_value(v[0], "NM:i") == "0");
    CHECK(tag_value(v[0], "MD:Z") == "40");
    CHECK(has_tag(v[0], "XM:i"));
    CHECK(v[0].tag.count("XR:Z") == 0);
    CHECK(v[0].tag.count("XG:Z") == 0);
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

TEST_CASE("aln_meth: identical reads share a search and still draw their own main hit")
{
    // a batch of mostly duplicated reads, so the duplicate pass runs, with repeats (several best
    // hits) so the draw matters
    std::vector<std::string> reads;
    for (int k = 0; k < 40; ++k) {
        reads.push_back(ref().read_r);                                   // a read with 31 repeat loci
        reads.push_back(ref().read_v); // the neutral conversion split must replay per read
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

TEST_CASE("aln_meth: aln_opt_init clears the --meth and --known-snvs pointers of a reused aln_opt_t")
{
    namespace ba = bwa_aln;
    ba::aln_opt_t o;
    ba::aln_meth_opt_t mo;
    ba::aln_meth_ref_t mr;
    ba::ks_sites_t sites;
    ba::ks_windows_t windows;
    o.meth = &mo;
    o.meth_ref = &mr;
    o.known_snvs = &sites;
    o.ks_windows = &windows;
    ba::aln_opt_init(&o);
    CHECK(o.meth == nullptr);
    CHECK(o.meth_ref == nullptr);
    CHECK(o.known_snvs == nullptr);
    CHECK(o.ks_windows == nullptr);
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
    // the full payload, computed on the contig the record names
    REQUIRE(v[0].flag == 0);
    REQUIRE((v[0].rname == "c1" || v[0].rname == "c2"));
    const std::string &own = v[0].rname == "c1" ? ref().c1 : ref().c2;
    CHECK(tag_value(v[0], "XM:Z") == expected_xm(own, v[0].pos - 1, conv(x, 'C', 'T'), true, false));
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
    // Independent oracle for the cheap placement: affine_cheap is the read without bases
    // [28, 33), so a 5I CIGAR is aM5I(55-a)M for each a where removing read[a, a+5) gives it.
    // The equivalent placements are equal-cost; the record takes the leftmost (the first a).
    const std::string &read = ref().affine_read;
    const std::string cheap = read.substr(0, 28) + read.substr(33);
    std::vector<std::string> cheap_cigars;
    for (int a = 1; a < 55; ++a)
        if (read.substr(0, a) + read.substr(a + 5) == cheap)
            cheap_cigars.push_back(std::to_string(a) + "M5I" + std::to_string(55 - a) + "M");
    REQUIRE(!cheap_cigars.empty());
    const std::string cheap_cigar = cheap_cigars.front();
    auto placed = [](const Rec &r, const char *contig, const std::string &cigar) {
        CHECK(r.flag == 0);
        CHECK(r.rname == contig);
        CHECK(r.pos == 1);
        CHECK(r.cigar == cigar);
    };
    for (bool nonstop : {false, true}) {
        std::vector<std::string> options = {"--meth=taps", "-n", "6", "-l", "1024",
                                            "-M", "6", "-O", "1", "-E", "1"};
        if (nonstop) options.push_back("-N");
        const auto hits = run_aln(options, {fq});
        REQUIRE(hits.size() == 1);
        placed(hits[0], "affine_cheap", cheap_cigar); // cost 5 vs 12, despite rewards 45 vs 52
        CHECK(tag_value(hits[0], "NM:i") == "5");
        CHECK(tag_value(hits[0], "MD:Z") == "55");
        CHECK(tag_value(hits[0], "X0:i") == "1");
        // Ordinary cost window ends at 11; -N also retains two worse-cost loci, the cost-12
        // mismatch placement and a shifted gapped placement, each scored along its search path.
        CHECK(tag_value(hits[0], "X1:i") == (nonstop ? "2" : "0"));
        CHECK(hits[0].mapq == (nonstop ? 20 : 37));
    }
    for (const auto &penalties : {std::vector<std::string>{"0","1","1"},
                                   std::vector<std::string>{"6","0","0"}}) {
        const auto zero = run_aln({"--meth=taps", "-n", "6", "-l", "1024",
                                   "-M", penalties[0], "-O", penalties[1], "-E", penalties[2]}, {fq});
        REQUIRE(zero.size()==1);
        if (penalties[0]=="0") {
            placed(zero[0], "affine_expensive", "60M");
            CHECK(tag_value(zero[0],"NM:i")=="2");
        } else {
            placed(zero[0], "affine_cheap", cheap_cigar);
            CHECK(tag_value(zero[0],"NM:i")=="5");
        }
        CHECK(zero[0].mapq==37);
    }
    // With expensive gaps the same refined candidates instead prefer two mismatches.
    const auto costly_gap = run_aln({"--meth=taps", "-n", "6", "-l", "1024",
                                    "-M", "3", "-O", "11", "-E", "4"}, {fq});
    REQUIRE(costly_gap.size() == 1);
    placed(costly_gap[0], "affine_expensive", "60M");
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

TEST_CASE("aln_meth: a refined window over the read's -n limit is not a locus")
{
    // The converted search frees a read C over a reference T (read 1, OT), so it proposes windows
    // within -n converted differences whose original-reference alignment has more. Stock aln never
    // reports such a hit, and neither does --meth: the window does not map the read, count toward
    // X0/X1 or set the genomic MAPQ. The conversion-blind cap still counts one that is within -n
    // with the mirror cell free (s, below).
    uint32_t st = 8675309;
    auto flank = [&]() { return aln_test::lcg_dna(st, 100); };
    auto swap_ag = [](char c) { return c == 'A' ? 'G' : 'A'; };
    // r (26 bp; a fractional -n 0.04 allows 2 differences): its only locus has 4 read-C/reference-T
    // mismatches, spread out so that the refinement stays ungapped.
    std::string r = aln_test::lcg_dna(st, 26), r_over;
    for (int i : {4, 10, 16, 22}) r[i] = 'C';
    r_over = r;
    for (int i : {4, 10, 16, 22}) r_over[i] = 'T';
    // s (30 bp): 3 read-C/reference-T mismatches (cost 9) at one locus, a one-base deletion
    // (cost 11, one difference) at another.
    std::string s = aln_test::lcg_dna(st, 30), s_over;
    for (int i : {5, 14, 23}) s[i] = 'C';
    s[15] = 'A';
    s_over = s;
    for (int i : {5, 14, 23}) s_over[i] = 'T';
    const std::string s_gap = s.substr(0, 15) + "G" + s.substr(15);
    // u (30 bp): 2 genuine A/G mismatches (cost 6) at one locus, 3 read-C/reference-T mismatches
    // (cost 9, within the retained best + M) at another.
    std::string u = aln_test::lcg_dna(st, 30), u_near, u_over;
    for (int i : {3, 13, 26}) u[i] = 'C';
    u[8] = 'A'; u[20] = 'G';
    u_near = u; u_near[8] = swap_ag(u[8]); u_near[20] = swap_ag(u[20]);
    u_over = u;
    for (int i : {3, 13, 26}) u_over[i] = 'T';
    bwa_tests::TempDir tmp{"aln_meth_strict_n"}; // this case's index, reads and output
    // a unique basename: shared memory is looked up by basename, which the directory does not isolate
    const std::string dir = tmp.path(), fa = tmp.unique_file();
    REQUIRE(aln_test::write_test_fasta(fa, {{"r_over", flank() + r_over + flank()},
                                            {"s_over", flank() + s_over + flank()},
                                            {"s_gap", flank() + s_gap + flank()},
                                            {"u_near", flank() + u_near + flank()},
                                            {"u_over", flank() + u_over + flank()}}));
    std::vector<std::string> index_args = {"index", "--meth", fa};
    std::vector<char *> argv;
    for (auto &arg : index_args) argv.push_back(&arg[0]);
    argv.push_back(nullptr);
    optind = 1;
#if defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || defined(__NetBSD__)
    optreset = 1;
#endif
    REQUIRE(bwa_index((int)index_args.size(), argv.data()) == 0);

    const std::string fq_r = write_fq_in(dir, "strict-n-r.fq", {r});
    // The fractional default (-n 0.04: 2 for 26 bp), an integer -n below 4, and -N: unmapped.
    for (const auto &options : {std::vector<std::string>{"--meth"},
                                std::vector<std::string>{"--meth", "-n", "3"},
                                std::vector<std::string>{"--meth", "-N"}}) {
        const auto v = run_aln(options, {fq_r}, fa, dir);
        REQUIRE(v.size() == 1);
        CHECK((v[0].flag & 4) != 0);
        CHECK(v[0].mapq == 0);
    }
    // Allowed 4 differences, the same window maps the read.
    const auto r4 = run_aln({"--meth", "-n", "4"}, {fq_r}, fa, dir);
    REQUIRE(r4.size() == 1);
    CHECK((r4[0].flag & 4) == 0);
    CHECK(r4[0].rname == "r_over");
    CHECK(r4[0].pos == 101);
    CHECK(r4[0].cigar == "26M");
    CHECK(tag_value(r4[0], "NM:i") == "4");

    // s: under -n 2 the cheaper 3-mismatch locus is dropped, so the gapped one is the unique hit;
    // under -n 3 the 3-mismatch locus is best and the gapped one is a retained competitor.
    const std::string fq_s = write_fq_in(dir, "strict-n-s.fq", {s});
    const auto s2 = run_aln({"--meth", "-n", "2"}, {fq_s}, fa, dir);
    REQUIRE(s2.size() == 1);
    CHECK(s2[0].rname == "s_gap");
    CHECK(s2[0].pos == 101); // the flank's 100 bp, then the planted G between s[14] (C) and s[15] (A)
    CHECK(s2[0].cigar == "15M1D15M");
    CHECK(tag_value(s2[0], "X0:i") == "1");
    CHECK(tag_value(s2[0], "X1:i") == "0");
    // the default --meth-mapq conversion-blind still counts the dropped locus: with the mirror
    // cell free its 3 mismatches cost nothing and it is within -n, so it beats the gapped hit (as
    // under collapsed) and the cap is 0; the genomic MAPQ is kept in YQ:i and under scoring
    CHECK(s2[0].mapq == 0);
    CHECK(tag_value(s2[0], "YQ:i") == "37");
    const auto s2_scoring = run_aln({"--meth", "--meth-mapq", "scoring", "-n", "2"}, {fq_s}, fa, dir);
    REQUIRE(s2_scoring.size() == 1);
    CHECK(s2_scoring[0].rname == "s_gap");
    CHECK(s2_scoring[0].mapq == 37);
    const auto s2_collapsed = run_aln({"--meth", "--meth-scoring", "collapsed", "-n", "2"}, {fq_s}, fa, dir);
    REQUIRE(s2_collapsed.size() == 1);
    CHECK(s2_collapsed[0].rname == "s_over"); // the oracle: collapsed keeps that locus, and prefers it
    const auto s3 = run_aln({"--meth", "-n", "3"}, {fq_s}, fa, dir);
    REQUIRE(s3.size() == 1);
    CHECK(s3[0].rname == "s_over");
    CHECK(tag_value(s3[0], "X0:i") == "1");
    CHECK(tag_value(s3[0], "X1:i") == "1");

    // u: the over-limit locus is within the retained cost window but does not count toward X1.
    const std::string fq_u = write_fq_in(dir, "strict-n-u.fq", {u});
    const auto u2 = run_aln({"--meth", "-n", "2"}, {fq_u}, fa, dir);
    REQUIRE(u2.size() == 1);
    CHECK(u2[0].rname == "u_near");
    CHECK(tag_value(u2[0], "NM:i") == "2");
    CHECK(tag_value(u2[0], "X0:i") == "1");
    CHECK(tag_value(u2[0], "X1:i") == "0");
    const auto u3 = run_aln({"--meth", "-n", "3"}, {fq_u}, fa, dir);
    REQUIRE(u3.size() == 1);
    CHECK(u3[0].rname == "u_near");
    CHECK(tag_value(u3[0], "X1:i") == "1");
}

TEST_CASE("aln_meth: the conversion-blind MAPQ takes each locus at its lowest mirror-free row")
{
    // A locus can hold several rows: search paths with a gap in different places. The genomic
    // search keeps the genomically cheapest (ties: the CIGAR); collapsed scoring, which the
    // conversion-blind cap matches, keeps the cheapest with the mirror cell (read C over reference
    // T) free. Read r (50 bp, read 1, unseeded) has p at offsets 16-18. Locus A is r with q in
    // place of p, one base longer, so deleting q[k] gives row Dk; B is the read's other locus.
    //   tie (p ACA, q ATAA, and a mismatch at offset 40): D1 and D2 cost 17 genomically, D1 with
    //     two mismatches, D2 with one and a mirror mismatch (mirror-free 14). The CIGAR tie-break
    //     keeps D1. B, a two-base deletion, costs 15: genomic MAPQ 23 (A within B + 3). With the
    //     mirror cell free A (14) beats B, so the cap is 0; collapsed places the read at A.
    //   replaced (p CAC, q TCTA): D0 costs 17 (two mismatches) and D3 20 (two mirror mismatches,
    //     mirror-free 14), so D0 replaces D3. B has four mismatches (12): genomic MAPQ 37 (A past
    //     B + 3). With the mirror cell free A (14) is within B + 3: the cap is collapsed's 23. At
    //     D0's mirror-free cost (17) A stays past B + 3 and the cap stays 37.
    // xa: A as the genomic search keeps it, in the scoring record's XA (none: past the cut)
    struct Case { const char *name, *p, *q; int mismatch_at; bool b_gapped; int scoring, blind; const char *at, *xa; };
    for (const Case &k : {Case{"tie", "ACA", "ATAA", 40, true, 23, 0, "A", "A,+101,17M1D33M,3;"},
                          Case{"replaced", "CAC", "TCTA", -1, false, 37, 23, "B", ""}}) {
        const std::string name = k.name;
        CAPTURE(name);
        uint32_t st = 20261123;
        auto flank = [&]() { return aln_test::lcg_dna(st, 100); };
        auto alt = [](char c) { return c == 'G' ? 'A' : 'G'; }; // never T: no mirror cell
        std::string r = aln_test::lcg_dna(st, 50);
        r.replace(16, 3, k.p);
        std::string a = r.substr(0, 16) + k.q + r.substr(19);
        if (k.mismatch_at >= 0) a[k.mismatch_at + 1] = alt(r[k.mismatch_at]); // A is one base longer
        std::string b = r;
        if (k.b_gapped) b.insert(30, std::string(1, alt(r[30])) + alt(alt(r[30])));
        else for (int i : {5, 25, 33, 44}) b[i] = alt(r[i]);
        bwa_tests::TempDir tmp{"aln_meth_blind_row"};
        const std::string dir = tmp.path(), fa = tmp.unique_file();
        REQUIRE(aln_test::write_test_fasta(fa, {{"A", flank() + a + flank()}, {"B", flank() + b + flank()}}));
        std::vector<std::string> index_args = {"index", "--meth", fa};
        std::vector<char *> argv;
        for (auto &arg : index_args) argv.push_back(&arg[0]);
        argv.push_back(nullptr);
        optind = 1;
#if defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || defined(__NetBSD__)
        optreset = 1;
#endif
        REQUIRE(bwa_index((int)index_args.size(), argv.data()) == 0);
        const std::string fq = write_fq_in(dir, "blind-row.fq", {r});
        const std::vector<std::string> base = {"--meth", "-n", "5", "-l", "1024"};
        auto run = [&](std::vector<std::string> extra) {
            std::vector<std::string> args = base;
            args.insert(args.end(), extra.begin(), extra.end());
            const auto v = run_aln(args, {fq}, fa, dir);
            REQUIRE(v.size() == 1);
            return v[0];
        };
        const Rec scoring = run({"--meth-mapq", "scoring"});
        CHECK(scoring.rname == "B");
        CHECK(scoring.mapq == k.scoring);
        CHECK((scoring.tag.count("XA:Z") ? scoring.tag.at("XA:Z") : std::string()) == k.xa); // D1, not D2
        const Rec blind = run({});
        CHECK(blind.rname == "B"); // placement is the genomic search's
        CHECK(blind.mapq == k.blind);
        CHECK(tag_value(blind, "YQ:i") == std::to_string(k.scoring));
        const Rec collapsed = run({"--meth-scoring", "collapsed"});
        CHECK(collapsed.rname == k.at);
        CHECK(collapsed.mapq == 23);
    }
}

TEST_CASE("aln_meth: aln_meth_blind_view takes a hit at its moved row, else at its own row mirror-free")
{
    // The single-end cap and the shadow pairing read every hit through aln_meth_blind_view: a hit
    // listed in blind_alt (sorted by position, then strand) is taken at that row whole (cost,
    // mismatches, gaps, span); any other at its own row less its mirror mismatches.
    namespace ba = bwa_aln;
    ba::aln_meth_read_t r;
    auto alt = [](bwtint_t pos, int strand, int32_t span, int32_t cost, int n_mm, int n_gapo, int n_gape) {
        ba::aln_meth_blind_row_t b;
        b.pos = pos; b.strand = (uint8_t)strand; b.span = span; b.cost = cost;
        b.n_mm = (uint8_t)n_mm; b.n_gapo = (uint8_t)n_gapo; b.n_gape = (uint8_t)n_gape;
        return b;
    };
    r.blind_alt = {alt(100, 0, 49, 14, 1, 1, 0), alt(100, 1, 52, 26, 1, 1, 1), alt(300, 0, 50, 6, 2, 0, 0)};
    ba::aln_meth_hit_t h = {};
    h.pos = 100; h.strand = 1; h.span = 50; h.cost = 30; h.n_mm = 3; h.n_gapo = 1; h.n_gape = 0; h.n_mirror = 1;
    const ba::aln_meth_blind_row_t moved = ba::aln_meth_blind_view(r, h, 3);
    CHECK(moved.pos == 100); CHECK((int)moved.strand == 1);
    CHECK(moved.span == 52); CHECK(moved.cost == 26);
    CHECK((int)moved.n_mm == 1); CHECK((int)moved.n_gapo == 1); CHECK((int)moved.n_gape == 1);
    h.pos = 200; // not listed: its own row, mirror-free
    const ba::aln_meth_blind_row_t own = ba::aln_meth_blind_view(r, h, 3);
    CHECK(own.pos == 200); CHECK((int)own.strand == 1);
    CHECK(own.span == 50); CHECK(own.cost == 27);
    CHECK((int)own.n_mm == 2); CHECK((int)own.n_gapo == 1); CHECK((int)own.n_gape == 0);
    h.pos = 300; h.strand = 1; // listed only on the other strand
    CHECK(ba::aln_meth_blind_view(r, h, 3).cost == 27);
}

TEST_CASE("aln_meth: the conversion-blind strata are collapsed's selection, locus by locus")
{
    // A property test of the conversion-blind strata: on a seeded synthetic reference, the strata
    // a genomic search builds with the mirror cell freed (blind_best, blind_c1, blind_c2, and each
    // locus within them: a hit at its mirror-free row, aln_meth_blind_view, or a blind_rows row)
    // must be what collapsed scoring, which frees the mirror cell itself, selects: its best, c1,
    // c2, each hit's cost, mismatches, gaps and span, and its main hit (blind_pick, the locus
    // blind_main draws among the mirror-free best loci with the read's seed). Each family is a 50-base
    // segment whose read 1 keeps some Cs; its copies carry mirror cells (a T under a kept C), real
    // mismatches and one-base gaps. Planted families: a one-gap best with a mirror-free best and a locus only
    // the mirror-free cut keeps, three best loci with a gapped mirror locus, and a copy past the
    // genomic cut. Reads where a search stops on -R are skipped: past such a stop the second
    // stratum is partial, as in stock aln. So are reads whose strata the stock cuts may have
    // pruned (blind_incomplete: their cap is 0). The hits compared are those the paired-end shadow
    // pairing keeps (aln_meth_blind_selected: within the mirror-free cut and, under the cuts, the
    // mirror-free ceiling). Every random value is drawn into a local first, so the
    // reference does not depend on the compiler's argument evaluation order.
    namespace ba = bwa_aln;
    bwa_tests::TempDir tmp{"aln_meth_blind_strata"};
    // a unique basename: shared memory is looked up by basename, which the directory does not isolate
    const std::string fa = tmp.unique_file();
    uint32_t st = 3;
    const int len = 50;
    auto rnd = [&](int n) { st = st * 1103515245u + 12345u; return (int)((st >> 8) % (uint32_t)n); };
    auto alt = [](char c) { return c == 'G' ? 'A' : 'G'; }; // never T: cannot match a converted C
    auto convert = [](std::string s, const std::set<int> &keep) {
        for (int i = 0; i < (int)s.size(); ++i)
            if (s[i] == 'C' && !keep.count(i)) s[i] = 'T';
        return s;
    };
    // A copy of seg: a T under each position of mirror, a mismatch against the read at each of
    // real (read is seg converted, so the substitute is never T and never C), then one gap at
    // gap_at (del: the copy lacks a base, so the read has an insertion; ins: the copy has one more).
    auto copy_of = [&](std::string seg, const std::set<int> &mirror, const std::set<int> &real, int gap_at, bool del) {
        const std::string read = convert(seg, {});
        for (int p : mirror) seg[p] = 'T';
        for (int p : real) seg[p] = alt(read[p]);
        if (gap_at >= 0) {
            if (del) seg.erase(gap_at, 1);
            else seg.insert(gap_at, 1, alt(seg[gap_at]));
        }
        return seg;
    };
    // A segment with Cs at the positions of meth (the methylated Cs a read can keep).
    auto segment = [&](const std::set<int> &meth) {
        std::string s = aln_test::lcg_dna(st, len);
        for (int p : meth) s[p] = 'C';
        return s;
    };
    // Distinct positions in [lo, hi) outside avoid.
    auto positions = [&](int n, int lo, int hi, const std::set<int> &avoid) {
        std::set<int> out;
        while ((int)out.size() < n) {
            const int p = lo + rnd(hi - lo);
            if (!avoid.count(p)) out.insert(p);
        }
        return out;
    };
    std::vector<std::pair<std::string, std::string>> contigs;
    std::vector<std::string> reads;
    auto family = [&](const std::vector<std::string> &copies) {
        std::string c = aln_test::lcg_dna(st, 80);
        for (const auto &k : copies) {
            const int spacer = 60 + rnd(60);
            c += k + aln_test::lcg_dna(st, spacer);
        }
        contigs.emplace_back("f" + std::to_string(contigs.size()), c);
    };
    {   // a one-gap best, a mirror-free best (three mismatches and a mirror), and a four-mismatch
        // locus only the mirror-free cut keeps
        const std::set<int> meth = {12};
        const std::string seg = segment(meth);
        const std::set<int> three = positions(3, 0, len, meth);
        const std::set<int> four = positions(4, 0, len, meth);
        family({copy_of(seg, {}, {}, 25, true), copy_of(seg, meth, three, -1, false), copy_of(seg, {}, four, -1, false)});
        reads.push_back(convert(seg, meth));
    }
    {   // three loci at 12, one at genomic 14 / mirror-free 11, one at 14
        const std::set<int> meth = {12};
        const std::string seg = segment(meth);
        std::vector<std::string> copies;
        for (int k = 0; k < 3; ++k) {
            const std::set<int> four = positions(4, 0, len, meth);
            copies.push_back(copy_of(seg, {}, four, -1, false));
        }
        copies.push_back(copy_of(seg, meth, {}, 30, true));
        const std::set<int> one = positions(1, 0, len, meth);
        copies.push_back(copy_of(seg, {}, one, 30, true));
        family(copies);
        reads.push_back(convert(seg, meth));
    }
    {   // a copy past the genomic cut, reached after the best
        const std::set<int> meth = {20};
        const std::string seg = segment(meth);
        const std::set<int> one = positions(1, 0, len, meth);
        family({seg, copy_of(seg, meth, one, -1, false)});
        reads.push_back(convert(seg, meth));
    }
    for (int f = 0; f < 24; ++f) {
        const int n_meth = 1 + rnd(3);
        const std::set<int> meth = positions(n_meth, 5, len - 5, {});
        const std::string seg = segment(meth);
        std::vector<std::string> copies;
        const int n_copies = 2 + rnd(6);
        for (int k = 0; k < n_copies; ++k) {
            std::set<int> mirror;
            for (int p : meth)
                if (rnd(3) == 0) mirror.insert(p);
            const int gap = rnd(10);
            const int n_real = rnd(5);
            const std::set<int> real = positions(n_real, 0, len, meth);
            const int gap_at = gap < 3 ? 15 + rnd(20) : -1;
            copies.push_back(copy_of(seg, mirror, real, gap_at, gap == 0));
        }
        family(copies);
        reads.push_back(convert(seg, meth));
        std::set<int> some;
        for (int p : meth)
            if (rnd(2)) some.insert(p);
        reads.push_back(convert(seg, some));
    }
    REQUIRE(aln_test::write_test_fasta(fa, contigs));
    std::vector<std::string> args = {"index", "--meth", fa};
    std::vector<char *> argv;
    for (auto &s : args) argv.push_back(&s[0]);
    argv.push_back(nullptr);
    optind = 1;
#if defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || defined(__NetBSD__)
    optreset = 1;
#endif
    REQUIRE(bwa_index((int)args.size(), argv.data()) == 0);
    FMI_search fmi((fa + ".meth").c_str());
    fmi.load_index(false, 1, true, BWAMEM_HUGETLB_OFF, false);
    auto *bwt = ba::bwt_from_fmi(&fmi);
    ba::aln_meth_ref_t mr;
    REQUIRE(ba::aln_meth_ref_load(&mr, fa.c_str(), fmi.idx->bns, 1) == 0);
    const auto fq = write_fq_in(tmp.path(), "blind-strata.fq", reads);
    const int n_reads = (int)reads.size();
    auto *opt = ba::gap_init_opt();
    opt->fnr = -1; opt->max_diff = 5;
    // the default penalties; M=1, O=3, E=1, where gaps are cheap against mismatches; and M=5, O=3,
    // E=1, where a mismatch costs more than a gap open
    long compared = 0, with_blind_rows = 0, with_moved_row = 0, with_other_main = 0;
    for (const auto &pen : {std::vector<int>{3, 11, 4}, std::vector<int>{1, 3, 1}, std::vector<int>{5, 3, 1}})
        for (int max_top2 : {30, 2})
            for (int seed_len : {1024, 20}) {
                CAPTURE(pen[0]); CAPTURE(max_top2); CAPTURE(seed_len);
                opt->s_mm = pen[0]; opt->s_gapo = pen[1]; opt->s_gape = pen[2];
                opt->max_top2 = max_top2; opt->seed_len = seed_len;
                // 0: genomic, conversion-blind; 1: collapsed; 2: genomic, scoring
                std::vector<ba::aln_meth_read_t> res[3];
                for (int mode = 0; mode < 3; ++mode) {
                    ba::aln_meth_opt_t mo;
                    mo.m.meth_mode = 1;
                    mo.m.meth_scoring = mode == 1 ? MEM_METH_SCORING_COLLAPSED : MEM_METH_SCORING_GENOMIC;
                    mo.m.meth_mapq = mode == 0 ? MEM_METH_MAPQ_CONVERSION_BLIND : MEM_METH_MAPQ_SCORING;
                    ba::aln_meth_opt_finish(&mo);
                    auto *input = ba::bwa_seq_open(fq.c_str());
                    int n = 0;
                    auto *seqs = ba::bwa_read_seq(input, n_reads, &n, BWA_MODE_COMPREAD, 0);
                    REQUIRE(n == n_reads);
                    res[mode].resize(n);
                    ba::aln_dup_state_t dup;
                    ba::aln_meth_search_batch(bwt, &mr, &mo, n, seqs, opt, 1, 1, 0, 1, 1000, &dup, res[mode].data());
                    ba::bwa_free_read_seq(n, seqs);
                    ba::bwa_seq_close(input);
                }
                for (int i = 0; i < n_reads; ++i) {
                    CAPTURE(i);
                    const ba::aln_meth_read_t &g = res[0][i], &c = res[1][i], &s = res[2][i];
                    // conversion-blind changes no genomic result: the same loci, counts and flags
                    // as --meth-mapq scoring
                    using Hit = std::tuple<int, int, int, bwtint_t, int32_t, int, int, int, int, int>;
                    auto hits = [](const ba::aln_meth_read_t &r) {
                        std::vector<Hit> v;
                        for (const auto &h : r.hits)
                            v.emplace_back(h.rid, h.hyp, h.strand, h.pos, h.cost, h.span, h.n_mm, h.n_gapo, h.n_gape, h.n_mirror);
                        return v;
                    };
                    CHECK(hits(g) == hits(s));
                    CHECK(g.c1 == s.c1);
                    CHECK(g.c2 == s.c2);
                    CHECK(g.flags == s.flags);
                    CHECK(g.main == s.main);
                    if (g.flags || c.flags || g.hits.empty() || c.hits.empty() || g.blind_incomplete) continue;
                    if (g.c1 > max_top2 || c.c1 > max_top2 || g.blind_c1 > max_top2) continue;
                    const bool gape = opt->mode & BWA_MODE_GAPE;
                    // each locus: (strand, position) -> (cost, mismatches, gap opens, extensions, span)
                    using Row = std::tuple<int64_t, int, int, int, int>;
                    std::map<std::pair<int, bwtint_t>, Row> cand, want;
                    for (const auto &h : g.hits) {
                        if (!ba::aln_meth_blind_selected(g, h, opt->s_mm, gape)) continue;
                        const ba::aln_meth_blind_row_t m = ba::aln_meth_blind_view(g, h, opt->s_mm);
                        cand[{m.strand, m.pos}] = Row{m.cost, m.n_mm, m.n_gapo, m.n_gape, m.span};
                        with_moved_row += m.cost != (int64_t)h.cost - (int64_t)opt->s_mm * h.n_mirror;
                    }
                    for (const auto &b : g.blind_rows) {
                        CHECK(cand.count({b.strand, b.pos}) == 0); // a hit is never a blind row
                        cand[{b.strand, b.pos}] = Row{b.cost, b.n_mm, b.n_gapo, b.n_gape, b.span};
                    }
                    const std::map<std::pair<int, bwtint_t>, Row> &blind = cand;
                    int64_t c_best = INT64_MAX;
                    for (const auto &h : c.hits) {
                        want[{h.strand, h.pos}] = Row{h.cost, h.n_mm, h.n_gapo, h.n_gape, h.span};
                        c_best = std::min(c_best, (int64_t)h.cost);
                    }
                    CHECK(g.blind_best == c_best);
                    CHECK(g.blind_c1 == c.c1);
                    CHECK(g.blind_c2 == c.c2);
                    CHECK(blind == want);
                    // collapsed's main hit is the conversion-blind draw: the loci at the best
                    // mirror-free cost, in collapsed's order, drawn with the read's seed
                    REQUIRE(g.blind_main >= 0);
                    CHECK(g.blind_main < c.c1);
                    const ba::aln_meth_blind_row_t &bm = g.blind_pick;
                    const ba::aln_meth_hit_t &cm = c.hits[c.main], &gm = g.hits[g.main];
                    CHECK(bm.rid == cm.rid);
                    CHECK((int)bm.hyp == (int)cm.hyp);
                    CHECK((int)bm.strand == (int)cm.strand);
                    CHECK(bm.pos == cm.pos);
                    with_other_main += bm.pos != gm.pos || bm.strand != gm.strand;
                    ++compared;
                    with_blind_rows += !g.blind_rows.empty();
                }
                for (auto &v : res) v.clear();
            }
    // the test compares most reads, some through rows only the mirror-free strata keep, some
    // through a hit whose lowest mirror-free row is not its genomic row, and some whose collapsed
    // main hit is not the genomic one
    CHECK(compared > 12 * 30);
    CHECK(with_blind_rows > 0);
    CHECK(with_moved_row > 0);
    CHECK(with_other_main > 0); // collapsed draws a main hit the genomic search does not
    free(opt);
    ba::aln_meth_ref_free(&mr);
    ba::bwt_destroy(bwt);
}

// --meth-mapq conversion-blind observes collapsed scoring: each mapped record has the lower of the
// scoring and collapsed MAPQs where collapsed places it at the same position and strand, else 0.
// The three runs are of the same input, so their records come in the same order.
static void check_blind_observes_collapsed(const std::vector<Rec> &blind, const std::vector<Rec> &scoring,
                                           const std::vector<Rec> &collapsed)
{
    REQUIRE(blind.size() == scoring.size());
    REQUIRE(blind.size() == collapsed.size());
    for (size_t i = 0; i < blind.size(); ++i) {
        const Rec &b = blind[i], &s = scoring[i], &c = collapsed[i];
        CAPTURE(b.name); CAPTURE(b.flag); CAPTURE(s.mapq); CAPTURE(c.rname); CAPTURE(c.pos); CAPTURE(c.mapq);
        REQUIRE(b.name == c.name);
        if (b.flag & 0x4) continue;
        const bool same = b.rname == c.rname && b.pos == c.pos && (b.flag & 0x10) == (c.flag & 0x10);
        CHECK(b.mapq == (same ? std::min(s.mapq, c.mapq) : 0));
    }
}

TEST_CASE("aln_meth: the conversion-blind pair MAPQ is collapsed's pairing's, capped at the genomic one")
{
    // Pairing's MAPQ depends on where the ends start: an end it moves gets at most the pair's MAPQ,
    // an end it leaves can get more (the sum of the ends', or the pair's plus 7 for an end of MAPQ
    // 0), and two moved ends get the pair's less 20. Collapsed scoring starts each end at its own
    // draw among the loci tied with the mirror cell (read 1 C over reference T, read 2 G over A)
    // free, which the genomic draw does not see, so the conversion-blind cap must pair from there.
    // Each pair has a home contig (read 1's segment, read 2's 300 bases on) and copies of a segment
    // on contigs of their own: exact (a genomic tie), or with a T under the read's kept C (one
    // mismatch genomically, none mirror-free). Sixteen copies of each pair draw with sixteen seeds.
    //   tie:  read 2 has an exact and a mirror copy (genomic MAPQ 0, collapsed 0). Where collapsed
    //         starts read 2 off home pairing moves it there: 29, where the genomic pair has 36.
    //   low:  read 2 has a mirror copy (genomic 23, collapsed 0): likewise 29 where the genomic has 60.
    //   both: each read has a mirror copy. Where collapsed starts both off home both move: 9, where
    //         the genomic pair has 46 (pairing from the genomic starts, both ends capped at 0, gave 0).
    // -A: no insert-size estimate (-a bounds the pairs), so no mate rescue either.
    bwa_tests::TempDir tmp{"aln_meth_blind_pair"};
    const std::string dir = tmp.path(), fa = tmp.unique_file();
    uint32_t st = 20261123;
    auto dna = [&](int n) { return aln_test::lcg_dna(st, n); };
    // read 1 keeps the C at offset 12 of its segment, read 2 the C at offset 20 of its segment
    // (read 2 reads the reverse strand: a kept G)
    std::vector<std::pair<std::string, std::string>> contigs;
    std::vector<std::string> r1, r2, names;
    auto add = [&](const std::string &name, bool dup2, bool mirror1, bool mirror2) {
        std::string s1 = dna(40), mid = dna(260), s2 = dna(40);
        s1[12] = 'C'; s2[20] = 'C';
        const std::string head = dna(100), tail = dna(100);
        contigs.emplace_back(name, head + s1 + mid + s2 + tail);
        auto copy = [&](const std::string &suffix, std::string seg, int t_at) {
            if (t_at >= 0) seg[t_at] = 'T';
            const std::string before = dna(100), after = dna(100);
            contigs.emplace_back(name + suffix, before + seg + after);
        };
        if (dup2) copy("_dup2", s2, -1);
        if (mirror1) copy("_mir1", s1, 12);
        if (mirror2) copy("_mir2", s2, 20);
        std::string a = conv(s1, 'C', 'T');
        a[12] = 'C';
        std::string b = conv(aln_test::revcomp(s2), 'G', 'A');
        b[19] = 'G'; // offset 20 of s2, read from the other end
        for (int k = 0; k < 16; ++k) {
            r1.push_back(a); r2.push_back(b);
            names.push_back(name + std::to_string(k));
        }
    };
    add("tie", true, false, true);
    add("low", false, false, true);
    add("both", false, true, true);
    REQUIRE(aln_test::write_test_fasta(fa, contigs));
    std::vector<std::string> index_args = {"index", "--meth", fa};
    std::vector<char *> argv;
    for (auto &arg : index_args) argv.push_back(&arg[0]);
    argv.push_back(nullptr);
    optind = 1;
#if defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || defined(__NetBSD__)
    optreset = 1;
#endif
    REQUIRE(bwa_index((int)index_args.size(), argv.data()) == 0);
    auto fq = [&](const std::string &file, const std::vector<std::string> &reads) {
        const std::string path = dir + "/" + file;
        std::ofstream out(path);
        for (size_t i = 0; i < reads.size(); ++i)
            out << '@' << names[i] << '\n' << reads[i] << "\n+\n" << std::string(reads[i].size(), 'I') << '\n';
        return path;
    };
    const std::string fq1 = fq("pair1.fq", r1), fq2 = fq("pair2.fq", r2);
    auto run = [&](std::vector<std::string> extra) {
        std::vector<std::string> args = {"--meth", "-A"};
        args.insert(args.end(), extra.begin(), extra.end());
        return run_aln(args, {fq1, fq2}, fa, dir);
    };
    const auto scoring = run({"--meth-mapq", "scoring"});
    const auto blind = run({});
    const auto collapsed = run({"--meth-scoring", "collapsed"});
    REQUIRE(blind.size() == 2 * names.size());
    for (size_t i = 0; i < blind.size(); ++i) { // placement is the genomic search's
        CAPTURE(blind[i].name);
        CHECK(blind[i].flag == scoring[i].flag);
        CHECK(blind[i].rname == scoring[i].rname);
        CHECK(blind[i].pos == scoring[i].pos);
        CHECK(blind[i].cigar == scoring[i].cigar);
    }
    check_blind_observes_collapsed(blind, scoring, collapsed);
    // each shape, as (pair, genomic MAPQ, conversion-blind MAPQ) of some record
    std::set<std::tuple<std::string, int, int>> seen;
    for (size_t i = 0; i < blind.size(); ++i)
        seen.emplace(blind[i].name.substr(0, blind[i].name.find_first_of("0123456789")), scoring[i].mapq, blind[i].mapq);
    CHECK(seen.count({"tie", 36, 29}) == 1);
    CHECK(seen.count({"low", 60, 29}) == 1);
    CHECK(seen.count({"both", 46, 9}) == 1);
}

TEST_CASE("aln_meth: the conversion-blind pair MAPQs are collapsed's, capped at the genomic ones, on random families")
{
    // A property test of the paired conversion-blind cap: on a seeded synthetic reference, each
    // record's conversion-blind MAPQ is the lower of its genomic and collapsed MAPQs where collapsed
    // places the end at the same position and strand, else 0 (check_blind_observes_collapsed), and
    // conversion-blind changes nothing else. Each family is a fragment on a home contig: read 1's
    // 40-base segment, a gap, read 2's segment. Its copies, on the same contig past a spacer, are
    // of either segment alone or of the whole fragment with a gap of another length (a suboptimal
    // pair), and carry mirror cells (a T under a C a read keeps), real mismatches or neither, as
    // can read 1's segment at home. Read 2's segment at home can carry two real mismatches, with a
    // mirror copy on a contig of its own (an improper pair under collapsed, which rescues read 2
    // home). Each family gives a few pairs, keeping some of the segments' Cs, some repeated (a
    // duplicate draws with its own seed). Plain pairs on contigs of their own give collapsed an
    // insert size too. Runs with an insert-size estimate
    // (collapsed's own) and mate rescue (collapsed's is rerun for the cap), and with -A. Every
    // random value is drawn into a local first, so the reference does not depend on the compiler's
    // argument evaluation order.
    bwa_tests::TempDir tmp{"aln_meth_blind_pairs"};
    const std::string dir = tmp.path(), fa = tmp.unique_file();
    uint32_t st = 3;
    const int len = 40;
    auto rnd = [&](int n) { st = st * 1103515245u + 12345u; return (int)((st >> 8) % (uint32_t)n); };
    auto alt = [](char c) { return c == 'G' ? 'A' : 'G'; }; // never T: cannot match a converted C
    // Distinct positions in [lo, hi) outside avoid.
    auto positions = [&](int n, int lo, int hi, const std::set<int> &avoid) {
        std::set<int> out;
        while ((int)out.size() < n) {
            const int p = lo + rnd(hi - lo);
            if (!avoid.count(p)) out.insert(p);
        }
        return out;
    };
    // seg with a T under each position of mirror and, at each of real, a base the read cannot match
    auto copy_of = [&](std::string seg, const std::set<int> &meth, const std::set<int> &mirror, const std::set<int> &real) {
        for (int p : mirror) seg[p] = 'T';
        for (int p : real) seg[p] = meth.count(p) ? 'A' : alt(seg[p] == 'C' ? 'T' : seg[p]);
        return seg;
    };
    std::vector<std::pair<std::string, std::string>> contigs;
    std::vector<std::string> r1, r2, names;
    for (int f = 0; f < 30; ++f) {
        // the Cs a read can keep: in read 1's segment (top strand), and in read 2's (a kept G)
        const int n_meth1 = 1 + rnd(2), n_meth2 = 1 + rnd(2);
        const std::set<int> meth1 = positions(n_meth1, 5, len - 5, {});
        const std::set<int> meth2 = positions(n_meth2, 5, len - 5, {});
        std::string s1 = aln_test::lcg_dna(st, len);
        std::string s2 = aln_test::lcg_dna(st, len);
        for (int p : meth1) s1[p] = 'C';
        for (int p : meth2) s2[p] = 'C';
        const int gap = 180 + rnd(120);
        const std::string mid = aln_test::lcg_dna(st, gap);
        // read 1's segment at home: as the read has it, with a T under a kept C, or with two real
        // mismatches (then a copy can be the genomic tie of home that is past collapsed's cut)
        const int home = rnd(3);
        std::string h1 = s1;
        if (home == 1) h1 = copy_of(s1, meth1, {*meth1.begin()}, {});
        else if (home == 2) {
            const std::set<int> two = positions(2, 0, len, meth1);
            h1 = copy_of(s1, meth1, {}, two);
        }
        // read 2's at home: as the read has it, or (in a third of the families) with two real
        // mismatches and a copy on a contig of its own with a T under each C read 2 can keep (its
        // collapsed best, with home past collapsed's cut: an improper pair, whose read 2 collapsed
        // rescues home)
        const int home2 = rnd(3);
        std::string h2 = s2;
        if (home2 == 2) {
            const std::set<int> two = positions(2, 0, len, meth2);
            h2 = copy_of(s2, meth2, {}, two);
            const std::string before = aln_test::lcg_dna(st, 100), after = aln_test::lcg_dna(st, 100);
            contigs.emplace_back("m" + std::to_string(f), before + copy_of(s2, meth2, meth2, {}) + after);
        }
        std::string c = aln_test::lcg_dna(st, 80);
        c += h1 + mid + h2;
        const int n_copies = 1 + rnd(5);
        for (int k = 0; k < n_copies; ++k) {
            const int spacer = 600 + rnd(200);
            c += aln_test::lcg_dna(st, spacer);
            const int kind = rnd(3); // 0: read 1's segment, 1: read 2's, 2: the fragment
            const int n_mirror1 = rnd(3), n_mirror2 = rnd(3), n_real1 = rnd(3), n_real2 = rnd(3);
            std::set<int> mirror1, mirror2;
            for (int p : meth1)
                if ((int)mirror1.size() < n_mirror1) mirror1.insert(p);
            for (int p : meth2)
                if ((int)mirror2.size() < n_mirror2) mirror2.insert(p);
            const std::set<int> real1 = positions(n_real1 == 2 ? 1 : 0, 0, len, meth1);
            const std::set<int> real2 = positions(n_real2 == 2 ? 1 : 0, 0, len, meth2);
            const std::string c1 = copy_of(s1, meth1, mirror1, real1), c2 = copy_of(s2, meth2, mirror2, real2);
            if (kind == 0) c += c1;
            else if (kind == 1) c += c2;
            else {
                const int other_gap = gap + 20 + rnd(60);
                c += c1 + aln_test::lcg_dna(st, other_gap) + c2;
            }
        }
        c += aln_test::lcg_dna(st, 80);
        const std::string name = "f" + std::to_string(f);
        contigs.emplace_back(name, c);
        const int n_pairs = 2 + rnd(3);
        for (int k = 0; k < n_pairs; ++k) {
            std::string a = conv(s1, 'C', 'T'), b = conv(aln_test::revcomp(s2), 'G', 'A');
            for (int p : meth1)
                if (rnd(4)) a[p] = 'C';
            for (int p : meth2)
                if (rnd(4)) b[len - 1 - p] = 'G';
            const int n_dup = 1 + (rnd(3) == 0 ? 2 : 0);
            for (int d = 0; d < n_dup; ++d) {
                r1.push_back(a); r2.push_back(b);
                names.push_back(name + "_" + std::to_string(k) + "_" + std::to_string(d));
            }
        }
    }
    // plain pairs on contigs of their own, at the families' inserts: an insert-size estimate for
    // collapsed (whose own selections are mostly repeats) as for the genomic run, so both rescue
    for (int k = 0; k < 24; ++k) {
        const int gap = 180 + rnd(120);
        const std::string s1 = aln_test::lcg_dna(st, len), mid = aln_test::lcg_dna(st, gap);
        const std::string s2 = aln_test::lcg_dna(st, len), tail = aln_test::lcg_dna(st, 80);
        const std::string head = aln_test::lcg_dna(st, 80);
        contigs.emplace_back("p" + std::to_string(k), head + s1 + mid + s2 + tail);
        r1.push_back(conv(s1, 'C', 'T')); r2.push_back(conv(aln_test::revcomp(s2), 'G', 'A'));
        names.push_back("p" + std::to_string(k));
    }
    REQUIRE(aln_test::write_test_fasta(fa, contigs));
    std::vector<std::string> index_args = {"index", "--meth", fa};
    std::vector<char *> argv;
    for (auto &arg : index_args) argv.push_back(&arg[0]);
    argv.push_back(nullptr);
    optind = 1;
#if defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || defined(__NetBSD__)
    optreset = 1;
#endif
    REQUIRE(bwa_index((int)index_args.size(), argv.data()) == 0);
    auto fq = [&](const std::string &file, const std::vector<std::string> &reads) {
        const std::string path = dir + "/" + file;
        std::ofstream out(path);
        for (size_t i = 0; i < reads.size(); ++i)
            out << '@' << names[i] << '\n' << reads[i] << "\n+\n" << std::string(reads[i].size(), 'I') << '\n';
        return path;
    };
    const std::string fq1 = fq("pairs1.fq", r1), fq2 = fq("pairs2.fq", r2);
    long moved = 0, lowered = 0, rescued = 0;
    // with the insert-size estimates (collapsed's and the genomic one) and mate rescue; -s drops the
    // rescue, -A both
    for (const std::vector<std::string> &extra :
         {std::vector<std::string>{}, std::vector<std::string>{"-s"}, std::vector<std::string>{"-A"}}) {
        CAPTURE(extra.size());
        auto run = [&](std::vector<std::string> mode) {
            std::vector<std::string> args = {"--meth"};
            args.insert(args.end(), extra.begin(), extra.end());
            args.insert(args.end(), mode.begin(), mode.end());
            return run_aln(args, {fq1, fq2}, fa, dir);
        };
        const auto scoring = run({"--meth-mapq", "scoring"});
        const auto blind = run({});
        const auto collapsed = run({"--meth-scoring", "collapsed"});
        REQUIRE(blind.size() == 2 * names.size());
        for (size_t i = 0; i < blind.size(); ++i) { // conversion-blind changes only MAPQ (and YQ:i)
            CAPTURE(blind[i].name);
            CHECK(blind[i].flag == scoring[i].flag);
            CHECK(blind[i].rname == scoring[i].rname);
            CHECK(blind[i].pos == scoring[i].pos);
            CHECK(blind[i].cigar == scoring[i].cigar);
            CHECK(blind[i].mapq <= scoring[i].mapq);
            lowered += blind[i].mapq < scoring[i].mapq;
            moved += blind[i].pos != collapsed[i].pos || blind[i].rname != collapsed[i].rname;
            const auto xt = collapsed[i].tag.find("XT:A");
            // collapsed rescues an end of the pair, lowering this end below the genomic MAPQ
            const size_t mate = i ^ 1;
            const auto xt_mate = collapsed[mate].tag.find("XT:A");
            const bool pair_rescued = (xt != collapsed[i].tag.end() && xt->second == "M")
                || (xt_mate != collapsed[mate].tag.end() && xt_mate->second == "M");
            rescued += pair_rescued && blind[i].mapq < scoring[i].mapq;
        }
        check_blind_observes_collapsed(blind, scoring, collapsed);
    }
    // the cap lowers some records, collapsed places some elsewhere, and collapsed's rescue of a pair
    // lowers some
    CHECK(lowered > 0);
    CHECK(moved > 0);
    CHECK(rescued > 0);
}

TEST_CASE("aln_meth: the conversion-blind pair MAPQ follows collapsed's mate rescue")
{
    // Collapsed scoring can leave a pair improper that the genomic run pairs, and then rescue a
    // mate: that lowers the anchor's MAPQ too (to at most its MAPQ less half the mate's, less 8),
    // so the conversion-blind cap must rerun collapsed's mate rescue. The pair's home contig has
    // read 1's segment and, 280 bases on, read 2's with two real mismatches; another contig has a
    // copy of read 2's segment with a T under each of the read's three kept Cs (three mismatches
    // genomically, in the seed, so -k drops it; none mirror-free). Genomically read 2 is home and
    // the pair proper (MAPQ 60). Collapsed places read 2 at the copy (home is past its cut), pairs
    // nothing, and rescues read 2 home from read 1: both ends get 11. A second pair, on contig plain,
    // has no kept C and read 2's segment there has four real mismatches (past -n): both runs leave
    // read 2 unmapped and rescue it, so the cap reruns collapsed's rescue of a pair it places one end
    // of. Pairs on contigs of their own, at inserts from 280 to 349, give both runs an insert size.
    // Also run with --exact-rescue-window and --short-insert, which the rerun rescue takes too.
    bwa_tests::TempDir tmp{"aln_meth_blind_rescue"};
    const std::string dir = tmp.path(), fa = tmp.unique_file();
    uint32_t st = 20261124;
    auto dna = [&](int n) { return aln_test::lcg_dna(st, n); };
    std::vector<std::pair<std::string, std::string>> contigs;
    std::vector<std::string> r1, r2, names;
    for (int k = 0; k < 24; ++k) {
        const std::string head = dna(100), s1 = dna(40), mid = dna(200 + 3 * k), s2 = dna(40), tail = dna(100);
        contigs.emplace_back("bg" + std::to_string(k), head + s1 + mid + s2 + tail);
        r1.push_back(conv(s1, 'C', 'T'));
        r2.push_back(conv(aln_test::revcomp(s2), 'G', 'A'));
        names.push_back("bg" + std::to_string(k));
    }
    std::string s1 = dna(40), s2 = dna(40);
    const int kept[3] = {10, 20, 30};
    for (int p : kept) s2[p] = 'C';
    // two real mismatches at home: a base that is neither the read's nor a conversion of it
    std::string h2 = s2;
    for (int p : {5, 35}) {
        const char read = s2[p] == 'C' ? 'T' : s2[p];
        h2[p] = read == 'G' ? 'A' : 'G';
    }
    std::string copy = s2;
    for (int p : kept) copy[p] = 'T';
    const std::string head = dna(100), mid = dna(240), tail = dna(100);
    contigs.emplace_back("home", head + s1 + mid + h2 + tail);
    const std::string before = dna(100), after = dna(100);
    contigs.emplace_back("copy", before + copy + after);
    std::string b = conv(aln_test::revcomp(s2), 'G', 'A');
    for (int p : kept) b[39 - p] = 'G';
    r1.push_back(conv(s1, 'C', 'T')); r2.push_back(b); names.push_back("target");
    {
        const std::string p1 = dna(40), p2 = dna(40);
        std::string m2 = p2; // four real mismatches, two at each end (local rescue clips them)
        for (int p : {1, 3, 36, 38}) {
            const char read = p2[p] == 'C' ? 'T' : p2[p];
            m2[p] = read == 'G' ? 'A' : 'G';
        }
        const std::string ph = dna(100), pm = dna(240), pt = dna(100);
        contigs.emplace_back("plain", ph + p1 + pm + m2 + pt);
        r1.push_back(conv(p1, 'C', 'T')); r2.push_back(conv(aln_test::revcomp(p2), 'G', 'A'));
        names.push_back("plain");
    }
    REQUIRE(aln_test::write_test_fasta(fa, contigs));
    std::vector<std::string> index_args = {"index", "--meth", fa};
    std::vector<char *> argv;
    for (auto &arg : index_args) argv.push_back(&arg[0]);
    argv.push_back(nullptr);
    optind = 1;
#if defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || defined(__NetBSD__)
    optreset = 1;
#endif
    REQUIRE(bwa_index((int)index_args.size(), argv.data()) == 0);
    auto fq = [&](const std::string &file, const std::vector<std::string> &reads) {
        const std::string path = dir + "/" + file;
        std::ofstream out(path);
        for (size_t i = 0; i < reads.size(); ++i)
            out << '@' << names[i] << '\n' << reads[i] << "\n+\n" << std::string(reads[i].size(), 'I') << '\n';
        return path;
    };
    const std::string fq1 = fq("rescue1.fq", r1), fq2 = fq("rescue2.fq", r2);
    auto xt = [](const Rec &r) { const auto it = r.tag.find("XT:A"); return it == r.tag.end() ? std::string() : it->second; };
    for (const std::vector<std::string> &extra :
         {std::vector<std::string>{}, std::vector<std::string>{"--exact-rescue-window", "--short-insert"}}) {
        CAPTURE(extra.size());
        auto run = [&](std::vector<std::string> mode) {
            std::vector<std::string> args = {"--meth"};
            args.insert(args.end(), extra.begin(), extra.end());
            args.insert(args.end(), mode.begin(), mode.end());
            return run_aln(args, {fq1, fq2}, fa, dir);
        };
        const auto scoring = run({"--meth-mapq", "scoring"});
        const auto blind = run({});
        const auto collapsed = run({"--meth-scoring", "collapsed"});
        REQUIRE(blind.size() == 2 * names.size());
        for (size_t i = 0; i < blind.size(); ++i) { // placement is the genomic search's
            CAPTURE(blind[i].name);
            CHECK(blind[i].flag == scoring[i].flag);
            CHECK(blind[i].rname == scoring[i].rname);
            CHECK(blind[i].pos == scoring[i].pos);
            CHECK(blind[i].cigar == scoring[i].cigar);
        }
        check_blind_observes_collapsed(blind, scoring, collapsed);
        // the target pair, as described: genomic proper at home; collapsed rescues read 2 home,
        // which lowers both ends below the genomic MAPQ
        const size_t t = blind.size() - 4;
        REQUIRE(blind[t].name == "target");
        REQUIRE(xt(collapsed[t + 1]) == "M");
        for (size_t e = t; e < t + 2; ++e) {
            CAPTURE(e - t);
            REQUIRE(scoring[e].rname == "home");
            REQUIRE((scoring[e].flag & 0x2) != 0);
            REQUIRE(collapsed[e].rname == "home");
            REQUIRE(collapsed[e].pos == scoring[e].pos);
            REQUIRE(collapsed[e].mapq < scoring[e].mapq);
            CHECK(blind[e].mapq == collapsed[e].mapq);
            CHECK(tag_value(blind[e], "YQ:i") == std::to_string(scoring[e].mapq));
        }
        if (extra.empty()) { // the shape described above
            CHECK(scoring[t].mapq == 60);
            CHECK(blind[t].mapq == 11);
            CHECK(blind[t + 1].mapq == 11);
        }
        // the plain pair: both runs rescue read 2, and the cap is collapsed's
        const size_t q = blind.size() - 2;
        REQUIRE(blind[q].name == "plain");
        REQUIRE(xt(scoring[q + 1]) == "M");
        REQUIRE(xt(collapsed[q + 1]) == "M");
        for (size_t e = q; e < q + 2; ++e) {
            CAPTURE(e - q);
            CHECK(collapsed[e].pos == scoring[e].pos);
            CHECK(blind[e].mapq == std::min(scoring[e].mapq, collapsed[e].mapq));
        }
    }
}

TEST_CASE("aln_meth: rows are scored along their search path, and -k limits the seed's strict differences"
          * doctest::test_suite("unit/aln_meth"))
{
    // w (40 bp, seeded under -l 20) has two loci. At k_over its reference holds a T under two read
    // Cs inside the seed: invisible to the converted search, two strict seed differences. At k_ok
    // it has three genuine mismatches past the seed. -k applies to strict differences, as -n does,
    // so under -k 1 only k_ok maps the read; under -k 2, or unseeded, k_over (cost 6) beats k_ok
    // (cost 9), which is then a retained competitor.
    uint32_t st = 4242;
    auto flank = [&]() { return aln_test::lcg_dna(st, 100); };
    std::string w = aln_test::lcg_dna(st, 40);
    for (int i : {3, 8}) w[i] = 'C';
    for (int i : {28, 32, 36}) w[i] = 'A';
    std::string k_over = w, k_ok = w;
    for (int i : {3, 8}) k_over[i] = 'T';
    for (int i : {28, 32, 36}) k_ok[i] = 'G';
    // g (40 bp): a one-base deletion at a locus otherwise exact; the row's alignment is the path's.
    std::string g = aln_test::lcg_dna(st, 40);
    g[19] = 'A'; g[20] = 'G'; g[21] = 'T';
    const std::string g_ref = g.substr(0, 20) + "C" + g.substr(20);
    bwa_tests::TempDir tmp{"aln_meth_strict_k"}; // this case's index, reads and output
    const std::string dir = tmp.path(), fa = dir + "/strict-k.fa";
    REQUIRE(aln_test::write_test_fasta(fa, {{"k_over", flank() + k_over + flank()},
                                            {"k_ok", flank() + k_ok + flank()},
                                            {"g_ref", flank() + g_ref + flank()}}));
    std::vector<std::string> index_args = {"index", "--meth", fa};
    std::vector<char *> argv;
    for (auto &arg : index_args) argv.push_back(&arg[0]);
    argv.push_back(nullptr);
    optind = 1;
#if defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || defined(__NetBSD__)
    optreset = 1;
#endif
    REQUIRE(bwa_index((int)index_args.size(), argv.data()) == 0);
    const std::string fq = write_fq_in(dir, "strict-k.fq", {w});
    const auto k1 = run_aln({"--meth", "-n", "3", "-l", "20", "-k", "1"}, {fq}, fa, dir);
    REQUIRE(k1.size() == 1);
    CHECK(k1[0].rname == "k_ok");
    CHECK(tag_value(k1[0], "NM:i") == "3");
    CHECK(tag_value(k1[0], "X0:i") == "1");
    CHECK(tag_value(k1[0], "X1:i") == "0");
    for (const auto &options : {std::vector<std::string>{"--meth", "-n", "3", "-l", "20", "-k", "2"},
                                std::vector<std::string>{"--meth", "-n", "3", "-l", "1024"}}) {
        const auto v = run_aln(options, {fq}, fa, dir);
        REQUIRE(v.size() == 1);
        CHECK(v[0].rname == "k_over");
        CHECK(tag_value(v[0], "NM:i") == "2");
        CHECK(tag_value(v[0], "X0:i") == "1");
        CHECK(tag_value(v[0], "X1:i") == "1");
    }
    const auto gv = run_aln({"--meth", "-n", "2", "-l", "1024"}, {write_fq_in(dir, "strict-k-g.fq", {g})}, fa, dir);
    REQUIRE(gv.size() == 1);
    CHECK(gv[0].rname == "g_ref");
    CHECK(gv[0].pos == 101);
    CHECK(gv[0].cigar == "20M1D20M"); // the deleted C sits between an A and a G: one placement
    CHECK(tag_value(gv[0], "NM:i") == "1");
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
            ba::aln_meth_occ_t o{};o.rid=0;o.hyp=1;o.strand=strand;
            o.pos=mr.bns->anns[0].offset+500;o.alignment.rb=o.pos;o.alignment.re=o.pos+57;
            o.alignment.cigar=strand ? std::vector<uint32_t>{3u<<4|4,57u<<4}
                                     : std::vector<uint32_t>{57u<<4,3u<<4|4};
            result.occ.push_back(std::move(o));
        }
        ba::aln_meth_store_loci(&result);
        ba::aln_meth_finish_se(&mr,&mo,opt,3,1,read,&result);
        auto span=[](const ba::bwa_cigar_t *c,int count) {
            int total=0;for(int k=0;k<count;++k) if(__cigar_op(c[k])!=2) total+=__cigar_len(c[k]);
            return total;
        };
        CHECK(span(read->cigar,read->n_cigar)==17060);
        REQUIRE(read->n_multi==1);
        CHECK(span(read->multi[0].cigar,read->multi[0].n_cigar)==17060);
        auto text=[](const ba::bwa_cigar_t *c,int count) {
            std::string s;
            for(int k=0;k<count;++k) s+=std::to_string(__cigar_len(c[k]))+"MIDS"[__cigar_op(c[k])];
            return s;
        };
        // The 17000 restored bases follow the read's 3' end: after the forward 3S, before the
        // reverse 3S. The primary is the first locus (strand `reverse`); XA holds the other.
        auto expected=[](int strand) { return std::string(strand ? "17003S57M" : "57M17003S"); };
        CHECK(text(read->cigar,read->n_cigar)==expected(reverse));
        CHECK(text(read->multi[0].cigar,read->multi[0].n_cigar)==expected(1-reverse));
        CHECK((int)read->len==17060);
        ba::bwa_free_read_seq(n,read);ba::bwa_seq_close(input);
    }
    free(opt);ba::aln_meth_ref_free(&mr);
}

TEST_CASE("aln_meth: a duplicate copies only the loci it keeps, and trim releases dropped loci")
{
    namespace ba=bwa_aln;
    ba::aln_meth_read_t src;src.main=0;src.c1=6;src.c2=4;src.hyp=1;src.n_lookups=99;
    // Every locus has 64 CIGAR operations (for the capacity checks), each list its own, so a copy
    // of the wrong locus's CIGAR shows.
    std::vector<std::vector<uint32_t>> input_cigar;
    for(int k=0;k<10;++k) {
        ba::aln_meth_occ_t o;o.rid=0;o.hyp=1;o.strand=k&1;o.pos=100*k;o.cost=k<6 ? 3 : 7;
        std::vector<uint32_t> cigar(64,50u<<4);cigar[k]=(uint32_t)(60+k)<<4;
        input_cigar.push_back(cigar);
        o.alignment.rb=o.pos;o.alignment.re=o.pos+50;o.alignment.cigar=std::move(cigar);
        src.occ.push_back(std::move(o));
    }
    // A finished representative: the search stored its loci as hits (64-op CIGARs in its table).
    ba::aln_meth_store_loci(&src);
    REQUIRE(src.hits.size()==10);REQUIRE(src.cig.size()==640);
    // An independent oracle for the draw: a seeded reservoir sample (splitmix64) over the six
    // best-cost loci (positions 0..500). The first tie is held without a draw; the k-th (k >= 2)
    // replaces the held one when the next draw % k is 0.
    auto oracle=[](uint64_t x) {
        int held=0;
        for(int k=2;k<=6;++k) {
            uint64_t z=(x+=0x9e3779b97f4a7c15ull);
            z=(z^(z>>30))*0xbf58476d1ce4e5b9ull;z=(z^(z>>27))*0x94d049bb133111ebull;z^=z>>31;
            if(z%(uint64_t)k==0) held=k-1;
        }
        return 100*held;
    };
    std::set<int> picked;
    for(uint64_t seed : {1ull,7ull,0x9e3779b97f4a7c15ull}) {
        // The reference behavior: copy everything, redraw, then trim.
        ba::aln_meth_read_t whole=src;ba::aln_meth_redraw_hits(&whole,seed);
        const auto drawn=whole.hits[whole.main].pos;
        CHECK(drawn==oracle(seed));
        picked.insert(oracle(seed));
        ba::aln_meth_trim(&whole,3);
        REQUIRE(whole.hits.size()==1);
        CHECK(whole.hits.capacity()==1);CHECK(whole.cig.size()==64);CHECK(whole.cig.capacity()==64);
        // Compares a locus of read x with whole's main hit, CIGAR included.
        auto same=[&whole](const ba::aln_meth_read_t &x,const ba::aln_meth_hit_t &a) {
            const ba::aln_meth_hit_t &b=whole.hits[0];
            return a.pos==b.pos&&a.strand==b.strand&&a.rid==b.rid&&a.hyp==b.hyp&&a.cost==b.cost
                &&cigar_of(x,a)==cigar_of(whole,b);
        };

        ba::aln_meth_read_t dup;dup.n_lookups=5;
        ba::aln_meth_copy_dup(&dup,&src,seed,0,3);
        REQUIRE(dup.hits.size()==1);
        CHECK(dup.hits.capacity()==1);CHECK(dup.cig.size()==64);CHECK(dup.cig.capacity()==64);
        CHECK(dup.main==0);CHECK(same(dup,dup.hits[0]));CHECK(dup.hits[0].pos==oracle(seed));
        CHECK(cigar_of(dup,dup.hits[0])==input_cigar[oracle(seed)/100]);
        CHECK(dup.c1==6);CHECK(dup.c2==4);CHECK(dup.hyp==1);CHECK(dup.flags==0);
        CHECK(dup.n_lookups==0);
        CHECK(src.hits.size()==10);CHECK(src.cig.size()==640);CHECK(src.n_lookups==99); // the representative is untouched

        // Keeping strata, or a list short enough to print as XA, copies every locus.
        ba::aln_meth_copy_dup(&dup,&src,seed,1,3);
        CHECK(dup.hits.size()==10);CHECK(dup.cig.size()==640);CHECK(same(dup,dup.hits[dup.main]));
        ba::aln_meth_copy_dup(&dup,&src,seed,0,9);
        CHECK(dup.hits.size()==10);CHECK(dup.cig.size()==640);CHECK(dup.hits[dup.main].pos==drawn);
        // A negative XA limit (--max-xa-hits -1) prints no XA: only the main hit, as aln_meth_trim.
        for(int limit : {-1,-5}) {
            CAPTURE(limit);
            ba::aln_meth_copy_dup(&dup,&src,seed,0,limit);
            REQUIRE(dup.hits.size()==1);CHECK(same(dup,dup.hits[0]));
            ba::aln_meth_read_t trimmed=src;ba::aln_meth_redraw_hits(&trimmed,seed);ba::aln_meth_trim(&trimmed,limit);
            REQUIRE(trimmed.hits.size()==1);CHECK(same(trimmed,trimmed.hits[0]));
        }
        // An incomplete (flagged) result keeps only its main hit however short the list.
        ba::aln_meth_read_t flagged=src;flagged.flags=ba::ALN_METH_BIND_BEST;
        ba::aln_meth_copy_dup(&dup,&flagged,seed,0,9);
        REQUIRE(dup.hits.size()==1);CHECK(same(dup,dup.hits[0]));
        CHECK(dup.flags==ba::ALN_METH_BIND_BEST);
    }
    CHECK(picked.size()>1); // the seeds select different loci, so a fixed choice cannot pass
    ba::aln_meth_read_t none;
    ba::aln_meth_read_t dup;ba::aln_meth_copy_dup(&dup,&none,1,0,3);
    CHECK(dup.main==-1);CHECK(dup.hits.empty());CHECK(dup.cig.empty());
}

TEST_CASE("aln_meth: duplicate reads under a negative XA limit keep only their main hit")
{
    namespace ba=bwa_aln;
    FMI_search fmi((ref().fa+".meth").c_str());fmi.load_index(false,1);
    ba::bwt_t *bwt=ba::bwt_from_fmi(&fmi);
    ba::aln_meth_ref_t mr;REQUIRE(ba::aln_meth_ref_load(&mr,ref().fa.c_str(),fmi.idx->bns,1)==0);
    ba::aln_meth_opt_t mo;ba::aln_meth_opt_finish(&mo);
    ba::gap_opt_t *opt=ba::gap_init_opt();
    // read_r has many equally good converted loci, so an untrimmed copy would keep all of them.
    const int n=8;
    const std::string fq=write_fq("negative-xa-dups.fq",std::vector<std::string>(n,ref().read_r));
    for(int limit : {-1,3}) {
        CAPTURE(limit);
        ba::bwa_seqio_t *input=ba::bwa_seq_open(fq.c_str());int got=0;
        ba::bwa_seq_t *reads=ba::bwa_read_seq(input,n,&got,BWA_MODE_COMPREAD,0);
        REQUIRE(got==n);
        std::vector<ba::aln_meth_read_t> batch(n);ba::aln_dup_state_t dup;
        REQUIRE(ba::aln_meth_search_batch(bwt,&mr,&mo,n,reads,opt,1,1,0,0,limit,&dup,batch.data())==n-1);
        // Every copy, like its representative, keeps only its main hit: there is no XA to print.
        for(int i=0;i<n;++i) {
            CAPTURE(i);
            REQUIRE(batch[i].main==0);CHECK(batch[i].hits.size()==1u);
            CHECK(batch[i].c1+batch[i].c2>1); // the trim dropped loci the copy would otherwise hold
        }
        ba::bwa_free_read_seq(n,reads);ba::bwa_seq_close(input);
    }
    free(opt);ba::aln_meth_ref_free(&mr);ba::bwt_destroy(bwt);
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
    ba::aln_meth_occ_t primary{};primary.rid=0;primary.hyp=1;primary.strand=0;
    primary.pos=c1_off+100;primary.alignment.rb=primary.pos;primary.alignment.re=primary.pos+400;
    primary.alignment.cigar={400u<<4};
    result.occ.push_back(std::move(primary));
    // An internal 256-base insertion: one gap open and 255 extensions fit a record, but the
    // alternative's edit distance (at least 256) does not fit XA's 8-bit field.
    ba::aln_meth_occ_t alt{};alt.rid=0;alt.hyp=1;alt.strand=0;alt.cost=1;
    alt.pos=c1_off+1000;alt.alignment.rb=alt.pos;alt.alignment.re=alt.pos+144;
    alt.alignment.cigar={72u<<4,256u<<4|1,72u<<4};
    result.occ.push_back(std::move(alt));
    ba::aln_meth_store_loci(&result);
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
            ba::aln_meth_occ_t o{};
            o.rid=0; o.hyp=1; o.strand=0; o.pos=100+k;
            o.cost=best+k;
            o.alignment.reward=100+k; // reward deliberately points in the opposite direction
            result.occ.push_back(std::move(o));
        }
    };
    // the retained loci as (pos, cost) pairs, in position order
    auto kept=[&]() {
        std::vector<std::pair<int64_t,int64_t>> v;
        for(const auto &o : result.hits) v.emplace_back(o.pos,(int64_t)o.cost);
        std::sort(v.begin(),v.end());
        return v;
    };
    auto first=[](int64_t best,int n) {
        std::vector<std::pair<int64_t,int64_t>> v;
        for(int k=0;k<n;++k) v.emplace_back(100+k,best+k);
        return v;
    };
    populate(0);ba::aln_meth_slot_finish(&slot);
    CHECK(result.c1==1); CHECK(result.c2==3);
    CHECK(kept()==first(0,4)); // costs 0..3: the window ends at best + s_mm, cost 4 is beyond it
    CHECK(result.hits[result.main].pos==100);
    opt->mode|=BWA_MODE_NONSTOP;
    populate(0);ba::aln_meth_slot_finish(&slot);
    CHECK(result.c1==1); CHECK(result.c2==4);
    CHECK(kept()==first(0,5));
    opt->mode&=~BWA_MODE_NONSTOP;opt->s_mm=INT_MAX;
    populate(INT_MAX-4);ba::aln_meth_slot_finish(&slot);
    CHECK(result.c1==1); CHECK(result.c2==4);
    CHECK(kept()==first(INT_MAX-4,5));
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
            ba::aln_meth_occ_t o{}; o.rid=0; o.hyp=1; o.strand=0; o.pos=100+k; o.cost=10+k;
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
    ba::aln_meth_occ_t o{}; o.rid=0; o.hyp=1; o.strand=0; o.pos=100; result.occ.push_back(std::move(o));
    ba::aln_meth_slot_finish(&slot);
    CHECK(result.flags==0);
    // Every candidate dropped (e.g. an oversized window was the only one): any of them could have
    // been the best, so the read is incomplete with nothing retained; with none dropped, complete.
    for(const auto &c : std::vector<std::pair<long,int>>{{1,ba::ALN_METH_BIND_UNREP},{0,0}}) {
        CAPTURE(c.first);
        result=ba::aln_meth_read_t();
        result.n_unrepresentable=c.first; slot.unrep_cost=c.first ? 0 : INT64_MAX;
        ba::aln_meth_slot_finish(&slot);
        CHECK(result.flags==c.second);
        CHECK(result.occ.empty()); CHECK(result.c1==0); CHECK(result.c2==0);
    }
    free(opt);
}

TEST_CASE("aln_meth: finished reads store compact loci and drop an incomplete read's worse loci"
          * doctest::test_suite("unit/aln_meth"))
{
    namespace ba = bwa_aln;
    ba::gap_opt_t *opt=ba::gap_init_opt();
    opt->s_mm=3;
    ba::aln_meth_group_t group; group.opt=opt; group.max_xa_hits=3;
    ba::aln_meth_read_t result;
    ba::aln_meth_slot_t slot; slot.grp=&group; slot.res=&result;
    // Three best-cost loci (pos 100, 102, 104), three worse ones, in a vector with spare capacity
    // as the search leaves it.
    auto populate=[&](int flags) {
        result=ba::aln_meth_read_t();
        result.occ.reserve(1024);
        const int64_t costs[6]={0,3,0,3,0,3};
        for(int k=0;k<6;++k) {
            ba::aln_meth_occ_t o;
            o.rid=0; o.hyp=1; o.strand=0; o.pos=100+k; o.cost=costs[k];
            o.alignment.rb=o.pos; o.alignment.re=o.pos+32;
            o.alignment.raw_cigar={32u<<4};
            o.alignment.cigar={32u<<4};
            result.occ.push_back(std::move(o));
        }
        result.flags=flags;
    };
    auto released=[&]() {
        // The search's loci are freed; the stored ones are sized exactly, one-op CIGARs inline.
        CHECK(result.occ.empty()); CHECK(result.occ.capacity()==0);
        CHECK(result.hits.capacity()==result.hits.size());
        CHECK(result.cig.empty());
        for(const auto &o:result.hits) CHECK(cigar_of(result,o)==std::vector<uint32_t>{32u<<4});
    };
    SUBCASE("paired, complete: every retained locus stays for pairing") {
        slot.keep=ba::ALN_METH_KEEP_STRATA; populate(0); ba::aln_meth_slot_finish(&slot);
        REQUIRE(result.hits.size()==6); CHECK(result.c1==3); CHECK(result.c2==3);
        for(int k=0;k<6;++k) CHECK(result.hits[k].pos==(bwtint_t)(100+k));
        released();
    }
    // Oracle for an incomplete read: the same read finished complete keeps all six loci; its
    // main-hit draw over them is what the incomplete read printed before anything was dropped.
    for (int keep : {(int)ba::ALN_METH_KEEP_COPIES, (int)ba::ALN_METH_KEEP_STRATA}) {
        CAPTURE(keep);
        SUBCASE("paired, incomplete: the same main hit; ties only while duplicate copies need them") {
            slot.keep=keep;
            std::set<bwtint_t> drawn;
            for(uint64_t seed=0;seed<16;++seed) {
                CAPTURE(seed);
                slot.rng=seed;
                populate(0); ba::aln_meth_slot_finish(&slot);
                REQUIRE(result.main>=0);
                const bwtint_t expect=result.hits[result.main].pos;
                populate(ba::ALN_METH_BIND_BEST); ba::aln_meth_slot_finish(&slot);
                CHECK(result.c1==3); CHECK(result.c2==3);
                REQUIRE(result.main>=0);
                CHECK(result.hits[result.main].pos==expect);
                CHECK(result.hits[result.main].cost==0);
                if (keep==ba::ALN_METH_KEEP_COPIES) {
                    // Copies redraw over the best-cost loci, so all of them stay, in order.
                    REQUIRE(result.hits.size()==3);
                    CHECK(result.hits[0].pos==100); CHECK(result.hits[1].pos==102); CHECK(result.hits[2].pos==104);
                    for(const auto &o:result.hits) CHECK(o.cost==0);
                } else {
                    // Kept only for pairing, which never reads an incomplete read past its main hit.
                    REQUIRE(result.hits.size()==1); CHECK(result.main==0);
                }
                drawn.insert(expect);
                released();
            }
            CHECK(drawn.size()==3); // the seeds reach every best-cost locus, so the check discriminates
        }
    }
    SUBCASE("single-end: the trim to the main hit also drops the search's capacity") {
        slot.keep=0; populate(0); ba::aln_meth_slot_finish(&slot);
        REQUIRE(result.hits.size()==1); CHECK(result.hits[0].cost==0);
        released();
    }
    SUBCASE("aln_meth_trim releases capacity on its own") {
        populate(0); result.main=0; ba::aln_meth_trim(&result,1); // 6 loci > max_xa_hits+1
        REQUIRE(result.hits.size()==1); released();
    }
    free(opt);
}

TEST_CASE("aln_meth: stored loci keep multi-op CIGARs in the read's table, through copies and trims"
          * doctest::test_suite("unit/aln_meth"))
{
    namespace ba = bwa_aln;
    const std::vector<std::vector<uint32_t>> cigars={{50u<<4},{20u<<4,1u<<4|1,29u<<4},{50u<<4},
        {3u<<4|4,47u<<4},{10u<<4,2u<<4|2,40u<<4}};
    ba::aln_meth_read_t result;
    for(int k=0;k<(int)cigars.size();++k) {
        ba::aln_meth_occ_t o;o.rid=k%2;o.hyp=k%2;o.strand=1-k%2;o.pos=1000+100*k;o.cost=3*k;
        o.n_mm=k;o.n_gapo=k%2;o.n_gape=2*k;o.n_conv=7+k;o.n_ins=k;o.n_del=k+1;
        o.alignment.rb=o.pos;o.alignment.re=o.pos+48+k;o.alignment.cigar=cigars[k];
        o.alignment.raw_cigar={9u<<4};
        result.occ.push_back(o);
    }
    result.main=3;
    ba::aln_meth_store_loci(&result);
    CHECK(result.occ.empty());
    REQUIRE(result.hits.size()==cigars.size());
    CHECK(result.cig.size()==3+2+3); // only the multi-op CIGARs are in the table
    for(int k=0;k<(int)cigars.size();++k) {
        CAPTURE(k);
        const auto &h=result.hits[k];
        CHECK(std::make_tuple((int)h.rid,(int)h.hyp,(int)h.strand,h.pos,(int)h.cost,(int)h.span)==
              std::make_tuple(k%2,k%2,1-k%2,(bwtint_t)(1000+100*k),3*k,48+k));
        CHECK(std::make_tuple((int)h.n_mm,(int)h.n_gapo,(int)h.n_gape,(int)h.n_conv)==
              std::make_tuple(k,k%2,2*k,7+k));
        CHECK(cigar_of(result,h)==cigars[k]);
    }
    ba::aln_meth_store_loci(&result); // no-op once finished
    CHECK(result.hits.size()==cigars.size());
    // A copy owns its own table.
    ba::aln_meth_read_t copy(result);
    for(int k=0;k<(int)cigars.size();++k) CHECK(cigar_of(copy,copy.hits[k])==cigars[k]);
    // Trimming to a multi-op main hit keeps exactly its CIGAR in a fresh table.
    ba::aln_meth_trim(&copy,1);
    REQUIRE(copy.hits.size()==1);CHECK(copy.main==0);CHECK(copy.hits.capacity()==1);
    CHECK(copy.hits[0].pos==1300);CHECK(cigar_of(copy,copy.hits[0])==cigars[3]);
    CHECK(copy.cig.size()==2);CHECK(copy.cig.capacity()==2);
    // ... and to a one-op main hit leaves no table at all.
    result.main=2;ba::aln_meth_trim(&result,1);
    REQUIRE(result.hits.size()==1);CHECK(result.hits[0].pos==1200);
    CHECK(cigar_of(result,result.hits[0])==cigars[2]);CHECK(result.cig.capacity()==0);
}

TEST_CASE("aln_meth: trimming an unstored read stores only its main locus"
          * doctest::test_suite("unit/aln_meth"))
{
    namespace ba = bwa_aln;
    const std::vector<std::vector<uint32_t>> cigars={{20u<<4,1u<<4|1,29u<<4},{},{10u<<4,2u<<4|2,40u<<4},{50u<<4}};
    for(int main_k=0;main_k<(int)cigars.size();++main_k) {
        CAPTURE(main_k);
        ba::aln_meth_read_t result;
        for(int k=0;k<(int)cigars.size();++k) {
            ba::aln_meth_occ_t o;o.rid=0;o.hyp=0;o.strand=0;o.pos=1000+100*k;o.cost=0;
            o.alignment.rb=o.pos;o.alignment.re=o.pos+50;o.alignment.cigar=cigars[k];
            result.occ.push_back(o);
        }
        result.main=main_k;
        ba::aln_meth_trim(&result,1); // 4 loci > max_xa_hits+1
        CHECK(result.occ.empty());
        REQUIRE(result.hits.size()==1);CHECK(result.hits.capacity()==1);CHECK(result.main==0);
        CHECK(result.hits[0].pos==(bwtint_t)(1000+100*main_k));
        CHECK(cigar_of(result,result.hits[0])==cigars[main_k]);
        const size_t table=cigars[main_k].size()>1?cigars[main_k].size():0;
        CHECK(result.cig.size()==table);CHECK(result.cig.capacity()==table);
        if(cigars[main_k].size()!=1) CHECK(result.hits[0].cig==0); // no stale table offset
    }
}

TEST_CASE("aln_meth: a locus that does not fit its stored record is fatal, not truncated"
          * doctest::test_suite("unit/aln_meth"))
{
    namespace ba = bwa_aln;
    bwa_tests::TempDir dir("aln_meth_fit");
    const std::string err_path=dir.path()+"/stderr";
    // Each case breaks one guard in meth_hit; case 0 fits.
    for(int bad=0;bad<13;++bad) {
        CAPTURE(bad);
        ba::aln_meth_read_t result;
        ba::aln_meth_occ_t o;o.rid=0;o.hyp=0;o.strand=0;o.pos=1000;o.cost=3;
        o.alignment.rb=o.pos;o.alignment.re=o.pos+50;o.alignment.cigar={50u<<4};
        switch(bad) {
            case 1: o.cost=(int64_t)INT_MAX+1; break;
            case 2: o.cost=-1; break;
            case 3: o.n_mm=256; break;
            case 4: o.n_mm=-1; break;
            case 5: o.n_gapo=256; break;
            case 6: o.n_gapo=-1; break;
            case 7: o.n_gape=256; break;
            case 8: o.n_gape=-1; break;
            case 9: o.alignment.re=o.alignment.rb+(int64_t)INT32_MAX+1; break;
            case 10: o.alignment.re=o.alignment.rb+(int64_t)INT32_MIN-1; break;
            case 11: o.alignment.rb=o.pos+1; o.alignment.re=o.alignment.rb+50; break; // pos must be the start
            case 12: o.alignment.rb=o.pos-1; o.alignment.re=o.alignment.rb+50; break;
        }
        result.occ.push_back(o);result.main=0;
        // The guard's own exit status reaches the parent: a forked child's exit() leaves the
        // shared index fixture alone (TempDir removes its directory only in the creating process).
        const int status=bwa_tests::run_in_child([&result] {
            ba::aln_meth_store_loci(&result);
            _exit(result.hits.size()==1 && result.hits[0].cost==3 ? 0 : 2);
        },err_path.c_str());
        std::ifstream in(err_path);
        const std::string err((std::istreambuf_iterator<char>(in)),std::istreambuf_iterator<char>());
        REQUIRE(WIFEXITED(status));
        if(bad==0) { CHECK(WEXITSTATUS(status)==0); CHECK(err.empty()); }
        else {
            CHECK(WEXITSTATUS(status)==1);
            CHECK(err.find("ERROR: --meth locus does not fit its stored record")!=std::string::npos);
        }
    }
}

TEST_CASE("aln_meth: a finished read's redraw draws what the search-time redraw would"
          * doctest::test_suite("unit/aln_meth"))
{
    namespace ba = bwa_aln;
    for(int moved : {0,1}) {
        CAPTURE(moved);
        ba::aln_meth_read_t search;
        const int costs[7]={3,0,0,6,0,0,0}, convs[7]={0,2,1,0,1,3,1};
        for(int k=0;k<7;++k) {
            ba::aln_meth_occ_t o;o.rid=0;o.hyp=1;o.strand=0;o.pos=100+k;o.cost=costs[k];o.n_conv=convs[k];
            o.alignment.rb=o.pos;o.alignment.re=o.pos+50;o.alignment.cigar={50u<<4};
            search.occ.push_back(o);
        }
        search.main=0;search.n_conv_moved=moved;
        ba::aln_meth_read_t stored(search);ba::aln_meth_store_loci(&stored);
        std::set<int> drawn;
        for(uint64_t seed=0;seed<64;++seed) {
            CAPTURE(seed);
            search.main=stored.main=0;
            ba::aln_meth_redraw(&search,seed);ba::aln_meth_redraw_hits(&stored,seed);
            CHECK(stored.main==search.main);
            drawn.insert(stored.main);
        }
        CHECK(drawn.size()==(moved?3u:5u)); // the draw ranges over every eligible tie
    }
}

TEST_CASE("aln_meth: capped duplicate copies draw the main hit a search of their own would"
          * doctest::test_suite("unit/aln_meth"))
{
    namespace ba=bwa_aln;
    FMI_search fmi((ref().fa+".meth").c_str()); fmi.load_index(false,1);
    auto *bwt=ba::bwt_from_fmi(&fmi);
    ba::aln_meth_ref_t mr;
    REQUIRE(ba::aln_meth_ref_load(&mr,ref().fa.c_str(),fmi.idx->bns,1)==0);
    auto *opt=ba::gap_init_opt(); opt->fnr=-1; opt->max_diff=3;
    opt->max_top2=1000; // so that stock's -R stop (more than 30 best loci) does not end it first
    // read_r has 31 equally good converted loci plus worse ones; a small lookup budget caps it.
    const int n=64;
    const auto fq=write_fq("capped-dups.fq",std::vector<std::string>(n,ref().read_r));
    // Budget 32 stops the search after all 31 ties and 1 worse locus. (The read's certified stop
    // completes it at 33 lookups, so a larger budget no longer caps it.)
    for (int keep_strata : {1,0}) for (int budget : {32}) {
        CAPTURE(keep_strata); CAPTURE(budget);
        ba::aln_meth_opt_t mo; mo.p_best=budget; mo.p_second=0;
        // The 31 ties are counted under collapsed scoring (the mirror cell is free too).
        mo.m.meth_scoring=MEM_METH_SCORING_COLLAPSED;
        ba::aln_meth_opt_finish(&mo);
        auto *input=ba::bwa_seq_open(fq.c_str()); int got=0;
        auto *reads=ba::bwa_read_seq(input,n,&got,BWA_MODE_COMPREAD,0);
        REQUIRE(got==n);
        std::vector<ba::aln_meth_read_t> batch(n);
        ba::aln_dup_state_t dup;
        const int n_dups=ba::aln_meth_search_batch(bwt,&mr,&mo,n,reads,opt,2,1,0,keep_strata,3,&dup,batch.data());
        REQUIRE(n_dups==n-1); // every copy took the duplicate path
        REQUIRE(batch[0].flags!=0);
        REQUIRE(batch[0].c1==31); REQUIRE(batch[0].c2>=1); // the draw matters and worse loci were dropped
        CHECK(batch[0].hits.size()==1u); // an incomplete read keeps only its main hit once its copies are drawn
        std::set<bwtint_t> drawn;
        for (int i=0;i<n;++i) {
            CAPTURE(i);
            // Independent reference: the copy searched on its own, at its own ordinal.
            ba::aln_meth_read_t single; ba::aln_dup_state_t single_dup;
            ba::aln_meth_search_batch(bwt,&mr,&mo,1,reads+i,opt,1,1,i,keep_strata,3,&single_dup,&single);
            REQUIRE(batch[i].main>=0); REQUIRE(single.main>=0);
            // The whole stored record: placement, counts, span and CIGAR.
            auto same=[&](const ba::aln_meth_hit_t &x,const ba::aln_meth_hit_t &y) {
                CHECK(std::make_tuple(x.rid,x.hyp,x.strand,x.pos,x.cost,x.span)==
                      std::make_tuple(y.rid,y.hyp,y.strand,y.pos,y.cost,y.span));
                CHECK(std::make_tuple(x.n_mm,x.n_gapo,x.n_gape,x.n_conv)==
                      std::make_tuple(y.n_mm,y.n_gapo,y.n_gape,y.n_conv));
                CHECK(cigar_of(batch[i],x)==cigar_of(single,y));
            };
            const auto &a=batch[i].hits[batch[i].main], &b=single.hits[single.main];
            same(a,b);
            CHECK(batch[i].c1==single.c1); CHECK(batch[i].c2==single.c2);
            REQUIRE(batch[i].hits.size()==single.hits.size());
            CHECK(batch[i].main==single.main);
            for (size_t k=0;k<batch[i].hits.size();++k) {
                const auto &x=batch[i].hits[k], &y=single.hits[k];
                same(x,y);
                CHECK(x.cost==a.cost);
            }
            drawn.insert(a.pos);
        }
        CHECK(drawn.size()>=2); // copies draw different loci, so a shared main would fail
        ba::bwa_free_read_seq(n,reads); ba::bwa_seq_close(input);
    }
    {
        // A complete paired-end copy is neither trimmed nor shrunk, so only a fresh copy frees the
        // buffer a reused result carries over from the previous batch.
        ba::aln_meth_opt_t mo; ba::aln_meth_opt_finish(&mo);
        auto *input=ba::bwa_seq_open(fq.c_str()); int got=0;
        auto *reads=ba::bwa_read_seq(input,n,&got,BWA_MODE_COMPREAD,0);
        REQUIRE(got==n);
        std::vector<ba::aln_meth_read_t> batch(n);
        for (auto &b:batch) b.hits.reserve(256);
        ba::aln_dup_state_t dup;
        REQUIRE(ba::aln_meth_search_batch(bwt,&mr,&mo,n,reads,opt,2,1,0,1,3,&dup,batch.data())==n-1);
        REQUIRE(batch[0].flags==0); REQUIRE(batch[0].hits.size()>1);
        for (int i=1;i<n;++i) {
            CAPTURE(i);
            CHECK(batch[i].hits.size()==batch[0].hits.size());
            CHECK(batch[i].hits.capacity()==batch[i].hits.size());
        }
        ba::bwa_free_read_seq(n,reads); ba::bwa_seq_close(input);
    }
    free(opt); ba::aln_meth_ref_free(&mr); ba::bwt_destroy(bwt);
}

TEST_CASE("aln_meth: identical reads with different --known-snvs candidates do not share a search")
{
    namespace ba=bwa_aln;
    FMI_search fmi((ref().fa+".meth").c_str()); fmi.load_index(false,1);
    auto *bwt=ba::bwt_from_fmi(&fmi);
    ba::aln_meth_ref_t mr;
    REQUIRE(ba::aln_meth_ref_load(&mr,ref().fa.c_str(),fmi.idx->bns,1)==0);
    auto *opt=ba::gap_init_opt(); opt->fnr=-1; opt->max_diff=2;
    // A site at c1[q] (A or T) with ALT G. The read, read 1 unmethylated, is c1[q - 10, q + 15) with
    // G at the site and two A<->G mismatches: three differences, over the budget of 2, so its own
    // search finds nothing; the site's window as a candidate places it with G (two).
    int q=1005;
    while (ref().c1[q]!='A' && ref().c1[q]!='T') ++q;
    std::string read=ref().c1.substr(q-10,25);
    read[10]='G';
    int n_mm=0;
    for (int k : {3,17,6,20,2,22,4,18})
        if (n_mm<2 && (read[k]=='A' || read[k]=='G')) read[k]=read[k]=='A'?'G':'A', ++n_mm;
    REQUIRE(n_mm==2);
    read=conv(read,'C','T');
    ba::ks_sites_t sites;
    sites.by_contig.assign(mr.bns->n_seqs,{});
    sites.by_contig[0]={{q,(uint8_t)(ref().c1[q]=='A'?0:3),2}};
    sites.n_sites=1;
    mr.known_snvs=&sites;
    const int n=8;
    const auto fq=write_fq("ks-dups.fq",std::vector<std::string>(n,read));
    // Only the odd copies have the window (as an incomplete window search gives one copy and not another).
    const int64_t rb=mr.bns->anns[0].offset+q-10;
    std::vector<std::vector<ba::aln_meth_candidate_t>> extra(n);
    for (int i=1;i<n;i+=2) extra[i]={{0,0,1,rb,rb+25}};
    ba::aln_meth_opt_t mo; ba::aln_meth_opt_finish(&mo);
    auto *input=ba::bwa_seq_open(fq.c_str()); int got=0;
    auto *reads=ba::bwa_read_seq(input,n,&got,BWA_MODE_COMPREAD,0);
    REQUIRE(got==n);
    std::vector<ba::aln_meth_read_t> batch(n);
    ba::aln_dup_state_t dup;
    // every copy but the first of each candidate set shares its search
    CHECK(ba::aln_meth_search_batch(bwt,&mr,&mo,n,reads,opt,1,1,0,0,3,&dup,batch.data(),extra.data())==n-2);
    for (int i=0;i<n;++i) {
        CAPTURE(i);
        ba::aln_meth_read_t single; ba::aln_dup_state_t single_dup;
        ba::aln_meth_search_batch(bwt,&mr,&mo,1,reads+i,opt,1,1,i,0,3,&single_dup,&single,&extra[i]);
        REQUIRE((single.main>=0)==(i%2==1)); // the test's premise: only the candidate places it
        REQUIRE(batch[i].main==single.main);
        CHECK(batch[i].c1==single.c1); CHECK(batch[i].c2==single.c2);
        REQUIRE(batch[i].hits.size()==single.hits.size());
        for (size_t k=0;k<single.hits.size();++k) {
            const auto &x=batch[i].hits[k], &y=single.hits[k];
            CHECK(std::make_tuple(x.rid,x.hyp,x.strand,x.pos,x.cost,x.span,x.alt,
                                  x.n_cig,x.n_conv,x.n_mm,x.n_gapo,x.n_gape)==
                  std::make_tuple(y.rid,y.hyp,y.strand,y.pos,y.cost,y.span,y.alt,
                                  y.n_cig,y.n_conv,y.n_mm,y.n_gapo,y.n_gape));
            if (x.n_cig==y.n_cig) { // the operations, not the cig field (an offset into r.cig)
                const uint32_t *cx=ba::aln_meth_hit_cigar(batch[i],x), *cy=ba::aln_meth_hit_cigar(single,y);
                CHECK(std::vector<uint32_t>(cx,cx+x.n_cig)==std::vector<uint32_t>(cy,cy+y.n_cig));
            }
        }
    }
    ba::bwa_free_read_seq(n,reads); ba::bwa_seq_close(input);
    free(opt); ba::aln_meth_ref_free(&mr); ba::bwt_destroy(bwt);
}

TEST_CASE("aln_meth: the -R stop counts only best loci within the final difference ceiling")
{
    // M=1, O=3, E=1: 31 best loci (cost 3) with three strict mismatches each, and one with a single
    // gap open (cost 3, one difference). The best row's differences set the ceiling to 2, which drops
    // the 31 at the finish, so they must not count toward -R: one best locus is not a repeat.
    namespace ba = bwa_aln;
    ba::gap_opt_t *opt = ba::gap_init_opt();
    opt->s_mm = 1; opt->s_gapo = 3; opt->s_gape = 1; opt->max_diff = 3; opt->max_top2 = 30;
    ba::aln_meth_opt_t mo; ba::aln_meth_opt_finish(&mo);
    ba::aln_meth_ref_t mr;
    ba::aln_meth_group_t group; group.ref = &mr; group.mo = &mo; group.opt = opt;
    ba::aln_meth_slot_t slot; slot.grp = &group;
    ba::gap_job_t job = {}; job.opt = *opt; job.max_diff = 3; job.score_cut = INT_MAX;
    slot.job = &job;
    ba::aln_meth_read_t result; slot.res = &result;
    for (int i = 0; i < 31; ++i) {
        ba::aln_meth_occ_t o{}; o.cost = 3; o.n_mm = 3; o.pos = 1000 * (i + 1);
        result.occ.push_back(o);
    }
    ba::aln_meth_occ_t gapped{}; gapped.cost = 3; gapped.n_gapo = 1; gapped.pos = 50000;
    result.occ.push_back(gapped);
    slot.best = 3;
    // While the next entry is at the best cost the stratum is incomplete: no ceiling, no -R stop.
    CHECK(ba::aln_meth_boundary(&slot, 3) == 0);
    CHECK(job.max_diff == 3);
    CHECK(ba::aln_meth_boundary(&slot, 4) == 0);
    CHECK(job.max_diff == 2);
    // 31 best loci within the ceiling are a repeat: 30 more one-gap rows (cost 3, one difference)
    for (int i = 0; i < 30; ++i) {
        ba::aln_meth_occ_t o{}; o.cost = 3; o.n_mm = 0; o.n_gapo = 1; o.pos = 100000 + 1000 * i;
        result.occ.push_back(o);
    }
    CHECK(ba::aln_meth_boundary(&slot, 4) == 1);
    free(opt);
}

TEST_CASE("aln_meth: the difference ceiling comes from the stored best rows")
{
    // M=1, O=3, E=1, -n 4. The best locus holds a three-mismatch row at cost 3 (a one-gap row at
    // the same locus and cost may have been resolved and replaced by it). Only stored rows set the
    // ceiling, as at the finish: 3 + 1 = 4, so the search keeps paths a four-mismatch row needs.
    namespace ba = bwa_aln;
    ba::gap_opt_t *opt = ba::gap_init_opt();
    opt->s_mm = 1; opt->s_gapo = 3; opt->s_gape = 1; opt->max_diff = 4;
    ba::aln_meth_opt_t mo; ba::aln_meth_opt_finish(&mo);
    ba::aln_meth_ref_t mr;
    ba::aln_meth_group_t group; group.ref = &mr; group.mo = &mo; group.opt = opt;
    ba::aln_meth_slot_t slot; slot.grp = &group;
    ba::gap_job_t job = {}; job.opt = *opt; job.max_diff = 4; job.score_cut = INT_MAX;
    slot.job = &job;
    ba::aln_meth_read_t result; slot.res = &result;
    ba::aln_meth_occ_t o{}; o.cost = 3; o.n_mm = 3; o.pos = 1000;
    result.occ.push_back(o);
    slot.best = 3;
    CHECK(ba::aln_meth_boundary(&slot, 4) == 0);
    CHECK(job.max_diff == 4);
    free(opt);
}

TEST_CASE("aln_meth: under conversion-blind the search ceiling is stock's, and a pruned mirror-free row is noted")
{
    // Default penalties (M=3, O=11, E=4), -n 5. The best locus is one gap open (cost 11, one
    // difference), so stock's ceiling is 2, and conversion-blind keeps it: the search is that of
    // --meth-mapq scoring. With the mirror cell free a locus can cost less than 11 with more
    // differences (X: three real mismatches and one mirror, 12, mirror-free 9 with three), and the
    // mirror-free strata then keep rows of up to four differences up to 9 + 3: at next 12 the
    // ceiling can prune one, so the strata are unknown (blind_need). They are complete with X1 in
    // place of X (a gap open and a mirror mismatch: 14, mirror-free 11 with one difference, so
    // the mirror-free ceiling is 2), at next 13 (above 9 + 3), or without X.
    namespace ba = bwa_aln;
    struct Case { int blind, x, next, need; }; // x: 0 none, 1 X, 2 X1
    for (const Case &k : {Case{0, 1, 12, 0}, Case{1, 1, 12, 1}, Case{1, 2, 12, 0}, Case{1, 1, 13, 0}, Case{1, 0, 12, 0}}) {
        CAPTURE(k.blind); CAPTURE(k.x); CAPTURE(k.next);
        ba::gap_opt_t *opt = ba::gap_init_opt();
        opt->max_diff = 5;
        ba::aln_meth_opt_t mo; ba::aln_meth_opt_finish(&mo);
        ba::aln_meth_ref_t mr;
        ba::aln_meth_group_t group; group.ref = &mr; group.mo = &mo; group.opt = opt;
        ba::aln_meth_slot_t slot; slot.grp = &group; slot.blind = k.blind; slot.max_diff = 5;
        ba::gap_job_t job = {}; job.opt = *opt; job.max_diff = 5; job.score_cut = INT_MAX;
        slot.job = &job;
        ba::aln_meth_read_t result; slot.res = &result;
        ba::aln_meth_occ_t gapped{}; gapped.cost = 11; gapped.n_gapo = 1; gapped.pos = 1000;
        gapped.alignment.rb = 1000;
        result.occ.push_back(gapped);
        if (k.x) {
            ba::aln_meth_occ_t x{}; x.n_mirror = 1; x.pos = 2000; x.alignment.rb = 2000;
            if (k.x == 1) { x.cost = 12; x.n_mm = 4; }
            else { x.cost = 14; x.n_mm = 1; x.n_gapo = 1; }
            result.occ.push_back(x);
        }
        slot.best = 11;
        CHECK(ba::aln_meth_boundary(&slot, k.next) == 0);
        CHECK(job.max_diff == 2);
        CHECK(slot.blind_need == k.need);
        CHECK(slot.need_next == (k.blind ? k.next : -1));
        // decided once, where the cuts first apply
        CHECK(ba::aln_meth_boundary(&slot, k.next + 1) == 0);
        CHECK(slot.blind_need == k.need);
        free(opt);
    }
}

TEST_CASE("aln_meth: under conversion-blind the -R stop is stock's, and notes unknown mirror-free strata")
{
    // Default penalties, -n 5, -R 30: 31 best loci at cost 12 (four mismatches) are a repeat, and
    // the search stops at next 13 with or without conversion-blind. A locus with one gap open, one
    // mismatch and that mismatch a mirror costs 14 but 11 with the mirror cell free, so the
    // mirror-free best is 11 and unique: its strata run to 14 and the entries from 13 on were not
    // expanded, so they are unknown (blind_need). Without it the mirror-free best loci are the 31,
    // over -R themselves, so a mirror-free run stops there too. Stopping at a later next (15, past
    // 11 + 3) leaves nothing unexpanded.
    namespace ba = bwa_aln;
    for (int blind = 0; blind < 2; ++blind)
        for (int mirror_row = 0; mirror_row < 2; ++mirror_row)
            for (int next : {13, 15}) {
                CAPTURE(blind); CAPTURE(mirror_row); CAPTURE(next);
                ba::gap_opt_t *opt = ba::gap_init_opt();
                opt->max_diff = 5; opt->max_top2 = 30;
                ba::aln_meth_opt_t mo; ba::aln_meth_opt_finish(&mo);
                ba::aln_meth_ref_t mr;
                ba::aln_meth_group_t group; group.ref = &mr; group.mo = &mo; group.opt = opt;
                ba::aln_meth_slot_t slot; slot.grp = &group; slot.blind = blind; slot.max_diff = 5;
                ba::gap_job_t job = {}; job.opt = *opt; job.max_diff = 5; job.score_cut = INT_MAX;
                slot.job = &job;
                ba::aln_meth_read_t result; slot.res = &result;
                for (int i = 0; i < 31; ++i) {
                    ba::aln_meth_occ_t o{}; o.cost = 12; o.n_mm = 4; o.pos = 1000 * (i + 1);
                    o.alignment.rb = (int64_t)o.pos;
                    result.occ.push_back(o);
                }
                if (mirror_row) {
                    ba::aln_meth_occ_t o{}; o.cost = 14; o.n_mm = 1; o.n_mirror = 1; o.n_gapo = 1; o.pos = 50000;
                    o.alignment.rb = 50000;
                    result.occ.push_back(o);
                }
                slot.best = 12;
                CHECK(ba::aln_meth_boundary(&slot, next) == 1);
                CHECK(slot.blind_need == (blind && mirror_row && next == 13));
                free(opt);
            }
}

TEST_CASE("aln_meth: the finish applies the difference ceiling only to a search that takes the cuts")
{
    // M=1, O=3: a one-gap best (cost 3, one difference) and a four-mismatch row (cost 4, within
    // best + M). Under the cuts the ceiling (1 + 1) drops the second; a search without them
    // (--known-snvs, -N, a negative penalty) keeps the cost window alone.
    namespace ba = bwa_aln;
    ba::gap_opt_t *opt = ba::gap_init_opt();
    opt->s_mm = 1; opt->s_gapo = 3; opt->s_gape = 1;
    ba::aln_meth_group_t group; group.opt = opt;
    ba::aln_meth_read_t result;
    ba::aln_meth_slot_t slot; slot.grp = &group; slot.res = &result; slot.keep = 1;
    for (int cut = 0; cut < 2; ++cut) {
        CAPTURE(cut);
        result = ba::aln_meth_read_t();
        ba::aln_meth_occ_t gapped{}; gapped.rid = 0; gapped.hyp = 1; gapped.pos = 100; gapped.cost = 3; gapped.n_gapo = 1;
        ba::aln_meth_occ_t mm{}; mm.rid = 0; mm.hyp = 1; mm.pos = 200; mm.cost = 4; mm.n_mm = 4;
        result.occ.push_back(gapped); result.occ.push_back(mm);
        slot.cut = cut;
        ba::aln_meth_slot_finish(&slot);
        CHECK(result.c1 == 1);
        CHECK(result.c2 == (cut ? 0 : 1));
        // the retained loci, as stored for pairing and output: the gapped best, then (without the
        // cuts) the four-mismatch row
        CHECK(result.occ.empty());
        REQUIRE(result.hits.size() == (cut ? 1u : 2u));
        CHECK(result.main == 0);
        for (const auto &h : result.hits) CHECK((h.rid == 0 && h.hyp == 1 && h.strand == 0));
        CHECK(result.hits[0].pos == 100);
        CHECK(result.hits[0].cost == 3);
        CHECK(result.hits[0].n_mm == 0);
        CHECK(result.hits[0].n_gapo == 1);
        if (!cut) {
            CHECK(result.hits[1].pos == 200);
            CHECK(result.hits[1].cost == 4);
            CHECK(result.hits[1].n_mm == 4);
            CHECK(result.hits[1].n_gapo == 0);
        }
    }
    free(opt);
}

namespace {
// A finish-only slot (no search) for the conversion-blind finish tests, default penalties.
struct BlindFinish {
    bwa_aln::gap_opt_t *opt = bwa_aln::gap_init_opt();
    bwa_aln::aln_meth_opt_t mo;
    bwa_aln::aln_meth_group_t group;
    bwa_aln::aln_meth_read_t result;
    bwa_aln::aln_meth_slot_t slot;
    explicit BlindFinish(int cut)
    {
        mo.m.meth_mode = 1; mo.m.meth_mapq = MEM_METH_MAPQ_CONVERSION_BLIND;
        bwa_aln::aln_meth_opt_finish(&mo);
        group.opt = opt; group.mo = &mo;
        slot.grp = &group; slot.res = &result; slot.keep = 1; slot.cut = cut; slot.blind = 1;
    }
    ~BlindFinish() { free(opt); }
    void add(bwtint_t pos, int64_t cost, int n_mm, int n_mirror, int n_gapo)
    {
        bwa_aln::aln_meth_occ_t o{};
        o.rid = 0; o.hyp = 1; o.pos = pos; o.cost = cost; o.n_mm = n_mm; o.n_mirror = n_mirror; o.n_gapo = n_gapo;
        o.alignment.rb = (int64_t)pos; o.alignment.re = (int64_t)pos + 40;
        result.occ.push_back(o);
    }
    std::vector<std::pair<bwtint_t, int32_t>> blind_rows() const
    {
        std::vector<std::pair<bwtint_t, int32_t>> v;
        for (const auto &b : result.blind_rows) v.emplace_back(b.pos, b.cost);
        std::sort(v.begin(), v.end());
        return v;
    }
};
} // namespace

TEST_CASE("aln_meth: a locus the finish's ceiling drops is a conversion-blind row")
{
    // Default penalties. G: one gap open (cost 11, one difference). X: four mismatches, one a
    // mirror (12, mirror-free 9 with three). X2: four mismatches (12). Under the cuts stock keeps
    // G alone (ceiling 2); the mirror-free strata are X (best), then G and X2 (within 9 + 3, and
    // four differences within X's 3 + 1). X and X2 are within the genomic cost window but over
    // its ceiling, so they are neither hits nor counted genomically: the shadow pairing can see
    // them only as blind rows. Without the cuts all three are hits.
    for (int cut = 0; cut < 2; ++cut) {
        CAPTURE(cut);
        BlindFinish f(cut);
        f.add(100, 11, 0, 0, 1);
        f.add(200, 12, 4, 1, 0);
        f.add(300, 12, 4, 0, 0);
        bwa_aln::aln_meth_slot_finish(&f.slot);
        CHECK(f.result.c1 == 1);
        CHECK(f.result.c2 == (cut ? 0 : 2));
        CHECK(f.result.hits.size() == (cut ? 1u : 3u));
        CHECK(f.result.blind_best == 9);
        CHECK(f.result.blind_c1 == 1);
        CHECK(f.result.blind_c2 == 2);
        using Rows = std::vector<std::pair<bwtint_t, int32_t>>;
        CHECK(f.blind_rows() == (cut ? Rows{{200, 9}, {300, 12}} : Rows{}));
    }
}

TEST_CASE("aln_meth: under the cuts the conversion-blind strata take the mirror-free ceiling")
{
    // Default penalties. G: one gap open (cost 11, one difference); Y: four mismatches (12). Under
    // the cuts stock's ceiling (2) drops Y, and so does collapsed's, which has the same rows: the
    // mirror-free strata are stock's selection on mirror-free costs, so Y is not counted there
    // either. Without the cuts both keep the cost window alone and count Y.
    for (int cut = 0; cut < 2; ++cut) {
        CAPTURE(cut);
        BlindFinish f(cut);
        f.add(100, 11, 0, 0, 1);
        f.add(200, 12, 4, 0, 0);
        bwa_aln::aln_meth_slot_finish(&f.slot);
        CHECK(f.result.c1 == 1);
        CHECK(f.result.c2 == (cut ? 0 : 1));
        CHECK(f.result.blind_best == 11);
        CHECK(f.result.blind_c1 == 1);
        CHECK(f.result.blind_c2 == (cut ? 0 : 1));
        CHECK(f.result.blind_rows.empty());
    }
}

TEST_CASE("aln_meth: under the cuts the conversion-blind strata are collapsed's selection")
{
    // A property test of the conversion-blind cuts: on a seeded synthetic reference, the strata a
    // genomic search builds with the mirror cell freed (blind_best, blind_c1, blind_c2, and the
    // loci: the hits within them plus blind_rows) must be what collapsed scoring, which frees the
    // mirror cell itself, selects (its best, c1, c2 and hits). Each family is a 50-base segment
    // whose read 1 keeps some Cs; its copies carry mirror cells (a T under a kept C), real
    // mismatches and one-base gaps. Planted families: the gapped-best triple (a one-gap best, a
    // mirror-free best with three mismatches and a mirror, a four-mismatch locus only the
    // mirror-free ceiling keeps), the -R case (three best loci at 12, a gapped mirror locus at
    // mirror-free 11, and a locus at 14 past the stop), and a copy past the genomic cut that the
    // search reaches after the best. Reads where either side stops on -R are skipped: past such a
    // stop the second stratum is partial, as in stock aln. So are reads whose strata the cuts may
    // have pruned (blind_incomplete, cap 0), which the test counts.
    // The loci compared are the ones the paired-end shadow pairing consumes: the hits it keeps
    // (aln_meth_blind_selected) and blind_rows. conversion-blind must leave the genomic search as
    // --meth-mapq scoring runs it: loci, counts, flags, the main hit, and the search's own work
    // (lookups, flushes, rows scored), including where -m binds.
    namespace ba = bwa_aln;
    bwa_tests::TempDir tmp{"aln_meth_blind_property"};
    // a unique basename: the index is looked up by basename, which the directory does not isolate
    const std::string fa = tmp.unique_file();
    uint32_t st = 20261123;
    const int len = 50;
    auto rnd = [&](int n) { st = st * 1103515245u + 12345u; return (int)((st >> 8) % (uint32_t)n); };
    auto alt = [](char c) { return c == 'G' ? 'A' : 'G'; }; // never T: cannot match a converted C
    auto convert = [](std::string s, const std::set<int> &keep) {
        for (int i = 0; i < (int)s.size(); ++i)
            if (s[i] == 'C' && !keep.count(i)) s[i] = 'T';
        return s;
    };
    // A copy of seg: a T under each position of mirror, a mismatch against the read at each of
    // real (read is seg converted, so the substitute is never T and never C), then one gap at
    // gap_at (del: the copy lacks a base, so the read has an insertion; ins: the copy has one more).
    auto copy_of = [&](std::string seg, const std::set<int> &mirror, const std::set<int> &real, int gap_at, bool del) {
        const std::string read = convert(seg, {});
        for (int p : mirror) seg[p] = 'T';
        for (int p : real) seg[p] = alt(read[p]);
        if (gap_at >= 0) {
            if (del) seg.erase(gap_at, 1);
            else seg.insert(gap_at, 1, alt(seg[gap_at]));
        }
        return seg;
    };
    // A segment with Cs at the positions of meth (the methylated Cs a read can keep).
    auto segment = [&](const std::set<int> &meth) {
        std::string s = aln_test::lcg_dna(st, len);
        for (int p : meth) s[p] = 'C';
        return s;
    };
    // Distinct positions in [lo, hi) outside avoid.
    auto positions = [&](int n, int lo, int hi, const std::set<int> &avoid) {
        std::set<int> out;
        while ((int)out.size() < n) {
            const int p = lo + rnd(hi - lo);
            if (!avoid.count(p)) out.insert(p);
        }
        return out;
    };
    std::vector<std::pair<std::string, std::string>> contigs;
    std::vector<std::string> reads;
    auto family = [&](const std::vector<std::string> &copies) {
        std::string c = aln_test::lcg_dna(st, 80);
        for (const auto &k : copies) c += k + aln_test::lcg_dna(st, 60 + rnd(60));
        contigs.emplace_back("f" + std::to_string(contigs.size()), c);
    };
    {   // the gapped-best triple
        const std::set<int> meth = {12};
        const std::string seg = segment(meth);
        family({copy_of(seg, {}, {}, 25, true), copy_of(seg, meth, positions(3, 0, len, meth), -1, false),
                copy_of(seg, {}, positions(4, 0, len, meth), -1, false)});
        reads.push_back(convert(seg, meth));
    }
    {   // -R: three loci at 12, one at genomic 14 / mirror-free 11, one at 14 past the stop
        const std::set<int> meth = {12};
        const std::string seg = segment(meth);
        std::vector<std::string> copies;
        for (int k = 0; k < 3; ++k) copies.push_back(copy_of(seg, {}, positions(4, 0, len, meth), -1, false));
        copies.push_back(copy_of(seg, meth, {}, 30, true));
        copies.push_back(copy_of(seg, {}, positions(1, 0, len, meth), 30, true));
        family(copies);
        reads.push_back(convert(seg, meth));
    }
    {   // a copy past the genomic cut, reached after the best
        const std::set<int> meth = {20};
        const std::string seg = segment(meth);
        family({seg, copy_of(seg, meth, positions(1, 0, len, meth), -1, false)});
        reads.push_back(convert(seg, meth));
    }
    for (int f = 0; f < 80; ++f) {
        const std::set<int> meth = positions(1 + rnd(3), 5, len - 5, {});
        const std::string seg = segment(meth);
        std::vector<std::string> copies;
        const int n_copies = 2 + rnd(6);
        for (int k = 0; k < n_copies; ++k) {
            std::set<int> mirror;
            for (int p : meth)
                if (rnd(3) == 0) mirror.insert(p);
            // every random value into a local first: the reference must not depend on the
            // compiler's argument evaluation order
            const int gap = rnd(10);
            const int n_real = rnd(5);
            const std::set<int> real = positions(n_real, 0, len, meth);
            const int gap_at = gap < 3 ? 15 + rnd(20) : -1;
            copies.push_back(copy_of(seg, mirror, real, gap_at, gap == 0));
        }
        family(copies);
        reads.push_back(convert(seg, meth));
        std::set<int> some;
        for (int p : meth)
            if (rnd(2)) some.insert(p);
        reads.push_back(convert(seg, some));
    }
    REQUIRE(aln_test::write_test_fasta(fa, contigs));
    std::vector<std::string> args = {"index", "--meth", fa};
    std::vector<char *> argv;
    for (auto &s : args) argv.push_back(&s[0]);
    argv.push_back(nullptr);
    optind = 1;
#if defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || defined(__NetBSD__)
    optreset = 1;
#endif
    REQUIRE(bwa_index((int)args.size(), argv.data()) == 0);
    FMI_search fmi((fa + ".meth").c_str());
    fmi.load_index(false, 1, true, BWAMEM_HUGETLB_OFF, false);
    auto *bwt = ba::bwt_from_fmi(&fmi);
    ba::aln_meth_ref_t mr;
    REQUIRE(ba::aln_meth_ref_load(&mr, fa.c_str(), fmi.idx->bns, 1) == 0);
    const auto fq = write_fq_in(tmp.path(), "blind-property.fq", reads);
    const int n_reads = (int)reads.size();
    auto *opt = ba::gap_init_opt();
    opt->fnr = -1; opt->max_diff = 5;
    // the planted cases' arithmetic is at the default penalties; the random families also run at
    // M=1, O=3, E=1, where gaps are cheap against mismatches, at M=5, O=3, E=1, where a gap open
    // costs less than a mismatch, and at -e 3 -E 1 (extensions not counted as differences); each
    // also (once per penalty set) with a -m queue limit that binds on some reads
    using Key = std::tuple<int, bwtint_t, int64_t, int, int, int>; // strand, position, mirror-free cost, counts
    long compared = 0, with_blind_rows = 0, incomplete = 0, se_zero = 0, bound = 0;
    const int max_entries0 = opt->max_entries, max_gape0 = opt->max_gape, mode0 = opt->mode;
    for (const auto &pen : {std::vector<int>{3, 11, 4, 0}, std::vector<int>{1, 3, 1, 0}, std::vector<int>{5, 3, 1, 0},
                            std::vector<int>{3, 11, 1, 3}})
        for (int max_top2 : {30, 2})
            for (int seed_len : {1024, 20})
            for (int max_entries : {max_entries0, 200}) {
                // the binding -m once per penalty set, to keep the test short
                if (max_entries != max_entries0 && (max_top2 != 30 || seed_len != 1024)) continue;
                CAPTURE(pen[0]); CAPTURE(pen[2]); CAPTURE(pen[3]); CAPTURE(max_top2); CAPTURE(seed_len); CAPTURE(max_entries);
                opt->s_mm = pen[0]; opt->s_gapo = pen[1]; opt->s_gape = pen[2];
                opt->max_gape = pen[3] ? pen[3] : max_gape0;
                opt->mode = pen[3] ? (mode0 & ~BWA_MODE_GAPE) : mode0;
                opt->max_top2 = max_top2; opt->seed_len = seed_len; opt->max_entries = max_entries;
                // 0: genomic, conversion-blind; 1: collapsed; 2: genomic, scoring
                std::vector<ba::aln_meth_read_t> res[3];
                for (int mode = 0; mode < 3; ++mode) {
                    ba::aln_meth_opt_t mo;
                    mo.m.meth_mode = 1;
                    mo.m.meth_scoring = mode == 1 ? MEM_METH_SCORING_COLLAPSED : MEM_METH_SCORING_GENOMIC;
                    mo.m.meth_mapq = mode == 0 ? MEM_METH_MAPQ_CONVERSION_BLIND : MEM_METH_MAPQ_SCORING;
                    ba::aln_meth_opt_finish(&mo);
                    auto *input = ba::bwa_seq_open(fq.c_str());
                    int n = 0;
                    auto *seqs = ba::bwa_read_seq(input, n_reads, &n, BWA_MODE_COMPREAD, 0);
                    REQUIRE(n == n_reads);
                    res[mode].resize(n);
                    ba::aln_dup_state_t dup;
                    ba::aln_meth_search_batch(bwt, &mr, &mo, n, seqs, opt, 1, 1, 0, 1, 1000, &dup, res[mode].data());
                    ba::bwa_free_read_seq(n, seqs);
                    ba::bwa_seq_close(input);
                }
                for (int i = 0; i < n_reads; ++i) {
                    CAPTURE(i);
                    const ba::aln_meth_read_t &g = res[0][i], &c = res[1][i], &s = res[2][i];
                    // conversion-blind changes no genomic result: the same loci, counts and flags
                    // as --meth-mapq scoring
                    using Hit = std::tuple<int, int, int, bwtint_t, int32_t, int, int, int, int, int>;
                    auto hits = [](const ba::aln_meth_read_t &r) {
                        std::vector<Hit> v;
                        for (const auto &h : r.hits)
                            v.emplace_back(h.rid, h.hyp, h.strand, h.pos, h.cost, h.span, h.n_mm, h.n_gapo, h.n_gape, h.n_mirror);
                        return v;
                    };
                    CHECK(hits(g) == hits(s));
                    CHECK(g.c1 == s.c1);
                    CHECK(g.c2 == s.c2);
                    CHECK(g.flags == s.flags);
                    CHECK(g.main == s.main);
                    CHECK(g.n_lookups == s.n_lookups);
                    CHECK(g.n_flushes == s.n_flushes);
                    CHECK(g.n_rescored == s.n_rescored);
                    bound += (s.flags & ba::ALN_METH_BIND_M) != 0;
                    if (g.flags || c.flags || g.hits.empty() || c.hits.empty()) continue;
                    const bool gape = opt->mode & BWA_MODE_GAPE;
                    if (g.blind_incomplete) {
                        // Counted: whether the single-end MAPQ is 0 without that rule too (the
                        // genomic one, more than one best locus; or the cap, the main hit not
                        // alone at the mirror-free best within its ceiling). Not guaranteed: the
                        // main locus's mirror-free row can be a moved row with more differences.
                        ++incomplete;
                        const ba::aln_meth_hit_t &h = g.hits[g.main];
                        const ba::aln_meth_blind_row_t m = ba::aln_meth_blind_view(g, h, opt->s_mm);
                        se_zero += g.c1 > 1 || !ba::aln_meth_blind_selected(g, h, opt->s_mm, gape)
                            || m.cost > g.blind_best || g.blind_c1 > 1;
                        continue;
                    }
                    if (c.c1 > max_top2 || g.blind_c1 > max_top2) continue;
                    // the loci the shadow pairing consumes
                    std::set<Key> blind, want;
                    for (const auto &h : g.hits) {
                        if (!ba::aln_meth_blind_selected(g, h, opt->s_mm, gape)) continue;
                        const ba::aln_meth_blind_row_t m = ba::aln_meth_blind_view(g, h, opt->s_mm);
                        blind.emplace(m.strand, m.pos, m.cost, m.n_mm, m.n_gapo, m.n_gape);
                    }
                    for (const auto &b : g.blind_rows) {
                        CHECK(blind.count(Key{b.strand, b.pos, b.cost, b.n_mm, b.n_gapo, b.n_gape}) == 0);
                        blind.emplace(b.strand, b.pos, b.cost, b.n_mm, b.n_gapo, b.n_gape);
                    }
                    int64_t c_best = INT64_MAX;
                    for (const auto &h : c.hits) {
                        want.emplace(h.strand, h.pos, h.cost, h.n_mm, h.n_gapo, h.n_gape);
                        c_best = std::min(c_best, (int64_t)h.cost);
                    }
                    CHECK(g.blind_best == c_best);
                    CHECK(g.blind_c1 == c.c1);
                    CHECK(g.blind_c2 == c.c2);
                    CHECK(blind == want);
                    ++compared;
                    with_blind_rows += !g.blind_rows.empty();
                }
                for (auto &v : res) v.clear();
            }
    // the test compares most reads, some of them through rows only the mirror-free strata keep,
    // and -m binds on some
    CHECK(compared > 24 * 100);
    CHECK(with_blind_rows > 0);
    CHECK(bound > 0);
    CHECK(incomplete > 0); // the need test runs end to end
    MESSAGE("compared " << compared << ", incomplete " << incomplete << " (single-end MAPQ 0 regardless: "
            << se_zero << "), -m bound " << bound);
    free(opt);
    ba::aln_meth_ref_free(&mr);
    ba::bwt_destroy(bwt);
}

TEST_CASE("aln_meth: the flush that binds the lookup budget resolves only a few more rows")
{
    namespace ba = bwa_aln;
    const std::string fq = write_fq("bind.fq", {conv(ref().c1.substr(500,50),'C','T')});
    ba::bwa_seqio_t *input = ba::bwa_seq_open(fq.c_str());
    REQUIRE(input != nullptr);
    int n=0;ba::bwa_seq_t *reads=ba::bwa_read_seq(input,1,&n,BWA_MODE_COMPREAD,0);
    REQUIRE(n==1);REQUIRE(reads!=nullptr);
    const std::string prefix=ref().fa+".meth";
    FMI_search fmi(prefix.c_str());fmi.load_index(false,1);
    ba::bwt_t *bwt=ba::bwt_from_fmi(&fmi);
    ba::aln_meth_ref_t mr;
    REQUIRE(ba::aln_meth_ref_load(&mr,ref().fa.c_str(),fmi.idx->bns,1)==0);
    ba::gap_opt_t *opt=ba::gap_init_opt();
    ba::aln_meth_opt_t mo;
    ba::aln_meth_group_t group;group.bwt=bwt;group.ref=&mr;group.mo=&mo;group.opt=opt;
    group.n_slots=1;group.slot.resize(1);
    ba::gap_slot_core_t core={};group.cores=&core;
    auto &slot=group.slot[0];slot.grp=&group;slot.job=&core.job;
    ba::aln_meth_read_t result;
    const ba::gap_entry_t exact={};
    // Each flush holds the intervals [k, k + width) of its widths, consecutive over the SA rows;
    // the rows' windows are whatever they are (the lookup count does not depend on them).
    REQUIRE(bwt->seq_len>4000);
    auto run=[&](int budget,const std::vector<std::vector<int>> &flushes) {
        mo.p_best=budget;mo.p_second=0;ba::aln_meth_opt_finish(&mo);
        ba::aln_meth_slot_start(&slot,reads,&result,0,3,0);
        std::vector<long> lookups;bwtint_t k=0;
        for (const auto &widths : flushes) {
            for (int w : widths) { ba::aln_meth_hit(&slot,&exact,k,k+w-1); k+=w; }
            core.state=ba::GAP_SLOT_PARKED;ba::aln_meth_flush_parked(&group);
            lookups.push_back(result.n_lookups);
        }
        return lookups;
    };
    SUBCASE("flushes within the budget resolve every row") {
        CHECK(run(1000,{{100,200},{300}})==std::vector<long>{300,600});
        CHECK(result.flags==0);
    }
    SUBCASE("a flush that fills the budget exactly does not bind it") {
        CHECK(run(1000,{{400},{600}})==std::vector<long>{400,1000});
        CHECK(result.flags==0);
    }
    SUBCASE("the binding flush resolves ALN_METH_BIND_ROWS more rows and marks the read") {
        CHECK(run(1000,{{100},{50,2000,10}})==std::vector<long>{100,100+ba::ALN_METH_BIND_ROWS});
        CHECK((result.flags&ba::ALN_METH_BIND_BEST)!=0);
    }
    SUBCASE("the binding rows are spread over the pending intervals in order") {
        // The windows of the rows that were looked up, and the rows that were a wrong copy:
        // together they tell which SA rows the flush resolved.
        using seen_t=std::pair<decltype(slot.windows.key),long>;
        // as sets: a flush may resolve a batch's rows in any order
        auto same=[](const seen_t &a,const seen_t &b) {
            return a.second==b.second && a.first.size()==b.first.size() &&
                   std::is_permutation(a.first.begin(),a.first.end(),b.first.begin());
        };
        auto rows_of=[&](int budget,const std::vector<std::pair<bwtint_t,int>> &intervals) {
            mo.p_best=budget;mo.p_second=0;ba::aln_meth_opt_finish(&mo);
            ba::aln_meth_slot_start(&slot,reads,&result,0,3,0);
            for (const auto &iv : intervals) ba::aln_meth_hit(&slot,&exact,iv.first,iv.first+iv.second-1);
            core.state=ba::GAP_SLOT_PARKED;ba::aln_meth_flush_parked(&group);
            return seen_t{slot.windows.key,result.n_wrong_copy};
        };
        // The oracles resolve their rows in full, within the budget: the first 40 rows plus 24
        // of the second interval (discovery order), and the third interval plus 24 of the second
        // (the reverse).
        const seen_t in_order=rows_of(1000,{{0,40},{40,24}}), reversed=rows_of(1000,{{80,40},{40,24}});
        REQUIRE(!same(in_order,reversed)); // the windows tell the two orders apart
        CHECK(run(100,{{40,40,40}})==std::vector<long>{ba::ALN_METH_BIND_ROWS});
        CHECK((result.flags&ba::ALN_METH_BIND_BEST)!=0);
        CHECK(same(seen_t{slot.windows.key,result.n_wrong_copy},in_order));
    }
    SUBCASE("an interval two paths reach counts once toward the budget") {
        // A one-base insertion's path (gap word: S index 25, length 1) to the same interval as the
        // exact one: its rows reuse the first path's SA rows, so 100 new rows fit the 150 left.
        mo.p_best=1000;mo.p_second=0;ba::aln_meth_opt_finish(&mo);
        ba::aln_meth_slot_start(&slot,reads,&result,0,3,0);
        ba::aln_meth_hit(&slot,&exact,0,849);
        core.state=ba::GAP_SLOT_PARKED;ba::aln_meth_flush_parked(&group);
        REQUIRE(result.n_lookups==850);
        ba::gap_entry_t ins={};ins.n_gapo=1;ins.n_ins=1;ins.last_diff_pos=25|(1u<<16);
        ba::aln_meth_hit(&slot,&exact,850,949);ba::aln_meth_hit(&slot,&ins,850,949);
        core.state=ba::GAP_SLOT_PARKED;ba::aln_meth_flush_parked(&group);
        CHECK(result.n_lookups==950);
        CHECK(result.flags==0);
    }
    SUBCASE("the binding flush never resolves past the budget") {
        CHECK(run(1000,{{980},{30}})==std::vector<long>{980,1000});
        CHECK((result.flags&ba::ALN_METH_BIND_BEST)!=0);
    }
    SUBCASE("a binding first flush") {
        CHECK(run(1000,{{1001}})==std::vector<long>{ba::ALN_METH_BIND_ROWS});
        CHECK((result.flags&ba::ALN_METH_BIND_BEST)!=0);
    }
    ba::aln_meth_slot_finish(&slot);
    ba::bwa_free_read_seq(n,reads);ba::bwa_seq_close(input);
    ba::aln_meth_ref_free(&mr);ba::bwt_destroy(bwt);free(opt);
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
    // A path with a deletion before the read's first base (strand 0: the window's first base), as
    // its gap event records it (ALN_METH_GAP: the complemented read's last index, 49).
    ba::gap_entry_t exact={},longer={};longer.n_del=1;longer.n_gapo=1;longer.last_diff_pos=(int)ALN_METH_GAP(49,1,1);
    for(int duplicate=0;duplicate<601;++duplicate) ba::aln_meth_hit(&slot,&exact,row,row);
    ba::aln_meth_hit(&slot,&longer,row,row);
    REQUIRE(slot.pend.size()==2);
    core.state=ba::GAP_SLOT_PARKED;ba::aln_meth_flush_parked(&group);
    CHECK(result.n_rescored==2);CHECK(slot.windows.size()==2);
    // The two spans have the same normalized start; their minimum-cost alignment wins.
    REQUIRE(result.occ.size()==1);
    CHECK(result.occ[0].pos==500);CHECK(result.occ[0].strand==0); // the row's own locus
    CHECK(result.occ[0].cost==0);
    CHECK(result.occ[0].n_mm==0);
    CHECK(result.occ[0].n_gapo==0);
    CHECK(result.occ[0].n_gape==0);
    CHECK(result.occ[0].n_conv==(int)std::count(ref().c1.begin()+500,ref().c1.begin()+550,'C'));
    ba::aln_meth_hit(&slot,&exact,row,row); // repeated interval across flushes
    CHECK(slot.pend.empty());
    ba::aln_meth_slot_finish(&slot);
    REQUIRE(result.hits.size()==1);CHECK(result.c1==1);
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
        // The read is converted at every C, so under cpg its non-CpG conversions are strict
        // mismatches, more than -n 3 allows: a limit of its length keeps the window to inspect.
        ba::aln_meth_slot_start(&slot,reads,&result,0,50,1);
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
    // An interval reached again along a path with other gaps reuses its SA rows: no lookup, so it
    // cannot bind the budget.
    ba::aln_meth_slot_start(&slot,reads,&result,2,3,1);
    ba::aln_meth_hit(&slot,&longer,row,row);
    ba::aln_meth_hit(&slot,&exact,row,row);
    core.state=ba::GAP_SLOT_PARKED;
    ba::aln_meth_flush_parked(&group);
    CHECK(result.n_lookups==1);CHECK(result.flags==0);
    REQUIRE(result.occ.size()==1);CHECK(result.occ[0].cost==0);
    ba::aln_meth_slot_finish(&slot);
    ba::aln_meth_slot_start(&slot,reads,&result,2,3,1);
    ba::aln_meth_hit(&slot,&longer,row,row); // the costlier path evaluated first
    ba::aln_meth_hit(&slot,&exact,row+1,row+1); // another interval hidden behind the cap
    core.state=ba::GAP_SLOT_PARKED;
    ba::aln_meth_flush_parked(&group);
    REQUIRE(result.occ.size()==1);CHECK(result.occ[0].cost==opt->s_gapo); // the gapped path's row
    ba::aln_meth_slot_finish(&slot);
    REQUIRE(result.main>=0);CHECK(result.flags!=0);CHECK(result.n_rescored==1);
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

TEST_CASE("aln_meth: stopped searches still score known-SNV extra windows unless CIGAR-capped"
          * doctest::test_suite("unit/aln_meth"))
{
    namespace ba = bwa_aln;
    const std::string fq = write_fq("stopped-extra.fq", {conv(ref().c1.substr(500, 50), 'C', 'T')});
    ba::bwa_seqio_t *input = ba::bwa_seq_open(fq.c_str());
    REQUIRE(input != nullptr);
    int n = 0;
    ba::bwa_seq_t *reads = ba::bwa_read_seq(input, 1, &n, BWA_MODE_COMPREAD, 0);
    REQUIRE(n == 1);
    FMI_search fmi((ref().fa + ".meth").c_str());
    fmi.load_index(false, 1);
    ba::aln_meth_ref_t mr;
    REQUIRE(ba::aln_meth_ref_load(&mr, ref().fa.c_str(), fmi.idx->bns, 1) == 0);
    auto *opt = ba::gap_init_opt();
    ba::aln_meth_opt_t mo;
    ba::aln_meth_opt_finish(&mo);
    ba::aln_meth_group_t group;
    group.ref = &mr; group.mo = &mo; group.opt = opt;
    ba::aln_meth_slot_t slot;
    slot.grp = &group;
    const int64_t rb = mr.bns->anns[0].offset + 500;
    // One known SNV inside the window, so it is scored as the reference and as that site's ALT.
    const int ref_code = (int)std::string("ACGT").find(ref().c1[520]);
    REQUIRE(ref_code >= 0);
    ba::ks_sites_t sites;
    sites.by_contig.assign(mr.bns->n_seqs, {});
    sites.by_contig[0] = {{520, (uint8_t)ref_code, (uint8_t)((ref_code + 1) & 3)}};
    sites.n_sites = 1;
    mr.known_snvs = &sites;
    for (bool path : {false, true}) for (int flags : std::vector<int>{0, ba::ALN_METH_BIND_BEST,
        ba::ALN_METH_BIND_CIGAR, ba::ALN_METH_BIND_CIGAR | ba::ALN_METH_BIND_BEST}) {
        CAPTURE(path); CAPTURE(flags);
        ba::aln_meth_read_t result;
        ba::aln_meth_slot_start(&slot, reads, &result, 0, 3, 1);
        std::vector<ba::aln_meth_candidate_t> extra{{0, 0, 1, rb, rb + 50, {}}};
        if (path) extra[0].cig = {50u << 4};
        slot.extra = &extra;
        slot.stop = 1; // lookup budget, certificate or search-cut stop is not a CIGAR cap
        result.flags = flags;
        ba::aln_meth_slot_finish(&slot);
        const bool capped = flags & ba::ALN_METH_BIND_CIGAR;
        CHECK(result.flags == flags);
        CHECK(result.n_rescored == (capped ? 0 : 2)); // the reference and the site allele must still be considered
        if (capped) {
            CHECK(result.main == -1);
            CHECK(result.hits.empty());
        } else {
            REQUIRE(result.main >= 0);
            CHECK(result.hits[result.main].pos == (bwtint_t)rb);
            CHECK(result.hits[result.main].cost == 0);
        }
        slot.extra = nullptr;
    }
    { // A read carrying the site's ALT: the site allele's row is cheaper than the reference's (one
      // mismatch there), so the hit credits it, which only scoring that allele can do.
        std::string alt_read = ref().c1.substr(500, 50);
        alt_read[20] = "ACGT"[(ref_code + 1) & 3];
        const std::string alt_fq = write_fq("stopped-extra-alt.fq", {conv(alt_read, 'C', 'T')});
        ba::bwa_seqio_t *alt_input = ba::bwa_seq_open(alt_fq.c_str());
        REQUIRE(alt_input != nullptr);
        int n_alt = 0;
        ba::bwa_seq_t *alt_reads = ba::bwa_read_seq(alt_input, 1, &n_alt, BWA_MODE_COMPREAD, 0);
        REQUIRE(n_alt == 1);
        for (bool path : {false, true}) {
            CAPTURE(path);
            ba::aln_meth_read_t result;
            ba::aln_meth_slot_start(&slot, alt_reads, &result, 0, 3, 1);
            std::vector<ba::aln_meth_candidate_t> extra{{0, 0, 1, rb, rb + 50, {}}};
            if (path) extra[0].cig = {50u << 4};
            slot.extra = &extra;
            slot.stop = 1;
            ba::aln_meth_slot_finish(&slot);
            CHECK(result.n_rescored == 2);
            REQUIRE(result.main >= 0);
            const auto &h = result.hits[result.main];
            CHECK(h.pos == (bwtint_t)rb);
            CHECK(h.alt == 1); // the window's first (only) site allele
            CHECK(h.cost == 0);
            CHECK(h.n_mm == 0);
            REQUIRE(h.n_cig == 1);
            CHECK(ba::aln_meth_hit_cigar(result, h)[0] == (50u << 4));
            slot.extra = nullptr;
        }
        ba::bwa_free_read_seq(n_alt, alt_reads);
        CHECK(ba::bwa_seq_close(alt_input) == 0);
    }
    free(opt);
    mr.known_snvs = nullptr;
    ba::aln_meth_ref_free(&mr);
    ba::bwa_free_read_seq(n, reads);
    CHECK(ba::bwa_seq_close(input) == 0);
}

TEST_CASE("aln_meth: past the mirror table's CIGAR budget the conversion-blind cap is 0")
{
    // The read is c1[500, 550) converted, with a C kept over one of its Ts (a mirror mismatch):
    // genomically one mismatch at a unique locus, MAPQ 37. The mirror table keeps that row with
    // its CIGARs (8 bytes: one raw and one display operation) under a 64 MiB budget of its own.
    // Seeding the table's counter rather than allocating 64 MiB: with room for the row the cap
    // leaves 37; without it the mirror-free strata are unknown, so the cap is 0, and the genomic
    // result is unchanged.
    namespace ba = bwa_aln;
    const std::string seg = ref().c1.substr(500, 50);
    std::string read = conv(seg, 'C', 'T');
    int at = -1;
    for (int i = 10; i < 40 && at < 0; ++i)
        if (seg[i] == 'T') at = i;
    REQUIRE(at >= 0);
    read[at] = 'C';
    const std::string fq = write_fq("mirror-budget.fq", {read});
    ba::bwa_seqio_t *input = ba::bwa_seq_open(fq.c_str());
    REQUIRE(input != nullptr);
    int n = 0;
    ba::bwa_seq_t *reads = ba::bwa_read_seq(input, 1, &n, BWA_MODE_COMPREAD, 0);
    REQUIRE(reads != nullptr);
    REQUIRE(n == 1);
    const std::string prefix = ref().fa + ".meth";
    FMI_search fmi(prefix.c_str());
    fmi.load_index(false, 1);
    ba::bwt_t *bwt = ba::bwt_from_fmi(&fmi);
    ba::aln_meth_ref_t mr;
    REQUIRE(ba::aln_meth_ref_load(&mr, ref().fa.c_str(), fmi.idx->bns, 1) == 0);
    ba::aln_meth_opt_t mo;
    mo.m.meth_mode = 1;
    mo.m.meth_mapq = MEM_METH_MAPQ_CONVERSION_BLIND;
    ba::aln_meth_opt_finish(&mo);
    ba::gap_opt_t *opt = ba::gap_init_opt();
    ba::aln_meth_group_t group;
    group.bwt = bwt; group.ref = &mr; group.mo = &mo; group.opt = opt;
    group.n_slots = 1; group.slot.resize(1);
    ba::gap_slot_core_t core = {};
    group.cores = &core;
    auto &slot = group.slot[0];
    slot.grp = &group; slot.job = &core.job;
    bwtint_t row = (bwtint_t)-1;
    for (bwtint_t k = 0; k <= bwt->seq_len; ++k) {
        int strand = 0;
        const auto pos = ba::bwa_pos_from_sa(mr.conv_bns, ba::bwt_sa(bwt, k), 50, &strand);
        if (pos != (bwtint_t)-1 && strand == 0 && bns_pos2rid(mr.conv_bns, pos) == 1 &&
            pos == (bwtint_t)(mr.conv_bns->anns[1].offset + 500)) { row = k; break; }
    }
    REQUIRE(row < bwt->seq_len);
    uint64_t remaining = 8;
    SUBCASE("room for the row") { remaining = 8; }
    SUBCASE("one byte short") { remaining = 7; }
    ba::aln_meth_read_t result;
    ba::aln_meth_slot_start(&slot, reads, &result, 0, 3, 0);
    slot.mirror_bytes = (64ULL << 20) - remaining;
    ba::gap_entry_t exact = {};
    ba::aln_meth_hit(&slot, &exact, row, row);
    core.state = ba::GAP_SLOT_PARKED;
    ba::aln_meth_flush_parked(&group);
    ba::aln_meth_slot_finish(&slot);
    const bool over = remaining < 8;
    CHECK(result.blind_partial == (over ? 1 : 0));
    CHECK(result.mirror_rows.empty()); // released at the finish either way
    REQUIRE(result.main >= 0);
    CHECK(result.flags == 0); // the genomic result does not see the table's budget
    REQUIRE(result.hits.size() == 1);
    CHECK(result.hits[0].pos == 500);
    CHECK(result.hits[0].n_mm == 1);
    CHECK(result.hits[0].n_mirror == 1);
    ba::bwase_initialize();
    ba::aln_meth_finish_se(&mr, &mo, opt, 3, 1, reads, &result);
    CHECK((int)reads[0].type == BWA_TYPE_UNIQUE);
    CHECK((int)reads[0].mapQ == (over ? 0 : 37));
    CHECK(result.yq == (over ? 37 : -1));
    free(opt);
    ba::bwt_destroy(bwt);
    ba::aln_meth_ref_free(&mr);
    ba::bwa_free_read_seq(n, reads);
    CHECK(ba::bwa_seq_close(input) == 0);
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
    opt->mode |= BWA_MODE_NONSTOP; // keep the gapped competitor eligible past best + M
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
    longer.n_del = longer.n_gapo = 1;
    longer.last_diff_pos = (int)ALN_METH_GAP(49, 1, 1); // leading reference deletion
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
    CHECK(slot.stop == ba::ALN_METH_STOP_CIGAR);
    CHECK(results[0].n_lookups == (retained ? 2 : 1));
    CHECK(slot.pend.empty());
    // A read stopped by the CIGAR budget must not even visit later ALT-window candidates.
    std::vector<ba::aln_meth_candidate_t> extra = {{0, 0, 1, 600, 600 + ALN_MAX_READ_LEN + 1}};
    slot.extra = &extra;
    ba::aln_meth_slot_finish(&slot);
    CHECK(results[0].n_unrepresentable == 0);
    CHECK((results[0].main >= 0) == retained);
    if (retained) {
        CHECK(results[0].hits[results[0].main].pos == 500);
        CHECK(results[0].hits[results[0].main].strand == 0);
    }
    // Reuse the same slot for a later read: the cap and stop state belong to one read only.
    slot.extra = nullptr;
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

TEST_CASE("aln_meth: a CIGAR stop survives the lookup budget ending the same flush")
{
    namespace ba = bwa_aln;
    const std::string fq = write_fq("cigar-lookup-cap.fq", {conv(ref().c1.substr(500, 50), 'C', 'T')});
    ba::bwa_seqio_t *input = ba::bwa_seq_open(fq.c_str());
    REQUIRE(input != nullptr);
    int n = 0;
    ba::bwa_seq_t *reads = ba::bwa_read_seq(input, 1, &n, BWA_MODE_COMPREAD, 0);
    REQUIRE(reads != nullptr);
    REQUIRE(n == 1);
    const std::string prefix = ref().fa + ".meth";
    FMI_search fmi(prefix.c_str());
    fmi.load_index(false, 1);
    ba::bwt_t *bwt = ba::bwt_from_fmi(&fmi);
    ba::aln_meth_ref_t mr;
    REQUIRE(ba::aln_meth_ref_load(&mr, ref().fa.c_str(), fmi.idx->bns, 1) == 0);
    ba::aln_meth_opt_t mo;
    ba::aln_meth_opt_finish(&mo);
    mo.p_best = 1; mo.p_second = 0; // one genomic lookup fits the budget
    ba::gap_opt_t *opt = ba::gap_init_opt();
    ba::aln_meth_group_t group;
    group.bwt = bwt; group.ref = &mr; group.mo = &mo; group.opt = opt;
    group.n_slots = 1; group.slot.resize(1);
    ba::gap_slot_core_t core = {};
    group.cores = &core;
    auto &slot = group.slot[0];
    slot.grp = &group; slot.job = &core.job;
    ba::aln_meth_read_t result;
    ba::aln_meth_slot_start(&slot, reads, &result, 0, 3, 1);
    bwtint_t row = (bwtint_t)-1;
    for (bwtint_t k = 0; k <= bwt->seq_len; ++k) {
        int strand = 0;
        const auto pos = ba::bwa_pos_from_sa(mr.conv_bns, ba::bwt_sa(bwt, k), 50, &strand);
        if (pos != (bwtint_t)-1 && strand == 0 && bns_pos2rid(mr.conv_bns, pos) == 1 &&
            pos == (bwtint_t)(mr.conv_bns->anns[1].offset + 500)) { row = k; break; }
    }
    REQUIRE(row < bwt->seq_len);
    // The one looked-up row's 8-byte alignment overflows the CIGAR budget, and the interval's
    // second row is past the lookup budget: both budgets end this flush.
    slot.alignment_bytes = (64ULL << 20) - 7;
    ba::gap_entry_t exact = {};
    ba::aln_meth_hit(&slot, &exact, row, row + 1);
    core.state = ba::GAP_SLOT_PARKED;
    ba::aln_meth_flush_parked(&group);
    CHECK(result.n_lookups == 1);
    CHECK(result.flags == (ba::ALN_METH_BIND_CIGAR | ba::ALN_METH_BIND_BEST));
    // The CIGAR stop is final: the lookup stop must not replace it.
    CHECK(slot.stop == ba::ALN_METH_STOP_CIGAR);
    // So a pending ALT-window candidate is not resolved past the CIGAR budget.
    std::vector<ba::aln_meth_candidate_t> extra = {{0, 0, 1, 600, 600 + ALN_MAX_READ_LEN + 1}};
    slot.extra = &extra;
    ba::aln_meth_slot_finish(&slot);
    CHECK(result.n_unrepresentable == 0);
    CHECK(slot.windows.size() == 1); // only the looked-up row's window
    CHECK(result.main == -1);
    ba::aln_meth_ref_free(&mr);
    ba::bwt_destroy(bwt);
    free(opt);
    ba::bwa_free_read_seq(n, reads);
    CHECK(ba::bwa_seq_close(input) == 0);
}

TEST_CASE("aln_meth: a genomic lookup-budget stop still resolves the ALT-window candidates")
{
    namespace ba = bwa_aln;
    const std::string fq = write_fq("lookup-cap.fq", {conv(ref().c1.substr(500, 50), 'C', 'T')});
    ba::bwa_seqio_t *input = ba::bwa_seq_open(fq.c_str());
    REQUIRE(input != nullptr);
    int n = 0;
    ba::bwa_seq_t *reads = ba::bwa_read_seq(input, 1, &n, BWA_MODE_COMPREAD, 0);
    REQUIRE(reads != nullptr);
    REQUIRE(n == 1);
    const std::string prefix = ref().fa + ".meth";
    FMI_search fmi(prefix.c_str());
    fmi.load_index(false, 1);
    ba::bwt_t *bwt = ba::bwt_from_fmi(&fmi);
    ba::aln_meth_ref_t mr;
    REQUIRE(ba::aln_meth_ref_load(&mr, ref().fa.c_str(), fmi.idx->bns, 1) == 0);
    ba::aln_meth_opt_t mo;
    ba::aln_meth_opt_finish(&mo);
    mo.p_best = 0; mo.p_second = 0; // no genomic lookup fits the budget
    ba::gap_opt_t *opt = ba::gap_init_opt();
    ba::aln_meth_group_t group;
    group.bwt = bwt; group.ref = &mr; group.mo = &mo; group.opt = opt;
    group.n_slots = 1; group.slot.resize(1);
    ba::gap_slot_core_t core = {};
    group.cores = &core;
    auto &slot = group.slot[0];
    slot.grp = &group; slot.job = &core.job;
    ba::aln_meth_read_t result;
    ba::aln_meth_slot_start(&slot, reads, &result, 0, 3, 1);
    ba::gap_entry_t exact = {};
    ba::aln_meth_hit(&slot, &exact, 0, 0);
    core.state = ba::GAP_SLOT_PARKED;
    ba::aln_meth_flush_parked(&group);
    CHECK(result.flags == ba::ALN_METH_BIND_BEST);
    CHECK(slot.stop == ba::ALN_METH_STOP_LOOKUP);
    CHECK(result.n_lookups == 0);
    CHECK(result.occ.empty());
    // The ALT-window search already found the read's exact window: it is still resolved.
    const int64_t rb = (int64_t)mr.bns->anns[0].offset + 500;
    std::vector<ba::aln_meth_candidate_t> extra = {{0, 0, 1, rb, rb + 50}};
    slot.extra = &extra;
    ba::aln_meth_slot_finish(&slot);
    CHECK(result.flags == ba::ALN_METH_BIND_BEST);
    REQUIRE(result.main >= 0);
    CHECK(result.hits[result.main].pos == 500);
    CHECK(result.hits[result.main].strand == 0);
    ba::bwase_initialize();
    ba::aln_meth_finish_se(&mr, &mo, opt, 3, 1, reads, &result);
    // Still incomplete: MAPQ 0, no exact X0/X1.
    CHECK((int)reads[0].mapQ == 0);
    CHECK((int)reads[0].c1 == 0);
    CHECK((int)reads[0].c2 == 0);
    kstring_t text = {0, 0, nullptr};
    ba::aln_meth_print_sam1(&text, &mr, &mo, reads, &result, 1, nullptr, BWA_MODE_COMPREAD, 30);
    const std::string record(text.s, text.l);
    free(text.s);
    CHECK(record.find("r0\t0\tc1\t501\t0\t50M\t") == 0);
    CHECK(record.find("X0:i:") == std::string::npos);
    CHECK(record.find("X1:i:") == std::string::npos);
    ba::aln_meth_ref_free(&mr);
    ba::bwt_destroy(bwt);
    free(opt);
    ba::bwa_free_read_seq(n, reads);
    CHECK(ba::bwa_seq_close(input) == 0);
}

TEST_CASE("aln_meth: ineligible searches batch pending hits instead of flushing score boundaries")
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
        if (mode==1) job.seed_width=&seed, job.opt.max_seed_diff=0; // seeded, no seed differences
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

TEST_CASE("aln_meth: suppressed boundaries do not move lookup-cap detection into a stratum")
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

TEST_CASE("aln_meth: mixed seed eligibility and caps survive real slot reuse and final flushing")
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
            // Finishing frees the search's loci; the stored ones keep what pairing and output read.
            CHECK(actual.occ.empty()); CHECK(single.occ.empty());
            REQUIRE(actual.hits.size()==single.hits.size());
            for (size_t k=0;k<actual.hits.size();++k) {
                const auto &a=actual.hits[k], &b=single.hits[k];
                CHECK(std::make_tuple(a.rid,a.hyp,a.strand,a.pos,a.cost,a.span)==
                      std::make_tuple(b.rid,b.hyp,b.strand,b.pos,b.cost,b.span));
                CHECK(std::make_tuple(a.n_mm,a.n_gapo,a.n_gape,a.n_conv)==
                      std::make_tuple(b.n_mm,b.n_gapo,b.n_gape,b.n_conv));
                CHECK(cigar_of(actual,a)==cigar_of(single,b));
            }
            if (budget==4596 && max_entries==2000000) {
                CHECK(actual.flags==0);
                CHECK_FALSE(actual.hits.empty()); // Pending exact hits must reach final refinement.
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
    ba::aln_meth_occ_t exact{}; exact.cost=0; result.occ.push_back(exact);
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
    // A seeded read stops only on a stratum within the seed's difference limit (-k, 2 here).
    ba::bwt_width_t seed={}; job.seed_width=&seed; REQUIRE(job.opt.max_seed_diff==2);
    CHECK(ba::aln_meth_can_stop(&slot,opt->s_mm+1));     // best 0: one difference
    result.occ[0].cost=3; CHECK_FALSE(ba::aln_meth_can_stop(&slot,6));
    CHECK(ba::aln_meth_can_stop(&slot,7));                 // best 3: two differences
    result.occ[0].cost=6; CHECK_FALSE(ba::aln_meth_can_stop(&slot,10)); // three: past -k
    job.seed_width=nullptr; CHECK(ba::aln_meth_can_stop(&slot,10)); job.seed_width=&seed;
    result.occ[0].cost=0; job.opt.max_seed_diff=1;
    CHECK(ba::aln_meth_can_stop(&slot,opt->s_mm+1));
    result.occ[0].cost=3; CHECK_FALSE(ba::aln_meth_can_stop(&slot,7));
    result.occ[0].cost=0; job.opt.max_seed_diff=0;
    CHECK_FALSE(ba::aln_meth_stop_eligible(&slot)); CHECK_FALSE(ba::aln_meth_can_stop(&slot,opt->s_mm+1));
    job.opt.max_seed_diff=2; job.seed_width=nullptr;
    job.opt.mode|=BWA_MODE_NONSTOP; CHECK_FALSE(ba::aln_meth_can_stop(&slot,opt->s_mm+1));
    job.opt.mode&=~BWA_MODE_NONSTOP;
    job.opt.max_diff=0; CHECK_FALSE(ba::aln_meth_can_stop(&slot,opt->s_mm+1)); job.opt.max_diff=4;
    job.opt.s_gapo=job.opt.s_mm; CHECK_FALSE(ba::aln_meth_can_stop(&slot,opt->s_mm+1));
    job.opt.s_gapo=11; job.opt.s_mm=0; CHECK_FALSE(ba::aln_meth_can_stop(&slot,1));
    job.opt.s_mm=3; job.opt.s_gape=-1; CHECK_FALSE(ba::aln_meth_can_stop(&slot,4));
    free(opt);
}

TEST_CASE("aln_meth: ALT width capture survives slot reuse reset and duplicate compaction")
{
    namespace ba = bwa_aln;
    FMI_search fmi((ref().fa + ".meth").c_str()); fmi.load_index(false, 1);
    auto *bwt = ba::bwt_from_fmi(&fmi); ba::aln_meth_ref_t mr;
    REQUIRE(ba::aln_meth_ref_load(&mr, ref().fa.c_str(), fmi.idx->bns, 1) == 0);
    auto *opt = ba::gap_init_opt(); opt->fnr = -1; opt->max_diff = 2;
    ba::aln_meth_opt_t mo; ba::aln_meth_opt_finish(&mo);
    std::vector<std::string> queries;
    for (int i = 0; i < 36; ++i) queries.push_back(ref().c1.substr(400 + (i % 19) * 3, i % 2 ? 35 : 17));
    queries.push_back(queries[0]); queries.push_back(queries[4]);
    const auto fq = write_fq("width-capture.fq", queries);
    for (int read_no : {0, 1}) {
        opt->seed_len = 18; // seeded; the round's unseeded compaction check below sets 1024
        auto *input = ba::bwa_seq_open(fq.c_str()); int n = 0;
        auto *reads = ba::bwa_read_seq(input, 100, &n, BWA_MODE_COMPREAD, 0); REQUIRE(n == (int)queries.size());
        std::vector<int> expected(n);
        for (int i = 0; i < n; ++i) {
            std::vector<uint8_t> converted(reads[i].len); std::vector<ba::bwt_width_t> width(reads[i].len + 1);
            ba::aln_meth_convert(read_no, reads[i].seq, reads[i].len, converted.data());
            ba::bwt_cal_width(bwt, reads[i].len, converted.data(), width.data());
            expected[i] = width[reads[i].len - 1].bid;
            CHECK((int)width[reads[i].len].bid == expected[i] + 1); // sentinel is NOT the certificate
        }
        std::vector<ba::aln_meth_read_t> batch(n); ba::aln_dup_state_t dup;
        CHECK(ba::aln_meth_search_batch(bwt, &mr, &mo, n, reads, opt, 2, read_no, 0, 1, 1000,
                                      &dup, batch.data(), nullptr, nullptr, true) >= 2);
        for (int i = 0; i < n; ++i) {
            CHECK(batch[i].full_width_bound.valid); CHECK(batch[i].full_width_bound.value == expected[i]);
        }
        // Reusing the result with capture disabled must clear old validity even on duplicate copies.
        dup = ba::aln_dup_state_t();
        ba::aln_meth_search_batch(bwt, &mr, &mo, n, reads, opt, 1, read_no, 0, 1, 1000, &dup, batch.data());
        for (const auto &r : batch) CHECK_FALSE(r.full_width_bound.valid);
        // Differing opt-in metadata prevents inappropriate sharing, and compaction uses uniq[k].
        ba::ks_sites_t no_sites; no_sites.by_contig.resize(mr.bns->n_seqs); mr.known_snvs = &no_sites;
        std::vector<ba::aln_meth_alt_bound_t> bounds(n);
        for (int i = 0; i < n; ++i) { bounds[i].valid = true; bounds[i].value = i; }
        dup = ba::aln_dup_state_t();
        CHECK(ba::aln_meth_search_batch(bwt, &mr, &mo, n, reads, opt, 1, read_no, 0, 1, 1000,
                                      &dup, batch.data(), nullptr, bounds.data()) == 0);
        for (int i = 0; i < n; ++i) {
            ba::aln_meth_read_t single;
            ba::aln_meth_search_batch(bwt, &mr, &mo, 1, reads + i, opt, 1, read_no, i, 1, 1000,
                                      nullptr, &single, nullptr, &bounds[i]);
            CHECK(batch[i].flags == single.flags); CHECK(batch[i].c1 == single.c1); CHECK(batch[i].c2 == single.c2);
            CHECK(batch[i].n_lookups == single.n_lookups);
        }
        // A duplicate before later distinct reads makes uniq[k] nonidentity; mixed bounds must
        // retain each later read's own stopping policy rather than the compacted position's.
        opt->seed_len = 1024;
        const auto swap_read = reads[1]; reads[1] = reads[n - 2]; reads[n - 2] = swap_read;
        std::vector<std::string> bases = queries; std::swap(bases[1], bases[n - 2]);
        for (int i = 0; i < n; ++i) { bounds[i].valid = true; bounds[i].value = reads[i].len == 17 ? 0 : 100; }
        // A copy repeats an earlier read's bases and bound; every other read is a representative.
        std::vector<bool> is_copy(n, false); int n_copies = 0;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < i && !is_copy[i]; ++j)
                is_copy[i] = bases[j] == bases[i] && bounds[j].value == bounds[i].value;
            n_copies += is_copy[i];
        }
        REQUIRE(n_copies == 2); // the two repeated queries, one now before later distinct reads
        CHECK(is_copy[1]);
        dup = ba::aln_dup_state_t();
        CHECK(ba::aln_meth_search_batch(bwt, &mr, &mo, n, reads, opt, 1, read_no, 0, 1, 1000,
                                      &dup, batch.data(), nullptr, bounds.data()) == n_copies);
        for (int i = 0; i < n; ++i) {
            CAPTURE(i);
            ba::aln_meth_read_t single;
            ba::aln_meth_search_batch(bwt, &mr, &mo, 1, reads + i, opt, 1, read_no, i, 1, 1000,
                                      nullptr, &single, nullptr, &bounds[i]);
            CHECK(batch[i].flags == single.flags); CHECK(batch[i].c1 == single.c1); CHECK(batch[i].c2 == single.c2);
            // Duplicate copies intentionally have zero work counters; a representative's are its own.
            CHECK(batch[i].n_lookups == (is_copy[i] ? 0 : single.n_lookups));
        }
        mr.known_snvs = nullptr;
        ba::bwa_free_read_seq(n, reads); ba::bwa_seq_close(input);
    }
    free(opt); ba::aln_meth_ref_free(&mr); ba::bwt_destroy(bwt);
}

TEST_CASE("aln_meth: known-site stop requires a strict valid same-read ALT certificate")
{
    namespace ba = bwa_aln;
    auto *opt = ba::gap_init_opt();
    ba::ks_sites_t sites; ba::aln_meth_ref_t ref; ref.known_snvs = &sites;
    ba::aln_meth_group_t group; group.ref = &ref; group.opt = opt;
    ba::gap_job_t job = {}; job.opt = *opt; job.opt.max_diff = 3;
    ba::aln_meth_read_t res; ba::aln_meth_occ_t best = {}; best.cost = 0; res.occ.push_back(best);
    ba::aln_meth_slot_t slot; slot.grp = &group; slot.job = &job; slot.res = &res; slot.len = 30;
    CHECK_FALSE(ba::aln_meth_can_stop(&slot, 4)); // unavailable is not zero candidates
    slot.alt_bound.valid = true; slot.alt_bound.value = 1;
    CHECK_FALSE(ba::aln_meth_can_stop(&slot, 4)); // LB == d cannot certify
    slot.alt_bound.value = 2;
    CHECK_FALSE(ba::aln_meth_can_stop(&slot, 3)); // next == M*d cannot certify
    CHECK(ba::aln_meth_can_stop(&slot, 4));
    // No extras list is consulted: a competitive ALT can lie beyond the reference horizon.
    std::vector<ba::aln_meth_candidate_t> empty; slot.extra = &empty;
    slot.alt_bound.value = 0; CHECK_FALSE(ba::aln_meth_can_stop(&slot, 100));
    slot.alt_bound.value = 2;
    for (int len : {0, 160, 161, 350}) {
        slot.len = len; CHECK(ba::aln_meth_can_stop(&slot, 4) == (len == 160));
    }
    slot.len = 30;
    for (int flag : {ba::ALN_METH_BIND_BEST, ba::ALN_METH_BIND_UNREP, ba::ALN_METH_BIND_CIGAR, ba::ALN_METH_BIND_M}) {
        res.flags = flag; CHECK_FALSE(ba::aln_meth_can_stop(&slot, 100));
    }
    res.flags = 0;
    ba::bwt_width_t seed = {}; job.seed_width = &seed;
    CHECK(ba::aln_meth_can_stop(&slot, 4)); // d=1 fits the seed and LB=2 excludes ALT competitors
    job.opt.max_seed_diff = 0; CHECK_FALSE(ba::aln_meth_can_stop(&slot, 100));
    job.opt.max_seed_diff = 1; CHECK(ba::aln_meth_can_stop(&slot, 4)); // d == -k
    slot.alt_bound.valid = false; CHECK_FALSE(ba::aln_meth_can_stop(&slot, 100));
    slot.alt_bound.valid = true; slot.alt_bound.value = 1; CHECK_FALSE(ba::aln_meth_can_stop(&slot, 100));
    slot.alt_bound.value = 3; res.occ[0].cost = 3; CHECK_FALSE(ba::aln_meth_can_stop(&slot, 100)); // d > -k
    job.opt.max_seed_diff = 2; CHECK(ba::aln_meth_can_stop(&slot, 7));
    res.occ[0].cost = 0; slot.alt_bound.value = 2;
    for (int len : {0, 160, 161, 350}) {
        slot.len = len; CHECK(ba::aln_meth_can_stop(&slot, 4) == (len == 160));
    }
    slot.len = 30; job.seed_width = nullptr;
    job.opt.mode |= BWA_MODE_NONSTOP; CHECK_FALSE(ba::aln_meth_can_stop(&slot, 100)); job.opt.mode &= ~BWA_MODE_NONSTOP;
    job.opt.s_gapo = 3; CHECK_FALSE(ba::aln_meth_can_stop(&slot, 100)); job.opt.s_gapo = 11;
    for (int m : {0, -1}) { job.opt.s_mm = m; CHECK_FALSE(ba::aln_meth_can_stop(&slot, 100)); }
    job.opt.s_mm = 3; job.opt.s_gape = -1; CHECK_FALSE(ba::aln_meth_can_stop(&slot, 100)); job.opt.s_gape = 0;
    job.opt.max_diff = 0; CHECK_FALSE(ba::aln_meth_can_stop(&slot, 100)); job.opt.max_diff = 3;
    job.opt.mode &= ~BWA_MODE_GAPE; CHECK(ba::aln_meth_can_stop(&slot, 4));
    res.occ[0].cost = 3; slot.alt_bound.value = 2; CHECK_FALSE(ba::aln_meth_can_stop(&slot, 7));
    slot.alt_bound.value = 3; CHECK(ba::aln_meth_can_stop(&slot, 7));
    free(opt);
}

TEST_CASE("aln_meth: known-site certificate preserves independent REF plus ALT window strata")
{
    namespace ba = bwa_aln;
    bwa_tests::TempDir tmp{"aln_meth_alt_stop"}; const std::string fa = tmp.path() + "/ref.fa";
    const std::string query = "ACGTCGATCGTA";
    std::vector<std::pair<std::string, std::string>> contigs{{"exact", query}};
    for (int i = 0; i < (int)query.size(); ++i) for (char b : std::string("ACGTN")) {
        std::string target = query; target[i] = b;
        contigs.push_back({"v" + std::to_string(contigs.size()), target});
    }
    // REF cost 2*M exceeds T=M; the single ALT at base zero restores cost M.
    // This competitor must survive even though genomic safe stopping alone could omit it.
    std::string improved = query; improved.front() = 'C'; improved.back() = 'C';
    contigs.push_back({"alt_beyond_ref_horizon", improved});
    REQUIRE(aln_test::write_test_fasta(fa, contigs));
    std::vector<std::string> args{"index", "--meth", fa}; std::vector<char *> argv;
    for (auto &a : args) argv.push_back(&a[0]); argv.push_back(nullptr); optind = 1;
#if defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || defined(__NetBSD__)
    optreset = 1;
#endif
    REQUIRE(bwa_index((int)args.size(), argv.data()) == 0);
    FMI_search fmi((fa + ".meth").c_str()); fmi.load_index(false, 1, true, BWAMEM_HUGETLB_OFF, false);
    auto *bwt = ba::bwt_from_fmi(&fmi); ba::aln_meth_ref_t mr;
    REQUIRE(ba::aln_meth_ref_load(&mr, fa.c_str(), fmi.idx->bns, 1) == 0);
    ba::ks_sites_t sites; sites.by_contig.resize(contigs.size());
    for (size_t i = 0; i < contigs.size(); ++i) {
        int p = i == contigs.size() - 1 ? 0 : i % 2 ? 0 : query.size() - 1; // both edges
        if (contigs[i].second[p] == 'N') p = p == 0 ? query.size() - 1 : 0; // sites require ACGT REF
        const char rb = contigs[i].second[p];
        sites.by_contig[i].push_back({p, (uint8_t)std::string("ACGT").find(rb),
                                     (uint8_t)(i == contigs.size() - 1 ? 0 : rb == 'A' ? 2 : 0)});
        ++sites.n_sites;
    }
    mr.known_snvs = &sites;
    ba::ks_windows_t windows; REQUIRE(ba::ks_windows_build(&sites, mr.bns, mr.pac, 1, 1, &windows) == 0);
    auto *opt = ba::gap_init_opt(); opt->fnr = -1; opt->max_diff = 2;
    using Key = std::tuple<int,int,int,int64_t,int64_t>;
    // The single-mismatch contigs put differences inside, outside and exactly at the seed edge.
    for (int seed_len : {1024, 6})
    for (int scoring : {MEM_METH_SCORING_NEUTRAL, MEM_METH_SCORING_COLLAPSED})
    for (int context : {METH_CTX_ALL, METH_CTX_CPG}) for (int read_no : {0, 1}) for (int hyp : {0, 1}) {
        ba::aln_meth_opt_t mo; mo.m.meth_scoring = scoring; mo.m.meth_context = context;
        mo.m.meth_chem = METH_CHEM_TAPS; ba::aln_meth_opt_finish(&mo);
        opt->seed_len = seed_len; CAPTURE(seed_len);
        const int strand = hyp ^ read_no;
        const auto fq = write_fq_in(tmp.path(), "r.fq", {strand ? aln_test::revcomp(query) : query});
        auto *input = ba::bwa_seq_open(fq.c_str()); int n = 0;
        auto *reads = ba::bwa_read_seq(input, 1, &n, BWA_MODE_COMPREAD, 0); REQUIRE(n == 1);
        ba::gap_opt_t discovery = *opt; discovery.mode |= BWA_MODE_NONSTOP; discovery.max_top2 = INT_MAX;
        ba::aln_meth_read_t wr;
        ba::aln_meth_search_batch(windows.bwt, windows.mref, &mo, 1, reads, &discovery, 1,
                                  read_no, 0, 1, 1000, nullptr, &wr, nullptr, nullptr, true);
        std::vector<ba::aln_meth_candidate_t> extra; ba::ks_meth_candidates(&windows, mr.bns, wr, &extra);
        ba::aln_meth_read_t result;
        // As ks_meth_search composes them: the production certificate and the window flags.
        const ba::aln_meth_alt_bound_t certificate = ba::aln_meth_alt_certificate(wr);
        const int wflags = wr.flags;
        ba::aln_meth_search_batch(bwt, &mr, &mo, 1, reads, opt, 1, read_no, 0, 1, 1000,
                                  nullptr, &result, &extra, &certificate, false, &wflags);
        CHECK(result.flags == 0); CHECK(wr.flags == 0);
        REQUIRE(wr.full_width_bound.valid);CHECK(wr.full_width_bound.value<=1);
        std::set<Key> expected, actual;
        for (int rid = 0; rid < (int)contigs.size(); ++rid) {
            int best = INT_MAX;
            for (bool alt : {false, true}) {
                std::string target = contigs[rid].second;
                if (alt) { const auto &s = sites.by_contig[rid][0]; target[s.pos] = "ACGT"[s.alt]; }
                int mm = 0;
                for (int k = 0; k < (int)query.size(); ++k) {
                    const char r = target[k], q = query[k], from = hyp ? 'C' : 'G', to = hyp ? 'T' : 'A';
                    const bool cpg = hyp ? (r == 'C' && k + 1 < (int)target.size() && target[k + 1] == 'G')
                                         : (r == 'G' && k > 0 && target[k - 1] == 'C');
                    const bool cv = r == from && q == to && (context == METH_CTX_ALL || cpg);
                    const bool mirror = context == METH_CTX_ALL && scoring == MEM_METH_SCORING_COLLAPSED && r == to && q == from;
                    mm += r == 'N' || (r != q && !cv && !mirror);
                }
                best = std::min(best, mm);
            }
            if (best <= 1) expected.emplace(rid, hyp, strand, mr.bns->anns[rid].offset, best * opt->s_mm);
        }
        for (const auto &h : result.hits) {
            INFO("duplicate " << h.rid << " hyp=" << h.hyp << " strand=" << h.strand << " pos=" << h.pos);
            CHECK(actual.emplace(h.rid, h.hyp, h.strand, h.pos, h.cost).second);
        }
        CAPTURE(scoring); CAPTURE(context); CAPTURE(read_no); CAPTURE(hyp);
        for (const auto &k : actual) if (!expected.count(k)) {
            INFO("unexpected " << std::get<0>(k) << " hyp=" << std::get<1>(k) << " strand=" << std::get<2>(k)
                 << " pos=" << std::get<3>(k) << " cost=" << std::get<4>(k)); CHECK(expected.count(k) == 1);
        }
        for (const auto &k : expected) if (!actual.count(k)) {
            INFO("missing " << std::get<0>(k) << " hyp=" << std::get<1>(k) << " strand=" << std::get<2>(k)
                 << " pos=" << std::get<3>(k) << " cost=" << std::get<4>(k)); CHECK(actual.count(k) == 1);
        }
        CHECK(actual == expected);
        int best_count = 0, second_count = 0;
        for (const auto &k : expected) { if (std::get<4>(k) == 0) ++best_count; else ++second_count; }
        CHECK(result.c1 + result.n_conv_moved == best_count);
        CHECK(result.c2 == second_count + result.n_conv_moved);
        ba::bwa_free_read_seq(n, reads); ba::bwa_seq_close(input);
    }
    free(opt); ba::ks_windows_free(&windows); ba::aln_meth_ref_free(&mr); ba::bwt_destroy(bwt);
}

TEST_CASE("aln_meth: seeded known-site searches use an actual ALT bound without losing retained loci")
{
    namespace ba = bwa_aln;
    bwa_tests::TempDir tmp{"aln_meth_seed_alt"};
    const std::string fa=tmp.file(tmp.path().substr(tmp.path().find_last_of('/')+1)+".fa");
    uint32_t state=314159; const std::string query=aln_test::lcg_dna(state,30);
    std::vector<std::pair<std::string,std::string>> contigs{{"exact",query},{"unrelated",std::string(321,'G')}};
    for(int i=0;i<30;++i) {
        auto target=query;
        for(int p:{i,(i+13)%30}) target[p]="CGTA"[std::string("ACGT").find(target[p])];
        contigs.push_back({"two_mm_"+std::to_string(i),target});
    }
    REQUIRE(aln_test::write_test_fasta(fa,contigs));
    std::vector<std::string> args{"index","--meth",fa};std::vector<char *> argv;
    for(auto &arg:args) argv.push_back(&arg[0]);argv.push_back(nullptr);optind=1;
#if defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || defined(__NetBSD__)
    optreset=1;
#endif
    REQUIRE(bwa_index((int)args.size(),argv.data())==0);
    FMI_search fmi((fa+".meth").c_str());fmi.load_index(false,1,true,BWAMEM_HUGETLB_OFF,false);
    auto *bwt=ba::bwt_from_fmi(&fmi);ba::aln_meth_ref_t ref;
    REQUIRE(ba::aln_meth_ref_load(&ref,fa.c_str(),fmi.idx->bns,1)==0);
    ba::ks_sites_t sites;sites.by_contig.resize(contigs.size());
    sites.by_contig[1].push_back({160,2,0});sites.n_sites=1;ref.known_snvs=&sites;
    ba::ks_windows_t windows;REQUIRE(ba::ks_windows_build(&sites,ref.bns,ref.pac,1,1,&windows)==0);
    auto *opt=ba::gap_init_opt();opt->fnr=-1;opt->max_diff=2;opt->seed_len=15;
    using Key=std::tuple<int,int,int,int64_t,int>;
    for(int context:{METH_CTX_ALL,METH_CTX_CPG}) for(int read_no:{0,1}) for(int hyp:{0,1}) {
        ba::aln_meth_opt_t mo;mo.m.meth_context=context;mo.m.meth_chem=METH_CHEM_TAPS;
        mo.m.meth_scoring=MEM_METH_SCORING_NEUTRAL;ba::aln_meth_opt_finish(&mo);
        const int strand=hyp^read_no;
        const auto fq=write_fq_in(tmp.path(),"r.fq",{strand?aln_test::revcomp(query):query});
        auto *input=ba::bwa_seq_open(fq.c_str());int n=0;
        auto *reads=ba::bwa_read_seq(input,1,&n,BWA_MODE_COMPREAD,0);REQUIRE(n==1);
        ba::gap_opt_t discovery=*opt;discovery.mode|=BWA_MODE_NONSTOP;discovery.max_top2=INT_MAX;
        ba::aln_meth_read_t wr;
        ba::aln_meth_search_batch(windows.bwt,windows.mref,&mo,1,reads,&discovery,1,read_no,0,1,1000,nullptr,&wr,nullptr,nullptr,true);
        REQUIRE(wr.flags==0);REQUIRE(wr.full_width_bound.valid);REQUIRE(wr.full_width_bound.value>1);
        std::vector<ba::aln_meth_candidate_t> extra;ba::ks_meth_candidates(&windows,ref.bns,wr,&extra);
        ba::aln_meth_read_t certified,uncertified;ba::aln_meth_alt_bound_t unavailable;
        ba::aln_meth_search_batch(bwt,&ref,&mo,1,reads,opt,1,read_no,0,1,1000,nullptr,&certified,&extra,&wr.full_width_bound);
        ba::aln_meth_search_batch(bwt,&ref,&mo,1,reads,opt,1,read_no,0,1,1000,nullptr,&uncertified,&extra,&unavailable);
        CHECK(certified.flags==0);CHECK(uncertified.flags==0);
        CHECK(certified.n_lookups<uncertified.n_lookups);
        std::set<Key> expected{{0,hyp,strand,0,0}},actual,control;
        for(const auto &hit:certified.hits) actual.emplace(hit.rid,hit.hyp,hit.strand,hit.pos,hit.cost);
        for(const auto &hit:uncertified.hits) control.emplace(hit.rid,hit.hyp,hit.strand,hit.pos,hit.cost);
        CHECK(actual==expected);CHECK(control==expected);
        for(const auto *result:{&certified,&uncertified}) {
            REQUIRE(result->hits.size()==1);
            const auto &hit=result->hits.front();
            CHECK(hit.rid==0);CHECK(hit.hyp==hyp);CHECK(hit.strand==strand);
            CHECK(hit.pos==0);CHECK(hit.span==30);CHECK(hit.pos+hit.span==30);
            CHECK(cigar_of(*result,hit)==std::vector<uint32_t>{30u<<4});
            CHECK(hit.n_cig==1);CHECK(hit.cost==0);CHECK(hit.n_mm==0);
            CHECK(hit.n_gapo==0);CHECK(hit.n_gape==0);CHECK(hit.n_conv==0);CHECK(hit.alt==0);
            CHECK(result->c1==1);CHECK(result->c2==0);CHECK(result->n_conv_moved==0);
            CHECK(result->main==0);CHECK(result->flags==0);
        }
        ba::bwa_free_read_seq(n,reads);ba::bwa_seq_close(input);
    }
    free(opt);ba::ks_windows_free(&windows);ba::aln_meth_ref_free(&ref);ba::bwt_destroy(bwt);
}

TEST_CASE("aln_meth: zero-cost stopping preserves an independently enumerated tiny candidate set")
{
    namespace ba=bwa_aln;
    bwa_tests::TempDir tmp{"aln_meth_zero_bound"}; // this case's index and reads, removed when it ends
    const std::string fa=tmp.path()+"/zero-bound.fa";
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
    FMI_search fmi((fa+".meth").c_str());fmi.load_index(false,1,true,BWAMEM_HUGETLB_OFF,false);
    auto *bwt=ba::bwt_from_fmi(&fmi);ba::aln_meth_ref_t mr;
    REQUIRE(ba::aln_meth_ref_load(&mr,fa.c_str(),fmi.idx->bns,1)==0);
    auto *opt=ba::gap_init_opt();opt->fnr=-1;opt->max_diff=2;
    using Key=std::tuple<int,int,int,int64_t,int64_t>;
    // Unseeded, and seeded (the 12-base read is longer than -l 6): the seeded search stops on
    // the same stratum, since one difference is within the seed's -k 2.
    for(int seed_len : {1024,6})
    for(int scoring : {MEM_METH_SCORING_COLLAPSED,MEM_METH_SCORING_GENOMIC,MEM_METH_SCORING_NEUTRAL})
        for(int context : {METH_CTX_ALL,METH_CTX_CPG}) for(int read_no : {0,1}) for(int hyp : {0,1}) {
            ba::aln_meth_opt_t mo;mo.m.meth_scoring=scoring;mo.m.meth_context=context;
            mo.m.meth_chem=METH_CHEM_TAPS;ba::aln_meth_opt_finish(&mo);
            opt->seed_len=seed_len;CAPTURE(seed_len);
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
            const auto fq=write_fq_in(tmp.path(),"zero-bound.fq",{input_read});auto *input=ba::bwa_seq_open(fq.c_str());
            int n=0;auto *reads=ba::bwa_read_seq(input,1,&n,BWA_MODE_COMPREAD,0);REQUIRE(n==1);
            ba::aln_meth_read_t result;ba::aln_dup_state_t dup;
            REQUIRE(ba::aln_meth_search_batch(bwt,&mr,&mo,1,reads,opt,2,read_no,0,1,1000,&dup,&result)==0);
            CHECK(result.flags==0);
            for(const auto &o:result.hits) actual.emplace(o.rid,o.hyp,o.strand,o.pos,o.cost);
            CHECK(actual==expected);
            ba::bwa_free_read_seq(n,reads);ba::bwa_seq_close(input);
        }
    free(opt);ba::aln_meth_ref_free(&mr);ba::bwt_destroy(bwt);
}

TEST_CASE("aln_meth: nonzero stopping preserves independent candidate strata with custom penalties")
{
    namespace ba=bwa_aln;
    bwa_tests::TempDir tmp{"aln_meth_nonzero_bound"}; // this case's index and reads, removed when it ends
    const std::string fa=tmp.path()+"/nonzero-bound.fa";
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
    FMI_search fmi((fa+".meth").c_str());fmi.load_index(false,1,true,BWAMEM_HUGETLB_OFF,false);
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
            const auto fq=write_fq_in(tmp.path(),"nonzero-bound.fq",{input_read});auto *input=ba::bwa_seq_open(fq.c_str());
            int n=0;auto *reads=ba::bwa_read_seq(input,1,&n,BWA_MODE_COMPREAD,0);REQUIRE(n==1);
            ba::aln_meth_read_t result;ba::aln_dup_state_t dup;
            REQUIRE(ba::aln_meth_search_batch(bwt,&mr,&mo,1,reads,opt,2,read_no,0,1,1000,&dup,&result)==0);
            CHECK(result.flags==0);
            for(const auto &o:result.hits) actual.emplace(o.rid,o.hyp,o.strand,o.pos,o.cost);
            CHECK(actual==expected);
            const int complete_c1=result.c1,complete_c2=result.c2;
            ba::bwa_free_read_seq(n,reads);ba::bwa_seq_close(input);
            // -N disables stopping but must produce exactly the same retained <=best+M strata.
            opt->mode|=BWA_MODE_NONSTOP;
            input=ba::bwa_seq_open(fq.c_str());reads=ba::bwa_read_seq(input,1,&n,BWA_MODE_COMPREAD,0);
            ba::aln_meth_read_t nonstop;ba::aln_dup_state_t dup_nonstop;
            REQUIRE(ba::aln_meth_search_batch(bwt,&mr,&mo,1,reads,opt,2,read_no,0,1,1000,&dup_nonstop,&nonstop)==0);
            CHECK(nonstop.flags==0);actual.clear();int within_c2=0;
            for(const auto &o:nonstop.hits) if(o.cost<=2*opt->s_mm) {
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
                    actual.clear();for(const auto &o:capped.hits)actual.emplace(o.rid,o.hyp,o.strand,o.pos,o.cost);
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
    // Every pair must be written, in input order, read 1 then read 2: matching truncated output
    // must not pass the comparisons below.
    auto all_pairs = [&](const std::vector<Rec> &recs) {
        REQUIRE(recs.size() == 2 * r1.size());
        for (size_t i = 0; i < r1.size(); ++i) {
            CAPTURE(i);
            CHECK(recs[2 * i].name == "r" + std::to_string(i));
            CHECK(recs[2 * i + 1].name == "r" + std::to_string(i));
            CHECK((recs[2 * i].flag & 0xc1) == 0x41);
            CHECK((recs[2 * i + 1].flag & 0xc1) == 0x81);
        }
    };
    all_pairs(run_aln({"--meth", "-t", "1"}, {f1, f2}));
    const std::string one = body(slurp_out());
    auto v = run_aln({"--meth", "-t", "3"}, {f1, f2});
    all_pairs(v);
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
    all_pairs(a);
    all_pairs(b);
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
        ba::aln_meth_occ_t o{};o.rid=0;o.hyp=1;o.strand=j;o.pos=base+(j?700:100);
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
    store_pair(r);ba::aln_meth_pe_stage1(opt,1,seqs,r);
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
            CHECK(rec[1][i].mapq==rec[0][i].mapq);
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

TEST_CASE("aln_meth: short-insert restores short overlap calibration in ordinary and meth paired output"
          * doctest::test_suite("unit/aln_meth"))
{
    bool chemistry = false, converted = false, rescue = false;
    SUBCASE("ordinary") {}
    SUBCASE("meth exact") { chemistry = true; }
    SUBCASE("meth converted") { chemistry = converted = true; }
    SUBCASE("ordinary with rescue enabled") { rescue = true; }
    SUBCASE("meth exact with rescue enabled") { chemistry = rescue = true; }
    SUBCASE("meth converted with rescue enabled") { chemistry = converted = rescue = true; }
    std::vector<std::string> one, two;
    for (int i = 0; i < 41; ++i) {
        const int start = 100 + 50 * i, len = i == 40 ? 150 : 30 + i % 16;
        std::string r1 = ref().c1.substr(start, len);
        std::string r2 = aln_test::revcomp(ref().c1.substr(start + 5, len));
        if (converted) { r1 = conv(r1, 'C', 'T'); r2 = conv(r2, 'G', 'A'); }
        one.push_back(r1); two.push_back(r2);
    }
    const std::vector<std::string> inputs = {write_fq("short-one.fq", one), write_fq("short-two.fq", two)};
    std::vector<std::string> args = {"-n", "0", "-l", "1024"};
    if (!rescue) args.push_back("-s");
    if (chemistry) args.push_back("--meth");
    const auto legacy = run_aln(args, inputs);
    REQUIRE(legacy.size() == 82);
    for (const auto &r : legacy) CHECK((r.flag & 2) == 0);
    args.push_back("--short-insert");
    const auto short_insert = run_aln(args, inputs);
    REQUIRE(short_insert.size() == 82);
    for (int i = 0; i < 41; ++i) for (int j = 0; j < 2; ++j) {
        CAPTURE(i); CAPTURE(j);
        const Rec &r = short_insert[2 * i + j];
        const Rec &old = legacy[2 * i + j];
        const int len = i == 40 ? 150 : 30 + i % 16;
        CHECK(r.name == "r" + std::to_string(i));
        const int start = 100 + 50 * i + 5 * j; // the fixture's mate start, 0-based
        CHECK(r.rname == "c1"); CHECK(r.pos == start + 1);
        CHECK(r.rnext == "="); CHECK(r.pnext == 101 + 50 * i + 5 * (1 - j));
        CHECK(r.cigar == std::to_string(len) + "M");
        CHECK(r.tlen == (j ? -1 : 1) * (len + 5));
        CHECK((r.flag & 16) == (j ? 16 : 0));
        CHECK((r.flag & 2) == (i < 40 ? 2 : 0));
        CHECK(old.mapq == 25); // -n 0 has no suboptimal mismatch stratum.
        CHECK(r.mapq == (i < 40 ? 50 : 25));
        CHECK(old.pos == r.pos); CHECK(old.cigar == r.cigar);
        CHECK(old.pnext == r.pnext); CHECK(old.tlen == r.tlen);
        CHECK(tag_value(r, "NM:i") == "0");
        const std::string aligned = j ? aln_test::revcomp(two[i]) : one[i];
        CHECK(tag_value(r, "MD:Z") == expected_md(ref().c1, start, aligned));
        if (chemistry) {
            CHECK(tag_value(r, "XR:Z") == (j ? "GA" : "CT"));
            CHECK(tag_value(r, "XG:Z") == "CT");
            CHECK(tag_value(r, "XM:Z") == expected_xm(ref().c1, start, aligned, true, false));
            CHECK(tag_value(old, "XM:Z") == tag_value(r, "XM:Z"));
        }
    }
    args.push_back("-A");
    const auto forced = run_aln(args, inputs);
    REQUIRE(forced.size() == 82);
    for (const auto &r : forced) CHECK((r.flag & 2) != 0);
}

TEST_CASE("aln_meth: short-insert installs a mate rescued at reference offset 0"
          * doctest::test_suite("unit/aln_meth"))
{
    // c1 is the first contig, so its first base is reference offset 0. The last pair's forward
    // read starts there and its mate covers the same 30 bases; with one strict mismatch the mate
    // is found only by rescue, at offset 0. An unmapped mate also has pos 0, so the rescue must
    // not be mistaken for "no change".
    bool meth = false;
    SUBCASE("ordinary") {}
    SUBCASE("meth") { meth = true; }
    std::vector<std::string> reads[2];
    for (int i = 0; i < 41; ++i) {
        const int start = 100 + 50 * i, len = 30 + i % 16;
        reads[0].push_back(ref().c1.substr(start, len));
        reads[1].push_back(aln_test::revcomp(ref().c1.substr(start + 5, len)));
    }
    reads[0].push_back(ref().c1.substr(0, 30));
    std::string mate = aln_test::revcomp(ref().c1.substr(0, 30));
    mate[15] = mate[15] == 'T' ? 'C' : 'T'; // strict under both ordinary and CT chemistry
    reads[1].push_back(mate);
    const std::vector<std::string> input = {write_fq("offset0-1.fq", reads[0]), write_fq("offset0-2.fq", reads[1])};
    std::vector<std::string> args = {"--short-insert", "-n", "0", "-l", "1024"};
    if (meth) args.push_back("--meth");
    const auto v = run_aln(args, input);
    REQUIRE(v.size() == 84);
    const auto &a = v[82], &b = v[83];
    CHECK(a.rname == "c1"); CHECK(a.pos == 1); CHECK(a.cigar == "30M"); CHECK((a.flag & 4) == 0);
    CHECK(b.rname == "c1"); CHECK(b.pos == 1); CHECK(b.cigar == "30M");
    CHECK((b.flag & 4) == 0); CHECK((b.flag & 16) != 0);
    CHECK(tag_value(b, "XT:A") == "M");
}

TEST_CASE("aln_meth: short-insert successfully rescues a singleton within unchanged rescue geometry"
          * doctest::test_suite("unit/aln_meth"))
{
    bool meth = false, converted = false;
    SUBCASE("ordinary") {}
    SUBCASE("meth exact") { meth = true; }
    SUBCASE("meth converted") { meth = converted = true; }
    std::vector<std::string> reads[2];
    for (int i = 0; i < 42; ++i) {
        const int start = 100 + 50 * i, len = i == 40 ? 150 : 30 + i % 16;
        std::string one = ref().c1.substr(start, i == 41 ? 20 : len);
        std::string two = aln_test::revcomp(ref().c1.substr(start + (i == 41 ? 20 : 5), i == 41 ? 30 : len));
        if (converted) { one = conv(one, 'C', 'T'); two = conv(two, 'G', 'A'); }
        if (i == 41) two[15] = two[15] == 'T' ? 'C' : 'T'; // strict under both ordinary and CT chemistry
        reads[0].push_back(one); reads[1].push_back(two);
    }
    const std::vector<std::string> input = {write_fq("short-rescue1.fq", reads[0]), write_fq("short-rescue2.fq", reads[1])};
    std::vector<std::string> args = {"--short-insert", "-n", "0", "-l", "1024"};
    if (meth) args.push_back("--meth");
    const auto rescued = run_aln(args, input);
    REQUIRE(rescued.size() == 84);
    const auto &a = rescued[82], &b = rescued[83];
    CHECK(a.rname == "c1"); CHECK(a.pos == 2151); CHECK(a.cigar == "20M");
    CHECK(b.rname == "c1"); CHECK(b.pos == 2171); CHECK(b.cigar == "30M");
    CHECK((a.flag & 4) == 0); CHECK((b.flag & 4) == 0);
    CHECK((a.flag & 2) != 0); CHECK((b.flag & 2) != 0); CHECK((b.flag & 16) != 0);
    CHECK(a.rnext == "="); CHECK(b.rnext == "="); CHECK(a.pnext == 2171); CHECK(b.pnext == 2151);
    CHECK(a.tlen == 50); CHECK(b.tlen == -50);
    CHECK(tag_value(b, "XT:A") == "M"); CHECK(tag_value(b, "NM:i") == "1");
    const std::string aligned = aln_test::revcomp(reads[1].back());
    CHECK(tag_value(b, "MD:Z") == expected_md(ref().c1, 2170, aligned));
    if (meth) {
        CHECK(tag_value(b, "XR:Z") == "GA"); CHECK(tag_value(b, "XG:Z") == "CT");
        CHECK(tag_value(b, "XM:Z") == expected_xm(ref().c1, 2170, aligned, true, false));
    }
    args.push_back("-s");
    const auto disabled = run_aln(args, input);
    REQUIRE(disabled.size() == 84);
    CHECK((disabled[83].flag & 4) != 0); CHECK((disabled[83].flag & 2) == 0);
    CHECK(disabled[83].cigar == "*"); CHECK(disabled[83].mapq == 0);
}

TEST_CASE("aln_meth: short-insert rescues a reverse mate whose overhang runs past the reference end"
          * doctest::test_suite("unit/aln_meth"))
{
    // The reverse-mate counterpart of the 5' overhang case: the last pair's fragment is the final
    // 35 bases of the reference (one contig, so its end is the reference end), and read 2 starts
    // with adapter, so in reference orientation its overhang runs past that end. Read 1 anchors
    // there; the rescue window from the anchor holds only the 35 fragment bases, fewer than read 2,
    // and local alignment clips the adapter.
    bool meth = false;
    SUBCASE("ordinary") {}
    SUBCASE("meth") { meth = true; }
    bwa_tests::TempDir tmp{"aln_meth_ref_end"};
    const std::string fa = tmp.unique_file();
    uint32_t seed = 752;
    const std::string contig = aln_test::lcg_dna(seed, 3000);
    REQUIRE(aln_test::write_test_fasta(fa, {{"end", contig}}));
    {
        std::vector<std::string> cmd = {"index", "--meth", fa};
        std::vector<char *> argv;
        for (auto &x : cmd) argv.push_back(&x[0]);
        argv.push_back(nullptr);
        optind = 1;
#if defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || defined(__NetBSD__)
        optreset = 1;
#endif
        REQUIRE(bwa_index((int)cmd.size(), argv.data()) == 0);
    }
    const std::string adapter = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"; // 33 bases
    const int frag = 35, pos = (int)contig.size() - frag + 1;
    std::vector<std::string> reads[2];
    for (int i = 0; i < 40; ++i) {
        const std::string insert = contig.substr(100 + 60 * i, 30 + i % 11);
        reads[0].push_back(insert); reads[1].push_back(aln_test::revcomp(insert));
    }
    const std::string fragment = contig.substr(pos - 1, frag);
    reads[0].push_back(fragment);
    reads[1].push_back(adapter + aln_test::revcomp(fragment)); // 68 bases
    const std::vector<std::string> input = {write_fq_in(tmp.path(), "refend1.fq", reads[0]),
                                            write_fq_in(tmp.path(), "refend2.fq", reads[1])};
    std::vector<std::string> args = {"-n", "0", "-l", "1024"};
    if (meth) args.push_back("--meth");
    const auto stock = run_aln(args, input, fa, tmp.path());
    REQUIRE(stock.size() == 82);
    CHECK((stock[81].flag & 4) != 0); // read 2 of the end pair stays unmapped without the option
    args.push_back("--short-insert");
    const auto rescued = run_aln(args, input, fa, tmp.path());
    REQUIRE(rescued.size() == 82);
    const auto &a = rescued[80], &b = rescued[81];
    CHECK(a.rname == "end"); CHECK(a.pos == pos); CHECK(a.cigar == std::to_string(frag) + "M");
    CHECK((b.flag & 4) == 0); CHECK((b.flag & 16) != 0);
    CHECK(b.rname == "end"); CHECK(b.pos == pos);
    CHECK(b.cigar == std::to_string(frag) + "M" + std::to_string(adapter.size()) + "S");
    CHECK((a.flag & 2) != 0); CHECK((b.flag & 2) != 0);
    CHECK(tag_value(b, "XT:A") == "M");
    // Mate fields: both start at pos and span the fragment, so the insert is frag.
    CHECK(a.rnext == "="); CHECK(b.rnext == "="); CHECK(a.pnext == pos); CHECK(b.pnext == pos);
    CHECK(a.tlen == frag); CHECK(b.tlen == -frag);
    CHECK(a.mapq > 0); CHECK(b.mapq > 0); // a unique fragment: the rescue keeps both ends confident
    // Tags: the aligned bases match the reference, so NM is 0 and MD/XM follow from the reference.
    CHECK(tag_value(a, "NM:i") == "0"); CHECK(tag_value(b, "NM:i") == "0");
    CHECK(tag_value(a, "MD:Z") == expected_md(contig, pos - 1, fragment));
    CHECK(tag_value(b, "MD:Z") == expected_md(contig, pos - 1, fragment));
    if (meth) {
        CHECK(tag_value(a, "XR:Z") == "CT"); CHECK(tag_value(a, "XG:Z") == "CT");
        CHECK(tag_value(b, "XR:Z") == "GA"); CHECK(tag_value(b, "XG:Z") == "CT");
        CHECK(tag_value(a, "XM:Z") == expected_xm(contig, pos - 1, fragment, true, false));
        const std::string xm_b = tag_value(b, "XM:Z"); // the clipped adapter makes no calls
        REQUIRE(xm_b.size() == fragment.size() + adapter.size());
        CHECK(xm_b.substr(0, fragment.size()) == expected_xm(contig, pos - 1, fragment, true, false));
        CHECK(xm_b.substr(fragment.size()) == std::string(adapter.size(), '.'));
    }
}

TEST_CASE("aln_meth: short-insert rescues a mate whose 5' overhang runs past its anchor"
          * doctest::test_suite("unit/aln_meth"))
{
    // Read 1 starts with adapter and then covers the whole fragment, which read 2 also covers.
    // End-to-end search cannot place read 1; rescue can, once its window may cover the anchor's
    // own span, and the clipped adapter leaves the alignment inside the fragment.
    bool meth = false, converted = false;
    SUBCASE("ordinary") {}
    SUBCASE("meth exact") { meth = true; }
    SUBCASE("meth converted") { meth = converted = true; }
    const std::string adapter = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT";
    std::vector<std::string> reads[2];
    for (int i = 0; i < 41; ++i) {
        const int start = 100 + 50 * i, len = 30 + i % 11;
        std::string insert = ref().c1.substr(start, len);
        std::string one = converted ? conv(insert, 'C', 'T') : insert;
        std::string two = aln_test::revcomp(insert);
        if (converted) two = conv(two, 'G', 'A');
        if (i == 40) one = adapter + one;
        reads[0].push_back(one); reads[1].push_back(two);
    }
    const std::vector<std::string> input = {write_fq("overhang1.fq", reads[0]), write_fq("overhang2.fq", reads[1])};
    std::vector<std::string> args = {"-n", "0", "-l", "1024"};
    if (meth) args.push_back("--meth");
    const auto stock = run_aln(args, input);
    REQUIRE(stock.size() == 82);
    CHECK((stock[80].flag & 4) != 0); // read 1 of the adapter pair stays unmapped without the option
    args.push_back("--short-insert");
    const auto rescued = run_aln(args, input);
    REQUIRE(rescued.size() == 82);
    const auto &a = rescued[80], &b = rescued[81];
    const int pos = 100 + 50 * 40 + 1, len = 30 + 40 % 11;
    CHECK(b.rname == "c1"); CHECK(b.pos == pos); CHECK(b.cigar == std::to_string(len) + "M");
    CHECK((a.flag & 4) == 0); CHECK(a.rname == "c1"); CHECK(a.pos == pos);
    CHECK(a.cigar == std::to_string(adapter.size()) + "S" + std::to_string(len) + "M");
    CHECK((a.flag & 2) != 0); CHECK((b.flag & 2) != 0);
    CHECK(tag_value(a, "XT:A") == "M");
    // Mate fields: both start at pos and span len reference bases, so the insert is len.
    CHECK(a.rnext == "="); CHECK(b.rnext == "="); CHECK(a.pnext == pos); CHECK(b.pnext == pos);
    CHECK(a.tlen == len); CHECK(b.tlen == -len);
    CHECK(a.mapq > 0); CHECK(b.mapq > 0); // a unique fragment: the rescue keeps both ends confident
    // Tags: the aligned bases match the reference (conversions are free under --meth), so NM is 0
    // and MD/XM follow from the reference and the aligned segment of each read.
    const std::string aligned_a = reads[0].back().substr(adapter.size()), aligned_b = aln_test::revcomp(reads[1].back());
    CHECK(tag_value(a, "NM:i") == "0"); CHECK(tag_value(b, "NM:i") == "0");
    CHECK(tag_value(a, "MD:Z") == expected_md(ref().c1, pos - 1, aligned_a));
    CHECK(tag_value(b, "MD:Z") == expected_md(ref().c1, pos - 1, aligned_b));
    if (meth) {
        CHECK(tag_value(a, "XR:Z") == "CT"); CHECK(tag_value(a, "XG:Z") == "CT");
        CHECK(tag_value(b, "XR:Z") == "GA"); CHECK(tag_value(b, "XG:Z") == "CT");
        CHECK(tag_value(b, "XM:Z") == expected_xm(ref().c1, pos - 1, aligned_b, true, false));
        const std::string xm_a = tag_value(a, "XM:Z"); // the clipped adapter makes no calls
        REQUIRE(xm_a.size() == adapter.size() + aligned_a.size());
        CHECK(xm_a.substr(0, adapter.size()) == std::string(adapter.size(), '.'));
        CHECK(xm_a.substr(adapter.size()) == expected_xm(ref().c1, pos - 1, aligned_a, true, false));
    }
    for (int i = 0; i < 40; ++i) { // the coincident starts of the other pairs are proper too
        CHECK((rescued[2 * i].flag & 2) != 0); CHECK((rescued[2 * i + 1].flag & 2) != 0);
    }
}

TEST_CASE("aln_meth: short-insert cap rescue skips nonrepresentable likelihoods without losing the batch"
          * doctest::test_suite("unit/aln_meth"))
{
    // Narrow calibration plus a lookup-capped 70bp mate outside its likelihood support. The mate
    // is poly-T, so a 4,700-base poly-A run gives it more placements than the default --meth lookup
    // budget (p_best + p_second = 4,596 rows) while the genome stays a few kilobases.
    // The control calibrates on 90-91bp inserts instead, so the same anchor and capped mate form a
    // pair with a finite penalty: cap rescue must then install the mate, which shows the rejected
    // case evaluates the candidate rather than skipping it.
    bool finite = false;
    SUBCASE("unrepresentable likelihood") {}
    SUBCASE("finite-penalty control") { finite = true; }
    bwa_tests::TempDir tmp{"aln_meth_short_cap"}; // this case's index, reads and output, removed when it ends
    const std::string fa = tmp.path() + "/short-cap.fa";
    uint32_t state = 77;
    std::string genome = aln_test::lcg_dna(state, 5000);
    const std::string repeat(70, 'A');
    genome.replace(120, 70, repeat);
    genome += std::string(4700, 'A') + aln_test::lcg_dna(state, 30);
    REQUIRE(aln_test::write_test_fasta(fa, {{"cap", genome}}));
    std::vector<std::string> index_args = {"index", "--meth", fa};
    std::vector<char *> argv;
    for (auto &arg : index_args) argv.push_back(&arg[0]);
    argv.push_back(nullptr);
    optind = 1;
#if defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || defined(__NetBSD__)
    optreset = 1;
#endif
    REQUIRE(bwa_index((int)index_args.size(), argv.data()) == 0);
    std::vector<std::string> reads[2];
    for (int i = 0; i < 22; ++i) {
        std::string one, two;
        if (i < 20) {
            const int start = 1000 + 200 * i;
            const int shift = (finite ? 40 : 5) + i % 2;
            one = genome.substr(start, 50); two = aln_test::revcomp(genome.substr(start + shift, 50));
        } else if (i == 20) { one = genome.substr(100, 20); two = aln_test::revcomp(repeat); }
        else one = two = std::string(150, 'N');
        reads[0].push_back(one); reads[1].push_back(two);
    }
    const std::vector<std::string> input = {write_fq_in(tmp.path(), "short-cap1.fq", reads[0]),
                                            write_fq_in(tmp.path(), "short-cap2.fq", reads[1])};
    const auto output = run_aln({"--meth", "--short-insert", "-n", "0", "-l", "1024"}, input, fa, tmp.path());
    REQUIRE(output.size() == 44);
    for (int i = 0; i < 20; ++i) for (int j = 0; j < 2; ++j) {
        const auto &r = output[2*i+j];
        CHECK(r.rname == "cap"); CHECK(r.pos == 1001 + 200*i + (j ? (finite ? 40 : 5)+i%2 : 0));
        CHECK(r.cigar == "50M"); CHECK((r.flag & 2) != 0);
    }
    const auto &anchor = output[40], &mate = output[41];
    CHECK(anchor.pos == 101); CHECK(anchor.cigar == "20M"); CHECK((anchor.flag & 2) == 0);
    CHECK((mate.flag & 2) == 0); CHECK(mate.mapq == 0);
    if (finite) {
        // The poly-A run at offset 120 is the only exact placement inside the fragment.
        CHECK(tag_value(mate, "XT:A") == "M");
        CHECK(mate.rname == "cap"); CHECK(mate.pos == 121);
        CHECK(mate.cigar == "70M"); CHECK((mate.flag & 16) != 0);
    } else {
        CHECK((mate.tag.count("XT:A") == 0 || mate.tag.at("XT:A") != "M"));
    }
    CHECK((output[42].flag & 4) != 0); CHECK((output[43].flag & 4) != 0);
}

TEST_CASE("aln_meth: short-insert active estimates retain neutral ceilings and incomplete restrictions"
          * doctest::test_suite("unit/aln_meth"))
{
    namespace ba = bwa_aln;
    FMI_search fmi((ref().fa + ".meth").c_str()); fmi.load_index(false, 1);
    ba::aln_meth_ref_t mr;
    REQUIRE(ba::aln_meth_ref_load(&mr, ref().fa.c_str(), fmi.idx->bns, 1) == 0);
    ba::aln_meth_opt_t mo;
    mo.m.meth_chem = METH_CHEM_TAPS; mo.m.meth_scoring = MEM_METH_SCORING_NEUTRAL;
    ba::aln_meth_opt_finish(&mo);
    auto *gap = ba::gap_init_opt(); auto *pe = ba::bwa_init_pe_opt();
    pe->short_insert = 1;
    ba::bwase_initialize();
    ba::bwa_seq_t calibration[2][41] = {};
    for (int j = 0; j < 2; ++j) for (int i = 0; i < 41; ++i) {
        auto &p = calibration[j][i]; p.mapQ = 37;
        p.len = p.full_len = i == 40 ? 150 : 30 + i % 16;
        p.pos = 100 + 50 * i + j * 5;
    }
    ba::bwa_seq_t *cal[2] = {calibration[0], calibration[1]};
    ba::isize_info_t ii = {}, none = {}; none.avg = -1;
    ba::bwa_pe_isize_ref(mr.bns->l_pac, 41, cal, &ii, pe, &none);
    REQUIRE(ii.avg > 0); REQUIRE(ii.high > 0); REQUIRE(ii.high_bayesian > 45);
    for (int incomplete : {0, (int)ba::ALN_METH_BIND_M, (int)ba::ALN_METH_BIND_BEST}) {
        CAPTURE(incomplete);
        ba::bwa_seq_t *seqs[2]; ba::bwa_seqio_t *input[2]; int n[2];
        for (int j = 0; j < 2; ++j) {
            const std::string read = j ? aln_test::revcomp(ref().c1.substr(105, 35)) : ref().c1.substr(100, 35);
            input[j] = ba::bwa_seq_open(write_fq(j ? "short-guard2.fq" : "short-guard1.fq", {read}).c_str());
            seqs[j] = ba::bwa_read_seq(input[j], 1, &n[j], BWA_MODE_COMPREAD, 0);
            REQUIRE(n[j] == 1);
        }
        ba::aln_meth_read_t result[2]; ba::aln_meth_read_t *r[2] = {result, result + 1};
        for (int j = 0; j < 2; ++j) {
            ba::aln_meth_occ_t o; o.rid = 0; o.hyp = 1; o.strand = j; o.pos = 100 + j * 5;
            o.alignment.rb = o.pos; o.alignment.re = o.pos + 35; o.alignment.cigar = {35u << 4};
            result[j].occ.push_back(o); result[j].main = 0; result[j].c1 = 1;
        }
        result[0].flags = incomplete;
        store_pair(r); ba::aln_meth_pe_stage1(gap, 1, seqs, r);
        // A neutral conversion-selected main locus must not gain independent evidence from its mate.
        result[0].n_conv_moved = 1; result[0].conv_gap = 2;
        if (!incomplete) seqs[0]->mapQ = 12;
        uint64_t total[2] = {}, mapped[2] = {};
        ba::aln_meth_pe_stage2(&mr, &mo, gap, pe, &ii, 1, seqs, r, total, mapped);
        CHECK(seqs[0]->pos == 100); CHECK(seqs[1]->pos == 105);
        CHECK(result[0].main == 0); CHECK(result[0].hyp == 1); CHECK(result[1].hyp == 1);
        if (incomplete) {
            CHECK((int)seqs[0]->mapQ == 0); CHECK((int)seqs[0]->c1 == 0);
            CHECK((seqs[0]->extra_flag & SAM_FPP) == 0);
            CHECK((seqs[1]->extra_flag & SAM_FPP) == 0);
            // Cap rescue may now search the overlapping mate's span and re-place a best-capped end
            // there, still at MAPQ 0 and unpaired; an end capped by -m is never rescued.
            CHECK(((int)seqs[0]->type == BWA_TYPE_MATESW) == (incomplete == (int)ba::ALN_METH_BIND_BEST));
        } else {
            CHECK((int)seqs[0]->mapQ == 12);
            CHECK((seqs[0]->extra_flag & SAM_FPP) != 0);
            CHECK((seqs[1]->extra_flag & SAM_FPP) != 0);
        }
        for (int j = 0; j < 2; ++j) { ba::bwa_free_read_seq(1, seqs[j]); ba::bwa_seq_close(input[j]); }
    }
    free(gap); free(pe); ba::aln_meth_ref_free(&mr);
}

TEST_CASE("aln_meth: rescue windows belong to the anchor contig at both boundaries")
{
    namespace ba=bwa_aln;
    FMI_search fmi((ref().fa+".meth").c_str());fmi.load_index(false,1);
    ba::aln_meth_ref_t mr;REQUIRE(ba::aln_meth_ref_load(&mr,ref().fa.c_str(),fmi.idx->bns,1)==0);
    ba::aln_meth_opt_t mo;ba::aln_meth_opt_finish(&mo);
    auto *opt=ba::gap_init_opt();auto *pe=ba::bwa_init_pe_opt();ba::bwase_initialize();
    bool conversion_guard=false,cap_bound=false,rescued_moved=false;
    SUBCASE("ordinary anchors") {}
    SUBCASE("cap-bound no-hit ends use bounded rescue without proper-pair claims") {cap_bound=true;}
    SUBCASE("conversion-selected anchors retain their ceiling") {conversion_guard=true;}
    SUBCASE("a rescued end leaves its conversion-selected locus and ceiling") {rescued_moved=true;}
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
        ba::aln_meth_occ_t o{};o.rid=c.rid;o.strand=c.strand;o.hyp=1-c.strand;
        o.pos=mr.bns->anns[c.rid].offset+c.anchor;o.alignment.rb=o.pos;o.alignment.re=o.pos+50;
        o.alignment.cigar={50u<<4};result[0].occ.push_back(o);result[0].main=0;result[0].c1=1;
        if(cap_bound) {
            result[1].flags=ba::ALN_METH_BIND_BEST;
            result[1].hits.reserve(64); // spare capacity
        }
        if(conversion_guard) {result[0].n_conv_moved=1;result[0].conv_gap=3;result[0].c2=1;}
        if(rescued_moved) {
            // The mate's own conversion-selected locus (ceiling min(60,6*1) = 6) is unpaired:
            // same strand as the anchor, well away from the rescue window.
            ba::aln_meth_occ_t m=o;m.rid=0;m.pos=mr.bns->anns[0].offset+1500;
            m.alignment.rb=m.pos;m.alignment.re=m.pos+50;
            result[1].occ.push_back(m);result[1].main=0;result[1].c1=1;result[1].c2=1;
            result[1].n_conv_moved=1;result[1].conv_gap=1;
        }
        store_pair(r);ba::aln_meth_pe_stage1(opt,1,seqs,r);
        if(conversion_guard) CHECK((int)seqs[0]->mapQ==18);
        if(rescued_moved) {CHECK((int)seqs[0]->mapQ==37);CHECK((int)seqs[1]->mapQ==6);}
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
                REQUIRE(result[1].hits.size()==1);CHECK(result[1].hits[0].cost==0);
                CHECK(result[1].hits.capacity()==1); // the rescue releases the old capacity
                CHECK(cigar_of(result[1],result[1].hits[0])==std::vector<uint32_t>{50u<<4});
            }
            if(conversion_guard) {CHECK((int)seqs[0]->mapQ==10);CHECK((int)seqs[1]->mapQ==10);}
            // Rescue moved the mate off that locus, so its MAPQ is the pair's, 37 - 6/2 - 8 = 26,
            // not the old locus's ceiling of 6.
            if(rescued_moved) {CHECK((int)seqs[0]->mapQ==26);CHECK((int)seqs[1]->mapQ==26);}
            if(seqs[1]->md) CHECK(std::string(seqs[1]->md)=="50");
            else CHECK(seqs[1]->md!=nullptr);
        } else if(!rescued_moved) {
            CHECK((int)seqs[1]->type==BWA_TYPE_NO_MATCH);
            CHECK((seqs[0]->extra_flag&SAM_FPP)==0);CHECK((seqs[1]->extra_flag&SAM_FPP)==0);
        } else {
            CHECK((int)seqs[1]->type!=BWA_TYPE_MATESW);CHECK((int)seqs[1]->mapQ<=6);
        }
        for(int j=0;j<2;++j){ba::bwa_free_read_seq(1,seqs[j]);ba::bwa_seq_close(inputs[j]);}
    }
    free(pe);free(opt);ba::aln_meth_ref_free(&mr);
}

TEST_CASE("aln_meth: an incomplete end's worse-cost loci never reach paired SAM output"
          * doctest::test_suite("unit/aln_meth"))
{
    // An incomplete end kept for pairing keeps only its main hit (and, while duplicate copies are
    // left to draw, its other best-cost loci). The records must be the same as when every locus was
    // kept: pairing, mate rescue, cap rescue and XA must never read past the main hit. The worse
    // locus shares the anchor's hypothesis at a proper-pair distance, so pairing takes it when the
    // end is complete (the positive control below). Cap rescue reads only the end's main-hit cost
    // as a bound and searches the window, where the mate aligns at 1010, not at 940.
    namespace ba=bwa_aln;
    bwa_tests::TempDir tmp{"aln_meth_worse_loci"}; // this case's reads, removed when it ends
    FMI_search fmi((ref().fa+".meth").c_str());fmi.load_index(false,1);
    ba::aln_meth_ref_t mr;REQUIRE(ba::aln_meth_ref_load(&mr,ref().fa.c_str(),fmi.idx->bns,1)==0);
    ba::aln_meth_opt_t mo;ba::aln_meth_opt_finish(&mo);
    auto *opt=ba::gap_init_opt();auto *pe=ba::bwa_init_pe_opt();ba::bwase_initialize();
    // Read 1 anchors on c2 (reverse, at 1110). Read 2 is c2[1010, 1060), forward: its true
    // locus, which cap rescue finds. The worse locus at c2+940 is also a proper-pair distance
    // away, but not where the read aligns, so only a reader of the dropped locus can report it.
    const std::string mate=ref().c2.substr(1010,50);
    const int64_t c1_off=mr.bns->anns[0].offset,c2_off=mr.bns->anns[1].offset;
    auto locus=[](int rid,int strand,int64_t pos,int64_t cost) {
        ba::aln_meth_occ_t o;o.rid=rid;o.strand=strand;o.hyp=1-strand;o.pos=pos;o.cost=cost;
        o.alignment.rb=pos;o.alignment.re=pos+50;o.alignment.cigar={50u<<4};
        o.n_mm=(int)(cost/3);
        return o;
    };
    for(int flags:{(int)ba::ALN_METH_BIND_BEST,(int)ba::ALN_METH_BIND_M}) for(int drawn:{0,1}) {
        CAPTURE(flags);CAPTURE(drawn);
        // filtered: 0 keeps every locus, 1 only the best-cost ones (a read with duplicate copies
        // left to draw), 2 only the main hit (what pairing gets).
        std::string sam[3];
        for(int filtered=0;filtered<3;++filtered) {
            ba::bwa_seq_t *seqs[2];ba::bwa_seqio_t *inputs[2];int n[2];
            for(int j=0;j<2;++j) {
                inputs[j]=ba::bwa_seq_open(write_fq_in(tmp.path(),j?"worse2.fq":"worse1.fq",{j?mate:std::string(50,'A')}).c_str());
                seqs[j]=ba::bwa_read_seq(inputs[j],1,&n[j],BWA_MODE_COMPREAD,0);REQUIRE(n[j]==1);
            }
            ba::aln_meth_read_t result[2];ba::aln_meth_read_t *r[2]={result,result+1};
            result[0].occ.push_back(locus(1,1,c2_off+1110,0));result[0].main=0;result[0].c1=1;
            // The incomplete end, in canonical order: two best-cost loci on c1 and, unless
            // filtered, a worse one on c2 beside the anchor. c1/c2 count all three either way.
            auto &rm=result[1];
            if(filtered<2 || !drawn) rm.occ.push_back(locus(0,0,c1_off+100,0));
            if(filtered<2 || drawn) rm.occ.push_back(locus(0,0,c1_off+600,0));
            if(!filtered) rm.occ.push_back(locus(1,0,c2_off+940,3));
            for(auto &o:rm.occ) o.hyp=result[0].occ[0].hyp; // pairable with the anchor
            rm.main=filtered==2?0:drawn;rm.c1=2;rm.c2=1;rm.flags=flags;
            store_pair(r);ba::aln_meth_pe_stage1(opt,1,seqs,r);
            ba::isize_info_t ii={};ii.avg=190;ii.std=10;ii.ap_prior=1e-5;ii.low=100;ii.high=250;ii.high_bayesian=250;
            uint64_t total[2]={},mapped[2]={};ba::aln_meth_pe_stage2(&mr,&mo,opt,pe,&ii,1,seqs,r,total,mapped);
            kstring_t text={0,0,nullptr};
            ba::aln_meth_print_sam1(&text,&mr,&mo,seqs[0],result,1,seqs[1],BWA_MODE_COMPREAD,30);
            ba::aln_meth_print_sam1(&text,&mr,&mo,seqs[1],result+1,0,seqs[0],BWA_MODE_COMPREAD,30);
            sam[filtered].assign(text.s,text.l);free(text.s);
            // An incomplete end is never a proper pair, and never placed at the dropped locus.
            CHECK((seqs[1]->extra_flag&SAM_FPP)==0);
            CHECK(seqs[1]->pos!=(bwtint_t)(c2_off+940));
            // The test exercises what it claims: cap rescue places a BIND_BEST end at its true
            // locus; a BIND_M end is not cap-rescued and keeps its drawn best-cost locus.
            if(flags==ba::ALN_METH_BIND_BEST) {
                CHECK(seqs[1]->pos==(bwtint_t)(c2_off+1010));CHECK((int)seqs[1]->type==BWA_TYPE_MATESW);
            } else {
                CHECK(seqs[1]->pos==(bwtint_t)(c1_off+(drawn?600:100)));
            }
            for(int j=0;j<2;++j){ba::bwa_free_read_seq(1,seqs[j]);ba::bwa_seq_close(inputs[j]);}
        }
        CHECK(sam[0]==sam[1]); CHECK(sam[0]==sam[2]);
        CHECK(sam[0].find("XA:Z:")==std::string::npos);
        // Independent oracle for both records. Read 1 is poly-A placed reverse at c2:1110, so it
        // prints as poly-T; both ends take the anchor's GA genome strand, where a reference G read
        // as A is a conversion, not an edit. A BIND_BEST end is cap-rescued to its true locus
        // (XT:A:M); a BIND_M end keeps its drawn locus as a repeat (XT:A:R). Neither incomplete
        // end has a unique-hit count, so only the anchor carries X0/X1.
        const bool rescued=flags==ba::ALN_METH_BIND_BEST;
        const std::string &mate_genome=rescued?ref().c2:ref().c1;
        const int mate_start=rescued?1010:(drawn?600:100);
        const std::string anchor_read(50,'T'),quals(50,'I');
        auto edits=[](const std::string &g,int at,const std::string &a) {
            int nm=0;
            for(size_t i=0;i<a.size();++i) nm+=a[i]!=g[at+i] && !(g[at+i]=='G' && a[i]=='A');
            return std::to_string(nm);
        };
        const std::string mate_rname=rescued?"c2":"c1",anchor_next=rescued?"=":"c1";
        const std::string want=
            "r0\t81\tc2\t1111\t37\t50M\t"+anchor_next+"\t"+std::to_string(mate_start+1)+"\t"+(rescued?"-150":"0")+
            "\t"+anchor_read+"\t"+quals+"\tXT:A:U\tNM:i:"+edits(ref().c2,1110,anchor_read)+
            "\tSM:i:37\tAM:i:0\tX0:i:1\tX1:i:0\tMD:Z:"+expected_md(ref().c2,1110,anchor_read)+
            "\tXR:Z:CT\tXG:Z:GA\tXM:Z:"+expected_xm(ref().c2,1110,anchor_read,false,false)+"\n"+
            "r0\t161\t"+mate_rname+"\t"+std::to_string(mate_start+1)+"\t0\t50M\t"+(rescued?"=":"c2")+"\t1111\t"+
            (rescued?"150":"0")+"\t"+mate+"\t"+quals+"\tXT:A:"+(rescued?"M":"R")+
            "\tNM:i:"+edits(mate_genome,mate_start,mate)+"\tSM:i:0\tAM:i:0\tMD:Z:"+
            expected_md(mate_genome,mate_start,mate)+"\tXR:Z:GA\tXG:Z:GA\tXM:Z:"+
            expected_xm(mate_genome,mate_start,mate,false,false)+"\n";
        CHECK(sam[0]==want);
    }
    {
        // Positive control: the same end, complete, with every locus kept, pairs at the worse locus,
        // so the incomplete cases above would report it if anything read past the main hit.
        ba::bwa_seq_t *seqs[2];ba::bwa_seqio_t *inputs[2];int n[2];
        for(int j=0;j<2;++j) {
            inputs[j]=ba::bwa_seq_open(write_fq_in(tmp.path(),j?"worse2.fq":"worse1.fq",{j?mate:std::string(50,'A')}).c_str());
            seqs[j]=ba::bwa_read_seq(inputs[j],1,&n[j],BWA_MODE_COMPREAD,0);REQUIRE(n[j]==1);
        }
        ba::aln_meth_read_t result[2];ba::aln_meth_read_t *r[2]={result,result+1};
        result[0].occ.push_back(locus(1,1,c2_off+1110,0));result[0].main=0;result[0].c1=1;
        auto &rm=result[1];
        rm.occ.push_back(locus(0,0,c1_off+100,0));
        rm.occ.push_back(locus(0,0,c1_off+600,0));
        rm.occ.push_back(locus(1,0,c2_off+940,3));
        for(auto &o:rm.occ) o.hyp=result[0].occ[0].hyp;
        rm.main=0;rm.c1=2;rm.c2=1;rm.flags=0;
        store_pair(r);ba::aln_meth_pe_stage1(opt,1,seqs,r);
        ba::isize_info_t ii={};ii.avg=190;ii.std=10;ii.ap_prior=1e-5;ii.low=100;ii.high=250;ii.high_bayesian=250;
        uint64_t total[2]={},mapped[2]={};ba::aln_meth_pe_stage2(&mr,&mo,opt,pe,&ii,1,seqs,r,total,mapped);
        CAPTURE(seqs[1]->pos-c2_off);
        CHECK(seqs[1]->pos==(bwtint_t)(c2_off+940));
        CHECK((seqs[1]->extra_flag&SAM_FPP)!=0);
        for(int j=0;j<2;++j){ba::bwa_free_read_seq(1,seqs[j]);ba::bwa_seq_close(inputs[j]);}
    }
    // Cap rescue bounds the candidate by the main hit's cost, not the first locus's. The mate now
    // has one edit that is no conversion on either strand (a reference A or T read as C), so its
    // window candidate at c2:1010 costs one mismatch. With a cost-0 main hit after a cost-3 locus
    // that candidate is worse than the bound and must be rejected; with the main hit at the
    // candidate's own cost (the control) it is installed, which pins the candidate's cost.
    std::string edited=mate;
    const size_t edit=edited.find_first_of("AT",10);
    REQUIRE(edit!=std::string::npos);edited[edit]='C';
    for(bool control:{false,true}) {
        CAPTURE(control);
        ba::bwa_seq_t *seqs[2];ba::bwa_seqio_t *inputs[2];int n[2];
        for(int j=0;j<2;++j) {
            inputs[j]=ba::bwa_seq_open(write_fq_in(tmp.path(),j?"bound2.fq":"bound1.fq",{j?edited:std::string(50,'A')}).c_str());
            seqs[j]=ba::bwa_read_seq(inputs[j],1,&n[j],BWA_MODE_COMPREAD,0);REQUIRE(n[j]==1);
        }
        ba::aln_meth_read_t result[2];ba::aln_meth_read_t *r[2]={result,result+1};
        result[0].occ.push_back(locus(1,1,c2_off+1110,0));result[0].main=0;result[0].c1=1;
        auto &rm=result[1];
        if(control) {rm.occ.push_back(locus(0,0,c1_off+100,3));rm.main=0;rm.c1=1;rm.c2=0;}
        else {
            rm.occ.push_back(locus(0,0,c1_off+600,3));rm.occ.push_back(locus(0,0,c1_off+100,0));
            rm.main=1;rm.c1=1;rm.c2=1;
        }
        for(auto &o:rm.occ) o.hyp=result[0].occ[0].hyp;
        rm.flags=ba::ALN_METH_BIND_BEST;
        store_pair(r);ba::aln_meth_pe_stage1(opt,1,seqs,r);
        ba::isize_info_t ii={};ii.avg=190;ii.std=10;ii.ap_prior=1e-5;ii.low=100;ii.high=250;ii.high_bayesian=250;
        uint64_t total[2]={},mapped[2]={};ba::aln_meth_pe_stage2(&mr,&mo,opt,pe,&ii,1,seqs,r,total,mapped);
        if(control) {
            CHECK(seqs[1]->pos==(bwtint_t)(c2_off+1010));CHECK((int)seqs[1]->type==BWA_TYPE_MATESW);
            CHECK(seqs[1]->score==3);CHECK((int)seqs[1]->n_mm==1);
        } else {
            CHECK(seqs[1]->pos==(bwtint_t)(c1_off+100));CHECK((int)seqs[1]->type!=BWA_TYPE_MATESW);
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

TEST_CASE("aln_meth: a pair whose insert penalty underflows is not scored")
{
    // A very narrow estimate (std sqrt(.2)) puts a 100 bp insert about 222 standard deviations
    // below the mean: erfc underflows to 0 and the penalty is +inf, which cannot be converted to
    // an int. The candidate is skipped, as --meth cap rescue skips it, so the reads stay unpaired.
    namespace ba=bwa_aln;
    FMI_search fmi((ref().fa+".meth").c_str());fmi.load_index(false,1);
    ba::aln_meth_ref_t mr;REQUIRE(ba::aln_meth_ref_load(&mr,ref().fa.c_str(),fmi.idx->bns,1)==0);
    ba::bwase_initialize();auto *pe=ba::bwa_init_pe_opt();
    for(const bool tail:{true,false}) {
        CAPTURE(tail);
        ba::bwa_seq_t reads[2]={};ba::bwa_seq_t *p[2]={reads,reads+1};
        ba::bwt_aln1_t hits[2][1]={};
        for(int j=0;j<2;++j){p[j]->len=p[j]->full_len=50;p[j]->type=BWA_TYPE_UNIQUE;
            p[j]->mapQ=p[j]->seQ=37;p[j]->aln=hits[j];}
        p[0]->pos=100;p[0]->strand=0;p[1]->pos=150;p[1]->strand=1;
        std::vector<pair64_t> rows={{100,0},{150,3}};
        pair64_v v;v.a=rows.data();v.n=v.m=rows.size();
        ba::isize_info_t ii={};ii.avg=tail?199.5:100.0;ii.std=std::sqrt(.2);
        ii.low=1;ii.high=300;ii.high_bayesian=300;
        ba::bwa_pe_pairing_same_contig(mr.bns,p,&v,pe,3,&ii);
        // The control (insert at the mean) pairs, so the skip is what leaves the tail case unpaired.
        CHECK(((p[0]->extra_flag&SAM_FPP)!=0)==!tail);
        CHECK(((p[1]->extra_flag&SAM_FPP)!=0)==!tail);
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
        for(int j=0;j<2;++j){ba::aln_meth_occ_t o{};o.rid=j;o.hyp=1;o.strand=j;o.pos=j?3100:2800;
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
        store_pair(r);ba::aln_meth_pe_stage1(opt,1,seqs,r);ba::isize_info_t ii={};ii.avg=ii.std=-1;
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
    const std::string planted=query;
    query.insert(30,"AAA");
    // The full planted alignment, no clipping: aM3IbM with a+b==60, where removing the three bases
    // at a gives back the planted read (a shifted but equivalent insertion is allowed).
    std::set<std::string> expected_cigars;
    for(int a=0;a<=60;++a)
        if(query.substr(0,a)+query.substr(a+3)==planted)
            expected_cigars.insert(std::to_string(a)+"M3I"+std::to_string(60-a)+"M");
    REQUIRE(expected_cigars.count("30M3I30M")==1);
    const auto hits=run_aln({"--meth=taps","--meth-context","cpg","-n","5","-l","1024","-o","1"},
                            {write_fq("cpg-insertion.fq",{query})});
    REQUIRE(hits.size()==1);CHECK(hits[0].rname=="c1");CHECK(hits[0].pos==start+1);
    CAPTURE(hits[0].cigar);
    REQUIRE(expected_cigars.count(hits[0].cigar)==1);
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
    bwa_tests::TempDir tmp{"aln_meth_context"}; // this case's reference, index and reads
    const std::string fa=tmp.path()+"/context-boundaries.fa";
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
                            {write_fq_in(tmp.path(),"context-boundaries.fq",{l,r})},fa);
    REQUIRE(edge.size()==2);
    CHECK(edge[0].rname=="left");CHECK(edge[1].rname=="right");
    CHECK(tag_value(edge[0],"NM:i")=="1");CHECK(tag_value(edge[1],"NM:i")=="1");
    CHECK(tag_value(edge[0],"MD:Z")=="59C0");CHECK(tag_value(edge[1],"MD:Z")=="0G59");
    CHECK(tag_value(edge[0],"XM:Z").back()=='U');
    CHECK(tag_value(edge[1],"XM:Z").front()=='U');
    const auto n=run_aln({"--meth=taps","--meth-context","cpg","-n","4","-l","1024"},
                         {write_fq_in(tmp.path(),"context-ambiguity.fq",{ambiguous})},fa);
    REQUIRE(n.size()==1);CHECK(n[0].rname=="ambiguous");
    CHECK(tag_value(n[0],"NM:i")=="1");
    const auto &xm=tag_value(n[0],"XM:Z");REQUIRE(xm.size()==80);
    CHECK(xm[29]=='u');CHECK(xm[30]=='.');
    CHECK(xm==expected_xm(ambiguous,0,ambiguous,true,true));
    // Meth MD reports original ambiguity, independently of the index's random PAC fill.
    CHECK(tag_value(n[0],"MD:Z")=="30N49");
    // A deletion over that base also shows it as N, not as the packed fill.
    std::string deleted=ambiguous;deleted.erase(30,1);
    const auto d=run_aln({"--meth=taps","--meth-context","cpg","-n","4","-l","1024"},
                         {write_fq("context-ambiguity-deletion.fq",{deleted})},fa);
    REQUIRE(d.size()==1);CHECK(d[0].rname=="ambiguous");CHECK(d[0].pos==1);
    CHECK(d[0].cigar=="30M1D49M");
    CHECK(tag_value(d[0],"MD:Z")=="30^N49");CHECK(tag_value(d[0],"NM:i")=="1");
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


TEST_CASE("aln_meth: incoming window flags precede neutral selection and duplicate sharing")
{
    namespace ba = bwa_aln;
    FMI_search fmi((ref().fa + ".meth").c_str()); fmi.load_index(false, 1, false, BWAMEM_HUGETLB_OFF, false);
    auto *bwt = ba::bwt_from_fmi(&fmi); ba::aln_meth_ref_t mr;
    REQUIRE(ba::aln_meth_ref_load(&mr, ref().fa.c_str(), fmi.idx->bns, 1) == 0);
    auto *opt = ba::gap_init_opt(); opt->fnr = -1; opt->max_diff = 2; opt->seed_len = 1024;
    ba::aln_meth_opt_t mo; mo.m.meth_chem = METH_CHEM_TAPS; mo.m.meth_scoring = MEM_METH_SCORING_NEUTRAL;
    ba::aln_meth_opt_finish(&mo);
    const auto fq = write_fq("incoming-flags.fq", std::vector<std::string>(8, ref().read_u));
    auto *input = ba::bwa_seq_open(fq.c_str()); int n = 0;
    auto *reads = ba::bwa_read_seq(input, 20, &n, BWA_MODE_COMPREAD, 0); REQUIRE(n == 8);
    std::vector<int> flags{0, ba::ALN_METH_BIND_M, 0, ba::ALN_METH_BIND_BEST,
        ba::ALN_METH_BIND_M, ba::ALN_METH_BIND_CIGAR, ba::ALN_METH_BIND_UNREP, 0};
    std::vector<ba::aln_meth_read_t> results(n); ba::aln_dup_state_t dup;
    CHECK(ba::aln_meth_search_batch(bwt, &mr, &mo, n, reads, opt, 1, 1, 0, 1, 1000,
          &dup, results.data(), nullptr, nullptr, false, flags.data()) == 3);
    for (int i = 0; i < n; ++i) {
        CHECK(results[i].flags == flags[i]);
        if (flags[i]) { CHECK(results[i].n_conv_moved == 0); CHECK(results[i].conv_gap == 0); }
        else CHECK(results[i].n_conv_moved == 3);
        ba::aln_meth_read_t single;
        ba::aln_meth_search_batch(bwt, &mr, &mo, 1, reads + i, opt, 1, 1, i, 1, 1000,
            nullptr, &single, nullptr, nullptr, false, &flags[i]);
        CHECK(results[i].main == single.main); CHECK(results[i].c1 == single.c1); CHECK(results[i].c2 == single.c2);
    }
    // Reset with no incoming flags restores complete-state tie-breaking for every reused result.
    dup = ba::aln_dup_state_t();
    ba::aln_meth_search_batch(bwt, &mr, &mo, n, reads, opt, 1, 1, 0, 1, 1000, &dup, results.data());
    for (const auto &r : results) { CHECK(r.flags == 0); CHECK(r.n_conv_moved == 3); }
    ba::bwa_free_read_seq(n, reads); ba::bwa_seq_close(input);
    free(opt); ba::aln_meth_ref_free(&mr); ba::bwt_destroy(bwt);
}

TEST_CASE("aln_meth: duplicate reads share a search only with the same candidate paths"
          * doctest::test_suite("unit/aln_meth"))
{
    namespace ba = bwa_aln;
    FMI_search fmi((ref().fa + ".meth").c_str()); fmi.load_index(false, 1, false, BWAMEM_HUGETLB_OFF, false);
    auto *bwt = ba::bwt_from_fmi(&fmi); ba::aln_meth_ref_t mr;
    REQUIRE(ba::aln_meth_ref_load(&mr, ref().fa.c_str(), fmi.idx->bns, 1) == 0);
    int c4 = -1;
    for (int i = 0; i < mr.bns->n_seqs; ++i) if (std::strcmp(mr.bns->anns[i].name, "c4") == 0) c4 = i;
    REQUIRE(c4 >= 0);
    auto *opt = ba::gap_init_opt(); opt->fnr = -1; opt->max_diff = 2; opt->seed_len = 1024;
    ba::aln_meth_opt_t mo; mo.m.meth_chem = METH_CHEM_TAPS; ba::aln_meth_opt_finish(&mo);
    const auto fq = write_fq("dup-candidate-paths.fq", std::vector<std::string>(2, ref().read_u));
    auto *input = ba::bwa_seq_open(fq.c_str()); int n = 0;
    auto *reads = ba::bwa_read_seq(input, 2, &n, BWA_MODE_COMPREAD, 0); REQUIRE(n == 2);
    // Both copies have a candidate window at read_u's locus; the second reaches it along another
    // path (an insertion and a deletion, the same span), which scores differently.
    const int len = (int)ref().read_u.size();
    const int64_t rb = mr.bns->anns[c4].offset + ref().u_pos - 1;
    const std::vector<uint32_t> ungapped{(uint32_t)len << 4};
    const std::vector<uint32_t> gapped{10u << 4, 1u << 4 | 1, 10u << 4, 1u << 4 | 2, (uint32_t)(len - 21) << 4};
    // Independently: every c4 window that read_u matches base for base, a read T over a reference C
    // (a TAPS conversion) included, is a tied ungapped locus with nothing to pay. The candidate's
    // (u_pos) is one; its gapped path costs more there and does not replace it.
    std::vector<std::pair<int64_t, int>> expected; // (locus, conversions)
    for (size_t p = 0; p + len <= ref().c4.size(); ++p) {
        int n_conv = 0, k = 0;
        for (; k < len; ++k) {
            const char a = ref().read_u[k], b = ref().c4[p + k];
            if (a == b) continue;
            if (a != 'T' || b != 'C') break;
            ++n_conv;
        }
        if (k == len) expected.emplace_back(mr.bns->anns[c4].offset + (int64_t)p, n_conv);
    }
    REQUIRE(expected.size() == 4u);
    CHECK(expected[0] == std::make_pair(rb, 0));
    for (bool same : {true, false}) {
        CAPTURE(same);
        std::vector<std::vector<ba::aln_meth_candidate_t>> extra(2);
        extra[0].push_back({c4, 0, 1, rb, rb + len, ungapped});
        extra[1].push_back({c4, 0, 1, rb, rb + len, same ? ungapped : gapped});
        std::vector<ba::aln_meth_read_t> results(n); ba::aln_dup_state_t dup;
        CHECK(ba::aln_meth_search_batch(bwt, &mr, &mo, n, reads, opt, 1, 1, 0, 1, 1000, &dup, results.data(),
                                        extra.data()) == (same ? 1 : 0));
        for (int i = 0; i < n; ++i) {
            CAPTURE(i);
            // The read searched alone, with no duplicate to share: what a shared or unshared
            // result must equal.
            ba::aln_meth_read_t single; ba::aln_dup_state_t single_dup;
            ba::aln_meth_search_batch(bwt, &mr, &mo, 1, reads + i, opt, 1, 1, i, 1, 1000, &single_dup, &single,
                                      &extra[i]);
            const ba::aln_meth_read_t &r = results[i];
            CHECK(r.flags == single.flags);
            CHECK(r.c1 == single.c1); CHECK(r.c2 == single.c2);
            REQUIRE(r.main == single.main);
            REQUIRE(r.hits.size() == single.hits.size());
            for (size_t k = 0; k < r.hits.size(); ++k) {
                const auto &x = r.hits[k], &y = single.hits[k];
                CHECK(std::make_tuple(x.rid, x.hyp, x.strand, x.pos, x.cost, x.span, x.alt) ==
                      std::make_tuple(y.rid, y.hyp, y.strand, y.pos, y.cost, y.span, y.alt));
                CHECK(std::make_tuple(x.n_mm, x.n_gapo, x.n_gape, x.n_conv) ==
                      std::make_tuple(y.n_mm, y.n_gapo, y.n_gape, y.n_conv));
                CHECK(cigar_of(r, x) == cigar_of(single, y));
            }
            // Independently: the expected loci, each ungapped with nothing to pay.
            CHECK(r.c1 == (int)expected.size()); CHECK(r.c2 == 0); CHECK(r.flags == 0);
            REQUIRE(r.hits.size() == expected.size());
            for (size_t k = 0; k < expected.size(); ++k) {
                CAPTURE(k);
                const auto &h = r.hits[k];
                CHECK(h.rid == c4); CHECK(h.strand == 0); CHECK(h.span == len);
                CHECK(h.pos == (bwtint_t)expected[k].first); CHECK(h.n_conv == expected[k].second);
                CHECK(h.cost == 0); CHECK(h.n_mm == 0); CHECK(h.n_gapo == 0); CHECK(h.n_gape == 0);
                CHECK(cigar_of(r, h) == ungapped);
            }
            CHECK((r.main >= 0 && (size_t)r.main < r.hits.size()));
        }
    }
    ba::bwa_free_read_seq(n, reads); ba::bwa_seq_close(input);
    free(opt); ba::aln_meth_ref_free(&mr); ba::bwt_destroy(bwt);
}

TEST_CASE("aln_meth: ks_meth_search applies incomplete window flags before the genomic search")
{
    namespace ba = bwa_aln;
    FMI_search fmi((ref().fa + ".meth").c_str()); fmi.load_index(false, 1, false, BWAMEM_HUGETLB_OFF, false);
    auto *bwt = ba::bwt_from_fmi(&fmi); ba::aln_meth_ref_t mr;
    REQUIRE(ba::aln_meth_ref_load(&mr, ref().fa.c_str(), fmi.idx->bns, 1) == 0);
    int c4 = -1;
    for (int i = 0; i < mr.bns->n_seqs; ++i) if (std::strcmp(mr.bns->anns[i].name, "c4") == 0) c4 = i;
    REQUIRE(c4 >= 0);
    // Eight sites just past read_u's exact locus, none under the read: every site's ALT window
    // holds the read (and its first one-conversion copy), so the window search has more rows to
    // look up than the genomic search.
    ba::ks_sites_t sites; sites.by_contig.resize(mr.bns->n_seqs);
    const int read_end = ref().u_pos - 1 + (int)ref().read_u.size();
    for (int k = 0; k < 8; ++k) {
        const int p = read_end + 10 + 5 * k;
        const uint8_t code = (uint8_t)std::string("ACGT").find(ref().c4[p]);
        REQUIRE(code < 4);
        sites.by_contig[c4].push_back({p, code, (uint8_t)(code ^ 1)}); // a transversion
        ++sites.n_sites;
    }
    mr.known_snvs = &sites;
    ba::ks_windows_t windows;
    {
        bwa_tests::StderrCapture quiet;
        REQUIRE(ba::ks_windows_build(&sites, mr.bns, mr.pac, 1, 1, &windows) == 0);
    }
    ba::aln_opt_t o; ba::aln_opt_init(&o);
    o.gap.fnr = -1; o.gap.max_diff = 2; o.gap.seed_len = 1024; o.gap.n_threads = 1;
    ba::aln_meth_opt_t mo; mo.m.meth_chem = METH_CHEM_TAPS; mo.m.meth_scoring = MEM_METH_SCORING_NEUTRAL;
    ba::aln_meth_opt_finish(&mo);
    o.meth = &mo; o.meth_ref = &mr; o.known_snvs = &sites; o.ks_windows = &windows;
    const auto fq = write_fq("ks-incoming-flags.fq", {ref().read_u});
    auto *input = ba::bwa_seq_open(fq.c_str()); int n = 0;
    auto *reads = ba::bwa_read_seq(input, 1, &n, BWA_MODE_COMPREAD, 0); REQUIRE(n == 1);
    const ba::gap_opt_t discovery = ba::ks_discovery_opt(o.gap);
    auto window_search = [&]() {
        ba::aln_meth_read_t w;
        ba::aln_meth_search_batch(windows.bwt, windows.mref, &mo, 1, reads, &discovery, 1, 1, 0, 1, o.max_xa_hits,
                                  nullptr, &w, nullptr, nullptr, true);
        return w;
    };
    // The genomic search given the window search's candidates and certificate but not its flags:
    // what the read would get if the flags were only OR-ed into the result afterwards.
    auto genomic_without_flags = [&](const ba::aln_meth_read_t &w) {
        std::vector<ba::aln_meth_candidate_t> extra; ba::ks_meth_candidates(&windows, mr.bns, w, &extra);
        const ba::aln_meth_alt_bound_t certificate = ba::aln_meth_alt_certificate(w);
        ba::aln_meth_read_t r;
        ba::aln_meth_search_batch(bwt, &mr, &mo, 1, reads, &o.gap, 1, 1, 0, 0, o.max_xa_hits, nullptr, &r, &extra,
                                  &certificate);
        return r;
    };
    // Calibrate: with the budget set to the rows the genomic search needs without stopping, the
    // window search (more rows) is bound while the genomic search alone completes.
    const auto full_window = window_search();
    REQUIRE(full_window.flags == 0);
    ba::aln_meth_read_t unstopped = full_window; unstopped.site_bound.valid = false; unstopped.full_width_bound.value = 0;
    const auto full_genomic = genomic_without_flags(unstopped);
    REQUIRE(full_genomic.flags == 0);
    REQUIRE(full_window.n_lookups > full_genomic.n_lookups);
    mo.p_best = (int)full_genomic.n_lookups; mo.p_second = 0;
    const auto bound_window = window_search();
    REQUIRE(bound_window.flags == ba::ALN_METH_BIND_BEST);
    // Flags OR-ed in after the search would leave neutral scoring's conversion tie-break in place.
    const auto or_after = genomic_without_flags(bound_window);
    REQUIRE(or_after.flags == 0); REQUIRE(or_after.n_conv_moved == 3);
    // ks_meth_search hands them to the genomic slot first: the read is incomplete from the start,
    // so the tie-break (which assumes every equal-cost locus was seen) is skipped.
    ba::aln_meth_read_t res;
    ba::ks_meth_search(&o, bwt, 1, reads, 1, 0, 0, nullptr, nullptr, &res);
    CHECK(res.flags == ba::ALN_METH_BIND_BEST);
    CHECK(res.n_conv_moved == 0); CHECK(res.conv_gap == 0);
    ba::bwa_free_read_seq(n, reads); ba::bwa_seq_close(input);
    mr.known_snvs = nullptr; ba::ks_windows_free(&windows); ba::aln_meth_ref_free(&mr); ba::bwt_destroy(bwt);
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
            ba::aln_meth_occ_t o{};o.rid=0;o.hyp=1;o.strand=0;o.pos=100+100*k;
            o.cost=k==3?12:6;o.n_conv=k==0?2:k==1?2+gap:1;
            o.n_mm=k==0?1:0;o.n_gapo=k==1?1:0;
            o.alignment.rb=o.pos;o.alignment.re=o.pos+40;
            o.alignment.cigar={40u<<4};o.alignment.raw_cigar=o.alignment.cigar;
            result.occ.push_back(o);
        }
        result.flags=incomplete?ba::ALN_METH_BIND_M:0;
        ba::aln_meth_slot_t slot;slot.grp=&group;slot.res=&result;slot.keep=1;
        ba::aln_meth_slot_finish(&slot);
        CHECK(result.hits[result.main].cost==6);
        if(incomplete) { CHECK(result.n_conv_moved==0);CHECK(result.conv_gap==0); }
        else {
            // The third cost-6 locus has fewer conversions than both other compositions.
            CHECK(result.hits[result.main].pos==300);
            CHECK(result.c1==1); CHECK(result.c2==3);
            CHECK(result.n_conv_moved==2); CHECK(result.conv_gap==1);
            for(int ordinal=0;ordinal<40;++ordinal) {
                ba::aln_meth_redraw_hits(&result,ordinal);
                CHECK(result.hits[result.main].pos==300);
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
            ba::aln_meth_occ_t o{};o.rid=0;o.hyp=1;o.strand=j;o.pos=j?500:350;
            o.alignment.rb=o.pos;o.alignment.re=o.pos+40;o.alignment.cigar={40u<<4};
            result[j].occ.push_back(o);result[j].main=0;result[j].c1=1;result[j].c2=1;
            result[j].n_conv_moved=1;result[j].conv_gap=j?2:1;
            if(j==0 && moved) {
                auto &moved_from=result[j].occ[0];
                moved_from.pos=20;moved_from.alignment.rb=20;moved_from.alignment.re=20+40;
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
        store_pair(r);ba::aln_meth_pe_stage1(opt,1,seqs,r);
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
            ba::aln_meth_occ_t o{};o.rid=3;o.hyp=1;o.strand=0;
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
            // The span-41 path: a deletion after query base 20 (strand 0: the complemented read's
            // index 19 is the base the deletion follows).
            ba::gap_entry_t entry={};entry.n_del=delta;
            if(delta) {entry.n_gapo=1;entry.last_diff_pos=(int)ALN_METH_GAP(19,1,1);}
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
        for(const auto &o:result[0].hits) keys.emplace(o.rid,o.hyp,o.strand,o.pos);
        CHECK(keys.size()==result[0].hits.size());
        const bool complete=budget==4596;
        if(complete) {CHECK(result[0].flags==0);CHECK(result[0].hits.size()==601);CHECK(result[0].c1==601);}
        else {CHECK(result[0].flags!=0);CHECK(result[0].hits.size()<=500);}
        // Force both single-end choices away from the sole compatible pair.
        for(int j=0;j<2;++j) for(size_t k=0;k<result[j].hits.size();++k) {
            const auto &o=result[j].hits[k];
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
            CHECK(result[0].hits[result[0].main].pos==seqs[0]->pos);
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

TEST_CASE("aln_meth: CIGAR-capped mates retain their placement even when the lookup cap also binds")
{
    namespace ba = bwa_aln;
    FMI_search fmi((ref().fa + ".meth").c_str());
    fmi.load_index(false, 1);
    ba::aln_meth_ref_t mr;
    REQUIRE(ba::aln_meth_ref_load(&mr, ref().fa.c_str(), fmi.idx->bns, 1) == 0);
    auto *opt = ba::gap_init_opt();
    auto *pe = ba::bwa_init_pe_opt();
    ba::aln_meth_opt_t mo;
    ba::aln_meth_opt_finish(&mo);
    ba::bwase_initialize();
    for (int end : {0, 1}) for (bool retained : {false, true})
        for (int flags : std::vector<int>{ba::ALN_METH_BIND_CIGAR, ba::ALN_METH_BIND_CIGAR | ba::ALN_METH_BIND_BEST,
                                         ba::ALN_METH_BIND_BEST}) {
            CAPTURE(end); CAPTURE(retained); CAPTURE(flags);
            ba::bwa_seq_t *seqs[2];
            ba::bwa_seqio_t *inputs[2];
            for (int j = 0; j < 2; ++j) {
                const std::string read = j ? aln_test::revcomp(ref().c1.substr(700, 50)) : ref().c1.substr(550, 50);
                inputs[j] = ba::bwa_seq_open(write_fq(j ? "cigar-pe-2.fq" : "cigar-pe-1.fq", {read}).c_str());
                REQUIRE(inputs[j] != nullptr);
                int n = 0;
                seqs[j] = ba::bwa_read_seq(inputs[j], 1, &n, BWA_MODE_COMPREAD, 0);
                REQUIRE(n == 1);
            }
            ba::aln_meth_read_t result[2];
            ba::aln_meth_read_t *r[2] = {result, result + 1};
            for (int j = 0; j < 2; ++j) {
                if (j == end && !retained) continue;
                ba::aln_meth_occ_t o{};
                o.rid = 0; o.hyp = 1; o.strand = j;
                o.pos = j == end ? 900 : (j ? 700 : 550);
                o.alignment.rb = o.pos; o.alignment.re = o.pos + 50;
                o.alignment.cigar = {50u << 4};
                result[j].occ.push_back(o);
                result[j].main = 0; result[j].c1 = 1;
            }
            result[end].flags = flags;
            store_pair(r);
            ba::aln_meth_pe_stage1(opt, 1, seqs, r);
            ba::isize_info_t ii = {};
            ii.avg = 200; ii.std = 20; ii.ap_prior = 1e-5;
            uint64_t total[2] = {}, mapped[2] = {};
            ba::aln_meth_pe_stage2(&mr, &mo, opt, pe, &ii, 1, seqs, r, total, mapped);
            if (flags & ba::ALN_METH_BIND_CIGAR) {
                CHECK((int)seqs[end]->type == (retained ? BWA_TYPE_REPEAT : BWA_TYPE_NO_MATCH));
                CHECK(result[end].main == (retained ? 0 : -1));
                if (retained) CHECK(seqs[end]->pos == 900);
            } else {
                CHECK((int)seqs[end]->type == BWA_TYPE_MATESW);
                CHECK(seqs[end]->pos == (end ? 700 : 550));
            }
            CHECK(result[end].flags == flags);
            CHECK((int)seqs[end]->mapQ == 0);
            CHECK((int)seqs[end]->c1 == 0);
            CHECK((int)seqs[end]->c2 == 0);
            CHECK(seqs[end]->n_multi == 0);
            CHECK((seqs[end]->extra_flag & SAM_FPP) == 0);
            kstring_t text = {0, 0, nullptr};
            ba::aln_meth_print_sam1(&text, &mr, &mo, seqs[end], r[end], end == 0,
                                   seqs[1 - end], BWA_MODE_COMPREAD, 30);
            const std::string record(text.s, text.l);
            free(text.s);
            CHECK(record.find("\tX0:i:") == std::string::npos);
            CHECK(record.find("\tX1:i:") == std::string::npos);
            CHECK(record.find("\tXA:Z:") == std::string::npos);
            for (int j = 0; j < 2; ++j) {
                ba::bwa_free_read_seq(1, seqs[j]);
                CHECK(ba::bwa_seq_close(inputs[j]) == 0);
            }
        }
    free(pe); free(opt); ba::aln_meth_ref_free(&mr);
}

TEST_CASE("aln_meth: cap rescue skips nonfinite insert penalties without losing safe placements")
{
    namespace ba = bwa_aln;
    FMI_search fmi((ref().fa + ".meth").c_str()); fmi.load_index(false, 1);
    ba::aln_meth_ref_t mr;
    REQUIRE(ba::aln_meth_ref_load(&mr, ref().fa.c_str(), fmi.idx->bns, 1) == 0);
    auto *opt = ba::gap_init_opt();
    auto *pe = ba::bwa_init_pe_opt();
    ba::aln_meth_opt_t mo; ba::aln_meth_opt_finish(&mo); ba::bwase_initialize();
    for (bool extreme : {false, true}) for (bool safe_anchor : {false, true}) {
        CAPTURE(extreme); CAPTURE(safe_anchor);
        ba::bwa_seq_t *seqs[2]; ba::bwa_seqio_t *inputs[2];
        for (int j = 0; j < 2; ++j) {
            const auto read = j ? aln_test::revcomp(ref().c1.substr(700, 50)) : ref().c1.substr(575, 50);
            inputs[j] = ba::bwa_seq_open(write_fq(j ? "penalty2.fq" : "penalty1.fq", {read}).c_str());
            REQUIRE(inputs[j] != nullptr);
            int n = 0; seqs[j] = ba::bwa_read_seq(inputs[j], 1, &n, BWA_MODE_COMPREAD, 0);
            REQUIRE(n == 1);
        }
        ba::aln_meth_read_t result[2]; ba::aln_meth_read_t *r[2] = {result, result + 1};
        for (int pos : safe_anchor ? std::vector<int>{575, 550} : std::vector<int>{575}) {
            ba::aln_meth_occ_t a{}; a.rid = 0; a.hyp = 1; a.pos = pos;
            a.alignment.rb = pos; a.alignment.re = pos + 50; a.alignment.cigar = {50u << 4};
            result[0].occ.push_back(a);
        }
        result[0].main = 0; result[0].c1 = result[0].occ.size();
        ba::aln_meth_occ_t retained{}; retained.rid = 0; retained.hyp = 1; retained.strand = 1;
        retained.pos = 900; retained.alignment.rb = 900; retained.alignment.re = 950;
        retained.alignment.cigar = {50u << 4};
        result[1].occ.push_back(retained); result[1].main = 0; result[1].flags = ba::ALN_METH_BIND_BEST;
        store_pair(r);
        ba::aln_meth_pe_stage1(opt, 1, seqs, r);
        ba::isize_info_t ii{};
        ii.avg = extreme ? 200 : 199.5;
        ii.std = extreme ? std::numeric_limits<double>::min() : std::sqrt(.2);
        ii.high = 202; ii.ap_prior = 1e-5;
        // The 175-base candidate underflows the Gaussian tail; the 200-base pair is safe.
        uint64_t total[2] = {}, mapped[2] = {};
        ba::aln_meth_pe_stage2(&mr, &mo, opt, pe, &ii, 1, seqs, r, total, mapped);
        REQUIRE(result[1].hits.size() == 1);
        CHECK(seqs[1]->pos == (safe_anchor ? 700 : 900));
        CHECK(result[1].hits[0].pos == seqs[1]->pos);
        CHECK(cigar_of(result[1], result[1].hits[0]) == std::vector<uint32_t>{50u << 4});
        CHECK(result[1].hits[0].cost == 0);
        REQUIRE(seqs[1]->n_cigar == 1);
        CHECK(__cigar_op(seqs[1]->cigar[0]) == 0); CHECK(__cigar_len(seqs[1]->cigar[0]) == 50);
        CHECK((int)seqs[1]->type == (safe_anchor ? BWA_TYPE_MATESW : BWA_TYPE_REPEAT));
        CHECK(result[0].main == (safe_anchor ? 1 : 0));
        CHECK(result[1].flags == ba::ALN_METH_BIND_BEST);
        CHECK((int)seqs[1]->mapQ == 0); CHECK((int)seqs[1]->seQ == 0);
        CHECK((int)seqs[1]->c1 == 0); CHECK((int)seqs[1]->c2 == 0);
        CHECK(seqs[1]->n_multi == 0); CHECK((seqs[1]->extra_flag & SAM_FPP) == 0);
        for (int j = 0; j < 2; ++j) {
            ba::bwa_free_read_seq(1, seqs[j]); CHECK(ba::bwa_seq_close(inputs[j]) == 0);
        }
    }
    free(pe); free(opt); ba::aln_meth_ref_free(&mr);
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
            ba::aln_meth_occ_t a{};a.rid=0;a.hyp=1;a.strand=0;a.pos=550+k;
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
        store_pair(r);ba::aln_meth_pe_stage1(opt,1,seqs,r);
        ba::isize_info_t ii={};ii.avg=200;ii.std=20;ii.high=0;ii.ap_prior=1e-5;
        if(bad_range && anchors<=255) {
            if(bad_range==1) ii.std=std::numeric_limits<double>::infinity();
            if(bad_range==2) ii.std=1000000;
            if(bad_range==3) result[0].hits[0].span=(int32_t)(mr.bns->anns[0].len+1-(int64_t)result[0].hits[0].pos);
            // The fatal path must exit 1; returning ends the child with 0 and fails the check.
            bwa_tests::TempDir err_dir{"aln_meth_cap_fatal"};
            const std::string err=err_dir.path()+"/stderr";
            const int status=bwa_tests::run_in_child([&] {
                // A window clipped to its contig exceeds 65535 bases only on a contig longer than
                // that: lengthen the anchors' contig (c1) in this child alone. The window size
                // check exits before the window is fetched.
                if(bad_range==2) mr.bns->anns[0].len=1<<20;
                uint64_t total[2]={},mapped[2]={};
                ba::aln_meth_pe_stage2(&mr,&mo,opt,pe,&ii,1,seqs,r,total,mapped);
            },err.c_str());
            CHECK(WIFEXITED(status));CHECK(WEXITSTATUS(status)==1);
            const char *expected[]={"","ERROR: --meth invalid cap-rescue insert-size range",
                                    "ERROR: --meth rescue window exceeds allocation budget",
                                    "ERROR: --meth cap-rescue anchor leaves contig"};
            std::ifstream in(err);std::stringstream text;text<<in.rdbuf();
            CHECK(text.str()==std::string(expected[bad_range])+"\n");
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
            REQUIRE(result[1].hits.size()==1);CHECK(result[1].hits[0].cost==0);
            CHECK(result[1].hits[0].n_conv==0);CHECK(result[1].hits[0].n_mm==0);
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
            ba::aln_meth_occ_t a{};a.rid=rid;a.hyp=1;a.strand=0;a.pos=base+100+400*(k^(int)order);
            a.alignment.rb=a.pos;a.alignment.re=a.pos+50;a.alignment.cigar={50u<<4};result[0].occ.push_back(a);
        }
        result[0].main=0;result[0].c1=2;result[1].flags=ba::ALN_METH_BIND_BEST;
        // Existing best-known cost zero still accepts neutral conversions at zero cost.
        ba::aln_meth_occ_t known{};known.rid=rid;known.hyp=1;known.strand=1;known.pos=base;
        known.alignment.rb=base;known.alignment.re=base+50;known.alignment.cigar={50u<<4};
        result[1].occ.push_back(known);result[1].main=0;
        store_pair(r);ba::aln_meth_pe_stage1(opt,1,seqs,r);
        ba::isize_info_t ii={};ii.avg=200;ii.std=5;ii.ap_prior=1e-5;
        uint64_t total[2]={},mapped[2]={};ba::aln_meth_pe_stage2(&mr,&mo,opt,pe,&ii,1,seqs,r,total,mapped);
        REQUIRE(result[1].hits.size()==1);const auto &o=result[1].hits[0];
        CHECK(o.pos==(bwtint_t)(base+(neutral?650:250)));
        CHECK(o.n_conv==(neutral?0:1));CHECK(o.n_mm==0);CHECK(o.cost==0);
        CHECK(cigar_of(result[1],o)==std::vector<uint32_t>{50u<<4});
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
    // Beside the shared MethRef index (ref().fa), under a name of its own.
    const std::string fa = ref().tmp.unique_file("_homopolymer.fa");
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
        ba::aln_meth_occ_t o{};o.rid=0;o.strand=0;o.hyp=1;
        o.pos=mr.bns->anns[0].offset+anchor;o.alignment.rb=o.pos;o.alignment.re=o.pos+50;
        o.alignment.cigar={50u<<4};result[0].occ.push_back(o);result[0].main=0;result[0].c1=1;
        store_pair(r);ba::aln_meth_pe_stage1(opt,1,seqs,r);
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

TEST_CASE("aln_meth: a read with a '-' base places on both strands instead of stopping the run")
{
    // aln reads '-' as base code 5 (bwa's own table), which the search accepts as an ambiguous
    // base. The refinement must take it as N on both strands, not reject the window and stop.
    const std::string x = ref().c1.substr(800, 50);
    std::string fwd = conv(x, 'C', 'T');
    std::string rev = conv(aln_test::revcomp(x), 'C', 'T');
    fwd[25] = '-';
    rev[25] = '-';
    const auto v = run_aln({"--meth"}, {write_fq("dash.fq", {fwd, rev})});
    REQUIRE(v.size() == 2);
    CHECK(v[0].flag == 0);
    CHECK(v[0].rname == "c1");
    CHECK(v[0].pos == 801);
    CHECK(v[0].cigar == "50M");
    CHECK(v[1].flag == 16);
    CHECK(v[1].rname == "c1");
    CHECK(v[1].pos == 801);
    CHECK(v[1].cigar == "50M");
    // SEQ prints the '-' as N (stock aln prints a NUL byte), on both strands.
    const std::string sam = slurp_out();
    CHECK(sam.find('\0') == std::string::npos);
    std::vector<std::string> seqs;
    std::stringstream ss(sam);
    for (std::string l; std::getline(ss, l);) {
        if (l.empty() || l[0] == '@') continue;
        std::stringstream fs(l);
        std::string f;
        for (int k = 0; k < 10 && std::getline(fs, f, '\t'); ++k) {}
        seqs.push_back(f);
    }
    REQUIRE(seqs.size() == 2);
    std::string want_fwd = fwd, want_rev = aln_test::revcomp(rev.substr(0, 25) + "N" + rev.substr(26));
    want_fwd[25] = 'N';
    CHECK(seqs[0] == want_fwd);
    CHECK(seqs[1] == want_rev);
}

TEST_CASE("aln_meth: a window longer than the realigner accepts is dropped, not fatal")
{
    // A 65535-base read (the longest aln accepts) with one deleted base spans 65536 reference
    // bases, one more than the realigner evaluates. That candidate cannot be scored: the read is
    // reported without it (unmapped here, as it is the only candidate) and the next read still maps.
    bwa_tests::TempDir tmp{"aln_meth_long"}; // this case's index and reads, removed when it ends
    const std::string fa = tmp.unique_file();
    uint32_t st = 777;
    const std::string contig = aln_test::lcg_dna(st, 70000);
    REQUIRE(aln_test::write_test_fasta(fa, {{"long", contig}}));
    std::vector<std::string> args = {"index", "--meth", fa};
    std::vector<char *> argv;
    for (auto &s : args) argv.push_back(&s[0]);
    argv.push_back(nullptr);
    optind = 1;
#if defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || defined(__NetBSD__)
    optreset = 1;
#endif
    REQUIRE(bwa_index((int)args.size(), argv.data()) == 0);
    std::string longest = contig.substr(1000, 65536);
    longest.erase(30000, 1); // one deletion in the read: 65535 bases over 65536
    REQUIRE(longest.size() == 65535);
    const std::string next = contig.substr(200, 50);
    const std::string out = tmp.path() + "/long.sam", fq = tmp.path() + "/long.fq";
    {
        std::ofstream reads(fq);
        const std::string rs[2] = {conv(longest, 'C', 'T'), conv(next, 'C', 'T')};
        for (int i = 0; i < 2; ++i)
            reads << "@r" << i << "\n" << rs[i] << "\n+\n" << std::string(rs[i].size(), 'I') << "\n";
        REQUIRE(reads.good());
    }
    std::vector<std::string> a = {"aln", "-f", out, "--meth", "-n", "1", fa, fq};
    std::vector<char *> v;
    for (auto &s : a) v.push_back(&s[0]);
    v.push_back(nullptr);
    REQUIRE(main_aln((int)a.size(), v.data()) == 0);
    std::ifstream in(out);
    const auto recs = parse_records(in);
    REQUIRE(recs.size() == 2);
    CHECK((recs[0].flag & 4) != 0);
    CHECK(recs[1].flag == 0);
    CHECK(recs[1].rname == "long");
    CHECK(recs[1].pos == 201);
    CHECK(recs[1].cigar == "50M");
}

TEST_CASE("aln_meth: mate rescue takes a read '-' base as N instead of stopping the run")
{
    std::vector<std::string> one,two;
    for(int i=0;i<45;++i) {
        const int start=100+50*i,span=180+i%21;
        one.push_back(ref().c1.substr(start,50));
        two.push_back(aln_test::revcomp(ref().c1.substr(start+span-50,50)));
    }
    // As in the rescue test above, search cannot propose the last mate; a '-' (code 5) in it
    // must reach the local rescue and the realigner as N.
    for(int at : {15,25,35}) two.back()[at]=two.back()[at]=='A'?'C':'A';
    two.back()[5]='-';
    const auto v=run_aln({"--meth","-n","0","-l","1024"},
                        {write_fq("rescue-dash-one.fq",one),write_fq("rescue-dash-two.fq",two)});
    REQUIRE(v.size()==90);
    CHECK((v[89].flag&4)==0);CHECK((v[89].flag&2)!=0);
    CHECK(v[89].pos==2433);CHECK(v[89].cigar=="50M");
    CHECK(tag_value(v[89],"XT:A")=="M");
    const std::string sam=slurp_out();
    CHECK(sam.find('\0')==std::string::npos);
}

namespace {
// Read 2 of a planted pair: the reverse complement of the fragment's last reference bases, with
// `shift` bases deleted from the read (shift > 0, so it spans 50 + shift reference bases) or
// -shift bases inserted (it spans 50 + shift), always 50 bases long. The fragment, from read 1's
// start to read 2's reference end, is `frag` bases.
std::string gapped_read2(const std::string &genome, int start, int frag, int shift)
{
    const int span = 50 + shift;
    std::string g = genome.substr(start + frag - span, span);
    if (shift > 0) g.erase(25, shift);
    else g.insert(23, std::string("GATTACA").substr(0, -shift));
    REQUIRE(g.size() == 50);
    return aln_test::revcomp(g);
}

// The mean the "[infer_isize] inferred external isize" report prints (3 decimals), or -1.
double reported_isize(const std::string &err)
{
    const std::string key = "inferred external isize from ";
    const size_t at = err.find(key);
    if (at == std::string::npos) return -1;
    const size_t colon = err.find(": ", at);
    return colon == std::string::npos ? -1 : std::strtod(err.c_str() + colon + 2, nullptr);
}
} // namespace

TEST_CASE("aln_meth: ref-span-insert holds a gapped mate's reference span, not its read length, to -a"
          * doctest::test_suite("unit/aln_meth"))
{
    // With -A -a 200, a pair is proper when its insert is at most 200. Read 2 of pair 0 has a
    // 4-base deletion: its insert is 202 by reference span, 198 by read length. Read 2 of pair 1
    // has a 4-base insertion: 198 by span, 202 by read length. TLEN is from the CIGAR either way.
    bool meth = false, converted = false;
    SUBCASE("ordinary") {}
    SUBCASE("meth exact") { meth = true; }
    SUBCASE("meth converted") { meth = converted = true; }
    const int start[2] = {500, 1500}, frag[2] = {202, 198}, shift[2] = {4, -4};
    std::vector<std::string> reads[2];
    for (int i = 0; i < 2; ++i) {
        std::string one = ref().c1.substr(start[i], 50), two = gapped_read2(ref().c1, start[i], frag[i], shift[i]);
        if (converted) { one = conv(one, 'C', 'T'); two = conv(two, 'G', 'A'); }
        reads[0].push_back(one); reads[1].push_back(two);
    }
    const std::vector<std::string> input = {write_fq("span-a1.fq", reads[0]), write_fq("span-a2.fq", reads[1])};
    std::vector<std::string> args = {"-A", "-a", "200", "-n", "3", "-e", "10", "-l", "1024"};
    if (meth) args.push_back("--meth");
    for (bool option : {false, true}) {
        CAPTURE(option);
        std::vector<std::string> a = args;
        if (option) a.push_back("--ref-span-insert");
        const auto out = run_aln(a, input);
        REQUIRE(out.size() == 4);
        for (int i = 0; i < 2; ++i) {
            CAPTURE(i);
            const Rec &r1 = out[2 * i], &r2 = out[2 * i + 1];
            const bool proper = option ? frag[i] <= 200 : frag[i] - shift[i] <= 200;
            CHECK(r1.flag == (proper ? 99 : 97)); CHECK(r2.flag == (proper ? 147 : 145));
            CHECK(r1.rname == "c1"); CHECK(r2.rname == "c1");
            CHECK(r1.pos == start[i] + 1); CHECK(r1.cigar == "50M");
            CHECK(r2.pos == start[i] + frag[i] - (50 + shift[i]) + 1);
            // The whole CIGAR is one indel between two matches, at a placement that turns the planted
            // read back into the reference: the aligner may pick any of several equally scored ones,
            // and under --meth a C/T difference on this strand is a free conversion.
            auto fold = [&](std::string x) { return meth ? conv(x, 'C', 'T') : x; };
            const std::string window = fold(ref().c1.substr(start[i] + frag[i] - (50 + shift[i]), 50 + shift[i]));
            const std::string planted = fold(aln_test::revcomp(gapped_read2(ref().c1, start[i], frag[i], shift[i])));
            const char op = shift[i] > 0 ? 'D' : 'I';
            int left = -1, right = -1;
            CHECK(std::sscanf(r2.cigar.c_str(), shift[i] > 0 ? "%dM4D%dM" : "%dM4I%dM", &left, &right) == 2);
            CHECK(r2.cigar == std::to_string(left) + "M4" + op + std::to_string(right) + "M");
            REQUIRE(left > 0); REQUIRE(right > 0); REQUIRE(left + right == 50 - (shift[i] > 0 ? 0 : 4));
            if (shift[i] > 0) CHECK(window.substr(0, left) + window.substr(left + 4) == planted);
            else CHECK(planted.substr(0, left) + planted.substr(left + 4) == window);
            CHECK(r1.tlen == frag[i]); CHECK(r2.tlen == -frag[i]);
        }
    }
}

TEST_CASE("aln_meth: ref-span-insert estimates the insert size from gapped mates' reference spans"
          * doctest::test_suite("unit/aln_meth"))
{
    // 30 planted pairs whose read 2 has a 4-base deletion: each fragment (read 1's start to read
    // 2's reference end) is 4 bases longer than the read-length measurement of it.
    bool meth = false, converted = false;
    SUBCASE("ordinary") {}
    SUBCASE("meth exact") { meth = true; }
    SUBCASE("meth converted") { meth = converted = true; }
    std::vector<std::string> reads[2];
    double sum = 0;
    for (int i = 0; i < 30; ++i) {
        const int start = 100 + 80 * i, frag = 180 + i % 21;
        std::string one = ref().c1.substr(start, 50), two = gapped_read2(ref().c1, start, frag, 4);
        if (converted) { one = conv(one, 'C', 'T'); two = conv(two, 'G', 'A'); }
        reads[0].push_back(one); reads[1].push_back(two);
        sum += frag;
    }
    const std::vector<std::string> input = {write_fq("span-ii1.fq", reads[0]), write_fq("span-ii2.fq", reads[1])};
    std::vector<std::string> args = {"-n", "2", "-e", "10", "-l", "1024"};
    if (meth) args.push_back("--meth");
    const int verbose = bwa_verbose;
    for (bool option : {false, true}) {
        CAPTURE(option);
        std::vector<std::string> a = args;
        if (option) a.push_back("--ref-span-insert");
        std::string err;
        {
            bwa_tests::StderrCapture capture;
            const auto out = run_aln(a, input);
            err = capture.text();
            REQUIRE(out.size() == 60);
            for (const auto &r : out) CHECK((r.flag & 2) != 0);
        }
        bwa_verbose = verbose;
        CHECK(reported_isize(err) == doctest::Approx(sum / 30 - (option ? 0 : 4)).epsilon(1e-5));
    }
}

TEST_CASE("aln_meth: exact-rescue-window rescues a mate at the stated window's end; plain aln keeps sampe's window"
          * doctest::test_suite("unit/aln_meth"))
{
    // 40 calibration pairs of 30 bp reads (inserts 96-104), then a singleton whose reverse mate,
    // which carries a transversion so -n 0 cannot place it, ends at the last base of the stated
    // window (avg + 3 std + 15 past the anchor's start, rounded out). sampe's window ends one or
    // two bases short, so stock rescue clips the mate there; the exact window holds all of it.
    namespace ba = bwa_aln;
    bool meth = false;
    SUBCASE("ordinary") {}
    SUBCASE("meth") { meth = true; }
    const int len = 30, n_cal = 40, anchor = 2300;
    std::vector<std::string> reads[2];
    ba::bwa_seq_t ends[2][n_cal] = {};
    for (int i = 0; i < n_cal; ++i) {
        const int start = 100 + 50 * i, insert = 100 + i % 9 - 4;
        reads[0].push_back(ref().c1.substr(start, len));
        reads[1].push_back(aln_test::revcomp(ref().c1.substr(start + insert - len, len)));
        ends[0][i].pos = start; ends[1][i].pos = start + insert - len; ends[1][i].strand = 1;
        ends[0][i].len = ends[1][i].len = len; ends[0][i].mapQ = ends[1][i].mapQ = 25;
    }
    ba::bwa_seq_t *seqs[2] = {ends[0], ends[1]};
    ba::isize_info_t ii = {};
    REQUIRE(ba::infer_isize(n_cal, seqs, &ii, 1e-5, 100000, nullptr, 0) == 0);
    int64_t sb, se, eb, ee;
    ba::aln_rescue_window(anchor, true, len, &ii, false, &sb, &se);
    ba::aln_rescue_window(anchor, true, len, &ii, true, &eb, &ee);
    // Independently of aln_rescue_window: the calibration inserts' mean, the stated window's end
    // (avg + 3 std + len/2 past the anchor, rounded up) and sampe's truncated end.
    double mean = 0;
    for (int i = 0; i < n_cal; ++i) mean += 100 + i % 9 - 4;
    mean /= n_cal;
    CHECK(ii.avg == doctest::Approx(mean));
    CHECK(ee == anchor + (int64_t)std::ceil(ii.avg + 3 * ii.std + len * 0.5));
    CHECK(se == (int64_t)((double)(int64_t)(anchor + ii.avg - 3 * ii.std - len * 1.5) + 6 * ii.std + 2 * len));
    const int lost = (int)(ee - se);
    REQUIRE(lost >= 1); REQUIRE(lost <= 2); REQUIRE(sb == eb);
    std::string mate = ref().c1.substr(ee - len, len);
    const char tv[] = {'C', 'A', 'T', 'G'}; // A>C, C>A, G>T, T>G: never a conversion
    mate[10] = tv[std::string("ACGT").find(mate[10])];
    reads[0].push_back(ref().c1.substr(anchor, len));
    reads[1].push_back(aln_test::revcomp(mate));
    const std::vector<std::string> input = {write_fq("exact-window1.fq", reads[0]), write_fq("exact-window2.fq", reads[1])};
    auto run = [&](std::vector<std::string> extra) {
        std::vector<std::string> args = {"-n", "0", "-l", "1024"};
        if (meth) args.push_back("--meth");
        args.insert(args.end(), extra.begin(), extra.end());
        const auto recs = run_aln(args, input);
        REQUIRE(recs.size() == 2 * (n_cal + 1));
        std::stringstream in(slurp_out()); // the records' text, without the header (@PG has the flags)
        std::string body;
        for (std::string l; std::getline(in, l);) if (l.empty() || l[0] != '@') body += l + "\n";
        return std::make_pair(recs, body);
    };
    const auto by_default = run({}), off = run({"--no-exact-rescue-window"}), on = run({"--exact-rescue-window"});
    const auto &stock = off.first[2 * n_cal + 1], &exact = on.first[2 * n_cal + 1];
    CHECK((stock.flag & 4) == 0); CHECK((exact.flag & 4) == 0);
    CHECK(stock.pos == ee - len + 1); CHECK(exact.pos == ee - len + 1);
    CHECK(stock.cigar == std::to_string(len - lost) + "M" + std::to_string(lost) + "S");
    CHECK(exact.cigar == std::to_string(len) + "M");
    CHECK(tag_value(exact, "NM:i") == "1");
    // Plain aln keeps sampe's window unless asked; --meth takes the exact one unless told not to.
    CHECK(by_default.second == (meth ? on.second : off.second));
    // The calibration pairs are proper and need no rescue: every field of their records, tags
    // included, is the same in both modes.
    std::stringstream off_body(off.second), on_body(on.second);
    for (int i = 0; i < 2 * n_cal; ++i) {
        CAPTURE(i);
        std::string off_line, on_line;
        REQUIRE(static_cast<bool>(std::getline(off_body, off_line)));
        REQUIRE(static_cast<bool>(std::getline(on_body, on_line)));
        CHECK(off_line == on_line);
    }
}

TEST_CASE("aln_meth: cap rescue takes a reverse anchor's end as mate rescue does")
{
    // A reverse anchor with a one-base deletion: 50 read bases over 51 reference bases. Mate rescue
    // places a forward mate's window from the anchor's read-length end (pos + len), or its reference
    // end under --ref-span-insert, and the cap-rescue window, exact or not, must start at the same
    // base. The window's offset is an integer here, so both modes start it at the same base. The
    // lookup-capped mate lies exactly at that start: from the read-length end its first base is in
    // the window and it is placed end to end; from the reference end the window starts a base
    // later, so the only candidate clips the first base and the mate keeps its retained locus.
    namespace ba = bwa_aln;
    bool ref_span = false, exact = true;
    SUBCASE("read-length end") {}
    SUBCASE("no-exact-rescue-window: read-length end") { exact = false; }
    SUBCASE("ref-span-insert: reference end") { ref_span = true; }
    FMI_search fmi((ref().fa + ".meth").c_str()); fmi.load_index(false, 1);
    ba::aln_meth_ref_t mr;
    REQUIRE(ba::aln_meth_ref_load(&mr, ref().fa.c_str(), fmi.idx->bns, 1) == 0);
    auto *opt = ba::gap_init_opt();
    auto *pe = ba::bwa_init_pe_opt();
    pe->exact_rescue_window = exact ? 1 : 0; pe->ref_span_insert = ref_span ? 1 : 0;
    ba::aln_meth_opt_t mo; ba::aln_meth_opt_finish(&mo); ba::bwase_initialize();
    const int len = 50, anchor = 1500, retained_pos = 400;
    ba::isize_info_t ii{};
    ii.avg = 200; ii.std = 10; ii.high = 300; ii.ap_prior = 1e-5;
    // The window's start from the read-length end: floor(-avg - 3 std - len / 2) before it.
    const int mate = anchor + len + (int)std::floor(-ii.avg - 3 * ii.std - len * 0.5);
    REQUIRE(mate == anchor - 205);
    ba::bwa_seq_t *seqs[2]; ba::bwa_seqio_t *inputs[2];
    for (int j = 0; j < 2; ++j) {
        const auto read = j ? aln_test::revcomp(ref().c1.substr(anchor, 25) + ref().c1.substr(anchor + 26, 25))
                            : ref().c1.substr(mate, len);
        inputs[j] = ba::bwa_seq_open(write_fq(j ? "gapped-anchor2.fq" : "gapped-anchor1.fq", {read}).c_str());
        REQUIRE(inputs[j] != nullptr);
        int n = 0; seqs[j] = ba::bwa_read_seq(inputs[j], 1, &n, BWA_MODE_COMPREAD, 0);
        REQUIRE(n == 1);
    }
    ba::aln_meth_read_t result[2]; ba::aln_meth_read_t *r[2] = {result, result + 1};
    ba::aln_meth_occ_t retained{}; retained.rid = 0; retained.hyp = 1;
    retained.pos = retained_pos; retained.alignment.rb = retained_pos; retained.alignment.re = retained_pos + len;
    retained.alignment.cigar = {(uint32_t)len << 4};
    result[0].occ.push_back(retained); result[0].main = 0; result[0].flags = ba::ALN_METH_BIND_BEST;
    ba::aln_meth_occ_t a{}; a.rid = 0; a.hyp = 1; a.strand = 1; a.pos = anchor;
    a.alignment.rb = anchor; a.alignment.re = anchor + len + 1;
    a.alignment.cigar = {25u << 4, 1u << 4 | 2, 25u << 4};
    result[1].occ.push_back(a); result[1].main = 0; result[1].c1 = 1;
    store_pair(r);
    ba::aln_meth_pe_stage1(opt, 1, seqs, r);
    uint64_t total[2] = {}, mapped[2] = {};
    ba::aln_meth_pe_stage2(&mr, &mo, opt, pe, &ii, 1, seqs, r, total, mapped);
    REQUIRE(result[0].hits.size() == 1);
    CHECK(seqs[0]->pos == (ref_span ? retained_pos : mate));
    CHECK(result[0].hits[0].pos == seqs[0]->pos);
    CHECK(cigar_of(result[0], result[0].hits[0]) == std::vector<uint32_t>{(uint32_t)len << 4});
    CHECK((int)seqs[0]->type == (ref_span ? BWA_TYPE_REPEAT : BWA_TYPE_MATESW));
    CHECK(seqs[1]->pos == (uint64_t)anchor); // the anchor stays where it was
    for (int j = 0; j < 2; ++j) {
        ba::bwa_free_read_seq(1, seqs[j]); CHECK(ba::bwa_seq_close(inputs[j]) == 0);
    }
    free(pe); free(opt); ba::aln_meth_ref_free(&mr);
}

TEST_CASE("aln_meth: cap rescue floors a forward anchor's mate at the anchor's end as mate rescue does")
{
    // A forward anchor with a one-base deletion: 50 read bases over 51 reference bases. Without
    // --short-insert, mate rescue starts a reverse mate's window no earlier than the anchor's end:
    // its read-length end (pos + len), or its reference end under --ref-span-insert. The capped
    // mate starts at pos + len, a base before the reference end: from the read-length floor it is
    // placed end to end, in either window mode; from the reference-end floor its first base lies
    // outside the window, so the only candidate clips it and the mate keeps its retained locus.
    namespace ba = bwa_aln;
    bool ref_span = false, exact = true;
    SUBCASE("read-length floor") {}
    SUBCASE("no-exact-rescue-window: read-length floor") { exact = false; }
    SUBCASE("ref-span-insert: reference-end floor") { ref_span = true; }
    FMI_search fmi((ref().fa + ".meth").c_str()); fmi.load_index(false, 1);
    ba::aln_meth_ref_t mr;
    REQUIRE(ba::aln_meth_ref_load(&mr, ref().fa.c_str(), fmi.idx->bns, 1) == 0);
    auto *opt = ba::gap_init_opt();
    auto *pe = ba::bwa_init_pe_opt();
    pe->exact_rescue_window = exact ? 1 : 0; pe->ref_span_insert = ref_span ? 1 : 0;
    ba::aln_meth_opt_t mo; ba::aln_meth_opt_finish(&mo); ba::bwase_initialize();
    const int len = 50, anchor = 1500, retained_pos = 2400, mate = anchor + len;
    ba::isize_info_t ii{};
    ii.avg = 100; ii.std = 10; ii.high = 200; ii.ap_prior = 1e-5;
    ba::bwa_seq_t *seqs[2]; ba::bwa_seqio_t *inputs[2];
    for (int j = 0; j < 2; ++j) {
        const auto read = j ? aln_test::revcomp(ref().c1.substr(mate, len))
                            : ref().c1.substr(anchor, 25) + ref().c1.substr(anchor + 26, 25);
        inputs[j] = ba::bwa_seq_open(write_fq(j ? "floor-anchor2.fq" : "floor-anchor1.fq", {read}).c_str());
        REQUIRE(inputs[j] != nullptr);
        int n = 0; seqs[j] = ba::bwa_read_seq(inputs[j], 1, &n, BWA_MODE_COMPREAD, 0);
        REQUIRE(n == 1);
    }
    ba::aln_meth_read_t result[2]; ba::aln_meth_read_t *r[2] = {result, result + 1};
    ba::aln_meth_occ_t a{}; a.rid = 0; a.hyp = 1; a.pos = anchor;
    a.alignment.rb = anchor; a.alignment.re = anchor + len + 1;
    a.alignment.cigar = {25u << 4, 1u << 4 | 2, 25u << 4};
    result[0].occ.push_back(a); result[0].main = 0; result[0].c1 = 1;
    ba::aln_meth_occ_t retained{}; retained.rid = 0; retained.hyp = 1; retained.strand = 1;
    retained.pos = retained_pos; retained.alignment.rb = retained_pos; retained.alignment.re = retained_pos + len;
    retained.alignment.cigar = {(uint32_t)len << 4};
    result[1].occ.push_back(retained); result[1].main = 0; result[1].flags = ba::ALN_METH_BIND_BEST;
    store_pair(r);
    ba::aln_meth_pe_stage1(opt, 1, seqs, r);
    uint64_t total[2] = {}, mapped[2] = {};
    ba::aln_meth_pe_stage2(&mr, &mo, opt, pe, &ii, 1, seqs, r, total, mapped);
    REQUIRE(result[1].hits.size() == 1);
    CHECK(seqs[1]->pos == (uint64_t)(ref_span ? retained_pos : mate));
    CHECK(result[1].hits[0].pos == seqs[1]->pos);
    CHECK(cigar_of(result[1], result[1].hits[0]) == std::vector<uint32_t>{(uint32_t)len << 4});
    CHECK((int)seqs[1]->type == (ref_span ? BWA_TYPE_REPEAT : BWA_TYPE_MATESW));
    CHECK(seqs[0]->pos == (uint64_t)anchor); // the anchor stays where it was
    for (int j = 0; j < 2; ++j) {
        ba::bwa_free_read_seq(1, seqs[j]); CHECK(ba::bwa_seq_close(inputs[j]) == 0);
    }
    free(pe); free(opt); ba::aln_meth_ref_free(&mr);
}

TEST_SUITE_END;

TEST_CASE("aln_meth: known alleles stop at the CIGAR cap and the next read resets"
          * doctest::test_suite("unit/aln_meth"))
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

    ba::gap_entry_t exact = {};
    auto flush = [&]() {
        ba::aln_meth_hit(&slot, &exact, row, row);
        core.state = ba::GAP_SLOT_PARKED;
        ba::aln_meth_flush_parked(&group);
    };
    flush();
    REQUIRE(results[0].occ.size() == 1);
    const auto before = results[0].occ[0];
    ba::ks_sites_t sites;
    sites.by_contig.assign(mr.bns->n_seqs, {});
    for (int p=505; p<545 && sites.n_sites<3; ++p) {
        if (ref().c1[p]!='A') continue;
        sites.by_contig[0].push_back({p,0,2});
        ++sites.n_sites;
    }
    REQUIRE(sites.n_sites == 3);
    mr.known_snvs = &sites;
    uint64_t remaining = 7;
    int expected = 1;
    SUBCASE("reference evaluation binds") { remaining = 7; expected = 1; }
    SUBCASE("first ALT evaluation binds") { remaining = 15; expected = 2; }
    SUBCASE("second ALT evaluation binds") { remaining = 23; expected = 3; }
    slot.alignment_bytes = (64ULL << 20) - remaining;
    slot.windows.clear(); slot.intervals.clear(); // re-evaluate the same locus with the seeded counter
    const auto rescored = results[0].n_rescored;
    flush();
    CHECK(results[0].n_rescored - rescored == expected);
    CHECK(results[0].flags == ba::ALN_METH_BIND_CIGAR);
    CHECK(slot.stop == 1);
    CHECK(slot.alignment_bytes == (64ULL << 20) - remaining + 8 * (expected - 1));
    CHECK(slot.pend.empty());
    REQUIRE(results[0].occ.size() == 1);
    const auto &after = results[0].occ[0];
    CHECK(std::make_tuple(after.pos, after.rid, after.hyp, after.strand, after.cost, after.ks_site,
                          after.alt, after.n_mm, after.n_gapo, after.n_gape, after.n_conv, after.n_ins, after.n_del,
                          after.alignment.rb, after.alignment.re, after.alignment.reward, after.alignment.band) ==
          std::make_tuple(before.pos, before.rid, before.hyp, before.strand, before.cost, before.ks_site,
                          before.alt, before.n_mm, before.n_gapo, before.n_gape, before.n_conv, before.n_ins, before.n_del,
                          before.alignment.rb, before.alignment.re, before.alignment.reward, before.alignment.band));
    CHECK(after.alignment.raw_cigar == before.alignment.raw_cigar);
    CHECK(after.alignment.cigar == before.alignment.cigar);
    ba::aln_meth_slot_finish(&slot);
    REQUIRE(results[0].main >= 0);
    CHECK(results[0].hits[results[0].main].pos == 500);
    ba::aln_meth_slot_start(&slot, reads+1, results+1, 1, 3, 1);
    CHECK(slot.stop == 0);
    CHECK(slot.alignment_bytes == 0);
    flush();
    CHECK(results[1].n_rescored == 4);
    CHECK(results[1].flags == 0);
    ba::aln_meth_slot_finish(&slot);
    REQUIRE(results[1].main >= 0);
    CHECK(results[1].hits[results[1].main].pos == 500);
    ba::aln_meth_ref_free(&mr);
    ba::bwt_destroy(bwt);
    free(opt);
    ba::bwa_free_read_seq(n, reads);
    CHECK(ba::bwa_seq_close(input) == 0);
}
