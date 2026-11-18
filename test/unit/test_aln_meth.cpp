// Unit tests for `aln --meth`: orientation, chemistry and tags, stock cost ranking,
// bounded candidate discovery, conservative caps, duplicate/thread invariance,
// and stock alignment against a --meth-built index.

#include "doctest/doctest.h"

#include "aln_test_ref.h"
#include "../framework/stderr_capture.h"
#include "../framework/run_in_child.h"
#include "../framework/temp_dir.h"
#include "../../src/aln/aln_known_snvs.h"
#include "../../src/aln/aln_main.h"
#include "../../src/aln/aln_meth.h"
#include "../../src/aln/aln_meth_search.h"
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
#include <functional>
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

// Turns off what slot_start set from the run's options: the stock cuts, and with them the
// certificates and the no-hit proof. For a test that drives a slot's flushes by hand (its slot has
// no job options), as the search without the cuts makes them.
void no_cuts(bwa_aln::aln_meth_slot_t &slot)
{
    slot.cut = slot.cert = slot.prove = 0;
    slot.prove2_at = INT64_MAX;
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

TEST_CASE("aln_meth: a repeat certificate reports the full X0 and draws its primary uniformly")
{
    // Under collapsed scoring read_r has 31 exact loci. With -R 5 the search ends once the best
    // stratum is complete with more than 5 of them (a repeat certificate): MAPQ 0, X0 the full count
    // (31), no X1. The primary is drawn per read: over 1,240 copies of the read every locus is drawn
    // about 40 times (chi-square, 30 degrees of freedom, p < 0.001 at 59.7).
    const int n = 31 * 40;
    const std::string fq = write_fq("cert-draw.fq", std::vector<std::string>(n, ref().read_r));
    const auto v = run_aln({"--meth=taps", "--meth-scoring", "collapsed", "-n", "2", "-l", "1024", "-R", "5"}, {fq});
    REQUIRE((int)v.size() == n);
    std::map<std::string, int> drawn;
    for (const auto &r : v) {
        CHECK(r.mapq == 0);
        CHECK(tag_value(r, "XT:A") == "R");
        CHECK(tag_value(r, "X0:i") == "31");
        CHECK(r.tag.count("X1:i") == 0);
        ++drawn[r.rname + ":" + std::to_string(r.pos) + (r.flag & 16 ? "-" : "+")];
    }
    CHECK(drawn.size() == 31u);
    double chi2 = 0;
    for (const auto &d : drawn) chi2 += (d.second - 40.0) * (d.second - 40.0) / 40.0;
    CAPTURE(chi2);
    CHECK(chi2 < 59.7);
}

TEST_CASE("aln_meth: a repeat certificate's draws do not depend on the thread count")
{
    // 96 copies of read_r (a repeat certificate under -R 5) from the command line: each copy's drawn
    // primary and X0 must match the single-threaded run. Duplicate sharing (the default) searches
    // the copies once; the per-copy searches split over threads are tested through the library
    // ("... match per-copy searches at any thread count").
    const std::string fq = write_fq("cert-threads.fq", std::vector<std::string>(96, ref().read_r));
    const std::vector<std::string> opts = {"--meth=taps", "--meth-scoring", "collapsed", "-n", "2", "-l", "1024", "-R", "5"};
    auto with = [&](const char *t) { auto o = opts; o.push_back("-t"); o.push_back(t); return run_aln(o, {fq}); };
    const auto one = with("1"), four = with("4");
    REQUIRE(one.size() == 96); REQUIRE(four.size() == one.size());
    for (size_t i = 0; i < one.size(); ++i) {
        CAPTURE(i);
        CHECK(one[i].mapq == 0);
        CHECK(std::make_tuple(one[i].flag, one[i].rname, one[i].pos, one[i].mapq, tag_value(one[i], "X0:i")) ==
              std::make_tuple(four[i].flag, four[i].rname, four[i].pos, four[i].mapq, tag_value(four[i], "X0:i")));
    }
}

TEST_CASE("aln_meth: a gapped best row's ceiling drops equal-cost mismatch loci before any certificate")
{
    // Default penalties (M 3, O 11, E 4): 31 loci where the read has five transversions (cost 15,
    // five differences) and one where it lacks two bases (cost 11 + 4 = 15, two differences). With
    // -R 5 a certificate waits for the cost-15 stratum to be complete, so the gapped row is verified
    // and its top-2 ceiling (three differences) drops the 31 mismatch loci: every copy is a unique
    // hit at the gapped locus, as stock's ceiling gives. A certificate inside the stratum would have
    // counted the mismatch loci (X0 > 5, MAPQ 0).
    bwa_tests::TempDir tmp{"aln_meth_cert_gap"};
    const std::string dir = tmp.path(), fa = dir + "/cert-gap.fa";
    uint32_t st = 2718;
    const std::string read = aln_test::lcg_dna(st, 40);
    auto tv = [](char b) { return b == 'A' ? 'T' : b == 'T' ? 'A' : b == 'C' ? 'G' : 'C'; };
    std::vector<std::pair<std::string, std::string>> contigs;
    for (int i = 0; i < 31; ++i) {
        std::string locus = read;
        for (int k = 0; k < 5; ++k) { const int at = 3 + 8 * k + i % 5; locus[at] = tv(locus[at]); }
        contigs.push_back({"mm" + std::to_string(i), aln_test::lcg_dna(st, 100) + locus + aln_test::lcg_dna(st, 100)});
    }
    contigs.push_back({"gap", aln_test::lcg_dna(st, 100) + read.substr(0, 20) + "GT" + read.substr(20) +
                              aln_test::lcg_dna(st, 100)});
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
    const int n = 256; // copies with their own draw seeds: none may draw a dropped locus
    const auto v = run_aln({"--meth", "-n", "5", "-l", "1024", "-R", "5"},
                           {write_fq_in(dir, "cert-gap.fq", std::vector<std::string>(n, read))}, fa, dir);
    REQUIRE((int)v.size() == n);
    for (const auto &r : v) {
        CAPTURE(r.rname);
        CHECK(r.rname == "gap");
        CHECK(tag_value(r, "XT:A") == "U");
        CHECK(tag_value(r, "X0:i") == "1");
        CHECK(r.mapq > 0);
    }
}

TEST_CASE("aln_meth: a certificate counts only the best rows within the top-2 ceiling")
{
    // Default penalties (M 3, O 11, E 4); every locus costs 18. gapA: the read lacks two bases and
    // has a read C over a reference T, invisible in converted space, so it is verified first, in
    // stratum 15, and sets the ceiling at four differences. Six loci with six transversions each
    // (six differences, over the ceiling: the finish drops them) are verified next, then the two
    // copies of gapB (two missing bases and a transversion: three differences), the widest hit. With
    // -R 5, counting the dropped loci would fire a repeat certificate before gapB is verified: X0 1
    // and a unique-looking MAPQ. Counting only the kept rows verifies gapB: three best loci.
    bwa_tests::TempDir tmp{"aln_meth_cert_ceiling"};
    const std::string dir = tmp.path(), fa = dir + "/cert-ceiling.fa";
    uint32_t st = 3141;
    std::string read = aln_test::lcg_dna(st, 40);
    read[30] = 'C';
    auto tv = [](char b) { return b == 'A' ? 'T' : b == 'T' ? 'A' : b == 'C' ? 'G' : 'C'; };
    std::vector<std::pair<std::string, std::string>> contigs;
    for (int i = 0; i < 6; ++i) {
        std::string mm = read;
        for (int at : {1, 7, 13, 19, 25, 33}) mm[at + i] = tv(mm[at + i]);
        contigs.push_back({"mm" + std::to_string(i), aln_test::lcg_dna(st, 100) + mm + aln_test::lcg_dna(st, 100)});
    }
    std::string gap_a = read.substr(0, 20) + "GT" + read.substr(20);
    gap_a[32] = 'T'; // read position 30
    contigs.push_back({"gapA", aln_test::lcg_dna(st, 100) + gap_a + aln_test::lcg_dna(st, 100)});
    std::string gap_b = read;
    gap_b[5] = tv(gap_b[5]);
    gap_b = gap_b.substr(0, 12) + "CA" + gap_b.substr(12);
    for (int k = 0; k < 2; ++k)
        contigs.push_back({"gapB" + std::to_string(k), aln_test::lcg_dna(st, 100) + gap_b + aln_test::lcg_dna(st, 100)});
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
    const auto v = run_aln({"--meth", "-n", "6", "-l", "1024", "-R", "5"},
                           {write_fq_in(dir, "cert-ceiling.fq", {read})}, fa, dir);
    REQUIRE(v.size() == 1);
    CHECK(v[0].rname.rfind("gap", 0) == 0);
    CHECK(tag_value(v[0], "XT:A") == "R");
    CHECK(tag_value(v[0], "X0:i") == "3");
    CHECK(v[0].mapq == 0);
}

TEST_CASE("aln_meth: a certificate waits for the best stratum to be complete")
{
    // Default penalties (M 3, O 11, E 4); every locus costs 18 and is found in stratum 18. Six
    // distinct loci with six transversions (six differences) are verified first, then one where the
    // read lacks two bases and has a transversion (three differences). The gapped row sets the top-2
    // ceiling at four differences, which drops the six mismatch loci: the read is a unique hit. With
    // -R 5, a repeat certificate taken inside stratum 18, before the gapped row is verified, would
    // count the six mismatch loci against a ceiling that is not final (X0:6, MAPQ 0).
    bwa_tests::TempDir tmp{"aln_meth_cert_stratum"};
    const std::string dir = tmp.path(), fa = dir + "/cert-stratum.fa";
    uint32_t st = 3141;
    const std::string read = aln_test::lcg_dna(st, 40);
    auto tv = [](char b) { return b == 'A' ? 'T' : b == 'T' ? 'A' : b == 'C' ? 'G' : 'C'; };
    std::vector<std::pair<std::string, std::string>> contigs;
    for (int i = 0; i < 6; ++i) {
        std::string mm = read;
        for (int at : {1, 7, 13, 19, 25, 33}) mm[at + i] = tv(mm[at + i]);
        contigs.push_back({"mm" + std::to_string(i), aln_test::lcg_dna(st, 100) + mm + aln_test::lcg_dna(st, 100)});
    }
    std::string gap = read;
    gap[5] = tv(gap[5]);
    gap = gap.substr(0, 12) + "CA" + gap.substr(12);
    contigs.push_back({"gap", aln_test::lcg_dna(st, 100) + gap + aln_test::lcg_dna(st, 100)});
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
    const auto v = run_aln({"--meth", "-n", "6", "-l", "1024", "-R", "5"},
                           {write_fq_in(dir, "cert-stratum.fq", {read})}, fa, dir);
    REQUIRE(v.size() == 1);
    CHECK(v[0].rname == "gap");
    CHECK(tag_value(v[0], "XT:A") == "U");
    CHECK(tag_value(v[0], "X0:i") == "1");
    CHECK(v[0].mapq > 0);
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

// A single-end --meth run through the library, so the test hooks of aln_opt_t can be set. tweak,
// when given, adjusts the options last, before the --meth matrices are built.
std::string run_se_api(const std::string &fq, int share_dup_search, int batch_size, const char *chem_taps,
                       const std::string &mate = "", int barcode_len = 0,
                       const std::function<void(bwa_aln::aln_opt_t &, bwa_aln::aln_meth_opt_t &)> &tweak = {})
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
    if (tweak) tweak(o, mo);
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

TEST_CASE("aln_meth: a repeat certificate's draws match per-copy searches at any thread count")
{
    // 96 copies of read_r (a repeat certificate under -R 5 and collapsed TAPS scoring). Without
    // duplicate sharing each copy is searched on its own and -t 4 splits them over threads: every
    // copy's primary and X0 must match the single-threaded run, and the shared search, whose first
    // copy makes the other copies' draws, must match both.
    const std::string fq = write_fq("cert-threads-api.fq", std::vector<std::string>(96, ref().read_r));
    auto run = [&](int share, int threads) {
        return run_se_api(fq, share, 0x40000, "taps", "", 0, [&](bwa_aln::aln_opt_t &o, bwa_aln::aln_meth_opt_t &mo) {
            o.gap.n_threads = threads;
            o.gap.max_top2 = 5;
            mo.m.meth_scoring = MEM_METH_SCORING_COLLAPSED;
        });
    };
    const std::string alone1 = run(0, 1);
    int n_rec = 0;
    std::istringstream sam(alone1);
    for (std::string line; std::getline(sam, line);) {
        if (line.empty() || line[0] == '@') continue;
        ++n_rec;
        CHECK(line.find("\tXT:A:R\t") != std::string::npos);
    }
    CHECK(n_rec == 96);
    CHECK(run(0, 4) == alone1);
    CHECK(run(1, 4) == alone1);
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
    // X0/X1 or set the MAPQ.
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
    CHECK(s2[0].mapq == 37);
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
        ba::aln_meth_opt_t mo; mo.p_best=budget; mo.p_second=0; mo.cert_budget=budget;
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
    std::vector<ba::aln_meth_extra_t> extra(n);
    for (int i=1;i<n;i+=2) extra[i].rows={{0,0,1,rb,rb+25,{}}};
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
    // 31 best loci within the ceiling are a repeat
    for (int i = 0; i < 30; ++i) {
        ba::aln_meth_occ_t o{}; o.cost = 3; o.n_mm = 1; o.n_gapo = 0; o.pos = 100000 + 1000 * i;
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
    // The budget's binding flush is the search's without stock cuts (its slot has no job options),
    // so each slot_start below turns the stock cuts off.
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
        ba::aln_meth_slot_start(&slot,reads,&result,0,3,0);no_cuts(slot);
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
    // The flushes here are the search's without stock cuts (the slot has no job options), so each
    // slot_start below turns the stock cuts off.
    ba::aln_meth_opt_t mo;ba::aln_meth_opt_finish(&mo);
    ba::gap_opt_t *opt=ba::gap_init_opt();
    ba::aln_meth_group_t group;group.bwt=bwt;group.ref=&mr;group.mo=&mo;group.opt=opt;
    group.n_slots=1;group.slot.resize(1);
    ba::gap_slot_core_t core={};group.cores=&core;
    auto &slot=group.slot[0];slot.grp=&group;slot.job=&core.job;
    ba::aln_meth_read_t result;
    ba::aln_meth_slot_start(&slot,reads,&result,0,3,1);no_cuts(slot);
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
        ba::aln_meth_slot_start(&slot,reads,&result,3,3,1);no_cuts(slot);
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
        ba::aln_meth_slot_start(&slot,reads,&result,0,50,1);no_cuts(slot);
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
    ba::aln_meth_slot_start(&slot,reads,&result,1,3,1);no_cuts(slot);
    ba::aln_meth_hit(&slot,&exact,row,row);core.state=ba::GAP_SLOT_PARKED;
    ba::aln_meth_flush_parked(&group);ba::aln_meth_slot_finish(&slot);
    CHECK(result.flags!=0);CHECK(result.main==-1);CHECK(result.c1==0);
    mo.p_best=1;mo.p_second=0;
    // An interval reached again along a path with other gaps reuses its SA rows: no lookup, so it
    // cannot bind the budget.
    ba::aln_meth_slot_start(&slot,reads,&result,2,3,1);no_cuts(slot);
    ba::aln_meth_hit(&slot,&longer,row,row);
    ba::aln_meth_hit(&slot,&exact,row,row);
    core.state=ba::GAP_SLOT_PARKED;
    ba::aln_meth_flush_parked(&group);
    CHECK(result.n_lookups==1);CHECK(result.flags==0);
    REQUIRE(result.occ.size()==1);CHECK(result.occ[0].cost==0);
    ba::aln_meth_slot_finish(&slot);
    ba::aln_meth_slot_start(&slot,reads,&result,2,3,1);no_cuts(slot);
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

TEST_CASE("aln_meth: an ALT-window candidate's path counts its seed gaps toward -k")
{
    // A 51-base read: c1[500, 550) with one extra base inserted near one end or the other. The
    // candidate gives that insertion as its path. Under -k 0 a seeded read drops the row whose
    // insertion is in the seed and keeps the other; an unseeded read keeps both.
    namespace ba = bwa_aln;
    FMI_search fmi((ref().fa + ".meth").c_str());
    fmi.load_index(false, 1);
    ba::aln_meth_ref_t mr;
    REQUIRE(ba::aln_meth_ref_load(&mr, ref().fa.c_str(), fmi.idx->bns, 1) == 0);
    auto *opt = ba::gap_init_opt();
    opt->max_seed_diff = 0;
    ba::aln_meth_opt_t mo;
    ba::aln_meth_opt_finish(&mo);
    ba::aln_meth_group_t group;
    group.ref = &mr; group.mo = &mo; group.opt = opt;
    ba::gap_slot_core_t core = {};
    ba::bwt_width_t seed_width = {};
    const int64_t rb = mr.bns->anns[0].offset + 500;
    const std::string window = ref().c1.substr(500, 50);
    for (bool seeded : {false, true}) {
        long dropped = 0;
        for (int at : {3, 47}) {
            CAPTURE(seeded); CAPTURE(at);
            std::string read = window.substr(0, at) + "A" + window.substr(at);
            const std::string fq = write_fq("seed-gap.fq", {conv(read, 'C', 'T')});
            ba::bwa_seqio_t *input = ba::bwa_seq_open(fq.c_str());
            int n = 0;
            ba::bwa_seq_t *reads = ba::bwa_read_seq(input, 1, &n, BWA_MODE_COMPREAD, 0);
            REQUIRE(n == 1);
            ba::aln_meth_slot_t slot;
            slot.grp = &group;
            ba::aln_meth_read_t result;
            ba::aln_meth_slot_start(&slot, reads, &result, 0, 3, 0);
            slot.job = &core.job;
            core.job.seed_width = seeded ? &seed_width : nullptr;
            core.job.c.seed_off = (int)reads->len - 32;
            ba::aln_meth_extra_t extra;
            extra.rows = {{0, 0, 1, rb, rb + 50, {(uint32_t)at << 4, 1u << 4 | 1, (uint32_t)(50 - at) << 4}}};
            slot.extra = &extra;
            ba::aln_meth_score_extras(&slot);
            CHECK(result.n_over_seed + (long)result.occ.size() == 1);
            dropped += result.n_over_seed;
            slot.extra = nullptr;
            slot.job = nullptr;
            ba::bwa_free_read_seq(n, reads);
            CHECK(ba::bwa_seq_close(input) == 0);
        }
        CHECK(dropped == (seeded ? 1 : 0));
    }
    free(opt);
    ba::aln_meth_ref_free(&mr);
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
        ba::aln_meth_extra_t extra;
        extra.rows = {{0, 0, 1, rb, rb + 50, {}}};
        if (path) extra.rows[0].cig = {50u << 4};
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
            ba::aln_meth_extra_t extra;
            extra.rows = {{0, 0, 1, rb, rb + 50, {}}};
            if (path) extra.rows[0].cig = {50u << 4};
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

TEST_CASE("aln_meth: the finish audits the cuts a search took on a best that was not final")
{
    namespace ba = bwa_aln;
    const std::string fq = write_fq("audit.fq", {conv(ref().c1.substr(500, 50), 'C', 'T')});
    ba::bwa_seqio_t *input = ba::bwa_seq_open(fq.c_str());
    REQUIRE(input != nullptr);
    int n = 0;
    ba::bwa_seq_t *reads = ba::bwa_read_seq(input, 1, &n, BWA_MODE_COMPREAD, 0);
    REQUIRE(n == 1);
    FMI_search fmi((ref().fa + ".meth").c_str());
    fmi.load_index(false, 1);
    ba::aln_meth_ref_t mr;
    REQUIRE(ba::aln_meth_ref_load(&mr, ref().fa.c_str(), fmi.idx->bns, 1) == 0);
    ba::ks_sites_t no_sites; // a --known-snvs run whose sites are elsewhere
    no_sites.by_contig.assign(mr.bns->n_seqs, {});
    auto *opt = ba::gap_init_opt();
    ba::aln_meth_opt_t mo;
    ba::aln_meth_opt_finish(&mo);
    ba::aln_meth_group_t group;
    group.ref = &mr; group.mo = &mo; group.opt = opt;
    ba::aln_meth_slot_t slot;
    slot.grp = &group;
    const int64_t rb = mr.bns->anns[0].offset + 500;
    // The read's one row (cost 0, no differences) comes from the ALT-window search, which stopped
    // with every window row below `horizon` found; the read retains rows up to 0 + M = 3.
    struct case_t { int64_t horizon; int ceil_applied, r_stop, want; };
    for (const case_t &c : {case_t{INT64_MAX, INT_MAX, 0, 0}, case_t{4, INT_MAX, 0, 0},
                            case_t{3, INT_MAX, 0, ba::ALN_METH_AUDIT_HORIZON},
                            case_t{4, 0, 0, ba::ALN_METH_AUDIT_CEILING}, // tighter than the final 0 + 1
                            case_t{4, 1, 0, 0},
                            case_t{4, INT_MAX, 1, ba::ALN_METH_AUDIT_REPEAT}}) // -R, yet one best locus
        for (const ba::ks_sites_t *sites : {(const ba::ks_sites_t *)nullptr, (const ba::ks_sites_t *)&no_sites}) {
            CAPTURE(c.horizon); CAPTURE(c.ceil_applied); CAPTURE(c.r_stop); CAPTURE(sites != nullptr);
            mr.known_snvs = sites;
            ba::aln_meth_read_t result;
            ba::aln_meth_slot_start(&slot, reads, &result, 0, 3, 1);
            ba::aln_meth_extra_t extra;
            extra.rows = {{0, 0, 1, rb, rb + 50, {50u << 4}}};
            extra.horizon = c.horizon;
            extra.flags = c.horizon == INT64_MAX ? 0 : ba::ALN_METH_BIND_BEST;
            slot.extra = &extra;
            slot.ceil_applied = c.ceil_applied; slot.r_stop = c.r_stop;
            ba::aln_meth_slot_finish(&slot);
            CHECK(result.audit == c.want);
            // under --known-snvs a failed audit makes the read incomplete; otherwise it is counted
            CHECK(result.flags == (c.want && sites ? ba::ALN_METH_BIND_BEST : 0));
            REQUIRE(result.main >= 0);
            CHECK(result.hits[result.main].pos == (bwtint_t)rb);
            slot.extra = nullptr;
        }
    mr.known_snvs = nullptr;
    free(opt);
    ba::aln_meth_ref_free(&mr);
    ba::bwa_free_read_seq(n, reads);
    CHECK(ba::bwa_seq_close(input) == 0);
}

TEST_CASE("aln_meth: --known-snvs takes the stock cuts in both searches but no certificates")
{
    namespace ba = bwa_aln;
    const std::string fq = write_fq("ks-cuts.fq", {conv(ref().c1.substr(500, 50), 'C', 'T')});
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
    // Only whether the run has sites matters to the slot's search settings.
    ba::ks_sites_t site_set;
    site_set.by_contig.assign(mr.bns->n_seqs, {});
    site_set.by_contig[0].push_back({510, 0, 1});
    site_set.n_sites = 1;
    const ba::ks_sites_t *sites = &site_set;
    struct expect_t { int cut, cert, prove; };
    auto start = [&](const std::string &search, const ba::ks_sites_t *known, int window, int keep_strata, expect_t want) {
        CAPTURE(search); CAPTURE(keep_strata);
        mr.known_snvs = known; mr.ks_window = window;
        group.keep_strata = keep_strata;
        ba::aln_meth_read_t result;
        ba::aln_meth_slot_start(&slot, reads, &result, 0, 3, 0);
        CHECK(slot.cut == want.cut);
        CHECK(slot.cert == want.cert);
        CHECK(slot.prove == want.prove);
        CHECK(slot.gen == window); // the ALT-window search is a candidate generator
    };
    for (int keep_strata : {0, 1}) {
        // Without --known-snvs: cuts, single-end certificates and the no-hit proof.
        start("plain", nullptr, 0, keep_strata, {1, !keep_strata, 1});
        // --known-snvs main search: the cuts and the no-hit proof, no certificates (its best is
        // final only through the ALT-window search's horizon, which the finish audits).
        start("known-snvs main", sites, 0, keep_strata, {1, 0, 1});
        // --known-snvs ALT-window search: the cuts, but no certificates and no no-hit proof.
        start("known-snvs window", nullptr, 1, keep_strata, {1, 0, 0});
    }
    // -N and a negative penalty still turn the cuts (and so the certificates) off.
    opt->mode |= BWA_MODE_NONSTOP;
    start("-N plain", nullptr, 0, 0, {0, 0, 1});
    start("-N known-snvs main", sites, 0, 0, {0, 0, 1});
    start("-N known-snvs window", nullptr, 1, 0, {0, 0, 0});
    opt->mode &= ~BWA_MODE_NONSTOP;
    opt->s_mm = -1;
    start("negative penalty known-snvs main", sites, 0, 0, {0, 0, 1});
    opt->s_mm = 3;
    // A --known-snvs read takes the cuts only when it fits an ALT window: its longest placement
    // (its length plus a deleted base per gap open and extension allowed) at most KS_WINDOW_HALF + 1
    // (161) bases. At -o 1 -e 6 (the defaults) that is 154 bases; at -o 2 -e 6, 153; at -o 0 -e 0,
    // 160. The ALT-window search keeps its cuts at any length.
    {
        const std::vector<int> lens{153, 154, 155, 160, 161};
        std::vector<std::string> seqs;
        for (int len : lens) seqs.push_back(conv(ref().c1.substr(500, len), 'C', 'T'));
        const std::string long_fq = write_fq("ks-cuts-long.fq", seqs);
        ba::bwa_seqio_t *long_input = ba::bwa_seq_open(long_fq.c_str());
        int n_long = 0;
        ba::bwa_seq_t *long_reads = ba::bwa_read_seq(long_input, (int)lens.size(), &n_long, BWA_MODE_COMPREAD, 0);
        REQUIRE(n_long == (int)lens.size());
        group.keep_strata = 0;
        struct gaps_t { int o, e, longest; };
        for (const gaps_t g : {gaps_t{1, 6, 154}, gaps_t{2, 6, 153}, gaps_t{0, 0, 160}}) {
            opt->max_gapo = g.o; opt->max_gape = g.e;
            for (int i = 0; i < n_long; ++i) {
                CAPTURE(g.o); CAPTURE(g.e); CAPTURE(lens[i]);
                const bool fits = lens[i] <= g.longest;
                CHECK(ba::ks_fits_window(lens[i], opt) == fits);
                ba::aln_meth_read_t result;
                mr.known_snvs = sites; mr.ks_window = 0;
                ba::aln_meth_slot_start(&slot, long_reads + i, &result, 0, 3, 0);
                CHECK(slot.cut == (int)fits);
                mr.known_snvs = nullptr; mr.ks_window = 1;
                ba::aln_meth_slot_start(&slot, long_reads + i, &result, 0, 3, 0);
                CHECK(slot.cut == 1);
                mr.ks_window = 0;
                ba::aln_meth_slot_start(&slot, long_reads + i, &result, 0, 3, 0);
                CHECK(slot.cut == 1); // without --known-snvs: stock's cuts at any length
            }
        }
        opt->max_gapo = 1; opt->max_gape = 6;
        ba::bwa_free_read_seq(n_long, long_reads);
        CHECK(ba::bwa_seq_close(long_input) == 0);
    }
    mr.known_snvs = nullptr; mr.ks_window = 0;
    free(opt);
    ba::aln_meth_ref_free(&mr);
    ba::bwa_free_read_seq(n, reads);
    CHECK(ba::bwa_seq_close(input) == 0);
}

TEST_CASE("aln_meth: the ALT-window generator takes only the cuts a lower final best keeps valid")
{
    namespace ba = bwa_aln;
    auto *opt = ba::gap_init_opt(); // M 3, O 11, E 4, -R 30, gap extensions count as differences
    ba::aln_meth_opt_t mo;
    ba::aln_meth_opt_finish(&mo);
    // 31 windows of nearby sites, each holding genomic position 5000 of contig 0 at its offset 100
    std::vector<ba::ks_window_t> win;
    std::vector<bntann1_t> anns(31);
    for (int k = 0; k < 31; ++k) {
        win.push_back({0, 4900, 321, 4950 - k, 0});
        anns[k] = bntann1_t();
        anns[k].offset = 321 * k; anns[k].len = 321;
    }
    bntseq_t wbns = bntseq_t();
    wbns.n_seqs = 31; wbns.anns = anns.data(); wbns.l_pac = 321 * 31;
    ba::aln_meth_ref_t mr;
    mr.bns = &wbns; mr.ks_window = 1; mr.ks_win = win.data();
    ba::aln_meth_group_t group;
    group.ref = &mr; group.mo = &mo; group.opt = opt;
    ba::gap_job_t job = ba::gap_job_t();
    job.opt = *opt;
    ba::aln_meth_slot_t slot;
    slot.grp = &group; slot.job = &job;
    ba::aln_meth_read_t result;
    slot.res = &result;
    auto reset = [&](int64_t best) {
        result = ba::aln_meth_read_t();
        job.max_diff = 5; job.score_cut = ba::GAP_SCORE_LIMIT - 1;
        slot.gen = 1; slot.cut = 1; slot.best = best;
    };
    auto row = [&](int window, int64_t cost) {
        ba::aln_meth_occ_t o;
        o.rid = window; o.hyp = 1; o.strand = 0; o.cost = cost;
        o.alignment.rb = anns[window].offset + 100;
        result.occ.push_back(o);
    };
    // A gapped best (cost 11, one difference): the score cut 14 and the ceiling it implies,
    // (11 + 3) / 3 = 4, not stock's 2, which a lower final best with more differences loosens.
    reset(11);
    row(0, 11);
    CHECK(ba::aln_meth_boundary(&slot, 12) == 0);
    CHECK(job.score_cut == 14);
    CHECK(job.max_diff == 4);
    // 31 rows at the best, all at one genomic locus: one locus for -R, so no stop.
    reset(0);
    for (int k = 0; k < 31; ++k) row(k, 0);
    CHECK(ba::aln_meth_boundary(&slot, 3) == 0);
    CHECK(job.max_diff == 1);
    CHECK(result.horizon == INT64_MAX);
    // 31 distinct loci: -R stops, but only once every row below B + M is resolved, and it records
    // where the window search stopped.
    reset(0);
    for (int k = 0; k < 31; ++k) { row(k, 0); result.occ.back().alignment.rb += k; }
    CHECK(ba::aln_meth_boundary(&slot, 2) == 0);
    CHECK(ba::aln_meth_boundary(&slot, 3) == 1);
    CHECK(result.horizon == 3);
    // Without each window's place on the genome, each window's row is a locus of its own.
    mr.ks_win = nullptr;
    reset(0);
    for (int k = 0; k < 31; ++k) row(k, 0);
    CHECK(ba::aln_meth_boundary(&slot, 3) == 1);
    mr.ks_win = win.data();
    // Its finish keeps every row within B + M (no top-2 ceiling, no trim, even when it is
    // incomplete) and states how far they are complete.
    for (int bound : {0, 1}) {
        CAPTURE(bound);
        reset(11);
        row(0, 11); row(1, 14); row(2, 15);
        result.occ[1].alignment.rb += 7; result.occ[1].n_mm = 4; // four differences, over stock's 2
        if (bound) { result.flags = ba::ALN_METH_BIND_BEST; result.horizon = 13; }
        for (auto &o : result.occ) o.pos = o.alignment.rb;
        ba::aln_meth_slot_finish(&slot);
        REQUIRE(result.hits.size() == 2);
        CHECK(result.hits[0].cost == 11);
        CHECK(result.hits[1].cost == 14);
        CHECK(result.horizon == (bound ? 13 : 15));
        CHECK(result.flags == (bound ? ba::ALN_METH_BIND_BEST : 0));
    }
    // The -m work bound ends the search between two flushes: complete only below the last
    // boundary it flushed, whatever its score cut.
    reset(11);
    row(0, 11);
    result.flags = ba::ALN_METH_BIND_M;
    slot.flushed_score = 12;
    for (auto &o : result.occ) o.pos = o.alignment.rb;
    ba::aln_meth_slot_finish(&slot);
    CHECK(result.horizon == 12);
    CHECK(result.flags == ba::ALN_METH_BIND_M);
    slot.flushed_score = -1;
    // A window row dropped as unrepresentable within the score cut: complete only below its cost,
    // and the generator is incomplete (its flags reach the main search's audit).
    reset(11);
    row(0, 11);
    result.n_unrepresentable = 1;
    slot.unrep_cost = 13;
    for (auto &o : result.occ) o.pos = o.alignment.rb;
    ba::aln_meth_slot_finish(&slot);
    CHECK(result.horizon == 13);
    CHECK(result.flags == ba::ALN_METH_BIND_UNREP);
    // Past the score cut (14) it leaves the horizon at the cut.
    reset(11);
    row(0, 11);
    result.n_unrepresentable = 1;
    slot.unrep_cost = 20;
    for (auto &o : result.occ) o.pos = o.alignment.rb;
    ba::aln_meth_slot_finish(&slot);
    CHECK(result.horizon == 15);
    CHECK(result.flags == 0);
    slot.unrep_cost = INT64_MAX;
    free(opt);
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
    ba::aln_meth_extra_t extra;
    extra.rows = {{0, 0, 1, 600, 600 + ALN_MAX_READ_LEN + 1}};
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
    slot.cert = 0; // a --known-snvs read (ALT-window candidates) takes no certificates: no cut
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
    ba::aln_meth_extra_t extra;
    extra.rows = {{0, 0, 1, 600, 600 + ALN_MAX_READ_LEN + 1}};
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
    slot.cert = 0; // a --known-snvs read (ALT-window candidates) takes no certificates: no cut
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
    ba::aln_meth_extra_t extra;
    extra.rows = {{0, 0, 1, rb, rb + 50}};
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

TEST_CASE("aln_meth: per-read ALT-window candidates survive duplicate compaction")
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
    const auto fq = write_fq("extra-compaction.fq", queries);
    for (int read_no : {0, 1}) {
        opt->seed_len = 18; // seeded; the round's unseeded compaction check below sets 1024
        auto *input = ba::bwa_seq_open(fq.c_str()); int n = 0;
        auto *reads = ba::bwa_read_seq(input, 100, &n, BWA_MODE_COMPREAD, 0); REQUIRE(n == (int)queries.size());
        std::vector<ba::aln_meth_read_t> batch(n); ba::aln_dup_state_t dup;
        // Differing candidates (here, their horizons) prevent sharing. A horizon of 0 fails the
        // finish-time audit under --known-snvs, so each read's own extra decides its flags.
        ba::ks_sites_t no_sites; no_sites.by_contig.resize(mr.bns->n_seqs); mr.known_snvs = &no_sites;
        std::vector<ba::aln_meth_extra_t> extra(n);
        for (int i = 0; i < n; ++i) extra[i].horizon = i;
        CHECK(ba::aln_meth_search_batch(bwt, &mr, &mo, n, reads, opt, 1, read_no, 0, 1, 1000,
                                      &dup, batch.data(), extra.data()) == 0);
        for (int i = 0; i < n; ++i) {
            ba::aln_meth_read_t single;
            ba::aln_meth_search_batch(bwt, &mr, &mo, 1, reads + i, opt, 1, read_no, i, 1, 1000,
                                      nullptr, &single, &extra[i]);
            CHECK(batch[i].flags == single.flags); CHECK(batch[i].c1 == single.c1); CHECK(batch[i].c2 == single.c2);
            CHECK(batch[i].n_lookups == single.n_lookups);
        }
        // A duplicate before later distinct reads makes uniq[k] nonidentity; mixed extras must
        // reach each later read as its own rather than the compacted position's.
        opt->seed_len = 1024;
        const auto swap_read = reads[1]; reads[1] = reads[n - 2]; reads[n - 2] = swap_read;
        std::vector<std::string> bases = queries; std::swap(bases[1], bases[n - 2]);
        for (int i = 0; i < n; ++i) extra[i].horizon = reads[i].len == 17 ? 0 : INT64_MAX;
        // A copy repeats an earlier read's bases and extra; every other read is a representative.
        std::vector<bool> is_copy(n, false); int n_copies = 0;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < i && !is_copy[i]; ++j)
                is_copy[i] = bases[j] == bases[i] && extra[j].horizon == extra[i].horizon;
            n_copies += is_copy[i];
        }
        REQUIRE(n_copies == 2); // the two repeated queries, one now before later distinct reads
        CHECK(is_copy[1]);
        dup = ba::aln_dup_state_t();
        CHECK(ba::aln_meth_search_batch(bwt, &mr, &mo, n, reads, opt, 1, read_no, 0, 1, 1000,
                                      &dup, batch.data(), extra.data()) == n_copies);
        int flagged = 0;
        for (int i = 0; i < n; ++i) {
            CAPTURE(i);
            ba::aln_meth_read_t single;
            ba::aln_meth_search_batch(bwt, &mr, &mo, 1, reads + i, opt, 1, read_no, i, 1, 1000,
                                      nullptr, &single, &extra[i]);
            CHECK(batch[i].flags == single.flags); CHECK(batch[i].c1 == single.c1); CHECK(batch[i].c2 == single.c2);
            flagged += batch[i].flags != 0;
            // Duplicate copies intentionally have zero work counters; a representative's are its own.
            CHECK(batch[i].n_lookups == (is_copy[i] ? 0 : single.n_lookups));
        }
        CHECK(flagged > 0); CHECK(flagged < n); // the premise: the extras make a difference
        mr.known_snvs = nullptr;
        ba::bwa_free_read_seq(n, reads); ba::bwa_seq_close(input);
        opt->seed_len = 18;
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
    ba::aln_meth_extra_t empty; slot.extra = &empty;
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

TEST_CASE("aln_meth: the known-site search preserves independent REF plus ALT window strata")
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
    // This competitor must survive even though the genomic search's cuts alone could omit it.
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
        ba::aln_opt_t o; ba::aln_opt_init(&o); o.gap = *opt; o.gap.n_threads = 1; o.max_xa_hits = 1000;
        o.meth = &mo; o.meth_ref = &mr; o.known_snvs = &sites; o.ks_windows = &windows;
        ba::aln_meth_read_t result;
        ba::ks_meth_search(&o, bwt, 1, reads, read_no, 0, 1, nullptr, nullptr, &result);
        CHECK(result.flags == 0);
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

TEST_CASE("aln_meth: seeded known-site searches keep their retained loci")
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
        // The --known-snvs search, and as its control the search without the sites.
        ba::aln_opt_t o;ba::aln_opt_init(&o);o.gap=*opt;o.gap.n_threads=1;o.max_xa_hits=1000;
        o.meth=&mo;o.meth_ref=&ref;o.known_snvs=&sites;o.ks_windows=&windows;
        ba::aln_meth_read_t certified,uncertified;
        ba::ks_meth_search(&o,bwt,1,reads,read_no,0,1,nullptr,nullptr,&certified);
        ref.known_snvs=nullptr;
        ba::aln_meth_search_batch(bwt,&ref,&mo,1,reads,opt,1,read_no,0,1,1000,nullptr,&uncertified);
        ref.known_snvs=&sites;
        CHECK(certified.flags==0);CHECK(uncertified.flags==0);
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


TEST_CASE("aln_meth: audit-flagged window searches undo neutral selection and split duplicate sharing")
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
    // Each read's ALT-window search: complete, or stopped with these flags below every cost (a
    // horizon of 0), which the finish-time audit hands on to the read under --known-snvs.
    std::vector<int> flags{0, ba::ALN_METH_BIND_M, 0, ba::ALN_METH_BIND_BEST,
        ba::ALN_METH_BIND_M, ba::ALN_METH_BIND_CIGAR, ba::ALN_METH_BIND_UNREP, 0};
    std::vector<ba::aln_meth_extra_t> extra(n);
    for (int i = 0; i < n; ++i) if (flags[i]) { extra[i].flags = flags[i]; extra[i].horizon = 0; }
    ba::ks_sites_t no_sites; no_sites.by_contig.resize(mr.bns->n_seqs); mr.known_snvs = &no_sites;
    std::vector<ba::aln_meth_read_t> results(n); ba::aln_dup_state_t dup;
    CHECK(ba::aln_meth_search_batch(bwt, &mr, &mo, n, reads, opt, 1, 1, 0, 1, 1000,
          &dup, results.data(), extra.data()) == 3);
    for (int i = 0; i < n; ++i) {
        CAPTURE(i);
        CHECK(results[i].flags == flags[i]);
        // (a duplicate copy ran no search: its audit, like its work counters, stays with its first copy)
        if (results[i].n_lookups) CHECK(results[i].audit == (flags[i] ? ba::ALN_METH_AUDIT_HORIZON : 0));
        if (flags[i]) { CHECK(results[i].n_conv_moved == 0); CHECK(results[i].conv_gap == 0); }
        else CHECK(results[i].n_conv_moved == 3);
        ba::aln_meth_read_t single;
        ba::aln_meth_search_batch(bwt, &mr, &mo, 1, reads + i, opt, 1, 1, i, 1, 1000,
            nullptr, &single, &extra[i]);
        CHECK(results[i].main == single.main); CHECK(results[i].c1 == single.c1); CHECK(results[i].c2 == single.c2);
        CHECK(results[i].n_conv_moved == single.n_conv_moved);
    }
    // Reset with no window candidates restores complete-state tie-breaking for every reused result.
    dup = ba::aln_dup_state_t();
    ba::aln_meth_search_batch(bwt, &mr, &mo, n, reads, opt, 1, 1, 0, 1, 1000, &dup, results.data());
    for (const auto &r : results) { CHECK(r.flags == 0); CHECK(r.n_conv_moved == 3); }
    mr.known_snvs = nullptr;
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
        std::vector<ba::aln_meth_extra_t> extra(2);
        extra[0].rows.push_back({c4, 0, 1, rb, rb + len, ungapped});
        extra[1].rows.push_back({c4, 0, 1, rb, rb + len, same ? ungapped : gapped});
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

TEST_CASE("aln_meth: ks_meth_search undoes the neutral tie-break of a read its audit flags")
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
    auto window_search = [&]() {
        ba::aln_meth_read_t w;
        ba::aln_meth_search_batch(windows.bwt, windows.mref, &mo, 1, reads, &o.gap, 1, 1, 0, 1, o.max_xa_hits,
                                  nullptr, &w);
        return w;
    };
    // The genomic search given the window search's candidates as complete (no flags, no horizon):
    // what the read would get if the window flags were only OR-ed into the result afterwards.
    auto genomic_without_flags = [&](const ba::aln_meth_read_t &w) {
        ba::aln_meth_extra_t extra; ba::ks_meth_candidates(&windows, mr.bns, w, &extra);
        extra.flags = 0; extra.horizon = INT64_MAX;
        ba::aln_meth_read_t r;
        ba::aln_meth_search_batch(bwt, &mr, &mo, 1, reads, &o.gap, 1, 1, 0, 0, o.max_xa_hits, nullptr, &r, &extra);
        return r;
    };
    // Calibrate: with the budget set to the rows the genomic search needs, the window search (more
    // rows) is bound while the genomic search alone completes.
    const auto full_window = window_search();
    REQUIRE(full_window.flags == 0);
    const auto full_genomic = genomic_without_flags(full_window);
    REQUIRE(full_genomic.flags == 0);
    REQUIRE(full_window.n_lookups > full_genomic.n_lookups);
    mo.p_best = (int)full_genomic.n_lookups; mo.p_second = 0;
    const auto bound_window = window_search();
    REQUIRE(bound_window.flags == ba::ALN_METH_BIND_BEST);
    REQUIRE(bound_window.horizon == 0); // bound in the read's best stratum
    // Flags OR-ed in after the search would leave neutral scoring's conversion tie-break in place.
    const auto or_after = genomic_without_flags(bound_window);
    REQUIRE(or_after.flags == 0); REQUIRE(or_after.n_conv_moved == 3);
    // ks_meth_search's audit finds the window search stopped below the read's costs: the read is
    // incomplete, so the tie-break (which assumes every equal-cost locus was seen) is undone.
    ba::aln_meth_read_t res;
    ba::ks_meth_search(&o, bwt, 1, reads, 1, 0, 0, nullptr, nullptr, &res);
    CHECK(res.flags == ba::ALN_METH_BIND_BEST);
    CHECK(res.audit == ba::ALN_METH_AUDIT_HORIZON);
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

TEST_CASE("aln_meth: a path's seed gap differences are the same from its gap events and from its layout")
{
    namespace ba = bwa_aln;
    uint32_t st = 4242;
    auto next = [&](int n) { st = st * 1103515245u + 12345u; return (int)((st >> 8) % (uint32_t)n); };
    int checked = 0;
    for (int trial = 0; trial < 2000; ++trial) {
        const int len = 20 + next(140), n_gapo = 1 + next(ba::ALN_METH_MAX_GAPS);
        uint32_t gap[ba::ALN_METH_MAX_GAPS];
        for (int g = 0; g < n_gapo; ++g) {
            const int del = next(2), gl = 1 + next(4);
            // an insertion's S index i covers S bases [i - gl + 1, i]; a deletion follows S base i
            const int i = del ? next(len - 1) : gl - 1 + next(len - gl + 1);
            gap[g] = ALN_METH_GAP(i, gl, del);
        }
        const int seed_off = next(len);
        for (int strand = 0; strand < 2; ++strand)
            for (int gape = 0; gape < 2; ++gape) {
                CAPTURE(trial); CAPTURE(strand); CAPTURE(gape);
                CHECK(ba::aln_meth_layout_seed_diffs(gap, n_gapo, len, strand, seed_off, gape) ==
                      ba::aln_meth_seed_gap_diffs(gap, n_gapo, seed_off, gape));
                ++checked;
            }
    }
    CHECK(checked == 8000);
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
    slot.cert = 0; // the lookup-budget flush: a repeat certificate would end this read's row work
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
