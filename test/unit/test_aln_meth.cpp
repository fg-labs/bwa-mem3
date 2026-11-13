// Unit tests for `aln --meth`: orientation, chemistry and tags, reward ranking,
// bounded candidate discovery, conservative caps, duplicate/thread invariance,
// and stock alignment against a --meth-built index.

#include "doctest/doctest.h"

#include "aln_test_ref.h"
#include "../framework/scoped_env.h"
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
    int flag = 0, pos = 0, mapq = 0;
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

std::vector<Rec> run_aln(std::vector<std::string> args, const std::vector<std::string> &inputs)
{
    const std::string out = ref().dir + "/out.sam";
    std::vector<std::string> a = {"aln", "-f", out};
    a.insert(a.end(), args.begin(), args.end());
    a.push_back(ref().fa);
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
        const char context = genome[pos + step] == partner ? 'z' : genome[pos + 2 * step] == partner ? 'x' : 'h';
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
std::string run_se_api(const std::string &fq, int share_dup_search, int batch_size, const char *chem_taps)
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
    REQUIRE(ba::aln_run_se_meth(&o, bwt, fq.c_str(), out) == 0);
    REQUIRE(ba::aln_out_close(out) == 0);
    ba::aln_meth_ref_free(&mr);
    ba::bwt_destroy(bwt);
    std::ifstream in(out_path);
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

} // namespace

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
    // the stats line of a run (BWA_ALN_METH_STATS): SA lookups, windows scored, flushes
    auto work = [](const std::string &fq) {
        bwa_tests::ScopedEnv stats("BWA_ALN_METH_STATS", "1");
        bwa_tests::StderrCapture cap;
        run_se_api(fq, 1, 0x40000, nullptr);
        const std::string err = cap.text();
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
    const std::string out = ref().dir + "/homopolymer.sam";
    auto run = [&](std::vector<std::string> mode) {
        std::vector<std::string> a = {"aln", "-f", out};
        a.insert(a.end(), mode.begin(), mode.end());
        a.insert(a.end(), {fa, write_fq("polyc.fq", {polyc})});
        std::vector<char *> v;
        for (auto &s : a) v.push_back(&s[0]);
        v.push_back(nullptr);
        REQUIRE(main_aln((int)a.size(), v.data()) == 0);
        std::ifstream in(out);
        std::vector<int> flags;
        for (std::string l; std::getline(in, l);)
            if (!l.empty() && l[0] != '@') flags.push_back(atoi(l.substr(l.find('\t') + 1).c_str()));
        return flags;
    };
    for (const char *mode : {"--meth=taps", "--meth=5base"}) {
        CAPTURE(mode);
        // No representable candidate, so the read is reported unmapped.
        const auto flags = run({mode});
        REQUIRE(flags.size() == 1);
        CHECK((flags[0] & 4) != 0);
    }
    // EM-seq's collapsed scoring frees read C over reference T, so the same read still places
    // (genomic scoring, the default, counts it as TAPS does).
    const auto em = run({"--meth=emseq", "--meth-scoring=collapsed"});
    REQUIRE(em.size() == 1);
    CHECK((em[0] & 4) == 0);
}
