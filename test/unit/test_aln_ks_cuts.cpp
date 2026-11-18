// `bwa-mem3 aln --known-snvs`: the search's cuts against an oracle. Each case builds a small
// reference with planted loci and known sites, and compares the read's record (placement, MAPQ,
// X0, X1 and XA) with the same aligner run without --known-snvs on the reference with the read's
// ALT allele substituted. For a locus that needs at most one ALT credited, the two must agree,
// whichever search (main or ALT window) proposes the locus and in whichever order.

#include "doctest/doctest.h"

#include "aln_test_ref.h"
#include "../framework/stderr_capture.h"
#include "../framework/temp_dir.h"
#include "../../src/aln/aln_main.h"
#include "../../src/bwa.h"

#include <getopt.h>
#include <algorithm>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

int bwa_index(int argc, char *argv[]); // src/bwtindex.cpp, the `index` command

namespace {

// A random sequence without C: under --meth (read 1, C->T) such a read converts to itself, so its
// costs are the same with and without --meth.
std::string no_c(uint32_t &st, int n)
{
    std::string s = aln_test::lcg_dna(st, n);
    for (char &c : s)
        if (c == 'C') c = "AGT"[(st = st * 1103515245u + 12345u) >> 16 & 1 ? 1 : 0];
    return s;
}

// A reference base that is a mismatch for read base b, with and without --meth: not b and not C
// (a read T over a reference C is a conversion under --meth).
char mm_base(char b) { return b == 'A' ? 'G' : b == 'G' ? 'T' : 'A'; }

// The read with a mismatch planted at each of `at` (read positions).
std::string with_mm(std::string s, const std::vector<int> &at)
{
    for (int i : at) s[i] = mm_base(s[i]);
    return s;
}

struct site_t {
    int pos;           // 0-based on c1
    char ref;          // the reference base
    std::string alts;  // comma-separated, as in the VCF
    char oracle;       // the base the oracle reference holds there (ref: unchanged)
};

/* One case: c1 (the reference), its sites, and the reads. Builds the reference and the oracle
 * reference (c1 with each site's `oracle` base in place), each with a plain and a --meth index. */
struct KsCase {
    bwa_tests::TempDir tmp{"aln_ks_cuts_unit"};
    std::string dir = tmp.path(), ref, oracle, vcf, fq;

    KsCase(const std::string &c1, const std::vector<site_t> &sites, const std::vector<std::string> &reads)
    {
        std::string alt = c1;
        for (const site_t &s : sites) {
            REQUIRE(c1[s.pos] == s.ref);
            alt[s.pos] = s.oracle;
        }
        ref = build(dir + "/ref.fa", c1);
        oracle = build(dir + "/oracle.fa", alt);
        vcf = dir + "/sites.vcf";
        FILE *fp = fopen(vcf.c_str(), "w");
        REQUIRE(fp != nullptr);
        fprintf(fp, "##fileformat=VCFv4.2\n##contig=<ID=c1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        std::vector<site_t> sorted = sites;
        std::sort(sorted.begin(), sorted.end(), [](const site_t &a, const site_t &b) { return a.pos < b.pos; });
        for (const site_t &s : sorted) fprintf(fp, "c1\t%d\t.\t%c\t%s\t.\tPASS\t.\n", s.pos + 1, s.ref, s.alts.c_str());
        REQUIRE(fclose(fp) == 0);
        fq = dir + "/r.fq";
        fp = fopen(fq.c_str(), "w");
        REQUIRE(fp != nullptr);
        for (size_t i = 0; i < reads.size(); ++i)
            fprintf(fp, "@r%zu\n%s\n+\n%s\n", i, reads[i].c_str(), std::string(reads[i].size(), 'I').c_str());
        REQUIRE(fclose(fp) == 0);
    }

    static std::string build(const std::string &fa, const std::string &c1)
    {
        FILE *fp = fopen(fa.c_str(), "w");
        REQUIRE(fp != nullptr);
        fprintf(fp, ">c1\n");
        for (size_t i = 0; i < c1.size(); i += 60) fprintf(fp, "%s\n", c1.substr(i, 60).c_str());
        REQUIRE(fclose(fp) == 0);
        bwa_tests::StderrCapture quiet;
        REQUIRE(bwa_idx_build(fa.c_str(), fa.c_str(), 0, 3) == 0);
        std::vector<std::string> a = {"index", "--meth", fa};
        std::vector<char *> argv;
        for (std::string &s : a) argv.push_back(&s[0]);
        argv.push_back(nullptr);
        optind = 1;
#if defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || defined(__NetBSD__)
        optreset = 1;
#endif
        REQUIRE(bwa_index((int)a.size(), argv.data()) == 0);
        return fa;
    }

    // aln with `args` on `reference` (the case's reads); the records' fields, header lines skipped
    std::vector<std::vector<std::string>> run(std::vector<std::string> args, const std::string &reference,
                                              const std::string &tag) const
    {
        const std::string out = dir + "/" + tag + ".sam";
        std::vector<std::string> a = {"aln", "-f", out};
        a.insert(a.end(), args.begin(), args.end());
        a.push_back(reference);
        a.push_back(fq);
        std::vector<char *> argv;
        for (std::string &s : a) argv.push_back(&s[0]);
        argv.push_back(nullptr);
        bwa_tests::StderrCapture quiet;
        REQUIRE(main_aln((int)a.size(), argv.data()) == 0);
        std::vector<std::vector<std::string>> recs;
        FILE *fp = fopen(out.c_str(), "r");
        REQUIRE(fp != nullptr);
        std::string line;
        for (int c; (c = fgetc(fp)) != EOF;) {
            if (c != '\n') { line += (char)c; continue; }
            if (!line.empty() && line[0] != '@') {
                std::vector<std::string> f;
                size_t a0 = 0;
                for (size_t b; (b = line.find('\t', a0)) != std::string::npos; a0 = b + 1) f.push_back(line.substr(a0, b - a0));
                f.push_back(line.substr(a0));
                recs.push_back(f);
            }
            line.clear();
        }
        fclose(fp);
        return recs;
    }
};

std::string tag(const std::vector<std::string> &r, const std::string &key) // e.g. "X1:i"
{
    for (size_t i = 11; i < r.size(); ++i)
        if (r[i].compare(0, key.size() + 1, key + ":") == 0) return r[i].substr(key.size() + 1);
    return "-";
}

// What the cases compare: strand and mapped flags, contig, position, CIGAR, MAPQ, X0, X1 and the
// XA alternatives (without their NM) as a set.
std::string summary(const std::vector<std::string> &r)
{
    std::string xa = tag(r, "XA:Z");
    std::vector<std::string> alts;
    for (size_t a = 0, b; xa != "-" && a < xa.size(); a = b + 1) {
        b = xa.find(';', a);
        if (b == std::string::npos) b = xa.size();
        // contig, strand and position, CIGAR; not NM, which counts the site's ALT as an edit
        // against the reference but not against the oracle
        if (b > a) alts.push_back(xa.substr(a, xa.rfind(',', b - 1) - a));
    }
    std::sort(alts.begin(), alts.end());
    std::string s = std::to_string(std::stoi(r[1]) & 0x14) + " " + r[2] + ":" + r[3] + " " + r[5] + " MAPQ " + r[4] +
                    " X0 " + tag(r, "X0:i") + " X1 " + tag(r, "X1:i") + " XA";
    for (const std::string &a : alts) s += " " + a;
    return s;
}

// The reads' summaries under --known-snvs on the case's reference, and the oracle's (no
// --known-snvs, the oracle reference), with `opts` (plus --meth when meth).
void check_against_oracle(const KsCase &k, bool meth, std::vector<std::string> opts, const std::vector<std::string> &want)
{
    if (meth) opts.insert(opts.begin(), "--meth");
    std::vector<std::string> ks_opts = opts;
    ks_opts.push_back("--known-snvs");
    ks_opts.push_back(k.vcf);
    const auto got = k.run(ks_opts, k.ref, meth ? "ks.meth" : "ks");
    const auto oracle = k.run(opts, k.oracle, meth ? "oracle.meth" : "oracle");
    REQUIRE(got.size() == oracle.size());
    REQUIRE(got.size() == want.size());
    for (size_t i = 0; i < got.size(); ++i) {
        CAPTURE(i);
        CHECK(summary(got[i]) == summary(oracle[i]));
        CHECK(summary(got[i]) == want[i]);
    }
}

std::string pos1(int p) { return std::to_string(p + 1); }

/* The reviewer's case: a 100-base read R with no C, and three loci within -n 5:
 *  A: R with one deleted base (50M1D50M): cost 11, one difference
 *  B: R with four mismatches, one of them a known site whose ALT is R's base: cost 9 with the ALT
 *     (three differences), 12 without
 *  C: R with four mismatches, no site: cost 12, four differences
 * The best is B with the ALT (9, three differences), so stock's ceiling is 4 and A and C are both
 * retained: X0 1, X1 2, MAPQ 20, XA A and C. A search that takes its cuts from A (its own best)
 * sets the ceiling to 2 and never proposes C. copies_of_a: how many copies of A (31 > -R 30). */
struct Repro {
    std::string c1, read;
    int a = 6000, b = 1000, c = 3000, b_site = 0;
    char b_ref = 0;
    explicit Repro(int copies_of_a)
    {
        uint32_t st = 20261112;
        read = no_c(st, 100);
        c1 = aln_test::lcg_dna(st, 6000 + 400 * copies_of_a + 2000);
        // B: mismatches at 50, 60, 70, 80 (outside the 32-base seed); the site is at 60
        const std::string seg_b = with_mm(read, {50, 60, 70, 80});
        c1.replace(b, 100, seg_b);
        b_site = b + 60;
        b_ref = seg_b[60];
        c1.replace(c, 100, with_mm(read, {45, 55, 65, 75}));
        // A: one reference base after read base 49 that neither neighbour equals (no shifted copy)
        char extra = 'A';
        while (extra == read[49] || extra == read[50]) extra = extra == 'A' ? 'G' : 'T';
        const std::string seg_a = read.substr(0, 50) + extra + read.substr(50);
        for (int k = 0; k < copies_of_a; ++k) c1.replace(a + 400 * k, 101, seg_a);
    }
    site_t site() const { return {b_site, b_ref, std::string(1, read[60]), read[60]}; }
    std::string want(int copies_of_a) const
    {
        std::string s = "0 c1:" + pos1(b) + " 100M MAPQ " + (copies_of_a == 1 ? "20" : "8") + " X0 1 X1 " +
                        std::to_string(copies_of_a + 1) + " XA";
        if (copies_of_a == 1) s += " c1,+" + pos1(c) + ",100M c1,+" + pos1(a) + ",50M1D50M"; // as sorted
        return s;
    }
};

void repro_case(int copies_of_a, bool meth)
{
    const Repro r(copies_of_a);
    const KsCase k(r.c1, {r.site()}, {r.read});
    check_against_oracle(k, meth, {"-n", "5"}, {r.want(copies_of_a)});
}

/* The ALT-window search's own best, a gapped ALT row, must not set its ceiling:
 * G: R with a deleted base and, at a site, a mismatch whose ALT is R's base: 11 with the ALT (one
 *    difference), 14 without. It is the ALT-window search's best.
 * H: R with three mismatches, no site: 9, the read's best (ceiling 4, retained cost 12).
 * I: R with five mismatches, one at a site whose ALT is R's base: 12 with the ALT (four
 *    differences), 15 without. Only the ALT-window search can propose it, and only if it does not
 *    cut at its own best's ceiling (2). */
void window_ceiling_case(bool meth)
{
    uint32_t st = 77;
    const std::string read = no_c(st, 100);
    std::string c1 = aln_test::lcg_dna(st, 12000);
    const int g = 1000, h = 4000, i = 7000;
    std::string seg_g = with_mm(read, {70});
    char extra = 'A';
    while (extra == read[49] || extra == read[50]) extra = extra == 'A' ? 'G' : 'T';
    seg_g = seg_g.substr(0, 50) + extra + seg_g.substr(50);
    c1.replace(g, 101, seg_g);
    c1.replace(h, 100, with_mm(read, {45, 60, 75}));
    const std::string seg_i = with_mm(read, {40, 50, 60, 70, 80});
    c1.replace(i, 100, seg_i);
    const int g_site = g + 71, i_site = i + 60; // seg_g's base 71 is read base 70
    const KsCase k(c1, {{g_site, seg_g[71], std::string(1, read[70]), read[70]},
                        {i_site, seg_i[60], std::string(1, read[60]), read[60]}}, {read});
    const std::string want = "0 c1:" + pos1(h) + " 100M MAPQ 20 X0 1 X1 2 XA c1,+" + pos1(g) +
                             ",50M1D50M c1,+" + pos1(i) + ",100M";
    check_against_oracle(k, meth, {"-n", "5"}, {want});
}

/* Overlapping ALT windows count one locus once toward -R:
 * P: R exactly, with 31 sites just before it, so each of their ALT windows holds P with the
 *    reference bases: 31 window rows at cost 0, all one locus.
 * Q: R with two mismatches, one at a site whose ALT is R's base: 3 with the ALT.
 * A window search that counts window rows stops at cost 0 on -R and never proposes Q. */
void window_overlap_case(bool meth)
{
    uint32_t st = 4242;
    const std::string read = no_c(st, 100);
    std::string c1 = aln_test::lcg_dna(st, 8000);
    const int p = 2000, q = 5000;
    c1.replace(p, 100, read);
    const std::string seg_q = with_mm(read, {50, 70});
    c1.replace(q, 100, seg_q);
    std::vector<site_t> sites;
    for (int s = p - 61; s < p; s += 2) sites.push_back({s, c1[s], std::string(1, mm_base(c1[s] == 'C' ? 'A' : c1[s])), c1[s]});
    REQUIRE(sites.size() == 31);
    sites.push_back({q + 50, seg_q[50], std::string(1, read[50]), read[50]});
    const KsCase k(c1, sites, {read});
    const std::string want = "0 c1:" + pos1(p) + " 100M MAPQ 23 X0 1 X1 1 XA c1,+" + pos1(q) + ",100M";
    check_against_oracle(k, meth, {"-n", "5"}, {want});
}

} // namespace

TEST_SUITE_BEGIN("unit/aln_ks_cuts");

TEST_CASE("aln_ks_cuts: --meth: the ALT-window search does not take a ceiling from its own gapped best")
{
    window_ceiling_case(true);
}

TEST_CASE("aln_ks_cuts: --meth: overlapping ALT windows count one locus once toward -R")
{
    window_overlap_case(true);
}

TEST_CASE("aln_ks_cuts: --meth: a cheaper ALT-window row loosens the ceiling the main search took from its own best")
{
    repro_case(1, true);
}

TEST_CASE("aln_ks_cuts: --meth: -R counts the final best, not the main search's own best stratum")
{
    // 31 copies of A: the main search's own best stratum holds 31 loci, more than -R 30, but the
    // read's best is B (one locus), so stock reports all 32 worse loci.
    repro_case(31, true);
}

TEST_CASE("aln_ks_cuts: a CpG allele preserves an eligible seed at equal whole-read cost")
{
    uint32_t st = 20261114;
    std::string read = no_c(st, 100); read[30] = 'T'; read[31] = 'A';
    std::string c1 = aln_test::lcg_dna(st, 5000);
    const int a = 2000;
    c1.replace(a, read.size(), read); c1[a + 30] = 'C';
    const KsCase k(c1, {{a + 31, 'A', "G", 'G'}}, {read});
    std::vector<std::string> opts = {"--meth=taps", "--meth-context", "cpg", "--meth-scoring", "genomic", "-n", "3", "-l", "32", "-k", "0"};
    auto ksopts = opts; ksopts.insert(ksopts.end(), {"--known-snvs", k.vcf});
    const auto got = k.run(ksopts, k.ref, "cpg.seed.ks");
    const auto oracle = k.run(opts, k.oracle, "cpg.seed.oracle");
    REQUIRE(got.size() == 1); REQUIRE(oracle.size() == 1);
    CAPTURE(summary(got[0])); CAPTURE(summary(oracle[0]));
    CHECK((std::stoi(oracle[0][1]) & 4) == 0);
    CHECK(summary(got[0]) == summary(oracle[0]));
}

TEST_CASE("aln_ks_cuts: --meth credits an ALT at one site and the reference at another")
{
    // M: R with four mismatches and two sites: at s1 the reference differs from R and the ALT is
    // R's base; at s2 the reference is R's base (the ALT is not). Crediting s1's ALT and keeping
    // s2's reference leaves four mismatches (cost 12, MAPQ 37 under -n 5); neither the reference
    // nor every ALT in place does (five mismatches, MAPQ 25).
    uint32_t st = 31337;
    const std::string read = no_c(st, 100);
    std::string c1 = aln_test::lcg_dna(st, 4000);
    const int m = 1500;
    const std::string seg = with_mm(read, {40, 50, 70, 80, 60});
    c1.replace(m, 100, seg);
    const int s1 = m + 60, s2 = m + 65;
    const KsCase k(c1, {{s1, seg[60], std::string(1, read[60]), read[60]},
                        {s2, seg[65], std::string(1, mm_base(seg[65])), seg[65]}}, {read});
    const std::string want = "0 c1:" + pos1(m) + " 100M MAPQ 37 X0 1 X1 0 XA";
    check_against_oracle(k, true, {"-n", "5"}, {want});
    // XV: neither site counts against the placement (each column matches one of its alleles)
    const auto got = k.run({"--meth", "-n", "5", "--known-snvs", k.vcf}, k.ref, "xv");
    REQUIRE(got.size() == 1);
    const std::string xv = tag(got[0], "XV:Z");
    CAPTURE(xv);
    CHECK(xv.find("c1:" + pos1(s1) + ":" + seg[60] + ">" + read[60] + ",A,4,") != std::string::npos);
    CHECK(xv.find("c1:" + pos1(s2) + ":" + seg[65] + ">" + mm_base(seg[65]) + ",R,4,") != std::string::npos);
}

TEST_CASE("aln_ks_cuts: --meth credits the second ALT of a multi-allelic site")
{
    // R carries the second ALT (in VCF order and in code order) of a two-ALT site, plus four
    // mismatches: five differences against the reference and against the first ALT, four with
    // the second ALT, which -n 4 allows (MAPQ 25: all of -n used).
    uint32_t st = 99;
    std::string read;
    int at = -1;
    while (at < 0) {
        read = no_c(st, 100);
        for (int q = 55; q < 66 && at < 0; ++q)
            if (read[q] == 'T') at = q;
    }
    std::string c1 = aln_test::lcg_dna(st, 4000);
    const int m = 1200;
    std::string seg = with_mm(read, {40, 48, 72, 80});
    seg[at] = 'A'; // the reference: A; the ALTs G and T (R's base)
    c1.replace(m, 100, seg);
    const KsCase k(c1, {{m + at, 'A', "G,T", 'T'}}, {read});
    const std::string want = "0 c1:" + pos1(m) + " 100M MAPQ 25 X0 1 X1 0 XA";
    check_against_oracle(k, true, {"-n", "4"}, {want});
}

TEST_SUITE_END();
