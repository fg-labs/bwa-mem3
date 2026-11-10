// `bwa-mem3 aln --known-snvs` (src/aln/aln_known_snvs.cpp): loading the sites, the errors that
// name a bad record, and output that does not change when no read is affected.

#include "doctest/doctest.h"

#include "aln_test_ref.h"
#include "../framework/stderr_capture.h"
#include "../../src/aln/aln_main.h"
#include "../../src/bwa.h"

#include <cstdio>
#include <string>
#include <vector>
#include <unistd.h>
#include <getopt.h>

int bwa_index(int argc, char *argv[]); // src/bwtindex.cpp, the `index` command

namespace {

struct KsIndex {
    bwa_tests::TempDir tmp_dir{"aln_known_snvs_unit"};
    std::string dir = tmp_dir.path(), prefix, fq;
    std::string c1, c2;

    KsIndex()
    {
        prefix = dir + "/ref.fa";
        REQUIRE(aln_test::write_test_fasta(prefix));
        REQUIRE(bwa_idx_build(prefix.c_str(), prefix.c_str(), 0, 3) == 0);
        const auto contigs = aln_test::test_contigs();
        c1 = contigs[0].second;
        c2 = contigs[1].second;
        fq = dir + "/r.fq";
        FILE *fp = fopen(fq.c_str(), "w");
        REQUIRE(fp != nullptr);
        for (int i = 0; i < 6; ++i)
            fprintf(fp, "@r%d\n%s\n+\n%s\n", i, c1.substr(80 + 30 * i, 32).c_str(), std::string(32, 'I').c_str());
        REQUIRE(fclose(fp) == 0);
    }
};

KsIndex &fixture()
{
    static KsIndex idx;
    return idx;
}

int run(std::vector<std::string> args)
{
    std::vector<char *> argv;
    std::string name = "aln";
    argv.push_back(&name[0]);
    for (std::string &s : args) argv.push_back(&s[0]);
    argv.push_back(nullptr);
    return main_aln((int)argv.size() - 1, argv.data());
}

std::string slurp(const std::string &path)
{
    std::string s;
    FILE *fp = fopen(path.c_str(), "r");
    if (!fp) return s;
    char buf[4096];
    size_t n;
    while ((n = fread(buf, 1, sizeof buf, fp)) > 0) s.append(buf, n);
    fclose(fp);
    return s;
}

// The records and every header line but @PG (it holds the command line).
std::string body(const std::string &sam)
{
    std::string out;
    size_t i = 0;
    while (i < sam.size()) {
        size_t j = sam.find('\n', i);
        if (j == std::string::npos) j = sam.size();
        if (sam.compare(i, 3, "@PG") != 0) out.append(sam, i, j - i + 1);
        i = j + 1;
    }
    return out;
}

std::string write_vcf(const std::string &path, const std::vector<std::string> &rows)
{
    FILE *fp = fopen(path.c_str(), "w");
    REQUIRE(fp != nullptr);
    fprintf(fp, "##fileformat=VCFv4.2\n##contig=<ID=c1>\n##contig=<ID=c2>\n##contig=<ID=c3>\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
    for (const std::string &r : rows) fprintf(fp, "%s\n", r.c_str());
    REQUIRE(fclose(fp) == 0);
    return path;
}

char other(char b, char not_this = 0)
{
    for (char c : {'A', 'C', 'G', 'T'})
        if (c != b && c != not_this) return c;
    return 'N';
}

std::string revcomp(const std::string &s)
{
    std::string r(s.rbegin(), s.rend());
    for (char &c : r) c = c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : c == 'T' ? 'A' : 'N';
    return r;
}

// Each SAM record's tab-separated fields (header lines skipped).
std::vector<std::vector<std::string>> records(const std::string &sam)
{
    std::vector<std::vector<std::string>> out;
    size_t i = 0;
    while (i < sam.size()) {
        size_t j = sam.find('\n', i);
        if (j == std::string::npos) j = sam.size();
        if (sam[i] != '@') {
            std::vector<std::string> f;
            size_t a = i;
            while (a <= j) {
                size_t b = sam.find('\t', a);
                if (b == std::string::npos || b > j) b = j;
                f.push_back(sam.substr(a, b - a));
                a = b + 1;
            }
            out.push_back(f);
        }
        i = j + 1;
    }
    return out;
}

std::string tags_of(const std::vector<std::string> &r)
{
    std::string t;
    for (size_t i = 11; i < r.size(); ++i) t += r[i] + " ";
    return t;
}

std::string row(const char *chrom, int pos1, char ref, const std::string &alt, const char *filter = "PASS")
{
    return std::string(chrom) + "\t" + std::to_string(pos1) + "\t.\t" + ref + "\t" + alt + "\t.\t" + filter + "\t.";
}

} // namespace

TEST_CASE("aln_known_snvs: sites load, non-SNVs and duplicates are counted, FILTER is ignored"
          * doctest::test_suite("unit/aln_known_snvs")) {
    KsIndex &x = fixture();
    const char r100 = x.c1[99], r200 = x.c1[199];
    const std::string vcf = write_vcf(x.dir + "/ok.vcf", {
        row("c1", 100, r100, std::string(1, other(r100))),
        row("c1", 200, r200, std::string(1, other(r200)) + "," + other(r200, other(r200)), "LowQual"),
        row("c1", 300, x.c1[299], std::string(1, x.c1[299]) + "T"),  // insertion: skipped
        row("c1", 100, r100, std::string(1, other(r100))),           // duplicate allele
    });
    const std::string out = x.dir + "/ok.sam";
    bwa_tests::StderrCapture capture;
    CHECK(run({"--known-snvs", vcf, "-f", out, x.prefix, x.fq}) == 0);
    CHECK(capture.text().find("3 SNV alleles from 4 records (1 alleles not SNVs, 1 duplicates skipped)")
          != std::string::npos);
}

TEST_CASE("aln_known_snvs: a contig, REF or position the reference does not have is an error naming the record"
          * doctest::test_suite("unit/aln_known_snvs")) {
    KsIndex &x = fixture();
    struct bad_t { std::string row, msg; };
    // c2 is revcomp(c1[100, 400)) + 50 Ns + 200 bases: positions 301-350 are N.
    const std::vector<bad_t> bads = {
        {row("chr9", 5, 'A', "C"), "chr9:5: contig 'chr9' is not in the reference"},
        {row("c1", 100, other(x.c1[99]), "T"), "c1:100: REF " + std::string(1, other(x.c1[99]))
             + " does not match the reference base " + x.c1[99]},
        {row("c2", 320, 'A', "C"), "c2:320: REF A, but the reference base is N"},
        {row("c1", 601, 'A', "C"), "c1:601 is past the end of the contig (600 bp)"},
    };
    for (size_t i = 0; i < bads.size(); ++i) {
        CAPTURE(bads[i].row);
        const std::string vcf = write_vcf(x.dir + "/bad" + std::to_string(i) + ".vcf", {bads[i].row});
        const std::string out = x.dir + "/bad.sam";
        bwa_tests::StderrCapture capture;
        CHECK(run({"--known-snvs", vcf, "-f", out, x.prefix, x.fq}) == 1);
        CHECK(capture.text().find(bads[i].msg) != std::string::npos);
    }
    bwa_tests::StderrCapture capture;
    CHECK(run({"--known-snvs", x.dir + "/missing.vcf", "-f", x.dir + "/m.sam", x.prefix, x.fq}) == 1);
    CHECK(capture.text().find("--known-snvs: cannot open") != std::string::npos);
}

TEST_CASE("aln_known_snvs: an empty site set leaves the output unchanged"
          * doctest::test_suite("unit/aln_known_snvs")) {
    KsIndex &x = fixture();
    const std::string base = x.dir + "/base.sam", with = x.dir + "/empty.sam";
    const std::string vcf = write_vcf(x.dir + "/empty.vcf", {row("c1", 300, x.c1[299], std::string(1, x.c1[299]) + "T")});
    bwa_tests::StderrCapture capture;
    REQUIRE(run({"-f", base, x.prefix, x.fq}) == 0);
    REQUIRE(run({"--known-snvs", vcf, "-f", with, x.prefix, x.fq}) == 0);
    CHECK(body(slurp(base)) == body(slurp(with)));
    REQUIRE(run({"-f", base + ".pe", x.prefix, x.fq, x.fq}) == 0);
    REQUIRE(run({"--known-snvs", vcf, "-f", with + ".pe", x.prefix, x.fq, x.fq}) == 0);
    CHECK(body(slurp(base + ".pe")) == body(slurp(with + ".pe")));
}

TEST_CASE("aln_known_snvs: an ALT read over the mismatch budget is placed; its REF twin's record is unchanged"
          * doctest::test_suite("unit/aln_known_snvs")) {
    KsIndex &x = fixture();
    // Site c1:500 (0-based 499), outside the part of c1 that c2 copies. A 20-base read has a
    // mismatch budget of 2 (aln -n 0.04): two errors plus the ALT is 3 against the reference, so
    // stock leaves it unmapped; crediting the ALT makes it 2, at MAPQ 25 (n_mm == max_diff).
    const char ref = x.c1[499], alt = other(ref);
    std::string ref_read = x.c1.substr(489, 20), alt_read = ref_read;
    alt_read[10] = alt;
    alt_read[3] = other(alt_read[3]);
    alt_read[16] = other(alt_read[16]);
    const std::string fq = x.dir + "/alt.fq";
    FILE *fp = fopen(fq.c_str(), "w");
    REQUIRE(fp != nullptr);
    fprintf(fp, "@ref\n%s\n+\n%s\n@alt\n%s\n+\n%s\n", ref_read.c_str(), std::string(20, 'I').c_str(),
            alt_read.c_str(), std::string(20, 'I').c_str());
    REQUIRE(fclose(fp) == 0);
    const std::string vcf = write_vcf(x.dir + "/site.vcf", {row("c1", 500, ref, std::string(1, alt))});
    const std::string stock = x.dir + "/stock.sam", ks = x.dir + "/ks.sam";
    bwa_tests::StderrCapture capture;
    REQUIRE(run({"-f", stock, x.prefix, fq}) == 0);
    REQUIRE(run({"--known-snvs", vcf, "-f", ks, x.prefix, fq}) == 0);
    const auto s = records(slurp(stock)), k = records(slurp(ks));
    REQUIRE(s.size() == 2);
    REQUIRE(k.size() == 2);
    // the REF read: the same record plus XV (REF; no differences away from the site; budget 2;
    // its ALT twin lands on the same locus)
    REQUIRE(k[0].size() == s[0].size() + 1);
    CHECK(std::vector<std::string>(k[0].begin(), k[0].end() - 1) == s[0]);
    const std::string site = std::string("c1:500:") + ref + ">" + alt;
    CHECK(k[0].back().rfind("XV:Z:" + site + ",R,0,2,1,", 0) == 0);
    CHECK(s[1][1] == "4");                     // stock: the ALT read is unmapped
    CHECK(k[1][1] == "0");                     // credited: forward, at the site
    CHECK(k[1][2] == "c1");
    CHECK(k[1][3] == "490");
    CHECK(k[1][4] == "25");
    CHECK(k[1][5] == "20M");
    const std::string tags = tags_of(k[1]);
    CHECK(tags.find("NM:i:3") != std::string::npos); // reported against the reference
    CHECK(tags.find("XM:i:2") != std::string::npos); // the search's mismatches, the ALT credited
    // ALT, 2 differences away from the site with a budget of 2: the budget rule leaves it out; its
    // REF twin (2 mismatches) lands on the same locus at MAPQ 25
    CHECK(tags.find("XV:Z:" + site + ",A,2,2,1,25") != std::string::npos);
}

TEST_CASE("aln_known_snvs: paired end, an ALT mate is placed by its own search or by mate rescue"
          * doctest::test_suite("unit/aln_known_snvs")) {
    KsIndex &x = fixture();
    // Site c1:500 (0-based 499). Pairs of 25-base reads with 55- to 65-base inserts from
    // c1[400, 600) (outside the part c2 copies) give the batch an insert size; two pairs with
    // 60-base inserts test the site:
    //  - "own": read 1 is the 20-base over-budget ALT read of the single-end test; its own search
    //    places it once the ALT is credited, and the pairing keeps it.
    //  - "resc": read 1 is c1[480, 512) with the ALT and three more mismatches: three differences
    //    with the ALT credited, over a 32-base read's budget of 2, so no search places it; mate
    //    rescue scores 32 - 4 * 3 = 20 against the window with the ALT (the minimum) and 16 against
    //    the reference, so only the ALT-aware rescue places it.
    const char ref = x.c1[499], alt = other(ref);
    std::string own = x.c1.substr(489, 20);
    own[10] = alt;
    own[3] = other(own[3]);
    own[16] = other(own[16]);
    std::string resc = x.c1.substr(480, 32);
    resc[19] = alt;
    for (int q : {4, 12, 27}) resc[q] = other(resc[q]);
    const std::string fq1 = x.dir + "/pe1.fq", fq2 = x.dir + "/pe2.fq";
    FILE *f1 = fopen(fq1.c_str(), "w"), *f2 = fopen(fq2.c_str(), "w");
    REQUIRE(f1 != nullptr);
    REQUIRE(f2 != nullptr);
    auto pair = [&](const std::string &name, const std::string &r1, int r2_end) {
        const std::string r2 = revcomp(x.c1.substr(r2_end - 25, 25));
        fprintf(f1, "@%s\n%s\n+\n%s\n", name.c_str(), r1.c_str(), std::string(r1.size(), 'I').c_str());
        fprintf(f2, "@%s\n%s\n+\n%s\n", name.c_str(), r2.c_str(), std::string(25, 'I').c_str());
    };
    for (int p = 400, i = 0; p + 65 <= 600; p += 5, ++i)
        pair("bg" + std::to_string(i), x.c1.substr(p, 25), p + 55 + i % 11);
    pair("own", own, 489 + 60);
    pair("resc", resc, 480 + 60);
    REQUIRE(fclose(f1) == 0);
    REQUIRE(fclose(f2) == 0);
    const std::string vcf = write_vcf(x.dir + "/pe.vcf", {row("c1", 500, ref, std::string(1, alt))});
    const std::string stock = x.dir + "/pe.stock.sam", ks = x.dir + "/pe.ks.sam";
    bwa_tests::StderrCapture capture;
    REQUIRE(run({"-f", stock, x.prefix, fq1, fq2}) == 0);
    REQUIRE(run({"--known-snvs", vcf, "-f", ks, x.prefix, fq1, fq2}) == 0);
    const auto s = records(slurp(stock)), k = records(slurp(ks));
    REQUIRE(s.size() == k.size());
    REQUIRE(s.size() >= 4);
    const size_t n = s.size();
    // the background pairs are unchanged; a read over the site (all REF) gains XV
    for (size_t i = 0; i + 4 < n; ++i) {
        std::vector<std::string> ki = k[i];
        if (ki.back().rfind("XV:Z:", 0) == 0) {
            CHECK(ki.back().rfind("XV:Z:c1:500:", 0) == 0);
            CHECK(ki.back().find(",R,0,") != std::string::npos);
            ki.pop_back();
        }
        CHECK(ki == s[i]);
    }
    const std::vector<std::string> &s_own = s[n - 4], &k_own = k[n - 4], &s_resc = s[n - 2], &k_resc = k[n - 2];
    REQUIRE(k_own[0] == "own");
    REQUIRE(k_resc[0] == "resc");
    const std::string site = std::string("c1:500:") + ref + ">" + alt;
    CHECK((std::stoi(s_own[1]) & 4) != 0);      // stock: unmapped
    CHECK((std::stoi(k_own[1]) & 0x12) == 0x2); // credited: forward, in a proper pair
    CHECK(k_own[3] == "490");
    CHECK(tags_of(k_own).find("XV:Z:" + site + ",A,2,2,") != std::string::npos);
    CHECK((std::stoi(s_resc[1]) & 4) != 0);
    CHECK((std::stoi(k_resc[1]) & 0x16) == 0x2); // rescued: mapped, forward, in a proper pair
    CHECK(k_resc[3] == "481");
    CHECK(k_resc[5] == "32M");
    CHECK(tags_of(k_resc).find("XV:Z:" + site + ",A,3,2,") != std::string::npos);
    // the mates are unchanged apart from their mate fields
    for (size_t r : {n - 3, n - 1}) {
        CHECK(k[r][2] == s[r][2]);
        CHECK(k[r][3] == s[r][3]);
        CHECK(k[r][5] == s[r][5]);
    }
}

TEST_CASE("aln_known_snvs: paired end under --meth, an over-budget ALT mate is placed and pairs"
          * doctest::test_suite("unit/aln_known_snvs")) {
    KsIndex &x = fixture();
    {
        std::vector<std::string> a = {"index", "--meth", x.prefix};
        std::vector<char *> argv;
        for (std::string &s : a) argv.push_back(&s[0]);
        argv.push_back(nullptr);
        bwa_tests::StderrCapture quiet;
        optind = 1;
#if defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || defined(__NetBSD__)
        optreset = 1;
#endif
        REQUIRE(bwa_index((int)a.size(), argv.data()) == 0);
    }
    // Unmethylated directional fragments from c1[400, 600): read 1 is the fragment's start with
    // every C read as T, read 2 the reverse complement of its end. The site is the first A or T in
    // c1[495, 505) with G as its ALT (no conversion explains either base on read 1's strand).
    int q = 495;
    while (q < 505 && x.c1[q] != 'A' && x.c1[q] != 'T') ++q;
    REQUIRE(q < 505);
    const char ref = x.c1[q], alt = 'G';
    auto conv = [](std::string s) { for (char &c : s) if (c == 'C') c = 'T'; return s; };
    const std::string fq1 = x.dir + "/m1.fq", fq2 = x.dir + "/m2.fq";
    FILE *f1 = fopen(fq1.c_str(), "w"), *f2 = fopen(fq2.c_str(), "w");
    REQUIRE(f1 != nullptr);
    REQUIRE(f2 != nullptr);
    auto pair = [&](const std::string &name, std::string frag) {
        const std::string r1 = conv(frag).substr(0, 25), r2 = revcomp(conv(frag).substr(frag.size() - 25));
        fprintf(f1, "@%s\n%s\n+\n%s\n", name.c_str(), r1.c_str(), std::string(25, 'I').c_str());
        fprintf(f2, "@%s\n%s\n+\n%s\n", name.c_str(), r2.c_str(), std::string(25, 'I').c_str());
    };
    for (int p = 400, i = 0; p + 65 <= 600; p += 5, ++i) pair("bg" + std::to_string(i), x.c1.substr(p, 55 + i % 11));
    // "own": the ALT plus two mismatches (bases no conversion explains) in read 1: three
    // differences against the reference, over a 25-base read's budget of 2; two with the ALT.
    const int beg = q - 10;
    std::string frag = x.c1.substr(beg, 60);
    frag[10] = alt;
    int n_mm = 0;
    for (int k = 0; k < 25 && n_mm < 2; ++k)
        if (k != 10 && k % 7 == 3 && (frag[k] == 'A' || frag[k] == 'G')) frag[k] = frag[k] == 'A' ? 'G' : 'A', ++n_mm;
    REQUIRE(n_mm == 2);
    pair("own", frag);
    REQUIRE(fclose(f1) == 0);
    REQUIRE(fclose(f2) == 0);
    const std::string vcf = write_vcf(x.dir + "/m.vcf", {row("c1", q + 1, ref, std::string(1, alt))});
    const std::string stock = x.dir + "/m.stock.sam", ks = x.dir + "/m.ks.sam";
    bwa_tests::StderrCapture capture;
    REQUIRE(run({"--meth", "-f", stock, x.prefix, fq1, fq2}) == 0);
    const int ks_rc = run({"--meth", "--known-snvs", vcf, "-f", ks, x.prefix, fq1, fq2});
    REQUIRE(ks_rc == 0);
    const auto s = records(slurp(stock)), k = records(slurp(ks));
    REQUIRE(s.size() == k.size());
    REQUIRE(s.size() >= 2);
    const size_t n = s.size();
    for (size_t i = 0; i + 2 < n; ++i) { // background: unchanged apart from XV on reads over the site
        std::vector<std::string> ki = k[i];
        if (ki.back().rfind("XV:Z:", 0) == 0) ki.pop_back();
        CHECK(ki == s[i]);
    }
    REQUIRE(k[n - 2][0] == "own");
    CHECK((std::stoi(s[n - 2][1]) & 4) != 0);       // stock: unmapped
    CHECK((std::stoi(k[n - 2][1]) & 0x16) == 0x2);  // credited: forward, in a proper pair
    CHECK(k[n - 2][3] == std::to_string(beg + 1));
    const std::string site = "c1:" + std::to_string(q + 1) + ":" + ref + ">" + alt;
    CHECK(tags_of(k[n - 2]).find("XV:Z:" + site + ",A,2,2,") != std::string::npos);
}

// ---- the window reference built in memory ----------------------------------------------------

#include "../../src/aln/aln_known_snvs.h"
#include "../../src/aln/aln_meth.h"
#include "../../src/meth_ref.h"

namespace {

struct KsWinFixture {
    bwa_tests::TempDir tmp_dir{"aln_ks_winref_unit"};
    std::string dir = tmp_dir.path(), ref, fa;
    bntseq_t *bns = nullptr;
    uint8_t *pac = nullptr;
    std::vector<bwa_aln::ks_window_t> win;

    KsWinFixture()
    {
        ref = dir + "/ref.fa";
        uint32_t st = 7;
        const std::string s1 = aln_test::lcg_dna(st, 400), s2 = aln_test::lcg_dna(st, 900);
        FILE *fp = fopen(ref.c_str(), "w");
        REQUIRE(fp != nullptr);
        fprintf(fp, ">chrA\n%s\n>chrB\n%s\n", s1.c_str(), s2.c_str());
        REQUIRE(fclose(fp) == 0);
        REQUIRE(bwa_idx_build(ref.c_str(), ref.c_str(), 0, 3) == 0);
        REQUIRE(meth_orig_ref_load(ref.c_str(), &bns, &pac, 1) == 0);
        const std::vector<std::vector<int64_t>> sites = {{5, 200, 395}, {2, 450, 897}};
        for (int rid = 0; rid < 2; ++rid)
            for (int64_t p : sites[rid]) {
                const int64_t beg = std::max<int64_t>(0, p - KS_WINDOW_HALF);
                const int64_t end = std::min<int64_t>(bns->anns[rid].len, p + KS_WINDOW_HALF + 1);
                const uint8_t refc = (pac[(bns->anns[rid].offset + p) >> 2] >> ((~(bns->anns[rid].offset + p) & 3) << 1)) & 3;
                win.push_back({rid, beg, (int32_t)(end - beg), p, (uint8_t)((refc + 1) & 3)});
            }
        fa = dir + "/windows.fa";
        fp = fopen(fa.c_str(), "w");
        REQUIRE(fp != nullptr);
        for (size_t i = 0; i < win.size(); ++i) {
            const auto &w = win[i];
            std::string seq((size_t)w.len, 'N');
            for (int32_t j = 0; j < w.len; ++j) {
                const int64_t k = bns->anns[w.rid].offset + w.beg + j;
                seq[(size_t)j] = "ACGT"[(pac[k >> 2] >> ((~k & 3) << 1)) & 3];
            }
            seq[(size_t)(w.site_pos - w.beg)] = "ACGT"[w.alt];
            fprintf(fp, ">w%zu\n%s\n", i, seq.c_str());
        }
        REQUIRE(fclose(fp) == 0);
    }
    ~KsWinFixture() { meth_orig_ref_free(&bns, &pac); }
};

void expect_bns_equal(const bntseq_t *a, const bntseq_t *b)
{
    REQUIRE(a->l_pac == b->l_pac);
    REQUIRE(a->n_seqs == b->n_seqs);
    CHECK(a->seed == b->seed);
    CHECK(a->n_holes == b->n_holes);
    CHECK((a->ambs == nullptr) == (b->ambs == nullptr));
    CHECK(a->fp_pac == b->fp_pac);
    CHECK((a->pos2rid_bucket == nullptr) == (b->pos2rid_bucket == nullptr));
    if (a->pos2rid_bucket && b->pos2rid_bucket) {
        const int64_t nb = (a->l_pac >> BNS_POS2RID_SHIFT) + 2;
        for (int64_t i = 0; i < nb; ++i) CHECK(a->pos2rid_bucket[i] == b->pos2rid_bucket[i]);
    }
    for (int i = 0; i < a->n_seqs; ++i) {
        CHECK(std::string(a->anns[i].name) == b->anns[i].name);
        CHECK(std::string(a->anns[i].anno) == b->anns[i].anno);
        CHECK(a->anns[i].offset == b->anns[i].offset);
        CHECK(a->anns[i].len == b->anns[i].len);
        CHECK(a->anns[i].n_ambs == b->anns[i].n_ambs);
        CHECK(a->anns[i].gi == b->anns[i].gi);
        CHECK(a->anns[i].is_alt == b->anns[i].is_alt);
    }
}

} // namespace

TEST_CASE("aln_known_snvs: window reference equals the indexed windows.fa")
{
    KsWinFixture f;
    REQUIRE(bwa_idx_build(f.fa.c_str(), f.fa.c_str(), 0, 3) == 0);
    bntseq_t *want_bns = nullptr;
    uint8_t *want_pac = nullptr;
    REQUIRE(meth_orig_ref_load(f.fa.c_str(), &want_bns, &want_pac, 1) == 0);
    bntseq_t *bns = nullptr;
    uint8_t *pac = nullptr;
    bwa_aln::ks_windows_ref(f.win, f.bns, f.pac, &bns, &pac);
    expect_bns_equal(want_bns, bns);
    CHECK(bns->n_holes == 0);
    CHECK(memcmp(want_pac, pac, (size_t)(want_bns->l_pac / 4 + 1)) == 0);
    bns_destroy(bns);
    free(pac);
    meth_orig_ref_free(&want_bns, &want_pac);
}

TEST_CASE("aln_known_snvs: --meth window conversion equals index --meth")
{
    KsWinFixture f;
    std::string a0 = "index", a1 = "--meth", a2 = f.fa;
    char *argv[] = {&a0[0], &a1[0], &a2[0], nullptr};
#if defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || defined(__NetBSD__)
    optind = 1;
    optreset = 1;
#else
    optind = 0;
#endif
    REQUIRE(bwa_index(3, argv) == 0);
    bntseq_t *want_bns = nullptr;
    uint8_t *want_pac = nullptr;
    REQUIRE(meth_orig_ref_load((f.fa + ".meth").c_str(), &want_bns, &want_pac, 1) == 0);
    bntseq_t *wbns = nullptr, *cbns = nullptr;
    uint8_t *wpac = nullptr, *cpac = nullptr;
    bwa_aln::ks_windows_ref(f.win, f.bns, f.pac, &wbns, &wpac);
    bwa_aln::ks_windows_meth_convert(wbns, wpac, &cbns, &cpac);
    expect_bns_equal(want_bns, cbns);
    CHECK(memcmp(want_pac, cpac, (size_t)(want_bns->l_pac / 4 + 1)) == 0);

    bwa_aln::aln_meth_ref_t ok;
    CHECK(bwa_aln::aln_meth_ref_adopt(&ok, wbns, wpac, cbns) == 0);
    bwa_aln::aln_meth_ref_free(&ok);

    bntseq_t *wbns2 = nullptr, *cbns2 = nullptr;
    uint8_t *wpac2 = nullptr, *cpac2 = nullptr;
    bwa_aln::ks_windows_ref(f.win, f.bns, f.pac, &wbns2, &wpac2);
    bwa_aln::ks_windows_meth_convert(wbns2, wpac2, &cbns2, &cpac2);
    cbns2->anns[3].name[1] = 'X';
    bwa_tests::StderrCapture cap;
    bwa_aln::aln_meth_ref_t bad;
    CHECK(bwa_aln::aln_meth_ref_adopt(&bad, wbns2, wpac2, cbns2) == -1);
    CHECK(bad.bns == nullptr);

    bns_destroy(cbns); free(cpac);
    bns_destroy(cbns2); free(cpac2);
    meth_orig_ref_free(&want_bns, &want_pac);
}
