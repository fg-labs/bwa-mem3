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
    // paired input with sites is refused until pairing uses them
    const std::string site = write_vcf(x.dir + "/one.vcf", {row("c1", 100, x.c1[99], std::string(1, other(x.c1[99])))});
    CHECK(run({"--known-snvs", site, "-f", with + ".pe1", x.prefix, x.fq, x.fq}) == 1);
    CHECK(capture.text().find("--known-snvs with paired-end input is not supported yet") != std::string::npos);
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
    auto records = [](const std::string &sam) {
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
    };
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
    const std::string tags = [&] { std::string t; for (size_t i = 11; i < k[1].size(); ++i) t += k[1][i] + " "; return t; }();
    CHECK(tags.find("NM:i:3") != std::string::npos); // reported against the reference
    CHECK(tags.find("XM:i:2") != std::string::npos); // the search's mismatches, the ALT credited
    // ALT, 2 differences away from the site with a budget of 2: the budget rule leaves it out; its
    // REF twin (2 mismatches) lands on the same locus at MAPQ 25
    CHECK(tags.find("XV:Z:" + site + ",A,2,2,1,25") != std::string::npos);
}
