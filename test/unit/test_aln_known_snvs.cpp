// `bwa-mem3 aln --known-snvs` (src/aln/aln_known_snvs.cpp): loading the sites, the errors that
// name a bad record, and output that does not change when no read is affected.

// HTSlib first: its kstring API must precede the aligner's same-guard header.
#include "htslib/vcf.h"
#include <zlib.h>
#include "doctest/doctest.h"

#include "aln_test_ref.h"
#include "../framework/run_in_child.h"
#include "../framework/stderr_capture.h"
#include "../../src/aln/aln_main.h"
#include "../../src/aln/aln_bwt.h"
#include "../../src/aln/aln_known_snvs.h"
#include "../../src/aln/aln_meth_realign.h"
#include "../../src/aln/bwape.h"
#include "../../src/aln/bwaseqio.h"
#include "../../src/aln/bwtgap.h"
#include "../../src/FMI_search.h"
#include "../../src/bntseq.h"
#include "../../src/aln/bwase.h"
#include "../../src/bwa.h"

#include <algorithm>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <unistd.h>
#include <dirent.h>
#include <getopt.h>
#include <map>
#include <memory>
#include <set>
#include <sys/stat.h>

#if defined(BWA_TESTS_CHILD_SKIP_LEAK_CHECK)
#include <sanitizer/lsan_interface.h>
#endif

int bwa_index(int argc, char *argv[]); // src/bwtindex.cpp, the `index` command

namespace {

// The index the cases share. Its files stay until the process exits rather than being dropped
// (TempIndex::drop_files) once loaded: most cases run the aln command, which loads the index from
// disk itself. A case's own scratch files (VCF, reads, output) go in its own TempDir, and a case
// that needs the converted index builds its own (KsMethIndex).
struct KsIndex {
    bwa_tests::TempDir tmp_dir{"aln_known_snvs_unit"};
    std::string dir = tmp_dir.path(), prefix, fq;
    std::string c1, c2;

    KsIndex()
    {
        prefix = tmp_dir.unique_file();
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

// Exact name/mate multiset, independent of output ordering and of a second execution.
void require_records(const std::string &sam, const std::vector<std::string> &names, bool paired)
{
    std::map<std::pair<std::string, int>, int> expected, actual;
    for (const auto &name : names) {
        ++expected[{name, paired ? 64 : 0}];
        if (paired) ++expected[{name, 128}];
    }
    const auto rec = records(sam);
    REQUIRE(rec.size() == names.size() * (paired ? 2 : 1));
    for (const auto &r : rec) {
        REQUIRE(r.size() >= 11);
        const int flag = std::stoi(r[1]);
        const int mate = paired ? flag & (64 | 128) : 0;
        ++actual[{r[0], mate}];
    }
    CHECK(actual == expected);
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

// Builds `index --meth` for the reference at prefix.
void build_meth_index(const std::string &prefix)
{
    std::vector<std::string> a = {"index", "--meth", prefix};
    std::vector<char *> argv;
    for (std::string &s : a) argv.push_back(&s[0]);
    argv.push_back(nullptr);
    bwa_tests::StderrCapture quiet;
    // BSD getopt resets with optreset; glibc needs optind = 0 for a full reset.
#if defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || defined(__NetBSD__)
    optind = 1;
    optreset = 1;
#else
    optind = 0;
#endif
    REQUIRE(bwa_index((int)a.size(), argv.data()) == 0);
}

// A case's own copy of the fixture with its `index --meth` converted index too, removed with the
// case: the shared fixture never holds converted index files.
struct KsMethIndex : KsIndex {
    KsMethIndex() { build_meth_index(prefix); }
};

std::string write_fq(const std::string &path, const std::vector<std::pair<std::string, std::string>> &reads)
{
    FILE *fp = fopen(path.c_str(), "w");
    REQUIRE(fp != nullptr);
    for (const auto &r : reads)
        fprintf(fp, "@%s\n%s\n+\n%s\n", r.first.c_str(), r.second.c_str(), std::string(r.second.size(), 'I').c_str());
    REQUIRE(fclose(fp) == 0);
    return path;
}

// MD:Z and the mismatch count of an ungapped forward read against the reference segment under it.
std::string md_of(const std::string &ref, const std::string &read, int *n_mm)
{
    std::string md;
    int run = 0;
    *n_mm = 0;
    for (size_t i = 0; i < read.size(); ++i) {
        if (read[i] == ref[i]) { ++run; continue; }
        md += std::to_string(run) + ref[i];
        run = 0;
        ++*n_mm;
    }
    return md + std::to_string(run);
}

// A forward, ungapped read's expected record at 1-based pos on c1, from the reference itself: SEQ,
// QUAL, CIGAR, NM and MD. The caller fills flag, MAPQ, mate fields and the remaining tags.
struct expect_t {
    std::vector<std::string> f;
    int nm = 0;
};
expect_t expect_fwd(const std::string &name, const std::string &c1, int pos, const std::string &read)
{
    expect_t e;
    const std::string md = md_of(c1.substr(pos - 1, read.size()), read, &e.nm);
    e.f = {name, "", "c1", std::to_string(pos), "", std::to_string(read.size()) + "M", "", "", "", read,
           std::string(read.size(), 'I'), "NM:i:" + std::to_string(e.nm), "MD:Z:" + md};
    return e;
}

// The tag with this two-letter name and type (e.g. "XV:Z:"), or "" when absent.
std::string tag(const std::vector<std::string> &r, const std::string &prefix)
{
    for (size_t i = 11; i < r.size(); ++i)
        if (r[i].rfind(prefix, 0) == 0) return r[i];
    return "";
}

} // namespace

TEST_CASE("aln_known_snvs: merge guards cover large options and distinct combined loci"
          * doctest::test_suite("unit/aln_known_snvs")) {
    namespace ba = bwa_aln;
    KsIndex &x = fixture();
    FMI_search fmi(x.prefix.c_str()); fmi.load_index(true, 1, false);
    ba::bwt_t *bwt = ba::bwt_from_fmi(&fmi);
    // Synthetic window hits use the main index's coordinate system; merging only reads rows.
    ba::ks_windows_t w; w.bns = fmi.idx->bns; w.bwt = bwt;
    w.win.resize(w.bns->n_seqs);
    w.win[0].rid = 0; w.win[0].beg = 0; w.win[0].len = x.c1.size();
    w.win[0].site_pos = 499; w.win[0].alt = 0;
    std::vector<ba::bwtint_t> rows(x.c1.size(), (ba::bwtint_t)-1);
    for (ba::bwtint_t r = 0; r <= bwt->seq_len; ++r) {
        int strand;
        const auto p = ba::bwa_sa2pos(w.bns, bwt, r, 120, &strand);
        if (strand == 0 && p < rows.size()) rows[p] = r;
    }
    auto hit = [&](int pos) {
        REQUIRE(rows[pos] != (ba::bwtint_t)-1);
        ba::bwt_aln1_t h = {}; h.k = h.l = rows[pos]; return h;
    };
    ba::gap_opt_t *opt = ba::gap_init_opt(); opt->fnr = -1; opt->max_diff = 2;
    std::vector<ba::bwt_aln1_t> main, win{hit(400)};
    ba::ks_outcome_t expected = ba::KS_CREDITED;
    SUBCASE("large gap extension bound falls back") {
        opt->max_gape = INT_MAX; expected = ba::KS_INCOMPLETE;
    }
    SUBCASE("large mismatch penalty retains the best gapped hit") {
        opt->s_mm = INT_MAX; win[0].score = 11; win[0].n_gapo = 1;
    }
    SUBCASE("exactly 256 combined loci are allowed") {
        for (int p = 0; p < 156; ++p) main.push_back(hit(p));
        win.clear(); for (int p = 380; p < 480; ++p) win.push_back(hit(p));
    }
    SUBCASE("257 combined loci fall back") {
        for (int p = 0; p < 157; ++p) main.push_back(hit(p));
        win.clear(); for (int p = 380; p < 480; ++p) win.push_back(hit(p));
        expected = ba::KS_INCOMPLETE;
    }
    SUBCASE("duplicate main rows are not distinct loci") {
        for (int p = 0; p < 200; ++p) main.push_back(hit(p / 2));
        win.clear(); for (int p = 380; p < 480; ++p) win.push_back(hit(p));
    }
    SUBCASE("unchanged window hits preserve the stock result") {
        for (int p = 340; p < 596; ++p) main.push_back(hit(p));
        expected = ba::KS_SAME;
    }
    uint8_t seq[120] = {};
    ba::bwa_seq_t s = {}; s.len = 120; s.seq = seq;
    s.aln = main.data(); s.n_aln = main.size(); s.pos = 7;
    int64_t alt = -1;
    // A search that stopped at -m may lack competitors: a read that would be credited keeps its
    // stock result (and no listed hits); the other outcomes are unchanged.
    std::vector<ba::pe_hit_t> stopped_hits(1);
    CHECK(ba::ks_override(&w, w.bns, bwt, opt, 0, &s, win.data(), win.size(), &alt, &stopped_hits, true)
          == (expected == ba::KS_CREDITED ? ba::KS_INCOMPLETE : expected));
    CHECK(s.pos == 7); CHECK(alt == -1); CHECK(stopped_hits.empty());
    CHECK(ba::ks_override(&w, w.bns, bwt, opt, 0, &s, win.data(), win.size(), &alt) == expected);
    if (expected != ba::KS_CREDITED) { CHECK(s.pos == 7); CHECK(alt == -1); }
    ba::bwt_destroy(bwt); free(opt);
}

TEST_CASE("aln_known_snvs: merged histories count each locus once and restore -R"
          * doctest::test_suite("unit/aln_known_snvs")) {
    namespace ba = bwa_aln;
    KsIndex &x = fixture();
    FMI_search fmi(x.prefix.c_str()); fmi.load_index(true, 1, false);
    ba::bwt_t *bwt = ba::bwt_from_fmi(&fmi);
    ba::ks_windows_t w; w.bns = fmi.idx->bns; w.bwt = bwt;
    w.win.resize(w.bns->n_seqs);
    w.win[0].rid = 0; w.win[0].len = x.c1.size(); w.win[0].site_pos = 410; w.win[0].alt = 0;
    auto hit = [&](int pos, int score) {
        ba::bwt_aln1_t h = {}; h.score = score; h.n_mm = score / 3;
        bool found = false;
        for (ba::bwtint_t r = 0; r <= bwt->seq_len; ++r) {
            int strand;
            if (ba::bwa_sa2pos(w.bns, bwt, r, 32, &strand) == (ba::bwtint_t)pos && strand == 0) {
                h.k = h.l = r; found = true; break;
            }
        }
        REQUIRE(found); return h;
    };
    ba::gap_opt_t *opt = ba::gap_init_opt(); opt->fnr = -1; opt->max_diff = 3;
    std::vector<ba::bwt_aln1_t> main{hit(400, 3), hit(400, 3), hit(100, 3)}, win{hit(400, 0)};
    int want_c2 = 1;
    SUBCASE("duplicate main histories are collapsed") {}
    SUBCASE("-R zero excludes the worse stratum") { opt->max_top2 = 0; want_c2 = 0; }
    SUBCASE("-R equality retains the worse stratum") { opt->max_top2 = 1; }
    SUBCASE("exhaustive mode ignores -R") { opt->mode |= BWA_MODE_NONSTOP; opt->max_top2 = 0; }
    ba::bwa_seq_t s = {}; s.len = 32; s.aln = main.data(); s.n_aln = main.size();
    uint8_t seq[32] = {}; s.seq = seq;
    int64_t alt;
    std::vector<ba::pe_hit_t> hits;
    REQUIRE(ba::ks_override(&w, w.bns, bwt, opt, 10, &s, win.data(), win.size(), &alt, &hits) == ba::KS_CREDITED);
    CHECK((int)s.c1 == 1); CHECK((int)s.c2 == want_c2);
    CHECK(hits.size() == (size_t)(1 + want_c2));
    CHECK(s.n_multi == want_c2);
    free(s.multi); free(opt); ba::bwt_destroy(bwt);
}

TEST_CASE("aln_known_snvs: contig-crossing refined alignments have no votes"
          * doctest::test_suite("unit/aln_known_snvs")) {
    namespace ba = bwa_aln;
    bntann1_t anns[2] = {}; anns[0].len = 20; anns[1].offset = 20; anns[1].len = 20;
    bntseq_t bns = {}; bns.n_seqs = 2; bns.anns = anns; bns.l_pac = 40;
    ba::ks_sites_t sites; sites.by_contig.resize(2); sites.by_contig[0] = {{15, 0, 1}};
    ba::gap_opt_t *opt = ba::gap_init_opt();
    uint8_t seq[25] = {}; ba::bwa_seq_t s = {};
    s.seq = seq; s.pos = 8; s.len = 25; s.type = BWA_TYPE_UNIQUE;
    ba::bwa_cigar_t cigar[4] = {__cigar_create(FROM_M, 12), __cigar_create(FROM_S, 13)};
    bool crosses = true;
    SUBCASE("ungapped crosses") {}
    SUBCASE("deletion crosses") {
        // Without the deletion the span (5 + 7) ends exactly at the contig end; with it, one past.
        cigar[0] = __cigar_create(FROM_M, 5); cigar[1] = __cigar_create(FROM_D, 1);
        cigar[2] = __cigar_create(FROM_M, 7); cigar[3] = __cigar_create(FROM_S, 13);
        s.cigar = cigar; s.n_cigar = 4;
    }
    SUBCASE("soft clip ends at the boundary") { s.cigar = cigar; s.n_cigar = 2; crosses = false; }
    std::vector<ba::ks_vote_t> votes;
    CHECK(ba::ks_site_votes(&sites, &bns, opt, &s, -1, -1, &votes) == (crosses ? -1 : 0));
    CHECK(votes.size() == (crosses ? 0 : 1));
    free(opt);
}

TEST_CASE("aln_known_snvs: XA uses selected-allele CIGAR and reference edit distance"
          * doctest::test_suite("unit/aln_known_snvs")) {
    bwa_tests::TempDir tmp{"ks_xa"}; const std::string dir = tmp.path();
    const std::string ref = "GCTAGTCAGTATGCATCCCATCTCGGCTCTCTGACTCCGG";
    std::string alt_ref = ref; alt_ref[25] = 'A';
    const std::string prefix = tmp.unique_file();
    REQUIRE(aln_test::write_test_fasta(prefix, {{"c1", alt_ref}, {"c2", ref}}));
    REQUIRE(bwa_idx_build(prefix.c_str(), prefix.c_str(), 0, 3) == 0);
    const std::string vcf = write_vcf(dir + "/s.vcf", {row("c2", 26, 'G', "A")});
    for (bool gapped : {false, true}) for (bool reverse : {false, true})
    for (bool paired : {false, true}) for (bool trimmed : {false, true}) {
        CAPTURE(gapped); CAPTURE(reverse); CAPTURE(paired); CAPTURE(trimmed);
        const std::string read = gapped ? "GCTAGTCAGTATGCATCCCATCTCGATCTCTGACTCCGG" : alt_ref;
        const std::string fq = dir + "/r.fq";
        const std::string bases = (reverse ? revcomp(read) : read) + (trimmed ? "NNNNN" : "");
        const std::string qual = std::string(read.size(), 'I') + (trimmed ? "!!!!!" : "");
        FILE *fp = fopen(fq.c_str(), "w"); REQUIRE(fp != nullptr);
        fprintf(fp, "@r\n%s\n+\n%s\n", bases.c_str(), qual.c_str()); REQUIRE(fclose(fp) == 0);
        const std::string mate = write_fq(dir + "/m.fq", {{"r", std::string(30, 'N')}});
        const std::string sam = dir + "/out.sam";
        std::vector<std::string> args = {"--known-snvs", vcf, "-q", "20", "-n", "2", "-f", sam, prefix, fq};
        if (paired) args.push_back(mate);
        bwa_tests::StderrCapture capture;
        REQUIRE(run(args) == 0);
        const auto rec = records(slurp(sam)); REQUIRE(rec.size() == (paired ? 2 : 1));
        CHECK(rec[0][2] == "c1");
        const std::string cigar = std::string(trimmed && reverse ? "5S" : "")
            + (gapped ? "26M1D13M" : "40M") + (trimmed && !reverse ? "5S" : "");
        const std::string expected = std::string("XA:Z:c2,") + (reverse ? "-" : "+") + "1,"
            + cigar + (gapped ? ",2;" : ",1;");
        CHECK(tag(rec[0], "XA:Z:") == expected);
    }
}

TEST_CASE("aln_known_snvs: sites load, non-SNVs and duplicates are counted, FILTER is ignored"
          * doctest::test_suite("unit/aln_known_snvs")) {
    KsIndex &x = fixture();
    bwa_tests::TempDir scratch{"ks_case"}; // this case's VCF, reads and output, removed when it ends
    const std::string dir = scratch.path();
    const char r100 = x.c1[99], r200 = x.c1[199];
    const std::string vcf = write_vcf(dir + "/ok.vcf", {
        row("c1", 100, r100, std::string(1, other(r100))),
        row("c1", 200, r200, std::string(1, other(r200)) + "," + other(r200, other(r200)), "LowQual"),
        row("c1", 300, x.c1[299], std::string(1, x.c1[299]) + "T"),  // insertion: skipped
        row("c1", 100, r100, std::string(1, other(r100))),           // duplicate allele
    });
    const std::string out = dir + "/ok.sam";
    bwa_tests::StderrCapture capture;
    CHECK(run({"--known-snvs", vcf, "-f", out, x.prefix, x.fq}) == 0);
    CHECK(capture.text().find("3 SNV alleles from 4 records (1 alleles not SNVs, 1 duplicates skipped)")
          != std::string::npos);
}

TEST_CASE("aln_known_snvs: sampled VCF compressed VCF and BCF load identical sites and reject truncation"
          * doctest::test_suite("unit/aln_known_snvs")) {
    namespace ba = bwa_aln;
    KsIndex &x = fixture();
    bwa_tests::TempDir tmp{"ks_formats"}; const std::string dir = tmp.path();
    const std::string text = dir + "/s.vcf";
    FILE *fp = fopen(text.c_str(), "w"); REQUIRE(fp != nullptr);
    fprintf(fp, "##fileformat=VCFv4.2\n##contig=<ID=c1>\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n");
    for (int p = 200; p < 500; ++p)
        fprintf(fp, "c1\t%d\t.\t%c\t%c\t.\tPASS\t.\tGT\t0/1\n", p + 1, x.c1[p], other(x.c1[p]));
    REQUIRE(fclose(fp) == 0);
    // Flush the header into its own BGZF block so truncation is necessarily in record data.
    for (bool binary : {false, true}) {
        htsFile *in = hts_open(text.c_str(), "r"); REQUIRE(in != nullptr);
        bcf_hdr_t *hdr = bcf_hdr_read(in); REQUIRE(hdr != nullptr);
        const std::string path = dir + (binary ? "/s.bcf" : "/s.vcf.gz");
        htsFile *out = hts_open(path.c_str(), binary ? "wb" : "wz"); REQUIRE(out != nullptr);
        REQUIRE(bcf_hdr_write(out, hdr) == 0); REQUIRE(hts_flush(out) == 0);
        bcf1_t *r = bcf_init(); REQUIRE(r != nullptr);
        int rc;
        while ((rc = bcf_read(in, hdr, r)) == 0) REQUIRE(bcf_write(out, hdr, r) == 0);
        REQUIRE(rc == -1);
        bcf_destroy(r); bcf_hdr_destroy(hdr);
        REQUIRE(hts_close(in) == 0); REQUIRE(hts_close(out) == 0);
        std::string bytes = slurp(path);
        REQUIRE(bytes.size() > 100);
        // A partial last data block, not just a missing optional EOF marker.
        bytes.resize(bytes.size() - 40);
        const std::string truncated = path + ".truncated";
        fp = fopen(truncated.c_str(), "wb"); REQUIRE(fp != nullptr);
        REQUIRE(fwrite(bytes.data(), 1, bytes.size(), fp) == bytes.size()); REQUIRE(fclose(fp) == 0);
    }
    // Ordinary gzip is a distinct supported transport from BGZF-compressed VCF.
    const std::string gzip = dir + "/ordinary.vcf.gz", contents = slurp(text);
    gzFile gz = gzopen(gzip.c_str(), "wb"); REQUIRE(gz != nullptr);
    REQUIRE(gzwrite(gz, contents.data(), contents.size()) == (int)contents.size()); REQUIRE(gzclose(gz) == Z_OK);
    std::string gzip_bytes = slurp(gzip); REQUIRE(gzip_bytes.size() > 40); gzip_bytes.resize(gzip_bytes.size() - 20);
    fp = fopen((gzip + ".truncated").c_str(), "wb"); REQUIRE(fp != nullptr);
    REQUIRE(fwrite(gzip_bytes.data(), 1, gzip_bytes.size(), fp) == gzip_bytes.size()); REQUIRE(fclose(fp) == 0);
    FMI_search fmi(x.prefix.c_str()); fmi.load_index(true, 1, false);
    ba::ks_sites_t baseline;
    REQUIRE(ba::ks_load(text.c_str(), fmi.idx->bns, fmi.idx->pac, &baseline) == 0);
    REQUIRE(baseline.n_sites == 300);
    for (const std::string suffix : {"/s.vcf.gz", "/s.bcf", "/ordinary.vcf.gz"}) {
        ba::ks_sites_t sites;
        REQUIRE(ba::ks_load((dir + suffix).c_str(), fmi.idx->bns, fmi.idx->pac, &sites) == 0);
        CHECK(sites.n_records == baseline.n_records); CHECK(sites.n_sites == baseline.n_sites);
        REQUIRE(sites.by_contig[0].size() == baseline.by_contig[0].size());
        for (size_t i = 0; i < baseline.by_contig[0].size(); ++i) {
            CHECK(sites.by_contig[0][i].pos == baseline.by_contig[0][i].pos);
            CHECK(sites.by_contig[0][i].ref == baseline.by_contig[0][i].ref);
            CHECK(sites.by_contig[0][i].alt == baseline.by_contig[0][i].alt);
        }
        bwa_tests::StderrCapture quiet;
        CHECK(ba::ks_load((dir + suffix + ".truncated").c_str(), fmi.idx->bns, fmi.idx->pac, &sites) == -1);
        const std::string error = quiet.text();
        if (suffix == "/ordinary.vcf.gz")
            CHECK(error.find("not a VCF or BCF file") != std::string::npos);
        else CHECK(error.find("read error after record") != std::string::npos);
    }
}

TEST_CASE("aln_known_snvs: a contig, REF or position the reference does not have is an error naming the record"
          * doctest::test_suite("unit/aln_known_snvs")) {
    KsIndex &x = fixture();
    bwa_tests::TempDir scratch{"ks_case"}; // this case's VCF, reads and output, removed when it ends
    const std::string dir = scratch.path();
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
        const std::string vcf = write_vcf(dir + "/bad" + std::to_string(i) + ".vcf", {bads[i].row});
        const std::string out = dir + "/bad.sam";
        bwa_tests::StderrCapture capture;
        CHECK(run({"--known-snvs", vcf, "-f", out, x.prefix, x.fq}) == 1);
        CHECK(capture.text().find(bads[i].msg) != std::string::npos);
    }
    bwa_tests::StderrCapture capture;
    CHECK(run({"--known-snvs", dir + "/missing.vcf", "-f", dir + "/m.sam", x.prefix, x.fq}) == 1);
    CHECK(capture.text().find("--known-snvs: cannot open") != std::string::npos);
}

TEST_CASE("aln_known_snvs: an empty site set leaves the output unchanged"
          * doctest::test_suite("unit/aln_known_snvs")) {
    KsIndex &x = fixture();
    bwa_tests::TempDir scratch{"ks_case"}; // this case's VCF, reads and output, removed when it ends
    const std::string dir = scratch.path();
    const std::string base = dir + "/base.sam", with = dir + "/empty.sam";
    const std::string vcf = write_vcf(dir + "/empty.vcf", {row("c1", 300, x.c1[299], std::string(1, x.c1[299]) + "T")});
    bwa_tests::StderrCapture capture;
    REQUIRE(run({"-f", base, x.prefix, x.fq}) == 0);
    REQUIRE(run({"--known-snvs", vcf, "-f", with, x.prefix, x.fq}) == 0);
    const std::vector<std::string> names{"r0", "r1", "r2", "r3", "r4", "r5"};
    require_records(slurp(base), names, false); require_records(slurp(with), names, false);
    CHECK(body(slurp(base)) == body(slurp(with)));
    REQUIRE(run({"-f", base + ".pe", x.prefix, x.fq, x.fq}) == 0);
    REQUIRE(run({"--known-snvs", vcf, "-f", with + ".pe", x.prefix, x.fq, x.fq}) == 0);
    require_records(slurp(base + ".pe"), names, true); require_records(slurp(with + ".pe"), names, true);
    CHECK(body(slurp(base + ".pe")) == body(slurp(with + ".pe")));
}

TEST_CASE("aln_known_snvs: an ALT read over the mismatch budget is placed; its REF twin's record is unchanged"
          * doctest::test_suite("unit/aln_known_snvs")) {
    KsIndex &x = fixture();
    bwa_tests::TempDir scratch{"ks_case"}; // this case's VCF, reads and output, removed when it ends
    const std::string dir = scratch.path();
    // Site c1:500 (0-based 499), outside the part of c1 that c2 copies. A 20-base read has a
    // mismatch budget of 2 (aln -n 0.04): two errors plus the ALT is 3 against the reference, so
    // stock leaves it unmapped; crediting the ALT makes it 2, at MAPQ 25 (n_mm == max_diff).
    const char ref = x.c1[499], alt = other(ref);
    std::string ref_read = x.c1.substr(489, 20), alt_read = ref_read;
    alt_read[10] = alt;
    alt_read[3] = other(alt_read[3]);
    alt_read[16] = other(alt_read[16]);
    const std::string fq = dir + "/alt.fq";
    FILE *fp = fopen(fq.c_str(), "w");
    REQUIRE(fp != nullptr);
    fprintf(fp, "@ref\n%s\n+\n%s\n@alt\n%s\n+\n%s\n", ref_read.c_str(), std::string(20, 'I').c_str(),
            alt_read.c_str(), std::string(20, 'I').c_str());
    REQUIRE(fclose(fp) == 0);
    const std::string vcf = write_vcf(dir + "/site.vcf", {row("c1", 500, ref, std::string(1, alt))});
    const std::string stock = dir + "/stock.sam", ks = dir + "/ks.sam";
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
    CHECK(k[0].back() == "XV:Z:" + site + ",R,0,2,1,37");
    CHECK(s[1][1] == "4");                     // stock: the ALT read is unmapped
    // Credited: forward at the read's true locus; NM and MD against the reference (three
    // mismatches, the ALT among them); XM the search's two, with the ALT credited; MAPQ 25
    // (n_mm == max_diff). XV: ALT, 2 differences away from the site with a budget of 2 (the budget
    // rule leaves it out); its REF twin (2 mismatches) lands on the same locus at MAPQ 25.
    const expect_t e = expect_fwd("alt", x.c1, 490, alt_read);
    REQUIRE(e.nm == 3);
    const std::vector<std::string> want = {"alt", "0", "c1", "490", "25", "20M", "*", "0", "0", alt_read,
        std::string(20, 'I'), "XT:A:U", e.f[11], "X0:i:1", "X1:i:0", "XM:i:" + std::to_string(e.nm - 1),
        "XO:i:0", "XG:i:0", e.f[12], "XV:Z:" + site + ",A," + std::to_string(e.nm - 1) + ",2,1,25"};
    CHECK(k[1] == want);
}

TEST_CASE("aln_known_snvs: a read whose search stops at -m keeps its stock record"
          * doctest::test_suite("unit/aln_known_snvs")) {
    // A 20-base read over a site with the ALT and two errors: three differences against the
    // reference, over its budget of 2, so stock leaves it unmapped; crediting places it. On a
    // 20 kb reference, -m 16 stops its genomic search but not its ALT-window search: the window
    // hit cannot be credited against a genomic search that may have missed competitors.
    bwa_tests::TempDir tmp{"ks_stopped"}; const std::string dir = tmp.path();
    uint32_t state = 911; const std::string c1 = aln_test::lcg_dna(state, 20000);
    const std::string prefix = tmp.unique_file();
    REQUIRE(aln_test::write_test_fasta(prefix, {{"c1", c1}}));
    REQUIRE(bwa_idx_build(prefix.c_str(), prefix.c_str(), 0, 3) == 0);
    const int site = 10000;
    const char ref = c1[site], alt = other(ref);
    std::string alt_read = c1.substr(site - 10, 20);
    alt_read[10] = alt;
    alt_read[3] = other(alt_read[3]);
    alt_read[16] = other(alt_read[16]);
    const std::string mate = revcomp(c1.substr(site + 30, 20));
    const std::string fq = dir + "/alt.fq", fq2 = dir + "/mate.fq";
    for (const auto &f : {std::make_pair(fq, alt_read), std::make_pair(fq2, mate)}) {
        FILE *fp = fopen(f.first.c_str(), "w");
        REQUIRE(fp != nullptr);
        fprintf(fp, "@alt\n%s\n+\n%s\n", f.second.c_str(), std::string(20, 'I').c_str());
        REQUIRE(fclose(fp) == 0);
    }
    const std::string vcf = write_vcf(dir + "/site.vcf", {row("c1", site + 1, ref, std::string(1, alt))});
    const std::string stock = dir + "/stock.sam", ks = dir + "/ks.sam";
    for (bool paired : {false, true}) {
        CAPTURE(paired);
        std::vector<std::string> in = {prefix, fq};
        if (paired) in.push_back(fq2);
        auto aln = [&](std::vector<std::string> a) {
            a.insert(a.end(), in.begin(), in.end());
            bwa_tests::StderrCapture capture;
            REQUIRE(run(a) == 0);
            return capture.text();
        };
        const std::string what = paired ? "ends" : "reads";
        // control: the default -m completes both searches and the read is credited
        std::string err = aln({"-v", "3", "-f", ks, "--known-snvs", vcf});
        CHECK(err.find(" 1 " + what + " credited, 0 at a site with no change, 0 kept stock") != std::string::npos);
        CHECK((std::stoi(records(slurp(ks))[0][1]) & 4) == 0);
        aln({"-m", "16", "-f", stock});
        err = aln({"-m", "16", "-v", "3", "-f", ks, "--known-snvs", vcf});
        INFO(err);
        CHECK(err.find(" 0 " + what + " credited, 0 at a site with no change, 1 kept stock") != std::string::npos);
        const auto s = records(slurp(stock)), k = records(slurp(ks));
        REQUIRE(s.size() == (paired ? 2u : 1u));
        REQUIRE(k.size() == s.size());
        CHECK((std::stoi(s[0][1]) & 4) != 0); // stock: unmapped
        if (!paired) CHECK(k == s);           // paired: ALT-aware mate rescue may still place it
    }
}

TEST_CASE("aln_known_snvs: paired end, an ALT mate is placed by its own search or by mate rescue"
          * doctest::test_suite("unit/aln_known_snvs")) {
    KsIndex &x = fixture();
    bwa_tests::TempDir scratch{"ks_case"}; // this case's VCF, reads and output, removed when it ends
    const std::string dir = scratch.path();
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
    const std::string fq1 = dir + "/pe1.fq", fq2 = dir + "/pe2.fq";
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
    const std::string vcf = write_vcf(dir + "/pe.vcf", {row("c1", 500, ref, std::string(1, alt))});
    const std::string stock = dir + "/pe.stock.sam", ks = dir + "/pe.ks.sam";
    bwa_tests::StderrCapture capture;
    REQUIRE(run({"-f", stock, x.prefix, fq1, fq2}) == 0);
    REQUIRE(run({"--known-snvs", vcf, "-f", ks, x.prefix, fq1, fq2}) == 0);
    const auto s = records(slurp(stock)), k = records(slurp(ks));
    REQUIRE(s.size() == k.size());
    std::vector<std::string> names;
    for (int i = 0; i < 28; ++i) names.push_back("bg" + std::to_string(i));
    names.push_back("own"); names.push_back("resc");
    require_records(slurp(stock), names, true); require_records(slurp(ks), names, true);
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
    const std::vector<std::string> &s_own = s[n - 4], &s_resc = s[n - 2];
    REQUIRE(k[n - 4][0] == "own");
    REQUIRE(k[n - 2][0] == "resc");
    const std::string site = std::string("c1:500:") + ref + ">" + alt;
    CHECK((std::stoi(s_own[1]) & 4) != 0); // stock: unmapped
    CHECK((std::stoi(s_resc[1]) & 4) != 0);
    // Credited and paired, both ends forward/reverse in a proper pair (99/147), 60-base inserts.
    // NM and MD come from the reference; XM is the search's count with the ALT credited; XV's
    // `away` leaves out the site. MAPQ, SM and AM are the pairing's.
    const expect_t eo = expect_fwd("own", x.c1, 490, own), er = expect_fwd("resc", x.c1, 481, resc);
    REQUIRE(eo.nm == 3);
    REQUIRE(er.nm == 4);
    const std::string q20(20, 'I'), q25(25, 'I'), q32(32, 'I');
    const std::vector<std::string> want_own = {"own", "99", "c1", "490", "60", "20M", "=", "525", "60", own, q20,
        "XT:A:U", eo.f[11], "SM:i:25", "AM:i:25", "X0:i:1", "X1:i:0", "XM:i:2", "XO:i:0", "XG:i:0", eo.f[12],
        "XV:Z:" + site + ",A,2,2,1,60"};
    // rescued: XT:A:M, no X0/X1; three differences away from the site, over the budget of 2
    const std::vector<std::string> want_resc = {"resc", "99", "c1", "481", "29", "32M", "=", "516", "60", resc, q32,
        "XT:A:M", er.f[11], "SM:i:29", "AM:i:29", "XM:i:3", "XO:i:0", "XG:i:0", er.f[12],
        "XV:Z:" + site + ",A,3,2,1,29"};
    CHECK(k[n - 4] == want_own);
    CHECK(k[n - 2] == want_resc);
    // The mates: stock's records with the pair's flags, mate fields and MAPQs (their own XT, NM,
    // X0/X1, XM/XO/XG and MD unchanged).
    auto mate = [&](std::vector<std::string> r, const char *mapq, const char *pnext, const char *sm, const char *am) {
        r[1] = "147", r[4] = mapq, r[6] = "=", r[7] = pnext, r[8] = "-60";
        for (std::string &t : r)
            if (t.rfind("SM:i:", 0) == 0) t = std::string("SM:i:") + sm;
            else if (t.rfind("AM:i:", 0) == 0) t = std::string("AM:i:") + am;
        return r;
    };
    CHECK(k[n - 3] == mate(s[n - 3], "60", "490", "37", "25"));
    CHECK(k[n - 1] == mate(s[n - 1], "29", "481", "29", "29"));
    CHECK(s[n - 3][5] == "25M");
    CHECK(s[n - 1][5] == "25M");
}

TEST_CASE("aln_known_snvs: unrelated sites retain safe methylation stopping and confidence"
          * doctest::test_suite("unit/aln_known_snvs")) {
    bwa_tests::TempDir tmp{"ks_stop_unrelated"}; const std::string dir = tmp.path();
    uint32_t state = 15652; const std::string read = aln_test::lcg_dna(state, 30);
    std::string copy = read;
    for (int p : {7, 19}) copy[p] = read[p] == 'A' ? 'G' : 'A';
    std::string genome = aln_test::lcg_dna(state, 100) + read + aln_test::lcg_dna(state, 100);
    for (int i = 0; i < 5000; ++i) genome += copy + "GAGAG";
    genome += aln_test::lcg_dna(state, 1000);
    const int site = genome.size() - 50;
    const std::string prefix = tmp.unique_file();
    REQUIRE(aln_test::write_test_fasta(prefix, {{"c1", genome}}));
    build_meth_index(prefix);
    std::string converted = read; for (char &b : converted) if (b == 'C') b = 'T';
    const std::string fq = write_fq(dir + "/r.fq", {{"r", converted}});
    const std::string vcf = write_vcf(dir + "/s.vcf", {row("c1", site + 1, genome[site], std::string(1, other(genome[site])))});
    bwa_tests::StderrCapture quiet;
    std::vector<std::string> stock;
    for (int mode = 0; mode < 3; ++mode) {
        const std::string sam = dir + "/out.sam";
        std::vector<std::string> args = {"--meth", "-f", sam};
        if (mode) args.insert(args.end(), {"--known-snvs", vcf});
        if (mode == 2) args.push_back("-N");
        args.insert(args.end(), {prefix, fq}); REQUIRE(run(args) == 0);
        require_records(slurp(sam), {"r"}, false); const auto rec = records(slurp(sam));
        if (!mode) stock = rec[0];
        else if (mode == 1) CHECK(rec[0] == stock);
        if (mode < 2) { CHECK(rec[0][4] == "37"); CHECK(tag(rec[0], "X0:i:") == "X0:i:1"); }
        else { CHECK(rec[0][4] == "0"); CHECK(tag(rec[0], "X0:i:").empty()); }
    }
}

namespace {

// One exact locus for `read` on c1 and 5000 two-difference copies on c2: without an early stop the
// genomic search spends its whole lookup budget on the copies and the read loses its confidence.
// The copies differ by transversions, which no conversion model scores as a match. The cases below
// share it (offsite()); each case's own files go in its own TempDir.
struct OffsiteStop {
    bwa_tests::TempDir tmp{"ks_stop_offsite"};
    std::string prefix, read, c1;
    static constexpr int locus = 100;

    OffsiteStop()
    {
        uint32_t state = 15652; read = aln_test::lcg_dna(state, 30);
        std::string copy = read;
        for (int p : {7, 19}) copy[p] = std::string("CATG")[std::string("ACGT").find(read[p])];
        c1 = aln_test::lcg_dna(state, locus) + read + aln_test::lcg_dna(state, 1000);
        std::string c2;
        for (int i = 0; i < 5000; ++i) c2 += copy + "GAGAG";
        c2 += aln_test::lcg_dna(state, 1000);
        prefix = tmp.unique_file();
        REQUIRE(aln_test::write_test_fasta(prefix, {{"c1", c1}, {"c2", c2}}));
        build_meth_index(prefix);
    }
    // The read's bases from either strand, with `base` at c1 position `site` when one is given.
    std::string oriented(bool reverse, int site = -1, char base = 0) const
    {
        std::string r = read;
        if (site >= 0) r[site - locus] = base;
        return reverse ? revcomp(r) : r;
    }
    // A single-end read (C->T converted, as read 1) from either strand of the exact locus.
    std::string fragment(bool reverse, int site = -1, char base = 0) const
    {
        std::string r = oriented(reverse, site, base);
        for (char &b : r) if (b == 'C') b = 'T';
        return r;
    }
};

OffsiteStop &offsite()
{
    static OffsiteStop x;
    return x;
}

// A transversion: no conversion model scores it as a match, on either converted copy.
char transversion(char b) { return std::string("CATG")[std::string("ACGT").find(b)]; }

// The default --meth lookup budget, which the genomic search spends in full when it cannot stop.
int full_lookup_budget()
{
    const bwa_aln::aln_meth_opt_t mo;
    return mo.p_best + mo.p_second;
}

// Direct searches against offsite(): its converted index and original reference, and the ALT
// windows of each site set (c1 positions, each with its transversion as ALT), built once per case.
struct OffsiteSearch {
    struct site_set_t {
        bwa_aln::ks_sites_t sites;
        bwa_aln::ks_windows_t windows;
    };
    struct outcome_t {
        bwa_aln::aln_meth_read_t window, genomic;
        bwa_aln::aln_meth_alt_bound_t certificate;
    };
    bwa_tests::TempDir tmp{"ks_stop_search"};
    FMI_search fmi{(offsite().prefix + ".meth").c_str()};
    bwa_aln::bwt_t *bwt = nullptr;
    bwa_aln::aln_meth_ref_t mr;
    bwa_aln::gap_opt_t opt; // -n 2 -o 0, unseeded
    std::map<std::vector<int>, std::unique_ptr<site_set_t>> built;

    OffsiteSearch()
    {
        fmi.load_index(/*load_pac=*/false, 1, /*build_seed_jump=*/false, BWAMEM_HUGETLB_OFF, /*allow_shm=*/false);
        bwt = bwa_aln::bwt_from_fmi(&fmi);
        REQUIRE(bwa_aln::aln_meth_ref_load(&mr, offsite().prefix.c_str(), fmi.idx->bns, 1) == 0);
        bwa_aln::gap_opt_t *g = bwa_aln::gap_init_opt();
        opt = *g; free(g);
        opt.fnr = -1; opt.max_diff = 2; opt.max_gapo = 0; opt.seed_len = 1024;
    }
    ~OffsiteSearch()
    {
        for (auto &b : built) bwa_aln::ks_windows_free(&b.second->windows);
        bwa_aln::aln_meth_ref_free(&mr);
        bwa_aln::bwt_destroy(bwt);
    }
    site_set_t &site_set(const std::vector<int> &site_pos)
    {
        std::unique_ptr<site_set_t> &w = built[site_pos];
        if (w) return *w;
        w.reset(new site_set_t);
        const std::string &c1 = offsite().c1;
        std::vector<std::string> rows;
        for (int p : site_pos) rows.push_back(row("c1", p + 1, c1[p], std::string(1, transversion(c1[p]))));
        const std::string vcf = write_vcf(tmp.file("s" + std::to_string(built.size()) + ".vcf"), rows);
        REQUIRE(bwa_aln::ks_load(vcf.c_str(), mr.bns, mr.pac, &w->sites) == 0);
        bwa_tests::StderrCapture quiet;
        REQUIRE(bwa_aln::ks_windows_build(&w->sites, mr.bns, mr.pac, 1, 1, &w->windows) == 0);
        return *w;
    }
    // `bases` (reference orientation of its strand, unconverted) as a one-read batch, converted
    // as read read_no (read 1: C->T, read 2: G->A). Free with free_read.
    bwa_aln::bwa_seq_t *read(const std::string &bases, int read_no)
    {
        std::string converted = bases;
        for (char &b : converted) if (b == (read_no ? 'C' : 'G')) b = read_no ? 'T' : 'A';
        const std::string fq = write_fq(tmp.file("r.fq"), {{"r", converted}});
        bwa_aln::bwa_seqio_t *input = bwa_aln::bwa_seq_open(fq.c_str());
        int n = 0;
        bwa_aln::bwa_seq_t *reads = bwa_aln::bwa_read_seq(input, 1, &n, BWA_MODE_COMPREAD, 0);
        bwa_aln::bwa_seq_close(input);
        REQUIRE(n == 1);
        return reads;
    }
    // The window search, then (with `genomic`) the genomic search with its candidates, certificate
    // and flags, as ks_meth_search composes them.
    outcome_t search(const std::vector<int> &site_pos, int read_no, const std::string &bases, int budget,
                     const bwa_aln::gap_opt_t &gap, bool genomic = true)
    {
        site_set_t &s = site_set(site_pos);
        mr.known_snvs = &s.sites;
        bwa_aln::aln_meth_opt_t mo; mo.p_best = budget; mo.p_second = 0; bwa_aln::aln_meth_opt_finish(&mo);
        bwa_aln::bwa_seq_t *reads = read(bases, read_no);
        const bwa_aln::gap_opt_t discovery = bwa_aln::ks_discovery_opt(gap);
        outcome_t out;
        bwa_aln::aln_meth_search_batch(s.windows.bwt, s.windows.mref, &mo, 1, reads, &discovery, 1, read_no, 0, 1,
                                       1000, nullptr, &out.window, nullptr, nullptr, true);
        out.certificate = bwa_aln::aln_meth_alt_certificate(out.window);
        if (genomic) {
            std::vector<bwa_aln::aln_meth_candidate_t> extra;
            bwa_aln::ks_meth_candidates(&s.windows, mr.bns, out.window, &extra);
            const int flags = out.window.flags;
            bwa_aln::aln_meth_search_batch(bwt, &mr, &mo, 1, reads, &gap, 1, read_no, 0, 1, 1000, nullptr,
                                           &out.genomic, &extra, &out.certificate, false, &flags);
        }
        bwa_aln::bwa_free_read_seq(1, reads);
        mr.known_snvs = nullptr;
        return out;
    }
};

} // namespace

TEST_CASE("aln_known_snvs: a site its ALT window holds beside the read keeps methylation stopping and confidence"
          * doctest::test_suite("unit/aln_known_snvs")) {
    const OffsiteStop &x = offsite();
    bwa_tests::TempDir tmp{"ks_stop_beside"}; const std::string dir = tmp.path();
    // The site is 50 bases past the read: its ALT window holds an exact copy of the read that does
    // not cover the site, so no candidate placement covers a site and the record must not change.
    const int site = OffsiteStop::locus + 79;
    const std::string fq = write_fq(dir + "/r.fq", {{"fwd", x.fragment(false)}, {"rev", x.fragment(true)}});
    const std::string vcf = write_vcf(dir + "/s.vcf", {row("c1", site + 1, x.c1[site], std::string(1, other(x.c1[site])))});
    bwa_tests::StderrCapture quiet;
    std::vector<std::vector<std::string>> stock;
    for (bool known : {false, true}) {
        const std::string sam = dir + "/out.sam";
        std::vector<std::string> args{"--meth", "-n", "2", "-o", "0", "-l", "1024", "-f", sam};
        if (known) args.insert(args.end(), {"--known-snvs", vcf});
        args.insert(args.end(), {x.prefix, fq}); REQUIRE(run(args) == 0);
        require_records(slurp(sam), {"fwd", "rev"}, false); const auto rec = records(slurp(sam));
        if (!known) stock = rec;
        for (size_t i = 0; i < rec.size(); ++i) {
            CAPTURE(known); CAPTURE(rec[i][0]);
            CHECK(rec[i] == stock[i]);
            CHECK(rec[i][2] == "c1"); CHECK(rec[i][3] == std::to_string(OffsiteStop::locus + 1));
            CHECK(rec[i][4] == "37"); CHECK(tag(rec[i], "X0:i:") == "X0:i:1"); CHECK(tag(rec[i], "X1:i:") == "X1:i:0");
            CHECK(tag(rec[i], "XV:Z:").empty());
        }
    }
}

TEST_CASE("aln_known_snvs: an ALT placement covering its site within the retained differences still blocks stopping"
          * doctest::test_suite("unit/aln_known_snvs")) {
    const OffsiteStop &x = offsite();
    bwa_tests::TempDir tmp{"ks_stop_under"}; const std::string dir = tmp.path();
    // The read covers the site with the reference base, so the site's ALT window places it with one
    // difference: inside the retained window (best cost plus one mismatch). The stop stays refused,
    // the copies exhaust the lookup budget and the placement is reported without confidence. (The
    // run without --known-snvs, MAPQ 37 and X0 1, is the previous case's.)
    const int site = OffsiteStop::locus + 10;
    const std::string fq = write_fq(dir + "/r.fq", {{"fwd", x.fragment(false)}, {"rev", x.fragment(true)}});
    const std::string vcf = write_vcf(dir + "/s.vcf", {row("c1", site + 1, x.c1[site], std::string(1, other(x.c1[site])))});
    bwa_tests::StderrCapture quiet;
    const std::string sam = dir + "/out.sam";
    REQUIRE(run({"--meth", "-n", "2", "-o", "0", "-l", "1024", "-f", sam, "--known-snvs", vcf, x.prefix, fq}) == 0);
    require_records(slurp(sam), {"fwd", "rev"}, false);
    for (const auto &r : records(slurp(sam))) {
        CAPTURE(r[0]);
        CHECK(r[2] == "c1"); CHECK(r[3] == std::to_string(OffsiteStop::locus + 1)); CHECK(r[5] == "30M");
        CHECK(r[4] == "0"); CHECK(tag(r, "X0:i:").empty()); CHECK(tag(r, "X1:i:").empty());
        CHECK(tag(r, "XV:Z:").find("c1:" + std::to_string(site + 1) + ":") != std::string::npos);
    }
}

TEST_CASE("aln_known_snvs: the ALT-window stop certificate bounds complete site-covering placements only"
          * doctest::test_suite("unit/aln_known_snvs")) {
    namespace ba = bwa_aln;
    const OffsiteStop &x = offsite();
    OffsiteSearch s;
    const int full_budget = full_lookup_budget(), locus = OffsiteStop::locus;
    for (int read_no : {0, 1}) for (bool reverse : {false, true}) {
        CAPTURE(read_no); CAPTURE(reverse);
        // Beside the read: the window index holds the read exactly, but no placement covers the
        // site. The width bound is 0; the site bound is one past the difference limit.
        auto beside = s.search({locus + 79}, read_no, x.oriented(reverse), full_budget, s.opt);
        CHECK(beside.window.flags == 0);
        CHECK(beside.window.full_width_bound.valid); CHECK(beside.window.full_width_bound.value == 0);
        CHECK(beside.window.site_bound.valid); CHECK(beside.window.site_bound.value == 3);
        CHECK(beside.certificate.valid); CHECK(beside.certificate.value == 3);
        CHECK(beside.genomic.flags == 0); CHECK(beside.genomic.n_lookups == 1);
        CHECK(beside.genomic.c1 == 1); CHECK(beside.genomic.c2 == 0);
        // Under the read, which carries the ALT: the ALT placement is exact. (The window search
        // alone: this site's REF read is the case after the loop.)
        auto under_alt = s.search({locus + 10}, read_no, x.oriented(reverse, locus + 10, transversion(x.c1[locus + 10])),
                                  full_budget, s.opt, false);
        CHECK(under_alt.window.flags == 0); CHECK(under_alt.window.site_bound.value == 0);
        CHECK(under_alt.certificate.value == 0);
        // Under the read at either end: the placement still covers its site. One base outside
        // either end: no placement of the read covers the site.
        for (int edge : {locus, locus + 29}) CHECK(s.search({edge}, read_no, x.oriented(reverse), full_budget, s.opt, false).certificate.value == 1);
        for (int edge : {locus - 1, locus + 30}) CHECK(s.search({edge}, read_no, x.oriented(reverse), full_budget, s.opt, false).certificate.value == 3);
    }
    // Under the read, which carries the reference base: the ALT placement has one mismatch. It could
    // be retained (best cost plus one mismatch), so the certificate refuses the stop and the genomic
    // search spends its whole budget on the copies.
    auto under_ref = s.search({locus + 10}, 1, x.oriented(false), full_budget, s.opt);
    CHECK(under_ref.window.flags == 0); CHECK(under_ref.window.site_bound.value == 1);
    CHECK(under_ref.certificate.value == 1);
    CHECK(under_ref.genomic.flags == ba::ALN_METH_BIND_BEST); CHECK(under_ref.genomic.n_lookups == full_budget);
    // Two windows hold the read beside their sites. A lookup budget of one row binds the window
    // search before it has seen both: it proved nothing about site-covering placements, so the
    // certificate falls back to the width bound and the genomic result stays incomplete.
    auto bound = s.search({locus + 79, locus + 85}, 1, x.oriented(false), 1, s.opt);
    CHECK(bound.window.flags == ba::ALN_METH_BIND_BEST);
    CHECK(bound.window.site_bound.valid); CHECK(bound.window.site_bound.value == 3);
    CHECK(bound.certificate.valid); CHECK(bound.certificate.value == bound.window.full_width_bound.value);
    CHECK(bound.certificate.value == 0);
    CHECK((bound.genomic.flags & ba::ALN_METH_BIND_BEST) != 0);
    // The same windows searched to completion certify the stop.
    auto unbound = s.search({locus + 79, locus + 85}, 1, x.oriented(false), full_budget, s.opt);
    CHECK(unbound.window.flags == 0); CHECK(unbound.certificate.value == 3);
    CHECK(unbound.genomic.flags == 0); CHECK(unbound.genomic.n_lookups == 1);
    // A seeded search proposes every ungapped window only up to the seed's difference limit.
    for (int seed_diff : {0, 1, 2, 5}) {
        ba::gap_opt_t seeded = s.opt; seeded.seed_len = 20; seeded.max_seed_diff = seed_diff;
        CAPTURE(seed_diff);
        CHECK(s.search({locus + 79}, 1, x.oriented(false), full_budget, seeded, false).window.site_bound.value
              == std::min(seed_diff, 2) + 1);
    }
    // -n and -k accept INT_MAX: the horizon saturates rather than overflowing. (Window search only,
    // under a one-row budget: the bound is set before the search runs.)
    ba::gap_opt_t unlimited = s.opt; unlimited.max_diff = INT_MAX; unlimited.max_entries = 1000;
    for (bool seeded : {false, true}) {
        CAPTURE(seeded);
        ba::gap_opt_t g = unlimited;
        if (seeded) { g.seed_len = 20; g.max_seed_diff = INT_MAX; }
        const auto o = s.search({locus + 79}, 1, x.oriented(false), 1, g, false);
        CHECK(o.window.site_bound.valid); CHECK(o.window.site_bound.value == INT_MAX);
    }
}

TEST_CASE("aln_known_snvs: the stop certificate counts only equal-span gapped window histories, at zero"
          * doctest::test_suite("unit/aln_known_snvs")) {
    namespace ba = bwa_aln;
    const OffsiteStop &x = offsite();
    OffsiteSearch s;
    // A 30-base read from the unique part of c1, past the read the other cases use.
    const int beg = 400, site = beg + 3;
    const std::string ref = x.c1.substr(beg, 31);
    ba::gap_opt_t gapped = s.opt; gapped.max_diff = 3; gapped.max_gapo = 2;
    // Equal span: one base deleted after read base 8 and one inserted after read base 23, and the
    // reference base at the site (one mismatch against its ALT window). Ungapped, the 15 shifted
    // bases between the gaps mismatch: the only histories within -n 3 covering the site are this
    // one (two gap opens, one mismatch) and its equivalent gap placements, so the bound is 0.
    const std::string equal = ref.substr(0, 8) + ref.substr(9, 15) + "A" + ref.substr(24, 6);
    REQUIRE(equal.size() == 30);
    auto e = s.search({site}, 1, equal, full_lookup_budget(), gapped, false);
    CHECK(e.window.flags == 0); CHECK(e.window.site_bound.valid); CHECK(e.window.site_bound.value == 0);
    // Unequal span: one base deleted, and the read carries the ALT. That history covers the site
    // with no mismatch, but its window is a base longer than the read: the bound ignores it. (One
    // gap open at most, so no second gap can restore the read's length.)
    std::string alt = x.c1.substr(beg, 31); alt[site - beg] = transversion(alt[site - beg]);
    const std::string longer = alt.substr(0, 15) + alt.substr(16, 15);
    REQUIRE(longer.size() == 30);
    ba::gap_opt_t one_gap = gapped; one_gap.max_gapo = 1;
    auto l = s.search({site}, 1, longer, full_lookup_budget(), one_gap, false);
    CHECK(l.window.flags == 0); CHECK(l.window.site_bound.valid); CHECK(l.window.site_bound.value == 4);
    CHECK(l.window.n_rescored > 0);
}

TEST_CASE("aln_known_snvs: ks_meth_search passes the site certificate and the window flags to the genomic search"
          * doctest::test_suite("unit/aln_known_snvs")) {
    namespace ba = bwa_aln;
    const OffsiteStop &x = offsite();
    OffsiteSearch s;
    const int locus = OffsiteStop::locus;
    ba::aln_opt_t o; ba::aln_opt_init(&o); o.gap = s.opt; o.gap.n_threads = 1;
    ba::aln_meth_opt_t mo; ba::aln_meth_opt_finish(&mo);
    o.meth = &mo; o.meth_ref = &s.mr;
    auto run_search = [&](const std::vector<int> &site_pos, const std::string &bases) {
        OffsiteSearch::site_set_t &w = s.site_set(site_pos);
        s.mr.known_snvs = o.known_snvs = &w.sites; o.ks_windows = &w.windows;
        ba::bwa_seq_t *reads = s.read(bases, 1);
        ba::aln_meth_read_t res;
        ba::ks_meth_search(&o, s.bwt, 1, reads, 1, 0, 0, nullptr, nullptr, &res);
        ba::bwa_free_read_seq(1, reads);
        s.mr.known_snvs = nullptr;
        return res;
    };
    // The site beside the read: the site bound (3), not the width bound (0), lets the read stop at
    // its exact locus after one lookup instead of spending the budget on the copies.
    const auto beside = run_search({locus + 79}, x.oriented(false));
    CHECK(beside.flags == 0); CHECK(beside.n_lookups == 1); CHECK(beside.c1 == 1); CHECK(beside.c2 == 0);
    // (That the window flags reach the genomic slot before it searches is pinned in test_aln_meth.cpp,
    // where TAPS neutral scoring makes the difference visible.)
}

TEST_CASE("aln_known_snvs: conservative downstream counting excludes both conversion-source alleles"
          * doctest::test_suite("unit/aln_known_snvs")) {
    namespace ba = bwa_aln;
    // This recipe belongs to downstream counting, not the aligner's R/A/H or twin semantics.
    auto old_recipe = [](const ba::ks_vote_t &v, int mapq, int same, int twin_mapq) {
        return (v.allele == 'R' || v.allele == 'A') && mapq >= 20 && same == 1 && twin_mapq >= 20
            && v.away <= v.max_diff - 1;
    };
    for (int hyp : {0, 1}) for (bool cpg : {false, true})
        for (int context : {0, 1, 2, 3}) { // active, wrong neighbor, ambiguous neighbor, boundary
        const uint8_t source = hyp ? 1 : 2, target = hyp ? 3 : 0;
        for (int kind : {0, 1, 2}) {
            const uint8_t refbase = kind == 2 ? (hyp ? 0 : 1) : kind == 1 ? target : source;
            const uint8_t altbase = kind == 2 ? (hyp ? 2 : 3) : kind == 1 ? source : (hyp ? 0 : 3);
            const int p = context == 3 ? (hyp ? 19 : 0) : 10;
            const int neighbor = hyp ? p + 1 : p - 1;
            uint8_t pac[5] = {}; const uint8_t context_base = hyp ? 2 : 1;
            auto set = [&](int i, uint8_t b) { pac[i / 4] |= b << ((3 - i % 4) * 2); };
            set(p, refbase);
            if (neighbor >= 0 && neighbor < 20) set(neighbor, context == 1 ? 0 : context_base);
            bntann1_t ann = {}; ann.name = const_cast<char *>("c1"); ann.len = 20;
            bntamb1_t amb = {}; amb.offset = neighbor; amb.len = 1;
            bntseq_t bns = {}; bns.n_seqs = 1; bns.anns = &ann; bns.l_pac = 20;
            if (context == 2) { bns.n_holes = 1; bns.ambs = &amb; }
            ba::ks_sites_t sites; sites.by_contig = {{{p, refbase, altbase}}};
            auto *opt = ba::gap_init_opt(); opt->fnr = -1; opt->max_diff = 2;
            // Independent downstream context rule: same contig, no N in either CpG base.
            const bool active = !cpg || (context == 0);
            const bool exclude = active && (refbase == source || altbase == source);
            for (uint8_t observed : {refbase, altbase}) {
                uint8_t seq[20] = {}; seq[p] = observed;
                ba::bwa_seq_t s = {}; s.len = s.clip_len = 20; s.type = BWA_TYPE_UNIQUE; s.seq = seq;
                std::vector<ba::ks_vote_t> votes;
                ba::ks_site_votes(&sites, &bns, opt, &s, -1, hyp, &votes, pac, cpg);
                REQUIRE(votes.size() == 1);
                CHECK(votes[0].conv == active);
                // Placeable literal source can be informative even when its counterpart is hidden.
                const bool old_count = old_recipe(votes[0], 37, 1, 37);
                const bool count = !exclude && old_count;
                if (exclude) CHECK_FALSE(count);
                else CHECK(count); // unchanged inactive-context and nonsource R/A controls
                if (active && kind == 1) {
                    CHECK(votes[0].allele == (observed == target ? 'H' : 'A'));
                    CHECK(old_count == (observed == source));
                }
            }
            free(opt);
        }
    }
    // Actual alignments reproduce the C>A paralog and complementary G>T observation asymmetry.
    bwa_tests::TempDir tmp{"ks_counting"}; const std::string dir = tmp.path();
    uint32_t state = 7731; std::string genome = aln_test::lcg_dna(state, 1000); genome[210] = 'C';
    std::string paralog = genome; paralog[210] = 'T';
    bwa_tests::StderrCapture quiet;
    for (bool ga : {false, true}) {
        const std::string prefix = tmp.unique_file(ga ? ".ga.fa" : ".ct.fa");
        REQUIRE(aln_test::write_test_fasta(prefix, {{"c1", ga ? revcomp(genome) : genome},
                                                 {"c2", ga ? revcomp(paralog) : paralog}}));
        build_meth_index(prefix);
        std::string a = genome.substr(200, 35); a[10] = 'A';
        const std::string fq = write_fq(dir + "/r.fq", {{"ref", genome.substr(200, 35)}, {"alt", a}});
        const std::string vcf = write_vcf(dir + "/s.vcf", {row("c1", ga ? 790 : 211, ga ? 'G' : 'C', ga ? "T" : "A")});
        const std::string sam = dir + "/out.sam";
        REQUIRE(run({"--meth", "--known-snvs", vcf, "-n", "2", "-o", "0", "-f", sam, prefix, fq}) == 0);
        require_records(slurp(sam), {"ref", "alt"}, false);
        const auto rec = records(slurp(sam));
        for (int i = 0; i < 2; ++i) {
            CHECK(rec[i][2] == "c1"); CHECK(rec[i][3] == (ga ? "766" : "201"));
            CHECK(rec[i][4] == "23"); CHECK(rec[i][5] == "35M");
            CHECK(tag(rec[i], "XG:Z:") == (ga ? "XG:Z:GA" : "XG:Z:CT"));
            const std::string site = ga ? "c1:790:G>T" : "c1:211:C>A";
            CHECK(tag(rec[i], "XV:Z:") == "XV:Z:" + site + (i ? ",A,0,2,0,0" : ",R,0,2,1,23"));
            // An independent consumer interprets the numeric twin fields as placement only.
            const std::string xv = tag(rec[i], "XV:Z:");
            const size_t last = xv.rfind(','), prior = xv.rfind(',', last - 1);
            const bool old_count = std::stoi(rec[i][4]) >= 20 && xv.substr(prior + 1, last - prior - 1) == "1"
                && std::stoi(xv.substr(last + 1)) >= 20;
            CHECK(old_count == (i == 0));
            const bool source_present = ga ? site.find("G>") != std::string::npos : site.find("C>") != std::string::npos;
            CHECK_FALSE((old_count && !source_present)); // both alleles excluded, without altering tags
        }
    }
    // TAPS/CpG: the literal source is identifiable ALT, but the target observation is hidden.
    genome[210] = 'T'; genome[211] = 'G';
    for (bool ga : {false, true}) {
        const std::string prefix = tmp.unique_file(ga ? ".tapsga.fa" : ".tapsct.fa");
        REQUIRE(aln_test::write_test_fasta(prefix, {{"c1", ga ? revcomp(genome) : genome}}));
        build_meth_index(prefix);
        std::string alt = genome.substr(200, 35); alt[10] = 'C';
        const std::string fq = write_fq(dir + "/taps.fq", {{"ref", genome.substr(200, 35)}, {"alt", alt}});
        const std::string vcf = write_vcf(dir + "/taps.vcf", {row("c1", ga ? 790 : 211, ga ? 'A' : 'T', ga ? "G" : "C")});
        const std::string sam = dir + "/taps.sam";
        REQUIRE(run({"--meth=taps", "--meth-context", "cpg", "--known-snvs", vcf, "-n", "2", "-o", "0", "-f", sam, prefix, fq}) == 0);
        require_records(slurp(sam), {"ref", "alt"}, false);
        const auto rec = records(slurp(sam));
        for (int i = 0; i < 2; ++i) {
            CHECK(rec[i][3] == (ga ? "766" : "201")); CHECK(rec[i][4] == "37"); CHECK(rec[i][5] == "35M");
            CHECK(tag(rec[i], "XG:Z:") == (ga ? "XG:Z:GA" : "XG:Z:CT"));
            const std::string site = ga ? "c1:790:A>G" : "c1:211:T>C";
            CHECK(tag(rec[i], "XV:Z:") == "XV:Z:" + site + (i ? ",A,0,2,1,37" : ",H,0,2,NA,NA"));
            const std::string xv = tag(rec[i], "XV:Z:");
            const bool old_count = xv.find(",A,0,2,1,37") != std::string::npos;
            CHECK(old_count == (i == 1));
            const char source = ga ? 'G' : 'C';
            const bool excluded = site.back() == source; // active CpG in either allele context
            CHECK_FALSE((old_count && !excluded));
        }
    }
    // Saved C>A/paralog observation pattern: both actual reads MAPQ23, but only C's A twin passes.
    ba::ks_vote_t literal_c = {}, literal_a = {}; literal_c.allele = 'R'; literal_a.allele = 'A';
    literal_c.max_diff = literal_a.max_diff = 2;
    CHECK(old_recipe(literal_c, 23, 1, 23)); CHECK_FALSE(old_recipe(literal_a, 23, 0, 0));
    const bool ct_source_exclusion = true; // C>A on active CT: exclude BOTH, not just the bad twin.
    CHECK_FALSE((!ct_source_exclusion && old_recipe(literal_c, 23, 1, 23)));
    CHECK_FALSE((!ct_source_exclusion && old_recipe(literal_a, 23, 0, 0)));
}

TEST_CASE("aln_known_snvs: overlapping mates swap the same site in plain and methylation twins"
          * doctest::test_suite("unit/aln_known_snvs")) {
    bwa_tests::TempDir tmp{"ks_overlap"}; const std::string dir = tmp.path();
    uint32_t state = 76323;
    std::string c1 = aln_test::lcg_dna(state, 1000); c1[500] = 'A';
    std::string c2 = aln_test::lcg_dna(state, 100) + c1.substr(475, 55) + aln_test::lcg_dna(state, 100);
    c2[125] = 'G';
    const std::string prefix = tmp.unique_file();
    FILE *fp = fopen(prefix.c_str(), "w"); REQUIRE(fp != nullptr);
    fprintf(fp, ">c1\n%s\n>c2\n%s\n", c1.c_str(), c2.c_str()); REQUIRE(fclose(fp) == 0);
    REQUIRE(bwa_idx_build(prefix.c_str(), prefix.c_str(), 0, 3) == 0);
    build_meth_index(prefix);
    const std::string vcf = write_vcf(dir + "/s.vcf", {row("c1", 501, 'A', "G")});
    bwa_tests::StderrCapture quiet;
    for (bool meth : {false, true}) {
        std::vector<std::vector<std::string>> result[2];
        for (int swapped = 0; swapped < 2; ++swapped) {
            std::string fq[2];
            for (int end = 0; end < 2; ++end) {
                const int beg = end ? 490 : 480;
                std::string read = c1.substr(beg, 35);
                if (swapped) read[500 - beg] = 'G'; // independent oracle swaps BOTH input mates
                if (meth) for (char &b : read) if (b == 'C') b = 'T';
                if (end) read = revcomp(read);
                fq[end] = write_fq(dir + "/r" + std::to_string(end) + ".fq", {{"r", read}});
            }
            const std::string sam = dir + "/" + std::to_string(swapped) + ".sam";
            std::vector<std::string> args{"--known-snvs", vcf, "-n", "2", "-l", "1024", "-f", sam};
            if (meth) args.push_back("--meth");
            args.insert(args.end(), {prefix, fq[0], fq[1]}); REQUIRE(run(args) == 0);
            require_records(slurp(sam), {"r"}, true); result[swapped] = records(slurp(sam));
        }
        for (int end = 0; end < 2; ++end) {
            const auto &original = result[0][end], &oracle = result[1][end];
            REQUIRE(original[2] == "c1"); REQUIRE(std::stoi(original[4]) > 0);
            const bool same = !(std::stoi(oracle[1]) & 4) && oracle[2] == original[2]
                && ((std::stoi(oracle[1]) ^ std::stoi(original[1])) & 16) == 0
                && std::abs(std::stoi(oracle[3]) - std::stoi(original[3])) <= 5;
            CHECK_FALSE(same); CHECK(oracle[4] == "0");
            const std::string xv = tag(original, "XV:Z:");
            CHECK(xv == "XV:Z:c1:501:A>G,R,0,2," + std::to_string(same) + "," + oracle[4]);
        }
    }
}

TEST_CASE("aln_known_snvs: paired end under --meth, an over-budget ALT mate is placed and pairs"
          * doctest::test_suite("unit/aln_known_snvs")) {
    KsMethIndex x;
    bwa_tests::TempDir scratch{"ks_case"}; // this case's VCF, reads and output, removed when it ends
    const std::string dir = scratch.path();
    // Unmethylated directional fragments from c1[400, 600): read 1 is the fragment's start with
    // every C read as T, read 2 the reverse complement of its end. The site is the first A or T in
    // c1[495, 505) with G as its ALT (no conversion explains either base on read 1's strand).
    int q = 495;
    while (q < 505 && x.c1[q] != 'A' && x.c1[q] != 'T') ++q;
    REQUIRE(q < 505);
    const char ref = x.c1[q], alt = 'G';
    auto conv = [](std::string s) { for (char &c : s) if (c == 'C') c = 'T'; return s; };
    const std::string fq1 = dir + "/m1.fq", fq2 = dir + "/m2.fq";
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
    const std::string vcf = write_vcf(dir + "/m.vcf", {row("c1", q + 1, ref, std::string(1, alt))});
    const std::string stock = dir + "/m.stock.sam", ks = dir + "/m.ks.sam";
    bwa_tests::StderrCapture capture;
    REQUIRE(run({"--meth", "-f", stock, x.prefix, fq1, fq2}) == 0);
    const int ks_rc = run({"--meth", "--known-snvs", vcf, "-f", ks, x.prefix, fq1, fq2});
    REQUIRE(ks_rc == 0);
    const auto s = records(slurp(stock)), k = records(slurp(ks));
    REQUIRE(s.size() == k.size());
    std::vector<std::string> names;
    for (int i = 0; i < 28; ++i) names.push_back("bg" + std::to_string(i));
    names.push_back("own");
    require_records(slurp(stock), names, true); require_records(slurp(ks), names, true);
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
    // Independent ungapped CT oracle: literal MD includes conversions, NM does not.
    // Both mates see the same converted fragment in reference orientation.
    auto expected = [&](int end) {
        const int start = beg + (end ? 35 : 0);
        const std::string aligned = conv(frag).substr(end ? 35 : 0, 25);
        const std::string genome = x.c1.substr(start, 25);
        int literal_nm = 0, nm = 0;
        const std::string md = md_of(genome, aligned, &literal_nm);
        std::string xm(25, '.');
        for (int i = 0; i < 25; ++i) {
            const char rb = genome[i], qb = aligned[i];
            if (rb != qb && !(rb == 'C' && qb == 'T')) ++nm;
            if (rb != 'C' || (qb != 'C' && qb != 'T')) continue;
            const char context = x.c1[start + i + 1] == 'G' ? 'z'
                : x.c1[start + i + 2] == 'G' ? 'x' : 'h';
            xm[i] = qb == 'C' ? context - 'a' + 'A' : context;
        }
        CHECK(nm == (end ? 0 : 3));
        std::vector<std::string> record = {"own", end ? "147" : "99", "c1", std::to_string(start + 1), "60",
            "25M", "=", std::to_string(beg + (end ? 1 : 36)), end ? "-60" : "60", aligned,
            std::string(25, 'I'), "XT:A:U", "NM:i:" + std::to_string(nm), end ? "SM:i:37" : "SM:i:25",
            "AM:i:25", "X0:i:1", "X1:i:0", "MD:Z:" + md, end ? "XR:Z:GA" : "XR:Z:CT", "XG:Z:CT", "XM:Z:" + xm};
        if (!end) record.push_back("XV:Z:" + site + ",A,2,2,1,60");
        return record;
    };
    CHECK(k[n - 2] == expected(0));
    CHECK(k[n - 1] == expected(1));
}

TEST_CASE("aln_known_snvs: with -e, a gap's extensions are not differences, as in the gap search"
          * doctest::test_suite("unit/aln_known_snvs")) {
    // -e clears BWA_MODE_GAPE: the gap search then counts a gap as one difference whatever its
    // length, and the ALT-window cuts and XV:Z's `away` must count it the same way.
    bwa_tests::TempDir tmp{"ks_gape"};
    const std::string dir = tmp.path(), prefix = tmp.unique_file();
    uint32_t st = 4242;
    const std::string a = aln_test::lcg_dna(st, 600), b0 = aln_test::lcg_dna(st, 100);
    // "tie": its exact copy is a[100, 140); c2 holds it again with 5 extra reference bases after
    // its 20th base and a site (REF not the read's base, ALT the read's base) under its 31st base.
    const std::string tie = a.substr(100, 40), extra = "GATCA";
    std::string c2 = b0 + tie.substr(0, 20) + extra + tie.substr(20) + aln_test::lcg_dna(st, 100);
    const int site2 = 100 + 20 + 5 + 10; // 0-based, under tie[30]
    c2[site2] = other(tie[30]);
    // "del": a[450, 493) with a[470, 473) deleted (one gap open, two extensions), over a site at
    // a[460] where the read carries REF. (Its window, a[300, 600), does not hold "tie"'s exact copy,
    // which would end the window search at that copy before it reached c2's gapped hit.)
    const std::string del = a.substr(450, 20) + a.substr(473, 20);
    REQUIRE(aln_test::write_test_fasta(prefix, {{"c1", a}, {"c2", c2}}));
    {
        bwa_tests::StderrCapture quiet;
        REQUIRE(bwa_idx_build(prefix.c_str(), prefix.c_str(), 0, 3) == 0);
    }
    const std::string vcf = write_vcf(dir + "/gape.vcf", {row("c1", 461, a[460], std::string(1, other(a[460]))),
                                                          row("c2", site2 + 1, c2[site2], std::string(1, tie[30]))});
    const std::string fq = write_fq(dir + "/gape.fq", {{"tie", tie}, {"del", del}});
    SUBCASE("the ALT-window cuts") {
        // -O 0 -E 0: the deletion is free, so c2's ALT-credited placement (one gap: one
        // difference) ties c1's exact one in score, and the difference cut (the exact best's 0, plus
        // one) keeps it. Counting the gap's four extensions too would drop it.
        const std::string out = dir + "/tie.sam";
        bwa_tests::StderrCapture capture;
        REQUIRE(run({"--known-snvs", vcf, "-n", "2", "-e", "5", "-O", "0", "-E", "0", "-l", "1024", "-f", out,
                     prefix, fq}) == 0);
        const auto r = records(slurp(out));
        REQUIRE(r.size() == 2);
        REQUIRE(r[0][0] == "tie");
        CHECK(tag(r[0], "X0:i:") == "X0:i:2"); // was 1: the 5-base gap counted as five differences
        CHECK(r[0][4] == "0");
    }
    SUBCASE("XV:Z's away") {
        const std::string out = dir + "/del.sam";
        bwa_tests::StderrCapture capture;
        REQUIRE(run({"--known-snvs", vcf, "-n", "2", "-e", "3", "-l", "1024", "-f", out, prefix, fq}) == 0);
        const auto r = records(slurp(out));
        REQUIRE(r.size() == 2);
        REQUIRE(r[1][0] == "del");
        REQUIRE(r[1][5].find("3D") != std::string::npos); // placed with the deletion
        const std::string site = std::string("c1:461:") + a[460] + ">" + other(a[460]);
        CHECK(tag(r[1], "XV:Z:").rfind("XV:Z:" + site + ",R,1,2,", 0) == 0); // was away 3
    }
}

TEST_CASE("aln_known_snvs: under --meth-context cpg a converted base counts as its allele only at a CpG"
          * doctest::test_suite("unit/aln_known_snvs")) {
    KsIndex &x = fixture();
    FMI_search fmi(x.prefix.c_str());
    fmi.load_index(/*load_pac=*/true, 1, /*build_seed_jump=*/false);
    const bntseq_t *bns = fmi.idx->bns;
    const uint8_t *pac = fmi.idx->pac;
    auto *opt = bwa_aln::gap_init_opt();
    // A C>A site at c1[p]; a 25-base forward read on the C->T strand (hyp 1) over it, from c1[p - 10].
    auto votes = [&](int p, char read_base, int cpg, std::vector<uint8_t> *twin) {
        bwa_aln::ks_sites_t sites;
        sites.by_contig.assign(bns->n_seqs, {});
        sites.by_contig[0] = {{p, 1, 0}};
        sites.n_sites = 1;
        std::vector<uint8_t> seq(25);
        for (int i = 0; i < 25; ++i) seq[i] = nst_nt4_table[(int)x.c1[p - 10 + i]];
        seq[10] = nst_nt4_table[(int)read_base];
        bwa_aln::bwa_seq_t s{};
        s.seq = seq.data();
        s.len = 25;
        s.type = BWA_TYPE_UNIQUE;
        s.pos = (bwtint_t)(p - 10);
        std::vector<bwa_aln::ks_vote_t> v;
        const int first = bwa_aln::ks_site_votes(&sites, bns, opt, &s, -1, /*hyp=*/1, &v, pac, cpg);
        REQUIRE(v.size() == 1);
        if (twin && first == 0) bwa_aln::ks_twin_seq(&s, v[0], 1, twin);
        return v[0].allele;
    };
    int p_cpg = -1, p_plain = -1; // a C followed by G, and a C followed by another base
    for (int p = 20; p < 500 && (p_cpg < 0 || p_plain < 0); ++p) {
        if (x.c1[p] != 'C') continue;
        if (x.c1[p + 1] == 'G' && p_cpg < 0) p_cpg = p;
        else if (x.c1[p + 1] != 'G' && x.c1[p - 1] != 'C' && p_plain < 0) p_plain = p;
    }
    REQUIRE(p_cpg > 0);
    REQUIRE(p_plain > 0);
    CHECK(votes(p_plain, 'T', 0, nullptr) == 'R'); // --meth-context all: the conversion of a REF C
    CHECK(votes(p_plain, 'T', 1, nullptr) == 'O'); // cpg: a mismatch, as the alignment scores it
    CHECK(votes(p_cpg, 'T', 1, nullptr) == 'R');
    CHECK(votes(p_plain, 'C', 1, nullptr) == 'R');
    // The twin of an ALT read swaps the REF C in; it reads as T only where it would convert.
    std::vector<uint8_t> twin;
    CHECK(votes(p_cpg, 'A', 1, &twin) == 'A');
    REQUIRE(twin.size() == 25);
    CHECK(twin[10] == 3);
    twin.clear();
    CHECK(votes(p_plain, 'A', 1, &twin) == 'A');
    REQUIRE(twin.size() == 25);
    CHECK(twin[10] == 1);
    twin.clear();
    CHECK(votes(p_plain, 'A', 0, &twin) == 'A');
    REQUIRE(twin.size() == 25);
    CHECK(twin[10] == 3);
    free(opt);
}

TEST_CASE("aln_known_snvs: trimmed votes retain the search budget and twins retain only searched bases"
          * doctest::test_suite("unit/aln_known_snvs")) {
    namespace ba = bwa_aln;
    ba::gap_opt_t *opt = ba::gap_init_opt();
    opt->fnr = 0.04;
    bntann1_t ann = {}; ann.name = const_cast<char *>("c;:,%1"); ann.len = 100;
    bntseq_t bns = {}; bns.n_seqs = 1; bns.anns = &ann; bns.l_pac = 100;
    ba::ks_sites_t sites; sites.by_contig = {{{10, 0, 1}}};
    ba::bwa_seq_t s = {}; s.len = s.full_len = 80; s.clip_len = 25; s.type = BWA_TYPE_UNIQUE;
    std::vector<uint8_t> seq(80, 0); s.seq = seq.data();
    ba::bwa_cigar_t cigar[2] = {__cigar_create(FROM_M, 25), __cigar_create(FROM_S, 55)};
    s.cigar = cigar; s.n_cigar = 2;
    for (int strand : {0, 1}) {
        s.strand = strand;
        if (strand) { std::fill(seq.begin(), seq.end(), 3); cigar[0] = __cigar_create(FROM_S, 55); cigar[1] = __cigar_create(FROM_M, 25); }
        std::vector<ba::ks_vote_t> votes;
        REQUIRE(ba::ks_site_votes(&sites, &bns, opt, &s, -1, -1, &votes) == 0);
        REQUIRE(votes.size() == 1);
        CHECK(votes[0].max_diff == ba::bwa_cal_maxdiff(25, BWA_AVG_ERR, opt->fnr));
        std::vector<uint8_t> twin;
        ba::ks_twin_seq(&s, votes[0], -1, &twin);
        REQUIRE(twin.size() == 25);
        CHECK(twin[strand ? 14 : 10] == (strand ? 2 : 1));
        std::string xv;
        ba::ks_format_xv(&bns, votes, 0, 1, 37, &xv);
        CHECK(xv.find("c%3B%3A%2C%251:") == 0);
    }
    free(opt);
}

TEST_CASE("aln_known_snvs: listed pairing retains the selected duplicate's reference span"
          * doctest::test_suite("unit/aln_known_snvs")) {
    namespace ba = bwa_aln;
    KsIndex &x = fixture();
    FMI_search fmi(x.prefix.c_str());
    fmi.load_index(/*load_pac=*/true, 1, /*build_seed_jump=*/false);
    ba::bwt_t *bwt = ba::bwt_from_fmi(&fmi);
    ba::bwa_seq_t a = {}, b = {};
    a.type = b.type = BWA_TYPE_UNIQUE;
    a.full_len = a.clip_len = a.len = b.full_len = b.clip_len = b.len = 50;
    a.c1 = b.c1 = 1; a.mapQ = b.mapQ = a.seQ = b.seQ = 37;
    a.pos = 1000; b.pos = 1100; b.strand = 1;
    ba::bwa_seq_t *seqs[2] = {&a, &b};
    std::vector<ba::pe_hit_t> h[2] = {{{1000, 0, 0, 0, 0, 0, 0}},
        {{1100, 1, 3, 1, 0, 0, 0}, {1100, 1, 0, 0, 0, 0, 20}}};
    const std::vector<ba::pe_hit_t> *hits[2] = {h, h + 1};
    auto *pe = ba::bwa_init_pe_opt(); auto *gap = ba::gap_init_opt();
    pe->ref_span_insert = 1; pe->n_multi = pe->N_multi = 0;
    ba::isize_info_t ii = {}; ii.avg = -1;
    // Grouping resolves an empty SA batch even when both ends use listed hits.
    ba::bwa_pe_pair(fmi.idx->bns, bwt, 1, seqs, &ii, pe, gap, hits);
    CHECK(b.ref_shift == 20);
    ba::bwt_destroy(bwt);
    free(pe); free(gap);
}

TEST_CASE("aln_known_snvs: listed extra hits keep the primary's start on the other strand"
          * doctest::test_suite("unit/aln_known_snvs")) {
    namespace ba = bwa_aln;
    KsIndex &x = fixture();
    FMI_search fmi(x.prefix.c_str());
    fmi.load_index(/*load_pac=*/true, 1, /*build_seed_jump=*/false);
    ba::bwt_t *bwt = ba::bwt_from_fmi(&fmi);
    ba::bwa_seq_t a = {}, b = {};
    a.type = b.type = BWA_TYPE_UNIQUE;
    a.full_len = a.clip_len = a.len = b.full_len = b.clip_len = b.len = 50;
    a.c1 = b.c1 = 1; a.c2 = 1; a.mapQ = b.mapQ = a.seQ = b.seQ = 37;
    a.pos = 1000; b.pos = 1100; b.strand = 1;
    ba::bwa_seq_t *seqs[2] = {&a, &b};
    // End 0's primary is (1000, +); (1000, -) is an alternative at the same start.
    std::vector<ba::pe_hit_t> h[2] = {{{1000, 0, 0, 0, 0, 0, 0}, {1000, 1, 3, 1, 0, 0, 0}},
        {{1100, 1, 0, 0, 0, 0, 0}}};
    const std::vector<ba::pe_hit_t> *hits[2] = {h, h + 1};
    auto *pe = ba::bwa_init_pe_opt(); auto *gap = ba::gap_init_opt();
    pe->n_multi = pe->N_multi = 3;
    ba::isize_info_t ii = {}; ii.avg = -1;
    ba::bwa_pe_pair(fmi.idx->bns, bwt, 1, seqs, &ii, pe, gap, hits);
    CHECK(a.pos == 1000); CHECK((int)a.strand == 0);
    REQUIRE(a.n_multi == 1);
    CHECK(a.multi[0].pos == 1000); CHECK((int)a.multi[0].strand == 1); CHECK((int)a.multi[0].mm == 1);
    free(a.multi); free(b.multi);
    ba::bwt_destroy(bwt);
    free(pe); free(gap);
}

TEST_CASE("aln_known_snvs: twins use trimmed bases for single ends and paired mate copies"
          * doctest::test_suite("unit/aln_known_snvs")) {
    KsMethIndex x;
    bwa_tests::TempDir tmp{"ks_trim"}; const std::string dir = tmp.path();
    int p = 480; while (x.c1[p] == 'C') ++p;
    const char ref = x.c1[p], alt = other(ref, 'C');
    const std::string vcf = write_vcf(dir + "/sites.vcf", {row("c1", p + 1, ref, std::string(1, alt))});
    auto ct = [](std::string s) { for (char &b : s) if (b == 'C') b = 'T'; return s; };
    // A 70-base fragment: read 1 covers the site with REF, read 2 is the fragment's last 35 bases.
    const std::string f1 = x.c1.substr(p - 10, 35), f2 = x.c1.substr(p + 25, 35);
    auto fq = [&](const std::string &path, const std::string &read, bool trim) {
        FILE *fp = fopen(path.c_str(), "w"); REQUIRE(fp != nullptr);
        const std::string bases = read + (trim ? std::string(55, 'N') : "");
        const std::string qual = std::string(35, 'I') + (trim ? std::string(55, '!') : "");
        fprintf(fp, "@r\n%s\n+\n%s\n", bases.c_str(), qual.c_str()); REQUIRE(fclose(fp) == 0);
    };
    const std::string site = "XV:Z:c1:" + std::to_string(p + 1) + ":" + ref + ">" + alt + ",R,0,2,1,";
    for (bool meth : {false, true}) for (bool paired : {false, true}) for (bool trim : {false, true}) {
        CAPTURE(meth); CAPTURE(paired); CAPTURE(trim);
        const std::string a = dir + "/a.fq", b = dir + "/b.fq", sam = dir + "/out.sam";
        fq(a, meth ? ct(f1) : f1, trim); fq(b, revcomp(meth ? ct(f2) : f2), trim);
        std::vector<std::string> args = {"-q", "20", "-f", sam, x.prefix, a};
        if (meth) args.insert(args.begin(), "--meth");
        args.insert(args.begin(), {"--known-snvs", vcf});
        if (paired) args.push_back(b);
        REQUIRE(run(args) == 0);
        const auto rec = records(slurp(sam)); REQUIRE(rec.size() == (paired ? 2 : 1));
        INFO(slurp(sam));
        // Both reads are exact and unique. The 70-base fragment pairs at the searched length (35)
        // but not at the untrimmed one (90), which pairing's minimum outer span uses: a trimmed
        // pair stays unpaired, single-end MAPQ 37, and its ALT twin pair must pair alike.
        const int mapq = paired && !trim ? 60 : 37;
        CHECK(rec[0][3] == std::to_string(p - 9));
        CHECK(rec[0][4] == std::to_string(mapq));
        CHECK(tag(rec[0], "XV:Z:") == site + std::to_string(mapq));
    }
}

TEST_CASE("aln_known_snvs: excluded ALT spans cannot raise the insert-size floor"
          * doctest::test_suite("unit/aln_known_snvs")) {
    namespace ba = bwa_aln;
    ba::bwa_seq_t ends[2][31] = {};
    uint8_t excluded[31] = {};
    for (int i = 0; i < 31; ++i) {
        ends[0][i].pos = i * 1000;
        ends[1][i].pos = i * 1000 + 130 + i % 11;
        for (int j = 0; j < 2; ++j) ends[j][i].len = 50, ends[j][i].mapQ = 37;
    }
    ends[1][30].ref_shift = 10000;
    ends[0][30].mapQ = ends[1][30].mapQ = 0;
    excluded[30] = 1;
    ba::bwa_seq_t *seqs[2] = {ends[0], ends[1]};
    ba::isize_info_t ii = {};
    REQUIRE(ba::infer_isize(31, seqs, &ii, 1e-5, 100000, nullptr, 0, false, true, excluded) == 0);
    CHECK(ii.avg > 0);
    CHECK(ba::infer_isize(31, seqs, &ii, 1e-5, 100000, nullptr, 0, false, true) == -1);
}

TEST_CASE("aln_known_snvs: XV escapes delimiters in contig names"
          * doctest::test_suite("unit/aln_known_snvs")) {
    bntann1_t ann = {}; ann.name = const_cast<char *>("c;:,%1");
    bntseq_t bns = {}; bns.anns = &ann;
    bwa_aln::ks_site_t site = {10, 0, 1};
    bwa_aln::ks_vote_t vote = {}; vote.site = &site; vote.rid = 0; vote.allele = 'R'; vote.max_diff = 2;
    std::string xv;
    bwa_aln::ks_format_xv(&bns, {vote}, 0, 1, 37, &xv);
    CHECK(xv == "c%3B%3A%2C%251:11:A>C,R,0,2,1,37");
}

TEST_CASE("aln_known_snvs: rejected junction rows cannot produce a credited merge with no best hit"
          * doctest::test_suite("unit/aln_known_snvs")) {
    namespace ba = bwa_aln;
    KsIndex &x = fixture();
    FMI_search fmi(x.prefix.c_str());
    fmi.load_index(true, 1, false);
    ba::bwt_t *bwt = ba::bwt_from_fmi(&fmi);
    ba::ks_sites_t sites; sites.by_contig.resize(fmi.idx->bns->n_seqs);
    sites.by_contig[0] = {{499, (uint8_t)nst_nt4_table[(int)x.c1[499]], 0}};
    sites.n_sites = 1;
    ba::ks_windows_t w;
    REQUIRE(ba::ks_windows_build(&sites, fmi.idx->bns, fmi.idx->pac, 1, 0, &w) == 0);
    ba::bwt_aln1_t main = {}, win = {};
    auto row_at = [](const ba::bwt_t *b, const bntseq_t *bns, int64_t pos, int span) {
        for (ba::bwtint_t r = 0; r <= b->seq_len; ++r) {
            int strand;
            if ((int64_t)ba::bwa_sa2pos(bns, b, r, span, &strand) == pos && strand == 0) return r;
        }
        return (ba::bwtint_t)-1;
    };
    main.k = (ba::bwtint_t)-1;
    for (ba::bwtint_t r = 0; r <= bwt->seq_len; ++r)
        if (ba::bwt_sa(bwt, r) == (ba::bwtint_t)fmi.idx->bns->l_pac - 10) main.k = r;
    main.l = main.k;
    win.k = win.l = row_at(w.bwt, w.bns, 150, 20);
    REQUIRE(main.k != (ba::bwtint_t)-1); REQUIRE(win.k != (ba::bwtint_t)-1);
    win.score = 3; win.n_mm = 1;
    ba::bwa_seq_t s = {}; s.len = 20; s.n_aln = 1; s.aln = &main;
    uint8_t seq[20] = {}; s.seq = seq;
    ba::gap_opt_t *opt = ba::gap_init_opt();
    int64_t alt = -1;
    const int status = bwa_tests::run_in_child([&] {
        const auto result = ba::ks_override(&w, fmi.idx->bns, bwt, opt, 0, &s, &win, 1, &alt);
        _exit(result == ba::KS_INCOMPLETE ? 0 : 2);
    });
    REQUIRE(WIFEXITED(status)); CHECK(WEXITSTATUS(status) == 0);
    ba::ks_windows_free(&w); ba::bwt_destroy(bwt); free(opt);
}

TEST_CASE("aln_known_snvs: -f naming the --known-snvs file, by path or hard link, is refused before writing"
          * doctest::test_suite("unit/aln_known_snvs")) {
    KsIndex &x = fixture();
    bwa_tests::TempDir scratch{"ks_case"}; // this case's VCF, reads and output, removed when it ends
    const std::string dir = scratch.path();
    const std::string vcf = write_vcf(dir + "/alias.vcf", {row("c1", 500, x.c1[499], std::string(1, other(x.c1[499])))});
    const std::string link = dir + "/alias.link.vcf";
    unlink(link.c_str());
    REQUIRE(::link(vcf.c_str(), link.c_str()) == 0);
    const std::string before = slurp(vcf);
    for (const std::string &out : {vcf, link}) {
        CAPTURE(out);
        bwa_tests::StderrCapture capture;
        CHECK(run({"--known-snvs", vcf, "-f", out, x.prefix, x.fq}) == 1);
        CHECK(capture.text().find("-f '" + out + "' is the --known-snvs file") != std::string::npos);
        CHECK(slurp(vcf) == before);
    }
}

TEST_CASE("aln_known_snvs: --known-snvs and a read file cannot both be standard input"
          * doctest::test_suite("unit/aln_known_snvs")) {
    KsIndex &x = fixture();
    bwa_tests::TempDir scratch{"ks_case"}; // this case's VCF, reads and output, removed when it ends
    const std::string dir = scratch.path();
    // The sites are read to the end before the reads, which would find standard input empty.
    for (const std::vector<std::string> &reads : {std::vector<std::string>{"-"},
                                                  std::vector<std::string>{"-", x.fq},
                                                  std::vector<std::string>{x.fq, "-"}}) {
        CAPTURE(reads.size());
        std::vector<std::string> args = {"--known-snvs", "-", "-f", dir + "/stdin.sam", x.prefix};
        args.insert(args.end(), reads.begin(), reads.end());
        bwa_tests::StderrCapture capture;
        CHECK(run(args) == 1);
        CHECK(capture.text().find("can stand for only one of --known-snvs and the read files") != std::string::npos);
    }
}

TEST_CASE("aln_known_snvs: a run that ends with exit() removes its window directory"
          * doctest::test_suite("unit/aln_known_snvs")) {
    KsIndex &x = fixture();
    bwa_tests::TempDir tmp{"ks_exit"};
    const std::string root = tmp.path();
    // Use a plain fork: run_in_child's glibc ASAN on_exit handler calls _exit and can bypass
    // the window cleanup handler inherited from earlier tests. Only this child's deliberately
    // abandoned allocations are exempt from leak checking; exit() must run the real cleanup.
    fflush(nullptr);
    const pid_t child = fork();
    REQUIRE(child >= 0);
    if (child == 0) {
#if defined(BWA_TESTS_CHILD_SKIP_LEAK_CHECK)
        __lsan_disable();
#endif
        setenv("TMPDIR", root.c_str(), 1);
        FMI_search fmi(x.prefix.c_str());
        fmi.load_index(/*load_pac=*/true, 1, /*build_seed_jump=*/false);
        bwa_aln::ks_sites_t sites;
        sites.by_contig.assign(fmi.idx->bns->n_seqs, {});
        sites.by_contig[0] = {{300, 0, 1}};
        sites.n_sites = 1;
        bwa_aln::ks_windows_t w;
        if (bwa_aln::ks_windows_build(&sites, fmi.idx->bns, fmi.idx->pac, 1, 0, &w) != 0) _exit(2);
        struct stat st;
        if (w.dir.rfind(root + "/", 0) != 0 || stat(w.dir.c_str(), &st) != 0) _exit(3);
        exit(1);
    }
    int status = 0;
    REQUIRE(waitpid(child, &status, 0) == child);
    REQUIRE(WIFEXITED(status));
    REQUIRE(WEXITSTATUS(status) == 1); // 2 or 3: the windows were not built where expected
    std::vector<std::string> left;
    DIR *d = opendir(root.c_str());
    REQUIRE(d != nullptr);
    for (struct dirent *e; (e = readdir(d)) != nullptr;)
        if (strcmp(e->d_name, ".") != 0 && strcmp(e->d_name, "..") != 0) left.push_back(e->d_name);
    closedir(d);
    CHECK(left.empty());
}

TEST_CASE("aln_known_snvs: mate rescue lists the sites of every contig its range crosses"
          * doctest::test_suite("unit/aln_known_snvs")) {
    KsIndex &x = fixture();
    bntseq_t *bns = bns_restore(x.prefix.c_str());
    REQUIRE(bns != nullptr);
    REQUIRE(bns->n_seqs == 3);
    bwa_aln::ks_sites_t sites;
    sites.by_contig.assign(3, {});
    // c1 is 600 bp; c2 starts at packed 600 and c3 after c2.
    sites.by_contig[0] = {{590, 0, 1}, {599, 2, 3}};
    sites.by_contig[1] = {{10, 1, 2}, {(int64_t)bns->anns[1].len - 1, 3, 0}};
    sites.by_contig[2] = {{5, 0, 1}};
    sites.n_sites = 5;
    const bwa_aln::ks_rescue_ctx_t ctx = {&sites, bns};
    const int64_t c2 = bns->anns[1].offset, c3 = bns->anns[2].offset;
    std::vector<std::pair<int64_t, uint8_t>> got;
    bwa_aln::ks_rescue_sites(&ctx, 580, 620, &got);
    CHECK(got == std::vector<std::pair<int64_t, uint8_t>>{{590, 1}, {599, 3}, {c2 + 10, 2}});
    got.clear();
    bwa_aln::ks_rescue_sites(&ctx, 595, c3 + 6, &got); // c1's tail, all of c2, c3's head
    CHECK(got == std::vector<std::pair<int64_t, uint8_t>>{{599, 3}, {c2 + 10, 2}, {c3 - 1, 0}, {c3 + 5, 1}});
    got.clear();
    bwa_aln::ks_rescue_sites(&ctx, 591, 599, &got); // between sites
    CHECK(got.empty());
    bns_destroy(bns);
}

TEST_CASE("aln_known_snvs: --meth mate rescue lists only its anchor contig's sites"
          * doctest::test_suite("unit/aln_known_snvs")) {
    KsIndex &x = fixture();
    bntseq_t *bns = bns_restore(x.prefix.c_str());
    REQUIRE(bns != nullptr);
    REQUIRE(bns->n_seqs == 3);
    // c1's last 40 bases hold more site alleles than the cap; c2 has one site near its start.
    bwa_aln::ks_sites_t sites;
    sites.by_contig.assign(3, {});
    for (int64_t p = 560; p < 600; ++p) sites.by_contig[0].push_back({p, 0, 1});
    REQUIRE(sites.by_contig[0].size() > KS_RESCUE_MAX_SITES);
    sites.by_contig[1] = {{10, 1, 2}};
    sites.n_sites = (int64_t)sites.by_contig[0].size() + 1;
    bwa_aln::aln_meth_ref_t ref;
    ref.bns = bns;
    ref.known_snvs = &sites;
    const int64_t c2 = bns->anns[1].offset;
    // A window from c1's tail into c2's head (a reverse anchor near c2's start): c1's sites would
    // fill the cap and c2's site, the only one rescue can use, would never be tried.
    std::vector<int64_t> keys;
    bwa_aln::aln_meth_window_sites(ref, 550, (int)(c2 + 50 - 550), &keys);
    CHECK(keys.size() == (size_t)KS_RESCUE_MAX_SITES);
    CHECK(std::find(keys.begin(), keys.end(), bwa_aln::aln_site_key(c2 + 10, 2)) == keys.end());
    bwa_aln::aln_meth_window_sites(ref, 550, (int)(c2 + 50 - 550), &keys, /*rid=*/1);
    CHECK(keys == std::vector<int64_t>{bwa_aln::aln_site_key(c2 + 10, 2)});
    // and the forward case: a window from c1 running past its end lists c1's sites only
    bwa_aln::aln_meth_window_sites(ref, 595, 30, &keys, /*rid=*/0);
    CHECK(keys == std::vector<int64_t>{bwa_aln::aln_site_key(595, 1), bwa_aln::aln_site_key(596, 1), bwa_aln::aln_site_key(597, 1),
                                       bwa_aln::aln_site_key(598, 1), bwa_aln::aln_site_key(599, 1)});
    bwa_aln::aln_meth_window_sites(ref, c2 + 20, 30, &keys, /*rid=*/0); // wholly outside the contig
    CHECK(keys.empty());
    bns_destroy(bns);
}

TEST_CASE("aln_known_snvs: at a multi-allelic site, the other ALT's vote counts against the credited ALT"
          * doctest::test_suite("unit/aln_known_snvs")) {
    KsIndex &x = fixture();
    bwa_tests::TempDir scratch{"ks_case"}; // this case's VCF, reads and output, removed when it ends
    const std::string dir = scratch.path();
    // Site c1:500 with both non-REF bases a1 < a2 as ALTs. The read is c1[484, 516) carrying a2 and
    // nothing else: credited with a2 it has no differences, so every vote's `away` is 0 (the site's
    // own base is never counted), whichever ALT the vote names.
    const char ref = x.c1[499];
    std::string alts;
    for (char c : {'A', 'C', 'G', 'T'})
        if (c != ref) alts += c;
    const char a1 = alts[0], a2 = alts[1];
    std::string read = x.c1.substr(484, 32);
    read[15] = a2;
    const std::string fq = write_fq(dir + "/multi.fq", {{"m", read}});
    const std::string vcf = write_vcf(dir + "/multi.vcf", {row("c1", 500, ref, std::string(1, a1) + "," + a2)});
    const std::string out = dir + "/multi.sam";
    bwa_tests::StderrCapture capture;
    REQUIRE(run({"--known-snvs", vcf, "-f", out, x.prefix, fq}) == 0);
    const auto k = records(slurp(out));
    REQUIRE(k.size() == 1);
    // Forward at its locus, unique: NM 1 and MD against the reference, XM 0 with a2 credited.
    const expect_t e = expect_fwd("m", x.c1, 485, read);
    REQUIRE(e.nm == 1);
    const std::vector<std::string> want = {"m", "0", "c1", "485", "37", "32M", "*", "0", "0", read,
        std::string(32, 'I'), "XT:A:U", e.f[11], "X0:i:1", "X1:i:0", "XM:i:0", "XO:i:0", "XG:i:0", e.f[12],
        std::string("XV:Z:c1:500:") + ref + ">" + a1 + ",O,0,2,NA,NA;c1:500:" + ref + ">" + a2 + ",A,0,2,1,37"};
    CHECK(k[0] == want);
}

TEST_CASE("aln_known_snvs: --meth evaluates every ALT of a multi-allelic site"
          * doctest::test_suite("unit/aln_known_snvs")) {
    KsMethIndex x;
    bwa_tests::TempDir scratch{"ks_case"}; // this case's VCF, reads and output, removed when it ends
    const std::string dir = scratch.path();
    // An unmethylated read 1 (every C read as T) from c1[q - 10, q + 15), where c1[q] is A or T,
    // with ALTs C and G there. The read carries G, the second ALT, plus two A<->G mismatches: three
    // differences against the reference and with C credited, over a 25-base read's budget
    // of 2; two with G credited.
    int q = 495;
    while (q < 505 && x.c1[q] != 'A' && x.c1[q] != 'T') ++q;
    REQUIRE(q < 505);
    const char ref = x.c1[q];
    const int beg = q - 10;
    std::string read = x.c1.substr(beg, 25);
    read[10] = 'G';
    int n_mm = 0;
    for (int k = 0; k < 25 && n_mm < 2; ++k)
        if (k != 10 && k % 7 == 3 && (read[k] == 'A' || read[k] == 'G')) read[k] = read[k] == 'A' ? 'G' : 'A', ++n_mm;
    REQUIRE(n_mm == 2);
    for (char &c : read)
        if (c == 'C') c = 'T';
    const std::string fq = write_fq(dir + "/mm.fq", {{"g", read}});
    const std::string vcf = write_vcf(dir + "/mm.vcf", {row("c1", q + 1, ref, "C,G")});
    const std::string stock = dir + "/mm.stock.sam", out = dir + "/mm.sam";
    bwa_tests::StderrCapture capture;
    REQUIRE(run({"--meth", "-f", stock, x.prefix, fq}) == 0);
    REQUIRE(run({"--meth", "--known-snvs", vcf, "-f", out, x.prefix, fq}) == 0);
    const auto s = records(slurp(stock)), k = records(slurp(out));
    REQUIRE(s.size() == 1);
    REQUIRE(k.size() == 1);
    CHECK((std::stoi(s[0][1]) & 4) != 0); // stock: unmapped
    CHECK(k[0][1] == "0");
    CHECK(k[0][2] == "c1");
    CHECK(k[0][3] == std::to_string(beg + 1));
    const std::string site = "c1:" + std::to_string(q + 1) + ":" + ref + ">";
    CHECK(k[0][4] == "25"); // two mismatches with G credited, the budget (with C, or no ALT, three)
    // The C vote is another base, two differences away with G credited; the REF twin (two
    // mismatches) lands on the same locus at MAPQ 25.
    CHECK(tag(k[0], "XV:Z:") == "XV:Z:" + site + "C,O,2,2,NA,NA;" + site + "G,A,2,2,1,25");
}

TEST_CASE("aln_known_snvs: ALT-aware mate rescue keeps a trial only if it fits the fragment and is not tied elsewhere"
          * doctest::test_suite("unit/aln_known_snvs")) {
    namespace ba = bwa_aln;
    // A forward anchor at [1000, 1030): under --short-insert its reverse mate must start at or
    // after 1000 and end at or after 1030.
    ba::bwa_seq_t anchor;
    memset(&anchor, 0, sizeof anchor);
    anchor.pos = 1000, anchor.len = 30, anchor.strand = 0;
    ba::pe_opt_t *popt = ba::bwa_init_pe_opt();
    ba::bwa_cigar_t m30 = __cigar_create(FROM_M, 30);
    auto trial = [&](int64_t beg, int score) { return ba::pe_rescue_trial_t{&m30, 1, 0, beg, score, -1, 0}; };
    // a trial whose window's best score ties at two loci: no alignment, only that score
    auto tied = [&](int score) { return ba::pe_rescue_trial_t{nullptr, 0, 0, 0, 0, -1, score}; };
    for (int short_insert : {0, 1}) {
        CAPTURE(short_insert);
        popt->short_insert = short_insert;
        // the reference trial fits; a better ALT trial starts before the anchor
        std::vector<ba::pe_rescue_trial_t> t = {trial(1010, 20), trial(990, 25)};
        CHECK(ba::aln_rescue_pick(t.data(), (int)t.size(), &anchor, popt) == (short_insert ? 0 : 1));
        // nothing fits
        t = {trial(990, 20), trial(995, 25)};
        CHECK(ba::aln_rescue_pick(t.data(), (int)t.size(), &anchor, popt) == (short_insert ? -1 : 1));
    }
    popt->short_insert = 0;
    // equal best at overlapping loci: the first stands (the same placement found again)
    std::vector<ba::pe_rescue_trial_t> t = {trial(1010, 20), trial(1010, 20), trial(1015, 20)};
    CHECK(ba::aln_rescue_pick(t.data(), (int)t.size(), &anchor, popt) == 0);
    // trial order, not position, picks among them: a later trial that starts earlier does not win
    t = {trial(1015, 20), trial(1010, 20)};
    CHECK(ba::aln_rescue_pick(t.data(), (int)t.size(), &anchor, popt) == 0);
    // equal best at loci that do not overlap: ambiguous, as one window's equal best is
    t = {trial(1010, 20), trial(1100, 20)};
    CHECK(ba::aln_rescue_pick(t.data(), (int)t.size(), &anchor, popt) == -1);
    t = {trial(1010, 18), trial(1100, 20), trial(1200, 20)};
    CHECK(ba::aln_rescue_pick(t.data(), (int)t.size(), &anchor, popt) == -1);
    // Each later interval overlaps the first, but the later intervals are disjoint.
    t = {trial(100, 20), trial(75, 20), trial(125, 20)};
    CHECK(ba::aln_rescue_pick(t.data(), (int)t.size(), &anchor, popt) == -1);
    // a better trial clears an earlier tie
    t = {trial(1010, 20), trial(1100, 20), trial(1200, 22)};
    CHECK(ba::aln_rescue_pick(t.data(), (int)t.size(), &anchor, popt) == 2);
    // trials without an alignment are skipped
    t = {{nullptr, 0, 0, 0, 0, -1, 0}, trial(1100, 20)};
    CHECK(ba::aln_rescue_pick(t.data(), (int)t.size(), &anchor, popt) == 1);
    // a trial tied in its own window at the pick's score or above leaves the mate ambiguous, in
    // either order (the reference window tied and an ALT one unique, or the reverse)
    t = {tied(20), trial(1100, 20)};
    CHECK(ba::aln_rescue_pick(t.data(), (int)t.size(), &anchor, popt) == -1);
    t = {trial(1100, 20), tied(20)};
    CHECK(ba::aln_rescue_pick(t.data(), (int)t.size(), &anchor, popt) == -1);
    t = {trial(1100, 20), tied(25)};
    CHECK(ba::aln_rescue_pick(t.data(), (int)t.size(), &anchor, popt) == -1);
    // a uniquely better trial stands
    t = {tied(20), trial(1100, 21)};
    CHECK(ba::aln_rescue_pick(t.data(), (int)t.size(), &anchor, popt) == 1);
    // and a tie alone rescues nothing
    t = {tied(20)};
    CHECK(ba::aln_rescue_pick(t.data(), (int)t.size(), &anchor, popt) == -1);
    free(popt);
}

namespace {

// Two sites 5 bases apart on c1, at positions p1 and p1 + 5 (0-based) whose REF is A, G or T,
// for --meth read 1 (C->T): no conversion explains REF or ALT there.
struct TwoSites {
    int p1 = -1, p2 = -1;
    char r1 = 0, r2 = 0;
    std::string not_c(char ref) const // the bases other than ref and C, in order
    {
        std::string s;
        for (char c : {'A', 'G', 'T'})
            if (c != ref) s += c;
        return s;
    }
    explicit TwoSites(const std::string &c1)
    {
        for (int i = 460; i < 520 && p1 < 0; ++i)
            if (c1[i] != 'C' && c1[i + 5] != 'C') p1 = i;
        REQUIRE(p1 >= 0);
        p2 = p1 + 5, r1 = c1[p1], r2 = c1[p2];
    }
};

std::string conv_ct(std::string s)
{
    for (char &c : s)
        if (c == 'C') c = 'T';
    return s;
}

std::string site_id(int p0, char ref, char alt)
{
    return "c1:" + std::to_string(p0 + 1) + ":" + ref + ">" + alt;
}

} // namespace

TEST_CASE("aln_known_snvs: --meth credits one site allele at a time, the other sites at REF"
          * doctest::test_suite("unit/aln_known_snvs")) {
    KsMethIndex x;
    bwa_tests::TempDir scratch{"ks_case"}; // this case's VCF, reads and output, removed when it ends
    const std::string dir = scratch.path();
    const TwoSites t(x.c1);
    const char a1 = t.not_c(t.r1)[0], a2 = t.not_c(t.r2)[0];
    // 25-base unmethylated read 1 from c1[p1 - 10, p1 + 15): site 1 is base 10, site 2 base 15.
    const int beg = t.p1 - 10;
    auto read = [&](char b1, char b2) {
        std::string r = x.c1.substr(beg, 25);
        r[10] = b1, r[15] = b2;
        return conv_ct(r);
    };
    const std::string fq = write_fq(dir + "/two.fq", {{"rr", read(t.r1, t.r2)}, {"ar", read(a1, t.r2)},
                                                        {"ra", read(t.r1, a2)}, {"aa", read(a1, a2)}});
    const std::string vcf = write_vcf(dir + "/two.vcf", {row("c1", t.p1 + 1, t.r1, std::string(1, a1)),
                                                           row("c1", t.p2 + 1, t.r2, std::string(1, a2))});
    const std::string out = dir + "/two.sam";
    bwa_tests::StderrCapture capture;
    REQUIRE(run({"--meth", "--known-snvs", vcf, "-f", out, x.prefix, fq}) == 0);
    const auto k = records(slurp(out));
    REQUIRE(k.size() == 4);
    const std::string s1 = site_id(t.p1, t.r1, a1), s2 = site_id(t.p2, t.r2, a2);
    for (const auto &r : k) {
        CAPTURE(r[0]);
        CHECK(r[1] == "0");
        CHECK(r[2] == "c1");
        CHECK(r[3] == std::to_string(beg + 1));
        CHECK(r[4] == "37");
    }
    // REF/REF: no differences; its twin (ALT/REF) lands on the same locus at the same MAPQ.
    CHECK(tag(k[0], "XV:Z:") == "XV:Z:" + s1 + ",R,0,2,1,37;" + s2 + ",R,0,2,NA,NA");
    // ALT/REF and REF/ALT: the one ALT credited, the other site at REF, no differences.
    CHECK(tag(k[1], "XV:Z:") == "XV:Z:" + s1 + ",A,0,2,1,37;" + s2 + ",R,0,2,NA,NA");
    CHECK(tag(k[2], "XV:Z:") == "XV:Z:" + s1 + ",R,0,2,1,37;" + s2 + ",A,0,2,NA,NA");
    // ALT/ALT: credited at one site (the lower position on a tie), one difference: site 2's ALT
    // away from site 1, none away from site 2. Its twin (REF/ALT) is credited at site 2.
    CHECK(tag(k[3], "XV:Z:") == "XV:Z:" + s1 + ",A,1,2,1,37;" + s2 + ",A,0,2,NA,NA");
}

TEST_CASE("aln_known_snvs: --meth credits each ALT of a tri-allelic site with its neighbour at REF"
          * doctest::test_suite("unit/aln_known_snvs")) {
    KsMethIndex x;
    bwa_tests::TempDir scratch{"ks_case"}; // this case's VCF, reads and output, removed when it ends
    const std::string dir = scratch.path();
    const TwoSites t(x.c1);
    const std::string alts = t.not_c(t.r1);
    const char a2 = t.not_c(t.r2)[0];
    const int beg = t.p1 - 10;
    auto read = [&](char b1) {
        std::string r = x.c1.substr(beg, 25);
        r[10] = b1;
        return conv_ct(r);
    };
    const std::string fq = write_fq(dir + "/tri.fq", {{"x", read(alts[0])}, {"y", read(alts[1])}});
    const std::string vcf = write_vcf(dir + "/tri.vcf", {row("c1", t.p1 + 1, t.r1, std::string(1, alts[0]) + "," + alts[1]),
                                                           row("c1", t.p2 + 1, t.r2, std::string(1, a2))});
    const std::string out = dir + "/tri.sam";
    bwa_tests::StderrCapture capture;
    REQUIRE(run({"--meth", "--known-snvs", vcf, "-f", out, x.prefix, fq}) == 0);
    const auto k = records(slurp(out));
    REQUIRE(k.size() == 2);
    const std::string sx = site_id(t.p1, t.r1, alts[0]), sy = site_id(t.p1, t.r1, alts[1]), s2 = site_id(t.p2, t.r2, a2);
    for (const auto &r : k) {
        CAPTURE(r[0]);
        CHECK(r[3] == std::to_string(beg + 1));
        CHECK(r[4] == "37");
    }
    CHECK(tag(k[0], "XV:Z:") == "XV:Z:" + sx + ",A,0,2,1,37;" + sy + ",O,0,2,NA,NA;" + s2 + ",R,0,2,NA,NA");
    CHECK(tag(k[1], "XV:Z:") == "XV:Z:" + sx + ",O,0,2,NA,NA;" + sy + ",A,0,2,1,37;" + s2 + ",R,0,2,NA,NA");
}

TEST_CASE("aln_known_snvs: --meth mate rescue credits one site allele with the other site at REF"
          * doctest::test_suite("unit/aln_known_snvs")) {
    KsMethIndex x;
    bwa_tests::TempDir scratch{"ks_case"}; // this case's VCF, reads and output, removed when it ends
    const std::string dir = scratch.path();
    const TwoSites t(x.c1);
    const char a1 = t.not_c(t.r1)[0], a2 = t.not_c(t.r2)[0];
    const std::string fq1 = dir + "/mr1.fq", fq2 = dir + "/mr2.fq";
    FILE *f1 = fopen(fq1.c_str(), "w"), *f2 = fopen(fq2.c_str(), "w");
    REQUIRE(f1 != nullptr);
    REQUIRE(f2 != nullptr);
    auto put = [&](const std::string &name, const std::string &r1, const std::string &frag) {
        const std::string r2 = revcomp(conv_ct(frag).substr(frag.size() - 25));
        fprintf(f1, "@%s\n%s\n+\n%s\n", name.c_str(), r1.c_str(), std::string(r1.size(), 'I').c_str());
        fprintf(f2, "@%s\n%s\n+\n%s\n", name.c_str(), r2.c_str(), std::string(25, 'I').c_str());
    };
    for (int p = 400, i = 0; p + 65 <= 600; p += 5, ++i) {
        const std::string frag = x.c1.substr(p, 55 + i % 11);
        put("bg" + std::to_string(i), conv_ct(frag).substr(0, 25), frag);
    }
    // A 23-base read 1 from c1[p1 - 10, p1 + 13): ALT at site 1 (base 10), REF at site 2 (base
    // 15), and its first three bases wrong. Its own search fails (three differences with the ALT
    // credited, over the budget of 2); rescue's local alignment clips the three and scores the
    // other 20 bases, enough to rescue only against the window with site 1's ALT and site 2 REF.
    const int beg = t.p1 - 10;
    const std::string frag = x.c1.substr(beg, 60);
    std::string r1 = frag.substr(0, 23);
    r1[10] = a1;
    for (int q = 0; q < 3; ++q) r1[q] = r1[q] == 'A' ? 'T' : r1[q] == 'G' ? 'T' : r1[q] == 'T' ? 'A' : 'A';
    put("resc", conv_ct(r1), frag);
    REQUIRE(fclose(f1) == 0);
    REQUIRE(fclose(f2) == 0);
    const std::string vcf = write_vcf(dir + "/mr.vcf", {row("c1", t.p1 + 1, t.r1, std::string(1, a1)),
                                                          row("c1", t.p2 + 1, t.r2, std::string(1, a2))});
    const std::string out = dir + "/mr.sam";
    bwa_tests::StderrCapture capture;
    REQUIRE(run({"--meth", "--known-snvs", vcf, "-f", out, x.prefix, fq1, fq2}) == 0);
    const auto k = records(slurp(out));
    std::vector<std::string> names;
    for (int i = 0; i < 28; ++i) names.push_back("bg" + std::to_string(i));
    names.push_back("resc");
    require_records(slurp(out), names, true);
    const auto &r = k[k.size() - 2];
    REQUIRE(r[0] == "resc");
    CHECK((std::stoi(r[1]) & 0x16) == 0x2); // rescued: mapped, forward, in a proper pair
    CHECK(r[3] == std::to_string(beg + 4));
    CHECK(r[5] == "3S20M");
    CHECK(tag(r, "XT:A:") == "XT:A:M"); // placed by rescue, not its own search
    const std::string s1 = site_id(t.p1, t.r1, a1), s2 = site_id(t.p2, t.r2, a2);
    CHECK(tag(r, "XV:Z:").rfind("XV:Z:" + s1 + ",A,0,2,", 0) == 0);
    CHECK(tag(r, "XV:Z:").find(";" + s2 + ",R,0,2,NA,NA") != std::string::npos);
}

TEST_CASE("aln_known_snvs: --meth ALT-window candidates obey the read's difference limit"
          * doctest::test_suite("unit/aln_known_snvs")) {
    KsMethIndex x;
    bwa_tests::TempDir scratch{"ks_case"}; // this case's VCF, reads and output, removed when it ends
    const std::string dir = scratch.path();
    // An unmethylated read 1 of 25 bases from c1[beg, beg + 25): the ALT (not C) of a site at base
    // 12, whose REF is not C, plus two A<->G mismatches, so three differences against the
    // reference and two with the ALT credited (the limit is 2).
    int beg = -1, q = 12;
    std::vector<int> mm;
    for (int b = 440; b < 540 && beg < 0; ++b) {
        if (x.c1[b + q] == 'C') continue;
        mm.clear();
        for (int k : {17, 20, 22, 23, 19, 21})
            if (mm.size() < 2 && (x.c1[b + k] == 'A' || x.c1[b + k] == 'G')) mm.push_back(k);
        if (mm.size() == 2) beg = b;
    }
    REQUIRE(beg >= 0);
    const char ref = x.c1[beg + q], alt = TwoSites(x.c1).not_c(ref)[0];
    std::string r = x.c1.substr(beg, 25);
    r[q] = alt;
    for (int k : mm) r[k] = r[k] == 'A' ? 'G' : 'A';
    const std::string fq = write_fq(dir + "/lim.fq", {{"r", conv_ct(r)}});
    const std::string site = row("c1", beg + q + 1, ref, std::string(1, alt));
    // `dense` adds three ALTs at each of bases 0-10: 33 site alleles before the read's, so the
    // read's ALT is past the 32 a candidate window tries.
    std::vector<std::string> dense;
    for (int k = 0; k <= 10; ++k) {
        const char b = x.c1[beg + k];
        std::string alts;
        for (char c : {'A', 'C', 'G', 'T'})
            if (c != b) alts += alts.empty() ? std::string(1, c) : std::string(",") + c;
        dense.push_back(row("c1", beg + k + 1, b, alts));
    }
    dense.push_back(site);
    const std::string one = write_vcf(dir + "/lim.one.vcf", {site}), many = write_vcf(dir + "/lim.dense.vcf", dense);
    bwa_tests::StderrCapture capture;
    const std::string so = dir + "/lim.stock.sam", oo = dir + "/lim.one.sam", dout = dir + "/lim.dense.sam";
    REQUIRE(run({"--meth", "-f", so, x.prefix, fq}) == 0);
    REQUIRE(run({"--meth", "--known-snvs", one, "-f", oo, x.prefix, fq}) == 0);
    REQUIRE(run({"--meth", "--known-snvs", many, "-f", dout, x.prefix, fq}) == 0);
    const auto s = records(slurp(so)), k1 = records(slurp(oo)), kd = records(slurp(dout));
    REQUIRE(s.size() == 1);
    REQUIRE(k1.size() == 1);
    REQUIRE(kd.size() == 1);
    CHECK(s[0][1] == "4"); // stock: three differences, over the limit
    // One site: the ALT credited, two differences, at the limit: MAPQ 25.
    CHECK(k1[0][1] == "0");
    CHECK(k1[0][3] == std::to_string(beg + 1));
    CHECK(k1[0][4] == "25");
    const std::string sid = "c1:" + std::to_string(beg + q + 1) + ":" + ref + ">" + alt;
    CHECK(tag(k1[0], "XV:Z:") == "XV:Z:" + sid + ",A,2,2,1,25");
    // Dense: the read's ALT is never tried, and every allele that is leaves three or more
    // differences, over the limit, so the read stays unmapped as in stock.
    CHECK(kd[0] == s[0]);
}

TEST_CASE("aln_known_snvs: --meth ALT-window placements over the limits are realigned within the window's extent"
          * doctest::test_suite("unit/aln_known_snvs")) {
    bwa_tests::TempDir tmp{"ks_meth_narrow"}; const std::string dir = tmp.path();
    // The read is the window's own length, ungapped against REF GAGGTTAGGTATATAAGGGTTAGGA with
    // the ALT (T>G at its base 21) and four mismatches (ATAT read as TATA). Its best-scoring
    // realignment trades those four mismatches for an insertion and a deletion (9M1I4M1D11M).
    uint32_t state = 4242;
    const std::string core = "GAGGTTAGGTATATAAGGGTTAGGA", read = "GAGGTTAGGTTATAAAGGGTGAGGA";
    const std::string c1 = aln_test::lcg_dna(state, 300) + core + aln_test::lcg_dna(state, 300);
    const std::string prefix = tmp.unique_file();
    REQUIRE(aln_test::write_test_fasta(prefix, {{"c1", c1}}));
    REQUIRE(bwa_idx_build(prefix.c_str(), prefix.c_str(), 0, 3) == 0);
    build_meth_index(prefix);
    const std::string fq = write_fq(dir + "/r.fq", {{"r", read}});
    const std::string vcf = write_vcf(dir + "/s.vcf", {row("c1", 321, 'T', "G")});
    const std::string sam = dir + "/out.sam";
    bwa_tests::StderrCapture quiet;
    auto one = [&](std::vector<std::string> args) {
        args.insert(args.end(), {"-l", "1024", "-f", sam, prefix, fq});
        REQUIRE(run(args) == 0);
        const auto rec = records(slurp(sam)); REQUIRE(rec.size() == 1);
        return rec[0];
    };
    SUBCASE("over the limits: the ungapped ALT placement is kept") {
        // -n 4 -o 0: the gapped realignment has two gap opens; ungapped, the read has four
        // mismatches with the ALT credited, within -n 4. Stock (no sites) leaves it unmapped.
        CHECK(one({"--meth", "-n", "4", "-o", "0"})[1] == "4");
        const auto r = one({"--meth", "--known-snvs", vcf, "-n", "4", "-o", "0"});
        INFO(slurp(sam));
        CHECK(r[1] == "0"); CHECK(r[2] == "c1"); CHECK(r[3] == "301"); CHECK(r[5] == "25M");
        CHECK(r[4] == "25"); // unique, at its difference limit
        CHECK(tag(r, "NM:i:") == "NM:i:5");
        CHECK(tag(r, "MD:Z:") == "MD:Z:10A0T0A0T6T4");
        // The REF twin has the same four mismatches and no ALT: it places as the genomic search
        // places it, at the same start.
        CHECK(tag(r, "XV:Z:") == "XV:Z:c1:321:T>G,A,4,4,1,37");
    }
    SUBCASE("within the limits: the best-scoring realignment is unchanged") {
        // -n 4 -o 2: the gapped realignment (two gap opens, one mismatch) fits and is kept.
        const auto r = one({"--meth", "--known-snvs", vcf, "-n", "4", "-o", "2"});
        INFO(slurp(sam));
        CHECK(r[1] == "0"); CHECK(r[2] == "c1"); CHECK(r[3] == "301"); CHECK(r[5] == "9M1I4M1D11M");
        CHECK(r[4] == "37");
        CHECK(tag(r, "NM:i:") == "NM:i:3");
        CHECK(tag(r, "XV:Z:") == "XV:Z:c1:321:T>G,A,2,4,1,37");
    }
}

TEST_CASE("aln_known_snvs: negative gap-open limits allow ungapped methylation ALT reads only"
          * doctest::test_suite("unit/aln_known_snvs")) {
    KsMethIndex x;
    bwa_tests::TempDir tmp{"ks_negative_gap"}; const std::string dir = tmp.path();
    int beg = 440, q = 12;
    while (x.c1[beg + q] == 'C') ++beg;
    const char ref = x.c1[beg + q], alt = TwoSites(x.c1).not_c(ref)[0];
    const std::string vcf = write_vcf(dir + "/s.vcf", {row("c1", beg + q + 1, ref, std::string(1, alt))});
    std::string a = x.c1.substr(beg, 30); a[q] = alt;
    std::string gapped = a; gapped.erase(20, 1);
    const std::string fq = write_fq(dir + "/r.fq", {{"ref", conv_ct(x.c1.substr(beg, 30))},
        {"alt", conv_ct(a)}, {"gap", conv_ct(gapped)}});
    bwa_tests::StderrCapture quiet;
    std::vector<std::vector<std::string>> zero;
    for (const std::string limit : {"0", "-1"}) {
        const std::string sam = dir + "/" + limit + ".sam";
        REQUIRE(run({"--meth", "--known-snvs", vcf, "-n", "0", "-l", "1024", "-o", limit,
                     "-f", sam, x.prefix, fq}) == 0);
        const auto rec = records(slurp(sam)); REQUIRE(rec.size() == 3);
        for (int i = 0; i < 2; ++i) {
            CHECK(rec[i][3] == std::to_string(beg + 1)); CHECK(rec[i][5] == "30M");
        }
        CHECK((std::stoi(rec[2][1]) & 4) != 0);
        if (limit == "0") zero = rec; else CHECK(rec == zero);
    }
    // With one difference available the deletion can fit, but only when gaps are enabled.
    const std::string gapfq = write_fq(dir + "/gap.fq", {{"gap", conv_ct(gapped)}});
    for (const std::string limit : {"1", "0", "-1"}) {
        const std::string sam = dir + "/gap" + limit + ".sam";
        REQUIRE(run({"--meth", "--known-snvs", vcf, "-n", "1", "-l", "1024", "-o", limit,
                     "-f", sam, x.prefix, gapfq}) == 0);
        const auto rec = records(slurp(sam)); REQUIRE(rec.size() == 1);
        if (limit == "1") {
            CHECK((std::stoi(rec[0][1]) & 4) == 0); CHECK(rec[0][5].find('D') != std::string::npos);
        } else CHECK((std::stoi(rec[0][1]) & 4) != 0);
    }
}

TEST_CASE("aln_known_snvs: hidden conversion votes never select twins"
          * doctest::test_suite("unit/aln_known_snvs")) {
    namespace ba = bwa_aln;
    bntann1_t ann = {}; ann.name = const_cast<char *>("c1"); ann.len = 30;
    bntseq_t bns = {}; bns.n_seqs = 1; bns.anns = &ann; bns.l_pac = 30;
    auto *opt = ba::gap_init_opt();
    uint8_t pac[8] = {};
    for (int hyp : {0, 1}) for (bool cpg : {false, true}) for (bool context : {false, true}) {
        const uint8_t ref = hyp ? 1 : 2, alt = hyp ? 3 : 0;
        std::fill(pac, pac + 8, 0);
        auto put = [&](int pos, int base) { pac[pos >> 2] |= base << ((~pos & 3) << 1); };
        put(10, ref); put(hyp ? 11 : 9, context ? (hyp ? 2 : 1) : (hyp ? 0 : 3));
        ba::ks_sites_t sites; sites.by_contig = {{{10, ref, alt}}};
        uint8_t seq[20] = {}; seq[10] = alt;
        ba::bwa_seq_t s = {}; s.len = 20; s.type = BWA_TYPE_UNIQUE; s.seq = seq;
        std::vector<ba::ks_vote_t> votes;
        const int first = ba::ks_site_votes(&sites, &bns, opt, &s, -1, hyp, &votes, pac, cpg);
        REQUIRE(votes.size() == 1);
        const bool hidden = !cpg || context;
        CHECK(votes[0].allele == (hidden ? 'H' : 'A'));
        CHECK(first == (hidden ? -1 : 0));
        std::string xv; ba::ks_format_xv(&bns, votes, first, -1, 0, &xv);
        CHECK(xv.substr(xv.size() - 5) == "NA,NA");
    }
    free(opt);
}

TEST_CASE("aln_known_snvs: ALT gapped refinement keeps the credited allele and CIGAR"
          * doctest::test_suite("unit/aln_known_snvs")) {
    bwa_tests::TempDir tmp{"ks_indel"}; const std::string dir = tmp.path();
    const std::string ref = "ATGCATCCCATCTCGGCTCTCTGACTCCGG";
    const std::string read = "ATGCATCCCATCTCGATCTCTGACTCCGG";
    const std::string prefix = tmp.unique_file(), fq = dir + "/r.fq", sam = dir + "/out.sam";
    REQUIRE(aln_test::write_test_fasta(prefix, {{"c1", ref}}));
    REQUIRE(bwa_idx_build(prefix.c_str(), prefix.c_str(), 0, 3) == 0);
    FILE *fp = fopen(fq.c_str(), "w"); REQUIRE(fp != nullptr);
    fprintf(fp, "@r\n%s\n+\n%s\n", read.c_str(), std::string(read.size(), 'I').c_str()); REQUIRE(fclose(fp) == 0);
    const std::string vcf = write_vcf(dir + "/s.vcf", {row("c1", 16, 'G', "A")});
    bwa_tests::StderrCapture capture;
    REQUIRE(run({"--known-snvs", vcf, "-n", "2", "-f", sam, prefix, fq}) == 0);
    const auto rec = records(slurp(sam)); REQUIRE(rec.size() == 1);
    CHECK(rec[0][5] == "16M1D13M");
    CHECK(tag(rec[0], "XV:Z:").find("c1:16:G>A,A,1,2,") != std::string::npos);
    CHECK(tag(rec[0], "NM:i:") == "NM:i:2");
}

TEST_CASE("aln_known_snvs: help excludes credited records from the stock identity guarantee"
          * doctest::test_suite("unit/aln_known_snvs")) {
    FILE *fp = tmpfile(); REQUIRE(fp != nullptr);
    fflush(stdout); const int saved = dup(fileno(stdout)); REQUIRE(saved >= 0);
    REQUIRE(dup2(fileno(fp), fileno(stdout)) >= 0);
    const int status = run({"--help"}); fflush(stdout);
    REQUIRE(dup2(saved, fileno(stdout)) >= 0); close(saved);
    REQUIRE(status == 0); rewind(fp);
    std::string text; char buf[4096]; size_t n;
    while ((n = fread(buf, 1, sizeof buf, fp))) text.append(buf, n);
    fclose(fp);
    CHECK(text.find("--ref-span-insert or --known-snvs") != std::string::npos);
}

TEST_CASE("aln_known_snvs: collapsed mirror cells are free for away but not biological votes"
          * doctest::test_suite("unit/aln_known_snvs")) {
    namespace ba = bwa_aln;
    bntann1_t ann = {}; ann.name = const_cast<char *>("c1"); ann.len = 40;
    bntseq_t bns = {}; bns.n_seqs = 1; bns.anns = &ann; bns.l_pac = 40;
    auto *opt = ba::gap_init_opt(); opt->fnr = 0; opt->max_diff = 0;
    for (int hyp : {0, 1}) {
        const uint8_t ref = hyp ? 3 : 0, observed = hyp ? 1 : 2, alt = hyp ? 0 : 3;
        ba::ks_sites_t sites; sites.by_contig = {{{10, ref, alt}}};
        uint8_t seq[20] = {}; seq[10] = observed;
        ba::bwa_seq_t read = {}; read.len = 20; read.type = BWA_TYPE_UNIQUE; read.seq = seq;
        std::vector<ba::ks_vote_t> votes;
        CHECK(ba::ks_site_votes(&sites, &bns, opt, &read, -1, hyp, &votes, nullptr, 0, 1) == -1);
        REQUIRE(votes.size() == 1); CHECK(votes[0].allele == 'O');
        CHECK(votes[0].away == 0); CHECK(votes[0].max_diff == 0);
        CHECK(!(votes[0].away <= votes[0].max_diff - 1));
    }
    free(opt);
}

TEST_CASE("aln_known_snvs: work-limited plain searches keep stock records, with duplicate sharing"
          * doctest::test_suite("unit/aln_known_snvs")) {
    bwa_tests::TempDir tmp{"ks_twin_limit"}; const std::string dir = tmp.path();
    uint32_t state = 177;
    std::string ref = aln_test::lcg_dna(state, 3000);
    const std::string bases = "CAGGGGTTTTAATGAATGTCGGATCCTCAT";
    ref.replace(1000, bases.size(), bases);
    const std::string prefix = tmp.unique_file();
    FILE *fp = fopen(prefix.c_str(), "w"); REQUIRE(fp != nullptr);
    fprintf(fp, ">c1\n%s\n", ref.c_str()); REQUIRE(fclose(fp) == 0);
    REQUIRE(bwa_idx_build(prefix.c_str(), prefix.c_str(), 0, 3) == 0);
    const std::string vcf = write_vcf(dir + "/s.vcf", {row("c1", 1011, 'A', "C")});
    const std::string fq = write_fq(dir + "/r.fq", {{"r0", bases}, {"r1", bases}});
    const std::string mate = write_fq(dir + "/mate.fq", {{"r0", revcomp(ref.substr(1040, 30))},
        {"r1", revcomp(ref.substr(1040, 30))}});
    bwa_tests::StderrCapture quiet;
    // At -m 2 the reads' own searches stop: stock leaves them unmapped, and a window hit is not
    // credited against an incomplete search, so they stay unmapped with no vote. At -m 100 every
    // search completes.
    for (const std::string limit : {"2", "100"}) {
        const std::string sam = dir + "/" + limit + ".sam";
        REQUIRE(run({"--known-snvs", vcf, "-n", "2", "-o", "0", "-l", "1024", "-m", limit,
                     "-f", sam, prefix, fq}) == 0);
        const auto rec = records(slurp(sam)); REQUIRE(rec.size() == 2);
        for (const auto &r : rec) {
            if (limit == "2") {
                CHECK(r[1] == "4"); CHECK(r[3] == "0"); CHECK(tag(r, "XV:Z:").empty());
                continue;
            }
            REQUIRE(r[3] == "1001"); REQUIRE(std::stoi(r[4]) > 0);
            CHECK(tag(r, "XV:Z:") == "XV:Z:c1:1011:A>C,R,0,2,1,37");
        }
    }
    // Duplicate search results must own their intervals and retain the representative's cutoff.
    FMI_search fmi(prefix.c_str()); fmi.load_index(true, 1, false);
    bwa_aln::bwt_t *bwt = bwa_aln::bwt_from_fmi(&fmi);
    auto *opt = bwa_aln::gap_init_opt(); opt->fnr = -1; opt->max_diff = 2; opt->max_gapo = 0;
    opt->max_entries = 2; opt->seed_len = INT_MAX;
    std::string swapped = bases; swapped[10] = 'C';
    uint8_t codes[30]; for (int i = 0; i < 30; ++i) codes[i] = std::string("ACGT").find(swapped[29 - i]);
    bwa_aln::bwa_seq_t searches[2] = {};
    for (auto &s : searches) { s.len = 30; s.seq = codes; }
    bwa_aln::aln_dup_state_t dup_state;
    CHECK(bwa_aln::aln_search_batch(bwt, 2, searches, opt, 1, &dup_state) == 1);
    CHECK(searches[0].search_incomplete == 1); CHECK(searches[1].search_incomplete == 1);
    CHECK(searches[0].aln != searches[1].aln);
    for (auto &s : searches) free(s.aln);
    free(opt); bwa_aln::bwt_destroy(bwt);
    // Paired, in either order: the work-limited ends keep stock's records (no window credit).
    for (bool reverse : {false, true}) {
        const std::string sam = dir + "/paired.sam", stock = dir + "/paired.stock.sam";
        const std::string r1 = reverse ? mate : fq, r2 = reverse ? fq : mate;
        REQUIRE(run({"-n", "2", "-o", "0", "-l", "1024", "-m", "2", "-f", stock, prefix, r1, r2}) == 0);
        REQUIRE(run({"--known-snvs", vcf, "-n", "2", "-o", "0", "-l", "1024", "-m", "2", "-f", sam,
                     prefix, r1, r2}) == 0);
        require_records(slurp(sam), {"r0", "r1"}, true);
        CHECK(records(slurp(sam)) == records(slurp(stock)));
    }
}

TEST_CASE("aln_known_snvs: twin spans use the original contig and preserve valid slop"
          * doctest::test_suite("unit/aln_known_snvs")) {
    bntann1_t ann[2] = {}; ann[0].offset = 0; ann[0].len = 100; ann[1].offset = 100; ann[1].len = 100;
    bntseq_t bns = {}; bns.anns = ann; bns.n_seqs = 2; bns.l_pac = 200;
    for (int64_t shift : {-5, 0, 5}) CHECK(bwa_aln::ks_twin_same_contig(&bns, 120 + shift, 20, 120));
    CHECK_FALSE(bwa_aln::ks_twin_same_contig(&bns, 114, 20, 120));
    CHECK_FALSE(bwa_aln::ks_twin_same_contig(&bns, 126, 20, 120));
    CHECK_FALSE(bwa_aln::ks_twin_same_contig(&bns, 99, 20, 100)); // neighboring contig, within slop
    CHECK_FALSE(bwa_aln::ks_twin_same_contig(&bns, 180, 21, 180)); // deletion-derived span past end
    CHECK(bwa_aln::ks_twin_same_contig(&bns, 180, 20, 180));
    CHECK_FALSE(bwa_aln::ks_twin_same_contig(&bns, 120, 0, 120));
}

TEST_CASE("aln_known_snvs: paired twin spans honor search shifts and refined CIGARs"
          * doctest::test_suite("unit/aln_known_snvs")) {
    namespace ba = bwa_aln;
    bntann1_t ann = {}; ann.len = 100;
    bntseq_t bns = {}; bns.anns = &ann; bns.n_seqs = 1; bns.l_pac = 100;
    ba::bwa_seq_t q = {}; q.len = 25; q.type = BWA_TYPE_UNIQUE;
    q.pos = 75; q.ref_shift = 1;
    CHECK_FALSE(ba::ks_twin_read_same_contig(&bns, &q, q.pos));
    q.pos = 76; q.ref_shift = -1;
    CHECK(ba::ks_twin_read_same_contig(&bns, &q, q.pos));
    ba::bwa_seq_t mate = {}; mate.len = 25; mate.type = BWA_TYPE_UNIQUE; mate.pos = 75; mate.ref_shift = 1;
    CHECK_FALSE(ba::ks_twin_pair_same_contig(&bns, &q, &mate, q.pos));
    CHECK_FALSE(ba::ks_twin_pair_same_contig(&bns, &mate, &q, mate.pos));
    mate.pos = 76; mate.ref_shift = -1;
    CHECK(ba::ks_twin_pair_same_contig(&bns, &q, &mate, q.pos));
    CHECK(ba::ks_twin_pair_same_contig(&bns, &mate, &q, mate.pos));
    mate.type = BWA_TYPE_NO_MATCH; mate.pos = 100; mate.ref_shift = INT_MAX;
    CHECK(ba::ks_twin_pair_same_contig(&bns, &q, &mate, q.pos));
    q.pos = 75; q.ref_shift = 0;
    CHECK(ba::ks_twin_read_same_contig(&bns, &q, q.pos));
    q.ref_shift = -25;
    CHECK_FALSE(ba::ks_twin_read_same_contig(&bns, &q, q.pos));
    q.ref_shift = INT_MIN;
    CHECK_FALSE(ba::ks_twin_read_same_contig(&bns, &q, q.pos));
    q.ref_shift = INT_MAX;
    CHECK_FALSE(ba::ks_twin_read_same_contig(&bns, &q, q.pos));
    ba::bwa_cigar_t cigar = __cigar_create(0, 25); q.cigar = &cigar; q.n_cigar = 1;
    CHECK(ba::ks_twin_read_same_contig(&bns, &q, q.pos)); // CIGAR overrides stale shift
    cigar = __cigar_create(0, 26);
    CHECK_FALSE(ba::ks_twin_read_same_contig(&bns, &q, q.pos));
    cigar = __cigar_create(1, 25);
    CHECK_FALSE(ba::ks_twin_read_same_contig(&bns, &q, q.pos)); // zero reference consumption
}

TEST_CASE("aln_known_snvs: plain twins cannot certify a contig bridge"
          * doctest::test_suite("unit/aln_known_snvs")) {
    bwa_tests::TempDir tmp{"ks_twin_bridge"}; const std::string dir = tmp.path();
    const std::string prefix = tmp.unique_file();
    REQUIRE(aln_test::write_test_fasta(prefix, {{"c1", "GCTA"}, {"c2", std::string(169, 'A') + "C"}}));
    REQUIRE(bwa_idx_build(prefix.c_str(), prefix.c_str(), 0, 3) == 0);
    const std::string fq = write_fq(dir + "/r.fq", {{"r", std::string(169, 'A') + "C"}});
    const std::string vcf = write_vcf(dir + "/s.vcf", {row("c2", 170, 'C', "A")});
    const std::string sam = dir + "/out.sam";
    bwa_tests::StderrCapture quiet;
    for (bool paired : {false, true}) {
        std::vector<std::string> args{"--known-snvs", vcf, "-n", "0", "-o", "0", "-l", "1024", "-f", sam, prefix, fq};
        if (paired) args.push_back(fq);
        REQUIRE(run(args) == 0); require_records(slurp(sam), {"r"}, paired);
        const auto rec = records(slurp(sam));
        for (const auto &r : rec) {
            CHECK(r[2] == "c2"); CHECK(r[3] == "1"); CHECK(r[5] == "170M");
            const std::string xv = tag(r, "XV:Z:");
            REQUIRE(xv.find("c2:170:C>A,R,0,0,") != std::string::npos);
            const size_t last = xv.rfind(','), before = xv.rfind(',', last - 1);
            CHECK(xv.substr(before + 1, last - before - 1) == "0");
        }
    }
}

TEST_CASE("aln_known_snvs: a plain paired twin's supporting mate cannot cross a contig end"
          * doctest::test_suite("unit/aln_known_snvs")) {
    bwa_tests::TempDir tmp{"ks_twin_mate_bridge"}; const std::string dir = tmp.path();
    const std::string prefix = tmp.unique_file();
    uint32_t state = 90210;
    const std::string c1 = aln_test::lcg_dna(state, 400), c2 = aln_test::lcg_dna(state, 400);
    REQUIRE(aln_test::write_test_fasta(prefix, {{"c1", c1}, {"c2", c2}}));
    REQUIRE(bwa_idx_build(prefix.c_str(), prefix.c_str(), 0, 3) == 0);
    // Read 1 lies inside c1 over the site (with REF), so its twin (the ALT) lands at the same place
    // inside c1. Read 2, and the twin's copy of it, reads across the c1/c2 junction: the twin is
    // not supported by a placement inside one contig.
    const int site = 220;
    const std::string r1 = c1.substr(200, 50), r2 = revcomp(c1.substr(370) + c2.substr(0, 20));
    const std::string fq1 = write_fq(dir + "/r1.fq", {{"p", r1}}), fq2 = write_fq(dir + "/r2.fq", {{"p", r2}});
    const std::string vcf = write_vcf(dir + "/s.vcf", {row("c1", site + 1, c1[site], std::string(1, other(c1[site])))});
    const std::string sam = dir + "/out.sam";
    bwa_tests::StderrCapture quiet;
    REQUIRE(run({"--known-snvs", vcf, "-n", "2", "-o", "0", "-l", "1024", "-f", sam, prefix, fq1, fq2}) == 0);
    require_records(slurp(sam), {"p"}, true);
    for (const auto &r : records(slurp(sam))) {
        if (!(std::stoi(r[1]) & 64)) continue;
        CHECK(r[2] == "c1"); CHECK(r[3] == "201"); CHECK(r[5] == "50M");
        const std::string xv = tag(r, "XV:Z:");
        REQUIRE(xv.rfind("XV:Z:c1:221:" + std::string(1, c1[site]) + ">", 0) == 0);
        CHECK(xv.find(",R,0,2,") != std::string::npos);
        const size_t last = xv.rfind(','), before = xv.rfind(',', last - 1);
        CHECK(xv.substr(before + 1, last - before - 1) == "0");
    }
}

TEST_CASE("aln_known_snvs: a paired twin's mate moved by pairing is checked with its new hit's span"
          * doctest::test_suite("unit/aln_known_snvs")) {
    namespace ba = bwa_aln;
    KsIndex &x = fixture();
    FMI_search fmi(x.prefix.c_str()); fmi.load_index(true, 1, false);
    ba::bwt_t *bwt = ba::bwt_from_fmi(&fmi);
    const bntseq_t *bns = fmi.idx->bns;
    const int64_t c1_len = (int64_t)x.c1.size();
    std::vector<ba::bwtint_t> rows(x.c1.size(), (ba::bwtint_t)-1);
    for (ba::bwtint_t r = 0; r <= bwt->seq_len; ++r) {
        int strand;
        const auto p = ba::bwa_sa2pos(bns, bwt, r, 50, &strand);
        if (strand == 0 && p < rows.size()) rows[p] = r;
    }
    // The mate (end 0, no hit list) is placed at its best hit, an insertion (ref_shift -1). Its
    // ungapped second hit is c1's last 50 bases, and only that one pairs with the twin (end 1,
    // listed, reverse, inside c2), so pairing moves the mate there without updating ref_shift. As in
    // stock, the forward start comes from the hit's end and the stale span: one base to the right,
    // so the mate's ungapped record (50M) would end one base past c1.
    const int64_t far = 100, near_end = c1_len - 50;
    REQUIRE(rows[far] != (ba::bwtint_t)-1); REQUIRE(rows[near_end] != (ba::bwtint_t)-1);
    std::vector<ba::bwt_aln1_t> aln(2);
    aln[0].k = aln[0].l = rows[far]; aln[0].n_gapo = 1; aln[0].n_ins = 1; aln[0].score = 11;
    aln[1].k = aln[1].l = rows[near_end]; aln[1].n_mm = 4; aln[1].score = 12;
    ba::bwa_seq_t mate = {}, twin = {};
    mate.type = twin.type = BWA_TYPE_UNIQUE;
    mate.full_len = mate.clip_len = mate.len = twin.full_len = twin.clip_len = twin.len = 50;
    mate.aln = aln.data(); mate.n_aln = 2; mate.c1 = 1; mate.c2 = 1;
    mate.pos = far; mate.ref_shift = -1; mate.n_gapo = 1; mate.score = 11;
    mate.mapQ = mate.seQ = twin.mapQ = twin.seQ = 37;
    twin.pos = c1_len + 10; twin.strand = 1; twin.c1 = 1;
    ba::bwa_seq_t *seqs[2] = {&mate, &twin};
    std::vector<ba::pe_hit_t> h[2] = {{}, {{c1_len + 10, 1, 0, 0, 0, 0, 0}}};
    const std::vector<ba::pe_hit_t> *hits[2] = {h, h + 1};
    auto *pe = ba::bwa_init_pe_opt(); auto *gap = ba::gap_init_opt();
    pe->n_multi = pe->N_multi = 0;
    ba::isize_info_t ii = {}; ii.avg = -1;
    int32_t moved[2] = {-2, -2};
    int32_t *moved_span[2] = {moved, moved + 1};
    ba::bwa_pe_pair(bns, bwt, 1, seqs, &ii, pe, gap, hits, moved_span);
    REQUIRE((int64_t)mate.pos == near_end + 1);
    CHECK(mate.ref_shift == -1); // stock leaves the previous hit's
    CHECK(moved[0] == 50);
    CHECK(moved[1] == -1);
    // With the stale shift the mate's span (49) fits inside c1; its hit's span (50) does not.
    CHECK(ba::ks_twin_pair_same_contig(bns, &twin, &mate, twin.pos));
    CHECK_FALSE(ba::ks_twin_pair_same_contig(bns, &twin, &mate, twin.pos, moved[1], moved[0]));
    // An end with a CIGAR (a rescued mate) is checked with its CIGAR, whatever pairing recorded.
    ba::bwa_cigar_t cigar = __cigar_create(FROM_M, 49); mate.cigar = &cigar; mate.n_cigar = 1;
    CHECK(ba::ks_twin_pair_same_contig(bns, &twin, &mate, twin.pos, moved[1], moved[0]));
    ba::bwt_destroy(bwt);
    free(pe); free(gap);
}

TEST_CASE("aln_known_snvs: negative plain twin placements reject moved tied and strand-inverted hits"
          * doctest::test_suite("unit/aln_known_snvs")) {
    namespace ba = bwa_aln;
    KsIndex &x = fixture(); FMI_search fmi(x.prefix.c_str()); fmi.load_index(true, 1, false);
    ba::bwt_t *bwt = ba::bwt_from_fmi(&fmi);
    auto *opt = ba::gap_init_opt(); ba::bwt_aln1_t hit = {};
    int strand = 0; ba::bwtint_t pos = (ba::bwtint_t)-1;
    for (ba::bwtint_t row = 0; row < bwt->seq_len; ++row) {
        pos = ba::bwa_sa2pos(fmi.idx->bns, bwt, row, 25, &strand);
        if (pos != (ba::bwtint_t)-1 && pos > 10) { hit.k = hit.l = row; break; }
    }
    REQUIRE(pos != (ba::bwtint_t)-1);
    ba::bwa_seq_t read = {}; read.len = 25; read.n_aln = 1; read.aln = &hit;
    ba::ks_windows_t windows; int mapq = -1;
    CHECK(ba::ks_twin_place(&windows, fmi.idx->bns, bwt, opt, &read, nullptr, 0, pos + 6, strand, &mapq) == 0);
    CHECK(mapq > 0);
    CHECK(ba::ks_twin_place(&windows, fmi.idx->bns, bwt, opt, &read, nullptr, 0, pos + 5, strand, &mapq) == 1);
    CHECK(ba::ks_twin_place(&windows, fmi.idx->bns, bwt, opt, &read, nullptr, 0, pos, 1 - strand, &mapq) == 0);
    hit.l = hit.k + 1;
    CHECK(ba::ks_twin_place(&windows, fmi.idx->bns, bwt, opt, &read, nullptr, 0, pos, strand, &mapq) == 0);
    ba::ks_site_t site = {10, 0, 1}; ba::ks_vote_t vote = {}; vote.rid = 0; vote.site = &site; vote.allele = 'A';
    std::string xv; ba::ks_format_xv(fmi.idx->bns, {vote}, 0, 0, 37, &xv);
    CHECK(xv.substr(xv.size() - 5) == ",0,37");
    free(opt); ba::bwt_destroy(bwt);
}

TEST_CASE("aln_known_snvs: credited twin checks the deletion-derived full span"
          * doctest::test_suite("unit/aln_known_snvs")) {
    namespace ba = bwa_aln;
    KsIndex &x = fixture(); FMI_search fmi(x.prefix.c_str()); fmi.load_index(true, 1, false);
    ba::bwt_t *bwt = ba::bwt_from_fmi(&fmi);
    ba::bwt_aln1_t hit = {}; bool found = false;
    for (ba::bwtint_t row = 0; row < bwt->seq_len; ++row) {
        int strand;
        if (ba::bwa_sa2pos(fmi.idx->bns, bwt, row, 26, &strand) == 400 && strand == 0) {
            hit.k = hit.l = row; found = true; break;
        }
    }
    REQUIRE(found); hit.n_del = 1; hit.n_gapo = 1; hit.score = 11;
    ba::ks_windows_t windows; windows.bns = fmi.idx->bns; windows.bwt = bwt;
    windows.win.resize(windows.bns->n_seqs);
    windows.win[0].rid = 0; windows.win[0].len = x.c1.size(); windows.win[0].site_pos = 410;
    ba::bwa_seq_t read = {}; read.len = 25;
    auto *opt = ba::gap_init_opt(); opt->fnr = -1; opt->max_diff = 1;
    int mapq = -1;
    CHECK(ba::ks_twin_place(&windows, fmi.idx->bns, bwt, opt, &read, &hit, 1, 400, 0, &mapq) == 1);
    CHECK(mapq > 0);
    // Keep the read's 25 bases inside the contig, but its one deleted base outside it.
    bntann1_t ann = fmi.idx->bns->anns[0]; ann.len = 425;
    bntseq_t shorter = *fmi.idx->bns; shorter.anns = &ann; shorter.n_seqs = 1;
    shorter.l_pac = ann.len; shorter.pos2rid_bucket = nullptr;
    CHECK(ba::ks_twin_place(&windows, &shorter, bwt, opt, &read, &hit, 1, 400, 0, &mapq) == 0);
    CHECK(mapq > 0);
    free(opt); ba::bwt_destroy(bwt);
}

TEST_CASE("aln_known_snvs: stock twin placement checks the deletion-derived full span"
          * doctest::test_suite("unit/aln_known_snvs")) {
    namespace ba = bwa_aln;
    KsIndex &x = fixture(); FMI_search fmi(x.prefix.c_str()); fmi.load_index(true, 1, false);
    ba::bwt_t *bwt = ba::bwt_from_fmi(&fmi);
    ba::bwt_aln1_t hit = {}; bool found = false;
    for (ba::bwtint_t row = 0; row < bwt->seq_len; ++row) {
        int strand;
        if (ba::bwa_sa2pos(fmi.idx->bns, bwt, row, 26, &strand) == 400 && strand == 0) {
            hit.k = hit.l = row; found = true; break;
        }
    }
    REQUIRE(found); hit.n_del = 1; hit.n_gapo = 1; hit.score = 11;
    // No window hits: the read's own search result (one row with one deletion) is the twin's.
    ba::bwa_seq_t read = {}; read.len = 25; read.n_aln = 1; read.aln = &hit;
    ba::ks_windows_t windows;
    auto *opt = ba::gap_init_opt(); opt->fnr = -1; opt->max_diff = 1;
    int mapq = -1;
    CHECK(ba::ks_twin_place(&windows, fmi.idx->bns, bwt, opt, &read, nullptr, 0, 400, 0, &mapq) == 1);
    CHECK(mapq > 0);
    // c1 ends after the read's 25 bases but before its deleted base (l_pac unchanged, so the
    // hit's position is as before).
    bntann1_t ann[2] = {fmi.idx->bns->anns[0], fmi.idx->bns->anns[0]};
    ann[0].len = 425; ann[1].offset = 425; ann[1].len = fmi.idx->bns->l_pac - 425;
    bntseq_t split = *fmi.idx->bns; split.anns = ann; split.n_seqs = 2; split.pos2rid_bucket = nullptr;
    mapq = -1;
    CHECK(ba::ks_twin_place(&windows, &split, bwt, opt, &read, nullptr, 0, 400, 0, &mapq) == 0);
    CHECK(mapq > 0);
    free(opt); ba::bwt_destroy(bwt);
}

TEST_CASE("aln_known_snvs: additional window searches honor duplicate sharing in every mode"
          * doctest::test_suite("unit/aln_known_snvs")) {
    KsMethIndex x;
    bwa_tests::TempDir tmp{"ks_dup"}; const std::string dir = tmp.path();
    const int site = 480;
    const char allele = other(x.c1[site], 'C');
    const std::string vcf = write_vcf(dir + "/s.vcf", {row("c1", site + 1, x.c1[site], std::string(1, allele))});
    const std::string a = dir + "/a.fq", b = dir + "/b.fq";
    for (int j = 0; j < 2; ++j) {
        FILE *fp = fopen((j ? b : a).c_str(), "w"); REQUIRE(fp != nullptr);
        std::string seq = j ? revcomp(x.c1.substr(510, 35)) : x.c1.substr(470, 35);
        if (!j) seq[10] = allele;
        for (int k = 0; k < 12; ++k) fprintf(fp, "@r%d\n%s\n+\n%s\n", k, seq.c_str(), std::string(35, 'I').c_str());
        REQUIRE(fclose(fp) == 0);
    }
    for (bool meth : {false, true}) for (bool paired : {false, true}) {
        std::string want;
        for (bool share : {false, true}) {
            const std::string sam = dir + "/out.sam";
            namespace ba = bwa_aln;
            ba::aln_opt_t o; ba::aln_opt_init(&o); o.share_dup_search = share;
            o.gap.fnr = -1; o.gap.max_diff = 0;
            ba::aln_meth_opt_t mo; mo.m.meth_mode = meth; ba::aln_meth_opt_finish(&mo);
            FMI_search fmi((meth ? x.prefix + ".meth" : x.prefix).c_str()); fmi.load_index(!meth, 1, false);
            ba::bwt_t *bwt = ba::bwt_from_fmi(&fmi); ba::aln_meth_ref_t mr;
            if (meth) { REQUIRE(ba::aln_meth_ref_load(&mr, x.prefix.c_str(), fmi.idx->bns, 1) == 0); o.meth = &mo; o.meth_ref = &mr; }
            const bntseq_t *bns = meth ? mr.bns : fmi.idx->bns;
            const uint8_t *pac = meth ? mr.pac : fmi.idx->pac;
            ba::ks_sites_t sites; REQUIRE(ba::ks_load(vcf.c_str(), bns, pac, &sites) == 0);
            ba::ks_windows_t windows; REQUIRE(ba::ks_windows_build(&sites, bns, pac, 1, meth, &windows) == 0);
            o.known_snvs = &sites; o.ks_windows = &windows;
            if (meth) mr.known_snvs = &sites;
            ba::aln_out_t *out = ba::aln_out_open(sam.c_str(), 0, 0, "", 1); REQUIRE(out != nullptr);
            const int result = meth ? (paired ? ba::aln_run_pe_meth(&o, bwt, a.c_str(), b.c_str(), out) : ba::aln_run_se_meth(&o, bwt, a.c_str(), out))
                                    : (paired ? ba::aln_run_pe(&o, bwt, bns, fmi.idx->pac, a.c_str(), b.c_str(), out) : ba::aln_run_se(&o, bwt, bns, fmi.idx->pac, a.c_str(), out));
            REQUIRE(result == 0); REQUIRE(ba::aln_out_close(out) == 0);
            ba::ks_windows_free(&windows); ba::bwt_destroy(bwt); if (meth) ba::aln_meth_ref_free(&mr);
            std::vector<std::string> names;
            for (int i = 0; i < 12; ++i) names.push_back("r" + std::to_string(i));
            require_records(slurp(sam), names, paired);
            const auto rec = records(slurp(sam));
            for (int i = 0; i < 12; ++i) {
                const auto &r = rec[paired ? i * 2 : i];
                CHECK((std::stoi(r[1]) & 4) == 0); CHECK(r[2] == "c1"); CHECK(r[3] == "471");
                CHECK(r[5] == "35M");
                const std::string sid = "c1:481:" + std::string(1, x.c1[site]) + ">" + allele;
                CHECK(tag(r, "XV:Z:").find("XV:Z:" + sid + ",A,0,0,") == 0);
            }
            if (share) CHECK(body(slurp(sam)) == want); else want = body(slurp(sam));
        }
    }
}

TEST_CASE("aln_known_snvs: exhaustive merged hits retain higher-difference alternatives"
          * doctest::test_suite("unit/aln_known_snvs")) {
    bwa_tests::TempDir tmp{"ks_nonstop"}; const std::string dir = tmp.path();
    uint32_t state = 773; const std::string read = aln_test::lcg_dna(state, 30);
    std::vector<std::pair<std::string, std::string>> contigs;
    for (int j = 0; j < 3; ++j) {
        std::string seq = read;
        if (j == 1) for (int p : {10, 20}) seq[p] = other(seq[p]);
        if (j == 2) for (int p : {7, 15, 23}) seq[p] = other(seq[p]);
        contigs.push_back({"c" + std::to_string(j + 1), aln_test::lcg_dna(state, 100) + seq + aln_test::lcg_dna(state, 100)});
    }
    const std::string prefix = tmp.unique_file(), fq = dir + "/r.fq", sam = dir + "/out.sam";
    REQUIRE(aln_test::write_test_fasta(prefix, contigs)); REQUIRE(bwa_idx_build(prefix.c_str(), prefix.c_str(), 0, 3) == 0);
    FILE *fp = fopen(fq.c_str(), "w"); REQUIRE(fp != nullptr); fprintf(fp, "@r\n%s\n+\n%s\n", read.c_str(), std::string(30, 'I').c_str()); REQUIRE(fclose(fp) == 0);
    const auto vcf = write_vcf(dir + "/s.vcf", {row("c2", 111, contigs[1].second[110], std::string(1, read[10]))});
    bwa_tests::StderrCapture capture;
    REQUIRE(run({"--known-snvs", vcf, "-N", "-n", "4", "-o", "0", "--max-xa-hits", "20", "-l", "1024", "-f", sam, prefix, fq}) == 0);
    const auto rec = records(slurp(sam)); REQUIRE(rec.size() == 1);
    INFO(slurp(sam));
    // The primary is the exact copy on c1; the two alternatives are exactly the fixture's other
    // copies: c2 (two mismatches, one of them the credited ALT) and c3 (three mismatches), each
    // forward at 101 with the edit distance against the original reference.
    CHECK(rec[0][2] == "c1"); CHECK(rec[0][3] == "101"); CHECK(rec[0][5] == "30M");
    CHECK(tag(rec[0], "X0:i:") == "X0:i:1");
    CHECK(tag(rec[0], "X1:i:") == "X1:i:2");
    const std::string xa = tag(rec[0], "XA:Z:");
    REQUIRE(xa.rfind("XA:Z:", 0) == 0);
    std::multiset<std::string> alts;
    for (size_t a = 5, b; a < xa.size(); a = b + 1) {
        b = xa.find(';', a);
        REQUIRE(b != std::string::npos);
        alts.insert(xa.substr(a, b - a));
    }
    CHECK(alts == std::multiset<std::string>{"c2,+101,30M,2", "c3,+101,30M,3"});
}

TEST_CASE("aln_known_snvs: extra window copies cannot hide a genomic alternative"
          * doctest::test_suite("unit/aln_known_snvs")) {
    bwa_tests::TempDir tmp{"ks_copies"}; const std::string dir = tmp.path();
    uint32_t state = 478; const std::string read = aln_test::lcg_dna(state, 30);
    std::string a = aln_test::lcg_dna(state, 250) + read + aln_test::lcg_dna(state, 250);
    std::string b = a; b[260] = other(read[10]); b[270] = other(read[20]);
    const std::string prefix = tmp.unique_file(), fq = dir + "/r.fq", sam = dir + "/out.sam";
    REQUIRE(aln_test::write_test_fasta(prefix, {{"c1", a}, {"c2", b}})); REQUIRE(bwa_idx_build(prefix.c_str(), prefix.c_str(), 0, 3) == 0);
    FILE *fp = fopen(fq.c_str(), "w"); REQUIRE(fp != nullptr); fprintf(fp, "@r\n%s\n+\n%s\n", read.c_str(), std::string(30, 'I').c_str()); REQUIRE(fclose(fp) == 0);
    std::vector<std::string> rows = {row("c2", 261, b[260], std::string(1, read[10]))};
    std::string want;
    for (bool many : {false, true}) {
        if (many) for (int p = 150; p < 190; ++p) rows.push_back(row("c1", p + 1, a[p], std::string(1, other(a[p]))));
        const auto vcf = write_vcf(dir + "/s.vcf", rows);
        bwa_tests::StderrCapture capture;
        REQUIRE(run({"--known-snvs", vcf, "-n", "2", "-l", "1024", "-f", sam, prefix, fq}) == 0);
        const auto rec = records(slurp(sam)); REQUIRE(rec.size() == 1);
        CHECK(tag(rec[0], "X1:i:") == "X1:i:1");
        // the c2 copy (one difference once its site is credited) is the one alternative
        const std::string xa = tag(rec[0], "XA:Z:");
        CHECK(xa.find("XA:Z:c2,+251,30M,") == 0);
        if (!many) want = xa; else CHECK(xa == want);
    }
}

TEST_CASE("aln_known_snvs: full-read ALT-window discovery cannot discover over-budget long meth reads"
          * doctest::test_suite("unit/aln_known_snvs")) {
    bwa_tests::TempDir tmp{"ks_long_discovery"}; const std::string dir = tmp.path();
    uint32_t state = 315;
    std::string ref = aln_test::lcg_dna(state, 1000); ref[500] = 'T';
    std::string altref = ref; altref[500] = 'A';
    const std::string fq = write_fq(dir + "/r.fq", {{"r", conv_ct(altref.substr(330, 350))}});
    const std::string shortfq = write_fq(dir + "/short.fq", {{"r", conv_ct(altref.substr(490, 35))}});
    const std::string vcf = write_vcf(dir + "/s.vcf", {row("c1", 501, 'T', "A")});
    std::string prefixes[2] = {tmp.unique_file(), tmp.unique_file(".alt.fa")};
    bwa_tests::StderrCapture quiet;
    for (int i = 0; i < 2; ++i) {
        REQUIRE(aln_test::write_test_fasta(prefixes[i], {{"c1", i ? altref : ref}}));
        build_meth_index(prefixes[i]);
    }
    const std::string sam = dir + "/out.sam";
    auto align = [&](int budget, bool explicit_alt, bool short_read) {
        std::vector<std::string> args{"--meth", "-n", std::to_string(budget), "-o", "0", "-l", "1024", "-f", sam};
        if (!explicit_alt) args.insert(args.end(), {"--known-snvs", vcf});
        args.insert(args.end(), {prefixes[explicit_alt], short_read ? shortfq : fq});
        REQUIRE(run(args) == 0); const auto rec = records(slurp(sam)); REQUIRE(rec.size() == 1); return rec[0];
    };
    CHECK((std::stoi(align(0, false, false)[1]) & 4) != 0);
    const auto oracle = align(0, true, false);
    CHECK(oracle[3] == "331"); CHECK(oracle[5] == "350M"); CHECK(std::stoi(oracle[4]) > 0);
    const auto genomic = align(1, false, false);
    CHECK(genomic[3] == oracle[3]); CHECK(genomic[5] == oracle[5]);
    CHECK(tag(genomic, "XV:Z:").find("c1:501:T>A,A,0,1,") != std::string::npos);
    const auto shortrec = align(0, false, true);
    CHECK(shortrec[3] == "491"); CHECK(shortrec[5] == "35M");
}

TEST_CASE("aln_known_snvs: long meth reads retain a credited second locus beyond the window length"
          * doctest::test_suite("unit/aln_known_snvs")) {
    bwa_tests::TempDir tmp{"ks_long"}; const std::string dir = tmp.path();
    uint32_t state = 110; const std::string read = aln_test::lcg_dna(state, 350);
    std::string second = read; std::vector<int> changes;
    for (int p = 30; p < 300 && changes.size() < 2; ++p) if (read[p] != 'C') { changes.push_back(p); second[p] = other(read[p], 'C'); }
    const std::string prefix = tmp.unique_file(), fq = dir + "/r.fq", sam = dir + "/out.sam";
    const std::string flank = aln_test::lcg_dna(state, 100);
    REQUIRE(aln_test::write_test_fasta(prefix, {{"c1", flank + read + flank}, {"c2", flank + second + flank}}));
    std::string converted = read; for (char &c : converted) if (c == 'C') c = 'T';
    FILE *fp = fopen(fq.c_str(), "w"); REQUIRE(fp != nullptr); fprintf(fp, "@r\n%s\n+\n%s\n", converted.c_str(), std::string(350, 'I').c_str()); REQUIRE(fclose(fp) == 0);
    build_meth_index(prefix);
    const auto vcf = write_vcf(dir + "/s.vcf", {row("c2", 101 + changes[0], second[changes[0]], std::string(1, read[changes[0]]))});
    bwa_tests::StderrCapture capture;
    REQUIRE(run({"--meth", "--known-snvs", vcf, "-n", "2", "-o", "0", "-l", "1024", "-f", sam, prefix, fq}) == 0);
    const auto rec = records(slurp(sam)); REQUIRE(rec.size() == 1);
    CHECK(tag(rec[0], "X1:i:") == "X1:i:1");
    INFO(slurp(sam));
    CHECK(tag(rec[0], "XA:Z:").find("c2,+101,350M,") != std::string::npos);
}

TEST_CASE("aln_known_snvs: incomplete meth window discovery cannot become a confident primary"
          * doctest::test_suite("unit/aln_known_snvs")) {
    namespace ba = bwa_aln;
    bwa_tests::TempDir tmp{"ks_cap"}; const std::string dir = tmp.path();
    uint32_t state = 119; const std::string unit = aln_test::lcg_dna(state, 100);
    std::string reference; std::vector<std::string> rows;
    for (int i = 0; i < 8; ++i) { reference += unit; rows.push_back(row("c1", 100 * i + 41, unit[40], std::string(1, other(unit[40])))); }
    const std::string prefix = tmp.unique_file(), fq = dir + "/r.fq", sam = dir + "/out.sam";
    REQUIRE(aln_test::write_test_fasta(prefix, {{"c1", reference}}));
    build_meth_index(prefix);
    std::string read = unit.substr(30, 25); read[10] = other(unit[40]);
    for (char &c : read) if (c == 'C') c = 'T';
    FILE *fp = fopen(fq.c_str(), "w"); REQUIRE(fp != nullptr); fprintf(fp, "@r\n%s\n+\n%s\n", read.c_str(), std::string(25, 'I').c_str()); REQUIRE(fclose(fp) == 0);
    FMI_search fmi((prefix + ".meth").c_str()); fmi.load_index(false, 1, false);
    ba::bwt_t *bwt = ba::bwt_from_fmi(&fmi); ba::aln_meth_ref_t mr;
    REQUIRE(ba::aln_meth_ref_load(&mr, prefix.c_str(), fmi.idx->bns, 1) == 0);
    ba::ks_sites_t sites; const auto vcf = write_vcf(dir + "/s.vcf", rows);
    REQUIRE(ba::ks_load(vcf.c_str(), mr.bns, mr.pac, &sites) == 0);
    ba::ks_windows_t windows; REQUIRE(ba::ks_windows_build(&sites, mr.bns, mr.pac, 1, 1, &windows) == 0);
    mr.known_snvs = &sites;
    ba::aln_opt_t o; ba::aln_opt_init(&o); o.gap.fnr = -1; o.gap.max_diff = 0;
    ba::aln_meth_opt_t mo; mo.m.meth_mode = 1; mo.p_best = 1; mo.p_second = 0; ba::aln_meth_opt_finish(&mo);
    o.meth = &mo; o.meth_ref = &mr; o.known_snvs = &sites; o.ks_windows = &windows;
    for (bool paired : {false, true}) {
        ba::aln_out_t *out = ba::aln_out_open(sam.c_str(), 0, 0, "", 1); REQUIRE(out != nullptr);
        REQUIRE((paired ? ba::aln_run_pe_meth(&o, bwt, fq.c_str(), fq.c_str(), out) : ba::aln_run_se_meth(&o, bwt, fq.c_str(), out)) == 0);
        REQUIRE(ba::aln_out_close(out) == 0);
        const auto rec = records(slurp(sam)); REQUIRE(rec.size() == (paired ? 2 : 1));
        for (size_t end = 0; end < rec.size(); ++end) {
            const auto &r = rec[end]; INFO(slurp(sam));
            if (end == 0) { // CT end has an exact retained ALT hit; GA mate has no such hit.
                CHECK((std::stoi(r[1]) & 4) == 0); CHECK(r[2] == "c1"); CHECK(r[5] == "25M");
                const int pos = std::stoi(r[3]);
                CHECK(pos >= 31); CHECK(pos <= 731); CHECK((pos - 31) % 100 == 0);
                CHECK(tag(r, "XV:Z:").find("c1:" + std::to_string(pos + 10) + ":") != std::string::npos);
            } else CHECK((std::stoi(r[1]) & 4) != 0);
            CHECK(r[4] == "0"); CHECK(tag(r, "X0:i:").empty());
            CHECK(tag(r, "X1:i:").empty()); CHECK(tag(r, "XA:Z:").empty());
        }
    }
    ba::ks_windows_free(&windows); ba::aln_meth_ref_free(&mr); ba::bwt_destroy(bwt);
}

TEST_CASE("aln_known_snvs: single-end twins tied elsewhere are not same-locus votes"
          * doctest::test_suite("unit/aln_known_snvs")) {
    bwa_tests::TempDir tmp{"ks_twin_tie"}; const std::string dir = tmp.path();
    uint32_t state = 961; const std::string a = aln_test::lcg_dna(state, 500);
    std::string b = a; b[210] = other(a[210], 'C');
    const std::string prefix = tmp.unique_file(), fq = dir + "/r.fq", sam = dir + "/out.sam";
    REQUIRE(aln_test::write_test_fasta(prefix, {{"c1", a}, {"c2", b}}));
    build_meth_index(prefix);
    const auto vcf = write_vcf(dir + "/s.vcf", {row("c1", 211, a[210], std::string(1, b[210]))});
    for (bool meth : {false, true}) {
        std::string read = a.substr(200, 35); if (meth) for (char &c : read) if (c == 'C') c = 'T';
        FILE *fp = fopen(fq.c_str(), "w"); REQUIRE(fp != nullptr); fprintf(fp, "@r\n%s\n+\n%s\n", read.c_str(), std::string(35, 'I').c_str()); REQUIRE(fclose(fp) == 0);
        std::vector<std::string> args = {"--known-snvs", vcf, "-f", sam}; if (meth) args.push_back("--meth"); args.push_back(prefix); args.push_back(fq);
        bwa_tests::StderrCapture capture; REQUIRE(run(args) == 0);
        const auto rec = records(slurp(sam)); REQUIRE(rec.size() == 1);
        INFO(slurp(sam)); CHECK(tag(rec[0], "XV:Z:").find(",R,0,2,0,0") != std::string::npos);
    }
}

TEST_CASE("aln_known_snvs: artificial cross-window matches cannot hide a valid placement"
          * doctest::test_suite("unit/aln_known_snvs")) {
    bwa_tests::TempDir tmp{"ks_junction"}; const std::string dir = tmp.path();
    uint32_t state = 893; const std::string read = aln_test::lcg_dna(state, 30);
    std::vector<std::pair<std::string, std::string>> contigs;
    for (int i = 0; i < 3; ++i) contigs.push_back({"c" + std::to_string(i + 1), aln_test::lcg_dna(state, 700)});
    // Windows centred at200 end at361; the concatenation is read[0:15]+read[15:30].
    contigs[0].second.replace(346, 15, read.substr(0, 15));
    contigs[1].second.replace(40, 15, read.substr(15));
    contigs[2].second.replace(200, 30, read);
    for (int p : {3, 10, 19}) contigs[2].second[200 + p] = other(read[p]);
    const std::string prefix = tmp.unique_file(), fq = dir + "/r.fq", sam = dir + "/out.sam";
    REQUIRE(aln_test::write_test_fasta(prefix, contigs)); REQUIRE(bwa_idx_build(prefix.c_str(), prefix.c_str(), 0, 3) == 0);
    FILE *fp = fopen(fq.c_str(), "w"); REQUIRE(fp != nullptr); fprintf(fp, "@r\n%s\n+\n%s\n", read.c_str(), std::string(30, 'I').c_str()); REQUIRE(fclose(fp) == 0);
    std::vector<std::string> rows = {row("c3", 211, contigs[2].second[210], std::string(1, read[10]))};
    std::string want;
    for (bool joined : {false, true}) {
        if (joined) for (int j = 0; j < 2; ++j) rows.push_back(row(contigs[j].first.c_str(), 201, contigs[j].second[200], std::string(1, other(contigs[j].second[200]))));
        const auto vcf = write_vcf(dir + "/s.vcf", rows);
        bwa_tests::StderrCapture capture;
        REQUIRE(run({"--known-snvs", vcf, "-n", "2", "-o", "0", "-l", "1024", "-f", sam, prefix, fq}) == 0);
        const auto rec = records(slurp(sam)); REQUIRE(rec.size() == 1);
        CHECK(rec[0][2] == "c3"); CHECK(rec[0][3] == "201");
        if (!joined) want = body(slurp(sam)); else CHECK(body(slurp(sam)) == want);
    }
}
