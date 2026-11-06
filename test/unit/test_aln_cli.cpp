// Command-line checks of `bwa-mem3 aln` (main_aln in src/aln/aln_main.cpp): values that are
// meaningless and on which stock crashes, aborts or silently drops every read are refused with an
// error and exit 1, as are malformed numbers; inputs that would lose data (standard input for both
// files, -f naming an input) are refused before any output exists. None of this changes the
// output of a run stock completes; aln_parity.sh holds that output to stock's.

#include "doctest/doctest.h"

#include "aln_test_ref.h"
#include "../framework/stderr_capture.h"
#include "../../src/aln/aln_main.h"
#include "../../src/aln/aln_meth.h"
#include "../../src/compat_target.h"
#include "../../src/bwa.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <initializer_list>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <vector>
#include <unistd.h>

namespace {

struct CliIndex {
    std::string dir, prefix, fq;

    CliIndex()
    {
        std::string tmpl_s = aln_test::tmp_template("aln_cli_unit");
        char *tmpl = &tmpl_s[0];
        REQUIRE(mkdtemp(tmpl) != nullptr);
        dir = tmpl;
        prefix = dir + "/ref.fa";
        REQUIRE(aln_test::write_test_fasta(prefix));
        REQUIRE(bwa_idx_build(prefix.c_str(), prefix.c_str(), 0, 3) == 0);
        const std::string c1 = aln_test::test_contigs()[0].second;
        fq = dir + "/r.fq";
        FILE *fp = fopen(fq.c_str(), "w");
        REQUIRE(fp != nullptr);
        for (int i = 0; i < 4; ++i)
            fprintf(fp, "@r%d\n%s\n+\n%s\n", i, c1.substr(420 + 20 * i, 36).c_str(), std::string(36, 'I').c_str());
        REQUIRE(fclose(fp) == 0);
    }
    ~CliIndex()
    {
        aln_test::remove_dir(dir);
    }
};

CliIndex &fixture()
{
    static CliIndex idx;
    return idx;
}

int run(std::initializer_list<std::string> args)
{
    std::vector<std::string> a(args);
    a.insert(a.begin(), "aln");
    std::vector<char *> argv;
    for (std::string &s : a) argv.push_back(&s[0]);
    argv.push_back(nullptr);
    return main_aln((int)a.size(), argv.data());
}

std::string slurp(const std::string &path)
{
    std::ifstream in(path);
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

int n_records(const std::string &sam)
{
    std::istringstream in(sam);
    int n = 0;
    for (std::string l; std::getline(in, l);) n += !l.empty() && l[0] != '@';
    return n;
}

} // namespace

TEST_CASE("aln_cli: repeated parsing resets state when option storage is reused"
          * doctest::test_suite("unit/aln_cli")) {
    char name[] = "aln", prefix[] = "unused-prefix", option[16] = "-s";
    char *first[] = {name, prefix, option, nullptr};
    {
        bwa_tests::StderrCapture capture;
        CHECK(main_aln(3, first) == 1);
        CHECK(capture.text().find("Usage:   bwa-mem3 aln") != std::string::npos);
    }

    // glibc permutes the positional argument past -s and retains a pointer at option + 2.
    // Reuse that storage so an incomplete reset would parse --meth as the short option -m eth.
    std::strcpy(option, "--meth");
    char *second[] = {name, option, nullptr};
    bwa_tests::StderrCapture capture;
    CHECK(main_aln(2, second) == 1);
    const std::string error = capture.text();
    CHECK(error.find("Usage:   bwa-mem3 aln") != std::string::npos);
    CHECK(error.find("-m wants") == std::string::npos);
}

TEST_CASE("aln_cli: short-insert is documented and refused before single-end I/O"
          * doctest::test_suite("unit/aln_cli")) {
    CliIndex &x = fixture();
    for (bool meth : {false, true}) {
        bwa_tests::StderrCapture capture;
        const std::string out = x.dir + "/short-single.sam";
        const int status = meth
            ? run({"--meth", "--short-insert", "-f", out, "/missing-index", "/missing-reads"})
            : run({"--short-insert", "-f", out, "/missing-index", "/missing-reads"});
        const std::string error = capture.text();
        CHECK(status == 1);
        CHECK(error.find("--short-insert requires paired-end input") != std::string::npos);
        CHECK(error.find("index") == std::string::npos);
        CHECK(access(out.c_str(), F_OK) != 0);
    }
    FILE *capture = tmpfile();
    REQUIRE(capture != nullptr);
    fflush(stdout);
    const int saved = dup(fileno(stdout));
    REQUIRE(saved >= 0);
    REQUIRE(dup2(fileno(capture), fileno(stdout)) >= 0);
    const int status = run({"--help"});
    fflush(stdout);
    REQUIRE(dup2(saved, fileno(stdout)) >= 0);
    close(saved);
    rewind(capture);
    std::string help;
    char buf[1024];
    while (fgets(buf, sizeof buf, capture)) help += buf;
    fclose(capture);
    CHECK(status == 0);
    CHECK(help.find("--short-insert") != std::string::npos);
    CHECK(help.find("except with --meth or --short-insert") != std::string::npos);
}

TEST_CASE("aln_cli: meaningless negative values and malformed numbers are refused"
          * doctest::test_suite("unit/aln_cli")) {
    CliIndex &x = fixture();
    const std::string out = x.dir + "/v.sam";
    // a vector of vectors, not nested initializer_lists: an inner list's backing array
    // does not outlive the full-expression that builds the outer one
    const std::vector<std::vector<std::string>> bads = {
             {"-l", "-1"}, {"-l", "-5"}, {"-B", "-1"}, {"-M", "-1"}, {"-O", "-1"},
             {"-E", "-2"}, {"-m", "-1"}, {"-n", "-2"}, {"-n", "-0.5"}, {"-M", "0", "-O", "0", "-E", "0"},
             {"-o", "-5"}, {"-t", "abc"}, {"-q", "20x"}, {"-k", ""}, {"--bam-level", "10"}, {"--max-occ", "1e3"},
             {"-c", "x"}, {"-c", "inf"}, {"-c", "nan"}, {"-n", "0.0"}, {"-n", "0e0"}, {"-n", "1e-60"},
             {"-r", "@RG\\tID:x\\n"}};
    for (const std::vector<std::string> &a : bads) {
        CAPTURE(a[0]);
        CAPTURE(a[1]);
        std::vector<std::string> full(a);
        for (const std::string &s : {std::string("-f"), out, x.prefix, x.fq}) full.push_back(s);
        std::vector<char *> argv;
        std::string name = "aln";
        argv.push_back(&name[0]);
        for (std::string &s : full) argv.push_back(&s[0]);
        argv.push_back(nullptr);
        CHECK(main_aln((int)argv.size() - 1, argv.data()) == 1);
        CHECK(access(out.c_str(), F_OK) != 0); // refused before any output exists
    }
}

TEST_CASE("aln_cli: values stock handles still run"
          * doctest::test_suite("unit/aln_cli")) {
    CliIndex &x = fixture();
    const std::string base = x.dir + "/base.sam";
    REQUIRE(run({"-f", base, x.prefix, x.fq}) == 0);
    CHECK(n_records(slurp(base)) == 4);
    // zero penalties that still leave room, -o -1 (no gap opens), -l 0, -m 0, -e -1, -B 0
    const std::vector<std::vector<std::string>> oks = {
             {"-E", "0"}, {"-M", "0", "-O", "0"}, {"-o", "-1"}, {"-l", "0"},
             {"-m", "0"}, {"-e", "-1"}, {"-B", "0"}, {"-n", "0"}, {"-t", "0"}, {"-i", "-1"}, {"--bam-level", "0"}};
    for (const std::vector<std::string> &a : oks) {
        CAPTURE(a[0]);
        const std::string out = x.dir + "/ok.sam";
        std::vector<std::string> full(a);
        for (const std::string &s : {std::string("-f"), out, x.prefix, x.fq}) full.push_back(s);
        std::vector<char *> argv;
        std::string name = "aln";
        argv.push_back(&name[0]);
        for (std::string &s : full) argv.push_back(&s[0]);
        argv.push_back(nullptr);
        CHECK(main_aln((int)argv.size() - 1, argv.data()) == 0);
        CHECK(n_records(slurp(out)) == 4);
    }
}

TEST_CASE("aln_cli: a run without -r does not reuse an earlier run's read group"
          * doctest::test_suite("unit/aln_cli")) {
    CliIndex &x = fixture();
    const std::string with = x.dir + "/rg.with.sam", without = x.dir + "/rg.without.sam";
    REQUIRE(run({"-r", "@RG\\tID:grp1\\tSM:s", "-f", with, x.prefix, x.fq}) == 0);
    CHECK(slurp(with).find("\tRG:Z:grp1") != std::string::npos);
    REQUIRE(run({"-f", without, x.prefix, x.fq}) == 0);
    CHECK(n_records(slurp(without)) == 4);
    CHECK(slurp(without).find("RG:Z:") == std::string::npos);
    CHECK(slurp(without).find("@RG") == std::string::npos);
}

TEST_CASE("aln_cli: standard input for both read files is refused"
          * doctest::test_suite("unit/aln_cli")) {
    CliIndex &x = fixture();
    const std::string out = x.dir + "/stdin2.sam";
    CHECK(run({"-f", out, x.prefix, "-", "-"}) == 1);
    CHECK(access(out.c_str(), F_OK) != 0);
}

TEST_CASE("aln_cli: -f naming an input is refused and the input is left intact"
          * doctest::test_suite("unit/aln_cli")) {
    CliIndex &x = fixture();
    const std::string copy = x.dir + "/copy.fq";
    {
        std::ofstream o(copy);
        o << slurp(x.fq);
    }
    const std::string before = slurp(copy);
    CHECK(run({"-f", copy, x.prefix, copy}) == 1);
    CHECK(slurp(copy) == before);
    CHECK(run({"-f", copy, x.prefix, x.fq, copy}) == 1);
    CHECK(slurp(copy) == before);
}

TEST_CASE("aln_cli: a malformed record exits 1 after writing what stock writes"
          * doctest::test_suite("unit/aln_cli")) {
    CliIndex &x = fixture();
    const std::string bad = x.dir + "/bad.fq";
    {
        std::ofstream o(bad);
        o << slurp(x.fq) << "@trunc\nACGTACGTAC\n+\nIII\n";
    }
    const std::string good = x.dir + "/good.sam", out = x.dir + "/bad.sam";
    REQUIRE(run({"-f", good, x.prefix, x.fq}) == 0);
    CHECK(run({"-f", out, x.prefix, bad}) == 1);
    CHECK(slurp(out) == slurp(good)); // the four good reads, as stock writes them
}

TEST_CASE("aln_cli: --bam with an -r @SQ that hides the index's contigs is an error"
          * doctest::test_suite("unit/aln_cli")) {
    CliIndex &x = fixture();
    const std::string out = x.dir + "/sq.bam";
    CHECK(run({"--bam", "-f", out, "-r", "@RG\\tID:x\\n@SQ\\tSN:foo\\tLN:10", x.prefix, x.fq}) == 1);
    // SAM keeps stock's output (the names stay text)
    CHECK(run({"-f", x.dir + "/sq.sam", "-r", "@RG\\tID:x\\n@SQ\\tSN:foo\\tLN:10", x.prefix, x.fq}) == 0);
    std::memset(bwa_rg_id, 0, sizeof bwa_rg_id); // -r sets this global; later tests expect no RG tag
}

TEST_CASE("aln_cli: -n with an exponent is a probability, as with a decimal point"
          * doctest::test_suite("unit/aln_cli")) {
    CliIndex &x = fixture();
    const std::string dec = x.dir + "/n_dec.sam", exp = x.dir + "/n_exp.sam", one = x.dir + "/n_one.sam";
    REQUIRE(run({"-f", dec, "-n", "0.01", x.prefix, x.fq}) == 0);
    REQUIRE(run({"-f", exp, "-n", "1e-2", x.prefix, x.fq}) == 0);
    REQUIRE(run({"-f", one, "-n", "1", x.prefix, x.fq}) == 0);
    CHECK(n_records(slurp(exp)) == 4);
    // the @PG line names the command line, so compare from the first record on
    const std::string d = slurp(dec), e = slurp(exp);
    CHECK(d.substr(d.find("\nr0\t")) == e.substr(e.find("\nr0\t")));
}

TEST_CASE("aln_cli: -f naming an index file is refused and the index is left intact"
          * doctest::test_suite("unit/aln_cli")) {
    CliIndex &x = fixture();
    for (const char *sfx : {".pac", ".ann", ".bwt.2bit.64"}) {
        CAPTURE(sfx);
        const std::string comp = x.prefix + sfx;
        const std::string before = slurp(comp);
        REQUIRE(!before.empty());
        CHECK(run({"-f", comp, x.prefix, x.fq}) == 1);
        CHECK(slurp(comp) == before);
    }
    CHECK(run({"-f", x.dir + "/after.sam", x.prefix, x.fq}) == 0); // the index still works
}

TEST_CASE("aln_cli: meth uses the original header and rejects user SQ records"
          * doctest::test_suite("unit/aln_cli")) {
    bwa_aln::aln_meth_ref_t ref;
    bntseq_t bns = {};
    bntann1_t ann = {};
    char name[] = "chr1";
    ann.name = name; ann.len = 100;
    bns.n_seqs = 1; bns.anns = &ann;
    ref.bns = &bns;
    kstring_t hdr = {0, 0, nullptr};
    REQUIRE(bwa_aln::aln_meth_header(&hdr, &ref, "@RG\tID:x", "@PG\tID:aln", nullptr) == 0);
    const std::string text(hdr.s, hdr.l);
    CHECK(text.find(BWAMEM3_DEFAULT_HD_LINE) == 0);
    CHECK(text.find("@SQ\tSN:chr1\tLN:100\n") != std::string::npos);
    CHECK(text.find("@RG\tID:x") != std::string::npos);
    free(hdr.s);
    hdr = {0, 0, nullptr};
    CHECK(bwa_aln::aln_meth_header(&hdr, &ref, "@RG\tID:x\n@SQ\tSN:chr1\tLN:100", nullptr, nullptr) == -1);
    CHECK(hdr.l == 0);
    free(hdr.s);
    CliIndex &x = fixture();
    const std::string out = x.dir + "/meth-sq.sam";
    CHECK(run({"--meth", "-r", "@RG\\tID:x\\n@SQ\\tSN:chr1\\tLN:100", "-f", out, x.prefix, x.fq}) == 1);
    CHECK(access(out.c_str(), F_OK) != 0);
    std::memset(bwa_rg_id, 0, sizeof bwa_rg_id);
}

TEST_CASE("aln_cli: meth reference checks contig names lengths and offsets"
          * doctest::test_suite("unit/aln_cli")) {
    CliIndex &x = fixture();
    bntseq_t *orig = bns_restore(x.prefix.c_str());
    REQUIRE(orig != nullptr);
    std::vector<bntann1_t> anns(2 * orig->n_seqs);
    std::vector<std::string> names(anns.size());
    bntseq_t conv = {};
    conv.n_seqs = (int)anns.size(); conv.l_pac = 2 * orig->l_pac; conv.anns = anns.data();
    for (int i = 0; i < orig->n_seqs; ++i) {
        for (int h = 0; h < 2; ++h) {
            const int j = 2 * i + h;
            names[j] = std::string(h ? "f" : "r") + orig->anns[i].name;
            anns[j].name = &names[j][0];
            anns[j].len = orig->anns[i].len;
            anns[j].offset = 2 * orig->anns[i].offset + h * (int64_t)orig->anns[i].len;
        }
    }
    bwa_aln::aln_meth_ref_t ref;
    REQUIRE(bwa_aln::aln_meth_ref_load(&ref, x.prefix.c_str(), &conv, 1) == 0);
    bwa_aln::aln_meth_ref_free(&ref);
    for (int mutation = 0; mutation < 5; ++mutation) {
        const bntann1_t saved = anns[0];
        if (mutation == 0) anns[0].name = const_cast<char *>("wrong");
        if (mutation == 1) ++anns[0].len;
        if (mutation == 2) ++anns[0].offset;
        if (mutation == 3) --conv.n_seqs;
        if (mutation == 4) --conv.l_pac;
        CHECK(bwa_aln::aln_meth_ref_load(&ref, x.prefix.c_str(), &conv, 1) == -1);
        CHECK(ref.bns == nullptr);
        CHECK(ref.pac == nullptr);
        anns[0] = saved; conv.n_seqs = (int)anns.size(); conv.l_pac = 2 * orig->l_pac;
    }
    bns_destroy(orig);
}

TEST_CASE("aln_cli: meth record writer shares sequence helpers and QC-fail options"
          * doctest::test_suite("unit/aln_cli")) {
    bntann1_t ann = {};
    char contig[] = "chr1", name[] = "read";
    ann.name = contig; ann.len = 100;
    bntseq_t bns = {};
    bns.n_seqs = 1; bns.l_pac = 100; bns.anns = &ann;
    bwa_aln::aln_meth_ref_t ref;
    ref.bns = &bns;
    bwa_aln::aln_meth_opt_t mo;
    mo.m.meth_tags = MEM_METH_TAG_XR | MEM_METH_TAG_XG; // XM needs sequence context, tested separately
    bwa_aln::aln_meth_read_t result;
    result.hyp = 1;
    unsigned char seq[] = {0, 1, 2, 3}, qual[] = "IIII";
    for (int strand = 0; strand < 2; ++strand) {
        for (char failed : {char(0), 'f', 'r'}) {
            bwa_aln::bwa_seq_t read = {};
            read.name = name; read.seq = seq; read.qual = qual;
            read.len = read.full_len = read.clip_len = 4;
            read.type = BWA_TYPE_UNIQUE; read.strand = strand;
            read.pos = 10; read.c1 = 1; read.mapQ = 37;
            mo.m.meth_set_as_failed = failed;
            kstring_t out = {0, 0, nullptr};
            bwa_aln::aln_meth_print_sam1(&out, &ref, &mo, &read, &result, 1, nullptr, BWA_MODE_COMPREAD, 30);
            const std::string text(out.s, out.l);
            const int flag = (strand ? 16 : 0) | (failed == 'f' ? 512 : 0);
            CHECK(text.find("read\t" + std::to_string(flag) + "\tchr1\t11\t37\t4M") == 0);
            CHECK(text.find("\tACGT\tIIII") != std::string::npos);
            CHECK(text.find("\tXR:Z:CT\tXG:Z:CT\n") != std::string::npos);
            CHECK(text.find("XM:i:") == std::string::npos);
            CHECK(text.find("XO:i:") == std::string::npos);
            CHECK(text.find("XG:i:") == std::string::npos);
            free(out.s);
        }
    }
}
