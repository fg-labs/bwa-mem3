// Command-line checks of `bwa-mem3 aln` (main_aln in src/aln/aln_main.cpp): values that are
// meaningless and on which stock crashes, aborts or silently drops every read are refused with an
// error and exit 1, as are malformed numbers; inputs that would lose data (standard input for both
// files, -f naming an input) are refused before any output exists. None of this changes the
// output of a run stock completes; aln_parity.sh holds that output to stock's.

#include "doctest/doctest.h"

#include "aln_test_ref.h"
#include "../../src/aln/aln_main.h"
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
        std::string cmd = "rm -rf '" + dir + "'";
        if (system(cmd.c_str()) != 0) std::fprintf(stderr, "warning: could not remove %s\n", dir.c_str());
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
             {"-c", "x"}};
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
