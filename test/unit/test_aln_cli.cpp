// Command-line checks of `bwa-mem3 aln` (main_aln in src/aln/aln_main.cpp): values that are
// meaningless and on which stock crashes, aborts or silently drops every read are refused with an
// error and exit 1, as are malformed numbers; inputs that would lose data (standard input for both
// files, -f naming an input) are refused before any output exists. None of this changes the
// output of a run stock completes; aln_parity.sh holds that output to stock's.

#include "doctest/doctest.h"

#include "aln_test_ref.h"
#include "aln_bam_decode.h"
#include "scoped_env.h"
#include "stderr_capture.h"
#include "../../src/FMI_search.h"
#include "../framework/stderr_capture.h"
#include "../../src/aln/aln_main.h"
#include "../../src/aln/aln_meth.h"
#include "../../src/compat_target.h"
#include "../../src/bwa.h"

#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <initializer_list>
#include <sstream>
#include <string>
#include <thread>
#include <sys/stat.h>
#include <vector>
#include <unistd.h>
#include <utime.h>

namespace {

struct CliIndex {
    bwa_tests::TempDir tmp_dir{"aln_cli_unit"}; // removed only by the process that made it
    std::string dir = tmp_dir.path(), prefix, fq;

    CliIndex()
    {
        prefix = dir + "/" + dir.substr(dir.rfind('/') + 1) + ".fa";
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
};

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

std::string expected_cli_records(int mapq, const std::string &read_group = "")
{
    const std::string c1 = aln_test::test_contigs()[0].second;
    std::string records;
    for (int i = 0; i < 4; ++i)
        records += "r" + std::to_string(i) + "\t0\tc1\t" + std::to_string(421 + 20 * i)
                + "\t" + std::to_string(mapq) + "\t36M\t*\t0\t0\t"
                + c1.substr(420 + 20 * i, 36) + "\t" + std::string(36, 'I')
                + (read_group.empty() ? "" : "\tRG:Z:" + read_group)
                + "\tXT:A:U\tNM:i:0\tX0:i:1\tX1:i:0\tXM:i:0\tXO:i:0\tXG:i:0\tMD:Z:36\n";
    return records;
}

int n_records(const std::string &sam)
{
    std::istringstream in(sam);
    int n = 0;
    for (std::string l; std::getline(in, l);) n += !l.empty() && l[0] != '@';
    return n;
}

} // namespace

TEST_CASE("aln_cli: synthetic index basename is unique to its temporary directory"
          * doctest::test_suite("unit/aln_cli")) {
    CliIndex x;
    const std::string leaf = x.dir.substr(x.dir.rfind('/') + 1);
    CHECK(x.prefix == x.dir + "/" + leaf + ".fa");
}

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
    CliIndex x;
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
    CHECK(help.find("unless --meth, --short-insert, --ref-span-insert or --exact-rescue-window") != std::string::npos);
}

TEST_CASE("aln_cli: exact-rescue-window is paired only, its --no- form is accepted anywhere, and both are documented"
          * doctest::test_suite("unit/aln_cli")) {
    CliIndex x;
    for (bool meth : {false, true}) {
        bwa_tests::StderrCapture capture;
        const std::string out = x.dir + "/exact-single.sam";
        const int status = meth
            ? run({"--meth", "--exact-rescue-window", "-f", out, "/missing-index", "/missing-reads"})
            : run({"--exact-rescue-window", "-f", out, "/missing-index", "/missing-reads"});
        const std::string error = capture.text();
        CHECK(status == 1);
        CHECK(error.find("--exact-rescue-window requires paired-end input") != std::string::npos);
        CHECK(error.find("index") == std::string::npos);
        CHECK(access(out.c_str(), F_OK) != 0);
    }
    { // the --no- form only restores sampe's window, so one FASTQ is fine; flags take no argument
        bwa_tests::StderrCapture capture;
        const std::string se = x.dir + "/exact-off-se.sam", pe = x.dir + "/exact-pe.sam";
        CHECK(run({"--no-exact-rescue-window", "-f", se, x.prefix, x.fq}) == 0);
        CHECK(n_records(slurp(se)) == 4);
        CHECK(run({"--exact-rescue-window", "--no-exact-rescue-window", "--exact-rescue-window", "-f", pe,
                   x.prefix, x.fq, x.fq}) == 0);
        CHECK(n_records(slurp(pe)) == 8);
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
    CHECK(help.find("  --exact-rescue-window") != std::string::npos);
    CHECK(help.find("[off; on under --meth]") != std::string::npos);
    CHECK(help.find("--no-exact-rescue-window keeps sampe's window under --meth") != std::string::npos);
}

TEST_CASE("aln_cli: ref-span-insert parses, is documented and is refused before single-end I/O"
          * doctest::test_suite("unit/aln_cli")) {
    CliIndex x;
    for (bool meth : {false, true}) {
        bwa_tests::StderrCapture capture;
        const std::string out = x.dir + "/span-single.sam";
        const int status = meth
            ? run({"--meth", "--ref-span-insert", "-f", out, "/missing-index", "/missing-reads"})
            : run({"--ref-span-insert", "-f", out, "/missing-index", "/missing-reads"});
        const std::string error = capture.text();
        CHECK(status == 1);
        CHECK(error.find("--ref-span-insert requires paired-end input") != std::string::npos);
        CHECK(error.find("index") == std::string::npos);
        CHECK(access(out.c_str(), F_OK) != 0);
    }
    { // parsed as a flag, not as taking the next argument: paired input runs to completion
        const std::string out = x.dir + "/span-pe.sam";
        bwa_tests::StderrCapture capture;
        CHECK(run({"--ref-span-insert", "-f", out, x.prefix, x.fq, x.fq}) == 0);
        CHECK(n_records(slurp(out)) == 8);
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
    CHECK(help.find("  --ref-span-insert") != std::string::npos);
}

TEST_CASE("aln_cli: meaningless negative values and malformed numbers are refused"
          * doctest::test_suite("unit/aln_cli")) {
    CliIndex x;
    bwa_tests::TempDir scratch{"aln_cli_case"}; // this case's files, removed when it ends
    const std::string out = scratch.path() + "/v.sam";
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

TEST_CASE("aln_cli: thread counts above the supported maximum are refused before work"
          * doctest::test_suite("unit/aln_cli")) {
    bwa_tests::TempDir scratch{"aln_cli_threads"};
    const std::string out = scratch.file("out.sam");
    for (const std::string &count : {std::string("257"), std::string("16000"), std::string("2147483647")}) {
        CAPTURE(count);
        bwa_tests::StderrCapture err;
        CHECK(run({"-t", count, "-f", out, scratch.file("missing.fa"), scratch.file("missing.fq")}) == 1);
        CHECK(err.text().find("-t wants an integer from") != std::string::npos);
        CHECK(err.text().find("Index file found") == std::string::npos);
        CHECK(access(out.c_str(), F_OK) != 0);
        CHECK(run({"-t", count, "--help"}) == 1);
    }
    CHECK(run({"-t", "256", "--help"}) == 0); // parse the boundary without starting 256 workers
}

TEST_CASE("aln_cli: small thread counts and the lower clamp still align reads"
          * doctest::test_suite("unit/aln_cli")) {
    CliIndex x;
    bwa_tests::TempDir scratch{"aln_cli_threads"};
    const std::string want = expected_cli_records(37);
    for (const std::string &count : {std::string("-1"), std::string("0"), std::string("1"), std::string("2")}) {
        CAPTURE(count);
        const std::string out = scratch.file("threads.sam");
        REQUIRE(run({"-t", count, "-f", out, x.prefix, x.fq}) == 0);
        const std::string sam = slurp(out);
        CHECK(n_records(sam) == 4);
        std::istringstream in(sam);
        std::string got;
        for (std::string line; std::getline(in, line);)
            if (!line.empty() && line[0] != '@') got += line + "\n";
        CHECK(got == want);
    }
}

TEST_CASE("aln_cli: values stock handles still run"
          * doctest::test_suite("unit/aln_cli")) {
    CliIndex x;
    bwa_tests::TempDir scratch{"aln_cli_case"}; // this case's files, removed when it ends
    const std::string base = scratch.path() + "/base.sam";
    REQUIRE(run({"-f", base, x.prefix, x.fq}) == 0);
    CHECK(n_records(slurp(base)) == 4);
    // zero penalties that still leave room, -o -1 (no gap opens), -l 0, -m 0, -e -1, -B 0
    const std::vector<std::vector<std::string>> oks = {
             {"-E", "0"}, {"-M", "0", "-O", "0"}, {"-o", "-1"}, {"-l", "0"},
             {"-m", "0"}, {"-e", "-1"}, {"-B", "0"}, {"-n", "0"}, {"-t", "0"}, {"-i", "-1"}, {"--bam-level", "0"}};
    for (const std::vector<std::string> &a : oks) {
        CAPTURE(a[0]);
        const std::string out = scratch.path() + "/ok.sam";
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
    // -n 0 allows no differences, so the four exact reads still map, uniquely, at their planted
    // c1 positions; a hit that uses the whole difference budget (here 0) gets samse's MAPQ 25,
    // not the 37 of a unique hit with differences to spare. The records are built from c1.
    const std::string n0 = scratch.path() + "/n0.sam";
    REQUIRE(run({"-n", "0", "-f", n0, x.prefix, x.fq}) == 0);
    const std::string want = expected_cli_records(25);
    std::string got;
    std::istringstream in(slurp(n0));
    for (std::string l; std::getline(in, l);)
        if (!l.empty() && l[0] != '@') got += l + "\n";
    CHECK(got == want);
}

TEST_CASE("aln_cli: read-group blank lines are refused before output"
          * doctest::test_suite("unit/aln_cli")) {
    CliIndex x;
    bwa_tests::TempDir scratch{"aln_cli_blank_rg"};
    const std::vector<std::string> invalid = {
        "@RG\\tID:x\\n", "@RG\\tID:x\\n\\nSM:y", "@RG\\tID:x\n\nSM:y"};
    for (size_t i = 0; i < invalid.size(); ++i) {
        CAPTURE(i);
        const std::string out = scratch.file("bad" + std::to_string(i) + ".sam");
        bwa_tests::StderrCapture err;
        CHECK(run({"-r", invalid[i], "-f", out, x.prefix, x.fq}) == 1);
        CHECK(err.text().find("empty SAM header line") != std::string::npos);
        CHECK(access(out.c_str(), F_OK) != 0);
    }
    const std::string out = scratch.file("valid.sam");
    REQUIRE(run({"-r", "@RG\\tID:x\\n@CO\\tcomment", "-f", out, x.prefix, x.fq}) == 0);
    CHECK(slurp(out).find("@RG\tID:x\n@CO\tcomment\n") != std::string::npos);
    CHECK(n_records(slurp(out)) == 4);
}

TEST_CASE("aln_cli: a run without -r does not reuse an earlier run's read group"
          * doctest::test_suite("unit/aln_cli")) {
    CliIndex x;
    bwa_tests::TempDir scratch{"aln_cli_case"}; // this case's files, removed when it ends
    const std::string with = scratch.path() + "/rg.with.sam", without = scratch.path() + "/rg.without.sam";
    REQUIRE(run({"-r", "@RG\\tID:grp1\\tSM:s", "-f", with, x.prefix, x.fq}) == 0);
    CHECK(slurp(with).find("\tRG:Z:grp1") != std::string::npos);
    REQUIRE(run({"-f", without, x.prefix, x.fq}) == 0);
    CHECK(n_records(slurp(without)) == 4);
    CHECK(slurp(without).find("RG:Z:") == std::string::npos);
    CHECK(slurp(without).find("@RG") == std::string::npos);
}

TEST_CASE("aln_cli: standard input for both read files is refused"
          * doctest::test_suite("unit/aln_cli")) {
    CliIndex x;
    bwa_tests::TempDir scratch{"aln_cli_case"}; // this case's files, removed when it ends
    const std::string out = scratch.path() + "/stdin2.sam";
    CHECK(run({"-f", out, x.prefix, "-", "-"}) == 1);
    CHECK(access(out.c_str(), F_OK) != 0);
}

TEST_CASE("aln_cli: -f naming an input is refused and the input is left intact"
          * doctest::test_suite("unit/aln_cli")) {
    CliIndex x;
    bwa_tests::TempDir scratch{"aln_cli_case"}; // this case's files, removed when it ends
    const std::string copy = scratch.path() + "/copy.fq";
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
    CliIndex x;
    bwa_tests::TempDir scratch{"aln_cli_case"}; // this case's files, removed when it ends
    const std::string bad = scratch.path() + "/bad.fq";
    {
        std::ofstream o(bad);
        o << slurp(x.fq) << "@trunc\nACGTACGTAC\n+\nIII\n";
    }
    const std::string good = scratch.path() + "/good.sam", out = scratch.path() + "/bad.sam";
    REQUIRE(run({"-f", good, x.prefix, x.fq}) == 0);
    CHECK(run({"-f", out, x.prefix, bad}) == 1);
    CHECK(slurp(out) == slurp(good)); // the four good reads, as stock writes them
}

TEST_CASE("aln_cli: --bam with an -r @SQ that hides the index's contigs is an error"
          * doctest::test_suite("unit/aln_cli")) {
    CliIndex x;
    bwa_tests::TempDir scratch{"aln_cli_case"}; // this case's files, removed when it ends
    const std::string out = scratch.path() + "/sq.bam";
    CHECK(run({"--bam", "-f", out, "-r", "@RG\\tID:x\\n@SQ\\tSN:foo\\tLN:10", x.prefix, x.fq}) == 1);
    // SAM keeps stock's output (the names stay text)
    CHECK(run({"-f", scratch.path() + "/sq.sam", "-r", "@RG\\tID:x\\n@SQ\\tSN:foo\\tLN:10", x.prefix, x.fq}) == 0);
    std::memset(bwa_rg_id, 0, sizeof bwa_rg_id); // -r sets this global; later tests expect no RG tag
}

TEST_CASE("aln_cli: output preflight does not create and delete absent paths"
          * doctest::test_suite("unit/aln_cli")) {
    CliIndex x;
    bwa_tests::TempDir scratch{"aln_cli_case"};
    const std::string target = scratch.path() + "/absent.sam";
    const std::string link = scratch.path() + "/out.sam";
    REQUIRE(symlink(target.c_str(), link.c_str()) == 0);
    const struct utimbuf old_times = {1000000000, 1000000000};
    for (const std::string &out : {target, link}) {
        REQUIRE(utime(scratch.path().c_str(), &old_times) == 0);
        // Invalid BAM reference lengths stop after preflight but before the real output open.
        CHECK(run({"--bam", "-f", out, "-r", "@RG\\tID:x\\n@SQ\\tSN:c1\\tLN:10", x.prefix, x.fq}) == 1);
        struct stat st;
        REQUIRE(stat(scratch.path().c_str(), &st) == 0);
        CHECK(st.st_mtime == old_times.modtime);
        CHECK(access(target.c_str(), F_OK) != 0);
        REQUIRE(lstat(link.c_str(), &st) == 0);
        CHECK(S_ISLNK(st.st_mode));
    }
    std::memset(bwa_rg_id, 0, sizeof bwa_rg_id);
}

TEST_CASE("aln_cli: --bam with an -r @SQ whose length differs from the index's is an error"
          * doctest::test_suite("unit/aln_cli")) {
    CliIndex x;
    bwa_tests::TempDir scratch{"aln_cli_case"}; // this case's files, removed when it ends
    // c1 is 600 bp in the index and the reads map at 421-497: an LN:10 header would put them
    // past the end of the declared contig
    const std::string sq = "@RG\\tID:x\\n@SQ\\tSN:c1\\tLN:10\\n@SQ\\tSN:c2\\tLN:550\\n@SQ\\tSN:c3\\tLN:128";
    const std::string out = scratch.path() + "/sq_len.bam";
    CHECK(run({"--bam", "-f", out, "-r", sq, x.prefix, x.fq}) == 1);
    // the same lengths as the index are accepted
    const std::string ok = "@RG\\tID:x\\n@SQ\\tSN:c1\\tLN:600\\n@SQ\\tSN:c2\\tLN:550\\n@SQ\\tSN:c3\\tLN:128";
    const std::string accepted = scratch.path() + "/sq_ok.bam";
    REQUIRE(run({"--bam", "-f", accepted, "-r", ok, x.prefix, x.fq}) == 0);
    CHECK(bwa_tests::decoded_aln_records(accepted) == expected_cli_records(37, "x"));
    // SAM keeps stock's output, which writes the header as given
    CHECK(run({"-f", scratch.path() + "/sq_len.sam", "-r", sq, x.prefix, x.fq}) == 0);
    std::memset(bwa_rg_id, 0, sizeof bwa_rg_id); // -r sets this global; later tests expect no RG tag
}

TEST_CASE("aln_cli: --bam rejects an empty final duplicate LN before opening output"
          * doctest::test_suite("unit/aln_cli")) {
    CliIndex x;
    bwa_tests::TempDir scratch{"aln_cli_case"};
    const std::string out = scratch.path() + "/empty_ln.bam";
    const std::string sq = "@RG\\tID:x\\n@SQ\\tSN:c1\\tLN:600\\tLN:";
    CHECK(run({"--bam", "-f", out, "-r", sq, x.prefix, x.fq}) == 1);
    CHECK(access(out.c_str(), F_OK) != 0);
    std::memset(bwa_rg_id, 0, sizeof bwa_rg_id);
}

TEST_CASE("aln_cli: alternative-hit limits reject INT_MAX before opening output"
          * doctest::test_suite("unit/aln_cli")) {
    CliIndex x;
    bwa_tests::TempDir scratch{"aln_cli_case"};
    for (const char *flag : {"--max-xa-hits", "--max-xa-hits-discordant"}) {
        CAPTURE(flag);
        const std::string out = scratch.path() + "/xa.sam";
        CHECK(run({flag, "2147483647", "-f", out, x.prefix, x.fq, x.fq}) == 1);
        CHECK(access(out.c_str(), F_OK) != 0);
        std::remove(out.c_str());
        CHECK(run({flag, "2147483646", "-f", out, x.prefix, x.fq, x.fq}) == 0);
        CHECK(n_records(slurp(out)) == 8);
        std::remove(out.c_str());
    }
}

TEST_CASE("aln_cli: -n with an exponent is a probability, as with a decimal point"
          * doctest::test_suite("unit/aln_cli")) {
    CliIndex x;
    bwa_tests::TempDir scratch{"aln_cli_case"}; // this case's files, removed when it ends
    const std::string dec = scratch.path() + "/n_dec.sam", exp = scratch.path() + "/n_exp.sam";
    REQUIRE(run({"-f", dec, "-n", "0.01", x.prefix, x.fq}) == 0);
    // the integer parse would refuse "1e-2", so this exit status is what shows it was read as a
    // probability (on these exact reads -n 1 and -n 0.01 give the same records)
    REQUIRE(run({"-f", exp, "-n", "1e-2", x.prefix, x.fq}) == 0);
    CHECK(n_records(slurp(exp)) == 4);
    // the @PG line names the command line, so compare from the first record on
    const std::string d = slurp(dec), e = slurp(exp);
    CHECK(d.substr(d.find("\nr0\t")) == e.substr(e.find("\nr0\t")));

    // A read with two mismatches and no closer hit tells the two budgets apart: -n 0.01 allows
    // three differences on 36 bp, while an integer-prefix parse of "1e-2" (stock's atoi) gives -n 1
    // and leaves it unmapped. The expected records are stock bwa 0.7.19's aln + samse output.
    std::string r = aln_test::test_contigs()[0].second.substr(450, 36);
    for (int i : {10, 25}) r[i] = r[i] == 'A' ? 'C' : r[i] == 'C' ? 'G' : r[i] == 'G' ? 'T' : 'A';
    const std::string fq2 = scratch.path() + "/two_mm.fq", exp2 = scratch.path() + "/n_exp2.sam", int1 = scratch.path() + "/n_int1.sam";
    {
        std::ofstream o(fq2);
        o << "@d2\n" << r << "\n+\n" << std::string(36, 'I') << "\n";
    }
    const std::string q(36, 'I');
    REQUIRE(run({"-f", exp2, "-n", "1e-2", x.prefix, fq2}) == 0);
    const std::string e2 = slurp(exp2);
    CHECK(e2.substr(e2.find("d2\t")) == "d2\t0\tc1\t451\t37\t36M\t*\t0\t0\t" + r + "\t" + q
                                           + "\tXT:A:U\tNM:i:2\tX0:i:1\tX1:i:0\tXM:i:2\tXO:i:0\tXG:i:0\tMD:Z:10C14A10\n");
    REQUIRE(run({"-f", int1, "-n", "1", x.prefix, fq2}) == 0);
    const std::string i1 = slurp(int1);
    CHECK(i1.substr(i1.find("d2\t")) == "d2\t4\t*\t0\t0\t*\t*\t0\t0\t" + r + "\t" + q + "\n");
}

TEST_CASE("aln_cli: -f naming an index file is refused and the index is left intact"
          * doctest::test_suite("unit/aln_cli")) {
    CliIndex x;
    bwa_tests::TempDir scratch{"aln_cli_case"}; // this case's files, removed when it ends
    for (const char *sfx : {".pac", ".ann", ".bwt.2bit.64"}) {
        CAPTURE(sfx);
        const std::string comp = x.prefix + sfx;
        const std::string before = slurp(comp);
        REQUIRE(!before.empty());
        CHECK(run({"-f", comp, x.prefix, x.fq}) == 1);
        CHECK(slurp(comp) == before);
    }
    CHECK(run({"-f", scratch.path() + "/after.sam", x.prefix, x.fq}) == 0); // the index still works
}

TEST_CASE("aln_cli: -f naming a file of an index --meth seed index is refused and left intact"
          * doctest::test_suite("unit/aln_cli")) {
    CliIndex x;
    // `index --meth` writes <prefix>.meth.{pac,ann,amb,bwt.2bit.64} beside the original index;
    // aln reads only the original, but overwriting the seed index would break `mem --meth`.
    for (const char *sfx : {".meth.pac", ".meth.ann", ".meth.amb", ".meth.bwt.2bit.64"}) {
        CAPTURE(sfx);
        const std::string comp = x.prefix + sfx;
        const std::string before = std::string("seed index stand-in ") + sfx + "\n";
        {
            std::ofstream o(comp);
            o << before;
        }
        CHECK(run({"-f", comp, x.prefix, x.fq}) == 1);
        CHECK(slurp(comp) == before);
        std::remove(comp.c_str());
    }
}

TEST_CASE("aln_cli: -f naming the file on standard input is refused and the input is left intact"
          * doctest::test_suite("unit/aln_cli")) {
    CliIndex x;
    bwa_tests::TempDir scratch{"aln_cli_case"}; // this case's files, removed when it ends
    const std::string copy = scratch.path() + "/stdin_copy.fq";
    {
        std::ofstream o(copy);
        o << slurp(x.fq);
    }
    const std::string before = slurp(copy);
    const int saved = dup(STDIN_FILENO);
    REQUIRE(saved >= 0);
    const int fd = open(copy.c_str(), O_RDONLY);
    REQUIRE(fd >= 0);
    REQUIRE(dup2(fd, STDIN_FILENO) == STDIN_FILENO);
    close(fd);
    const int se = run({"-f", copy, x.prefix, "-"});
    const int pe = run({"-f", copy, x.prefix, x.fq, "-"});
    REQUIRE(dup2(saved, STDIN_FILENO) == STDIN_FILENO);
    close(saved);
    CHECK(se == 1);
    CHECK(pe == 1);
    CHECK(slurp(copy) == before);
}

TEST_CASE("aln_cli: -f naming a dangling symlink writes through the link"
          * doctest::test_suite("unit/aln_cli")) {
    CliIndex x;
    bwa_tests::TempDir scratch{"aln_cli_case"}; // this case's files, removed when it ends
    const std::string target = scratch.path() + "/link_target.sam", link = scratch.path() + "/link.sam";
    REQUIRE(symlink(target.c_str(), link.c_str()) == 0);
    REQUIRE(run({"-f", link, x.prefix, x.fq}) == 0);
    struct stat st;
    REQUIRE(lstat(link.c_str(), &st) == 0);
    CHECK(S_ISLNK(st.st_mode));
    CHECK(n_records(slurp(target)) == 4);
}

TEST_CASE("aln_cli: -f naming a FIFO gives its reader the whole output"
          * doctest::test_suite("unit/aln_cli")) {
    CliIndex x;
    bwa_tests::TempDir scratch{"aln_cli_case"}; // this case's files, removed when it ends
    const std::string fifo = scratch.path() + "/out.fifo";
    REQUIRE(mkfifo(fifo.c_str(), 0600) == 0);
    std::string first;
    std::atomic<bool> reader_done{false};
    std::thread reader([&] {
        auto drain = [&](std::string *dst) {
            const int fd = open(fifo.c_str(), O_RDONLY);
            if (fd < 0) return;
            char buf[4096];
            for (ssize_t n; (n = read(fd, buf, sizeof buf)) > 0;) if (dst) dst->append(buf, (size_t)n);
            close(fd);
        };
        drain(&first);
        // A writer that opened and closed the FIFO early (a probe) ended the read above with no
        // records; open it again so the real write does not block forever, and fail below.
        if (n_records(first) == 0) drain(nullptr);
        reader_done.store(true);
    });
    const int rc = run({"-f", fifo, x.prefix, x.fq});
    // run() may have failed before opening the FIFO, leaving the reader blocked in open(): open and
    // close it as a writer until the reader is done (non-blocking, so it fails rather than waits
    // while the reader is not yet, or no longer, in open()).
    while (!reader_done.load()) {
        const int fd = open(fifo.c_str(), O_WRONLY | O_NONBLOCK);
        if (fd >= 0) close(fd);
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    reader.join();
    CHECK(rc == 0);
    // These exact reads are unique forward matches at their planted c1 positions.
    const std::string c1 = aln_test::test_contigs()[0].second;
    std::string want, got;
    for (int i = 0; i < 4; ++i)
        want += "r" + std::to_string(i) + "\t0\tc1\t" + std::to_string(421 + 20 * i) + "\t37\t36M\t*\t0\t0\t"
                + c1.substr(420 + 20 * i, 36) + "\t" + std::string(36, 'I')
                + "\tXT:A:U\tNM:i:0\tX0:i:1\tX1:i:0\tXM:i:0\tXO:i:0\tXG:i:0\tMD:Z:36\n";
    std::istringstream in(first);
    for (std::string l; std::getline(in, l);)
        if (!l.empty() && l[0] != '@') got += l + "\n";
    CHECK(got == want);
}

TEST_CASE("aln_cli: aln loads the index without building mem's seed jump table"
          * doctest::test_suite("unit/aln_cli")) {
    // The table builder is the only reader of BWA3_BWTSEED_JUMP, and it reports an invalid
    // value on stderr; so an invalid value shows whether a load built the table.
    CliIndex x;
    bwa_tests::TempDir scratch{"aln_cli_case"}; // this case's files, removed when it ends
    bwa_tests::ScopedEnv env("BWA3_BWTSEED_JUMP", "not-a-depth");
    const std::string diag = "BWA3_BWTSEED_JUMP=\"not-a-depth\"";
    SUBCASE("a building load reports the value (positive control)") {
        bwa_tests::StderrCapture cap;
        {
            FMI_search fmi(x.prefix.c_str());
            fmi.load_index(true, 1);
        }
        CHECK(cap.text().find(diag) != std::string::npos);
    }
    SUBCASE("main_aln does not") {
        const std::string sai = scratch.path() + "/r.sai";
        bwa_tests::StderrCapture cap;
        const int rc = run({"-f", sai, x.prefix, x.fq});
        const std::string err = cap.text();
        CHECK(rc == 0);
        CHECK(err.find(diag) == std::string::npos);
    }
}

TEST_CASE("aln_cli: a batch size below 1 is an error in every driver"
          * doctest::test_suite("unit/aln_cli")) {
    CliIndex x;
    bwa_tests::TempDir scratch{"aln_cli_case"}; // this case's files, removed when it ends
    FMI_search fmi(x.prefix.c_str());
    fmi.load_index(true, 1, /*build_seed_jump=*/false);  // as main_aln loads it
    bwa_aln::bwt_t *bwt = bwa_aln::bwt_from_fmi(&fmi);
    bwa_aln::aln_opt_t o;
    bwa_aln::aln_opt_init(&o);
    for (int bs : {0, -1}) {
        CAPTURE(bs);
        o.batch_size = bs;
        const std::string out_path = scratch.path() + "/batch" + std::to_string(-bs) + ".sam";
        bwa_aln::aln_out_t *out = bwa_aln::aln_out_open(out_path.c_str(), 0, 0, "");
        REQUIRE(out != nullptr);
        CHECK(bwa_aln::aln_run_se(&o, bwt, fmi.idx->bns, fmi.idx->pac, x.fq.c_str(), out) == 1);
        CHECK(bwa_aln::aln_run_pe(&o, bwt, fmi.idx->bns, fmi.idx->pac, x.fq.c_str(), x.fq.c_str(), out) == 1);
        CHECK(bwa_aln::aln_run_se_meth(&o, bwt, x.fq.c_str(), out) == 1);
        CHECK(bwa_aln::aln_run_pe_meth(&o, bwt, x.fq.c_str(), x.fq.c_str(), out) == 1);
        CHECK(bwa_aln::aln_out_close(out) == 0);
    }
    bwa_aln::bwt_destroy(bwt);
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
    CliIndex x;
    const std::string out = x.dir + "/meth-sq.sam";
    CHECK(run({"--meth", "-r", "@RG\\tID:x\\n@SQ\\tSN:chr1\\tLN:100", "-f", out, x.prefix, x.fq}) == 1);
    CHECK(access(out.c_str(), F_OK) != 0);
    std::memset(bwa_rg_id, 0, sizeof bwa_rg_id);
}

TEST_CASE("aln_cli: meth reference checks contig names lengths and offsets"
          * doctest::test_suite("unit/aln_cli")) {
    CliIndex x;
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
    // Asymmetric sequence and qualities: the reverse-strand record must reverse-complement the
    // sequence and reverse the qualities.
    unsigned char seq[] = {0, 0, 1, 2}, qual[5];
    for (int strand = 0; strand < 2; ++strand) {
        for (char failed : {char(0), 'f', 'r'}) {
            std::memcpy(qual, "ABCD", sizeof qual); // the writer reverses a reverse-strand quality in place
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
            CHECK(text.find(strand ? "\tCGTT\tDCBA\t" : "\tAACG\tABCD\t") != std::string::npos);
            CHECK(text.find("\tXR:Z:CT\tXG:Z:CT\n") != std::string::npos);
            CHECK(text.find("XM:i:") == std::string::npos);
            CHECK(text.find("XO:i:") == std::string::npos);
            CHECK(text.find("XG:i:") == std::string::npos);
            free(out.s);
        }
    }
}
