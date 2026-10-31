// Unit tests for the single-end half of `bwa-mem3 aln` (src/aln/bwase.cpp, aln_out.cpp,
// aln_main.cpp). Byte-for-byte agreement with stock bwa aln + samse is
// test/regression/aln_parity.sh; these cases pin the header layout and the once-per-run
// drand48 seeding on the shared synthetic reference (the BAM sink is test_aln_se_bam.cpp).

#include "doctest/doctest.h"

#include "aln_test_ref.h"
#include "../../src/FMI_search.h"
#include "../../src/aln/aln_main.h"
#include "../../src/aln/bwase.h"
#include "../../src/bwa.h"
#include "../../src/read_index_ele.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include <unistd.h>

namespace {

namespace ba = bwa_aln;

struct SeIndex {
    std::string dir, prefix;
    std::unique_ptr<FMI_search> fmi;
    ba::bwt_t *bwt = nullptr;

    SeIndex()
    {
        char tmpl[] = "/tmp/aln_se_unit_XXXXXX";
        REQUIRE(mkdtemp(tmpl) != nullptr);
        dir = tmpl;
        prefix = dir + "/ref.fa";
        REQUIRE(aln_test::write_test_fasta(prefix));
        REQUIRE(bwa_idx_build(prefix.c_str(), prefix.c_str(), 0, 3) == 0);
        fmi.reset(new FMI_search(prefix.c_str()));
        fmi->load_index(true, 1);
        bwt = ba::bwt_from_fmi(fmi.get());
    }
    ~SeIndex()
    {
        ba::bwt_destroy(bwt);
        fmi.reset();
        std::string cmd = "rm -rf '" + dir + "'";
        if (system(cmd.c_str()) != 0) std::fprintf(stderr, "warning: could not remove %s\n", dir.c_str());
    }
    const bntseq_t *bns() const { return fmi->idx->bns; }
};

SeIndex &index_fixture()
{
    static SeIndex idx;
    return idx;
}

std::string header(const char *hdr_line, const char *pg_line)
{
    kstring_t ks = {0, 0, 0};
    ba::aln_print_sam_hdr(&ks, index_fixture().bns(), hdr_line, pg_line);
    std::string s(ks.s ? ks.s : "", ks.l);
    free(ks.s);
    return s;
}

std::vector<std::string> lines(const std::string &s)
{
    std::vector<std::string> v;
    std::istringstream in(s);
    for (std::string l; std::getline(in, l);) v.push_back(l);
    return v;
}

std::string write_fastq(const std::string &dir, const std::string &name, const std::vector<std::string> &reads)
{
    std::string path = dir + "/" + name;
    FILE *fp = fopen(path.c_str(), "w");
    REQUIRE(fp != nullptr);
    for (size_t i = 0; i < reads.size(); ++i)
        fprintf(fp, "@r%zu\n%s\n+\n%s\n", i, reads[i].c_str(), std::string(reads[i].size(), 'I').c_str());
    fclose(fp);
    return path;
}

// Runs the SE pipeline over `fq` and returns the SAM text it wrote.
std::string run_se(int batch_size, const std::string &fq, const std::string &out_path)
{
    SeIndex &ix = index_fixture();
    ba::aln_opt_t o;
    ba::aln_opt_init(&o);
    o.batch_size = batch_size;
    o.gap.max_diff = 2;
    o.gap.fnr = -1.0;
    ba::aln_out_t *out = ba::aln_out_open(out_path.c_str(), 0, 0, ix.bns(), "");
    REQUIRE(out != nullptr);
    REQUIRE(ba::aln_run_se(&o, ix.bwt, ix.bns(), ix.fmi->idx->pac, fq.c_str(), out) == 0);
    REQUIRE(ba::aln_out_close(out) == 0);
    std::ifstream in(out_path);
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

} // namespace

TEST_CASE("aln_se: header is @HD, @SQ lines, then RG, then PG"
          * doctest::test_suite("unit/aln_se")) {
    std::vector<std::string> l = lines(header("@RG\tID:x\tSM:y", "@PG\tID:bwa-mem3\n"));
    REQUIRE(l.size() == 6);
    CHECK(l[0] == "@HD\tVN:1.5\tSO:unsorted\tGO:query");
    CHECK(l[1] == "@SQ\tSN:c1\tLN:600");
    CHECK(l[2] == "@SQ\tSN:c2\tLN:550");
    CHECK(l[3] == "@SQ\tSN:c3\tLN:128");
    CHECK(l[4] == "@RG\tID:x\tSM:y");
    CHECK(l[5] == "@PG\tID:bwa-mem3");
    // no -r and a @PG line without its newline
    l = lines(header(nullptr, "@PG\tID:p"));
    REQUIRE(l.size() == 5);
    CHECK(l[4] == "@PG\tID:p");
}

TEST_CASE("aln_se: @HD in -r text suppresses the default @HD"
          * doctest::test_suite("unit/aln_se")) {
    std::vector<std::string> l = lines(header("@HD\tVN:1.6\tSO:coordinate\n@RG\tID:x", nullptr));
    REQUIRE(l.size() == 5);
    CHECK(l[0] == "@SQ\tSN:c1\tLN:600");
    CHECK(l[3] == "@HD\tVN:1.6\tSO:coordinate");
    CHECK(l[4] == "@RG\tID:x");
    int n_hd = 0;
    for (const std::string &s : l) n_hd += s.compare(0, 3, "@HD") == 0;
    CHECK(n_hd == 1);
}

TEST_CASE("aln_se: drand48 is seeded once per run, not per batch"
          * doctest::test_suite("unit/aln_se")) {
    SeIndex &ix = index_fixture();
    // c1[150:186] lies in the c1/c2 reverse-complement copy: two equally good hits, so every
    // read draws from drand48. Four copies of it, as one batch of 4 and as two batches of 2.
    std::string rep = aln_test::test_contigs()[0].second.substr(150, 36);
    std::string fq = write_fastq(ix.dir, "rep.fq", {rep, rep, rep, rep});
    std::string one = run_se(4, fq, ix.dir + "/one.sam");
    std::string two = run_se(2, fq, ix.dir + "/two.sam");
    CHECK(one == two);
    CHECK(one.find("X0:i:2") != std::string::npos);
    // A re-run in the same process reseeds, so it reproduces the first run.
    CHECK(run_se(4, fq, ix.dir + "/again.sam") == one);
}
