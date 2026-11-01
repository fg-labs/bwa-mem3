// BAM sink of `bwa-mem3 aln` (src/aln/aln_out.cpp): reads the file back with htslib.
// A separate file from test_aln_se.cpp because htslib's kstring_t and bwa-mem3's have
// different C++ names, so this file must not call the functions that take one.

// htslib before any bwa-mem3 header that pulls in kstring.h (they share the KSTRING_H guard).
#include "htslib/sam.h"
#include "htslib/kstring.h"

#include "doctest/doctest.h"

#include "aln_test_ref.h"
#include "../../src/FMI_search.h"
#include "../../src/aln/aln_main.h"
#include "../../src/bwa.h"
#include "../../src/read_index_ele.h"

#include <cstdio>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>
#include <unistd.h>

namespace {

namespace ba = bwa_aln;

struct BamIndex {
    std::string dir, prefix;
    std::unique_ptr<FMI_search> fmi;
    ba::bwt_t *bwt = nullptr;

    BamIndex()
    {
        std::string tmpl_s = aln_test::tmp_template("aln_se_bam");
        char *tmpl = &tmpl_s[0];
        REQUIRE(mkdtemp(tmpl) != nullptr);
        dir = tmpl;
        prefix = dir + "/ref.fa";
        REQUIRE(aln_test::write_test_fasta(prefix));
        REQUIRE(bwa_idx_build(prefix.c_str(), prefix.c_str(), 0, 3) == 0);
        fmi.reset(new FMI_search(prefix.c_str()));
        fmi->load_index(true, 1);
        bwt = ba::bwt_from_fmi(fmi.get());
    }
    ~BamIndex()
    {
        ba::bwt_destroy(bwt);
        fmi.reset();
        std::string cmd = "rm -rf '" + dir + "'";
        if (system(cmd.c_str()) != 0) std::fprintf(stderr, "warning: could not remove %s\n", dir.c_str());
    }
    const bntseq_t *bns() const { return fmi->idx->bns; }
};

BamIndex &index_fixture()
{
    static BamIndex idx;
    return idx;
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

} // namespace

TEST_CASE("aln_se: BAM output holds the reads, mapped where the SAM puts them"
          * doctest::test_suite("unit/aln_se")) {
    BamIndex &ix = index_fixture();
    const std::string c1 = aln_test::test_contigs()[0].second;
    std::string fq = write_fastq(ix.dir, "bam.fq", {c1.substr(450, 32), c1.substr(10, 30), "ACGTACGTACGTACGTACGT"});
    std::string bam = ix.dir + "/out.bam";

    ba::aln_opt_t o;
    ba::aln_opt_init(&o);
    const char *hdr = "@HD\tVN:1.5\tSO:unsorted\tGO:query\n@SQ\tSN:c1\tLN:600\n@SQ\tSN:c2\tLN:550\n@SQ\tSN:c3\tLN:128\n";
    ba::aln_out_t *out = ba::aln_out_open(bam.c_str(), 1, 1, hdr);
    REQUIRE(out != nullptr);
    REQUIRE(ba::aln_run_se(&o, ix.bwt, ix.bns(), ix.fmi->idx->pac, fq.c_str(), out) == 0);
    REQUIRE(ba::aln_out_close(out) == 0);

    htsFile *fp = hts_open(bam.c_str(), "r");
    REQUIRE(fp != nullptr);
    sam_hdr_t *h = sam_hdr_read(fp);
    REQUIRE(h != nullptr);
    CHECK(sam_hdr_nref(h) == 3);
    bam1_t *b = bam_init1();
    int n = 0;
    std::vector<std::string> names;
    while (sam_read1(fp, h, b) >= 0) {
        names.push_back(bam_get_qname(b));
        if (n == 0) {
            CHECK(b->core.tid == 0);
            CHECK(b->core.pos == 450);
            CHECK(b->core.l_qseq == 32);
        }
        ++n;
    }
    CHECK(n == 3);
    CHECK(names == std::vector<std::string>({"r0", "r1", "r2"}));
    bam_destroy1(b);
    sam_hdr_destroy(h);
    hts_close(fp);
}
