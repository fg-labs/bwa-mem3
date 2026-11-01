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
#include <cstring>
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

namespace {

// Every record of `path` (SAM or BAM), as htslib reads it.
std::vector<bam1_t *> read_all(const std::string &path)
{
    std::vector<bam1_t *> v;
    htsFile *fp = hts_open(path.c_str(), "r");
    REQUIRE(fp != nullptr);
    sam_hdr_t *h = sam_hdr_read(fp);
    REQUIRE(h != nullptr);
    for (bam1_t *b = bam_init1(); sam_read1(fp, h, b) >= 0; b = bam_init1()) v.push_back(b);
    sam_hdr_destroy(h);
    hts_close(fp);
    return v;
}

void free_all(std::vector<bam1_t *> &v)
{
    for (bam1_t *b : v) bam_destroy1(b);
    v.clear();
}

const char *kHdr = "@HD\tVN:1.5\tSO:unsorted\tGO:query\n@SQ\tSN:c1\tLN:600\n@SQ\tSN:c2\tLN:550\n@SQ\tSN:c3\tLN:128\n";

int run_to(const std::string &path, int bam, const char *hdr, const std::string &fq)
{
    BamIndex &ix = index_fixture();
    ba::aln_opt_t o;
    ba::aln_opt_init(&o);
    ba::aln_out_t *out = ba::aln_out_open(path.c_str(), bam, 1, hdr);
    REQUIRE(out != nullptr);
    int ret = ba::aln_run_se(&o, ix.bwt, ix.bns(), ix.fmi->idx->pac, fq.c_str(), out);
    if (ba::aln_out_close(out) != 0) ret = 1;
    return ret;
}

} // namespace

TEST_CASE("aln_se: BAM output holds the same records as the SAM output"
          * doctest::test_suite("unit/aln_se")) {
    BamIndex &ix = index_fixture();
    const std::string c1 = aln_test::test_contigs()[0].second;
    // mapped forward, mapped reverse with a mismatch (MD, NM), in the c1/c2 repeat (XA), unmapped
    std::string rev = aln_test::revcomp(c1.substr(470, 36));
    rev[10] = rev[10] == 'A' ? 'C' : 'A';
    std::string fq = write_fastq(ix.dir, "bam.fq",
                                 {c1.substr(450, 32), c1.substr(10, 30), "ACGTACGTACGTACGTACGT", rev, c1.substr(150, 36)});
    const std::string bam = ix.dir + "/out.bam", sam = ix.dir + "/out.sam";
    REQUIRE(run_to(bam, 1, kHdr, fq) == 0);
    REQUIRE(run_to(sam, 0, kHdr, fq) == 0);

    std::vector<bam1_t *> got = read_all(bam), want = read_all(sam);
    REQUIRE(got.size() == 5);
    REQUIRE(want.size() == got.size());
    for (size_t i = 0; i < got.size(); ++i) {
        CAPTURE(i);
        const bam1_core_t &g = got[i]->core, &w = want[i]->core;
        CHECK(g.tid == w.tid);
        CHECK(g.pos == w.pos);
        CHECK(g.flag == w.flag);
        CHECK(g.qual == w.qual);
        CHECK(g.mtid == w.mtid);
        CHECK(g.mpos == w.mpos);
        CHECK(g.isize == w.isize);
        CHECK(g.n_cigar == w.n_cigar);
        CHECK(g.l_qseq == w.l_qseq);
        // qname, CIGAR, seq, qual and every aux tag, byte for byte
        REQUIRE(got[i]->l_data == want[i]->l_data);
        CHECK(std::memcmp(got[i]->data, want[i]->data, got[i]->l_data) == 0);
    }
    // and the records are where they should be
    CHECK(got[0]->core.tid == 0);
    CHECK(got[0]->core.pos == 450);
    CHECK(got[1]->core.tid == 0);
    CHECK(got[1]->core.pos == 10);
    CHECK((got[2]->core.flag & BAM_FUNMAP) != 0);
    CHECK((got[3]->core.flag & BAM_FREVERSE) != 0);
    CHECK(bam_aux_get(got[3], "MD") != nullptr);
    CHECK(bam_aux_get(got[4], "XA") != nullptr);
    free_all(got);
    free_all(want);
}

TEST_CASE("aln_se: BAM output refuses a record whose reference is not in the header"
          * doctest::test_suite("unit/aln_se")) {
    // An -r text with its own @SQ lines replaces the index's (stock samse does the same); SAM keeps
    // the reference name as text, but sam_parse1 would turn it into tid -1 and mark the read
    // unmapped, so the BAM path stops with an error instead of writing a corrupt record.
    BamIndex &ix = index_fixture();
    const std::string c1 = aln_test::test_contigs()[0].second;
    std::string fq = write_fastq(ix.dir, "badsq.fq", {c1.substr(450, 32)});
    const char *hdr = "@HD\tVN:1.5\tSO:unsorted\tGO:query\n@SQ\tSN:foo\tLN:10\n";
    CHECK(run_to(ix.dir + "/badsq.bam", 1, hdr, fq) == 1);
    CHECK(run_to(ix.dir + "/badsq.sam", 0, hdr, fq) == 0); // SAM is stock's output
}
