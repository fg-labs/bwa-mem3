// Unit tests for the paired-end half of `bwa-mem3 aln` (src/aln/bwape.cpp and aln_run_pe in
// aln_main.cpp). Byte-for-byte agreement with stock bwa aln + sampe is
// test/regression/aln_parity.sh; these cases pin two behaviours the parity script reaches only
// indirectly: the carry-over of the insert-size estimate across batches (in upstream's
// bwa_cal_pac_pos_pe, which the stages are checked against, and in the shipped driver
// aln_run_pe over several batches), and the per-file batch used for the gap-open cap.

#include "doctest/doctest.h"

#include "aln_test_ref.h"
#include "../../src/FMI_search.h"
#include "../../src/aln/aln_main.h"
#include "../../src/aln/bwape.h"
#include "../../src/aln/bwase.h"
#include "../../src/aln/bwaseqio.h"
#include "../../src/aln/bwtgap.h"
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

// One 6,000 bp random contig: every read is unique, so well-placed pairs get mapQ >= 20.
struct PeIndex {
    bwa_tests::TempDir tmp_dir{"aln_pe_unit"}; // removed only by the process that made it
    std::string dir = tmp_dir.path(), prefix, contig;
    std::unique_ptr<FMI_search> fmi;
    ba::bwt_t *bwt = nullptr;

    PeIndex()
    {
        prefix = dir + "/" + dir.substr(dir.rfind('/') + 1) + ".fa";
        uint32_t state = 77;
        contig = aln_test::lcg_dna(state, 6000);
        FILE *fp = fopen(prefix.c_str(), "w");
        REQUIRE(fp != nullptr);
        bool ok = fprintf(fp, ">u1\n") >= 0;
        for (size_t i = 0; i < contig.size(); i += 60) ok = ok && fprintf(fp, "%s\n", contig.substr(i, 60).c_str()) >= 0;
        REQUIRE((fclose(fp) == 0 && ok));
        REQUIRE(bwa_idx_build(prefix.c_str(), prefix.c_str(), 0, 3) == 0);
        fmi.reset(new FMI_search(prefix.c_str()));
        fmi->load_index(true, 1, /*build_seed_jump=*/false);  // as main_aln loads it
        bwt = ba::bwt_from_fmi(fmi.get());
        tmp_dir.remove(); // loaded: the cases need only the in-memory index
    }
    ~PeIndex()
    {
        ba::bwt_destroy(bwt);
        fmi.reset();
    }
    const bntseq_t *bns() const { return fmi->idx->bns; }
};

PeIndex &index_fixture()
{
    static PeIndex idx;
    return idx;
}

std::string write_fastq(const std::string &path, const std::vector<std::string> &reads)
{
    FILE *fp = fopen(path.c_str(), "w");
    REQUIRE(fp != nullptr);
    bool ok = true;
    for (size_t i = 0; i < reads.size(); ++i)
        ok = ok && fprintf(fp, "@p%zu\n%s\n+\n%s\n", i, reads[i].c_str(), std::string(reads[i].size(), 'I').c_str()) >= 0;
    REQUIRE((fclose(fp) == 0 && ok));
    return path;
}

// n pairs of `len` bp from fragments of 280-320 bp: read 1 forward at the start, read 2
// reverse-complemented at the end.
void well_placed(const std::string &contig, int n, int len, int first, std::vector<std::string> *r1,
                 std::vector<std::string> *r2)
{
    for (int i = 0; i < n; ++i) {
        int frag = 280 + (i * 7) % 41;
        int pos = first + i * 20;
        r1->push_back(contig.substr(pos, len));
        r2->push_back(aln_test::revcomp(contig.substr(pos + frag - len, len)));
    }
}

void random_pairs(uint32_t seed, int n, int len, std::vector<std::string> *r1, std::vector<std::string> *r2)
{
    uint32_t state = seed;
    for (int i = 0; i < n; ++i) {
        r1->push_back(aln_test::lcg_dna(state, len));
        r2->push_back(aln_test::lcg_dna(state, len));
    }
}

// Reads one batch of each file, searches, and runs the pairing stage.
struct Batch {
    ba::bwa_seq_t *seqs[2] = {nullptr, nullptr};
    int n = 0;
    ~Batch() { for (int j = 0; j < 2; ++j) if (seqs[j]) ba::bwa_free_read_seq(n, seqs[j]); }
};

void load_batch(const ba::gap_opt_t &gap, const std::string &f1, const std::string &f2, Batch *b)
{
    ba::bwa_seqio_t *ks1 = ba::bwa_seq_open(f1.c_str()), *ks2 = ba::bwa_seq_open(f2.c_str());
    REQUIRE(ks1 != nullptr);
    REQUIRE(ks2 != nullptr);
    int n1 = 0, n2 = 0;
    b->seqs[0] = ba::bwa_read_seq(ks1, 0x40000, &n1, gap.mode, gap.trim_qual);
    b->seqs[1] = ba::bwa_read_seq(ks2, 0x40000, &n2, gap.mode, gap.trim_qual);
    REQUIRE(n1 == n2);
    b->n = n1;
    ba::bwa_seq_close(ks1);
    ba::bwa_seq_close(ks2);
    PeIndex &ix = index_fixture();
    ba::aln_search_batch(ix.bwt, b->n, b->seqs[0], &gap, 1);
    ba::aln_search_batch(ix.bwt, b->n, b->seqs[1], &gap, 1);
}

void check_same_ii(const ba::isize_info_t &a, const ba::isize_info_t &b)
{
    CHECK(a.avg == b.avg);
    CHECK(a.std == b.std);
    CHECK(a.low == b.low);
    CHECK(a.high == b.high);
    CHECK(a.high_bayesian == b.high_bayesian);
}

} // namespace

TEST_CASE("aln_pe: infer_isize failure falls back to the previous batch's estimate"
          * doctest::test_suite("unit/aln_pe")) {
    PeIndex &ix = index_fixture();
    bwa_tests::TempDir scratch{"aln_pe_case"}; // this case's files, removed when it ends
    ba::bwase_initialize();
    srand48(ix.bns()->seed);
    ba::gap_opt_t *gap = ba::gap_init_opt();
    ba::pe_opt_t *pe = ba::bwa_init_pe_opt();
    std::vector<std::string> a1, a2, b1, b2, c1, c2;
    well_placed(ix.contig, 200, 36, 100, &a1, &a2);   // batch 1
    random_pairs(5, 200, 36, &b1, &b2);               // batch 2: unmappable
    well_placed(ix.contig, 200, 36, 1000, &c1, &c2);  // batch 3
    ba::isize_info_t last;
    last.avg = -1.0;

    Batch bt1;
    load_batch(*gap, write_fastq(scratch.path() + "/a1.fq", a1), write_fastq(scratch.path() + "/a2.fq", a2), &bt1);
    ba::isize_info_t ii1;
    ba::bwa_cal_pac_pos_pe(ix.bns(), ix.bwt, bt1.n, bt1.seqs, &ii1, pe, gap, &last);
    REQUIRE(ii1.avg > 0.0);
    CHECK(ii1.avg > 280.0);
    CHECK(ii1.avg < 320.0);
    CHECK(ii1.high_bayesian > ii1.low);

    // batch 2: its own estimate fails, so the pairing uses batch 1's
    Batch bt2;
    load_batch(*gap, write_fastq(scratch.path() + "/b1.fq", b1), write_fastq(scratch.path() + "/b2.fq", b2), &bt2);
    ba::isize_info_t ii2;
    ba::bwa_cal_pac_pos_pe(ix.bns(), ix.bwt, bt2.n, bt2.seqs, &ii2, pe, gap, &ii1);
    {
        ba::isize_info_t own;
        CHECK(ba::infer_isize(bt2.n, bt2.seqs, &own, pe->ap_prior, ix.bwt->seq_len / 2) == -1);
        CHECK(own.avg < 0.0);
    }
    check_same_ii(ii2, ii1);
    CHECK(ii2.ap_prior == ii1.ap_prior);

    // batch 3: a fresh estimate again, from its own pairs
    Batch bt3;
    load_batch(*gap, write_fastq(scratch.path() + "/c1.fq", c1), write_fastq(scratch.path() + "/c2.fq", c2), &bt3);
    ba::isize_info_t ii3;
    ba::bwa_cal_pac_pos_pe(ix.bns(), ix.bwt, bt3.n, bt3.seqs, &ii3, pe, gap, &ii2);
    CHECK(ii3.avg > 280.0);
    CHECK(ii3.avg < 320.0);

    // without a previous estimate batch 2 stays at "no estimate"
    Batch bt2b;
    load_batch(*gap, scratch.path() + "/b1.fq", scratch.path() + "/b2.fq", &bt2b);
    ba::isize_info_t none, ii2b;
    none.avg = -1.0;
    ba::bwa_cal_pac_pos_pe(ix.bns(), ix.bwt, bt2b.n, bt2b.seqs, &ii2b, pe, gap, &none);
    CHECK(ii2b.avg < 0.0);
    free(gap);
    free(pe);
}

TEST_CASE("aln_pe: the gap-open cap is per file (helper values) and mixed read lengths run end to end"
          * doctest::test_suite("unit/aln_pe")) {
    PeIndex &ix = index_fixture();
    bwa_tests::TempDir scratch{"aln_pe_case"}; // this case's files, removed when it ends
    ba::gap_opt_t *gap = ba::gap_init_opt();
    gap->max_gapo = 5; // the cap is min(max_gapo, differences allowed for the longest read)
    std::vector<std::string> r1, r2, long1, long2;
    well_placed(ix.contig, 30, 36, 100, &r1, &r2);          // R1: 36 bp reads
    well_placed(ix.contig, 30, 150, 100, &long1, &long2);   // R2: 150 bp reads
    const std::string f1 = write_fastq(scratch.path() + "/m1.fq", r1);
    const std::string f2 = write_fastq(scratch.path() + "/m2.fq", long2);
    Batch bt;
    {
        ba::bwa_seqio_t *ks1 = ba::bwa_seq_open(f1.c_str()), *ks2 = ba::bwa_seq_open(f2.c_str());
        REQUIRE(ks1 != nullptr);
        REQUIRE(ks2 != nullptr);
        int n1 = 0, n2 = 0;
        bt.seqs[0] = ba::bwa_read_seq(ks1, 0x40000, &n1, gap->mode, gap->trim_qual);
        bt.seqs[1] = ba::bwa_read_seq(ks2, 0x40000, &n2, gap->mode, gap->trim_qual);
        REQUIRE(n1 == 30);
        REQUIRE(n2 == 30);
        bt.n = 30;
        ba::bwa_seq_close(ks1);
        ba::bwa_seq_close(ks2);
    }
    const int cap1 = ba::aln_batch_max_gapo(bt.n, bt.seqs[0], gap);
    const int cap2 = ba::aln_batch_max_gapo(bt.n, bt.seqs[1], gap);
    CHECK(cap1 == ba::bwa_cal_maxdiff(36, BWA_AVG_ERR, gap->fnr));
    CHECK(cap2 == 5);
    CHECK(cap1 < cap2); // searching the two files as one batch would give R1 the larger cap
    // NOTE: this checks the helper's values and that the driver runs on files of different read
    // lengths. It cannot tell a per-file search from a combined one: the cap never binds in
    // the search (a gap open is only tried from an entry with differences still allowed, so the
    // entry's gap opens stay below the file's own cap), so the output is the same either way.
    // Output parity for this shape is the pe_mixlen case in aln_parity.sh.

    // end to end: the driver runs and writes both mates of every pair
    ba::aln_opt_t o;
    ba::aln_opt_init(&o);
    o.gap.max_gapo = 5;
    const std::string out_path = scratch.path() + "/pe.sam";
    ba::aln_out_t *out = ba::aln_out_open(out_path.c_str(), 0, 0, "");
    REQUIRE(out != nullptr);
    REQUIRE(ba::aln_run_pe(&o, ix.bwt, ix.bns(), ix.fmi->idx->pac, f1.c_str(), f2.c_str(), out) == 0);
    REQUIRE(ba::aln_out_close(out) == 0);
    std::ifstream in(out_path);
    int n_rec = 0;
    std::vector<std::string> names;
    std::vector<int> flags;
    for (std::string l; std::getline(in, l);) {
        if (l.empty() || l[0] == '@') continue;
        ++n_rec;
        std::istringstream f(l);
        std::string name;
        int flag = 0;
        f >> name >> flag;
        names.push_back(name);
        flags.push_back(flag);
    }
    REQUIRE(n_rec == 60);
    for (int i = 0; i < 30; ++i) { // pair i is p<i> as read 1, then p<i> as read 2
        CAPTURE(i);
        const std::string want = "p" + std::to_string(i);
        CHECK(names[2 * i] == want);
        CHECK(names[2 * i + 1] == want);
        CHECK((flags[2 * i] & 0xc1) == 0x41);     // paired, first in pair
        CHECK((flags[2 * i + 1] & 0xc1) == 0x81); // paired, second in pair
    }
    free(gap);
}

namespace {

// Runs aln_run_pe (the shipped driver) and returns the FLAG of every record, in output order.
std::vector<int> run_pe_flags(int batch_size, const std::string &f1, const std::string &f2, const std::string &out_path)
{
    PeIndex &ix = index_fixture();
    ba::aln_opt_t o;
    ba::aln_opt_init(&o);
    o.batch_size = batch_size;
    ba::aln_out_t *out = ba::aln_out_open(out_path.c_str(), 0, 0, "");
    REQUIRE(out != nullptr);
    REQUIRE(ba::aln_run_pe(&o, ix.bwt, ix.bns(), ix.fmi->idx->pac, f1.c_str(), f2.c_str(), out) == 0);
    REQUIRE(ba::aln_out_close(out) == 0);
    std::ifstream in(out_path);
    std::vector<int> flags;
    for (std::string l; std::getline(in, l);) {
        if (l.empty() || l[0] == '@') continue;
        size_t t = l.find('\t');
        REQUIRE(t != std::string::npos);
        flags.push_back(std::atoi(l.c_str() + t + 1));
    }
    return flags;
}

} // namespace

TEST_CASE("aln_pe: aln_run_pe carries the insert-size estimate into a batch too small for its own"
          * doctest::test_suite("unit/aln_pe")) {
    // 200 well-placed pairs (fragments of 280-320 bp), then 5 pairs from 470 bp fragments, in
    // batches of 200. The 5-pair batch is below infer_isize's minimum of 20 pairs, so it has no
    // estimate of its own. With the first batch's estimate carried over (aln_run_pe's last_ii),
    // 470 bp is far above its upper bound and the pairs are not proper (no 0x2). Without an
    // estimate the bound is -a (500), so the same 5 pairs on their own are proper: the flag can
    // only come from the carry-over.
    PeIndex &ix = index_fixture();
    bwa_tests::TempDir scratch{"aln_pe_case"}; // this case's files, removed when it ends
    std::vector<std::string> r1, r2;
    well_placed(ix.contig, 200, 36, 100, &r1, &r2);
    std::vector<std::string> t1, t2;
    for (int i = 0; i < 5; ++i) {
        const int pos = 4500 + i * 40, frag = 470;
        t1.push_back(ix.contig.substr(pos, 36));
        t2.push_back(aln_test::revcomp(ix.contig.substr(pos + frag - 36, 36)));
    }
    std::vector<std::string> a1(r1), a2(r2);
    a1.insert(a1.end(), t1.begin(), t1.end());
    a2.insert(a2.end(), t2.begin(), t2.end());
    const std::string f1 = write_fastq(scratch.path() + "/carry1.fq", a1), f2 = write_fastq(scratch.path() + "/carry2.fq", a2);

    std::vector<int> alone = run_pe_flags(200, write_fastq(scratch.path() + "/tail1.fq", t1), write_fastq(scratch.path() + "/tail2.fq", t2),
                                          scratch.path() + "/tail.sam");
    REQUIRE(alone.size() == 10);
    for (int f : alone) CHECK((f & 0x2) != 0);

    std::vector<int> two = run_pe_flags(200, f1, f2, scratch.path() + "/carry_200.sam");
    REQUIRE(two.size() == 410);
    for (size_t i = 400; i < 410; ++i) {
        CAPTURE(i);
        CHECK((two[i] & 0x2) == 0);
    }
    // the first batch is flagged as in a single-batch run of the same pairs
    std::vector<int> one = run_pe_flags(0x40000, f1, f2, scratch.path() + "/carry_one.sam");
    REQUIRE(one.size() == 410);
    for (size_t i = 0; i < 400; ++i) CHECK((one[i] & 0x2) == (two[i] & 0x2));

    // and each first-batch record is what its construction implies, independent of either run:
    // the reads are exact copies of a unique contig, so R1 is forward at pos, R2 reverse at
    // pos + frag - 36, both 36M, a proper pair with TLEN +-frag (well_placed's arithmetic)
    std::ifstream in(scratch.path() + "/carry_200.sam");
    std::vector<std::vector<std::string>> recs;
    for (std::string l; std::getline(in, l);) {
        if (l.empty() || l[0] == '@') continue;
        std::vector<std::string> f;
        std::istringstream ls(l);
        for (std::string x; std::getline(ls, x, '\t');) f.push_back(x);
        REQUIRE(f.size() >= 11);
        recs.push_back(f);
    }
    REQUIRE(recs.size() == 410);
    for (int i = 0; i < 200; ++i) {
        const int frag = 280 + (i * 7) % 41, pos = 100 + i * 20; // well_placed(..., 100, ...)
        const int p1 = pos + 1, p2 = pos + frag - 36 + 1;
        const std::vector<std::string> &a = recs[2 * i], &b = recs[2 * i + 1];
        CAPTURE(i);
        CHECK(a[1] == "99");
        CHECK(b[1] == "147");
        CHECK(a[2] == "u1");
        CHECK(b[2] == "u1");
        CHECK(a[3] == std::to_string(p1));
        CHECK(b[3] == std::to_string(p2));
        CHECK(std::atoi(a[4].c_str()) >= 20);
        CHECK(std::atoi(b[4].c_str()) >= 20);
        CHECK(a[5] == "36M");
        CHECK(b[5] == "36M");
        CHECK(a[6] == "=");
        CHECK(b[6] == "=");
        CHECK(a[7] == std::to_string(p2));
        CHECK(b[7] == std::to_string(p1));
        CHECK(a[8] == std::to_string(frag));
        CHECK(b[8] == std::to_string(-frag));
        CHECK(a[9] == r1[i]);
        CHECK(b[9] == ix.contig.substr(pos + frag - 36, 36)); // reverse strand: SEQ is the reference strand
    }
    // the tail pairs are still mapped where they were built, not lost (an unmapped pair would also
    // lack 0x2): R1 forward at pos, R2 reverse at pos + 470 - 36, both 36M, mates on the same
    // contig with TLEN +-470, flagged paired and mate-strand but not proper (97 and 145)
    for (int i = 0; i < 5; ++i) {
        const int pos = 4500 + i * 40, frag = 470; // the tail loop above
        const int p1 = pos + 1, p2 = pos + frag - 36 + 1;
        const std::vector<std::string> &a = recs[400 + 2 * i], &b = recs[401 + 2 * i];
        CAPTURE(i);
        CHECK(a[1] == "97");
        CHECK(b[1] == "145");
        CHECK(a[2] == "u1");
        CHECK(b[2] == "u1");
        CHECK(a[3] == std::to_string(p1));
        CHECK(b[3] == std::to_string(p2));
        CHECK(std::atoi(a[4].c_str()) >= 20);
        CHECK(std::atoi(b[4].c_str()) >= 20);
        CHECK(a[5] == "36M");
        CHECK(b[5] == "36M");
        CHECK(a[6] == "=");
        CHECK(b[6] == "=");
        CHECK(a[7] == std::to_string(p2));
        CHECK(b[7] == std::to_string(p1));
        CHECK(a[8] == std::to_string(frag));
        CHECK(b[8] == std::to_string(-frag));
        CHECK(a[9] == t1[i]);
        CHECK(b[9] == ix.contig.substr(pos + frag - 36, 36));
    }
}

TEST_CASE("aln_pe: aln_run_pe bounds the joined barcode to BWA_MAX_BCLEN"
          * doctest::test_suite("unit/aln_pe")) {
    // main_aln limits -B to BWA_MAX_BCLEN / 2 with two read files, but aln_run_pe is a public
    // entry point and bwa_read_seq accepts barcodes up to BWA_MAX_BCLEN. Two 40 bp barcodes
    // join to 80 bases, more than the 63 the buffer holds: the join keeps the first 63.
    PeIndex &ix = index_fixture();
    bwa_tests::TempDir scratch{"aln_pe_case"}; // this case's files, removed when it ends
    std::vector<std::string> r1, r2;
    well_placed(ix.contig, 5, 120, 100, &r1, &r2);
    const std::string f1 = write_fastq(scratch.path() + "/bc1.fq", r1);
    const std::string f2 = write_fastq(scratch.path() + "/bc2.fq", r2);
    ba::aln_opt_t o;
    ba::aln_opt_init(&o);
    const int l_bc = 40;
    REQUIRE(2 * l_bc > BWA_MAX_BCLEN);
    o.gap.mode |= l_bc << 24;
    const std::string out_path = scratch.path() + "/bc.sam";
    ba::aln_out_t *out = ba::aln_out_open(out_path.c_str(), 0, 0, "");
    REQUIRE(out != nullptr);
    REQUIRE(ba::aln_run_pe(&o, ix.bwt, ix.bns(), ix.fmi->idx->pac, f1.c_str(), f2.c_str(), out) == 0);
    REQUIRE(ba::aln_out_close(out) == 0);
    std::ifstream in(out_path);
    int n_rec = 0;
    for (std::string l; std::getline(in, l);) {
        if (l.empty() || l[0] == '@') continue;
        const int pair = n_rec / 2;
        ++n_rec;
        CAPTURE(pair);
        const size_t at = l.find("\tBC:Z:");
        REQUIRE(at != std::string::npos);
        const size_t st = at + 6, en = l.find('\t', st);
        const std::string bc = l.substr(st, en == std::string::npos ? std::string::npos : en - st);
        const std::string want = (r1[pair].substr(0, l_bc) + r2[pair].substr(0, l_bc)).substr(0, BWA_MAX_BCLEN);
        CHECK(bc == want); // both mates carry the same joined barcode
    }
    CHECK(n_rec == 10);
}
