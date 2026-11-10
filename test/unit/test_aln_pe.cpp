// Unit tests for the paired-end half of `bwa-mem3 aln` (src/aln/bwape.cpp and aln_run_pe in
// aln_main.cpp). Byte-for-byte agreement with stock bwa aln + sampe is
// test/regression/aln_parity.sh; these cases pin two behaviours the parity script reaches only
// indirectly: the carry-over of the insert-size estimate across batches (in upstream's
// bwa_cal_pac_pos_pe, which the stages are checked against, and in the shipped driver
// aln_run_pe over several batches), and the per-file batch used for the gap-open cap.

#include "doctest/doctest.h"

#include "aln_test_gz.h"
#include "aln_test_ref.h"
#include "temp_index.h"
#include "../../src/FMI_search.h"
#include "../../src/aln/aln_main.h"
#include "../../src/aln/bwape.h"
#include "../../src/aln/bwase.h"
#include "../../src/aln/bwaseqio.h"
#include "../../src/aln/bwtgap.h"
#include "../../src/bwa.h"
#include "../../src/kvec.h"
#include "../../src/read_index_ele.h"

#include <cmath>
#include <cstddef>
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

// Reads are named p0, p1, ...; read `renamed` (if >= 0) is named x<i> instead.
std::string write_fastq(const std::string &path, const std::vector<std::string> &reads, int renamed = -1)
{
    FILE *fp = fopen(path.c_str(), "w");
    REQUIRE(fp != nullptr);
    bool ok = true;
    for (size_t i = 0; i < reads.size(); ++i)
        ok = ok && fprintf(fp, "@%s%zu\n%s\n+\n%s\n", (int)i == renamed ? "x" : "p", i, reads[i].c_str(),
                           std::string(reads[i].size(), 'I').c_str()) >= 0;
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

TEST_CASE("aln_pe: short-insert preserves the stock pairing option prefix"
          * doctest::test_suite("unit/aln_pe")) {
    // Existing callers' field displacements must not change when an opt-in field is added.
    struct StockOptions {
        int max_isize, force_isize;
        int max_occ;
        int n_multi, N_multi;
        int type, is_sw, is_preload;
        double ap_prior;
    };
    CHECK(offsetof(ba::pe_opt_t,max_isize)==offsetof(StockOptions,max_isize));
    CHECK(offsetof(ba::pe_opt_t,force_isize)==offsetof(StockOptions,force_isize));
    CHECK(offsetof(ba::pe_opt_t,max_occ)==offsetof(StockOptions,max_occ));
    CHECK(offsetof(ba::pe_opt_t,n_multi)==offsetof(StockOptions,n_multi));
    CHECK(offsetof(ba::pe_opt_t,N_multi)==offsetof(StockOptions,N_multi));
    CHECK(offsetof(ba::pe_opt_t,type)==offsetof(StockOptions,type));
    CHECK(offsetof(ba::pe_opt_t,is_sw)==offsetof(StockOptions,is_sw));
    CHECK(offsetof(ba::pe_opt_t,is_preload)==offsetof(StockOptions,is_preload));
    CHECK(offsetof(ba::pe_opt_t,ap_prior)==offsetof(StockOptions,ap_prior));
    CHECK(offsetof(ba::pe_opt_t,short_insert)>=sizeof(StockOptions));
    auto *opt=ba::bwa_init_pe_opt();
    REQUIRE(opt!=nullptr);
    CHECK(opt->short_insert==0);
    free(opt);
}

TEST_CASE("aln_pe: short-insert estimates all short spans despite a long mapped or unmapped outlier"
          * doctest::test_suite("unit/aln_pe")) {
    std::vector<ba::bwa_seq_t> ends[2];
    for (auto &end : ends) end.resize(41);
    for (int i = 0; i < 41; ++i) for (int j = 0; j < 2; ++j) {
        auto &p = ends[j][i];
        p.mapQ = 37;
        p.pos = 100 + 200 * i + 5 * j;
        p.len = p.full_len = i == 40 ? 150 : 20 + i % 16;
    }
    ba::bwa_seq_t *seqs[2] = {ends[0].data(), ends[1].data()};
    for (int n : {40, 41}) for (bool mapped : {false, true}) {
        CAPTURE(n); CAPTURE(mapped);
        ends[0][40].mapQ = ends[1][40].mapQ = mapped ? 37 : 0;
        ba::isize_info_t ii = {};
        CHECK(ba::infer_isize(n, seqs, &ii, 1e-5, 100000, nullptr, 1, true) == 0);
        CHECK(ii.low == 12);
        CHECK(ii.high == 52);
        CHECK(ii.avg == doctest::Approx(31.7));
        CHECK(ii.std == doctest::Approx(std::sqrt(823.4 / 40)));
        CHECK(ii.high_bayesian > 40);
        CHECK(ii.high_bayesian < 155);
    }
}

TEST_CASE("aln_pe: short-insert clamps the raw lower fence and normalizes failed estimates"
          * doctest::test_suite("unit/aln_pe")) {
    std::vector<ba::bwa_seq_t> ends[2];
    for (auto &end : ends) end.resize(20);
    ba::bwa_seq_t *seqs[2] = {ends[0].data(), ends[1].data()};
    for (int i = 0; i < 20; ++i) for (int j = 0; j < 2; ++j) {
        auto &p = ends[j][i];
        p.mapQ = 37; p.len = p.full_len = 20;
        p.pos = 100 + 100 * i + j * (i < 10 ? 1 : 30);
    }
    ba::isize_info_t ii = {};
    CHECK(ba::infer_isize(20, seqs, &ii, 1e-5, 100000, nullptr, 1, true) == 0);
    CHECK(ii.low == 1);
    CHECK(ii.high == 108);
    CHECK(ii.avg == doctest::Approx(35.5));
    CHECK(ba::infer_isize(19, seqs, &ii, 1e-5, 100000, nullptr, 1, true) == -1);
    CHECK(ii.avg == -1); CHECK(ii.std == -1);
    CHECK(ii.low == 0); CHECK(ii.high == 0); CHECK(ii.high_bayesian == 0);
    // Equal fences and constant spans retain the inherited invalid-variance failure.
    for (int i = 0; i < 20; ++i) ends[1][i].pos = ends[0][i].pos + 5;
    CHECK(ba::infer_isize(20, seqs, &ii, 1e-5, 100000, nullptr, 1, true) == -1);
    CHECK(ii.avg == -1); CHECK(ii.std == -1);
    CHECK(ii.low == 0); CHECK(ii.high == 0); CHECK(ii.high_bayesian == 0);
}

TEST_CASE("aln_pe: short-insert constant spans fail safely with an unmapped long outlier"
          * doctest::test_suite("unit/aln_pe")) {
    ba::bwa_seq_t ends[2][21] = {};
    for (int j = 0; j < 2; ++j) for (int i = 0; i < 21; ++i) {
        auto &p = ends[j][i];
        p.mapQ = i < 20 ? 37 : 0;
        p.len = p.full_len = i < 20 ? 20 : 150;
        p.pos = 100 + 200 * i + 5 * j;
    }
    ba::bwa_seq_t *seqs[2] = {ends[0], ends[1]};
    ba::isize_info_t ii = {};
    CHECK(ba::infer_isize(21, seqs, &ii, 1e-5, 100000) == -1);
    CHECK(ii.low == 150); CHECK(ii.high == 25);
    CHECK(ba::infer_isize(21, seqs, &ii, 1e-5, 100000, nullptr, 1, true) == -1);
    CHECK(ii.avg == -1); CHECK(ii.std == -1);
    CHECK(ii.low == 0); CHECK(ii.high == 0); CHECK(ii.high_bayesian == 0);
}

TEST_CASE("aln_pe: short-insert zero-span failure leaves no active estimate"
          * doctest::test_suite("unit/aln_pe")) {
    ba::bwa_seq_t ends[2][20] = {};
    for (auto &end : ends) for (auto &p : end) p.mapQ = 37;
    ba::bwa_seq_t *seqs[2] = {ends[0], ends[1]};
    ba::isize_info_t ii = {};
    CHECK(ba::infer_isize(20, seqs, &ii, 1e-5, 100000, nullptr, 1, true) == -1);
    CHECK(ii.avg == -1); CHECK(ii.std == -1);
    CHECK(ii.low == 0); CHECK(ii.high == 0); CHECK(ii.high_bayesian == 0);
}

TEST_CASE("aln_pe: short-insert retains previous-estimate fallback and explicit -A precedence"
          * doctest::test_suite("unit/aln_pe")) {
    auto *pe = ba::bwa_init_pe_opt();
    CHECK(pe->short_insert == 0);
    pe->short_insert = 1;
    ba::bwa_seq_t ends[2][20] = {};
    for (int i = 0; i < 20; ++i) for (int j = 0; j < 2; ++j) {
        auto &p = ends[j][i];
        p.mapQ = 37; p.len = p.full_len = 20;
        p.pos = 100 + 100 * i + j * (5 + i % 10);
    }
    ba::bwa_seq_t *seqs[2] = {ends[0], ends[1]};
    ba::isize_info_t none = {}, prior = {}, ii = {};
    none.avg = -1;
    ba::bwa_pe_isize_ref(100000, 20, seqs, &prior, pe, &none);
    REQUIRE(prior.avg > 0);
    // One fewer qualifying pair must reuse the complete prior estimate.
    ends[0][19].mapQ = 19;
    ba::bwa_pe_isize_ref(100000, 20, seqs, &ii, pe, &prior);
    check_same_ii(ii, prior); CHECK(ii.ap_prior == prior.ap_prior);
    ba::bwa_pe_isize_ref(100000, 20, seqs, &ii, pe, &none);
    CHECK(ii.avg == -1); CHECK(ii.std == -1);
    CHECK(ii.low == 0); CHECK(ii.high == 0); CHECK(ii.high_bayesian == 0);
    pe->force_isize = 1;
    ba::bwa_pe_isize_ref(100000, 20, seqs, &ii, pe, &prior);
    CHECK(ii.avg == -1); CHECK(ii.std == -1); CHECK(ii.low == 0); CHECK(ii.high == 0);
    ends[0][19].mapQ = 37;
    ba::bwa_pe_isize_ref(100000, 20, seqs, &ii, pe, &prior);
    CHECK(ii.avg == -1); CHECK(ii.std == -1); CHECK(ii.low == 0); CHECK(ii.high == 0);
    free(pe);
}

TEST_CASE("aln_pe: short-insert rejects underflowed pairing likelihoods without changing finite rankings"
          * doctest::test_suite("unit/aln_pe")) {
    ba::bwa_seq_t calibration[2][22] = {};
    for (int j = 0; j < 2; ++j) for (int i = 0; i < 22; ++i) {
        auto &p = calibration[j][i]; p.mapQ = i == 21 ? 0 : 37;
        p.len = p.full_len = i == 21 ? 150 : i == 20 ? 20 : 50;
        p.pos = 100 + 200 * i + (j ? (i == 20 ? 1 : 5 + i % 2) : 0);
    }
    ba::bwa_seq_t *cal[2] = {calibration[0], calibration[1]};
    ba::isize_info_t ii = {};
    REQUIRE(ba::infer_isize(22, cal, &ii, 1e-5, 100000, nullptr, 1, true) == 0);
    CHECK(ii.avg == doctest::Approx(55.5)); CHECK(ii.std == doctest::Approx(std::sqrt(.2)));
    CHECK(ii.high_bayesian == 58);
    // Independent log-tail reference: a 21bp span is >77 standard deviations away.
    const long double z = (ii.avg - 21) / ii.std;
    const long double tail_score = 4.343L * (z*z/2 + std::log(z) + .5L*std::log(2*acosl(-1.L)));
    CHECK(tail_score > 12000); CHECK(tail_score < 14000);
    CHECK(std::erfc((double)(z / std::sqrt(2.L))) == 0);
    auto *pe = ba::bwa_init_pe_opt(); pe->short_insert = 1; ba::bwase_initialize();
    for (bool alternative : {false, true}) {
        ba::bwa_seq_t ends[2] = {}; ba::bwt_aln1_t hits[2][2] = {};
        for (int j = 0; j < 2; ++j) {
            ends[j].aln = hits[j]; ends[j].pos = 10000 + j; ends[j].strand = j;
            ends[j].len = ends[j].full_len = 20; ends[j].mapQ = 25;
        }
        ba::bwa_seq_t *pair[2] = {ends, ends + 1};
        pair64_t rows[3] = {{10000, 0}, {10001, 3}, {10036, 7}};
        pair64_v arr = {}; arr.a = rows; arr.n = arr.m = alternative ? 3 : 2;
        ba::bwa_pe_pairing(pair, &arr, pe, 3, &ii);
        CHECK(ends[0].pos == 10000); CHECK(ends[1].pos == (alternative ? 10036 : 10001));
        CHECK(((ends[0].extra_flag & SAM_FPP) != 0) == alternative);
        CHECK(((ends[1].extra_flag & SAM_FPP) != 0) == alternative);
    }
    free(pe);
}

TEST_CASE("aln_pe: short-insert pairs a coincident start whatever the two hits' indices, unless -A"
          * doctest::test_suite("unit/aln_pe")) {
    // Rows sort by (position, hit index << 2 | strand << 1 | mate), so at one position the
    // reverse mate's hit 0 sorts before the forward mate's hit 1. The pair must still be found.
    auto *pe = ba::bwa_init_pe_opt();
    ba::bwase_initialize();
    ba::isize_info_t ii = {}; ii.low = 1; ii.high = ii.high_bayesian = 69; ii.avg = 32; ii.std = 7;
    for (int option : {0, 1, 2}) { // 2: --short-insert with -A, which overrides it
        CAPTURE(option);
        pe->short_insert = option != 0; pe->force_isize = option == 2;
        ba::bwa_seq_t ends[2] = {};
        ba::bwt_aln1_t hits0[2] = {}, hits1[1] = {};
        ends[0].aln = hits0; ends[1].aln = hits1;
        for (int j = 0; j < 2; ++j) {
            ends[j].type = BWA_TYPE_REPEAT; ends[j].strand = j;
            ends[j].len = ends[j].full_len = 20; ends[j].mapQ = ends[j].seQ = 25;
        }
        ends[0].pos = 5000; ends[1].pos = 100;
        ba::bwa_seq_t *pair[2] = {ends, ends + 1};
        pair64_t rows[3] = {{5000, 0u << 2 | 0 << 1 | 0}, {100, 1u << 2 | 0 << 1 | 0}, {100, 0u << 2 | 1 << 1 | 1}};
        pair64_v arr = {}; arr.a = rows; arr.n = arr.m = 3;
        ba::bwa_pe_pairing(pair, &arr, pe, 3, &ii);
        CHECK(ends[0].pos == (option == 1 ? 100 : 5000));
        for (const auto &p : ends) CHECK(((p.extra_flag & SAM_FPP) != 0) == (option == 1));
    }
    free(pe);
}

TEST_CASE("aln_pe: short-insert pairs coincident starts but keeps the other geometry and same-contig requirements"
          * doctest::test_suite("unit/aln_pe")) {
    auto *pe = ba::bwa_init_pe_opt();
    ba::bwase_initialize();
    ba::isize_info_t ii = {}; ii.avg = ii.std = -1;
    bntann1_t ann[2] = {}; ann[0].len = ann[1].len = 200; ann[1].offset = 200;
    bntseq_t bns = {}; bns.n_seqs = 2; bns.l_pac = 400; bns.anns = ann;
    // stock_proper: without the option. proper: with it, which only adds coincident starts.
    struct Geometry {int start, end, len0, len1, full0, full1; bool proper, stock_proper;};
    const Geometry cases[] = {
        {100, 100, 20, 20, 20, 20, true, false}, // coincident starts: both reads cover the fragment
        {100, 100, 20, 19, 20, 20, false, false}, // coincident, but the forward read runs past the reverse end
        {100, 100, 19, 20, 20, 20, true, false}, // coincident, the reverse read is the longer
        {100, 99, 20, 20, 20, 20, false, false}, // the reverse read starts first
        {100, 101, 20, 20, 20, 20, true, true},
        {100, 101, 10, 19, 20, 20, true, true}, // outer span == max full_len
        {100, 101, 10, 18, 20, 20, false, false},
        {100, 110, 30, 20, 30, 20, true, true}, // asymmetric full lengths
        {100, 109, 30, 20, 30, 20, false, false},
        {190, 205, 20, 20, 20, 20, true, true} // only ordinary path may pair across contigs
    };
    for (int option : {0, 1}) for (bool active : {false, true}) for (bool same_contig : {false, true})
    for (const auto &c : cases) {
        pe->short_insert = option;
        CAPTURE(option);
        ii = {}; ii.avg = ii.std = -1;
        if (active) { ii.low = 1; ii.high = ii.high_bayesian = 69; ii.avg = 32; ii.std = 7; }
        CAPTURE(active); CAPTURE(same_contig); CAPTURE(c.start); CAPTURE(c.end); CAPTURE(c.len1);
        ba::bwa_seq_t ends[2] = {};
        ba::bwt_aln1_t hits[2] = {};
        for (int j = 0; j < 2; ++j) {
            ends[j].aln = hits + j; ends[j].type = BWA_TYPE_UNIQUE;
            ends[j].pos = j ? c.end : c.start; ends[j].strand = j;
            ends[j].len = j ? c.len1 : c.len0; ends[j].full_len = j ? c.full1 : c.full0;
            ends[j].mapQ = ends[j].seQ = 25;
        }
        ba::bwa_seq_t *pair[2] = {ends, ends + 1};
        pair64_t rows[2] = {{(uint64_t)c.start, 0}, {(uint64_t)c.end, 3}};
        pair64_v arr = {}; arr.a = rows; arr.n = arr.m = 2;
        if (same_contig) ba::bwa_pe_pairing_same_contig(&bns, pair, &arr, pe, 3, &ii);
        else ba::bwa_pe_pairing(pair, &arr, pe, 3, &ii);
        const bool want = option ? c.proper : c.stock_proper;
        const bool proper = want && !(same_contig && c.start < 200 && c.end >= 200);
        for (const auto &p : ends) CHECK(((p.extra_flag & SAM_FPP) != 0) == proper);
    }
    free(pe);
}

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

// ---- SA lookups made in batches (bwt_sa_batch in bwa_pe_se_pos, bwa_pe_pair, bwa_cal_pac_pos) ----

namespace {

// A contig built to give reads many hits: random sequence with a 300 bp segment copied 20
// times forward and 10 times reverse-complemented, and a 160 bp "AC" tandem repeat. The index
// files are deleted as soon as it is loaded.
struct RepIndex {
    std::string contig;
    std::unique_ptr<bwa_tests::TempIndex> idx;
    ba::bwt_t *bwt = nullptr;

    RepIndex()
    {
        uint32_t state = 91;
        const std::string seg = aln_test::lcg_dna(state, 300);
        contig = aln_test::lcg_dna(state, 1500);
        for (int i = 0; i < 20; ++i) contig += seg + aln_test::lcg_dna(state, 400 + 37 * (i % 6));
        for (int i = 0; i < 10; ++i) contig += aln_test::revcomp(seg) + aln_test::lcg_dna(state, 500);
        for (int i = 0; i < 80; ++i) contig += "AC";
        contig += aln_test::lcg_dna(state, 1500);
        idx.reset(new bwa_tests::TempIndex({{"r1", contig}}));
        idx->drop_files();
        bwt = ba::bwt_from_fmi(idx->fmi());
    }
    ~RepIndex() { ba::bwt_destroy(bwt); }
    const bntseq_t *bns() const { return idx->bns(); }
};

RepIndex &rep_fixture()
{
    static RepIndex idx;
    return idx;
}

// Pairs of 36 bp (read 1) and 40 bp (read 2) from 260-340 bp fragments starting every 25 bp along
// the whole contig (so many land in the repeats). Every fifth read 2 has a mismatch, every third
// read 1 a 1 bp deletion and every fourth read 2 a 1 bp insertion, so the ends' reference lengths
// (len + ref_shift) differ between the ends and between hits; every seventh pair, the first
// included, is random (unmapped), and every eleventh has only read 1 placed. Written as
// rep1.fq and rep2.fq in `tmp`.
void rep_pairs(const bwa_tests::TempDir &tmp, std::string *f1, std::string *f2)
{
    RepIndex &ix = rep_fixture();
    std::vector<std::string> r1, r2;
    const int len1 = 36, len2 = 40;
    uint32_t state = 3;
    for (int pos = 0, i = 0; pos + 342 < (int)ix.contig.size(); pos += 25, ++i) {
        if (i % 7 == 0) {
            r1.push_back(aln_test::lcg_dna(state, len1));
            r2.push_back(aln_test::lcg_dna(state, len2));
        }
        const int frag = 260 + (i * 11) % 81;
        // read 1: with a deletion, 37 reference bases give 36 read bases (ref_shift +1)
        std::string a = ix.contig.substr(pos, len1 + (i % 3 == 0));
        if (i % 3 == 0) a.erase(len1 / 2, 1);
        r1.push_back(a);
        // read 2: with an insertion, 39 reference bases give 40 read bases (ref_shift -1)
        const int span2 = len2 - (i % 4 == 0);
        std::string m = ix.contig.substr(pos + frag - span2, span2);
        if (i % 4 == 0) m.insert((size_t)span2 / 2, 1, m[span2 / 2] == 'A' ? 'G' : 'A');
        m = aln_test::revcomp(m);
        if (i % 5 == 0) m[len2 / 4] = m[len2 / 4] == 'A' ? 'C' : 'A';
        r2.push_back(i % 11 == 0 ? aln_test::lcg_dna(state, len2) : m);
    }
    *f1 = write_fastq(tmp.file("rep1.fq"), r1);
    *f2 = write_fastq(tmp.file("rep2.fq"), r2);
}

// Reads the two files into *b; with `searched` (the same files, already loaded and searched with
// the same options) its search results are copied in, else the reads are searched. The search is
// deterministic, and repeating it would dominate the tests' run time.
void load_rep_batch(const ba::gap_opt_t &gap, const std::string &f1, const std::string &f2, Batch *b,
                    const Batch *searched = nullptr)
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
    if (searched) {
        REQUIRE(searched->n == b->n);
        for (int j = 0; j < 2; ++j)
            for (int i = 0; i < b->n; ++i) { // what aln_search_batch sets on each read
                ba::bwa_seq_t *q = b->seqs[j] + i;
                const ba::bwa_seq_t *s = searched->seqs[j] + i;
                q->sa = 0, q->type = BWA_TYPE_NO_MATCH, q->c1 = q->c2 = 0;
                q->n_aln = s->n_aln;
                q->aln = nullptr;
                if (s->n_aln) {
                    q->aln = (ba::bwt_aln1_t *)malloc((size_t)s->n_aln * sizeof(ba::bwt_aln1_t));
                    if (q->aln == nullptr) FAIL("out of memory copying the search results");
                    memcpy(q->aln, s->aln, (size_t)s->n_aln * sizeof(ba::bwt_aln1_t));
                }
            }
        return;
    }
    RepIndex &ix = rep_fixture();
    ba::aln_search_batch(ix.bwt, b->n, b->seqs[0], &gap, 1);
    ba::aln_search_batch(ix.bwt, b->n, b->seqs[1], &gap, 1);
}

// The rep_pairs reads searched with the default gap options, once for every case that needs them.
const Batch &rep_searched()
{
    static Batch searched;
    static const bool loaded = [] {
        bwa_tests::TempDir tmp("aln_pe_rep_unit");
        std::string f1, f2;
        rep_pairs(tmp, &f1, &f2);
        ba::gap_opt_t *gap = ba::gap_init_opt();
        load_rep_batch(*gap, f1, f2, &searched);
        free(gap);
        return true;
    }();
    (void)loaded;
    return searched;
}

// Every field the position, pairing and extra-hit stages write, for one read.
typedef std::vector<long long> ReadState;
ReadState read_state(const ba::bwa_seq_t *p)
{
    ReadState v = {p->type, (long long)p->pos, p->strand, p->mapQ, (long long)p->seQ, p->extra_flag, p->n_mm,
                   p->n_gapo, p->n_gape, p->score, (long long)p->c1, (long long)p->c2, (long long)p->sa, p->n_multi};
    for (int k = 0; k < p->n_multi; ++k) {
        const ba::bwt_multi1_t *q = p->multi + k;
        v.insert(v.end(), {(long long)q->pos, q->strand, q->gap, q->mm, q->ref_shift});
    }
    return v;
}

// read_state of every read, read 1 and read 2 of each pair in turn.
std::vector<ReadState> batch_state(const Batch &b)
{
    std::vector<ReadState> v;
    v.reserve((size_t)b.n * 2);
    for (int i = 0; i < b.n; ++i)
        for (int j = 0; j < 2; ++j) v.push_back(read_state(b.seqs[j] + i));
    return v;
}

// The index of the first read whose state differs between a and b, or -1 if none does (one
// assertion per comparison rather than one per read keeps the cases fast).
long long first_diff(const std::vector<ReadState> &a, const std::vector<ReadState> &b)
{
    for (size_t r = 0; r < a.size() || r < b.size(); ++r)
        if (r >= a.size() || r >= b.size() || a[r] != b[r]) return (long long)r;
    return -1;
}

// Reference: upstream's bwa_sa2pos body as one function (bwt_sa, then its position tail), written
// out here so the references do not share bwa_pos_from_sa with the code under test.
ba::bwtint_t ref_sa2pos(const bntseq_t *bns, const ba::bwt_t *bwt, ba::bwtint_t sapos, int ref_len, int *strand)
{
    ba::bwtint_t pos_f = ba::bwt_sa(bwt, sapos); // position on the forward-reverse coordinate
    int is_rev;
    *strand = 0;
    if (pos_f < (ba::bwtint_t)bns->l_pac && (ba::bwtint_t)bns->l_pac < pos_f + ref_len) return (ba::bwtint_t)-1;
    pos_f = bns_depos(bns, pos_f, &is_rev);
    *strand = !is_rev;
    if (is_rev) pos_f = pos_f + 1 < (ba::bwtint_t)ref_len ? 0 : pos_f - ref_len + 1;
    return pos_f;
}

// Reference: bwa_pe_se_pos as upstream wrote it, one bwa_sa2pos per read.
void ref_pe_se_pos(const bntseq_t *bns, const ba::bwt_t *bwt, int n_seqs, ba::bwa_seq_t *seqs[2], const ba::gap_opt_t *gopt)
{
    for (int i = 0; i != n_seqs; ++i)
        for (int j = 0; j < 2; ++j) {
            ba::bwa_seq_t *p = seqs[j] + i;
            if (p->type == BWA_TYPE_UNIQUE || p->type == BWA_TYPE_REPEAT) {
                int strand;
                int max_diff = gopt->fnr > 0.0 ? ba::bwa_cal_maxdiff(p->len, BWA_AVG_ERR, gopt->fnr) : gopt->max_diff;
                p->seQ = p->mapQ = ba::bwa_approx_mapQ(p, max_diff);
                p->pos = ref_sa2pos(bns, bwt, p->sa, p->len + p->ref_shift, &strand);
                p->strand = strand;
                if (p->pos == (ba::bwtint_t)-1) p->type = BWA_TYPE_NO_MATCH;
            }
        }
}

// Reference: bwa_cal_pac_pos as upstream wrote it.
void ref_cal_pac_pos(const bntseq_t *bns, const ba::bwt_t *bwt, int n_seqs, ba::bwa_seq_t *seqs, int max_mm, float fnr)
{
    for (int i = 0; i != n_seqs; ++i) {
        ba::bwa_seq_t *p = &seqs[i];
        if (p->type == BWA_TYPE_UNIQUE || p->type == BWA_TYPE_REPEAT) { // upstream's bwa_cal_pac_pos_core
            int strand;
            const int max_diff = fnr > 0.0 ? ba::bwa_cal_maxdiff(p->len, BWA_AVG_ERR, fnr) : max_mm;
            p->seQ = p->mapQ = ba::bwa_approx_mapQ(p, max_diff);
            p->pos = ref_sa2pos(bns, bwt, p->sa, p->len + p->ref_shift, &strand);
            p->strand = strand;
            p->seQ = p->mapQ = ba::bwa_approx_mapQ(p, max_diff);
            if (p->pos == (ba::bwtint_t)-1) p->type = BWA_TYPE_NO_MATCH;
        }
        int n_multi = 0;
        for (int j = 0; j < p->n_multi; ++j) {
            int strand;
            ba::bwt_multi1_t *q = p->multi + j;
            q->pos = ref_sa2pos(bns, bwt, q->pos, p->len + q->ref_shift, &strand);
            q->strand = strand;
            if (q->pos != p->pos && q->pos != (ba::bwtint_t)-1) p->multi[n_multi++] = *q;
        }
        p->n_multi = n_multi;
    }
}

// Reference: bwa_pe_pair as upstream wrote it (its loop over the pairs, one bwa_sa2pos per row).
int ref_pe_pair(const bntseq_t *bns, const ba::bwt_t *bwt, int n_seqs, ba::bwa_seq_t *seqs[2],
                const ba::isize_info_t *ii, const ba::pe_opt_t *opt, const ba::gap_opt_t *gopt)
{
    int cnt_chg = 0;
    pair64_v arr = {0, 0, 0};
    for (int i = 0; i != n_seqs; ++i) {
        ba::bwa_seq_t *p[2] = {seqs[0] + i, seqs[1] + i};
        if ((p[0]->type == BWA_TYPE_UNIQUE || p[0]->type == BWA_TYPE_REPEAT)
            && (p[1]->type == BWA_TYPE_UNIQUE || p[1]->type == BWA_TYPE_REPEAT)) {
            long long n_occ[2];
            for (int j = 0; j < 2; ++j) {
                n_occ[j] = 0;
                for (int k = 0; k < p[j]->n_aln; ++k) n_occ[j] += p[j]->aln[k].l - p[j]->aln[k].k + 1;
            }
            if (n_occ[0] > opt->max_occ || n_occ[1] > opt->max_occ) continue;
            arr.n = 0;
            for (int j = 0; j < 2; ++j)
                for (int k = 0; k < p[j]->n_aln; ++k) {
                    const ba::bwt_aln1_t *r = p[j]->aln + k;
                    for (ba::bwtint_t l = r->k; l <= r->l; ++l) {
                        int strand;
                        pair64_t x;
                        x.x = ref_sa2pos(bns, bwt, l, p[j]->len + p[j]->ref_shift, &strand);
                        x.y = k << 2 | strand << 1 | j;
                        kv_push(pair64_t, arr, x);
                    }
                }
            // bwa_pe_pairing sorts with pdqsort_128, not upstream's ks_introsort_128, so this
            // reference shares that sort; test_pair64_sort.cpp checks the two sorts agree.
            cnt_chg += ba::bwa_pe_pairing(p, &arr, opt, gopt->s_mm, ii);
        }
        if (opt->N_multi || opt->n_multi) {
            for (int j = 0; j < 2; ++j) {
                if (p[j]->type != BWA_TYPE_NO_MATCH) {
                    if (!(p[j]->extra_flag & SAM_FPP) && p[1 - j]->type != BWA_TYPE_NO_MATCH)
                        ba::bwa_aln2seq_core(p[j]->n_aln, p[j]->aln, p[j], 0, p[j]->c1 + p[j]->c2 - 1 > opt->N_multi ? opt->n_multi : opt->N_multi);
                    else ba::bwa_aln2seq_core(p[j]->n_aln, p[j]->aln, p[j], 0, opt->n_multi);
                    int n_multi = 0;
                    for (int k = 0; k < p[j]->n_multi; ++k) {
                        int strand;
                        ba::bwt_multi1_t *q = p[j]->multi + k;
                        q->pos = ref_sa2pos(bns, bwt, q->pos, p[j]->len + q->ref_shift, &strand);
                        q->strand = strand;
                        if (q->pos != p[j]->pos && q->pos != (ba::bwtint_t)-1) p[j]->multi[n_multi++] = *q;
                    }
                    p[j]->n_multi = n_multi;
                }
            }
        }
    }
    kv_destroy(arr);
    return cnt_chg;
}

} // namespace

TEST_CASE("aln_pe: batched single-end positions equal one bwa_sa2pos per read"
          * doctest::test_suite("unit/aln_pe")) {
    RepIndex &ix = rep_fixture();
    ba::bwase_initialize();
    ba::gap_opt_t *gap = ba::gap_init_opt();
    bwa_tests::TempDir tmp("aln_pe_rep_unit");
    std::string f1, f2;
    rep_pairs(tmp, &f1, &f2);
    Batch a, b;
    load_rep_batch(*gap, f1, f2, &a, &rep_searched());
    load_rep_batch(*gap, f1, f2, &b, &rep_searched());
    srand48(ix.bns()->seed);
    ba::bwa_pe_aln2seq(a.n, a.seqs);
    srand48(ix.bns()->seed);
    ba::bwa_pe_aln2seq(b.n, b.seqs);
    REQUIRE(first_diff(batch_state(a), batch_state(b)) == -1);
    ref_pe_se_pos(ix.bns(), ix.bwt, a.n, a.seqs, gap);
    ba::bwa_pe_se_pos(ix.bns(), ix.bwt, b.n, b.seqs, gap);
    int n_mapped = 0;
    for (int i = 0; i < a.n; ++i) n_mapped += a.seqs[0][i].type != BWA_TYPE_NO_MATCH;
    CHECK(n_mapped > a.n / 2);
    CHECK(first_diff(batch_state(a), batch_state(b)) == -1);
    free(gap);
}

TEST_CASE("aln_se: batched bwa_cal_pac_pos equals upstream's loop"
          * doctest::test_suite("unit/aln_pe")) {
    RepIndex &ix = rep_fixture();
    ba::bwase_initialize();
    ba::gap_opt_t *gap = ba::gap_init_opt();
    bwa_tests::TempDir tmp("aln_pe_rep_unit");
    std::string f1, f2;
    rep_pairs(tmp, &f1, &f2);
    const Batch &searched = rep_searched();
    for (int max_xa : {3, 0, 50}) {
        Batch a, b; // fresh reads for each max_xa, so no pass starts from another's extra hits
        load_rep_batch(*gap, f1, f2, &a, &searched);
        load_rep_batch(*gap, f1, f2, &b, &searched);
        for (int j = 0; j < 2; ++j) { // each file as a single-end batch
            srand48(ix.bns()->seed);
            for (int i = 0; i < a.n; ++i) ba::bwa_aln2seq_core(a.seqs[j][i].n_aln, a.seqs[j][i].aln, &a.seqs[j][i], 1, max_xa);
            srand48(ix.bns()->seed);
            for (int i = 0; i < b.n; ++i) ba::bwa_aln2seq_core(b.seqs[j][i].n_aln, b.seqs[j][i].aln, &b.seqs[j][i], 1, max_xa);
            int n_multi = 0;
            for (int i = 0; i < a.n; ++i) n_multi += a.seqs[j][i].n_multi;
            if (max_xa == 50) CHECK(n_multi > a.n); // repeats give extra hits
            if (max_xa == 0) CHECK(n_multi == 0);
            ref_cal_pac_pos(ix.bns(), ix.bwt, a.n, a.seqs[j], gap->max_diff, gap->fnr);
            ba::bwa_cal_pac_pos(ix.bns(), ix.bwt, b.n, b.seqs[j], gap->max_diff, gap->fnr);
            int n_shift = 0; // reads and extra hits whose reference length differs from the read's
            for (int i = 0; i < a.n; ++i) {
                n_shift += a.seqs[j][i].ref_shift != 0;
                for (int k = 0; k < a.seqs[j][i].n_multi; ++k) n_shift += a.seqs[j][i].multi[k].ref_shift != 0;
            }
            CHECK(n_shift > 0);
            long long diff = -1;
            for (int i = 0; i < a.n && diff < 0; ++i)
                if (read_state(a.seqs[j] + i) != read_state(b.seqs[j] + i)) diff = i;
            REQUIRE_MESSAGE(diff == -1, "max_xa " << max_xa << " file " << j << " read " << diff);
        }
    }
    free(gap);
}

namespace {

struct GroupOpts { int max_occ, n_multi, N_multi; };

// bwa_pe_pair, and bwa_pe_pair_grouped at group weights from one pair per group to the whole
// batch, against upstream's loop (ref_pe_pair) under the pairing options `op`. One case per
// option set keeps each case fast.
void check_grouping(const GroupOpts &op)
{
    RepIndex &ix = rep_fixture();
    ba::bwase_initialize();
    ba::gap_opt_t *gap = ba::gap_init_opt();
    bwa_tests::TempDir tmp("aln_pe_rep_unit");
    std::string f1, f2;
    rep_pairs(tmp, &f1, &f2);
    const Batch &searched = rep_searched(); // every run below starts from a fresh copy of its results
    ba::pe_opt_t *pe = ba::bwa_init_pe_opt();
    pe->max_occ = op.max_occ, pe->n_multi = op.n_multi, pe->N_multi = op.N_multi;
    std::vector<ReadState> want;
    int want_chg = 0;
    // rows 0: upstream's loop (ref_pe_pair), the reference; -1: bwa_pe_pair itself (group weight
    // ALN_PE_GROUP_ROWS); then group weights from one (every pair its own group) to the batch
    for (long long rows : {0LL, -1LL, 1LL, 333LL, 1LL << 20}) {
        Batch b;
        load_rep_batch(*gap, f1, f2, &b, &searched);
        srand48(ix.bns()->seed);
        ba::isize_info_t ii, last;
        last.avg = -1.0;
        ba::bwa_pe_aln2seq(b.n, b.seqs);
        ba::bwa_pe_se_pos(ix.bns(), ix.bwt, b.n, b.seqs, gap);
        ba::bwa_pe_isize(ix.bwt, b.n, b.seqs, &ii, pe, &last);
        // pairing rows (both ends mapped, within max_occ), as the pairing stage sees them
        long long n_pair_rows = 0;
        for (int i = 0; i < b.n; ++i) {
            long long occ[2] = {0, 0};
            bool mapped = true;
            for (int j = 0; j < 2; ++j) {
                const ba::bwa_seq_t *q = b.seqs[j] + i;
                mapped = mapped && (q->type == BWA_TYPE_UNIQUE || q->type == BWA_TYPE_REPEAT);
                for (int k = 0; k < q->n_aln; ++k) occ[j] += q->aln[k].l - q->aln[k].k + 1;
            }
            if (mapped && occ[0] <= op.max_occ && occ[1] <= op.max_occ) n_pair_rows += occ[0] + occ[1];
        }
        // with the default max_occ, the default group size (ALN_PE_GROUP_ROWS) splits the batch
        if (op.max_occ == 100000) CHECK_MESSAGE(n_pair_rows > ba::ALN_PE_GROUP_ROWS, n_pair_rows << " pairing rows");
        const int chg = rows == 0 ? ref_pe_pair(ix.bns(), ix.bwt, b.n, b.seqs, &ii, pe, gap)
                                  : rows == -1 ? ba::bwa_pe_pair(ix.bns(), ix.bwt, b.n, b.seqs, &ii, pe, gap)
                                                 : ba::bwa_pe_pair_grouped(ix.bns(), ix.bwt, b.n, b.seqs, &ii, pe, gap, rows);
        std::vector<ReadState> got = batch_state(b);
        if (rows == 0) {
            want = got, want_chg = chg;
            // the input exercises what the grouping must preserve
            int n_big = 0, n_fpp = 0, n_multi = 0;
            for (int i = 0; i < b.n; ++i) {
                long long occ = 0;
                for (int j = 0; j < 2; ++j)
                    for (int k = 0; k < b.seqs[j][i].n_aln; ++k) occ += b.seqs[j][i].aln[k].l - b.seqs[j][i].aln[k].k + 1;
                n_big += occ > 50;
                n_fpp += (b.seqs[0][i].extra_flag & SAM_FPP) != 0;
                n_multi += b.seqs[0][i].n_multi + b.seqs[1][i].n_multi;
            }
            CHECK(ii.avg > 0.0);
            CHECK(n_big > 0);
            CHECK(n_fpp > b.n / 4);
            if (op.n_multi + op.N_multi) CHECK(n_multi > 0);
            continue;
        }
        CHECK_MESSAGE(chg == want_chg, "group rows " << rows << " max_occ " << op.max_occ);
        const long long diff = first_diff(got, want);
        REQUIRE_MESSAGE(diff == -1, "group rows " << rows << " max_occ " << op.max_occ << " read " << diff);
    }
    free(pe);
    free(gap);
}

} // namespace

TEST_CASE("aln_pe: batched pairing equals upstream's loop for every grouping (default options)"
          * doctest::test_suite("unit/aln_pe")) {
    check_grouping({100000, 3, 10});
}

// a max_occ that sends the repeat pairs down upstream's `continue` (which also skips their
// extra hits), with few and with many extra hits allowed
TEST_CASE("aln_pe: batched pairing equals upstream's loop for every grouping (low max_occ)"
          * doctest::test_suite("unit/aln_pe")) {
    check_grouping({20, 3, 10});
}

TEST_CASE("aln_pe: batched pairing equals upstream's loop for every grouping (low max_occ, many extra hits)"
          * doctest::test_suite("unit/aln_pe")) {
    check_grouping({20, 30, 200});
}

TEST_CASE("aln_pe: batched pairing equals upstream's loop for every grouping (no extra hits)"
          * doctest::test_suite("unit/aln_pe")) {
    check_grouping({100000, 0, 0});
}

TEST_CASE("aln_pe: batched pairing equals upstream's loop for every grouping (many extra hits)"
          * doctest::test_suite("unit/aln_pe")) {
    check_grouping({100000, 30, 200});
}

namespace {

// Runs aln_run_pe (the shipped driver) with batches of batch_size pairs. Returns its exit code and
// the SAM record lines (no header) in output order; *err (if given) gets what it printed on stderr.
int run_pe_lines(int batch_size, const std::string &f1, const std::string &f2, const std::string &out_path,
                 std::vector<std::string> *lines, std::string *err = nullptr, int share_dup_search = 1)
{
    PeIndex &ix = index_fixture();
    ba::aln_opt_t o;
    ba::aln_opt_init(&o);
    o.batch_size = batch_size;
    o.share_dup_search = share_dup_search;
    ba::aln_out_t *out = ba::aln_out_open(out_path.c_str(), 0, 0, "");
    REQUIRE(out != nullptr);
    const std::string err_path = out_path + ".err";
    fflush(stderr);
    const int saved = dup(fileno(stderr));
    FILE *ef = fopen(err_path.c_str(), "w");
    REQUIRE(saved >= 0);
    REQUIRE(ef != nullptr);
    REQUIRE(dup2(fileno(ef), fileno(stderr)) >= 0);
    const int ret = ba::aln_run_pe(&o, ix.bwt, ix.bns(), ix.fmi->idx->pac, f1.c_str(), f2.c_str(), out);
    fflush(stderr);
    REQUIRE(dup2(saved, fileno(stderr)) >= 0);
    close(saved);
    fclose(ef);
    REQUIRE(ba::aln_out_close(out) == 0);
    if (err) {
        std::ifstream e(err_path);
        std::stringstream ss;
        ss << e.rdbuf();
        *err = ss.str();
    }
    std::ifstream in(out_path);
    lines->clear();
    for (std::string l; std::getline(in, l);)
        if (!l.empty() && l[0] != '@') lines->push_back(l);
    return ret;
}

// Column `col` (0-based) of every SAM line.
std::vector<std::string> sam_column(const std::vector<std::string> &lines, int col)
{
    std::vector<std::string> v;
    for (const std::string &l : lines) {
        size_t st = 0;
        for (int c = 0; c < col; ++c) st = l.find('\t', st) + 1;
        v.push_back(l.substr(st, l.find('\t', st) - st));
    }
    return v;
}

// Runs aln_run_pe and returns the FLAG of every record, in output order.
std::vector<int> run_pe_flags(int batch_size, const std::string &f1, const std::string &f2, const std::string &out_path)
{
    std::vector<std::string> lines;
    REQUIRE(run_pe_lines(batch_size, f1, f2, out_path, &lines) == 0);
    std::vector<int> flags;
    for (const std::string &f : sam_column(lines, 1)) flags.push_back(std::atoi(f.c_str()));
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

TEST_CASE("aln_pe: duplicated pairs give the same records whether or not they share a search"
          * doctest::test_suite("unit/aln_pe")) {
    // 200 pairs: R1 is mostly copies of a few reads, R2 has no two alike (as when one mate is
    // more error-prone), in one batch and in batches of 7 (each file keeps its own duplicate check).
    PeIndex &ix = index_fixture();
    std::vector<std::string> r1, r2;
    well_placed(ix.contig, 200, 36, 100, &r1, &r2);
    for (size_t i = 0; i < r1.size(); ++i)
        if (i % 4 != 0) r1[i] = r1[(i * 7) % 12];
    bwa_tests::TempDir tmp{"aln_pe_dupl"}; // this case's scratch files, removed when it ends
    const std::string f1 = write_fastq(tmp.file("dupl1.fq"), r1), f2 = write_fastq(tmp.file("dupl2.fq"), r2);
    struct RestoreVerbose { int v = bwa_verbose; ~RestoreVerbose() { bwa_verbose = v; } } restore;
    bwa_verbose = 4; // the driver then reports how many reads of each file shared a search
    auto shared = [](const std::string &err) { // R1's reads that shared a search, over all batches
        long n = 0;
        const std::string tag = "[M::aln] read 1: ";
        for (size_t at = err.find(tag); at != std::string::npos; at = err.find(tag, at + 1))
            n += std::atol(err.c_str() + at + tag.size());
        return n;
    };
    for (int bs : {0x40000, 7}) {
        INFO("batch size " << bs);
        std::vector<std::string> on, off;
        std::string err_on, err_off;
        REQUIRE(run_pe_lines(bs, f1, f2, tmp.file("dupl_on.sam"), &on, &err_on, 1) == 0);
        REQUIRE(run_pe_lines(bs, f1, f2, tmp.file("dupl_off.sam"), &off, &err_off, 0) == 0);
        REQUIRE(on.size() == 400);
        CHECK(on == off);
        // not vacuous: R1's copies shared a search with the switch on, and nothing did with it off
        if (bs > 400) CHECK(shared(err_on) > 100); // batches of 7 hold few copies of each other
        CHECK(shared(err_off) == 0);
        CHECK(err_off.find("shared an identical read's search") != std::string::npos);
    }
    // R1 a single read copied 200 times, in batches of 7: every R1 batch has duplicates while no R2
    // batch has any, so a check R2 turned off for both files would leave R1's copies searched. With
    // a check per file, every R1 batch shares: 6 reads in each of 28 batches, 3 in the last.
    const std::string f1_same = write_fastq(tmp.file("dupl1_same.fq"), std::vector<std::string>(r1.size(), r1[0]));
    std::vector<std::string> on, off;
    std::string err_on;
    REQUIRE(run_pe_lines(7, f1_same, f2, tmp.file("same_on.sam"), &on, &err_on, 1) == 0);
    REQUIRE(run_pe_lines(7, f1_same, f2, tmp.file("same_off.sam"), &off, nullptr, 0) == 0);
    REQUIRE(on.size() == 400);
    CHECK(on == off);
    CHECK(shared(err_on) == 28 * 6 + 3);
}

namespace {

// The record names of pairs p0 .. p(n - 1), two records each.
std::vector<std::string> pair_names(int n)
{
    std::vector<std::string> v;
    for (int i = 0; i < n; ++i) v.insert(v.end(), 2, "p" + std::to_string(i));
    return v;
}

} // namespace

// The driver reads the next batch (or pair of batches) while the current one is aligned. A run that
// stops early (pair names that differ, one file ending first) must still stop at the same record
// as when each batch was read only after the previous one was written, with batches of 3 pairs so
// that every stop falls after a batch has already been read ahead.
TEST_CASE("aln_pe: aln_run_pe stops at the same pair with the next batch read ahead"
          * doctest::test_suite("unit/aln_pe")) {
    PeIndex &ix = index_fixture();
    bwa_tests::TempDir scratch{"aln_pe_case"}; // this case's files, removed when it ends
    std::vector<std::string> r1, r2;
    well_placed(ix.contig, 10, 36, 100, &r1, &r2);
    std::vector<std::string> lines;
    std::string err;
    SUBCASE("pair names differ in the second batch: output ends with that pair, exit 1") {
        CHECK(run_pe_lines(3, write_fastq(scratch.path() + "/ra_n1.fq", r1), write_fastq(scratch.path() + "/ra_n2.fq", r2, 4),
                           scratch.path() + "/ra_n.sam", &lines, &err) == 1);
        std::vector<std::string> want = pair_names(4);
        want.push_back("p4"), want.push_back("x4");
        CHECK(sam_column(lines, 0) == want);
        CHECK(err.find("paired reads have different names: \"p4\", \"x4\"") != std::string::npos);
    }
    SUBCASE("read 2 file ends first: the pairs up to its last read, exit 0") {
        std::vector<std::string> s2(r2.begin(), r2.begin() + 7);
        CHECK(run_pe_lines(3, write_fastq(scratch.path() + "/ra_a1.fq", r1), write_fastq(scratch.path() + "/ra_a2.fq", s2),
                           scratch.path() + "/ra_a.sam", &lines, &err) == 0);
        CHECK(sam_column(lines, 0) == pair_names(7));
        CHECK(err.find("read counts differ in this batch: R1 3, R2 1") != std::string::npos);
        CHECK(err.find("read 2 file ended before read 1 file") != std::string::npos);
        CHECK(err.find("read 1 file ended before read 2 file") == std::string::npos);
    }
    SUBCASE("read 1 file ends first: the pairs up to its last read, exit 0, and R2's rest is noticed") {
        std::vector<std::string> s1(r1.begin(), r1.begin() + 7);
        CHECK(run_pe_lines(3, write_fastq(scratch.path() + "/ra_b1.fq", s1), write_fastq(scratch.path() + "/ra_b2.fq", r2),
                           scratch.path() + "/ra_b.sam", &lines, &err) == 0);
        CHECK(sam_column(lines, 0) == pair_names(7));
        // R2 is read only after a non-empty R1 batch, so p9 is still unread when R1 ends: the
        // driver reads it and warns (a read-ahead that read R2 past R1's end would not)
        CHECK(err.find("read 1 file ended before read 2 file") != std::string::npos);
    }
    SUBCASE("equal files in batches of 3 and in one batch: the same records") {
        const std::string f1 = write_fastq(scratch.path() + "/ra_e1.fq", r1), f2 = write_fastq(scratch.path() + "/ra_e2.fq", r2);
        std::vector<std::string> one;
        CHECK(run_pe_lines(3, f1, f2, scratch.path() + "/ra_e3.sam", &lines) == 0);
        CHECK(run_pe_lines(0x40000, f1, f2, scratch.path() + "/ra_e.sam", &one) == 0);
        CHECK(sam_column(lines, 0) == pair_names(10));
        // the fields that do not depend on the per-batch insert-size estimate: QNAME, RNAME, POS,
        // CIGAR, SEQ, QUAL (FLAG, MAPQ and the mate fields can follow the estimate)
        for (int col : {0, 2, 3, 5, 9, 10}) {
            CAPTURE(col);
            CHECK(sam_column(lines, col) == sam_column(one, col));
        }
    }
}

// A gzip read error in R2's next batch (read ahead while the current pairs are aligned) ends the
// run after the current batch is written, with stock's message and exit status 1; the pairs
// written are the first pairs, those of the batches read before the failing read (counted from
// zlib directly, not with the reader under test).
TEST_CASE("aln_pe: a gzip read error in the next R2 batch still lets the current batch be written"
          * doctest::test_suite("unit/aln_pe")) {
    PeIndex &ix = index_fixture();
    bwa_tests::TempDir scratch{"aln_pe_case"}; // this case's files, removed when it ends
    std::vector<std::string> r1, r2;
    well_placed(ix.contig, 200, 36, 100, &r1, &r2);
    std::vector<std::string> a1, a2;
    for (int k = 0; k < 30; ++k) a1.insert(a1.end(), r1.begin(), r1.end()), a2.insert(a2.end(), r2.begin(), r2.end());
    const std::string f1 = write_fastq(scratch.path() + "/gz1.fq", a1);
    const std::string f2 = aln_test::write_fastq_gz_corrupt_tail(scratch.path() + "/gz2.fq.gz", a2, "p");
    const int batch = 1000;
    const int expected_pairs = aln_test::gz_reads_before_error(f2, batch);
    REQUIRE(expected_pairs >= 2 * batch); // the error falls in a batch read ahead while an earlier one is aligned
    std::vector<std::string> lines;
    std::string err;
    CHECK(run_pe_lines(batch, f1, f2, scratch.path() + "/gz.sam", &lines, &err) == 1);
    CHECK((int)lines.size() == 2 * expected_pairs);
    CHECK(sam_column(lines, 0) == pair_names(expected_pairs)); // the first pairs, in order
    CHECK(err.find("[gzread] " + f2 + ": invalid block type\n") != std::string::npos);
}

// The same error in R2's first batch: R1's first batch was read and is freed (LeakSanitizer checks
// it under ASAN), nothing is written, and the run returns 1 with the message.
TEST_CASE("aln_pe: a gzip read error in the first R2 batch writes nothing and returns 1"
          * doctest::test_suite("unit/aln_pe")) {
    PeIndex &ix = index_fixture();
    bwa_tests::TempDir scratch{"aln_pe_case"}; // this case's files, removed when it ends
    std::vector<std::string> r1, r2;
    well_placed(ix.contig, 50, 36, 100, &r1, &r2);
    const std::string f1 = write_fastq(scratch.path() + "/gzfirst1.fq", r1);
    const std::string f2 = aln_test::write_fastq_gz_corrupt_tail(scratch.path() + "/gzfirst2.fq.gz", r2, "p");
    const int batch = 1000; // more than the file holds: the first read of R2 fails
    REQUIRE(aln_test::gz_reads_before_error(f2, batch) == 0);
    std::vector<std::string> lines;
    std::string err;
    CHECK(run_pe_lines(batch, f1, f2, scratch.path() + "/gzfirst.sam", &lines, &err) == 1);
    CHECK(lines.empty());
    CHECK(err.find("[gzread] " + f2 + ": invalid block type\n") != std::string::npos);
}

TEST_CASE("aln_pe: short-insert rescue windows may cover the anchor and keep the mate inside the fragment"
          * doctest::test_suite("unit/aln_pe")) {
    auto *pe = ba::bwa_init_pe_opt();
    ba::bwa_seq_t anchor = {};
    anchor.pos = 1000; anchor.len = 30; anchor.ref_shift = 2; // 32 reference bases (a deletion)
    auto cigar = [](std::initializer_list<std::pair<int, int>> ops) {
        std::vector<ba::bwa_cigar_t> c;
        for (auto [op, n] : ops) c.push_back(__cigar_create(op, n));
        return c;
    };
    for (int option : {0, 1}) {
        CAPTURE(option);
        pe->short_insert = option;
        CHECK(ba::aln_rescue_window_floor(&anchor, pe) == (option ? 1000 : 1030));
        CHECK(ba::aln_rescue_window_ceiling(&anchor, pe) == (option ? 1032 : 1000));
    }
    const auto m32 = cigar({{FROM_M, 32}}), clipped = cigar({{FROM_S, 25}, {FROM_M, 32}});
    const auto m20 = cigar({{FROM_M, 20}}), m40 = cigar({{FROM_M, 40}});
    anchor.strand = 0; // forward anchor [1000, 1032): a reverse mate may start at 1000 and must end at >= 1032
    CHECK(ba::aln_rescue_within_fragment(&anchor, 1000, m32.data(), (int)m32.size()));
    CHECK(ba::aln_rescue_within_fragment(&anchor, 1000, clipped.data(), (int)clipped.size()));
    CHECK(ba::aln_rescue_within_fragment(&anchor, 1010, m40.data(), (int)m40.size()));
    CHECK_FALSE(ba::aln_rescue_within_fragment(&anchor, 999, m40.data(), (int)m40.size()));  // starts first
    CHECK_FALSE(ba::aln_rescue_within_fragment(&anchor, 1005, m20.data(), (int)m20.size())); // ends inside
    anchor.strand = 1; // reverse anchor [1000, 1032): a forward mate must start at <= 1000 and end by 1032
    CHECK(ba::aln_rescue_within_fragment(&anchor, 1000, m32.data(), (int)m32.size()));
    CHECK(ba::aln_rescue_within_fragment(&anchor, 1000, clipped.data(), (int)clipped.size()));
    CHECK(ba::aln_rescue_within_fragment(&anchor, 980, m40.data(), (int)m40.size()));
    CHECK_FALSE(ba::aln_rescue_within_fragment(&anchor, 1001, m20.data(), (int)m20.size())); // starts after
    CHECK_FALSE(ba::aln_rescue_within_fragment(&anchor, 995, m40.data(), (int)m40.size()));  // ends after
    free(pe);
}

TEST_CASE("aln_pe: ref-span-insert is an appended option that defaults off"
          * doctest::test_suite("unit/aln_pe")) {
    CHECK(offsetof(ba::pe_opt_t, ref_span_insert) >= offsetof(ba::pe_opt_t, short_insert) + sizeof(int));
    auto *opt = ba::bwa_init_pe_opt();
    REQUIRE(opt != nullptr);
    CHECK(opt->ref_span_insert == 0);
    free(opt);
}

TEST_CASE("aln_pe: ref-span-insert measures each insert-size observation from the alignments' reference spans"
          * doctest::test_suite("unit/aln_pe")) {
    // 30 planted pairs: read 1 forward, 50 bp, ungapped; read 2 reverse, 50 bp, with a 4-base
    // deletion (spans 54 reference bases) or, every third pair, a 3-base insertion (spans 47).
    // The fragment, from read 1's start to read 2's reference end, is frag_i bases. Read-length
    // measurement sees frag_i - ref_shift_i instead.
    ba::bwa_seq_t ends[2][30] = {};
    double sum_span = 0, sum_len = 0;
    for (int i = 0; i < 30; ++i) {
        const int frag = 180 + i % 21, shift = i % 3 == 0 ? -3 : 4;
        auto &a = ends[0][i], &b = ends[1][i];
        a.mapQ = b.mapQ = 37;
        a.len = a.full_len = b.len = b.full_len = 50;
        a.pos = 100 + 300 * i;
        b.ref_shift = shift; b.strand = 1;
        b.pos = a.pos + frag - (50 + shift);
        sum_span += frag; sum_len += frag - shift;
    }
    ba::bwa_seq_t *seqs[2] = {ends[0], ends[1]};
    ba::isize_info_t stock = {}, span = {};
    REQUIRE(ba::infer_isize(30, seqs, &stock, 1e-5, 100000) == 0);
    REQUIRE(ba::infer_isize(30, seqs, &span, 1e-5, 100000, nullptr, 0, false, true) == 0);
    CHECK(stock.avg == doctest::Approx(sum_len / 30));
    CHECK(span.avg == doctest::Approx(sum_span / 30));
    // The option flag on pe_opt_t reaches infer_isize through the shared driver.
    auto *pe = ba::bwa_init_pe_opt();
    ba::isize_info_t none = {}, ii = {};
    none.avg = -1;
    ba::bwa_pe_isize_ref(100000, 30, seqs, &ii, pe, &none, nullptr, 0);
    check_same_ii(ii, stock);
    pe->ref_span_insert = 1;
    ba::bwa_pe_isize_ref(100000, 30, seqs, &ii, pe, &none, nullptr, 0);
    check_same_ii(ii, span);
    free(pe);
}

TEST_CASE("aln_pe: ref-span-insert pairs on the reverse hit's reference span against -a"
          * doctest::test_suite("unit/aln_pe")) {
    // No estimate (as with -A): a pair is proper when its insert is at most -a 200. Read 1 is
    // forward at 1000, 50 bp. Read 2 is reverse, 50 bp, and its hit has a 20-base deletion (70
    // reference bases) or a 20-base insertion (30 reference bases).
    auto *pe = ba::bwa_init_pe_opt();
    pe->max_isize = 200;
    ba::bwase_initialize();
    struct Planted { int n_del, n_ins, pos; bool stock_proper, span_proper; };
    const Planted cases[] = {
        {20, 0, 1150, true, false}, // read length: 1150 + 50 - 1000 = 200; span: 220
        {0, 20, 1165, false, true}, // read length: 215; span: 1165 + 30 - 1000 = 195
        {0, 0, 1150, true, true},   // ungapped: both 200
        {0, 0, 1151, false, false}, // ungapped: both 201
    };
    for (int option : {0, 1}) for (const auto &c : cases) {
        CAPTURE(option); CAPTURE(c.n_del); CAPTURE(c.n_ins); CAPTURE(c.pos);
        pe->ref_span_insert = option;
        ba::isize_info_t ii = {}; ii.avg = ii.std = -1;
        ba::bwa_seq_t ends[2] = {};
        ba::bwt_aln1_t hits[2] = {};
        hits[1].n_del = c.n_del; hits[1].n_ins = c.n_ins; hits[1].n_gapo = c.n_del || c.n_ins;
        for (int j = 0; j < 2; ++j) {
            ends[j].aln = hits + j; ends[j].type = BWA_TYPE_UNIQUE; ends[j].strand = j;
            ends[j].len = ends[j].full_len = 50; ends[j].mapQ = ends[j].seQ = 25;
        }
        ends[0].pos = 1000; ends[1].pos = c.pos; ends[1].ref_shift = c.n_del - c.n_ins;
        ba::bwa_seq_t *pair[2] = {ends, ends + 1};
        pair64_t rows[2] = {{1000, 0}, {(uint64_t)c.pos, 3}};
        pair64_v arr = {}; arr.a = rows; arr.n = arr.m = 2;
        ba::bwa_pe_pairing(pair, &arr, pe, 3, &ii);
        const bool want = option ? c.span_proper : c.stock_proper;
        for (const auto &p : ends) CHECK(((p.extra_flag & SAM_FPP) != 0) == want);
    }
    free(pe);
}

TEST_CASE("aln_pe: ref-span-insert gives a read that pairing moves its new hit's reference span"
          * doctest::test_suite("unit/aln_pe")) {
    // Read 2's best hit is ungapped and far away; its other hit (index 1) has a 4-base deletion
    // and pairs with read 1. Stock moves the read but keeps the first hit's ref_shift (0); with
    // the option the moved read takes its new hit's span, which refinement and output then use.
    auto *pe = ba::bwa_init_pe_opt();
    ba::bwase_initialize();
    for (int option : {0, 1}) {
        CAPTURE(option);
        pe->ref_span_insert = option;
        ba::isize_info_t ii = {}; ii.avg = ii.std = -1;
        ba::bwa_seq_t ends[2] = {};
        ba::bwt_aln1_t hits0[1] = {}, hits1[2] = {};
        hits1[1].n_del = 4; hits1[1].n_gapo = 1; hits1[1].n_gape = 3;
        ends[0].aln = hits0; ends[1].aln = hits1;
        for (int j = 0; j < 2; ++j) {
            ends[j].type = BWA_TYPE_REPEAT; ends[j].strand = j;
            ends[j].len = ends[j].full_len = 50; ends[j].mapQ = ends[j].seQ = 25;
        }
        ends[0].pos = 1000; ends[1].pos = 5000;
        ba::bwa_seq_t *pair[2] = {ends, ends + 1};
        pair64_t rows[3] = {{1000, 0u << 2 | 0 << 1 | 0}, {5000, 0u << 2 | 1 << 1 | 1}, {1150, 1u << 2 | 1 << 1 | 1}};
        pair64_v arr = {}; arr.a = rows; arr.n = arr.m = 3;
        ba::bwa_pe_pairing(pair, &arr, pe, 3, &ii);
        CHECK(ends[1].pos == 1150);
        CHECK((ends[1].extra_flag & SAM_FPP) != 0);
        CHECK((int)ends[1].n_gapo == 1);
        CHECK(ends[1].ref_shift == (option ? 4 : 0));
    }
    free(pe);
}

TEST_CASE("aln_pe: ref-span-insert places the rescue window from the anchor's reference end"
          * doctest::test_suite("unit/aln_pe")) {
    auto *pe = ba::bwa_init_pe_opt();
    ba::bwa_seq_t anchor = {};
    anchor.pos = 1000; anchor.len = 30; anchor.ref_shift = 2; // 32 reference bases (a deletion)
    for (int short_insert : {0, 1}) for (int option : {0, 1}) {
        CAPTURE(short_insert); CAPTURE(option);
        pe->short_insert = short_insert; pe->ref_span_insert = option;
        CHECK(ba::aln_rescue_anchor_end(&anchor, pe) == (option ? 1032 : 1030));
        // A forward anchor's mate starts after its end (stock) or at its start (--short-insert).
        CHECK(ba::aln_rescue_window_floor(&anchor, pe) == (short_insert ? 1000 : option ? 1032 : 1030));
        // --short-insert already ends a reverse anchor's window at its reference end.
        CHECK(ba::aln_rescue_window_ceiling(&anchor, pe) == (short_insert ? 1032 : 1000));
    }
    free(pe);
}

TEST_CASE("aln_pe: exact-rescue-window is an appended option that defaults off"
          * doctest::test_suite("unit/aln_pe")) {
    CHECK(offsetof(ba::pe_opt_t, exact_rescue_window) >= offsetof(ba::pe_opt_t, ref_span_insert) + sizeof(int));
    auto *opt = ba::bwa_init_pe_opt();
    REQUIRE(opt != nullptr);
    CHECK(opt->exact_rescue_window == 0);
    free(opt);
}

namespace {
// sampe's rescue-window arithmetic as bwa_paired_sw wrote it (each assignment truncates).
void sampe_window(int64_t base, bool mate_reverse, int len, const ba::isize_info_t &ii, int64_t *beg, int64_t *end)
{
    if (mate_reverse) *beg = base + ii.avg - 3 * ii.std - len * 1.5;
    else *beg = base - ii.avg - 3 * ii.std - len * 0.5;
    *end = *beg + 6 * ii.std + 2 * len;
}
} // namespace

TEST_CASE("aln_pe: the stock rescue window loses a base when 6 std crosses an integer; the exact one does not"
          * doctest::test_suite("unit/aln_pe")) {
    // A forward 150 bp anchor and its reverse 150 bp mate. std 12.6685 and 12.6646 differ by 0.004,
    // and 6 std crosses 76 between them (at 12.6667). Stock's end is floor(start) + floor(6 std) +
    // 2 len, so every window of this geometry loses its last base at once; the stated end, avg +
    // 3 std + len / 2 past the anchor, moves by 0.012 of a base and stays inside the same base.
    const int64_t base = 1000000;
    ba::isize_info_t wide = {}, narrow = {};
    wide.avg = narrow.avg = 300.3;
    wide.std = 12.6685; narrow.std = 12.6646;
    int64_t sb[2], se[2], eb[2], ee[2];
    ba::aln_rescue_window(base, true, 150, &wide, false, &sb[0], &se[0]);
    ba::aln_rescue_window(base, true, 150, &narrow, false, &sb[1], &se[1]);
    CHECK(sb[0] == base + 37); CHECK(sb[1] == base + 37);
    CHECK(se[0] == base + 413); CHECK(se[1] == base + 412); // the cohort-wide step
    ba::aln_rescue_window(base, true, 150, &wide, true, &eb[0], &ee[0]);
    ba::aln_rescue_window(base, true, 150, &narrow, true, &eb[1], &ee[1]);
    CHECK(eb[0] == base + 37); CHECK(eb[1] == base + 37);
    CHECK(ee[0] == base + 414); CHECK(ee[1] == base + 414); // ceil(300.3 + 3 std + 75), 413.31 and 413.29
    // The same holds for a reverse anchor's forward mate placed from the anchor's end.
    ba::aln_rescue_window(base, false, 150, &wide, false, &sb[0], &se[0]);
    ba::aln_rescue_window(base, false, 150, &narrow, false, &sb[1], &se[1]);
    CHECK(se[0] - sb[0] == 376); CHECK(se[1] - sb[1] == 375);
    ba::aln_rescue_window(base, false, 150, &wide, true, &eb[0], &ee[0]);
    ba::aln_rescue_window(base, false, 150, &narrow, true, &eb[1], &ee[1]);
    CHECK(eb[0] == eb[1]); CHECK(ee[0] == ee[1]);
}

TEST_CASE("aln_pe: the exact rescue window rounds the insert range outward and is monotone in std"
          * doctest::test_suite("unit/aln_pe")) {
    const int64_t base = 1000000;
    // (Not avg 250.5 with these even lengths: a start whose fraction is exactly 1/2 makes stock's two
    // floors cancel, and its end then steps once per third of std as well, apart from rounding.)
    for (double avg : {99.75, 250.2, 301.37}) for (int len : {30, 100, 150}) for (bool rev : {true, false}) {
        CAPTURE(avg); CAPTURE(len); CAPTURE(rev);
        ba::isize_info_t ii = {};
        ii.avg = avg;
        int64_t pb = 0, pe = 0, pse = 0;
        int stock_moves = 0, exact_moves = 0, stock_drops = 0;
        for (int step = 0; step <= 20000; ++step) {
            ii.std = 5.0 + step * 0.0007; // 5 to 19
            int64_t b, e, sb, se, rb, re;
            ba::aln_rescue_window(base, rev, len, &ii, false, &sb, &se);
            sampe_window(base, rev, len, ii, &rb, &re);
            REQUIRE(sb == rb); REQUIRE(se == re); // without the option: sampe's window
            ba::aln_rescue_window(base, rev, len, &ii, true, &b, &e);
            const double lo = rev ? avg - 3 * ii.std - len * 1.5 : -avg - 3 * ii.std - len * 0.5;
            const double hi = lo + 6 * ii.std + 2 * len;
            // the smallest integer window holding [base + lo, base + hi)
            REQUIRE((double)(b - base) <= lo); REQUIRE(lo < (double)(b - base) + 1);
            REQUIRE((double)(e - base) >= hi - 1e-9); REQUIRE((double)(e - base) - 1 < hi + 1e-9);
            // it holds stock's window, and adds at most two bases at the end
            REQUIRE(b <= sb); REQUIRE(sb - b <= 1); REQUIRE(e >= se); REQUIRE(e - se <= 2);
            if (step) {
                REQUIRE(b <= pb); REQUIRE(e >= pe); // monotone: a wider std never narrows it
                exact_moves += e != pe;
                stock_moves += se != pse;
                stock_drops += se < pse;
            }
            pb = b; pe = e; pse = se;
        }
        // The stated end crosses a base each third of a unit of std (42 times here). Stock's end
        // steps up with floor(6 std) (84 times) and back with floor(start) (42 times).
        CHECK(exact_moves == 42);
        CHECK(stock_moves == 3 * exact_moves);
        CHECK(stock_drops == exact_moves);
    }
}
