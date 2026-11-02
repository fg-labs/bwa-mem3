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
#include "../../src/kvec.h"
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
#include <zlib.h>

namespace {

namespace ba = bwa_aln;

// One 6,000 bp random contig: every read is unique, so well-placed pairs get mapQ >= 20.
struct PeIndex {
    std::string dir, prefix, contig;
    std::unique_ptr<FMI_search> fmi;
    ba::bwt_t *bwt = nullptr;

    PeIndex()
    {
        std::string tmpl_s = aln_test::tmp_template("aln_pe_unit");
        char *tmpl = &tmpl_s[0];
        REQUIRE(mkdtemp(tmpl) != nullptr);
        dir = tmpl;
        prefix = dir + "/ref.fa";
        uint32_t state = 77;
        contig = aln_test::lcg_dna(state, 6000);
        FILE *fp = fopen(prefix.c_str(), "w");
        REQUIRE(fp != nullptr);
        fprintf(fp, ">u1\n");
        for (size_t i = 0; i < contig.size(); i += 60) fprintf(fp, "%s\n", contig.substr(i, 60).c_str());
        REQUIRE(fclose(fp) == 0);
        REQUIRE(bwa_idx_build(prefix.c_str(), prefix.c_str(), 0, 3) == 0);
        fmi.reset(new FMI_search(prefix.c_str()));
        fmi->load_index(true, 1);
        bwt = ba::bwt_from_fmi(fmi.get());
    }
    ~PeIndex()
    {
        ba::bwt_destroy(bwt);
        fmi.reset();
        std::string cmd = "rm -rf '" + dir + "'";
        if (system(cmd.c_str()) != 0) std::fprintf(stderr, "warning: could not remove %s\n", dir.c_str());
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
    for (size_t i = 0; i < reads.size(); ++i)
        fprintf(fp, "@%s%zu\n%s\n+\n%s\n", (int)i == renamed ? "x" : "p", i, reads[i].c_str(),
                std::string(reads[i].size(), 'I').c_str());
    fclose(fp);
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
    load_batch(*gap, write_fastq(ix.dir + "/a1.fq", a1), write_fastq(ix.dir + "/a2.fq", a2), &bt1);
    ba::isize_info_t ii1;
    ba::bwa_cal_pac_pos_pe(ix.bns(), ix.bwt, bt1.n, bt1.seqs, &ii1, pe, gap, &last);
    REQUIRE(ii1.avg > 0.0);
    CHECK(ii1.avg > 280.0);
    CHECK(ii1.avg < 320.0);
    CHECK(ii1.high_bayesian > ii1.low);

    // batch 2: its own estimate fails, so the pairing uses batch 1's
    Batch bt2;
    load_batch(*gap, write_fastq(ix.dir + "/b1.fq", b1), write_fastq(ix.dir + "/b2.fq", b2), &bt2);
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
    load_batch(*gap, write_fastq(ix.dir + "/c1.fq", c1), write_fastq(ix.dir + "/c2.fq", c2), &bt3);
    ba::isize_info_t ii3;
    ba::bwa_cal_pac_pos_pe(ix.bns(), ix.bwt, bt3.n, bt3.seqs, &ii3, pe, gap, &ii2);
    CHECK(ii3.avg > 280.0);
    CHECK(ii3.avg < 320.0);

    // without a previous estimate batch 2 stays at "no estimate"
    Batch bt2b;
    load_batch(*gap, ix.dir + "/b1.fq", ix.dir + "/b2.fq", &bt2b);
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
    ba::gap_opt_t *gap = ba::gap_init_opt();
    gap->max_gapo = 5; // the cap is min(max_gapo, differences allowed for the longest read)
    std::vector<std::string> r1, r2, long1, long2;
    well_placed(ix.contig, 30, 36, 100, &r1, &r2);          // R1: 36 bp reads
    well_placed(ix.contig, 30, 150, 100, &long1, &long2);   // R2: 150 bp reads
    const std::string f1 = write_fastq(ix.dir + "/m1.fq", r1);
    const std::string f2 = write_fastq(ix.dir + "/m2.fq", long2);
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
    const std::string out_path = ix.dir + "/pe.sam";
    ba::aln_out_t *out = ba::aln_out_open(out_path.c_str(), 0, 0, "");
    REQUIRE(out != nullptr);
    REQUIRE(ba::aln_run_pe(&o, ix.bwt, ix.bns(), ix.fmi->idx->pac, f1.c_str(), f2.c_str(), out) == 0);
    REQUIRE(ba::aln_out_close(out) == 0);
    std::ifstream in(out_path);
    int n_rec = 0;
    for (std::string l; std::getline(in, l);) n_rec += !l.empty() && l[0] != '@';
    CHECK(n_rec == 60);
    free(gap);
}

// ---- SA lookups made in batches (bwt_sa_batch in bwa_pe_se_pos, bwa_pe_pair, bwa_cal_pac_pos) ----

namespace {

// A contig built to give reads many hits: random sequence with a 300 bp segment copied 20
// times forward and 10 times reverse-complemented, and a 160 bp "AC" tandem repeat.
struct RepIndex {
    std::string dir, prefix, contig;
    std::unique_ptr<FMI_search> fmi;
    ba::bwt_t *bwt = nullptr;

    RepIndex()
    {
        std::string tmpl_s = aln_test::tmp_template("aln_pe_rep_unit");
        char *tmpl = &tmpl_s[0];
        REQUIRE(mkdtemp(tmpl) != nullptr);
        dir = tmpl;
        prefix = dir + "/ref.fa";
        uint32_t state = 91;
        const std::string seg = aln_test::lcg_dna(state, 300);
        contig = aln_test::lcg_dna(state, 1500);
        for (int i = 0; i < 20; ++i) contig += seg + aln_test::lcg_dna(state, 400 + 37 * (i % 6));
        for (int i = 0; i < 10; ++i) contig += aln_test::revcomp(seg) + aln_test::lcg_dna(state, 500);
        for (int i = 0; i < 80; ++i) contig += "AC";
        contig += aln_test::lcg_dna(state, 1500);
        FILE *fp = fopen(prefix.c_str(), "w");
        REQUIRE(fp != nullptr);
        fprintf(fp, ">r1\n");
        for (size_t i = 0; i < contig.size(); i += 60) fprintf(fp, "%s\n", contig.substr(i, 60).c_str());
        REQUIRE(fclose(fp) == 0);
        REQUIRE(bwa_idx_build(prefix.c_str(), prefix.c_str(), 0, 3) == 0);
        fmi.reset(new FMI_search(prefix.c_str()));
        fmi->load_index(true, 1);
        bwt = ba::bwt_from_fmi(fmi.get());
    }
    ~RepIndex()
    {
        ba::bwt_destroy(bwt);
        fmi.reset();
        std::string cmd = "rm -rf '" + dir + "'";
        if (system(cmd.c_str()) != 0) std::fprintf(stderr, "warning: could not remove %s\n", dir.c_str());
    }
    const bntseq_t *bns() const { return fmi->idx->bns; }
};

RepIndex &rep_fixture()
{
    static RepIndex idx;
    return idx;
}

// Pairs of 36 bp (read 1) and 40 bp (read 2) from 260-340 bp fragments starting every 13 bp along
// the whole contig (so many land in the repeats). Every fifth read 2 has a mismatch, every third
// read 1 a 1 bp deletion and every fourth read 2 a 1 bp insertion, so the ends' reference lengths
// (len + ref_shift) differ between the ends and between hits; every seventh pair, the first
// included, is random (unmapped), and every eleventh has only read 1 placed.
void rep_pairs(std::string *f1, std::string *f2)
{
    RepIndex &ix = rep_fixture();
    std::vector<std::string> r1, r2;
    const int len1 = 36, len2 = 40;
    uint32_t state = 3;
    for (int pos = 0, i = 0; pos + 342 < (int)ix.contig.size(); pos += 13, ++i) {
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
    *f1 = write_fastq(ix.dir + "/rep1.fq", r1);
    *f2 = write_fastq(ix.dir + "/rep2.fq", r2);
}

// Reads the two files into *b and searches them; with `searched` (the same files, already
// loaded and searched with the same options) its search results are copied instead, since the
// search is deterministic and repeating it would dominate the tests' run time.
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
                    REQUIRE(q->aln != nullptr);
                    memcpy(q->aln, s->aln, (size_t)s->n_aln * sizeof(ba::bwt_aln1_t));
                }
            }
        return;
    }
    RepIndex &ix = rep_fixture();
    ba::aln_search_batch(ix.bwt, b->n, b->seqs[0], &gap, 1);
    ba::aln_search_batch(ix.bwt, b->n, b->seqs[1], &gap, 1);
}

// Every field the position, pairing and extra-hit stages write, for one read.
std::string read_state(const ba::bwa_seq_t *p)
{
    std::ostringstream o;
    o << p->type << ' ' << p->pos << ' ' << p->strand << ' ' << p->mapQ << ' ' << p->seQ << ' '
      << p->extra_flag << ' ' << p->n_mm << ' ' << p->n_gapo << ' ' << p->n_gape << ' ' << p->score
      << ' ' << p->c1 << ' ' << p->c2 << ' ' << p->sa << " m" << p->n_multi;
    for (int k = 0; k < p->n_multi; ++k) {
        const ba::bwt_multi1_t *q = p->multi + k;
        o << " [" << q->pos << ' ' << q->strand << ' ' << q->gap << ' ' << q->mm << ' ' << q->ref_shift << ']';
    }
    return o.str();
}

std::vector<std::string> batch_state(const Batch &b)
{
    std::vector<std::string> v;
    for (int i = 0; i < b.n; ++i)
        for (int j = 0; j < 2; ++j) v.push_back(read_state(b.seqs[j] + i));
    return v;
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
    std::string f1, f2;
    rep_pairs(&f1, &f2);
    Batch a, b;
    load_rep_batch(*gap, f1, f2, &a);
    load_rep_batch(*gap, f1, f2, &b, &a);
    srand48(ix.bns()->seed);
    ba::bwa_pe_aln2seq(a.n, a.seqs);
    srand48(ix.bns()->seed);
    ba::bwa_pe_aln2seq(b.n, b.seqs);
    REQUIRE(batch_state(a) == batch_state(b));
    ref_pe_se_pos(ix.bns(), ix.bwt, a.n, a.seqs, gap);
    ba::bwa_pe_se_pos(ix.bns(), ix.bwt, b.n, b.seqs, gap);
    int n_mapped = 0;
    for (int i = 0; i < a.n; ++i) n_mapped += a.seqs[0][i].type != BWA_TYPE_NO_MATCH;
    CHECK(n_mapped > a.n / 2);
    CHECK(batch_state(a) == batch_state(b));
    free(gap);
}

TEST_CASE("aln_se: batched bwa_cal_pac_pos equals upstream's loop"
          * doctest::test_suite("unit/aln_pe")) {
    RepIndex &ix = rep_fixture();
    ba::bwase_initialize();
    ba::gap_opt_t *gap = ba::gap_init_opt();
    std::string f1, f2;
    rep_pairs(&f1, &f2);
    Batch searched;
    load_rep_batch(*gap, f1, f2, &searched);
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
            for (int i = 0; i < a.n; ++i)
                REQUIRE_MESSAGE(read_state(a.seqs[j] + i) == read_state(b.seqs[j] + i), "max_xa " << max_xa << " file " << j << " read " << i);
        }
    }
    free(gap);
}

TEST_CASE("aln_pe: batched pairing equals upstream's loop for every grouping of the pairs"
          * doctest::test_suite("unit/aln_pe")) {
    RepIndex &ix = rep_fixture();
    ba::bwase_initialize();
    ba::gap_opt_t *gap = ba::gap_init_opt();
    std::string f1, f2;
    rep_pairs(&f1, &f2);
    Batch searched; // searched once; every run below starts from a fresh copy of its results
    load_rep_batch(*gap, f1, f2, &searched);
    struct Opts { int max_occ, n_multi, N_multi; };
    // default; a max_occ that sends the repeat pairs down upstream's `continue` (which also
    // skips their extra hits), with few and with many extra hits allowed; no extra hits; many
    const Opts opts[] = {{100000, 3, 10}, {20, 3, 10}, {20, 30, 200}, {100000, 0, 0}, {100000, 30, 200}};
    for (const Opts &op : opts) {
        ba::pe_opt_t *pe = ba::bwa_init_pe_opt();
        pe->max_occ = op.max_occ, pe->n_multi = op.n_multi, pe->N_multi = op.N_multi;
        std::vector<std::string> want;
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
            std::vector<std::string> got = batch_state(b);
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
            for (size_t r = 0; r < want.size(); ++r)
                REQUIRE_MESSAGE(got[r] == want[r], "group rows " << rows << " max_occ " << op.max_occ << " read " << r);
        }
        free(pe);
    }
    free(gap);
}

namespace {

// Runs aln_run_pe (the shipped driver) with batches of batch_size pairs. Returns its exit code and
// the SAM record lines (no header) in output order; *err (if given) gets what it printed on stderr.
int run_pe_lines(int batch_size, const std::string &f1, const std::string &f2, const std::string &out_path,
                 std::vector<std::string> *lines, std::string *err = nullptr)
{
    PeIndex &ix = index_fixture();
    ba::aln_opt_t o;
    ba::aln_opt_init(&o);
    o.batch_size = batch_size;
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
    const std::string f1 = write_fastq(ix.dir + "/carry1.fq", a1), f2 = write_fastq(ix.dir + "/carry2.fq", a2);

    std::vector<int> alone = run_pe_flags(200, write_fastq(ix.dir + "/tail1.fq", t1), write_fastq(ix.dir + "/tail2.fq", t2),
                                          ix.dir + "/tail.sam");
    REQUIRE(alone.size() == 10);
    for (int f : alone) CHECK((f & 0x2) != 0);

    std::vector<int> two = run_pe_flags(200, f1, f2, ix.dir + "/carry_200.sam");
    REQUIRE(two.size() == 410);
    for (size_t i = 400; i < 410; ++i) {
        CAPTURE(i);
        CHECK((two[i] & 0x2) == 0);
    }
    // the first batch is flagged as in a single-batch run of the same pairs
    std::vector<int> one = run_pe_flags(0x40000, f1, f2, ix.dir + "/carry_one.sam");
    REQUIRE(one.size() == 410);
    for (size_t i = 0; i < 400; ++i) CHECK((one[i] & 0x2) == (two[i] & 0x2));
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
    std::vector<std::string> r1, r2;
    well_placed(ix.contig, 10, 36, 100, &r1, &r2);
    std::vector<std::string> lines;
    std::string err;
    SUBCASE("pair names differ in the second batch: output ends with that pair, exit 1") {
        CHECK(run_pe_lines(3, write_fastq(ix.dir + "/ra_n1.fq", r1), write_fastq(ix.dir + "/ra_n2.fq", r2, 4),
                           ix.dir + "/ra_n.sam", &lines, &err) == 1);
        std::vector<std::string> want = pair_names(4);
        want.push_back("p4"), want.push_back("x4");
        CHECK(sam_column(lines, 0) == want);
        CHECK(err.find("paired reads have different names: \"p4\", \"x4\"") != std::string::npos);
    }
    SUBCASE("read 2 file ends first: the pairs up to its last read, exit 0") {
        std::vector<std::string> s2(r2.begin(), r2.begin() + 7);
        CHECK(run_pe_lines(3, write_fastq(ix.dir + "/ra_a1.fq", r1), write_fastq(ix.dir + "/ra_a2.fq", s2),
                           ix.dir + "/ra_a.sam", &lines, &err) == 0);
        CHECK(sam_column(lines, 0) == pair_names(7));
        CHECK(err.find("read counts differ in this batch: R1 3, R2 1") != std::string::npos);
        CHECK(err.find("read 2 file ended before read 1 file") != std::string::npos);
        CHECK(err.find("read 1 file ended before read 2 file") == std::string::npos);
    }
    SUBCASE("read 1 file ends first: the pairs up to its last read, exit 0, and R2's rest is noticed") {
        std::vector<std::string> s1(r1.begin(), r1.begin() + 7);
        CHECK(run_pe_lines(3, write_fastq(ix.dir + "/ra_b1.fq", s1), write_fastq(ix.dir + "/ra_b2.fq", r2),
                           ix.dir + "/ra_b.sam", &lines, &err) == 0);
        CHECK(sam_column(lines, 0) == pair_names(7));
        // R2 is read only after a non-empty R1 batch, so p9 is still unread when R1 ends: the
        // driver reads it and warns (a read-ahead that read R2 past R1's end would not)
        CHECK(err.find("read 1 file ended before read 2 file") != std::string::npos);
    }
    SUBCASE("equal files in batches of 3 and in one batch: the same records") {
        const std::string f1 = write_fastq(ix.dir + "/ra_e1.fq", r1), f2 = write_fastq(ix.dir + "/ra_e2.fq", r2);
        std::vector<std::string> one;
        CHECK(run_pe_lines(3, f1, f2, ix.dir + "/ra_e3.sam", &lines) == 0);
        CHECK(run_pe_lines(0x40000, f1, f2, ix.dir + "/ra_e.sam", &one) == 0);
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
// written are those of the batches read before the failing read (counted with the reader).
TEST_CASE("aln_pe: a gzip read error in the next R2 batch still lets the current batch be written"
          * doctest::test_suite("unit/aln_pe")) {
    PeIndex &ix = index_fixture();
    std::vector<std::string> r1, r2;
    well_placed(ix.contig, 200, 36, 100, &r1, &r2);
    std::vector<std::string> a1, a2;
    for (int k = 0; k < 30; ++k) a1.insert(a1.end(), r1.begin(), r1.end()), a2.insert(a2.end(), r2.begin(), r2.end());
    const std::string f1 = write_fastq(ix.dir + "/gz1.fq", a1);
    const std::string f2 = ix.dir + "/gz2.fq.gz";
    {
        gzFile gz = gzopen(f2.c_str(), "wb");
        REQUIRE(gz != nullptr);
        for (size_t i = 0; i < a2.size(); ++i) {
            std::string rec = "@p" + std::to_string(i) + "\n" + a2[i] + "\n+\n" + std::string(a2[i].size(), 'I') + "\n";
            REQUIRE(gzwrite(gz, rec.data(), (unsigned)rec.size()) == (int)rec.size());
        }
        REQUIRE(gzclose(gz) == Z_OK);
        const unsigned char bad[] = {0x1f, 0x8b, 0x08, 0, 0, 0, 0, 0, 0, 0x03, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff};
        FILE *fp = fopen(f2.c_str(), "ab");
        REQUIRE(fp != nullptr);
        REQUIRE(fwrite(bad, 1, sizeof bad, fp) == sizeof bad);
        fclose(fp);
    }
    const int batch = 1000;
    int expected_pairs = 0, n_batches = 0;
    {
        ba::bwa_seqio_t *bs = ba::bwa_seq_open(f2.c_str());
        REQUIRE(bs != nullptr);
        ba::bwa_seq_defer_read_errors(bs);
        int n;
        ba::bwa_seq_t *seqs;
        while ((seqs = ba::bwa_read_seq(bs, batch, &n, 0, 0)) != nullptr) {
            expected_pairs += n, ++n_batches;
            ba::bwa_free_read_seq(n, seqs);
        }
        REQUIRE(ba::bwa_seq_read_failed(bs));
        ba::bwa_seq_close(bs);
    }
    REQUIRE(n_batches >= 2);
    std::vector<std::string> lines;
    std::string err;
    CHECK(run_pe_lines(batch, f1, f2, ix.dir + "/gz.sam", &lines, &err) == 1);
    CHECK((int)lines.size() == 2 * expected_pairs);
    CHECK(err.find("[gzread] ") != std::string::npos);
}
