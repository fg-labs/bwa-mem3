// Unit tests for the bwa backtracking search port (src/aln/bwtgap.{h,cpp}, bwtaln.cpp).
//
// The reads come from the shared synthetic reference (aln_test_ref.h) through the
// ported FASTQ reader, so they have bwa's in-memory form (seq reversed, rseq reverse
// complemented). The byte-for-byte comparison against stock `bwa aln` .sai output is
// test/regression/aln_sai_parity.sh; these cases pin the search's behaviour on a
// reference small enough to reason about.

#include "doctest/doctest.h"

#include "aln_test_ref.h"
#include "../../src/FMI_search.h"
#include "../../src/aln/aln_bwt.h"
#include "../../src/aln/bwaseqio.h"
#include "../../src/aln/bwtgap.h"
#include "../../src/bwa.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>
#include <vector>
#include <unistd.h>

namespace {

namespace ba = bwa_aln;

struct SearchIndex {
    std::string dir, prefix;
    std::unique_ptr<FMI_search> fmi;
    ba::bwt_t *bwt = nullptr;

    SearchIndex()
    {
        char tmpl[] = "/tmp/aln_search_unit_XXXXXX";
        REQUIRE(mkdtemp(tmpl) != nullptr);
        dir = tmpl;
        prefix = dir + "/ref.fa";
        REQUIRE(aln_test::write_test_fasta(prefix));
        REQUIRE(bwa_idx_build(prefix.c_str(), prefix.c_str(), 0, 3) == 0);
        fmi.reset(new FMI_search(prefix.c_str()));
        fmi->load_index(false, 1);
        bwt = ba::bwt_from_fmi(fmi.get());
    }
    ~SearchIndex()
    {
        ba::bwt_destroy(bwt);
        fmi.reset();
        std::string cmd = "rm -rf '" + dir + "'";
        if (system(cmd.c_str()) != 0) std::fprintf(stderr, "warning: could not remove %s\n", dir.c_str());
    }
};

SearchIndex &index_fixture()
{
    static SearchIndex idx;
    return idx;
}

// Reads for a batch, loaded through bwa_read_seq.
struct Batch {
    ba::bwa_seq_t *seqs = nullptr;
    int n = 0;
    explicit Batch(const std::vector<std::string> &reads)
    {
        char tmpl[] = "/tmp/aln_search_fq_XXXXXX";
        int fd = mkstemp(tmpl);
        REQUIRE(fd >= 0);
        close(fd);
        FILE *fp = fopen(tmpl, "w");
        REQUIRE(fp != nullptr);
        for (size_t i = 0; i < reads.size(); ++i)
            fprintf(fp, "@r%zu\n%s\n+\n%s\n", i, reads[i].c_str(), std::string(reads[i].size(), 'I').c_str());
        fclose(fp);
        ba::bwa_seqio_t *bs = ba::bwa_seq_open(tmpl);
        REQUIRE(bs != nullptr);
        seqs = ba::bwa_read_seq(bs, 0x40000, &n, 0, 0);
        ba::bwa_seq_close(bs);
        remove(tmpl);
        REQUIRE(n == (int)reads.size());
    }
    ~Batch() { ba::bwa_free_read_seq(n, seqs); }
    Batch(const Batch &) = delete;
    Batch &operator=(const Batch &) = delete;
};

ba::gap_opt_t make_opt(int max_diff, int max_gapo)
{
    ba::gap_opt_t *o = ba::gap_init_opt();
    ba::gap_opt_t opt = *o;
    free(o);
    opt.max_diff = max_diff;
    opt.fnr = -1.0; // fixed number of differences, as `-n INT`
    opt.max_gapo = max_gapo;
    return opt;
}

const std::string &c1() { static const std::string s = aln_test::test_contigs()[0].second; return s; }

// 32-mer well outside the c1/c2 reverse-complement copy (c1[100:400]).
std::string exact32() { return c1().substr(450, 32); }

// 40-mer with c1[480] deleted: a 1-bp deletion in the read relative to the reference.
std::string deletion40() { return c1().substr(460, 20) + c1().substr(481, 20); }

bool same_aln(const ba::bwt_aln1_t &a, const ba::bwt_aln1_t &b)
{
    return a.n_mm == b.n_mm && a.n_gapo == b.n_gapo && a.n_gape == b.n_gape && a.score == b.score
        && a.n_ins == b.n_ins && a.n_del == b.n_del && a.k == b.k && a.l == b.l;
}

} // namespace

TEST_CASE("aln_search: exact 32-mer from c1 gives one hit, score 0, n_mm 0"
          * doctest::test_suite("unit/aln_search")) {
    SearchIndex &ix = index_fixture();
    Batch b({exact32()});
    ba::gap_opt_t opt = make_opt(2, 1);
    ba::aln_search_batch(ix.bwt, b.n, b.seqs, &opt, 1);
    REQUIRE(b.seqs[0].n_aln == 1);
    const ba::bwt_aln1_t &a = b.seqs[0].aln[0];
    CHECK((int)a.score == 0);
    CHECK((int)a.n_mm == 0);
    CHECK((int)a.n_gapo == 0);
    CHECK(a.l - a.k + 1 == 1);
    CHECK((int)b.seqs[0].type == BWA_TYPE_NO_MATCH);
    CHECK(b.seqs[0].sa == 0);
    CHECK((uint64_t)b.seqs[0].c1 == 0);
    CHECK((uint64_t)b.seqs[0].c2 == 0);
}

TEST_CASE("aln_search: one planted mismatch gives n_mm 1 at -n 2"
          * doctest::test_suite("unit/aln_search")) {
    SearchIndex &ix = index_fixture();
    std::string r = exact32();
    r[16] = r[16] == 'A' ? 'C' : 'A';
    Batch b({r});
    ba::gap_opt_t opt = make_opt(2, 1);
    ba::aln_search_batch(ix.bwt, b.n, b.seqs, &opt, 1);
    REQUIRE(b.seqs[0].n_aln >= 1);
    CHECK((int)b.seqs[0].aln[0].n_mm == 1);
    CHECK((int)b.seqs[0].aln[0].n_gapo == 0);
    CHECK((int)b.seqs[0].aln[0].score == opt.s_mm);
}

TEST_CASE("aln_search: one planted 1-bp deletion found with -o 1, not with -o 0"
          * doctest::test_suite("unit/aln_search")) {
    SearchIndex &ix = index_fixture();
    Batch with_gap({deletion40()});
    ba::gap_opt_t o1 = make_opt(2, 1);
    ba::aln_search_batch(ix.bwt, with_gap.n, with_gap.seqs, &o1, 1);
    REQUIRE(with_gap.seqs[0].n_aln >= 1);
    CHECK((int)with_gap.seqs[0].aln[0].n_gapo == 1);
    CHECK((int)with_gap.seqs[0].aln[0].n_del == 1);
    CHECK((int)with_gap.seqs[0].aln[0].n_mm == 0);

    Batch no_gap({deletion40()});
    ba::gap_opt_t o0 = make_opt(2, 0);
    ba::aln_search_batch(ix.bwt, no_gap.n, no_gap.seqs, &o0, 1);
    CHECK(no_gap.seqs[0].n_aln == 0);
}

TEST_CASE("aln_search: AC-repeat read gives c1 > 1 (the repeat type is set later, in aln2seq)"
          * doctest::test_suite("unit/aln_search")) {
    SearchIndex &ix = index_fixture();
    std::string r;
    for (int i = 0; i < 16; ++i) r += "AC";
    Batch b({r});
    ba::gap_opt_t opt = make_opt(2, 1);
    ba::aln_search_batch(ix.bwt, b.n, b.seqs, &opt, 1);
    REQUIRE(b.seqs[0].n_aln >= 1);
    // The search records the number of top hits as the SA interval's width; bwa's
    // c1 is that sum over the best-scoring alignments (computed after the search).
    uint64_t top = 0;
    for (int i = 0; i < b.seqs[0].n_aln; ++i)
        if (b.seqs[0].aln[i].score == b.seqs[0].aln[0].score)
            top += b.seqs[0].aln[i].l - b.seqs[0].aln[i].k + 1;
    CHECK(top > 1);
}

TEST_CASE("aln_search: max_gapo cap uses the batch's longest read"
          * doctest::test_suite("unit/aln_search")) {
    // Stock formula (bwtaln.c): the differences allowed for an l-bp read at error rate
    // 0.02 and miss probability 0.04. Evaluated on paper: 20 bp -> 2, 100 bp -> 5.
    CHECK(ba::bwa_cal_maxdiff(20, BWA_AVG_ERR, 0.04) == 2);
    CHECK(ba::bwa_cal_maxdiff(100, BWA_AVG_ERR, 0.04) == 5);

    ba::gap_opt_t opt = ba::gap_init_opt()[0];
    opt.fnr = 0.04;
    opt.max_diff = -1;
    opt.max_gapo = 4; // `-n 0.04 -o 4`

    const std::string r20 = c1().substr(450, 20), r100 = c1().substr(480, 100);
    Batch alone(std::vector<std::string>{r20});
    Batch mixed(std::vector<std::string>{r20, r100});
    Batch longer(std::vector<std::string>{r100});

    // The cap follows the longest read of the batch handed to the search, not the read
    // being searched: 2 for the lone 20 bp read, 4 once a 100 bp read shares its batch.
    CHECK(ba::aln_batch_max_gapo(alone.n, alone.seqs, &opt) == 2);
    CHECK(ba::aln_batch_max_gapo(mixed.n, mixed.seqs, &opt) == 4);
    CHECK(ba::aln_batch_max_gapo(longer.n, longer.seqs, &opt) == 4);
    // A smaller -o is never raised, and a fixed -n bounds it independent of the batch.
    opt.max_gapo = 1;
    CHECK(ba::aln_batch_max_gapo(mixed.n, mixed.seqs, &opt) == 1);
    opt.fnr = -1.0; opt.max_diff = 3; opt.max_gapo = 4;
    CHECK(ba::aln_batch_max_gapo(alone.n, alone.seqs, &opt) == 3);
    CHECK(ba::aln_batch_max_gapo(mixed.n, mixed.seqs, &opt) == 3);

    SearchIndex &ix = index_fixture();
    // The cap is applied through the search: with the same options the short read's
    // result does not depend on who shares its batch (its own differences allowed,
    // bwa_cal_maxdiff(20) == 2, bound it first), and the 20 bp read with two
    // 1-bp deletions is found in both.
    std::string two_del = c1().substr(450, 8) + c1().substr(459, 10) + c1().substr(470, 12);
    opt = ba::gap_init_opt()[0];
    opt.fnr = 0.04; opt.max_diff = -1; opt.max_gapo = 4; opt.indel_end_skip = 5;
    Batch a2(std::vector<std::string>{two_del});
    Batch m2(std::vector<std::string>{two_del, r100});
    ba::aln_search_batch(ix.bwt, a2.n, a2.seqs, &opt, 1);
    ba::aln_search_batch(ix.bwt, m2.n, m2.seqs, &opt, 1);
    CHECK(a2.seqs[0].n_aln >= 1); // found (with two gap opens) in both
    REQUIRE(a2.seqs[0].n_aln == m2.seqs[0].n_aln);
    for (int i = 0; i < a2.seqs[0].n_aln; ++i) CHECK(same_aln(a2.seqs[0].aln[i], m2.seqs[0].aln[i]));
}

TEST_CASE("aln_search: name/seq/rseq/qual survive the search unchanged"
          * doctest::test_suite("unit/aln_search")) {
    SearchIndex &ix = index_fixture();
    std::string mm = exact32();
    mm[10] = mm[10] == 'G' ? 'T' : 'G';
    Batch b({exact32(), mm, deletion40()});
    struct Snap { std::string name, seq, rseq, qual; };
    std::vector<Snap> before;
    for (int i = 0; i < b.n; ++i) {
        const ba::bwa_seq_t &s = b.seqs[i];
        before.push_back({s.name, std::string((const char *)s.seq, s.len),
                          std::string((const char *)s.rseq, s.len),
                          std::string((const char *)s.qual, s.len)});
    }
    ba::gap_opt_t opt = make_opt(2, 1);
    ba::aln_search_batch(ix.bwt, b.n, b.seqs, &opt, 2);
    for (int i = 0; i < b.n; ++i) {
        const ba::bwa_seq_t &s = b.seqs[i];
        REQUIRE(s.name != nullptr);
        REQUIRE(s.seq != nullptr);
        REQUIRE(s.rseq != nullptr);
        REQUIRE(s.qual != nullptr);
        CHECK(before[i].name == s.name);
        CHECK(before[i].seq == std::string((const char *)s.seq, s.len));
        CHECK(before[i].rseq == std::string((const char *)s.rseq, s.len));
        CHECK(before[i].qual == std::string((const char *)s.qual, s.len));
    }
}

TEST_CASE("aln_search: n_threads 1 and 4 give identical aln arrays"
          * doctest::test_suite("unit/aln_search")) {
    SearchIndex &ix = index_fixture();
    std::vector<std::string> reads;
    for (int i = 0; i < 40; ++i) {
        std::string r = c1().substr(i * 13, 32 + (i % 5));
        if (i % 3 == 1) r[12] = r[12] == 'A' ? 'T' : 'A';
        if (i % 7 == 2) r = deletion40();
        if (i % 11 == 3) { r.clear(); for (int j = 0; j < 16; ++j) r += "AC"; }
        reads.push_back(r);
    }
    Batch one(reads), four(reads);
    ba::gap_opt_t opt = make_opt(2, 1);
    ba::aln_search_batch(ix.bwt, one.n, one.seqs, &opt, 1);
    ba::aln_search_batch(ix.bwt, four.n, four.seqs, &opt, 4);
    int found = 0;
    for (int i = 0; i < one.n; ++i) {
        REQUIRE(one.seqs[i].n_aln == four.seqs[i].n_aln);
        found += one.seqs[i].n_aln > 0;
        for (int j = 0; j < one.seqs[i].n_aln; ++j) CHECK(same_aln(one.seqs[i].aln[j], four.seqs[i].aln[j]));
    }
    CHECK(found > 20); // the comparison is not vacuous
}

// ---------------------------------------------------------------------------------------
// bwt_match_gap (upstream's loop) and gap_job_start/gap_job_run (the resumable copy the
// aligner runs) must give the same hits for the same read and options.

namespace {

struct GapResult {
    std::vector<ba::bwt_aln1_t> aln;
};

// Prepares one read as bwa_cal_sa_reg_gap does, then runs it through both loops.
// `codes` are nt4 codes (4 = N); the search takes the complemented read.
// `batch_max_gapo` >= 0 is the driver's batch-wide cap on max_gapo (aln_batch_max_gapo);
// -1 caps it with this read's own max_diff.
void both_loops(const ba::bwt_t *bwt, const std::vector<uint8_t> &codes, const ba::gap_opt_t &base,
                GapResult &ref, GapResult &job, int batch_max_gapo = -1)
{
    const int len = (int)codes.size();
    ba::gap_opt_t opt = base;
    if (opt.fnr > 0.0) opt.max_diff = ba::bwa_cal_maxdiff(len, BWA_AVG_ERR, opt.fnr);
    if (batch_max_gapo >= 0) opt.max_gapo = batch_max_gapo;
    else if (opt.max_diff < opt.max_gapo) opt.max_gapo = opt.max_diff;
    ba::gap_opt_t ropt = opt;
    ropt.seed_len = opt.seed_len < len ? opt.seed_len : 0x7fffffff;

    std::vector<ba::bwt_width_t> w(len + 1), seed_w(opt.seed_len + 1);
    std::vector<uint8_t> cseq(len + 1);
    ba::bwt_cal_width(bwt, len, codes.data(), w.data());
    if (len > opt.seed_len)
        ba::bwt_cal_width(bwt, opt.seed_len, codes.data() + (len - opt.seed_len), seed_w.data());
    for (int k = 0; k < len; ++k) cseq[k] = codes[k] > 3 ? 4 : 3 - codes[k];
    ba::bwt_width_t *sw = len <= opt.seed_len ? nullptr : seed_w.data();

    ba::gap_stack_t *st1 = ba::gap_init_stack(opt.max_diff, opt.max_gapo, opt.max_gape, &opt);
    int n1 = 0;
    ba::bwt_aln1_t *a1 = ba::bwt_match_gap(const_cast<ba::bwt_t *>(bwt), len, cseq.data(), w.data(), sw, &ropt, &n1, st1);
    ref.aln.assign(a1, a1 + n1);
    free(a1);
    ba::gap_destroy_stack(st1);

    // gap_shadow narrows the widths in place while a search runs, so the job gets fresh ones
    ba::bwt_cal_width(bwt, len, codes.data(), w.data());
    if (len > opt.seed_len)
        ba::bwt_cal_width(bwt, opt.seed_len, codes.data() + (len - opt.seed_len), seed_w.data());
    ba::gap_job_t j;
    memset(&j, 0, sizeof j);
    j.opt = ropt;
    j.len = len;
    j.seq = cseq.data();
    j.width = w.data();
    j.seed_width = sw;
    j.stack = ba::gap_init_stack(opt.max_diff, opt.max_gapo, opt.max_gape, &opt);
    if (!ba::gap_job_start(&j, bwt)) {
        int guard = 0;
        while (!ba::gap_job_run(&j, bwt)) REQUIRE(++guard < 100000000);
    }
    job.aln.assign(j.aln, j.aln + j.n_aln);
    free(j.aln);
    ba::gap_destroy_stack(j.stack);
}

std::vector<uint8_t> to_codes(const std::string &s)
{
    // bwa_read_seq keeps the read reversed in `seq`; the search complements that copy
    std::vector<uint8_t> v;
    for (auto it = s.rbegin(); it != s.rend(); ++it)
        v.push_back(*it == 'A' ? 0 : *it == 'C' ? 1 : *it == 'G' ? 2 : *it == 'T' ? 3 : 4);
    return v;
}

} // namespace

TEST_CASE("aln_search: gap_job_run gives the same hits as bwt_match_gap"
          * doctest::test_suite("unit/aln_search")) {
    SearchIndex &ix = index_fixture();
    auto contigs = aln_test::test_contigs();
    const std::string &s1 = contigs[0].second, &s2 = contigs[1].second;
    auto flip = [](std::string r, size_t i) { r[i] = r[i] == 'A' ? 'T' : 'A'; return r; };

    std::vector<std::string> reads;
    reads.push_back("");                                    // length 0
    reads.push_back("A");                                   // length 1
    reads.push_back(s1.substr(10, 5));                      // short
    reads.push_back(s1.substr(450, 32));                    // exact
    reads.push_back(s1.substr(200, 36));                    // in the reverse-complement copy
    for (int i = 0; i < 12; ++i) {                          // mismatches
        std::string r = s1.substr(i * 40, 30 + i);
        reads.push_back(flip(r, 5));
        reads.push_back(flip(flip(r, 3), r.size() - 4));
        reads.push_back(flip(flip(flip(r, 2), 11), r.size() - 3));
    }
    reads.push_back(s1.substr(30, 20) + s1.substr(52, 20));           // 2 bp deletion
    reads.push_back(s1.substr(30, 20) + "GG" + s1.substr(50, 20));    // insertion
    reads.push_back(s1.substr(300, 25) + s1.substr(327, 30));         // 2 bp deletion, longer
    reads.push_back(s1.substr(100, 24) + "T" + s1.substr(125, 24));   // deletion plus context
    reads.push_back(s1.substr(400, 18) + "ACA" + s1.substr(418, 20)); // insertion plus mismatch below
    reads.push_back(flip(s1.substr(400, 18) + "ACA" + s1.substr(418, 20), 30));
    {   // N-heavy
        std::string r = s1.substr(60, 40);
        for (size_t i : {3, 9, 10, 17, 25, 33}) r[i] = 'N';
        reads.push_back(r);
        reads.push_back(std::string(30, 'N'));
        reads.push_back("N" + s1.substr(70, 30));
        reads.push_back(s1.substr(70, 30) + "N");
        reads.push_back(s2.substr(290, 40));                // spans the N run of the reference
    }
    {   // tandem repeats: many equal suffixes
        std::string r;
        for (int i = 0; i < 16; ++i) r += "AC";
        reads.push_back(r);
        reads.push_back(flip(r, 9));
        reads.push_back(r.substr(0, 11));
    }
    {   // the exact-match tail (no differences left): reads with an N (which force a width restart,
        // so none reaches the tail with an N still ahead), a miss part-way, long tails
        std::string r = s1.substr(500, 40);
        for (size_t n : {1, 6, 20, 33, 38}) {
            std::string a = r; a[n] = 'N';
            reads.push_back(a);
            reads.push_back(flip(a, n < 20 ? 30 : 4));
            reads.push_back(flip(flip(a, n < 20 ? 30 : 4), 15));
        }
        uint32_t x = 12345;
        for (int t = 0; t < 6; ++t) { // not in the reference: every tail misses
            std::string q;
            for (int i = 0; i < 28 + 4 * t; ++i) { x = x * 1103515245u + 12345u; q += "ACGT"[(x >> 16) & 3]; }
            reads.push_back(q);
            reads.push_back(s1.substr(410 + 7 * t, 18) + q.substr(0, 10 + t)); // matches, then misses
        }
        reads.push_back(flip(s1.substr(420, 150), 75));     // a long tail
        reads.push_back(flip(flip(s1.substr(470, 90), 3), 85));
    }

    struct Case { const char *name; int max_diff; int max_gapo; int max_entries; int max_top2; int mode_or; int seed_len; int max_seed_diff; double fnr; };
    const Case cases[] = {
        {"defaults",       -1, 1, 2000000, 30,         0,                 32, 2, 0.04},
        {"fixed -n 3 -o 2", 3, 2, 2000000, 30,         0,                 32, 2, -1.0},
        {"tiny -m",        -1, 1, 5,       30,         0,                 32, 2, 0.04},
        {"-m 1",            3, 1, 1,       30,         0,                 32, 2, -1.0},
        {"-m 3",            3, 2, 3,       30,         0,                 32, 2, -1.0},
        {"-m 10",           3, 2, 10,      30,         0,                 32, 2, -1.0},
        {"-m 40",           4, 2, 40,      30,         0,                 32, 2, -1.0},
        {"small -R",       -1, 1, 2000000, 1,          0,                 32, 2, 0.04},
        {"-R 0",            3, 1, 2000000, 0,          0,                 32, 2, -1.0},
        {"NONSTOP",         3, 2, 2000000, 0x7fffffff, BWA_MODE_NONSTOP,  32, 2, -1.0},
        {"NONSTOP -n 0.04",-1, 1, 2000000, 0x7fffffff, BWA_MODE_NONSTOP,  32, 2, 0.04},
        {"short seed",      3, 1, 2000000, 30,         0,                 12, 1, -1.0},
        {"aDNA-like",      -1, 2, 2000000, 30,         0,                 1024, 2, 0.04},
        {"-L log gap",      3, 2, 2000000, 30,         BWA_MODE_LOGGAP,   32, 2, -1.0},
    };
    long compared = 0, with_hits = 0;
    for (const Case &c : cases) {
        ba::gap_opt_t *o = ba::gap_init_opt();
        ba::gap_opt_t opt = *o;
        free(o);
        opt.max_diff = c.max_diff;
        opt.fnr = (float)c.fnr;
        opt.max_gapo = c.max_gapo;
        opt.max_entries = c.max_entries;
        opt.max_top2 = c.max_top2;
        opt.mode |= c.mode_or;
        opt.seed_len = c.seed_len;
        opt.max_seed_diff = c.max_seed_diff;
        for (size_t i = 0; i < reads.size(); ++i) {
            INFO("case " << std::string(c.name) << ", read " << i << " (" << reads[i] << ")");
            GapResult ref, job;
            both_loops(ix.bwt, to_codes(reads[i]), opt, ref, job);
            REQUIRE(ref.aln.size() == job.aln.size());
            for (size_t k = 0; k < ref.aln.size(); ++k) CHECK(same_aln(ref.aln[k], job.aln[k]));
            ++compared;
            with_hits += !ref.aln.empty();
        }
    }
    CHECK(compared > 400);
    CHECK(with_hits > compared / 4); // not vacuous
}

// The aligner's driver (aln_search_batch: ALN_GROUP reads in flight per thread, each slot
// refilled as its read finishes) gives every read exactly bwt_match_gap's hits.
TEST_CASE("aln_search: aln_search_batch gives the same hits as bwt_match_gap for each read"
          * doctest::test_suite("unit/aln_search")) {
    SearchIndex &ix = index_fixture();
    auto contigs = aln_test::test_contigs();
    const std::string &s1 = contigs[0].second, &s2 = contigs[1].second;
    auto flip = [](std::string r, size_t i) { r[i] = r[i] == 'A' ? 'T' : 'A'; return r; };
    std::vector<std::string> reads;
    uint32_t x = 777;
    auto rnd = [&x](uint32_t n) { x = x * 1103515245u + 12345u; return (x >> 16) % n; };
    for (int i = 0; i < 150; ++i) { // mixed lengths and kinds, so slots finish at different times
        const std::string &src = (i % 4 == 3) ? s2 : s1;
        int len = 20 + (int)rnd(60);
        std::string r = src.substr(rnd((uint32_t)(src.size() - len)), len);
        switch (i % 6) {
        case 1: r = flip(r, rnd(len)); break;
        case 2: r = flip(flip(r, rnd(len)), rnd(len)); break;
        case 3: r[rnd(len)] = 'N'; r = flip(r, rnd(len)); break;
        case 4: r = r.substr(0, len / 2) + r.substr(len / 2 + 1 + rnd(2)); break; // deletion
        case 5: { std::string q; for (int k = 0; k < len; ++k) q += "ACGT"[rnd(4)]; r = (i % 12 == 5) ? q : r.substr(0, len / 2) + q.substr(len / 2); } break;
        default: break;
        }
        reads.push_back(r);
    }
    for (int i = 0; i < 8; ++i) { std::string r; for (int j = 0; j < 12 + i; ++j) r += "AC"; reads.push_back(r); }
    // edge lengths and contents: a single base, all-N, exactly seed_len (32) and seed_len + 1
    reads.push_back("A");
    reads.push_back("N");
    reads.push_back(std::string(25, 'N'));
    reads.push_back(s1.substr(200, 32));
    reads.push_back(s1.substr(200, 33));
    reads.push_back(flip(s1.substr(240, 32), 10));
    reads.push_back(flip(s1.substr(240, 33), 10));

    struct Case { const char *name; int max_diff; int max_gapo; int max_entries; double fnr; };
    const Case cases[] = {
        {"-n 2 -o 1", 2, 1, 2000000, -1.0},
        {"-n 3 -o 2", 3, 2, 2000000, -1.0},
        {"defaults", -1, 1, 2000000, 0.04},
        {"-n 3 -m 5", 3, 1, 5, -1.0},
    };
    long with_hits = 0, compared = 0;
    for (const Case &c : cases) {
        ba::gap_opt_t opt = make_opt(c.max_diff, c.max_gapo);
        opt.fnr = (float)c.fnr;
        opt.max_entries = c.max_entries;
        for (int nt : {1, 3}) {
            Batch b(reads);
            ba::aln_search_batch(ix.bwt, b.n, b.seqs, &opt, nt);
            const int cap = ba::aln_batch_max_gapo(b.n, b.seqs, &opt); // what the driver applies
            for (int i = 0; i < b.n; ++i) {
                INFO("case " << std::string(c.name) << ", -t " << nt << ", read " << i << " (" << reads[i] << ")");
                GapResult ref, job;
                both_loops(ix.bwt, to_codes(reads[i]), opt, ref, job, cap);
                REQUIRE(b.seqs[i].n_aln == (int)ref.aln.size());
                for (int k = 0; k < b.seqs[i].n_aln; ++k) CHECK(same_aln(ref.aln[k], b.seqs[i].aln[k]));
                ++compared;
                with_hits += !ref.aln.empty();
            }
        }
    }
    CHECK(compared == 8L * (long)reads.size());
    CHECK(with_hits > compared / 3); // not vacuous
}

// The resumable width walk (width_walk_init/width_walk_run, which the driver interleaves) must
// write exactly what bwt_cal_width writes, and end on the same bid.
TEST_CASE("aln_search: width_walk_run matches bwt_cal_width"
          * doctest::test_suite("unit/aln_search")) {
    SearchIndex &ix = index_fixture();
    auto contigs = aln_test::test_contigs();
    const std::string &s1 = contigs[0].second;
    auto flip = [](std::string r, size_t i) { r[i] = r[i] == 'A' ? 'T' : 'A'; return r; };
    std::vector<std::string> reads;
    reads.push_back("");                                  // len 0
    reads.push_back("A");                                 // len 1
    reads.push_back(s1.substr(10, 1));
    reads.push_back(std::string(40, 'N'));                // all N
    { std::string r = s1.substr(60, 50); r[20] = r[21] = r[22] = 'N'; reads.push_back(r); } // N run mid-read
    { std::string r = s1.substr(60, 50); r[7] = 'N'; reads.push_back(r); }
    reads.push_back(flip(s1.substr(300, 60), 30));        // mid-read miss (k > l restart)
    reads.push_back(s1.substr(300, 28) + "ACGTTGCAACGTTGCAGGCCTTAAGGCC"); // match, then off-reference
    reads.push_back(s1.substr(500, 32));                  // seed_len
    reads.push_back(s1.substr(500, 33));                  // seed_len + 1
    reads.push_back(s1.substr(400, 400));                 // long
    reads.push_back(flip(flip(s1.substr(400, 400), 100), 300));
    long restarts = 0;
    for (const std::string &r : reads) {
        INFO("read " << r);
        std::vector<uint8_t> codes = to_codes(r);
        const int len = (int)codes.size();
        std::vector<ba::bwt_width_t> a(len + 1), b(len + 1);
        for (auto &x : a) { x.w = 0; x.bid = 0; }
        for (auto &x : b) { x.w = 0; x.bid = 0; }
        int bid = ba::bwt_cal_width(ix.bwt, len, codes.data(), a.data());
        ba::width_walk_t ww;
        ba::width_walk_init(&ww, ix.bwt, len, codes.data(), b.data());
        int guard = 0;
        while (!ba::width_walk_run(&ww, ix.bwt)) REQUIRE(++guard < 10000000);
        CHECK(ww.bid == bid);
        for (int i = 0; i <= len; ++i) {
            CHECK(a[i].w == b[i].w);
            CHECK(a[i].bid == b[i].bid);
        }
        restarts += bid > 1;
    }
    CHECK(restarts >= 6); // restarts are exercised
}
