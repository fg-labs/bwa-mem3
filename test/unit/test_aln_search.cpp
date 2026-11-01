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

    SearchIndex() : SearchIndex(aln_test::test_contigs()) {}
    explicit SearchIndex(const std::vector<std::pair<std::string, std::string>> &contigs)
    {
        std::string tmpl_s = aln_test::tmp_template("aln_search_unit");
        char *tmpl = &tmpl_s[0];
        REQUIRE(mkdtemp(tmpl) != nullptr);
        dir = tmpl;
        prefix = dir + "/ref.fa";
        REQUIRE(aln_test::write_test_fasta(prefix, contigs));
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

// A reference with a diverged copy: p2 is p1[50:350] with a 1-bp deletion (at p1 150), a 1-bp
// insertion (after p1 250) and two mismatches (p1 100 and 300), so one read can align exactly or
// with mismatches to p1 and with a gap to p2.
std::vector<std::pair<std::string, std::string>> paralog_contigs()
{
    uint32_t st = 77;
    std::string p1 = aln_test::lcg_dna(st, 400), p3 = aln_test::lcg_dna(st, 300);
    std::string c = p1.substr(50, 300); // c[i] = p1[50 + i]
    auto other = [](char b) { return b == 'A' ? 'C' : 'A'; };
    c[250] = other(c[250]);                  // p1 300
    c.insert(201, 1, other(c[200]));         // after p1 250
    c.erase(100, 1);                         // p1 150
    c[50] = other(c[50]);                    // p1 100
    return {{"p1", p1}, {"p2", c}, {"p3", p3}};
}

SearchIndex &paralog_fixture()
{
    static SearchIndex idx(paralog_contigs());
    return idx;
}

// Reads for a batch, loaded through bwa_read_seq.
struct Batch {
    ba::bwa_seq_t *seqs = nullptr;
    int n = 0;
    explicit Batch(const std::vector<std::string> &reads)
    {
        std::string tmpl_s = aln_test::tmp_template("aln_search_fq");
        char *tmpl = &tmpl_s[0];
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
    // c3 is 64 x "AC" (128 bp), so a 32-mer "ACAC..." starts at every even offset 0..96: 49
    // places, all with zero differences; no other contig holds it.
    CHECK(top == 49);
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
    long exact_yields = 0; // gap_job_run calls that stopped inside an exact-match tail
    int max_entries = 0; // job only: the most entries its stack held
};

// Prepares one read as bwa_cal_sa_reg_gap does, then runs it through both loops.
// `codes` are nt4 codes (4 = N); the search takes the complemented read.
// `batch_max_gapo` >= 0 is the driver's batch-wide cap on max_gapo (aln_batch_max_gapo);
// -1 caps it with this read's own max_diff.
// `pool_keep` is the job stack's pool cap (free chunks kept); production uses ALN_POOL_KEEP.
void both_loops(const ba::bwt_t *bwt, const std::vector<uint8_t> &codes, const ba::gap_opt_t &base,
                GapResult &ref, GapResult &job, int batch_max_gapo = -1, int pool_keep = ba::ALN_POOL_KEEP)
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
    ba::gap_pool_t pool;
    ba::gap_init_pool(&pool, pool_keep);
    j.stack = ba::gap_init_cstack(opt.max_diff, opt.max_gapo, opt.max_gape, &opt, &pool);
    if (!ba::gap_job_start(&j, bwt)) {
        int guard = 0;
        while (!ba::gap_job_run(&j, bwt)) {
            REQUIRE(++guard < 100000000);
            job.exact_yields += j.pending == ba::JOB_EXACT;
        }
    }
    // a finished job has given every chunk back: none in use, at most pool_keep kept
    CHECK(j.stack->n_entries == 0);
    CHECK(pool.n_alloc == pool.n_free);
    CHECK(pool.n_free <= pool_keep);
    job.aln.assign(j.aln, j.aln + j.n_aln);
    job.max_entries = j.max_entries;
    free(j.aln);
    ba::gap_destroy_cstack(j.stack);
    ba::gap_destroy_pool(&pool);
    CHECK(pool.n_alloc == 0);
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
        {"-n 0",            0, 1, 2000000, 30,         0,                 32, 2, -1.0}, // root entry runs as an exact tail
    };
    long compared = 0, with_hits = 0, exact_yields = 0;
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
            exact_yields += job.exact_yields;
        }
    }
    CHECK(compared > 400);
    CHECK(with_hits > compared / 4); // not vacuous
    CHECK(exact_yields > compared); // the exact-match tails yield (one stop per base), not run in one call
}

// Stack-heavy searches: many differences allowed (NONSTOP keeps every branch alive), so the job's
// stack holds thousands of entries spread over its buckets, buckets drain and refill, chunks
// cross between buckets through the pool, and -m caps fall on and around the chunk size. Small
// pool caps make the pool free chunks it cannot keep.
TEST_CASE("aln_search: gap_job_run matches bwt_match_gap on large stacks and around the chunk size"
          * doctest::test_suite("unit/aln_search")) {
    SearchIndex &ix = index_fixture();
    auto contigs = aln_test::test_contigs();
    const std::string &s1 = contigs[0].second, &s2 = contigs[1].second;
    auto flip = [](std::string r, size_t i) { r[i] = r[i] == 'A' ? 'T' : 'A'; return r; };
    std::vector<std::string> reads;
    uint32_t x = 4242;
    auto rnd = [&x](uint32_t n) { x = x * 1103515245u + 12345u; return (x >> 16) % n; };
    for (int t = 0; t < 16; ++t) {
        int len = 40 + (int)rnd(60);
        const std::string &src = t & 1 ? s2 : s1;
        std::string r = src.substr(rnd((uint32_t)(src.size() - len)), len);
        for (int f = 0; f < t % 4; ++f) r = flip(r, rnd(len));
        reads.push_back(r);
        std::string q;
        for (int k = 0; k < len; ++k) q += "ACGT"[rnd(4)];
        reads.push_back(q);
    }
    { std::string r; for (int i = 0; i < 25; ++i) r += "AC"; reads.push_back(flip(r, 20)); }

    const int C = ba::GAP_CHUNK_ENTRIES;
    struct Case { const char *name; int max_diff; int max_gapo; int max_gape; int max_entries; int mode_or; int pool_keep; };
    const Case cases[] = {
        {"-n 5 -o 2 NONSTOP",           5, 2, 6, 2000000, BWA_MODE_NONSTOP, ba::ALN_POOL_KEEP},
        {"-n 6 -o 2 -e 3 NONSTOP keep 0", 6, 2, 3, 2000000, BWA_MODE_NONSTOP, 0},
        {"-n 6 -o 2 keep 2",            6, 2, 6, 2000000, 0,                2},
        {"-n 6 -o 2 -m 1",              6, 2, 6, 1,       BWA_MODE_NONSTOP, ba::ALN_POOL_KEEP},
        {"-m chunk-1",                  6, 2, 6, C - 1,   BWA_MODE_NONSTOP, 1},
        {"-m chunk",                    6, 2, 6, C,       BWA_MODE_NONSTOP, 1},
        {"-m chunk+1",                  6, 2, 6, C + 1,   BWA_MODE_NONSTOP, 1},
        {"-m 2chunk+1",                 6, 2, 6, 2 * C + 1, BWA_MODE_NONSTOP, 3},
        {"-m 4chunk+1",                 6, 2, 6, 4 * C + 1, BWA_MODE_NONSTOP, ba::ALN_POOL_KEEP},
        {"-n 8 -o 2 NONSTOP keep 4",    8, 2, 6, 2000000, BWA_MODE_NONSTOP, 4},
    };
    long compared = 0, with_hits = 0, big = 0, capped = 0, max_seen = 0;
    for (const Case &c : cases) {
        long case_capped = 0, case_big = 0;
        ba::gap_opt_t opt = make_opt(c.max_diff, c.max_gapo);
        opt.max_gape = c.max_gape;
        opt.max_entries = c.max_entries;
        opt.mode |= c.mode_or;
        opt.max_top2 = c.mode_or & BWA_MODE_NONSTOP ? 0x7fffffff : opt.max_top2;
        for (size_t i = 0; i < reads.size(); ++i) {
            INFO("case " << std::string(c.name) << ", read " << i << " (" << reads[i] << ")");
            GapResult ref, job;
            both_loops(ix.bwt, to_codes(reads[i]), opt, ref, job, -1, c.pool_keep);
            REQUIRE(ref.aln.size() == job.aln.size());
            for (size_t k = 0; k < ref.aln.size(); ++k) CHECK(same_aln(ref.aln[k], job.aln[k]));
            ++compared;
            with_hits += !ref.aln.empty();
            case_big += job.max_entries > 4 * C;
            if (job.max_entries > max_seen) max_seen = job.max_entries;
            case_capped += job.max_entries > c.max_entries;
        }
        // each case reaches what it is there for: an -m case stops on its cap at least once, and
        // an uncapped NONSTOP case builds a stack of more than 4 chunks at least once
        INFO("case " << std::string(c.name) << ": stopped by -m " << case_capped << ", over 4 chunks " << case_big);
        if (c.max_entries < 2000000) CHECK(case_capped > 0);
        else if (c.mode_or & BWA_MODE_NONSTOP) CHECK(case_big > 0);
        big += case_big;
        capped += case_capped;
    }
    MESSAGE("compared " << compared << ", with hits " << with_hits << ", stacks over 4 chunks " << big
            << ", stopped by -m " << capped << ", largest stack " << max_seen);
    // floors for this fixed-seed data set (the reads come from a seeded LCG)
    CHECK(with_hits > compared / 4);
    CHECK(big > 10);     // large stacks are exercised
    CHECK(max_seen > 16 * C);
    CHECK(capped > 20);  // and the -m stop at, below and above the chunk size
}

// Targets gap_job_run's top-of-stack forwarding (the last child of an expansion used without a
// store and a pop when it is the entry the pop would return, counted as pushed). Forwarding is
// common (the match child usually lands in the lowest bucket), so a -m cap swept over every small
// value makes the -m break fall on a forwarded entry (counted with the stack) as well as on stored
// ones; -R swept over small values does the same for the top2 break. Scoring variants change which
// child is pushed last and which bucket is lowest: with -O below -M an insertion child can sit below
// the last (mismatch) child, which must then be stored and the older entry popped; reads with an N
// give four children in the same bucket. Forwarding is output-neutral by design, so this checks
// that forwarding is correct, not that it happens: a build that never forwards passes too (and
// only runs slower).
TEST_CASE("aln_search: gap_job_run matches bwt_match_gap at every small -m and -R and under odd scoring"
          * doctest::test_suite("unit/aln_search")) {
    SearchIndex &ix = index_fixture();
    auto contigs = aln_test::test_contigs();
    const std::string &s1 = contigs[0].second, &s2 = contigs[1].second;
    auto flip = [](std::string r, size_t i) { r[i] = r[i] == 'A' ? 'T' : 'A'; return r; };
    std::vector<std::string> reads;
    uint32_t x = 31337;
    auto rnd = [&x](uint32_t n) { x = x * 1103515245u + 12345u; return (x >> 16) % n; };
    for (int t = 0; t < 10; ++t) {
        int len = 25 + (int)rnd(40);
        const std::string &src = t & 1 ? s2 : s1;
        std::string r = src.substr(rnd((uint32_t)(src.size() - len)), len);
        for (int f = 0; f < 1 + t % 3; ++f) r = flip(r, rnd(len));
        if (t % 4 == 3) r[rnd(len)] = 'N';
        if (t % 5 == 4) r = r.substr(0, len / 2) + r.substr(len / 2 + 1); // deletion
        reads.push_back(r);
    }
    { std::string r; for (int i = 0; i < 14; ++i) r += "AC"; reads.push_back(r); reads.push_back(flip(r, 7)); }
    { std::string r = s1.substr(120, 36); r[10] = 'N'; r[25] = 'N'; reads.push_back(r); }

    struct Score { const char *name; int s_mm, s_gapo, s_gape; };
    const Score scores[] = {{"stock", 3, 11, 4}, {"-M 5 -O 2 -E 1", 5, 2, 1}, {"-M 1 -O 1 -E 1", 1, 1, 1}};
    long compared = 0, with_hits = 0, capped = 0, m_compared = 0;
    for (const Score &sc : scores) {
        // -m 1..48 at -n 3 -o 2, NONSTOP off and on
        for (int nonstop = 0; nonstop < 2; ++nonstop) {
            for (int m = 1; m <= 48; ++m) {
                ba::gap_opt_t opt = make_opt(3, 2);
                opt.s_mm = sc.s_mm; opt.s_gapo = sc.s_gapo; opt.s_gape = sc.s_gape;
                opt.max_entries = m;
                if (nonstop) { opt.mode |= BWA_MODE_NONSTOP; opt.max_top2 = 0x7fffffff; }
                for (size_t i = 0; i < reads.size(); ++i) {
                    INFO("scoring " << sc.name << ", NONSTOP " << nonstop << ", -m " << m << ", read " << i << " (" << reads[i] << ")");
                    GapResult ref, job;
                    both_loops(ix.bwt, to_codes(reads[i]), opt, ref, job);
                    REQUIRE(ref.aln.size() == job.aln.size());
                    for (size_t k = 0; k < ref.aln.size(); ++k) CHECK(same_aln(ref.aln[k], job.aln[k]));
                    ++compared;
                    ++m_compared;
                    with_hits += !ref.aln.empty();
                    capped += job.max_entries > m;
                }
            }
        }
        // -R 0..12 (top2 break), at -n 2 and -n 4
        for (int top2 = 0; top2 <= 12; ++top2) {
            for (int nd : {2, 4}) {
                ba::gap_opt_t opt = make_opt(nd, 1);
                opt.s_mm = sc.s_mm; opt.s_gapo = sc.s_gapo; opt.s_gape = sc.s_gape;
                opt.max_top2 = top2;
                for (size_t i = 0; i < reads.size(); ++i) {
                    INFO("scoring " << sc.name << ", -R " << top2 << ", -n " << nd << ", read " << i << " (" << reads[i] << ")");
                    GapResult ref, job;
                    both_loops(ix.bwt, to_codes(reads[i]), opt, ref, job);
                    REQUIRE(ref.aln.size() == job.aln.size());
                    for (size_t k = 0; k < ref.aln.size(); ++k) CHECK(same_aln(ref.aln[k], job.aln[k]));
                    ++compared;
                    with_hits += !ref.aln.empty();
                }
            }
        }
    }
    MESSAGE("compared " << compared << ", with hits " << with_hits << ", stopped by -m " << capped << " of "
            << m_compared);
    CHECK(with_hits > compared / 4);
    CHECK(capped > m_compared / 2); // the -m sweep does stop most of its searches
}

// The score break (stop once the entry popped scores above the first hit's score + s_mm) only
// changes the hits when an entry above that score could still become a hit: one more difference
// than the best hit, where that difference is a gap (it costs more than s_mm). Reads with a
// mismatch or an indel next to a short repeat give such second alignments; -M 1 makes every gap
// cost more than the cut allows.
TEST_CASE("aln_search: gap_job_run matches bwt_match_gap where the score break decides the hits"
          * doctest::test_suite("unit/aln_search")) {
    SearchIndex &ix = paralog_fixture();
    auto contigs = paralog_contigs();
    const std::string &s1 = contigs[0].second, &s2 = contigs[1].second;
    uint32_t x = 2718;
    auto rnd = [&x](uint32_t n) { x = x * 1103515245u + 12345u; return (x >> 16) % n; };
    auto mutate = [&](std::string r, size_t i) { r[i] = "ACGT"[(std::string("ACGT").find(r[i]) + 1 + rnd(3)) % 4]; return r; };
    std::vector<std::string> reads;
    // 300 reads: an off-by-one score_cut (best_score + s_mm + 1) changes the hits of 12 read/option pairs
    for (int t = 0; t < 300; ++t) {
        int len = 30 + (int)rnd(30);
        const std::string &src = t % 3 == 2 ? s2 : s1;
        // mostly from the diverged span p1[50:350] (or its copy), so most reads have a second locus
        size_t start = t % 6 == 5 ? rnd((uint32_t)(src.size() - len - 4)) : (&src == &s1 ? 50 : 0) + rnd(300 - len - 4);
        std::string r = src.substr(start, len + 2);
        switch (t % 5) {
        case 4: break; // unchanged
        case 0: r = mutate(r, 8 + rnd(len - 16)); break;
        case 1: { size_t p = 8 + rnd(len - 16); r = r.substr(0, p) + r.substr(p + 1); } break; // deletion
        case 2: { size_t p = 8 + rnd(len - 16); r = r.substr(0, p) + r[p] + r.substr(p); } break; // duplicated base
        default: r = mutate(mutate(r, 6 + rnd(len / 2 - 6)), len / 2 + rnd(len / 2 - 6)); break;
        }
        reads.push_back(r.substr(0, len));
    }
    struct Case { const char *name; int max_diff, max_gapo, s_mm, s_gapo, s_gape; };
    const Case cases[] = {
        {"-n 2 -o 1", 2, 1, 3, 11, 4}, {"-n 3 -o 1", 3, 1, 3, 11, 4}, {"-n 3 -o 2", 3, 2, 3, 11, 4},
        {"-n 2 -o 1 -M 1", 2, 1, 1, 11, 4}, {"-n 3 -o 2 -M 1", 3, 2, 1, 11, 4}, {"-n 3 -o 2 -M 2 -O 5 -E 3", 3, 2, 2, 5, 3},
    };
    long compared = 0, with_hits = 0, gapped = 0;
    for (const Case &c : cases) {
        ba::gap_opt_t opt = make_opt(c.max_diff, c.max_gapo);
        opt.s_mm = c.s_mm; opt.s_gapo = c.s_gapo; opt.s_gape = c.s_gape;
        for (size_t i = 0; i < reads.size(); ++i) {
            INFO("case " << std::string(c.name) << ", read " << i << " (" << reads[i] << ")");
            GapResult ref, job;
            both_loops(ix.bwt, to_codes(reads[i]), opt, ref, job);
            REQUIRE(ref.aln.size() == job.aln.size());
            for (size_t k = 0; k < ref.aln.size(); ++k) CHECK(same_aln(ref.aln[k], job.aln[k]));
            ++compared;
            with_hits += !ref.aln.empty();
            bool g = false;
            for (const auto &a : ref.aln) g = g || a.n_gapo;
            gapped += g;
        }
    }
    MESSAGE("compared " << compared << ", with hits " << with_hits << ", with a gapped hit " << gapped);
    CHECK(with_hits > compared / 2);
    CHECK(gapped > compared / 10);
}

// gap_cstack_t (the job's stack) against gap_stack_t (upstream's): random push/pop sequences
// must pop the same entries in the same order, with the same counts and best bucket after
// every step; and the pool must hold exactly the chunks the buckets need.
TEST_CASE("aln_search: gap_cstack_t pops in gap_stack_t's order"
          * doctest::test_suite("unit/aln_search")) {
    ba::gap_opt_t opt = make_opt(4, 2);
    opt.max_gape = 3;
    const int C = ba::GAP_CHUNK_ENTRIES;
    REQUIRE(C > 8);
    uint32_t x = 99;
    auto rnd = [&x](uint32_t n) { x = x * 1103515245u + 12345u; return (x >> 8) % n; };
    for (int keep : {0, 1, ba::ALN_POOL_KEEP}) {
        ba::gap_pool_t pool;
        ba::gap_init_pool(&pool, keep);
        ba::gap_stack_t *a = ba::gap_init_stack(opt.max_diff, opt.max_gapo, opt.max_gape, &opt);
        ba::gap_cstack_t *b = ba::gap_init_cstack(opt.max_diff, opt.max_gapo, opt.max_gape, &opt, &pool);
        REQUIRE(a->n_stacks == b->n_stacks);
        ba::gap_stack_reset_ext(a); // as bwt_match_gap does before its first push
        long n_pops = 0, seq_id = 0, peak_in_use = 0;
        // compared with plain C++ and asserted once per call: a doctest assertion per bucket per
        // step would cost ~0.5 s
        auto check_state = [&]() {
            const char *bad = nullptr;
            int bad_bucket = -1;
            long chunks = 0;
            if (a->n_entries != b->n_entries) bad = "n_entries";
            else if (a->best != b->best) bad = "best";
            for (int s = 0; s < b->n_stacks && !bad; ++s) {
                if (a->stacks[s].n_entries != b->stacks[s].n_entries) { bad = "bucket n_entries"; bad_bucket = s; }
                chunks += (b->stacks[s].n_entries + C - 1) / C;
            }
            if (!bad && pool.n_alloc - pool.n_free != chunks) bad = "in-use chunks != entries' chunks";
            if (!bad && pool.n_free > keep) bad = "pool keeps more than max_free";
            // free chunks are reused before new ones are allocated: the pool never holds more
            // chunks than were in use at once
            if (pool.n_alloc - pool.n_free > peak_in_use) peak_in_use = pool.n_alloc - pool.n_free;
            if (!bad && pool.n_alloc > peak_in_use) bad = "pool allocated a chunk while holding a free one";
            INFO("sequence " << seq_id << ", pop " << n_pops << ": " << (bad ? bad : "") << " (bucket " << bad_bucket << ")");
            REQUIRE(bad == nullptr);
        };
        // one push into the bucket of (mm, go, ge), with the other fields drawn at random
        auto push = [&](int mm, int go, int ge) {
            int i = (int)rnd(200), state = (int)rnd(3), is_diff = (int)rnd(2);
            ba::bwtint_t k = rnd(100000), l = k + rnd(1000);
            int ins = (int)rnd(5), del = (int)rnd(5);
            ba::gap_stack_push_ext(a, i, k, l, mm, go, ge, ins, del, state, is_diff, &opt);
            ba::gap_cstack_push_ext(b, i, k, l, mm, go, ge, ins, del, state, is_diff, &opt);
        };
        auto pop = [&]() {
            ba::gap_entry_t ea, eb;
            ba::gap_stack_pop_ext(a, &ea);
            ba::gap_cstack_pop_ext(b, &eb);
            const bool same = ea.info == eb.info && ea.n_mm == eb.n_mm && ea.n_gapo == eb.n_gapo
                && ea.n_gape == eb.n_gape && ea.state == eb.state && ea.n_ins == eb.n_ins
                && ea.n_del == eb.n_del && ea.last_diff_pos == eb.last_diff_pos && ea.k == eb.k && ea.l == eb.l;
            INFO("sequence " << seq_id << ", pop " << n_pops << ": info " << ea.info << "/" << eb.info
                 << " k " << ea.k << "/" << eb.k << " l " << ea.l << "/" << eb.l);
            REQUIRE(same);
            ++n_pops;
        };
        auto rand_push = [&]() {
            push((int)rnd(opt.max_diff + 1), (int)rnd(opt.max_gapo + 1), (int)rnd(opt.max_gape + 1));
        };
        // (1) one bucket filled to 1, C-1, C, C+1, 2C, 2C+1 entries, drained, refilled
        for (int n : {1, C - 1, C, C + 1, 2 * C, 2 * C + 1}) {
            ++seq_id;
            for (int r = 0; r < 2; ++r) {
                for (int t = 0; t < n; ++t) { push(1, 0, 0); check_state(); }
                for (int t = 0; t < n; ++t) { pop(); check_state(); }
            }
            // oscillate across the chunk boundary
            for (int t = 0; t < n; ++t) push(1, 0, 0);
            for (int r = 0; r < 5; ++r) { push(1, 0, 0); check_state(); pop(); check_state(); pop(); check_state(); push(1, 0, 0); }
            while (a->n_entries) { pop(); check_state(); }
        }
        // (2) random sequences: bursts of pushes into random buckets, then runs of pops, so
        // the best bucket moves up and down and buckets drain and refill
        for (int sq = 0; sq < 30; ++sq) {
            ++seq_id;
            int steps = 50 + (int)rnd(400);
            for (int st = 0; st < steps; ++st) {
                if (a->n_entries == 0 || rnd(100) < 55) {
                    int burst = rnd(10) == 0 ? (int)rnd(3 * C) : 1 + (int)rnd(8);
                    if (rnd(4) == 0) { // all in one bucket, to cross chunk boundaries
                        int mm = (int)rnd(opt.max_diff + 1), go = (int)rnd(opt.max_gapo + 1), ge = (int)rnd(opt.max_gape + 1);
                        for (int t = 0; t < burst; ++t) push(mm, go, ge);
                    } else {
                        for (int t = 0; t < burst; ++t) rand_push();
                    }
                } else {
                    int run = rnd(10) == 0 ? (int)rnd(3 * C) : 1 + (int)rnd(8);
                    for (int t = 0; t < run && a->n_entries; ++t) pop();
                }
                check_state();
            }
            if (sq % 2) { // empty both the way the searches do at the start of a read
                ba::gap_stack_reset_ext(a);
                ba::gap_release_cstack(b);
                check_state();
                REQUIRE(pool.n_alloc == pool.n_free);
            } else {
                while (a->n_entries) pop();
                check_state();
            }
        }
        // (3) two stacks sharing the pool, as a thread's slots do: interleaved pushes and pops on
        // both, each checked against its own gap_stack_t, and the pool's in-use chunks equal to
        // the two stacks' chunks together
        {
            ++seq_id;
            ba::gap_stack_t *a2 = ba::gap_init_stack(opt.max_diff, opt.max_gapo, opt.max_gape, &opt);
            ba::gap_cstack_t *b2 = ba::gap_init_cstack(opt.max_diff, opt.max_gapo, opt.max_gape, &opt, &pool);
            ba::gap_stack_reset_ext(a2);
            auto chunks_of = [&](const ba::gap_cstack_t *s) {
                long n = 0;
                for (int q = 0; q < s->n_stacks; ++q) n += (s->stacks[q].n_entries + C - 1) / C;
                return n;
            };
            bool ok = true;
            for (int st = 0; st < 3000 && ok; ++st) {
                const bool second = rnd(2) != 0;
                ba::gap_stack_t *ra = second ? a2 : a;
                ba::gap_cstack_t *rb = second ? b2 : b;
                if (ra->n_entries == 0 || rnd(100) < 52) {
                    const int burst = rnd(20) == 0 ? (int)rnd(2 * C) : 1 + (int)rnd(6);
                    const int mm = (int)rnd(opt.max_diff + 1), go = (int)rnd(opt.max_gapo + 1), ge = (int)rnd(opt.max_gape + 1);
                    for (int t = 0; t < burst; ++t) {
                        int i = (int)rnd(200), state = (int)rnd(3), is_diff = (int)rnd(2);
                        ba::bwtint_t k = rnd(100000), l = k + rnd(1000);
                        ba::gap_stack_push_ext(ra, i, k, l, mm, go, ge, 0, 0, state, is_diff, &opt);
                        ba::gap_cstack_push_ext(rb, i, k, l, mm, go, ge, 0, 0, state, is_diff, &opt);
                    }
                } else {
                    const int run = rnd(20) == 0 ? (int)rnd(2 * C) : 1 + (int)rnd(6);
                    for (int t = 0; t < run && ra->n_entries && ok; ++t) {
                        ba::gap_entry_t ea, eb;
                        ba::gap_stack_pop_ext(ra, &ea);
                        ba::gap_cstack_pop_ext(rb, &eb);
                        ok = ea.info == eb.info && ea.k == eb.k && ea.l == eb.l && ea.last_diff_pos == eb.last_diff_pos
                            && ea.n_mm == eb.n_mm && ea.n_gapo == eb.n_gapo && ea.n_gape == eb.n_gape && ea.state == eb.state;
                    }
                }
                ok = ok && a->n_entries == b->n_entries && a2->n_entries == b2->n_entries && a->best == b->best
                    && a2->best == b2->best && pool.n_alloc - pool.n_free == chunks_of(b) + chunks_of(b2)
                    && pool.n_free <= keep;
            }
            INFO("keep " << keep << ": two stacks on one pool");
            REQUIRE(ok);
            ba::gap_destroy_stack(a2);
            ba::gap_destroy_cstack(b2);
            ba::gap_stack_reset_ext(a);
            ba::gap_release_cstack(b);
            REQUIRE(pool.n_alloc == pool.n_free);
        }
        CHECK(n_pops > 15000); // a floor for this fixed seed
        ba::gap_destroy_stack(a);
        ba::gap_destroy_cstack(b);
        CHECK(pool.n_alloc == pool.n_free);
        CHECK(pool.n_free <= keep);
        ba::gap_destroy_pool(&pool);
        CHECK(pool.n_alloc == 0);
    }
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
        {"-n 0", 0, 1, 2000000, -1.0}, // the root entry runs as an exact tail
        // more differences, so the slots' stacks hold more chunks and draw on the thread's pool;
        // the pool cap (ALN_POOL_KEEP) is not reached here: the single-stack tests (pool caps 0
        // and 1) and "two stacks sharing the pool" cover the path where the pool frees chunks
        {"-n 5 -o 2", 5, 2, 2000000, -1.0},
        {"-n 5 -o 2 -m chunk+1", 5, 2, ba::GAP_CHUNK_ENTRIES + 1, -1.0},
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
    CHECK(compared == 2L * (long)(sizeof cases / sizeof cases[0]) * (long)reads.size());
    CHECK(with_hits > compared / 3); // not vacuous
}

// gap_slots_run's slot states: batches smaller than, equal to and just over the group size
// (ALN_GROUP: 8, or 12 on Apple arm64) and run on more threads than reads, so some threads get no
// read and some slots never fill or go empty early; all-N and one-base reads end in gap_job_start,
// before any search step, so their slot is refilled straight from GAP_SLOT_WIDTHS.
TEST_CASE("aln_search: aln_search_batch gives bwt_match_gap's hits for batches around the group size"
          * doctest::test_suite("unit/aln_search")) {
    SearchIndex &ix = index_fixture();
    auto contigs = aln_test::test_contigs();
    const std::string &s1 = contigs[0].second;
    uint32_t x = 4242;
    auto rnd = [&x](uint32_t n) { x = x * 1103515245u + 12345u; return (x >> 16) % n; };
    std::vector<std::string> pool;
    for (int i = 0; i < 30; ++i) {
        if (i % 5 == 1) { pool.push_back(std::string(20 + i, 'N')); continue; }
        if (i % 5 == 3) { pool.push_back("ACGT"[i % 4] + std::string()); continue; }
        int len = 25 + (int)rnd(30);
        std::string r = s1.substr(rnd((uint32_t)(s1.size() - len)), len);
        if (i % 2) r[rnd(len)] = r[0] == 'A' ? 'C' : 'A';
        pool.push_back(r);
    }
    ba::gap_opt_t opt = make_opt(2, 1);
    long compared = 0, with_hits = 0;
    for (int n : {1, 2, 7, 8, 9, 11, 12, 13, 25}) {
        std::vector<std::string> reads(pool.begin(), pool.begin() + n);
        for (int nt : {1, 3, 16}) {
            Batch b(reads);
            ba::aln_search_batch(ix.bwt, b.n, b.seqs, &opt, nt);
            const int cap = ba::aln_batch_max_gapo(b.n, b.seqs, &opt);
            for (int i = 0; i < b.n; ++i) {
                INFO("batch " << n << ", -t " << nt << ", read " << i << " (" << reads[i] << ")");
                GapResult ref, job;
                both_loops(ix.bwt, to_codes(reads[i]), opt, ref, job, cap);
                REQUIRE(b.seqs[i].n_aln == (int)ref.aln.size());
                for (int k = 0; k < b.seqs[i].n_aln; ++k) CHECK(same_aln(ref.aln[k], b.seqs[i].aln[k]));
                ++compared;
                with_hits += !ref.aln.empty();
            }
        }
    }
    CHECK(compared == 3L * (1 + 2 + 7 + 8 + 9 + 11 + 12 + 13 + 25));
    CHECK(with_hits > compared / 3);
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
    reads.push_back(s1.substr(150, 400));                 // long (c1 is 600 bp)
    reads.push_back(flip(flip(s1.substr(150, 400), 100), 300));
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
        int yields = 0;
        while (!ba::width_walk_run(&ww, ix.bwt))
            if (++yields >= 10000000) break;
        // one yield before each Occ lookup, i.e. per non-N base, so the walk interleaves
        int lookups = 0;
        for (uint8_t c : codes) lookups += c < 4;
        CHECK(yields == lookups);
        CHECK(ww.bid == bid);
        for (int i = 0; i <= len; ++i) {
            CHECK(a[i].w == b[i].w);
            CHECK(a[i].bid == b[i].bid);
        }
        restarts += bid > 1;
    }
    CHECK(restarts >= 6); // restarts are exercised
}
