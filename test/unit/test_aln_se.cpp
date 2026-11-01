// Unit tests for the single-end half of `bwa-mem3 aln` (src/aln/bwase.cpp, aln_out.cpp,
// aln_main.cpp). Byte-for-byte agreement with stock bwa aln + samse is
// test/regression/aln_parity.sh; these cases pin the header layout and the once-per-run
// drand48 seeding on the shared synthetic reference (the BAM sink is test_aln_se_bam.cpp).

#include "doctest/doctest.h"

#include "aln_test_ref.h"
#include "../../src/FMI_search.h"
#include "../../src/aln/aln_main.h"
#include "../../src/aln/bwase.h"
#include "../../src/aln/bwaseqio.h"
#include "../../src/bwa.h"
#include "../../src/read_index_ele.h"

#include <climits>
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
        std::string tmpl_s = aln_test::tmp_template("aln_se_unit");
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
    ba::aln_out_t *out = ba::aln_out_open(out_path.c_str(), 0, 0, "");
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

TEST_CASE("aln_se: a huge -l allocates no seed widths for shorter reads and matches -l 1024"
          * doctest::test_suite("unit/aln_se")) {
    // The seed-width buffer is sized only for a read longer than -l. Sized per search slot up
    // front, -l 2^28 would zero-fill 8 x 4 GiB per thread here; with INT_MAX, -l + 1 overflowed.
    SeIndex &ix = index_fixture();
    const std::string c1 = aln_test::test_contigs()[0].second;
    std::string fq = write_fastq(ix.dir, "seedw.fq", {c1.substr(40, 36), c1.substr(300, 36), c1.substr(77, 30)});
    auto run = [&](int seed_len, const std::string &out_path) {
        ba::aln_opt_t o;
        ba::aln_opt_init(&o);
        o.gap.seed_len = seed_len;
        ba::aln_out_t *out = ba::aln_out_open(out_path.c_str(), 0, 0, "");
        REQUIRE(out != nullptr);
        REQUIRE(ba::aln_run_se(&o, ix.bwt, ix.bns(), ix.fmi->idx->pac, fq.c_str(), out) == 0);
        REQUIRE(ba::aln_out_close(out) == 0);
        std::ifstream in(out_path);
        std::stringstream ss;
        ss << in.rdbuf();
        return ss.str();
    };
    const std::string base = run(1024, ix.dir + "/seedw_1024.sam");
    CHECK(base.find("\tXT:A:U") != std::string::npos);
    CHECK(run(1 << 28, ix.dir + "/seedw_big.sam") == base);
    CHECK(run(INT_MAX, ix.dir + "/seedw_max.sam") == base);
    // and a read longer than -l still gets its seed widths
    CHECK(run(20, ix.dir + "/seedw_20.sam").find("\tXT:A:U") != std::string::npos);
}

TEST_CASE("aln_pe: drand48 is seeded once per run, not per batch, in aln_run_pe too"
          * doctest::test_suite("unit/aln_se")) {
    // Both ends of each pair lie in c1[100:400], which c2 holds reverse-complemented, so each
    // end has two equally good hits and draws from drand48. Eight copies of the pair, in batches
    // of 4 and of 2: the draws stay in one serial pass seeded once per run, so the records match.
    SeIndex &ix = index_fixture();
    const std::string c1 = aln_test::test_contigs()[0].second;
    const std::string m1 = c1.substr(150, 36), m2 = aln_test::revcomp(c1.substr(150 + 200 - 36, 36));
    std::vector<std::string> v1(8, m1), v2(8, m2);
    const std::string f1 = write_fastq(ix.dir, "rep_pe1.fq", v1), f2 = write_fastq(ix.dir, "rep_pe2.fq", v2);
    auto run = [&](int bs, const std::string &path) {
        ba::aln_opt_t o;
        ba::aln_opt_init(&o);
        o.batch_size = bs;
        ba::aln_out_t *out = ba::aln_out_open(path.c_str(), 0, 0, "");
        REQUIRE(out != nullptr);
        REQUIRE(ba::aln_run_pe(&o, ix.bwt, ix.bns(), ix.fmi->idx->pac, f1.c_str(), f2.c_str(), out) == 0);
        REQUIRE(ba::aln_out_close(out) == 0);
        std::ifstream in(path);
        std::stringstream ss;
        ss << in.rdbuf();
        return ss.str();
    };
    const std::string four = run(4, ix.dir + "/rep_pe_4.sam");
    CHECK(four.find("X0:i:2") != std::string::npos);
    CHECK(run(2, ix.dir + "/rep_pe_2.sam") == four);
    CHECK(run(4, ix.dir + "/rep_pe_again.sam") == four);
}

// --- SAM formatting without printf -------------------------------------------------------
// bwa_print_sam1 formats with kputs/kputc/kputl so the formatting threads do not contend on
// libc's locale lock. These cases hold that text equal to what ksprintf produced.

namespace {

// The integer helpers behind bwa_print_sam1 are kputl((long)v); check them against "%d" /
// "%lld" at the edges, so a change to either side shows up.
std::string via_kputl(long long v)
{
    kstring_t ks = {0, 0, 0};
    kputl((long)v, &ks);
    std::string s(ks.s, ks.l);
    free(ks.s);
    return s;
}

std::string via_ksprintf_ll(long long v)
{
    kstring_t ks = {0, 0, 0};
    ksprintf(&ks, "%lld", v);
    std::string s(ks.s, ks.l);
    free(ks.s);
    return s;
}

// The pre-change bwa_print_sam1, verbatim apart from the helpers it needs (ksprintf everywhere).
struct RefSam {
    static int64_t pos_5(const ba::bwa_seq_t *p)
    {
        if (p->type != BWA_TYPE_NO_MATCH) return p->strand ? ba::pos_end(p) : p->pos;
        return -1;
    }
    static int64_t pos_end_multi(const ba::bwt_multi1_t *p, int len)
    {
        if (p->cigar) {
            int64_t x = p->pos;
            for (int j = 0; j != (int)p->n_cigar; ++j) {
                int op = __cigar_op(p->cigar[j]);
                if (op == 0 || op == 2) x += __cigar_len(p->cigar[j]);
            }
            return x;
        }
        return p->pos + len;
    }
    static void print_seq(kstring_t *out, const ba::bwa_seq_t *seq)
    {
        if (seq->strand == 0) {
            for (int i = 0; i < (int)seq->full_len; ++i) kputc("ACGTN"[seq->seq[i]], out);
        } else {
            for (int i = (int)seq->full_len - 1; i >= 0; --i) kputc("TGCAN"[seq->seq[i]], out);
        }
    }
    static void print(kstring_t *out, const bntseq_t *bns, ba::bwa_seq_t *p, const ba::bwa_seq_t *mate, int mode,
                      int max_top2)
    {
        int j;
        if (p->type != BWA_TYPE_NO_MATCH || (mate && mate->type != BWA_TYPE_NO_MATCH)) {
            int seqid, nn, am = 0, flag = p->extra_flag;
            char XT;
            if (p->type == BWA_TYPE_NO_MATCH) {
                p->pos = mate->pos;
                p->strand = mate->strand;
                flag |= SAM_FSU;
                j = 1;
            } else j = ba::pos_end(p) - p->pos;
            nn = bns_cnt_ambi(bns, p->pos, j, &seqid);
            if (p->type != BWA_TYPE_NO_MATCH &&
                (uint64_t)p->pos + (uint64_t)j - (uint64_t)bns->anns[seqid].offset > (uint64_t)bns->anns[seqid].len)
                flag |= SAM_FSU;
            if (p->strand) flag |= SAM_FSR;
            if (mate) {
                if (mate->type != BWA_TYPE_NO_MATCH) {
                    if (mate->strand) flag |= SAM_FMR;
                } else flag |= SAM_FMU;
            }
            ksprintf(out, "%s\t%d\t%s\t", p->name, flag, bns->anns[seqid].name);
            ksprintf(out, "%d\t%d\t", (int)(p->pos - bns->anns[seqid].offset + 1), p->mapQ);
            if (p->cigar) {
                for (j = 0; j != p->n_cigar; ++j)
                    ksprintf(out, "%d%c", __cigar_len(p->cigar[j]), "MIDS"[__cigar_op(p->cigar[j])]);
            } else if (p->type == BWA_TYPE_NO_MATCH) ksprintf(out, "*");
            else ksprintf(out, "%dM", p->len);
            if (mate && mate->type != BWA_TYPE_NO_MATCH) {
                int m_seqid;
                long long isize;
                am = mate->seQ < p->seQ ? mate->seQ : p->seQ;
                bns_cnt_ambi(bns, mate->pos, mate->len, &m_seqid);
                ksprintf(out, "\t%s\t", (seqid == m_seqid) ? "=" : bns->anns[m_seqid].name);
                isize = (seqid == m_seqid) ? pos_5(mate) - pos_5(p) : 0;
                if (p->type == BWA_TYPE_NO_MATCH) isize = 0;
                ksprintf(out, "%d\t%lld\t", (int)(mate->pos - bns->anns[m_seqid].offset + 1), isize);
            } else if (mate) ksprintf(out, "\t=\t%d\t0\t", (int)(p->pos - bns->anns[seqid].offset + 1));
            else ksprintf(out, "\t*\t0\t0\t");
            print_seq(out, p);
            kputc('\t', out);
            if (p->qual) {
                if (p->strand) ba::seq_reverse(p->len, p->qual, 0);
                ksprintf(out, "%s", p->qual);
            } else ksprintf(out, "*");
            if (bwa_rg_id[0]) ksprintf(out, "\tRG:Z:%s", bwa_rg_id);
            if (p->bc[0]) ksprintf(out, "\tBC:Z:%s", p->bc);
            if (p->clip_len < (int)p->full_len) ksprintf(out, "\tXC:i:%d", p->clip_len);
            if (p->type != BWA_TYPE_NO_MATCH) {
                XT = "NURM"[p->type];
                if (nn > 10) XT = 'N';
                ksprintf(out, "\tXT:A:%c\t%s:i:%d", XT, (mode & BWA_MODE_COMPREAD) ? "NM" : "CM", p->nm);
                if (nn) ksprintf(out, "\tXN:i:%d", nn);
                if (mate) ksprintf(out, "\tSM:i:%d\tAM:i:%d", (int)p->seQ, am);
                if (p->type != BWA_TYPE_MATESW) {
                    ksprintf(out, "\tX0:i:%d", (int)p->c1);
                    if ((int)p->c1 <= max_top2) ksprintf(out, "\tX1:i:%d", (int)p->c2);
                }
                ksprintf(out, "\tXM:i:%d\tXO:i:%d\tXG:i:%d", p->n_mm, p->n_gapo, p->n_gapo + p->n_gape);
                if (p->md) ksprintf(out, "\tMD:Z:%s", p->md);
                if (p->n_multi) {
                    ksprintf(out, "\tXA:Z:");
                    for (int i = 0; i < p->n_multi; ++i) {
                        ba::bwt_multi1_t *q = p->multi + i;
                        j = pos_end_multi(q, p->len) - q->pos;
                        nn = bns_cnt_ambi(bns, q->pos, j, &seqid);
                        ksprintf(out, "%s,%c%d,", bns->anns[seqid].name, q->strand ? '-' : '+',
                                 (int)(q->pos - bns->anns[seqid].offset + 1));
                        if (q->cigar) {
                            for (int k = 0; k < (int)q->n_cigar; ++k)
                                ksprintf(out, "%d%c", __cigar_len(q->cigar[k]), "MIDS"[__cigar_op(q->cigar[k])]);
                        } else ksprintf(out, "%dM", p->len);
                        ksprintf(out, ",%d;", q->gap + q->mm);
                    }
                }
            }
            kputc('\n', out);
        } else {
            int flag = p->extra_flag | SAM_FSU;
            if (mate && mate->type == BWA_TYPE_NO_MATCH) flag |= SAM_FMU;
            ksprintf(out, "%s\t%d\t*\t0\t0\t*\t*\t0\t0\t", p->name, flag);
            print_seq(out, p);
            kputc('\t', out);
            if (p->qual) {
                if (p->strand) ba::seq_reverse(p->len, p->qual, 0);
                ksprintf(out, "%s", p->qual);
            } else ksprintf(out, "*");
            if (bwa_rg_id[0]) ksprintf(out, "\tRG:Z:%s", bwa_rg_id);
            if (p->bc[0]) ksprintf(out, "\tBC:Z:%s", p->bc);
            if (p->clip_len < (int)p->full_len) ksprintf(out, "\tXC:i:%d", p->clip_len);
            kputc('\n', out);
        }
    }
};

struct ReadBuf {
    std::vector<ubyte_t> seq, qual;
    std::string name;
    std::vector<ba::bwa_cigar_t> cigar, mcigar;
    std::vector<ba::bwt_multi1_t> multi;
    std::string md;
    ba::bwa_seq_t s;
    ReadBuf(int len, uint32_t type, bwtint_t pos, int strand)
    {
        std::memset(&s, 0, sizeof s);
        for (int i = 0; i < len; ++i) { seq.push_back((ubyte_t)(i % 5)); qual.push_back((ubyte_t)('!' + i % 40)); }
        qual.push_back(0);
        name = "read/1";
        s.name = &name[0];
        s.seq = seq.data();
        s.qual = qual.data();
        s.len = s.full_len = len;
        s.clip_len = len;
        s.type = type;
        s.pos = pos;
        s.strand = strand;
    }
    void finish() { s.cigar = cigar.empty() ? nullptr : cigar.data(); s.n_cigar = (int)cigar.size();
                    s.multi = multi.empty() ? nullptr : multi.data(); s.n_multi = (int)multi.size();
                    s.md = md.empty() ? nullptr : &md[0]; }
};

// Prints (p, mate) with the shipped printer and the ksprintf reference on separate copies
// (the printers reverse .qual in place on the reverse strand) and compares the text.
void check_same(ReadBuf &a, ReadBuf *mate, int mode, int max_top2)
{
    SeIndex &ix = index_fixture();
    a.finish();
    if (mate) mate->finish();
    ReadBuf b = a;            // vectors are copied; repoint the raw pointers into the copy
    b.s.name = &b.name[0]; b.s.seq = b.seq.data(); b.s.qual = b.qual.data();
    b.s.cigar = b.cigar.empty() ? nullptr : b.cigar.data();
    b.s.multi = b.multi.empty() ? nullptr : b.multi.data();
    for (size_t i = 0; i < b.multi.size(); ++i) b.multi[i].cigar = a.multi[i].cigar ? b.mcigar.data() : nullptr;
    b.s.md = b.md.empty() ? nullptr : &b.md[0];
    kstring_t k1 = {0, 0, 0}, k2 = {0, 0, 0};
    ba::bwa_print_sam1(&k1, ix.bns(), &a.s, mate ? &mate->s : nullptr, mode, max_top2);
    RefSam::print(&k2, ix.bns(), &b.s, mate ? &mate->s : nullptr, mode, max_top2);
    std::string s1(k1.s, k1.l), s2(k2.s, k2.l);
    free(k1.s); free(k2.s);
    CHECK(s1 == s2);
    CHECK(!s1.empty());
}

} // namespace

TEST_CASE("aln_se: kputl formats integers like ksprintf at the edges"
          * doctest::test_suite("unit/aln_se")) {
    static_assert(sizeof(long) == sizeof(long long), "kputl must hold 64 bits");
    const long long vals[] = {0, 1, -1, 9, 10, -10, 99, 100, 2147483647LL, -2147483647LL - 1, 2147483648LL,
                              4294967295LL, 4294967296LL, 123456789012345LL, -123456789012345LL,
                              9223372036854775807LL, -9223372036854775807LL};
    for (long long v : vals) CHECK(via_kputl(v) == via_ksprintf_ll(v));
    // The int paths print through (long)int: INT_MIN must not wrap (kputw's negation does).
    kstring_t ks = {0, 0, 0};
    ksprintf(&ks, "%d", INT_MIN);
    CHECK(via_kputl((long)INT_MIN) == std::string(ks.s, ks.l));
    free(ks.s);
}

TEST_CASE("aln_se: bwa_print_sam1 matches the ksprintf printer"
          * doctest::test_suite("unit/aln_se")) {
    const int mode = 0;
    std::strcpy(bwa_rg_id, "");
    {   // unique forward hit on c1, ungapped, with MD, a barcode and a clipped read
        ReadBuf r(36, BWA_TYPE_UNIQUE, 100, 0);
        r.s.mapQ = 37; r.s.seQ = 37; r.s.c1 = 1; r.s.c2 = 0; r.s.nm = 2; r.s.n_mm = 2;
        r.s.clip_len = 30; r.md = "10A5C19";
        std::strcpy(r.s.bc, "ACGTAC");
        check_same(r, nullptr, mode, 1);
    }
    {   // reverse hit on c2 (offset 600), gapped CIGAR, read group, repeat with XA hits
        std::strcpy(bwa_rg_id, "grp1");
        ReadBuf r(40, BWA_TYPE_REPEAT, 650, 1);
        r.cigar = {(ba::bwa_cigar_t)__cigar_create(0, 10), (ba::bwa_cigar_t)__cigar_create(1, 2),
                   (ba::bwa_cigar_t)__cigar_create(0, 20), (ba::bwa_cigar_t)__cigar_create(2, 3),
                   (ba::bwa_cigar_t)__cigar_create(3, 8)};
        r.s.mapQ = 0; r.s.seQ = 0; r.s.c1 = 5; r.s.c2 = 3; r.s.nm = 4095; r.s.n_mm = 255; r.s.n_gapo = 255;
        r.s.n_gape = 255; r.md = "10^AC20";
        r.mcigar = {(ba::bwa_cigar_t)__cigar_create(0, 40)};
        ba::bwt_multi1_t m0; std::memset(&m0, 0, sizeof m0);
        m0.pos = 20; m0.strand = 1; m0.gap = 1; m0.mm = 2; m0.n_cigar = 1;
        ba::bwt_multi1_t m1; std::memset(&m1, 0, sizeof m1);
        m1.pos = 610; m1.strand = 0; m1.gap = 0; m1.mm = 255;
        r.multi = {m0, m1};
        r.multi[0].cigar = r.mcigar.data();
        r.multi[1].cigar = nullptr;
        check_same(r, nullptr, mode | BWA_MODE_COMPREAD, 100);
        std::strcpy(bwa_rg_id, "");
    }
    {   // a pair on one contig (negative and positive insert size), and across contigs
        ReadBuf r1(30, BWA_TYPE_UNIQUE, 200, 0), r2(30, BWA_TYPE_UNIQUE, 150, 1);
        r1.s.mapQ = 60; r1.s.seQ = 40; r1.s.c1 = 1; r1.s.extra_flag = SAM_FPD | SAM_FR1;
        r2.s.mapQ = 55; r2.s.seQ = 30; r2.s.c1 = 1; r2.s.extra_flag = SAM_FPD | SAM_FR2;
        check_same(r1, &r2, mode, 1);
        check_same(r2, &r1, mode, 1);
        ReadBuf r3(30, BWA_TYPE_MATESW, 700, 0);
        r3.s.mapQ = 20; r3.s.seQ = 10; r3.s.extra_flag = SAM_FPD | SAM_FR2;
        check_same(r1, &r3, mode, 1);
        check_same(r3, &r1, mode, 1);
    }
    {   // unmapped, alone and next to a mapped mate and an unmapped mate
        ReadBuf u(25, BWA_TYPE_NO_MATCH, 0, 0), v(25, BWA_TYPE_NO_MATCH, 0, 1), m(25, BWA_TYPE_UNIQUE, 120, 1);
        m.s.mapQ = 25; m.s.seQ = 25; m.s.c1 = 1;
        u.s.clip_len = 20;
        check_same(u, nullptr, mode, 1);
        check_same(u, &v, mode, 1);
        check_same(u, &m, mode, 1);
        u.qual.clear(); u.s.qual = nullptr;
        check_same(u, nullptr, mode, 1);
    }
}
