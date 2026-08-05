// test/unit/test_bam_pa_tag_parity.cpp — the `pa` tag must mean the same thing
// in SAM text and in BAM (fg-labs/bwa-mem3#365).
//
// `pa` is the ratio of a hit's score to the score of the better overlapping ALT
// hit. bwa-mem3 emits it from three writers — the SAM-text path
// (mem_aln2sam), the --bam path (mem_aln_to_bam) and the --meth path
// (meth_mem_aln_to_bam) — and `--bam` is documented as a container choice, not
// a content change. So a `--bam` record must carry exactly the value a
// `samtools view -b` of the same run's SAM text would store: float32 of the
// three-decimal SAM token.
//
// The BAM writers used to store the raw unrounded quotient instead, so
// `pa:f:0.806723` came out of --bam where the SAM path (and both upstreams,
// which are SAM-text-only and both round) emit `pa:f:0.807`. These tests are
// against the SAM token rather than against any particular rounding
// expression, because the two only agree if the BAM value is derived from the
// same decimal rendering: `round(x * 1000) / 1000` is NOT equivalent (it
// rounds halves away from zero where printf rounds them to even, so 39/48
// gives 0.813 instead of the SAM path's 0.812).
//
// htslib headers first: they share the KSTRING_H include guard with
// bwa-mem3's kstring.h, and the BAM writers are compiled against htslib's.
// Matching that order here keeps `kstring_t` a single type.
#include "htslib/sam.h"
#include "htslib/kstring.h"

#include "doctest/doctest.h"

#include "bam_writer.h"
#include "bwa.h"
#include "bwamem.h"

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

// bwamem.h declares mem_aln2sam against bwa-mem3's own kstring_t, but htslib's
// kstring.h won the shared KSTRING_H guard above, so the declaration bwamem.h
// produces in this translation unit mangles to a symbol libbwa.a does not
// define. Re-declare it against bwa-mem3's struct tag (same layout, different
// name) so the call resolves to the definition in bwamem.cpp.
struct __kstring_t { size_t l, m; char *s; };
void mem_aln2sam(const mem_opt_t *opt, const bntseq_t *bns, __kstring_t *str,
                 bseq1_t *s, int n, const mem_aln_t *list, int which,
                 const mem_aln_t *m_);

namespace {

// A one-contig reference. The writers only read anns[rid].{name,anno,is_alt}
// and n_seqs/l_pac for a mapped record, so a stack fixture is enough; `name`
// and `anno` are not owned by bntann1_t in bwa, so string literals are correct.
struct RefFixture {
    bntann1_t ann{};
    bntseq_t  bns{};

    RefFixture() {
        ann.offset = 0;
        ann.len    = 1000;
        ann.n_ambs = 0;
        ann.gi     = 0;
        ann.is_alt = 0;
        ann.name   = const_cast<char *>("chr1");
        ann.anno   = const_cast<char *>("");
        bns.l_pac  = 1000;
        bns.n_seqs = 1;
        bns.anns   = &ann;
    }
};

// One mapped read. `seq` is 2-bit (nt4) encoded, which is the form both
// writers expect by the time they run; `qual` is ASCII phred+33.
struct ReadFixture {
    static const int kLen = 10;

    std::string        name{"read1"};
    std::vector<char>  seq;
    std::string        qual;
    bseq1_t            s{};

    ReadFixture() : seq(kLen), qual(kLen, 'I') {
        for (int i = 0; i < kLen; ++i) seq[i] = (char)(i & 3);
        s.l_seq = kLen;
        s.id    = 0;
        s.name  = const_cast<char *>(name.c_str());
        s.seq   = seq.data();
        s.qual  = const_cast<char *>(qual.c_str());
    }
};

// A CIGAR buffer in mem_reg2aln's layout: `n_cigar` ops followed by the MD:Z
// string, which both writers read as `(char *)(cigar + n_cigar)`.
struct CigarFixture {
    uint32_t buf[4]{};

    CigarFixture() {
        buf[0] = ((uint32_t)ReadFixture::kLen << 4) | 0;  // 10M
        memcpy(&buf[1], "10", 3);                         // MD:Z:10
    }
};

// A primary, forward-strand alignment scoring `score` against an ALT hit
// scoring `alt_sc`.
mem_aln_t make_aln(uint32_t *cigar, int score, int alt_sc, int flag = 0) {
    mem_aln_t a{};
    a.pos     = 100;
    a.rid     = 0;
    a.flag    = flag;
    a.is_rev  = 0;
    a.is_alt  = 0;
    a.mapq    = 60;
    a.NM      = 0;
    a.n_cigar = 1;
    a.cigar   = cigar;
    a.XA      = NULL;
    a.HN      = -1;
    a.score   = score;
    a.sub     = -1;
    a.alt_sc  = alt_sc;
    a.meth_hypothesis = -1;
    return a;
}

// The `pa` token of a SAM record, or "" when the record carries no pa tag.
// mem_aln2sam terminates the record with '\n', so the token ends at whichever
// of the tag separator or the line terminator comes first.
std::string sam_pa(const std::string &record) {
    const size_t at = record.find("\tpa:f:");
    if (at == std::string::npos) return "";
    const size_t begin = at + 6;
    const size_t end   = record.find_first_of("\t\n", begin);
    return record.substr(begin, end == std::string::npos ? std::string::npos : end - begin);
}

// Render one alignment through the SAM-text writer.
std::string emit_sam(const mem_opt_t *opt, const bntseq_t *bns, bseq1_t *s,
                     const mem_aln_t *aln) {
    __kstring_t str = {0, 0, NULL};
    mem_aln2sam(opt, bns, &str, s, 1, aln, 0, NULL);
    std::string out = str.s ? std::string(str.s, str.l) : std::string();
    free(str.s);
    return out;
}

// The float32 a `samtools view -b` of `token` would store, i.e. what the BAM
// writers have to reproduce. strtod-then-narrow mirrors htslib's own SAM aux
// parser, which reads an 'f' field as a double and stores it as a float.
float sam_token_as_float(const std::string &token) {
    return (float)strtod(token.c_str(), NULL);
}

// The token the SAM path is contractually required to produce, computed here
// rather than read back out of the writer. bwa (bwamem.c) and bwa-mem2
// (bwamem.cpp) both render `pa` with "%.3f"; asserting against that keeps the
// oracle independent of bwa_format_pa_value, which this suite also exercises
// through the BAM side. Without it, a regression that reverted BOTH writers to
// the raw quotient would leave every cross-container assertion satisfied.
std::string upstream_pa_token(int score, int alt_sc) {
    char buf[64];
    snprintf(buf, sizeof buf, "%.3f", (double)score / alt_sc);
    return buf;
}

// The aux tag names of a SAM record, in emission order.
std::vector<std::string> sam_tag_names(const std::string &record) {
    std::vector<std::string> names;
    size_t begin = 0;
    for (int field = 0; begin != std::string::npos; ++field) {
        const size_t end = record.find_first_of("\t\n", begin);
        if (field >= 11 && end - begin >= 5)
            names.push_back(record.substr(begin, 2));
        if (end == std::string::npos) break;
        begin = end + 1;
        if (begin >= record.size() || record[end] == '\n') break;
    }
    return names;
}

// The aux tag names of a bam1_t, in stored order. Walks the aux block the way
// htslib's own iterator does rather than probing for known tags, so a tag
// emitted out of order is visible.
std::vector<std::string> bam_tag_names(bam1_t *b) {
    std::vector<std::string> names;
    uint8_t *aux = bam_get_aux(b);
    uint8_t *end = b->data + b->l_data;
    while (aux != NULL && aux + 4 <= end) {
        names.push_back(std::string((const char *)aux, 2));
        uint8_t *next = bam_aux_next(b, aux + 2);
        if (next == NULL) break;
        aux = next - 2;
    }
    return names;
}

// Index of `name` in `names`, or -1.
int index_of(const std::vector<std::string> &names, const std::string &name) {
    for (size_t i = 0; i < names.size(); ++i)
        if (names[i] == name) return (int)i;
    return -1;
}

}  // namespace

TEST_CASE("pa: the --bam value is the float32 of the SAM text token"
          * doctest::test_suite("unit/bam_pa_tag_parity")) {
    // Ratios chosen to cover the three ways the two renderings can disagree:
    // a quotient with more than three decimals (the original defect), exact
    // halves at the fourth decimal (where printf rounds to even and round()
    // does not), and quotients that need no rounding at all.
    struct Case { int score, alt_sc; } cases[] = {
        {100, 124}, {19, 128}, {150, 151}, {39, 48}, {13, 16}, {5, 16},
        {1, 16},    {45, 48},  {7, 80},    {9, 80},  {60, 60}, {3, 16},
        {1, 3},     {2, 3},    {7, 9},     {123, 200},
    };

    mem_opt_t *opt = mem_opt_init();
    RefFixture   ref;
    ReadFixture  read;
    CigarFixture cig;

    for (const Case &c : cases) {
        CAPTURE(c.score);
        CAPTURE(c.alt_sc);

        mem_aln_t aln = make_aln(cig.buf, c.score, c.alt_sc);

        const std::string pa_text = sam_pa(emit_sam(opt, &ref.bns, &read.s, &aln));
        REQUIRE(pa_text != "");
        // The SAM side against its own oracle first, so the cross-container
        // assertion below cannot be satisfied by both writers being wrong.
        CHECK(pa_text == upstream_pa_token(c.score, c.alt_sc));

        bam1_t *b = bam_writer_alloc();
        REQUIRE(b != nullptr);
        REQUIRE(mem_aln_to_bam(b, opt, &ref.bns, &read.s, 1, &aln, 0, NULL) == 0);

        uint8_t *tag = bam_aux_get(b, "pa");
        REQUIRE(tag != nullptr);
        CHECK((float)bam_aux2f(tag) == sam_token_as_float(pa_text));

        bam_writer_free(b);
    }

    free(opt);
}

TEST_CASE("pa: ties at the fourth decimal round the way the SAM text does"
          * doctest::test_suite("unit/bam_pa_tag_parity")) {
    // 39/48 == 0.8125 exactly, so this is a genuine halfway case. printf's
    // "%.3f" rounds it to even ("0.812"); `round(x * 1000) / 1000` rounds it
    // away from zero (0.813). Pinning the literal keeps a future "just use
    // round()" simplification from silently re-splitting the two writers.
    mem_opt_t *opt = mem_opt_init();
    RefFixture   ref;
    ReadFixture  read;
    CigarFixture cig;

    mem_aln_t aln = make_aln(cig.buf, /*score=*/39, /*alt_sc=*/48);

    CHECK(sam_pa(emit_sam(opt, &ref.bns, &read.s, &aln)) == "0.812");

    bam1_t *b = bam_writer_alloc();
    REQUIRE(b != nullptr);
    REQUIRE(mem_aln_to_bam(b, opt, &ref.bns, &read.s, 1, &aln, 0, NULL) == 0);

    uint8_t *tag = bam_aux_get(b, "pa");
    REQUIRE(tag != nullptr);
    CHECK((float)bam_aux2f(tag) == 0.812f);

    bam_writer_free(b);
    free(opt);
}

TEST_CASE("pa: the BAM writer emits pa after SA:Z, as the SAM writer does"
          * doctest::test_suite("unit/bam_pa_tag_parity")) {
    // Both writers emit SA:Z and pa from the same non-secondary block, and the
    // BAM writer used to emit pa several tags earlier — before SA:Z, and before
    // the block was even entered. Order is observable: `samtools view` prints
    // aux fields in stored order, so a --bam run and a SAM run of the same
    // alignment produced differently ordered records. Every other case here
    // passes n_alns == 1, which emits no SA:Z and so cannot see this.
    mem_opt_t *opt = mem_opt_init();
    RefFixture   ref;
    ReadFixture  read;
    CigarFixture cig;

    // Two primary hits: `which = 0` is the record under test, and the second
    // (supplementary, 0x800 — not secondary) is what puts it in SA:Z.
    mem_aln_t alns[2] = {
        make_aln(cig.buf, /*score=*/100, /*alt_sc=*/124),
        make_aln(cig.buf, /*score=*/40,  /*alt_sc=*/0, /*flag=*/0x800),
    };
    alns[1].pos = 500;

    __kstring_t str = {0, 0, NULL};
    mem_aln2sam(opt, &ref.bns, &str, &read.s, 2, alns, 0, NULL);
    const std::string record = str.s ? std::string(str.s, str.l) : std::string();
    free(str.s);

    const std::vector<std::string> sam_names = sam_tag_names(record);
    REQUIRE(index_of(sam_names, "SA") >= 0);
    REQUIRE(index_of(sam_names, "pa") >= 0);
    CHECK(index_of(sam_names, "SA") < index_of(sam_names, "pa"));

    bam1_t *b = bam_writer_alloc();
    REQUIRE(b != nullptr);
    REQUIRE(mem_aln_to_bam(b, opt, &ref.bns, &read.s, 2, alns, 0, NULL) == 0);

    const std::vector<std::string> bam_names = bam_tag_names(b);
    REQUIRE(index_of(bam_names, "SA") >= 0);
    REQUIRE(index_of(bam_names, "pa") >= 0);
    CHECK(index_of(bam_names, "SA") < index_of(bam_names, "pa"));
    // Stronger than the pairwise order above: the whole tag sequence agrees.
    CHECK(bam_names == sam_names);

    bam_writer_free(b);
    free(opt);
}

TEST_CASE("pa: secondary records carry no pa tag in either writer"
          * doctest::test_suite("unit/bam_pa_tag_parity")) {
    // mem_aln2sam emits pa only inside its `!(flag & 0x100)` block, so a
    // secondary record has no pa token to match. The BAM writer emitted one
    // regardless, which is the same SAM/BAM disagreement in a different form.
    mem_opt_t *opt = mem_opt_init();
    RefFixture   ref;
    ReadFixture  read;
    CigarFixture cig;

    mem_aln_t aln = make_aln(cig.buf, /*score=*/100, /*alt_sc=*/124, /*flag=*/0x100);

    CHECK(sam_pa(emit_sam(opt, &ref.bns, &read.s, &aln)) == "");

    bam1_t *b = bam_writer_alloc();
    REQUIRE(b != nullptr);
    REQUIRE(mem_aln_to_bam(b, opt, &ref.bns, &read.s, 1, &aln, 0, NULL) == 0);

    CHECK(bam_aux_get(b, "pa") == nullptr);

    bam_writer_free(b);
    free(opt);
}
