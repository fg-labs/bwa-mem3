// Unit tests for how bwa_gen_cigar3 derives NM and MD under --meth.
//
// NM counts an aligned column only when the scoring matrix penalises it, so a
// bisulfite conversion (the cell the per-hypothesis matrix frees) is not an
// edit; inserted and deleted bases always count. MD is always literal: it lists
// every reference base the read differs from, conversions included, so CIGAR +
// SEQ + MD rebuilds the real reference and the conversions stay readable from
// the record. NM is therefore smaller than the edit count MD and CIGAR imply by
// exactly the number of matrix-freed columns.

#include "doctest/doctest.h"
#include "bwa.h"
#include "bwamem.h"
#include "meth_model.h"
#include "meth_orig_ref_fixture.h"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <random>
#include <string>
#include <utility>
#include <vector>

namespace {

const uint8_t kNt4Ambig = 4;

uint8_t nt4(char c)
{
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default:  return kNt4Ambig;
    }
}

char complement(char c)
{
    switch (c) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default:  return 'N';
    }
}

std::string reverse_complement(const std::string &s)
{
    std::string rc(s.rbegin(), s.rend());
    for (char &c : rc) c = complement(c);
    return rc;
}

std::string random_bases(size_t len, uint32_t seed)
{
    std::mt19937 rng(seed);
    std::string s(len, 'A');
    for (char &c : s) c = "ACGT"[rng() & 3];
    return s;
}

// Replace every `from` with `to`, returning how many bases changed.
int convert_all(std::string &s, char from, char to)
{
    int n_converted = 0;
    for (char &c : s)
        if (c == from) c = to, ++n_converted;
    return n_converted;
}

// The output of one bwa_gen_cigar3 call, decoded.
struct GenCigarResult {
    std::vector<uint32_t> cigar;
    std::string md;
    int nm = -1;
};

std::vector<uint8_t> nt4_codes(const std::string &bases)
{
    std::vector<uint8_t> codes(bases.size());
    for (size_t i = 0; i < codes.size(); ++i) codes[i] = nt4(bases[i]);
    return codes;
}

// Decode (and free) a generator's buffer: the n_cigar ops, then the MD.
GenCigarResult take_cigar(uint32_t *cigar, int n_cigar, int nm)
{
    GenCigarResult r;
    REQUIRE(cigar != nullptr);
    r.cigar.assign(cigar, cigar + n_cigar);
    r.md = (const char *)(cigar + n_cigar);
    r.nm = nm;
    free(cigar);
    return r;
}

GenCigarResult gen_cigar(const meth_test::OrigRefFixture &ref, const int8_t *mat,
                         const std::string &query_in_ref_frame, int64_t rb, int64_t re,
                         int band, int nm_from_mat)
{
    std::vector<uint8_t> query = nt4_codes(query_in_ref_frame);
    int score = 0, n_cigar = 0, nm = -1;
    uint32_t *cigar = bwa_gen_cigar3(mat, 6, 1, 6, 1, band, ref.bns.l_pac, ref.orig_pac,
                                     (int)query.size(), query.data(), rb, re, &score,
                                     &n_cigar, &nm, nm_from_mat);
    return take_cigar(cigar, n_cigar, nm);
}

// As gen_cigar, through bwa_gen_cigar4 with an explicit `want_md` (1 = the
// literal MD, 2 = --md short) and alphabet size `m`.
GenCigarResult gen_cigar_md(const meth_test::OrigRefFixture &ref, const int8_t *mat,
                            const std::string &query_in_ref_frame, int64_t rb, int64_t re,
                            int band, int nm_from_mat, int want_md, int m = 5)
{
    std::vector<uint8_t> query = nt4_codes(query_in_ref_frame);
    int score = 0, n_cigar = 0, nm = -1;
    uint32_t *cigar = bwa_gen_cigar4(mat, m, &ref.bns, 6, 1, 6, 1, band, ref.bns.l_pac,
                                     ref.orig_pac, (int)query.size(), query.data(), rb, re,
                                     &score, &n_cigar, &nm, nm_from_mat, want_md);
    return take_cigar(cigar, n_cigar, nm);
}

// The MD of an ungapped alignment of `len` columns that lists exactly the
// columns in `listed` (column -> reference base, in forward-genome order).
std::string md_listing(size_t len, const std::vector<std::pair<size_t, char>> &listed)
{
    std::string md;
    size_t next = 0;
    for (const auto &col : listed) {
        md += std::to_string(col.first - next) + col.second;
        next = col.first + 1;
    }
    return md + std::to_string(len - next);
}

// The reference an alignment claims, rebuilt from SEQ (forward-genome
// orientation), CIGAR and MD as `samtools calmd` would, plus how many columns
// MD lists although SEQ already has that base there. A literal MD lists exactly
// the columns where SEQ differs from the reference, so a correct MD has none of
// those and rebuilds the real reference.
struct MdRebuild {
    std::string reference;
    int n_listed_where_seq_matches = 0;
};

MdRebuild rebuild_reference(const std::string &seq, const std::vector<uint32_t> &cigar,
                            const std::string &md)
{
    MdRebuild r;
    std::string &reference = r.reference;
    size_t q = 0, m = 0;
    int pending_matches = 0;
    auto next_md_token = [&]() {
        while (m < md.size() && isdigit((unsigned char)md[m]))
            pending_matches = pending_matches * 10 + (md[m++] - '0');
    };
    for (uint32_t c : cigar) {
        const int op = c & 0xf, len = (int)(c >> 4);
        if (op == 0) {
            for (int i = 0; i < len; ++i, ++q) {
                next_md_token();
                if (pending_matches > 0) {
                    reference += seq[q];
                    --pending_matches;
                } else {
                    if (md[m] == seq[q]) ++r.n_listed_where_seq_matches;
                    reference += md[m++];
                }
            }
        } else if (op == 1) {
            q += len;
        } else if (op == 2) {
            next_md_token();
            REQUIRE(pending_matches == 0);
            REQUIRE(md[m] == '^');
            ++m;
            for (int i = 0; i < len; ++i) reference += md[m++];
        }
    }
    next_md_token();
    CHECK(pending_matches == 0);
    CHECK(m == md.size());
    return r;
}

// MD must list exactly the columns where `seq` differs from `expected_reference`.
void check_md_is_literal(const std::string &seq, const GenCigarResult &r,
                         const std::string &expected_reference)
{
    const MdRebuild rebuilt = rebuild_reference(seq, r.cigar, r.md);
    CHECK(rebuilt.reference == expected_reference);
    CHECK(rebuilt.n_listed_where_seq_matches == 0);
}

// The reference bases MD lists at mismatched columns (deletions excluded).
std::string md_mismatch_bases(const std::string &md)
{
    std::string bases;
    bool in_deletion = false;
    for (char c : md) {
        if (c == '^') in_deletion = true;
        else if (isdigit((unsigned char)c)) in_deletion = false;
        else if (!in_deletion) bases += c;
    }
    return bases;
}

int count_md_mismatches(const std::string &md)
{
    return (int)md_mismatch_bases(md).size();
}

int count_cigar_gaps(const std::vector<uint32_t> &cigar)
{
    int n = 0;
    for (uint32_t c : cigar)
        if ((c & 0xf) == 1 || (c & 0xf) == 2) n += (int)(c >> 4);
    return n;
}

mem_opt_t *meth_opt(int scoring)
{
    mem_opt_t *o = mem_opt_init();
    o->meth.meth_scoring = scoring;
    mem_opt_fill_meth_mat(o);
    return o;
}

// One forward contig, with an OT read over [kBeg, kEnd) whose every C is a
// conversion (C->T) and whose first A carries a real A->G variant.
const int64_t kContigLen = 160, kBeg = 20, kEnd = 140;

struct OtReadCase {
    std::string genome = random_bases(kContigLen, 20260930);
    meth_test::OrigRefFixture ref{genome};
    std::string region = genome.substr(kBeg, kEnd - kBeg);
    std::string read = region;
    int n_conversions = convert_all(read, 'C', 'T');
    size_t variant_at = region.find('A');

    OtReadCase()
    {
        REQUIRE(n_conversions > 10);
        REQUIRE(variant_at != std::string::npos);
        read[variant_at] = 'G';
    }
};

}  // namespace

TEST_CASE("--meth: NM skips conversions while MD lists them (ungapped, forward)"
          * doctest::test_suite("unit/meth_nm_md"))
{
    OtReadCase c;
    mem_opt_t *o = meth_opt(MEM_METH_SCORING_GENOMIC);

    const GenCigarResult r = gen_cigar(c.ref, o->mat_ot, c.read, kBeg, kEnd, 0, 1);

    CHECK(r.nm == 1);
    CHECK(count_md_mismatches(r.md) == c.n_conversions + 1);
    check_md_is_literal(c.read, r, c.region);
    free(o);
}

TEST_CASE("--meth: MD stays literal on the reverse strand"
          * doctest::test_suite("unit/meth_nm_md"))
{
    // A region on the reverse strand is fetched as the reverse complement of the
    // forward bases, so the read is built in that frame. MD is written in
    // forward-genome orientation, against SEQ = reverse complement of the read.
    OtReadCase c;
    mem_opt_t *o = meth_opt(MEM_METH_SCORING_GENOMIC);
    const int64_t l_pac = c.ref.bns.l_pac;
    std::string read_in_ref_frame = reverse_complement(c.region);
    const int n_conversions = convert_all(read_in_ref_frame, 'C', 'T');

    const GenCigarResult r = gen_cigar(c.ref, o->mat_ot, read_in_ref_frame,
                                       2 * l_pac - kEnd, 2 * l_pac - kBeg, 0, 1);

    CHECK(r.nm == 0);
    CHECK(count_md_mismatches(r.md) == n_conversions);
    check_md_is_literal(reverse_complement(read_in_ref_frame), r, c.region);
    free(o);
}

TEST_CASE("--meth: a gapped alignment counts the gap in NM and keeps conversions in MD"
          * doctest::test_suite("unit/meth_nm_md"))
{
    OtReadCase c;
    mem_opt_t *o = meth_opt(MEM_METH_SCORING_GENOMIC);
    const int deletion_len = 3, deletion_at = 50;
    std::string read =
        c.region.substr(0, deletion_at) + c.region.substr(deletion_at + deletion_len);
    const int n_conversions = convert_all(read, 'C', 'T');

    const GenCigarResult r = gen_cigar(c.ref, o->mat_ot, read, kBeg, kEnd, 10, 1);

    REQUIRE(count_cigar_gaps(r.cigar) == deletion_len);
    CHECK(r.nm == deletion_len);
    // Where the aligner places the deletion decides which bases it swallows, so
    // the number of listed conversions is not pinned; that MD lists exactly the
    // differing columns, and that each is a converted ref C, is.
    check_md_is_literal(read, r, c.region);
    const std::string listed = md_mismatch_bases(r.md);
    CHECK(listed.find_first_not_of('C') == std::string::npos);
    CHECK((int)listed.size() <= n_conversions);
    free(o);
}

TEST_CASE("--meth-scoring neutral: a zero-scored conversion is not an edit but is listed"
          * doctest::test_suite("unit/meth_nm_md"))
{
    OtReadCase c;
    mem_opt_t *o = meth_opt(MEM_METH_SCORING_NEUTRAL);

    const GenCigarResult r = gen_cigar(c.ref, o->mat_ot, c.read, kBeg, kEnd, 0, 1);

    CHECK(r.nm == 1);
    CHECK(count_md_mismatches(r.md) == c.n_conversions + 1);
    check_md_is_literal(c.read, r, c.region);
    free(o);
}

TEST_CASE("--meth-scoring collapsed: a mirror-cell variant is hidden from NM, not from MD"
          * doctest::test_suite("unit/meth_nm_md"))
{
    OtReadCase c;
    mem_opt_t *o = meth_opt(MEM_METH_SCORING_COLLAPSED);
    std::string read = c.region;
    const size_t t_at = read.find('T');
    REQUIRE(t_at != std::string::npos);
    read[t_at] = 'C';   // ref T, read C: the mirror cell collapsed frees

    const GenCigarResult r = gen_cigar(c.ref, o->mat_ot, read, kBeg, kEnd, 0, 1);

    CHECK(r.nm == 0);
    CHECK(count_md_mismatches(r.md) == 1);
    check_md_is_literal(read, r, c.region);
    free(o);
}

TEST_CASE("--meth-scoring genomic: a mirror-cell variant is an edit and is listed"
          * doctest::test_suite("unit/meth_nm_md"))
{
    // genomic frees only the conversion cell, so ref T x read C stays penalised:
    // NM follows the matrix, not a fixed C/T rule.
    OtReadCase c;
    mem_opt_t *o = meth_opt(MEM_METH_SCORING_GENOMIC);
    std::string read = c.read;
    const size_t t_at = c.region.find('T');
    REQUIRE(t_at != std::string::npos);
    REQUIRE(t_at != c.variant_at);
    read[t_at] = 'C';

    const GenCigarResult r = gen_cigar(c.ref, o->mat_ot, read, kBeg, kEnd, 0, 1);

    CHECK(r.nm == 2);
    CHECK(count_md_mismatches(r.md) == c.n_conversions + 2);
    check_md_is_literal(read, r, c.region);
    free(o);
}

TEST_CASE("--meth OB: G->A conversions are not edits but are listed"
          * doctest::test_suite("unit/meth_nm_md"))
{
    OtReadCase c;
    mem_opt_t *o = meth_opt(MEM_METH_SCORING_GENOMIC);
    std::string read = c.region;
    const int n_conversions = convert_all(read, 'G', 'A');
    REQUIRE(n_conversions > 10);

    const GenCigarResult r = gen_cigar(c.ref, o->mat_ob, read, kBeg, kEnd, 0, 1);

    CHECK(r.nm == 0);
    CHECK(md_mismatch_bases(r.md) == std::string(n_conversions, 'G'));
    check_md_is_literal(read, r, c.region);
    free(o);
}

TEST_CASE("--meth: a read N is an edit and is listed"
          * doctest::test_suite("unit/meth_nm_md"))
{
    OtReadCase c;
    mem_opt_t *o = meth_opt(MEM_METH_SCORING_GENOMIC);
    std::string read = c.read;
    read[c.variant_at] = 'N';   // in place of the A->G variant

    const GenCigarResult r = gen_cigar(c.ref, o->mat_ot, read, kBeg, kEnd, 0, 1);

    CHECK(r.nm == 1);
    CHECK(count_md_mismatches(r.md) == c.n_conversions + 1);
    check_md_is_literal(read, r, c.region);
    free(o);
}

TEST_CASE("without --meth NM and MD both count every literal difference"
          * doctest::test_suite("unit/meth_nm_md"))
{
    OtReadCase c;
    mem_opt_t *o = mem_opt_init();

    const GenCigarResult r = gen_cigar(c.ref, o->mat, c.read, kBeg, kEnd, 0, 0);

    CHECK(r.nm == c.n_conversions + 1);
    CHECK(count_md_mismatches(r.md) == c.n_conversions + 1);
    check_md_is_literal(c.read, r, c.region);
    free(o);
}

// --md short (bwa_gen_cigar4 want_md 2): MD lists only the mismatches NM
// counts, so under --meth a matrix-freed conversion reads as a match and
// CIGAR + SEQ + MD rebuilds the read's own bases at those columns. Its
// mismatch count is therefore NM minus the gap bases, and with nm_from_mat 0
// it is the literal MD.

TEST_CASE("--md short: a converted read lists only its real variant (forward)"
          * doctest::test_suite("unit/meth_nm_md"))
{
    OtReadCase c;
    mem_opt_t *o = meth_opt(MEM_METH_SCORING_GENOMIC);

    const GenCigarResult full = gen_cigar_md(c.ref, o->mat_ot, c.read, kBeg, kEnd, 0, 1, 1);
    const GenCigarResult shrt = gen_cigar_md(c.ref, o->mat_ot, c.read, kBeg, kEnd, 0, 1, 2);

    CHECK(shrt.cigar == full.cigar);
    CHECK(shrt.nm == full.nm);
    CHECK(shrt.nm == 1);
    CHECK(shrt.md == std::to_string(c.variant_at) + "A"
                     + std::to_string(c.region.size() - c.variant_at - 1));
    // The rebuilt reference is the read with only the variant column restored.
    std::string expected = c.read;
    expected[c.variant_at] = 'A';
    CHECK(rebuild_reference(c.read, shrt.cigar, shrt.md).reference == expected);
    CHECK(shrt.md.size() < full.md.size());
    free(o);
}

TEST_CASE("--md short: a fully converted read is all matches (reverse strand)"
          * doctest::test_suite("unit/meth_nm_md"))
{
    OtReadCase c;
    mem_opt_t *o = meth_opt(MEM_METH_SCORING_GENOMIC);
    const int64_t l_pac = c.ref.bns.l_pac;
    std::string read_in_ref_frame = reverse_complement(c.region);
    REQUIRE(convert_all(read_in_ref_frame, 'C', 'T') > 10);

    const GenCigarResult r = gen_cigar_md(c.ref, o->mat_ot, read_in_ref_frame,
                                          2 * l_pac - kEnd, 2 * l_pac - kBeg, 0, 1, 2);

    CHECK(r.nm == 0);
    CHECK(r.md == std::to_string(kEnd - kBeg));
    free(o);
}

TEST_CASE("--md short: a real variant on the reverse strand is listed at its forward column"
          * doctest::test_suite("unit/meth_nm_md"))
{
    // The fetch frame is the reverse complement, so a ref-frame A->G at frame
    // column a is forward column len-1-a, listed as the forward base T.
    OtReadCase c;
    mem_opt_t *o = meth_opt(MEM_METH_SCORING_GENOMIC);
    const int64_t l_pac = c.ref.bns.l_pac;
    const size_t len = c.region.size();
    std::string read_in_ref_frame = reverse_complement(c.region);
    const size_t a_at = read_in_ref_frame.find('A');
    REQUIRE(a_at != std::string::npos);
    REQUIRE(convert_all(read_in_ref_frame, 'C', 'T') > 10);
    read_in_ref_frame[a_at] = 'G';

    const GenCigarResult r = gen_cigar_md(c.ref, o->mat_ot, read_in_ref_frame,
                                          2 * l_pac - kEnd, 2 * l_pac - kBeg, 0, 1, 2);

    CHECK(r.nm == 1);
    CHECK(r.md == md_listing(len, {{len - 1 - a_at, 'T'}}));
    free(o);
}

TEST_CASE("--md short: OB G->A conversions are not listed, a C->T variant is"
          * doctest::test_suite("unit/meth_nm_md"))
{
    // genomic OB frees only ref G x read A, so a C->T at a ref C stays an edit.
    OtReadCase c;
    mem_opt_t *o = meth_opt(MEM_METH_SCORING_GENOMIC);
    std::string read = c.region;
    REQUIRE(convert_all(read, 'G', 'A') > 10);

    const GenCigarResult clean = gen_cigar_md(c.ref, o->mat_ob, read, kBeg, kEnd, 0, 1, 2);
    CHECK(clean.nm == 0);
    CHECK(clean.md == std::to_string(kEnd - kBeg));

    const size_t c_at = c.region.find('C');
    REQUIRE(c_at != std::string::npos);
    read[c_at] = 'T';
    const GenCigarResult r = gen_cigar_md(c.ref, o->mat_ob, read, kBeg, kEnd, 0, 1, 2);
    CHECK(r.nm == 1);
    CHECK(r.md == md_listing(c.region.size(), {{c_at, 'C'}}));
    free(o);
}

TEST_CASE("--md short: a gapped alignment keeps the deletion and drops the conversions"
          * doctest::test_suite("unit/meth_nm_md"))
{
    OtReadCase c;
    mem_opt_t *o = meth_opt(MEM_METH_SCORING_GENOMIC);
    const int deletion_len = 3, deletion_at = 50;
    std::string read =
        c.region.substr(0, deletion_at) + c.region.substr(deletion_at + deletion_len);
    REQUIRE(convert_all(read, 'C', 'T') > 10);

    const GenCigarResult full = gen_cigar_md(c.ref, o->mat_ot, read, kBeg, kEnd, 10, 1, 1);
    const GenCigarResult shrt = gen_cigar_md(c.ref, o->mat_ot, read, kBeg, kEnd, 10, 1, 2);

    REQUIRE(count_cigar_gaps(shrt.cigar) == deletion_len);
    CHECK(shrt.cigar == full.cigar);
    CHECK(shrt.nm == deletion_len);
    CHECK(count_md_mismatches(shrt.md) == 0);
    CHECK(std::count(shrt.md.begin(), shrt.md.end(), '^') == 1);
    // The deleted reference bases are still listed, so the rebuilt reference
    // has the region's length.
    CHECK(rebuild_reference(read, shrt.cigar, shrt.md).reference.size() == c.region.size());
    free(o);
}

TEST_CASE("--md short: in every scoring mode its mismatches are NM's"
          * doctest::test_suite("unit/meth_nm_md"))
{
    // collapsed also frees the mirror cell, so a ref-T read-C variant drops out
    // of MD as it does from NM; genomic and neutral keep it.
    const int modes[] = {MEM_METH_SCORING_COLLAPSED, MEM_METH_SCORING_GENOMIC,
                         MEM_METH_SCORING_NEUTRAL};
    for (int mode : modes) {
        CAPTURE(mode);
        OtReadCase c;
        mem_opt_t *o = meth_opt(mode);
        std::string read = c.read;
        const size_t t_at = c.region.find('T');
        REQUIRE(t_at != std::string::npos);
        REQUIRE(t_at != c.variant_at);
        read[t_at] = 'C';

        const GenCigarResult r = gen_cigar_md(c.ref, o->mat_ot, read, kBeg, kEnd, 0, 1, 2);

        CHECK(count_md_mismatches(r.md) == r.nm);
        CHECK(r.nm == (mode == MEM_METH_SCORING_COLLAPSED ? 1 : 2));
        // The exact columns: the A->G variant always; the mirror T->C only
        // where the matrix penalises it.
        std::vector<std::pair<size_t, char>> listed{{c.variant_at, 'A'}};
        if (mode != MEM_METH_SCORING_COLLAPSED) listed.push_back({t_at, 'T'});
        std::sort(listed.begin(), listed.end());
        CHECK(r.md == md_listing(c.region.size(), listed));
        free(o);
    }
}

TEST_CASE("--md short: without --meth it is the literal MD"
          * doctest::test_suite("unit/meth_nm_md"))
{
    OtReadCase c;
    mem_opt_t *o = mem_opt_init();

    const GenCigarResult full = gen_cigar_md(c.ref, o->mat, c.read, kBeg, kEnd, 0, 0, 1);
    const GenCigarResult shrt = gen_cigar_md(c.ref, o->mat, c.read, kBeg, kEnd, 0, 0, 2);

    CHECK(shrt.md == full.md);
    CHECK(shrt.nm == full.nm);
    CHECK(shrt.cigar == full.cigar);
    free(o);
}

TEST_CASE("--md short: under --meth-context cpg only a CpG conversion drops out"
          * doctest::test_suite("unit/meth_nm_md"))
{
    // The 7-code alphabet labels a C followed by G as C5, whose conversion the
    // OT matrix frees; a conversion at any other C stays an edit, so short MD
    // keeps it, and its mismatch count is NM.
    OtReadCase c;
    mem_opt_t *o = mem_opt_init();
    o->meth.meth_mode = 1;
    o->meth.meth_context = METH_CTX_CPG;
    o->meth.meth_scoring = MEM_METH_SCORING_GENOMIC;
    mem_opt_fill_meth_mat(o);
    const size_t cpg = c.region.find("CG");
    REQUIRE(cpg != std::string::npos);
    size_t plain = std::string::npos;
    for (size_t i = 0; i + 1 < c.region.size(); ++i)
        if (c.region[i] == 'C' && c.region[i + 1] != 'G') { plain = i; break; }
    REQUIRE(plain != std::string::npos);
    std::string read = c.region;
    read[cpg] = 'T';
    read[plain] = 'T';

    const GenCigarResult full = gen_cigar_md(c.ref, o->mat7_ot, read, kBeg, kEnd, 0, 1, 1, METH_M);
    const GenCigarResult shrt = gen_cigar_md(c.ref, o->mat7_ot, read, kBeg, kEnd, 0, 1, 2, METH_M);

    CHECK(full.nm == 1);
    CHECK(shrt.nm == 1);
    CHECK(count_md_mismatches(full.md) == 2);
    CHECK(shrt.md == std::to_string(plain) + "C" + std::to_string(c.region.size() - plain - 1));
    free(o);
}
