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
#include "meth_orig_ref_fixture.h"

#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <random>
#include <string>
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

GenCigarResult gen_cigar(const meth_test::OrigRefFixture &ref, const int8_t *mat,
                         const std::string &query_in_ref_frame, int64_t rb, int64_t re,
                         int band, int nm_from_mat)
{
    std::vector<uint8_t> query(query_in_ref_frame.size());
    for (size_t i = 0; i < query.size(); ++i) query[i] = nt4(query_in_ref_frame[i]);
    int score = 0, n_cigar = 0, nm = -1;
    uint32_t *cigar = bwa_gen_cigar3(mat, 6, 1, 6, 1, band, ref.bns.l_pac, ref.orig_pac,
                                     (int)query.size(), query.data(), rb, re, &score,
                                     &n_cigar, &nm, nm_from_mat);
    GenCigarResult r;
    REQUIRE(cigar != nullptr);
    r.cigar.assign(cigar, cigar + n_cigar);
    r.md = (const char *)(cigar + n_cigar);
    r.nm = nm;
    free(cigar);
    return r;
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
    o->meth_scoring = scoring;
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
