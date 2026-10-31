// Unit tests for the bwa FM-index API adapter (src/aln/aln_bwt.{h,cpp}).
//
// The adapter answers bwa's occ/SA queries from bwa-mem3's own index. These
// tests build a small index from the shared synthetic reference and check the
// adapter's bwa-semantics contract directly: exact-match recovery of known
// positions, the `$`-removed occ identity, 2occ/2occ4 agreement with occ, and
// L2. The row-for-row comparison against a stock `bwa index` is
// test/aln_index_parity_test.

#include "doctest/doctest.h"

#include "aln_test_ref.h"
#include "../../src/FMI_search.h"
#include "../../src/aln/aln_bwt.h"
#include "../../src/bwa.h"

#include <cstdio>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>
#include <unistd.h>

namespace {

namespace ba = bwa_aln;

// Builds and loads the index once for all cases; removes the files at exit.
struct AlnIndex {
    std::string dir, prefix;
    int64_t l_pac = 0;
    std::unique_ptr<FMI_search> fmi;
    ba::bwt_t *bwt = nullptr;

    AlnIndex()
    {
        char tmpl[] = "/tmp/aln_bwt_unit_XXXXXX";
        REQUIRE(mkdtemp(tmpl) != nullptr);
        dir = tmpl;
        prefix = dir + "/ref.fa";
        REQUIRE(aln_test::write_test_fasta(prefix));
        REQUIRE(bwa_idx_build(prefix.c_str(), prefix.c_str(), 0, 3) == 0);
        fmi.reset(new FMI_search(prefix.c_str()));
        fmi->load_index(false, 1);
        bwt = ba::bwt_from_fmi(fmi.get());
        l_pac = (int64_t)bwt->seq_len / 2;
    }
    ~AlnIndex()
    {
        ba::bwt_destroy(bwt);
        fmi.reset();
        std::string cmd = "rm -rf '" + dir + "'";
        if (system(cmd.c_str()) != 0) std::fprintf(stderr, "warning: could not remove %s\n", dir.c_str());
    }
};

AlnIndex &index_fixture()
{
    static AlnIndex idx;
    return idx;
}

} // namespace

TEST_CASE("aln_bwt: exact match of every reference 20-mer recovers its position"
          * doctest::test_suite("unit/aln_bwt")) {
    AlnIndex &ix = index_fixture();
    const std::string c1 = aln_test::test_contigs()[0].second;
    const int len = 20;
    for (size_t i = 0; i + len <= c1.size(); ++i) {
        // bwt_match_exact_alt consumes str from its last base to its first (a
        // standard backward search), so the pattern's own 0..3 codes find the
        // pattern. (bwa aln stores reads reversed and walks them from the front,
        // which is the same backward search of the read as sequenced.)
        std::string pat = c1.substr(i, len);
        ba::ubyte_t str[len];
        for (int j = 0; j < len; ++j)
            str[j] = (ba::ubyte_t)(pat[j] == 'A' ? 0 : pat[j] == 'C' ? 1 : pat[j] == 'G' ? 2 : 3);
        ba::bwtint_t k = 0, l = ix.bwt->seq_len;
        int n = ba::bwt_match_exact_alt(ix.bwt, len, str, &k, &l);
        REQUIRE(n > 0);
        bool found = false;
        for (ba::bwtint_t r = k; r <= l && !found; ++r) found = (ba::bwt_sa(ix.bwt, r) == (ba::bwtint_t)i);
        CHECK_MESSAGE(found, "position " << i << " not among the " << n << " hits");
    }
}

TEST_CASE("aln_bwt: occ is monotone and sums to row count"
          * doctest::test_suite("unit/aln_bwt")) {
    AlnIndex &ix = index_fixture();
    for (int c = 0; c < 4; ++c) CHECK(ba::bwt_occ(ix.bwt, (ba::bwtint_t)(-1), (ba::ubyte_t)c) == 0);
    ba::bwtint_t prev[4] = {0, 0, 0, 0};
    for (ba::bwtint_t k = 0; k <= ix.bwt->seq_len; ++k) {
        ba::bwtint_t sum = 0;
        for (int c = 0; c < 4; ++c) {
            ba::bwtint_t o = ba::bwt_occ(ix.bwt, k, (ba::ubyte_t)c);
            REQUIRE(o >= prev[c]);
            REQUIRE(o <= prev[c] + 1);
            prev[c] = o;
            sum += o;
        }
        // Rows 0..k hold k+1 characters, one of which is the `$` once k reaches primary.
        REQUIRE(sum == k + 1 - (k >= ix.bwt->primary ? 1 : 0));
    }
}

TEST_CASE("aln_bwt: ba::bwt_2occ / ba::bwt_2occ4 agree with ba::bwt_occ"
          * doctest::test_suite("unit/aln_bwt")) {
    AlnIndex &ix = index_fixture();
    const ba::bwtint_t n = ix.bwt->seq_len, p = ix.bwt->primary;
    std::vector<ba::bwtint_t> ks = {(ba::bwtint_t)(-1), p - 1, p, p + 1, 0, n};
    uint32_t state = 5;
    auto rnd = [&]() { state = state * 1103515245u + 12345u; return (ba::bwtint_t)((state >> 8) % (n + 1)); };
    int checked = 0;
    for (int it = 0; it < 10000; ++it) {
        ba::bwtint_t k, l;
        if (it < 36) { k = ks[it / 6]; l = ks[it % 6]; }
        else { k = rnd(); l = rnd(); }
        if (k != (ba::bwtint_t)(-1) && (l == (ba::bwtint_t)(-1) || k > l)) std::swap(k, l); // bwa requires k <= l
        if (k != (ba::bwtint_t)(-1) && l == (ba::bwtint_t)(-1)) continue;
        ba::bwtint_t ok, ol, ck[4], cl[4];
        ba::bwt_2occ4(ix.bwt, k, l, ck, cl);
        for (int c = 0; c < 4; ++c) {
            ba::bwt_2occ(ix.bwt, k, l, (ba::ubyte_t)c, &ok, &ol);
            REQUIRE(ok == ba::bwt_occ(ix.bwt, k, (ba::ubyte_t)c));
            REQUIRE(ol == ba::bwt_occ(ix.bwt, l, (ba::ubyte_t)c));
            REQUIRE(ck[c] == ok);
            REQUIRE(cl[c] == ol);
        }
        ++checked;
    }
    CHECK(checked > 9000);
}

TEST_CASE("aln_bwt: L2 is cumulative base counts of the fwd+revcomp text"
          * doctest::test_suite("unit/aln_bwt")) {
    AlnIndex &ix = index_fixture();
    // Count bases straight from the packed forward strand (.pac, 2 bits per
    // base, N positions already replaced by bwa's pseudo-random bases). Each
    // forward base c also appears as 3 - c on the reverse strand.
    FILE *fp = std::fopen((ix.prefix + ".pac").c_str(), "rb");
    REQUIRE(fp != nullptr);
    std::vector<uint8_t> pac((size_t)(ix.l_pac / 4 + 1));
    REQUIRE(std::fread(pac.data(), 1, pac.size(), fp) == pac.size());
    std::fclose(fp);
    ba::bwtint_t fwd[4] = {0, 0, 0, 0};
    for (int64_t j = 0; j < ix.l_pac; ++j) ++fwd[pac[j >> 2] >> ((~j & 3) << 1) & 3];

    CHECK(ix.bwt->L2[0] == 0);
    CHECK(ix.bwt->L2[4] == ix.bwt->seq_len);
    for (int c = 0; c < 4; ++c)
        CHECK(ix.bwt->L2[c + 1] - ix.bwt->L2[c] == fwd[c] + fwd[3 - c]);
}
