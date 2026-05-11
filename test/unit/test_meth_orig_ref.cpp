// Unit tests for src/meth_orig_ref.{h,cpp} and the meth_chrom_map_t
// basics. The orig-ref module is a lazy view over the doubled-c2t bns +
// pac that recovers original ref bases per slice via dual 2-bit decode +
// the 5-row (f, r) -> original table. Tests build a tiny in-memory
// doubled bns + pac fixture and exercise the slice paths.

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "doctest/doctest.h"

#include "bntseq.h"
#include "meth_bam.h"
#include "meth_orig_ref.h"
#include "meth_orig_ref_fixture.h"

namespace {

struct synth_bns_t {
    bntseq_t bns;
    std::vector<bntann1_t> anns;
    std::vector<bntamb1_t> ambs;
    std::vector<std::string> name_storage;
};

void build_synth_doubled_bns(synth_bns_t &s, int n_real, int per_len) {
    s.bns = bntseq_t{};
    s.bns.l_pac = (int64_t)n_real * 2 * per_len;
    s.bns.n_seqs = n_real * 2;
    s.anns.resize(s.bns.n_seqs);
    s.name_storage.reserve((size_t)s.bns.n_seqs);
    for (int i = 0; i < n_real; ++i) {
        for (int half = 0; half < 2; ++half) {
            int rid = i * 2 + half;
            char prefix = (half == 0) ? 'r' : 'f';
            char buf[32];
            std::snprintf(buf, sizeof(buf), "%cchr%d", prefix, i);
            s.name_storage.push_back(std::string(buf));
            bntann1_t &a = s.anns[rid];
            a.offset = (int64_t)rid * per_len;
            a.len    = per_len;
            a.n_ambs = 0;
            a.gi     = 0;
            a.is_alt = 0;
            a.name   = const_cast<char *>(s.name_storage.back().c_str());
            a.anno   = const_cast<char *>("");
        }
    }
    s.bns.anns = s.anns.data();
    s.bns.ambs = s.ambs.data();
    s.bns.n_holes = 0;
}

}  // namespace

TEST_CASE("meth_chrom_map: collapses f/r doubled contigs to one output per real chrom") {
    synth_bns_t s;
    build_synth_doubled_bns(s, /*n_real=*/3, /*per_len=*/16);

    meth_chrom_map_t *cmap = meth_chrom_map_build_from_bns(&s.bns);
    REQUIRE(cmap != nullptr);
    CHECK(cmap->n_internal == 6);
    CHECK(cmap->n_output == 3);

    CHECK(cmap->out_tid[0] == cmap->out_tid[1]);
    CHECK(cmap->out_tid[2] == cmap->out_tid[3]);
    CHECK(cmap->out_tid[4] == cmap->out_tid[5]);

    CHECK(cmap->direction[0] == 'r');
    CHECK(cmap->direction[1] == 'f');
    CHECK(cmap->direction[2] == 'r');
    CHECK(cmap->direction[3] == 'f');

    CHECK(std::string(cmap->output_names[0]) == "chr0");
    CHECK(std::string(cmap->output_names[1]) == "chr1");
    CHECK(std::string(cmap->output_names[2]) == "chr2");

    meth_chrom_map_free(cmap);
}

TEST_CASE("meth_orig_ref: 5-row recovery returns original bases through slice") {
    /* All four bases per the 5-row table:
     *   pos 0: A -> (f=A, r=A)
     *   pos 1: C -> (f=T, r=C)
     *   pos 2: G -> (f=G, r=A)
     *   pos 3: T -> (f=T, r=T)
     */
    meth_test::OrigRefFixture f("ACGT");
    uint8_t out[4];
    std::memset(out, 0, sizeof(out));
    meth_orig_ref_slice(f.orig, f.real_tid, 0, 4, out);
    CHECK(std::string((char *)out, 4) == "ACGT");
}

TEST_CASE("meth_orig_ref: OOB slice positions emit 'N'") {
    meth_test::OrigRefFixture f("ACGT");
    uint8_t out[8];
    std::memset(out, 0, sizeof(out));
    meth_orig_ref_slice(f.orig, f.real_tid, -2, 6, out);
    /* expect: N N A C G T N N */
    CHECK(out[0] == 'N');
    CHECK(out[1] == 'N');
    CHECK(out[2] == 'A');
    CHECK(out[3] == 'C');
    CHECK(out[4] == 'G');
    CHECK(out[5] == 'T');
    CHECK(out[6] == 'N');
    CHECK(out[7] == 'N');
}

TEST_CASE("meth_orig_ref: ambs intervals fill with 'N'") {
    /* "ACNTACNT" — N at pos 2 and pos 6. The fixture marks both via
     * bns->ambs on each half of the doubled pac, just like bwa indexer
     * does at index time. */
    meth_test::OrigRefFixture f("ACNTACNT");
    uint8_t out[8];
    std::memset(out, 0, sizeof(out));
    meth_orig_ref_slice(f.orig, f.real_tid, 0, 8, out);
    CHECK(std::string((char *)out, 8) == "ACNTACNT");
}

TEST_CASE("meth_orig_ref_load: NULL bns/pac returns NULL") {
    meth_chrom_map_t cmap{};
    cmap.n_internal = 0;
    cmap.n_output = 0;
    CHECK(meth_orig_ref_load(nullptr, nullptr, &cmap) == nullptr);
}

TEST_CASE("OrigRefFixture: bns.ambs is sorted by offset (bns_iter_ambi precondition)") {
    /* bns_iter_ambi binary-searches bns->ambs assuming non-decreasing
     * offset+len. With 2+ 'N's in the source sequence the fixture must
     * lay out all r-X (low-offset) intervals before all f-X intervals
     * to keep that invariant. Regression for the interleaved-push bug. */
    meth_test::OrigRefFixture f("NNANNA");
    REQUIRE(f.bns.n_holes >= 2);
    for (int i = 1; i < f.bns.n_holes; ++i) {
        CHECK(f.bns.ambs[i - 1].offset <= f.bns.ambs[i].offset);
    }
}

TEST_CASE("meth_orig_ref: multi-N forward sequence masks every N via slice") {
    /* End-to-end check that exercises bns_iter_ambi over an ambs array
     * with multiple entries per side. The previous interleaved fixture
     * ordering caused bns_iter_ambi's binary search to skip valid
     * overlapping intervals (and, with 3+ N's, walk into the wrong half
     * and translate to negative dst[] indices). */
    meth_test::OrigRefFixture f("NNNAN");
    uint8_t out[5];
    std::memset(out, 0, sizeof(out));
    meth_orig_ref_slice(f.orig, f.real_tid, 0, 5, out);
    CHECK(std::string((char *)out, 5) == "NNNAN");
}
