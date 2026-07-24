// Unit tests for bns_pos2rid and its bns_build_pos2rid acceleration table.
//
// bns_pos2rid answers "which contig covers packed forward position pos_f?".
// Two code paths implement the same predicate — the largest rid with
// anns[rid].offset <= pos_f — and which one runs depends only on whether the
// bucket table was built:
//
//   * table built   -> bucket bracket + bounded binary search
//   * table NULL    -> the original plain binary search over anns[]
//
// These tests pin the contract rather than either implementation: every query
// is checked against an independent linear reference, and the two paths are
// checked against each other on the same bns. The dense layout is the one that
// matters for the bucket path's cost bound — thousands of contig starts inside
// a single bucket window.

#include "doctest/doctest.h"

#include "bntseq.h"

#include <cstdlib>
#include <cstring>
#include <vector>

namespace {

/* A minimal bntseq_t owning only what bns_pos2rid reads: l_pac, n_seqs and
 * anns[].offset. Deliberately not built via bns_restore — no files, no pac,
 * and no name/anno strings to free — so the table logic can be exercised over
 * contig layouts that no real index on disk would give us. */
struct FakeBns {
    bntseq_t bns;
    std::vector<bntann1_t> anns;

    /* Lay out `n_seqs` contigs of the given lengths back to back from offset 0. */
    explicit FakeBns(const std::vector<int64_t> &contig_lengths) {
        std::memset(&bns, 0, sizeof(bns));
        anns.assign(contig_lengths.size(), bntann1_t());
        int64_t offset = 0;
        for (size_t i = 0; i < contig_lengths.size(); ++i) {
            std::memset(&anns[i], 0, sizeof(anns[i]));
            anns[i].offset = offset;
            anns[i].len    = (int32_t)contig_lengths[i];
            offset += contig_lengths[i];
        }
        bns.l_pac  = offset;
        bns.n_seqs = (int32_t)contig_lengths.size();
        bns.anns   = anns.data();
        bns.pos2rid_bucket = NULL;
    }

    ~FakeBns() { free(bns.pos2rid_bucket); }

    FakeBns(const FakeBns &)            = delete;
    FakeBns &operator=(const FakeBns &) = delete;

    /* Independent reference: scan every contig and keep the last one that
     * starts at or before pos_f. O(n_seqs) and obviously correct. */
    int expected_rid(int64_t pos_f) const {
        int rid = 0;
        for (int i = 0; i < bns.n_seqs; ++i)
            if (anns[i].offset <= pos_f) rid = i;
        return rid;
    }
};

/* Positions worth probing for a given layout: every contig boundary and its
 * immediate neighbours, both ends of the genome, and a coarse stride sweep
 * that lands inside contig bodies. */
std::vector<int64_t> probe_positions(const FakeBns &fake) {
    std::vector<int64_t> positions;
    for (int i = 0; i < fake.bns.n_seqs; ++i) {
        int64_t start = fake.anns[i].offset;
        positions.push_back(start - 1);
        positions.push_back(start);
        positions.push_back(start + 1);
    }
    positions.push_back(0);
    positions.push_back(fake.bns.l_pac - 1);
    const int64_t stride = fake.bns.l_pac / 97 + 1;
    for (int64_t p = 0; p < fake.bns.l_pac; p += stride) positions.push_back(p);
    return positions;
}

/* Assert both code paths satisfy the contract on every probe position: build
 * the table and query, then drop the table and query again, comparing each
 * against the linear reference (and hence against each other). */
void check_both_paths(FakeBns &fake) {
    const std::vector<int64_t> positions = probe_positions(fake);

    bns_build_pos2rid(&fake.bns);
    REQUIRE(fake.bns.pos2rid_bucket != NULL);
    std::vector<int> via_bucket;
    via_bucket.reserve(positions.size());
    for (size_t i = 0; i < positions.size(); ++i)
        via_bucket.push_back(bns_pos2rid(&fake.bns, positions[i]));

    free(fake.bns.pos2rid_bucket);
    fake.bns.pos2rid_bucket = NULL;
    for (size_t i = 0; i < positions.size(); ++i) {
        const int64_t pos = positions[i];
        const int expected = pos < 0 ? 0 : fake.expected_rid(pos);
        CAPTURE(pos);
        CHECK(via_bucket[i] == expected);
        CHECK(bns_pos2rid(&fake.bns, pos) == expected);
    }
}

/* Contig widths relative to BNS_POS2RID_SHIFT drive which regime the bucket
 * path lands in, so express the layouts in terms of the bucket width itself
 * instead of hard-coding a stride that a future shift change would silently
 * invalidate. */
const int64_t kBucketWidth = (int64_t)1 << BNS_POS2RID_SHIFT;

std::vector<int64_t> uniform_contigs(int n_seqs, int64_t len) {
    return std::vector<int64_t>((size_t)n_seqs, len);
}

} // namespace

TEST_CASE("bns_pos2rid: sparse layout (contigs much wider than a bucket)") {
    // The hg38-like regime: most buckets contain no contig start at all.
    FakeBns fake(uniform_contigs(64, kBucketWidth * 7 + 123));
    check_both_paths(fake);
}

TEST_CASE("bns_pos2rid: dense layout (thousands of contig starts per bucket)") {
    // The regime a panel/transcriptome/adapter reference produces: contigs far
    // narrower than a bucket, so one bucket window brackets many contig starts.
    // This is what a linear scan from the bucket entry would walk.
    FakeBns fake(uniform_contigs(4096, 8));
    REQUIRE(fake.bns.l_pac / kBucketWidth < 4); // genuinely dense, not spread out
    check_both_paths(fake);
}

TEST_CASE("bns_pos2rid: every contig start lands on a bucket boundary") {
    // Offsets exactly equal to a bucket's first position exercise the <= in
    // both the table build and the lookup.
    FakeBns fake(uniform_contigs(32, kBucketWidth));
    check_both_paths(fake);
}

TEST_CASE("bns_pos2rid: ragged contig widths straddling bucket boundaries") {
    std::vector<int64_t> lengths;
    for (int i = 0; i < 200; ++i)
        lengths.push_back(1 + (int64_t)((i * 2654435761u) % (kBucketWidth * 2)));
    FakeBns fake(lengths);
    check_both_paths(fake);
}

TEST_CASE("bns_pos2rid: single contig") {
    FakeBns fake(uniform_contigs(1, kBucketWidth * 3));
    check_both_paths(fake);
}

TEST_CASE("bns_pos2rid: positions at or past l_pac return -1 on both paths") {
    FakeBns fake(uniform_contigs(16, 1000));
    const int64_t l_pac = fake.bns.l_pac;

    bns_build_pos2rid(&fake.bns);
    REQUIRE(fake.bns.pos2rid_bucket != NULL);
    CHECK(bns_pos2rid(&fake.bns, l_pac) == -1);
    CHECK(bns_pos2rid(&fake.bns, l_pac + kBucketWidth * 4) == -1);

    free(fake.bns.pos2rid_bucket);
    fake.bns.pos2rid_bucket = NULL;
    CHECK(bns_pos2rid(&fake.bns, l_pac) == -1);
    CHECK(bns_pos2rid(&fake.bns, l_pac + kBucketWidth * 4) == -1);
}

TEST_CASE("bns_build_pos2rid: degenerate inputs leave the table unbuilt") {
    // A NULL table is the documented "fall back to binary search" signal, so
    // each of these must leave it NULL rather than allocate something empty.
    SUBCASE("NULL bns") {
        bns_build_pos2rid(NULL); // must not crash
    }
    SUBCASE("empty genome") {
        const std::vector<int64_t> no_contigs;
        FakeBns fake(no_contigs);
        bns_build_pos2rid(&fake.bns);
        CHECK(fake.bns.pos2rid_bucket == NULL);
    }
    SUBCASE("l_pac of zero") {
        FakeBns fake(uniform_contigs(4, 0));
        bns_build_pos2rid(&fake.bns);
        CHECK(fake.bns.pos2rid_bucket == NULL);
    }
}

TEST_CASE("bns_build_pos2rid: is idempotent") {
    FakeBns fake(uniform_contigs(8, kBucketWidth * 2));
    bns_build_pos2rid(&fake.bns);
    REQUIRE(fake.bns.pos2rid_bucket != NULL);
    int32_t *first = fake.bns.pos2rid_bucket;
    bns_build_pos2rid(&fake.bns);
    CHECK(fake.bns.pos2rid_bucket == first); // second call must not reallocate
}
