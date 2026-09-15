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
#include <utility>
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

// --- bns_intv2rid --------------------------------------------------------
//
// bns_intv2rid maps a packed interval [rb, re) to the single contig covering
// it, or a sentinel: -2 when the interval straddles the forward/reverse split,
// -1 when the two ends fall in different contigs. It has a single-bucket fast
// path that, when both depos'd ends land in one bucket whose bracket is empty,
// returns that rid without a second bns_pos2rid call. That fast path must be
// byte-identical to the two-call form; these tests pin it against an INDEPENDENT
// linear oracle (never the bucket table) and cross-check the table-built path
// against the table-NULL path on the same bns.

namespace {

/* bns_pos2rid's contract expressed via the linear reference: out-of-range -> -1,
 * negative clamps to rid 0 (which expected_rid already returns). */
int oracle_pos2rid(const FakeBns &fake, int64_t pos_f) {
    if (pos_f >= fake.bns.l_pac) return -1;
    return fake.expected_rid(pos_f);
}

/* Independent bns_intv2rid oracle: mirrors the documented contract using only
 * bns_depos (a pure coordinate flip) and the linear expected_rid, so it shares
 * no code path -- and hence no bug -- with the bucket fast path under test. */
int oracle_intv2rid(const FakeBns &fake, int64_t rb, int64_t re) {
    const int64_t l_pac = fake.bns.l_pac;
    if (rb < l_pac && re > l_pac) return -2; // straddles the forward/reverse split
    int is_rev;
    const int rid_b = oracle_pos2rid(fake, bns_depos(&fake.bns, rb, &is_rev));
    if (rb >= re) return rid_b;              // degenerate empty interval
    const int rid_e = oracle_pos2rid(fake, bns_depos(&fake.bns, re - 1, &is_rev));
    return rid_b == rid_e ? rid_b : -1;
}

/* rb anchors worth probing: a capped set of contig boundaries and bucket
 * boundaries (and their neighbours), the genome ends, and the reverse-strand
 * mirror of each so both the forward and the depos'd reverse path are hit. The
 * cap keeps the sweep under the ~100 ms unit budget even on the dense layout. */
std::vector<int64_t> interval_anchors(const FakeBns &fake) {
    const int64_t l_pac = fake.bns.l_pac;
    std::vector<int64_t> a;
    const int n = fake.bns.n_seqs;
    const int step = n > 64 ? n / 64 : 1; // sample ~64 contigs regardless of density
    for (int i = 0; i < n; i += step) {
        const int64_t s = fake.anns[i].offset;
        for (int64_t p : {s - 1, s, s + 1}) if (p >= 0 && p < l_pac) a.push_back(p);
    }
    for (int64_t b = 0; b * kBucketWidth < l_pac && b < 64; ++b) {
        const int64_t p = b * kBucketWidth;
        for (int64_t q : {p - 1, p, p + 1}) if (q >= 0 && q < l_pac) a.push_back(q);
    }
    a.push_back(0);
    a.push_back(l_pac - 1);
    const size_t fwd = a.size();
    for (size_t i = 0; i < fwd; ++i) a.push_back(2 * l_pac - 1 - a[i]); // reverse mirror
    return a;
}

/* For every (rb, re) in the anchor x delta grid, assert the table-NULL fallback
 * and the table-built fast path both equal the independent oracle AND each
 * other -- the operative byte-identity bar, checked directly instead of via
 * end-to-end alignment output. */
void check_intervals(FakeBns &fake) {
    const int64_t l_pac = fake.bns.l_pac;
    const std::vector<int64_t> anchors = interval_anchors(fake);
    const int64_t deltas[] = {0, 1, 2, 7, kBucketWidth - 1, kBucketWidth,
                              kBucketWidth + 1, 2 * kBucketWidth, l_pac};

    std::vector<std::pair<int64_t, int64_t>> intervals;
    std::vector<int> expected;
    for (int64_t rb : anchors) {
        for (int64_t d : deltas) {
            const int64_t re = rb + d;    // d >= 0, so rb <= re (impl asserts this)
            if (re > 2 * l_pac) continue; // stay within the packed range
            intervals.emplace_back(rb, re);
            expected.push_back(oracle_intv2rid(fake, rb, re));
        }
    }

    REQUIRE(fake.bns.pos2rid_bucket == NULL); // FakeBns starts with no table
    std::vector<int> via_fallback;
    via_fallback.reserve(intervals.size());
    for (size_t i = 0; i < intervals.size(); ++i) {
        CAPTURE(intervals[i].first);
        CAPTURE(intervals[i].second);
        const int got = bns_intv2rid(&fake.bns, intervals[i].first, intervals[i].second);
        CHECK(got == expected[i]);
        via_fallback.push_back(got);
    }

    bns_build_pos2rid(&fake.bns);
    REQUIRE(fake.bns.pos2rid_bucket != NULL);
    for (size_t i = 0; i < intervals.size(); ++i) {
        CAPTURE(intervals[i].first);
        CAPTURE(intervals[i].second);
        const int got = bns_intv2rid(&fake.bns, intervals[i].first, intervals[i].second);
        CHECK(got == expected[i]);       // fast path == independent oracle
        CHECK(got == via_fallback[i]);   // fast path == original two-call form
    }
}

} // namespace

TEST_CASE("bns_intv2rid: fast path matches the oracle across layouts"
          * doctest::test_suite("unit/bns_pos2rid")) {
    // Each layout only changes the FakeBns contig geometry before running the
    // same anchor x delta sweep in check_intervals(), so they are SUBCASEs of a
    // single parameterized case rather than separate TEST_CASEs.
    SUBCASE("sparse layout (empty-bracket fast path fires)") {
        // Contigs much wider than a bucket: most buckets have an empty bracket,
        // so the single-bucket short-circuit is the path actually exercised.
        FakeBns fake(uniform_contigs(64, kBucketWidth * 7 + 123));
        check_intervals(fake);
    }

    SUBCASE("dense layout (many contig starts per bucket)") {
        FakeBns fake(uniform_contigs(4096, 8));
        check_intervals(fake);
    }

    SUBCASE("contig starts on bucket boundaries") {
        FakeBns fake(uniform_contigs(32, kBucketWidth));
        check_intervals(fake);
    }

    SUBCASE("ragged contig widths straddling bucket boundaries") {
        std::vector<int64_t> lengths;
        for (int i = 0; i < 200; ++i)
            lengths.push_back(1 + (int64_t)((i * 2654435761u) % (kBucketWidth * 2)));
        FakeBns fake(lengths);
        check_intervals(fake);
    }

    SUBCASE("single contig") {
        FakeBns fake(uniform_contigs(1, kBucketWidth * 3));
        check_intervals(fake);
    }
}

TEST_CASE("bns_intv2rid: named edge cases on a small multi-contig layout"
          * doctest::test_suite("unit/bns_pos2rid")) {
    // rid0 [0,100), rid1 [100,250), rid2 [250,300); l_pac = 300.
    FakeBns fake(std::vector<int64_t>{100, 150, 50});
    const int64_t l_pac = fake.bns.l_pac;
    REQUIRE(l_pac == 300);
    bns_build_pos2rid(&fake.bns);

    CHECK(bns_intv2rid(&fake.bns, 10, 50) == 0);   // wholly inside one contig
    CHECK(bns_intv2rid(&fake.bns, 120, 200) == 1);
    CHECK(bns_intv2rid(&fake.bns, 90, 110) == -1); // spans a contig boundary
    CHECK(bns_intv2rid(&fake.bns, 120, 120) == 1); // degenerate empty interval
    CHECK(bns_intv2rid(&fake.bns, 100, 101) == 1); // single base at a contig start
    CHECK(bns_intv2rid(&fake.bns, 50, l_pac + 10) == -2); // straddles f/r split
    // Reverse-strand interval inside rid2's forward span [250,300).
    CHECK(bns_intv2rid(&fake.bns, 2 * l_pac - 300, 2 * l_pac - 250) == 2);
}
