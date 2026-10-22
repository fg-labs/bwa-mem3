// Unit tests for the known-contig window clamp: bns_pos_in_rid and the
// bns_fetch_seq_v2_rid / bns_fetch_bounds_rid variants.
//
// The extension and mate-rescue fetches pass the contig they expect as a hint,
// and the clamp tests it with bns_pos_in_rid instead of calling bns_pos2rid.
// That is only byte-identical if bns_pos_in_rid(pos_f, rid) is exactly
// `bns_pos2rid(pos_f) == rid`, and if the hinted fetches return what the
// unhinted ones do for ANY hint -- right, wrong, or out of range. These tests
// pin both, on layouts with zero-length contigs, gaps between contigs (offsets
// from an .ann are not required to be back to back) and the last contig.

#include "doctest/doctest.h"

#include "bntseq.h"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <vector>

namespace {

/* What the clamp reads: l_pac, n_seqs, anns[].offset/len, plus the pos2rid
 * table. Contig i starts gaps[i] bases after the end of contig i - 1, so a
 * non-zero gap leaves bases that belong to no contig's [offset, offset + len)
 * but that bns_pos2rid still assigns to the contig before them. */
struct LayoutBns {
    bntseq_t bns;
    std::vector<bntann1_t> anns;

    LayoutBns(const std::vector<int64_t> &lens, const std::vector<int64_t> &gaps, int64_t tail) {
        std::memset(&bns, 0, sizeof(bns));
        anns.assign(lens.size(), bntann1_t());
        int64_t off = 0;
        for (size_t i = 0; i < lens.size(); ++i) {
            off += gaps[i];
            std::memset(&anns[i], 0, sizeof(anns[i]));
            anns[i].offset = off;
            anns[i].len = (int32_t)lens[i];
            off += lens[i];
        }
        bns.l_pac = off + tail;
        bns.n_seqs = (int32_t)lens.size();
        bns.anns = anns.data();
        bns.pos2rid_bucket = NULL;
    }
    ~LayoutBns() { free(bns.pos2rid_bucket); }
    LayoutBns(const LayoutBns &) = delete;
    LayoutBns &operator=(const LayoutBns &) = delete;

    /* Every contig boundary and its neighbours, both genome ends, and a stride sweep. */
    std::vector<int64_t> positions() const {
        std::vector<int64_t> p;
        for (const bntann1_t &a : anns)
            for (int64_t d = -2; d <= 2; ++d) {
                p.push_back(a.offset + d);
                p.push_back(a.offset + a.len + d);
            }
        for (int64_t d = -2; d <= 2; ++d) p.push_back(bns.l_pac + d);
        for (int64_t x = 0; x < bns.l_pac; x += bns.l_pac / 61 + 1) p.push_back(x);
        // bns_pos_in_rid requires pos_f >= 0 (a bns_depos result always is).
        p.erase(std::remove_if(p.begin(), p.end(), [](int64_t x) { return x < 0; }), p.end());
        std::sort(p.begin(), p.end());
        p.erase(std::unique(p.begin(), p.end()), p.end());
        return p;
    }
};

/* lens, gaps and tail for the layouts under test. */
struct Layout { const char *name; std::vector<int64_t> lens, gaps; int64_t tail; };

std::vector<Layout> layouts() {
    const int64_t W = (int64_t)1 << BNS_POS2RID_SHIFT;
    return {
        {"back to back", {100, 150, 50}, {0, 0, 0}, 0},
        {"zero-length contigs, first, middle and last",
         {0, 40, 0, 0, 25, 60, 0}, {0, 0, 0, 0, 0, 0, 0}, 0},
        {"gaps between contigs and after the last", {30, 20, 45, 10}, {0, 7, 1, 13}, 9},
        {"gap before the first contig", {30, 20}, {6, 2}, 0},
        {"single contig", {500}, {0}, 0},
        {"contigs wider than a pos2rid bucket", {3 * W + 5, W, 2 * W - 1}, {0, 3, 0}, 0},
        {"many contig starts per bucket", std::vector<int64_t>(40, 9), std::vector<int64_t>(40, 1), 4},
    };
}

void for_each_layout_and_table(void (*check)(LayoutBns &)) {
    for (const Layout &l : layouts()) {
        LayoutBns fake(l.lens, l.gaps, l.tail);
        CAPTURE(l.name);
        check(fake);                      // binary-search fallback
        bns_build_pos2rid(&fake.bns);
        REQUIRE(fake.bns.pos2rid_bucket != NULL);
        check(fake);                      // bucket table
    }
}

/* Independent of bns_pos2rid: scan every contig and keep the last one starting
 * at or before pos_f (rid 0 when none does); -1 at or past l_pac. */
int linear_rid(const LayoutBns &fake, int64_t pos_f) {
    if (pos_f >= fake.bns.l_pac) return -1;
    int rid = 0;
    for (int i = 0; i < fake.bns.n_seqs; ++i)
        if (fake.anns[i].offset <= pos_f) rid = i;
    return rid;
}

void check_predicate(LayoutBns &fake) {
    int n_true = 0;
    for (int64_t pos : fake.positions()) {
        const int want = bns_pos2rid(&fake.bns, pos);
        CAPTURE(pos);
        REQUIRE(want == linear_rid(fake, pos));  // the reference itself agrees with the oracle
        // bns_pos_in_rid reads only anns[rid] and anns[rid + 1], so the rids
        // around the answer and at both ends of anns[] cover every outcome.
        const int n = fake.bns.n_seqs;
        for (int rid : {0, 1, want - 2, want - 1, want, want + 1, want + 2, n - 2, n - 1}) {
            if (rid < 0 || rid >= n) continue;
            CAPTURE(rid);
            const int got = bns_pos_in_rid(&fake.bns, pos, rid);
            CHECK(got == (want == rid));
            n_true += got;
        }
    }
    CHECK(n_true > 0);
}

/* Windows of several widths centred on each probe position, on both strands,
 * with the clamp run unhinted and with hints that are right, adjacent, at the
 * ends of anns[] and out of range. */
void check_bounds(LayoutBns &fake) {
    const int64_t l_pac = fake.bns.l_pac, L2 = l_pac << 1;
    const int64_t halves[] = {0, 1, 5, 37, 1000};
    int n_hint_hit = 0;
    for (int64_t f : fake.positions()) {
        if (f >= l_pac) continue;
        for (int64_t mid : {f, L2 - 1 - f}) {
            for (int64_t h : halves) {
                const int64_t b0 = mid - h < 0 ? 0 : mid - h;
                const int64_t e0 = mid + h + 1 > L2 ? L2 : mid + h + 1;
                int64_t rb = b0, re = e0;
                int rid = -9;
                bns_fetch_bounds(&fake.bns, &rb, mid, &re, &rid);
                const int n = fake.bns.n_seqs;
                for (int hint : {-2, -1, 0, 1, rid - 1, rid, rid + 1, n - 1, n, n + 1}) {
                    CAPTURE(mid);
                    CAPTURE(h);
                    CAPTURE(hint);
                    // Reversed order exercises the clamp's swap too.
                    int64_t hb = (hint & 1) ? e0 : b0, he = (hint & 1) ? b0 : e0;
                    int hrid = -9;
                    bns_fetch_bounds_rid(&fake.bns, &hb, mid, &he, hint, &hrid);
                    CHECK(hb == rb);
                    CHECK(he == re);
                    CHECK(hrid == rid);
                    n_hint_hit += (hint == rid);
                }
            }
        }
    }
    CHECK(n_hint_hit > 0);
}

} // namespace

TEST_CASE("bns_pos_in_rid: equals bns_pos2rid(pos) == rid"
          * doctest::test_suite("unit/bns_pos2rid")) {
    for_each_layout_and_table(check_predicate);
}

TEST_CASE("bns_fetch_bounds_rid: same bounds and rid as bns_fetch_bounds for any hint"
          * doctest::test_suite("unit/bns_pos2rid")) {
    for_each_layout_and_table(check_bounds);
}

TEST_CASE("bns_fetch_seq_v2_rid: same window as bns_fetch_seq_v2, .0123 and pac-fetch"
          * doctest::test_suite("unit/bns_pos2rid")) {
    // Gapped layout with a zero-length contig; every base gets a distinct-ish
    // value so a shifted window would not compare equal.
    LayoutBns fake({40, 0, 25, 31}, {0, 3, 0, 2}, 5);
    bns_build_pos2rid(&fake.bns);
    const int64_t l_pac = fake.bns.l_pac, L2 = l_pac << 1;
    std::vector<uint8_t> pac((size_t)((l_pac + 3) >> 2), 0);
    for (int64_t i = 0; i < l_pac; ++i) _set_pac(pac.data(), i, (i * 7 + i / 5) & 3);
    std::vector<uint8_t> ref_string((size_t)L2);   // the .0123 layout: forward then reverse complement
    for (int64_t i = 0; i < l_pac; ++i) {
        ref_string[(size_t)i] = _get_pac(pac.data(), i);
        ref_string[(size_t)(L2 - 1 - i)] = 3 - _get_pac(pac.data(), i);
    }

    for (int64_t mid = 0; mid < L2; ++mid) {
        // A window that stays on mid's strand, so the clamp never bridges.
        const int64_t lo = mid < l_pac ? 0 : l_pac, hi = mid < l_pac ? l_pac : L2;
        const int64_t b0 = mid - 6 < lo ? lo : mid - 6, e0 = mid + 9 > hi ? hi : mid + 9;
        int rid = -9;
        int64_t rb = b0, re = e0;
        uint8_t *want = bns_fetch_seq_v2(&fake.bns, pac.data(), &rb, mid, &re, &rid,
                                         ref_string.data(), NULL);
        std::vector<uint8_t> want_bytes(want, want + (re - rb));
        for (int hint = -1; hint <= fake.bns.n_seqs; ++hint) {
            CAPTURE(mid);
            CAPTURE(hint);
            int hrid = -9;
            int64_t hb = b0, he = e0;
            uint8_t *got = bns_fetch_seq_v2_rid(&fake.bns, pac.data(), &hb, mid, &he, hint, &hrid,
                                                ref_string.data(), NULL);
            CHECK(got == want);           // .0123: the same pointer into ref_string
            CHECK(hb == rb);
            CHECK(he == re);
            CHECK(hrid == rid);
            hb = b0, he = e0, hrid = -9;
            got = bns_fetch_seq_v2_rid(&fake.bns, pac.data(), &hb, mid, &he, hint, &hrid, NULL, NULL);
            REQUIRE(got != NULL);         // pac-fetch: the same bases, rebuilt from pac
            CHECK(hb == rb);
            CHECK(he == re);
            CHECK(hrid == rid);
            CHECK(std::vector<uint8_t>(got, got + (he - hb)) == want_bytes);
        }
    }
}
