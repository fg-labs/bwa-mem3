// test/unit/test_mem_seed_bns.cpp
//
// The reference handles `mem` gets from the index loader (bns_restore of the
// .ann/.amb files and the 2-bit .pac), checked against the FASTA the index was
// built from. The other bns unit tests run on hand-built bntseq_t structures;
// these build a real index with `bwa_idx_build` and load it back the way `mem`
// does (FMI_search::load_index), so the indexer's contig table and hole list,
// the loader's parsing, and the lookups `mem` makes on them are checked end to
// end. The oracle is the FASTA, read by the bwa conventions documented in
// src/bntseq.h and the bwa index format:
//
//   - contigs are concatenated in FASTA order (offset = sum of the previous
//     lengths); the name is the header's first word, the annotation the rest
//     (the .ann file says "(null)" for none, which the loader reads back as an
//     empty string);
//   - a hole is a maximal run of one repeated non-ACGT character (case
//     sensitive, never crossing a contig end), recorded with its offset,
//     length and that character; lower-case acgt are bases;
//   - the .pac holds every base 2-bit encoded, and some base for a hole;
//   - positions at or past l_pac are on the reverse strand, position p there
//     mirroring forward position 2*l_pac - 1 - p.

#include "doctest/doctest.h"
#include "mem_fixture.h"

#include <cstdint>
#include <cstring>
#include <memory>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "FMI_search.h"
#include "bntseq.h"

namespace {

using bwa_tests::Contig;

// Three contigs with every kind of hole: an N run inside, a contig starting
// with N, an `n` run right after an `N` run (two holes: the characters
// differ), an IUPAC R, a run ending a contig followed by a contig starting
// with the same character (two holes: a run never crosses contigs), a
// lower-case stretch (bases, not holes), and a header with an annotation.
std::vector<Contig> bns_reference() {
    std::mt19937 rng(424242u);
    auto bases = [&](size_t n) -> std::string {
        std::string s(n, 'A');
        for (size_t i = 0; i < n; ++i) s[i] = "ACGT"[rng() % 4];
        return s;
    };
    std::string a = bases(400), b = bases(300), c = bases(250);
    a.replace(100, 20, std::string(20, 'N'));
    a.replace(200, 1, "R");
    a.replace(395, 5, std::string(5, 'N'));
    b.replace(0, 7, std::string(7, 'N'));
    b.replace(150, 6, std::string(6, 'N'));
    b.replace(156, 4, std::string(4, 'n'));
    for (size_t i = 50; i < 80; ++i) c[i] = (char)(c[i] - 'A' + 'a');
    c.replace(249, 1, "N");
    return {{"ctgA", a}, {"ctgB second contig, annotated", b}, {"ctgC", c}};
}

struct Hole {
    int64_t offset;
    int32_t len;
    char amb;
};

int base_code(char ch) {
    switch (ch) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 4;
    }
}

// The holes of the FASTA, by the rule in the header.
std::vector<Hole> holes_of(const std::vector<Contig> &ref) {
    std::vector<Hole> h;
    int64_t off = 0;
    for (const Contig &c : ref) {
        for (size_t i = 0; i < c.seq.size(); ++i) {
            if (base_code(c.seq[i]) < 4) continue;
            if (i > 0 && c.seq[i - 1] == c.seq[i]) ++h.back().len;
            else h.push_back({off + (int64_t)i, 1, c.seq[i]});
        }
        off += (int64_t)c.seq.size();
    }
    return h;
}

struct BnsFixture {
    std::vector<Contig> ref;
    std::string F;                  // the contigs concatenated
    std::vector<int64_t> start;     // contig offsets, plus l_pac at the end
    std::unique_ptr<bwa_tests::TempIndex> ix;
    BnsFixture() : ref(bns_reference()) {
        start.push_back(0);
        for (const Contig &c : ref) {
            F += c.seq;
            start.push_back((int64_t)F.size());
        }
        ix.reset(new bwa_tests::TempIndex(ref));
    }
    int64_t l_pac() const { return (int64_t)F.size(); }
    // The contig holding forward position p.
    int rid_of(int64_t p) const {
        int r = 0;
        while (start[(size_t)r + 1] <= p) ++r;
        return r;
    }
    // Base at doubled-text position p, as `mem` sees it: the FASTA base, or
    // the .pac's base at a hole; complemented on the reverse strand.
    int text(int64_t p) const {
        const bool rev = p >= l_pac();
        const int64_t f = rev ? 2 * l_pac() - 1 - p : p;
        int c = base_code(F[(size_t)f]);
        if (c == 4) c = _get_pac(ix->pac(), f);
        return rev ? 3 - c : c;
    }
};

BnsFixture &bns_fixture() {
    static std::unique_ptr<BnsFixture> f;
    if (!f) f.reset(new BnsFixture);
    return *f;
}

int collect_visit(int64_t st, int64_t en, void *ctx) {
    std::vector<std::pair<int64_t, int64_t>> *v = (std::vector<std::pair<int64_t, int64_t>> *)ctx;
    v->push_back(std::make_pair(st, en));
    return 0;
}

int stop_after_first(int64_t st, int64_t en, void *ctx) {
    (void)en;
    *(int64_t *)ctx = st;
    return 7;
}

}  // namespace

TEST_CASE("bns_restore: a built index's contig table and holes and pac match its FASTA"
          * doctest::test_suite("unit/mem_seed")) {
    BnsFixture &f = bns_fixture();
    const bntseq_t *bns = f.ix->bns();
    REQUIRE(bns != nullptr);
    CHECK(bns->l_pac == f.l_pac());
    REQUIRE(bns->n_seqs == 3);
    const char *names[] = {"ctgA", "ctgB", "ctgC"};
    const char *annos[] = {"", "second contig, annotated", ""};
    const std::vector<Hole> holes = holes_of(f.ref);
    for (int r = 0; r < 3; ++r) {
        CAPTURE(r);
        CHECK(std::string(bns->anns[r].name) == names[r]);
        CHECK(std::string(bns->anns[r].anno) == annos[r]);
        CHECK(bns->anns[r].offset == f.start[(size_t)r]);
        CHECK(bns->anns[r].len == f.start[(size_t)r + 1] - f.start[(size_t)r]);
        CHECK(bns->anns[r].is_alt == 0);
        int n_ambs = 0;
        for (const Hole &h : holes) n_ambs += (h.offset >= f.start[(size_t)r] && h.offset < f.start[(size_t)r + 1]);
        CHECK(bns->anns[r].n_ambs == n_ambs);
    }
    // 20 N, R, 5 N | 7 N, 6 N, 4 n | 1 N.
    REQUIRE(holes.size() == 7);
    REQUIRE(bns->n_holes == (int32_t)holes.size());
    for (size_t i = 0; i < holes.size(); ++i) {
        CAPTURE(i);
        CHECK(bns->ambs[i].offset == holes[i].offset);
        CHECK(bns->ambs[i].len == holes[i].len);
        CHECK(bns->ambs[i].amb == holes[i].amb);
    }
    for (int64_t p = 0; p < f.l_pac(); ++p) {
        const int c = base_code(f.F[(size_t)p]);
        if (c < 4 && _get_pac(f.ix->pac(), p) != c) {
            CAPTURE(p);
            FAIL("pac base differs from the FASTA");
        }
    }
}

TEST_CASE("bns_pos2rid and bns_intv2rid: positions and intervals on both strands resolve to the FASTA's contigs"
          * doctest::test_suite("unit/mem_seed")) {
    BnsFixture &f = bns_fixture();
    const bntseq_t *bns = f.ix->bns();
    const int64_t L = f.l_pac();
    for (int64_t p = 0; p < L; ++p) {
        if (bns_pos2rid(bns, p) != f.rid_of(p)) {
            CAPTURE(p);
            CHECK(bns_pos2rid(bns, p) == f.rid_of(p));
        }
    }
    CHECK(bns_pos2rid(bns, L) == -1);
    // [rb, re) in doubled coordinates: -2 across the strand boundary, else the
    // contig of both ends (mirrored on the reverse strand) or -1 when they differ.
    int64_t n_bridge = 0, n_split = 0, n_bad = 0;
    for (int64_t rb = 0; rb < 2 * L; ++rb) {
        for (int64_t len = 1; len <= 30 && rb + len <= 2 * L; len += 7) {
            const int64_t re = rb + len;
            int want;
            if (rb < L && re > L) {
                want = -2;
                ++n_bridge;
            } else {
                const int64_t fb = rb < L ? rb : 2 * L - 1 - rb, fe = re - 1 < L ? re - 1 : 2 * L - 1 - (re - 1);
                want = f.rid_of(fb) == f.rid_of(fe) ? f.rid_of(fb) : -1;
                n_split += (want == -1);
            }
            if (bns_intv2rid(bns, rb, re) != want && n_bad++ < 5) {
                CAPTURE(rb);
                CAPTURE(re);
                CHECK(bns_intv2rid(bns, rb, re) == want);
            }
        }
    }
    CHECK(n_bad == 0);
    CHECK(n_bridge > 0);
    CHECK(n_split > 0);
}

TEST_CASE("bns_fetch_seq_into and bns_fetch_seq_v2_rid: windows clamp to mid's contig and strand and hold its bases"
          * doctest::test_suite("unit/mem_seed")) {
    BnsFixture &f = bns_fixture();
    const bntseq_t *bns = f.ix->bns();
    const int64_t L = f.l_pac();
    std::vector<uint8_t> dst((size_t)(2 * L));
    int64_t n_clamped = 0, n_windows = 0;
    for (int64_t mid = 0; mid < 2 * L; mid += 3) {
        const int64_t spans[][2] = {{0, 1}, {5, 9}, {40, 40}, {120, 7}, {500, 500}};
        for (const int64_t *sp : spans) {
            CAPTURE(mid);
            CAPTURE(sp[0]);
            CAPTURE(sp[1]);
            // Expected clamp: the contig of mid, mirrored on the reverse strand.
            const bool rev = mid >= L;
            const int rid = f.rid_of(rev ? 2 * L - 1 - mid : mid);
            int64_t lo = f.start[(size_t)rid], hi = f.start[(size_t)rid + 1];
            if (rev) {
                const int64_t t = lo;
                lo = 2 * L - hi;
                hi = 2 * L - t;
            }
            const int64_t beg0 = mid - sp[0], end0 = mid + sp[1];
            const int64_t want_beg = beg0 > lo ? beg0 : lo, want_end = end0 < hi ? end0 : hi;
            n_clamped += (want_beg != beg0 || want_end != end0);
            ++n_windows;

            int64_t beg = beg0, end = end0, len = -1;
            int got_rid = -9;
            bns_fetch_seq_into(bns, f.ix->pac(), &beg, mid, &end, &got_rid, dst.data(), &len);
            CHECK(got_rid == rid);
            REQUIRE(beg == want_beg);
            REQUIRE(end == want_end);
            REQUIRE(len == want_end - want_beg);
            bool same = true;
            for (int64_t i = 0; i < len; ++i) same = same && dst[(size_t)i] == f.text(want_beg + i);
            CHECK(same);

            // The pac-fetch v2 path (no .0123 loaded, as `mem` runs), with no
            // hint, the right hint and a wrong one.
            const int hints[] = {-1, rid, (rid + 1) % 3};
            for (int hint : hints) {
                CAPTURE(hint);
                int64_t b2 = end0, e2 = beg0;   // swapped on purpose: the clamp orders them
                int r2 = -9;
                const uint8_t *w = bns_fetch_seq_v2_rid(bns, f.ix->pac(), &b2, mid, &e2, hint, &r2, nullptr, nullptr);
                CHECK(r2 == rid);
                REQUIRE(b2 == want_beg);
                REQUIRE(e2 == want_end);
                REQUIRE(w != nullptr);
                bool same2 = true;
                for (int64_t i = 0; i < e2 - b2; ++i) same2 = same2 && w[i] == f.text(want_beg + i);
                CHECK(same2);
            }
        }
    }
    CHECK(n_clamped > 0);
    CHECK(n_clamped < n_windows);
    // A window across the strand boundary yields nothing.
    int64_t len = -1;
    bns_get_seq_into(L, f.ix->pac(), L - 3, L + 3, dst.data(), &len);
    CHECK(len == 0);
}

TEST_CASE("bns_iter_ambi: visits exactly the holes overlapping a window in order and clipped to it"
          * doctest::test_suite("unit/mem_seed")) {
    BnsFixture &f = bns_fixture();
    const bntseq_t *bns = f.ix->bns();
    const std::vector<Hole> holes = holes_of(f.ref);
    const int64_t L = f.l_pac();
    for (int64_t pos = 0; pos < L; ++pos) {
        const int lens[] = {1, 4, 25, 160};
        for (int len : lens) {
            std::vector<std::pair<int64_t, int64_t>> want, got;
            for (const Hole &h : holes) {
                const int64_t st = h.offset > pos ? h.offset : pos;
                const int64_t en = h.offset + h.len < pos + len ? h.offset + h.len : pos + len;
                if (st < en) want.push_back(std::make_pair(st, en));
            }
            const int rc = bns_iter_ambi(bns, pos, len, collect_visit, &got);
            if (got != want || rc != 0) {
                CAPTURE(pos);
                CAPTURE(len);
                CHECK(got.size() == want.size());
                CHECK(rc == 0);
                FAIL("bns_iter_ambi visits differ from the FASTA's holes");
            }
            if (want.size() > 0) {
                // A nonzero return stops the walk and is passed back.
                int64_t first = -1;
                CHECK(bns_iter_ambi(bns, pos, len, stop_after_first, &first) == 7);
                CHECK(first == want[0].first);
            }
        }
    }
}
