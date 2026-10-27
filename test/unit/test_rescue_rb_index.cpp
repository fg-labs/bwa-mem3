// test/unit/test_rescue_rb_index.cpp — rescue_rb_index (src/rescue_rb_index.h) and
// mem_rescue_pushed_rb_t.
//
// Once a mem_matesw_batch_post call has changed a mate list, mem_pair_resolve_batch_post answers
// the later anchors' skip sets from a sorted copy of the list's rb values. After each call it adds
// the rb values the call pushed, and it discards the copy when the list lost a region
// (n_after != n_before + pushed_rb.n). These tests drive rescue_rb_index through synthetic call
// sequences. The independent oracle is the per-region mem_infer_dir loop that
// mem_matesw_batch_pre / _post ran before the interval form (the same oracle as
// test_matesw_skip.cpp). The sequences cover:
//   - pushes of 0 to 4 regions, at the ends and between existing values, with duplicate rb values;
//   - removals, which must discard the copy;
//   - a removal masked by a push (n unchanged, so not a count match), which must also discard it,
//     and whose stale copy would have answered wrongly.
// The last cases run mem_matesw_batch_post itself on a one-contig fixture. They check the rb values
// it reports in mem_rescue_pushed_rb_t against the regions it inserted:
//   - one rescue, and two rescues in one call;
//   - existing regions outside the skip windows, which must be kept;
//   - an existing region outside the window that the rescued region makes redundant, which the
//     call's dedup drops, so the count check must discard the copy.
// Fixtures are built in memory.
#include "doctest/doctest.h"

extern "C" {
#include "bntseq.h"
}
#include "bwamem.h"
#include "rescue_rb_index.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <random>
#include <vector>

#define _set_pac(pac, l, c) ((pac)[(l) >> 2] |= (c) << ((~(l) & 3) << 1))

namespace {

// The orientation check as mem_matesw_batch_pre / _post wrote it before mem_matesw_skip.
unsigned oracle_bits(int64_t l_pac, int64_t rb, const std::vector<mem_alnreg_t> &regs,
                     const mem_pestat_t pes[4])
{
    unsigned bits = 0;
    for (int r = 0; r < 4; ++r) bits |= (pes[r].failed ? 1u : 0u) << r;
    for (const mem_alnreg_t &m : regs) {
        int64_t dist;
        const int r = mem_infer_dir(l_pac, rb, m.rb, &dist);
        if (dist >= pes[r].low && dist <= pes[r].high) bits |= 1u << r;
    }
    return bits;
}

mem_alnreg_t region_at(int64_t rb)
{
    mem_alnreg_t a;
    std::memset(&a, 0, sizeof(a));
    a.rb = rb;
    a.re = rb + 100;
    return a;
}

void random_pes(std::mt19937_64 &rng, mem_pestat_t pes[4])
{
    std::memset(pes, 0, 4 * sizeof(*pes));
    for (int r = 0; r < 4; ++r) {
        pes[r].failed = (rng() % 5 == 0);
        pes[r].low = (int)(rng() % 400) - 50;
        pes[r].high = pes[r].low + (int)(rng() % 900);
    }
}

unsigned probe_at(const rescue_rb_index &idx, int64_t l_pac, int64_t rb, const mem_pestat_t pes[4])
{
    mem_matesw_skip_query_t q;
    mem_matesw_skip_query_init(l_pac, rb, pes, &q);
    return idx.probe(&q);
}

// Every probe of the live copy equals the oracle over the list, for anchors spread over both
// strands and anchors at existing rb values; returns the number of mismatches.
int probe_mismatches(std::mt19937_64 &rng, int64_t l_pac, const rescue_rb_index &idx,
                     const std::vector<mem_alnreg_t> &regs, const mem_pestat_t pes[4])
{
    int bad = 0;
    for (int t = 0; t < 24; ++t) {
        int64_t rb = (int64_t)(rng() % (uint64_t)(2 * l_pac));
        if (t % 3 == 0 && !regs.empty()) rb = regs[rng() % regs.size()].rb + (int64_t)(rng() % 7) - 3;
        if (rb < 0) rb = 0;
        bad += probe_at(idx, l_pac, rb, pes) != oracle_bits(l_pac, rb, regs, pes);
    }
    return bad;
}

std::vector<int64_t> sorted_rbs(const mem_alnreg_t *a, size_t n)
{
    std::vector<int64_t> v(n);
    for (size_t k = 0; k < n; ++k) v[k] = a[k].rb;
    std::sort(v.begin(), v.end());
    return v;
}

}  // namespace

TEST_CASE("rescue_rb_index: build, probe and pushed inserts match a scan of the list"
          * doctest::test_suite("unit/rescue_rb_index"))
{
    std::mt19937_64 rng(20261025);
    int bad = 0, builds = 0, updates = 0;
    for (int trial = 0; trial < 200; ++trial) {
        const int64_t l_pac = trial % 4 == 0 ? 3000 : 1000000;
        mem_pestat_t pes[4];
        random_pes(rng, pes);
        std::vector<mem_alnreg_t> regs;
        const size_t n0 = rng() % 80;
        for (size_t k = 0; k < n0; ++k) regs.push_back(region_at((int64_t)(rng() % (uint64_t)(2 * l_pac))));
        rescue_rb_index idx;
        idx.build(regs.data(), regs.size());
        ++builds;
        REQUIRE(idx.live);
        CHECK(idx.matches(regs.data(), regs.size()));
        bad += probe_mismatches(rng, l_pac, idx, regs, pes);
        for (int call = 0; call < 12; ++call) {
            // One call: push 0 to 4 regions anywhere in the list (some at an existing rb value).
            const size_t n_before = regs.size();
            mem_rescue_pushed_rb_t pushed_rb;
            pushed_rb.n = (int)(rng() % 5);
            for (int k = 0; k < pushed_rb.n; ++k) {
                int64_t rb = (int64_t)(rng() % (uint64_t)(2 * l_pac));
                if (k == 1 && !regs.empty()) rb = regs[rng() % regs.size()].rb;  // duplicate value
                pushed_rb.rb[k] = rb;
                regs.insert(regs.begin() + (long)(rng() % (regs.size() + 1)), region_at(rb));
            }
            CHECK(idx.update(n_before, regs.size(), pushed_rb));
            ++updates;
            REQUIRE(idx.live);
            CHECK(idx.matches(regs.data(), regs.size()));
            bad += probe_mismatches(rng, l_pac, idx, regs, pes);
        }
    }
    CAPTURE(builds);
    CAPTURE(updates);
    CHECK(bad == 0);
}

TEST_CASE("rescue_rb_index: a call that drops a region discards the copy"
          * doctest::test_suite("unit/rescue_rb_index"))
{
    const int64_t l_pac = 100000;
    std::vector<mem_alnreg_t> regs;
    for (int64_t rb : {500, 1200, 1200, 40000, 150000}) regs.push_back(region_at(rb));
    mem_rescue_pushed_rb_t none;
    none.n = 0;

    SUBCASE("a pure removal")
    {
        rescue_rb_index idx;
        idx.build(regs.data(), regs.size());
        const size_t n_before = regs.size();
        regs.erase(regs.begin() + 1);
        CHECK_FALSE(idx.update(n_before, regs.size(), none));
        CHECK_FALSE(idx.live);
        // A discarded copy stays discarded, even across a later call that changes nothing.
        CHECK_FALSE(idx.update(regs.size(), regs.size(), none));
        CHECK_FALSE(idx.live);
    }
    SUBCASE("a push and a removal in the same call (same n, not a count match)")
    {
        rescue_rb_index idx;
        idx.build(regs.data(), regs.size());
        const std::vector<mem_alnreg_t> before = regs;
        const size_t n_before = regs.size();
        mem_rescue_pushed_rb_t pushed_rb;
        pushed_rb.n = 1;
        pushed_rb.rb[0] = 70000;
        regs.push_back(region_at(70000));
        regs.erase(regs.begin());  // the dedup dropped the old region at rb 500
        REQUIRE(regs.size() == n_before);
        CHECK_FALSE(idx.update(n_before, regs.size(), pushed_rb));
        CHECK_FALSE(idx.live);

        // The copy the count check prevents: the pre-call list plus the pushed values, with the
        // dropped region still in it. At rb 520 it still sees rb 500 (orientation 3 at distance
        // 20), which the list no longer holds, so it would answer a different skip set.
        rescue_rb_index stale;
        stale.build(before.data(), before.size());
        REQUIRE(stale.update(n_before, n_before + 1, pushed_rb));
        mem_pestat_t pes[4];
        std::memset(pes, 0, sizeof(pes));
        for (int r = 0; r < 4; ++r) pes[r].low = 0, pes[r].high = 50;
        const unsigned want = oracle_bits(l_pac, 520, regs, pes);
        CHECK((want & 8u) == 0);
        CHECK(probe_at(stale, l_pac, 520, pes) != want);
    }
    SUBCASE("a call that pushed nothing and dropped nothing keeps the copy")
    {
        rescue_rb_index idx;
        idx.build(regs.data(), regs.size());
        CHECK(idx.update(regs.size(), regs.size(), none));
        CHECK(idx.live);
        CHECK(idx.matches(regs.data(), regs.size()));
    }
}

namespace {

// A one-contig reference of L random bases with its reverse complement in the second half, as
// mem_matesw_batch_post's scalar path (gar == -1) reads it, and the mate's bases ms.
struct PostFixture {
    static constexpr int64_t L = 1024;
    static constexpr int l_ms = 100;
    std::vector<uint8_t> ref, pac;
    bntann1_t anns[1];
    bntseq_t bns;
    mem_opt_t *opt = nullptr;
    mem_pestat_t pes[4];
    mem_cache mmc;
    kswr_t *aln_storage = nullptr;

    PostFixture()
    {
        ref.resize(2 * L);
        std::mt19937_64 rng(7);
        for (int64_t i = 0; i < L; ++i) ref[i] = (uint8_t)(rng() & 3);
        std::memset(anns, 0, sizeof(anns));
        anns[0].len = (int32_t)L;
        anns[0].name = const_cast<char *>("chr_test");
        anns[0].anno = const_cast<char *>("");
        std::memset(&bns, 0, sizeof(bns));
        bns.l_pac = L;
        bns.n_seqs = 1;
        bns.anns = anns;
        opt = mem_opt_init();
        std::memset(pes, 0, sizeof(pes));
        for (int r = 0; r < 4; ++r) {
            pes[r].failed = 1;
            pes[r].low = 50, pes[r].high = 250, pes[r].avg = 150.0, pes[r].std = 30.0;
        }
        std::memset(&mmc, 0, sizeof(mmc));
        finish();
    }
    ~PostFixture() { free(opt); }

    // Copy ref[from, from + n) to ref[to, ...), then rebuild the reverse half and the pac.
    void copy_block(int64_t from, int64_t to, int64_t n)
    {
        for (int64_t i = 0; i < n; ++i) ref[to + i] = ref[from + i];
        finish();
    }
    void finish()
    {
        for (int64_t i = 0; i < L; ++i) ref[2 * L - 1 - i] = (uint8_t)(3 ^ ref[i]);
        pac.assign((size_t)((L + 3) / 4 + 1), 0);
        for (int64_t i = 0; i < L; ++i) _set_pac(pac.data(), i, ref[i]);
        mmc.ref_string = ref.data();
    }
    std::vector<uint8_t> mate_at(int64_t pos) const
    {
        return std::vector<uint8_t>(ref.begin() + pos, ref.begin() + pos + l_ms);
    }
    static mem_alnreg_t anchor_at(int64_t rb)
    {
        mem_alnreg_t a;
        std::memset(&a, 0, sizeof(a));
        a.rb = rb, a.re = rb + 100, a.qe = 100, a.score = 100, a.secondary = -1, a.secondary_all = -1;
        return a;
    }
    // A forward-strand mate region already in the list.
    static mem_alnreg_t mate_region(int64_t rb, int score)
    {
        mem_alnreg_t m;
        std::memset(&m, 0, sizeof(m));
        m.rb = rb, m.re = rb + 100, m.qb = 0, m.qe = 100, m.score = score, m.secondary = -1;
        m.secondary_all = -1;
        return m;
    }
    int post(const mem_alnreg_t &a, const std::vector<uint8_t> &ms, mem_alnreg_v *ma,
             mem_rescue_pushed_rb_t *pushed_rb)
    {
        int32_t gar[4] = {-1, -1, -1, -1};
        return mem_matesw_batch_post(opt, &bns, pac.data(), pes, &a, l_ms, ms.data(), ma, &aln_storage, 0,
                                     gar, &mmc, 0, nullptr, nullptr, -1, nullptr, -1, pushed_rb);
    }
};

mem_alnreg_v vec_of(const std::vector<mem_alnreg_t> &regs)
{
    mem_alnreg_v ma;
    std::memset(&ma, 0, sizeof(ma));
    for (const mem_alnreg_t &m : regs) kv_push(mem_alnreg_t, ma, m);
    return ma;
}

}  // namespace

TEST_CASE("mem_matesw_batch_post reports the rb of each region it pushes (mem_rescue_pushed_rb_t)"
          * doctest::test_suite("unit/rescue_rb_index"))
{
    PostFixture f;
    REQUIRE(f.opt != nullptr);

    SUBCASE("one rescued region")
    {
        // FF (orientation 0) only: the anchor at 200 searches [250, 550); the mate is ref[250, 350).
        f.pes[0].failed = 0;
        mem_alnreg_v ma = vec_of({});
        mem_rescue_pushed_rb_t pushed_rb;
        pushed_rb.n = -7;
        CHECK(f.post(PostFixture::anchor_at(200), f.mate_at(250), &ma, &pushed_rb) > 0);
        REQUIRE(ma.n == 1);
        REQUIRE(pushed_rb.n == 1);
        CHECK(pushed_rb.rb[0] == ma.a[0].rb);
        CHECK(ma.a[0].rb == 250);
        free(ma.a);
    }
    SUBCASE("two orientations rescue in one call, next to regions outside the skip windows")
    {
        // The anchor at 500: FF searches [550, 850), RR [250, 550). The mate is ref[300, 400),
        // copied to ref[600, 700), so both orientations find it. The existing regions (at 10, and
        // on the reverse strand at 1500) are outside every orientation's window, so neither skips
        // an orientation and neither is redundant with a rescued region.
        f.copy_block(300, 600, PostFixture::l_ms);
        f.pes[0].failed = 0;
        f.pes[3].failed = 0;
        const std::vector<mem_alnreg_t> old = {PostFixture::mate_region(10, 60),
                                               PostFixture::mate_region(1500, 60)};
        mem_alnreg_v ma = vec_of(old);
        const std::vector<int64_t> old_rbs = sorted_rbs(ma.a, ma.n);
        const size_t n_before = ma.n;
        mem_rescue_pushed_rb_t pushed_rb;
        pushed_rb.n = -7;
        CHECK(f.post(PostFixture::anchor_at(500), f.mate_at(300), &ma, &pushed_rb) > 0);
        REQUIRE(pushed_rb.n == 2);
        CHECK(ma.n == n_before + 2);
        std::vector<int64_t> want = old_rbs;
        for (int k = 0; k < pushed_rb.n; ++k) want.push_back(pushed_rb.rb[k]);
        std::sort(want.begin(), want.end());
        CHECK(sorted_rbs(ma.a, ma.n) == want);
        std::vector<int64_t> got_pushed(pushed_rb.rb, pushed_rb.rb + 2);
        std::sort(got_pushed.begin(), got_pushed.end());
        CHECK(got_pushed == std::vector<int64_t>{300, 600});

        // The copy the caller keeps follows: built before the call, updated after it.
        rescue_rb_index idx;
        idx.build(old.data(), old.size());
        CHECK(idx.update(n_before, ma.n, pushed_rb));
        CHECK(idx.matches(ma.a, ma.n));
        free(ma.a);
    }
    SUBCASE("a rescued region makes an existing one redundant: the dedup drops it, the copy is discarded")
    {
        // FF from the anchor at 200 rescues the mate at [250, 350). An existing region at
        // [248, 348) lies at distance 48, below pes[0].low = 50, so it does not skip FF; it overlaps
        // the rescued region on 98 of 100 reference bases and on the whole query (above
        // mask_level_redun), and scores lower, so the call's dedup drops it.
        f.pes[0].failed = 0;
        const std::vector<mem_alnreg_t> old = {PostFixture::mate_region(248, 40)};
        mem_alnreg_v ma = vec_of(old);
        const size_t n_before = ma.n;
        mem_rescue_pushed_rb_t pushed_rb;
        pushed_rb.n = -7;
        CHECK(f.post(PostFixture::anchor_at(200), f.mate_at(250), &ma, &pushed_rb) > 0);
        REQUIRE(pushed_rb.n == 1);
        CHECK(pushed_rb.rb[0] == 250);
        CHECK(ma.n != n_before + (size_t)pushed_rb.n);
        REQUIRE(ma.n == 1);
        CHECK(ma.a[0].rb == 250);

        rescue_rb_index idx;
        idx.build(old.data(), old.size());
        CHECK_FALSE(idx.update(n_before, ma.n, pushed_rb));
        CHECK_FALSE(idx.live);
        free(ma.a);
    }
    SUBCASE("every orientation failed: the skip-all return reports nothing")
    {
        mem_alnreg_v ma = vec_of({});
        mem_rescue_pushed_rb_t pushed_rb;
        pushed_rb.n = -7;
        CHECK(f.post(PostFixture::anchor_at(200), f.mate_at(250), &ma, &pushed_rb) == 0);
        CHECK(pushed_rb.n == 0);
        CHECK(ma.n == 0);
        free(ma.a);
    }
}
