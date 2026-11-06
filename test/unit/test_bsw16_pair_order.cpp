// test/unit/test_bsw16_pair_order.cpp — the 16-bit kernel's per-call pair order
// (bsw16_len1_major_order, sortPairsLen16).
//
// Pinned here, on synthetic pairs:
//   - the key choice: len1-major only when strictly more than half of the pairs have
//     len1 < len2; ties, all-equal-length pairs, n = 0 and n = 1 go the stated way;
//   - the order sortPairsLen16 produces matches an independent oracle, a std::stable_sort of
//     the input on (major key, minor key), for both choices and for the diverted call, whose
//     len1 key is max(len1, len2); equal pairs keep input order (stable) and no pair is
//     dropped or duplicated (every slot's full record is checked);
//   - n = 0 leaves the array untouched;
//   - through bsw_extend_tier_call and a recording IBandedPairWiseSW, the kernel entry, band
//     and full pair order production hands the kernel: the 8-bit tier's len2 order at
//     w <= its int8 band limit (BSW_EXT_CALL_8_MAX_W), its divert to getScores16 with the
//     diverted two-key order past the limit (both majority outcomes), the 16-bit tier's
//     two-key order, and the scalar tier's input order.

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

#include "doctest/doctest.h"
#include "bandedSWA.h"
#include "bwamem.h"
#include "run_in_child.h"        // forked death test
#include "temp_dir.h"
#include "scoped_rescue_env.h"   // ScopedDedupMode

namespace {

SeqPair make_pair(int len1, int len2, int id)
{
    SeqPair sp = {};
    sp.len1 = len1; sp.len2 = len2; sp.id = id;
    sp.idr = 100 + id; sp.idq = 200 + id; sp.h0 = 1000 + id;
    sp.seqid = 2000 + id; sp.regid = 3000 + id;
    sp.score = 4000 + id; sp.tle = 5000 + id; sp.gtle = 6000 + id; sp.qle = 7000 + id;
    sp.gscore = 8000 + id; sp.max_off = 9000 + id;
    sp.tight_band = 10 + id; sp.chain_band = 20 + id; sp.ugp_walk_score = 30 + id;
    sp.ugp_r_attempted = id % 2; sp.meth_hyp = id % 3 - 1;
    return sp;
}

// The order sortPairsLen16 must produce: stable on (major, minor).
std::vector<SeqPair> oracle(std::vector<SeqPair> v, bool len1_major, bool diverted)
{
    auto l1 = [&](const SeqPair &p) { return diverted ? std::max(p.len1, p.len2) : p.len1; };
    std::stable_sort(v.begin(), v.end(), [&](const SeqPair &a, const SeqPair &b) {
        const int am = len1_major ? l1(a) : a.len2, bm = len1_major ? l1(b) : b.len2;
        if (am != bm) return am < bm;
        const int an = len1_major ? a.len2 : l1(a), bn = len1_major ? b.len2 : l1(b);
        return an < bn;
    });
    return v;
}

void check_pair(const SeqPair &got, const SeqPair &want)
{
    CHECK(got.id == want.id);
    CHECK(got.len1 == want.len1);
    CHECK(got.len2 == want.len2);
    CHECK(got.idr == want.idr);
    CHECK(got.idq == want.idq);
    CHECK(got.h0 == want.h0);
    CHECK(got.seqid == want.seqid);
    CHECK(got.regid == want.regid);
    CHECK(got.score == want.score);
    CHECK(got.tle == want.tle);
    CHECK(got.gtle == want.gtle);
    CHECK(got.qle == want.qle);
    CHECK(got.gscore == want.gscore);
    CHECK(got.max_off == want.max_off);
    CHECK(got.tight_band == want.tight_band);
    CHECK(got.chain_band == want.chain_band);
    CHECK(got.ugp_walk_score == want.ugp_walk_score);
    CHECK(got.ugp_r_attempted == want.ugp_r_attempted);
    CHECK(got.meth_hyp == want.meth_hyp);
}

void check_sorted(std::vector<SeqPair> in, bool diverted)
{
    const bool len1_major = bsw16_len1_major_order(in.data(), (int32_t) in.size());
    const std::vector<SeqPair> want = oracle(in, len1_major, diverted);
    std::vector<SeqPair> tmp(in.size() + 1);
    std::vector<int32_t> hist(MAX_SEQ_LEN16 + 1, 0x5a5a5a5a);
    sortPairsLen16(in.data(), (int32_t) in.size(), tmp.data(), hist.data(), diverted);
    REQUIRE(in.size() == want.size());
    for (size_t i = 0; i < in.size(); i++) {
        CAPTURE(i);
        check_pair(in[i], want[i]);
    }
}

}  // namespace

TEST_CASE("bsw16 pair order: the key choice follows the len1 < len2 majority"
          * doctest::test_suite("unit/bwamem"))
{
    SUBCASE("n = 0 goes len2-major") {
        CHECK_FALSE(bsw16_len1_major_order(nullptr, 0));
    }
    SUBCASE("n = 1") {
        SeqPair shorter_ref = make_pair(30, 80, 0), longer_ref = make_pair(80, 30, 0), equal = make_pair(50, 50, 0);
        CHECK(bsw16_len1_major_order(&shorter_ref, 1));
        CHECK_FALSE(bsw16_len1_major_order(&longer_ref, 1));
        CHECK_FALSE(bsw16_len1_major_order(&equal, 1));   // len1 == len2 is not len1 < len2
    }
    SUBCASE("a tie goes len2-major") {
        std::vector<SeqPair> v = {make_pair(30, 80, 0), make_pair(80, 30, 1)};
        CHECK_FALSE(bsw16_len1_major_order(v.data(), 2));
        v = {make_pair(10, 20, 0), make_pair(20, 10, 1), make_pair(15, 25, 2), make_pair(25, 25, 3)};
        CHECK_FALSE(bsw16_len1_major_order(v.data(), 4));
    }
    SUBCASE("a strict majority of len1 < len2 goes len1-major") {
        std::vector<SeqPair> v = {make_pair(10, 20, 0), make_pair(20, 10, 1), make_pair(15, 25, 2)};
        CHECK(bsw16_len1_major_order(v.data(), 3));
        v.push_back(make_pair(40, 40, 3));                 // 2 of 4: a tie again
        CHECK_FALSE(bsw16_len1_major_order(v.data(), 4));
    }
    SUBCASE("nine pairs: 5 of 9 shorter-reference goes len1-major, 4 of 9 does not") {
        // Literal expectations, independent of the selector the sort oracle reuses. Five of
        // nine pairs have len1 < len2 (a strict majority past n = 4); the other four have
        // len1 >= len2, two of them equal.
        std::vector<SeqPair> v = {make_pair(10, 60, 0), make_pair(90, 60, 1), make_pair(20, 70, 2),
                                  make_pair(60, 60, 3), make_pair(30, 80, 4), make_pair(85, 40, 5),
                                  make_pair(40, 90, 6), make_pair(70, 70, 7), make_pair(50, 95, 8)};
        CHECK(bsw16_len1_major_order(v.data(), 9));
        // Turn pair 8 into len1 == len2: 4 of 9 shorter-reference, a strict minority.
        v[8] = make_pair(95, 95, 8);
        CHECK_FALSE(bsw16_len1_major_order(v.data(), 9));
        // And 4 of 8 (drop the last pair) is a tie, len2-major.
        CHECK_FALSE(bsw16_len1_major_order(v.data(), 8));
    }
    SUBCASE("a strict majority of len1 >= len2 goes len2-major") {
        std::vector<SeqPair> v = {make_pair(300, 200, 0), make_pair(250, 250, 1), make_pair(20, 90, 2)};
        CHECK_FALSE(bsw16_len1_major_order(v.data(), 3));
    }
}

TEST_CASE("bsw16 pair order: sortPairsLen16 is the stable two-key order"
          * doctest::test_suite("unit/bwamem"))
{
    std::mt19937 rng(20261101);
    SUBCASE("n = 0 leaves the array untouched") {
        SeqPair guard = make_pair(7, 9, 42);
        std::vector<SeqPair> tmp(1);
        std::vector<int32_t> hist(MAX_SEQ_LEN16 + 1, 0);
        sortPairsLen16(&guard, 0, tmp.data(), hist.data(), false);
        check_pair(guard, make_pair(7, 9, 42));
    }
    SUBCASE("n = 1") {
        check_sorted({make_pair(30, 80, 0)}, false);
        check_sorted({make_pair(80, 30, 0)}, true);
    }
    SUBCASE("len1 < len2 majority (short-read shape), with repeated keys") {
        for (int trial = 0; trial < 20; trial++) {
            std::vector<SeqPair> v;
            const int n = 1 + (int) (rng() % 200);
            for (int k = 0; k < n; k++) {
                int len2 = 60 + (int) (rng() % 4) * 10;          // few distinct len2: ties
                int len1 = (rng() % 8 == 0) ? len2 + (int) (rng() % 20) : (int) (rng() % len2);
                v.push_back(make_pair(len1, len2, k));
            }
            check_sorted(v, false);
        }
    }
    SUBCASE("len1 >= len2 majority (long-extension shape), with repeated keys") {
        for (int trial = 0; trial < 20; trial++) {
            std::vector<SeqPair> v;
            const int n = 1 + (int) (rng() % 300);
            for (int k = 0; k < n; k++) {
                int len2 = 100 + (int) (rng() % 5) * 50;
                int len1 = (rng() % 8 == 0) ? (int) (rng() % len2) : len2 + (int) (rng() % 4) * 25;
                v.push_back(make_pair(len1, len2, k));
            }
            check_sorted(v, false);
            check_sorted(v, true);   // the diverted call's max(len1, len2) key
        }
    }
    SUBCASE("an exact tie sorts len2-major") {
        std::vector<SeqPair> v = {make_pair(90, 10, 0), make_pair(10, 90, 1), make_pair(50, 10, 2),
                                  make_pair(10, 50, 3), make_pair(50, 50, 4), make_pair(20, 40, 5)};
        // 3 of 6 have len1 < len2: not a strict majority
        REQUIRE_FALSE(bsw16_len1_major_order(v.data(), (int32_t) v.size()));
        check_sorted(v, false);
    }
}

#if HAVE_BSW_VECTOR_8_16
namespace {

// Records which kernel entry the production tier call reached, at what band, and the
// order of the pairs it was handed. Writes a content-derived score so the scatter back to
// each pair can be checked too.
class RecordingBsw final : public IBandedPairWiseSW {
public:
    enum Entry { NONE, SCALAR, GS8, GS16 };
    Entry entry = NONE;
    int calls = 0;
    int32_t w = -1;
    std::vector<int> ids;

    int scalarBandedSWA(int, const uint8_t *, int, const uint8_t *, int32_t, int, int *, int *, int *, int *,
                        int *) override { return 0; }
    void scalarBandedSWAWrapper(SeqPair *pa, uint8_t *, uint8_t *, int n, int, int32_t w_) override {
        record(SCALAR, pa, n, w_);
    }
    void getScores8(SeqPair *pa, uint8_t *, uint8_t *, int32_t n, uint16_t, int32_t w_) override {
        record(GS8, pa, n, w_);
    }
    void getScores16(SeqPair *pa, uint8_t *, uint8_t *, int32_t n, uint16_t, int32_t w_) override {
        record(GS16, pa, n, w_);
    }
    uint64_t sw_cells() const override { return 0; }
    void set_guard_overshoot(bool) override {}
    void set_meth_cpg(int, int, int) override {}
    void set_lane_compaction(int, int) override {}
    BswCompactCounts lane_compaction_counts() const override { return {}; }
    BswScoreDomainCounts bsw8_score_domain_counts() const override { return {}; }
    void set_lane_compaction16(int, int) override {}
    BswCompactCounts lane_compaction16_counts() const override { return {}; }

private:
    void record(Entry e, SeqPair *pa, int n, int32_t w_) {
        entry = e; calls++; w = w_;
        ids.clear();
        for (int i = 0; i < n; i++) {
            ids.push_back(pa[i].id);
            pa[i].score = 1000 + pa[i].len1;
            pa[i].tle = pa[i].len2; pa[i].qle = pa[i].h0;
            // h0 is distinct per pair (40 + i), so a result scattered to the wrong pair fails.
            pa[i].gscore = 2000 + pa[i].h0; pa[i].gtle = 3000 + pa[i].h0; pa[i].max_off = pa[i].h0 - 40;
        }
    }
};

struct Dims { int len1, len2; };

// Drive bsw_extend_tier_call on pairs with these dimensions (pair i: id i, a distinct h0
// and distinct bases, so dedup finds no repeat and scores every pair) and return the recording
// double. With dedup on, these batches therefore never build a repeat group: the
// copy-to-repeat path is covered by the 8-bit dispatch case in test_bandedswa_lane_indep.cpp.
RecordingBsw run_tier_call(BswExtCall tier, int w, const std::vector<Dims> &dims)
{
    const int n = (int) dims.size();
    std::vector<uint8_t> ref, qer;
    std::vector<SeqPair> pa(n + 128), scratch(n + 128);
    std::vector<int32_t> hist(MAX_SEQ_LEN16 + 1, 0);
    for (int i = 0; i < n; i++) {
        SeqPair sp = make_pair(dims[i].len1, dims[i].len2, i);
        sp.h0 = 40 + i;
        sp.seqid = 500 + i; sp.regid = 9 * i;
        sp.idr = (int) ref.size(); sp.idq = (int) qer.size();
        for (int k = 0; k < sp.len1; k++) ref.push_back((uint8_t) ((k + i) % 4));
        for (int k = 0; k < sp.len2; k++) qer.push_back((uint8_t) ((k * 3 + i) % 4));
        pa[i] = sp;
    }
    ref.resize(ref.size() + 4096, 0);
    qer.resize(qer.size() + 4096, 0);
    RecordingBsw rec;
    bsw_extend_tier_call(tier, &rec, pa.data(), ref.data(), qer.data(), n, w, scratch.data(), hist.data());
    // Every pair comes back once, all six result fields the double's for its own dimensions.
    std::vector<char> seen(n, 0);
    for (int s = 0; s < n; s++) {
        const int i = pa[s].seqid - 500;
        REQUIRE(i >= 0);
        REQUIRE(i < n);
        REQUIRE(pa[s].regid == 9 * i);
        REQUIRE_FALSE(seen[i]);
        seen[i] = 1;
        CHECK(pa[s].score == 1000 + dims[i].len1);
        CHECK(pa[s].tle == dims[i].len2);
        CHECK(pa[s].qle == 40 + i);
        CHECK(pa[s].gscore == 2040 + i);
        CHECK(pa[s].gtle == 3040 + i);
        CHECK(pa[s].max_off == i);
    }
    return rec;
}

// len1 >= len2 for every pair (the shape the 8-bit envelope admits). Sorted by len2 alone:
// 4,1,2,0,3; len2-major with len1 ties: 4,1,2,3,0; len1 alone: 1,4,3,0,2.
const std::vector<Dims> kLongRef = {{50, 30}, {40, 20}, {60, 20}, {45, 30}, {40, 10}};
// Three of five pairs have len1 < len2. len1-major (len2 ties): 2,3,0,4,1; keyed on
// max(len1, len2): 2,4,1,0,3; len2-major (len1 ties): 1,2,4,0,3.
const std::vector<Dims> kShortRef = {{30, 80}, {70, 20}, {10, 50}, {20, 90}, {60, 50}};

}  // namespace

TEST_CASE("extension tier call: the kernel, band and pair order production hands the kernel"
          * doctest::test_suite("unit/bwamem"))
{
    const int lim = BSW_EXT_CALL_8_MAX_W;
    for (const char *dedup : {"off", "on"}) {
        const std::string dedup_name(dedup);
        CAPTURE(dedup_name);
        const bwa_tests::ScopedDedupMode dedup_mode(dedup);
        SUBCASE("8-bit tier at w <= the int8 band limit: getScores8, len2 order") {
            for (int w : {100, lim}) {
                const RecordingBsw r = run_tier_call(BSW_EXT_CALL_8, w, kLongRef);
                CHECK(r.calls == 1);
                CHECK(r.entry == RecordingBsw::GS8);
                CHECK(r.w == w);
                CHECK(r.ids == std::vector<int>({4, 1, 2, 0, 3}));
            }
        }
        SUBCASE("8-bit tier past the int8 band limit: diverted getScores16, len2-major when len1 >= len2") {
            for (int w : {lim + 1, lim + 4, 200}) {
                const RecordingBsw r = run_tier_call(BSW_EXT_CALL_8, w, kLongRef);
                CHECK(r.calls == 1);
                CHECK(r.entry == RecordingBsw::GS16);
                CHECK(r.w == w);
                CHECK(r.ids == std::vector<int>({4, 1, 2, 3, 0}));
            }
        }
        SUBCASE("8-bit tier past the int8 band limit: diverted getScores16, max(len1, len2)-major on a len1 < len2 majority") {
            const RecordingBsw r = run_tier_call(BSW_EXT_CALL_8, lim + 3, kShortRef);
            CHECK(r.calls == 1);
            CHECK(r.entry == RecordingBsw::GS16);
            CHECK(r.w == lim + 3);
            CHECK(r.ids == std::vector<int>({2, 4, 1, 0, 3}));
        }
        SUBCASE("16-bit tier: len2-major when len1 >= len2") {
            const RecordingBsw r = run_tier_call(BSW_EXT_CALL_16, 100, kLongRef);
            CHECK(r.calls == 1);
            CHECK(r.entry == RecordingBsw::GS16);
            CHECK(r.w == 100);
            CHECK(r.ids == std::vector<int>({4, 1, 2, 3, 0}));
        }
        SUBCASE("16-bit tier: len1-major on a len1 < len2 majority") {
            const RecordingBsw r = run_tier_call(BSW_EXT_CALL_16, 100, kShortRef);
            CHECK(r.calls == 1);
            CHECK(r.entry == RecordingBsw::GS16);
            CHECK(r.ids == std::vector<int>({2, 3, 0, 4, 1}));
        }
        SUBCASE("scalar tier: the wrapper, input order") {
            const RecordingBsw r = run_tier_call(BSW_EXT_CALL_SCALAR, 100, kShortRef);
            CHECK(r.calls == 1);
            CHECK(r.entry == RecordingBsw::SCALAR);
            CHECK(r.ids == std::vector<int>({0, 1, 2, 3, 4}));
        }
    }
}

TEST_CASE("extension tier call: an unknown tier value is fatal, not a silent scalar call"
          * doctest::test_suite("unit/bwamem"))
{
    // bsw_extend_tier_call xasserts its tier (an xassert aborts with a message, in all builds).
    // Without the check a stray value would run the scalar wrapper, the very oracle the 8-bit
    // dispatch test compares against, and pass.
    for (int bad : {-1, 3}) {
        CAPTURE(bad);
        const bwa_tests::TempDir dir("bsw_tier_call_death");
        const std::string err = dir.path() + "/child.err";
        const int status = bwa_tests::run_in_child([bad] {
            run_tier_call(static_cast<BswExtCall>(bad), 100, kShortRef);
        }, err.c_str());
        CHECK(WIFSIGNALED(status));
        CHECK(WTERMSIG(status) == SIGABRT);
        FILE *f = fopen(err.c_str(), "r");
        REQUIRE(f != nullptr);
        char buf[512] = {0};
        const size_t got = fread(buf, 1, sizeof buf - 1, f);
        fclose(f);
        CHECK(std::string(buf, got).find("bsw_extend_tier_call: unknown tier") != std::string::npos);
    }
}
#endif  // HAVE_BSW_VECTOR_8_16
