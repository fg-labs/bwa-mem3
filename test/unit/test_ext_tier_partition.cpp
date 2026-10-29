// test/unit/test_ext_tier_partition.cpp — the extension tier partition (sortPairsLenExt).
//
// sortPairsLenExt splits a side's extension pairs into the three kernel tiers and lays them
// out [8-bit | 16-bit | scalar], each tier in input order. The 8-bit pairs are compacted in
// place; the 16-bit and scalar pairs are staged in tempArray (from its front and its back)
// and copied back. Pinned here, on synthetic pairs whose tier is chosen by construction:
//
//   - the layout and every slot's full record match an independent oracle, a std::stable_sort
//     of the input by the tier each pair was built for (so in-tier order is input order and a
//     reversed scalar run, a dropped pair or a duplicated pair all fail);
//   - the counts match the oracle's;
//   - tempArray of exactly count slots is enough: a guard slot past it stays untouched;
//   - an all-8-bit batch neither moves a pair nor writes tempArray;
//   - count == 0 leaves both arrays untouched and returns zero counts.
//
// The 8-bit/16-bit boundary is bsw8_envelope_ok's (pinned by its own tests); each template
// REQUIREs that the gate agrees with the tier it was built for. The 16-bit/scalar bounds
// (both lengths and minval in [0, MAX_SEQ_LEN16)) are pinned from both sides.
//
// BWAMEM3_DISABLE_BSW8 (read once per process by bsw8_envelope_ok) sends every pair off the
// 8-bit path, so the templates built for it would land in 16-bit. When it is set, the cases
// that depend on the 8-bit gate are skipped with a message; the empty-side case still runs.

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <random>
#include <vector>

#include "doctest/doctest.h"
#include "bandedSWA.h"
#include "bwamem.h"

namespace {

// Production defaults (-A 1, -w 100, -d 100).
constexpr int kA = 1, kW = 100, kZdrop = 100;

enum Tier { T8 = 0, T16 = 1, T1 = 2 };

struct Template { int len1, len2, h0; Tier tier; };

// Whether the 8-bit path is disabled for this process, parsed as bsw8_envelope_ok does.
bool bsw8_disabled()
{
    const char *e = std::getenv("BWAMEM3_DISABLE_BSW8");
    return e && std::atoi(e);
}

// One or more ways into each tier at the production scoring. The subcases below refer to
// templates by index, so append new ones at the end.
const Template kTemplates[] = {
    {100, 100,  50, T8},                    // max score 150 < 254
    {150, 120,   0, T8},                    // h0 == 0 stays 8-bit
    { 60,  60, 190, T8},                    // max score 250, just under 254
    {100, 100, 200, T16},                   // max score 300: over the byte, 16-bit range
    { 50,  60,  10, T16},                   // len1 < len2: 8-bit gate rejects
    { 20,  20,  -1, T16},                   // negative h0 with minval 19 >= 0
    {2000, 1500, 30, T16},                  // longer than MAX_SEQ_LEN8
    { 20,  20, -30, T1},                    // minval -10 < 0
    {MAX_SEQ_LEN16, 100, 0, T1},            // len1 at the 16-bit limit
    {100, MAX_SEQ_LEN16, 0, T1},            // len2 at the 16-bit limit
    {MAX_SEQ_LEN16 - 1, 200, MAX_SEQ_LEN16 - 200, T1},  // minval == MAX_SEQ_LEN16
    // The accepting side of each 16-bit/scalar bound (none qualifies for 8-bit: h0 < 0, or a
    // length >= MAX_SEQ_LEN8).
    { 20,  20, -20, T16},                   // minval == 0
    {MAX_SEQ_LEN16 - 1, 100, 0, T16},       // len1 just under the 16-bit limit
    {100, MAX_SEQ_LEN16 - 1, 0, T16},       // len2 just under the 16-bit limit
    {MAX_SEQ_LEN16 - 1, 200, MAX_SEQ_LEN16 - 201, T16},  // minval == MAX_SEQ_LEN16 - 1
};

// A pair built from template t, with every field distinct per index so a record copied to the
// wrong slot, or half-copied, is caught by the field-wise compare.
SeqPair make_pair(const Template &t, int idx)
{
    SeqPair p = {};
    p.idr = 1000 * idx + 1; p.idq = 1000 * idx + 2; p.id = idx;
    p.len1 = t.len1; p.len2 = t.len2; p.h0 = t.h0;
    p.seqid = idx; p.regid = 3 * idx + 1;
    p.score = idx + 11; p.tle = idx + 12; p.gtle = idx + 13; p.qle = idx + 14;
    p.gscore = idx + 15; p.max_off = idx + 16;
    p.tight_band = idx % 7; p.chain_band = idx % 5; p.ugp_walk_score = idx + 17;
    p.ugp_r_attempted = (uint8_t)(idx & 1); p.meth_hyp = (int8_t)(idx % 3 - 1);
    return p;
}

bool same_record(const SeqPair &a, const SeqPair &b)
{
    return a.idr == b.idr && a.idq == b.idq && a.id == b.id && a.len1 == b.len1 &&
           a.len2 == b.len2 && a.h0 == b.h0 && a.seqid == b.seqid && a.regid == b.regid &&
           a.score == b.score && a.tle == b.tle && a.gtle == b.gtle && a.qle == b.qle &&
           a.gscore == b.gscore && a.max_off == b.max_off && a.tight_band == b.tight_band &&
           a.chain_band == b.chain_band && a.ugp_walk_score == b.ugp_walk_score &&
           a.ugp_r_attempted == b.ugp_r_attempted && a.meth_hyp == b.meth_hyp;
}

// Partition `in` (templates by index into kTemplates) and check it against the oracle.
void check_partition(const std::vector<int> &tmpl_idx)
{
    const int n = (int)tmpl_idx.size();
    std::vector<SeqPair> pairs(n);
    std::vector<Tier> tier(n);
    for (int i = 0; i < n; i++) {
        pairs[i] = make_pair(kTemplates[tmpl_idx[i]], i);
        tier[i] = kTemplates[tmpl_idx[i]].tier;
    }

    // Oracle: stable sort of the input by its built-for tier.
    std::vector<int> order(n);
    for (int i = 0; i < n; i++) order[i] = i;
    std::stable_sort(order.begin(), order.end(), [&](int x, int y) { return tier[x] < tier[y]; });
    const int want8  = (int)std::count(tier.begin(), tier.end(), T8);
    const int want16 = (int)std::count(tier.begin(), tier.end(), T16);
    const int want1  = n - want8 - want16;

    // tempArray of exactly n slots, followed by a guard slot that must stay untouched.
    const SeqPair guard = make_pair(kTemplates[0], -7);
    std::vector<SeqPair> tmp(n + 1, guard);
    std::vector<SeqPair> got = pairs;
    int n8 = -1, n16 = -1, n1 = -1;
    sortPairsLenExt(got.data(), n, tmp.data(), n8, n16, n1, kA, kW, kZdrop);

    CHECK(n8 == want8);
    CHECK(n16 == want16);
    CHECK(n1 == want1);
    CHECK(same_record(tmp[n], guard));
    int bad = 0;
    for (int k = 0; k < n; k++)
        if (!same_record(got[k], pairs[order[k]])) {
            if (++bad <= 5) FAIL_CHECK("slot " << k << " holds id " << got[k].id
                                       << ", expected id " << pairs[order[k]].id);
        }
    CHECK(bad == 0);
}

}  // namespace

TEST_CASE("ext tier partition: each template takes the tier it was built for"
          * doctest::test_suite("unit/ext-tier-partition")) {
    if (bsw8_disabled()) { MESSAGE("BWAMEM3_DISABLE_BSW8 is set: skipped"); return; }
    for (const Template &t : kTemplates) {
        CAPTURE(t.len1); CAPTURE(t.len2); CAPTURE(t.h0);
        REQUIRE(bsw8_envelope_ok(t.len1, t.len2, kW, kA, kZdrop, t.h0) == (t.tier == T8));
        SeqPair p = make_pair(t, 0), tmp[1];
        int n8 = -1, n16 = -1, n1 = -1;
        sortPairsLenExt(&p, 1, tmp, n8, n16, n1, kA, kW, kZdrop);
        CHECK(n8  == (t.tier == T8  ? 1 : 0));
        CHECK(n16 == (t.tier == T16 ? 1 : 0));
        CHECK(n1  == (t.tier == T1  ? 1 : 0));
    }
}

TEST_CASE("ext tier partition: [8-bit | 16-bit | scalar], each tier in input order"
          * doctest::test_suite("unit/ext-tier-partition")) {
    if (bsw8_disabled()) { MESSAGE("BWAMEM3_DISABLE_BSW8 is set: skipped"); return; }
    const int nt = (int)(sizeof(kTemplates) / sizeof(kTemplates[0]));
    SUBCASE("tiers interleaved, scalar first") {
        std::vector<int> v;
        for (int r = 0; r < 4; r++) for (int i = nt - 1; i >= 0; i--) v.push_back(i);
        check_partition(v);
    }
    SUBCASE("only 16-bit and scalar pairs (no in-place run)") {
        check_partition({3, 7, 4, 8, 5, 9, 6, 10, 11, 12, 13, 14, 3, 7});
    }
    SUBCASE("only scalar pairs") { check_partition({7, 8, 9, 10, 7, 8}); }
    SUBCASE("only 16-bit pairs") { check_partition({3, 4, 5, 6, 11, 12, 13, 14, 3, 4}); }
    SUBCASE("8-bit run then the rest") { check_partition({0, 1, 2, 0, 1, 3, 7, 4, 8}); }
    SUBCASE("random batches") {
        std::mt19937 rng(20261030);
        for (int rep = 0; rep < 200; rep++) {
            const int n = 1 + (int)(rng() % 300);
            std::vector<int> v(n);
            for (int &x : v) x = (int)(rng() % nt);
            check_partition(v);
        }
    }
}

TEST_CASE("ext tier partition: an all-8-bit batch moves nothing and leaves tempArray alone"
          * doctest::test_suite("unit/ext-tier-partition")) {
    if (bsw8_disabled()) { MESSAGE("BWAMEM3_DISABLE_BSW8 is set: skipped"); return; }
    const int n = 257;
    std::vector<SeqPair> pairs(n);
    for (int i = 0; i < n; i++) pairs[i] = make_pair(kTemplates[i % 3], i);   // templates 0-2: 8-bit
    const SeqPair poison = make_pair(kTemplates[7], -1);
    std::vector<SeqPair> tmp(n, poison);
    std::vector<SeqPair> got = pairs;
    int n8 = -1, n16 = -1, n1 = -1;
    sortPairsLenExt(got.data(), n, tmp.data(), n8, n16, n1, kA, kW, kZdrop);
    CHECK(n8 == n);
    CHECK(n16 == 0);
    CHECK(n1 == 0);
    int moved = 0, touched = 0;
    for (int i = 0; i < n; i++) {
        moved   += !same_record(got[i], pairs[i]);
        touched += !same_record(tmp[i], poison);
    }
    CHECK(moved == 0);
    CHECK(touched == 0);
}

TEST_CASE("ext tier partition: an empty side returns zero counts and touches nothing"
          * doctest::test_suite("unit/ext-tier-partition")) {
    const SeqPair poison = make_pair(kTemplates[7], -1);
    SeqPair pa = poison, tmp = poison;
    int n8 = -1, n16 = -1, n1 = -1;
    sortPairsLenExt(&pa, 0, &tmp, n8, n16, n1, kA, kW, kZdrop);
    CHECK(n8 == 0);
    CHECK(n16 == 0);
    CHECK(n1 == 0);
    CHECK(same_record(pa, poison));
    CHECK(same_record(tmp, poison));
}
