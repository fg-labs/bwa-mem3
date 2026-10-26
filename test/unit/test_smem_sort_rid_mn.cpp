// Unit tests for smem_sort_rid_mn: the fused SMEM sort used by mem_collect_smem.
//
// Oracle: the two-step order it replaces -- a stable sort by rid (std::
// stable_sort, independent of the counting-sort helpers that smem_sort_rid_mn
// shares with FMI_search::sortSMEMs), then the shipped ks_introsort_mem_intv1
// (bwamem.cpp's KSORT_INIT(mem_intv1, intv_lt1)) on each rid block. Every
// oracle call also checks that FMI_search::sortSMEMs produces the same rid
// order as the std::stable_sort step, so sortSMEMs stays pinned too. ks_introsort is unstable, so on (rid, m, n) ties the output order is
// whatever permutation the introsort happens to apply. smem_sort_rid_mn must
// reproduce that exact permutation, so every test gives tied entries distinct
// payloads (info, k, l, s) and compares the full arrays field by field.

#include "doctest/doctest.h"
#include "FMI_search.h"
#include "ksort.h"

#include <algorithm>
#include <cstring>
#include <random>
#include <vector>

// Defined by KSORT_INIT(mem_intv1, SMEM, intv_lt1) in bwamem.cpp.
void ks_introsort_mem_intv1(size_t n, SMEM a[]);

namespace {

SMEM make_smem(uint32_t rid, uint32_t m, uint32_t n, int64_t tag)
{
    SMEM s;
    memset(&s, 0, sizeof(s));
    s.info = (uint64_t)tag * 7u + 1u;
    s.rid  = rid;
    s.m    = m;
    s.n    = n;
    s.k    = tag;
    s.l    = -tag;
    s.s    = tag * 3 + 1;
    return s;
}

bool same_smem(const SMEM &x, const SMEM &y);

// FMI_search::sortSMEMs on one block. It reads no index state, so an unloaded
// FMI_search is enough.
std::vector<SMEM> sort_smems_rid(std::vector<SMEM> a)
{
    static FMI_search fmi("unused-no-index");
    SmemSortScratch os = {NULL, 0, NULL, 0, NULL, 0};
    int64_t total[1] = {(int64_t)a.size()};
    fmi.sortSMEMs(a.data(), total, 1, 0, 1, os);
    _mm_free(os.cnt);
    _mm_free(os.tmp);
    return a;
}

// The pre-fusion order: stable by rid, then per-rid ks_introsort(mem_intv1).
std::vector<SMEM> oracle(std::vector<SMEM> a)
{
    std::vector<SMEM> by_rid = sort_smems_rid(a);
    std::stable_sort(a.begin(), a.end(),
                     [](const SMEM &x, const SMEM &y) { return x.rid < y.rid; });
    bool rid_step_same = by_rid.size() == a.size();
    for (size_t k = 0; rid_step_same && k < a.size(); k++)
        rid_step_same = same_smem(by_rid[k], a[k]);
    CHECK_MESSAGE(rid_step_same, "FMI_search::sortSMEMs differs from a stable sort by rid");
    size_t i = 0;
    while (i < a.size()) {
        size_t j = i + 1;
        while (j < a.size() && a[j].rid == a[i].rid) j++;
        ks_introsort_mem_intv1(j - i, a.data() + i);
        i = j;
    }
    return a;
}

std::vector<SMEM> fused(std::vector<SMEM> a, SmemSortScratch &scratch)
{
    smem_sort_rid_mn(a.data(), (int64_t)a.size(), scratch);
    return a;
}

bool same_smem(const SMEM &x, const SMEM &y)
{
    return x.info == y.info && x.rid == y.rid && x.m == y.m && x.n == y.n &&
           x.k == y.k && x.l == y.l && x.s == y.s;
}

// One assertion per array (not per field) keeps the large sweeps fast; on a
// mismatch the first differing index is reported.
// Elements compared so far; sweeps assert on it so an empty sweep cannot pass.
size_t g_checked = 0;

void check_same(const std::vector<SMEM> &got, const std::vector<SMEM> &want)
{
    REQUIRE(got.size() == want.size());
    g_checked += got.size();
    long first_diff = -1;
    for (size_t i = 0; i < got.size(); i++)
        if (!same_smem(got[i], want[i])) { first_diff = (long)i; break; }
    CHECK(first_diff == -1);
}

void free_scratch(SmemSortScratch &s)
{
    _mm_free(s.cnt);
    _mm_free(s.tmp);
    _mm_free(s.pairs);
    s = SmemSortScratch{NULL, 0, NULL, 0, NULL, 0};
}

// Random array: nreads reads starting at rid0, three seeding passes appended in
// read order (as mem_collect_smem lays them out), keys drawn from a small space
// so (m, n) ties are frequent.
std::vector<SMEM> random_batch(std::mt19937_64 &rng, uint32_t rid0, uint32_t nreads,
                               int max_per_read, uint32_t key_space)
{
    std::vector<SMEM> a;
    int64_t tag = 0;
    for (int pass = 0; pass < 3; pass++) {
        for (uint32_t r = 0; r < nreads; r++) {
            int cnt = (int)(rng() % (uint64_t)(max_per_read + 1));
            for (int c = 0; c < cnt; c++) {
                uint32_t m = (uint32_t)(rng() % key_space);
                uint32_t n = m + (uint32_t)(rng() % key_space);
                a.push_back(make_smem(rid0 + r, m, n, tag++));
            }
        }
    }
    return a;
}

// ---- Adversarial input for ks_introsort's combsort fallback ----
// ks_introsort switches a partition to ks_combsort once its depth budget
// (2 * ceil(log2 n)) runs out. A median-of-3 killer drives it there. The input
// is built with McIlroy's "antiqsort" adversary run against the same ksort.h
// template: values start as "gas" and are frozen lazily so that every pivot
// choice is as bad as possible.
struct AdvItem { uint32_t id; };
std::vector<uint32_t> g_adv_val;
uint32_t g_adv_gas, g_adv_nsolid, g_adv_candidate;
bool adv_gas_lt(uint32_t x, uint32_t y)
{
    std::vector<uint32_t> &v = g_adv_val;
    if (v[x] == g_adv_gas && v[y] == g_adv_gas) {
        if (x == g_adv_candidate) v[x] = g_adv_nsolid++;
        else                      v[y] = g_adv_nsolid++;
    }
    if (v[x] == g_adv_gas)      g_adv_candidate = x;
    else if (v[y] == g_adv_gas) g_adv_candidate = y;
    return v[x] < v[y];
}

// Probe instantiation on (key, id): counts comparisons between two in-array
// elements a fixed distance d > 1 apart. Only combsort makes many of those
// (a full pass at one gap); partitioning compares against a pivot copy, the
// median-of-3 makes three, and insertion sort compares neighbours.
struct AdvProbe { uint64_t key; uint32_t id; };
const AdvProbe *g_probe_base = nullptr;
size_t g_probe_n = 0;
std::vector<size_t> g_probe_gap_hits;
bool adv_probe_lt(const AdvProbe *x, const AdvProbe *y)
{
    if (g_probe_base && x >= g_probe_base && x < g_probe_base + g_probe_n &&
        y >= g_probe_base && y < g_probe_base + g_probe_n) {
        size_t d = (size_t)(x > y ? x - y : y - x);
        if (d > 1) g_probe_gap_hits[d]++;
    }
    return x->key < y->key;
}

}  // namespace

#define adv_gas_lt_m(a, b) adv_gas_lt((a).id, (b).id)
KSORT_INIT(ssf_adv_gas, AdvItem, adv_gas_lt_m)
#define adv_probe_lt_m(a, b) adv_probe_lt(&(a), &(b))
KSORT_INIT(ssf_adv_probe, AdvProbe, adv_probe_lt_m)

namespace {

// Keys (one per position) that push ks_introsort into its combsort fallback.
std::vector<uint64_t> combsort_killer_keys(size_t n)
{
    g_adv_val.assign(n, (uint32_t)n);
    g_adv_gas = (uint32_t)n; g_adv_nsolid = 0; g_adv_candidate = 0;
    std::vector<AdvItem> items(n);
    for (size_t i = 0; i < n; i++) items[i].id = (uint32_t)i;
    ks_introsort_ssf_adv_gas(n, items.data());
    std::vector<uint64_t> keys(n);
    for (size_t i = 0; i < n; i++) {
        if (g_adv_val[i] == g_adv_gas) g_adv_val[i] = g_adv_nsolid++;
        keys[i] = g_adv_val[i];
    }
    return keys;
}

// Largest number of same-gap (d > 1) in-array comparisons made while sorting keys.
size_t max_same_gap_compares(const std::vector<uint64_t> &keys)
{
    std::vector<AdvProbe> a(keys.size());
    for (size_t i = 0; i < keys.size(); i++) { a[i].key = keys[i]; a[i].id = (uint32_t)i; }
    g_probe_base = a.data(); g_probe_n = a.size();
    g_probe_gap_hits.assign(a.size() + 1, 0);
    ks_introsort_ssf_adv_probe(a.size(), a.data());
    g_probe_base = nullptr;
    return *std::max_element(g_probe_gap_hits.begin(), g_probe_gap_hits.end());
}

}  // namespace

TEST_CASE("smem_sort_rid_mn: empty batch, a single SMEM, and two elements"
          * doctest::test_suite("unit/smem"))
{
    SmemSortScratch scratch = {NULL, 0, NULL, 0, NULL, 0};
    std::vector<SMEM> none;
    check_same(fused(none, scratch), oracle(none));

    std::vector<SMEM> one = {make_smem(5, 3, 9, 1)};
    check_same(fused(one, scratch), oracle(one));

    std::vector<SMEM> two_swap = {make_smem(5, 4, 9, 1), make_smem(5, 3, 9, 2)};
    check_same(fused(two_swap, scratch), oracle(two_swap));

    std::vector<SMEM> two_tie = {make_smem(5, 3, 9, 1), make_smem(5, 3, 9, 2)};
    check_same(fused(two_tie, scratch), oracle(two_tie));

    std::vector<SMEM> two_rids = {make_smem(6, 1, 2, 1), make_smem(5, 3, 9, 2)};
    check_same(fused(two_rids, scratch), oracle(two_rids));
    free_scratch(scratch);
}

TEST_CASE("smem_sort_rid_mn: key is m then n, compared as one 64-bit value"
          * doctest::test_suite("unit/smem"))
{
    SmemSortScratch scratch = {NULL, 0, NULL, 0, NULL, 0};
    std::vector<SMEM> a = {
        make_smem(0, 2, 0xFFFFFFFFu, 1), make_smem(0, 3, 0, 2),
        make_smem(0, 2, 1, 3),           make_smem(0, 0xFFFFFFFFu, 0xFFFFFFFFu, 4),
        make_smem(0, 0, 0, 5),           make_smem(0, 2, 1, 6),
    };
    check_same(fused(a, scratch), oracle(a));
    free_scratch(scratch);
}

TEST_CASE("smem_sort_rid_mn: ties keep the introsort's permutation (small and large blocks)"
          * doctest::test_suite("unit/smem"))
{
    SmemSortScratch scratch = {NULL, 0, NULL, 0, NULL, 0};
    // All-equal keys in one read: every pair is a tie, the permutation is all
    // that distinguishes outputs. Sizes straddle the 16-element insertion-sort
    // cut-off and reach the partitioning path.
    size_t checked0 = g_checked;
    for (int sz : {3, 16, 17, 33, 200, 1000}) {
        std::vector<SMEM> a;
        for (int i = 0; i < sz; i++) a.push_back(make_smem(7, 10, 20, i));
        CAPTURE(sz);
        check_same(fused(a, scratch), oracle(a));
    }
    REQUIRE(g_checked - checked0 == (size_t)(3 + 16 + 17 + 33 + 200 + 1000));
    // Two keys interleaved, many ties each.
    std::vector<SMEM> b;
    for (int i = 0; i < 500; i++) b.push_back(make_smem(1, (i % 2) ? 4 : 9, 30, i));
    check_same(fused(b, scratch), oracle(b));
    free_scratch(scratch);
}

TEST_CASE("smem_sort_rid_mn: sorted, reversed and sawtooth large blocks"
          * doctest::test_suite("unit/smem"))
{
    SmemSortScratch scratch = {NULL, 0, NULL, 0, NULL, 0};
    std::vector<SMEM> asc, desc, saw;
    for (int i = 0; i < 3000; i++) {
        asc.push_back(make_smem(2, (uint32_t)i, (uint32_t)i + 5, i));
        desc.push_back(make_smem(2, (uint32_t)(3000 - i), 1, i));
        saw.push_back(make_smem(2, (uint32_t)(i % 37), (uint32_t)(i % 5), i));
    }
    check_same(fused(asc, scratch), oracle(asc));
    check_same(fused(desc, scratch), oracle(desc));
    check_same(fused(saw, scratch), oracle(saw));
    free_scratch(scratch);
}

TEST_CASE("smem_sort_rid_mn: random batches match the two-step sort, scratch reused and grown"
          * doctest::test_suite("unit/smem"))
{
    std::mt19937_64 rng(0x5eed5eedULL);
    SmemSortScratch scratch = {NULL, 0, NULL, 0, NULL, 0};
    // Mixed sizes so the reused scratch is grown, then reused smaller, then
    // grown again; nonzero rid0 exercises the minRid offset.
    struct { uint32_t rid0, nreads; int max_per_read; uint32_t key_space; } cfg[] = {
        {0, 10, 4, 3},     {1000, 200, 12, 8},  {0, 3, 60, 4},
        {77, 50, 2, 2},    {0, 1, 400, 6},      {5, 2000, 20, 50},
        {123456, 64, 30, 1}, {0, 500, 8, 1000},
    };
    int batches = 0, nonempty = 0;
    size_t elems = 0, checked0 = g_checked;
    for (int rep = 0; rep < 4; rep++) {
        for (auto &c : cfg) {
            std::vector<SMEM> a = random_batch(rng, c.rid0, c.nreads, c.max_per_read, c.key_space);
            batches++; elems += a.size(); if (a.size() > 1) nonempty++;
            CAPTURE(rep); CAPTURE(c.rid0); CAPTURE(c.nreads); CAPTURE(a.size());
            check_same(fused(a, scratch), oracle(a));
        }
    }
    // An empty sweep must not pass vacuously.
    REQUIRE(batches == 4 * (int)(sizeof(cfg) / sizeof(cfg[0])));
    REQUIRE(nonempty > batches / 2);
    REQUIRE(elems > 10000);
    REQUIRE(g_checked - checked0 == elems);
    free_scratch(scratch);
}

TEST_CASE("smem_sort_rid_mn: reads with no SMEMs between, and a single rid"
          * doctest::test_suite("unit/smem"))
{
    std::mt19937_64 rng(42);
    SmemSortScratch scratch = {NULL, 0, NULL, 0, NULL, 0};
    std::vector<SMEM> a;
    int64_t tag = 0;
    for (uint32_t r : {3u, 9u, 9u, 40u, 3u, 40u, 9u}) {   // sparse rids, empty buckets between
        for (int c = 0; c < 20; c++)
            a.push_back(make_smem(r, (uint32_t)(rng() % 5), (uint32_t)(rng() % 3), tag++));
    }
    check_same(fused(a, scratch), oracle(a));

    std::vector<SMEM> single;
    for (int c = 0; c < 64; c++)
        single.push_back(make_smem(11, (uint32_t)(rng() % 6), (uint32_t)(rng() % 6), c));
    check_same(fused(single, scratch), oracle(single));
    free_scratch(scratch);
}

TEST_CASE("smem_sort_rid_mn: tied keys through the introsort's combsort fallback"
          * doctest::test_suite("unit/smem"))
{
    const size_t n = 4000;
    // Control: a shuffled input never runs out of depth, so the probe stays low.
    std::vector<uint64_t> shuffled(n);
    for (size_t i = 0; i < n; i++) shuffled[i] = i;
    std::shuffle(shuffled.begin(), shuffled.end(), std::mt19937_64(7));
    REQUIRE(max_same_gap_compares(shuffled) < n / 32);

    std::vector<uint64_t> keys = combsort_killer_keys(n);
    // The killer must actually reach combsort: a full pass at one gap compares
    // far more same-distance pairs than partitioning ever does.
    REQUIRE(max_same_gap_compares(keys) > n / 8);

    // Halve the keys so every value appears twice (tied (m, n) with distinct
    // payloads), and confirm the tied input still reaches combsort.
    std::vector<uint64_t> tied(n);
    for (size_t i = 0; i < n; i++) tied[i] = keys[i] / 2;
    REQUIRE(max_same_gap_compares(tied) > n / 8);

    SmemSortScratch scratch = {NULL, 0, NULL, 0, NULL, 0};
    size_t checked0 = g_checked;
    for (const std::vector<uint64_t> *ks : {&keys, &tied}) {
        std::vector<SMEM> a;
        for (size_t i = 0; i < n; i++) {
            uint64_t k = (*ks)[i];
            // key = m << 32 | n: put the high half of the value in m.
            a.push_back(make_smem(9, (uint32_t)(k >> 6), (uint32_t)(k & 63), (int64_t)i));
        }
        // Precede the bucket with another read so the killer bucket starts at
        // a nonzero offset in the pairs array.
        std::vector<SMEM> b = {make_smem(8, 1, 1, -1), make_smem(8, 0, 1, -2)};
        b.insert(b.end(), a.begin(), a.end());
        check_same(fused(b, scratch), oracle(b));
    }
    REQUIRE(g_checked - checked0 == 2 * (n + 2));
    free_scratch(scratch);
}
