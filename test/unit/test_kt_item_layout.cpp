// test/unit/test_kt_item_layout.cpp — how a kt_for() pass is cut into work items.
//
// kt_for() hands each worker items of BATCH_SIZE reads, except about the last
// n_threads * BATCH_SIZE reads of the pass (rounded up to whole batches), which
// are cut into quarter-size items so the workers reach the end-of-pass barrier
// closer together. That is output-neutral only if the items still tile the
// pass exactly: every read in exactly one item, in order, no item larger than
// BATCH_SIZE (per-thread scratch is sized for one), and no item starting on an
// odd read (a read pair must never straddle two items -- the paired-end
// callbacks index pairs as seq_id >> 1).
//
// The first group checks kt_item_layout / kt_item_range directly over a grid of
// pass sizes, thread counts and batch sizes; the second drives real kt_for()
// passes and records what the callbacks were actually handed, which also covers
// the strided dispatch and the work-stealing loop. All inputs are synthetic.

#include "doctest/doctest.h"

#include "kthread.h"
#include "bwamem.h"   // worker_t
#include "macro.h"    // BATCH_SIZE

#include <cstdlib>
#include <vector>

namespace {

struct Item { long st, ed; };

// Every item of a layout, in index order, via kt_item_range.
std::vector<Item> items_of(const kt_item_layout_t &L)
{
    std::vector<Item> out;
    for (long k = 0;; ++k) {
        int st = -1, ed = -1;
        if (!kt_item_range(&L, k, &st, &ed)) break;
        out.push_back({st, ed});
        if (out.size() > (size_t)L.n + 1) break;   // never loops forever (n_items CHECK fails)
    }
    return out;
}

// The invariants every layout must satisfy, independent of how it is computed.
void check_layout(long n, int n_threads, long batch)
{
    CAPTURE(n); CAPTURE(n_threads); CAPTURE(batch);
    kt_item_layout_t L;
    kt_item_layout(n, n_threads, batch, &L);
    const std::vector<Item> items = items_of(L);

    CHECK((long)items.size() == L.n_items);
    // Tiles [0, n): contiguous, in order, no gaps, no overlaps, nothing past n.
    // Violations are counted rather than CHECKed per item: a pass can have
    // thousands of items, and a per-item CHECK would dominate the runtime.
    long next = 0, gaps = 0, empty = 0, oversize = 0, odd_start = 0, odd_len = 0;
    for (size_t k = 0; k < items.size(); ++k) {
        const long len = items[k].ed - items[k].st;
        gaps      += items[k].st != next;
        empty     += len <= 0;
        oversize  += len > batch;                       // fits a full-batch scratch
        odd_start += items[k].st % 2 != 0;              // a pair never straddles items
        odd_len   += k + 1 < items.size() && len % 2 != 0;  // only the last may be odd
        next = items[k].ed;
    }
    CHECK(gaps == 0);
    CHECK(empty == 0);
    CHECK(oversize == 0);
    CHECK(odd_start == 0);
    CHECK(odd_len == 0);
    CHECK(next == n);

    // The head is the legacy layout exactly: item k is [k*batch, (k+1)*batch).
    long not_legacy = 0;
    for (long k = 0; k < L.n_full; ++k)
        not_legacy += items[k].st != k * batch || items[k].ed != (k + 1) * batch;
    CHECK(not_legacy == 0);
    // The finer items cover the last n_threads * batch reads (rounded up to a
    // whole batch so the head stays batch-aligned), or the whole pass if smaller.
    const long tail = (long)n_threads * batch;
    const long tail_reads = n - L.full_end;
    CHECK(L.full_end % batch == 0);
    if (n > tail) {
        CHECK(tail_reads >= tail);
        CHECK(tail_reads < tail + batch);
    } else {
        CHECK(L.full_end == 0);
    }
    long tail_oversize = 0;
    for (long k = L.n_full; k < L.n_items; ++k)
        tail_oversize += items[k].ed - items[k].st > L.tail_batch;
    CHECK(tail_oversize == 0);
}

} // namespace

TEST_CASE("kt_item_layout: items tile every pass exactly, in order, on even boundaries"
          * doctest::test_suite("unit/kt_for")) {
    const long batches[] = { 512, 1024, 16 };
    const int threads[] = { 1, 2, 3, 7, 16, 64 };
    for (long batch : batches) {
        for (int nt : threads) {
            const long t = (long)nt * batch;
            const long ns[] = { 0, 1, 2, 3, batch - 1, batch, batch + 1, batch + 2,
                                t - 2, t - 1, t, t + 1, t + 2, t + batch - 1, t + batch,
                                t + batch + 1, 3 * t + 5, 10 * t + batch / 2 + 2,
                                1000000, 1000001 };
            for (long n : ns) check_layout(n, nt, batch);
        }
    }
}

TEST_CASE("kt_item_layout: the tail is cut into quarter-batch items"
          * doctest::test_suite("unit/kt_for")) {
    // 16 threads, 1024-read batches, a 40-batch pass: the first 24 batches stay
    // whole and the last 16 become 64 items of 256 reads.
    kt_item_layout_t L;
    kt_item_layout(40 * 1024, 16, 1024, &L);
    CHECK(L.tail_batch == 256);
    CHECK(L.n_full == 24);
    CHECK(L.full_end == 24 * 1024);
    CHECK(L.n_items == 24 + 64);

    // A pass no bigger than one round is all quarter items.
    kt_item_layout(16 * 512, 16, 512, &L);
    CHECK(L.n_full == 0);
    CHECK(L.n_items == 64);

    // A ragged end: the last item is short, the rest are full quarter items.
    kt_item_layout(5000, 4, 512, &L);
    const std::vector<Item> items = items_of(L);
    REQUIRE(!items.empty());
    CHECK(items.back().ed == 5000);
    CHECK(items.back().ed - items.back().st == (5000 - L.full_end) % 128);
}

TEST_CASE("kt_item_range: indices outside the layout are rejected"
          * doctest::test_suite("unit/kt_for")) {
    kt_item_layout_t L;
    kt_item_layout(3000, 2, 512, &L);
    int st = -7, ed = -7;
    CHECK(kt_item_range(&L, L.n_items, &st, &ed) == 0);
    CHECK(kt_item_range(&L, L.n_items + 100, &st, &ed) == 0);
    CHECK(kt_item_range(&L, -1, &st, &ed) == 0);
    CHECK(st == -7);
    CHECK(ed == -7);

    kt_item_layout(0, 16, 512, &L);
    CHECK(L.n_items == 0);
    CHECK(kt_item_range(&L, 0, &st, &ed) == 0);
}

// ---- real kt_for() passes ----------------------------------------------------

namespace {

// Per read: how many items covered it. Per call: whether the item was legal.
std::vector<int> g_hits;
int g_bad_items;

void mark_item(void *data, int st, int len, int tid)
{
    worker_t *w = (worker_t *) data;
    if (st < 0 || st % 2 != 0 || len <= 0 || len > BATCH_SIZE || tid < 0 || tid >= w->nthreads
        || (size_t)st + (size_t)len > g_hits.size()) {
        __sync_fetch_and_add(&g_bad_items, 1);
        return;   // never index past g_hits
    }
    for (int i = st; i < st + len; ++i) __sync_fetch_and_add(&g_hits[i], 1);
}

} // namespace

TEST_CASE("kt_for: every read is handed to exactly one callback, in legal items"
          * doctest::test_suite("unit/kt_for")) {
    worker_t *w = (worker_t *) calloc(1, sizeof(worker_t));
    REQUIRE(w != NULL);
    const int threads[] = { 1, 3, 16 };
    for (int nt : threads) {
        w->nthreads = nt;
        const int ns[] = { 0, 1, 2, BATCH_SIZE, nt * BATCH_SIZE - 2, nt * BATCH_SIZE + 1,
                           7 * nt * BATCH_SIZE + 333 };
        for (int n : ns) {
            CAPTURE(nt); CAPTURE(n);
            g_hits.assign((size_t)n, 0);
            g_bad_items = 0;
            kt_for(mark_item, w, n);
            CHECK(g_bad_items == 0);
            int wrong = 0;
            for (int i = 0; i < n; ++i) wrong += g_hits[i] != 1;
            CHECK(wrong == 0);
        }
    }
    kt_pool_destroy();
    free(w);
}
