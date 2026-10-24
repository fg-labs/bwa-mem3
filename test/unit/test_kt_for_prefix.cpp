// test/unit/test_kt_for_prefix.cpp — kt_for_prefix hands its caller finished
// work items at the front of a pass, in read order, while the pass still runs.
//
// The pipeline's output streaming (stream_emit_prefix in fastmap.cpp) writes
// the reads it is offered straight to the output, so the contract is what keeps
// the output byte-identical: every offer is exactly the next work item, in item
// order, and only once func has returned for it (whatever order the workers
// finish in); a decline (0) is offered again after the next item finishes,
// including one that finishes while on_prefix is still running; and the pass
// itself still runs every item, as kt_for does.
//
// Items are made to finish out of order by sleeping in the first ones. Slow
// last items keep the pass open so the offers can be observed while they run:
// they are held until the caller has consumed a set number of items, not for a
// fixed time, so a slow or loaded runner cannot end the pass before the offers
// the checks depend on. The decline cases also hold the later fast items until
// the declined call has started (or has returned), so their completions land
// during (or after) it rather than all before. Each hold gives up after 10 s,
// and any hold that gave up fails the case, so a broken kt_for_prefix fails
// the checks instead of hanging. All inputs are synthetic.

#include "doctest/doctest.h"

#include "kthread.h"
#include "bwamem.h"   // worker_t
#include "macro.h"    // BATCH_SIZE

#include <atomic>
#include <chrono>
#include <cstdlib>
#include <vector>
#include <unistd.h>

namespace {

struct Item { int st, ed; };

std::vector<Item> layout_items(int n, int nt)
{
    kt_item_layout_t L;
    kt_item_layout(n, nt, BATCH_SIZE, &L);
    std::vector<Item> out;
    int st, ed;
    for (long k = 0; kt_item_range(&L, k, &st, &ed); ++k) out.push_back({st, ed});
    return out;
}

// Per read: set by func once its item is done (after any sleep).
std::vector<std::atomic<int>> *g_done;
std::vector<Item> g_items;
long g_slow_from = -1;        // items from here on are held, keeping the pass open
int g_slow_release_at = 1;    // ... until this many items have been consumed
std::atomic<int> g_consumed{0};
std::atomic<int> g_slow_started{0};    // slow items whose func is running
std::atomic<int> g_slow_finished{0};
std::atomic<int> g_slow_timed_out{0};  // a slow item gave up waiting
long g_gate_from = -1;        // fast items from here on wait for g_gate_open
std::atomic<int> g_gate_open{0};
std::atomic<int> g_gate_timed_out{0};  // a gated item gave up waiting
int g_early_sleep_items = 0;  // the first items sleep briefly, so later ones finish first
long g_trickle_from = -1;     // items from here on sleep a little each (completions keep coming)
useconds_t g_trickle_us = 1000;
std::atomic<int> g_runs{0};   // func calls

// Waits until `ready` holds, for at most 10 s. Returns whether it did.
template <class F> bool wait_until(F ready)
{
    const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(10);
    while (!ready()) {
        if (std::chrono::steady_clock::now() >= deadline) return false;
        usleep(200);
    }
    return true;
}

bool is_slow(long k) { return g_slow_from >= 0 && k >= g_slow_from; }

void run_item(void *data, int st, int len, int tid)
{
    (void) data; (void) tid;
    long k = -1;
    for (size_t i = 0; i < g_items.size(); ++i)
        if (g_items[i].st == st) { k = (long) i; break; }
    if (k >= 0 && k < g_early_sleep_items) usleep(500);
    if (g_gate_from >= 0 && k >= g_gate_from && !is_slow(k))
        if (!wait_until([] { return g_gate_open.load() != 0; }))
            g_gate_timed_out.store(1);
    if (g_trickle_from >= 0 && k >= g_trickle_from && !is_slow(k)) usleep(g_trickle_us);
    if (is_slow(k)) {
        g_slow_started.fetch_add(1);
        if (!wait_until([] { return g_consumed.load() >= g_slow_release_at; }))
            g_slow_timed_out.store(1);
    }
    for (int i = st; i < st + len; ++i) (*g_done)[i].store(1, std::memory_order_relaxed);
    if (is_slow(k)) g_slow_finished.store(1);
    g_runs.fetch_add(1);
}

struct Ctx {
    int next = 0;            // first read not yet consumed
    int calls = 0;
    int bad_order = 0;       // offer did not start at `next`
    int not_one_item = 0;    // offer was not exactly the next layout item
    int unfinished = 0;      // a read offered before its item had finished
    int consumed_while_slow = 0;  // items consumed while no slow item had finished
    int decline_first = 0;   // decline the first offer
    int wait_in_decline = 0; // ... after every fast item has finished meanwhile
    int open_gate_on_return = 0;  // ... and open g_gate as it returns, not as it starts
    int decline_wait_timed_out = 0;  // the wait_in_decline wait gave up
    size_t item = 0;
};

int on_prefix(void *ctx_, int st, int ed)
{
    Ctx *c = (Ctx *) ctx_;
    ++c->calls;
    if (c->decline_first && c->calls == 1) {
        if (!c->open_gate_on_return) g_gate_open.store(1);
        if (c->wait_in_decline) {
            // Let every fast item finish while this call runs, so no completion
            // arrives after it returns: only counting those that landed during
            // the call can prompt the retry. A worker records an item's
            // completion before it takes its next item, so once every fast item
            // has run and every worker is held in a slow item (one per worker),
            // every fast completion has been recorded.
            const int n_slow = (int) g_items.size() - (int) g_slow_from;
            if (!wait_until([n_slow] {
                    return g_runs.load() >= (int) g_items.size() - n_slow
                           && g_slow_started.load() >= n_slow;
                }))
                c->decline_wait_timed_out = 1;
        }
        g_gate_open.store(1);
        return 0;
    }
    c->bad_order += st != c->next;
    c->not_one_item += c->item >= g_items.size()
                       || g_items[c->item].st != st || g_items[c->item].ed != ed;
    for (int i = st; i < ed; ++i) c->unfinished += (*g_done)[i].load(std::memory_order_relaxed) != 1;
    c->consumed_while_slow += !g_slow_finished.load();
    c->next = ed;
    ++c->item;
    g_consumed.fetch_add(1);
    return 1;
}

void reset(int n, int nt, long slow_from, int early_sleep)
{
    delete g_done;
    g_done = new std::vector<std::atomic<int>>((size_t) n);
    g_items = layout_items(n, nt);
    g_slow_from = slow_from;
    g_slow_release_at = 1;
    g_consumed.store(0);
    g_slow_started.store(0);
    g_slow_finished.store(0);
    g_slow_timed_out.store(0);
    g_gate_from = -1;
    g_gate_open.store(0);
    g_gate_timed_out.store(0);
    g_early_sleep_items = early_sleep;
    g_runs.store(0);
    g_trickle_from = -1;
    g_trickle_us = 1000;
}

} // namespace

TEST_CASE("kt_for_prefix: offers are the finished items, one at a time, in read order"
          * doctest::test_suite("unit/kt_for")) {
    const int threads[] = { 1, 3, 8 };
    for (int nt : threads) {
        const int ns[] = { 2 * BATCH_SIZE + 6, 3 * nt * BATCH_SIZE + 2 };
        for (int n : ns) {
            CAPTURE(nt); CAPTURE(n);
            reset(n, nt, /*slow_from=*/-1, /*early_sleep=*/0);
            const long last = (long) g_items.size() - 1;
            g_slow_from = last;
            g_slow_release_at = (int) last;   // every earlier item consumed first
            g_early_sleep_items = (int) g_items.size() / 2;
            Ctx c;
            kt_for_prefix(run_item, &c, n, nt, on_prefix, &c);
            CHECK(g_runs.load() == (int) g_items.size());   // every item ran
            int not_done = 0;
            for (int i = 0; i < n; ++i) not_done += (*g_done)[i].load() != 1;
            CHECK(not_done == 0);
            CHECK(c.bad_order == 0);
            CHECK(c.not_one_item == 0);
            CHECK(c.unfinished == 0);
            // The slow last item holds the pass open, so every earlier item is
            // offered (and consumed) before the pass ends.
            CHECK(c.item >= g_items.size() - 1);
            CHECK(c.consumed_while_slow >= 1);
            CHECK(g_slow_timed_out.load() == 0);
        }
    }
    kt_pool_destroy();
}

TEST_CASE("kt_for_prefix: a declined offer is made again after the next item finishes"
          * doctest::test_suite("unit/kt_for")) {
    const int nt = 3, n = 4 * BATCH_SIZE;
    reset(n, nt, -1, 0);
    g_slow_from = (long) g_items.size() - 1;
    // Every item after the first waits until the declined call is returning
    // (the gate opens as it returns), then takes g_trickle_us, so the
    // completions that prompt the retry come after kt_for_prefix has recorded
    // the decline unless the caller stalls for longer than that between the
    // two. (One that lands during the call is the next case.)
    g_gate_from = 1;
    g_trickle_from = 1;
    Ctx c;
    c.decline_first = 1;
    c.open_gate_on_return = 1;
    kt_for_prefix(run_item, &c, n, nt, on_prefix, &c);
    CHECK(c.calls >= 2);
    CHECK(c.bad_order == 0);
    CHECK(c.not_one_item == 0);
    CHECK(c.unfinished == 0);
    CHECK(c.consumed_while_slow >= 1);
    CHECK(g_slow_timed_out.load() == 0);
    CHECK(g_gate_timed_out.load() == 0);
    CHECK(g_runs.load() == (int) g_items.size());
    kt_pool_destroy();
}

TEST_CASE("kt_for_prefix: an item that finishes during a declined offer counts"
          * doctest::test_suite("unit/kt_for")) {
    const int nt = 2, n = 3 * BATCH_SIZE;
    reset(n, nt, -1, 0);
    g_slow_from = (long) g_items.size() - nt;  // one per worker, held until item 0 is consumed
    g_gate_from = 1;     // the other items wait for the declined call to start,
    g_trickle_from = 1;  // so they finish during it, not before it
    Ctx c;
    c.decline_first = 1;
    c.wait_in_decline = 1;
    kt_for_prefix(run_item, &c, n, nt, on_prefix, &c);
    CHECK(c.calls >= 2);
    CHECK(c.bad_order == 0);
    CHECK(c.unfinished == 0);
    // Every fast item finished during the declined call; the retry came from
    // counting them, not from a later completion, so it ran while the slow
    // items still did.
    CHECK(c.consumed_while_slow >= 1);
    CHECK(c.decline_wait_timed_out == 0);
    CHECK(g_slow_timed_out.load() == 0);
    CHECK(g_gate_timed_out.load() == 0);
    kt_pool_destroy();
}

TEST_CASE("kt_for_prefix: an empty pass offers nothing, and kt_for still works after"
          * doctest::test_suite("unit/kt_for")) {
    reset(0, 4, -1, 0);
    Ctx c;
    kt_for_prefix(run_item, &c, 0, 4, on_prefix, &c);
    CHECK(c.calls == 0);
    CHECK(g_runs.load() == 0);
    // And the plain kt_for still works on the same pool afterwards.
    worker_t *w = (worker_t *) calloc(1, sizeof(worker_t));
    REQUIRE(w != NULL);
    w->nthreads = 4;
    reset(BATCH_SIZE, 4, -1, 0);
    kt_for(run_item, w, BATCH_SIZE);
    CHECK(g_runs.load() == (int) g_items.size());
    free(w);
    kt_pool_destroy();
    delete g_done; g_done = NULL;
}
