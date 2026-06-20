// Regression test for the kt_for() persistent worker pool.
//
// bwa-mem3 calls kt_for() three times per chunk (worker_bwt / worker_aln /
// worker_sam) across many chunks. Per-worker scratch buffers (e.g.
// mmc->enc_qdb[tid]) persist across chunks and are grown with realloc(). If
// kt_for() spawned a fresh set of pthreads on every call (the historical
// behaviour), the buffer for slot `tid` would be allocated by one chunk's
// worker thread and realloc'd by a *different* OS thread on the next chunk —
// after the first thread exited. mimalloc v3.x mishandles that cross-thread
// realloc of an abandoned-heap block and corrupts the heap (benign under glibc
// and mimalloc v2.x, so no sanitizer flags it), crashing several chunks later.
//
// The fix makes kt_for() drive a persistent worker pool, so worker slot `tid`
// is always the same OS thread: every realloc of that slot's buffers happens on
// the thread that allocated them. This test asserts exactly that invariant —
// the same OS thread services a given slot across successive kt_for() calls.
// It passes with the pool and fails with per-call threads.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <pthread.h>
#include <unistd.h>

#include "kthread.h"
#include "bwamem.h"   // worker_t
#include "macro.h"    // BATCH_SIZE

static const int NT = 4;
static pthread_t g_seen[2][NT];
static int       g_seen_ok[2][NT];  // 1 once slot t has been observed on call c
static int       g_call = 0;

// Sums the work items (len) handed to the callback, so the edge-case phase can
// check that a recreated single-thread pool runs all the work and that a
// zero-item call is a clean no-op rather than a hang.
static long g_processed = 0;
static void count_work(void *data, int st, int len, int tid) {
    (void)data; (void)st; (void)tid;
    __sync_fetch_and_add(&g_processed, (long)len);
}

// kt_for() invokes this as func(data, st, len, tid). `tid` is the worker slot
// index; record which OS thread is servicing that slot on this call. The brief
// sleep keeps all NT workers busy on their own strided batch simultaneously so
// none finishes early and steals another slot's work — that keeps the
// slot->thread mapping stable for this test (work-stealing is irrelevant to the
// invariant under test: thread persistence across kt_for() calls).
static void record_slot_thread(void *data, int st, int len, int tid) {
    (void)data; (void)st; (void)len;
    if (tid >= 0 && tid < NT) {
        g_seen[g_call][tid] = pthread_self();
        g_seen_ok[g_call][tid] = 1;
    }
    usleep(20000);  // 20 ms: long enough to overlap all workers, trivial total
}

int main(void) {
    // kt_for() reads only ->nthreads from its data argument. worker_t is large
    // (per-thread buffer arrays), so allocate it on the heap.
    worker_t *w = (worker_t *) calloc(1, sizeof(worker_t));
    if (w == NULL) { fprintf(stderr, "calloc failed\n"); return 2; }
    w->nthreads = NT;

    // n = NT * BATCH_SIZE gives exactly one strided batch per worker slot, so
    // every slot's thread runs record_slot_thread() at least once.
    const int n = NT * BATCH_SIZE;

    g_call = 0; kt_for(record_slot_thread, w, n);
    g_call = 1; kt_for(record_slot_thread, w, n);

    kt_pool_destroy();

    // Edge cases. kt_pool_destroy() above tore down the NT-thread pool, so the
    // next kt_for() lazily builds a fresh one — exercising destroy->recreate and,
    // with nthreads=1, the single-thread pool path.
    w->nthreads = 1;
    g_processed = 0;
    kt_for(count_work, w, BATCH_SIZE);
    if (g_processed != (long)BATCH_SIZE) {
        fprintf(stderr, "FAIL: recreated single-thread pool processed %ld items, expected %d\n",
                g_processed, BATCH_SIZE);
        return 1;
    }
    // Zero work items must be a clean no-op (no work run, no deadlock).
    g_processed = 0;
    kt_for(count_work, w, 0);
    if (g_processed != 0) {
        fprintf(stderr, "FAIL: kt_for(n=0) processed %ld items, expected 0 (no-op)\n", g_processed);
        return 1;
    }

    kt_pool_destroy();
    free(w);

    int fail = 0;
    for (int t = 0; t < NT; ++t) {
        if (!g_seen_ok[0][t] || !g_seen_ok[1][t]) {
            fprintf(stderr, "slot %d was not observed in both kt_for() calls "
                    "-> cannot confirm thread persistence\n", t);
            fail = 1;
            continue;
        }
        if (!pthread_equal(g_seen[0][t], g_seen[1][t])) {
            fprintf(stderr,
                    "slot %d ran on different OS threads across kt_for() calls "
                    "-> pool is not persistent (cross-thread realloc would occur)\n", t);
            fail = 1;
        }
    }
    if (fail) {
        fprintf(stderr, "FAIL: kt_for() is not reusing a persistent worker pool\n");
        return 1;
    }
    printf("OK: kt_for() reuses a persistent worker pool "
           "(same OS thread services each slot across calls)\n");
    return 0;
}
