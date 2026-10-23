/*************************************************************************************
                           The MIT License

   Copyright Attractive Chaos <attractor@live.co.uk>
   BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.

   Modified Copyright (C) 2019  Intel Corporation, Heng Li.
   Authors: Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>;
         Heng Li <hli@jimmy.harvard.edu>.
*****************************************************************************************/

#include "kthread.h"
#include "stage_prof.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#ifdef BWAMEM3_TESTING
#include <time.h>
#endif


/* Apple Silicon QoS (Quality of Service) support
 * This helps the scheduler preferentially place compute threads on P-cores
 * (performance cores) rather than E-cores (efficiency cores) */
#ifdef __APPLE__
#include <pthread/qos.h>
#include <sys/sysctl.h>

/* Get the number of performance cores on Apple Silicon
 * Returns -1 if unable to detect (e.g., on Intel Macs) */
static int get_pcore_count() {
    int pcore_count = 0;
    size_t size = sizeof(pcore_count);

    /* Try Apple Silicon specific sysctl first */
    if (sysctlbyname("hw.perflevel0.physicalcpu", &pcore_count, &size, NULL, 0) == 0) {
        return pcore_count;
    }

    /* Fallback: assume half of physical cores are P-cores on Apple Silicon,
     * or return -1 on Intel (no hybrid cores) */
    int total_cores = 0;
    size = sizeof(total_cores);
    if (sysctlbyname("hw.physicalcpu", &total_cores, &size, NULL, 0) == 0) {
        /* Check if this is Apple Silicon by looking for perflevel */
        int levels = 0;
        size = sizeof(levels);
        if (sysctlbyname("hw.nperflevels", &levels, &size, NULL, 0) == 0 && levels > 1) {
            /* Apple Silicon with hybrid cores - rough estimate */
            return total_cores / 2;
        }
    }
    return -1;  /* Intel Mac or detection failed */
}
#endif /* __APPLE__ */

/* Work-item layout of one kt_for() pass.
 *
 * A pass ends when its slowest worker finishes its last item, and every other
 * worker idles at that barrier meanwhile. With every item `batch` reads, a
 * worker that picks up a full item just before the others run dry holds the
 * whole pool for up to one item's time (tens of ms on hg38 reads). So about the
 * last n_threads * batch reads (rounded up to whole batches) -- roughly the
 * pool's final round -- are cut into items of batch / KT_TAIL_DIV reads, which
 * lets the workers finish closer together. Items before the tail are exactly
 * the legacy batch-sized items.
 *
 * This makes the item boundaries depend on n_threads (-t). Before it they did
 * not: every item was [k * batch, (k + 1) * batch) whatever -t, so at a fixed
 * -K the reads that shared an item were the same at every thread count. Output
 * stays identical because no callback's per-read result depends on which other
 * reads share its item:
 *  - every kt_for callback computes per-read results into per-read slots, and
 *    the per-tid scratch it borrows is sized for a full batch (tail items are
 *    smaller);
 *  - item boundaries stay even (batch and tail_batch are both even), so a read
 *    pair never straddles two items;
 *  - the batched SIMD kernels that do group reads within an item return each
 *    pair's result independently of its lane group: the banded extension
 *    kernels since PR #565 (test_bandedswa_lane_indep.cpp), and the kswv
 *    mate-rescue kernels, which match the scalar ksw_align2 on every field for
 *    every pair of mixed-length and partial batches (test_kswv_correctness.cpp,
 *    test_kswv_mixed_qlen.cpp) and whose banded rescue path is exact per job
 *    (rescue_band.h). The 8-bit vs 16-bit choice is made per pair in both.
 * Measured: unsorted output at -t 1, 7 and 16 with a fixed -K was identical on
 * 5M-pair WGS and WES slices; test/regression/thread_determinism.sh compares
 * -t 1, 3 and 4 at a fixed -K.
 *
 * Items are numbered in read order and tile [0, n) contiguously; the strided
 * dispatch and work stealing in ktf_run() are unchanged, they just count items
 * from the layout instead of dividing by BATCH_SIZE. */
#define KT_TAIL_DIV 4
static_assert(BATCH_SIZE % (2 * KT_TAIL_DIV) == 0,
              "kt_for tail items must be an even number of reads (whole pairs)");

void kt_item_layout(long n, int n_threads, long batch, kt_item_layout_t *L)
{
	xassert(batch > 0 && batch % (2 * KT_TAIL_DIV) == 0,
	        "kt_item_layout: batch must be a positive multiple of 2 * KT_TAIL_DIV");
	if (n < 0) n = 0;
	if (n_threads < 1) n_threads = 1;
	L->n = n;
	L->batch = batch;
	L->tail_batch = batch / KT_TAIL_DIV;
	const long tail = (long)n_threads * batch;
	L->full_end = n > tail ? ((n - tail) / batch) * batch : 0;
	L->n_full = L->full_end / batch;
	L->n_items = L->n_full + (n - L->full_end + L->tail_batch - 1) / L->tail_batch;
}

int kt_item_range(const kt_item_layout_t *L, long k, int *st, int *ed)
{
	if (k < 0 || k >= L->n_items) return 0;
	long s, e;
	if (k < L->n_full) { s = k * L->batch; e = s + L->batch; }
	else { s = L->full_end + (k - L->n_full) * L->tail_batch; e = s + L->tail_batch; }
	if (e > L->n) e = L->n;
	*st = (int)s; *ed = (int)e;
	return 1;
}

static inline long steal_work(kt_for_t *t)
{
	int i, min_i = -1;
	long k, min = LONG_MAX;
	for (i = 0; i < t->n_threads; ++i)
		if (min > t->w[i].i) min = t->w[i].i, min_i = i;
	k = __sync_fetch_and_add(&t->w[min_i].i, t->n_threads);
	return k >= t->layout.n_items? -1 : k;
}

static void kt_item_finished(kt_for_t *t, long k);
#ifdef BWAMEM3_TESTING
static int kt_testing_hold_last_item(void);
static void kt_testing_wait_for_release(void);
#endif

/* The per-thread work loop: strided dispatch followed by work-stealing.
 * Extracted from the old per-call thread body so it can be driven by a
 * persistent pool worker (see below). Items come from the pass's layout
 * (kt_item_layout); tid passed to func() is the worker's slot index
 * (w - w->t->w). */
static void ktf_run(ktf_worker_t *w)
{
	long i;
	double _c0 = sp_enabled() ? sp_thread_cpu() : 0.0;
	if (sp_enabled()) sp_encode_reset();


	int st, ed;
	for (;;) {
		i = __sync_fetch_and_add(&w->i, w->t->n_threads);
		if (!kt_item_range(&w->t->layout, i, &st, &ed)) break;
		w->t->func(w->t->data, st, ed-st, w - w->t->w);
		if (w->t->done) {
			kt_item_finished(w->t, i);
#ifdef BWAMEM3_TESTING
			if (i == w->t->layout.n_items - 1 && w->t->layout.n_items >= 2 && kt_testing_hold_last_item())
				kt_testing_wait_for_release();
#endif
		}
	}

	while ((i = steal_work(w->t)) >= 0) {
		if (!kt_item_range(&w->t->layout, i, &st, &ed)) break;
		w->t->func(w->t->data, st, ed-st, w - w->t->w);
		if (w->t->done) {
			kt_item_finished(w->t, i);
#ifdef BWAMEM3_TESTING
			if (i == w->t->layout.n_items - 1 && w->t->layout.n_items >= 2 && kt_testing_hold_last_item())
				kt_testing_wait_for_release();
#endif
		}
	}
	if (sp_enabled()) { w->cpu_busy = sp_thread_cpu() - _c0; w->encode = sp_encode_get(); }
}

/* ---------------------------------------------------------------------------
 * Persistent worker pool for kt_for().
 *
 * The original kt_for() spawned n worker pthreads on every call and joined
 * them at the end. bwa-mem3 calls kt_for() three times per chunk (worker_bwt /
 * worker_aln / worker_sam) across many chunks, so a per-thread scratch buffer
 * (e.g. mmc->enc_qdb[tid]) is malloc()'d by one chunk's worker thread and then
 * realloc()'d by a *different* OS thread on a later chunk — the first thread
 * has already exited. That cross-thread realloc of an abandoned-heap block is
 * mishandled by mimalloc v3.x and silently corrupts the heap, crashing with
 * SIGSEGV several chunks later (it is benign under glibc / mimalloc v2.x, which
 * is why no sanitizer flags it).
 *
 * A persistent pool keeps the same n worker threads alive for the whole run, so
 * worker slot `tid` is always the same OS thread: every realloc of that slot's
 * buffers happens on the thread that allocated them — no cross-thread realloc.
 * It also removes the per-chunk thread create/join overhead.
 * ------------------------------------------------------------------------- */
typedef struct {
	int started;
	int n_threads;
	pthread_t *threads;
	ktf_worker_t *w;          /* per-worker steal state (n_threads entries) */
	kt_for_t job;             /* current job: func / data / n / w           */
	pthread_mutex_t mtx;
	pthread_cond_t cv_go;     /* signalled when a new generation is posted  */
	pthread_cond_t cv_done;   /* signalled when the last worker finishes    */
	long generation;          /* bumped once per kt_for() call              */
	int n_left;               /* workers not yet done with `generation`     */
	int shutdown;
	long n_finished;          /* kt_for_prefix items finished, all passes   */
#ifdef BWAMEM3_TESTING
	int hold_released;        /* test hold (kt_testing_hold_last_item)      */
	pthread_cond_t cv_hold;   /* signalled when hold_released is set        */
#endif
} kt_pool_t;

static kt_pool_t g_kt_pool = {0};

#ifdef BWAMEM3_TESTING
/* Test-only (TESTING_BUILD=1): BWAMEM3_TESTING_HOLD_LAST_ITEM=1 keeps the worker
 * that finished a kt_for_prefix pass's last item from leaving the pass until the
 * caller has consumed a prefix, and has a declined caller re-offer on a 1 ms
 * timed wait. Every pass of two or more items then streams at least one, so a
 * test can assert streaming engages without depending on whether the caller
 * runs before every item has finished (test/regression/stream_output_order.sh).
 * A pass whose prefixes are never consumed is released after 30 s with a
 * warning instead of hanging. Read once. */
static int kt_testing_hold_last_item(void)
{
	static int hold = -1;
	if (hold < 0) {
		const char *e = getenv("BWAMEM3_TESTING_HOLD_LAST_ITEM");
		hold = (e && e[0] == '1' && e[1] == '\0');
	}
	return hold;
}

static void kt_testing_wait_for_release(void)
{
	struct timespec dl;
	clock_gettime(CLOCK_REALTIME, &dl);
	dl.tv_sec += 30;
	pthread_mutex_lock(&g_kt_pool.mtx);
	while (!g_kt_pool.hold_released)
		if (pthread_cond_timedwait(&g_kt_pool.cv_hold, &g_kt_pool.mtx, &dl) != 0) {
			fprintf(stderr, "WARNING: BWAMEM3_TESTING_HOLD_LAST_ITEM: no prefix consumed in 30 s; releasing\n");
			g_kt_pool.hold_released = 1;
		}
	pthread_mutex_unlock(&g_kt_pool.mtx);
}
#endif

/* kt_for_prefix: item k of the current pass has returned from func. The store
 * releases everything func wrote for the item's reads to the caller, which
 * reads done[k] with acquire before handing those reads to on_prefix. cv_done
 * has one waiter (the kt_for_prefix caller), so a signal wakes it. */
static void kt_item_finished(kt_for_t *t, long k)
{
	__atomic_store_n(&t->done[k], 1, __ATOMIC_RELEASE);
	pthread_mutex_lock(&g_kt_pool.mtx);
	++g_kt_pool.n_finished;
	pthread_cond_signal(&g_kt_pool.cv_done);
	pthread_mutex_unlock(&g_kt_pool.mtx);
}

static void *kt_pool_worker(void *arg)
{
	long k = (long)(intptr_t)arg;
	long seen = 0;
	for (;;) {
		pthread_mutex_lock(&g_kt_pool.mtx);
		while (g_kt_pool.generation == seen && !g_kt_pool.shutdown)
			pthread_cond_wait(&g_kt_pool.cv_go, &g_kt_pool.mtx);
		if (g_kt_pool.shutdown) { pthread_mutex_unlock(&g_kt_pool.mtx); break; }
		seen = g_kt_pool.generation;
		pthread_mutex_unlock(&g_kt_pool.mtx);

		ktf_run(&g_kt_pool.w[k]);

		pthread_mutex_lock(&g_kt_pool.mtx);
		if (--g_kt_pool.n_left == 0) pthread_cond_signal(&g_kt_pool.cv_done);
		pthread_mutex_unlock(&g_kt_pool.mtx);
	}
	return NULL;
}

static void kt_pool_init(int n_threads)
{
	g_kt_pool.n_threads = n_threads;
	g_kt_pool.threads = (pthread_t*) malloc(n_threads * sizeof(pthread_t));
	/* L27: cache-line-aligned so each worker's steal counter sits on its own
	 * line (sizeof(ktf_worker_t) is padded to 64, so n_threads*sizeof is a
	 * multiple of 64 as aligned_alloc requires). free() handles it. */
	g_kt_pool.w       = (ktf_worker_t*) aligned_alloc(64, n_threads * sizeof(ktf_worker_t));
	if (g_kt_pool.threads == NULL || g_kt_pool.w == NULL) {
		perror("Allocation of kt_for worker pool failed");
		exit(EXIT_FAILURE);
	}
	pthread_mutex_init(&g_kt_pool.mtx, NULL);
	pthread_cond_init(&g_kt_pool.cv_go, NULL);
	pthread_cond_init(&g_kt_pool.cv_done, NULL);
#ifdef BWAMEM3_TESTING
	pthread_cond_init(&g_kt_pool.cv_hold, NULL);
#endif
	g_kt_pool.generation = 0;
	g_kt_pool.n_left = 0;
	g_kt_pool.shutdown = 0;

	pthread_attr_t attr;
	pthread_attr_init(&attr);
#ifdef __APPLE__
	/* Prefer P-cores on Apple Silicon for compute-intensive alignment work. */
	int pcore_count = get_pcore_count();
	if (pcore_count > 0)
		pthread_attr_set_qos_class_np(&attr, QOS_CLASS_USER_INITIATED, 0);
#endif
	for (int i = 0; i < n_threads; ++i) {
		/* A failed worker would leave kt_for() waiting on cv_done for a
		 * completion that never arrives (deadlock), so fail loudly instead.
		 * pthread_create() returns the error code directly rather than via
		 * errno, so report it with strerror(rc), not perror(). */
		int rc = pthread_create(&g_kt_pool.threads[i], &attr, kt_pool_worker, (void*)(intptr_t)i);
		if (rc != 0) {
			fprintf(stderr, "ERROR: kt_for worker pool: pthread_create failed (worker %d): %s\n", i, strerror(rc));
			exit(EXIT_FAILURE);
		}
	}
	pthread_attr_destroy(&attr);
	g_kt_pool.started = 1;
}

/* Tear the pool down (joins the workers). Safe to call when no pool exists. */
void kt_pool_destroy(void)
{
	if (!g_kt_pool.started) return;
	pthread_mutex_lock(&g_kt_pool.mtx);
	g_kt_pool.shutdown = 1;
	pthread_cond_broadcast(&g_kt_pool.cv_go);
	pthread_mutex_unlock(&g_kt_pool.mtx);
	for (int i = 0; i < g_kt_pool.n_threads; ++i)
		pthread_join(g_kt_pool.threads[i], NULL);
	free(g_kt_pool.threads);
	free(g_kt_pool.w);
	pthread_mutex_destroy(&g_kt_pool.mtx);
	pthread_cond_destroy(&g_kt_pool.cv_go);
	pthread_cond_destroy(&g_kt_pool.cv_done);
#ifdef BWAMEM3_TESTING
	pthread_cond_destroy(&g_kt_pool.cv_hold);
#endif
	memset(&g_kt_pool, 0, sizeof(g_kt_pool));
}

void kt_for(void (*func)(void*, int, int, int), void *data, int n)
{
	kt_for_nt(func, data, n, ((worker_t*) data)->nthreads);
}

/* Create the pool on first use, or rebuild it for a new size. kt_for() is only
 * ever called from the serialized step-1 of kt_pipeline (one chunk in flight at
 * a time), so lazy init here is race-free, and n_threads (= opt->n_threads) is
 * constant for the whole run. Rebuild the pool defensively if that invariant
 * ever changes. */
static void kt_pool_ensure(int n_threads)
{
	if (!g_kt_pool.started)
		kt_pool_init(n_threads);
	else if (g_kt_pool.n_threads != n_threads) {
		kt_pool_destroy();
		kt_pool_init(n_threads);
	}
}

/* Post one pass, cut into items by L (kt_item_layout), to the pool. Caller
 * holds g_kt_pool.mtx and then waits on cv_done until n_left reaches 0. */
static void kt_pool_post(void (*func)(void*, int, int, int), void *data, const kt_item_layout_t *L,
                         int n_threads, uint8_t *done)
{
	g_kt_pool.job.func = func;
	g_kt_pool.job.data = data;
	g_kt_pool.job.n_threads = n_threads;
	g_kt_pool.job.w = g_kt_pool.w;
	g_kt_pool.job.done = done;
	g_kt_pool.job.layout = *L;
	for (int i = 0; i < n_threads; ++i) {
		g_kt_pool.w[i].t = &g_kt_pool.job;
		g_kt_pool.w[i].i = i;
	}
	g_kt_pool.n_left = n_threads;
	++g_kt_pool.generation;
	pthread_cond_broadcast(&g_kt_pool.cv_go);
}

/* After a pass: all workers are done with this generation (n_left == 0) and
 * will not touch their per-slot stats again until the next kt_for() call, so
 * the pool's worker slots can be read without the lock here. */
static void kt_pool_harvest(int n_threads)
{
	if (sp_enabled()) {
		double *busy = (double*) malloc(n_threads * sizeof(double));
		xassert(busy != NULL, "out of memory: busy");
		double sum = 0, esum = 0;
		for (int i = 0; i < n_threads; ++i) { busy[i] = g_kt_pool.w[i].cpu_busy; sum += busy[i]; esum += g_kt_pool.w[i].encode; }
		g_ktfor.proc_cpu += sum;                        /* accumulate across kt_for calls in a step */
		g_ktfor.encode   += esum;                       /* SAM/BAM-build CPU (only worker_sam adds) */
		sp_thread_stats(&g_ktfor, busy, n_threads);     /* balance stats from the most recent call */
		free(busy);
	}
}

void kt_for_nt(void (*func)(void*, int, int, int), void *data, int n, int n_threads)
{
	if (n_threads < 1) n_threads = 1;  /* nthreads is int16_t; never spin up a 0/negative pool */
	kt_pool_ensure(n_threads);

	kt_item_layout_t L;
	kt_item_layout(n, n_threads, BATCH_SIZE, &L);
	pthread_mutex_lock(&g_kt_pool.mtx);
	kt_pool_post(func, data, &L, n_threads, NULL);
	while (g_kt_pool.n_left > 0)
		pthread_cond_wait(&g_kt_pool.cv_done, &g_kt_pool.mtx);
	pthread_mutex_unlock(&g_kt_pool.mtx);

	kt_pool_harvest(n_threads);
}

void kt_for_prefix(void (*func)(void*, int, int, int), void *data, int n, int n_threads,
                   kt_prefix_f on_prefix, void *ctx)
{
	if (n_threads < 1) n_threads = 1;
	kt_pool_ensure(n_threads);

	kt_item_layout_t L;
	kt_item_layout(n, n_threads, BATCH_SIZE, &L);
	uint8_t *done = (uint8_t*) calloc(L.n_items > 0 ? L.n_items : 1, 1);
	xassert(done != NULL, "out of memory: kt_for_prefix done flags");

	pthread_mutex_lock(&g_kt_pool.mtx);
#ifdef BWAMEM3_TESTING
	const int hold = kt_testing_hold_last_item();
	g_kt_pool.hold_released = 0;
#endif
	kt_pool_post(func, data, &L, n_threads, done);
	/* next: the first item not yet consumed. One item per on_prefix call, so
	 * once the pass ends the caller is held up by at most the item it is
	 * handling. A declined item is offered again after the next item
	 * completion (or not at all, if none is left); one that lands while
	 * on_prefix runs counts, so the decline does not wait for another. */
	long next = 0;
	int declined = 0;
	while (g_kt_pool.n_left > 0) {
		if (!declined && next < L.n_items && __atomic_load_n(&done[next], __ATOMIC_ACQUIRE)) {
			const long finished = g_kt_pool.n_finished;
			int st, ed;
			kt_item_range(&L, next, &st, &ed);
			pthread_mutex_unlock(&g_kt_pool.mtx);
			const int consumed = on_prefix(ctx, st, ed);
			pthread_mutex_lock(&g_kt_pool.mtx);
			if (consumed) {
				++next;
#ifdef BWAMEM3_TESTING
				if (hold && !g_kt_pool.hold_released) {
					g_kt_pool.hold_released = 1;
					pthread_cond_broadcast(&g_kt_pool.cv_hold);
				}
#endif
			}
			else if (g_kt_pool.n_finished == finished) declined = 1;
			continue;
		}
#ifdef BWAMEM3_TESTING
		if (hold && declined) {
			/* No item completion may come while the last item is held, so
			 * re-offer on a timer: the decline ends when its turn comes. */
			struct timespec dl;
			clock_gettime(CLOCK_REALTIME, &dl);
			dl.tv_nsec += 1000000;
			if (dl.tv_nsec >= 1000000000) { dl.tv_sec += 1; dl.tv_nsec -= 1000000000; }
			pthread_cond_timedwait(&g_kt_pool.cv_done, &g_kt_pool.mtx, &dl);
			declined = 0;
			continue;
		}
#endif
		pthread_cond_wait(&g_kt_pool.cv_done, &g_kt_pool.mtx);
		declined = 0;
	}
	g_kt_pool.job.done = NULL;
	pthread_mutex_unlock(&g_kt_pool.mtx);
	free(done);

	kt_pool_harvest(n_threads);
}
