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
#include <stdio.h>

#if AFF && (__linux__)
extern int affy[256];
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

extern uint64_t tprof[LIM_R][LIM_C];

static inline long steal_work(kt_for_t *t)
{
	int i, min_i = -1;
	long k, min = LONG_MAX;
	for (i = 0; i < t->n_threads; ++i)
		if (min > t->w[i].i) min = t->w[i].i, min_i = i;
	k = __sync_fetch_and_add(&t->w[min_i].i, t->n_threads);
	// return k >= t->n? -1 : k;
	return k*BATCH_SIZE >= t->n? -1 : k;
}

/******** Current working code *********/
static void *ktf_worker(void *data)
{
	ktf_worker_t *w = (ktf_worker_t*)data;
	long i;
	int tid = w->i;

#if AFF && (__linux__)
	fprintf(stderr, "i: %d, CPU: %d\n", tid , sched_getcpu());
#endif
	
	for (;;) {
		i = __sync_fetch_and_add(&w->i, w->t->n_threads);
		int st = i * BATCH_SIZE;
		if (st >= w->t->n) break;
		int ed = (i + 1) * BATCH_SIZE < w->t->n? (i + 1) * BATCH_SIZE : w->t->n;
		// w->t->func(w->t->data, st, ed-st, tid);
        w->t->func(w->t->data, st, ed-st, w - w->t->w);
	}

	while ((i = steal_work(w->t)) >= 0) {
		int st = i * BATCH_SIZE;
		int ed = (i + 1) * BATCH_SIZE < w->t->n? (i + 1) * BATCH_SIZE : w->t->n;
		w->t->func(w->t->data, st, ed-st, w - w->t->w);
	}
	pthread_exit(0);
}

void kt_for(void (*func)(void*, int, int, int), void *data, int n)
{
	int i;
	kt_for_t t;
	pthread_t *tid;
	worker_t *w = (worker_t*) data;
	t.func = func, t.data = data, t.n_threads = w->nthreads, t.n = n;
	t.w = (ktf_worker_t*) malloc (t.n_threads * sizeof(ktf_worker_t));
    assert(t.w != NULL);
	tid = (pthread_t*) malloc (t.n_threads * sizeof(pthread_t));
    assert(tid != NULL);
	for (i = 0; i < t.n_threads; ++i)
		t.w[i].t = &t, t.w[i].i = i;

	pthread_attr_t attr;
    pthread_attr_init(&attr);

#ifdef __APPLE__
    /* Set QoS to USER_INITIATED for compute-intensive alignment work
     * This hints to the scheduler to prefer P-cores (performance cores)
     * over E-cores (efficiency cores) on Apple Silicon */
    static int pcore_count = -2;  /* -2 = not yet queried */
    if (pcore_count == -2) {
        pcore_count = get_pcore_count();
    }

    /* Only set QoS on Apple Silicon (pcore_count > 0) */
    if (pcore_count > 0) {
        pthread_attr_set_qos_class_np(&attr, QOS_CLASS_USER_INITIATED, 0);
    }
#endif

	// printf("getcpu: %d\n", sched_getcpu());
	for (i = 0; i < t.n_threads; ++i) {
#if AFF && (__linux__)
		cpu_set_t cpus;
		CPU_ZERO(&cpus);
		// CPU_SET(i, &cpus);
		CPU_SET(affy[i], &cpus);
		pthread_attr_setaffinity_np(&attr, sizeof(cpu_set_t), &cpus);
		pthread_create(&tid[i], &attr, ktf_worker, &t.w[i]);
#elif defined(__APPLE__)
        /* Use attr with QoS settings on Apple */
        pthread_create(&tid[i], &attr, ktf_worker, &t.w[i]);
#else
		pthread_create(&tid[i], NULL, ktf_worker, &t.w[i]);
#endif
	}
	for (i = 0; i < t.n_threads; ++i) pthread_join(tid[i], 0);

    pthread_attr_destroy(&attr);
    free(t.w);
	free(tid);
}
