/*************************************************************************************
                           The MIT License

   BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),
   Copyright (C) 2019  Intel Corporation, Heng Li.

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

Authors: Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>.
*****************************************************************************************/

#include <stdio.h>
#include "macro.h"
#include <stdint.h>
#include <assert.h>
#include <math.h>       /* NAN: proc_freq_slept_s */
#include <time.h>       /* clock_gettime, nanosleep: proc_freq */
#include "bwa.h"        /* bwa_verbose: gates the per-category UGP histograms */
#include "profiling.h"
#if defined(__ARM_NEON) || defined(__aarch64__)
#include "simd_compat.h" /* __rdtsc: sse2neon's _rdtsc (CNTVCT_EL0) */
#else
#include <immintrin.h>   /* __rdtsc (GCC < 11: the src/utils.h shim, via bwa.h) */
#endif

// The clock the tick counter is measured against when the rate has to be
// measured: on x86, and on aarch64 only if CNTFRQ_EL0 reads 0 (see below).
// CLOCK_MONOTONIC_RAW is preferred because it is not slewed by NTP.
#if defined(CLOCK_MONOTONIC_RAW)
#define PROC_FREQ_CLOCK CLOCK_MONOTONIC_RAW
#else
#define PROC_FREQ_CLOCK CLOCK_MONOTONIC
#endif

// Over a whole-run window the clock and the counter must also agree about
// system sleep, and on x86 they need not: what the TSC does across S3 or
// s2idle varies by platform (it may keep counting, stop, or reset), while on
// Linux CLOCK_MONOTONIC_RAW excludes suspend and on an Intel Mac it includes
// it. A window that spans a suspend is therefore not used. It is detected with
// a pair of clocks that differ only in whether they count suspend: on Linux
// CLOCK_BOOTTIME counts it and CLOCK_MONOTONIC does not; on macOS
// CLOCK_MONOTONIC_RAW counts it and CLOCK_UPTIME_RAW does not (clock_gettime(3):
// CLOCK_UPTIME_RAW "does not increment while the system is asleep"). If their
// difference grew across the window, the system slept and the rate is measured
// over the sleep calibration instead. Where neither pair exists the window is
// used unchecked; a counter that did not advance (c1 <= c0) still falls back.
// The cost of getting it wrong is only the timing lines on stderr, never the
// SAM output.
#if defined(__linux__) && defined(CLOCK_BOOTTIME)
#define PROC_FREQ_AWAKE_CLOCK CLOCK_MONOTONIC
#define PROC_FREQ_TOTAL_CLOCK CLOCK_BOOTTIME
#elif defined(__APPLE__) && defined(CLOCK_UPTIME_RAW) && defined(CLOCK_MONOTONIC_RAW)
#define PROC_FREQ_AWAKE_CLOCK CLOCK_UPTIME_RAW
#define PROC_FREQ_TOTAL_CLOCK CLOCK_MONOTONIC_RAW
#endif

// A window whose suspend clocks drifted apart by more than this slept. Far
// above the microseconds between the two clock reads; a preemption between
// them long enough to cross it costs only a 5 ms sleep calibration.
static const double proc_freq_suspend_floor_s = 1e-3;

// Seconds from a to b.
static double timespec_elapsed_s(const struct timespec *a, const struct timespec *b)
{
    return (double)(b->tv_sec - a->tv_sec) + (double)(b->tv_nsec - a->tv_nsec) * 1e-9;
}

// Profiling globals. Defined here rather than in main.cpp so they are part of
// libbwa.a and visible to consumers that link the library without pulling in
// main.o (e.g. language bindings).
uint64_t proc_freq = 0;
alignas(64) tprof_row_t tprof[LIM_R] = {};
uint64_t prof[LIM_R] = {0};

// On aarch64, __rdtsc() is sse2neon's _rdtsc(), which reads the generic
// timer's virtual count CNTVCT_EL0, and CNTFRQ_EL0 holds that counter's
// frequency. The rate is read rather than measured: exact, with no window
// and no sleep. Firmware programs CNTFRQ_EL0 and a nonzero value is trusted
// as programmed; if it reads 0 the rate is measured as on x86 instead.
//
// Elsewhere the rate is measured against PROC_FREQ_CLOCK over a window,
// dividing the ticks by the wall time that clock reports for the same window.
// The TSC is invariant on every x86-64 this targets, so the residual error is
// the nanosecond-scale bracket of each (clock, counter) sample (see
// proc_freq_sample()) spread over the window, far below the precision the
// timing reports claim.
//
// The window is the run itself: proc_freq_mark_start() opens it at the start
// of main() and proc_freq_resolve() closes it when a report first reads the
// rate.
// proc_freq is read by nothing during the run, only by those end-of-run
// reports, so no process sleeps before doing its work. The run's window is also
// the longer one: seconds for an alignment rather than milliseconds.
//
// History: this was `tim = __rdtsc(); sleep(1); proc_freq = __rdtsc() - tim;`
// in main(), a full second of wall time on EVERY invocation (run_unit_tests.sh
// spent ~60 s of its 74 s asleep), then a 5 ms nanosleep() window, still paid
// up front by every process. A sleeping core sits at its lowest P-state, so on
// parts whose TSC is *not* invariant a sleep measured precisely the frequency
// the aligner never runs at, where the run's own window averages over the
// frequencies the run ran at.

// The window's opening (clock, counter) pair. A zero proc_freq_t0 means the
// window was never opened: a monotonic clock never reads exactly zero after
// boot.
static struct timespec proc_freq_t0 = {0, 0};
static uint64_t proc_freq_c0 = 0;
// The suspend clocks' difference (total minus awake, in seconds) when the
// window opened; NaN when it could not be read, which leaves the window
// unchecked.
static double proc_freq_slept0 = 0.0;

#if defined(__aarch64__)
// 1: read CNTFRQ_EL0 (always, outside the unit test). Set only by the test
// hook below: 0 skips the read to reach the measured paths on aarch64, and -1
// treats the register as reading 0 to pin the fall-through.
static int proc_freq_use_cntfrq = 1;

static uint64_t read_cntfrq(void)
{
    uint64_t freq;
    __asm__ __volatile__("mrs %0, cntfrq_el0" : "=r"(freq));
    return freq;
}
#endif

// Number of times calibrate_proc_freq() has run; read by the test hook only.
static int proc_freq_sleep_calibrations = 0;
// The source proc_freq_resolved() last named, NULL until it names one, so a
// pre-set rate is named once; also read by the test hook.
static const char *proc_freq_source = NULL;
// Test hook state: when positive, the PROC_FREQ_CLOCK read with this 1-based
// index (counted from the last reset) fails as clock_gettime() would.
static int proc_freq_fail_read = 0;
static int proc_freq_reads = 0;

static int proc_freq_clock(struct timespec *t)
{
    if (++proc_freq_reads == proc_freq_fail_read) return -1;
    return clock_gettime(PROC_FREQ_CLOCK, t);
}

// Read the clock and the counter as one pair. The counter is read on both
// sides of the clock and the pair keeps the midpoint of the two counts, so the
// pair's error is at most half the bracket. A thread preempted between its
// reads widens the bracket, so of a few tries the tightest is kept: a split
// pair is discarded rather than biasing the rate, in either direction. 0 when
// a clock read failed, and t and c are then not to be used.
static int proc_freq_sample(struct timespec *t, uint64_t *c)
{
    const int tries = 4;
    uint64_t best = UINT64_MAX;
    for (int i = 0; i < tries; ++i) {
        struct timespec ti = {0, 0};
        uint64_t a = __rdtsc();
        if (proc_freq_clock(&ti) != 0) return 0;
        uint64_t b = __rdtsc();
        uint64_t span = b >= a ? b - a : UINT64_MAX - 1;
        if (span < best) {
            best = span;
            *t = ti;
            *c = a + (b >= a ? span / 2 : 0);
        }
    }
    return 1;
}

// Seconds the system has spent suspended since boot, up to a constant: the
// suspend clocks' difference. NaN where there is no such pair or a read failed.
static double proc_freq_slept_s(void)
{
#if defined(PROC_FREQ_AWAKE_CLOCK)
    struct timespec awake = {0, 0}, total = {0, 0};
    if (clock_gettime(PROC_FREQ_AWAKE_CLOCK, &awake) == 0
        && clock_gettime(PROC_FREQ_TOTAL_CLOCK, &total) == 0)
        return timespec_elapsed_s(&awake, &total);
#endif
    return NAN;
}

void proc_freq_mark_start(void)
{
    struct timespec t0 = {0, 0};
    uint64_t c0 = 0;
    if (!proc_freq_sample(&t0, &c0)) return;
    proc_freq_slept0 = proc_freq_slept_s();
    proc_freq_c0 = c0;
    proc_freq_t0 = t0;
}

// The fallback, taken by proc_freq_resolve() when the run's window was never
// opened (a library consumer that does not go through main(), or
// proc_freq_mark_start()'s clock_gettime() failed), is under 1 ms, cannot be
// closed (clock_gettime() failed), spans a system suspend, or saw a counter
// that did not advance (c1 <= c0); and on aarch64 only when CNTFRQ_EL0 also
// reads 0. Measure over a 5 ms sleep, doubling per retry. nanosleep() returns
// early when a signal arrives; we divide by the *measured* elapsed time rather
// than the requested one, so a short window is still usable -- but if it came
// back too short, or either clock read failed, try again on a longer one.
static uint64_t calibrate_proc_freq(void)
{
    ++proc_freq_sleep_calibrations;
    const long base_window_ns = 5L * 1000 * 1000;
    const int max_attempts = 4;
    for (int attempt = 0; attempt < max_attempts; ++attempt) {
        struct timespec window = {0, base_window_ns << attempt};
        struct timespec t0 = {0, 0}, t1 = {0, 0};
        uint64_t c0 = 0, c1 = 0;

        if (!proc_freq_sample(&t0, &c0)) continue;
        nanosleep(&window, NULL);
        if (!proc_freq_sample(&t1, &c1)) continue;

        double elapsed = timespec_elapsed_s(&t0, &t1);
        if (elapsed >= 1e-4 && c1 > c0)
            return (uint64_t)((double)(c1 - c0) / elapsed + 0.5);
    }
    // Not reachable in practice (it would take max_attempts consecutive
    // windows cut short to under 100 us, or clock reads failing). Say so on
    // stderr rather than silently returning a bogus rate: without this, a
    // pathological environment -- nanosleep() blocked by a seccomp filter, a
    // process under relentless signal pressure -- would print implausible
    // timings with nothing explaining why. stderr only, so the SAM on stdout
    // is unaffected.
    fprintf(stderr, "[W::%s] could not calibrate the processor tick rate in "
                    "%d attempts; timing reports will be meaningless.\n",
            __func__, max_attempts);
    // Return 1 rather than 0 so the reports print implausible seconds
    // instead of dividing by zero.
    return 1;
}

// At -v 4, name the source the rate came from, so a test can tell a run that
// measured over its window from one that slept for the rate. Named once: a
// later report that reuses the rate prints nothing.
static uint64_t proc_freq_resolved(const char *source)
{
    proc_freq_source = source;
    if (bwa_verbose >= 4)
        fprintf(stderr, "[D::proc_freq_resolve] tick rate %.6f MHz from %s\n",
                proc_freq * 1.0 / 1e6, source);
    return proc_freq;
}

uint64_t proc_freq_resolve(void)
{
    if (proc_freq != 0)
        return proc_freq_source ? proc_freq : proc_freq_resolved("the caller");
#if defined(__aarch64__)
    if (proc_freq_use_cntfrq) {
        proc_freq = proc_freq_use_cntfrq > 0 ? read_cntfrq() : 0;
        if (proc_freq != 0) return proc_freq_resolved("CNTFRQ_EL0");
    }
#endif
    if (proc_freq_t0.tv_sec != 0 || proc_freq_t0.tv_nsec != 0) {
        struct timespec t1 = {0, 0};
        uint64_t c1 = 0;
        if (proc_freq_sample(&t1, &c1)) {
            double elapsed = timespec_elapsed_s(&proc_freq_t0, &t1);
            // A NaN on either side (no suspend clocks) compares false: unchecked.
            bool slept = proc_freq_slept_s() - proc_freq_slept0 > proc_freq_suspend_floor_s;
            // The window is used once it reaches 1 ms, where the tens of
            // nanoseconds spent reading the clock and the counter are a few
            // hundredths of a percent of it. A shorter one (e.g. a usage
            // error) is measured over the sleep instead, which at 5 ms is
            // the more precise of the two. The floor is not raised to 5 ms
            // to match: every x86 run shorter than that (a tiny `index` or
            // `mem`) would then sleep for its rate again, which is what
            // measuring over the run exists to avoid.
            if (elapsed >= 1e-3 && c1 > proc_freq_c0 && !slept) {
                proc_freq = (uint64_t)((double)(c1 - proc_freq_c0) / elapsed + 0.5);
                return proc_freq_resolved("the run window");
            }
        }
    }
    proc_freq = calibrate_proc_freq();
    return proc_freq_resolved("a sleep calibration");
}

// Test hooks for test/unit/test_proc_freq.cpp, which declares them itself.
// Reset every piece of proc_freq state to how a fresh process starts, and say
// whether the aarch64 path reads CNTFRQ_EL0 (1), skips it (0), or treats it as
// reading 0 (-1); ignored elsewhere. Each branch of proc_freq_resolve() can
// then be reached in one process.
void bwamem3_proc_freq_reset_for_test(int use_counter_register)
{
    proc_freq = 0;
    proc_freq_t0.tv_sec = 0;
    proc_freq_t0.tv_nsec = 0;
    proc_freq_c0 = 0;
    proc_freq_slept0 = 0.0;
    proc_freq_sleep_calibrations = 0;
    proc_freq_source = NULL;
    proc_freq_fail_read = 0;
    proc_freq_reads = 0;
#if defined(__aarch64__)
    proc_freq_use_cntfrq = use_counter_register;
#else
    (void)use_counter_register;
#endif
}

// Overwrite the window's opening counter value, so a test can put it ahead
// of the counter and reach the c1 <= c0 fallback.
void bwamem3_proc_freq_set_c0_for_test(uint64_t c0)
{
    proc_freq_c0 = c0;
}

// How many times the sleep calibration has run since the last reset.
int bwamem3_proc_freq_sleep_calibrations(void)
{
    return proc_freq_sleep_calibrations;
}

// The source the last resolve named ("the caller", "CNTFRQ_EL0", "the run
// window", "a sleep calibration"), or NULL if none since the last reset.
const char *bwamem3_proc_freq_source_for_test(void)
{
    return proc_freq_source;
}

// Make the n-th PROC_FREQ_CLOCK read after the last reset fail (1-based; 0
// disables), so a test can reach the clock_gettime() failure paths. Each
// (clock, counter) sample makes four reads, fewer if one fails.
void bwamem3_proc_freq_fail_clock_read_for_test(int n)
{
    proc_freq_fail_read = n;
}

// Move the opened window's suspend reading back by s seconds, as if the
// system had slept that long since the window opened. 0 where this platform
// has no suspend clocks to detect a sleep with (the hook then does nothing).
int bwamem3_proc_freq_add_suspend_for_test(double s)
{
#if defined(PROC_FREQ_AWAKE_CLOCK)
    proc_freq_slept0 -= s;
    return 1;
#else
    (void)s;
    return 0;
#endif
}

int find_opt(const tprof_row_t &a, int len, uint64_t *max, uint64_t *min, double *avg)
{
    *max = 0;
    *min = 1e15;
    *avg = 0;

    int i=0;
    for (i=0; i<len; i++)
    {
        if (a[i] > *max) *max = a[i];
        if (a[i] < *min) *min = a[i];
        *avg += a[i];
    }
    *avg /= len;

    return 1;
}

int display_stats(int nthreads)
{
    uint64_t max, min;
    double avg;
    proc_freq_resolve();  // sets proc_freq, which every line below divides by
    fprintf(stderr, "No. of OMP threads: %d\n", nthreads);
    fprintf(stderr, "Processor is running @%lf MHz\n", proc_freq*1.0/1e6);
    fprintf(stderr, "Runtime profile:\n");

    /* tprof[][] has a fixed column dimension of LIM_C (= MAX_THREADS = 256).
     * Every per-thread aggregation below — both the explicit `for (t=0; t<nthreads;...)`
     * loops and the find_opt(tprof[X], nthreads, ...) calls — would index past
     * the column dimension when the user runs with -t > LIM_C, walking into
     * adjacent globals. Clamp the per-thread profile span to LIM_C; the actual
     * thread count is already logged above. */
    if (nthreads > LIM_C) {
        fprintf(stderr,
                "WARNING: per-thread profile truncated to %d threads (-t %d > LIM_C=%d)\n",
                LIM_C, nthreads, LIM_C);
        nthreads = LIM_C;
    }

    fprintf(stderr, "\n\tTime taken for main_mem function: %0.2lf sec\n\n",
            tprof[MEM][0]*1.0/proc_freq);

    fprintf(stderr, "\tIO times (sec) :\n");
    find_opt(tprof[READ_IO], 1, &max, &min, &avg);
    fprintf(stderr, "\tReading IO time (reads) avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);

    find_opt(tprof[SAM_IO], 1, &max, &min, &avg);
    fprintf(stderr, "\tWriting IO time (SAM) avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);

    find_opt(tprof[REF_IO], 1, &max, &min, &avg);
    fprintf(stderr, "\tReading IO time (Reference Genome) avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);

    find_opt(tprof[FMI], 1, &max, &min, &avg);
    fprintf(stderr, "\tIndex read time avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);

    fprintf(stderr, "\n\tOverall time (sec) (Excluding Index reading time):\n");
    // find_opt(tprof[PROCESS], 1, &max, &min, &avg);
    fprintf(stderr, "\tPROCESS() (Total compute time + (read + SAM) IO time) : %0.2lf\n",
            tprof[PROCESS][0]*1.0/proc_freq);

    find_opt(tprof[MEM_PROCESS2], 1, &max, &min, &avg);
    fprintf(stderr, "\tMEM_PROCESS_SEQ() (Total compute time (Kernel + SAM)), avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);

    fprintf(stderr, "\n\t SAM Processing time (sec):\n");
    find_opt(tprof[WORKER20], 1, &max, &min, &avg);
    fprintf(stderr, "\t--WORKER_SAM avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
    
#if HIDE
    find_opt(tprof[SAM1], 1, &max, &min, &avg);
    fprintf(stderr, "\t\tWORKER_SAM1 avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
    find_opt(tprof[SAM2], 1, &max, &min, &avg);
    fprintf(stderr, "\t\tWORKER_SAM2 avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
    
    find_opt(tprof[ALIGN1], 1, &max, &min, &avg);
    fprintf(stderr, "\t\t\tWORKER_ALIGN1 avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
    
    find_opt(tprof[SAM3], 1, &max, &min, &avg);
    fprintf(stderr, "\t\tWORKER_SAM3 avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
#endif
    
    fprintf(stderr, "\n\tKernels' compute time (sec):\n");
    find_opt(tprof[WORKER10], 1, &max, &min, &avg);
    fprintf(stderr, "\tTotal kernel (smem+sal+bsw) time avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
    
#if HIDE
    find_opt(tprof[MEM_ALN_M1], nthreads, &max, &min, &avg);
    fprintf(stderr, "\t\tMEM_ALN_CHAIN_FLT avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
    
    find_opt(tprof[MEM_ALN_M2], nthreads, &max, &min, &avg);
    fprintf(stderr, "\t\tMEM_ALN_CHAIN_SEED avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
#endif
    
    find_opt(tprof[MEM_COLLECT], nthreads, &max, &min, &avg);
    fprintf(stderr, "\t\tSMEM compute avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);

    find_opt(tprof[MEM_CHAIN], nthreads, &max, &min, &avg);
    fprintf(stderr, "\t\tMEM_CHAIN avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
    
    find_opt(tprof[MEM_SA_BLOCK], nthreads, &max, &min, &avg);
    fprintf(stderr, "\t\tSAL compute avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
    
    #if 1 //HIDE
    find_opt(tprof[MEM_SA], nthreads, &max, &min, &avg);
    fprintf(stderr, "\t\t\t\tMEM_SA avg: %0.2lf, (%0.2lf, %0.2lf)\n\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
    #endif
    
    // printf("\n\t BSW compute time (sec):\n");
    find_opt(tprof[MEM_ALN2], nthreads, &max, &min, &avg);
    fprintf(stderr, "\t\tBSW time, avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);

    // PR 26a/c.1/e/26c.2/Q1/Q2/Q3 diagnostics. Gated behind bwa_verbose >= 4
    // — these are histogram-heavy profile-build sections, useful for tuning
    // but noisy in production stderr. Match the convention used elsewhere
    // for high-verbosity diagnostics in bwa.cpp.
    if (bwa_verbose >= 4) {
    // PR 26a/c.1/e diagnostics.
    {
        uint64_t la=0, lh=0, lt=0, ltb1=0, ltb2=0, ltb3=0, ln=0;
        uint64_t ra=0, rh=0, rt=0, rtb1=0, rtb2=0, rtb3=0, rn=0;
        for (int t = 0; t < nthreads; t++) {
            la   += tprof[UGP_L_ATTEMPT][t];
            lh   += tprof[UGP_L_HIT][t];
            lt   += tprof[UGP_L_TIGHT][t];
            ltb1 += tprof[UGP_L_TB_1_8][t];
            ltb2 += tprof[UGP_L_TB_9_32][t];
            ltb3 += tprof[UGP_L_TB_33_MAX][t];
            ln   += tprof[UGP_L_DISP_NARROW][t];
            ra   += tprof[UGP_R_ATTEMPT][t];
            rh   += tprof[UGP_R_HIT][t];
            rt   += tprof[UGP_R_TIGHT][t];
            rtb1 += tprof[UGP_R_TB_1_8][t];
            rtb2 += tprof[UGP_R_TB_9_32][t];
            rtb3 += tprof[UGP_R_TB_33_MAX][t];
            rn   += tprof[UGP_R_DISP_NARROW][t];
        }
        if (la + ra > 0) {
            fprintf(stderr, "\n\tUngapped fast-path / tight_band diagnostics:\n");
            fprintf(stderr, "\t  LEFT   att=%lu  hit=%lu (%.2f%%)  tight=%lu  tb[1..8]=%lu  tb[9..32]=%lu  tb[33..]=%lu  disp_narrow=%lu\n",
                    (unsigned long)la, (unsigned long)lh,
                    la ? 100.0 * (double)lh / (double)la : 0.0,
                    (unsigned long)lt, (unsigned long)ltb1,
                    (unsigned long)ltb2, (unsigned long)ltb3,
                    (unsigned long)ln);
            fprintf(stderr, "\t  RIGHT  att=%lu  hit=%lu (%.2f%%)  tight=%lu  tb[1..8]=%lu  tb[9..32]=%lu  tb[33..]=%lu  disp_narrow=%lu\n",
                    (unsigned long)ra, (unsigned long)rh,
                    ra ? 100.0 * (double)rh / (double)ra : 0.0,
                    (unsigned long)rt, (unsigned long)rtb1,
                    (unsigned long)rtb2, (unsigned long)rtb3,
                    (unsigned long)rn);
        }
    }

    /* PR 26c.2 design instrumentation. */
    {
        const char *fine_lbl[UGP_FINE_NBINS] = {
            "0", "1-2", "3-4", "5-8", "9-16", "17-32", "33-48", "49-64", "65-80", "81-100"
        };
        const char *tier_lbl[2]   = { "8bit ", "16bit" };
        const char *side_lbl[2]   = { "LEFT ", "RIGHT" };
        const char *band_lbl[4]   = { "tb=0", "tb 1-8", "tb 9-32", "tb 33+" };
        const char *narrow_lbl[UGP_NARROW_SZ_NBINS] = {
            "0", "1-15", "16-31", "32-63", "64-127", "128-255", "256-511", "512+"
        };

        /* Group A: fine tight_band histogram per side. */
        for (int side = 0; side < 2; side++) {
            uint64_t total = 0;
            uint64_t bins[UGP_FINE_NBINS] = {0};
            for (int b = 0; b < UGP_FINE_NBINS; b++) {
                for (int t = 0; t < nthreads; t++)
                    bins[b] += tprof[UGP_FINE_BASE + side * UGP_FINE_NBINS + b][t];
                total += bins[b];
            }
            if (total == 0) continue;
            fprintf(stderr, "\t  %s tight_band fine: total=%lu  ",
                    side_lbl[side], (unsigned long)total);
            for (int b = 0; b < UGP_FINE_NBINS; b++)
                fprintf(stderr, "%s=%lu ", fine_lbl[b], (unsigned long)bins[b]);
            fprintf(stderr, "\n");
        }

        /* Group B: tier × band per side. The 2 tiers and 4 bands here mirror
         * the producer-side layout written into UGP_TIER_TB_BASE; keep these
         * named so any future bin/tier-count tweak stays in sync. */
        const int tier_n = 2;
        const int band_n = 4;
        for (int side = 0; side < 2; side++) {
            for (int tier = 0; tier < tier_n; tier++) {
                uint64_t bins[band_n] = {0};
                uint64_t total = 0;
                for (int b = 0; b < band_n; b++) {
                    for (int t = 0; t < nthreads; t++)
                        bins[b] += tprof[UGP_TIER_TB_BASE
                                         + side * (tier_n * band_n)
                                         + tier * band_n + b][t];
                    total += bins[b];
                }
                if (total == 0) continue;
                fprintf(stderr, "\t  %s %s tier:       total=%lu  ",
                        side_lbl[side], tier_lbl[tier], (unsigned long)total);
                for (int b = 0; b < band_n; b++)
                    fprintf(stderr, "%s=%lu ", band_lbl[b], (unsigned long)bins[b]);
                fprintf(stderr, "\n");
            }
        }

        /* Group C: per-batch narrow bucket size histogram per side. */
        for (int side = 0; side < 2; side++) {
            uint64_t bins[UGP_NARROW_SZ_NBINS] = {0};
            uint64_t total = 0;
            for (int b = 0; b < UGP_NARROW_SZ_NBINS; b++) {
                for (int t = 0; t < nthreads; t++)
                    bins[b] += tprof[UGP_NARROW_SZ_BASE + side * UGP_NARROW_SZ_NBINS + b][t];
                total += bins[b];
            }
            if (total == 0) continue;
            fprintf(stderr, "\t  %s narrow bucket sz/batch: total=%lu  ",
                    side_lbl[side], (unsigned long)total);
            for (int b = 0; b < UGP_NARROW_SZ_NBINS; b++)
                fprintf(stderr, "%s=%lu ", narrow_lbl[b], (unsigned long)bins[b]);
            fprintf(stderr, "\n");
        }
    }

    /* Q1+Q2 instrumentation: post-SW ungapped/gapped final outcome and
     * per-side alignment-score histogram. */
    {
        const char *side_lbl[2]  = { "LEFT ", "RIGHT" };
        const char *score_lbl[UGP_SCORE_HIST_NBINS] = {
            "0-10", "11-25", "26-50", "51-75", "76-100", "101-125", "126-150", "151+"
        };
        for (int side = 0; side < 2; side++) {
            uint64_t ung = 0, gap = 0;
            for (int t = 0; t < nthreads; t++) {
                ung += tprof[UGP_OUTCOME_BASE + side * 2 + 0][t];
                gap += tprof[UGP_OUTCOME_BASE + side * 2 + 1][t];
            }
            uint64_t total = ung + gap;
            if (total == 0) continue;
            fprintf(stderr, "\t  %s outcome: total=%lu  ungapped=%lu (%.2f%%)  gapped=%lu (%.2f%%)\n",
                    side_lbl[side], (unsigned long)total,
                    (unsigned long)ung, 100.0 * (double)ung / (double)total,
                    (unsigned long)gap, 100.0 * (double)gap / (double)total);
        }
        for (int side = 0; side < 2; side++) {
            uint64_t bins[UGP_SCORE_HIST_NBINS] = {0};
            uint64_t total = 0;
            for (int b = 0; b < UGP_SCORE_HIST_NBINS; b++) {
                for (int t = 0; t < nthreads; t++)
                    bins[b] += tprof[UGP_SCORE_HIST_BASE + side * UGP_SCORE_HIST_NBINS + b][t];
                total += bins[b];
            }
            if (total == 0) continue;
            fprintf(stderr, "\t  %s score hist:  total=%lu  ",
                    side_lbl[side], (unsigned long)total);
            for (int b = 0; b < UGP_SCORE_HIST_NBINS; b++)
                fprintf(stderr, "%s=%lu ", score_lbl[b], (unsigned long)bins[b]);
            fprintf(stderr, "\n");
        }
    }

    /* Q3 instrumentation: per-category delta-from-perfect histograms (LEFT).
     * delta = (h0 + a*len2) - aln_score; 0 = perfect ungapped extension. */
    {
        const char *cat_lbl[UGP_CAT_NCAT] = {
            "ALL                 ",
            "UNGAP_FINAL         ",
            "GAPPED_FINAL        ",
            "HIT                 ",
            "UNGAP_FINAL_NOT_HIT "
        };
        const char *cat_delta_lbl[UGP_CAT_NBINS] = {
            "0", "1-5", "6-10", "11-25", "26-50", "51-75", "76-100", "101+"
        };
        uint64_t any_ung = 0, any_fin = 0;
        for (int c = 0; c < UGP_CAT_NCAT; c++)
            for (int b = 0; b < UGP_CAT_NBINS; b++)
                for (int t = 0; t < nthreads; t++) {
                    any_ung += tprof[UGP_L_CAT_UNG_BASE + c * UGP_CAT_NBINS + b][t];
                    any_fin += tprof[UGP_L_CAT_FIN_BASE + c * UGP_CAT_NBINS + b][t];
                }
        if (any_ung > 0 || any_fin > 0) {
            fprintf(stderr, "\n\tQ3 LEFT per-category delta-from-perfect histograms:\n");
        }
        if (any_ung > 0) {
            fprintf(stderr, "\t  -- delta on would-be ungapped extension score --\n");
            for (int c = 0; c < UGP_CAT_NCAT; c++) {
                uint64_t bins[UGP_CAT_NBINS] = {0};
                uint64_t total = 0;
                for (int b = 0; b < UGP_CAT_NBINS; b++) {
                    for (int t = 0; t < nthreads; t++)
                        bins[b] += tprof[UGP_L_CAT_UNG_BASE + c * UGP_CAT_NBINS + b][t];
                    total += bins[b];
                }
                fprintf(stderr, "\t  cat%d %s ung total=%lu  ",
                        c, cat_lbl[c], (unsigned long)total);
                for (int b = 0; b < UGP_CAT_NBINS; b++)
                    fprintf(stderr, "%s=%lu ", cat_delta_lbl[b], (unsigned long)bins[b]);
                fprintf(stderr, "\n");
            }
        }
        if (any_fin > 0) {
            fprintf(stderr, "\t  -- delta on final committed alignment score --\n");
            for (int c = 0; c < UGP_CAT_NCAT; c++) {
                uint64_t bins[UGP_CAT_NBINS] = {0};
                uint64_t total = 0;
                for (int b = 0; b < UGP_CAT_NBINS; b++) {
                    for (int t = 0; t < nthreads; t++)
                        bins[b] += tprof[UGP_L_CAT_FIN_BASE + c * UGP_CAT_NBINS + b][t];
                    total += bins[b];
                }
                fprintf(stderr, "\t  cat%d %s fin total=%lu  ",
                        c, cat_lbl[c], (unsigned long)total);
                for (int b = 0; b < UGP_CAT_NBINS; b++)
                    fprintf(stderr, "%s=%lu ", cat_delta_lbl[b], (unsigned long)bins[b]);
                fprintf(stderr, "\n");
            }
        }
    }
    } /* end bwa_verbose >= 4 gate for diagnostic histograms */

    #if HIDE
    /* The seed-buffer re-alloc report that used to live here read tprof[PE11],
     * [PE12] and [PE13]. Those counters were removed from the chaining path (they
     * were never reported outside this disabled block and caused cross-thread
     * false sharing), so the report is gone with its writers rather than left to
     * print a permanent zero. */
    double res = 0.0, max_ = 0, min_=1e10;
    for (int i=0; i<nthreads; i++) {
        double val = (tprof[ALIGN1][i]*1.0) / tprof[MEM_CHAIN][i];
        res += val;
        if (max_ < val) max_ = val;
        if (min_ > val) min_ = val;
    }
    fprintf(stderr, "\tAvg. FM-index traversal per get_sa_entry(): avg: %lf, max: %lf, min: %lf\n",
            res/nthreads, max_, min_);

    int64_t tot_inst1 = 0, tot_inst2 = 0;
    for (int i=0; i<nthreads; i++) {
        tot_inst1 += tprof[SAM1][i];
        tot_inst2 += tprof[SAM2][i];
    }
    
    fprintf(stderr, "\ttot_inst1: %ld, tot_inst2: %ld, over %d threads\n",
            tot_inst1, tot_inst2, nthreads);
    #endif
    
#if HIDE
    fprintf(stderr, "\n BSW Perf.:\n");
    find_opt(tprof[MEM_ALN2_B], 1, &max, &min, &avg);
    fprintf(stderr, "\tLeft 16-bit time avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
    find_opt(tprof[MEM_ALN2_D], 1, &max, &min, &avg);
    fprintf(stderr, "\tLeft 8-bit time avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
    find_opt(tprof[MEM_ALN2_C], 1, &max, &min, &avg);
    fprintf(stderr, "\tRight 16-bit time avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
    find_opt(tprof[MEM_ALN2_E], 1, &max, &min, &avg);
    fprintf(stderr, "\tRight 8-bit time avg: %0.2lf, (%0.2lf, %0.2lf)\n",
            avg*1.0/proc_freq, max*1.0/proc_freq, min*1.0/proc_freq);
#endif

#if HIDE
    fprintf(stderr, "\nSTATSV\n");
    fprintf(stderr, "%0.2lf\n", tprof[PROCESS][0]*1.0/proc_freq);
    find_opt(tprof[MEM_PROCESS2], 1, &max, &min, &avg);
    fprintf(stderr, "%0.2lf\n", avg*1.0/proc_freq);
    find_opt(tprof[WORKER10], 1, &max, &min, &avg);
    fprintf(stderr, "%0.2lf\n", avg*1.0/proc_freq);
    find_opt(tprof[MEM_COLLECT], nthreads, &max, &min, &avg);
    fprintf(stderr, "%0.2lf\n", avg*1.0/proc_freq);
    find_opt(tprof[MEM_SA_BLOCK], nthreads, &max, &min, &avg);
    fprintf(stderr, "%0.2lf\n", avg*1.0/proc_freq);
    find_opt(tprof[MEM_SA], nthreads, &max, &min, &avg);
    fprintf(stderr, "%0.2lf\n", avg*1.0/proc_freq);
    double val = 0;
    find_opt(tprof[MEM_ALN2_UP], nthreads, &max, &min, &avg);
    fprintf(stderr, "%0.2lf\n", avg*1.0/proc_freq);
    val += avg;
    find_opt(tprof[CLEFT], nthreads, &max, &min, &avg);
    fprintf(stderr, "%0.2lf\n", avg*1.0/proc_freq);
    val += avg;
    find_opt(tprof[CRIGHT], nthreads, &max, &min, &avg);
    fprintf(stderr, "%0.2lf\n", avg*1.0/proc_freq);
    val += avg;
    find_opt(tprof[MEM_ALN2_DOWN], nthreads, &max, &min, &avg);
    fprintf(stderr, "%0.2lf\n", avg*1.0/proc_freq);
    val += avg;
    fprintf(stderr, "%0.2lf\n", val*1.0/proc_freq);
    fprintf(stderr, "%0.2lf\n", (tprof[REF_IO][0] + tprof[FMI][0])*1.0/proc_freq);
    fprintf(stderr, "%0.2lf\n", tprof[READ_IO][0]*1.0/proc_freq);
#endif
    // printf("\tMemory usage (GB):\n");
    // printf("\tAvg: %0.2lf, Peak: %0.2lf\n", tprof[PE21][0]*1.0/1e9, tprof[PE22][0]*1.0/1e9);

    return 1;
}

