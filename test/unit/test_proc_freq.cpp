// Unit tests for proc_freq_resolve() (src/profiling.cpp), the __rdtsc() tick
// rate every timing report divides by.
//
// The CLI test (test/proc_freq_calibration_test.sh) checks only the source a
// real run used (at -v 4), so each branch is pinned here directly:
//
//   * a nonzero proc_freq set by the caller is kept, and named as the source;
//   * aarch64: the rate is CNTFRQ_EL0, read without a window or a sleep;
//   * aarch64: a CNTFRQ_EL0 that reads 0 falls through to the measured paths;
//   * a normal window is used WITHOUT the sleep calibration -- the point of
//     measuring over the run is that no process sleeps to get the rate;
//   * the sleep calibration is the fallback when the window was never opened,
//     is under 1 ms, could not be opened or closed (clock_gettime() failed),
//     spans a system suspend, or saw a counter that did not advance (c1 <= c0);
//   * the sleep calibration retries an attempt whose opening or closing clock
//     read failed rather than dividing by a bogus elapsed time.
//
// Not reached here: calibrate_proc_freq()'s "could not calibrate" warning
// after four failed attempts, which would take faking nanosleep() as well,
// and which only prints.
//
// Every rate is also checked against one this test measures itself, with
// std::chrono::steady_clock around a sleep, so a scaling bug in the resolve
// math (a dropped divide, wrong units) fails on every architecture rather
// than only where it happens to leave a 1..100000 MHz band.
//
// The hooks below are defined in profiling.cpp for this test only. The reset
// hook's argument says whether aarch64 reads CNTFRQ_EL0 (1), skips it (0), or
// treats it as reading 0 (-1), so the measured branches are reachable there
// too; the c0 hook puts the window's opening count ahead of the counter; the
// clock-read hook makes one clock_gettime() fail; the suspend hook fakes a
// system sleep inside the window; the source hook reports what the last
// resolve named.

#include "doctest/doctest.h"

#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <ctime>
#include <string>

#include "utils.h"       /* __rdtsc shim for GCC < 11 on x86 */
#include "profiling.h"
#if defined(__ARM_NEON) || defined(__aarch64__)
#include "simd_compat.h" /* __rdtsc: sse2neon's _rdtsc (CNTVCT_EL0) */
#else
#include <immintrin.h>
#endif

extern void bwamem3_proc_freq_reset_for_test(int use_counter_register);
extern void bwamem3_proc_freq_set_c0_for_test(uint64_t c0);
extern int bwamem3_proc_freq_sleep_calibrations(void);
extern const char *bwamem3_proc_freq_source_for_test(void);
extern void bwamem3_proc_freq_fail_clock_read_for_test(int n);
extern int bwamem3_proc_freq_add_suspend_for_test(double s);

namespace {

void sleep_ms(long ms)
{
    struct timespec req = {ms / 1000, (ms % 1000) * 1000L * 1000L};
    while (nanosleep(&req, &req) != 0) {
        // Only a signal is retried; any other error would loop forever.
        if (errno != EINTR) FAIL("nanosleep failed: " << std::strerror(errno));
    }
}

std::string source() {
    const char *s = bwamem3_proc_freq_source_for_test();
    return s ? s : "<none>";
}

// Every rate here comes from a constant-rate counter, so the only error is a
// clock/counter read pair split by preemption. profiling.cpp brackets each
// clock read with two counter reads and keeps the tightest of a few tries, so
// a split pair is rare there; the reference below reads clock outside counter
// (clock, counter, ..., counter, clock), so a split stretches the clock's span
// past the counter's and biases it LOW: a thread preempted for 1 ms inside a
// 10 ms window reads 10% low. A loaded runner does that now and then, so each
// rate below is measured up to kAttempts times and the reference keeps the
// highest of its attempts. A scaling bug in the resolve math (a dropped
// divide, wrong units) is deterministic and fails every attempt.
constexpr int kAttempts = 3;

// The tick rate measured independently of profiling.cpp: counter ticks over
// a 10 ms sleep, divided by steady_clock's elapsed time for the same span,
// highest of kAttempts. Measured once per process and shared by the cases, so
// each stays inside the unit-test time budget.
double reference_rate()
{
    static const double rate = [] {
        double best = 0.0;
        for (int attempt = 0; attempt < kAttempts; ++attempt) {
            const auto t0 = std::chrono::steady_clock::now();
            const uint64_t c0 = __rdtsc();
            sleep_ms(10);
            const uint64_t c1 = __rdtsc();
            const auto t1 = std::chrono::steady_clock::now();
            const double r = (double)(c1 - c0) / std::chrono::duration<double>(t1 - t0).count();
            if (r > best) best = r;
        }
        return best;
    }();
    return rate;
}

// Within 5% of the reference: room for clock read jitter, and still far
// tighter than any scaling bug, which is off by orders of magnitude.
bool near(double rate, double ref) { return std::fabs(rate / ref - 1.0) < 0.05; }

// Leave the globals as a fresh process has them for whatever runs next.
struct ResetAfter {
    ~ResetAfter() { bwamem3_proc_freq_reset_for_test(1); }
};

} // namespace

TEST_CASE("proc_freq: a nonzero rate set by the caller is kept"
          * doctest::test_suite("unit/proc_freq")) {
    ResetAfter guard;
    bwamem3_proc_freq_reset_for_test(1);
    proc_freq = 12345;
    CHECK(proc_freq_resolve() == 12345);
    CHECK(proc_freq == 12345);
    CHECK(bwamem3_proc_freq_sleep_calibrations() == 0);
    // Named as the source like every other one, so -v 4 reports it.
    CHECK(source() == "the caller");
}

TEST_CASE("proc_freq: opening the window resolves nothing and does not sleep"
          * doctest::test_suite("unit/proc_freq")) {
    ResetAfter guard;
    bwamem3_proc_freq_reset_for_test(0);
    proc_freq_mark_start();
    CHECK(proc_freq == 0);
    CHECK(bwamem3_proc_freq_sleep_calibrations() == 0);
}

TEST_CASE("proc_freq: a window that was never opened falls back to the sleep"
          * doctest::test_suite("unit/proc_freq")) {
    ResetAfter guard;
    const double ref = reference_rate();
    uint64_t f = 0;
    for (int attempt = 0; attempt < kAttempts; ++attempt) {
        bwamem3_proc_freq_reset_for_test(0);
        f = proc_freq_resolve();
        if (near((double)f, ref)) break;
    }
    CHECK(bwamem3_proc_freq_sleep_calibrations() == 1);
    CHECK(f == proc_freq);
    CHECK_MESSAGE(near((double)f, ref), "fallback ", f, " vs reference ", ref);
}

TEST_CASE("proc_freq: a window under 1 ms falls back to the sleep"
          * doctest::test_suite("unit/proc_freq")) {
    ResetAfter guard;
    const double ref = reference_rate();
    // Back-to-back calls leave a window of microseconds. Preemption could
    // stretch one attempt past 1 ms, or split a read pair in the 5 ms sleep
    // it falls back to, so allow a few attempts; it would take every one of
    // them being preempted to fail.
    bool fell_back = false;
    uint64_t f = 0;
    for (int attempt = 0; attempt < 5 && !(fell_back && near((double)f, ref)); ++attempt) {
        bwamem3_proc_freq_reset_for_test(0);
        proc_freq_mark_start();
        f = proc_freq_resolve();
        fell_back = bwamem3_proc_freq_sleep_calibrations() == 1;
    }
    CHECK(fell_back);
    CHECK_MESSAGE(near((double)f, ref), "short-window fallback ", f, " vs reference ", ref);
}

TEST_CASE("proc_freq: a normal window is used without the sleep calibration"
          * doctest::test_suite("unit/proc_freq")) {
    ResetAfter guard;
    const double ref = reference_rate();

    uint64_t fallback = 0, f = 0;
    int window_sleeps = -1;
    for (int attempt = 0; attempt < kAttempts; ++attempt) {
        bwamem3_proc_freq_reset_for_test(0);
        fallback = proc_freq_resolve();
        REQUIRE(bwamem3_proc_freq_sleep_calibrations() == 1);

        bwamem3_proc_freq_reset_for_test(0);
        proc_freq_mark_start();
        sleep_ms(10);
        f = proc_freq_resolve();
        window_sleeps = bwamem3_proc_freq_sleep_calibrations();
        if (window_sleeps == 0 && near((double)f, ref) && near((double)f, (double)fallback)) break;
    }
    CHECK(window_sleeps == 0);
    CHECK(source() == "the run window");
    CHECK_MESSAGE(near((double)f, ref), "run window ", f, " vs reference ", ref);
    CHECK_MESSAGE(near((double)f, (double)fallback), "run window ", f, " vs fallback ", fallback);

    // Resolved once; a later report reuses the rate.
    CHECK(proc_freq_resolve() == f);
    CHECK(bwamem3_proc_freq_sleep_calibrations() == 0);
}

TEST_CASE("proc_freq: a counter that did not advance falls back to the sleep"
          * doctest::test_suite("unit/proc_freq")) {
    ResetAfter guard;
    const double ref = reference_rate();
    // A window long enough to be used, but with its opening count put ahead
    // of the counter (as a counter reset across a suspend leaves it): the
    // c1 <= c0 arm must fall back rather than divide a wrapped difference.
    bool fell_back = false;
    uint64_t f = 0;
    for (int attempt = 0; attempt < kAttempts && !(fell_back && near((double)f, ref)); ++attempt) {
        bwamem3_proc_freq_reset_for_test(0);
        proc_freq_mark_start();
        sleep_ms(2);
        bwamem3_proc_freq_set_c0_for_test(UINT64_MAX);
        f = proc_freq_resolve();
        fell_back = bwamem3_proc_freq_sleep_calibrations() == 1;
    }
    CHECK(fell_back);
    CHECK_MESSAGE(near((double)f, ref), "counter-reset fallback ", f, " vs reference ", ref);
}

// Each (clock, counter) sample makes four clock reads, so read 1 opens
// a window (or the sleep calibration's first attempt) and read 5 closes it.
// The opening read failing while the closing one succeeds would, unchecked,
// divide by the clock's whole uptime and resolve a tiny rate.
TEST_CASE("proc_freq: the sleep calibration retries a failed clock read"
          * doctest::test_suite("unit/proc_freq")) {
    ResetAfter guard;
    const double ref = reference_rate();
    for (int fail : {1, 5}) {
        CAPTURE(fail);
        uint64_t f = 0;
        for (int attempt = 0; attempt < kAttempts; ++attempt) {
            bwamem3_proc_freq_reset_for_test(0);
            bwamem3_proc_freq_fail_clock_read_for_test(fail);
            f = proc_freq_resolve();
            if (near((double)f, ref)) break;
        }
        CHECK(bwamem3_proc_freq_sleep_calibrations() == 1);
        CHECK(source() == "a sleep calibration");
        CHECK_MESSAGE(near((double)f, ref), "calibration after a failed read ", f, " vs reference ", ref);
    }
}

TEST_CASE("proc_freq: a window whose clock read failed falls back to the sleep"
          * doctest::test_suite("unit/proc_freq")) {
    ResetAfter guard;
    const double ref = reference_rate();
    // 1: proc_freq_mark_start() cannot open the window; 5: proc_freq_resolve()
    // cannot close it.
    for (int fail : {1, 5}) {
        CAPTURE(fail);
        bool fell_back = false;
        uint64_t f = 0;
        for (int attempt = 0; attempt < kAttempts && !(fell_back && near((double)f, ref)); ++attempt) {
            bwamem3_proc_freq_reset_for_test(0);
            bwamem3_proc_freq_fail_clock_read_for_test(fail);
            proc_freq_mark_start();
            sleep_ms(2);
            f = proc_freq_resolve();
            fell_back = bwamem3_proc_freq_sleep_calibrations() == 1;
        }
        CHECK(fell_back);
        CHECK_MESSAGE(near((double)f, ref), "failed-read fallback ", f, " vs reference ", ref);
    }
}

TEST_CASE("proc_freq: a window that spans a suspend falls back to the sleep"
          * doctest::test_suite("unit/proc_freq")) {
    ResetAfter guard;
    bwamem3_proc_freq_reset_for_test(0);
    proc_freq_mark_start();
    if (!bwamem3_proc_freq_add_suspend_for_test(60.0)) {
        MESSAGE("no suspend clocks on this platform; the window is used unchecked");
        return;
    }
    const double ref = reference_rate();
    bool fell_back = false;
    uint64_t f = 0;
    for (int attempt = 0; attempt < kAttempts && !(fell_back && near((double)f, ref)); ++attempt) {
        bwamem3_proc_freq_reset_for_test(0);
        proc_freq_mark_start();
        sleep_ms(2);
        bwamem3_proc_freq_add_suspend_for_test(60.0);
        f = proc_freq_resolve();
        fell_back = bwamem3_proc_freq_sleep_calibrations() == 1;
    }
    CHECK(fell_back);
    CHECK_MESSAGE(near((double)f, ref), "suspend fallback ", f, " vs reference ", ref);
}

#if defined(__aarch64__)
// Measures over a 10 ms window with CNTFRQ_EL0 read as 0 (use_counter_register
// as given), and checks the window was used without a sleep, near ref.
static void check_cntfrq_zero_falls_through(int use_counter_register, double ref)
{
    uint64_t f = 0;
    int sleeps = -1;
    for (int attempt = 0; attempt < kAttempts; ++attempt) {
        bwamem3_proc_freq_reset_for_test(use_counter_register);
        proc_freq_mark_start();
        sleep_ms(10);
        f = proc_freq_resolve();
        sleeps = bwamem3_proc_freq_sleep_calibrations();
        if (sleeps == 0 && near((double)f, ref)) break;
    }
    CHECK(sleeps == 0);
    CHECK_MESSAGE(near((double)f, ref), "CNTFRQ_EL0 == 0 fall-through ", f, " vs reference ", ref);
}

TEST_CASE("proc_freq: aarch64 falls through to the run window when CNTFRQ_EL0 reads 0"
          * doctest::test_suite("unit/proc_freq")) {
    ResetAfter guard;
    // -1: the register read is replaced by 0, as on firmware that left it so.
    check_cntfrq_zero_falls_through(-1, reference_rate());
}

TEST_CASE("proc_freq: aarch64 reads CNTFRQ_EL0 without a window or a sleep"
          * doctest::test_suite("unit/proc_freq")) {
    ResetAfter guard;
    uint64_t cntfrq;
    __asm__ __volatile__("mrs %0, cntfrq_el0" : "=r"(cntfrq));
    const double ref = reference_rate();
    if (cntfrq == 0) {
        // Firmware left the register at 0, which production handles by
        // measuring instead: check that on the real register read.
        MESSAGE("CNTFRQ_EL0 reads 0 on this host; checking the measured fall-through");
        check_cntfrq_zero_falls_through(1, ref);
        return;
    }

    // Never opened: still exact, no sleep.
    bwamem3_proc_freq_reset_for_test(1);
    CHECK(proc_freq_resolve() == cntfrq);
    CHECK(bwamem3_proc_freq_sleep_calibrations() == 0);

    // Opened and closed at once (under 1 ms): still exact, no sleep.
    bwamem3_proc_freq_reset_for_test(1);
    proc_freq_mark_start();
    CHECK(proc_freq_resolve() == cntfrq);
    CHECK(bwamem3_proc_freq_sleep_calibrations() == 0);

    // The register's rate is the rate __rdtsc() actually ticks at.
    CHECK_MESSAGE(near((double)cntfrq, ref), "CNTFRQ_EL0 ", cntfrq, " vs reference ", ref);
}
#endif
