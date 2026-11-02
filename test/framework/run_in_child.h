// test/framework/run_in_child.h — run a test body in a forked child process.
//
// For death tests: code under test that ends the process (an xassert's abort(),
// err_fatal's exit(EXIT_FAILURE), a guard's _exit) cannot run in the test
// runner itself, so it runs in a fork and the parent inspects the wait status.
//
// The child touches no doctest machinery. It restores the default disposition
// of every signal doctest's crash handler takes (kDoctestCrashSignals) before
// unblocking those signals and running the body. It inherits that handler,
// which would otherwise catch the signal, write a spurious "test case CRASHED"
// report and a failing summary to
// the runner's output (stdout, or the --out file of a JUnit run), and only then
// re-raise it. With the defaults, a body that aborts or crashes ends the child
// by that signal, which the parent sees in the wait status.
//
// Core dumps are disabled, and stderr goes to `stderr_path` (default
// /dev/null) so a test can assert on the diagnostic.
//
// The parent's stdio buffers are flushed before the fork: a body that ends in
// exit() flushes the child's inherited copy of them, which would otherwise
// print the runner's pending output a second time.
//
// Under AddressSanitizer (glibc) the child skips LeakSanitizer's check at its
// exit: a body that ends the process on purpose (often with reader threads
// still holding buffers) leaves its allocations behind by design, and a leak
// report would change the exit status the test checks. An on_exit handler,
// which runs before the sanitizer's own exit handler, flushes stdio and ends
// the child with the status it was exiting with. Memory errors are still
// reported where they happen, and they change the status (a CI run's
// ASAN_OPTIONS exitcode makes that status distinct). The parent's own exit
// check still covers every allocation the test runner makes.

#ifndef BWA_TESTS_RUN_IN_CHILD_H
#define BWA_TESTS_RUN_IN_CHILD_H

#include "doctest/doctest.h"

#include <csignal>
#include <cstdio>
#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>

#if defined(__SANITIZE_ADDRESS__)
#define BWA_TESTS_ASAN 1
#elif defined(__has_feature)
#if __has_feature(address_sanitizer)
#define BWA_TESTS_ASAN 1
#endif
#endif
#if defined(BWA_TESTS_ASAN) && defined(__GLIBC__)
#include <cstdlib>
#define BWA_TESTS_CHILD_SKIP_LEAK_CHECK 1
#endif

namespace bwa_tests {

/* The signals doctest's POSIX crash handler takes: `signalDefs` in
 * ext/doctest/doctest.h (2.4.11). Keep in step on a doctest upgrade. */
constexpr int kDoctestCrashSignals[] = {SIGINT, SIGILL, SIGFPE, SIGSEGV, SIGTERM, SIGABRT};

#if defined(BWA_TESTS_CHILD_SKIP_LEAK_CHECK)
/* Ends a death-test child before LeakSanitizer's exit check (see above). */
inline void child_exit_without_leak_check(int status, void *) {
    fflush(NULL);
    _exit(status);
}
#endif

/* Run `body` in a forked child and return its wait status. A body that returns
 * ends the child with _exit(0), so any other status means it aborted or exited
 * on its own. */
template <typename Body>
int run_in_child(Body &&body, const char *stderr_path = "/dev/null") {
    fflush(NULL);
    const pid_t pid = fork();
    if (pid == 0) {
#if defined(BWA_TESTS_CHILD_SKIP_LEAK_CHECK)
        on_exit(child_exit_without_leak_check, nullptr);
#endif
        for (const int sig : kDoctestCrashSignals) signal(sig, SIG_DFL);
        struct rlimit no_core = {0, 0};
        setrlimit(RLIMIT_CORE, &no_core);
        (void)freopen(stderr_path, "w", stderr);
        /* freopen may leave stderr fully buffered (it does on macOS), and a
         * body that ends in _exit would then drop what it wrote. */
        setvbuf(stderr, NULL, _IONBF, 0);
        // Reset dispositions and disable core dumps before delivering pending signals.
        sigset_t unblock;
        sigemptyset(&unblock);
        for (const int sig : kDoctestCrashSignals) sigaddset(&unblock, sig);
        sigprocmask(SIG_UNBLOCK, &unblock, nullptr);
        body();
        _exit(0);
    }
    REQUIRE(pid > 0);
    int status = 0;
    REQUIRE(waitpid(pid, &status, 0) == pid);
    return status;
}

}  // namespace bwa_tests

#endif  // BWA_TESTS_RUN_IN_CHILD_H
