// test/framework/run_in_child.h — run a test body in a forked child process.
//
// For death tests: code under test that ends the process (an xassert's abort(),
// err_fatal's exit(EXIT_FAILURE), a guard's _exit) cannot run in the test
// runner itself, so it runs in a fork and the parent inspects the wait status.
//
// The child touches no doctest machinery. It restores the default SIGABRT
// disposition before running the body: it inherits doctest's crash handler,
// which would otherwise catch the abort, report a spurious "test case
// CRASHED", and then resume the child's copy of the test runner, re-running the
// rest of the suite inside the fork. Core dumps are disabled, and stderr goes
// to `stderr_path` (default /dev/null) so a test can assert on the diagnostic.
//
// The parent's stdio buffers are flushed before the fork: a body that ends in
// exit() flushes the child's inherited copy of them, which would otherwise
// print the runner's pending output a second time.

#ifndef BWA_TESTS_RUN_IN_CHILD_H
#define BWA_TESTS_RUN_IN_CHILD_H

#include "doctest/doctest.h"

#include <csignal>
#include <cstdio>
#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>

namespace bwa_tests {

/* Run `body` in a forked child and return its wait status. A body that returns
 * ends the child with _exit(0), so any other status means it aborted or exited
 * on its own. */
template <typename Body>
int run_in_child(Body &&body, const char *stderr_path = "/dev/null") {
    fflush(NULL);
    const pid_t pid = fork();
    if (pid == 0) {
        signal(SIGABRT, SIG_DFL);
        struct rlimit no_core = {0, 0};
        setrlimit(RLIMIT_CORE, &no_core);
        (void)freopen(stderr_path, "w", stderr);
        /* freopen may leave stderr fully buffered (it does on macOS), and a
         * body that ends in _exit would then drop what it wrote. */
        setvbuf(stderr, NULL, _IONBF, 0);
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
