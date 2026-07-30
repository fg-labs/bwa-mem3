// Unit tests for the xassert() guard in utils.h.
//
// xassert is the project's release-safe fatal check. Two properties make it
// usable as the guard on every allocation site, and neither is obvious from a
// call site — a future edit to the macro can silently lose either one:
//
//   * It still fires under NDEBUG. That is the whole reason allocation checks
//     use xassert instead of assert(): a release build has to enforce them.
//   * It is a single statement, so a trailing `else` binds to the caller's
//     `if`, not to an `if` hidden inside the macro. The macro used to expand
//     to a bare `if`, which stole the caller's `else` branch.
//
// NDEBUG is defined below, ahead of every include, so this translation unit
// compiles in exactly the configuration the properties are about: the release
// one, where a plain assert() would vanish.

#ifndef NDEBUG
#define NDEBUG 1
#endif

#include "doctest/doctest.h"

#include "utils.h"

#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>

namespace {

/* Run `body` in a forked child with stderr and core dumps suppressed, and
 * report how it terminated. A child that runs to completion exits 0, so any
 * other status means the body aborted or exited on its own.
 *
 * The child touches no doctest machinery — it does not report assertions, and
 * it restores the default SIGABRT disposition before running `body`. It
 * inherits doctest's crash handler, which would otherwise catch the abort,
 * report a spurious "test case CRASHED", and then resume the child's copy of
 * the test runner, re-running the rest of the suite inside the fork. With
 * SIG_DFL the child just dies on the signal, which is what the parent waits
 * for. */
int run_in_child(void (*body)()) {
    const pid_t pid = fork();
    if (pid == 0) {
        signal(SIGABRT, SIG_DFL);
        struct rlimit no_core = {0, 0};
        setrlimit(RLIMIT_CORE, &no_core);
        (void)freopen("/dev/null", "w", stderr);
        body();
        _exit(0);
    }
    REQUIRE(pid > 0);
    int status = 0;
    REQUIRE(waitpid(pid, &status, 0) == pid);
    return status;
}

void fail_an_xassert() { xassert(1 == 2, "test: xassert must fire under NDEBUG"); }

void pass_an_xassert() { xassert(1 == 1, "test: xassert must not fire when true"); }

}  // namespace

TEST_CASE("xassert is not compiled out under NDEBUG") {
    const int status = run_in_child(fail_an_xassert);
    /* _err_fatal_simple_core ends in abort(), so the child dies on SIGABRT
     * rather than exiting. Either way it must not have reached _exit(0). */
    const bool exited_cleanly = WIFEXITED(status) && WEXITSTATUS(status) == 0;
    CHECK_FALSE(exited_cleanly);
    if (WIFSIGNALED(status)) CHECK(WTERMSIG(status) == SIGABRT);
}

TEST_CASE("xassert does nothing when its condition holds") {
    const int status = run_in_child(pass_an_xassert);
    /* REQUIRE, not CHECK: WEXITSTATUS is only defined when WIFEXITED holds, and
     * CHECK is non-fatal -- on a signal-terminated child it would record the
     * failure and then evaluate WEXITSTATUS anyway, which is undefined. The
     * fatal form makes the second line unreachable with an invalid status. */
    REQUIRE(WIFEXITED(status));
    CHECK(WEXITSTATUS(status) == 0);
}

TEST_CASE("a trailing else binds to the caller's if, not to xassert's") {
    /* With the old bare-`if` expansion the `else` below attached to the
     * macro's own `if`, so the whole statement was skipped and `branch`
     * stayed 0. As a single statement, the `else` is the caller's. */
    int branch = 0;
    if (false)
        xassert(1 == 1, "test: never evaluated");
    else
        branch = 1;
    CHECK(branch == 1);

    /* The mirror case: the caller's `if` is taken, so the `else` must not run
     * even though xassert's condition holds and its inner `if` is false. */
    branch = 0;
    if (true)
        xassert(1 == 1, "test: never fires");
    else
        branch = 2;
    CHECK(branch == 0);
}
