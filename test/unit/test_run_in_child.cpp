// Unit tests for bwa_tests::run_in_child (test/framework/run_in_child.h), the
// fork helper the death tests use.
//
// The child inherits doctest's crash handler across fork. Left in place, a
// child that crashes runs it and writes a "test case CRASHED" report and a
// failing summary into the runner's output, as if the suite had failed. The
// helper must hand the body every signal the handler takes at its default
// disposition. That is checked directly (sigaction in the child, whatever
// reporter the run uses), and end to end for two signals that leave no core
// or crash report: the child dies by the signal, and its stdout, where the
// console reporter writes, stays empty.
//
// Fixtures are generated in-memory; no test data files are read.

#include "doctest/doctest.h"
#include "run_in_child.h"
#include "temp_dir.h"

#include <csignal>
#include <cstdlib>
#include <initializer_list>
#include <string>

#if defined(BWA_TESTS_CHILD_SKIP_LEAK_CHECK)
#include <thread>

namespace {
// A terminated thread leaves no live stack roots for LSan's conservative scan.
// The volatile access and compiler barrier keep the allocation in optimized builds.
__attribute__((noinline)) void leak_unreachable_memory() {
    void *p = malloc(12345);
    if (p == nullptr) _exit(90);
    *static_cast<volatile char *>(p) = 42;
    asm volatile("" : : "r"(p) : "memory");
}

void make_child_leak() {
    std::thread worker(leak_unreachable_memory);
    worker.join();
}
}  // namespace
#endif

TEST_CASE("run_in_child: a leaking exit preserves 37 with glibc ASAN and enabled LSan"
          * doctest::test_suite("unit/run_in_child")) {
#if defined(BWA_TESTS_CHILD_SKIP_LEAK_CHECK)
    bwa_tests::TempDir dir{"child_lsan"};
    const std::string canary_log = dir.file("canary.err");
    // No helper here: the unprotected canary must reach LSan's real exit check.
    fflush(nullptr);
    const pid_t canary = fork();
    REQUIRE(canary >= 0);
    if (canary == 0) {
        if (freopen(canary_log.c_str(), "w", stderr) == nullptr) _exit(90);
        make_child_leak();
        exit(37);
    }
    int canary_status = 0;
    REQUIRE(waitpid(canary, &canary_status, 0) == canary);
    const std::string diagnostic = dir.read("canary.err");
    INFO(diagnostic);
    // LSan's exit status follows ASAN_OPTIONS=exitcode, which may itself be 37,
    // so the report text, not the status, proves the leak check ran.
    REQUIRE(WIFEXITED(canary_status));
    REQUIRE(diagnostic.find("LeakSanitizer: detected memory leaks") != std::string::npos);
    REQUIRE(diagnostic.find("12345 byte(s)") != std::string::npos);

    const std::string child_log = dir.file("child.err");
    const int status = bwa_tests::run_in_child([] {
        make_child_leak();
        exit(37);  // _exit would bypass the on_exit/LSan contract under test.
    }, child_log.c_str());
    INFO(dir.read("child.err"));
    REQUIRE(WIFEXITED(status));
    CHECK(WEXITSTATUS(status) == 37);
    CHECK(dir.read("child.err").find("LeakSanitizer") == std::string::npos);
#else
    MESSAGE("SKIP: leaking exit contract requires Linux glibc ASAN with LSan");
#endif
}

TEST_CASE("run_in_child: the body runs with every signal doctest handles at its default"
          * doctest::test_suite("unit/run_in_child")) {
    for (const int sig : {SIGINT, SIGILL, SIGFPE, SIGSEGV, SIGTERM, SIGABRT}) {
        CAPTURE(sig);
        const int status = bwa_tests::run_in_child([sig] {
            struct sigaction cur;
            if (sigaction(sig, nullptr, &cur) != 0) _exit(90);
            _exit(cur.sa_handler == SIG_DFL ? 0 : 91);
        });
        CHECK((WIFEXITED(status) && WEXITSTATUS(status) == 0));
    }
}

TEST_CASE("run_in_child: a child ended by such a signal dies by it, without doctest's report"
          * doctest::test_suite("unit/run_in_child")) {
    bwa_tests::TempDir dir{"run_in_child"};
    for (const int sig : {SIGTERM, SIGINT}) {
        CAPTURE(sig);
        const std::string out = dir.file("child_stdout.txt");
        const int status = bwa_tests::run_in_child([&out, sig] {
            if (!freopen(out.c_str(), "w", stdout)) _exit(90);
            raise(sig);
            _exit(91);  // not reached: the default action ends the child
        });
        CHECK((WIFSIGNALED(status) && WTERMSIG(status) == sig));
        CHECK(dir.read("child_stdout.txt").empty());
    }
}

TEST_CASE("run_in_child: inherited blocked crash signals are unblocked only in the child"
          * doctest::test_suite("unit/run_in_child")) {
    sigset_t blocked, original;
    REQUIRE(sigemptyset(&blocked) == 0);
    for (const int sig : {SIGINT, SIGILL, SIGFPE, SIGSEGV, SIGTERM, SIGABRT})
        REQUIRE(sigaddset(&blocked, sig) == 0);
    REQUIRE(sigprocmask(SIG_BLOCK, &blocked, &original) == 0);
    struct RestoreMask {
        sigset_t mask;
        ~RestoreMask() { sigprocmask(SIG_SETMASK, &mask, nullptr); }
    } restore{original};
    const int status = bwa_tests::run_in_child([] {
        sigset_t current;
        if (sigprocmask(SIG_SETMASK, nullptr, &current) != 0) _exit(90);
        for (const int sig : {SIGINT, SIGILL, SIGFPE, SIGSEGV, SIGTERM, SIGABRT})
            if (sigismember(&current, sig) != 0) _exit(91);
        raise(SIGTERM);
        _exit(92);
    });
    CHECK((WIFSIGNALED(status) && WTERMSIG(status) == SIGTERM));
    sigset_t parent;
    REQUIRE(sigprocmask(SIG_SETMASK, nullptr, &parent) == 0);
    for (const int sig : {SIGINT, SIGILL, SIGFPE, SIGSEGV, SIGTERM, SIGABRT})
        CHECK(sigismember(&parent, sig) == 1);
}

TEST_CASE("run_in_child: a body that returns exits 0, and one that exits keeps its code"
          * doctest::test_suite("unit/run_in_child")) {
    const int ok = bwa_tests::run_in_child([] {});
    REQUIRE(WIFEXITED(ok));
    CHECK(WEXITSTATUS(ok) == 0);
    const int coded = bwa_tests::run_in_child([] { exit(37); });
    REQUIRE(WIFEXITED(coded));
    CHECK(WEXITSTATUS(coded) == 37);
}
