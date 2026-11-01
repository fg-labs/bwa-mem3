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

TEST_CASE("run_in_child: the body runs with every signal doctest handles at its default"
          * doctest::test_suite("unit/run_in_child")) {
    for (const int sig : bwa_tests::kDoctestCrashSignals) {
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

TEST_CASE("run_in_child: a body that returns exits 0, and one that exits keeps its code"
          * doctest::test_suite("unit/run_in_child")) {
    const int ok = bwa_tests::run_in_child([] {});
    REQUIRE(WIFEXITED(ok));
    CHECK(WEXITSTATUS(ok) == 0);
    const int coded = bwa_tests::run_in_child([] { exit(37); });
    REQUIRE(WIFEXITED(coded));
    CHECK(WEXITSTATUS(coded) == 37);
}
