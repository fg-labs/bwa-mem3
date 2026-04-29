// test/framework/test_main.cpp
//
// Shared entry point for both bwa_mem3_tests_unit and
// bwa_mem3_tests_integration. This file is the ONE translation unit that
// defines DOCTEST_CONFIG_IMPLEMENT; it is linked into each test binary.

#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest/doctest.h"

#include "junit_reporter.h"

int main(int argc, char **argv) {
    bwa_tests::install_matrix_row_reporter_hook();

    doctest::Context ctx;
    ctx.applyCommandLine(argc, argv);
    return ctx.run();
}
