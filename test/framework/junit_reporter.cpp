// test/framework/junit_reporter.cpp
#include "junit_reporter.h"

#include <cstdlib>
#include <cstdio>
#include <string>

namespace bwa_tests {

void install_matrix_row_reporter_hook() {
    // doctest's built-in junit reporter writes a default testsuite name
    // ("doctest"). We rename that via the environment variable doctest
    // reads at report-open time: DOCTEST_CONFIG_* isn't the right
    // mechanism here; instead, we set a wrapper-recognized env var that
    // the downstream JUnit aggregator (GitHub Actions) can key on.
    //
    // For Task #1 we do the minimum: print the matrix row as a stderr
    // banner before the test run so the CI log itself disambiguates.
    // A richer reporter override is scheduled for a follow-up PR once
    // we see what GitHub's test-results action actually consumes.
    const char *row = std::getenv("BWA_MEM3_TEST_MATRIX_ROW");
    if (row && row[0] != '\0') {
        std::fprintf(stderr, "\n=== bwa-mem3 tests — matrix row: %s ===\n\n", row);
    }
}

} // namespace bwa_tests
