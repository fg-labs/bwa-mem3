// test/framework/junit_reporter.h
//
// Install a doctest reporter hook that decorates the built-in junit
// reporter output with the CI matrix row name (read from the
// BWA_MEM3_TEST_MATRIX_ROW env var; empty if unset). This keeps the
// per-row result artifacts disambiguable after actions/upload-artifact
// collects them.

#ifndef BWA_TESTS_JUNIT_REPORTER_H
#define BWA_TESTS_JUNIT_REPORTER_H

namespace bwa_tests {

// Must be called before doctest::Context::run(). Idempotent.
void install_matrix_row_reporter_hook();

} // namespace bwa_tests

#endif
