// test/unit/test_mem_io_utils.cpp
//
// The src/utils.cpp helpers on the `mem` path: the klib sorts of 64-bit keys
// (ks_introsort_64, which orders seeds for chaining, and ks_combsort_64, its
// fallback) and the checked I/O wrappers (err_fwrite, err_fputs, err_fputc,
// err_fclose, err_fread_noeof, err_gzread, err_xopen_core, err_xzopen_core),
// whose contract is to behave like the stdio/zlib call they wrap and to end the
// process (EXIT_FAILURE) instead of returning a short result.
//
// Oracles: std::sort for the sorts (any correct sort of plain integers gives
// the one ascending order); the bytes written, read back with plain stdio and
// zlib, for the I/O wrappers; a forked child's exit status for the failures.

#include "doctest/doctest.h"
#include "mem_fixture.h"

#include "utils.h"

#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>
#include <zlib.h>

#include <algorithm>
#include <csignal>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <vector>

// Defined by KSORT_INIT in src/utils.cpp; not declared in utils.h.
void ks_combsort_64(size_t n, uint64_t *a);

namespace {

std::vector<std::vector<uint64_t> > sort_inputs() {
    std::mt19937_64 rng(97531u);
    std::vector<std::vector<uint64_t> > in;
    const size_t sizes[] = {0, 1, 2, 3, 15, 16, 17, 31, 64, 100, 257, 1000, 4096};
    for (size_t n : sizes) {
        std::vector<uint64_t> v(n);
        for (uint64_t &x : v) x = rng();
        in.push_back(v);                            // random 64-bit
        for (uint64_t &x : v) x = rng() % 7;        // many duplicates
        in.push_back(v);
        for (size_t i = 0; i < n; ++i) v[i] = i;    // sorted
        in.push_back(v);
        for (size_t i = 0; i < n; ++i) v[i] = n - i;  // reversed
        in.push_back(v);
        for (size_t i = 0; i < n; ++i) v[i] = (i % 2) ? UINT64_MAX - i : i;  // organ pipe + extremes
        in.push_back(v);
    }
    return in;
}

// Runs `body` in a forked child with stderr silenced and no core file; returns
// the wait status.
template <class F>
int in_child(F body) {
    fflush(stdout);
    fflush(stderr);
    const pid_t pid = fork();
    if (pid == 0) {
        signal(SIGABRT, SIG_DFL);
        struct rlimit no_core = {0, 0};
        setrlimit(RLIMIT_CORE, &no_core);
        if (freopen("/dev/null", "w", stderr) == NULL) _exit(3);
        body();
        _exit(0);
    }
    REQUIRE(pid > 0);
    int status = 0;
    REQUIRE(waitpid(pid, &status, 0) == pid);
    return status;
}

bool exited_failure(int status) { return WIFEXITED(status) && WEXITSTATUS(status) == EXIT_FAILURE; }

std::string slurp(const std::string &path) {
    FILE *fp = fopen(path.c_str(), "rb");
    REQUIRE(fp != NULL);
    std::string s;
    char buf[1024];
    size_t n;
    while ((n = fread(buf, 1, sizeof(buf), fp)) > 0) s.append(buf, n);
    fclose(fp);
    return s;
}

}  // namespace

TEST_CASE("mem_io utils: ks_introsort_64 and ks_combsort_64 sort ascending"
          * doctest::test_suite("unit/util")) {
    const std::vector<std::vector<uint64_t> > inputs = sort_inputs();
    for (size_t k = 0; k < inputs.size(); ++k) {
        CAPTURE(k);
        CAPTURE(inputs[k].size());
        std::vector<uint64_t> want = inputs[k];
        std::sort(want.begin(), want.end());
        std::vector<uint64_t> intro = inputs[k], comb = inputs[k];
        ks_introsort_64(intro.size(), intro.empty() ? NULL : intro.data());
        CHECK(intro == want);
        if (!comb.empty()) ks_combsort_64(comb.size(), comb.data());
        CHECK(comb == want);
    }
}

TEST_CASE("mem_io utils: the checked writers write exactly what they are given"
          * doctest::test_suite("unit/util")) {
    bwa_tests::TempDir dir;
    const std::string path = dir.file("out.txt");
    FILE *fp = err_xopen_core("test", path.c_str(), "w");
    REQUIRE(fp != NULL);
    const std::string block(70000, 'x');  // larger than a stdio buffer
    CHECK(err_fwrite("ab\0c", 1, 4, fp) == 4);
    CHECK(err_fputs("line\ttab\n", fp) >= 0);
    CHECK(err_fputc('Z', fp) == 'Z');
    CHECK(err_fwrite(block.data(), 7, 10000, fp) == 10000);
    CHECK(err_fclose(fp) == 0);
    CHECK(slurp(path) == std::string("ab\0c", 4) + "line\ttab\n" + "Z" + block);

    // err_fread_noeof reads back exactly; err_xopen_core("-") is stdin/stdout.
    FILE *in = err_xopen_core("test", path.c_str(), "r");
    REQUIRE(in != NULL);
    char buf[13];
    CHECK(err_fread_noeof(buf, 1, 13, in) == 13);
    CHECK(std::string(buf, 13) == std::string("ab\0c", 4) + "line\ttab\n");
    CHECK(err_fclose(in) == 0);
    CHECK(err_xopen_core("test", "-", "r") == stdin);
    CHECK(err_xopen_core("test", "-", "w") == stdout);
}

TEST_CASE("mem_io utils: err_gzread returns the decompressed bytes of plain and gzip input"
          * doctest::test_suite("unit/util")) {
    bwa_tests::TempDir dir;
    std::string text;
    for (int i = 0; i < 5000; ++i) text += "@r" + std::to_string(i) + "\nACGT\n+\nIIII\n";
    dir.write("plain.txt", text);
    gzFile gz = gzopen(dir.file("z.gz").c_str(), "wb");
    REQUIRE(gz != NULL);
    REQUIRE(gzwrite(gz, text.data(), (unsigned)text.size()) == (int)text.size());
    REQUIRE(gzclose(gz) == Z_OK);
    const char *files[] = {"plain.txt", "z.gz"};
    for (const char *f : files) {
        CAPTURE(f);
        gzFile g = err_xzopen_core("test", dir.file(f).c_str(), "r");
        REQUIRE(g != NULL);
        std::string got;
        char buf[777];
        int n;
        while ((n = err_gzread(g, buf, sizeof(buf))) > 0) got.append(buf, (size_t)n);
        CHECK(n == 0);
        CHECK(err_gzclose(g) == Z_OK);
        CHECK(got == text);
    }
}

TEST_CASE("mem_io utils: a failed open or write ends the process"
          * doctest::test_suite("unit/util")) {
    bwa_tests::TempDir dir;
    dir.write("ro.txt", "abc");
    const std::string ro = dir.file("ro.txt"), missing = dir.file("no/such/file");

    CHECK(exited_failure(in_child([&]() { err_xopen_core("test", missing.c_str(), "r"); })));
    CHECK(exited_failure(in_child([&]() { err_xzopen_core("test", missing.c_str(), "r"); })));
    // Writing to a stream opened read-only fails at the call.
    CHECK(exited_failure(in_child([&]() {
        FILE *fp = fopen(ro.c_str(), "r");
        if (fp == NULL) _exit(5);
        err_fwrite("x", 1, 1, fp);
    })));
    // Asking for more bytes than remain is an error for err_fread_noeof.
    CHECK(exited_failure(in_child([&]() {
        FILE *fp = fopen(ro.c_str(), "r");
        if (fp == NULL) _exit(5);
        char buf[8];
        err_fread_noeof(buf, 1, sizeof(buf), fp);
    })));
    // The controls: the same calls on good input return normally.
    const int ok = in_child([&]() {
        FILE *fp = err_xopen_core("test", ro.c_str(), "r");
        char buf[3];
        err_fread_noeof(buf, 1, 3, fp);
        err_fclose(fp);
    });
    CHECK(WIFEXITED(ok));
    CHECK(WEXITSTATUS(ok) == 0);
}
