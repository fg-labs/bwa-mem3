// test/unit/test_index_load_threads.cpp
//
// Regression test for the fail-closed BWA3_LOAD_THREADS parse in
// index_load_threads_capped (src/FMI_search.cpp), and for the per-CPU default
// cap (4 on AMD Zen 3, 8 elsewhere). The parse cases call the capped form with
// an explicit cap so they hold on any host, Zen 3 CI runners included.
//
// The override used to go through atoi, which silently accepts garbage
// ("abc" -> 0, "0x8" -> 0) and any magnitude ("100000" -> ~1280 pthread_creates
// on a bandwidth-bound load). The hardened parser strtol's the value and, on
// anything non-numeric / with trailing garbage / non-positive, warns and keeps
// the computed default instead of running with a bad worker count. A valid
// override is clamped to a 64-thread ceiling. These are the "malformed
// BWA3_LOAD_THREADS" fail-closed cases; they need no index on disk, so they are
// a pure unit test of the parse contract.

#include "doctest/doctest.h"
#include "../../src/FMI_search.h"
#include "../../src/simd_dispatch.h"

#include <cstdlib>
#include <string>

namespace {

// RAII guard so each case sets BWA3_LOAD_THREADS in isolation and the process
// environment is restored to its prior state on scope exit, even if a CHECK
// throws. It captures whether the variable was originally set (and its value)
// before mutating it, so a BWA3_LOAD_THREADS already present in the test
// runner's environment survives the test instead of being unconditionally
// cleared -- keeping later cases (and any other test in the binary) isolated.
struct LoadThreadsEnv {
    explicit LoadThreadsEnv(const char *value)
    {
        if (const char *prev = getenv("BWA3_LOAD_THREADS")) {
            had_prev_ = true;
            prev_     = prev;  // copy: the setenv/unsetenv below may invalidate `prev`
        }
        if (value != nullptr)
            setenv("BWA3_LOAD_THREADS", value, /*overwrite=*/1);
        else
            unsetenv("BWA3_LOAD_THREADS");
    }
    ~LoadThreadsEnv()
    {
        if (had_prev_)
            setenv("BWA3_LOAD_THREADS", prev_.c_str(), /*overwrite=*/1);
        else
            unsetenv("BWA3_LOAD_THREADS");
    }

  private:
    bool had_prev_ = false;
    std::string prev_;
};

}  // namespace

TEST_CASE("index_load_threads: no override clamps the caller's request to [1,8]"
          * doctest::test_suite("unit/index_load_threads")) {
    LoadThreadsEnv env(nullptr);  // ensure a stray value from the shell can't leak in

    // A non-positive request clamps up to 1, never to 0 (callers spawn `t`
    // pread workers and would otherwise create none / divide by zero).
    CHECK(index_load_threads_capped(0, 8) == 1);
    CHECK(index_load_threads_capped(-4, 8) == 1);

    // In-range requests pass through; over the cap saturates at 8.
    CHECK(index_load_threads_capped(1, 8) == 1);
    CHECK(index_load_threads_capped(4, 8) == 4);
    CHECK(index_load_threads_capped(8, 8) == 8);
    CHECK(index_load_threads_capped(9, 8) == 8);
    CHECK(index_load_threads_capped(64, 8) == 8);

    // The Zen 3 cap of 4 saturates the same way.
    CHECK(index_load_threads_capped(0, 4) == 1);
    CHECK(index_load_threads_capped(3, 4) == 3);
    CHECK(index_load_threads_capped(4, 4) == 4);
    CHECK(index_load_threads_capped(16, 4) == 4);

    // A cap below 1 still yields one worker.
    CHECK(index_load_threads_capped(16, 0) == 1);
}

TEST_CASE("index_load_threads: the default cap is 4 on AMD Zen 3 and 8 elsewhere"
          * doctest::test_suite("unit/index_load_threads")) {
    CHECK(index_load_threads_default_cap(1, 0x19u, 0x01u) == 4);   // Milan (EPYC 7xx3)
    CHECK(index_load_threads_default_cap(1, 0x19u, 0x21u) == 4);   // Vermeer
    CHECK(index_load_threads_default_cap(1, 0x19u, 0x11u) == 8);   // Genoa (Zen 4)
    CHECK(index_load_threads_default_cap(1, 0x19u, 0xa0u) == 8);   // Bergamo (Zen 4c)
    CHECK(index_load_threads_default_cap(1, 0x17u, 0x31u) == 8);   // Zen 2
    CHECK(index_load_threads_default_cap(1, 0x1au, 0x02u) == 8);   // Zen 5
    CHECK(index_load_threads_default_cap(0, 0x6u, 0x8fu) == 8);    // Intel
    CHECK(index_load_threads_default_cap(0, 0u, 0u) == 8);         // arm64 / no CPUID

    // Independent expectations at every known range boundary and its neighbors.
    // In particular, 80h-9Fh are Zen 4, not the complement of a Zen 4 denylist;
    // unrecognized models must retain the generic eight-worker default.
    struct ModelCase { uint32_t model; int zen3; int cap; };
    for (const ModelCase c : {ModelCase{0x00u, 1, 4}, {0x01u, 1, 4},
                             {0x0fu, 1, 4}, {0x10u, 0, 8}, {0x1fu, 0, 8},
                             {0x20u, 1, 4}, {0x21u, 1, 4}, {0x5fu, 1, 4},
                             {0x60u, 0, 8}, {0x61u, 0, 8}, {0x7fu, 0, 8},
                             {0x80u, 0, 8}, {0x90u, 0, 8}, {0x9fu, 0, 8},
                             {0xa0u, 0, 8}, {0xafu, 0, 8}, {0xb0u, 0, 8},
                             {0xffu, 0, 8}, {0x100u, 0, 8}}) {
        CAPTURE(c.model);
        CHECK(bwa3_cpu_is_zen3(1, 0x19u, c.model) == c.zen3);
        CHECK(index_load_threads_default_cap(1, 0x19u, c.model) == c.cap);
        CHECK(bwa3_cpu_is_zen3(0, 0x19u, c.model) == 0);
        CHECK(index_load_threads_default_cap(0, 0x19u, c.model) == 8);
    }
    CHECK(bwa3_cpu_is_zen3(1, 0x18u, 0x01u) == 0);
    CHECK(bwa3_cpu_is_zen3(1, 0x1au, 0x01u) == 0);

    // index_load_threads is the capped form at this host's cap.
    LoadThreadsEnv env(nullptr);
    int is_amd = 0;
    uint32_t family = 0, model = 0;
    bwa3_host_cpu_identity(&is_amd, &family, &model);
    const int cap = index_load_threads_default_cap(is_amd, family, model);
    CHECK(index_load_threads_host_cap() == cap);
    for (int n : {0, 1, 4, 8, 16, 64}) {
        CAPTURE(n);
        CHECK(index_load_threads(n) == index_load_threads_capped(n, cap));
    }
}

TEST_CASE("index_load_threads: a malformed override is ignored, keeping the default"
          * doctest::test_suite("unit/index_load_threads")) {
    // Every malformed form must fall back to the computed default (here 4),
    // proving the parse fails closed rather than running with a bad count. Each
    // form is its own SUBCASE so a failure names the offending input directly
    // instead of hiding behind a shared loop iteration.
    const char *value = nullptr;
    SUBCASE("abc (non-numeric)")                      { value = "abc"; }
    SUBCASE("empty (treated as unset: default stands)") { value = ""; }
    SUBCASE("12x (trailing garbage)")                { value = "12x"; }
    SUBCASE("0x8 (hex-looking; base-10 strtol stops at 'x')") { value = "0x8"; }
    SUBCASE("whitespace only (no digits consumed)")  { value = "  "; }
    SUBCASE("0 (not positive)")                      { value = "0"; }
    SUBCASE("-1 (not positive)")                     { value = "-1"; }
    SUBCASE("3.5 (trailing '.5')")                   { value = "3.5"; }
    SUBCASE("overflows long (errno == ERANGE)")      { value = "9999999999999999999999"; }

    LoadThreadsEnv env(value);
    CHECK(index_load_threads_capped(4, 8) == 4);
    CHECK(index_load_threads_capped(16, 4) == 4);
}

TEST_CASE("index_load_threads: a valid override wins and is clamped to the 64 ceiling"
          * doctest::test_suite("unit/index_load_threads")) {
    {
        LoadThreadsEnv env("16");  // above the [1,8] default cap, below the ceiling
        CHECK(index_load_threads_capped(4, 8) == 16);
        CHECK(index_load_threads_capped(16, 4) == 16);  // overrides the Zen 3 cap too
        CHECK(index_load_threads(4) == 16);              // and this host's cap
    }
    {
        LoadThreadsEnv env("1");
        CHECK(index_load_threads_capped(8, 8) == 1);  // override can also lower the count
    }
    {
        LoadThreadsEnv env("64");
        CHECK(index_load_threads_capped(4, 8) == 64);
    }
    {
        LoadThreadsEnv env("100000");  // the runaway atoi used to accept
        CHECK(index_load_threads_capped(4, 8) == 64);
    }
    {
        LoadThreadsEnv env("+5");  // strtol accepts a leading '+'
        CHECK(index_load_threads_capped(4, 8) == 5);
    }
}

TEST_CASE("index_load_threads: override_applied is true only for a valid override"
          * doctest::test_suite("unit/index_load_threads")) {
    // The --huge-pages load reads with more workers unless the user capped them:
    // it must treat only a value index_load_threads accepted as a cap, not an
    // unset, empty or malformed one that index_load_threads ignored.
    struct Case { const char *value; bool applied; };
    for (const Case c : {Case{nullptr, false}, Case{"", false}, Case{"abc", false},
                         Case{"0", false}, Case{"-3", false}, Case{"4x", false},
                         Case{"3", true}, Case{"100000", true}}) {
        CAPTURE(c.value ? c.value : "(unset)");
        LoadThreadsEnv env(c.value);
        bool applied = !c.applied;  // poisoned: the call must always write it
        const int t = index_load_threads_capped(32, 8, &applied);
        CHECK(applied == c.applied);
        CHECK(t == (c.applied ? (c.value[0] == '3' ? 3 : 64) : 8));
    }
}
