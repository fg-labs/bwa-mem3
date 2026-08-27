// test/unit/test_lockstep_width.cpp
//
// Unit tests for the phase-2 SMEM lockstep width resolution: the pure policy
// that turns a raw memory-level-parallelism estimate (or a BWA3_SMEM_LOCKSTEP_N
// override) into the width the lockstep SMEM driver runs at. The measurement
// itself (the startup pointer-chase) and the one-shot installer
// (bwa3_init_smem_lockstep_width, idempotent + global) are validated by the
// whole-run byte-identity and throughput gates, not here; this file pins the
// clamp/override policy those two consult -- including the initializer's
// probe-vs-pin decision, which lives entirely in bwa3_lockstep_width_parse_env.
//
// Every value fed to the functions under test is derived from the configured
// [SMEM_LOCKSTEP_N, SMEM_LOCKSTEP_N_MAX] range, never the 16/64 shipped default,
// so the assertions follow the contract under any valid -D override of either
// width (the compile-time guard in lockstep_width.h bounds that range). Cases
// that only exist for a non-degenerate range (a strict interior value) or that
// would overflow at the INT32_MAX extreme are compile-guarded accordingly.

#include "doctest/doctest.h"
#include "../../src/lockstep_width.h"

#include <cstdint>  // INT32_MAX
#include <cstdio>   // snprintf

// Stringify a macro's *value* (two-level expansion) so the configured widths
// can be fed to parse_env as decimal string literals.
#define BWA3_STR2(x) #x
#define BWA3_STR(x) BWA3_STR2(x)

// --- bwa3_lockstep_width_from_probe: clamp raw MLP to [floor, ceiling] ---

TEST_CASE("lockstep width: from_probe clamps a raw MLP estimate to [floor, ceiling]"
          * doctest::test_suite("unit/smem")) {
    SUBCASE("a non-positive probe (measurement failure) falls back to the floor") {
        CHECK(bwa3_lockstep_width_from_probe(0)  == SMEM_LOCKSTEP_N);
        CHECK(bwa3_lockstep_width_from_probe(-3) == SMEM_LOCKSTEP_N);
    }
    SUBCASE("the floor and ceiling map to themselves") {
        CHECK(bwa3_lockstep_width_from_probe(SMEM_LOCKSTEP_N)     == SMEM_LOCKSTEP_N);
        CHECK(bwa3_lockstep_width_from_probe(SMEM_LOCKSTEP_N_MAX) == SMEM_LOCKSTEP_N_MAX);
    }
    SUBCASE("a probe above the ceiling clamps to the max") {
#if SMEM_LOCKSTEP_N_MAX < INT32_MAX
        // Guarded so SMEM_LOCKSTEP_N_MAX + 1 cannot overflow int32_t.
        CHECK(bwa3_lockstep_width_from_probe(SMEM_LOCKSTEP_N_MAX + 1) == SMEM_LOCKSTEP_N_MAX);
#endif
        CHECK(bwa3_lockstep_width_from_probe(INT32_MAX) == SMEM_LOCKSTEP_N_MAX);
    }
#if SMEM_LOCKSTEP_N > 1
    SUBCASE("a probe below the floor clamps up to the floor") {
        CHECK(bwa3_lockstep_width_from_probe(SMEM_LOCKSTEP_N - 1) == SMEM_LOCKSTEP_N);
        CHECK(bwa3_lockstep_width_from_probe(1)                   == SMEM_LOCKSTEP_N);
    }
#endif
#if SMEM_LOCKSTEP_N_MAX > SMEM_LOCKSTEP_N
    SUBCASE("a probe inside the range is used verbatim") {
        constexpr int32_t mid = SMEM_LOCKSTEP_N + (SMEM_LOCKSTEP_N_MAX - SMEM_LOCKSTEP_N) / 2;
        CHECK(bwa3_lockstep_width_from_probe(mid) == mid);
    }
#endif
}

// --- bwa3_lockstep_width_parse_env: classify the BWA3_SMEM_LOCKSTEP_N override ---
//
// Tri-state contract the initializer keys off: > 0 pins the width (skip probe),
// 0 means unset (no override -- keep the compile-time default, probe only if
// opted in), -1 means set-but-invalid (report, then resolve as if unset).
// The regression this guards: an invalid override must NOT resolve to a usable
// width -- previously a nonempty-but-garbage value silently disabled auto-tuning
// by flooring, and an overflowed value silently selected the max.

TEST_CASE("lockstep width: parse_env classifies the override into pin/unset/invalid"
          * doctest::test_suite("unit/smem")) {
    SUBCASE("unset or empty is 0 (no override; caller keeps default, probes only if opted in)") {
        CHECK(bwa3_lockstep_width_parse_env(nullptr) == 0);
        CHECK(bwa3_lockstep_width_parse_env("")      == 0);
    }
    SUBCASE("a valid integer at the floor or ceiling pins that width") {
        CHECK(bwa3_lockstep_width_parse_env(BWA3_STR(SMEM_LOCKSTEP_N))     == SMEM_LOCKSTEP_N);
        CHECK(bwa3_lockstep_width_parse_env(BWA3_STR(SMEM_LOCKSTEP_N_MAX)) == SMEM_LOCKSTEP_N_MAX);
        CHECK(bwa3_lockstep_width_parse_env("1") == 1);  // escape hatch: env may go below the floor
    }
#if SMEM_LOCKSTEP_N_MAX > SMEM_LOCKSTEP_N
    SUBCASE("a valid integer inside the range pins that width") {
        constexpr int32_t mid = SMEM_LOCKSTEP_N + (SMEM_LOCKSTEP_N_MAX - SMEM_LOCKSTEP_N) / 2;
        char buf[16];
        snprintf(buf, sizeof(buf), "%d", mid);
        CHECK(bwa3_lockstep_width_parse_env(buf) == mid);
    }
#endif
#if SMEM_LOCKSTEP_N_MAX < INT32_MAX
    SUBCASE("a value above the ceiling clamps to the max") {
        char buf[16];
        snprintf(buf, sizeof(buf), "%lld", (long long)SMEM_LOCKSTEP_N_MAX + 1);
        CHECK(bwa3_lockstep_width_parse_env(buf) == SMEM_LOCKSTEP_N_MAX);
    }
#endif
    SUBCASE("a malformed value is -1 (invalid: reported, not a silent floor)") {
        CHECK(bwa3_lockstep_width_parse_env("garbage") == -1);
        CHECK(bwa3_lockstep_width_parse_env("12x")     == -1);
    }
    SUBCASE("a non-positive value is -1 (invalid: reported, not a silent floor)") {
        CHECK(bwa3_lockstep_width_parse_env("0")  == -1);
        CHECK(bwa3_lockstep_width_parse_env("-5") == -1);
    }
    SUBCASE("an overflowed value is -1 (invalid: reported, not a silent clamp to max)") {
        CHECK(bwa3_lockstep_width_parse_env("99999999999999999999999999") == -1);
    }
}

// --- bwa3_lockstep_probe_enabled: gate the startup MLP probe opt-in ---
//
// BWA3_SMEM_LOCKSTEP_PROBE is a truthy opt-in: the startup sweep runs only for
// a value that is present, non-empty, and not "0". Unset, empty, and "0" all
// leave the probe off (the shipped default keeps the compile-time width). This
// pins the gating the initializer keys off so an empty or "0" value can never
// silently pay the probe's startup cost.

TEST_CASE("lockstep width: probe_enabled is a truthy opt-in (unset/empty/\"0\" disable)"
          * doctest::test_suite("unit/smem")) {
    SUBCASE("unset, empty, or \"0\" leaves the probe disabled") {
        CHECK(bwa3_lockstep_probe_enabled(nullptr) == 0);
        CHECK(bwa3_lockstep_probe_enabled("")      == 0);
        CHECK(bwa3_lockstep_probe_enabled("0")     == 0);
    }
    SUBCASE("any other value enables the probe") {
        CHECK(bwa3_lockstep_probe_enabled("1")   != 0);
        CHECK(bwa3_lockstep_probe_enabled("2")   != 0);
        CHECK(bwa3_lockstep_probe_enabled("yes") != 0);
    }
}
