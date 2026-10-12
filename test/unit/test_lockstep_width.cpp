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
#include "../../src/lockstep_width_internal.h"
#include "../../src/simd_dispatch.h"

#include <cstdint>  // INT32_MAX
#include <cstdio>   // snprintf
#include <cstdlib>  // getenv, setenv, unsetenv
#include <initializer_list>
#include <string>

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

// --- third-pass bwtseed lockstep: on/off policy ---
//
// On by default everywhere; the env pin is a strict 0/1 so a typo can never
// silently flip the driver. The per-run installer is exercised below; the two
// drivers' whole-run output identity is pinned by
// test/regression/bwtseed_lockstep_driver_identity.sh.

TEST_CASE("bwtseed lockstep: parse_env is a strict 0/1 pin with unset and invalid distinct"
          * doctest::test_suite("unit/smem")) {
    SUBCASE("unset or empty is -1 (the compiled default)") {
        CHECK(bwa3_bwtseed_lockstep_parse_env(nullptr) == -1);
        CHECK(bwa3_bwtseed_lockstep_parse_env("")      == -1);
    }
    SUBCASE("\"0\" and \"1\" pin off and on") {
        CHECK(bwa3_bwtseed_lockstep_parse_env("0") == 0);
        CHECK(bwa3_bwtseed_lockstep_parse_env("1") == 1);
    }
    SUBCASE("anything else is -2 (invalid: reported, then the default)") {
        CHECK(bwa3_bwtseed_lockstep_parse_env("2")    == -2);
        CHECK(bwa3_bwtseed_lockstep_parse_env("yes")  == -2);
        CHECK(bwa3_bwtseed_lockstep_parse_env("01")   == -2);
        CHECK(bwa3_bwtseed_lockstep_parse_env(" 1")   == -2);
    }
}

TEST_CASE("bwtseed lockstep: init resolves default from the compile-time default, "
          "not a leaked prior value"
          * doctest::test_suite("unit/smem")) {
    // Regression: bwa3_init_bwtseed_lockstep must fall back to the compile-time
    // default (on), never the mutable g_bwtseed_lockstep. Otherwise an explicit
    // pin leaks into a later run whose environment is unset or invalid, which must
    // revert to the compiled default.
    const int32_t saved_global = g_bwtseed_lockstep;
    const char *saved_env = getenv("BWA3_BWTSEED_LOCKSTEP");
    const bool had_env = saved_env != NULL;
    const std::string saved_env_val = had_env ? saved_env : "";

    setenv("BWA3_BWTSEED_LOCKSTEP", "0", 1);
    bwa3_init_bwtseed_lockstep();
    CHECK(g_bwtseed_lockstep == 0);   // explicit pin honored

    unsetenv("BWA3_BWTSEED_LOCKSTEP");
    bwa3_init_bwtseed_lockstep();
    CHECK(g_bwtseed_lockstep == 1);   // reverts to the compiled default, not the leaked pin

    setenv("BWA3_BWTSEED_LOCKSTEP", "0", 1);
    bwa3_init_bwtseed_lockstep();
    REQUIRE(g_bwtseed_lockstep == 0);
    setenv("BWA3_BWTSEED_LOCKSTEP", "2", 1);   // invalid: reported, then the default
    bwa3_init_bwtseed_lockstep();
    CHECK(g_bwtseed_lockstep == 1);

    if (had_env) setenv("BWA3_BWTSEED_LOCKSTEP", saved_env_val.c_str(), 1);
    else unsetenv("BWA3_BWTSEED_LOCKSTEP");
    g_bwtseed_lockstep = saved_global;
}

// --- bwa3_bwtseed_width_parse_env: classify the BWA3_BWTSEED_LOCKSTEP_N override ---
//
// The third-pass lockstep WIDTH override (distinct from the 0/1 on/off pin above).
// Values are derived from the configured [1, BWTSEED_LOCKSTEP_N_MAX] range, never the
// arch-specific default, so the assertions hold under any -D override of either bound
// (the compile guard in lockstep_width.h enforces 1 <= N <= N_MAX). The idempotent
// installer bwa3_init_bwtseed_lockstep_width is validated by the whole-run byte-identity
// gate, not here; this pins the pure policy it consults.

TEST_CASE("bwtseed width: parse_env classifies the override into pin/unset/invalid"
          * doctest::test_suite("unit/smem")) {
    SUBCASE("unset or empty is 0 (keep the arch default)") {
        CHECK(bwa3_bwtseed_width_parse_env(nullptr) == 0);
        CHECK(bwa3_bwtseed_width_parse_env("")      == 0);
    }
    SUBCASE("a positive integer in range is that width; 1 is the floor escape hatch") {
        CHECK(bwa3_bwtseed_width_parse_env(BWA3_STR(BWTSEED_LOCKSTEP_N))     == BWTSEED_LOCKSTEP_N);
        CHECK(bwa3_bwtseed_width_parse_env(BWA3_STR(BWTSEED_LOCKSTEP_N_MAX)) == BWTSEED_LOCKSTEP_N_MAX);
        CHECK(bwa3_bwtseed_width_parse_env("1") == 1);
    }
#if BWTSEED_LOCKSTEP_N_MAX > 2
    SUBCASE("a value strictly inside the range is used verbatim") {
        char buf[16];
        snprintf(buf, sizeof(buf), "%d", BWTSEED_LOCKSTEP_N_MAX - 1);
        CHECK(bwa3_bwtseed_width_parse_env(buf) == BWTSEED_LOCKSTEP_N_MAX - 1);
    }
#endif
    SUBCASE("a value above the ceiling clamps to the max") {
        char buf[16];
        snprintf(buf, sizeof(buf), "%d", BWTSEED_LOCKSTEP_N_MAX + 1);
        CHECK(bwa3_bwtseed_width_parse_env(buf) == BWTSEED_LOCKSTEP_N_MAX);
    }
    SUBCASE("a non-integer is invalid (-1)") {
        CHECK(bwa3_bwtseed_width_parse_env("garbage") == -1);
        CHECK(bwa3_bwtseed_width_parse_env("12x")     == -1);
    }
    SUBCASE("zero or negative is invalid (-1)") {
        CHECK(bwa3_bwtseed_width_parse_env("0")  == -1);
        CHECK(bwa3_bwtseed_width_parse_env("-5") == -1);
    }
    SUBCASE("an overflowing value is invalid (-1)") {
        CHECK(bwa3_bwtseed_width_parse_env("99999999999999999999999999") == -1);
    }
}

// --- host-tuned default widths (AMD Zen 5) ---
//
// The CPU identity is injected, so the policy is tested on any build host.
// Expected values are derived with the same clamp the policy documents
// ([compile-time default, *_LOCKSTEP_N_MAX]), so the assertions hold under any
// valid -D override of the widths or their caps.

namespace {
constexpr int32_t clamp_w(int32_t w, int32_t lo, int32_t hi) {
    return w < lo ? lo : (w > hi ? hi : w);
}
constexpr int32_t kZen5Smem =
    clamp_w(BWA3_ZEN5_SMEM_LOCKSTEP_N, SMEM_LOCKSTEP_N, SMEM_LOCKSTEP_N_MAX);
constexpr int32_t kZen5Bwtseed =
    clamp_w(BWA3_ZEN5_BWTSEED_LOCKSTEP_N, BWTSEED_LOCKSTEP_N, BWTSEED_LOCKSTEP_N_MAX);
// The source a host-default width must report: HOST only when it moved the width.
constexpr int kZen5SmemSrc = kZen5Smem != SMEM_LOCKSTEP_N
    ? BWA3_LOCKSTEP_WIDTH_SRC_HOST : BWA3_LOCKSTEP_WIDTH_SRC_DEFAULT;
constexpr int kZen5BwtseedSrc = kZen5Bwtseed != BWTSEED_LOCKSTEP_N
    ? BWA3_LOCKSTEP_WIDTH_SRC_HOST : BWA3_LOCKSTEP_WIDTH_SRC_DEFAULT;
constexpr int kSmemMaxSrc = SMEM_LOCKSTEP_N_MAX != SMEM_LOCKSTEP_N
    ? BWA3_LOCKSTEP_WIDTH_SRC_HOST : BWA3_LOCKSTEP_WIDTH_SRC_DEFAULT;
constexpr int kBwtseedMaxSrc = BWTSEED_LOCKSTEP_N_MAX != BWTSEED_LOCKSTEP_N
    ? BWA3_LOCKSTEP_WIDTH_SRC_HOST : BWA3_LOCKSTEP_WIDTH_SRC_DEFAULT;
// CPUID leaf 1 EAX from its fields.
constexpr uint32_t leaf1(uint32_t base_fam, uint32_t ext_fam,
                         uint32_t base_model = 0, uint32_t ext_model = 0) {
    return (ext_fam << 20) | (ext_model << 16) | (base_fam << 8) | (base_model << 4);
}
// A CPUID vendor register from its four characters (low byte first), independent
// of the build host's byte order.
constexpr uint32_t reg4(const char (&c)[5]) {
    return (uint32_t)(unsigned char)c[0] | ((uint32_t)(unsigned char)c[1] << 8) |
           ((uint32_t)(unsigned char)c[2] << 16) | ((uint32_t)(unsigned char)c[3] << 24);
}
// The Linux kernel's family-1Ah Zen 5 model ranges (bsp_init_amd), restated
// independently of the implementation.
bool kernel_zen5_model(uint32_t m) {
    return m <= 0x2f || (m >= 0x40 && m <= 0x4f) || (m >= 0x60 && m <= 0x7f) ||
           (m >= 0xd0 && m <= 0xd7);
}
bool kernel_zen6_model(uint32_t m) {
    return (m >= 0x50 && m <= 0x5f) || (m >= 0x80 && m <= 0xaf) ||
           (m >= 0xc0 && m <= 0xcf) || (m >= 0xd8 && m <= 0xef);
}
constexpr uint32_t kZen5Model = 0x02;  // EPYC 9R45 (Turin), the measured part
}  // namespace

TEST_CASE("host width: x86 vendor check matches only AuthenticAMD in EBX/EDX/ECX order"
          * doctest::test_suite("unit/smem")) {
    // "AuthenticAMD" = EBX "Auth", EDX "enti", ECX "cAMD".
    CHECK(bwa3_x86_vendor_is_amd(reg4("Auth"), reg4("cAMD"), reg4("enti")) != 0);
    CHECK(bwa3_x86_vendor_is_amd(0x68747541u, 0x444d4163u, 0x69746e65u) != 0);
    // "GenuineIntel" = EBX "Genu", EDX "ineI", ECX "ntel".
    CHECK(bwa3_x86_vendor_is_amd(reg4("Genu"), reg4("ntel"), reg4("ineI")) == 0);
    // "HygonGenuine" (family 18h) is not AMD for this policy.
    CHECK(bwa3_x86_vendor_is_amd(reg4("Hygo"), reg4("uine"), reg4("nGen")) == 0);
    // Swapped ECX/EDX (the string read in EBX, ECX, EDX order) is rejected.
    CHECK(bwa3_x86_vendor_is_amd(reg4("Auth"), reg4("enti"), reg4("cAMD")) == 0);
    CHECK(bwa3_x86_vendor_is_amd(0, 0, 0) == 0);
}

TEST_CASE("host width: x86 family and model decode"
          * doctest::test_suite("unit/smem")) {
    SUBCASE("family: base + extended family only when the base family is 0xF") {
        CHECK(bwa3_x86_cpu_family(leaf1(0xf, 0x0b)) == 0x1au);  // Zen 5 / Zen 6
        CHECK(bwa3_x86_cpu_family(leaf1(0xf, 0x0a)) == 0x19u);  // Zen 3 / Zen 4
        CHECK(bwa3_x86_cpu_family(leaf1(0xf, 0x08)) == 0x17u);  // Zen 1 / Zen 2
        CHECK(bwa3_x86_cpu_family(leaf1(0x6, 0x00, 0xf)) == 0x6u);
        CHECK(bwa3_x86_cpu_family(leaf1(0x6, 0x0b)) == 0x6u);
    }
    SUBCASE("model: extended model is the high nibble for base family 0x6 / 0xF") {
        CHECK(bwa3_x86_cpu_model(leaf1(0xf, 0x0b, 0x2, 0x0)) == 0x02u);
        CHECK(bwa3_x86_cpu_model(leaf1(0xf, 0x0b, 0x4, 0x4)) == 0x44u);
        CHECK(bwa3_x86_cpu_model(leaf1(0xf, 0x0b, 0x0, 0x5)) == 0x50u);
        CHECK(bwa3_x86_cpu_model(leaf1(0xf, 0x0b, 0x7, 0xd)) == 0xd7u);
        CHECK(bwa3_x86_cpu_model(leaf1(0x6, 0x00, 0xf, 0x8)) == 0x8fu);  // Intel SPR
        // Any other base family uses the base model alone.
        CHECK(bwa3_x86_cpu_model(leaf1(0x5, 0x00, 0x4, 0x3)) == 0x04u);
    }
    SUBCASE("real signatures") {
        // EPYC 9R45 (Zen 5, the measured part): family 1Ah, model 02h, stepping 1.
        CHECK(bwa3_x86_cpu_family(0x00b00f21u) == 0x1au);
        CHECK(bwa3_x86_cpu_model(0x00b00f21u)  == 0x02u);
        // Ryzen 9000 (Zen 5, Granite Ridge): family 1Ah, model 44h, stepping 0.
        CHECK(bwa3_x86_cpu_family(0x00b40f40u) == 0x1au);
        CHECK(bwa3_x86_cpu_model(0x00b40f40u)  == 0x44u);
        // EPYC 9004 (Zen 4, Genoa): family 19h, model 11h, stepping 1.
        CHECK(bwa3_x86_cpu_family(0x00a10f11u) == 0x19u);
        CHECK(bwa3_x86_cpu_model(0x00a10f11u)  == 0x11u);
    }
}

TEST_CASE("host width: AMD family 1Ah models match the kernel's Zen 5 ranges"
          * doctest::test_suite("unit/smem")) {
    SUBCASE("range boundaries") {
        for (uint32_t m : {0x00u, 0x02u, 0x2fu, 0x40u, 0x44u, 0x4fu, 0x60u, 0x7fu, 0xd0u, 0xd7u})
            CHECK_MESSAGE(bwa3_amd_family_1ah_model_is_zen5(m) != 0, "model " << m);
        for (uint32_t m : {0x30u, 0x3fu, 0x50u, 0x5fu, 0x80u, 0xafu, 0xb0u, 0xbfu,
                           0xc0u, 0xcfu, 0xd8u, 0xefu, 0xf0u, 0xffu, 0x100u})
            CHECK_MESSAGE(bwa3_amd_family_1ah_model_is_zen5(m) == 0, "model " << m);
    }
    SUBCASE("every model: Zen 5 iff the kernel says Zen 5; never Zen 6") {
        for (uint32_t m = 0; m < 0x100u; m++) {
            const bool z5 = bwa3_amd_family_1ah_model_is_zen5(m) != 0;
            CHECK_MESSAGE(z5 == kernel_zen5_model(m), "model " << m);
            if (kernel_zen6_model(m)) CHECK_MESSAGE(!z5, "Zen 6 model " << m);
        }
    }
    SUBCASE("vendor and family gate the model check") {
        CHECK(bwa3_cpu_is_amd_zen5(1, 0x1a, kZen5Model) != 0);
        CHECK(bwa3_cpu_is_amd_zen5(0, 0x1a, kZen5Model) == 0);  // non-AMD
        CHECK(bwa3_cpu_is_amd_zen5(1, 0x19, kZen5Model) == 0);  // Zen 3 / Zen 4 family
        CHECK(bwa3_cpu_is_amd_zen5(1, 0x1a, 0x50) == 0);        // Zen 6
        CHECK(bwa3_cpu_is_amd_zen5(1, 0x1a, 0x30) == 0);        // unlisted 1Ah model
    }
}

TEST_CASE("host width: AMD Zen 5 gets the Zen 5 widths, every other host the defaults"
          * doctest::test_suite("unit/smem")) {
    SUBCASE("AMD Zen 5 (family 1Ah, Zen 5 models)") {
        for (uint32_t m : {0x02u, 0x11u, 0x24u, 0x44u, 0x70u, 0xd0u}) {
            CHECK(bwa3_smem_host_default_width(1, 0x1a, m)    == kZen5Smem);
            CHECK(bwa3_bwtseed_host_default_width(1, 0x1a, m) == kZen5Bwtseed);
        }
    }
    SUBCASE("AMD family 1Ah Zen 6 and unlisted models keep the defaults") {
        for (uint32_t m : {0x50u, 0x80u, 0xafu, 0xc0u, 0xd8u, 0xefu, 0x30u, 0xb0u, 0xf0u}) {
            CHECK(bwa3_smem_host_default_width(1, 0x1a, m)    == SMEM_LOCKSTEP_N);
            CHECK(bwa3_bwtseed_host_default_width(1, 0x1a, m) == BWTSEED_LOCKSTEP_N);
        }
    }
    SUBCASE("AMD family 19h (Zen 3 / Zen 4) and 17h keep the defaults") {
        CHECK(bwa3_smem_host_default_width(1, 0x19, 0x11)    == SMEM_LOCKSTEP_N);
        CHECK(bwa3_bwtseed_host_default_width(1, 0x19, 0x11) == BWTSEED_LOCKSTEP_N);
        CHECK(bwa3_smem_host_default_width(1, 0x17, 0x31)    == SMEM_LOCKSTEP_N);
        CHECK(bwa3_bwtseed_host_default_width(1, 0x17, 0x31) == BWTSEED_LOCKSTEP_N);
    }
    SUBCASE("a non-AMD vendor keeps the defaults, even at family 1Ah") {
        CHECK(bwa3_smem_host_default_width(0, 0x6, 0x8f)    == SMEM_LOCKSTEP_N);  // Intel
        CHECK(bwa3_bwtseed_host_default_width(0, 0x6, 0x8f) == BWTSEED_LOCKSTEP_N);
        CHECK(bwa3_smem_host_default_width(0, 0x1a, kZen5Model)    == SMEM_LOCKSTEP_N);
        CHECK(bwa3_bwtseed_host_default_width(0, 0x1a, kZen5Model) == BWTSEED_LOCKSTEP_N);
    }
    SUBCASE("a non-x86 host (identity 0/0/0) keeps the defaults") {
        CHECK(bwa3_smem_host_default_width(0, 0, 0)    == SMEM_LOCKSTEP_N);
        CHECK(bwa3_bwtseed_host_default_width(0, 0, 0) == BWTSEED_LOCKSTEP_N);
    }
#if SMEM_LOCKSTEP_N <= 32 && SMEM_LOCKSTEP_N_MAX >= 32 && \
    BWTSEED_LOCKSTEP_N <= 32 && BWTSEED_LOCKSTEP_N_MAX >= 32 && \
    BWA3_ZEN5_SMEM_LOCKSTEP_N == 32 && BWA3_ZEN5_BWTSEED_LOCKSTEP_N == 32
    SUBCASE("shipped configuration: Zen 5 is 32/32") {
        CHECK(bwa3_smem_host_default_width(1, 0x1a, kZen5Model)    == 32);
        CHECK(bwa3_bwtseed_host_default_width(1, 0x1a, kZen5Model) == 32);
    }
#endif
    SUBCASE("the host default stays inside [compile-time default, cap]") {
        const int32_t s = bwa3_smem_host_default_width(1, 0x1a, kZen5Model);
        const int32_t b = bwa3_bwtseed_host_default_width(1, 0x1a, kZen5Model);
        CHECK(s >= SMEM_LOCKSTEP_N);
        CHECK(s <= SMEM_LOCKSTEP_N_MAX);
        CHECK(b >= BWTSEED_LOCKSTEP_N);
        CHECK(b <= BWTSEED_LOCKSTEP_N_MAX);
    }
}

// The resolvers return {width, source}; check both at once.
#define CHECK_RES(expr, want_width, want_src)        \
    do {                                             \
        const bwa3_lockstep_width_res_t r_ = (expr); \
        CHECK(r_.width == (want_width));             \
        CHECK(r_.src == (want_src));                 \
    } while (0)

TEST_CASE("host width: SMEM resolution precedence is env pin > probe > host > default"
          * doctest::test_suite("unit/smem")) {
    const int32_t zen5 = bwa3_smem_host_default_width(1, 0x1a, kZen5Model);
    const int32_t other = bwa3_smem_host_default_width(1, 0x19, 0x11);
    REQUIRE(other == SMEM_LOCKSTEP_N);
    SUBCASE("a valid pin is taken alone, over the probe and the host default") {
        CHECK_RES(bwa3_resolve_smem_width(8, 1, zen5), 8, BWA3_LOCKSTEP_WIDTH_SRC_ENV);
        CHECK_RES(bwa3_resolve_smem_width(8, 0, zen5), 8, BWA3_LOCKSTEP_WIDTH_SRC_ENV);
        CHECK_RES(bwa3_resolve_smem_width(1, 0, zen5), 1, BWA3_LOCKSTEP_WIDTH_SRC_ENV);
        CHECK_RES(bwa3_resolve_smem_width(SMEM_LOCKSTEP_N, 0, SMEM_LOCKSTEP_N),
                  SMEM_LOCKSTEP_N, BWA3_LOCKSTEP_WIDTH_SRC_ENV);
    }
    SUBCASE("the opted-in probe wins over the host default (unset or invalid pin)") {
        CHECK(bwa3_resolve_smem_width(0, 1, zen5).src  == BWA3_LOCKSTEP_WIDTH_SRC_PROBE);
        CHECK(bwa3_resolve_smem_width(0, 1, other).src == BWA3_LOCKSTEP_WIDTH_SRC_PROBE);
        CHECK(bwa3_resolve_smem_width(-1, 1, zen5).src == BWA3_LOCKSTEP_WIDTH_SRC_PROBE);
    }
    SUBCASE("no pin, no probe: the host default, HOST only when it moved the width") {
        CHECK_RES(bwa3_resolve_smem_width(0, 0, zen5),  kZen5Smem, kZen5SmemSrc);
        CHECK_RES(bwa3_resolve_smem_width(-1, 0, zen5), kZen5Smem, kZen5SmemSrc);
        CHECK_RES(bwa3_resolve_smem_width(0, 0, other),  SMEM_LOCKSTEP_N,
                  BWA3_LOCKSTEP_WIDTH_SRC_DEFAULT);
        CHECK_RES(bwa3_resolve_smem_width(-1, 0, other), SMEM_LOCKSTEP_N,
                  BWA3_LOCKSTEP_WIDTH_SRC_DEFAULT);
    }
    SUBCASE("an out-of-range host width is clamped into [default, cap]") {
        CHECK_RES(bwa3_resolve_smem_width(0, 0, SMEM_LOCKSTEP_N_MAX + 100),
                  SMEM_LOCKSTEP_N_MAX, kSmemMaxSrc);
        CHECK_RES(bwa3_resolve_smem_width(0, 0, 0),  SMEM_LOCKSTEP_N,
                  BWA3_LOCKSTEP_WIDTH_SRC_DEFAULT);
        CHECK_RES(bwa3_resolve_smem_width(0, 0, -5), SMEM_LOCKSTEP_N,
                  BWA3_LOCKSTEP_WIDTH_SRC_DEFAULT);
    }
    SUBCASE("a measured probe result is clamped and labelled PROBE") {
        CHECK_RES(bwa3_resolve_smem_width_from_probe(SMEM_LOCKSTEP_N, zen5),
                  SMEM_LOCKSTEP_N, BWA3_LOCKSTEP_WIDTH_SRC_PROBE);
        CHECK_RES(bwa3_resolve_smem_width_from_probe(1, zen5),
                  SMEM_LOCKSTEP_N, BWA3_LOCKSTEP_WIDTH_SRC_PROBE);
        CHECK_RES(bwa3_resolve_smem_width_from_probe(SMEM_LOCKSTEP_N_MAX + 10, other),
                  SMEM_LOCKSTEP_N_MAX, BWA3_LOCKSTEP_WIDTH_SRC_PROBE);
    }
    SUBCASE("a probe that cannot measure falls back to the host default") {
        CHECK_RES(bwa3_resolve_smem_width_from_probe(0, zen5),  kZen5Smem, kZen5SmemSrc);
        CHECK_RES(bwa3_resolve_smem_width_from_probe(-3, zen5), kZen5Smem, kZen5SmemSrc);
        CHECK_RES(bwa3_resolve_smem_width_from_probe(0, other), SMEM_LOCKSTEP_N,
                  BWA3_LOCKSTEP_WIDTH_SRC_DEFAULT);
    }
}

TEST_CASE("host width: bwtSeed resolution precedence is env pin > host > default"
          * doctest::test_suite("unit/smem")) {
    const int32_t zen5 = bwa3_bwtseed_host_default_width(1, 0x1a, kZen5Model);
    const int32_t other = bwa3_bwtseed_host_default_width(0, 0x6, 0x8f);
    REQUIRE(other == BWTSEED_LOCKSTEP_N);
    SUBCASE("a valid pin is taken alone") {
        CHECK_RES(bwa3_resolve_bwtseed_width(8, zen5), 8, BWA3_LOCKSTEP_WIDTH_SRC_ENV);
        CHECK_RES(bwa3_resolve_bwtseed_width(1, other), 1, BWA3_LOCKSTEP_WIDTH_SRC_ENV);
        CHECK_RES(bwa3_resolve_bwtseed_width(BWTSEED_LOCKSTEP_N, BWTSEED_LOCKSTEP_N),
                  BWTSEED_LOCKSTEP_N, BWA3_LOCKSTEP_WIDTH_SRC_ENV);
    }
    SUBCASE("no pin (unset or invalid): the host default, HOST only when it moved") {
        CHECK_RES(bwa3_resolve_bwtseed_width(0, zen5),  kZen5Bwtseed, kZen5BwtseedSrc);
        CHECK_RES(bwa3_resolve_bwtseed_width(-1, zen5), kZen5Bwtseed, kZen5BwtseedSrc);
        CHECK_RES(bwa3_resolve_bwtseed_width(0, other),  BWTSEED_LOCKSTEP_N,
                  BWA3_LOCKSTEP_WIDTH_SRC_DEFAULT);
        CHECK_RES(bwa3_resolve_bwtseed_width(-1, other), BWTSEED_LOCKSTEP_N,
                  BWA3_LOCKSTEP_WIDTH_SRC_DEFAULT);
    }
    SUBCASE("an out-of-range host width is clamped into [default, cap]") {
        CHECK_RES(bwa3_resolve_bwtseed_width(0, BWTSEED_LOCKSTEP_N_MAX + 100),
                  BWTSEED_LOCKSTEP_N_MAX, kBwtseedMaxSrc);
        CHECK_RES(bwa3_resolve_bwtseed_width(0, 0), BWTSEED_LOCKSTEP_N,
                  BWA3_LOCKSTEP_WIDTH_SRC_DEFAULT);
    }
}

TEST_CASE("host width: every source has its documented log name"
          * doctest::test_suite("unit/smem")) {
    CHECK(std::string(bwa3_lockstep_width_source_name(BWA3_LOCKSTEP_WIDTH_SRC_DEFAULT)) == "default");
    CHECK(std::string(bwa3_lockstep_width_source_name(BWA3_LOCKSTEP_WIDTH_SRC_HOST))
          == "host default: AMD Zen 5");
    CHECK(std::string(bwa3_lockstep_width_source_name(BWA3_LOCKSTEP_WIDTH_SRC_PROBE)) == "startup probe");
    CHECK(std::string(bwa3_lockstep_width_source_name(BWA3_LOCKSTEP_WIDTH_SRC_ENV)) == "env pin");
    // Out-of-range values fall back to "default".
    CHECK(std::string(bwa3_lockstep_width_source_name(99)) == "default");
    CHECK(std::string(bwa3_lockstep_width_source_name(-1)) == "default");
    // The installed-source accessors always report a known source.
    CHECK(bwa3_smem_lockstep_n_source()    >= BWA3_LOCKSTEP_WIDTH_SRC_DEFAULT);
    CHECK(bwa3_smem_lockstep_n_source()    <= BWA3_LOCKSTEP_WIDTH_SRC_ENV);
    CHECK(bwa3_bwtseed_lockstep_n_source() >= BWA3_LOCKSTEP_WIDTH_SRC_DEFAULT);
    CHECK(bwa3_bwtseed_lockstep_n_source() <= BWA3_LOCKSTEP_WIDTH_SRC_ENV);
}

TEST_CASE("host width: the live host identity is consistent with the build arch"
          * doctest::test_suite("unit/smem")) {
    int is_amd = -1;
    uint32_t family = 0xffffffffu, model = 0xffffffffu;
    bwa3_host_cpu_identity(&is_amd, &family, &model);
#if defined(__x86_64__) || defined(__i386__)
    CHECK((is_amd == 0 || is_amd == 1));
    CHECK(family != 0u);   // every x86 CPU reports a nonzero family
    CHECK(model <= 0xffu);
#else
    CHECK(is_amd == 0);    // no CPUID: no host-tuned default applies
    CHECK(family == 0u);
    CHECK(model == 0u);
#endif
    // Repeated reads return the cached identity.
    int is_amd2 = -1;
    uint32_t family2 = 0, model2 = 0;
    bwa3_host_cpu_identity(&is_amd2, &family2, &model2);
    CHECK(is_amd2 == is_amd);
    CHECK(family2 == family);
    CHECK(model2 == model);
}
