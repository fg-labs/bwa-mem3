#include "system.h"

#include <cstdint>
#include <cstdio>
#include <cstdlib>

using bwa::system_detail::parse_cgroup_memory_max;
using bwa::system_detail::parse_cgroup_cpu_max;
using bwa::system_detail::parse_cgroup_v1_memory_limit;
using bwa::system_detail::parse_cgroup_v1_cpu;

// Always-on check that doesn't get elided under -DNDEBUG (which would
// strip plain assert() and turn this test into a no-op).
#define CHECK(cond) do { \
    if (!(cond)) { \
        std::fprintf(stderr, "CHECK failed at %s:%d: %s\n", __FILE__, __LINE__, #cond); \
        std::fflush(stderr); \
        std::abort(); \
    } \
} while (0)

int main() {
    // cgroup v2 memory.max
    CHECK(parse_cgroup_memory_max("max\n") == -1);
    CHECK(parse_cgroup_memory_max("  max  \n") == -1);
    CHECK(parse_cgroup_memory_max("12345\n") == 12345);
    CHECK(parse_cgroup_memory_max("1073741824") == 1073741824);
    // kernel "unlimited" sentinel around 2^62-ish — anything >= 2^62 is -1.
    CHECK(parse_cgroup_memory_max("9223372036854771712") == -1);
    CHECK(parse_cgroup_memory_max("bogus") == -1);
    CHECK(parse_cgroup_memory_max("") == -1);

    // cgroup v2 cpu.max
    CHECK(parse_cgroup_cpu_max("max 100000\n") == -1);
    CHECK(parse_cgroup_cpu_max("max\n") == -1);
    CHECK(parse_cgroup_cpu_max("100000 100000\n") == 1);
    CHECK(parse_cgroup_cpu_max("200000 100000\n") == 2);
    CHECK(parse_cgroup_cpu_max("250000 100000\n") == 3);   // ceil(2.5)
    CHECK(parse_cgroup_cpu_max("50000 100000\n") == 1);    // ceil(0.5)
    CHECK(parse_cgroup_cpu_max("100000\n") == -1);         // period missing
    CHECK(parse_cgroup_cpu_max("0 100000\n") == -1);
    CHECK(parse_cgroup_cpu_max("100000 0\n") == -1);
    // Trailing non-whitespace garbage on the period field must be rejected
    // (mirrors parse_cgroup_v1_cpu). Whitespace after period is fine.
    CHECK(parse_cgroup_cpu_max("100000 100000xyz") == -1);
    CHECK(parse_cgroup_cpu_max("100000 100000 \n") == 1);
    // Overflow-safe ceil_div_to_int: a huge period with a small quota has a true
    // ceil of 1, not INT_MAX (regression: the earlier `quota + period - 1`
    // additive-overflow guard returned INT_MAX for this input). A huge quota
    // still clamps to INT_MAX.
    CHECK(parse_cgroup_cpu_max("5 9223372036854775807\n") == 1);
    CHECK(parse_cgroup_cpu_max("9223372036854775807 2\n") == 2147483647);

    // cgroup v1 memory.limit_in_bytes
    CHECK(parse_cgroup_v1_memory_limit("1000000000") == 1000000000);
    CHECK(parse_cgroup_v1_memory_limit("9223372036854771712") == -1);
    CHECK(parse_cgroup_v1_memory_limit("bogus") == -1);

    // cgroup v1 CFS CPU budget
    CHECK(parse_cgroup_v1_cpu("-1", "100000") == -1);
    CHECK(parse_cgroup_v1_cpu("200000", "100000") == 2);
    CHECK(parse_cgroup_v1_cpu("150000", "100000") == 2);   // ceil(1.5)
    CHECK(parse_cgroup_v1_cpu("0", "100000") == -1);
    CHECK(parse_cgroup_v1_cpu("100000", "0") == -1);
    CHECK(parse_cgroup_v1_cpu("abc", "100000") == -1);
    // Trailing whitespace (typical pseudo-file content) must still parse,
    // but trailing non-whitespace garbage must be rejected on both fields.
    CHECK(parse_cgroup_v1_cpu("200000\n", "100000\n") == 2);
    CHECK(parse_cgroup_v1_cpu("200000  ", "100000  ") == 2);
    CHECK(parse_cgroup_v1_cpu("200000abc", "100000") == -1);
    CHECK(parse_cgroup_v1_cpu("200000", "100000xyz") == -1);

    // Public detectors: sanity bounds on whatever host we're running on.
    int64_t mem = bwa::detect_total_memory_bytes();
    CHECK(mem > 0);
    // 1 PiB upper bound is far above any real machine yet still finite, so
    // the assertion catches accidental gigantic returns without flaking on
    // big memory-optimised cloud / HPC instances.
    CHECK(mem <= (int64_t)1 << 50);

    // Batch memory budget policy. The regression this guards: an hg38 build
    // needs ~72 GiB, and the previous `min(50% of RAM, 32 GiB)` default
    // refused it on every host — including a 256 GiB server.
    const int64_t kGiB      = 1LL << 30;
    const int64_t kHg38Need = 77216326020LL;   // (2 * 3217346917 + 1) * 12 B
    using bwa::resolve_batch_memory_budget;
    using bwa::required_total_for_batch_budget;

    CHECK(resolve_batch_memory_budget(0) == -1);
    CHECK(resolve_batch_memory_budget(-1) == -1);

    // Reserve floor dominates below 40 GiB: budget = total - 2 GiB.
    CHECK(resolve_batch_memory_budget(8 * kGiB) == 6 * kGiB);
    CHECK(resolve_batch_memory_budget(32 * kGiB) == 30 * kGiB);
    // Proportional reserve (5%) takes over above 40 GiB.
    CHECK(resolve_batch_memory_budget(40 * kGiB) == 38 * kGiB);
    CHECK(resolve_batch_memory_budget(80 * kGiB) == 76 * kGiB);
    CHECK(resolve_batch_memory_budget(256 * kGiB) == 256 * kGiB - 256 * kGiB / 20);
    // Small hosts still resolve to something usable, never <= 0.
    CHECK(resolve_batch_memory_budget(1 * kGiB) == 512 * 1024 * 1024);
    CHECK(resolve_batch_memory_budget(1) > 0);
    // Monotonic: more RAM never yields a smaller budget.
    for (int64_t t = 1; t <= 512 * kGiB; t = t * 3 / 2 + 1)
        CHECK(resolve_batch_memory_budget(t) <= resolve_batch_memory_budget(t * 3 / 2 + 1));

    // hg38 must build on hosts that can plainly hold it, and must not be
    // promised on hosts that cannot.
    CHECK(resolve_batch_memory_budget(96 * kGiB) >= kHg38Need);
    CHECK(resolve_batch_memory_budget(128 * kGiB) >= kHg38Need);
    CHECK(resolve_batch_memory_budget(256 * kGiB) >= kHg38Need);
    CHECK(resolve_batch_memory_budget(64 * kGiB) < kHg38Need);

    // The inverse is what the error message quotes, so it must be both
    // sufficient AND minimal: quoting a host size larger than necessary sends
    // someone shopping for RAM they do not need.
    CHECK(required_total_for_batch_budget(0) == -1);
    CHECK(required_total_for_batch_budget(-1) == -1);

    // Sufficiency and minimality across all three reserve regimes.
    for (int64_t b = 1; b <= 512 * kGiB; b = b * 3 / 2 + 1) {
        int64_t total = required_total_for_batch_budget(b);
        CHECK(total > 0);
        CHECK(resolve_batch_memory_budget(total) >= b);
        // Minimal: one byte less must not suffice.
        if (total > 1) CHECK(resolve_batch_memory_budget(total - 1) < b);
    }

    // Half-reserve regime (total <= 4 GiB): budget = ceil(total/2), so the
    // smallest sufficient total is 2b-1 -- NOT b + 2 GiB. A 1-byte host
    // resolves to a 1-byte budget.
    CHECK(required_total_for_batch_budget(1) == 1);
    CHECK(required_total_for_batch_budget(2) == 3);
    CHECK(required_total_for_batch_budget(2 * kGiB) == 4 * kGiB - 1);

    // Transition into the fixed 2 GiB reserve, then into the 5% reserve.
    CHECK(required_total_for_batch_budget(2 * kGiB + 1) == 4 * kGiB + 1);
    CHECK(required_total_for_batch_budget(38 * kGiB - 1) == 40 * kGiB - 1);
    CHECK(required_total_for_batch_budget(38 * kGiB) == 40 * kGiB);

    // Large budgets whose required total still fits int64_t must be answered,
    // not rejected: the answer is only ~1.053x the budget.
    const int64_t huge = INT64_MAX / 2;
    const int64_t huge_total = required_total_for_batch_budget(huge);
    CHECK(huge_total > 0);
    CHECK(resolve_batch_memory_budget(huge_total) >= huge);
    // Genuinely infeasible: no total resolves to a budget of INT64_MAX, since
    // the reserve is always positive at that scale.
    CHECK(required_total_for_batch_budget(INT64_MAX) == -1);

    CHECK(required_total_for_batch_budget(kHg38Need) < 80 * kGiB);

    // fmt_bytes: a --max-memory of 15M must not read as "0.0 GiB" (it did).
    CHECK(bwa::fmt_bytes(15 * 1024 * 1024) == "15.0 MiB");
    CHECK(bwa::fmt_bytes(1200 * 1024 * 1024) == "1.2 GiB");
    CHECK(bwa::fmt_bytes(64 * kGiB) == "64.0 GiB");
    CHECK(bwa::fmt_bytes(512) == "512 B");
    CHECK(bwa::fmt_bytes(0) == "0 B");
    // Boundaries: exactly 1 MiB / 1 GiB step up to the larger unit.
    CHECK(bwa::fmt_bytes(1024 * 1024) == "1.0 MiB");
    CHECK(bwa::fmt_bytes(1024 * 1024 - 1) == "1048575 B");
    CHECK(bwa::fmt_bytes(kGiB) == "1.0 GiB");

    int cpu = bwa::detect_cpu_count();
    CHECK(cpu >= 1);
    // 4096 was too tight for HPC / cloud nodes: AWS u-7i.metal-224xl already
    // exposes 896 vCPUs and the trajectory is upward. 65536 keeps the test
    // robust for the foreseeable future without giving up the upper sanity
    // check entirely.
    CHECK(cpu <= (1 << 16));

    std::printf("system_test OK  (mem=%.2f GiB, cpu=%d)\n",
                (double)mem / (1024.0 * 1024.0 * 1024.0), cpu);
    return 0;
}
