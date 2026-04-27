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
