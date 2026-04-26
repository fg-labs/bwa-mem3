#ifndef BWA_SYSTEM_H
#define BWA_SYSTEM_H

#include <cstdint>

namespace bwa {

// Total memory available to the process in bytes. On Linux returns
// min(cgroup v2 memory.max, cgroup v1 memory.limit_in_bytes, physical RAM).
// Cgroup "unlimited" sentinels are ignored and the physical RAM is used.
// On macOS returns sysctl hw.memsize. Returns -1 only on catastrophic
// failure (no physical RAM reading available).
int64_t detect_total_memory_bytes();

// Number of CPUs available to the process. On Linux returns
// min(cgroup v2 cpu.max quota/period, cgroup v1 cfs_quota/cfs_period,
// physical CPUs). On macOS returns sysconf(_SC_NPROCESSORS_ONLN).
// Always returns >= 1.
int detect_cpu_count();

namespace system_detail {

// Parse cgroup v2 memory.max content. Returns -1 for "max" or unparseable
// input, else the limit in bytes. Kernel "unlimited" sentinel (>= 1<<62)
// is collapsed to -1 for v1-compat semantics.
int64_t parse_cgroup_memory_max(const char* text);

// Parse cgroup v2 cpu.max ("<quota> <period>"). Returns -1 when unlimited
// or unparseable, else ceil(quota / period) CPUs (min 1 when positive).
int parse_cgroup_cpu_max(const char* text);

// Parse cgroup v1 memory.limit_in_bytes. Returns -1 when at the kernel
// unlimited sentinel, else the limit in bytes.
int64_t parse_cgroup_v1_memory_limit(const char* text);

// Parse cgroup v1 CFS CPU budget. quota_text may be "-1" (unlimited);
// period_text must be a positive integer. Returns -1 when unlimited or
// unparseable, else ceil(quota / period).
int parse_cgroup_v1_cpu(const char* quota_text, const char* period_text);

} // namespace system_detail

} // namespace bwa

#endif
