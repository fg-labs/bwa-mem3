#ifndef BWA_SYSTEM_H
#define BWA_SYSTEM_H

#include <cstdint>
#include <string>

namespace bwa {

// Format a byte count for human-facing output: "N.N GiB", "N.N MiB", or "N B".
// Sized so a --max-memory of 15M reads as "15.0 MiB" rather than "0.0 GiB".
std::string fmt_bytes(int64_t bytes);

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

// Peak-memory budget for a one-shot batch job, given the total memory
// available to the process (i.e. detect_total_memory_bytes()).
//
// The policy is "the machine, less a headroom reserve" rather than a fixed
// fraction: a batch build exists to consume the host it was given, and a
// fractional split refuses work that the host can plainly do. The reserve is
// max(2 GiB, 5% of total), itself clamped to half of total so that small hosts
// still resolve to a usable budget instead of zero.
//
// Returns -1 when total_bytes <= 0; callers must supply their own fallback.
int64_t resolve_batch_memory_budget(int64_t total_bytes);

// The SMALLEST total memory for which resolve_batch_memory_budget() yields at
// least `budget_bytes` -- exact in every reserve regime, so the "retry on a host
// with >= N" hint it feeds never overstates what a host needs. Returns -1 when
// `budget_bytes` is non-positive, or when no total is sufficient (the reserve is
// always positive, so budgets within ~5% of INT64_MAX are unreachable).
int64_t required_total_for_batch_budget(int64_t budget_bytes);

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
