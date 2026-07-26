#include "system.h"

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <climits>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <unistd.h>

#ifdef __APPLE__
#include <sys/sysctl.h>
#include <sys/types.h>
#endif

namespace bwa {

std::string fmt_bytes(int64_t b) {
    char out[32];
    if      (b >= (1LL << 30)) std::snprintf(out, sizeof(out), "%.1f GiB", (double)b / (double)(1LL << 30));
    else if (b >= (1LL << 20)) std::snprintf(out, sizeof(out), "%.1f MiB", (double)b / (double)(1LL << 20));
    else                       std::snprintf(out, sizeof(out), "%lld B",   (long long)b);
    return std::string(out);
}

namespace {

// Slurp up to 4 KiB from `path` into `out`; returns true on success.
//
// Intentionally bounded: this helper exists for cgroup pseudo-files
// (memory.max, cpu.max, cpu.cfs_quota_us, …), all of which fit comfortably
// under 4 KiB. Returns false if the file would be truncated so callers
// can't accidentally consume garbage from a larger file.
bool read_small_file(const char* path, std::string& out) {
    FILE* fp = std::fopen(path, "r");
    if (!fp) return false;
    char buf[4096];
    size_t n = std::fread(buf, 1, sizeof(buf) - 1, fp);
    bool truncated = (n == sizeof(buf) - 1) && !std::feof(fp);
    std::fclose(fp);
    if (truncated) return false;
    buf[n] = '\0';
    out.assign(buf, n);
    return true;
}

// Parse a base-10 non-negative integer from a possibly-padded string. Returns
// -1 on any parse failure (empty, negative, trailing garbage other than
// whitespace). Treats values >= 2^62 as "unlimited" (kernel sentinel).
int64_t parse_u64_or_unlimited(const char* s) {
    while (*s && std::isspace((unsigned char)*s)) ++s;
    if (!*s) return -1;
    char* end = nullptr;
    errno = 0;
    long long v = std::strtoll(s, &end, 10);
    if (errno != 0 || end == s || v < 0) return -1;
    while (*end && std::isspace((unsigned char)*end)) ++end;
    if (*end) return -1;
    if (v >= (1LL << 62)) return -1;
    return (int64_t)v;
}

// Ceil-div for positive inputs. Caller guards period > 0.
int ceil_div_to_int(int64_t quota, int64_t period) {
    int64_t c = (quota + period - 1) / period;
    if (c < 1) c = 1;
    if (c > INT_MAX) c = INT_MAX;
    return (int)c;
}

int64_t detect_physical_memory_bytes() {
#ifdef __APPLE__
    int64_t mem = 0;
    size_t len = sizeof(mem);
    int mib[2] = { CTL_HW, HW_MEMSIZE };
    if (sysctl(mib, 2, &mem, &len, nullptr, 0) == 0 && mem > 0) return mem;
    return -1;
#else
    long pages = sysconf(_SC_PHYS_PAGES);
    long page  = sysconf(_SC_PAGE_SIZE);
    if (pages > 0 && page > 0) return (int64_t)pages * (int64_t)page;
    return -1;
#endif
}

int detect_physical_cpu_count() {
    long n = sysconf(_SC_NPROCESSORS_ONLN);
    if (n < 1) return 1;
    if (n > INT_MAX) return INT_MAX;
    return (int)n;
}

// Linux-only: read the cgroup v2 and v1 memory caps, returning the smallest
// finite value, or -1 if unlimited / unavailable.
int64_t detect_cgroup_memory_bytes() {
    int64_t best = -1;
    std::string buf;
    if (read_small_file("/sys/fs/cgroup/memory.max", buf)) {
        int64_t v = system_detail::parse_cgroup_memory_max(buf.c_str());
        if (v > 0) best = v;
    }
    if (read_small_file("/sys/fs/cgroup/memory/memory.limit_in_bytes", buf)) {
        int64_t v = system_detail::parse_cgroup_v1_memory_limit(buf.c_str());
        if (v > 0 && (best < 0 || v < best)) best = v;
    }
    return best;
}

int detect_cgroup_cpu_count() {
    int best = -1;
    std::string buf;
    if (read_small_file("/sys/fs/cgroup/cpu.max", buf)) {
        int v = system_detail::parse_cgroup_cpu_max(buf.c_str());
        if (v > 0) best = v;
    }
    std::string q, p;
    bool have_v1_quota  = read_small_file("/sys/fs/cgroup/cpu/cpu.cfs_quota_us",  q);
    bool have_v1_period = read_small_file("/sys/fs/cgroup/cpu/cpu.cfs_period_us", p);
    if (have_v1_quota && have_v1_period) {
        int v = system_detail::parse_cgroup_v1_cpu(q.c_str(), p.c_str());
        if (v > 0 && (best < 0 || v < best)) best = v;
    }
    return best;
}

} // anonymous namespace

namespace system_detail {

int64_t parse_cgroup_memory_max(const char* text) {
    while (*text && std::isspace((unsigned char)*text)) ++text;
    if (std::strncmp(text, "max", 3) == 0) {
        const char* p = text + 3;
        while (*p && std::isspace((unsigned char)*p)) ++p;
        if (*p == '\0') return -1;
    }
    return parse_u64_or_unlimited(text);
}

int parse_cgroup_cpu_max(const char* text) {
    while (*text && std::isspace((unsigned char)*text)) ++text;
    if (std::strncmp(text, "max", 3) == 0 &&
        (text[3] == '\0' || std::isspace((unsigned char)text[3]))) {
        return -1;
    }
    char* end = nullptr;
    errno = 0;
    long long quota = std::strtoll(text, &end, 10);
    if (errno != 0 || end == text || quota <= 0) return -1;
    while (*end && std::isspace((unsigned char)*end)) ++end;
    if (!*end) return -1;
    const char* period_text = end;
    errno = 0;
    long long period = std::strtoll(period_text, &end, 10);
    if (errno != 0 || end == period_text || period <= 0) return -1;
    // Reject trailing non-whitespace garbage (mirrors parse_cgroup_v1_cpu).
    while (*end && std::isspace((unsigned char)*end)) ++end;
    if (*end) return -1;
    return ceil_div_to_int(quota, period);
}

int64_t parse_cgroup_v1_memory_limit(const char* text) {
    return parse_u64_or_unlimited(text);
}

int parse_cgroup_v1_cpu(const char* quota_text, const char* period_text) {
    while (*quota_text && std::isspace((unsigned char)*quota_text)) ++quota_text;
    if (*quota_text == '-') return -1; // "-1" or any negative means unlimited
    char* end = nullptr;
    errno = 0;
    long long quota = std::strtoll(quota_text, &end, 10);
    if (errno != 0 || end == quota_text || quota <= 0) return -1;
    while (*end && std::isspace((unsigned char)*end)) ++end;
    if (*end) return -1;
    errno = 0;
    long long period = std::strtoll(period_text, &end, 10);
    if (errno != 0 || end == period_text || period <= 0) return -1;
    while (*end && std::isspace((unsigned char)*end)) ++end;
    if (*end) return -1;
    return ceil_div_to_int(quota, period);
}

} // namespace system_detail

int64_t detect_total_memory_bytes() {
    int64_t phys   = detect_physical_memory_bytes();
    int64_t cgroup = detect_cgroup_memory_bytes();
    if (phys < 0 && cgroup < 0) return -1;
    if (phys < 0) return cgroup;
    if (cgroup < 0) return phys;
    return std::min(phys, cgroup);
}

int detect_cpu_count() {
    int phys   = detect_physical_cpu_count();
    int cgroup = detect_cgroup_cpu_count();
    if (cgroup < 1) return phys;
    int r = std::min(phys, cgroup);
    return r < 1 ? 1 : r;
}

namespace {

// Headroom left to the OS and to whatever else shares the host: the larger of
// an absolute floor and a proportional share.
constexpr int64_t kBatchReserveFloorBytes = 2LL << 30;   // 2 GiB
constexpr int64_t kBatchReserveDivisor    = 20;          // 5% of total

} // namespace

int64_t resolve_batch_memory_budget(int64_t total_bytes) {
    if (total_bytes <= 0) return -1;
    int64_t reserve = std::max(kBatchReserveFloorBytes,
                               total_bytes / kBatchReserveDivisor);
    // Never reserve more than half: on a host smaller than 2x the floor the
    // proportional rule would otherwise resolve to a zero or negative budget.
    reserve = std::min(reserve, total_bytes / 2);
    return total_bytes - reserve;
}

int64_t required_total_for_batch_budget(int64_t budget_bytes) {
    if (budget_bytes <= 0) return -1;
    // The reserve has three regimes -- half of total below 4 GiB, a flat 2 GiB
    // up to 40 GiB, then 5% -- so a closed-form inverse needs three branches
    // plus their boundary cases, and any of them can silently disagree with the
    // forward policy. Binary search the forward function instead: it is
    // monotonic non-decreasing, so the smallest sufficient total is well
    // defined, the answer is exact in every regime by construction, and the
    // inverse cannot drift if the policy above is ever retuned. ~63 iterations,
    // evaluated at most once per refused build.
    int64_t lo = 1, hi = INT64_MAX, best = -1;
    while (lo <= hi) {
        int64_t mid = lo + (hi - lo) / 2;          // no overflow
        if (resolve_batch_memory_budget(mid) >= budget_bytes) {
            best = mid;
            hi = mid - 1;                          // mid >= lo >= 1, so no underflow
        } else {
            if (mid == INT64_MAX) break;           // infeasible: nothing suffices
            lo = mid + 1;
        }
    }
    return best;
}

} // namespace bwa
