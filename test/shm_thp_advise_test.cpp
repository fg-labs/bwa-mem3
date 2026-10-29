/* Verifies that the shm index segment is advised MADV_HUGEPAGE, on both the
 * stager's mapping and the attached one.
 *
 * A disk-loaded index lives in buffers advised MADV_HUGEPAGE (FMI_search.cpp),
 * so its FM-index and SA lookups run on 2 MB pages. The shm segment is a tmpfs
 * file mapping; without the same advice a `huge=advise` /dev/shm mount backs it
 * with 4 KB pages, and the TLB-bound lookups cost about a third more CPU than a
 * disk load. bwa_shm_stage and bwa_shm_attach map the segment through one
 * helper that advises it. The stage side matters most: tmpfs fixes a page's
 * size when it is first written, which is during the stager's fill. The stager
 * unmaps its mapping before it returns, so the test inspects it through
 * bwa_shm_stage_map_hook; the attach side is the mapping the aligner reads.
 *
 * Linux records the advice as the `hg` flag in the mapping's VmFlags line in
 * /proc/self/smaps, whatever the /dev/shm mount option, so the check does not
 * depend on the runner's mount. Off Linux, or on a kernel built without
 * transparent huge pages (where madvise rejects the advice), it prints SKIP.
 *
 * Like the other shm tests (shm_lock_destroy_test, shm_round_trip_test.sh),
 * it calls bwa_shm_destroy before and after, which drops every index this user
 * has staged on the host, not just its own. Run it where nothing else needs a
 * staged index.
 *
 * Usage: shm_thp_advise_test <indexed-prefix>
 */

#include "bwa_shm.h"

#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <sys/mman.h>
#include <unistd.h>

#define CHECK(cond) do { \
    if (!(cond)) { \
        std::fprintf(stderr, "CHECK failed at %s:%d: %s (errno=%d %s)\n", \
                     __FILE__, __LINE__, #cond, errno, std::strerror(errno)); \
        std::fflush(stderr); \
        bwa_shm_destroy(); \
        std::abort(); \
    } \
} while (0)

#if defined(__linux__)
/* The VmFlags of the /proc/self/smaps mapping that contains `addr`, or "" if
 * there is none. */
static std::string vmflags_of(const void *addr)
{
    FILE *fp = std::fopen("/proc/self/smaps", "r");
    if (fp == NULL) return "";
    const uintptr_t a = (uintptr_t)addr;
    bool in_vma = false;
    std::string flags;
    char line[4096];
    while (std::fgets(line, sizeof(line), fp) != NULL) {
        /* A mapping's header line starts "lo-hi "; its field lines ("Size:",
         * "AnonHugePages:", ...) never parse as two hex numbers joined by '-'. */
        unsigned long lo, hi;
        if (std::sscanf(line, "%lx-%lx ", &lo, &hi) == 2) {
            in_vma = (a >= lo && a < hi);
        } else if (in_vma && std::strncmp(line, "VmFlags:", 8) == 0) {
            flags = line + 8;
            break;
        }
    }
    std::fclose(fp);
    return flags;
}

/* The address of the last whole page of a `len`-byte mapping at `base` (which
 * mmap returns page-aligned). Advice that covered only part of the segment
 * would split the mapping, so checking both ends catches it. */
static const void *last_page_of(const void *base, size_t len)
{
    const size_t page = (size_t)sysconf(_SC_PAGESIZE);
    return (const uint8_t *)base + (len / page - 1) * page;
}

/* The VmFlags of the stager's mapping at its first and last page, captured
 * while it is still mapped. */
static std::string g_stage_flags_first, g_stage_flags_last;
static bool g_stage_hook_ran = false;

static void stage_map_hook(const void *base, size_t len)
{
    g_stage_hook_ran = true;
    g_stage_flags_first = vmflags_of(base);
    if (len >= (size_t)sysconf(_SC_PAGESIZE))
        g_stage_flags_last = vmflags_of(last_page_of(base, len));
}

/* 0 if `flags` (the VmFlags at `what`) carries the hg advice, else prints FAIL. */
static int check_hg(const char *what, const std::string &flags)
{
    if (flags.empty()) {
        std::fprintf(stderr, "FAIL: no /proc/self/smaps entry for the %s\n", what);
        return 1;
    }
    if (flags.find(" hg") == std::string::npos) {
        std::fprintf(stderr, "FAIL: the %s is not advised MADV_HUGEPAGE "
                             "(VmFlags:%s)\n", what, flags.c_str());
        return 1;
    }
    return 0;
}
#endif

int main(int argc, char *argv[])
{
    if (argc != 2) {
        std::fprintf(stderr, "usage: %s <indexed-prefix>\n", argv[0]);
        return 2;
    }
#if !defined(__linux__)
    (void)argv;
    std::printf("SKIP: shm_thp_advise_test (not Linux: no MADV_HUGEPAGE or /proc/self/smaps)\n");
    return 0;
#else
    /* A kernel built without transparent huge pages rejects the advice with
     * EINVAL; probe it on an anonymous page rather than trusting sysfs. */
    const size_t page = (size_t)sysconf(_SC_PAGESIZE);
    void *probe = mmap(NULL, page, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    if (probe == MAP_FAILED) {
        std::fprintf(stderr, "FAIL: mmap of the madvise probe page: %s\n", std::strerror(errno));
        return 1;
    }
    const int probe_rc = madvise(probe, page, MADV_HUGEPAGE);
    const int probe_errno = errno;
    munmap(probe, page);
    if (probe_rc != 0) {
        if (probe_errno == EINVAL) {
            std::printf("SKIP: shm_thp_advise_test (kernel built without transparent huge pages)\n");
            return 0;
        }
        std::fprintf(stderr, "FAIL: madvise(MADV_HUGEPAGE) probe: %s\n", std::strerror(probe_errno));
        return 1;
    }
    const char *prefix = argv[1];
    CHECK(bwa_shm_destroy() == 0);
    bwa_shm_stage_map_hook = stage_map_hook;
    CHECK(bwa_shm_stage(prefix) == 0);
    bwa_shm_stage_map_hook = NULL;
    CHECK(g_stage_hook_ran);

    size_t len = 0;
    uint8_t *base = bwa_shm_attach(prefix, &len);
    CHECK(base != NULL);
    CHECK(len > 0);
    /* bwamem_madv_hugepage advises only whole pages, so a segment smaller than
     * one gets no advice at all. */
    if (len < page) {
        bwa_shm_destroy();
        std::printf("SKIP: shm_thp_advise_test (segment smaller than a page; no advice possible)\n");
        return 0;
    }

    /* The advice covers the page-aligned interior, which starts at the mapping
     * base (mmap returns a page-aligned address) and ends at its last whole
     * page. */
    const int failed = check_hg("staging mapping's first page", g_stage_flags_first)
                     + check_hg("staging mapping's last page", g_stage_flags_last)
                     + check_hg("attached mapping's first page", vmflags_of(base))
                     + check_hg("attached mapping's last page", vmflags_of(last_page_of(base, len)));
    if (failed) {
        bwa_shm_destroy();
        return 1;
    }
    CHECK(bwa_shm_destroy() == 0);
    std::printf("shm_thp_advise_test: OK (staging and attached segment VmFlags have hg at both ends)\n");
    return 0;
#endif
}
