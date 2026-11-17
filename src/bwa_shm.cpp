/*************************************************************************************
                           The MIT License

   BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),
   Copyright (C) 2026  Fulcrum Genomics, contributors of bwa-mem3.

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*****************************************************************************************/

#include "bwa_shm.h"
#include "FMI_search.h"    /* CP_FILENAME_SUFFIX, CP_OCC, CP_SHIFT */
#include "macro.h"         /* SA_COMPX */
#include "bntseq.h"        /* bntseq_t, bns_restore, bns_destroy */
#include "utils.h"         /* err_fread_noeof, err_fclose, xopen */
#include "bwa_madvise.h"   /* bwamem_madv_hugepage */
#include "simd_dispatch.h" /* bwa3_host_cpu_identity */
#include "bwa_thp.h"       /* bwa_smaps_vma_*, bwa_proc_read_stream: the shared smaps reader */

#include <stdio.h>
#include <cstdlib>
#include <cstring>
#include <sys/mman.h>      /* shm_open, mmap, munmap */
#include <sys/ipc.h>
#include <sys/shm.h>       /* shmget, shmat, shmctl: the 1 GB segment */
#include <sys/stat.h>      /* mode constants, stat */
#include <sys/statvfs.h>   /* statvfs — /dev/shm capacity preflight */
#include <fcntl.h>         /* O_* */
#include <unistd.h>        /* close, ftruncate */
#include <getopt.h>        /* getopt_long for `bwa-mem3 shm --meth ...` */
#include <semaphore.h>     /* sem_open, sem_trywait, sem_post, sem_close, sem_unlink */
#include <errno.h>
#include <limits.h>        /* PATH_MAX */
#include <signal.h>        /* kill(pid, 0) liveness probe */
#include <time.h>
#include <atomic>
#include <string>
#include <vector>

/* /bwactl layout: 16-bit n_entries, 16-bit next_write_offset, then packed
 * (int64_t l_mem, char[] basename) entries. The header is 4 bytes total
 * (two uint16_t fields). next_write_offset starts at 4 on a fresh segment. */
#define BWA_SHM_CTL_HEADER_BYTES 4
namespace {
struct ctl_header { uint16_t n_entries, next_write; };
static_assert(sizeof(ctl_header) == BWA_SHM_CTL_HEADER_BYTES, "registry header is 4 bytes");
}

/* The registry header is one 32-bit word, published (release) after the entry
 * bytes, so a lockless reader never counts bytes not yet written. */
static ctl_header ctl_header_load(const void *m)
{
    const uint32_t w = __atomic_load_n((const uint32_t *)m, __ATOMIC_ACQUIRE);
    ctl_header h;
    std::memcpy(&h, &w, sizeof(h));
    return h;
}

static void ctl_header_store(void *m, ctl_header h)
{
    uint32_t w;
    std::memcpy(&w, &h, sizeof(w));
    __atomic_store_n((uint32_t *)m, w, __ATOMIC_RELEASE);
}

/* The /dev/shm free-space probe; a test may substitute its own. */
#ifndef BWA_SHM_STATVFS
#define BWA_SHM_STATVFS statvfs
#endif

/* shmget flags and page size of a 1 GB stage (`shm --huge-pages`); a test may
 * substitute its own. Negative flags: not supported on this platform. */
#if defined(__linux__) && !defined(SHM_HUGE_1GB)
#define SHM_HUGE_1GB (30 << 26)   /* log2(1 GiB) << SHM_HUGE_SHIFT */
#endif
#ifndef BWA_SHM_HUGE_FLAGS
#ifdef __linux__
#define BWA_SHM_HUGE_FLAGS (SHM_HUGETLB | SHM_HUGE_1GB)
#else
#define BWA_SHM_HUGE_FLAGS (-1)
#endif
#endif
#ifndef BWA_SHM_HUGE_PAGE
#define BWA_SHM_HUGE_PAGE (1ull << 30)
#endif
#if defined(__linux__) && !defined(MADV_POPULATE_WRITE)
#define MADV_POPULATE_WRITE 23   /* Linux 5.14 */
#endif

/* Whether the host is an AMD x86 CPU, and whether its kernel removes a System V
 * segment when the last process detaches (kernel.shm_rmid_forced): a 1 GB stage
 * is refused on both (see huge_refused). A test may substitute its own. */
#ifndef BWA_SHM_HOST_IS_AMD
#define BWA_SHM_HOST_IS_AMD host_is_amd
#endif
#ifndef BWA_SHM_RMID_FORCED
#define BWA_SHM_RMID_FORCED rmid_forced
#endif

/* Unlink POSIX index segment `path`; a test may substitute an unlink that fails. */
#ifndef BWA_SHM_UNLINK
#define BWA_SHM_UNLINK(path) shm_unlink(path)
#endif

/* Remove System V segment `id`; a test may substitute a removal that fails. */
#ifndef BWA_SHM_RMID
#define BWA_SHM_RMID(id) shmctl((id), IPC_RMID, NULL)
#endif

/* The per-stage token; a test may substitute its own. */
#ifndef BWA_SHM_STAGE_TOKEN
#define BWA_SHM_STAGE_TOKEN stage_token
#endif

/* How long `shm -d` waits for the stager lock before proceeding without it. */
#ifndef BWA_SHM_DESTROY_LOCK_WAIT_MS
#define BWA_SHM_DESTROY_LOCK_WAIT_MS 5000
#endif

/* How long a process waiting for the stager lock sleeps between tries. */
#define BWA_SHM_LOCK_RETRY_MS 10

/* Test hooks, no-ops by default: the first runs in `shm -d` once it holds the
 * stager lock, before it walks the registry; the second runs each time a try
 * for the stager lock finds it taken. */
#ifndef BWA_SHM_DESTROY_LOCKED_HOOK
#define BWA_SHM_DESTROY_LOCKED_HOOK() ((void)0)
#endif
#ifndef BWA_SHM_LOCK_BUSY_HOOK
#define BWA_SHM_LOCK_BUSY_HOOK() ((void)0)
#endif

/* How long a stage waits for another live stager of the same index. */
#ifndef BWA_SHM_STAGE_WAIT_S
#define BWA_SHM_STAGE_WAIT_S 3600
#endif

/* .bwt.2bit.64 prefix: int64_t reference_seq_len, int64_t count[5]. */
#define BWA_BWT_2BIT_HEADER_BYTES (sizeof(int64_t) * 6)

/* Sanity bound on reference_seq_len read from disk. The bwa-mem2 FM-index
 * has long enforced reference_seq_len <= 0x7fffffffff (about 549 Gbp); see
 * FMI_search::load_index. Anything beyond this is a malformed file. */
#define BWA_SHM_MAX_REF_SEQ_LEN 0x7fffffffffLL

/* Sanity bound on the total packed segment size we'll create. 1 TiB is far
 * larger than any realistic reference genome's index and well below what
 * 64-bit arithmetic can express, so overflow checks against this max are
 * effectively impossible to trip in practice but cheap to enforce. */
#define BWA_SHM_MAX_TOTAL_SIZE (1ull << 40)

/* Return a pointer to the filename component of `hint` — everything after
 * the last '/'. If `hint` has no '/', returns `hint`. */
static const char *prefix_basename(const char *hint)
{
    if (hint == NULL || hint[0] == '\0') return hint;
    const char *p = hint + std::strlen(hint) - 1;
    while (p >= hint && *p != '/') --p;
    return p + 1;
}

/* Build "<a><b>" into `out` (PATH_MAX-bounded), matching the path-building
 * convention used elsewhere in the codebase (see bntseq.cpp, FMI_search.cpp,
 * fastmap.cpp). Aborts via err_fatal if the result would overflow PATH_MAX —
 * every caller passes a PATH_MAX buffer and treats truncation as a hard error. */
static void path_concat2(char out[PATH_MAX], const char *a, const char *b)
{
    int n = std::snprintf(out, PATH_MAX, "%s%s", a, b);
    if (n < 0 || n >= PATH_MAX)
        err_fatal(__func__, "path too long: '%s%s'", a, b);
}

/* The cross-process lock that serializes read-modify-write sequences on
 * /bwactl. POSIX named semaphores work on both Linux and macOS (unlike
 * flock(2) on POSIX shm fds, which macOS rejects with EOPNOTSUPP).
 *
 * The lock is taken by name on every try, not by waiting on one handle:
 * `shm -d` unlinks the semaphore and leaves the unlinked one taken, so a
 * waiter that kept its handle would wait forever while a new stager creates
 * and takes a fresh one under the same name. Reopening on each try moves a
 * waiter to whatever semaphore the name refers to now, so only one is ever
 * taken. */
enum ctl_lock_result { CTL_LOCK_HELD, CTL_LOCK_GONE, CTL_LOCK_TIMEOUT, CTL_LOCK_ERROR };

/* Take the lock, trying every BWA_SHM_LOCK_RETRY_MS. With create, a missing
 * lock is created; without it, a missing lock is CTL_LOCK_GONE (errno ENOENT).
 * wait_ms < 0 waits for as long as it takes; otherwise CTL_LOCK_TIMEOUT once
 * wait_ms have passed. On CTL_LOCK_HELD *out is the held semaphore (release
 * it with ctl_lock_release); otherwise *out is NULL, and CTL_LOCK_ERROR has
 * printed the failing call and left its errno. */
static ctl_lock_result ctl_lock_take(sem_t **out, bool create, int wait_ms)
{
    *out = NULL;
    for (int waited_ms = 0;; waited_ms += BWA_SHM_LOCK_RETRY_MS) {
        sem_t *s = create ? sem_open(BWA_SHM_LOCK_NAME, O_CREAT, 0644, 1)
                          : sem_open(BWA_SHM_LOCK_NAME, 0);
        if (s == SEM_FAILED && !create && errno == ENOENT) return CTL_LOCK_GONE;
        if (s == SEM_FAILED) {
            const int err = errno;
            std::fprintf(stderr, "[E::%s] sem_open(%s) failed: %s\n",
                         __func__, BWA_SHM_LOCK_NAME, std::strerror(err));
            errno = err;
            return CTL_LOCK_ERROR;
        }
        int rc;
        while ((rc = sem_trywait(s)) < 0 && errno == EINTR) {}
        if (rc == 0) { *out = s; return CTL_LOCK_HELD; }
        const int err = errno;
        sem_close(s);
        if (err != EAGAIN) {
            std::fprintf(stderr, "[E::%s] sem_trywait(%s) failed: %s\n",
                         __func__, BWA_SHM_LOCK_NAME, std::strerror(err));
            errno = err;
            return CTL_LOCK_ERROR;
        }
        BWA_SHM_LOCK_BUSY_HOOK();
        if (wait_ms >= 0 && waited_ms >= wait_ms) return CTL_LOCK_TIMEOUT;
        usleep(BWA_SHM_LOCK_RETRY_MS * 1000);
    }
}

/* Take the lock for a stager, waiting for as long as it takes. Returns the
 * held semaphore, or NULL: errno ENOENT when, without create, the lock is
 * gone (`shm -d` removed it), otherwise the failing call's. */
static sem_t *ctl_lock_acquire(bool create = true)
{
    sem_t *s = NULL;
    return ctl_lock_take(&s, create, -1) == CTL_LOCK_HELD ? s : NULL;
}

static void ctl_lock_release(sem_t *s)
{
    if (s == NULL) return;
    sem_post(s);
    sem_close(s);
}

/* Map `/bwactl` read-write (callers hold the lock); `*fd_out` is closed with
 * ctl_close. With create, a missing or short registry is created empty; without
 * it, either fails with ENOENT, so a registry `shm -d` removed stays removed.
 * On failure errno is the failing call's. */
static void *ctl_open_rw(int *fd_out, bool create = true)
{
    int created = 0;
    int fd = shm_open(BWA_SHM_CTL_NAME, O_RDWR | (create ? O_CREAT : 0), 0644);
    if (fd < 0) return NULL;
    /* close() may clobber errno; keep the failing call's. */
    auto fail = [fd](int err) -> void * { close(fd); errno = err; return NULL; };
    struct stat st;
    if (fstat(fd, &st) < 0) return fail(errno);
    if ((size_t)st.st_size < BWA_SHM_CTL_SIZE) {   /* macOS sizes an object once */
        if (!create) return fail(ENOENT);
        if (ftruncate(fd, BWA_SHM_CTL_SIZE) < 0) return fail(errno);
        created = 1;
    }
    void *m = mmap(NULL, BWA_SHM_CTL_SIZE, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    if (m == MAP_FAILED) return fail(errno);
    if (created) {
        std::memset(m, 0, BWA_SHM_CTL_SIZE);
        ctl_header_store(m, ctl_header{ 0, BWA_SHM_CTL_HEADER_BYTES });
    }
    *fd_out = fd;
    return m;
}

static void ctl_close(void *m, int fd)
{
    if (m != NULL) munmap(m, BWA_SHM_CTL_SIZE);
    if (fd >= 0)   close(fd);
}

/* Map an index segment and ask for transparent huge pages on it, as the
 * disk-load path does for its index buffers (FMI_search.cpp). On 4 KB pages
 * the FM-index and SA lookups miss the dTLB about 12x as often, and an
 * shm-attached alignment costs about a third more CPU than a disk load.
 *
 * tmpfs fixes a page's size when it is first written, so the stager must
 * advise its mapping before the fill; the attach advises its own mapping too.
 * The advice takes effect only if /dev/shm is mounted huge=advise (or
 * huge=within_size, which needs no advice); the default huge=never ignores it.
 * Returns MAP_FAILED on failure, like mmap. */
static void *idx_segment_map(int fd, size_t len, int prot)
{
    void *m = mmap(NULL, len, prot, MAP_SHARED, fd, 0);
    if (m != MAP_FAILED) bwamem_madv_hugepage(m, len);
    return m;
}

int bwa_shm_hugepage_warn_wanted(int on_huge, uint64_t seg_len, uint64_t pmd_bytes,
                                 const char *shmem_enabled)
{
    if (on_huge != 0) return 0;                 /* on huge pages, or cannot tell */
    if (seg_len < pmd_bytes) return 0;          /* never PMD-mapped: no huge page fits */
    if (shmem_enabled == NULL) return 0;        /* no shmem THP in the kernel */
    if (std::strstr(shmem_enabled, "[deny]") != NULL) return 0;   /* off for every mount */
    return 1;
}

/* "4 KB", "64 KB", "2 MB", "512 MB", "1 GB": a page size in the message. */
static std::string shm_page_size_text(uint64_t bytes)
{
    static const char *const unit[] = { "B", "KB", "MB", "GB" };
    int u = 0;
    while (u < 3 && bytes >= 1024 && bytes % 1024 == 0) { bytes /= 1024; ++u; }
    return std::to_string((unsigned long long)bytes) + " " + unit[u];
}

int bwa_shm_format_hugepage_warning(char *buf, size_t n, const char *func, const char *name,
                                    uint64_t page_bytes, uint64_t pmd_bytes)
{
    return std::snprintf(buf, n,
        "[W::%s] the shm segment for '%s' is on %s pages, not %s huge pages. "
        "Fix: `sudo mount -o remount,huge=advise /dev/shm` (add huge=advise "
        "to the /dev/shm options in /etc/fstab to keep it across reboots), then restage, since "
        "page size is fixed at staging: `bwa-mem3 shm -d`, then your `bwa-mem3 shm` command "
        "again.\n",
        func, name, shm_page_size_text(page_bytes).c_str(), shm_page_size_text(pmd_bytes).c_str());
}

/* Warn, once per process, when the segment mapped at `base` is not on huge
 * pages. It checks what the kernel did, not the mount options: the smaps entry
 * for the mapping, whose pages the caller has touched (the stager's pack, the
 * attach's magic load). smaps is streamed (bwa_proc_read_stream) and the read
 * stops after that one VMA. Silent at -v < 2, when it cannot tell (smaps
 * unreadable), and when the remount it advises cannot help
 * (bwa_shm_hugepage_warn_wanted: a segment smaller than one PMD page, no shmem
 * THP in the kernel, or shmem_enabled set to deny). The caller never passes a
 * 1 GB (hugetlb) segment, which smaps would report as huge anyway. Linux only;
 * a no-op elsewhere. Never changes the exit status or output. Called with no
 * stager lock held: it reads procfs and sysfs and writes stderr. */
#if defined(__linux__)
static int shm_smaps_vma_feed(void *p, const char *data, size_t n)
{
    bwa_smaps_vma_feed((bwa_smaps_vma_parser_t *)p, data, n);
    return bwa_smaps_vma_done((const bwa_smaps_vma_parser_t *)p);
}

/* The PMD (transparent huge page) size: 2 MB on x86-64 and on arm64 with 4 KB
 * base pages, 512 MB on arm64 with 64 KB ones. */
static uint64_t shm_pmd_bytes(void)
{
    char *text = NULL;
    size_t cap = 0;
    unsigned long long v = 0;
    if (bwa_proc_read_text("/sys/kernel/mm/transparent_hugepage/hpage_pmd_size", &text, &cap, 64) > 0)
        v = std::strtoull(text, NULL, 10);
    std::free(text);
    return v > 0 ? (uint64_t)v : 2ULL << 20;
}
#endif

static std::atomic<bool> g_small_pages_warned{false};

static void shm_warn_if_small_pages(const char *func, const void *base, uint64_t seg_len,
                                    const char *name)
{
#if defined(__linux__)
    if (bwa_verbose < 2 || g_small_pages_warned.load(std::memory_order_relaxed)) return;
    const uint64_t pmd = shm_pmd_bytes();
    if (seg_len < pmd) return;   /* no huge page fits: skip the smaps walk */
    char *thp = NULL;
    size_t cap = 0;
    const bool have_thp =
        bwa_proc_read_text("/sys/kernel/mm/transparent_hugepage/shmem_enabled", &thp, &cap, 4096) >= 0;
    bwa_smaps_vma_parser_t vma;
    bwa_smaps_vma_init(&vma, (unsigned long long)(uintptr_t)base);
    const bool want = bwa_proc_read_stream("/proc/self/smaps", shm_smaps_vma_feed, &vma) == 0 &&
                      bwa_shm_hugepage_warn_wanted(bwa_smaps_vma_finish(&vma), seg_len, pmd,
                                                   have_thp ? thp : NULL);
    std::free(thp);
    if (!want || g_small_pages_warned.exchange(true)) return;
    const uint64_t page = (uint64_t)sysconf(_SC_PAGESIZE);
    const int n = bwa_shm_format_hugepage_warning(NULL, 0, func, name, page, pmd);
    if (n <= 0) return;
    std::string msg((size_t)n + 1, '\0');
    bwa_shm_format_hugepage_warning(&msg[0], msg.size(), func, name, page, pmd);
    std::fputs(msg.c_str(), stderr);
#else
    (void)func; (void)base; (void)seg_len; (void)name;
#endif
}

/* Test-only hook, called after the lock is released and before the pack
 * (shm_thp_advise_test, shm_stage_test). NULL in production. */
void (*bwa_shm_stage_map_hook)(const void *base, size_t len) = NULL;

/* The registry's one parser: cb(l_mem, name-in-m, ctx) per entry of mapped
 * registry `m` until cb returns non-zero; returns cb's last value (0 if none). */
typedef int (*bwa_shm_ctl_cb_t)(int64_t l_mem, const char *name, void *ctx);

static int ctl_walk_map(const void *m, bwa_shm_ctl_cb_t cb, void *ctx)
{
    const uint8_t *base = (const uint8_t *)m;
    ctl_header h = ctl_header_load(base);
    const uint8_t *p   = base + BWA_SHM_CTL_HEADER_BYTES;
    const uint8_t *end = base + BWA_SHM_CTL_SIZE;
    int rc = 0;
    for (uint16_t i = 0; i < h.n_entries; ++i) {
        if (p + sizeof(int64_t) >= end) break;
        int64_t l_mem = 0;
        std::memcpy(&l_mem, p, sizeof(int64_t));
        const char *name = (const char *)(p + sizeof(int64_t));
        size_t name_max  = (size_t)(end - (const uint8_t *)name);
        size_t namelen   = strnlen(name, name_max);
        if (namelen == name_max) break;
        rc = cb(l_mem, name, ctx);
        if (rc != 0) break;
        p = (const uint8_t *)name + namelen + 1;
    }
    return rc;
}

/* ctl_walk_map over the live registry, mapped read-only without taking the
 * lock. A walk racing a locked ctl_remove can misread entries (a shifted name's
 * suffix can read as an entry), so a caller that removes what it finds, as
 * `-d` does, holds the lock itself. Returns 0 if the registry is absent, -1 on
 * open/fstat/mmap error (errno kept). */
static int ctl_walk(bwa_shm_ctl_cb_t cb, void *ctx)
{
    int fd = shm_open(BWA_SHM_CTL_NAME, O_RDONLY, 0);
    if (fd < 0) return errno == ENOENT ? 0 : -1;
    struct stat st;
    if (fstat(fd, &st) < 0) {
        const int err = errno;
        close(fd);
        errno = err;
        return -1;
    }
    /* Short (aborted creation; macOS unlinked at size 0) means empty. */
    if ((size_t)st.st_size < BWA_SHM_CTL_SIZE) {
        close(fd);
        return 0;
    }
    void *m = mmap(NULL, BWA_SHM_CTL_SIZE, PROT_READ, MAP_SHARED, fd, 0);
    close(fd);
    if (m == MAP_FAILED) return -1;
    int rc = ctl_walk_map(m, cb, ctx);
    munmap(m, BWA_SHM_CTL_SIZE);
    return rc;
}

/* The name of `name`'s entry inside the mapped registry `m`, or NULL. */
static char *ctl_find(const void *m, const char *name)
{
    const char *want[2] = { name, NULL };   /* [0] sought, [1] found */
    ctl_walk_map(m, [](int64_t, const char *entry, void *ctx) {
        const char **w = (const char **)ctx;
        return std::strcmp(entry, w[0]) == 0 ? (w[1] = entry, 1) : 0;
    }, want);
    return (char *)want[1];
}

/* Delete `name`'s entry from the mapped registry `m` (held under the lock):
 * move the later entries down and fix n_entries / next_write. */
static void ctl_remove(void *m, const char *name)
{
    char *found = ctl_find(m, name);
    if (found == NULL) return;
    uint8_t *base = (uint8_t *)m;
    uint8_t *entry = (uint8_t *)found - sizeof(int64_t);
    size_t len = sizeof(int64_t) + std::strlen(found) + 1;
    ctl_header h = ctl_header_load(base);
    uint8_t *tail = base + h.next_write;
    if (tail < entry + len) tail = entry + len;   /* inconsistent cursor */
    std::memmove(entry, entry + len, (size_t)(tail - entry - len));
    std::memset(tail - len, 0, len);
    h.n_entries  = (uint16_t)(h.n_entries - 1);
    h.next_write = (uint16_t)(tail - len - base);
    ctl_header_store(base, h);
}

/* A complete index header this release reads: its version, and a checkpoint
 * layout it knows. The one test every reader of a segment header makes. */
static bool shm_header_known(const bwa_shm_header_t *h)
{
    return h->version == BWA_SHM_VERSION && h->cp_layout <= BWA_SHM_CP_LAYOUT_COMPACT;
}

/* Segment states, from the header (see bwa_shm_header_t). */
namespace { enum seg_state { SEG_ABSENT, SEG_COMPLETE, SEG_STAGING, SEG_DEAD, SEG_BAD }; }

/* `h` is the stub of a 1 GB stage (see bwa_shm.h), staging or complete. */
static bool is_huge_stub(const bwa_shm_header_t *h)
{
    return h->reserved[2] != 0 && (h->magic == 0 || h->magic == BWA_SHM_HUGE_MAGIC);
}

/* The System V key of the 1 GB segment staged with `token`. Never IPC_PRIVATE. */
static key_t huge_key(uint64_t token)
{
    const key_t k = (key_t)((token >> 1) & 0x7fffffff);
    return k != IPC_PRIVATE ? k : 1;
}

/* The namespace that /proc/self/ns/`kind` names, as its inode number; 0 where
 * there is none to read (not Linux). The kernel numbers namespaces with 32-bit
 * inodes; a wider one is taken as unreadable. */
static uint64_t self_namespace(const char *kind)
{
    char path[64];
    std::snprintf(path, sizeof path, "/proc/self/ns/%s", kind);
    struct stat st;
    if (stat(path, &st) != 0) return 0;
    if ((uint64_t)st.st_ino > 0xffffffffull) {
        errno = ERANGE;
        return 0;
    }
    return (uint64_t)st.st_ino;
}

/* This process's IPC namespace. System V keys are per IPC namespace, and
 * /dev/shm can be shared across namespaces (a bind mount into a container). */
static uint64_t ipc_namespace(void) { return self_namespace("ipc"); }

/* This process's PID namespace. The kernel reports a segment's creator PID as
 * seen from the asking process's PID namespace. */
static uint64_t pid_namespace(void) { return self_namespace("pid"); }

/* A 1 GB stub's reserved[4]: the stager's IPC namespace in the low 32 bits and
 * its PID namespace in the high 32 bits. 0 records neither (not Linux). */
static uint64_t stub_namespaces(void) { return ipc_namespace() | (pid_namespace() << 32); }
static uint64_t stub_ipc_namespace(const bwa_shm_header_t *stub)
{
    return stub->reserved[4] & 0xffffffffull;
}
static uint64_t stub_pid_namespace(const bwa_shm_header_t *stub)
{
    return stub->reserved[4] >> 32;
}

/* The id of the System V segment that the stage recorded in `stub` created,
 * or one of:
 *   HUGE_ABSENT     no segment has the stub's key;
 *   HUGE_FOREIGN    one does, but this stage did not create it;
 *   HUGE_ERROR      it could not be inspected (errno kept), which includes a
 *                   segment this process may not read: that does not show
 *                   whose it is;
 *   HUGE_ELSEWHERE  the stub was staged in another IPC namespace
 *                   (reserved[4]), so its segment cannot be seen from here
 *                   and nothing about it is concluded;
 *   HUGE_UNPROVEN   the stage's creator PID is needed to identify the
 *                   segment, and the stager was in another PID namespace, so
 *                   the PID the kernel reports here cannot be compared with it.
 * The stage's own segment is the one whose id the stub records. A stager
 * killed between creating the segment and recording its id leaves reserved[3]
 * zero, and then it is the one that the kernel says the stub's PID created, with
 * the stub's size. The stub holds the PID, token and size before the segment
 * exists, so another program's segment on the same key never matches, and only
 * a match is ever removed. */
namespace {
enum { HUGE_ABSENT = -1, HUGE_FOREIGN = -2, HUGE_ERROR = -3, HUGE_ELSEWHERE = -4,
       HUGE_UNPROVEN = -5 };
}
static int huge_segment(const bwa_shm_header_t *stub)
{
    if (stub->reserved[4] != 0 && stub_ipc_namespace(stub) != ipc_namespace())
        return HUGE_ELSEWHERE;
    const int id = shmget(huge_key(stub->reserved[1]), 0, 0);
    if (id < 0) return errno == ENOENT ? HUGE_ABSENT : HUGE_ERROR;
    if (stub->reserved[3] != 0) return stub->reserved[3] == (uint64_t)id + 1 ? id : HUGE_FOREIGN;
    struct shmid_ds ds;
    if (shmctl(id, IPC_STAT, &ds) < 0)   /* EIDRM, EINVAL: removed meanwhile */
        return errno == EIDRM || errno == EINVAL ? HUGE_ABSENT : HUGE_ERROR;
    if (stub_pid_namespace(stub) != 0 && stub_pid_namespace(stub) != pid_namespace())
        return HUGE_UNPROVEN;
    return (uint64_t)ds.shm_cpid == stub->reserved[0] &&
           (uint64_t)ds.shm_segsz == stub->reserved[2] ? id : HUGE_FOREIGN;
}

/* Why huge_segment's result `id` (< 0, errno as it left it) is not a usable
 * segment. */
static const char *huge_trouble(int id)
{
    return id == HUGE_ABSENT ? "it is gone" :
           id == HUGE_FOREIGN ? "it is another program's, not the one the stage created" :
           id == HUGE_ELSEWHERE ? (ipc_namespace() == 0 ?
                                   "this process cannot read its IPC namespace "
                                   "(/proc/self/ns/ipc) to compare with the stage's" :
                                   "it is in another IPC namespace (the stage, every `mem` "
                                   "and `shm -d` must share one)") :
           id == HUGE_UNPROVEN ? "the stage was interrupted before it recorded the segment, "
                                 "and its stager was in another PID namespace, so the "
                                 "segment's creator cannot be checked from here" :
                                 std::strerror(errno);
}

/* `h` is an incomplete 1 GB stub whose stager was in another PID namespace
 * than this process. */
static bool stager_elsewhere(const bwa_shm_header_t *h)
{
    return h->magic == 0 && is_huge_stub(h) && stub_pid_namespace(h) != 0 &&
           stub_pid_namespace(h) != pid_namespace();
}

/* Whether the complete index at `base`, whose header's total_size the caller
 * has checked against what is mapped, has the sections that every shm load
 * reads, within the index, aligned for the types the loads read them as and
 * of the sizes those loads require, so that a
 * stage `shm -l` calls usable is one `mem` can load rather than exit on. The
 * checks are the loaders' own: FMI_search::shm_scalars_ok for the FM-index
 * scalars and section sizes, indexEle::shm_bns_ok for BNS_STRUCT and the BNS
 * sections it sizes (BNS_STRUCT itself must be sizeof(bntseq_t)), and
 * indexEle::shm_names_walk for the names and annotations in BNS_NAMES. PAC and
 * REF_STRING need not be present: a seed-only (--meth) stage packs them
 * empty, and the loads that read them check them. A present PAC must be empty
 * (seed-only) or of the size the PAC load requires. */
static bool index_sections_ok(const uint8_t *base)
{
    uint64_t off = 0, size = 0;
    if (bwa_shm_section_find(base, BWA_SHM_SEC_FMI_SCALARS, &off, &size) != 0 ||
        size != BWA_SHM_FMI_SCALARS_BYTES)
        return false;
    int64_t fmi[3];
    if (!FMI_search::shm_scalars_ok(base + off, fmi, NULL, 0)) return false;
    /* The loads alias CP_OCC, SA_MS (packed 64-bit words), SA_LS and AMBS in
     * place as typed arrays, so each must start at an offset aligned for its
     * element type (the packer's 64-byte layout always is). */
    static const struct { uint32_t kind; int fmi; size_t align; } sized[] = {
        {BWA_SHM_SEC_FMI_CP_OCC, 0, alignof(CP_OCC)}, {BWA_SHM_SEC_FMI_SA_MS, 1, alignof(uint64_t)},
        {BWA_SHM_SEC_FMI_SA_LS, 2, alignof(uint32_t)} };
    for (const auto &r : sized)
        if (bwa_shm_section_find(base, r.kind, &off, &size) != 0 ||
            (int64_t)size != fmi[r.fmi] || off % r.align != 0)
            return false;
    uint64_t off_ambs = 0, sz_ambs = 0, sz_anns = 0, off_names = 0, sz_names = 0;
    if (bwa_shm_section_find(base, BWA_SHM_SEC_BNS_STRUCT, &off, &size) != 0 ||
        size != sizeof(bntseq_t) ||
        bwa_shm_section_find(base, BWA_SHM_SEC_BNS_AMBS, &off_ambs, &sz_ambs) != 0 ||
        off_ambs % alignof(bntamb1_t) != 0 ||
        bwa_shm_section_find(base, BWA_SHM_SEC_BNS_ANNS, NULL, &sz_anns) != 0 ||
        bwa_shm_section_find(base, BWA_SHM_SEC_BNS_NAMES, &off_names, &sz_names) != 0)
        return false;
    bntseq_t bns;
    std::memcpy(&bns, base + off, sizeof(bns));
    if (!indexEle::shm_bns_ok(&bns, sz_ambs, sz_anns, NULL, 0)) return false;
    if (!indexEle::shm_names_walk((const char *)(base + off_names), sz_names, bns.n_seqs,
                                  NULL, NULL, 0))
        return false;
    if (bwa_shm_section_find(base, BWA_SHM_SEC_PAC, NULL, &size) == 0 && size != 0 &&
        (int64_t)size != indexEle::shm_pac_bytes_for(bns.l_pac))
        return false;
    return true;
}

/* Whether the index at `base` (checked as bwa_shm_attach checks it) holds the
 * packed reference a load with PAC reads: a PAC section of the size
 * indexEle::shm_pac_bytes_for gives for its BNS_STRUCT's l_pac. A seed-only
 * (`shm --meth`) stage packs PAC empty, which only a seed load can use. */
static bool index_pac_ok(const uint8_t *base)
{
    uint64_t off = 0, size = 0;
    if (bwa_shm_section_find(base, BWA_SHM_SEC_BNS_STRUCT, &off, &size) != 0 ||
        size != sizeof(bntseq_t))
        return false;
    bntseq_t bns;
    std::memcpy(&bns, base + off, sizeof(bns));
    return bns.l_pac > 0 && bwa_shm_section_find(base, BWA_SHM_SEC_PAC, NULL, &size) == 0 &&
           (int64_t)size == indexEle::shm_pac_bytes_for(bns.l_pac);
}

/* Attach the System V segment of complete 1 GB stub `stub` read-only, after
 * checking that it is usable: the segment that the stub's stage created
 * (huge_segment), of the size the stub records, and holding the complete index
 * that stage packed (a header with the stub's token, and the sections every
 * load reads: index_sections_ok). Returns the attachment,
 * or (void *)-1 with *why saying why not. *id_out (if given) gets
 * huge_segment's result. */
static void *huge_attach(const bwa_shm_header_t *stub, int *id_out, const char **why)
{
    const int id = huge_segment(stub);
    if (id_out != NULL) *id_out = id;
    *why = id < 0 ? huge_trouble(id) : NULL;
    struct shmid_ds ds;
    void *m = (void *)-1;
    if (*why != NULL) return m;
    if (shmctl(id, IPC_STAT, &ds) < 0) *why = std::strerror(errno);
    else if ((uint64_t)ds.shm_segsz != stub->reserved[2])
        *why = "it is not the size its stub records";
    else if ((m = shmat(id, NULL, SHM_RDONLY)) == (void *)-1)
        *why = std::strerror(errno);
    else {
        const uint64_t magic = __atomic_load_n((const uint64_t *)m, __ATOMIC_ACQUIRE);
        bwa_shm_header_t hdr;
        std::memcpy(&hdr, m, sizeof(hdr));
        if (magic != BWA_SHM_MAGIC || hdr.reserved[1] != stub->reserved[1])
            *why = "it is not the index that was staged";
        else if (!shm_header_known(&hdr) || hdr.cp_layout != BWA_SHM_CP_LAYOUT_DEFAULT ||
                 hdr.total_size > stub->reserved[2] || hdr.total_size < sizeof(bwa_shm_header_t))
            *why = "its index header has another version or checkpoint layout, or a size past "
                   "the segment or shorter than the header";
        else if (!index_sections_ok((const uint8_t *)m))
            *why = "its index lacks a section that every load reads, or has one past its "
                   "end, misaligned for its type or of a size the load cannot use, or "
                   "holds names or FM-index scalars a load rejects";
        if (*why != NULL) {
            shmdt(m);
            m = (void *)-1;
        }
    }
    return m;
}

/* Probe `idx_path`; *hdr gets its header (zeroed if none could be read).
 * SEG_BAD from a failed call sets *failed (if given) to its name, errno kept;
 * otherwise *failed is NULL and SEG_BAD means a foreign magic or version, or a
 * complete 1 GB stub whose System V segment is not usable from here, by the
 * checks bwa_shm_attach makes (huge_attach says why). judge_huge=false skips
 * that check (a complete 1 GB stub is then SEG_COMPLETE) for a caller that
 * only needs the header. */
static seg_state segment_state(const char *idx_path, bwa_shm_header_t *hdr,
                               const char **failed = NULL, bool judge_huge = true)
{
    const char *dummy;
    if (failed == NULL) failed = &dummy;
    *failed = NULL;
    std::memset(hdr, 0, sizeof(*hdr));
    int fd = shm_open(idx_path, O_RDONLY, 0);
    if (fd < 0 && errno != ENOENT) *failed = "shm_open";
    if (fd < 0) return errno == ENOENT ? SEG_ABSENT : SEG_BAD;   /* keeps errno */
    struct stat st;
    const bool sized = fstat(fd, &st) == 0;
    if (!sized) *failed = "fstat";
    const bool small = sized && st.st_size < (off_t)sizeof(*hdr);   /* never written */
    void *m = sized && !small ? mmap(NULL, sizeof(*hdr), PROT_READ, MAP_SHARED, fd, 0)
                              : MAP_FAILED;
    const int err = errno;
    close(fd);
    errno = err;
    if (m == MAP_FAILED && sized && !small) *failed = "mmap";
    if (m == MAP_FAILED) return small ? SEG_DEAD : SEG_BAD;
    const uint64_t magic = __atomic_load_n((const uint64_t *)m, __ATOMIC_ACQUIRE);
    std::memcpy(hdr, m, sizeof(*hdr));
    hdr->magic = magic;   /* the value judged below, not a later one */
    munmap(m, sizeof(*hdr));
    errno = 0;
    /* A complete 1 GB stub is as good as its segment: complete only if attach
     * would take it. One whose segment is gone (removed by hand, or by the
     * host) is reclaimed like a dead stage. */
    if (magic == BWA_SHM_HUGE_MAGIC && hdr->version == BWA_SHM_VERSION && is_huge_stub(hdr)) {
        if (!judge_huge) return SEG_COMPLETE;   /* the caller wants the header only */
        int id;
        const char *why;
        void *seg = huge_attach(hdr, &id, &why);
        if (seg != (void *)-1) shmdt(seg);
        errno = 0;
        return seg != (void *)-1 ? SEG_COMPLETE : id == HUGE_ABSENT ? SEG_DEAD : SEG_BAD;
    }
    /* Shorter than its header declares: bwa_shm_attach refuses it, so it is bad. */
    if (magic == BWA_SHM_MAGIC && shm_header_known(hdr))
        return hdr->total_size <= (uint64_t)st.st_size ? SEG_COMPLETE : SEG_BAD;
    if (magic != 0) return SEG_BAD;
    /* A 1 GB stager in another PID namespace cannot be checked for life from
     * here (its PID names another process here, or none): taken as staging,
     * never as dead, so its live stage is not reclaimed. */
    if (stager_elsewhere(hdr)) return SEG_STAGING;
    const uint64_t pid = hdr->reserved[0];
    const bool alive = pid != 0 && pid <= (uint64_t)INT_MAX &&
                       (kill((pid_t)pid, 0) == 0 || errno == EPERM);
    return alive ? SEG_STAGING : SEG_DEAD;
}

/* segment_state of registry entry `name`'s segment. */
static seg_state name_state(const char *name, bwa_shm_header_t *hdr = NULL)
{
    char idx_path[PATH_MAX];
    path_concat2(idx_path, BWA_SHM_IDX_PREFIX, name);
    bwa_shm_header_t own;
    return segment_state(idx_path, hdr != NULL ? hdr : &own);
}

/* For each entry whose segment is still being staged, append to ctx (a
 * std::vector<staging_need>) its size and the bytes it has already allocated
 * (st_blocks; where that is unavailable, as macOS reports 0, none count). A
 * 1 GB stage takes nothing from /dev/shm. */
namespace { struct staging_need { uint64_t l_mem, allocated; }; }
static int staging_bytes_cb(int64_t l_mem, const char *name, void *ctx)
{
    if (l_mem <= 0) return 0;
    char idx_path[PATH_MAX];
    path_concat2(idx_path, BWA_SHM_IDX_PREFIX, name);
    /* Header only: a complete 1 GB stub is never SEG_STAGING however it is
     * judged, so attaching its segment to check its sections would be wasted
     * work under the stager lock. */
    bwa_shm_header_t hdr;
    if (segment_state(idx_path, &hdr, NULL, /*judge_huge=*/false) != SEG_STAGING ||
        is_huge_stub(&hdr))
        return 0;
    uint64_t allocated = 0;
    const int fd = shm_open(idx_path, O_RDONLY, 0);
    if (fd >= 0) {
        struct stat st;
        if (fstat(fd, &st) == 0 && st.st_blocks > 0) allocated = (uint64_t)st.st_blocks * 512u;
        close(fd);
    }
    ((std::vector<staging_need> *)ctx)->push_back(staging_need{ (uint64_t)l_mem, allocated });
    return 0;
}

/* `bytes` rounded up to whole allocation units of `unit` bytes (0: 1). tmpfs
 * allocates whole pages, so a segment's last, partial page takes a page. */
static uint64_t round_to_unit(uint64_t bytes, uint64_t unit)
{
    return unit <= 1 ? bytes : (bytes + unit - 1) / unit * unit;
}

/* Remove registry entry `name`'s segment at `idx_path`. Every removal of a
 * registered segment, by `shm -d`, by the reclaim of an interrupted stage or by
 * a stage's rollback, goes through here (a stage also removes, by id, the
 * System V segment it created itself: stage_guard, packing_shmid_remove). For
 * a 1 GB stub the System V segment goes
 * first: if that fails, the stub stays to name it and the result is -1, so the
 * caller keeps the registry entry and a later `shm -d` can retry. So does a
 * header that cannot be read, since it may be such a stub (unless
 * own_sysv_settled: the stub's segment is then known to be gone), and a
 * shm_unlink that fails other than for a segment already gone. Otherwise 0.
 *
 * own_sysv_settled: the caller is the stage that made this stub, and it knows
 * that no System V segment of its own is left (it never created one, or it has
 * removed it by id), so none is looked up by key: another segment on that key
 * is never this stage's to remove. */
static int segment_remove(const char *idx_path, const char *name,
                          bool own_sysv_settled = false)
{
    bwa_shm_header_t hdr;
    const char *failed;
    segment_state(idx_path, &hdr, &failed, /*judge_huge=*/false);
    if (failed != NULL && !own_sysv_settled) {
        std::fprintf(stderr, "[E::%s] cannot read the header of '%s' (%s(%s): %s), so it may "
                     "name a 1 GB segment; keeping it and its entry\n", __func__, name,
                     failed, idx_path, std::strerror(errno));
        return -1;
    }
    if (is_huge_stub(&hdr) && !own_sysv_settled) {
        const int id = huge_segment(&hdr);
        if (id == HUGE_ERROR || id == HUGE_ELSEWHERE || id == HUGE_UNPROVEN ||
            (id >= 0 && BWA_SHM_RMID(id) < 0 && errno != EIDRM && errno != EINVAL)) {
            const int err = errno;
            std::fprintf(stderr, "[E::%s] cannot remove the 1 GB segment of '%s' (System V "
                         "key 0x%x): %s; keeping its entry\n", __func__, name,
                         (unsigned)huge_key(hdr.reserved[1]),
                         id == HUGE_ELSEWHERE ? "it is in another IPC namespace; run "
                                                "bwa-mem3 shm -d there" :
                         id == HUGE_UNPROVEN ? "the stage was interrupted before it recorded "
                                               "the segment, and its stager was in another PID "
                                               "namespace, so the segment's creator cannot be "
                                               "checked from here" :
                         id >= 0 && err == EPERM ? "Operation not permitted (it is another "
                                                   "user's; run bwa-mem3 shm -d as that user "
                                                   "or as root)" :
                         id == HUGE_ERROR && err == EACCES ?
                                         "Permission denied (this user may not inspect the "
                                         "segment on that key, which may be another user's; "
                                         "run bwa-mem3 shm -d as root)" : std::strerror(err));
            /* The end of a PID namespace frees nothing: only the IPC namespace
             * holds the segment. Without the creator check nothing shows that
             * the segment on that key is the stage's, so none is suggested for
             * removal. */
            if (id == HUGE_UNPROVEN)
                std::fprintf(stderr, "[E::%s] run bwa-mem3 shm -d in the stager's PID namespace "
                             "(the container that staged it), which can check the segment. "
                             "From here nothing shows whether the segment on key 0x%x is the "
                             "stage's, so it and the stub are left as they are\n", __func__,
                             (unsigned)huge_key(hdr.reserved[1]));
            /* A namespace that has ended took its segments with it, and then
             * nothing but the stub is left and no namespace can remove it. A
             * container's exit does not end its namespace while another
             * process, or another container, is still in it. */
            /* A stage interrupted before it recorded its segment's id is
             * identified by the segment's creator PID, which only the
             * stager's PID namespace can compare (HUGE_UNPROVEN), so its
             * removal needs that namespace too. */
            if (id == HUGE_ELSEWHERE && hdr.reserved[3] == 0 && stub_pid_namespace(&hdr) != 0)
                std::fprintf(stderr, "[E::%s] run it in IPC namespace ipc:[%llu] and in the "
                             "stager's PID namespace pid:[%llu] (sudo nsenter -t <pid of a "
                             "process in both> --ipc --pid bwa-mem3 shm -d): the stage was "
                             "interrupted before it recorded the segment's id, so the segment "
                             "is identified by its creator PID, which only that PID namespace "
                             "can check. Run from the IPC namespace alone, shm -d cannot check "
                             "it and keeps it. While it cannot be checked, leave the stub "
                             "/dev/shm%s: deleting it while the segment lives leaves its pages "
                             "pinned with nothing that names them\n", __func__,
                             (unsigned long long)stub_ipc_namespace(&hdr),
                             (unsigned long long)stub_pid_namespace(&hdr), idx_path);
            else if (id == HUGE_ELSEWHERE)
                std::fprintf(stderr, "[E::%s] run it in IPC namespace ipc:[%llu] (sudo nsenter "
                             "-t <pid of a process in it> --ipc bwa-mem3 shm -d), which removes "
                             "the segment and the stub. The "
                             "segment is gone only once that namespace is: no process in it, "
                             "and nothing holding it open or bind-mounted (a container's exit "
                             "is not enough, and lsns may not list a namespace held only by "
                             "an open file). Unless that is certain, leave the stub /dev/shm%s: "
                             "deleting it while the segment lives leaves its pages pinned "
                             "with nothing that names them\n", __func__,
                             (unsigned long long)stub_ipc_namespace(&hdr), idx_path);
            return -1;
        }
        if (id == HUGE_FOREIGN)
            std::fprintf(stderr, "[W::%s] leaving the System V segment with key 0x%x: the "
                         "stage of '%s' did not create it\n", __func__,
                         (unsigned)huge_key(hdr.reserved[1]), name);
    }
    if (BWA_SHM_UNLINK(idx_path) != 0 && errno != ENOENT) {
        std::fprintf(stderr, "[E::%s] cannot remove '%s' (shm_unlink(%s): %s); keeping its "
                     "entry\n", __func__, name, idx_path, std::strerror(errno));
        return -1;
    }
    return 0;
}

/* A per-stage token (header reserved[1]); it only has to differ between
 * stagers, so time, PID and a splitmix finalizer suffice. Never 0. */
static uint64_t stage_token(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    uint64_t x = ((uint64_t)ts.tv_sec * 1000000000ull + (uint64_t)ts.tv_nsec) ^
                 ((uint64_t)getpid() << 40);
    x ^= x >> 33; x *= 0xff51afd7ed558ccdull; x ^= x >> 33;
    return x | 1;
}

extern "C" {

int bwa_shm_compute(const char *prefix, bwa_shm_layout_t *layout, bool bns_only,
                    int target_sa_compx, int n_threads, bool compact)
{
    if (prefix == NULL || layout == NULL) return -1;

    memset(layout, 0, sizeof(*layout));
    layout->densify_threads = (n_threads > 1) ? n_threads : 1;
    layout->compact = compact;
    size_t plen = strlen(prefix);
    if (plen + 1 > sizeof(layout->prefix)) {
        fprintf(stderr, "[E::%s] prefix too long\n", __func__);
        return -1;
    }
    memcpy(layout->prefix, prefix, plen + 1);

    /* 1. Load BNS (.ann + .amb + name/anno strings; opens .pac fp which we
     *    don't need here). */
    bntseq_t *bns = bns_restore(prefix);
    if (bns == NULL) {
        fprintf(stderr, "[E::%s] bns_restore(%s) failed\n", __func__, prefix);
        return -1;
    }
    /* Close the .pac stream that bns_restore left open; pack_into re-opens
     * it via xopen on the .pac path. Zero the field so bns_destroy doesn't
     * double-close. */
    if (bns->fp_pac != NULL) {
        err_fclose(bns->fp_pac);
        bns->fp_pac = NULL;
    }
    layout->bns = bns;

    if (bns->n_seqs < 0 || bns->n_holes < 0 || bns->l_pac <= 0) {
        fprintf(stderr,
            "[E::%s] malformed BNS in %s (n_seqs=%d n_holes=%d l_pac=%lld)\n",
            __func__, prefix, bns->n_seqs, bns->n_holes, (long long)bns->l_pac);
        bwa_shm_layout_free(layout);
        return -1;
    }

    /* 2. Peek FMI scalars + sentinel_index out of <prefix>.bwt.2bit.64. */
    char cp_path[PATH_MAX];
    path_concat2(cp_path, prefix, CP_FILENAME_SUFFIX);
    FILE *cp = xopen(cp_path, "rb");
    err_fread_noeof(&layout->reference_seq_len, sizeof(int64_t), 1, cp);
    err_fread_noeof(layout->count, sizeof(int64_t), 5, cp);

    if (layout->reference_seq_len <= 0 ||
        layout->reference_seq_len > BWA_SHM_MAX_REF_SEQ_LEN) {
        fprintf(stderr,
            "[E::%s] %s: reference_seq_len=%lld out of bounds (1..%lld)\n",
            __func__, cp_path, (long long)layout->reference_seq_len,
            (long long)BWA_SHM_MAX_REF_SEQ_LEN);
        err_fclose(cp);
        bwa_shm_layout_free(layout);
        return -1;
    }
    /* count[i] read from disk is the cumulative occurrence; +1 mirrors the
     * adjustment in FMI_search::load_index. A raw value is in
     * [0, reference_seq_len - 1]: the reference holds one sentinel besides the
     * counted bases. A larger one would let an interval reach row
     * reference_seq_len + 1, whose occurrence block lies past cp_occ. The
     * counts are cumulative, so none may be below the one before it (equal is
     * a base with no occurrences): FMI_search::shm_scalars_ok refuses a
     * decreasing pair, so a stage of one would succeed and be unusable. */
    for (int i = 0; i < 5; ++i) {
        if (layout->count[i] < 0 || layout->count[i] > layout->reference_seq_len - 1) {
            fprintf(stderr,
                "[E::%s] %s: count[%d]=%lld out of bounds (0..%lld)\n",
                __func__, cp_path, i, (long long)layout->count[i],
                (long long)layout->reference_seq_len - 1);
            err_fclose(cp);
            bwa_shm_layout_free(layout);
            return -1;
        }
        if (i > 0 && layout->count[i] < layout->count[i - 1]) {
            /* Both raw: the +1 adjustment is applied after this loop (the compact
             * layout's count form reads the raw counts). */
            fprintf(stderr,
                "[E::%s] %s: count[%d]=%lld < count[%d]=%lld (the counts must not "
                "decrease)\n",
                __func__, cp_path, i, (long long)layout->count[i], i - 1,
                (long long)layout->count[i - 1]);
            err_fclose(cp);
            bwa_shm_layout_free(layout);
            return -1;
        }
    }
    /* Every row but the sentinel's holds a base, so the last raw count is
     * reference_seq_len - 1 (reference_seq_len once +1-adjusted):
     * FMI_search::shm_scalars_ok refuses any other, so a stage of one would
     * succeed and be unusable. */
    if (layout->count[4] != layout->reference_seq_len - 1) {
        fprintf(stderr, "[E::%s] %s: count[4]=%lld != reference_seq_len - 1 (%lld)\n", __func__, cp_path,
                (long long)layout->count[4], (long long)layout->reference_seq_len - 1);
        err_fclose(cp);
        bwa_shm_layout_free(layout);
        return -1;
    }
    /* The compact layout's span size and count form, from the RAW counts as a
     * compact disk load decides them. */
    if (compact) {
        const FMI_search::CpCompactForm form = FMI_search::cp_compact_form(layout->count);
        layout->cp_span_shift = form.span_shift;
        layout->cp_wide       = form.wide ? 1 : 0;
    }
    for (int i = 0; i < 5; ++i) layout->count[i] += 1;

    /* The disk SA rate, by the same detect_sa_compx() the disk loader uses, so
     * the segment is sized for what that loader would read. */
    {
        struct stat st;
        /* A failed fstat must not pass size 0: that reads as a legacy rate. */
        if (fstat(fileno(cp), &st) != 0) {
            fprintf(stderr, "[E::%s] fstat(%s) failed: %s\n",
                    __func__, cp_path, strerror(errno));
            err_fclose(cp);
            bwa_shm_layout_free(layout);
            return -1;
        }
        layout->disk_sa_compx = detect_sa_compx(fileno(cp), (int64_t)st.st_size, layout->reference_seq_len, SA_COMPX);
    }

    /* `shm -u`: a target strictly denser than disk (in [0,6]) sizes the SA
     * sections for it, and pack_into synthesizes the extra samples. Any other
     * target warns and stages the disk rate. */
    layout->sa_compx = layout->disk_sa_compx;
    if (target_sa_compx >= 0) {
        if (compact) {
            /* The densify walk resolves rows against the file's 64-byte records,
             * which a compact segment does not hold; refused for any -u value,
             * a no-op one included, so the rule does not depend on the index. */
            fprintf(stderr,
                "[E::%s] -u (SA densify) cannot be combined with the compact checkpoint layout "
                "(--compact-occ); stage one or the other\n", __func__);
            err_fclose(cp);
            bwa_shm_layout_free(layout);
            return -1;
        }
        if (target_sa_compx < layout->disk_sa_compx && target_sa_compx <= 6) {
            layout->sa_compx = target_sa_compx;
        } else {
            fprintf(stderr,
                "[W::%s] -u %d ignored: not strictly denser than the on-disk SA "
                "rate (shift %lld) or out of [0,6]; staging the disk rate\n",
                __func__, target_sa_compx, (long long)layout->disk_sa_compx);
        }
    }

    int64_t cp_occ_count = (layout->reference_seq_len >> CP_SHIFT) + 1;
    int64_t cp_occ_bytes = cp_occ_count * (int64_t)sizeof(CP_OCC);   /* on disk */
    /* What the segment holds: the file's records, or the compact ones. */
    int64_t cp_sec_bytes = compact ? cp_occ_count * (int64_t)sizeof(CP_OCC_PLANES) : cp_occ_bytes;
    /* The span bases, held only when the counts are span-relative (as a disk load). */
    int64_t cp_base_bytes = compact && layout->cp_wide
        ? ((layout->reference_seq_len >> layout->cp_span_shift) + 1) * (int64_t)(4 * sizeof(int64_t))
        : 0;

    /* SA section sizing uses the STAGED stride (layout->sa_compx), which may be
     * denser than disk when -u densifies. */
    int64_t sa_count     = (layout->reference_seq_len >> layout->sa_compx) + 1;
    /* SA_MS holds the packed high parts, not the disk bytes: the same size the attach
     * check (FMI_search::shm_scalars_ok) requires. */
    int64_t sa_ms_bytes  = FMI_search::sa_hi_bytes_for(layout->reference_seq_len, layout->sa_compx);
    int64_t sa_ls_bytes  = sa_count     * (int64_t)sizeof(uint32_t);

    /* The sentinel_index int64 sits AFTER the on-disk SA arrays in
     * <prefix>.bwt.2bit.64, so its file offset must be computed from the DISK
     * stride (layout->disk_sa_compx) regardless of what we stage. */
    int64_t disk_sa_count = (layout->reference_seq_len >> layout->disk_sa_compx) + 1;
    int64_t sentinel_off  = (int64_t)BWA_BWT_2BIT_HEADER_BYTES + cp_occ_bytes
                          + disk_sa_count * (int64_t)sizeof(int8_t)
                          + disk_sa_count * (int64_t)sizeof(uint32_t);
    if (fseek(cp, (long)sentinel_off, SEEK_SET) != 0) {
        fprintf(stderr, "[E::%s] fseek(%lld) failed in %s\n",
                __func__, (long long)sentinel_off, cp_path);
        err_fclose(cp);
        bwa_shm_layout_free(layout);
        return -1;
    }
    err_fread_noeof(&layout->sentinel_index, sizeof(int64_t), 1, cp);
    err_fclose(cp);

    if (layout->sentinel_index < 0 ||
        layout->sentinel_index >= layout->reference_seq_len) {
        fprintf(stderr,
            "[E::%s] %s: sentinel_index=%lld out of bounds (0..%lld)\n",
            __func__, cp_path, (long long)layout->sentinel_index,
            (long long)layout->reference_seq_len);
        bwa_shm_layout_free(layout);
        return -1;
    }

    /* 3. REF_STRING (.0123) is never staged or stat'd: `mem` fetches bases
     *    from `.pac`. (bns_only, below, separately drops PAC.) */
    layout->ref_string_len = 0;

    /* 4. BNS_NAMES section size: walk anns[]. */
    int64_t bns_names_bytes = 0;
    for (int i = 0; i < bns->n_seqs; ++i) {
        bns_names_bytes += (int64_t)strlen(bns->anns[i].name) + 1
                         + (int64_t)strlen(bns->anns[i].anno) + 1;
    }

    /* 5. Compute section sizes (matching the historical layout). */
    int64_t fmi_scalars_bytes = (int64_t)BWA_SHM_FMI_SCALARS_BYTES;
    int64_t bns_struct_bytes  = (int64_t)sizeof(bntseq_t);
    int64_t bns_ambs_bytes    = (int64_t)bns->n_holes * (int64_t)sizeof(bntamb1_t);
    int64_t bns_anns_bytes    = (int64_t)bns->n_seqs  * (int64_t)sizeof(bntann1_t);
    /* PAC is omitted from a bns_only (D3 --meth seed) segment: `mem --meth`
     * extends against the ORIGINAL pac, never the seed's. Zero-size section. */
    int64_t pac_bytes         = bns_only ? 0 : (bns->l_pac / 4 + 1);

    /* 6. Lay out sections at 64-byte aligned offsets. */
    const uint32_t N_SECTIONS = compact ? 12 : 10;
    uint64_t cursor = sizeof(bwa_shm_header_t)
                    + (uint64_t)N_SECTIONS * sizeof(bwa_shm_section_t);
    cursor = (cursor + 63ull) & ~63ull;

    bwa_shm_section_t *sec = layout->sections;
    auto place = [&](uint32_t kind, int64_t size, uint32_t i) {
        sec[i].kind   = kind;
        sec[i].offset = cursor;
        sec[i].size   = (uint64_t)size;
        cursor = (cursor + (uint64_t)size + 63ull) & ~63ull;
    };

    place(BWA_SHM_SEC_FMI_SCALARS, fmi_scalars_bytes, 0);
    place(compact ? BWA_SHM_SEC_FMI_CP_PLANES : BWA_SHM_SEC_FMI_CP_OCC, cp_sec_bytes, 1);
    place(BWA_SHM_SEC_FMI_SA_MS,   sa_ms_bytes,       2);
    place(BWA_SHM_SEC_FMI_SA_LS,   sa_ls_bytes,       3);
    place(BWA_SHM_SEC_BNS_STRUCT,  bns_struct_bytes,  4);
    place(BWA_SHM_SEC_BNS_AMBS,    bns_ambs_bytes,    5);
    place(BWA_SHM_SEC_BNS_ANNS,    bns_anns_bytes,    6);
    place(BWA_SHM_SEC_BNS_NAMES,   bns_names_bytes,   7);
    place(BWA_SHM_SEC_PAC,         pac_bytes,         8);
    place(BWA_SHM_SEC_REF_STRING,  layout->ref_string_len, 9);
    if (compact) {
        place(BWA_SHM_SEC_FMI_CP_BASE,   cp_base_bytes,                          10);
        place(BWA_SHM_SEC_FMI_CP_LAYOUT, (int64_t)BWA_SHM_FMI_CP_LAYOUT_BYTES,   11);
    }

    if (cursor > BWA_SHM_MAX_TOTAL_SIZE) {
        fprintf(stderr,
            "[E::%s] computed segment size %llu exceeds max %llu\n",
            __func__, (unsigned long long)cursor,
            (unsigned long long)BWA_SHM_MAX_TOTAL_SIZE);
        bwa_shm_layout_free(layout);
        return -1;
    }

    layout->n_sections = N_SECTIONS;
    layout->total_size = cursor;
    return 0;
}

int bwa_shm_pack_into(const bwa_shm_layout_t *layout, uint8_t *dest)
{
    if (layout == NULL || dest == NULL || layout->bns == NULL) return -1;
    if (layout->compact && (uintptr_t)dest % alignof(CP_OCC_PLANES) != 0) {
        fprintf(stderr, "[E::%s] compact destination must be %zu-byte aligned\n",
                __func__, alignof(CP_OCC_PLANES));
        return -1;
    }
    const bntseq_t *bns = (const bntseq_t *)layout->bns;
    const bwa_shm_section_t *sec = layout->sections;

    /* 1. Header: last (step 9). 2. Section table. */
    memcpy(dest + sizeof(bwa_shm_header_t),
           layout->sections,
           (size_t)layout->n_sections * sizeof(bwa_shm_section_t));

    /* 3. FMI_SCALARS (assembled on stack). */
    {
        uint8_t scratch[BWA_SHM_FMI_SCALARS_BYTES];
        memcpy(scratch,        &layout->reference_seq_len, sizeof(int64_t));
        memcpy(scratch + 8,    layout->count,              sizeof(int64_t) * 5);
        memcpy(scratch + 48,   &layout->sentinel_index,    sizeof(int64_t));
        memcpy(scratch + 56,   &layout->sa_compx,          sizeof(int64_t));
        memcpy(dest + sec[0].offset, scratch, sizeof(scratch));
    }

    /* 4. CP_OCC + SA_MS + SA_LS from <prefix>.bwt.2bit.64.
     *
     * No densify (sa_compx == disk_sa_compx): the SA sample arrays on disk are
     * what we stage: cp_occ and sa_ls stream straight into their sections, and
     * the one-byte sa_ms is read in bounded chunks and packed into its section.
     *
     * Densify (sa_compx < disk_sa_compx): cp_occ still streams straight in, but
     * the disk SA arrays are at the coarser disk stride. Read them into heap
     * temporaries, then synthesize the denser sample table directly into the
     * (larger) sec[2] (packed) / sec[3] via FMI_search::densify_sa_into — an LF-walk over
     * the just-staged cp_occ + the disk samples. The result is byte-identical
     * to an index physically built at the denser stride. */
    const bool densify = (layout->sa_compx != layout->disk_sa_compx);

    char cp_path[PATH_MAX];
    path_concat2(cp_path, layout->prefix, CP_FILENAME_SUFFIX);
    FILE *cp = xopen(cp_path, "rb");
    if (fseek(cp, (long)BWA_BWT_2BIT_HEADER_BYTES, SEEK_SET) != 0) {
        fprintf(stderr, "[E::%s] fseek on %s failed\n", __func__, cp_path);
        err_fclose(cp);
        return -1;
    }
    if (!layout->compact) {
        /* CP_OCC first (identical in both modes). */
        err_fread_noeof(dest + sec[1].offset, 1, (size_t)sec[1].size, cp);
    } else {
        /* The compact layout, converted and checked as a compact disk load
         * converts and checks it (FMI_search::cp_compact_convert): the records,
         * then the span bases when the counts are wide; then skip the file's
         * records. */
        const int64_t n_cp = (layout->reference_seq_len >> CP_SHIFT) + 1;
        const FMI_search::CpCompactForm form = { (int)layout->cp_span_shift, layout->cp_wide != 0 };
        int64_t (*base)[4] = NULL;
        int64_t n_spans = 0;
        /* layout->count is +1-adjusted; the conversion checks the raw C[] array */
        int64_t raw_count[5];
        for (int i = 0; i < 5; ++i) raw_count[i] = layout->count[i] - 1;
        const bool sized = (uint64_t)n_cp * sizeof(CP_OCC_PLANES) == sec[1].size;
        const bool ok = sized &&
            FMI_search::cp_compact_convert(fileno(cp), cp_path, layout->reference_seq_len,
                                           raw_count, layout->sentinel_index, form,
                                           (CP_OCC_PLANES *)(dest + sec[1].offset), &base, &n_spans,
                                           layout->densify_threads);
        const bool base_sized = ok && (uint64_t)n_spans * sizeof base[0] == sec[10].size;
        if (base_sized) {
            if (sec[10].size != 0) memcpy(dest + sec[10].offset, base, (size_t)sec[10].size);
            const int64_t lay[2] = { layout->cp_span_shift, layout->cp_wide };
            memcpy(dest + sec[11].offset, lay, sizeof lay);
        }
        free(base);
        if (!base_sized) {
            fprintf(stderr, "[E::%s] %s: %s\n", __func__, cp_path,
                    sized && !ok ? "a checkpoint is not one-hot or out of range, or its counts "
                                   "are not the running totals (corrupt index?)"
                                 : "compact section sizes disagree with the layout");
            err_fclose(cp);
            return -1;
        }
        if (fseeko(cp, (off_t)BWA_BWT_2BIT_HEADER_BYTES + (off_t)n_cp * (off_t)sizeof(CP_OCC), SEEK_SET) != 0) {
            fprintf(stderr, "[E::%s] fseeko on %s failed\n", __func__, cp_path);
            err_fclose(cp);
            return -1;
        }
    }

    const int hi_log2 = FMI_search::sa_hi_log2_for(layout->reference_seq_len);
    /* SA_MS: the packed high parts. Both branches assign every word of the section. */
    uint64_t *hi_dst = (uint64_t *)(dest + sec[2].offset);
    if (!densify) {
        /* The disk's one byte per sample, read and packed into the section a small piece at
         * a time. One thread: -t only sets the -u densify pass's workers, and the pack is a
         * fraction of the stage's time. */
        const int64_t n_s = (layout->reference_seq_len >> layout->sa_compx) + 1;
        const off_t sa_off = ftello(cp);
        const int rc = sa_off < 0 ? FMI_search::SA_HI_READ_OOM
                     : FMI_search::sa_hi_pread_pack(fileno(cp), cp_path, sa_off, n_s, hi_log2,
                                                    hi_dst, /*nthreads=*/1);
        if (rc == 0 && fseeko(cp, sa_off + (off_t)n_s, SEEK_SET) != 0) {
            fprintf(stderr, "[E::%s] fseeko on %s failed\n", __func__, cp_path);
            err_fclose(cp);
            return -1;
        }
        if (rc != 0) {
            if (rc == FMI_search::SA_HI_READ_OOM)
                fprintf(stderr, "[E::%s] could not read the SA samples of %s (out of memory "
                        "or ftello failed)\n", __func__, cp_path);
            else
                fprintf(stderr, "[E::%s] SA sample high part out of range in %s "
                        "(corrupt or mismatched index?)\n", __func__, cp_path);
            err_fclose(cp);
            return -1;
        }
        err_fread_noeof(dest + sec[3].offset, 1, (size_t)sec[3].size, cp);
        err_fclose(cp);
    } else {
        const int64_t disk_count =
            (layout->reference_seq_len >> layout->disk_sa_compx) + 1;
        int8_t   *disk_ms = (int8_t   *)malloc((size_t)disk_count * sizeof(int8_t));
        uint32_t *disk_ls = (uint32_t *)malloc((size_t)disk_count * sizeof(uint32_t));
        if (disk_ms == NULL || disk_ls == NULL) {
            fprintf(stderr,
                "[E::%s] OOM allocating %lld disk SA samples for densify\n",
                __func__, (long long)disk_count);
            free(disk_ms); free(disk_ls);
            err_fclose(cp);
            return -1;
        }
        err_fread_noeof(disk_ms, sizeof(int8_t),   (size_t)disk_count, cp);
        err_fread_noeof(disk_ls, sizeof(uint32_t), (size_t)disk_count, cp);
        err_fclose(cp);

        /* Resolve the denser samples against the staged cp_occ + disk samples.
         * The FMI_search view borrows all pointers (cleared before its dtor). */
        FMI_search fmi(layout->prefix);
        const bool densified = fmi.densify_sa_into(
            (const CP_OCC *)(dest + sec[1].offset),
            layout->count,
            layout->reference_seq_len,
            layout->sentinel_index,
            disk_ms, disk_ls, layout->disk_sa_compx,
            /*dst_ms=*/NULL, hi_dst,
            (uint32_t *)(dest + sec[3].offset),
            layout->sa_compx,
            layout->densify_threads);
        if (!densified) {
            fprintf(stderr, "[E::%s] could not densify the SA samples of %s "
                    "(corrupt or mismatched index?)\n", __func__, cp_path);
            free(disk_ms); free(disk_ls);
            return -1;
        }

        free(disk_ms);
        free(disk_ls);

        fprintf(stderr,
            "[M::%s] densified SA samples %lld -> %lld (stride shift %lld -> %lld, %d thread%s)\n",
            __func__,
            (long long)disk_count,
            (long long)((layout->reference_seq_len >> layout->sa_compx) + 1),
            (long long)layout->disk_sa_compx, (long long)layout->sa_compx,
            layout->densify_threads, layout->densify_threads == 1 ? "" : "s");
    }

    /* 5. BNS_STRUCT, BNS_AMBS, BNS_ANNS — all small. SAM-A3: pos2rid_bucket is
     * a derived, process-private heap allocation that the attach path rebuilds
     * against its own anns[] (bwa_idx_load_ele_from_shm), so stage a NULL
     * rather than publishing this process's heap address to every reader of
     * the segment. Copy through a local so the caller's bns is untouched. */
    {
        bntseq_t bns_image = *bns;
        bns_image.pos2rid_bucket = NULL;
        memcpy(dest + sec[4].offset, &bns_image, sizeof(bns_image));
    }
    memcpy(dest + sec[5].offset, bns->ambs, (size_t)sec[5].size);
    memcpy(dest + sec[6].offset, bns->anns, (size_t)sec[6].size);

    /* 6. BNS_NAMES: concatenate name\0 anno\0 for each seq. The bns_anns
     *    payload's name/anno pointers are not patched here — the unpacker
     *    points them into BNS_NAMES at attach time once the segment base
     *    address is known. */
    {
        uint8_t *dst = dest + sec[7].offset;
        for (int i = 0; i < bns->n_seqs; ++i) {
            size_t n = strlen(bns->anns[i].name) + 1;
            memcpy(dst, bns->anns[i].name, n); dst += n;
            size_t a = strlen(bns->anns[i].anno) + 1;
            memcpy(dst, bns->anns[i].anno, a); dst += a;
        }
    }

    /* 7. PAC: stream <prefix>.pac. Skipped for a bns_only (D3 --meth seed)
     * segment, where the PAC section is zero-length. */
    if (sec[8].size > 0) {
        char pac_path[PATH_MAX];
        path_concat2(pac_path, layout->prefix, ".pac");
        FILE *fp = xopen(pac_path, "rb");
        err_fread_noeof(dest + sec[8].offset, 1, (size_t)sec[8].size, fp);
        err_fclose(fp);
    }

    /* 8. REF_STRING: stream <prefix>.0123. Skipped for a bns_only (D3 --meth
     * seed) segment — the REF_STRING section is zero-length and the seed
     * `.0123` need not exist on disk. */
    if (sec[9].size > 0) {
        char zer_path[PATH_MAX];
        path_concat2(zer_path, layout->prefix, ".0123");
        FILE *fp = xopen(zer_path, "rb");
        err_fread_noeof(dest + sec[9].offset, 1, (size_t)sec[9].size, fp);
        err_fclose(fp);
    }

    /* 9. Header: every field but magic, then magic with a release store, so a
     *    reader that acquire-loads BWA_SHM_MAGIC sees all the bytes above. */
    bwa_shm_header_t hdr;
    memset(&hdr, 0, sizeof(hdr));
    hdr.version     = BWA_SHM_VERSION;
    hdr.n_sections  = (uint16_t)layout->n_sections;
    hdr.cp_layout   = layout->compact ? BWA_SHM_CP_LAYOUT_COMPACT : BWA_SHM_CP_LAYOUT_DEFAULT;
    hdr.total_size  = layout->total_size;
    hdr.reserved[0] = layout->stager_pid;
    hdr.reserved[1] = layout->stage_token;
    memcpy(dest + sizeof(hdr.magic), (const uint8_t *)&hdr + sizeof(hdr.magic),
           sizeof(hdr) - sizeof(hdr.magic));
    __atomic_store_n((uint64_t *)dest, BWA_SHM_MAGIC, __ATOMIC_RELEASE);
    return 0;
}

void bwa_shm_layout_free(bwa_shm_layout_t *layout)
{
    if (layout == NULL) return;
    if (layout->bns != NULL) {
        bns_destroy((bntseq_t *)layout->bns);
        layout->bns = NULL;
    }
}

int bwa_shm_pack_from_disk(const char *prefix, uint8_t **buf_out, uint64_t *len_out)
{
    *buf_out = NULL;
    *len_out = 0;

    bwa_shm_layout_t layout;
    if (bwa_shm_compute(prefix, &layout) != 0) return -1;

    /* Zeroed, as a new segment is: bwa_shm_pack_into leaves the alignment
     * padding between sections unwritten, and a staged index is compared with
     * this buffer byte for byte. */
    uint8_t *buf = (uint8_t *)calloc(1, layout.total_size);
    if (buf == NULL) {
        bwa_shm_layout_free(&layout);
        return -1;
    }
    if (bwa_shm_pack_into(&layout, buf) != 0) {
        free(buf);
        bwa_shm_layout_free(&layout);
        return -1;
    }

    *buf_out = buf;
    *len_out = layout.total_size;
    bwa_shm_layout_free(&layout);
    return 0;
}

int bwa_shm_test(const char *prefix)
{
    if (prefix == NULL || prefix[0] == '\0') return -1;
    return ctl_walk([](int64_t, const char *entry, void *want) {
        return std::strcmp(entry, (const char *)want) == 0 ? 1 : 0;
    }, (void *)prefix_basename(prefix));
}

int bwa_shm_test_complete(const char *prefix)
{
    const int registered = bwa_shm_test(prefix);
    if (registered != 1) return registered;
    return name_state(prefix_basename(prefix)) == SEG_COMPLETE ? 1 : 0;
}

static int list_cb(int64_t l_mem, const char *name, void * /*ctx*/)
{
    bwa_shm_header_t hdr;
    const seg_state st = name_state(name, &hdr);
    const bool compact = st == SEG_COMPLETE && hdr.cp_layout == BWA_SHM_CP_LAYOUT_COMPACT;
    std::printf("%s\t%lld%s\n", name, (long long)l_mem,
                st == SEG_STAGING || st == SEG_DEAD ? "\tincomplete" :
                st == SEG_ABSENT ? "\tmissing" : st == SEG_BAD ? "\tunusable" :
                compact ? "\tcompact" : "");
    return 0;
}

int bwa_shm_list(void)
{
    int rc = ctl_walk(list_cb, NULL);
    return rc < 0 ? rc : 0;
}

/* ctx: destroy_result. */
namespace {
struct destroy_result {
    bool locked = false;   /* `shm -d` holds the stager lock */
    int kept = 0;
    int live_kept = 0;     /* kept for a stager that may hold the lock */
    std::vector<std::string> removed;
};
}
static int destroy_cb(int64_t /*l_mem*/, const char *name, void *ctx)
{
    destroy_result *r = (destroy_result *)ctx;
    char idx_path[PATH_MAX];
    path_concat2(idx_path, BWA_SHM_IDX_PREFIX, name);
    /* Without the lock, a live 1 GB stager that has not recorded its segment
     * may still be about to create it (it does so under the lock): removing
     * the stub now would leave that segment with nothing to name it if the
     * stager were then killed. Keep it until the stager is done or dead. One
     * in another PID namespace counts as live (segment_state). */
    bwa_shm_header_t hdr;
    const seg_state st = r->locked ? SEG_ABSENT : segment_state(idx_path, &hdr);
    if (st == SEG_STAGING && is_huge_stub(&hdr) && hdr.reserved[3] == 0) {
        std::fprintf(stderr, "[E::bwa_shm_destroy] '%s' is being staged on 1 GB pages by "
                     "PID %llu%s, which may still create its segment (the stager lock could "
                     "not be taken); keeping its entry. Run bwa-mem3 shm -d again once that "
                     "stage has finished or its stager has exited\n", name,
                     (unsigned long long)hdr.reserved[0],
                     stager_elsewhere(&hdr) ? " of another PID namespace" : "");
        ++r->kept;
        ++r->live_kept;
        return 0;
    }
    if (segment_remove(idx_path, name) != 0) ++r->kept;
    else r->removed.push_back(name);
    return 0;
}

int bwa_shm_destroy(void)
{
    /* Hold the lock across the walk, so no stager compacts the registry under
     * it and shifts a name's suffix into place as an entry. A missing lock is
     * created and taken: with no lock, no stager can be inside a registry
     * update. Registry updates take milliseconds; a lock held past the wait
     * was left by a stager killed inside one (no SEM_UNDO), so proceed
     * without it, and unlinking it below is the recovery. Each try reopens
     * the lock by name (ctl_lock_take), so a lock that a concurrent `shm -d`
     * unlinked and left taken is not waited on until the timeout. */
    sem_t *lock = NULL;
    const ctl_lock_result got = ctl_lock_take(&lock, /*create=*/true,
                                              BWA_SHM_DESTROY_LOCK_WAIT_MS);
    const bool held = got == CTL_LOCK_HELD;
    if (got == CTL_LOCK_TIMEOUT)
        std::fprintf(stderr, "[W::%s] stager lock %s still held after %d ms; "
                     "removing everything without it\n", __func__, BWA_SHM_LOCK_NAME,
                     (int)BWA_SHM_DESTROY_LOCK_WAIT_MS);
    else if (!held)
        std::fprintf(stderr, "[W::%s] could not take the stager lock %s (see above); "
                     "removing everything without it\n", __func__, BWA_SHM_LOCK_NAME);
    /* A segment that could not be removed keeps the registry and its own
     * entry, which still leads the next `shm -d` to it; the entries of the
     * segments that were removed go, but only under the lock: without it a
     * live stager may be rewriting the registry, so they stay, listed as
     * `missing`, and the next stage drops them. A lock that could not be
     * taken is unlinked (the recovery above), unless a 1 GB stager that may
     * still hold it was kept: a new lock would let the next `shm -d` take it
     * and remove that stager's stub. */
    if (held) BWA_SHM_DESTROY_LOCKED_HOOK();
    destroy_result r;
    r.locked = held;
    const int rc = ctl_walk(destroy_cb, &r);
    const int kept = r.kept;
    bool unlinked = false;
    if (rc >= 0 && kept == 0) {
        shm_unlink(BWA_SHM_CTL_NAME);
        /* ENOENT: the name is already gone, so this handle is orphaned too. */
        unlinked = sem_unlink(BWA_SHM_LOCK_NAME) == 0 || errno == ENOENT;
    } else if (rc >= 0) {
        int fd = -1;
        void *ctl = held ? ctl_open_rw(&fd, /*create=*/false) : NULL;
        if (ctl != NULL) {
            for (const std::string &name : r.removed) ctl_remove(ctl, name.c_str());
            ctl_close(ctl, fd);
        }
        if (!held && r.live_kept == 0) sem_unlink(BWA_SHM_LOCK_NAME);
        if (!held && r.live_kept != 0)
            std::fprintf(stderr, "[E::%s] left the stager lock %s in place, held by the stager "
                         "named above: new stages wait until it is released. If that "
                         "stager is dead (for one in another PID namespace, check there), "
                         "remove the lock (rm /dev/shm/sem.%s on Linux) and run "
                         "bwa-mem3 shm -d again\n", __func__, BWA_SHM_LOCK_NAME,
                         &BWA_SHM_LOCK_NAME[1]);
        std::fprintf(stderr, "[E::%s] %d segment%s could not be removed; kept the registry "
                     "and %s, so bwa-mem3 shm -d can be run again\n", __func__,
                     kept, kept == 1 ? "" : "s", kept == 1 ? "its entry" : "their entries");
    }
    /* An unlinked lock is left taken: a stager still waiting on it would
     * otherwise take it while a new stager takes the semaphore now under
     * that name. Stagers retry by name (ctl_lock_take), so they move on. */
    if (held) {
        if (!unlinked) sem_post(lock);
        sem_close(lock);
    }
    return rc < 0 ? rc : kept != 0 ? -1 : 0;
}

/* The System V segment of a 1 GB stage that this process is packing, or -1.
 * The pack's file readers exit the process on an I/O error (err_fatal), which
 * skips the stage_guard, so the segment is also removed at exit: a stage that
 * fails this way leaves an incomplete stub and no pinned pages. */
static volatile sig_atomic_t g_packing_shmid = -1;   /* read by packing_sigbus */
static void packing_shmid_remove(void)
{
    if (g_packing_shmid >= 0) BWA_SHM_RMID(g_packing_shmid);
}

/* SIGBUS while packing a 1 GB stage is a page the kernel would not fault in
 * (a hugetlb cgroup limit or NUMA policy, on a kernel that cannot prefault,
 * or a hardware error). Its default action skips atexit, so the handler
 * removes the segment, says why, and dies of the signal as before: the stage
 * leaves an incomplete stub and no pinned pages, as for an I/O error. */
static void packing_sigbus(int sig)
{
    /* shmctl is not on POSIX's async-signal-safe list, but it is a bare
     * system call with no user-space state, like those that are. */
    const int id = g_packing_shmid;
    const bool removed = id >= 0 && BWA_SHM_RMID(id) == 0;
    static const char msg[] = "[E::bwa_shm_stage] SIGBUS while packing on 1 GB pages: a "
        "page could not be faulted in (a hugetlb cgroup limit, a NUMA memory policy, or a "
        "hardware memory error). Nothing was staged";
    static const char gone[] = "; removed the segment\n";
    static const char kept[] = "; the segment could not be removed: bwa-mem3 shm -d "
        "removes it\n";
    if (write(STDERR_FILENO, msg, sizeof(msg) - 1) < 0 ||
        write(STDERR_FILENO, removed ? gone : kept, (removed ? sizeof(gone) : sizeof(kept)) - 1) < 0)
        { /* nothing to do */ }
    signal(sig, SIG_DFL);
    raise(sig);
}

/* Installs packing_sigbus for its lifetime when `on`, restoring the previous
 * action after. */
namespace {
struct packing_sigbus_guard {
    bool on;
    struct sigaction old;
    explicit packing_sigbus_guard(bool enable) : on(enable) {
        if (!on) return;
        struct sigaction sa;
        std::memset(&sa, 0, sizeof(sa));
        sa.sa_handler = packing_sigbus;
        sigemptyset(&sa.sa_mask);
        on = sigaction(SIGBUS, &sa, &old) == 0;
    }
    ~packing_sigbus_guard() { if (on) sigaction(SIGBUS, &old, NULL); }
};
} // namespace

/* Everything bwa_shm_stage holds, released on every return path; while
 * `created` is set, destruction also rolls the stage back. */
namespace {
struct stage_guard {
    const char *name, *idx_path;
    bwa_shm_layout_t *layout;
    sem_t *lock = NULL;
    void  *ctl = NULL;
    int    ctl_fd = -1;
    void  *idx = MAP_FAILED;
    size_t idx_len = 0;                  /* length of the `idx` mapping */
    int    idx_fd = -1;
    /* A 1 GB stage: `stub` maps the POSIX segment and `idx` is the attachment
     * of System V segment `shmid`. */
    void  *stub = MAP_FAILED;
    int    shmid = -1;
    bool   created = false, released = false;

    stage_guard(const char *n, const char *p, bwa_shm_layout_t *l)
        : name(n), idx_path(p), layout(l) {}

    /* create=false (after release) recreates nothing that `shm -d` removed;
     * a missing lock or registry then fails with errno ENOENT. */
    bool lock_registry(bool create = true) {
        if ((lock = ctl_lock_acquire(create)) == NULL) return false;
        if ((ctl = ctl_open_rw(&ctl_fd, create)) != NULL) return true;
        const int err = errno;
        if (create || err != ENOENT)
            std::fprintf(stderr, "[E::bwa_shm_stage] failed to open control segment %s: %s\n",
                         BWA_SHM_CTL_NAME, std::strerror(err));
        unlock();
        errno = err;
        return false;
    }
    void unlock() {
        ctl_close(ctl, ctl_fd);
        ctl_lock_release(lock);
        ctl = NULL, ctl_fd = -1, lock = NULL;
    }
    /* The segment at idx_path is the one this stage created. */
    bool still_ours() {
        bwa_shm_header_t hdr;
        segment_state(idx_path, &hdr);
        return !released || hdr.reserved[1] == layout->stage_token;
    }
    /* Roll back an unfinished stage, touching the segment only if it still
     * carries our token. Without the lock, leave the entry: `-l` shows it
     * `missing` and the next stage drops it. */
    ~stage_guard() {
        g_packing_shmid = -1;
        /* This process created `shmid`, so it is removed whatever has become
         * of the stub (`shm -d` may have taken it); and before the detach,
         * while this attachment still keeps the id this segment's. This
         * process knows what became of its own System V segment (none was
         * created, or it is removed here), so the stub's rollback looks up
         * none by key, unless that removal failed and the stub must keep
         * naming the segment: a segment found on the key would be another
         * stage's (two stages of one process can draw the same key). */
        bool sysv_settled = true;
        if (created && shmid >= 0 && BWA_SHM_RMID(shmid) < 0 && errno != EIDRM &&
            errno != EINVAL)
            sysv_settled = false;
        if (idx != MAP_FAILED && shmid >= 0) shmdt(idx);
        else if (idx != MAP_FAILED) munmap(idx, idx_len);
        if (stub != MAP_FAILED) munmap(stub, sizeof(bwa_shm_header_t));
        if (idx_fd >= 0) close(idx_fd);
        if (created) {
            const bool locked = lock != NULL || lock_registry(/*create=*/false);
            if (still_ours() && segment_remove(idx_path, name, sysv_settled) == 0 && locked)
                ctl_remove(ctl, name);
        }
        unlock();
        bwa_shm_layout_free(layout);
    }
};
} // namespace

/* errno `err` is EHWPOISON: a page with a memory error the hardware reported. */
static bool is_hwpoison(int err)
{
#ifdef EHWPOISON
    return err == EHWPOISON;
#else
    (void)err;
    return false;
#endif
}

static bool host_is_amd(void)
{
    int is_amd = 0;
    uint32_t family = 0, model = 0;
    bwa3_host_cpu_identity(&is_amd, &family, &model);
    return is_amd != 0;
}

/* kernel.shm_rmid_forced: 1 if set, 0 if read as 0, -1 if it cannot be read
 * (no such file, masked, or not a number; errno says why, EINVAL for the
 * last). */
static int rmid_forced(void)
{
    FILE *f = fopen("/proc/sys/kernel/shm_rmid_forced", "r");
    if (f == NULL) return -1;
    int v = 0;
    const bool got = fscanf(f, "%d", &v) == 1;
    const int err = ferror(f) != 0 ? errno : EINVAL;
    fclose(f);
    if (!got) {
        errno = err;
        return -1;
    }
    return v != 0 ? 1 : 0;
}

/* Whether this host cannot take a 1 GB stage at all, saying why if so. */
static bool huge_refused(void)
{
    if (BWA_SHM_HUGE_FLAGS < 0) {
        std::fprintf(stderr, "[E::bwa_shm_stage] --huge-pages needs Linux\n");
        return true;
    }
    /* 1 GB staging is a gain on Intel and on Graviton4. On AMD it measured no
     * gain (Zen 4: 0.4 to 0.6% slower than /dev/shm on 2 MB pages), so it is
     * refused there rather than left to make every later run slower. The
     * environment variable is for measuring another AMD CPU. */
    const char *amd_anyway = getenv("BWAMEM3_SHM_HUGE_PAGES_AMD");
    if (BWA_SHM_HOST_IS_AMD() && !(amd_anyway != NULL && std::strcmp(amd_anyway, "1") == 0)) {
        std::fprintf(stderr, "[E::bwa_shm_stage] --huge-pages is refused on AMD x86 hosts: "
                     "staging on 1 GB pages measured no gain there (AMD Zen 4: 0.4 to 0.6%% "
                     "slower than /dev/shm on 2 MB pages). Nothing was staged; stage without "
                     "--huge-pages. To measure it on this host anyway, set "
                     "BWAMEM3_SHM_HUGE_PAGES_AMD=1\n");
        return true;
    }
    const int forced = BWA_SHM_RMID_FORCED();
    if (forced != 0) {
        if (forced < 0)
            std::fprintf(stderr, "[E::bwa_shm_stage] --huge-pages cannot read "
                         "kernel.shm_rmid_forced (/proc/sys/kernel/shm_rmid_forced): %s. If it "
                         "were 1, the kernel would remove the 1 GB segment as soon as this "
                         "command exits. Nothing was staged; stage without --huge-pages, or "
                         "where that setting can be read\n", std::strerror(errno));
        else
            std::fprintf(stderr, "[E::bwa_shm_stage] --huge-pages cannot work with "
                         "kernel.shm_rmid_forced=1: the kernel would remove the 1 GB segment "
                         "as soon as this command exits. Nothing was staged; stage without "
                         "--huge-pages, or set it to 0 (sudo sysctl kernel.shm_rmid_forced=0)\n");
        return true;
    }
#ifdef __linux__
    /* The stub records the stager's IPC and PID namespaces: the first keeps a
     * process in another IPC namespace from taking the segment for gone, the
     * second one in another PID namespace from misreading its creator. */
    for (const char *kind : {"ipc", "pid"}) {
        if (self_namespace(kind) == 0) {
            std::fprintf(stderr, "[E::bwa_shm_stage] --huge-pages needs /proc/self/ns/%s, to "
                         "record the %s namespace of the 1 GB stage: %s. Nothing was staged; "
                         "stage without --huge-pages\n", kind, kind[0] == 'i' ? "IPC" : "PID",
                         errno == ERANGE ? "its inode number does not fit in 32 bits"
                                         : std::strerror(errno));
            return true;
        }
    }
#endif
    return false;
}

/* `prefix`, whose segment header is `hdr`, is already staged: say so and
 * return 0. A 1 GB stage of an index that is staged on other pages fails (-1)
 * instead, since it did not get what it asked for; the index is left as it is. */
static int report_already_staged(const char *prefix, const bwa_shm_header_t *hdr, bool huge)
{
    const bool on_huge = is_huge_stub(hdr);
    if (huge && !on_huge) {
        std::fprintf(stderr, "[E::bwa_shm_stage] index '%s' is already in shm, and not on "
                     "1 GB pages; it is left as it is. To move it, run bwa-mem3 shm -d (which "
                     "removes every staged index) and stage it again\n", prefix);
        return -1;
    }
    std::fprintf(stderr, "[M::bwa_shm_stage] index '%s' is already in shm%s\n", prefix,
                 on_huge && !huge ? " (on 1 GB pages)" : "");
    return 0;
}

/* Whether System V segment `id` is gone or marked for removal (`ipcrm` while
 * this process is attached leaves it so until the detach). */
static bool huge_removed(int id)
{
    struct shmid_ds ds;
    if (shmctl(id, IPC_STAT, &ds) < 0) return true;
#ifdef SHM_DEST
    if ((ds.shm_perm.mode & SHM_DEST) != 0) return true;
#endif
    return false;
}
/* A complete segment already staged for `name`: 0 if it holds the requested
 * checkpoint layout, -1 with an error if it holds the other one (one segment
 * per index; `shm -d` first to restage it). */
static int staged_layout_check(const char *func, const char *prefix,
                               const bwa_shm_header_t *hdr, bool compact)
{
    const bool staged_compact = hdr->cp_layout == BWA_SHM_CP_LAYOUT_COMPACT;
    if (staged_compact == compact) return 0;
    std::fprintf(stderr, "[E::%s] index '%s' is already in shm with the %s checkpoint layout; "
                 "run bwa-mem3 shm -d first (it drops every staged index) to stage it with the "
                 "%s one\n", func, prefix,
                 staged_compact ? "compact (--compact-occ)" : "default",
                 compact ? "compact (--compact-occ)" : "default");
    return -1;
}

int bwa_shm_stage(const char *prefix, bool bns_only, int target_sa_compx,
                  int n_threads, bool huge, bool compact)
{
    if (prefix == NULL || prefix[0] == '\0') {
        std::fprintf(stderr, "[E::%s] empty prefix\n", __func__);
        return -1;
    }
    const char *name = prefix_basename(prefix);
    if (name[0] == '\0') {
        std::fprintf(stderr, "[E::%s] prefix has no basename component\n", __func__);
        return -1;
    }
    char idx_path[PATH_MAX];
    path_concat2(idx_path, BWA_SHM_IDX_PREFIX, name);

    /* A 1 GB stage holds the default checkpoint layout only: its stub and
     * attach checks know that segment version alone. */
    if (huge && compact) {
        std::fprintf(stderr, "[E::%s] --huge-pages cannot be combined with the compact checkpoint "
                     "layout (--compact-occ); nothing was staged\n", __func__);
        return -1;
    }
    if (compact && target_sa_compx >= 0) {
        std::fprintf(stderr, "[E::%s] -u (SA densify) cannot be combined with the compact checkpoint "
                     "layout (--compact-occ); stage one or the other\n", __func__);
        return -1;
    }
    /* A host that cannot take a 1 GB stage refuses one whatever is staged. */
    if (huge && huge_refused()) return -1;

    /* 0. Already staged: skip the disk reads and the lock. Step 2 decides. */
    bwa_shm_header_t hdr;
    if (bwa_shm_test(prefix) == 1 && segment_state(idx_path, &hdr) == SEG_COMPLETE)
        return staged_layout_check(__func__, prefix, &hdr, compact) != 0
               ? -1 : report_already_staged(prefix, &hdr, huge);

    /* 1. Compute the layout (loads BNS; peeks scalars from the FMI file). */
    bwa_shm_layout_t layout;
    if (bwa_shm_compute(prefix, &layout, bns_only, target_sa_compx, n_threads, compact) != 0) {
        std::fprintf(stderr, "[E::%s] failed to compute layout for '%s'\n",
                     __func__, prefix);
        return -1;
    }
    layout.stager_pid  = (uint64_t)getpid();
    layout.stage_token = BWA_SHM_STAGE_TOKEN();
    /* A 1 GB stage puts the index in a System V segment of whole 1 GiB pages
     * and only a header-sized stub in the POSIX segment. (Rounded here, not by
     * bwamem_hugetlb_round_up, because a test substitutes the page size.) */
    const uint64_t huge_bytes = !huge ? 0 :
        (layout.total_size + BWA_SHM_HUGE_PAGE - 1) / BWA_SHM_HUGE_PAGE * BWA_SHM_HUGE_PAGE;
    const size_t posix_len = huge ? sizeof(bwa_shm_header_t) : (size_t)layout.total_size;
    stage_guard g(name, idx_path, &layout);

    /* 2. Under the lock, decide from the entry and the segment header, then
     *    create and register a staging segment. Behind a live stager of the
     *    same index, wait without the lock: the lock must never be held long. */
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    const time_t deadline = now.tv_sec + BWA_SHM_STAGE_WAIT_S;
    seg_state state;
    unsigned long long pid;
    const char *failed;
    for (bool announced = false;; announced = true) {
        if (!g.lock_registry()) return -1;
        state = segment_state(idx_path, &hdr, &failed);
        pid = (unsigned long long)hdr.reserved[0];
        if (state != SEG_STAGING) break;
        g.unlock();
        /* Waiting cannot end well for a stager whose life cannot be checked
         * from here: it would never be seen to exit, nor reclaimed if dead. */
        if (stager_elsewhere(&hdr)) {
            std::fprintf(stderr, "[E::%s] '%s' is being staged on 1 GB pages, or was left by "
                         "an interrupted stage, by PID %llu of another PID namespace, which "
                         "cannot be checked from here; not reclaiming it. Stage it, or run "
                         "bwa-mem3 shm -d, in that PID namespace\n",
                         __func__, name, pid);
            return -1;
        }
        if (!announced)
            std::fprintf(stderr, "[M::%s] '%s' is being staged by PID %llu; waiting\n",
                         __func__, name, pid);
        clock_gettime(CLOCK_MONOTONIC, &now);
        if (now.tv_sec >= deadline) {
            std::fprintf(stderr, "[E::%s] gave up after %d s waiting for PID %llu to stage "
                         "'%s'; if PID %llu is not a running bwa-mem3 shm, run bwa-mem3 shm -d\n",
                         __func__, (int)BWA_SHM_STAGE_WAIT_S, pid, name, pid);
            return -1;
        }
        sleep(1);
    }
    const int probe_errno = errno;
    if (state == SEG_COMPLETE && ctl_find(g.ctl, name) != NULL) {
        g.unlock();   /* the header is copied: report without holding the stager lock */
        return staged_layout_check(__func__, prefix, &hdr, compact) != 0
               ? -1 : report_already_staged(prefix, &hdr, huge);
    }
    if (state == SEG_COMPLETE || state == SEG_BAD) {   /* never adopt or replace it */
        int complete_id = -1;
        /* `shm -d` walks the registry, so it reaches only a registered stub. */
        const bool registered = ctl_find(g.ctl, name) != NULL;
        static const char unreached[] = "; it has no registry entry, so bwa-mem3 shm -d does "
                                        "not reach it";
        if (failed != NULL)
            std::fprintf(stderr, "[E::%s] %s(%s): %s\n", __func__, failed, idx_path,
                         std::strerror(probe_errno));
        else if (hdr.magic == BWA_SHM_HUGE_MAGIC && hdr.version != BWA_SHM_VERSION)
            std::fprintf(stderr, "[E::%s] '%s' has a 1 GB stub %s of version %u (expected "
                         "%u), staged by another bwa-mem3. Not replacing it%s%s\n", __func__,
                         name, idx_path, hdr.version, BWA_SHM_VERSION,
                         registered ? "; bwa-mem3 shm -d removes it" : unreached,
                         registered ? "" : ". Do not delete the stub while the segment it "
                                           "names may exist, or that segment's pages stay "
                                           "pinned with nothing that names them");
        else if (state == SEG_BAD && hdr.magic == BWA_SHM_HUGE_MAGIC && is_huge_stub(&hdr)) {
            int id;
            const char *why;
            void *seg = huge_attach(&hdr, &id, &why);
            const bool usable = seg != (void *)-1;
            if (usable) {   /* usable after all: it changed meanwhile */
                shmdt(seg);
                why = "it changed while it was being checked";
            }
            std::fprintf(stderr, "[E::%s] '%s' has a 1 GB stub %s, but its 1 GB segment "
                         "(System V key 0x%x) cannot be used: %s. Not replacing it%s",
                         __func__, name, idx_path, (unsigned)huge_key(hdr.reserved[1]),
                         why,
                         !registered ? unreached :
                         id == HUGE_FOREIGN ? "; bwa-mem3 shm -d removes the stub and leaves "
                                              "that segment" :
                         id == HUGE_ELSEWHERE ? "; run bwa-mem3 shm -d in the IPC namespace "
                                                "that staged it" :
                         id >= 0 && !usable ? "; bwa-mem3 shm -d removes it and its segment"
                                            : "");
            /* Unregistered: what to do by hand, as for a complete stub
             * without an entry below. Only the stage's own segment (id >= 0)
             * may be removed; a stub whose segment may still exist stays. */
            if (registered)
                std::fprintf(stderr, "\n");
            else if (id >= 0 && !usable)
                std::fprintf(stderr, ". It names System V segment id %d: after checking that "
                             "no stager is active, that `ipcs -m` still lists id %d with key "
                             "0x%x, and that it has no attached process (nattch 0 in `ipcs -m "
                             "-i %d`), remove it with `ipcrm -m %d`, then delete the stub (rm "
                             "/dev/shm%s)\n", id, id, (unsigned)huge_key(hdr.reserved[1]), id,
                             id, idx_path);
            else if (id == HUGE_FOREIGN)
                std::fprintf(stderr, ". The segment on its key is not the stage's, so delete "
                             "only the stub (rm /dev/shm%s)\n", idx_path);
            else if (id == HUGE_ABSENT)
                std::fprintf(stderr, ". Its segment is gone, so delete the stub (rm "
                             "/dev/shm%s)\n", idx_path);
            else
                std::fprintf(stderr, ". Do not delete the stub while its segment may exist, "
                             "or that segment's pages stay pinned with nothing that names "
                             "them\n");
        }
        else if (state == SEG_COMPLETE && is_huge_stub(&hdr) &&
                 (complete_id = huge_segment(&hdr)) >= 0) {
            /* Complete: huge_segment finds the stage's own segment (unless it
             * went meanwhile: then the generic message below). */
            const int id = complete_id;
            std::fprintf(stderr, "[E::%s] 1 GB stub %s exists without a registry entry; not "
                         "replacing it. It names System V segment id %d (key 0x%x, %llu "
                         "bytes). After checking that no stager is active, that `ipcs -m` "
                         "still lists id %d with that key and size, and that it has no "
                         "attached process (nattch 0 in `ipcs -m -i %d`), remove it with "
                         "`ipcrm -m %d`, delete the stub (rm /dev/shm%s), then retry; if "
                         "the key or size differ, the segment is not this stage's: delete "
                         "only the stub\n", __func__,
                         idx_path, id, (unsigned)huge_key(hdr.reserved[1]),
                         (unsigned long long)hdr.reserved[2], id, id, id, idx_path);
        }
        else if (is_huge_stub(&hdr))   /* its segment went between the two looks */
            std::fprintf(stderr, "[E::%s] 1 GB stub %s exists without a registry entry, and "
                         "its segment changed while it was being checked; not replacing it. "
                         "Retry; do not delete the stub while its segment may exist, or its "
                         "pages stay pinned with nothing that names them\n", __func__,
                         idx_path);
        else if (state == SEG_BAD && hdr.magic == BWA_SHM_MAGIC && !shm_header_known(&hdr) &&
                 registered)
            /* Ours, registered, but a segment version or checkpoint layout this release
             * does not read: staged by a different release, which `shm -d` removes like any
             * registered segment. A segment this release reads that is bad (shorter than
             * its header declares) gets the generic message below. Printed under the
             * stager lock like the other refusals in this block; the guard releases it on
             * the return just after. */
            std::fprintf(stderr, "[E::%s] '%s' is already in shm as segment version %u "
                         "(checkpoint layout %u), staged by a different bwa-mem3 release (this "
                         "one stages version %u); run bwa-mem3 shm -d first, then stage it "
                         "again\n", __func__, name, hdr.version, (unsigned)hdr.cp_layout,
                         BWA_SHM_VERSION);
        else
            std::fprintf(stderr, "[E::%s] segment %s exists %s; not replacing it. Remove "
                         "it with shm_unlink(\"%s\") after checking that no stager is "
                         "active, then retry\n", __func__, idx_path,
                         state == SEG_BAD ? "with a foreign magic or version, or "
                                            "shorter than its header declares"
                                          : "without a registry entry", idx_path);
        return -1;
    }
    if (state == SEG_DEAD) {
        if (hdr.magic == BWA_SHM_HUGE_MAGIC)
            std::fprintf(stderr, "[M::%s] reclaiming '%s': its 1 GB segment (System V key "
                         "0x%x) is gone\n", __func__, name, (unsigned)huge_key(hdr.reserved[1]));
        else
            std::fprintf(stderr, "[M::%s] reclaiming '%s' left by an interrupted stage (PID "
                         "%llu)\n", __func__, name, pid);
        if (segment_remove(idx_path, name) != 0) return -1;
    }
    ctl_remove(g.ctl, name);   /* a stale entry, if any */

    /* Bound-check the new entry. next_write is a uint16_t and BWA_SHM_CTL_SIZE
     * (65536) does not fit one, so reject the exact-fill case too: the cast
     * back must not wrap to 0 and let the next stage overwrite the header. */
    ctl_header h = ctl_header_load(g.ctl);
    const size_t name_bytes  = std::strlen(name) + 1;
    const size_t entry_bytes = sizeof(int64_t) + name_bytes;
    if ((size_t)h.next_write + entry_bytes >= BWA_SHM_CTL_SIZE) {
        std::fprintf(stderr, "[E::%s] control segment full (cannot stage '%s')\n",
                     __func__, name);
        return -1;
    }
    g.idx_fd = shm_open(idx_path, O_CREAT | O_RDWR | O_EXCL, 0644);
    if (g.idx_fd < 0) {
        std::fprintf(stderr, "[E::%s] shm_open(%s) failed: %s\n",
                     __func__, idx_path, std::strerror(errno));
        return -1;
    }
    g.created = true;

    /* /dev/shm preflight: tmpfs ftruncate succeeds past the limit and the pack
     * would later die of SIGBUS, so refuse up front with the fix. Stages in
     * progress allocate pages as they write, so their unallocated remainder is
     * reserved too. Sampled before free space, so pages written in between count
     * twice rather than not at all (conservative). f_bavail is in f_frsize units,
     * tmpfs's page, and every segment, this one and those in progress, takes
     * whole ones: a stub admitted into the bytes left in another stage's last
     * page would take that page from it. A statvfs failure is not fatal. A 1 GB
     * stage needs room for its stub only. */
    std::vector<staging_need> staging_needs;
    ctl_walk_map(g.ctl, staging_bytes_cb, &staging_needs);
    struct statvfs svfs;
    if (BWA_SHM_STATVFS("/dev/shm", &svfs) == 0) {
        const uint64_t unit = svfs.f_frsize != 0 ? (uint64_t)svfs.f_frsize
                                                 : (uint64_t)svfs.f_bsize;
        uint64_t staging = 0;
        for (const staging_need &s : staging_needs) {
            const uint64_t need = round_to_unit(s.l_mem, unit);
            if (s.allocated < need) staging += need - s.allocated;
        }
        uint64_t avail = unit * (uint64_t)svfs.f_bavail;
        avail = avail > staging ? avail - staging : 0;
        if (unit != 0 && round_to_unit((uint64_t)posix_len, unit) > avail) {
            if (huge)
                std::fprintf(stderr, "[E::%s] /dev/shm has no room for the %llu-byte stub of "
                             "'%s' (%llu bytes free)\n", __func__,
                             (unsigned long long)posix_len, name, (unsigned long long)avail);
            else std::fprintf(stderr,
                "[E::%s] /dev/shm has %llu bytes free but staging '%s' "
                "needs %llu bytes. Remount with larger size, e.g. "
                "`sudo mount -o remount,size=%lluG /dev/shm`.\n",
                __func__, (unsigned long long)avail, name,
                (unsigned long long)layout.total_size,
                (unsigned long long)((layout.total_size + (1ULL<<30) - 1) / (1ULL<<30) + 2));
            return -1;
        }
    }
    if (ftruncate(g.idx_fd, (off_t)posix_len) < 0) {
        std::fprintf(stderr, "[E::%s] ftruncate(%s, %llu) failed: %s\n",
                     __func__, idx_path,
                     (unsigned long long)posix_len, std::strerror(errno));
        return -1;
    }
    g.idx = idx_segment_map(g.idx_fd, posix_len, PROT_READ | PROT_WRITE);
    g.idx_len = posix_len;
    if (g.idx == MAP_FAILED) {
        std::fprintf(stderr, "[E::%s] mmap(%s) failed: %s\n",
                     __func__, idx_path, std::strerror(errno));
        return -1;
    }
    bwa_shm_header_t *first = (bwa_shm_header_t *)g.idx;   /* magic stays 0 */
    first->reserved[0] = layout.stager_pid;
    first->reserved[1] = layout.stage_token;
    if (huge) {
        first->version     = BWA_SHM_VERSION;
        first->total_size  = layout.total_size;
        first->reserved[2] = huge_bytes;
        first->reserved[4] = stub_namespaces();
    }

    uint8_t *entry = (uint8_t *)g.ctl + h.next_write;
    const int64_t l_mem = (int64_t)layout.total_size;
    std::memcpy(entry,                   &l_mem, sizeof(int64_t));
    std::memcpy(entry + sizeof(int64_t), name,   name_bytes);
    h.n_entries  = (uint16_t)(h.n_entries + 1);
    h.next_write = (uint16_t)(h.next_write + entry_bytes);
    ctl_header_store(g.ctl, h);   /* after the entry bytes */
    if (huge) {
        /* Created only now, with the stub registered and naming it, so a stage
         * killed at any point leaves nothing that `shm -d` or a restage cannot
         * find; and under the lock, so `shm -d` cannot remove the stub in
         * between. Creating reserves pages that the pool already holds and
         * writes none, which takes milliseconds. */
        g.stub = g.idx;
        g.idx = MAP_FAILED;
        g.shmid = shmget(huge_key(layout.stage_token), (size_t)huge_bytes,
                         IPC_CREAT | IPC_EXCL | BWA_SHM_HUGE_FLAGS | 0644);
        /* Recorded at once: until it is, only the kernel's record of the
         * creating PID ties the segment to this stub (huge_segment). */
        if (g.shmid >= 0) first->reserved[3] = (uint64_t)g.shmid + 1;
        void *m = g.shmid < 0 ? (void *)-1 : shmat(g.shmid, NULL, 0);
        if (m == (void *)-1) {
            const int err = errno;
            const bool made = g.shmid >= 0;
            std::fprintf(stderr, "[E::%s] cannot %s a System V segment of %llu x 1 GB pages "
                         "for '%s': %s (%s). Nothing was staged; rerun without --huge-pages "
                         "to stage in /dev/shm\n", __func__,
                         made ? "attach" : "create",
                         (unsigned long long)(huge_bytes / BWA_SHM_HUGE_PAGE), name,
                         std::strerror(err),
                         made ? "the address space or mapping limit of this process" :
                         err == ENOMEM ? "the 1 GB pool has too few unreserved pages "
                                         "(free_hugepages less resv_hugepages in /sys/kernel/"
                                         "mm/hugepages/hugepages-1048576kB; nr_hugepages there "
                                         "sets the pool), or a hugetlb cgroup limit is "
                                         "reached" :
                         err == EPERM ? "it needs CAP_IPC_LOCK or membership of the group "
                                        "named by vm.hugetlb_shm_group (sudo sysctl "
                                        "vm.hugetlb_shm_group=<gid>)" :
                         err == ENODEV || err == ENOSYS ?
                                         "this kernel or CPU has no 1 GB huge page size" :
                         err == EINVAL ? "the segment exceeds kernel.shmmax, or this host has "
                                         "no 1 GB page size" :
                         err == ENOSPC ? "kernel.shmall or kernel.shmmni is exhausted" :
                         err == EEXIST ? "another segment already has this stage's System V "
                                         "key; stage again, each stage draws a new key" :
                                         "shmget");
            return -1;
        }
        g.idx = m;
        g.idx_len = (size_t)huge_bytes;
        static bool at_exit = false;
        if (!at_exit) at_exit = atexit(packing_shmid_remove) == 0;
        g_packing_shmid = g.shmid;
    }
    g.unlock();
    g.released = true;

#ifdef __linux__
    /* Fault the 1 GB pages in before the pack writes them. The pool reservation
     * made by shmget does not cover a hugetlb cgroup limit or a NUMA policy
     * that excludes the nodes holding the free pages; a page refused then is
     * SIGBUS in the pack, and here it is an error; so is a page that the
     * hardware reports poisoned (EHWPOISON). Linux before 5.14 has no
     * MADV_POPULATE_WRITE (EINVAL), and the pack faults the pages itself; a
     * page refused then is SIGBUS, which removes the segment (packing_sigbus). */
    if (huge && madvise(g.idx, g.idx_len, MADV_POPULATE_WRITE) != 0 &&
        (errno == ENOMEM || errno == EFAULT || is_hwpoison(errno))) {
        const int err = errno;
        std::fprintf(stderr, "[E::%s] cannot fault in the %llu x 1 GB pages of '%s': %s (%s). "
                     "Nothing was staged\n", __func__,
                     (unsigned long long)(huge_bytes / BWA_SHM_HUGE_PAGE), name,
                     std::strerror(err),
                     is_hwpoison(err) ? "a page of the pool has a hardware memory error" :
                                        "a hugetlb cgroup limit, or a NUMA memory policy that "
                                        "excludes the nodes holding the pool's free pages");
        return -1;
    }
#endif

    /* 3. Pack without the lock; the magic, written last, marks it complete.
     *    A 1 GB page refused at its first write (a kernel without
     *    MADV_POPULATE_WRITE) is SIGBUS: remove the segment as it kills us. */
    const packing_sigbus_guard sigbus(huge);
    if (bwa_shm_stage_map_hook) bwa_shm_stage_map_hook(g.idx, g.idx_len);
    if (bwa_shm_pack_into(&layout, (uint8_t *)g.idx) != 0) {
        std::fprintf(stderr, "[E::%s] pack_into failed for %s\n", __func__, prefix);
        return -1;
    }

    /* 4. `shm -d`, or a stager that misjudged us dead, may have removed our
     *    entry or segment meanwhile; check without recreating anything. */
    const bool locked = g.lock_registry(/*create=*/false);
    if (!locked && errno != ENOENT) {
        std::fprintf(stderr, "[E::%s] cannot retake the stager lock to finish '%s': %s\n",
                     __func__, name, std::strerror(errno));
        return -1;
    }
    if (!locked || ctl_find(g.ctl, name) == NULL || !g.still_ours() ||
        (huge && huge_removed(g.shmid))) {
        std::fprintf(stderr, "[E::%s] '%s' was removed during staging (by bwa-mem3 shm -d, "
                     "another stager, or ipcrm)\n", __func__, name);
        return -1;
    }
    g.created = false;
    g_packing_shmid = -1;
    if (huge) {   /* the index is complete; now so is the stub */
        __atomic_store_n(&((bwa_shm_header_t *)g.stub)->magic, BWA_SHM_HUGE_MAGIC,
                         __ATOMIC_RELEASE);
        std::fprintf(stderr, "[M::%s] staged '%s' (%llu bytes) on %llu x 1 GB pages "
                     "(System V segment id %d, key 0x%x)\n", __func__, name,
                     (unsigned long long)layout.total_size,
                     (unsigned long long)(huge_bytes / BWA_SHM_HUGE_PAGE), g.shmid,
                     (unsigned)huge_key(layout.stage_token));
        return 0;
    }
    /* The stage is complete. Warn with no lock held (the pack has written every
     * page); a 1 GB stage is on huge pages by construction and returned above. */
    g.unlock();
    shm_warn_if_small_pages(__func__, g.idx, layout.total_size, name);
    std::fprintf(stderr, "[M::%s] staged '%s' (%llu bytes) as %s\n",
                 __func__, name, (unsigned long long)layout.total_size, idx_path);
    return 0;
}

uint8_t *bwa_shm_attach(const char *prefix, size_t *len_out, size_t *map_len_out,
                        bool need_pac, int cp_layout)
{
    if (len_out != NULL) *len_out = 0;
    if (map_len_out != NULL) *map_len_out = 0;
    if (prefix == NULL || prefix[0] == '\0') return NULL;

    /* Cheap pre-check: skip the open if there's no registry entry. */
    if (bwa_shm_test(prefix) != 1) return NULL;

    const char *name = prefix_basename(prefix);
    char idx_path[PATH_MAX];
    path_concat2(idx_path, BWA_SHM_IDX_PREFIX, name);

    int fd = shm_open(idx_path, O_RDONLY, 0);
    if (fd < 0) {
        std::fprintf(stderr, "[E::%s] shm_open(%s) failed: %s\n",
                     __func__, idx_path, std::strerror(errno));
        return NULL;
    }

    struct stat st;
    if (fstat(fd, &st) < 0) {
        std::fprintf(stderr, "[E::%s] fstat(%s) failed: %s\n",
                     __func__, idx_path, std::strerror(errno));
        close(fd);
        return NULL;
    }
    if (st.st_size < (off_t)sizeof(bwa_shm_header_t)) {
        std::fprintf(stderr, "[E::%s] segment %s too small (%lld bytes)\n",
                     __func__, idx_path, (long long)st.st_size);
        close(fd);
        return NULL;
    }
    size_t len = (size_t)st.st_size;

    void *m = idx_segment_map(fd, len, PROT_READ);
    close(fd);
    if (m == MAP_FAILED) {
        std::fprintf(stderr, "[E::%s] mmap(%s, %llu) failed: %s\n",
                     __func__, idx_path, (unsigned long long)len, std::strerror(errno));
        return NULL;
    }

    /* Validate the on-segment header. Acquire-load the magic first: it is
     * published last, so once it reads BWA_SHM_MAGIC the rest is in place. A
     * zero magic is a stage still in progress or one that was interrupted. */
    uint64_t magic = __atomic_load_n((const uint64_t *)m, __ATOMIC_ACQUIRE);
    bwa_shm_header_t hdr;
    std::memcpy(&hdr, m, sizeof(hdr));
    hdr.magic = magic;
    /* A complete 1 GB stub: the index is the System V segment that its stage
     * created, attached by huge_attach after the checks that segment_state
     * makes too (the size the stub records, the length mapped here, and a
     * complete header that carries the stub's token). */
    const bool from_huge = magic == BWA_SHM_HUGE_MAGIC && is_huge_stub(&hdr);
    if (from_huge) {
        munmap(m, len);
        if (hdr.version != BWA_SHM_VERSION) {
            std::fprintf(stderr, "[E::%s] unsupported 1 GB stub version %u in %s (expected "
                         "%u): it was staged by a different bwa-mem3 release and is not attached; "
                         "run bwa-mem3 shm -d (which removes all staged indices), then stage it "
                         "again with this release\n", __func__, hdr.version, idx_path, BWA_SHM_VERSION);
            return NULL;
        }
        const bwa_shm_header_t stub = hdr;
        int id;
        const char *why;
        m = huge_attach(&stub, &id, &why);
        len = (size_t)stub.reserved[2];
        if (m == (void *)-1) {
            std::fprintf(stderr, "[M::%s] the 1 GB segment of '%s' (System V key 0x%x) "
                         "cannot be used: %s; loading from disk%s\n", __func__, name,
                         (unsigned)huge_key(stub.reserved[1]), why,
                         id == HUGE_ABSENT ? " (stage it again to use shared memory)" : "");
            return NULL;
        }
        magic = __atomic_load_n((const uint64_t *)m, __ATOMIC_ACQUIRE);
        std::memcpy(&hdr, m, sizeof(hdr));
    }
    if (magic != BWA_SHM_MAGIC) {
        if (magic == 0)
            std::fprintf(stderr, "[M::%s] '%s' is not fully staged (being staged or "
                         "interrupted); loading from disk\n", __func__, name);
        else
            std::fprintf(stderr, "[E::%s] bad magic in %s\n", __func__, idx_path);
        bwa_shm_detach((uint8_t *)m, len);
        return NULL;
    }
    /* For a 1 GB index huge_attach has already checked the version and size
     * below; they guard the POSIX path. */
    if (!shm_header_known(&hdr)) {
        std::fprintf(stderr, "[E::%s] unsupported segment version %u (checkpoint layout %u) in %s "
                     "(expected version %u): it was staged by a different bwa-mem3 release and "
                     "is not attached; run bwa-mem3 shm -d, then stage it again with this "
                     "release\n", __func__, hdr.version, (unsigned)hdr.cp_layout, idx_path,
                     BWA_SHM_VERSION);
        bwa_shm_detach((uint8_t *)m, len);
        return NULL;
    }
    /* total_size is the exact data written; st_size may be page-rounded up.
     * Require total_size <= segment size to tolerate ftruncate's rounding. */
    if (hdr.total_size > (uint64_t)len) {
        std::fprintf(stderr, "[E::%s] header total_size=%llu exceeds segment size %llu in %s\n",
                     __func__, (unsigned long long)hdr.total_size,
                     (unsigned long long)len, idx_path);
        bwa_shm_detach((uint8_t *)m, len);
        return NULL;
    }
    /* A caller bound to one checkpoint layout (the seed facade, the tests)
     * cannot use the other: decline it, so the load reads the disk. */
    if (cp_layout >= 0 && (unsigned)cp_layout != hdr.cp_layout) {
        std::fprintf(stderr, "[W::%s] '%s' is staged in shared memory with the %s checkpoint "
                     "layout and this load uses the %s one; not attaching it, so a private copy "
                     "loads from disk\n", __func__, name,
                     hdr.cp_layout == BWA_SHM_CP_LAYOUT_COMPACT ? "compact" : "default",
                     cp_layout == (int)BWA_SHM_CP_LAYOUT_COMPACT ? "compact" : "default");
        bwa_shm_detach((uint8_t *)m, len);
        return NULL;
    }
    /* A load that reads the packed reference cannot use a seed-only stage,
     * and would exit on it: decline it here, so the load reads the disk. */
    if (need_pac && !index_pac_ok((const uint8_t *)m)) {
        std::fprintf(stderr, "[M::%s] '%s' is staged without a usable packed reference, "
                     "which this load needs (a seed-only `shm --meth` stage, for one); "
                     "loading from disk\n", __func__, name);
        bwa_shm_detach((uint8_t *)m, len);
        return NULL;
    }

    /* The magic load above touched the first page. A 1 GB segment is on huge
     * pages by construction. */
    if (!from_huge) shm_warn_if_small_pages(__func__, m, hdr.total_size, name);

    if (len_out != NULL) *len_out = (size_t)hdr.total_size;
    if (map_len_out != NULL) *map_len_out = len;
    return (uint8_t *)m;
}

void bwa_shm_detach(uint8_t *base, size_t map_len)
{
    /* A 1 GB index is a System V attachment, any other a mapping; shmdt fails
     * (EINVAL) on an address that is not an attachment's. */
    if (base != NULL && shmdt(base) != 0) munmap(base, map_len);
}

int bwa_shm_section_find(const uint8_t *base, uint32_t kind,
                         uint64_t *offset_out, uint64_t *size_out)
{
    if (base == NULL) return -1;
    const bwa_shm_header_t *hdr = (const bwa_shm_header_t *)base;
    if (hdr->magic   != BWA_SHM_MAGIC)   return -1;
    if (!shm_header_known(hdr))          return -1;

    /* Bound-check the section table itself before walking it. A corrupt
     * n_sections could otherwise drive us past the mapped region. */
    if (hdr->total_size < sizeof(bwa_shm_header_t)) return -1;
    const uint64_t table_bytes =
        (uint64_t)hdr->n_sections * (uint64_t)sizeof(bwa_shm_section_t);
    if (table_bytes > hdr->total_size - sizeof(bwa_shm_header_t)) return -1;

    const bwa_shm_section_t *table =
        (const bwa_shm_section_t *)(base + sizeof(bwa_shm_header_t));
    for (uint32_t i = 0; i < hdr->n_sections; ++i) {
        if (table[i].kind == kind) {
            /* Reject sections whose extent escapes the segment. */
            if (table[i].size > hdr->total_size ||
                table[i].offset > hdr->total_size - table[i].size) {
                return -1;
            }
            if (offset_out != NULL) *offset_out = table[i].offset;
            if (size_out   != NULL) *size_out   = table[i].size;
            return 0;
        }
    }
    return -1;
}

int bwa_shm_segment_compact(const uint8_t *base)
{
    return base != NULL && ((const bwa_shm_header_t *)base)->cp_layout == BWA_SHM_CP_LAYOUT_COMPACT;
}

/* Resolve the on-disk prefix for `bwa-mem3 shm --meth <idxbase>` so that
 * the registry basename matches what `bwa-mem3 mem --meth <idxbase>` will
 * later look up. Mirrors fastmap.cpp's D3 seed-suffix resolution: append
 * `.meth` (the converted seed FM-index) unless `prefix` already ends with it. */
static int resolve_meth_prefix(const char *prefix, char out[PATH_MAX])
{
    static const char SUFFIX[] = ".meth";
    static const size_t SUFLEN = sizeof(SUFFIX) - 1;
    size_t plen = std::strlen(prefix);
    int already = (plen >= SUFLEN) &&
                  (std::strcmp(prefix + plen - SUFLEN, SUFFIX) == 0);
    if (already) {
        if (plen + 1 > (size_t)PATH_MAX) {
            std::fprintf(stderr, "[E::%s] prefix too long\n", __func__);
            return -1;
        }
        strcpy(out, prefix);
    } else {
        if (plen + SUFLEN + 1 > (size_t)PATH_MAX) {
            std::fprintf(stderr, "[E::%s] prefix too long for --meth\n", __func__);
            return -1;
        }
        path_concat2(out, prefix, SUFFIX);
    }
    return 0;
}

/* If the plain FM-index `<prefix>.bwt.2bit.64` is absent but the seed index
 * `<prefix>.meth.bwt.2bit.64` is present, the user likely meant `--meth`. Print
 * a hint and return 1. Returns 0 when no hint applies (let the regular stage
 * path produce its own error). Probes `.bwt.2bit.64`, not `.0123`: the unpacked
 * `.0123` is no longer built by default (mem pac-fetches from `.pac`), so a
 * complete plain index has no `.0123` — only `.bwt.2bit.64` is a reliable
 * "index present" sentinel. */
static int hint_missing_meth_flag(const char *prefix)
{
    char fmi_path[PATH_MAX];
    path_concat2(fmi_path, prefix, ".bwt.2bit.64");
    struct stat st;
    if (stat(fmi_path, &st) == 0) return 0;        /* plain index present */
    if (errno != ENOENT) return 0;                 /* I/O error: defer */

    char meth_fmi[PATH_MAX];
    path_concat2(meth_fmi, prefix, ".meth.bwt.2bit.64");
    if (stat(meth_fmi, &st) != 0) return 0;        /* no meth seed index either */

    std::fprintf(stderr,
        "[E::main_shm] no FMI at '%s.bwt.2bit.64' but a meth seed index is present at "
        "'%s.meth.*'\n"
        "             did you mean `bwa-mem3 shm --meth %s`?\n",
        prefix, prefix, prefix);
    return 1;
}

static void print_shm_usage(void)
{
    std::fprintf(stderr,
        "\nUsage: bwa-mem3 shm [-d|-l|--help] [--meth] [--compact-occ] [--huge-pages] [-u INT] [-t INT] [idxbase]\n\n"
        "Options:\n"
        "  -d        destroy all indices in shared memory (matches bwa v1 behavior)\n"
        "  -l        list indices in shared memory; a third column flags\n"
        "            incomplete, missing or unusable entries, and `compact` marks\n"
        "            one staged with --compact-occ\n"
        "  --meth    stage a `bwa-mem3 index --meth` index — auto-appends\n"
        "            `.meth` to <idxbase>, mirroring `mem --meth`\n"
        "  --compact-occ\n"
        "            stage the FM-index occurrence checkpoints in the 32-byte\n"
        "            compact layout (about 3 GB less shared memory on hg38); `mem`\n"
        "            attaches it with or without --compact-occ. One layout per\n"
        "            index: run `shm -d` to restage in the other. BWA3_COMPACT_OCC\n"
        "            does not apply to `shm`. Not with -u\n"
        "  -u INT    densify the staged SA sample table to sample-rate shift INT\n"
        "            (in [0,6]) when it is denser than the on-disk rate. Trades a\n"
        "            one-time staging cost + extra shared memory for faster SA\n"
        "            resolution at `mem` time, WITHOUT rebuilding the on-disk\n"
        "            index. e.g. `-u 2` stages a stride-4 table from a stock\n"
        "            stride-8 (default) index. A value >= the disk rate is a\n"
        "            no-op (warns). Output is byte-identical either way.\n"
        "  -t, --threads INT\n"
        "            worker threads for the -u densify pass and the --compact-occ\n"
        "            conversion (default 1). The per-sample walks are independent,\n"
        "            so this scales near-linearly; has no effect without -u\n"
        "            actually densifying (a target denser than the on-disk rate)\n"
        "            or --compact-occ.\n"
        "  --huge-pages\n"
        "            stage the index on 1 GB huge pages: a System V shared memory\n"
        "            segment from the host's reserved 1 GB pool (10 pages for a\n"
        "            default hg38 stage), which `mem` attaches as usual. Cuts dTLB\n"
        "            page walks on Intel x86 and Graviton4 (arm64). Refused on AMD\n"
        "            x86, where it measured no gain. Linux only; needs CAP_IPC_LOCK\n"
        "            or membership of the vm.hugetlb_shm_group group, and one IPC\n"
        "            namespace for the stage, `mem` and `shm -d`. Fails, staging\n"
        "            nothing, when the pages cannot be had or the index is already\n"
        "            staged on other pages.\n"
        "  -h --help print this help and exit\n\n"
        "Stage with no flags: `bwa-mem3 shm <idxbase>` loads the index into\n"
        "POSIX shared memory; subsequent `bwa-mem3 mem <idxbase> ...` runs\n"
        "auto-attach instead of re-reading from disk. For meth indices, pass\n"
        "the same plain `<idxbase>` to all three commands plus `--meth` on\n"
        "`index`, `shm`, and `mem` (the c2t suffix is auto-appended).\n\n"
        "Footgun: if you re-build the index, run `bwa-mem3 shm -d` first.\n"
        "There is no staleness check -- a stale segment will silently mis-align.\n\n"
        "Interrupted stage: `shm -l` marks it `incomplete` and `mem` loads\n"
        "       from disk; rerun `bwa-mem3 shm <idxbase>` to reclaim and restage.\n"
        "       Stagers hold a named POSIX semaphore only for a few ms of\n"
        "       registry update; if one is killed inside that window, later\n"
        "       stages wait, and `bwa-mem3 shm -d` unlinks the semaphore.\n\n"
        "macOS: POSIX shm has implementation-defined per-segment caps; large\n"
        "       indices may simply fail to stage. Prefer Linux for production.\n"
        "Linux: /dev/shm defaults to ~50%% of RAM on bare metal; in containers\n"
        "       it is often much smaller and may need raising via --shm-size\n"
        "       (Docker) or an emptyDir tmpfs (Kubernetes).\n\n");
}

int main_shm(int argc, char *argv[])
{
    int c, to_list = 0, to_drop = 0, meth = 0, ret = 0;
    int target_sa_compx = -1;   /* -1 = stage the on-disk SA rate unchanged */
    int n_threads = 1;          /* worker count for the -u densify pass / compact conversion */
    bool compact = false;       /* --compact-occ (only `mem` reads BWA3_COMPACT_OCC) */
    int huge = 0;               /* --huge-pages */

    /* Pre-scan for --help / -h. getopt_long would treat -h as unknown
     * (no short opt declared) and print a generic error, but the top-level
     * usage advertises `bwa-mem3 <command> --help`, so we honor both. */
    for (int i = 1; i < argc; ++i) {
        const char *a = argv[i];
        /* Skip the value token of an argument-taking option so `-u --help` /
         * `-t --help` is a bad-value error from getopt below, not a spurious
         * help request. */
        if (std::strcmp(a, "-u") == 0 || std::strcmp(a, "-t") == 0 ||
            std::strcmp(a, "--threads") == 0) {
            ++i;
            continue;
        }
        if (std::strcmp(a, "--help") == 0 || std::strcmp(a, "-h") == 0) {
            print_shm_usage();
            return 0;
        }
    }

    /* getopt's global state may carry over from prior subcommand parsing.
     * Reset for a clean reparse. POSIX requires optind=1; macOS/BSD also
     * needs optreset=1 to clear the static `place` pointer inside getopt. */
    optind = 1;
    opterr = 1;
#if defined(__APPLE__) || defined(__FreeBSD__) || defined(__OpenBSD__) || defined(__NetBSD__)
    extern int optreset;
    optreset = 1;
#endif
    static struct option long_opts[] = {
        {"meth",    no_argument,       0, 1000},
        {"compact-occ", no_argument,   0, 1002},
        {"threads", required_argument, 0, 't' },
        {"huge-pages", no_argument,    0, 1001},
        {0,         0,                 0, 0   }
    };
    while ((c = getopt_long(argc, argv, "ldu:t:", long_opts, NULL)) >= 0) {
        if      (c == 'l')   to_list = 1;
        else if (c == 'd')   to_drop = 1;
        else if (c == 1000)  meth    = 1;
        else if (c == 1002)  compact = true;
        else if (c == 1001)  huge    = 1;
        else if (c == 't') {
            /* Worker count for the -u densify pass. Bare positive integer;
             * clamped to a sane upper bound so a fat-fingered value can't spin
             * up an absurd team. Meaningful only with -u (warned below). */
            char *end = NULL;
            errno = 0;
            long v = strtol(optarg, &end, 10);
            if (errno != 0 || end == optarg || *end != '\0' || v < 1 || v > 1024) {
                std::fprintf(stderr,
                    "[E::%s] --threads expects an integer in [1,1024], got '%s'\n",
                    __func__, optarg);
                return 1;
            }
            n_threads = (int)v;
        }
        else if (c == 'u') {
            /* Target SA sample-rate shift to densify to at stage time. Validate
             * as a bare non-negative integer in [0,6] (the writer-enforced
             * range; the sample period 1<<u must divide CP_BLOCK_SIZE=64).
             * Whether it is actually denser than THIS index's disk rate is
             * checked in bwa_shm_compute, which has read the disk rate. */
            char *end = NULL;
            errno = 0;
            long v = strtol(optarg, &end, 10);
            if (errno != 0 || end == optarg || *end != '\0' || v < 0 || v > 6) {
                std::fprintf(stderr,
                    "[E::%s] -u expects an integer in [0,6], got '%s'\n",
                    __func__, optarg);
                return 1;
            }
            target_sa_compx = (int)v;
        }
        else                 return 1;   /* getopt printed the error */
    }
    if (optind == argc && !to_list && !to_drop) {
        print_shm_usage();
        return 1;
    }
    if (optind < argc && (to_list || to_drop)) {
        std::fprintf(stderr,
            "[E::%s] -l or -d cannot be combined with idxbase\n", __func__);
        return 1;
    }
    if (target_sa_compx >= 0 && (to_list || to_drop)) {
        std::fprintf(stderr,
            "[E::%s] -u cannot be combined with -l or -d\n", __func__);
        return 1;
    }
    if (n_threads > 1 && (to_list || to_drop)) {
        std::fprintf(stderr,
            "[E::%s] --threads cannot be combined with -l or -d\n", __func__);
        return 1;
    }
    if (compact && (to_list || to_drop)) {
        std::fprintf(stderr,
            "[E::%s] --compact-occ cannot be combined with -l or -d\n", __func__);
        return 1;
    }
    if (compact && target_sa_compx >= 0) {
        std::fprintf(stderr,
            "[E::%s] --compact-occ cannot be combined with -u: the densify walk reads the "
            "index file's checkpoint layout, which a compact segment does not hold\n", __func__);
        return 1;
    }
    if (n_threads > 1 && target_sa_compx < 0 && !compact) {
        std::fprintf(stderr,
            "[W::%s] --threads only affects the -u densify pass and the "
            "--compact-occ conversion; ignoring without either\n", __func__);
    }
    if (huge && (to_list || to_drop)) {
        std::fprintf(stderr,
            "[E::%s] --huge-pages cannot be combined with -l or -d\n", __func__);
        return 1;
    }
    if (meth && (to_list || to_drop)) {
        std::fprintf(stderr,
            "[E::%s] --meth cannot be combined with -l or -d\n", __func__);
        return 1;
    }
    if (optind + 1 < argc) {
        std::fprintf(stderr,
            "[E::%s] too many positional arguments; expected at most one idxbase\n",
            __func__);
        return 1;
    }

    /* Stage. */
    if (optind < argc) {
        const char *user_prefix = argv[optind];
        char       resolved[PATH_MAX];
        if (meth) {
            if (resolve_meth_prefix(user_prefix, resolved) != 0) return 1;
            user_prefix = resolved;
        } else if (hint_missing_meth_flag(user_prefix)) {
            return 1;
        }
        /* --meth stages a SEED-only segment: omit PAC + REF_STRING (.0123),
         * which `mem --meth` never reads. Saves ~14.5 GB of shm on hg38. */
        if (bwa_shm_stage(user_prefix, /*bns_only=*/meth != 0, target_sa_compx,
                          n_threads, huge != 0, compact) < 0) {
            std::fprintf(stderr,
                "[E::%s] failed to stage '%s' in shared memory\n",
                __func__, user_prefix);
            ret = 1;
        }
    }
    if (to_list)  ret |= (bwa_shm_list()    < 0) ? 1 : 0;
    if (to_drop)  ret |= (bwa_shm_destroy() < 0) ? 1 : 0;
    return ret;
}

} /* extern "C" */
