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
#include "bwa_hugepages.h" /* bwamem_hugetlb_round_up, bwamem_hugetlb_prefault, BWAMEM_HUGETLB_MIN_BYTES */

#include <stdio.h>
#include <cstdlib>
#include <cstring>
#include <sys/mman.h>      /* shm_open, mmap, munmap */
#include <sys/stat.h>      /* mode constants, stat */
#include <sys/statvfs.h>   /* statvfs — /dev/shm capacity preflight */
#include <fcntl.h>         /* O_* */
#include <unistd.h>        /* close, ftruncate */
#include <getopt.h>        /* getopt_long for `bwa-mem3 shm --meth ...` */
#include <semaphore.h>     /* sem_open, sem_wait, sem_post, sem_close, sem_unlink */
#include <errno.h>
#include <limits.h>        /* PATH_MAX */
#include <signal.h>        /* kill(pid, 0) liveness probe */
#include <strings.h>       /* strncasecmp */
#include <time.h>
#if defined(__linux__)
#include <sys/vfs.h>       /* fstatfs — hugetlbfs mount check */
#ifndef HUGETLBFS_MAGIC
#define HUGETLBFS_MAGIC 0x958458f6   /* <linux/magic.h> */
#endif
#endif

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

/* How long `shm -d` waits for the stager lock before proceeding without it. */
#ifndef BWA_SHM_DESTROY_LOCK_WAIT_MS
#define BWA_SHM_DESTROY_LOCK_WAIT_MS 5000
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

/* Acquire the cross-process lock that serializes read-modify-write
 * sequences on /bwactl. POSIX named semaphores work on both Linux and
 * macOS (unlike flock(2) on POSIX shm fds, which macOS rejects with
 * EOPNOTSUPP). Returns the open semaphore on success — held by us until
 * ctl_lock_release — or NULL on failure. */
static sem_t *ctl_lock_acquire(bool create = true)
{
    sem_t *s = create ? sem_open(BWA_SHM_LOCK_NAME, O_CREAT, 0644, 1)
                      : sem_open(BWA_SHM_LOCK_NAME, 0);
    if (s == SEM_FAILED && !create && errno == ENOENT) return NULL;   /* removed */
    if (s == SEM_FAILED) {
        const int err = errno;
        std::fprintf(stderr, "[E::%s] sem_open(%s) failed: %s\n",
                     __func__, BWA_SHM_LOCK_NAME, std::strerror(err));
        errno = err;
        return NULL;
    }
    while (sem_wait(s) < 0) {
        const int err = errno;
        if (err == EINTR) continue;
        std::fprintf(stderr, "[E::%s] sem_wait(%s) failed: %s\n",
                     __func__, BWA_SHM_LOCK_NAME, std::strerror(err));
        sem_close(s);
        errno = err;
        return NULL;
    }
    return s;
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

/* ---- 1 GB hugetlb staging: see bwa_shm.h ------------------------------ */

static const uint64_t BWA_SHM_ONE_GIB = 1ull << 30;

#if defined(__linux__)
/* Whether the open directory `fd` is on a hugetlbfs mount whose page size is
 * 1 GiB. statfs reports a hugetlbfs mount's page size as f_bsize. */
static bool hugetlbfs_1g_dir(int fd)
{
    struct statfs sf;
    return fstatfs(fd, &sf) == 0 &&
           (unsigned long)sf.f_type == (unsigned long)HUGETLBFS_MAGIC &&
           (uint64_t)sf.f_bsize == BWA_SHM_ONE_GIB;
}
#else
static bool hugetlbfs_1g_dir(int) { return false; }
#endif

/* Test seams, in the style of BWA_SHM_STATVFS: the hugetlbfs check on an open
 * directory, the page prefault, and the size below which an index stays in
 * /dev/shm. A test substitutes them to run the redirect paths on any OS
 * against an ordinary directory, without a 1 GB page pool. */
#ifndef BWA_SHM_HUGETLB_FS_OK
#define BWA_SHM_HUGETLB_FS_OK(dirfd) hugetlbfs_1g_dir(dirfd)
#endif
#ifndef BWA_SHM_HUGETLB_PREFAULT
#define BWA_SHM_HUGETLB_PREFAULT(buf, nthreads) bwamem_hugetlb_prefault(buf, nthreads)
#endif
#ifndef BWA_SHM_HUGETLB_MIN_BYTES
#define BWA_SHM_HUGETLB_MIN_BYTES ((uint64_t)BWAMEM_HUGETLB_MIN_BYTES)
#endif

/* A stub's path is absolute and NUL-terminated within the stub. */
static bool hugetlb_stub_path_ok(const char *path)
{
    return path != NULL && path[0] == '/' &&
           strnlen(path, BWA_SHM_HUGETLB_PATH_MAX) < BWA_SHM_HUGETLB_PATH_MAX;
}

/* The index fits its file, which is a whole number of 1 GiB pages. */
static bool hugetlb_stub_sizes_ok(uint64_t total_size, uint64_t file_size)
{
    return total_size != 0 && total_size <= file_size && file_size % BWA_SHM_ONE_GIB == 0;
}

extern "C" int bwa_shm_hugetlb_stub_encode(const char *path, uint64_t total_size,
                                           uint64_t file_size, bwa_shm_hugetlb_stub_t *out)
{
    if (out == NULL) return -1;
    std::memset(out, 0, sizeof(*out));
    if (!hugetlb_stub_path_ok(path) || !hugetlb_stub_sizes_ok(total_size, file_size))
        return -1;
    out->magic      = BWA_SHM_HUGETLB_STUB_MAGIC;
    out->version    = BWA_SHM_HUGETLB_STUB_VERSION;
    out->kind       = BWA_SHM_HUGETLB_KIND_1G;
    out->total_size = total_size;
    out->file_size  = file_size;
    std::memcpy(out->path, path, std::strlen(path) + 1);
    return 0;
}

extern "C" int bwa_shm_hugetlb_stub_decode(const void *buf, size_t len,
                                           bwa_shm_hugetlb_stub_t *out)
{
    if (buf == NULL || out == NULL || len != sizeof(bwa_shm_hugetlb_stub_t)) return -1;
    bwa_shm_hugetlb_stub_t stub;
    std::memcpy(&stub, buf, sizeof(stub));
    if (stub.magic != BWA_SHM_HUGETLB_STUB_MAGIC ||
        stub.version != BWA_SHM_HUGETLB_STUB_VERSION ||
        stub.kind != BWA_SHM_HUGETLB_KIND_1G) return -1;
    if (!hugetlb_stub_path_ok(stub.path) ||
        !hugetlb_stub_sizes_ok(stub.total_size, stub.file_size)) return -1;
    *out = stub;
    return 0;
}

/* Whether a hugetlbfs `pagesize=` value is 1 GiB. The kernel prints it in the
 * largest of K or M that divides it (1 GiB is "1024M"); also accept the "1G"
 * and "1048576K" spellings a mount command takes. */
static bool hugetlb_pagesize_is_1g(const char *v, size_t n)
{
    static const char *const ok[] = {"1024M", "1G", "1048576K"};
    for (const char *o : ok) {
        if (std::strlen(o) == n && strncasecmp(v, o, n) == 0) return true;
    }
    return false;
}

extern "C" int bwa_shm_parse_hugetlbfs_1g_mount(const char *line, char *dir_out, size_t dir_len)
{
    if (line == NULL || dir_out == NULL || dir_len == 0) return 0;
    /* Split the first four space-separated fields. */
    const char *f[4]; size_t fl[4];
    const char *p = line;
    for (int i = 0; i < 4; ++i) {
        while (*p == ' ' || *p == '\t') ++p;
        if (*p == '\0' || *p == '\n') return 0;
        f[i] = p;
        while (*p != '\0' && *p != ' ' && *p != '\t' && *p != '\n') ++p;
        fl[i] = (size_t)(p - f[i]);
    }
    if (fl[2] != 9 || std::strncmp(f[2], "hugetlbfs", 9) != 0) return 0;

    /* Options: comma-separated, find pagesize=. */
    bool is_1g = false;
    const char *o = f[3], *oend = f[3] + fl[3];
    while (o < oend) {
        const char *c = o;
        while (c < oend && *c != ',') ++c;
        static const char key[] = "pagesize=";
        const size_t klen = sizeof(key) - 1;
        if ((size_t)(c - o) > klen && std::strncmp(o, key, klen) == 0)
            is_1g = hugetlb_pagesize_is_1g(o + klen, (size_t)(c - o) - klen);
        o = c + 1;
    }
    if (!is_1g) return 0;

    /* Decode the directory: /proc/mounts escapes space, tab, newline and
     * backslash as a backslash and three octal digits. */
    size_t w = 0;
    for (size_t i = 0; i < fl[1]; ++i) {
        char ch = f[1][i];
        if (ch == '\\' && i + 3 < fl[1]) {
            const char *d = f[1] + i + 1;
            if (d[0] >= '0' && d[0] <= '3' && d[1] >= '0' && d[1] <= '7' &&
                d[2] >= '0' && d[2] <= '7') {
                ch = (char)(((d[0] - '0') << 6) | ((d[1] - '0') << 3) | (d[2] - '0'));
                i += 3;
            }
        }
        if (w + 1 >= dir_len) return 0;
        dir_out[w++] = ch;
    }
    dir_out[w] = '\0';
    return dir_out[0] == '/' ? 1 : 0;
}

#if defined(__linux__)
/* The first hugetlbfs mount in /proc/mounts with 1 GB pages that this process
 * can create files in. Returns 0 with dir_out filled, or ENOENT. */
static int hugetlbfs_find_1g_mount(char *dir_out, size_t dir_len)
{
    FILE *fp = std::fopen("/proc/mounts", "r");
    if (fp == NULL) return ENOENT;
    char  *line = NULL;
    size_t cap = 0;
    int    rc = ENOENT;
    while (rc != 0 && getline(&line, &cap, fp) >= 0) {
        if (bwa_shm_parse_hugetlbfs_1g_mount(line, dir_out, dir_len) != 1 ||
            access(dir_out, W_OK | X_OK) != 0) continue;
        const int dfd = open(dir_out, O_RDONLY | O_DIRECTORY | O_CLOEXEC);
        if (dfd < 0) continue;
        if (BWA_SHM_HUGETLB_FS_OK(dfd)) rc = 0;
        close(dfd);
    }
    std::free(line);
    std::fclose(fp);
    if (rc != 0 && dir_len > 0) dir_out[0] = '\0';
    return rc;
}
#endif

/* Split a stub's absolute `path` into its directory (into `dir`) and its last
 * component (returned, pointing into `path`). */
static const char *hugetlb_split_path(const char *path, char dir[BWA_SHM_HUGETLB_PATH_MAX])
{
    const char *base = prefix_basename(path);
    size_t dlen = (size_t)(base - path);
    if (dlen > 1) --dlen;                  /* drop the separator, but keep "/" */
    std::memcpy(dir, path, dlen);
    dir[dlen] = '\0';
    return base;
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

/* Segment states, from the header (see bwa_shm_header_t). A redirect stub
 * (see bwa_shm.h) is staging, dead or bad by the same header fields; once
 * complete, it is SEG_COMPLETE only if the file it names checks out, and
 * SEG_DANGLING if that file is gone. */
namespace {
enum seg_state { SEG_ABSENT, SEG_COMPLETE, SEG_STAGING, SEG_DEAD, SEG_BAD, SEG_DANGLING };

/* What segment_state learned about a redirect stub. */
struct seg_redirect {
    bool redirect = false;        /* the segment is a stub (staging or complete) */
    bwa_shm_hugetlb_stub_t stub;  /* its bytes, when `redirect` */
    uid_t owner = 0;              /* the segment's owner */
    const char *why = NULL;       /* complete stub, SEG_BAD or SEG_DANGLING: why */
};
} // namespace

/* Open the directory of the file a stub names and check the file in it, for
 * registry name `name`: it must be `bwaidx-<name>` (a stub for one index never
 * leads to another's file), a regular file and not a symlink, owned by `owner`
 * (the stub segment's owner, so a stub lets nobody touch a file its author
 * could not), in a directory on a hugetlbfs mount with 1 GB pages. Returns the
 * directory's fd, with *base the file's name in it and *st its lstat, so the
 * caller acts on the file through the directory it checked. Returns -1 with
 * errno set (ENOENT when the directory or file is gone) and *why saying which
 * check failed. */
static int hugetlb_target_open(const bwa_shm_hugetlb_stub_t *stub, const char *name,
                               uid_t owner, const char **base, struct stat *st,
                               const char **why)
{
    static const size_t plen = sizeof(BWA_SHM_HUGETLB_FILE_PREFIX) - 1;
    if (!hugetlb_stub_path_ok(stub->path)) {
        *why = "its path is not absolute";
        errno = EINVAL;
        return -1;
    }
    char dir[BWA_SHM_HUGETLB_PATH_MAX];
    const char *b = hugetlb_split_path(stub->path, dir);
    if (std::strncmp(b, BWA_SHM_HUGETLB_FILE_PREFIX, plen) != 0 ||
        std::strcmp(b + plen, name) != 0) {
        *why = "its file is not named " BWA_SHM_HUGETLB_FILE_PREFIX "<name>";
        errno = EINVAL;
        return -1;
    }
    const int dfd = open(dir, O_RDONLY | O_DIRECTORY | O_CLOEXEC);
    if (dfd < 0) {
        *why = errno == ENOENT ? "its directory is gone" : "its directory cannot be opened";
        return -1;
    }
    auto fail = [dfd, why](int err, const char *w) {
        close(dfd);
        *why = w;
        errno = err;
        return -1;
    };
    if (fstatat(dfd, b, st, AT_SYMLINK_NOFOLLOW) != 0) {
        const int err = errno;
        return fail(err, err == ENOENT ? "its file is gone" : "its file cannot be read");
    }
    if (!S_ISREG(st->st_mode))
        return fail(EINVAL, "its file is a symlink or not a regular file");
    if (st->st_uid != owner)
        return fail(EPERM, "its file is not owned by the redirect's owner");
    if (!BWA_SHM_HUGETLB_FS_OK(dfd))
        return fail(ENOTSUP, "its directory is not on a hugetlbfs mount with 1 GB pages");
    *base = b;
    *why = NULL;
    return dfd;
}

/* Whether the header `h` of a redirect's file (magic `magic`, read before the
 * rest) is the complete index its stub records: this version, and the stub's
 * stage token and size. That rejects a stale stub over a newer file, a
 * replaced file, and a file still being filled. */
static bool hugetlb_header_matches(const bwa_shm_header_t *h, uint64_t magic,
                                   const bwa_shm_hugetlb_stub_t *stub)
{
    return magic == BWA_SHM_MAGIC && h->version == BWA_SHM_VERSION &&
           h->reserved[1] == stub->stage_token && h->total_size == stub->total_size;
}

/* The state of a complete stub's file (see hugetlb_target_open): reads its
 * 64-byte header with pread, which on hugetlbfs instantiates no page. */
static seg_state hugetlb_target_state(const bwa_shm_hugetlb_stub_t *stub, const char *name,
                                      uid_t owner, const char **why)
{
    const char *base = NULL;
    struct stat st;
    const int dfd = hugetlb_target_open(stub, name, owner, &base, &st, why);
    if (dfd < 0) return errno == ENOENT ? SEG_DANGLING : SEG_BAD;
    seg_state state = SEG_BAD;
    if ((uint64_t)st.st_size != stub->file_size) {
        *why = "its file is not the size the redirect records";
    } else {
        const int fd = openat(dfd, base, O_RDONLY | O_NOFOLLOW | O_CLOEXEC);
        struct stat fst;
        bwa_shm_header_t h;
        if (fd >= 0 && fstat(fd, &fst) == 0 && fst.st_dev == st.st_dev &&
            fst.st_ino == st.st_ino &&
            pread(fd, &h, sizeof(h), 0) == (ssize_t)sizeof(h) &&
            hugetlb_header_matches(&h, h.magic, stub))
            state = SEG_COMPLETE;
        else
            *why = "its file does not hold the index the redirect records";
        if (fd >= 0) close(fd);
    }
    close(dfd);
    return state;
}

/* Probe `idx_path`, the segment of registry entry `name`; *hdr gets its header
 * (zeroed if none could be read; a stub's first 64 bytes otherwise). SEG_BAD
 * from a failed call sets *failed (if given) to its name, errno kept;
 * otherwise *failed is NULL and SEG_BAD means a foreign magic or version, or a
 * complete stub whose file does not check out (*rd says why). *rd, if given,
 * says whether the segment is a redirect stub and holds it. */
static seg_state segment_state(const char *idx_path, const char *name, bwa_shm_header_t *hdr,
                               const char **failed = NULL, seg_redirect *rd = NULL)
{
    const char *dummy;
    if (failed == NULL) failed = &dummy;
    *failed = NULL;
    seg_redirect local;
    if (rd == NULL) rd = &local;
    rd->redirect = false;
    rd->why = NULL;
    std::memset(hdr, 0, sizeof(*hdr));
    int fd = shm_open(idx_path, O_RDONLY, 0);
    if (fd < 0 && errno != ENOENT) *failed = "shm_open";
    if (fd < 0) return errno == ENOENT ? SEG_ABSENT : SEG_BAD;   /* keeps errno */
    struct stat st;
    const bool sized = fstat(fd, &st) == 0;
    if (!sized) *failed = "fstat";
    const bool small = sized && st.st_size < (off_t)sizeof(*hdr);   /* never written */
    /* Large enough to be a stub (macOS rounds a shm object up to whole pages,
     * so the size alone does not tell a stub from a small index). */
    const bool stub_sized = sized && st.st_size >= (off_t)sizeof(bwa_shm_hugetlb_stub_t);
    const size_t map_len = stub_sized ? sizeof(bwa_shm_hugetlb_stub_t) : sizeof(*hdr);
    void *m = sized && !small ? mmap(NULL, map_len, PROT_READ, MAP_SHARED, fd, 0)
                              : MAP_FAILED;
    const int err = errno;
    close(fd);
    errno = err;
    if (m == MAP_FAILED && sized && !small) *failed = "mmap";
    if (m == MAP_FAILED) return small ? SEG_DEAD : SEG_BAD;
    const uint64_t magic = __atomic_load_n((const uint64_t *)m, __ATOMIC_ACQUIRE);
    std::memcpy(hdr, m, sizeof(*hdr));
    if (stub_sized) std::memcpy(&rd->stub, m, sizeof(rd->stub));
    munmap(m, map_len);
    errno = 0;
    rd->owner = st.st_uid;
    /* Before its magic is published, a stub is told from an index segment by
     * its version and kind (which overlay the header's version and
     * n_sections): an index has version 0 until its header is written, then
     * BWA_SHM_VERSION, never the stub's. */
    rd->redirect = stub_sized &&
        (magic == BWA_SHM_HUGETLB_STUB_MAGIC ||
         (magic == 0 && hdr->version == BWA_SHM_HUGETLB_STUB_VERSION &&
          hdr->n_sections == BWA_SHM_HUGETLB_KIND_1G));
    if (rd->redirect && magic == BWA_SHM_HUGETLB_STUB_MAGIC) {
        bwa_shm_hugetlb_stub_t stub;
        if (bwa_shm_hugetlb_stub_decode(&rd->stub, sizeof(stub), &stub) != 0) {
            rd->why = "the redirect is malformed";
            return SEG_BAD;
        }
        const seg_state s = hugetlb_target_state(&rd->stub, name, rd->owner, &rd->why);
        errno = 0;
        return s;
    }
    /* Shorter than its header declares: bwa_shm_attach refuses it, so it is bad. */
    if (magic == BWA_SHM_MAGIC && hdr->version == BWA_SHM_VERSION)
        return hdr->total_size <= (uint64_t)st.st_size ? SEG_COMPLETE : SEG_BAD;
    if (magic != 0) return SEG_BAD;
    const uint64_t pid = hdr->reserved[0];
    const bool alive = pid != 0 && pid <= (uint64_t)INT_MAX &&
                       (kill((pid_t)pid, 0) == 0 || errno == EPERM);
    return alive ? SEG_STAGING : SEG_DEAD;
}

/* segment_state of registry entry `name`'s segment. */
static seg_state name_state(const char *name, seg_redirect *rd = NULL)
{
    char idx_path[PATH_MAX];
    path_concat2(idx_path, BWA_SHM_IDX_PREFIX, name);
    bwa_shm_header_t hdr;
    return segment_state(idx_path, name, &hdr, NULL, rd);
}

/* Sum, over entries whose segment is still being staged, the bytes each has yet
 * to allocate (ctx: uint64_t): l_mem less the segment's allocated blocks. Where
 * st_blocks is unavailable (macOS reports 0) the full l_mem counts. A 1 GB page
 * stage takes nothing from /dev/shm: its pages were reserved from the hugetlb
 * pool when it mapped its file. */
static int staging_bytes_cb(int64_t l_mem, const char *name, void *ctx)
{
    seg_redirect rd;
    if (l_mem <= 0 || name_state(name, &rd) != SEG_STAGING || rd.redirect) return 0;
    char idx_path[PATH_MAX];
    path_concat2(idx_path, BWA_SHM_IDX_PREFIX, name);
    uint64_t allocated = 0;
    const int fd = shm_open(idx_path, O_RDONLY, 0);
    if (fd >= 0) {
        struct stat st;
        if (fstat(fd, &st) == 0 && st.st_blocks > 0) allocated = (uint64_t)st.st_blocks * 512u;
        close(fd);
    }
    if (allocated < (uint64_t)l_mem) *(uint64_t *)ctx += (uint64_t)l_mem - allocated;
    return 0;
}

/* Remove the hugetlbfs file a stub names, through the directory that
 * hugetlb_target_open checked. Returns 0 when it is removed or already gone,
 * or when it fails those checks (it is not this stub's to remove; a warning
 * says so), and -1 when a removal failed. */
static int hugetlb_file_unlink_checked(const bwa_shm_hugetlb_stub_t *stub, const char *name,
                                       uid_t owner)
{
    const char *base = NULL, *why = NULL;
    struct stat st;
    const int dfd = hugetlb_target_open(stub, name, owner, &base, &st, &why);
    if (dfd < 0) {
        if (errno != ENOENT)
            std::fprintf(stderr, "[W::%s] not removing %.*s: %s\n", __func__,
                         (int)strnlen(stub->path, sizeof(stub->path)), stub->path, why);
        return 0;
    }
    int rc = 0;
    if (unlinkat(dfd, base, 0) != 0 && errno != ENOENT) {
        std::fprintf(stderr, "[W::%s] could not remove %s: %s; its 1 GB pages stay "
                     "reserved until it is removed\n", __func__, stub->path,
                     std::strerror(errno));
        rc = -1;
    }
    close(dfd);
    return rc;
}

/* Remove registry entry `name`'s segment at `idx_path`. Every removal of a
 * segment, by `shm -d`, by the reclaim of an interrupted stage or by a stage's
 * rollback, goes through here, and this is the only code that removes a
 * hugetlbfs file. A redirect's file goes first and the stub after it, so a
 * file is never left without the stub that names it. The file of a staging
 * stub is removed once the stub records that it was created (its header is
 * written last, so nothing more can be checked); the file of a complete stub
 * only if it is the index that stub records. Returns 0, or -1 when a file
 * could not be removed (its pages stay reserved); other errors are ignored,
 * since the segment may already be gone. */
static int segment_remove(const char *idx_path, const char *name)
{
    bwa_shm_header_t hdr;
    seg_redirect rd;
    const seg_state state = segment_state(idx_path, name, &hdr, NULL, &rd);
    int rc = 0;
    if (rd.redirect) {
        const bool staging = rd.stub.magic == 0;
        if (staging ? rd.stub.file_size != 0 : state == SEG_COMPLETE)
            rc = hugetlb_file_unlink_checked(&rd.stub, name, rd.owner);
        else if (!staging && state == SEG_BAD)
            std::fprintf(stderr, "[W::%s] removing the redirect %s but not the file it names, "
                         "%.*s: %s\n", __func__, idx_path,
                         (int)strnlen(rd.stub.path, sizeof(rd.stub.path)), rd.stub.path,
                         rd.why != NULL ? rd.why : "it does not check out");
    }
    shm_unlink(idx_path);
    return rc;
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
                    int target_sa_compx, int n_threads)
{
    if (prefix == NULL || layout == NULL) return -1;

    memset(layout, 0, sizeof(*layout));
    layout->densify_threads = (n_threads > 1) ? n_threads : 1;
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
     * adjustment in FMI_search::load_index. Each value must be in
     * [0, reference_seq_len] after adjustment. */
    for (int i = 0; i < 5; ++i) {
        if (layout->count[i] < 0 || layout->count[i] > layout->reference_seq_len) {
            fprintf(stderr,
                "[E::%s] %s: count[%d]=%lld out of bounds (0..%lld)\n",
                __func__, cp_path, i, (long long)layout->count[i],
                (long long)layout->reference_seq_len);
            err_fclose(cp);
            bwa_shm_layout_free(layout);
            return -1;
        }
        layout->count[i] += 1;
    }

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
    int64_t cp_occ_bytes = cp_occ_count * (int64_t)sizeof(CP_OCC);

    /* SA section sizing uses the STAGED stride (layout->sa_compx), which may be
     * denser than disk when -u densifies. */
    int64_t sa_count     = (layout->reference_seq_len >> layout->sa_compx) + 1;
    int64_t sa_ms_bytes  = sa_count     * (int64_t)sizeof(int8_t);
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
    static const uint32_t N_SECTIONS = 10;
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
    place(BWA_SHM_SEC_FMI_CP_OCC,  cp_occ_bytes,      1);
    place(BWA_SHM_SEC_FMI_SA_MS,   sa_ms_bytes,       2);
    place(BWA_SHM_SEC_FMI_SA_LS,   sa_ls_bytes,       3);
    place(BWA_SHM_SEC_BNS_STRUCT,  bns_struct_bytes,  4);
    place(BWA_SHM_SEC_BNS_AMBS,    bns_ambs_bytes,    5);
    place(BWA_SHM_SEC_BNS_ANNS,    bns_anns_bytes,    6);
    place(BWA_SHM_SEC_BNS_NAMES,   bns_names_bytes,   7);
    place(BWA_SHM_SEC_PAC,         pac_bytes,         8);
    place(BWA_SHM_SEC_REF_STRING,  layout->ref_string_len, 9);

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
     * exactly what we stage, so cp_occ / sa_ms / sa_ls stream contiguously
     * straight into their sections — one fp, three freads.
     *
     * Densify (sa_compx < disk_sa_compx): cp_occ still streams straight in, but
     * the disk SA arrays are at the coarser disk stride. Read them into heap
     * temporaries, then synthesize the denser sample table directly into the
     * (larger) sec[2]/sec[3] via FMI_search::densify_sa_into — an LF-walk over
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
    /* CP_OCC first (identical in both modes). */
    err_fread_noeof(dest + sec[1].offset, 1, (size_t)sec[1].size, cp);

    if (!densify) {
        err_fread_noeof(dest + sec[2].offset, 1, (size_t)sec[2].size, cp);
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
        fmi.densify_sa_into(
            (const CP_OCC *)(dest + sec[1].offset),
            layout->count,
            layout->reference_seq_len,
            layout->sentinel_index,
            disk_ms, disk_ls, layout->disk_sa_compx,
            (int8_t   *)(dest + sec[2].offset),
            (uint32_t *)(dest + sec[3].offset),
            layout->sa_compx,
            layout->densify_threads);

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
    hdr.n_sections  = layout->n_sections;
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

    uint8_t *buf = (uint8_t *)malloc(layout.total_size);
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

static int list_cb(int64_t l_mem, const char *name, void * /*ctx*/)
{
    const seg_state st = name_state(name);
    std::printf("%s\t%lld%s\n", name, (long long)l_mem,
                st == SEG_STAGING || st == SEG_DEAD ? "\tincomplete" :
                st == SEG_ABSENT || st == SEG_DANGLING ? "\tmissing" :
                st == SEG_BAD ? "\tunusable" : "");
    return 0;
}

int bwa_shm_list(void)
{
    int rc = ctl_walk(list_cb, NULL);
    return rc < 0 ? rc : 0;
}

/* ctx: int, set to 1 if a segment's hugetlbfs file could not be removed. */
static int destroy_cb(int64_t /*l_mem*/, const char *name, void *ctx)
{
    char idx_path[PATH_MAX];
    path_concat2(idx_path, BWA_SHM_IDX_PREFIX, name);
    if (segment_remove(idx_path, name) != 0) *(int *)ctx = 1;
    return 0;
}

int bwa_shm_destroy(void)
{
    /* Hold the lock across the walk, so no stager compacts the registry under
     * it and shifts a name's suffix into place as an entry. Registry updates
     * take milliseconds; a lock held past the wait was left by a stager killed
     * inside one (no SEM_UNDO), so proceed without it, and unlinking it below
     * is the recovery. A missing lock (ENOENT) is normal. */
    sem_t *lock = sem_open(BWA_SHM_LOCK_NAME, 0);
    bool held = false;
    for (int waited_ms = 0; lock != SEM_FAILED; waited_ms += 10) {
        if (sem_trywait(lock) == 0) { held = true; break; }
        if (errno != EAGAIN && errno != EINTR) break;
        if (waited_ms >= BWA_SHM_DESTROY_LOCK_WAIT_MS) {
            std::fprintf(stderr, "[W::%s] stager lock %s still held after %d ms; "
                         "removing everything without it\n", __func__, BWA_SHM_LOCK_NAME,
                         (int)BWA_SHM_DESTROY_LOCK_WAIT_MS);
            break;
        }
        usleep(10000);
    }
    int unremoved = 0;
    const int rc = ctl_walk(destroy_cb, &unremoved);
    if (rc >= 0) {
        shm_unlink(BWA_SHM_CTL_NAME);
        sem_unlink(BWA_SHM_LOCK_NAME);
    }
    if (lock != SEM_FAILED) {
        if (held) sem_post(lock);
        sem_close(lock);
    }
    /* A hugetlbfs file left behind pins pages nothing else reclaims. */
    return rc < 0 ? rc : unremoved ? -1 : 0;
}

/* Everything bwa_shm_stage holds, released on every return path; while
 * `created` is set, destruction also rolls the stage back. */
namespace {
struct stage_guard {
    const char *name, *idx_path;
    bwa_shm_layout_t *layout;
    sem_t *lock = NULL;
    void  *ctl = NULL;
    int    ctl_fd = -1;
    void  *idx = MAP_FAILED;             /* the index: the segment, or the hugetlbfs file */
    size_t idx_len = 0;                  /* length of the `idx` mapping */
    int    idx_fd = -1;
    bwa_shm_hugetlb_stub_t *stub = NULL; /* the segment's redirect stub, when mapped */
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
        segment_state(idx_path, name, &hdr);
        return !released || hdr.reserved[1] == layout->stage_token;
    }
    /* Roll back an unfinished stage, touching the segment only if it still
     * carries our token. Without the lock, leave the entry: `-l` shows it
     * `missing` and the next stage drops it. segment_remove follows our stub
     * to the hugetlbfs file only once the stub records that we created it. */
    ~stage_guard() {
        if (idx != MAP_FAILED) munmap(idx, idx_len);
        if (stub != NULL) munmap(stub, sizeof(*stub));
        if (idx_fd >= 0) close(idx_fd);
        if (created) {
            const bool locked = lock != NULL || lock_registry(/*create=*/false);
            if (still_ours()) {
                if (locked) ctl_remove(ctl, name);
                segment_remove(idx_path, name);
            }
        }
        unlock();
        bwa_shm_layout_free(layout);
    }
};
} // namespace

/* Explain why `shm --huge-pages` stages in POSIX shared memory instead: the
 * 1 GB page path failed with `err`. `dir` is what --huge-pages named ("" for
 * a search of /proc/mounts); `path` the file, once one was chosen. */
static void report_hugetlb_fallback(int err, const char *dir, const char *path,
                                    uint64_t total_size)
{
    const double gib = (double)total_size / (double)BWA_SHM_ONE_GIB;
    const size_t pages = bwamem_hugetlb_round_up((size_t)total_size) >> 30;
    const char *where = path != NULL && path[0] != '\0' ? path
                      : dir[0] != '\0' ? dir : "the hugetlbfs mount";
    std::fprintf(stderr, "[M::bwa_shm_stage] --huge-pages: ");
    if (err == ENOSYS) {
        std::fprintf(stderr, "1 GB pages are supported only on Linux");
    } else if (err == ENOENT && dir[0] == '\0') {
        std::fprintf(stderr, "found no hugetlbfs mount with pagesize=1G that this user can "
                     "write in /proc/mounts (mount one, e.g. `sudo mkdir -p /mnt/huge-1g && "
                     "sudo mount -t hugetlbfs -o pagesize=1G,uid=$(id -u) none /mnt/huge-1g`)");
    } else if (err == ENOTSUP) {
        std::fprintf(stderr, "'%s' is not a hugetlbfs mount with pagesize=1G",
                     dir[0] != '\0' ? dir : where);
    } else if (err == ENOMEM || err == ENOSPC) {
        std::fprintf(stderr, "could not reserve %zu x 1 GB pages for the index (%.2f GiB) on "
                     "%s (the host's 1 GB pool is missing or has too few free pages, see "
                     "free_hugepages minus resv_hugepages in "
                     "/sys/kernel/mm/hugepages/hugepages-1048576kB; or the mount's size= limit, "
                     "or a hugetlb cgroup reservation limit, was reached)",
                     pages, gib, where);
    } else if (err == EFAULT) {
        std::fprintf(stderr, "could not fault in the 1 GB pages for the index (%.2f GiB; "
                     "needs Linux 5.14+ for MADV_POPULATE_WRITE, or a hugetlb cgroup usage "
                     "limit was reached)", gib);
    } else {
        std::fprintf(stderr, "could not put the index (%.2f GiB) on 1 GB pages in %s (%s)",
                     gib, where, std::strerror(err));
    }
    std::fprintf(stderr, "; staging in POSIX shared memory instead\n");
}

/* Workers for faulting in the 1 GB pages of a hugetlb-staged index. The
 * stage is one-off, so use the host's CPUs (capped like the index loader)
 * unless -t asked for more. */
static int hugetlb_prefault_threads(int n_threads)
{
    long ncpu = sysconf(_SC_NPROCESSORS_ONLN);
    int nt = ncpu > 0 ? (int)(ncpu < 64 ? ncpu : 64) : 1;
    return n_threads > nt ? n_threads : nt;
}

/* Before the lock: whether an index of `total_size` bytes, asked onto 1 GB
 * pages with --huge-pages[=dir], goes there. Returns true with `path` set to
 * `<mount>/bwaidx-<name>`, or prints why not and returns false. */
static bool hugetlb_choose_path(const char *dir, const char *name, uint64_t total_size,
                                char path[BWA_SHM_HUGETLB_PATH_MAX])
{
    path[0] = '\0';
    if (total_size <= (uint64_t)BWA_SHM_HUGETLB_MIN_BYTES) {
        std::fprintf(stderr,
            "[M::bwa_shm_stage] --huge-pages: the index (%.2f GiB) is too small for 1 GB "
            "pages (%llu MiB or less); staging in POSIX shared memory\n",
            (double)total_size / (double)BWA_SHM_ONE_GIB,
            (unsigned long long)((uint64_t)BWA_SHM_HUGETLB_MIN_BYTES >> 20));
        return false;
    }
    char mount_dir[PATH_MAX];
    int err = 0;
    if (dir[0] == '\0') {
#if defined(__linux__)
        err = hugetlbfs_find_1g_mount(mount_dir, sizeof(mount_dir));
#else
        err = ENOSYS;
#endif
    } else if (realpath(dir, mount_dir) == NULL) {   /* the stub holds an absolute path */
        err = errno != 0 ? errno : ENOENT;
    } else {
        const int dfd = open(mount_dir, O_RDONLY | O_DIRECTORY | O_CLOEXEC);
        if (dfd < 0) {
            err = errno;
        } else {
#if defined(__linux__)
            if (!BWA_SHM_HUGETLB_FS_OK(dfd)) err = ENOTSUP;
#else
            if (!BWA_SHM_HUGETLB_FS_OK(dfd)) err = ENOSYS;
#endif
            close(dfd);
        }
    }
    if (err == 0) {
        const size_t mlen = std::strlen(mount_dir);
        const char *sep = mlen > 0 && mount_dir[mlen - 1] == '/' ? "" : "/";
        const int n = std::snprintf(path, BWA_SHM_HUGETLB_PATH_MAX, "%s%s%s%s", mount_dir, sep,
                                    BWA_SHM_HUGETLB_FILE_PREFIX, name);
        if (n < 0 || n >= PATH_MAX) {
            path[0] = '\0';
            err = ENAMETOOLONG;
        }
    }
    if (err != 0) {
        report_hugetlb_fallback(err, dir, NULL, total_size);
        return false;
    }
    return true;
}

/* Under the lock, with `g` holding the empty segment it just created: make the
 * segment a redirect stub for `path`, then create, size and map that file,
 * which reserves its 1 GB pages from the pool. Nothing here waits on the pool
 * being zeroed; that is the prefault, after the lock.
 *
 * Returns 0 with g.idx mapping the file. Returns 1 when the pages cannot be
 * had: the file and stub are removed, a note says why, and g holds a fresh
 * empty segment for a POSIX stage. Returns -1 on an error that ends the stage,
 * including a file already at `path` (never adopted or removed). */
static int hugetlb_reserve(stage_guard &g, const char *hugetlb_dir, const char *path)
{
    const bwa_shm_layout_t &layout = *g.layout;
    const size_t map_len = bwamem_hugetlb_round_up((size_t)layout.total_size);

    /* 1. The stub, while its magic and file_size stay 0. */
    bwa_shm_hugetlb_stub_t stub;
    if (bwa_shm_hugetlb_stub_encode(path, layout.total_size, map_len, &stub) != 0) {
        std::fprintf(stderr, "[E::bwa_shm_stage] cannot record %s in a redirect\n", path);
        return -1;
    }
    stub.magic       = 0;
    stub.file_size   = 0;
    stub.stager_pid  = layout.stager_pid;
    stub.stage_token = layout.stage_token;
    void *m = MAP_FAILED;
    if (ftruncate(g.idx_fd, (off_t)sizeof(stub)) != 0 ||
        (m = mmap(NULL, sizeof(stub), PROT_READ | PROT_WRITE, MAP_SHARED, g.idx_fd, 0))
            == MAP_FAILED) {
        std::fprintf(stderr, "[E::bwa_shm_stage] cannot write the 1 GB page redirect %s: %s\n",
                     g.idx_path, std::strerror(errno));
        return -1;
    }
    g.stub = (bwa_shm_hugetlb_stub_t *)m;
    std::memcpy(g.stub, &stub, sizeof(stub));

    /* 2. Create the file, through its directory. An existing file is no stage's
     *    of this registry (it would be named by a stub), so leave it alone. */
    char dir[BWA_SHM_HUGETLB_PATH_MAX];
    const char *base = hugetlb_split_path(path, dir);
    const int dfd = open(dir, O_RDONLY | O_DIRECTORY | O_CLOEXEC);
    int err = dfd < 0 ? errno : 0;
    int ffd = -1;
    if (err == 0) {
        ffd = openat(dfd, base, O_CREAT | O_EXCL | O_RDWR | O_CLOEXEC | O_NOFOLLOW, 0644);
        if (ffd < 0 && errno == EEXIST) {
            std::fprintf(stderr,
                "[E::bwa_shm_stage] %s already exists but no staged index names it: it may "
                "be left by a stage whose redirect in shared memory was removed (for "
                "example by an older bwa-mem3), or belong to a registry on another host or "
                "container. Remove it (`rm %s`) if no bwa-mem3 uses it, then retry\n",
                path, path);
            close(dfd);
            return -1;
        }
        if (ffd < 0) err = errno;
    }
    /* 3. Only now does the stub name the file, so a removal that follows the
     *    stub never touches a file this stage did not create. */
    if (ffd >= 0) g.stub->file_size = map_len;

    /* 4. Recheck the mount, size the file and map it shared. The mapping
     *    reserves its pages (no MAP_NORESERVE), so a short pool fails here
     *    with ENOMEM rather than with SIGBUS during the pack. */
    void *idx = MAP_FAILED;
    if (err == 0 && !BWA_SHM_HUGETLB_FS_OK(dfd)) err = ENOTSUP;
    if (err == 0 && ftruncate(ffd, (off_t)map_len) != 0) err = errno;
    if (err == 0) {
        idx = mmap(NULL, map_len, PROT_READ | PROT_WRITE, MAP_SHARED, ffd, 0);
        if (idx == MAP_FAILED) err = errno != 0 ? errno : ENOMEM;
    }
    if (ffd >= 0) close(ffd);
    if (err == 0) {
        close(dfd);
        g.idx = idx;
        g.idx_len = map_len;
        return 0;
    }

    /* 5. Fall back: remove our file, then the stub, and start over with a new
     *    segment (macOS sizes a shared memory object only once). */
    report_hugetlb_fallback(err, hugetlb_dir, path, layout.total_size);
    if (ffd >= 0 && unlinkat(dfd, base, 0) != 0 && errno != ENOENT)
        std::fprintf(stderr, "[W::bwa_shm_stage] could not remove %s: %s; its 1 GB pages "
                     "stay reserved until it is removed\n", path, std::strerror(errno));
    if (dfd >= 0) close(dfd);
    munmap(g.stub, sizeof(*g.stub));
    g.stub = NULL;
    close(g.idx_fd);
    shm_unlink(g.idx_path);
    g.idx_fd = shm_open(g.idx_path, O_CREAT | O_RDWR | O_EXCL, 0644);
    if (g.idx_fd < 0) {
        std::fprintf(stderr, "[E::bwa_shm_stage] shm_open(%s) failed: %s\n",
                     g.idx_path, std::strerror(errno));
        return -1;
    }
    return 1;
}

/* "Already staged" for bwa_shm_stage, with a note when it was asked for the
 * other backing: a staged index stays where it is. */
static void report_already_staged(const char *prefix, const char *hugetlb_dir, bool redirect)
{
    std::fprintf(stderr, "[M::bwa_shm_stage] index '%s' is already in shm\n", prefix);
    if (hugetlb_dir != NULL && !redirect)
        std::fprintf(stderr, "[M::bwa_shm_stage] --huge-pages: it stays in POSIX shared "
                     "memory; to restage it on 1 GB pages, run bwa-mem3 shm -d first\n");
    else if (hugetlb_dir == NULL && redirect)
        std::fprintf(stderr, "[M::bwa_shm_stage] it stays on 1 GB pages; to restage it in "
                     "POSIX shared memory, run bwa-mem3 shm -d first\n");
}

/* bwa_shm_stage's result when the 1 GB pages could not be faulted in after
 * the stage was registered: it has been rolled back, and the caller stages
 * again in POSIX shared memory. */
#define STAGE_RETRY_POSIX 1

static int stage_impl(const char *prefix, bool bns_only, int target_sa_compx,
                      int n_threads, const char *hugetlb_dir)
{
    if (prefix == NULL || prefix[0] == '\0') {
        std::fprintf(stderr, "[E::bwa_shm_stage] empty prefix\n");
        return -1;
    }
    const char *name = prefix_basename(prefix);
    if (name[0] == '\0') {
        std::fprintf(stderr, "[E::bwa_shm_stage] prefix has no basename component\n");
        return -1;
    }
    char idx_path[PATH_MAX];
    path_concat2(idx_path, BWA_SHM_IDX_PREFIX, name);

    /* 0. Already staged: skip the disk reads and the lock. Step 2 decides. */
    bwa_shm_header_t hdr;
    seg_redirect rd;
    if (bwa_shm_test(prefix) == 1 &&
        segment_state(idx_path, name, &hdr, NULL, &rd) == SEG_COMPLETE) {
        report_already_staged(prefix, hugetlb_dir, rd.redirect);
        return 0;
    }

    /* 1. Compute the layout (loads BNS; peeks scalars from the FMI file). */
    bwa_shm_layout_t layout;
    if (bwa_shm_compute(prefix, &layout, bns_only, target_sa_compx, n_threads) != 0) {
        std::fprintf(stderr, "[E::bwa_shm_stage] failed to compute layout for '%s'\n",
                     prefix);
        return -1;
    }
    layout.stager_pid  = (uint64_t)getpid();
    layout.stage_token = stage_token();
    stage_guard g(name, idx_path, &layout);

    /* 1b. --huge-pages: find the mount and the file's path, before the lock. */
    char huge_path[BWA_SHM_HUGETLB_PATH_MAX];
    huge_path[0] = '\0';
    const bool want_huge = hugetlb_dir != NULL &&
        hugetlb_choose_path(hugetlb_dir, name, layout.total_size, huge_path);

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
        state = segment_state(idx_path, name, &hdr, &failed, &rd);
        pid = (unsigned long long)hdr.reserved[0];
        if (state != SEG_STAGING) break;
        g.unlock();
        if (!announced)
            std::fprintf(stderr, "[M::bwa_shm_stage] '%s' is being staged by PID %llu; "
                         "waiting\n", name, pid);
        clock_gettime(CLOCK_MONOTONIC, &now);
        if (now.tv_sec >= deadline) {
            std::fprintf(stderr, "[E::bwa_shm_stage] gave up after %d s waiting for PID %llu "
                         "to stage '%s'; if PID %llu is not a running bwa-mem3 shm, run "
                         "bwa-mem3 shm -d\n", (int)BWA_SHM_STAGE_WAIT_S, pid, name, pid);
            return -1;
        }
        sleep(1);
    }
    const int probe_errno = errno;
    if (state == SEG_COMPLETE && ctl_find(g.ctl, name) != NULL) {
        report_already_staged(prefix, hugetlb_dir, rd.redirect);
        return 0;
    }
    if (state == SEG_COMPLETE || state == SEG_BAD) {   /* never adopt or replace it */
        if (failed != NULL)
            std::fprintf(stderr, "[E::bwa_shm_stage] %s(%s): %s\n", failed, idx_path,
                         std::strerror(probe_errno));
        else if (rd.redirect)
            std::fprintf(stderr, "[E::bwa_shm_stage] segment %s is a 1 GB page redirect to "
                         "%.*s %s%s%s; not replacing it. Remove it with shm_unlink(\"%s\") "
                         "after checking that no stager is active, and the file if no "
                         "bwa-mem3 uses it, then retry\n", idx_path,
                         (int)strnlen(rd.stub.path, sizeof(rd.stub.path)), rd.stub.path,
                         state == SEG_BAD ? "that does not check out ("
                                          : "without a registry entry",
                         state == SEG_BAD && rd.why != NULL ? rd.why : "",
                         state == SEG_BAD ? ")" : "", idx_path);
        else
            std::fprintf(stderr, "[E::bwa_shm_stage] segment %s exists %s; not replacing it. "
                         "Remove it with shm_unlink(\"%s\") after checking that no stager "
                         "is active, then retry\n", idx_path,
                         state == SEG_BAD ? "with a foreign magic or version, or "
                                            "shorter than its header declares"
                                          : "without a registry entry", idx_path);
        return -1;
    }
    if (state == SEG_DEAD) {
        std::fprintf(stderr, "[M::bwa_shm_stage] reclaiming '%s' left by an interrupted stage "
                     "(PID %llu)\n", name, pid);
        segment_remove(idx_path, name);
    } else if (state == SEG_DANGLING) {
        std::fprintf(stderr, "[M::bwa_shm_stage] reclaiming '%s': %s, the file its 1 GB page "
                     "redirect names, is gone\n", name, rd.stub.path);
        segment_remove(idx_path, name);
    }
    ctl_remove(g.ctl, name);   /* a stale entry, if any */

    /* Bound-check the new entry. next_write is a uint16_t and BWA_SHM_CTL_SIZE
     * (65536) does not fit one, so reject the exact-fill case too: the cast
     * back must not wrap to 0 and let the next stage overwrite the header. */
    ctl_header h = ctl_header_load(g.ctl);
    const size_t name_bytes  = std::strlen(name) + 1;
    const size_t entry_bytes = sizeof(int64_t) + name_bytes;
    if ((size_t)h.next_write + entry_bytes >= BWA_SHM_CTL_SIZE) {
        std::fprintf(stderr, "[E::bwa_shm_stage] control segment full (cannot stage '%s')\n",
                     name);
        return -1;
    }
    g.idx_fd = shm_open(idx_path, O_CREAT | O_RDWR | O_EXCL, 0644);
    if (g.idx_fd < 0) {
        std::fprintf(stderr, "[E::bwa_shm_stage] shm_open(%s) failed: %s\n",
                     idx_path, std::strerror(errno));
        return -1;
    }
    g.created = true;

    bool on_huge = false;
    if (want_huge) {
        const int r = hugetlb_reserve(g, hugetlb_dir, huge_path);
        if (r < 0) return -1;
        on_huge = r == 0;
    }
    if (!on_huge) {
        /* /dev/shm preflight: tmpfs ftruncate succeeds past the limit and the
         * pack would later die of SIGBUS, so refuse up front with the fix.
         * Stages in progress allocate pages as they write, so their unallocated
         * remainder is reserved too. Sampled before free space, so pages written
         * in between count twice rather than not at all (conservative).
         * f_bavail is in f_frsize units. A statvfs failure is not fatal. */
        uint64_t staging = 0;
        ctl_walk_map(g.ctl, staging_bytes_cb, &staging);
        struct statvfs svfs;
        if (BWA_SHM_STATVFS("/dev/shm", &svfs) == 0) {
            uint64_t avail = (uint64_t)svfs.f_frsize * (uint64_t)svfs.f_bavail;
            avail = avail > staging ? avail - staging : 0;
            if ((uint64_t)layout.total_size > avail) {
                std::fprintf(stderr,
                    "[E::bwa_shm_stage] /dev/shm has %llu bytes free but staging '%s' "
                    "needs %llu bytes. Remount with larger size, e.g. "
                    "`sudo mount -o remount,size=%lluG /dev/shm`.\n",
                    (unsigned long long)avail, name,
                    (unsigned long long)layout.total_size,
                    (unsigned long long)((layout.total_size + (1ULL<<30) - 1) / (1ULL<<30) + 2));
                return -1;
            }
        }
        if (ftruncate(g.idx_fd, (off_t)layout.total_size) < 0) {
            std::fprintf(stderr, "[E::bwa_shm_stage] ftruncate(%s, %llu) failed: %s\n",
                         idx_path, (unsigned long long)layout.total_size,
                         std::strerror(errno));
            return -1;
        }
        g.idx = idx_segment_map(g.idx_fd, (size_t)layout.total_size, PROT_READ | PROT_WRITE);
        g.idx_len = (size_t)layout.total_size;
        if (g.idx == MAP_FAILED) {
            std::fprintf(stderr, "[E::bwa_shm_stage] mmap(%s) failed: %s\n",
                         idx_path, std::strerror(errno));
            return -1;
        }
        ((bwa_shm_header_t *)g.idx)->reserved[0] = layout.stager_pid;   /* magic stays 0 */
        ((bwa_shm_header_t *)g.idx)->reserved[1] = layout.stage_token;
    }

    uint8_t *entry = (uint8_t *)g.ctl + h.next_write;
    const int64_t l_mem = (int64_t)layout.total_size;
    std::memcpy(entry,                   &l_mem, sizeof(int64_t));
    std::memcpy(entry + sizeof(int64_t), name,   name_bytes);
    h.n_entries  = (uint16_t)(h.n_entries + 1);
    h.next_write = (uint16_t)(h.next_write + entry_bytes);
    ctl_header_store(g.ctl, h);   /* after the entry bytes */
    g.unlock();
    g.released = true;

    /* 3. Without the lock: fault in the 1 GB pages (the kernel zeroes each in
     *    the fault that first touches it; over several threads, not one at a
     *    time in the pack), then pack. The magic, written last, marks the index
     *    complete; for a redirect, the stub's magic after the file's. */
    if (on_huge) {
        bwamem_hugetlb_buf_t buf = { g.idx, g.idx_len };
        if (BWA_SHM_HUGETLB_PREFAULT(&buf, hugetlb_prefault_threads(n_threads)) != 0) {
            g.idx = MAP_FAILED;   /* the prefault unmapped it */
            g.idx_len = 0;
            report_hugetlb_fallback(EFAULT, hugetlb_dir, huge_path, layout.total_size);
            return STAGE_RETRY_POSIX;
        }
    }
    if (bwa_shm_stage_map_hook) bwa_shm_stage_map_hook(g.idx, (size_t)layout.total_size);
    if (bwa_shm_pack_into(&layout, (uint8_t *)g.idx) != 0) {
        std::fprintf(stderr, "[E::bwa_shm_stage] pack_into failed for %s\n", prefix);
        return -1;
    }
    if (on_huge)
        __atomic_store_n(&g.stub->magic, BWA_SHM_HUGETLB_STUB_MAGIC, __ATOMIC_RELEASE);

    /* 4. `shm -d`, or a stager that misjudged us dead, may have removed our
     *    entry or segment meanwhile; check without recreating anything. */
    const bool locked = g.lock_registry(/*create=*/false);
    if (!locked && errno != ENOENT) {
        std::fprintf(stderr, "[E::bwa_shm_stage] cannot retake the stager lock to finish "
                     "'%s': %s\n", name, std::strerror(errno));
        return -1;
    }
    if (!locked || ctl_find(g.ctl, name) == NULL || !g.still_ours()) {
        std::fprintf(stderr, "[E::bwa_shm_stage] '%s' was removed during staging (by bwa-mem3 "
                     "shm -d or another stager)\n", name);
        return -1;
    }
    g.created = false;
    if (on_huge)
        std::fprintf(stderr, "[M::bwa_shm_stage] staged '%s' (%llu bytes) on %zu x 1 GB pages "
                     "as %s (redirect %s)\n", name, (unsigned long long)layout.total_size,
                     g.idx_len >> 30, huge_path, idx_path);
    else
        std::fprintf(stderr, "[M::bwa_shm_stage] staged '%s' (%llu bytes) as %s\n",
                     name, (unsigned long long)layout.total_size, idx_path);
    return 0;
}

int bwa_shm_stage(const char *prefix, bool bns_only, int target_sa_compx,
                  int n_threads, const char *hugetlb_dir)
{
    int rc = stage_impl(prefix, bns_only, target_sa_compx, n_threads, hugetlb_dir);
    if (rc == STAGE_RETRY_POSIX)
        rc = stage_impl(prefix, bns_only, target_sa_compx, n_threads, NULL);
    return rc;
}

/* Map the hugetlbfs file a complete stub names, read-only, after the checks of
 * hugetlb_target_open and a size check on the opened file. It is already on
 * 1 GB pages, so it gets no transparent huge page advice. Returns the mapping
 * with *len its length, or MAP_FAILED after saying why. The caller validates
 * the file's header. */
static void *hugetlb_map_target(const bwa_shm_hugetlb_stub_t *stub, const char *name,
                                uid_t owner, size_t *len)
{
    const char *base = NULL, *why = NULL;
    struct stat st;
    const int dfd = hugetlb_target_open(stub, name, owner, &base, &st, &why);
    if (dfd < 0) {
        std::fprintf(stderr, "[E::bwa_shm_attach] not following the 1 GB page redirect for "
                     "'%s' to %s: %s\n", name, stub->path, why);
        return MAP_FAILED;
    }
    void *m = MAP_FAILED;
    const int fd = openat(dfd, base, O_RDONLY | O_NOFOLLOW | O_CLOEXEC);
    close(dfd);
    struct stat fst;
    if (fd < 0 || fstat(fd, &fst) != 0 || fst.st_dev != st.st_dev || fst.st_ino != st.st_ino ||
        (uint64_t)fst.st_size != stub->file_size) {
        std::fprintf(stderr, "[E::bwa_shm_attach] %s changed, or is not the size its "
                     "redirect records\n", stub->path);
    } else {
        m = mmap(NULL, (size_t)stub->file_size, PROT_READ, MAP_SHARED, fd, 0);
        if (m == MAP_FAILED)
            std::fprintf(stderr, "[E::bwa_shm_attach] mmap(%s, %llu) failed: %s\n",
                         stub->path, (unsigned long long)stub->file_size,
                         std::strerror(errno));
    }
    if (fd >= 0) close(fd);
    if (m != MAP_FAILED) *len = (size_t)stub->file_size;
    return m;
}

uint8_t *bwa_shm_attach(const char *prefix, size_t *len_out, size_t *map_len_out)
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

    /* A complete redirect stub (`shm --huge-pages`): map the hugetlbfs file it
     * names instead, and validate that file's header below. A stub still being
     * staged has magic 0 and falls through to the "not fully staged" path. */
    bwa_shm_hugetlb_stub_t stub;
    bool redirect = false;
    if (len >= sizeof(stub) &&
        __atomic_load_n((const uint64_t *)m, __ATOMIC_ACQUIRE) == BWA_SHM_HUGETLB_STUB_MAGIC) {
        const int bad = bwa_shm_hugetlb_stub_decode(m, sizeof(stub), &stub);
        munmap(m, len);
        if (bad) {
            std::fprintf(stderr, "[E::%s] malformed 1 GB page redirect in %s\n",
                         __func__, idx_path);
            return NULL;
        }
        m = hugetlb_map_target(&stub, name, st.st_uid, &len);
        if (m == MAP_FAILED) return NULL;
        redirect = true;
    }

    /* Validate the on-segment header. Acquire-load the magic first: it is
     * published last, so once it reads BWA_SHM_MAGIC the rest is in place. A
     * zero magic is a stage still in progress or one that was interrupted. */
    const uint64_t magic = __atomic_load_n((const uint64_t *)m, __ATOMIC_ACQUIRE);
    bwa_shm_header_t hdr;
    std::memcpy(&hdr, m, sizeof(hdr));
    if (magic != BWA_SHM_MAGIC) {
        if (magic == 0)
            std::fprintf(stderr, "[M::%s] '%s' is not fully staged (being staged or "
                         "interrupted); loading from disk\n", __func__, name);
        else
            std::fprintf(stderr, "[E::%s] bad magic in %s\n", __func__,
                         redirect ? stub.path : idx_path);
        munmap(m, len);
        return NULL;
    }
    if (hdr.version != BWA_SHM_VERSION) {
        std::fprintf(stderr, "[E::%s] unsupported segment version %u in %s (expected %u)\n",
                     __func__, hdr.version, redirect ? stub.path : idx_path, BWA_SHM_VERSION);
        munmap(m, len);
        return NULL;
    }
    /* total_size is the exact data written; st_size may be page-rounded up.
     * Require total_size <= segment size to tolerate ftruncate's rounding. */
    if (hdr.total_size > (uint64_t)len) {
        std::fprintf(stderr, "[E::%s] header total_size=%llu exceeds segment size %llu in %s\n",
                     __func__, (unsigned long long)hdr.total_size,
                     (unsigned long long)len, redirect ? stub.path : idx_path);
        munmap(m, len);
        return NULL;
    }
    if (redirect && !hugetlb_header_matches(&hdr, magic, &stub)) {
        std::fprintf(stderr, "[E::%s] %s does not hold the index its redirect %s records\n",
                     __func__, stub.path, idx_path);
        munmap(m, len);
        return NULL;
    }

    if (redirect)
        std::fprintf(stderr, "[M::%s] index '%s' is on %zu x 1 GB pages (%s)\n",
                     __func__, name, len >> 30, stub.path);
    if (len_out != NULL) *len_out = (size_t)hdr.total_size;
    if (map_len_out != NULL) *map_len_out = len;
    return (uint8_t *)m;
}

int bwa_shm_section_find(const uint8_t *base, uint32_t kind,
                         uint64_t *offset_out, uint64_t *size_out)
{
    if (base == NULL) return -1;
    const bwa_shm_header_t *hdr = (const bwa_shm_header_t *)base;
    if (hdr->magic   != BWA_SHM_MAGIC)   return -1;
    if (hdr->version != BWA_SHM_VERSION) return -1;

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
        "\nUsage: bwa-mem3 shm [-d|-l|--help] [--meth] [-u INT] [-t INT]\n"
        "                    [--huge-pages[=DIR]] [idxbase]\n\n"
        "Options:\n"
        "  -d        destroy all indices in shared memory (matches bwa v1 behavior)\n"
        "  -l        list indices in shared memory; a third column flags\n"
        "            incomplete, missing or unusable entries\n"
        "  --meth    stage a `bwa-mem3 index --meth` index — auto-appends\n"
        "            `.meth` to <idxbase>, mirroring `mem --meth`\n"
        "  -u INT    densify the staged SA sample table to sample-rate shift INT\n"
        "            (in [0,6]) when it is denser than the on-disk rate. Trades a\n"
        "            one-time staging cost + extra shared memory for faster SA\n"
        "            resolution at `mem` time, WITHOUT rebuilding the on-disk\n"
        "            index. e.g. `-u 2` stages a stride-4 table from a stock\n"
        "            stride-8 (default) index. A value >= the disk rate is a\n"
        "            no-op (warns). Output is byte-identical either way.\n"
        "  -t, --threads INT\n"
        "            worker threads for the -u densify pass (default 1). The\n"
        "            per-sample walks are independent, so this scales near-\n"
        "            linearly; has no effect unless -u actually densifies (a\n"
        "            target denser than the on-disk rate).\n"
        "  --huge-pages[=DIR]\n"
        "            Linux: stage the index on 1 GB huge pages, as a file on the\n"
        "            hugetlbfs mount DIR (default: the first one with\n"
        "            pagesize=1G this user can write, from /proc/mounts), instead\n"
        "            of in /dev/shm. Needs 1 GB pages reserved for the whole\n"
        "            index (11 for hg38). `mem` attaches the same way, and -d\n"
        "            removes the file. Falls back to /dev/shm, with a note, when\n"
        "            that is not possible or the index is 512 MiB or less.\n"
        "            Output is byte-identical either way.\n"
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
        "       stages block, and `bwa-mem3 shm -d` unlinks the semaphore.\n\n"
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
    int n_threads = 1;          /* worker count for the -u densify pass */
    int huge_pages = 0;         /* --huge-pages: stage on 1 GB hugetlb pages */
    const char *hugetlb_dir = "";   /* --huge-pages=DIR; "" = find a mount */

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
        {"threads", required_argument, 0, 't' },
        {"huge-pages", optional_argument, 0, 1001},
        {0,         0,                 0, 0   }
    };
    while ((c = getopt_long(argc, argv, "ldu:t:", long_opts, NULL)) >= 0) {
        if      (c == 'l')   to_list = 1;
        else if (c == 'd')   to_drop = 1;
        else if (c == 1000)  meth    = 1;
        else if (c == 1001) {
            huge_pages = 1;
            if (optarg != NULL) {
                if (optarg[0] == '\0') {
                    std::fprintf(stderr, "[E::%s] --huge-pages= needs a directory\n",
                                 __func__);
                    return 1;
                }
                hugetlb_dir = optarg;
            }
        }
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
    if (n_threads > 1 && target_sa_compx < 0) {
        std::fprintf(stderr,
            "[W::%s] --threads only affects the -u densify pass; ignoring "
            "without -u\n", __func__);
    }
    if (huge_pages && (to_list || to_drop)) {
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
                          n_threads, huge_pages ? hugetlb_dir : NULL) < 0) {
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
