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
#include "bwa.h"           /* BWA_CTL_SIZE — see static_assert below */
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
#if defined(__linux__)
#include <sys/vfs.h>       /* fstatfs, statfs — hugetlbfs mount check */
#endif

#ifndef HUGETLBFS_MAGIC
#define HUGETLBFS_MAGIC 0x958458f6   /* <linux/magic.h> */
#endif

/* /bwactl layout: 16-bit n_entries, 16-bit next_write_offset, then packed
 * (int64_t l_mem, char[] basename) entries. The header is 4 bytes total
 * (two uint16_t fields). next_write_offset starts at 4 on a fresh segment. */
#define BWA_SHM_CTL_HEADER_BYTES 4

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

/* The two constants must agree. bwa.h::BWA_CTL_SIZE is a v1-era stub that
 * we'll remove in a separate cleanup commit; until then, this assert
 * prevents silent drift if either definition is edited. */
static_assert(BWA_SHM_CTL_SIZE == BWA_CTL_SIZE,
              "BWA_SHM_CTL_SIZE must equal bwa.h::BWA_CTL_SIZE");

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
static sem_t *ctl_lock_acquire(void)
{
    sem_t *s = sem_open(BWA_SHM_LOCK_NAME, O_CREAT, 0644, 1);
    if (s == SEM_FAILED) {
        std::fprintf(stderr, "[E::%s] sem_open(%s) failed: %s\n",
                     __func__, BWA_SHM_LOCK_NAME, std::strerror(errno));
        return NULL;
    }
    while (sem_wait(s) < 0) {
        if (errno == EINTR) continue;
        std::fprintf(stderr, "[E::%s] sem_wait(%s) failed: %s\n",
                     __func__, BWA_SHM_LOCK_NAME, std::strerror(errno));
        sem_close(s);
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

/* Open `/bwactl` read-write, creating it if necessary and zero-initializing
 * a fresh segment with next-write offset = BWA_SHM_CTL_HEADER_BYTES.
 *
 * On success returns the mmap'd address (a `BWA_SHM_CTL_SIZE`-byte writable
 * mapping) and `*fd_out` holds the still-open file descriptor — caller
 * closes it after the mapping is no longer needed (typical pattern: close
 * immediately, since the mapping survives close).
 *
 * Concurrent RMW callers must hold ctl_lock_acquire across the open / read
 * / write / close cycle so the n_entries / next_write fields stay coherent
 * across processes. */
static void *ctl_open_rw(int *fd_out)
{
    int created = 0;
    int fd = shm_open(BWA_SHM_CTL_NAME, O_RDWR, 0);
    if (fd < 0) {
        if (errno != ENOENT) return NULL;
        fd = shm_open(BWA_SHM_CTL_NAME, O_CREAT | O_RDWR | O_EXCL, 0644);
        if (fd < 0) {
            fd = shm_open(BWA_SHM_CTL_NAME, O_RDWR, 0);
            if (fd < 0) return NULL;
        } else {
            created = 1;
        }
    }
    if (ftruncate(fd, BWA_SHM_CTL_SIZE) < 0) { close(fd); return NULL; }
    void *m = mmap(NULL, BWA_SHM_CTL_SIZE, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    if (m == MAP_FAILED) { close(fd); return NULL; }
    if (created) {
        std::memset(m, 0, BWA_SHM_CTL_SIZE);
        uint16_t next_write = BWA_SHM_CTL_HEADER_BYTES;
        std::memcpy((uint8_t *)m + 2, &next_write, sizeof(uint16_t));
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

uint64_t bwa_shm_hugetlb_min_bytes = (uint64_t)BWAMEM_HUGETLB_MIN_BYTES;

/* Whether the last component of `path` starts with BWA_SHM_HUGETLB_FILE_PREFIX.
 * Both the attach and the destroy follow a stub only to such a file, so a
 * stub can never make `shm -d` remove some unrelated file. */
static bool hugetlb_path_is_index_file(const char *path)
{
    const char *base = prefix_basename(path);
    return base != NULL &&
           std::strncmp(base, BWA_SHM_HUGETLB_FILE_PREFIX,
                        sizeof(BWA_SHM_HUGETLB_FILE_PREFIX) - 1) == 0 &&
           base[sizeof(BWA_SHM_HUGETLB_FILE_PREFIX) - 1] != '\0';
}

/* Field checks shared by encode and decode. */
static bool hugetlb_stub_fields_ok(const char *path, size_t path_cap,
                                   uint64_t total_size, uint64_t file_size)
{
    if (path == NULL || path[0] != '/') return false;
    if (strnlen(path, path_cap) == path_cap) return false;      /* no NUL */
    if (!hugetlb_path_is_index_file(path)) return false;
    if (total_size == 0 || total_size > file_size) return false;
    if (file_size % BWA_SHM_ONE_GIB != 0) return false;
    return true;
}

extern "C" int bwa_shm_hugetlb_stub_encode(const char *path, uint64_t total_size,
                                           uint64_t file_size, bwa_shm_hugetlb_stub_t *out)
{
    if (out == NULL) return -1;
    std::memset(out, 0, sizeof(*out));
    if (!hugetlb_stub_fields_ok(path, PATH_MAX, total_size, file_size)) return -1;
    out->magic      = BWA_SHM_HUGETLB_STUB_MAGIC;
    out->version    = BWA_SHM_HUGETLB_STUB_VERSION;
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
    if (stub.magic != BWA_SHM_HUGETLB_STUB_MAGIC) return -1;
    if (stub.version != BWA_SHM_HUGETLB_STUB_VERSION) return -1;
    if (!hugetlb_stub_fields_ok(stub.path, sizeof(stub.path),
                                stub.total_size, stub.file_size)) return -1;
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

extern "C" void bwa_shm_hugetlb_close(bwa_shm_hugetlb_file_t *f, bool unlink_file)
{
    if (f == NULL) return;
    if (f->base != NULL) munmap(f->base, f->map_len);
    /* fd >= 0 means this process created the file (O_EXCL), so only then is
     * it ours to remove. */
    if (f->fd >= 0) {
        close(f->fd);
        if (unlink_file && f->path[0] != '\0') unlink(f->path);
    }
    f->fd = -1;
    f->base = NULL;
    f->map_len = 0;
    f->path[0] = '\0';
}

#if defined(__linux__)

/* Whether `fd` (or the directory at `path` when fd < 0) is on a hugetlbfs
 * mount whose page size is 1 GiB. statfs reports a hugetlbfs mount's page
 * size as f_bsize. */
static bool hugetlbfs_1g_at(int fd, const char *path)
{
    struct statfs sf;
    const int rc = fd >= 0 ? fstatfs(fd, &sf) : statfs(path, &sf);
    return rc == 0 && (unsigned long)sf.f_type == (unsigned long)HUGETLBFS_MAGIC &&
           (uint64_t)sf.f_bsize == BWA_SHM_ONE_GIB;
}

/* The first hugetlbfs mount in /proc/mounts with 1 GB pages that this process
 * can create files in. Returns 0 with dir_out filled, or ENOENT. */
static int hugetlbfs_find_1g_mount(char *dir_out, size_t dir_len)
{
    FILE *fp = std::fopen("/proc/mounts", "r");
    if (fp == NULL) return ENOENT;
    char  *line = NULL;
    size_t cap = 0;
    int    rc = ENOENT;
    while (getline(&line, &cap, fp) >= 0) {
        if (bwa_shm_parse_hugetlbfs_1g_mount(line, dir_out, dir_len) == 1 &&
            access(dir_out, W_OK | X_OK) == 0 && hugetlbfs_1g_at(-1, dir_out)) {
            rc = 0;
            break;
        }
    }
    std::free(line);
    std::fclose(fp);
    if (rc != 0 && dir_len > 0) dir_out[0] = '\0';
    return rc;
}

extern "C" int bwa_shm_hugetlb_create(const char *dir, const char *name, uint64_t bytes,
                                      int nthreads, bwa_shm_hugetlb_file_t *out)
{
    out->fd = -1; out->base = NULL; out->map_len = 0; out->path[0] = '\0';
    if (dir == NULL || name == NULL || name[0] == '\0' || std::strchr(name, '/') != NULL)
        return EINVAL;
    const size_t map_len = bwamem_hugetlb_round_up((size_t)bytes);
    if (map_len == 0) return EINVAL;

    char mount_dir[PATH_MAX];
    if (dir[0] == '\0') {
        const int rc = hugetlbfs_find_1g_mount(mount_dir, sizeof(mount_dir));
        if (rc != 0) return rc;
    } else {
        /* The stub must hold an absolute path; resolve a relative one. */
        if (realpath(dir, mount_dir) == NULL) return errno != 0 ? errno : ENOENT;
        if (!hugetlbfs_1g_at(-1, mount_dir)) return ENOTSUP;
    }
    const int n = std::snprintf(out->path, sizeof(out->path), "%s/%s%s", mount_dir,
                                BWA_SHM_HUGETLB_FILE_PREFIX, name);
    if (n < 0 || (size_t)n >= sizeof(out->path)) { out->path[0] = '\0'; return ENAMETOOLONG; }

    out->fd = open(out->path, O_CREAT | O_EXCL | O_RDWR | O_CLOEXEC, 0644);
    if (out->fd < 0) return errno;   /* path kept for the caller's message */
    int err = 0;
    /* Guard against a mount that changed under us since the directory check. */
    if (!hugetlbfs_1g_at(out->fd, NULL)) err = ENOTSUP;
    /* hugetlbfs takes a size only in whole huge pages. */
    else if (ftruncate(out->fd, (off_t)map_len) != 0) err = errno;
    if (err == 0) {
        /* A shared hugetlbfs mapping reserves all its pages here (no
         * MAP_NORESERVE), so a pool or size= limit that is too small fails
         * now with ENOMEM instead of with SIGBUS during the fill. */
        void *m = mmap(NULL, map_len, PROT_READ | PROT_WRITE, MAP_SHARED, out->fd, 0);
        if (m == MAP_FAILED) err = errno != 0 ? errno : ENOMEM;
        else { out->base = m; out->map_len = map_len; }
    }
    if (err == 0) {
        /* Fault the pages in over several threads before the (serial) fill:
         * the kernel zeroes each 1 GB page in the fault that first touches
         * it, so the fill alone would zero them one at a time. The prefault
         * also turns a fault-time failure (a hugetlb cgroup limit) into an
         * error rather than a SIGBUS. On failure it unmaps the buffer. */
        bwamem_hugetlb_buf_t buf = {out->base, out->map_len};
        if (bwamem_hugetlb_prefault(&buf, nthreads) != 0) {
            out->base = NULL; out->map_len = 0;
            err = EFAULT;
        }
    }
    if (err != 0) bwa_shm_hugetlb_close(out, /*unlink_file=*/true);
    return err;
}

#else /* not Linux */

extern "C" int bwa_shm_hugetlb_create(const char *dir, const char *name, uint64_t bytes,
                                      int nthreads, bwa_shm_hugetlb_file_t *out)
{
    (void)dir; (void)name; (void)bytes; (void)nthreads;
    out->fd = -1; out->base = NULL; out->map_len = 0; out->path[0] = '\0';
    return ENOSYS;
}

#endif

/* Test-only hook: test/shm_thp_advise_test.cpp checks the stager's mapping
 * for the huge-page advice through it, since bwa_shm_stage unmaps that mapping
 * before it returns. Not set on any production path. */
void (*bwa_shm_stage_map_hook)(const void *base, size_t len) = NULL;

/* Walk the registry, calling `cb(l_mem, name, ctx)` for each entry. cb may
 * return 0 to continue or non-zero to stop. Returns the cb's last return
 * value (or 0 if the registry is absent or empty), -1 on open/mmap error. */
typedef int (*bwa_shm_ctl_cb_t)(int64_t l_mem, const char *name, void *ctx);

static int ctl_walk(bwa_shm_ctl_cb_t cb, void *ctx)
{
    int fd = shm_open(BWA_SHM_CTL_NAME, O_RDONLY, 0);
    if (fd < 0) return errno == ENOENT ? 0 : -1;
    /* A previous failed stage can leave a zero-sized registry; on macOS an
     * unlinked-but-still-referenced segment also lingers at size 0. mmap
     * rejects 0-byte mappings with EINVAL. ctl_open_rw always ftruncates
     * the segment to the full BWA_SHM_CTL_SIZE on creation, so anything
     * smaller is junk left behind by an aborted stage — treat as empty. */
    struct stat st;
    if (fstat(fd, &st) < 0 || (size_t)st.st_size < BWA_SHM_CTL_SIZE) {
        close(fd);
        return 0;
    }
    void *m = mmap(NULL, BWA_SHM_CTL_SIZE, PROT_READ, MAP_SHARED, fd, 0);
    close(fd);
    if (m == MAP_FAILED) return -1;

    const uint8_t *base = (const uint8_t *)m;
    uint16_t n_entries = 0;
    std::memcpy(&n_entries, base, sizeof(uint16_t));
    const uint8_t *p   = base + BWA_SHM_CTL_HEADER_BYTES;
    const uint8_t *end = base + BWA_SHM_CTL_SIZE;
    int rc = 0;
    for (uint16_t i = 0; i < n_entries; ++i) {
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

    munmap(m, BWA_SHM_CTL_SIZE);
    return rc;
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

    /* Peek the trailing sa_compx tag (written by write_fm_index_streaming,
     * see fm_index_writer.cpp) before sizing the SA sample arrays, via the
     * shared detect_sa_compx() helper (FMI_search.h/.cpp). Mirrors
     * FMI_search::load_index's disk-path detection exactly (same helper),
     * so bwa_shm_compute (which sizes the shm segment) agrees with the
     * loader it stages for. See detect_sa_compx()'s doc comment for why the
     * tail-detection heuristic is safe against legacy (no-tail) indexes. */
    {
        struct stat st;
        /* A failed fstat() must not silently fall back to file_size=0: that
         * would make detect_sa_compx() treat a genuinely non-default-rate
         * index as legacy (rate 3), sizing this staged segment's SA-sample
         * sections for the wrong rate. Fail the load instead of guessing,
         * mirroring FMI_search::load_index's disk-path handling of the same
         * fstat() call. */
        if (fstat(fileno(cp), &st) != 0) {
            fprintf(stderr, "[E::%s] fstat(%s) failed: %s\n",
                    __func__, cp_path, strerror(errno));
            err_fclose(cp);
            bwa_shm_layout_free(layout);
            return -1;
        }
        layout->disk_sa_compx = detect_sa_compx(fileno(cp), (int64_t)st.st_size, layout->reference_seq_len, SA_COMPX);
    }

    /* On-the-fly SA densification (bwa-mem3 shm -u INT). The SA sample table is
     * staged at layout->sa_compx; by default that equals the disk rate. A
     * target strictly denser than disk (and within the writer-enforced [0,6])
     * sizes the SA sections for the denser rate here; bwa_shm_pack_into then
     * synthesizes the added samples via an LF-walk against the disk samples.
     * A target that is not strictly denser (>= disk, or out of range) is
     * ignored with a warning and the disk rate is staged unchanged — a request
     * to "densify" to an equal/coarser rate is a no-op, not an error. */
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

    /* 3. REF_STRING (.0123) is NEVER staged. `mem` pac-fetches the original
     * reference from `.pac` on demand (bns_get_seq_v2), so the unpacked `.0123`
     * is never read from shm — for meth OR plain. Zero-length REF_STRING section
     * for every index → −6.4 GB per staged segment on hg38, and the `.0123` need
     * not exist on disk (commit "stop building .0123"). This is DECOUPLED from
     * `bns_only`: that flag still gates PAC below — the seed --meth segment drops
     * PAC too, but plain/original segments KEEP PAC (pac-fetch reads it). The
     * `.0123` is not stat'd, so staging works on an index that doesn't have it. */
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

    /* 1. Header. */
    bwa_shm_header_t hdr;
    memset(&hdr, 0, sizeof(hdr));
    hdr.magic       = BWA_SHM_MAGIC;
    hdr.version     = BWA_SHM_VERSION;
    hdr.n_sections  = layout->n_sections;
    hdr.total_size  = layout->total_size;
    memcpy(dest, &hdr, sizeof(hdr));

    /* 2. Section table. */
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

static int test_cb(int64_t /*l_mem*/, const char *name, void *ctx)
{
    return std::strcmp(name, (const char *)ctx) == 0 ? 1 : 0;
}

int bwa_shm_test(const char *prefix)
{
    if (prefix == NULL || prefix[0] == '\0') return -1;
    return ctl_walk(test_cb, (void *)prefix_basename(prefix));
}

static int list_cb(int64_t l_mem, const char *name, void * /*ctx*/)
{
    std::printf("%s\t%lld\n", name, (long long)l_mem);
    return 0;
}

int bwa_shm_list(void)
{
    int rc = ctl_walk(list_cb, NULL);
    return rc < 0 ? rc : 0;
}

/* If the POSIX segment `idx_path` is a hugetlb redirect stub, remove the
 * hugetlbfs file it points to. Only a file named bwaidx-* on a hugetlbfs
 * mount with 1 GB pages is removed (the stub decode checks the name), so a
 * stub cannot direct `shm -d` at any other file. Errors are ignored, like the
 * shm_unlink beside it. */
static void hugetlb_unlink_redirect_target(const char *idx_path)
{
    int fd = shm_open(idx_path, O_RDONLY, 0);
    if (fd < 0) return;
    struct stat st;
    if (fstat(fd, &st) == 0 && (size_t)st.st_size == sizeof(bwa_shm_hugetlb_stub_t)) {
        void *m = mmap(NULL, sizeof(bwa_shm_hugetlb_stub_t), PROT_READ, MAP_SHARED, fd, 0);
        if (m != MAP_FAILED) {
            bwa_shm_hugetlb_stub_t stub;
            if (bwa_shm_hugetlb_stub_decode(m, sizeof(stub), &stub) == 0) {
#if defined(__linux__)
                struct statfs sf;
                if (statfs(stub.path, &sf) == 0 &&
                    (unsigned long)sf.f_type == (unsigned long)HUGETLBFS_MAGIC)
                    unlink(stub.path);
#endif
            }
            munmap(m, sizeof(bwa_shm_hugetlb_stub_t));
        }
    }
    close(fd);
}

static int destroy_cb(int64_t /*l_mem*/, const char *name, void * /*ctx*/)
{
    char idx_path[PATH_MAX];
    path_concat2(idx_path, BWA_SHM_IDX_PREFIX, name);
    hugetlb_unlink_redirect_target(idx_path);
    shm_unlink(idx_path);   /* ignore errors (entry may already be gone) */
    return 0;
}

int bwa_shm_destroy(void)
{
    int rc = ctl_walk(destroy_cb, NULL);
    if (rc < 0) return rc;
    shm_unlink(BWA_SHM_CTL_NAME);
    /* Recover from a stuck lock: a stager that segfaulted or was kill -9'd
     * mid-stage holds the named semaphore until it is unlinked, blocking
     * every subsequent stager forever. POSIX has no SEM_UNDO, so unlinking
     * alongside the registry is the operator-visible recovery path. Errors
     * are ignored (ENOENT is the common case on a clean drop). */
    sem_unlink(BWA_SHM_LOCK_NAME);
    return 0;
}

/* /dev/shm capacity preflight for a POSIX-staged index; see bwa_shm_stage. */
static int devshm_preflight(const char *name, uint64_t total_size)
{
    /* ftruncate on a tmpfs file lazily reserves space — it succeeds even
     * when the requested size exceeds the tmpfs limit. The actual ENOSPC
     * surfaces during pack_into's freads as SIGBUS / EFAULT when a page
     * fault on a write to the mmap'd region cannot allocate a backing
     * page. The user-visible failure is then a confusing
     * `[fread] Bad address` from err_fread_noeof, with no indication
     * that /dev/shm was the cause.
     *
     * statvfs("/dev/shm") tells us the tmpfs's available bytes up-front,
     * so we can refuse cleanly with a message that points at the actual
     * fix. Default tmpfs size on Linux is RAM/2 (e.g. 16 GB on a 32 GB
     * c7a.4xlarge / c7i.4xlarge), so the hg38 FMI index (~17 GB
     * total_size) doesn't fit out of the box on common AWS instance
     * sizes. The remediation is a one-liner — `mount -o remount,size=…`
     * — but only useful if the user knows that's the problem. */
    struct statvfs svfs;
    if (statvfs("/dev/shm", &svfs) == 0) {
        /* f_bavail is in units of f_frsize (fundamental block size),
         * not f_bsize (preferred I/O block size); on Linux tmpfs the
         * two are equal but POSIX permits them to differ (e.g. ZFS
         * reports f_frsize=1MiB), so use f_frsize for portability. */
        uint64_t avail = (uint64_t)svfs.f_frsize * (uint64_t)svfs.f_bavail;
        if (total_size > avail) {
            std::fprintf(stderr,
                "[E::bwa_shm_stage] /dev/shm has %llu bytes free but staging '%s' "
                "needs %llu bytes. Remount with larger size, e.g. "
                "`sudo mount -o remount,size=%lluG /dev/shm`.\n",
                (unsigned long long)avail,
                name,
                (unsigned long long)total_size,
                (unsigned long long)((total_size + (1ULL<<30) - 1) / (1ULL<<30) + 2));
            return -1;
        }
    }
    /* statvfs failure (e.g. /dev/shm not present, ENOSYS on exotic
     * kernels) is non-fatal — fall through and let the existing
     * shm_open / ftruncate / pack_into path run. The pre-existing
     * (less helpful) error message in pack_into is still better
     * than refusing to stage based on a missing diagnostic. */
    return 0;
}

/* Explain why `shm --huge-pages` is staging on default pages after
 * bwa_shm_hugetlb_create returned `err`. */
static void report_hugetlb_fallback(int err, const char *dir,
                                    const bwa_shm_hugetlb_file_t *hf,
                                    uint64_t total_size)
{
    const double gib = (double)total_size / (double)BWA_SHM_ONE_GIB;
    const size_t pages = bwamem_hugetlb_round_up((size_t)total_size) >> 30;
    const char *where = dir[0] != '\0' ? dir : "the hugetlbfs mount";
    std::fprintf(stderr, "[M::bwa_shm_stage] --huge-pages: ");
    if (err == ENOSYS) {
        std::fprintf(stderr, "1 GB pages are supported only on Linux");
    } else if (err == ENOENT && dir[0] == '\0') {
        std::fprintf(stderr, "found no hugetlbfs mount with pagesize=1G that this user can "
                     "write in /proc/mounts (mount one, e.g. `sudo mkdir -p /mnt/huge-1g && "
                     "sudo mount -t hugetlbfs -o pagesize=1G,uid=$(id -u) none /mnt/huge-1g`)");
    } else if (err == ENOTSUP) {
        std::fprintf(stderr, "'%s' is not a hugetlbfs mount with pagesize=1G", dir);
    } else if (err == EEXIST) {
        std::fprintf(stderr, "%s already exists, probably left by an interrupted stage "
                     "(remove it to free its pages)", hf->path);
    } else if (err == ENOMEM) {
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
                     gib, hf->path[0] != '\0' ? hf->path : where, std::strerror(err));
    }
    std::fprintf(stderr, "; staging in POSIX shared memory instead\n");
}

/* Workers for faulting in the 1 GB pages of a hugetlb-staged segment. The
 * stage is one-off, so use the host's CPUs (capped like the index loader)
 * unless -t asked for more. */
static int hugetlb_prefault_threads(int n_threads)
{
    long ncpu = sysconf(_SC_NPROCESSORS_ONLN);
    int nt = ncpu > 0 ? (int)(ncpu < 64 ? ncpu : 64) : 1;
    return n_threads > nt ? n_threads : nt;
}

int bwa_shm_stage(const char *prefix, bool bns_only, int target_sa_compx,
                  int n_threads, const char *hugetlb_dir)
{
    if (prefix == NULL || prefix[0] == '\0') {
        std::fprintf(stderr, "[E::%s] empty prefix\n", __func__);
        return -1;
    }
    int already = bwa_shm_test(prefix);
    if (already < 0) {
        std::fprintf(stderr, "[E::%s] failed to probe registry\n", __func__);
        return -1;
    }
    if (already == 1) {
        std::fprintf(stderr, "[M::%s] index '%s' is already in shm\n",
                     __func__, prefix);
        return 0;
    }

    const char *name = prefix_basename(prefix);
    if (name[0] == '\0') {
        std::fprintf(stderr, "[E::%s] prefix has no basename component\n", __func__);
        return -1;
    }

    /* 1. Compute the layout (loads BNS; peeks scalars from the FMI file). */
    bwa_shm_layout_t layout;
    if (bwa_shm_compute(prefix, &layout, bns_only, target_sa_compx, n_threads) != 0) {
        std::fprintf(stderr, "[E::%s] failed to compute layout for '%s'\n",
                     __func__, prefix);
        return -1;
    }

    /* 2. Acquire the cross-process lock, then open the control segment R/W.
     *    Holding the lock across the read-modify-write window prevents two
     *    concurrent stagers from racing each other to append the same
     *    basename or stomping each other's next_write cursor. The pre-stage
     *    bwa_shm_test() above is unlocked (a cheap fast path); we re-scan
     *    the live registry below now that the lock is held to catch the
     *    case where another stager appended this prefix between our test
     *    and our acquire. */
    sem_t *ctl_lock = ctl_lock_acquire();
    if (ctl_lock == NULL) {
        bwa_shm_layout_free(&layout);
        return -1;
    }

    int   ctl_fd  = -1;
    void *ctl_map = ctl_open_rw(&ctl_fd);
    if (ctl_map == NULL) {
        std::fprintf(stderr, "[E::%s] failed to open control segment %s: %s\n",
                     __func__, BWA_SHM_CTL_NAME, std::strerror(errno));
        ctl_lock_release(ctl_lock);
        bwa_shm_layout_free(&layout);
        return -1;
    }

    /* 3. Re-scan the live registry for this basename. */
    uint16_t n_entries = 0, next_write = 0;
    std::memcpy(&n_entries,  (uint8_t *)ctl_map,     sizeof(uint16_t));
    std::memcpy(&next_write, (uint8_t *)ctl_map + 2, sizeof(uint16_t));
    {
        const uint8_t *p   = (const uint8_t *)ctl_map + BWA_SHM_CTL_HEADER_BYTES;
        const uint8_t *end = (const uint8_t *)ctl_map + BWA_SHM_CTL_SIZE;
        for (uint16_t i = 0; i < n_entries; ++i) {
            if (p + sizeof(int64_t) >= end) break;
            const char *entry_name = (const char *)(p + sizeof(int64_t));
            size_t name_max = (size_t)(end - (const uint8_t *)entry_name);
            size_t namelen  = strnlen(entry_name, name_max);
            if (namelen == name_max) break;
            if (std::strcmp(entry_name, name) == 0) {
                std::fprintf(stderr,
                    "[M::%s] index '%s' was staged by another process\n",
                    __func__, prefix);
                ctl_close(ctl_map, ctl_fd);
                ctl_lock_release(ctl_lock);
                bwa_shm_layout_free(&layout);
                return 0;
            }
            p = (const uint8_t *)entry_name + namelen + 1;
        }
    }

    /* 4. Bound-check the new entry. next_write is a uint16_t and
     *    BWA_SHM_CTL_SIZE is 0x10000 (65536), which uint16_t cannot represent.
     *    Reject the exact-fill case so the post-write cast back to uint16_t
     *    cannot wrap to 0 and have the next stage overwrite the header. */
    size_t name_bytes  = std::strlen(name) + 1;
    size_t entry_bytes = sizeof(int64_t) + name_bytes;
    if ((size_t)next_write + entry_bytes >= BWA_SHM_CTL_SIZE) {
        std::fprintf(stderr,
                     "[E::%s] control segment full (cannot stage '%s')\n",
                     __func__, name);
        ctl_close(ctl_map, ctl_fd);
        ctl_lock_release(ctl_lock);
        bwa_shm_layout_free(&layout);
        return -1;
    }

    /* 4b. Decide where the index goes. With --huge-pages it goes to a
     *     hugetlbfs file (step 5b) unless the index is too small for 1 GB
     *     pages; otherwise, or when that fails, into the POSIX segment, after
     *     the /dev/shm capacity preflight. */
    bool want_huge = hugetlb_dir != NULL;
    if (want_huge && layout.total_size <= bwa_shm_hugetlb_min_bytes) {
        std::fprintf(stderr,
            "[M::%s] --huge-pages: the index (%.2f GiB) is too small for 1 GB pages "
            "(%llu MiB or less); staging in POSIX shared memory\n",
            __func__, (double)layout.total_size / (double)BWA_SHM_ONE_GIB,
            (unsigned long long)(bwa_shm_hugetlb_min_bytes >> 20));
        want_huge = false;
    }
    if (!want_huge && devshm_preflight(name, layout.total_size) != 0) {
        ctl_close(ctl_map, ctl_fd);
        ctl_lock_release(ctl_lock);
        bwa_shm_layout_free(&layout);
        return -1;
    }

    /* 5. Create the per-index segment. EEXIST means another stager raced us
     *    to a segment with the same basename — treat that as already-staged
     *    rather than unlinking and re-creating, which would invalidate any
     *    process currently attached to it. */
    char idx_path[PATH_MAX];
    path_concat2(idx_path, BWA_SHM_IDX_PREFIX, name);
    int idx_fd = shm_open(idx_path, O_CREAT | O_RDWR | O_EXCL, 0644);
    if (idx_fd < 0 && errno == EEXIST) {
        std::fprintf(stderr,
            "[M::%s] segment %s already exists; treating as already staged\n",
            __func__, idx_path);
        ctl_close(ctl_map, ctl_fd);
        ctl_lock_release(ctl_lock);
        bwa_shm_layout_free(&layout);
        return 0;
    }
    if (idx_fd < 0) {
        std::fprintf(stderr, "[E::%s] shm_open(%s) failed: %s\n",
                     __func__, idx_path, std::strerror(errno));
        ctl_close(ctl_map, ctl_fd);
        ctl_lock_release(ctl_lock);
        bwa_shm_layout_free(&layout);
        return -1;
    }
    /* 5b. --huge-pages: the index goes to a 1 GB page hugetlbfs file, and the
     *     POSIX segment just claimed will hold a redirect to it (step 6b). */
    bwa_shm_hugetlb_file_t hfile = {-1, NULL, 0, {0}};
    bool on_huge = false;
    if (want_huge) {
        const int err = bwa_shm_hugetlb_create(hugetlb_dir, name, layout.total_size,
                                               hugetlb_prefault_threads(n_threads), &hfile);
        if (err == 0) {
            on_huge = true;
        } else {
            report_hugetlb_fallback(err, hugetlb_dir, &hfile, layout.total_size);
            if (devshm_preflight(name, layout.total_size) != 0) {
                close(idx_fd);
                shm_unlink(idx_path);
                ctl_close(ctl_map, ctl_fd);
                ctl_lock_release(ctl_lock);
                bwa_shm_layout_free(&layout);
                return -1;
            }
        }
    }

    void *idx_map = NULL;
    if (on_huge) {
        idx_map = hfile.base;
    } else if (ftruncate(idx_fd, (off_t)layout.total_size) < 0) {
        std::fprintf(stderr, "[E::%s] ftruncate(%s, %llu) failed: %s\n",
                     __func__, idx_path,
                     (unsigned long long)layout.total_size, std::strerror(errno));
        close(idx_fd);
        shm_unlink(idx_path);
        ctl_close(ctl_map, ctl_fd);
        ctl_lock_release(ctl_lock);
        bwa_shm_layout_free(&layout);
        return -1;
    } else {
        idx_map = idx_segment_map(idx_fd, (size_t)layout.total_size,
                                  PROT_READ | PROT_WRITE);
    }
    if (idx_map == MAP_FAILED) {
        std::fprintf(stderr, "[E::%s] mmap(%s) failed: %s\n",
                     __func__, idx_path, std::strerror(errno));
        close(idx_fd);
        shm_unlink(idx_path);
        ctl_close(ctl_map, ctl_fd);
        ctl_lock_release(ctl_lock);
        bwa_shm_layout_free(&layout);
        return -1;
    }
    if (bwa_shm_stage_map_hook) bwa_shm_stage_map_hook(idx_map, (size_t)layout.total_size);

    /* 6. Stream-pack directly into the mmap'd region. No heap intermediate. */
    if (bwa_shm_pack_into(&layout, (uint8_t *)idx_map) != 0) {
        std::fprintf(stderr, "[E::%s] pack_into failed for %s\n", __func__, prefix);
        if (on_huge) bwa_shm_hugetlb_close(&hfile, /*unlink_file=*/true);
        else         munmap(idx_map, (size_t)layout.total_size);
        close(idx_fd);
        shm_unlink(idx_path);
        ctl_close(ctl_map, ctl_fd);
        ctl_lock_release(ctl_lock);
        bwa_shm_layout_free(&layout);
        return -1;
    }
    /* 6b. For a hugetlb-staged index, write the redirect into the POSIX
     *     segment. The file keeps its pages after this process unmaps it. */
    char   huge_path[PATH_MAX];
    size_t huge_len = 0;
    if (on_huge) {
        std::memcpy(huge_path, hfile.path, sizeof(huge_path));
        huge_len = hfile.map_len;
        bwa_shm_hugetlb_stub_t stub;
        void *stub_map = MAP_FAILED;
        const bool ok =
            bwa_shm_hugetlb_stub_encode(hfile.path, layout.total_size, hfile.map_len, &stub) == 0 &&
            ftruncate(idx_fd, (off_t)sizeof(stub)) == 0 &&
            (stub_map = mmap(NULL, sizeof(stub), PROT_READ | PROT_WRITE, MAP_SHARED,
                             idx_fd, 0)) != MAP_FAILED;
        if (ok) {
            std::memcpy(stub_map, &stub, sizeof(stub));
            munmap(stub_map, sizeof(stub));
        }
        bwa_shm_hugetlb_close(&hfile, /*unlink_file=*/!ok);
        if (!ok) {
            std::fprintf(stderr, "[E::%s] failed to write the 1 GB page redirect %s: %s\n",
                         __func__, idx_path, std::strerror(errno));
            close(idx_fd);
            shm_unlink(idx_path);
            ctl_close(ctl_map, ctl_fd);
            ctl_lock_release(ctl_lock);
            bwa_shm_layout_free(&layout);
            return -1;
        }
    } else {
        munmap(idx_map, (size_t)layout.total_size);
    }
    close(idx_fd);

    /* 7. Append the entry to the control segment. */
    uint8_t *entry = (uint8_t *)ctl_map + next_write;
    int64_t l_mem = (int64_t)layout.total_size;
    std::memcpy(entry,                       &l_mem, sizeof(int64_t));
    std::memcpy(entry + sizeof(int64_t),     name,   name_bytes);
    n_entries  = (uint16_t)(n_entries + 1);
    next_write = (uint16_t)(next_write + entry_bytes);
    std::memcpy((uint8_t *)ctl_map,     &n_entries,  sizeof(uint16_t));
    std::memcpy((uint8_t *)ctl_map + 2, &next_write, sizeof(uint16_t));

    ctl_close(ctl_map, ctl_fd);
    ctl_lock_release(ctl_lock);

    if (on_huge)
        std::fprintf(stderr, "[M::%s] staged '%s' (%llu bytes) on %zu x 1 GB pages as %s "
                     "(redirect %s)\n",
                     __func__, name, (unsigned long long)layout.total_size,
                     huge_len >> 30, huge_path, idx_path);
    else
        std::fprintf(stderr, "[M::%s] staged '%s' (%llu bytes) as %s\n",
                     __func__, name,
                     (unsigned long long)layout.total_size, idx_path);
    bwa_shm_layout_free(&layout);
    return 0;
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

    /* A segment staged with `shm --huge-pages` is a redirect stub naming a
     * hugetlbfs file that holds the index: map that file instead. It is
     * already on 1 GB pages, so it gets no THP advice. */
    uint64_t magic = 0;
    std::memcpy(&magic, m, sizeof(magic));
    if (magic == BWA_SHM_HUGETLB_STUB_MAGIC) {
        bwa_shm_hugetlb_stub_t stub;
        const int bad = bwa_shm_hugetlb_stub_decode(m, len, &stub);
        munmap(m, len);
        if (bad) {
            std::fprintf(stderr, "[E::%s] malformed 1 GB page redirect in %s\n",
                         __func__, idx_path);
            return NULL;
        }
        int hfd = open(stub.path, O_RDONLY | O_CLOEXEC);
        if (hfd < 0) {
            std::fprintf(stderr, "[E::%s] open(%s) failed: %s\n",
                         __func__, stub.path, std::strerror(errno));
            return NULL;
        }
        if (fstat(hfd, &st) < 0 || (uint64_t)st.st_size != stub.file_size) {
            std::fprintf(stderr, "[E::%s] %s is not the size its redirect in %s records\n",
                         __func__, stub.path, idx_path);
            close(hfd);
            return NULL;
        }
        len = (size_t)stub.file_size;
        m = mmap(NULL, len, PROT_READ, MAP_SHARED, hfd, 0);
        close(hfd);
        if (m == MAP_FAILED) {
            std::fprintf(stderr, "[E::%s] mmap(%s, %llu) failed: %s\n",
                         __func__, stub.path, (unsigned long long)len, std::strerror(errno));
            return NULL;
        }
        std::fprintf(stderr, "[M::%s] index '%s' is on %zu x 1 GB pages (%s)\n",
                     __func__, name, len >> 30, stub.path);
    }

    /* Validate the on-segment header. */
    bwa_shm_header_t hdr;
    std::memcpy(&hdr, m, sizeof(hdr));
    if (hdr.magic != BWA_SHM_MAGIC) {
        std::fprintf(stderr, "[E::%s] bad magic in %s\n", __func__, idx_path);
        munmap(m, len);
        return NULL;
    }
    if (hdr.version != BWA_SHM_VERSION) {
        std::fprintf(stderr, "[E::%s] unsupported segment version %u in %s (expected %u)\n",
                     __func__, hdr.version, idx_path, BWA_SHM_VERSION);
        munmap(m, len);
        return NULL;
    }
    /* total_size is the exact data written; st_size may be page-rounded up.
     * Require total_size <= segment size to tolerate ftruncate's rounding. */
    if (hdr.total_size > (uint64_t)len) {
        std::fprintf(stderr, "[E::%s] header total_size=%llu exceeds segment size %llu in %s\n",
                     __func__, (unsigned long long)hdr.total_size,
                     (unsigned long long)len, idx_path);
        munmap(m, len);
        return NULL;
    }

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
        "  -l        list names of indices in shared memory\n"
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
        "            segment (11 for hg38). `mem` attaches the same way and -d\n"
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
        "Stuck-lock recovery: concurrent stagers are serialized by a named\n"
        "       POSIX semaphore. If a stager is kill -9'd mid-stage, the lock\n"
        "       persists and subsequent stages block forever. `bwa-mem3 shm -d`\n"
        "       unlinks the semaphore alongside the registry; rerun afterwards.\n\n"
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
