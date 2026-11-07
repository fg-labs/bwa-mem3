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
#include <time.h>

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

/* Segment states, from the header (see bwa_shm_header_t). */
namespace { enum seg_state { SEG_ABSENT, SEG_COMPLETE, SEG_STAGING, SEG_DEAD, SEG_BAD }; }

/* Probe `idx_path`; *hdr gets its header (zeroed if none could be read).
 * SEG_BAD from a failed call sets *failed (if given) to its name, errno kept;
 * otherwise *failed is NULL and SEG_BAD means a foreign magic or version. */
static seg_state segment_state(const char *idx_path, bwa_shm_header_t *hdr,
                               const char **failed = NULL)
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
    munmap(m, sizeof(*hdr));
    errno = 0;
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
static seg_state name_state(const char *name)
{
    char idx_path[PATH_MAX];
    path_concat2(idx_path, BWA_SHM_IDX_PREFIX, name);
    bwa_shm_header_t hdr;
    return segment_state(idx_path, &hdr);
}

/* Sum, over entries whose segment is still being staged, the bytes each has yet
 * to allocate (ctx: uint64_t): l_mem less the segment's allocated blocks. Where
 * st_blocks is unavailable (macOS reports 0) the full l_mem counts. */
static int staging_bytes_cb(int64_t l_mem, const char *name, void *ctx)
{
    if (l_mem <= 0 || name_state(name) != SEG_STAGING) return 0;
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
                st == SEG_ABSENT ? "\tmissing" : st == SEG_BAD ? "\tunusable" : "");
    return 0;
}

int bwa_shm_list(void)
{
    int rc = ctl_walk(list_cb, NULL);
    return rc < 0 ? rc : 0;
}

static int destroy_cb(int64_t /*l_mem*/, const char *name, void * /*ctx*/)
{
    char idx_path[PATH_MAX];
    path_concat2(idx_path, BWA_SHM_IDX_PREFIX, name);
    shm_unlink(idx_path);   /* ignore errors (entry may already be gone) */
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
    const int rc = ctl_walk(destroy_cb, NULL);
    if (rc >= 0) {
        shm_unlink(BWA_SHM_CTL_NAME);
        sem_unlink(BWA_SHM_LOCK_NAME);
    }
    if (lock != SEM_FAILED) {
        if (held) sem_post(lock);
        sem_close(lock);
    }
    return rc < 0 ? rc : 0;
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
    void  *idx = MAP_FAILED;
    int    idx_fd = -1;
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
        if (idx != MAP_FAILED) munmap(idx, (size_t)layout->total_size);
        if (idx_fd >= 0) close(idx_fd);
        if (created) {
            const bool locked = lock != NULL || lock_registry(/*create=*/false);
            if (still_ours()) {
                if (locked) ctl_remove(ctl, name);
                shm_unlink(idx_path);
            }
        }
        unlock();
        bwa_shm_layout_free(layout);
    }
};
} // namespace

int bwa_shm_stage(const char *prefix, bool bns_only, int target_sa_compx,
                  int n_threads)
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

    /* 0. Already staged: skip the disk reads and the lock. Step 2 decides. */
    bwa_shm_header_t hdr;
    if (bwa_shm_test(prefix) == 1 && segment_state(idx_path, &hdr) == SEG_COMPLETE) {
        std::fprintf(stderr, "[M::%s] index '%s' is already in shm\n", __func__, prefix);
        return 0;
    }

    /* 1. Compute the layout (loads BNS; peeks scalars from the FMI file). */
    bwa_shm_layout_t layout;
    if (bwa_shm_compute(prefix, &layout, bns_only, target_sa_compx, n_threads) != 0) {
        std::fprintf(stderr, "[E::%s] failed to compute layout for '%s'\n",
                     __func__, prefix);
        return -1;
    }
    layout.stager_pid  = (uint64_t)getpid();
    layout.stage_token = stage_token();
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
        std::fprintf(stderr, "[M::%s] index '%s' is already in shm\n", __func__, prefix);
        return 0;
    }
    if (state == SEG_COMPLETE || state == SEG_BAD) {   /* never adopt or replace it */
        if (failed != NULL)
            std::fprintf(stderr, "[E::%s] %s(%s): %s\n", __func__, failed, idx_path,
                         std::strerror(probe_errno));
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
        std::fprintf(stderr, "[M::%s] reclaiming '%s' left by an interrupted stage (PID %llu)\n",
                     __func__, name, pid);
        shm_unlink(idx_path);
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
     * twice rather than not at all (conservative). f_bavail is in f_frsize units.
     * A statvfs failure is not fatal. */
    uint64_t staging = 0;
    ctl_walk_map(g.ctl, staging_bytes_cb, &staging);
    struct statvfs svfs;
    if (BWA_SHM_STATVFS("/dev/shm", &svfs) == 0) {
        uint64_t avail = (uint64_t)svfs.f_frsize * (uint64_t)svfs.f_bavail;
        avail = avail > staging ? avail - staging : 0;
        if ((uint64_t)layout.total_size > avail) {
            std::fprintf(stderr,
                "[E::%s] /dev/shm has %llu bytes free but staging '%s' "
                "needs %llu bytes. Remount with larger size, e.g. "
                "`sudo mount -o remount,size=%lluG /dev/shm`.\n",
                __func__, (unsigned long long)avail, name,
                (unsigned long long)layout.total_size,
                (unsigned long long)((layout.total_size + (1ULL<<30) - 1) / (1ULL<<30) + 2));
            return -1;
        }
    }
    if (ftruncate(g.idx_fd, (off_t)layout.total_size) < 0) {
        std::fprintf(stderr, "[E::%s] ftruncate(%s, %llu) failed: %s\n",
                     __func__, idx_path,
                     (unsigned long long)layout.total_size, std::strerror(errno));
        return -1;
    }
    g.idx = idx_segment_map(g.idx_fd, (size_t)layout.total_size, PROT_READ | PROT_WRITE);
    if (g.idx == MAP_FAILED) {
        std::fprintf(stderr, "[E::%s] mmap(%s) failed: %s\n",
                     __func__, idx_path, std::strerror(errno));
        return -1;
    }
    ((bwa_shm_header_t *)g.idx)->reserved[0] = layout.stager_pid;   /* magic stays 0 */
    ((bwa_shm_header_t *)g.idx)->reserved[1] = layout.stage_token;

    uint8_t *entry = (uint8_t *)g.ctl + h.next_write;
    const int64_t l_mem = (int64_t)layout.total_size;
    std::memcpy(entry,                   &l_mem, sizeof(int64_t));
    std::memcpy(entry + sizeof(int64_t), name,   name_bytes);
    h.n_entries  = (uint16_t)(h.n_entries + 1);
    h.next_write = (uint16_t)(h.next_write + entry_bytes);
    ctl_header_store(g.ctl, h);   /* after the entry bytes */
    g.unlock();
    g.released = true;

    /* 3. Pack without the lock; the magic, written last, marks it complete. */
    if (bwa_shm_stage_map_hook) bwa_shm_stage_map_hook(g.idx, (size_t)layout.total_size);
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
    if (!locked || ctl_find(g.ctl, name) == NULL || !g.still_ours()) {
        std::fprintf(stderr, "[E::%s] '%s' was removed during staging (by bwa-mem3 shm -d "
                     "or another stager)\n", __func__, name);
        return -1;
    }
    g.created = false;
    std::fprintf(stderr, "[M::%s] staged '%s' (%llu bytes) as %s\n",
                 __func__, name, (unsigned long long)layout.total_size, idx_path);
    return 0;
}

uint8_t *bwa_shm_attach(const char *prefix, size_t *len_out)
{
    if (len_out != NULL) *len_out = 0;
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
    const uint64_t magic = __atomic_load_n((const uint64_t *)m, __ATOMIC_ACQUIRE);
    bwa_shm_header_t hdr;
    std::memcpy(&hdr, m, sizeof(hdr));
    if (magic != BWA_SHM_MAGIC) {
        if (magic == 0)
            std::fprintf(stderr, "[M::%s] '%s' is not fully staged (being staged or "
                         "interrupted); loading from disk\n", __func__, name);
        else
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
        "\nUsage: bwa-mem3 shm [-d|-l|--help] [--meth] [-u INT] [-t INT] [idxbase]\n\n"
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
        {0,         0,                 0, 0   }
    };
    while ((c = getopt_long(argc, argv, "ldu:t:", long_opts, NULL)) >= 0) {
        if      (c == 'l')   to_list = 1;
        else if (c == 'd')   to_drop = 1;
        else if (c == 1000)  meth    = 1;
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
                          n_threads) < 0) {
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
