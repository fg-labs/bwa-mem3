/*************************************************************************************
                           The MIT License

   BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),
   Copyright (C) 2019  Intel Corporation, Heng Li.

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

Authors: Sanchit Misra <sanchit.misra@intel.com>; Vasimuddin Md <vasimuddin.md@intel.com>;
*****************************************************************************************/

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <climits>
#include <cstring>
#include <vector>
#include <atomic>         /* BWA3_KS_DEDUP stats counters + auto-controller state */
#include <algorithm>      /* std::fill for the ks-dedup slot table */
#include <chrono>         /* ks-dedup auto controller: per-chunk net timing */
#include <mutex>          /* ks-dedup auto controller: pooled Welford + once_flag */
#include <cmath>          /* sqrt/fabs/isfinite for the net-cycles z-test + parse */
#include <cstdarg>
#include <pthread.h>
#include <unistd.h>       /* pread, _exit */
#include <sys/mman.h>     /* munmap */
#include <sys/stat.h>     /* fstat */
#include <fcntl.h>        /* open (per-worker index fds), posix_fadvise */
#include "bwa_madvise.h"
#include "bwa_shm.h"
#include "utils.h"        /* ATTRIBUTE, err_fread_noeof */
#include "io_utils.h"     /* io_request_size, IO_MAX_ONCE */
#include "FMI_search.h"
#include "ksort.h"        /* smem_sort_rid_mn: per-rid introsort on (key, idx) pairs */
#include "smem_parked.h"
#include "profiling.h"
#include "libsais_build.h"

/*
 * Parallel index load.
 *
 * The bwa-mem3 FM-index (cp_occ + compressed SA samples, ~10 GB on hg38) is
 * otherwise slurped by a single-threaded fread whose warm-cache cost is one
 * core's page-fault + memcpy bandwidth out of the page cache (~0.8 s on hg38).
 * Splitting each big array across a few workers is memory-bandwidth bound and
 * cuts that ~4x. pread (not read) means no worker touches a shared file
 * offset, so no locking is needed; each worker also reads through its own
 * descriptor (see open_worker_fd).
 *
 * The destination stays the _mm_malloc'd + MADV_HUGEPAGE buffer, so the
 * transparent-hugepage coverage the hot Occ-lookup loop relies on is preserved
 * -- unlike an mmap alias of the file, which lands on non-THP pages and slows
 * alignment enough to erase the load saving. (An shm segment is advised the
 * same way, but gets huge pages only if /dev/shm is mounted huge=advise or
 * huge=within_size; see idx_segment_map in bwa_shm.cpp.)
 */

/* Declared in FMI_search.h; see there for the contract. Uses a floor division
 * so the chunk size never dips below FMI_PREAD_MIN_CHUNK: at `nbytes` a whole
 * multiple of the floor, `nbytes / min_chunk + 1` would hand out one chunk more
 * than the array can cover at that size (8 MB read across 2 workers = 4 MB
 * chunks). Only reachable for small arrays -- a caller asking for more workers
 * than a test-sized reference can feed, or a large BWA3_LOAD_THREADS. */
int fmi_pread_worker_count(size_t nbytes, int nthreads)
{
    if (nthreads < 1) nthreads = 1;
    size_t max_threads = nbytes / FMI_PREAD_MIN_CHUNK;
    if (max_threads < 1) max_threads = 1;
    if ((size_t)nthreads > max_threads) nthreads = (int)max_threads;
    return nthreads;
}

size_t fmi_pread_request_size(size_t remaining)
{
    // Delegate to the one shared clamp so read and write agree on the cap; see
    // io_utils.h. pread() has the same macOS >INT_MAX EINVAL cap as pwrite().
    return io_request_size(remaining, IO_MAX_ONCE);
}

namespace {

struct PreadChunk { int fd; char *dst; size_t nbytes; off_t off; };

/* Bail out of a failed chunk read.
 *
 * _exit, not exit: this runs on a worker thread while its siblings are still
 * inside pread(). exit() would run atexit handlers and static destructors
 * (mimalloc teardown, htslib) against live threads, and two workers reaching it
 * at once -- what a truncated index produces, since every chunk past the real
 * EOF fails together -- is undefined. stderr is unbuffered, but flush anyway so
 * the diagnostic survives if a caller made it buffered.
 *
 * noreturn keeps the caller's `if (r < 0) { ... }` terminating the way the
 * inline exit() it replaced did; format(printf) keeps -Wformat checking the
 * call sites, matching utils.h's err_fatal family. */
void pread_chunk_fail(const char *fmt, ...)
    ATTRIBUTE((noreturn)) ATTRIBUTE((format(printf, 1, 2)));

void pread_chunk_fail(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fflush(stderr);
    _exit(EXIT_FAILURE);
}

/* Open a private read descriptor on `path` for one pread worker, or return
 * `shared_fd` if that fails or `path` no longer names the file behind
 * `shared_fd` (replaced since it was opened -- reading it would splice two
 * files' bytes into one array). A private descriptor is what keeps the load
 * fast on S3 Mountpoint: its prefetcher keeps one sequential cursor per open
 * file and resets on every jump, so workers sharing a descriptor stall at
 * random-read throughput (awslabs/mountpoint-s3#1886). */
int open_worker_fd(const char *path, int shared_fd)
{
    int fd = open(path, O_RDONLY | O_CLOEXEC);
    if (fd < 0) return shared_fd;
    struct stat mine, shared;
    if (fstat(fd, &mine) != 0 || fstat(shared_fd, &shared) != 0
        || mine.st_dev != shared.st_dev || mine.st_ino != shared.st_ino) {
        close(fd);
        return shared_fd;
    }
    return fd;
}

void *pread_chunk_worker(void *arg)
{
    PreadChunk *c = static_cast<PreadChunk *>(arg);
    size_t done = 0;
    /* macOS pread() rejects any count > INT_MAX with EINVAL, so a chunk larger
     * than 2GiB fails outright -- with the 8-worker cap that is every index over
     * ~17.2GB (the hg38 --meth FM-index is 20.9GB). Request a bounded slice per
     * call; the loop already accumulates partial reads. Linux has no such limit,
     * which is why this only ever bit local macOS runs. */
    while (done < c->nbytes) {
        ssize_t r = pread(c->fd, c->dst + done,
                          fmi_pread_request_size(c->nbytes - done),
                          c->off + (off_t)done);
        if (r < 0) {
            if (errno == EINTR) continue;
            pread_chunk_fail("ERROR: pread failed during index load: %s\n", strerror(errno));
        }
        if (r == 0) {  // index file shorter than the header says it should be
            // Counts are for THIS worker's chunk, not the whole array: a
            // truncated index trips every chunk past the real EOF, so several
            // of these can interleave, each describing a different slice.
            pread_chunk_fail("ERROR: unexpected EOF during index load "
                             "(chunk at offset %lld: %zu of %zu bytes)\n",
                             (long long)c->off, done, c->nbytes);
        }
        done += (size_t)r;
    }
    return NULL;
}

// Read `nbytes` at file offset `off` of `path` into `dst` using up to
// `nthreads` workers. Chunk 0 reads through `shared_fd` (already open on
// `path`); every other chunk gets its own descriptor from open_worker_fd.
void parallel_pread(int shared_fd, const char *path, void *dst, size_t nbytes,
                    off_t off, int nthreads)
{
    if (nbytes == 0) return;
    nthreads = fmi_pread_worker_count(nbytes, nthreads);

    std::vector<PreadChunk> chunks((size_t)nthreads);
    std::vector<pthread_t>  tids((size_t)nthreads);
    std::vector<bool>       spawned((size_t)nthreads, false);
    size_t base = nbytes / (size_t)nthreads, rem = nbytes % (size_t)nthreads, cum = 0;
    int n_sharing = 0;
    for (int i = 0; i < nthreads; i++) {
        size_t len = base + ((size_t)i < rem ? 1 : 0);
        int fd = i == 0 ? shared_fd : open_worker_fd(path, shared_fd);
        if (i > 0 && fd == shared_fd) n_sharing++;
        chunks[i] = PreadChunk{ fd, (char *)dst + cum, len, off + (off_t)cum };
        cum += len;
    }
    // Still a correct load, but one that can stall on S3 Mountpoint, so say why.
    if (n_sharing > 0)
        fprintf(stderr, "[W::%s] could not reopen %s for %d of %d index-load workers; "
                "they share one descriptor, which is slow on filesystems with "
                "per-handle read-ahead (e.g. S3 Mountpoint)\n",
                __func__, path, n_sharing, nthreads - 1);
    // Main thread takes chunk 0; spawn workers for the rest. If a spawn fails
    // (e.g. thread limit), fall back to reading that chunk inline so the load
    // still completes correctly, just less parallel.
    for (int i = 1; i < nthreads; i++) {
        if (pthread_create(&tids[i], NULL, pread_chunk_worker, &chunks[i]) == 0)
            spawned[i] = true;
        else
            pread_chunk_worker(&chunks[i]);
    }
    pread_chunk_worker(&chunks[0]);
    for (int i = 1; i < nthreads; i++)
        if (spawned[i]) pthread_join(tids[i], NULL);
    for (int i = 1; i < nthreads; i++)
        if (chunks[i].fd != shared_fd) close(chunks[i].fd);
}

}  // namespace

// Worker count for the index load: the caller's -t, capped (bandwidth bound
// past ~8), with an explicit BWA3_LOAD_THREADS override for tuning.
// Exposed (declared in FMI_search.h) so the env-parse validation is
// unit-testable without loading a real index.
int index_load_threads(int n_threads)
{
    int t = n_threads > 0 ? n_threads : 1;
    if (t > 8) t = 8;
    const char *e = getenv("BWA3_LOAD_THREADS");
    if (e != NULL && *e != '\0') {
        // Strict parse: atoi silently accepts garbage (abc -> 0, 0x8 -> 0) and
        // any huge value (100000 -> ~1280 pthread_creates). On malformed input
        // warn and keep the computed default -- ignore-and-continue like
        // BWA3_SMEM_LOCKSTEP_N (bwa3_lockstep_width_parse_env), NOT the hard
        // exit(1) that BWA_INDEX_THREADS's parse_ll takes. Clamp to a 64-thread
        // ceiling (the load is bandwidth-bound past ~8).
        char *end = NULL;
        errno = 0;
        long v = strtol(e, &end, 10);
        if (errno != 0 || end == e || *end != '\0' || v <= 0) {
            fprintf(stderr, "ERROR: BWA3_LOAD_THREADS=\"%s\" is not a positive integer; ignoring it\n", e);
        } else {
            if (v > 64) v = 64;   // ceiling; the load is bandwidth-bound past ~8
            t = (int)v;
        }
    }
    return t;
}

/* Declared in FMI_search.h; see there for the contract. ftello reports the
 * stream's logical position (buffered bytes included), so it is the right
 * offset to hand pread; the fseeko afterwards drops the now-stale buffer and
 * re-anchors the stream past the array we just read behind its back. */
void fmi_pread_from_stream(FILE *fp, const char *path, void *dst, size_t nbytes, int nthreads)
{
    off_t off = ftello(fp);
    if (off < 0) {
        fprintf(stderr, "ERROR: ftello failed during index load: %s\n", strerror(errno));
        exit(EXIT_FAILURE);
    }
    parallel_pread(fileno(fp), path, dst, nbytes, off, nthreads);
    if (fseeko(fp, off + (off_t)nbytes, SEEK_SET) != 0) {
        fprintf(stderr, "ERROR: fseeko failed during index load: %s\n", strerror(errno));
        exit(EXIT_FAILURE);
    }
}

/* Declared in FMI_search.h; see there for the contract.
 *
 * The on-disk layout is: header(48B: reference_seq_len + count[5], i.e.
 * 6*sizeof(int64_t)) + cp_occ + SA-sample arrays (sized by sa_compx) +
 * sentinel_index(8B) [+ sa_compx tag(8B), new format only]. Since the
 * SA-sample-array size depends on sa_compx, the tag's own file offset can't
 * be computed without already knowing sa_compx -- so instead we read the
 * file's last 8 bytes as a `candidate` sa_compx and verify it reproduces the
 * file's actual size when plugged back into the layout-size formula
 * (`off_sent_for`). If it doesn't -- including on a legacy, pre-this-feature
 * index, which carries no tag at all -- we fall back to `default_compx`.
 *
 * Why this is safe without a version field: a legacy index (no tag, always
 * built at the fixed compile-time SA_COMPX rate) has file size
 *
 *   legacy_size = off_sent_for(SA_COMPX) + 8            (sentinel only)
 *
 * while a false-positive match against some candidate `c` would require
 *
 *   legacy_size == off_sent_for(c) + 16                  (sentinel + tag)
 *
 * Each unit of SA-sample count changes off_sent_for by exactly 5 bytes
 * (1-byte int8_t + 4-byte uint32_t per sample), so equating the two reduces
 * to 5*Delta == 8 for some integer Delta = sa_sample_cnt(SA_COMPX) -
 * sa_sample_cnt(c). Since 5 does not divide 8, no integer Delta -- and hence
 * no candidate `c`, including c == SA_COMPX itself (Delta == 0, 0 != 8) --
 * can satisfy it. A legacy index can therefore never be misdetected as
 * carrying a valid tail. This is a property of the byte accounting, not of
 * a chosen constant, so it must not be "shored up" with a magic-number or
 * version field: doing so would change the on-disk format and break this
 * branch's byte-identity guarantee against stock indexes for no benefit. */
int64_t detect_sa_compx(int fd, int64_t file_size, int64_t ref_seq_len, int64_t default_compx)
{
    int64_t candidate = -1;
    if (file_size >= (int64_t)(2 * sizeof(int64_t))) {
        /* Robustly read the trailing 8-byte candidate tag. A read error or a
         * short/EOF read here must NOT silently fall back to default_compx: a
         * genuinely non-default index would then be mis-sized (wrong SA-sample
         * offsets, no diagnostic) -- the same failure the fstat() guards at
         * both call sites already prevent. Retry EINTR and treat any other
         * failure as a fatal load error, mirroring pread_chunk_worker() above.
         * Only a fully-read candidate that fails the layout test below
         * legitimately falls back to default_compx. */
        char   *dst  = (char *)&candidate;
        off_t   off  = (off_t)(file_size - (int64_t)sizeof(int64_t));
        size_t  done = 0;
        while (done < sizeof(int64_t)) {
            ssize_t r = pread(fd, dst + done, sizeof(int64_t) - done,
                              off + (off_t)done);
            if (r < 0) {
                if (errno == EINTR) continue;
                fprintf(stderr, "ERROR: pread failed reading SA-rate tag "
                                "during index load: %s\n", strerror(errno));
                exit(EXIT_FAILURE);
            }
            if (r == 0) {
                fprintf(stderr, "ERROR: unexpected EOF reading SA-rate tag "
                                "during index load\n");
                exit(EXIT_FAILURE);
            }
            done += (size_t)r;
        }
    }

    auto off_sent_for = [ref_seq_len](int64_t compx) -> int64_t {
        // Mirrors HDR_BYTES in fm_index_writer.cpp / BWA_BWT_2BIT_HEADER_BYTES
        // in bwa_shm.cpp: ref_seq_len + count[5].
        const int64_t hdr_bytes      = 6 * (int64_t)sizeof(int64_t);
        const int64_t cp_occ_cnt     = (ref_seq_len >> CP_SHIFT) + 1;
        const int64_t sa_sample_cnt  = (ref_seq_len >> compx) + 1;
        return hdr_bytes + cp_occ_cnt * (int64_t)sizeof(CP_OCC)
                         + sa_sample_cnt * (int64_t)sizeof(int8_t)
                         + sa_sample_cnt * (int64_t)sizeof(uint32_t);
    };

    if (candidate >= 0 && candidate <= CP_SHIFT &&
        off_sent_for(candidate) + 2 * (int64_t)sizeof(int64_t) == file_size) {
        return candidate;
    }
    return default_compx;
}

/* Build "<prefix><suffix>" into `out` (sized `outsz`); aborts on overflow.
 * Replaces the prior strcpy_s/strcat_s pattern for assembling FMI sidecar
 * paths from a user-supplied prefix. */
static void fmi_build_path(char *out, size_t outsz, const char *prefix, const char *suffix)
{
    int n = snprintf(out, outsz, "%s%s", prefix, suffix);
    if (n < 0 || (size_t)n >= outsz) {
        fprintf(stderr, "ERROR: FMI path too long for prefix '%s' (suffix '%s')\n",
                prefix, suffix);
        exit(EXIT_FAILURE);
    }
}

FMI_search::FMI_search(const char *fname)
{
    fprintf(stderr, "* Entering FMI_search\n");
    fmi_build_path(file_name, sizeof(file_name), fname, "");
    reference_seq_len = 0;
    sentinel_index = 0;
    sa_compx = SA_COMPX;
    sa_compx_mask = SA_COMPX_MASK;
    sa_ls_word = NULL;
    sa_ms_byte = NULL;
    cp_occ = NULL;
    shm_base = NULL;
    shm_len = 0;
    seed_jump = NULL;
    seed_jump_depth = 0;

    /* one_hot_mask_array is now a file-scope compile-time constant table (see
     * FMI_search.h); nothing to allocate or initialize here. */
}

FMI_search::~FMI_search()
{
    /* When attached from shm, cp_occ/sa_*_byte/sa_*_word alias mmap'd pages
     * owned by the shm segment, so we don't _mm_free them. We do munmap the
     * mapping itself — leaving it would leak VA in long-lived processes that
     * construct and destroy FMI_search repeatedly. */
    if (shm_base != NULL) {
        munmap(shm_base, shm_len);
        shm_base = NULL;
        shm_len  = 0;
    } else {
        if (sa_ms_byte) _mm_free(sa_ms_byte);
        if (sa_ls_word) _mm_free(sa_ls_word);
        if (cp_occ)     _mm_free(cp_occ);
    }
    /* Process-private on both paths (never part of the shm segment). */
    if (seed_jump) _mm_free(seed_jump);
}

int64_t FMI_search::cp_occ_size_bytes() const {
    return ((reference_seq_len >> CP_SHIFT) + 1) * (int64_t)sizeof(CP_OCC);
}

int64_t FMI_search::sa_sample_count() const {
    return (reference_seq_len >> sa_compx) + 1;
}

/* ===== Seed jump table (see FMI_search.h, seed_jump) =====
 *
 * Used by the bwtSeed walks (bsd_plan_jump) and the SMEM forward walks
 * (ls_init_slot, ls_advance_forward_table), so despite its name
 * BWA3_BWTSEED_JUMP also sets the table for the first- and second-pass SMEM
 * walks.
 * Depth: BWA3_BWTSEED_JUMP (0 = no table, so every walk steps from its first
 * base as before; 2..SEED_JUMP_DEPTH_MAX = the table's depth), default
 * SEED_JUMP_DEPTH_DEFAULT. Unset or empty silently means the default, as for
 * the other integer knobs; any other value is reported and the default used. An internal tuning knob: output is byte-identical
 * at every value. At depth D the table holds (4^(D+1) - 4) / 3 entries of 24
 * bytes: 2 MB at 8, 33.6 MB at 10, 134 MB at 11, 537 MB at 12. */
#define SEED_JUMP_DEPTH_DEFAULT 10
#define SEED_JUMP_DEPTH_MAX     12
/* Levels 1..SEED_JUMP_PREFIX are built serially; deeper levels are built per
 * subtree of a SEED_JUMP_PREFIX-base prefix (64 of them), one subtree at a time
 * per worker, so no level needs a barrier. */
#define SEED_JUMP_PREFIX        3
/* Parents ahead whose two cp_occ blocks a builder prefetches. */
#define SEED_JUMP_PF            16

static int seed_jump_env_depth()
{
    const char *e = getenv("BWA3_BWTSEED_JUMP");
    if (e == NULL || *e == '\0') return SEED_JUMP_DEPTH_DEFAULT;
    /* Digits only: strtol alone would also take leading space and a sign. */
    char *end = NULL;
    errno = 0;
    long v = strtol(e, &end, 10);
    if (*e < '0' || *e > '9' || errno != 0 || end == e || *end != '\0'
        || v == 1 || v > SEED_JUMP_DEPTH_MAX) {
        fprintf(stderr, "ERROR: BWA3_BWTSEED_JUMP=\"%s\" is not 0 or an integer in [2, %d]; "
                "ignoring it and using the default depth %d\n", e, SEED_JUMP_DEPTH_MAX,
                SEED_JUMP_DEPTH_DEFAULT);
        return SEED_JUMP_DEPTH_DEFAULT;
    }
    return (int)v;
}

/* Level L's children of the level L - 1 strings [first, first + n): child c
 * of parent q is q * 4 + c, computed exactly as bsd_advance_step steps (swap
 * k <-> l, backwardExt with 3 - c, swap back). */
void FMI_search::seed_jump_children(int L, int64_t first, int64_t n)
{
    const SeedJumpEntry *par = seed_jump + seed_jump_level_off(L - 1) + first;
    SeedJumpEntry *ch = seed_jump + seed_jump_level_off(L) + first * 4;
    for (int64_t i = 0; i < n; i++) {
#ifdef ENABLE_PREFETCH
        if (i + SEED_JUMP_PF < n) {
            const SeedJumpEntry &q = par[i + SEED_JUMP_PF];
            _mm_prefetch((const char *)(&cp_occ[(q.l >> CP_SHIFT)]), _MM_HINT_T0);
            _mm_prefetch((const char *)(&cp_occ[((q.l + q.s) >> CP_SHIFT)]), _MM_HINT_T0);
        }
#endif
        SMEM smem_;
        memset(&smem_, 0, sizeof(smem_));
        smem_.k = par[i].l;
        smem_.l = par[i].k;
        smem_.s = par[i].s;
        /* All four before any store, so the compiler can share the two
         * blocks' loads and popcounts among them (a store to ch[] may alias
         * cp_occ). */
        SMEM r[4];
        for (int c = 0; c < 4; c++) r[c] = backwardExt(smem_, (uint8_t)(3 - c));
        for (int c = 0; c < 4; c++) {
            ch[i * 4 + c].k = r[c].l;
            ch[i * 4 + c].l = r[c].k;
            ch[i * 4 + c].s = r[c].s;
        }
    }
}

namespace {
struct SeedJumpBuild {
    FMI_search          *fmi;
    std::atomic<int64_t> next;
    int64_t              n_prefix;
};
}  // namespace

void *FMI_search::seed_jump_worker(void *arg)
{
    SeedJumpBuild *b = (SeedJumpBuild *)arg;
    for (;;) {
        const int64_t p = b->next.fetch_add(1, std::memory_order_relaxed);
        if (p >= b->n_prefix) break;
        /* Every level below prefix p, level by level. */
        for (int L = SEED_JUMP_PREFIX + 1; L <= b->fmi->seed_jump_depth; L++) {
            const int64_t n = (int64_t)1 << (2 * (L - 1 - SEED_JUMP_PREFIX));
            b->fmi->seed_jump_children(L, p * n, n);
        }
    }
    return NULL;
}

void FMI_search::build_seed_jump_table(int n_threads)
{
    if (seed_jump != NULL) { _mm_free(seed_jump); seed_jump = NULL; }
    seed_jump_depth = 0;
    const int D = seed_jump_env_depth();
    if (D < 2) return;

    const int64_t total = seed_jump_level_off(D + 1);
    const size_t bytes = (size_t)total * sizeof(SeedJumpEntry);
    seed_jump = (SeedJumpEntry *)_mm_malloc(bytes, 64);
    if (seed_jump == NULL) {
        /* The table only saves steps, so run without it rather than fail. */
        fprintf(stderr, "ERROR: could not allocate the %.1f MB seed jump table "
                "(BWA3_BWTSEED_JUMP depth %d); seeding without it\n", bytes / 1e6, D);
        return;
    }
    bwamem_madv_hugepage(seed_jump, bytes);
    seed_jump_depth = D;

    for (int a = 0; a < 4; a++) {   // bsd_init_slot's 1-base interval
        seed_jump[a].k = count[a];
        seed_jump[a].l = count[3 - a];
        seed_jump[a].s = count[a + 1] - count[a];
    }
    const int P = D < SEED_JUMP_PREFIX ? D : SEED_JUMP_PREFIX;
    for (int L = 2; L <= P; L++)
        seed_jump_children(L, 0, (int64_t)1 << (2 * (L - 1)));
    if (D == P) return;

    SeedJumpBuild b;
    b.fmi = this;
    b.next.store(0, std::memory_order_relaxed);
    b.n_prefix = (int64_t)1 << (2 * P);
    int nt = n_threads > 0 ? n_threads : 1;
    if (nt > b.n_prefix) nt = (int)b.n_prefix;
    std::vector<pthread_t> tids((size_t)nt);
    std::vector<bool>      spawned((size_t)nt, false);
    /* A failed spawn only makes the build less parallel: the main thread's
     * loop takes every prefix no worker took. */
    for (int i = 1; i < nt; i++)
        spawned[i] = pthread_create(&tids[i], NULL, seed_jump_worker, &b) == 0;
    seed_jump_worker(&b);
    for (int i = 1; i < nt; i++)
        if (spawned[i]) pthread_join(tids[i], NULL);
}

void FMI_search::load_index_from_shm(uint8_t *base, size_t len)
{
    if (base == NULL) {
        fprintf(stderr, "ERROR! load_index_from_shm called with NULL base\n");
        exit(EXIT_FAILURE);
    }
    uint64_t off = 0, sz = 0;

    if (bwa_shm_section_find(base, BWA_SHM_SEC_FMI_SCALARS, &off, &sz) != 0
        || sz != BWA_SHM_FMI_SCALARS_BYTES) {
        fprintf(stderr, "ERROR! shm segment missing or malformed FMI_SCALARS\n");
        exit(EXIT_FAILURE);
    }
    memcpy(&reference_seq_len, base + off,                                    sizeof(int64_t));
    memcpy(count,              base + off + sizeof(int64_t),                  sizeof(int64_t) * 5);
    memcpy(&sentinel_index,    base + off + sizeof(int64_t) * 6,              sizeof(int64_t));
    /* sa_compx: the SA sample-rate shift the staged index was built with.
     * IMPORTANT FIX: this used to stay at the constructor's compile-time
     * default (SA_COMPX) because nothing here overrode it, so a non-default
     * `-u`-built index attached via shm was mis-sized (sa_sample_count() and
     * every sa_ms_byte/sa_ls_word index below used the wrong shift). Reading
     * it from the packed scalars (see BWA_SHM_FMI_SCALARS_BYTES / the PACK
     * side in bwa_shm.cpp) makes shm-attach agree with the disk loader's
     * tail-detected rate. */
    memcpy(&sa_compx,          base + off + sizeof(int64_t) * 7,              sizeof(int64_t));

    /* Validate scalars before we use reference_seq_len in cp_occ_size_bytes()
     * and the SA size accessors. Bounds match the disk path's asserts in
     * load_index() and the writer-side checks in bwa_shm_compute(). A corrupt
     * segment would otherwise drive negative or overflowing section sizes. */
    if (reference_seq_len <= 0 || reference_seq_len > 0x7fffffffffLL) {
        fprintf(stderr,
            "ERROR! shm FMI_SCALARS: reference_seq_len=%lld out of bounds\n",
            (long long)reference_seq_len);
        exit(EXIT_FAILURE);
    }
    /* count[] is +1-adjusted by bwa_shm_compute; range [1, ref_seq_len+1], and
     * as the cumulative C[] array it must be non-decreasing. These bounds match
     * the disk path's post-adjustment asserts in load_index(): a count[i] == 0
     * drives smem.k out of cp_occ on the first backwardExt, and a non-monotone
     * pair yields a negative smem.s (interval width) during seed extraction, so
     * a corrupt segment must fail closed here rather than attach with an invalid
     * C array. */
    for (int i = 0; i < 5; ++i) {
        if (count[i] < 1 || count[i] > reference_seq_len + 1) {
            fprintf(stderr,
                "ERROR! shm FMI_SCALARS: count[%d]=%lld out of bounds (ref_seq_len=%lld)\n",
                i, (long long)count[i], (long long)reference_seq_len);
            exit(EXIT_FAILURE);
        }
        if (i > 0 && count[i] < count[i - 1]) {
            fprintf(stderr,
                "ERROR! shm FMI_SCALARS: count[%d]=%lld < count[%d]=%lld "
                "(not monotonically non-decreasing)\n",
                i, (long long)count[i], i - 1, (long long)count[i - 1]);
            exit(EXIT_FAILURE);
        }
    }
    if (sentinel_index < 0 || sentinel_index >= reference_seq_len) {
        fprintf(stderr,
            "ERROR! shm FMI_SCALARS: sentinel_index=%lld out of bounds (ref_seq_len=%lld)\n",
            (long long)sentinel_index, (long long)reference_seq_len);
        exit(EXIT_FAILURE);
    }
    /* [0,6]: same range the -u CLI validates and write_fm_index_streaming
     * enforces (the sample period 1<<sa_compx must divide CP_BLOCK_SIZE=64). */
    if (sa_compx < 0 || sa_compx > 6) {
        fprintf(stderr,
            "ERROR! shm FMI_SCALARS: sa_compx=%lld out of bounds\n",
            (long long)sa_compx);
        exit(EXIT_FAILURE);
    }
    /* Only compute the mask once sa_compx is known to be in [0,6]; shifting by
     * an out-of-range value (e.g. a corrupt -1 or 64) would be UB. */
    sa_compx_mask = (1LL << sa_compx) - 1;

    if (bwa_shm_section_find(base, BWA_SHM_SEC_FMI_CP_OCC, &off, &sz) != 0
        || (int64_t)sz != cp_occ_size_bytes()) {
        fprintf(stderr, "ERROR! shm segment missing or sized FMI_CP_OCC\n");
        exit(EXIT_FAILURE);
    }
    cp_occ = (CP_OCC *)(base + off);

    if (bwa_shm_section_find(base, BWA_SHM_SEC_FMI_SA_MS, &off, &sz) != 0
        || (int64_t)sz != sa_ms_byte_size_bytes()) {
        fprintf(stderr, "ERROR! shm segment missing or sized FMI_SA_MS\n");
        exit(EXIT_FAILURE);
    }
    sa_ms_byte = (int8_t *)(base + off);

    if (bwa_shm_section_find(base, BWA_SHM_SEC_FMI_SA_LS, &off, &sz) != 0
        || (int64_t)sz != sa_ls_word_size_bytes()) {
        fprintf(stderr, "ERROR! shm segment missing or sized FMI_SA_LS\n");
        exit(EXIT_FAILURE);
    }
    sa_ls_word = (uint32_t *)(base + off);

    shm_base = base;
    shm_len  = len;

    fprintf(stderr, "* FMI attached from shm: ref_seq_len=%ld sentinel_index=%ld\n",
            (long)reference_seq_len, (long)sentinel_index);
}

void FMI_search::densify_sa_into(const CP_OCC *cp_occ_src, const int64_t count_src[5],
                                 int64_t ref_seq_len, int64_t sentinel_idx,
                                 const int8_t *src_ms, const uint32_t *src_ls,
                                 int64_t src_compx,
                                 int8_t *dst_ms, uint32_t *dst_ls, int64_t dst_compx,
                                 int n_threads)
{
    /* Stride validity: the resolver reads a sample at (row >> sa_compx), so a
     * denser destination (smaller shift) is required for this to add samples,
     * and both shifts must be in the writer-enforced [0,6] range (the sample
     * period 1<<compx must divide CP_BLOCK_SIZE=64). */
    if (!(dst_compx >= 0 && dst_compx < src_compx && src_compx <= 6)) {
        fprintf(stderr,
            "ERROR! densify_sa_into: invalid strides (dst=%lld src=%lld)\n",
            (long long)dst_compx, (long long)src_compx);
        exit(EXIT_FAILURE);
    }

    /* This borrows the object's own index-state members (cp_occ/sa_*) and clears
     * them before return, so it must only be called on a fresh FMI_search that
     * owns no loaded buffers -- otherwise those buffers would leak. Both current
     * callers (bwa_shm_pack_into, main_resa) pass a freshly-constructed object;
     * guard the precondition so a future caller on a loaded index fails loudly. */
    xassert(cp_occ == NULL && sa_ms_byte == NULL && sa_ls_word == NULL,
            "densify_sa_into called on an FMI_search with loaded index buffers");

    /* Validate count_src (the +1-adjusted cumulative C[]) before it drives the
     * LF-walk: get_sa_entry_compressed() computes sp = count[b] + occ and indexes
     * cp_occ[sp >> CP_SHIFT] with no bound check, so a non-monotonic or
     * out-of-range C[] would read past cp_occ. The shm and re-sa callers both
     * range-check on their own, but only re-sa checks monotonicity, so guard it
     * here at the shared chokepoint. Valid adjusted range is [1, ref_seq_len+1]. */
    for (int i = 0; i < 5; ++i) {
        if (count_src[i] < 1 || count_src[i] > ref_seq_len + 1 ||
            (i > 0 && count_src[i] < count_src[i - 1])) {
            fprintf(stderr,
                "ERROR! densify_sa_into: invalid count[] (non-monotonic or out of "
                "[1, %lld])\n", (long long)(ref_seq_len + 1));
            exit(EXIT_FAILURE);
        }
    }

    /* Borrow the source index state so get_sa_entry_compressed() resolves any
     * BWT row against the SOURCE samples (its walk terminates on the source
     * mask). These members are cleared before return; nothing here is owned. */
    reference_seq_len = ref_seq_len;
    sentinel_index    = sentinel_idx;
    memcpy(count, count_src, sizeof(int64_t) * 5);
    cp_occ        = const_cast<CP_OCC *>(cp_occ_src);
    sa_ms_byte    = const_cast<int8_t *>(src_ms);
    sa_ls_word    = const_cast<uint32_t *>(src_ls);
    sa_compx      = src_compx;
    sa_compx_mask = (1LL << src_compx) - 1;

    const int64_t src_mask = (1LL << src_compx) - 1;
    const int64_t dst_count = (ref_seq_len >> dst_compx) + 1;
    const int T = (n_threads > 1) ? n_threads : 1;

    /* For each destination sample slot j, its BWT row is row = j << dst_compx.
     * If that row is already a SOURCE-sampled row (row & src_mask == 0), copy
     * the stored value verbatim — the source and destination tables must agree
     * bit-for-bit on shared rows. Otherwise recover it by the resolver's own
     * LF-walk against the source stride.
     *
     * Parallel-safe: each iteration reads only shared, immutable index state
     * (cp_occ / count / sa_*_byte / masks via get_sa_entry_compressed) and
     * writes disjoint dst_ms[j]/dst_ls[j] slots — no shared mutable state, so
     * no locking. schedule(dynamic) balances the uneven per-row walk lengths
     * (already-sampled rows are O(1); added rows walk up to 2^src_compx-1
     * LF steps). Output is independent of thread count and scheduling. */
#ifdef _OPENMP
    #pragma omp parallel for num_threads(T) schedule(dynamic, 4096)
#else
    (void)T;
#endif
    for (int64_t j = 0; j < dst_count; ++j) {
        const int64_t row = j << dst_compx;
        int64_t v;
        if ((row & src_mask) == 0) {
            const int64_t s = row >> src_compx;
            v = ((int64_t)src_ms[s] << 32) + (int64_t)src_ls[s];
        } else {
            v = get_sa_entry_compressed(row, /*tid=*/0);
        }
        dst_ms[j] = (int8_t)((v >> 32) & 0xff);
        dst_ls[j] = (uint32_t)(v & 0xffffffffULL);
    }

    /* Drop the borrowed pointers so ~FMI_search() frees nothing (they alias
     * the caller's shm buffer / heap temporaries). */
    cp_occ     = NULL;
    sa_ms_byte = NULL;
    sa_ls_word = NULL;
}

int FMI_search::build_index(bool emit_unpacked_ref, int sa_compx) {

    char *prefix = file_name;

    // Read single-strand length from .ann to compute doubled pac_len.
    char ann_path[PATH_MAX];
    fmi_build_path(ann_path, sizeof(ann_path), prefix, ".ann");
    FILE* fann = fopen(ann_path, "r");
    if (fann == NULL) {
        fprintf(stderr, "ERROR: cannot open '%s'\n", ann_path);
        return 1;
    }
    int64_t l_pac = 0;
    int n_seqs = 0, seed = 0;
    if (fscanf(fann, "%" SCNd64 " %d %d", &l_pac, &n_seqs, &seed) != 3) {
        fprintf(stderr, "ERROR: malformed '%s'\n", ann_path);
        fclose(fann);
        return 1;
    }
    fclose(fann);
    // Defensive: a non-positive l_pac (corrupt or zero-length .ann) would
    // pass a bad pac_len into libsais and only fail much later. Catch it
    // up front with an actionable message, mirroring bntseq_restore_core.
    if (l_pac <= 0 || n_seqs < 0) {
        fprintf(stderr, "ERROR: malformed '%s' (l_pac=%" PRId64 ", n_seqs=%d)\n",
                ann_path, l_pac, n_seqs);
        return 1;
    }
    int64_t pac_len = 2 * l_pac;

    auto parse_ll = [](const char* s, const char* name,
                       long long max_val = LLONG_MAX) -> long long {
        char* end = nullptr;
        errno = 0;
        long long v = strtoll(s, &end, 10);
        if (errno || end == s || *end != '\0' || v <= 0 || v > max_val) {
            fprintf(stderr, "ERROR: invalid %s='%s' (expected positive integer)\n",
                    name, s);
            exit(1);
        }
        return v;
    };
    LibsaisBuildOpts opts;
    if (const char* th = getenv("BWA_INDEX_THREADS"))
        opts.num_threads      = (int)parse_ll(th, "BWA_INDEX_THREADS", INT_MAX);
    if (const char* mm = getenv("BWA_INDEX_MAX_MEMORY"))
        opts.max_memory_bytes = parse_ll(mm, "BWA_INDEX_MAX_MEMORY");
    if (const char* td = getenv("BWA_INDEX_TMPDIR"))
        opts.tmpdir           = td;
    if (const char* mu = getenv("BWA_INDEX_MAX_MEMORY_USER"))
        opts.max_memory_user_specified = (mu[0] == '1');
    opts.emit_unpacked_ref = emit_unpacked_ref;
    opts.sa_compx          = sa_compx;
    return libsais_build_fm_index(prefix, pac_len, opts);
}

void FMI_search::load_index(bool load_pac, int n_threads, bool build_seed_jump)
{
    /* Try the staged shm segment first. On hit, both the FMI internals
     * (cp_occ / sa_*) and the BNS+PAC are attached as views into the
     * mapping; the segment lifetime belongs to the shm-loading process. */
    {
        size_t   shm_attach_len = 0;
        uint8_t *shm_base_local = bwa_shm_attach(file_name, &shm_attach_len);
        if (shm_base_local != NULL) {
            load_index_from_shm(shm_base_local, shm_attach_len);
            /* D3 --meth seed segment is staged bns_only: no PAC section, and
             * idx->pac stays NULL (extension uses meth_orig_pac). */
            bwa_idx_load_ele_from_shm(shm_base_local, shm_attach_len, load_pac);
            /* The seed jump table is process-private (not part of the staged
             * segment), so an attached index builds its own. */
            if (build_seed_jump) build_seed_jump_table(index_load_threads(n_threads));
            fprintf(stderr, "* FMI+BNS%s attached from shm; "
                    "skipping disk load.\n", load_pac ? "+PAC" : "");
            return;
        }
    }

    // Running total of index bytes allocated so far in this function.
    // Passed to assert_not_null so a failed allocation reports both the
    // attempted size and how much we'd already committed before failing.
    int64_t index_alloc = 0;

    // Worker count for the big-array reads below (read once, out of any loop).
    const int load_nt = index_load_threads(n_threads);

    char *ref_file_name = file_name;
    //beCalls = 0;
    char cp_file_name[PATH_MAX];
    fmi_build_path(cp_file_name, sizeof(cp_file_name), ref_file_name, CP_FILENAME_SUFFIX);

    // Read the BWT and FM index of the reference sequence
    FILE *cpstream = NULL;
    cpstream = fopen(cp_file_name,"rb");
    if (cpstream == NULL)
    {
        fprintf(stderr, "ERROR! Unable to open the file: %s\n", cp_file_name);
        exit(EXIT_FAILURE);
    }
    else
    {
        fprintf(stderr, "* Index file found. Loading index from %s\n", cp_file_name);
    }

#if defined(__linux__)
    // Kick readahead for the whole index so the parallel preads below hit warm
    // pages on a cold cache (no-op benefit when already cached).
    posix_fadvise(fileno(cpstream), 0, 0, POSIX_FADV_WILLNEED);
#endif

    err_fread_noeof(&reference_seq_len, sizeof(int64_t), 1, cpstream);
    // Fail closed on a truncated/corrupt index, mirroring load_index_from_shm's
    // scalar bounds checks (a bad reference_seq_len/count[]/sentinel_index drives
    // smem.k past cp_occ on the first backwardExt -> wrong seeds, not a crash).
    // xassert survives an EXTRA_CXXFLAGS=-DNDEBUG opt-out build; a plain assert would not.
    xassert(reference_seq_len > 0 && reference_seq_len <= 0x7fffffffffLL,
            "FMI index: reference_seq_len out of bounds (corrupt .bwt.2bit.64?)");

    fprintf(stderr, "* Reference seq len for bi-index = %lld\n", (long long)reference_seq_len);

    // Peek the trailing sa_compx field (if present) before sizing the SA
    // sample arrays below, via the shared detect_sa_compx() helper (see its
    // doc comment in FMI_search.h and implementation comment above for why
    // the tail-detection heuristic is safe). pread doesn't disturb
    // cpstream's buffered read position, so this doesn't affect the
    // sequential reads that follow.
    {
        struct stat st;
        /* A failed fstat() must not silently fall back to file_size=0: that
         * would make detect_sa_compx() treat a genuinely non-default-rate
         * index as legacy (rate 3), mis-sizing the SA-sample arrays below
         * and producing wrong seed coordinates with no diagnostic. Fail the
         * load instead -- consistent with the ftello/fseeko error handling
         * just above in fmi_pread_from_stream(). */
        if (fstat(fileno(cpstream), &st) != 0) {
            fprintf(stderr, "ERROR: fstat failed during index load: %s\n", strerror(errno));
            exit(EXIT_FAILURE);
        }
        sa_compx = detect_sa_compx(fileno(cpstream), (int64_t)st.st_size, reference_seq_len, SA_COMPX);
        sa_compx_mask = (1LL << sa_compx) - 1;
    }

    // create checkpointed occ
    int64_t cp_occ_size = (reference_seq_len >> CP_SHIFT) + 1;
    cp_occ = NULL;

    err_fread_noeof(&count[0], sizeof(int64_t), 5, cpstream);
    int64_t cp_occ_bytes = cp_occ_size * sizeof(CP_OCC);
    cp_occ = (CP_OCC *)_mm_malloc(cp_occ_bytes, 64);
    index_alloc += cp_occ_bytes;
    assert_not_null(cp_occ, cp_occ_bytes, index_alloc);
    bwamem_madv_hugepage(cp_occ, cp_occ_bytes);

    fmi_pread_from_stream(cpstream, cp_file_name, cp_occ, (size_t)cp_occ_bytes, load_nt);
    int64_t ii = 0;
    // Validate the RAW on-disk count[] before the +1 adjustment below. Doing the
    // bounds/monotonicity checks on the raw values (rather than after +1) is what
    // makes this overflow-safe: a corrupt count[ii] == INT64_MAX would trigger
    // signed-overflow UB in `count[ii] + 1` before any post-adjustment check
    // could fail closed. On disk count[] is the cumulative C[] array of the
    // reference, so a valid raw entry lies in [0, reference_seq_len]; the shm
    // producer (bwa_shm_compute) applies the same +1 and load_index_from_shm
    // then validates the adjusted [1, ref_seq_len+1] range.
    for (ii = 0; ii < 5; ii++) {
        // A raw count[ii] outside [0, ref_seq_len] (e.g. a disk -1, or a value
        // exceeding the reference length) drives smem.k out of cp_occ on the
        // first backwardExt once adjusted.
        xassert(count[ii] >= 0 && count[ii] <= reference_seq_len,
                "FMI index: count[] out of bounds (corrupt .bwt.2bit.64?)");
        // count[] is the cumulative C[] array, so it must be non-decreasing; a
        // non-monotone pair yields a negative smem.s (interval width) downstream.
        // The +1 is applied uniformly, so raw monotonicity == adjusted monotonicity.
        xassert(ii == 0 || count[ii] >= count[ii - 1],
                "FMI index: count[] not monotonically non-decreasing (corrupt .bwt.2bit.64?)");
    }
    for(ii = 0; ii < 5; ii++)// update read count structure
    {
        count[ii] = count[ii] + 1;
    }
    // Validate count[] (post +1-adjustment, range [1, ref_seq_len+1]) as
    // load_index_from_shm does; a count[a] exceeding the reference length drives
    // smem.k out of cp_occ on the first backwardExt.
    for (ii = 0; ii < 5; ii++) {
        // Lower bound is 1, not 0: the +1 adjustment above is applied before
        // this check (load_index_from_shm instead validates the raw disk value
        // in [0, ref_seq_len] and adjusts after), so a disk count[ii] == -1
        // becomes 0 here and would slip past a >= 0 test -- leaving the disk
        // loader open to an invalid FM interval that the shm producer rejects.
        xassert(count[ii] >= 1 && count[ii] <= reference_seq_len + 1,
                "FMI index: count[] out of bounds (corrupt .bwt.2bit.64?)");
        // count[] is the cumulative C[] array, so it must be non-decreasing; a
        // non-monotone pair yields a negative smem.s (interval width) downstream.
        xassert(ii == 0 || count[ii] >= count[ii - 1],
                "FMI index: count[] not monotonically non-decreasing (corrupt .bwt.2bit.64?)");
    }

    #if SA_COMPRESSION

    int64_t reference_seq_len_ = (reference_seq_len >> sa_compx) + 1;
    int64_t sa_ms_bytes = reference_seq_len_ * sizeof(int8_t);
    int64_t sa_ls_bytes = reference_seq_len_ * sizeof(uint32_t);
    sa_ms_byte = (int8_t *)_mm_malloc(sa_ms_bytes, 64);
    index_alloc += sa_ms_bytes;
    assert_not_null(sa_ms_byte, sa_ms_bytes, index_alloc);
    sa_ls_word = (uint32_t *)_mm_malloc(sa_ls_bytes, 64);
    index_alloc += sa_ls_bytes;
    assert_not_null(sa_ls_word, sa_ls_bytes, index_alloc);
    bwamem_madv_hugepage(sa_ms_byte, sa_ms_bytes);
    bwamem_madv_hugepage(sa_ls_word, sa_ls_bytes);
    fmi_pread_from_stream(cpstream, cp_file_name, sa_ms_byte, (size_t)sa_ms_bytes, load_nt);
    fmi_pread_from_stream(cpstream, cp_file_name, sa_ls_word, (size_t)sa_ls_bytes, load_nt);

    #else

    int64_t sa_ms_bytes = reference_seq_len * sizeof(int8_t);
    int64_t sa_ls_bytes = reference_seq_len * sizeof(uint32_t);
    sa_ms_byte = (int8_t *)_mm_malloc(sa_ms_bytes, 64);
    index_alloc += sa_ms_bytes;
    assert_not_null(sa_ms_byte, sa_ms_bytes, index_alloc);
    sa_ls_word = (uint32_t *)_mm_malloc(sa_ls_bytes, 64);
    index_alloc += sa_ls_bytes;
    assert_not_null(sa_ls_word, sa_ls_bytes, index_alloc);
    bwamem_madv_hugepage(sa_ms_byte, sa_ms_bytes);
    bwamem_madv_hugepage(sa_ls_word, sa_ls_bytes);
    fmi_pread_from_stream(cpstream, cp_file_name, sa_ms_byte, (size_t)sa_ms_bytes, load_nt);
    fmi_pread_from_stream(cpstream, cp_file_name, sa_ls_word, (size_t)sa_ls_bytes, load_nt);

    #endif

    sentinel_index = -1;
    #if SA_COMPRESSION
    err_fread_noeof(&sentinel_index, sizeof(int64_t), 1, cpstream);
    xassert(sentinel_index >= 0 && sentinel_index < reference_seq_len,
            "FMI index: sentinel_index out of bounds (corrupt .bwt.2bit.64?)");
    fprintf(stderr, "* sentinel-index: %lld\n", (long long)sentinel_index);
    #endif
    fclose(cpstream);

    int64_t x;
    #if !SA_COMPRESSION
    for(x = 0; x < reference_seq_len; x++)
    {
        // fprintf(stderr, "x: %ld\n", x);
        #if SA_COMPRESSION
        if(get_sa_entry_compressed(x) == 0) {
            sentinel_index = x;
            break;
        }
        #else
        if(get_sa_entry(x) == 0) {
            sentinel_index = x;
            break;
        }
        #endif
    }
    fprintf(stderr, "\nsentinel_index: %lld\n", (long long)x);
    #endif

    // cp_occ, count[] and sentinel_index are in: build the seed jump table.
    if (build_seed_jump) build_seed_jump_table(load_nt);

    fprintf(stderr, "* Count:\n");
    for(x = 0; x < 5; x++)
    {
        fprintf(stderr, "%lld,\t%lld\n", (long long)x, (long long)count[x]);
    }
    fprintf(stderr, "\n");  

    fprintf(stderr, "* Reading other elements of the index from files %s\n",
            ref_file_name);
    /* D3 --meth: BNS only for the seed index (skip the ~1.6 GB seed pac). The
     * seed bns drives the seed->original remap; extension uses meth_orig_pac. */
    /* Pass the worker count already resolved for this load (line above), so the
     * .pac slurp reuses it rather than re-parsing BWA3_LOAD_THREADS a second
     * time within the same load. */
    bwa_idx_load_ele(ref_file_name, load_pac ? BWA_IDX_ALL : BWA_IDX_BNS, load_nt);

    fprintf(stderr, "* Done reading Index!!\n");
}

void FMI_search::getSMEMsOnePosOneThread(uint8_t *enc_qdb,
                                         int32_t *query_pos_array,
                                         int32_t *min_intv_array,
                                         int32_t *rid_array,
                                         int32_t numReads,
                                         int32_t batch_size,
                                         const bseq1_t *seq_,
                                         int32_t *query_cum_len_ar,
                                         int32_t max_readlength,
                                         int32_t minSeedLen,
                                         SMEM *matchArray,
                                         int64_t *__numTotalSmem)
{
    int64_t numTotalSmem = *__numTotalSmem;
    // Heap-allocate to avoid stack overflow on long reads (PacBio HiFi /
    // ONT 1 Mbp+); mirrors the lockstep path's max_readlength sizing.
    int64_t prevArray_bytes = (int64_t)max_readlength * sizeof(SMEM);
    SMEM *prevArray = (SMEM *)_mm_malloc((size_t)prevArray_bytes, 64);
    assert_not_null(prevArray, prevArray_bytes, prevArray_bytes);

    // Hoist the FM-index cumulative-count table into a local for the batch (as
    // the bwa-mem2 reference does). count[] is a lifetime-constant member never
    // mutated during seeding, so this is a pure cache of identical values: it
    // trades the per-seed `this->count[..]` member loads for stack reads. Byte-
    // identical to reading the member directly.
    const int64_t counts[5] = {count[0], count[1], count[2], count[3], count[4]};

    uint32_t i;
    // Perform SMEM for original reads
    for(i = 0; i < numReads; i++)
    {
        int x = query_pos_array[i];
        int32_t rid = rid_array[i];
        int next_x = x + 1;

        int readlength = seq_[rid].l_seq;
        int offset = query_cum_len_ar[rid];
        // uint8_t a = enc_qdb[rid * readlength + x];
        uint8_t a = enc_qdb[offset + x];

        if(a < 4)
        {
            SMEM smem;
            smem.info = 0;   // carries no SA tag (ls_s1_sa_tag)
            smem.rid = rid;
            smem.m = x;
            smem.n = x;
            smem.k = counts[a];
            smem.l = counts[3 - a];
            smem.s = counts[a+1] - counts[a];
            int numPrev = 0;
            
            int j;
            for(j = x + 1; j < readlength; j++)
            {
                // a = enc_qdb[rid * readlength + j];
                a = enc_qdb[offset + j];
                next_x = j + 1;
                if(a < 4)
                {
                    SMEM smem_ = smem;

                    // Forward extension is backward extension with the BWT of reverse complement:
                    // swap k <-> l, backwardExt, swap back. `l` is LIVE here — after the swap it is
                    // the k that drives the next forward step, and the full backwardExt's
                    // l-cumulation is what produces it. So backwardExt_konly (which drops the
                    // backward-phase-dead l-chain) does NOT apply to the forward pass; the full
                    // backwardExt is required. konly is backward-only by construction.
                    smem_.k = smem.l;
                    smem_.l = smem.k;
                    SMEM newSmem_ = backwardExt(smem_, 3 - a);
                    //SMEM newSmem_ = forwardExt(smem_, 3 - a);
                    SMEM newSmem = newSmem_;
                    newSmem.k = newSmem_.l;
                    newSmem.l = newSmem_.k;
                    newSmem.n = j;

                    int32_t s_neq_mask = newSmem.s != smem.s;

                    prevArray[numPrev] = smem;
                    numPrev += s_neq_mask;
                    if(newSmem.s < min_intv_array[i])
                    {
                        next_x = j;
                        break;
                    }
                    smem = newSmem;
#ifdef ENABLE_PREFETCH
                    /* Next iteration swaps k<->l then backwardExt reads
                     * sp = smem.l, ep = smem.l + smem.s (see ls_advance_forward_step
                     * for the full derivation). Prefetching smem.k targeted an
                     * address never touched; ep was never prefetched. Pure hint. */
                    _mm_prefetch((const char *)(&cp_occ[(smem.l) >> CP_SHIFT]), _MM_HINT_T0);
                    _mm_prefetch((const char *)(&cp_occ[(smem.l + smem.s) >> CP_SHIFT]), _MM_HINT_T0);
#endif
                }
                else
                {
                    break;
                }
            }
            if(smem.s >= min_intv_array[i])
            {

                prevArray[numPrev] = smem;
                numPrev++;
            }

            SMEM *prev;
            prev = prevArray;

            int p;
            for(p = 0; p < (numPrev/2); p++)
            {
                SMEM temp = prev[p];
                prev[p] = prev[numPrev - p - 1];
                prev[numPrev - p - 1] = temp;
            }

            // Backward search.
            //
            // This scalar path deliberately keeps the FULL backwardExt (which
            // also maintains the reverse-complement start `l`), NOT the
            // backwardExt_konly used by the lockstep production path
            // (ls_advance_backward_step). Two reasons, both intentional — do not
            // "unify" this to konly:
            //   1. Oracle independence. This function is the reference the
            //      lockstep path is validated against in
            //      test/smem_lockstep_parity_test.cpp; that test cross-checks
            //      konly's k/s against this full implementation's k/s. Making
            //      both use konly would remove the independent check.
            //   2. It is compiled out of production at the shipping
            //      SMEM_LOCKSTEP_N > 1 (getSMEMsAllPosOneThread dispatches to the
            //      lockstep driver), so the extra l-chain work here is off every
            //      shipped path — there is nothing to optimize.
            // (A SMEM_LOCKSTEP_N <= 1 build would run this in production and pay
            // the full-backwardExt cost, but that config is unsupported and the
            // parity test #errors on it.)
            int cur_j = readlength;
            for(j = x - 1; j >= 0; j--)
            {
                int numCurr = 0;
                int curr_s = -1;
                // a = enc_qdb[rid * readlength + j];
                a = enc_qdb[offset + j];

                if(a > 3)
                {
                    break;
                }
                for(p = 0; p < numPrev; p++)
                {
                    SMEM smem = prev[p];
                    SMEM newSmem = backwardExt(smem, a);
                    newSmem.m = j;

                    if((newSmem.s < min_intv_array[i]) && ((smem.n - smem.m + 1) >= minSeedLen))
                    {
                        cur_j = j;

                        matchArray[numTotalSmem++] = smem;
                        break;
                    }
                    if((newSmem.s >= min_intv_array[i]) && (newSmem.s != curr_s))
                    {
                        curr_s = newSmem.s;
                        prev[numCurr++] = newSmem;
#ifdef ENABLE_PREFETCH
                        _mm_prefetch((const char *)(&cp_occ[(newSmem.k) >> CP_SHIFT]), _MM_HINT_T0);
                        _mm_prefetch((const char *)(&cp_occ[(newSmem.k + newSmem.s) >> CP_SHIFT]), _MM_HINT_T0);
#endif
                        break;
                    }
                }
                p++;
                for(; p < numPrev; p++)
                {
                    SMEM smem = prev[p];

                    SMEM newSmem = backwardExt(smem, a);
                    newSmem.m = j;


                    if((newSmem.s >= min_intv_array[i]) && (newSmem.s != curr_s))
                    {
                        curr_s = newSmem.s;
                        prev[numCurr++] = newSmem;
#ifdef ENABLE_PREFETCH
                        _mm_prefetch((const char *)(&cp_occ[(newSmem.k) >> CP_SHIFT]), _MM_HINT_T0);
                        _mm_prefetch((const char *)(&cp_occ[(newSmem.k + newSmem.s) >> CP_SHIFT]), _MM_HINT_T0);
#endif
                    }
                }
                numPrev = numCurr;
                if(numCurr == 0)
                {
                    break;
                }
            }
            if(numPrev != 0)
            {
                SMEM smem = prev[0];
                if(((smem.n - smem.m + 1) >= minSeedLen))
                {

                    matchArray[numTotalSmem++] = smem;
                }
                numPrev = 0;
            }
        }
        query_pos_array[i] = next_x;
    }
    (*__numTotalSmem) = numTotalSmem;
    _mm_free(prevArray);
}

// ===== Lockstep SMEM batching =====
// Per-slot state for the lockstep SMEM walk. One instance per in-flight read
// in the batch. Every field mirrors a per-read local in the scalar
// getSMEMsOnePosOneThread body, so parity with the scalar path follows from
// composing the existing primitives (backwardExt, count[], cp_occ) in the
// same sequence the scalar would.
//
// The prev[] and match_buf[] buffers are heap-allocated by the driver
// (getSMEMsOnePosOneThread_lockstep), sized from the batch's
// max_readlength. No compile-time cap — long reads (PacBio HiFi, ONT)
// fit cleanly.

enum LockstepPhase : uint8_t {
    PH_FWD      = 0,  // forward extension inner loop active
    PH_BWD_INIT = 1,  // between-phases housekeeping pending
    PH_BWD      = 2,  // backward search outer loop active
    PH_DONE     = 3,  // slot finished (ready: match_buf not yet retired) or idle
    PH_S1_LF    = 4,  // unique forward match: LF walk from k to a sampled SA row
    PH_S1_CMP   = 5,  // unique forward match: compare the rest of the read with the text
    PH_FWD_TABLE = 6  // forward walk's first visit: its first steps from seed_jump
};

/* A phase-1 forward walk that reaches a unique match (s == 1) goes on comparing
 * the read with the reference text only when at least this many read bases are
 * left: resolving the match's text position costs an LF walk to a sampled SA
 * row (the SA is sampled by row, so about 7 steps on average at the default
 * 1 sample in 8 rows; sa_compx is per index) and a pac read. */
static constexpr int32_t LS_S1_MIN_REST = 16;

/* The text position of a unique forward match, kept in SMEM.info (a field
 * nothing else reads) for SA resolution. When ls_s1_lf_step resolves
 * v = SA[k] for a walk that started at start_pos (the match's m), it stores
 * q = v - start_pos with bit 63 set. Every SMEM the walk emits from that
 * entry, the entry itself or a backward extension of it, has s == 1 and
 * SA[k] = q + m: a backward step that keeps s == 1 moves the one suffix's
 * start a base left (SA[LF(k)] = SA[k] - 1) as it moves m a base left.
 * get_sa_entries_prefetch then writes the coordinate q + m instead of walking
 * k to a sampled row a second time, when its caller trusts the tag (see
 * trust_s1_tag there). mem_collect_smem then passes the tag on by containment
 * (s1_tag_contained). Every other SMEM the seeding path makes has info 0. */
static constexpr uint64_t LS_S1_SA_TAG = 1ULL << 63;
static inline uint64_t ls_s1_sa_tag(int64_t q)
{
    return LS_S1_SA_TAG | ((uint64_t)q & ~LS_S1_SA_TAG);
}
static inline int64_t ls_s1_sa_q(uint64_t info)
{
    return (int64_t)(info << 1) >> 1;   // q, sign-extended from bit 62
}
/* A unique match carrying the tag, and its one SA coordinate (q + m). */
static inline bool ls_s1_tagged(const SMEM &smem)
{
    return smem.s == 1 && (smem.info & LS_S1_SA_TAG) != 0;
}
static inline int64_t ls_s1_sa_coord(const SMEM &smem)
{
    return ls_s1_sa_q(smem.info) + (int64_t)smem.m;
}

/* An untagged s == 1 SMEM U of a read that lies inside a tagged SMEM T of the
 * same read (m_T <= m_U and n_U <= n_T) takes T's info unchanged: read[m_U..n_U]
 * occurs once in the text, and T's one occurrence, at q_T + m_T, holds it at
 * q_T + m_U, so SA[k_U] = q_T + m_U, the relation the tag states. Any T that
 * contains U gives the same q. a[0..n) is one read's SMEMs sorted by (m, n)
 * (mem_collect_smem's per-read sort), so one forward pass keeps, of the tagged
 * SMEMs before i, the one with the largest n (all have m_T <= m_i). A U whose
 * only container has the same m sorts before it (n_U <= n_T), so a tagged T
 * also tags the untagged s == 1 SMEMs of its m run that precede it. Containment
 * is transitive, so SMEMs tagged here need not join the running maximum.
 * The common case (an untagged SMEM) takes no data-dependent branch.
 * Every step above relies on the one-read, (m, n)-sorted precondition: an
 * unsorted or mixed-read array would hand a tag to an SMEM that is not inside
 * its source and so write a wrong coordinate without any other symptom. The
 * precondition is therefore checked (one predictable compare per SMEM). */
void FMI_search::s1_tag_contained(SMEM *a, int64_t n)
{
    int64_t best_n = -1;      // largest n_T over the tagged s == 1 SMEMs before i
    uint64_t best_info = 0;
    int64_t run = 0;          // first index of i's m run
    uint32_t m_prev = UINT32_MAX;
    for (int64_t i = 0; i < n; i++) {
        const uint64_t info = a[i].info;
        const uint32_t m = a[i].m;
        const bool s1 = a[i].s == 1;
        const int64_t ni = (int64_t)a[i].n;
        if (__builtin_expect(i != 0 && !(a[i].rid == a[0].rid &&
                                         (m > m_prev || (m == m_prev && a[i].n >= a[i - 1].n))), 0))
            err_fatal(__func__, "SMEMs must be one read's, sorted by (m, n), but SMEM %lld "
                      "(rid %u, m %u, n %u) follows SMEM %lld (rid %u, m %u, n %u) and the "
                      "run started at rid %u", (long long)i, a[i].rid, a[i].m, a[i].n,
                      (long long)(i - 1), a[i - 1].rid, a[i - 1].m, a[i - 1].n, a[0].rid);
        run = (m == m_prev) ? run : i;
        m_prev = m;
        if (s1 & ((info & LS_S1_SA_TAG) != 0)) {
            for (int64_t j = run; j < i; j++)
                if (a[j].s == 1 && (a[j].info & LS_S1_SA_TAG) == 0) a[j].info = info;
            if (ni > best_n) {
                best_n = ni;
                best_info = info;
            }
        } else {
            a[i].info = (s1 & (ni <= best_n)) ? best_info : info;
        }
    }
}

/* The unique-match text compare (ls_s1_*), and with it every SA tag, needs the
 * reference this index was built from: its pac, loaded with the index unless
 * the caller skipped it (--meth loads the seed index without), and of the
 * length the index's text implies (forward + reverse complement + sentinel). */
bool FMI_search::s1_text_available() const
{
    return idx != nullptr && idx->pac != nullptr && idx->bns != nullptr &&
           2 * idx->bns->l_pac + 1 == reference_seq_len;
}


// Per-thread cache for the lockstep SMEM driver's prev[]/match_buf[]
// buffers. RAII wrapper exists so the buffers are released at thread
// exit instead of leaking until process exit (LSan complaint, issue
// #116). Holding the cache as `static thread_local LockstepSmemCache`
// inside getSMEMsOnePosOneThread_lockstep keeps the per-call reuse
// behavior unchanged; only the cleanup point moves.
struct LockstepSmemCache {
    SMEM  *prev     = nullptr;
    SMEM  *match    = nullptr;
    size_t per_slot = 0;
    int32_t n_slots = 0;   // lockstep width the buffers were sized for
    ~LockstepSmemCache() {
        if (prev  != nullptr) _mm_free(prev);
        if (match != nullptr) _mm_free(match);
    }
};

// LISA trick #4: hybrid SoA. The per-slot "hot" state stays compact
// (168 bytes) so the driver's pass over slots[0..N) hits a tight
// L1-resident array instead of jumping 32KB strides across embedded
// buffers. The bulk arrays (prev[], match_buf[]; ~32 KB each) live in
// separate per-slot allocations, referenced by pointer (prev and prev_top
// bound the slot's prev entries; the forward phase writes them downward
// from prev_top, see the members below).
struct FMI_search::BatchSlot {
    // Input identity — copied at init, never mutated thereafter.
    int32_t input_idx;           // index into the caller's input arrays
    int32_t rid;                 // rid_array[input_idx]
    int32_t start_pos;           // query_pos_array[input_idx] (saved for bwd init)
    int32_t min_intv;            // min_intv_array[input_idx]
    int32_t readlength;          // seq_[rid].l_seq
    int32_t offset;              // query_cum_len_ar[rid]

    // Output for query_pos_array[input_idx] write-back when the read is retired.
    int32_t next_x;

    // Walk state (mirrors scalar locals).
    SMEM smem;                   // current SA interval
    int32_t j;                   // current query position in the active phase's inner loop
    int32_t cur_j;               // backward phase bookkeeping (scalar's cur_j)
    LockstepPhase phase;

    // Per-slot bulk-buffer state.
    int32_t numPrev;
    int32_t match_count;
    bool    ready;

    // Per-slot bulk arrays (driver-allocated, >= max_readlength entries). The
    // forward phase writes its prev entries downward from prev_top (one past
    // the slot's prev buffer); a walk from start_pos makes at most readlength -
    // start_pos <= max_readlength of them, so the last, lowest one is at or
    // above the buffer's base. ls_prepare_backward points prev at it, leaving
    // prev[0..numPrev) in the reversed order the scalar path builds by
    // swapping, with no swap; the backward phase reads and compacts via prev.
    SMEM    *prev;
    SMEM    *prev_top;
    SMEM    *match_buf;

    // Unique-match text compare (PH_S1_LF, PH_S1_CMP): the LF walk's current
    // row and step count, then (PH_S1_CMP) the text position of read base j.
    int64_t lf_row;
    int64_t lf_off;
    // The walk's first steps from the jump table (seed_jump): when read[x ..
    // x + jdp) is jdp >= 2 ACGT bases (jdp <= seed_jump_depth), jcode is their
    // 2-bit code, first base most significant, and the slot's first forward
    // visit (PH_FWD_TABLE, ls_advance_forward_table) takes depths 2..jdp from
    // the table. jdp == 0: every step is a backwardExt.
    uint32_t jcode;
    int32_t  jdp;
    // The bwtSeed pass's first walk on this read (BsdFirstWalk), while the
    // forward phase still has to note where it stops: the record, the first
    // step at which it may stop (min_len - 1) and its max_intv. bsd_j is
    // INT32_MAX when nothing is pending, so the forward step's one test
    // (j >= bsd_j) is false.
    BsdFirstWalk *bsd;
    int32_t bsd_j;
    int32_t bsd_intv;
};

void FMI_search::ls_prefetch_cp_occ(const BatchSlot *s)
{
#ifdef ENABLE_PREFETCH
    _mm_prefetch((const char *)(&cp_occ[(s->smem.k) >> CP_SHIFT]), _MM_HINT_T0);
    _mm_prefetch((const char *)(&cp_occ[(s->smem.l) >> CP_SHIFT]), _MM_HINT_T0);
    _mm_prefetch((const char *)(&cp_occ[(s->smem.k + s->smem.s) >> CP_SHIFT]), _MM_HINT_T0);
    _mm_prefetch((const char *)(&cp_occ[(s->smem.l + s->smem.s) >> CP_SHIFT]), _MM_HINT_T0);
#else
    (void)s;
#endif
}

// The 2-bit code of q[0..n), first base most significant, for a jump-table
// lookup; *any is the OR of the n bases, so *any >= 4 when one is an N.
static inline uint32_t seed_jump_pack(const uint8_t *q, int32_t n, uint32_t *any)
{
    uint32_t code = 0, a = 0;
    for (int32_t i = 0; i < n; i++) {
        code = (code << 2) | (q[i] & 3u);
        a |= q[i];
    }
    *any = a;
    return code;
}

// Hint a jump-table entry: a 24-byte entry can straddle two lines.
inline void FMI_search::seed_jump_prefetch(const SeedJumpEntry *e)
{
#ifdef ENABLE_PREFETCH
    _mm_prefetch((const char *)e, _MM_HINT_T0);
    _mm_prefetch((const char *)e + sizeof(SeedJumpEntry) - 1, _MM_HINT_T0);
#else
    (void)e;
#endif
}

// Populate a slot from the caller's input arrays at the given input_idx.
// After this call either:
//   phase == PH_FWD  — slot is ready to step the forward-extension inner loop
//   phase == PH_FWD_TABLE — the same, with its first steps planned from the
//                      jump table (ls_advance_forward_table)
//   phase == PH_DONE — the first base is non-ACGT; scalar skips this read, so
//                      we match that (zero matches emitted, ready to retire).
void FMI_search::ls_init_slot(BatchSlot *s,
                              int32_t input_idx,
                              const int32_t *query_pos_array,
                              const int32_t *min_intv_array,
                              const int32_t *rid_array,
                              const bseq1_t *seq_,
                              const int32_t *query_cum_len_ar,
                              const uint8_t *enc_qdb,
                              BsdFirstWalkHint *bsd_hint)
{
    s->input_idx   = input_idx;
    s->rid         = rid_array[input_idx];
    s->bsd_j       = INT32_MAX;
    s->start_pos   = query_pos_array[input_idx];
    s->min_intv    = min_intv_array[input_idx];
    /* backwardExt_konly leaves the new k stale on its s'->0 exit and relies on
     * the caller discarding that seed via `newSmem.s < min_intv` — which only
     * holds when min_intv >= 1 (an s'==0 result must compare < min_intv). Every
     * current source of min_intv satisfies this (mem_collect_smem seeds 1 or
     * p->s+1 >= 2; max_mem_intv is guarded > 0 and never reaches this path), so
     * this guard documents and enforces the precondition rather than papering
     * over a live bug. Once per slot init — off the per-step hot path. */
    xassert(s->min_intv >= 1, "ls_init_slot: min_intv must be >= 1");
    s->readlength  = seq_[s->rid].l_seq;
    s->offset      = query_cum_len_ar[s->rid];
    s->next_x      = s->start_pos + 1;
    s->numPrev     = 0;
    s->match_count = 0;
    s->ready       = false;

    int32_t x = s->start_pos;
    uint8_t a = enc_qdb[s->offset + x];
    if (a < 4) {
        s->smem.info = 0;   // no text position yet (ls_s1_sa_tag)
        s->smem.rid = s->rid;
        s->smem.m   = x;
        s->smem.n   = x;
        s->smem.k   = count[a];
        s->smem.l   = count[3 - a];
        s->smem.s   = count[a + 1] - count[a];
        s->j        = x + 1;
        s->phase    = PH_FWD;
        /* Plan the walk's first steps from the jump table: the longest run
         * read[x .. x + jdp) of ACGT bases with jdp <= seed_jump_depth (0 when
         * there is no table), and hint the entries of depths 2..jdp a pass
         * ahead of the visit that reads them (between the cp_occ blocks the
         * other slots stream, even levels 5-7, 512 KB, miss L2 unhinted). */
        s->jdp = 0;
        int32_t lim = s->readlength - x;
        if (lim > seed_jump_depth) lim = seed_jump_depth;
        if (lim >= 2) {
            const uint8_t *q = enc_qdb + s->offset + x;
            uint32_t any;
            uint32_t code = seed_jump_pack(q, lim, &any);
            int32_t d = lim;
            if (any >= 4) {   // an N: the run ends at the first one
                d = 1;
                while (q[d] < 4) d++;
                code >>= 2 * (lim - d);
            }
            if (d >= 2) {
                s->jcode = code;
                s->jdp   = d;
                s->phase = PH_FWD_TABLE;
#ifdef ENABLE_PREFETCH
                for (int32_t L = 2; L <= d; L++) {
                    seed_jump_prefetch(seed_jump + seed_jump_level_off(L) + (code >> (2 * (d - L))));
                }
#endif
            }
        }
        /* The bwtSeed pass walks from x = 0 with the same steps as this walk
         * as long as this one goes on, which at min_intv 1 is while s >= 1. */
        if (bsd_hint != nullptr && x == 0 && s->min_intv == 1 && bsd_hint->max_intv >= 1) {
            s->bsd      = &bsd_hint->rec[s->rid];
            s->bsd_j    = bsd_hint->min_len - 1;
            s->bsd_intv = bsd_hint->max_intv;
        }
    } else {
        // Scalar path skips the whole read when the first base is N.
        // Match that by going straight to DONE with zero matches.
        s->phase = PH_DONE;
        s->ready = true;
    }
}
// The bwtSeed pass's first walk (BsdFirstWalk) at forward step j = s->j, which
// extended read[0..j-1] to the interval (k, l, sz): the walk stops here when
// sz < max_intv (the caller tests j >= min_len - 1), emitting the interval as
// its seed when sz > 0, and goes on at x = j + 1.
inline void FMI_search::ls_bsd_note(BatchSlot *s, int64_t k, int64_t l, int64_t sz)
{
    if (sz < s->bsd_intv) {
        BsdFirstWalk *r = s->bsd;
        r->k = k;
        r->l = l;
        r->s = sz;
        r->n = s->j;
        r->x = s->j + 1;
        s->bsd_j = INT32_MAX;
    }
}

// The forward walk stops before the bwtSeed walk has: at read base p, the read's
// end (p >= readlength), a non-ACGT base, or an ACGT base the interval does not
// extend with (s = 0 from step p on). From there the bwtSeed walk steps with
// s = 0, which ends it at the first non-ACGT base p (x = p + 1), at the first
// step p >= min_len - 1 (s = 0 < max_intv, no seed; x = p + 1), or at the
// read's end (x = readlength).
void FMI_search::ls_bsd_resolve(BatchSlot *s, const uint8_t *enc_qdb, int32_t p)
{
    const uint8_t *q = enc_qdb + s->offset;
    const int32_t rl = s->readlength;
    while (p < rl && q[p] < 4 && p < s->bsd_j) p++;
    BsdFirstWalk *r = s->bsd;
    r->s = 0;
    r->x = p < rl ? p + 1 : rl;
    s->bsd_j = INT32_MAX;
}

// Advance one slot through one step of forward extension.
// Mirrors the inner j-loop body of the scalar getSMEMsOnePosOneThread
// (src/FMI_search.cpp — the for(j = x+1; j < readlength; j++) loop).
// On the step that exits forward-ext (end-of-read, non-ACGT at j, or
// s < min_intv), transitions phase to PH_BWD_INIT. When s1_text is set and the
// step leaves a unique match (s == 1) with at least LS_S1_MIN_REST read bases
// to go, transitions to PH_S1_LF instead (ls_s1_lf_step, ls_s1_compare).
void FMI_search::ls_advance_forward_step(BatchSlot *s, const uint8_t *enc_qdb,
                                         bool s1_text)
{
    if (s->j >= s->readlength) {
        // Ran off the end of the read; keep the still-valid smem if any.
        if (s->smem.s >= s->min_intv) {
            s->prev_top[-1 - s->numPrev] = s->smem;
            s->numPrev++;
        }
        if (s->bsd_j != INT32_MAX) ls_bsd_resolve(s, enc_qdb, s->j);
        s->phase = PH_BWD_INIT;
        return;
    }

    uint8_t a = enc_qdb[s->offset + s->j];
    s->next_x = s->j + 1;
    if (a >= 4) {
        // Non-ACGT base terminates forward ext.
        if (s->smem.s >= s->min_intv) {
            s->prev_top[-1 - s->numPrev] = s->smem;
            s->numPrev++;
        }
        if (s->bsd_j != INT32_MAX) ls_bsd_resolve(s, enc_qdb, s->j);
        s->phase = PH_BWD_INIT;
        return;
    }

    SMEM smem_ = s->smem;
    smem_.k = s->smem.l;
    smem_.l = s->smem.k;
    SMEM newSmem_ = backwardExt(smem_, 3 - a);
    SMEM newSmem  = newSmem_;
    newSmem.k = newSmem_.l;
    newSmem.l = newSmem_.k;
    newSmem.n = s->j;

    int32_t s_neq_mask = (newSmem.s != s->smem.s);
    s->prev_top[-1 - s->numPrev] = s->smem;
    s->numPrev += s_neq_mask;

    if (newSmem.s < s->min_intv) {
        /* With a bwtSeed walk pending, min_intv is 1: s = 0 from here on. */
        if (s->bsd_j != INT32_MAX) ls_bsd_resolve(s, enc_qdb, s->j);
        s->next_x = s->j;
        s->phase = PH_BWD_INIT;
        return;
    }

    if (s->j >= s->bsd_j) ls_bsd_note(s, newSmem.k, newSmem.l, newSmem.s);
    s->smem = newSmem;
    s->j++;
    /* A unique match (s == 1; newSmem.s >= min_intv above, so min_intv <= 1)
     * keeps its k on every further forward step: the one suffix it names is
     * either extended or dropped. So instead of stepping the rest of the read
     * one base at a time, resolve that suffix's text position from k and
     * compare the read with the text (ls_s1_lf_step, ls_s1_compare). */
    if (s1_text && newSmem.s == 1 && s->readlength - s->j >= LS_S1_MIN_REST) {
        s->lf_row = newSmem.k;
        s->lf_off = 0;
        s->phase  = PH_S1_LF;
        ls_prefetch_sa_row(newSmem.k);
        return;
    }
#ifdef ENABLE_PREFETCH
    /* The NEXT forward step swaps k<->l (see smem_ above) before calling
     * backwardExt, which reads sp = smem.k and ep = smem.k + smem.s. After the
     * swap that is OUR smem.l and smem.l + smem.s -- so those are the two lines
     * to prefetch. The old code prefetched smem.k (an address the next step
     * never touches) and smem.l, leaving the ep line to always miss cold.
     * Same targets as bsd_prefetch_cp_occ(), which fixed this for round 3 in
     * #242; the fix was simply never propagated here. Prefetch is a pure hint,
     * so this cannot change output. */
    _mm_prefetch((const char *)(&cp_occ[(s->smem.l) >> CP_SHIFT]), _MM_HINT_T0);
    _mm_prefetch((const char *)(&cp_occ[(s->smem.l + s->smem.s) >> CP_SHIFT]), _MM_HINT_T0);
#endif
}

// Prefetch what the next call_one_step on row p reads: the SA sample lines
// for a sampled row, else the row's checkpoint block (the same selection as
// get_sa_entries_prefetch's prefetch_row). A hint only.
inline void FMI_search::ls_prefetch_sa_row(int64_t p) const
{
#ifdef ENABLE_PREFETCH
    uintptr_t m = (uintptr_t)0 - (uintptr_t)((p & sa_compx_mask) == 0);
    BWA_OPAQUE(m);
    const uintptr_t cp = (uintptr_t)&cp_occ[p >> CP_SHIFT];
    const uintptr_t ms = (uintptr_t)&sa_ms_byte[p >> sa_compx];
    const uintptr_t ls = (uintptr_t)&sa_ls_word[p >> sa_compx];
    _mm_prefetch((const char *)((ms & m) | (cp & ~m)), _MM_HINT_T0);
    _mm_prefetch((const char *)((ls & m) | (cp & ~m)), _MM_HINT_T0);
#else
    (void)p;
#endif
}

// The bidirectional index's text is the forward reference (pac, l_pac bases)
// followed by its reverse complement (libsais_build.cpp, write_doubled_pac).
// ls_text_match compares read bases q[j..rl) with text bases t.. and returns
// the first j at which they differ or the text ends (t reaches 2 * l_pac), or
// rl. A non-ACGT read base (> 3) never equals a text base. The aligned bodies
// compare 8 bases at a time: a pac byte expands through a table to its four
// bases, one per byte in text order (forward), or complemented in reverse
// order (the reverse-complement half, read from the end of pac backwards).
struct LsPacLut {
    uint32_t fwd[256];
    uint32_t rev[256];
    constexpr LsPacLut() : fwd(), rev() {
        for (int b = 0; b < 256; b++) {
            uint32_t f = 0, r = 0;
            for (int i = 0; i < 4; i++) {
                const uint32_t base = (uint32_t)(b >> (6 - 2 * i)) & 3;  // _get_pac order
                f |= base << (8 * i);
                r |= (3 - base) << (8 * (3 - i));
            }
            fwd[b] = f;
            rev[b] = r;
        }
    }
};
static constexpr LsPacLut ls_pac_lut{};

static inline int32_t ls_text_match(const uint8_t *q, int32_t j, int32_t rl,
                                    const uint8_t *pac, int64_t l_pac, int64_t t)
{
    const int64_t tend = 2 * l_pac;
    while (j < rl && t < l_pac) {
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
        if ((t & 3) == 0) {
            while (j + 8 <= rl && t + 8 <= l_pac) {
                uint64_t r;
                memcpy(&r, q + j, 8);
                const uint8_t *pb = pac + (t >> 2);
                const uint64_t x = r ^ ((uint64_t)ls_pac_lut.fwd[pb[0]] |
                                        ((uint64_t)ls_pac_lut.fwd[pb[1]] << 32));
                if (x != 0) return j + (int32_t)(__builtin_ctzll(x) >> 3);
                j += 8;
                t += 8;
            }
            if (!(j < rl && t < l_pac)) break;
        }
#endif
        if (q[j] != (uint8_t)_get_pac(pac, t)) return j;
        j++;
        t++;
    }
    while (j < rl && t < tend) {
        const int64_t u = tend - 1 - t;   // pac base whose complement is text base t
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
        if ((u & 3) == 3) {
            int64_t uu = u;
            while (j + 8 <= rl && uu >= 7) {
                uint64_t r;
                memcpy(&r, q + j, 8);
                const uint8_t *pb = pac + (uu >> 2);
                const uint64_t x = r ^ ((uint64_t)ls_pac_lut.rev[pb[0]] |
                                        ((uint64_t)ls_pac_lut.rev[pb[-1]] << 32));
                if (x != 0) return j + (int32_t)(__builtin_ctzll(x) >> 3);
                j += 8;
                uu -= 8;
            }
            t = tend - 1 - uu;
            if (!(j < rl && t < tend)) break;
            if (q[j] != (uint8_t)(3 - _get_pac(pac, uu))) return j;
            j++;
            t++;
            continue;
        }
#endif
        if (q[j] != (uint8_t)(3 - _get_pac(pac, u))) return j;
        j++;
        t++;
    }
    return j;
}

// One pass of a unique match's SA walk: one LF step from lf_row toward a
// sampled row, or, on a sampled (or the sentinel) row, the match's text
// position. call_one_step is the SA resolver's own step, with the exact
// sentinel offset (drop_sentinel_offset 0) under every --compat, since the
// position here locates text, not a reported coordinate.
void FMI_search::ls_s1_lf_step(BatchSlot *s, const uint8_t *pac, int64_t l_pac)
{
    int64_t v;
    if (call_one_step(s->lf_row, v, s->lf_off, 0) == 0) {
        s->lf_row = v;
        ls_prefetch_sa_row(v);
        return;
    }
    // v = SA[k], where read[m..n] starts in the text; read base j = n + 1
    // faces text base v + (j - m). Keep v for SA resolution (ls_s1_sa_tag);
    // m is still the walk's start_pos here.
    s->smem.info = ls_s1_sa_tag(v - (int64_t)s->smem.m);
    const int64_t t = v + (int64_t)(s->j - (int32_t)s->smem.m);
    s->lf_row = t;
    s->phase  = PH_S1_CMP;
#ifdef ENABLE_PREFETCH
    // The bases still to compare span at most a couple of pac lines: hint the
    // first and the last (clamped to the text; the walk stops at its end).
    const int64_t tend = 2 * l_pac;
    int64_t t_last = t + (int64_t)(s->readlength - s->j) - 1;
    if (t_last > tend - 1) t_last = tend - 1;
    const int64_t t0 = t < tend ? t : tend - 1;
    const int64_t b0 = (t0 < l_pac ? t0 : tend - 1 - t0) >> 2;
    const int64_t b1 = (t_last < l_pac ? t_last : tend - 1 - t_last) >> 2;
    _mm_prefetch((const char *)&pac[b0], _MM_HINT_T0);
    _mm_prefetch((const char *)&pac[b1], _MM_HINT_T0);
#else
    (void)l_pac;
#endif
}

// Finish a unique match's forward walk by comparing read bases j.. with text
// bases lf_row.., and take the exit the base-at-a-time walk would take at the
// first base that ends it: the read's end (next_x stays readlength), a
// non-ACGT read base (next_x = j + 1), or a base the text does not continue
// with, the text's end (the sentinel) included (next_x = j). Every exit keeps
// the match (s = 1 >= min_intv) as the walk's last prev entry, with n = j - 1;
// k and s are those the walk would end with. l is left stale, as the backward
// phase's backwardExt_konly leaves it: nothing reads l after seeding.
void FMI_search::ls_s1_compare(BatchSlot *s, const uint8_t *enc_qdb,
                               const uint8_t *pac, int64_t l_pac)
{
    const uint8_t *q = enc_qdb + s->offset;
    const int32_t rl = s->readlength;
    const int32_t j = ls_text_match(q, s->j, rl, pac, l_pac, s->lf_row);
    const int32_t next_x = j >= rl ? rl : q[j] > 3 ? j + 1 : j;
    /* A pending bwtSeed walk had s = 1 on read bases s->j .. j - 1 and stops
     * on the first of them at or past min_len - 1 when 1 < max_intv (step
     * s->j - 1 came before min_len - 1, or ls_bsd_note would have recorded it),
     * with this match's k as its seed; otherwise it goes on from base j. */
    if (s->bsd_j != INT32_MAX) {
        if (s->bsd_intv > 1 && s->bsd_j < j) {
            BsdFirstWalk *r = s->bsd;
            r->k = s->smem.k;
            r->l = s->smem.l;
            r->s = 1;
            r->n = s->bsd_j;
            r->x = s->bsd_j + 1;
            s->bsd_j = INT32_MAX;
        } else {
            ls_bsd_resolve(s, enc_qdb, j);
        }
    }
    s->smem.n = j - 1;
    s->next_x = next_x;
    s->prev_top[-1 - s->numPrev] = s->smem;
    s->numPrev++;
    s->phase = PH_BWD_INIT;
}

// A walk's first forward visit when ls_init_slot planned table steps (jdp >= 2):
// the steps to depths 2..jdp, read[x + 1 .. x + jdp), each one as
// ls_advance_forward_step takes it, with the new interval read from the jump
// table instead of a backwardExt. A table entry is the interval the walk's own
// forward step reaches from the 1-base interval (seed_jump_children), and
// read[x .. x + jdp) is ACGT inside the read, so only the step's s < min_intv
// exit (phase 2's min_intv = s + 1 included) and the unique-match hand-off to
// PH_S1_LF can end the run early. Otherwise the slot steps on from depth
// jdp + 1 with the same prefetch as a forward step.
void FMI_search::ls_advance_forward_table(BatchSlot *s, const uint8_t *enc_qdb,
                                          bool s1_text)
{
    const int32_t dp = s->jdp;
    const uint32_t code = s->jcode;
    const int32_t x = s->start_pos;
    const int32_t min_intv = s->min_intv;
    SMEM *const top = s->prev_top;
    int32_t numPrev = s->numPrev;
    SMEM cur = s->smem;
    s->jdp = 0;

    int64_t off = seed_jump_level_off(1);
    for (int32_t d = 2; d <= dp; d++) {
        off = off * 4 + 4;   // seed_jump_level_off(d)
        const SeedJumpEntry *e = seed_jump + off + (code >> (2 * (dp - d)));
        const int32_t j = x + d - 1;
        SMEM nw = cur;
        nw.k = e->k;
        nw.l = e->l;
        nw.s = e->s;
        nw.n = j;

        top[-1 - numPrev] = cur;
        numPrev += (nw.s != cur.s);

        if (nw.s < min_intv) {
            s->numPrev = numPrev;
            s->smem = cur;
            s->j = j;
            /* With a bwtSeed walk pending, min_intv is 1: s = 0 from here on. */
            if (s->bsd_j != INT32_MAX) ls_bsd_resolve(s, enc_qdb, j);
            s->next_x = j;
            s->phase = PH_BWD_INIT;
            return;
        }
        if (j >= s->bsd_j) {
            s->j = j;
            ls_bsd_note(s, nw.k, nw.l, nw.s);
        }
        cur = nw;
        if (s1_text && nw.s == 1 && s->readlength - (j + 1) >= LS_S1_MIN_REST) {
            s->numPrev = numPrev;
            s->smem = cur;
            s->j = j + 1;
            s->next_x = j + 1;
            s->lf_row = nw.k;
            s->lf_off = 0;
            s->phase  = PH_S1_LF;
            ls_prefetch_sa_row(nw.k);
            return;
        }
    }
    s->numPrev = numPrev;
    s->smem = cur;
    s->j = x + dp;
    s->next_x = x + dp;
    s->phase = PH_FWD;
#ifdef ENABLE_PREFETCH
    _mm_prefetch((const char *)(&cp_occ[(cur.l) >> CP_SHIFT]), _MM_HINT_T0);
    _mm_prefetch((const char *)(&cp_occ[(cur.l + cur.s) >> CP_SHIFT]), _MM_HINT_T0);
#endif
}

// Between-phases housekeeping: point prev at the forward phase's entries and
// set up backward-phase cursors. The slot leaves as PH_BWD (ready for
// ls_advance_backward_step) or PH_DONE (the backward outer loop would end on
// its first step, so go straight to the final prev[0] emit).
void FMI_search::ls_prepare_backward(BatchSlot *s)
{
    static_assert(sizeof(void *) != 8 || sizeof(BatchSlot) == 168,
                  "BatchSlot size changed: update the slot-size comments");
    // The forward phase wrote its numPrev entries downward from prev_top, so
    // from the last one written they run in the reversed order the scalar path
    // makes by swapping. The backward phase's compaction writes prev[numCurr]
    // with numCurr <= p, at or below the entry it reads: inside the entries.
    s->prev = s->prev_top - s->numPrev;
    s->cur_j = s->readlength;
    s->j = s->start_pos - 1;  // first position for the backward outer loop

    if (s->numPrev == 0) {
        // Nothing to emit; scalar's final `if (numPrev != 0)` block is a no-op.
        s->phase = PH_DONE;
        s->ready = true;
        return;
    }
    if (s->j < 0) {
        // Backward outer loop cannot execute (start_pos == 0). Transition to
        // PH_BWD so the first ls_advance_backward_step call sees j < 0, jumps
        // to DONE, and runs the final prev[0] emit with the correct minSeedLen.
        s->phase = PH_BWD;
        return;
    }
#ifdef ENABLE_PREFETCH
    // The first backward step extends every prev entry from its k and k + s
    // checkpoint blocks. The forward step prefetched only the l side, so
    // without this no prev entry's k-side blocks get a hint. Hint every entry
    // into L1 now, a full pass of the other slots ahead of that step. A hint
    // only; nothing reads its result.
    for (int p = 0; p < s->numPrev; p++) {
        _mm_prefetch((const char *)(&cp_occ[(s->prev[p].k) >> CP_SHIFT]), _MM_HINT_T0);
        _mm_prefetch((const char *)(&cp_occ[(s->prev[p].k + s->prev[p].s) >> CP_SHIFT]), _MM_HINT_T0);
    }
#endif
    s->phase = PH_BWD;
}

// Advance one slot through ONE outer-j iteration of backward search.
// Mirrors the body of the scalar `for (j = x-1; j >= 0; j--)` outer loop
// in getSMEMsOnePosOneThread. On terminating conditions (numCurr == 0,
// j < 0 after decrement, or non-ACGT at j), runs the scalar's final
// prev[0] emit and transitions to PH_DONE.
void FMI_search::ls_advance_backward_step(BatchSlot *s,
                                          const uint8_t *enc_qdb,
                                          int32_t minSeedLen)
{
    if (s->j < 0) goto DONE;

    {
        int32_t numCurr = 0;
        int32_t curr_s = -1;
        uint8_t a = enc_qdb[s->offset + s->j];
        if (a > 3) goto DONE;

        /* numPrev==1 straight-line path: the two-loop machinery and the curr_s
         * dedup collapse to a single element. Byte-identical. A hand-kept copy
         * of the general loops below (emit predicate, prefetch pair, numCurr/j--
         * bookkeeping): mirror any edit to that logic here, or numPrev==1 reads
         * (unique matches, the common case) silently diverge; the lockstep
         * parity test's Case 13 guards it. It keeps the member call and s->
         * reads on purpose: one step has no per-iteration reloads for the
         * general loops' hoisted locals to save. */
        if (s->numPrev == 1) {
            SMEM smem = s->prev[0];
            SMEM newSmem = backwardExt_konly(smem, a);
            newSmem.m = s->j;
            if (newSmem.s < s->min_intv) {
                if ((smem.n - smem.m + 1) >= minSeedLen) {
                    s->cur_j = s->j;
                    s->match_buf[s->match_count++] = smem;
                }
                /* numCurr stays 0 */
            } else {
                s->prev[0] = newSmem;
#ifdef ENABLE_PREFETCH
                _mm_prefetch((const char *)(&cp_occ[(newSmem.k) >> CP_SHIFT]), _MM_HINT_T0);
                _mm_prefetch((const char *)(&cp_occ[(newSmem.k + newSmem.s) >> CP_SHIFT]), _MM_HINT_T0);
#endif
                numCurr = 1;
            }
            s->numPrev = numCurr;
            s->j--;
            if (numCurr == 0) goto DONE;
            return;
        }

        /* The slot's fields, cp_occ and count[a] are read once into locals for
         * the two general loops, and numPrev / j written back once after them.
         * Read through s and this, clang reloads them after every prev[] store
         * (the SMEM stores may alias them as far as it can tell), an extra L1
         * load ahead of each prev[p] load, compare and prefetch. prev[] and
         * match_buf[] are driver-owned buffers that overlap neither the slot nor
         * the index, so this is byte-identical. */
        SMEM *const prev = s->prev;
        const int32_t numPrev = s->numPrev;
        const int32_t min_intv = s->min_intv;
        const int32_t j = s->j;
        const CP_OCC *const cpo = cp_occ;
        const int64_t count_a = count[a];
        int p;
        for (p = 0; p < numPrev; p++) {
            SMEM smem = prev[p];
            SMEM newSmem = backwardExt_konly_cp(cpo, count_a, smem, a);
            newSmem.m = j;
            if ((newSmem.s < min_intv) && ((smem.n - smem.m + 1) >= minSeedLen)) {
                s->cur_j = j;
                s->match_buf[s->match_count++] = smem;
                break;
            }
            if ((newSmem.s >= min_intv) && (newSmem.s != curr_s)) {
                curr_s = newSmem.s;
                prev[numCurr++] = newSmem;
#ifdef ENABLE_PREFETCH
                _mm_prefetch((const char *)(&cpo[(newSmem.k) >> CP_SHIFT]), _MM_HINT_T0);
                _mm_prefetch((const char *)(&cpo[(newSmem.k + newSmem.s) >> CP_SHIFT]), _MM_HINT_T0);
#endif
                break;
            }
        }
        p++;
        for (; p < numPrev; p++) {
            SMEM smem = prev[p];
            SMEM newSmem = backwardExt_konly_cp(cpo, count_a, smem, a);
            newSmem.m = j;
            if ((newSmem.s >= min_intv) && (newSmem.s != curr_s)) {
                curr_s = newSmem.s;
                prev[numCurr++] = newSmem;
#ifdef ENABLE_PREFETCH
                _mm_prefetch((const char *)(&cpo[(newSmem.k) >> CP_SHIFT]), _MM_HINT_T0);
                _mm_prefetch((const char *)(&cpo[(newSmem.k + newSmem.s) >> CP_SHIFT]), _MM_HINT_T0);
#endif
            }
        }
        s->numPrev = numCurr;
        s->j = j - 1;          // advance for the next outer step
        if (numCurr == 0) goto DONE;
        return;                // stay in PH_BWD; next call continues
    }

DONE:
    // Scalar end-of-function final prev[0] emit (lines ~650-659 of scalar).
    if (s->numPrev != 0) {
        SMEM smem = s->prev[0];
        if ((smem.n - smem.m + 1) >= minSeedLen) {
            s->match_buf[s->match_count++] = smem;
        }
        s->numPrev = 0;
    }
    s->phase = PH_DONE;
    s->ready = true;
}
// ===== End lockstep SMEM batching =====

// ===== Lockstep bwtSeed batching =====
// Per-slot state for the forward-only bwtSeed walk. Mirrors the locals of
// the scalar bwtSeedStrategyAllPosOneThread inner-j loop, plus an outer-x
// cursor. Each slot owns the entire scalar `while (x < x_end)` outer
// loop — outer transitions are handled inline in bsd_advance_step (Heng
// review action: don't burn a stepping-pass on outer-only housekeeping
// when the inner-j loop is shallow, which is exactly bwtSeed's regime).
//
// match_buf[] is a per-slot accumulator flushed in input order; sized
// by the driver from the batch's max_readlength.
enum BwtSeedPhase : uint8_t {
    BSD_FWD  = 0,  // in active inner-j extension (s->j is next j to step;
                   // when s->j < 0 the next step does outer seek first)
    BSD_DONE = 1,  // slot retired; match_buf ready for in-order flush
};

// Per-thread cache for the lockstep bwtSeed driver's match_buf[]. RAII
// wrapper parallels LockstepSmemCache above; released at thread exit so
// LSan stays quiet (issue #116 family). Tracked separately from the SMEM
// cache because the access pattern is different (single buffer, not
// prev/match pair) and the watermark can move independently.
struct LockstepBwtSeedCache {
    SMEM  *match    = nullptr;
    size_t per_slot = 0;
    int32_t n_slots = 0;   // lockstep width the buffer was sized for
    ~LockstepBwtSeedCache() {
        if (match != nullptr) _mm_free(match);
    }
};

struct FMI_search::BwtSeedSlot {
    int32_t input_idx;     // index into the caller's input arrays; doubles as smem.rid
    int32_t max_intv;      // max_intv_array[input_idx] cached at init
    int32_t readlength;    // seq_[input_idx].l_seq
    int32_t offset;        // query_cum_len_ar[input_idx]

    SMEM    smem;          // current SA interval; rid/m fixed at outer-x init
    int32_t x;             // outer cursor: position of current smem's seed base
    int32_t j;             // inner cursor: next j to step; j < 0 = need outer seek
    int32_t x_end;         // first x from which no walk can emit (bsd_init_slot)

    BwtSeedPhase phase;
    int32_t match_count;
    bool    ready;

    SMEM   *match_buf;     // -> thread-local cache slice, sized = max_readlength

    // Jump-table entry for read[x .. x + jd), planned (and prefetched) when the
    // walk started at x and taken on the slot's next visit; NULL = step.
    const SeedJumpEntry *jump;
};

// First x from which no bwtSeed walk can emit. A walk from x emits only at a
// step whose length n - m + 1 is at least minSeedLen, and its last possible
// step ends at the read's last base, so a walk starting at
// x >= readlength - minSeedLen + 1 never emits, nor does any walk after it.
// minSeedLen <= 1 (-k 0 or below) keeps the read end as the bound, so no walk
// starts past the read. Shared by the scalar and lockstep drivers.
static inline int32_t bwtseed_x_end(int32_t readlength, int32_t minSeedLen)
{
    return minSeedLen > 1 ? readlength - minSeedLen + 1 : readlength;
}

// Same-slot T0 prefetch for the next inner-j step. The bwtSeed walk is
// forward-only: bsd_advance_step swaps k<->l and calls backwardExt on the
// swapped interval, so the next occ read lands in cp_occ at (smem.l >> CP_SHIFT)
// and ((smem.l + smem.s) >> CP_SHIFT) — i.e. sp and ep of the swapped interval,
// NOT smem.k. Prefetching smem.k (the un-swapped field) targets the wrong lines,
// so every continuation-step prefetch was a no-op and each step paid full cp_occ
// latency. smem.s is unchanged by the k<->l swap.
void FMI_search::bsd_prefetch_cp_occ(const BwtSeedSlot *s)
{
#ifdef ENABLE_PREFETCH
    _mm_prefetch((const char *)(&cp_occ[(s->smem.l) >> CP_SHIFT]), _MM_HINT_T0);
    _mm_prefetch((const char *)(&cp_occ[(s->smem.l + s->smem.s) >> CP_SHIFT]), _MM_HINT_T0);
#else
    (void)s;
#endif
}

// Plan the jump for a walk that starts at s->x (smem already the 1-base
// interval, j = x + 1): when read[x .. x + jd) is jd ACGT bases, point s->jump
// at that string's entry in the level-jd table jlev and prefetch it; the next
// visit (bsd_advance_step) takes it. Otherwise s->jump stays NULL and the walk
// steps from x as before, reaching the N or the read end exactly as it did.
// The read-end check cannot fire for a walk that starts before x_end
// (bsd_init_slot), since jd <= minSeedLen - 1 (clamped in
// bwtSeedStrategyAllPosOneThread_lockstep); it keeps the read in bounds anyway.
void FMI_search::bsd_plan_jump(BwtSeedSlot *s, const uint8_t *enc_qdb,
                               const SeedJumpEntry *jlev, int32_t jd)
{
    if (s->x + jd > s->readlength) return;
    const uint8_t *q = enc_qdb + s->offset + s->x;
    uint32_t any;
    const uint32_t code = seed_jump_pack(q, jd, &any);
    if (any >= 4) return;   // some base >= 4 (N): the walk stops there
    s->jump = jlev + code;
    seed_jump_prefetch(s->jump);
}

// Initialize a fresh slot for the given input read. Scans forward from
// x=0 for the first ACGT base (matching the scalar's outer while-loop
// skipping non-ACGT) and either enters BSD_FWD with smem seeded (and its
// jump planned, bsd_plan_jump), or goes straight to BSD_DONE when no walk left
// in the read can emit. With a hint, a read whose first walk the SMEM pass
// recorded starts after that walk instead.
//
// No walk from x_end on can emit (bwtseed_x_end): it would only step to an N
// or the read end. The slot retires as soon as its cursor reaches x_end (here,
// at the outer seek, or where a walk ends in bsd_advance_step) instead of
// walking there.
void FMI_search::bsd_init_slot(BwtSeedSlot *s,
                               int32_t input_idx,
                               const int32_t *max_intv_array,
                               const bseq1_t *seq_,
                               const int32_t *query_cum_len_ar,
                               const uint8_t *enc_qdb,
                               const SeedJumpEntry *jlev, int32_t jd,
                               int32_t minSeedLen,
                               const BsdFirstWalkHint *bsd_hint)
{
    s->input_idx   = input_idx;
    s->max_intv    = max_intv_array[input_idx];
    s->readlength  = seq_[input_idx].l_seq;
    s->offset      = query_cum_len_ar[input_idx];
    s->match_count = 0;
    s->ready       = false;
    s->jump        = nullptr;
    s->x_end       = bwtseed_x_end(s->readlength, minSeedLen);

    int32_t x = 0;
    /* The SMEM pass recorded this read's first walk (BsdFirstWalk): emit its
     * seed and go on where it ended, as the walk itself would. */
    if (bsd_hint != nullptr) {
        const BsdFirstWalk *r = &bsd_hint->rec[input_idx];
        if (r->x > 0 && s->max_intv == bsd_hint->max_intv) {
            if (r->s > 0) {
                s->smem.info = 0;   // untagged: get_sa_entries_prefetch reads info
                s->smem.rid = (uint32_t)input_idx;
                s->smem.m   = 0;
                s->smem.n   = r->n;
                s->smem.k   = r->k;
                s->smem.l   = r->l;
                s->smem.s   = r->s;
                s->match_buf[s->match_count++] = s->smem;
            }
            x = r->x;
        }
    }
    while (x < s->x_end && enc_qdb[s->offset + x] >= 4) x++;
    s->x = x;
    if (x >= s->x_end) {
        s->phase = BSD_DONE;
        s->ready = true;
        return;
    }
    uint8_t a = enc_qdb[s->offset + x];
    s->smem.info = 0;   // untagged: get_sa_entries_prefetch reads info
    s->smem.rid = (uint32_t)input_idx;
    s->smem.m   = x;
    s->smem.n   = x;
    s->smem.k   = count[a];
    s->smem.l   = count[3 - a];
    s->smem.s   = count[a + 1] - count[a];
    s->j        = x + 1;
    s->phase    = BSD_FWD;
    if (jd > 0) bsd_plan_jump(s, enc_qdb, jlev, jd);
}

// Advance one slot by one unit of work. Does one cp_occ access per call when
// the slot is active, so the lockstep pipeline keeps uniform per-pass memory
// parallelism (Heng action 2a — collapse outer transition into the step
// function rather than burning a stepping-pass on housekeeping), except at a
// walk's start when a jump is planned: that visit only prefetches the jump
// table entry, and the next one loads it in place of the walk's first
// jd - 1 steps.
//
// The semantics are identical to the scalar bwtSeedStrategyAllPosOneThread
// inner-j body + outer-x advancement:
//   for(j = x+1; j < readlength; j++) {
//       next_x = j+1; a = enc_qdb[offset + j];
//       if (a < 4) {
//           SMEM newSmem = forward-ext(smem, 3-a); newSmem.n = j;
//           smem = newSmem;
//           if (s < max_intv && (n - m + 1) >= minSeedLen) {
//               if (s > 0) emit smem;
//               break;
//           }
//       } else { break; }
//   }
//   x = next_x;
void FMI_search::bsd_advance_step(BwtSeedSlot *s,
                                  const uint8_t *enc_qdb,
                                  int32_t minSeedLen,
                                  const SeedJumpEntry *jlev, int32_t jd)
{
    // Jump planned when this walk started at x (prefetched then): the interval
    // of read[x .. x + jd), the state the walk reaches after jd - 1 steps, none
    // of which can exit (length < minSeedLen, no N, inside the read).
    if (s->jump != nullptr) {
        const SeedJumpEntry e = *s->jump;
        s->jump   = nullptr;
        s->smem.k = e.k;
        s->smem.l = e.l;
        s->smem.s = e.s;
        s->smem.n = s->x + jd - 1;
        s->j      = s->x + jd;
        bsd_prefetch_cp_occ(s);
        return;
    }

    // Outer seek: if j < 0 the slot just broke out of an inner-j extension
    // and needs to (re-)init smem at the next ACGT x. Spin through non-ACGT
    // bases inline (scalar does the same — at most O(N_run) here, and N runs
    // are rare in Illumina data). Falls through to one inner-j step. From
    // x_end on no walk can emit (bsd_init_slot), so the slot retires there.
    if (s->j < 0) {
        while (s->x < s->x_end && enc_qdb[s->offset + s->x] >= 4) s->x++;
        if (s->x >= s->x_end) {
            s->phase = BSD_DONE;
            s->ready = true;
            return;
        }
        uint8_t a = enc_qdb[s->offset + s->x];
        s->smem.rid = (uint32_t)s->input_idx;
        s->smem.m   = s->x;
        s->smem.n   = s->x;
        s->smem.k   = count[a];
        s->smem.l   = count[3 - a];
        s->smem.s   = count[a + 1] - count[a];
        s->j        = s->x + 1;
        if (jd > 0) {
            bsd_plan_jump(s, enc_qdb, jlev, jd);
            if (s->jump != nullptr) return;   // jump on the next visit
        }
        // Fall through to inner-j step.
    }

    // Inner-j step.
    if (s->j >= s->readlength) {
        // Scalar: inner for-loop exhausted; next_x stays at last iteration's
        // j+1 = readlength, at or past x_end, so the slot retires now.
        s->x = s->readlength;
        s->j = -1;
        s->phase = BSD_DONE;
        s->ready = true;
        return;
    }

    uint8_t a = enc_qdb[s->offset + s->j];
    if (a >= 4) {
        // Scalar: non-ACGT in inner-j sets next_x = j+1 and breaks. A walk
        // ending at x_end or past it is the read's last (bsd_init_slot).
        s->x = s->j + 1;
        s->j = -1;
        if (s->x >= s->x_end) {
            s->phase = BSD_DONE;
            s->ready = true;
        }
        return;
    }

    // Forward extension = backward extension on the BWT of the reverse
    // complement: swap k <-> l around backwardExt(., 3-a) and swap back.
    SMEM smem_ = s->smem;
    smem_.k = s->smem.l;
    smem_.l = s->smem.k;
    SMEM newSmem_ = backwardExt(smem_, 3 - a);
    SMEM newSmem  = newSmem_;
    newSmem.k = newSmem_.l;
    newSmem.l = newSmem_.k;
    newSmem.n = s->j;

    if ((newSmem.s < s->max_intv) && ((newSmem.n - newSmem.m + 1) >= minSeedLen)) {
        if (newSmem.s > 0) {
            s->match_buf[s->match_count++] = newSmem;
        }
        s->x = s->j + 1;
        s->j = -1;
        if (s->x >= s->x_end) {
            s->phase = BSD_DONE;
            s->ready = true;
        }
        return;
    }

    s->smem = newSmem;
    s->j++;
    bsd_prefetch_cp_occ(s);
}
// ===== End lockstep bwtSeed batching =====

void FMI_search::getSMEMsAllPosOneThread(uint8_t *enc_qdb,
                                         int32_t *min_intv_array,
                                         int32_t *rid_array,
                                         int32_t numReads,
                                         int32_t batch_size,
                                         const bseq1_t *seq_,
                                         int32_t *query_cum_len_ar,
                                         int32_t max_readlength,
                                         int32_t minSeedLen,
                                         SMEM *matchArray,
                                         int64_t *__numTotalSmem,
                                         BsdFirstWalkHint *bsd_hint)
{
    size_t query_pos_bytes = (size_t)numReads * sizeof(int32_t);
    int32_t *query_pos_array = (int32_t *)_mm_malloc(query_pos_bytes, 64);
    assert_not_null(query_pos_array, query_pos_bytes, query_pos_bytes);

    int32_t i;
    for(i = 0; i < numReads; i++)
        query_pos_array[i] = 0;
    /* Mark every read's BsdFirstWalk unrecorded before round 1, so a record
     * the bwtSeed pass trusts (x > 0) is always this call's. Only the lockstep
     * forward walk records one; other builds leave them all unrecorded. */
    if (bsd_hint != nullptr)
        for(i = 0; i < numReads; i++)
            bsd_hint->rec[rid_array[i]].x = -1;

    int32_t numActive = numReads;
    (*__numTotalSmem) = 0;

    do
    {
        int32_t head = 0;
        int32_t tail = 0;
        for(head = 0; head < numActive; head++)
        {
            int readlength = seq_[rid_array[head]].l_seq;
            if(query_pos_array[head] < readlength)
            {
                rid_array[tail] = rid_array[head];
                query_pos_array[tail] = query_pos_array[head];
                min_intv_array[tail] = min_intv_array[head];
                tail++;
            }
        }
#if SMEM_LOCKSTEP_N > 1
        getSMEMsOnePosOneThread_lockstep(enc_qdb,
                                         query_pos_array,
                                         min_intv_array,
                                         rid_array,
                                         tail,
                                         batch_size,
                                         seq_,
                                         query_cum_len_ar,
                                         max_readlength,
                                         minSeedLen,
                                         matchArray,
                                         __numTotalSmem,
                                         bsd_hint);
        // Only the first round's walks start at query position 0.
        bsd_hint = nullptr;
#else
        getSMEMsOnePosOneThread(enc_qdb,
                                query_pos_array,
                                min_intv_array,
                                rid_array,
                                tail,
                                batch_size,
                                seq_,
                                query_cum_len_ar,
                                max_readlength,
                                minSeedLen,
                                matchArray,
                                __numTotalSmem);
#endif
        numActive = tail;
    } while(numActive > 0);

    _mm_free(query_pos_array);
}

void FMI_search::getSMEMsOnePosOneThread_lockstep(uint8_t *enc_qdb,
                                                   int32_t *query_pos_array,
                                                   int32_t *min_intv_array,
                                                   int32_t *rid_array,
                                                   int32_t numReads,
                                                   int32_t batch_size,
                                                   const bseq1_t *seq_,
                                                   int32_t *query_cum_len_ar,
                                                   int32_t max_readlength,
                                                   int32_t minSeedLen,
                                                   SMEM *matchArray,
                                                   int64_t *__numTotalSmem,
                                                   BsdFirstWalkHint *bsd_hint)
{
    (void)batch_size;

    if (numReads == 0) return;

    const int32_t N = g_smem_lockstep_n;
    // slots[] is a fixed SMEM_LOCKSTEP_N_MAX stack array indexed in [0, N);
    // the env parser and the host default clamp to that range, but a direct
    // g_smem_lockstep_n write (e.g. a test) could exceed it -- guard the bound.
    xassert(N >= 1 && N <= SMEM_LOCKSTEP_N_MAX,
            "g_smem_lockstep_n out of range [1, SMEM_LOCKSTEP_N_MAX]");
    // LISA trick #4: hybrid SoA layout. `slots[]` holds only the small hot
    // state (168 B per slot: ~2.6 KB at N=16, ~5.3 KB at N=32; see BatchSlot).
    // Bulk per-slot buffers (prev/match_buf) live separately and are reused
    // across calls via thread_local caches sized from the batch's
    // max_readlength so reads of any length fit (see issue #44 / PR #55).
    //
    // The outer driver getSMEMsAllPosOneThread runs this in a do/while
    // loop, so allocating per-call (~2*N*max_readlength*sizeof(SMEM) ≈
    // 384 KB at N=8, max_readlength=1500) imposed measurable allocator
    // pressure. Cache per-thread (so OMP workers don't share) and grow
    // monotonically — max_readlength is bounded by the driver batch and
    // increases rarely in practice.
    //
    // TODO(memory): the cache only grows. For mixed-length workloads (e.g.
    // a single ONT-class read with max_readlength≈1e6 mid-stream — sized
    // to ~640 MB per thread, ~10 GB across 16 OMP workers — followed by
    // short reads), the high-water mark is held until thread exit.
    // Acceptable for the smoke1M Illumina PE150 benchmark (max_readlength≈
    // 150 → ~38 KB/thread). Revisit if a streaming aligner or long-running
    // service pipeline appears: gate the realloc on a configured upper
    // bound (e.g. MAX_SMEM_PER_SLOT) or shrink when cache.per_slot greatly
    // exceeds the current batch.
    /* Fixed-size on the stack (compile-time MAX); only the first N are used,
     * where N = g_smem_lockstep_n is the runtime width resolved at startup
     * (env pin, opt-in probe, host-tuned default, or compile-time default). */
    BatchSlot slots[SMEM_LOCKSTEP_N_MAX] = {};
    static thread_local LockstepSmemCache cache;
    /* The unique-match text compare (ls_s1_*) runs only when the index's
     * reference is at hand (s1_text_available); otherwise the forward walk
     * steps every base. */
    const bool s1_text = s1_text_available();
    const uint8_t *s1_pac = s1_text ? idx->pac : nullptr;
    const int64_t s1_l_pac = s1_text ? idx->bns->l_pac : 0;
    const size_t per_slot_smems = (size_t)max_readlength;
    if (per_slot_smems > cache.per_slot || N > cache.n_slots) {
        if (cache.prev  != nullptr) _mm_free(cache.prev);
        if (cache.match != nullptr) _mm_free(cache.match);
        // Grow on EITHER dimension: g_smem_lockstep_n (N) can rise after first
        // use (BWA3_SMEM_LOCKSTEP_N / the MLP probe / the host default, or a
        // direct write in a test), and sizing only on per_slot
        // would leave slots[N-1] pointing past the allocation.
        const size_t alloc_slots = (size_t)(N > cache.n_slots ? N : cache.n_slots);
        const size_t alloc_per   = per_slot_smems > cache.per_slot ? per_slot_smems : cache.per_slot;
        const size_t total_slot_bytes = alloc_slots * alloc_per * sizeof(SMEM);
        cache.prev  = (SMEM *)_mm_malloc(total_slot_bytes, 64);
        assert_not_null(cache.prev, total_slot_bytes, total_slot_bytes);
        cache.match = (SMEM *)_mm_malloc(total_slot_bytes, 64);
        assert_not_null(cache.match, total_slot_bytes, (size_t)2 * total_slot_bytes);
        cache.per_slot = alloc_per;
        cache.n_slots  = (int32_t)alloc_slots;
    }
    for (int32_t s = 0; s < N; s++) {
        // The forward phase writes below prev_top; ls_prepare_backward moves
        // prev to its last entry, so the base set here is only a defined value.
        slots[s].prev      = cache.prev  + (size_t)s * cache.per_slot;
        slots[s].prev_top  = cache.prev  + (size_t)(s + 1) * cache.per_slot;
        slots[s].match_buf = cache.match + (size_t)s * cache.per_slot;
    }

    // Seed the first min(N, numReads) slots.
    int32_t initial = (numReads < N) ? numReads : N;
    for (int32_t s = 0; s < initial; s++) {
        ls_init_slot(&slots[s], s, query_pos_array, min_intv_array,
                     rid_array, seq_, query_cum_len_ar, enc_qdb, bsd_hint);
        if (slots[s].phase == PH_FWD) ls_prefetch_cp_occ(&slots[s]);
    }
    // Any unused slots (numReads < N) start idle: DONE with no input, so the
    // stepping pass does nothing for them.
    for (int32_t s = initial; s < N; s++) {
        slots[s].phase = PH_DONE;
        slots[s].ready = false;
        slots[s].input_idx = -1;
    }

    /* A slot is refilled as soon as its read finishes, not when the read is
     * emitted: a read that finishes behind a slower, earlier one would otherwise
     * hold its slot idle until that read is out. Finished reads go through
     * `parked`, which emits them in input order (smem_parked.h), so matchArray is
     * unchanged; only which reads share a stepping pass changes, and a read's walk
     * depends on nothing but its own slot. Inputs start at most `window` reads
     * ahead of the oldest unemitted one, which bounds how many reads are parked. */
    static thread_local LockstepParked parked;
    parked.reset(LockstepParked::window_for(N));
    int32_t next_input = initial;
    int32_t idle       = N - initial;   // slots with no input
    int64_t numTotalSmem = *__numTotalSmem;

    while (parked.cursor() < numReads) {
        // --- Stepping pass: advance each slot by one phase-step. ---
        for (int32_t s = 0; s < N; s++) {
            BatchSlot *sl = &slots[s];
            switch (sl->phase) {
                case PH_FWD:
                    ls_advance_forward_step(sl, enc_qdb, s1_text);
                    // If forward just transitioned to PH_BWD_INIT, run the
                    // between-phases housekeeping immediately. prepare_backward
                    // may transition straight to PH_DONE (e.g. empty prev[]).
                    if (sl->phase == PH_BWD_INIT) {
                        ls_prepare_backward(sl);
                    }
                    break;
                case PH_S1_LF:
                    ls_s1_lf_step(sl, s1_pac, s1_l_pac);
                    break;
                case PH_S1_CMP:
                    ls_s1_compare(sl, enc_qdb, s1_pac, s1_l_pac);
                    ls_prepare_backward(sl);
                    break;
                case PH_FWD_TABLE:
                    ls_advance_forward_table(sl, enc_qdb, s1_text);
                    if (sl->phase == PH_BWD_INIT) {
                        ls_prepare_backward(sl);
                    }
                    break;
                case PH_BWD_INIT:
                    ls_prepare_backward(sl);
                    break;
                case PH_BWD:
                    ls_advance_backward_step(sl, enc_qdb, minSeedLen);
                    break;
                case PH_DONE:
                    break;
            }
            // Retire a finished read and refill the slot here, at the slot that
            // just finished (a separate pass over every slot costs more than the
            // idle steps it saves). A fresh read can be done at once (first base
            // non-ACGT), so loop.
            while (sl->phase == PH_DONE && sl->ready) {
                query_pos_array[sl->input_idx] = sl->next_x;
                parked.finish(sl->input_idx, sl->match_buf, sl->match_count,
                              matchArray, &numTotalSmem);
                if (next_input < numReads && parked.admits(next_input)) {
                    ls_init_slot(sl, next_input, query_pos_array, min_intv_array,
                                 rid_array, seq_, query_cum_len_ar, enc_qdb, bsd_hint);
                    if (sl->phase == PH_FWD) ls_prefetch_cp_occ(sl);
                    next_input++;
                } else {
                    sl->input_idx = -1;
                    sl->ready = false;
                    idle++;
                }
            }
            /* No cross-slot prefetch for another slot's next visit: each slot
             * already hints what its next step reads, at the end of its last one
             * (ls_advance_forward_step, ls_advance_forward_table,
             * ls_advance_backward_step, ls_prepare_backward, ls_s1_lf_step) or when
             * filled (ls_init_slot's jump-table entry, ls_prefetch_cp_occ). */
        }

        // --- Refill the slots the window turned away, now that it may admit
        // them. A fill that is done at once is retired on its slot's next visit.
        if (idle != 0) {
            for (int32_t s = 0; s < N && next_input < numReads && parked.admits(next_input); s++) {
                BatchSlot *sl = &slots[s];
                if (sl->input_idx >= 0) continue;
                ls_init_slot(sl, next_input, query_pos_array, min_intv_array,
                             rid_array, seq_, query_cum_len_ar, enc_qdb, bsd_hint);
                if (sl->phase == PH_FWD) ls_prefetch_cp_occ(sl);
                next_input++;
                idle--;
            }
        }
    }

    *__numTotalSmem = numTotalSmem;
    /* prev/match buffers are owned by the thread_local cache above; intentionally
     * not freed here so the next call (same thread) reuses them without a
     * round-trip through the allocator. They are released by
     * LockstepSmemCache's destructor at thread exit. */
}

int64_t FMI_search::bwtSeedStrategyAllPosOneThread(uint8_t *enc_qdb,
                                                   int32_t *max_intv_array,
                                                   int32_t numReads,
                                                   const bseq1_t *seq_,
                                                   int32_t *query_cum_len_ar,
                                                   int32_t minSeedLen,
                                                   SMEM *matchArray)
{
    int32_t i;

    int64_t numTotalSeed = 0;

    for(i = 0; i < numReads; i++)
    {
        int readlength = seq_[i].l_seq;
        // No walk from x_end on can emit (bwtseed_x_end), so the read ends there.
        const int x_end = bwtseed_x_end(readlength, minSeedLen);
        int32_t x = 0;
        while(x < x_end)
        {
            int next_x = x + 1;

            // Forward search
            SMEM smem;
            smem.info = 0;   // carries no SA tag (ls_s1_sa_tag)
            smem.rid = i;
            smem.m = x;
            smem.n = x;
            
            int offset = query_cum_len_ar[i];
            uint8_t a = enc_qdb[offset + x];
            // uint8_t a = enc_qdb[i * readlength + x];

            if(a < 4)
            {
                smem.k = count[a];
                smem.l = count[3 - a];
                smem.s = count[a+1] - count[a];


                int j;
                for(j = x + 1; j < readlength; j++)
                {
                    next_x = j + 1;
                    // a = enc_qdb[i * readlength + j];
                    a = enc_qdb[offset + j];
                    if(a < 4)
                    {
                        SMEM smem_ = smem;

                        // Forward extension is backward extension with the BWT of reverse complement
                        smem_.k = smem.l;
                        smem_.l = smem.k;
                        SMEM newSmem_ = backwardExt(smem_, 3 - a);
                        //SMEM smem = backwardExt(smem, 3 - a);
                        //smem.n = j;
                        SMEM newSmem = newSmem_;
                        newSmem.k = newSmem_.l;
                        newSmem.l = newSmem_.k;
                        newSmem.n = j;
                        smem = newSmem;
#ifdef ENABLE_PREFETCH
                        /* Same correction as the lockstep path: the next step reads
                         * sp = smem.l, ep = smem.l + smem.s after the k<->l swap.
                         * This is the x86 default round-3 path (the lockstep variant
                         * is arm64-gated), so it was missing the #242 fix entirely. */
                        _mm_prefetch((const char *)(&cp_occ[(smem.l) >> CP_SHIFT]), _MM_HINT_T0);
                        _mm_prefetch((const char *)(&cp_occ[(smem.l + smem.s) >> CP_SHIFT]), _MM_HINT_T0);
#endif


                        if((smem.s < max_intv_array[i]) && ((smem.n - smem.m + 1) >= minSeedLen))
                        {

                            if(smem.s > 0)
                            {
                                matchArray[numTotalSeed++] = smem;
                            }
                            break;
                        }
                    }
                    else
                    {

                        break;
                    }
                }

            }
            x = next_x;
        }
    }
    return numTotalSeed;
}

// Lockstep-batched bwtSeed driver. Drop-in for bwtSeedStrategyAllPosOneThread
// — same input contract, same emission order (reads in input order; within
// a read, outer-x ascending; at most one SMEM per outer x). The interleaving
// runs BWTSEED_LOCKSTEP_N reads' forward-extension walks together so the
// independent cp_occ cache-line misses issue concurrently into Zen 4's
// load queue rather than serializing on per-read load latency.
//
// Structurally parallel to getSMEMsOnePosOneThread_lockstep above; see that
// function's comments for the hybrid SoA per-slot/bulk-buffer layout
// rationale. Prefetch model: same-slot T0 only, at the end of each step
// (bsd_advance_step) and at slot init; as in the SMEM driver, there is no
// cross-slot T1 lookahead (see the stepping pass below).
int64_t FMI_search::bwtSeedStrategyAllPosOneThread_lockstep(uint8_t *enc_qdb,
                                                            int32_t *max_intv_array,
                                                            int32_t numReads,
                                                            const bseq1_t *seq_,
                                                            int32_t *query_cum_len_ar,
                                                            int32_t minSeedLen,
                                                            SMEM *matchArray,
                                                            int32_t max_readlength,
                                                            const BsdFirstWalkHint *bsd_hint)
{
    if (numReads <= 0) return 0;
    if (bsd_hint != nullptr && bsd_hint->min_len != minSeedLen) bsd_hint = nullptr;

    // Runtime lockstep width (g_bwtseed_lockstep_n, resolved once at startup from
    // a BWA3_BWTSEED_LOCKSTEP_N pin, else the host default: 32 on AMD Zen 5
    // (bwa3_bwtseed_host_default_width), the compile-time default elsewhere). Byte-identical
    // across widths -- batching only. The on-stack slot array is sized to the
    // compile-time BWTSEED_LOCKSTEP_N_MAX so the runtime width can exceed the default.
    const int32_t N = g_bwtseed_lockstep_n;
    // slots[] is a fixed BWTSEED_LOCKSTEP_N_MAX stack array indexed in [0, N);
    // the env parser and the host default clamp to that range, but a direct
    // g_bwtseed_lockstep_n write (e.g. a test) could exceed it -- guard the bound.
    xassert(N >= 1 && N <= BWTSEED_LOCKSTEP_N_MAX,
            "g_bwtseed_lockstep_n out of range [1, BWTSEED_LOCKSTEP_N_MAX]");

    // Per-thread cache for match_buf[] slices. One pointer (no prev[]
    // needed — bwtSeed has only a forward pass). Sized from max_readlength;
    // scalar's worst-case emit is one SMEM per x ∈ [0, readlength) so the
    // bound holds. Grown monotonically across calls.
    BwtSeedSlot slots[BWTSEED_LOCKSTEP_N_MAX] = {};
    static thread_local LockstepBwtSeedCache cache;
    const size_t per_slot_smems = (size_t)max_readlength;
    if (per_slot_smems > cache.per_slot || N > cache.n_slots) {
        if (cache.match != nullptr) _mm_free(cache.match);
        // Grow on EITHER dimension: g_bwtseed_lockstep_n (N) can rise after first
        // use (a BWA3_BWTSEED_LOCKSTEP_N pin / the host default, or a direct
        // write in a test), and
        // sizing only on per_slot would leave slots[N-1] pointing past the
        // allocation. Mirrors the phase-2 SMEM cache (LockstepSmemCache) above.
        const size_t alloc_slots = (size_t)(N > cache.n_slots ? N : cache.n_slots);
        const size_t alloc_per   = per_slot_smems > cache.per_slot ? per_slot_smems : cache.per_slot;
        const size_t total_bytes = alloc_slots * alloc_per * sizeof(SMEM);
        cache.match = (SMEM *)_mm_malloc(total_bytes, 64);
        assert_not_null(cache.match, total_bytes, total_bytes);
        cache.per_slot = alloc_per;
        cache.n_slots  = (int32_t)alloc_slots;
    }
    for (int32_t s = 0; s < N; s++) {
        slots[s].match_buf = cache.match + (size_t)s * cache.per_slot;
    }

    // Jump depth for this call: walks jump from their first base to length
    // jd = min(table depth, minSeedLen - 1) (bsd_plan_jump); 0 = step only
    // (no table, BWA3_BWTSEED_JUMP=0, or jd < 2, where there is nothing to skip).
    int32_t jd = 0;
    const SeedJumpEntry *jlev = nullptr;
    if (seed_jump != nullptr) {
        jd = seed_jump_depth < minSeedLen - 1 ? seed_jump_depth : minSeedLen - 1;
        if (jd >= 2) jlev = seed_jump + seed_jump_level_off(jd);
        else jd = 0;
    }

    // Seed the first min(N, numReads) slots.
    const int32_t initial = (numReads < N) ? numReads : N;
    for (int32_t s = 0; s < initial; s++) {
        bsd_init_slot(&slots[s], s, max_intv_array, seq_, query_cum_len_ar, enc_qdb, jlev, jd,
                      minSeedLen, bsd_hint);
        if (slots[s].phase == BSD_FWD && slots[s].jump == nullptr) bsd_prefetch_cp_occ(&slots[s]);
    }
    // Unused slots get input_idx = -1 so the flush pass skips them.
    for (int32_t s = initial; s < N; s++) {
        slots[s].phase     = BSD_DONE;
        slots[s].ready     = false;
        slots[s].input_idx = -1;
    }

    int32_t next_input   = initial;
    int32_t flush_cursor = 0;
    int32_t flush_slot   = 0;   // slot holding input flush_cursor (see the flush pass)
    int64_t numTotalSeed = 0;

    while (flush_cursor < numReads) {
        // Stepping pass: advance each active slot by one unit of work.
        for (int32_t s = 0; s < N; s++) {
            if (slots[s].phase == BSD_FWD) {
                bsd_advance_step(&slots[s], enc_qdb, minSeedLen, jlev, jd);
            }
            // No cross-slot prefetch for another slot's next visit: a slot that
            // steps next already hinted its two cp_occ lines (smem.l and
            // smem.l + smem.s) at the end of its last step (bsd_advance_step) or
            // at init; a slot with a planned jump next reads only its jump-table
            // entry, which bsd_plan_jump prefetched; and a fresh seed (j < 0)
            // reads lines that are not known until its next start x is found.
            // A second hint would only repeat the first or miss.
        }

        // Flush pass: emit retired slots in input order; recycle into the
        // next pending read. As in the SMEM driver, input i always lives in
        // slot i mod N (slot s is seeded with input s, and the slot freed by
        // input k is refilled with input k + N), so the slot holding input
        // `flush_cursor` is a counter that wraps at N instead of a scan.
        while (flush_cursor < numReads) {
            const int32_t s = flush_slot;
            if (!(slots[s].phase == BSD_DONE && slots[s].ready)) break;
            for (int32_t m = 0; m < slots[s].match_count; m++) {
                matchArray[numTotalSeed++] = slots[s].match_buf[m];
            }
            flush_cursor++;
            if (++flush_slot == N) flush_slot = 0;

            if (next_input < numReads) {
                bsd_init_slot(&slots[s], next_input,
                              max_intv_array, seq_, query_cum_len_ar, enc_qdb, jlev, jd,
                              minSeedLen, bsd_hint);
                if (slots[s].phase == BSD_FWD && slots[s].jump == nullptr)
                    bsd_prefetch_cp_occ(&slots[s]);
                next_input++;
            } else {
                slots[s].input_idx = -1;
                slots[s].ready     = false;
            }
        }
    }

    return numTotalSeed;
}


/* ---- Shared pieces of the SMEM rid counting sort (sortSMEMs and
 * smem_sort_rid_mn). ---- */

/* Grow a reused scratch buffer to hold at least `need` elements: free + fresh
 * 64-byte-aligned allocation, null-checked, never shrunk (so a buffer is
 * bounded by the largest batch it has seen). Contents are not preserved. */
template <typename T>
static T *smem_scratch_reserve(T *&buf, int64_t &cap, int64_t need)
{
    if (need > cap) {
        _mm_free(buf);
        buf = (T *) _mm_malloc((size_t)need * sizeof(T), 64);
        if (buf == NULL) {
            fprintf(stderr, "ERROR: out of memory in %s\n", __func__);
            exit(EXIT_FAILURE);
        }
        cap = need;
    }
    return buf;
}

/* rid range present in a[0, n), n >= 1. A per-thread block holds only its own
 * reads, so the dense rid index is offset by the observed minimum. */
static int64_t smem_rid_range(const SMEM *a, int64_t n, uint32_t *minRidOut)
{
    uint32_t minRid = a[0].rid, maxRid = a[0].rid;
    for (int64_t i = 1; i < n; i++) {
        uint32_t r = a[i].rid;
        if (r < minRid) minRid = r;
        if (r > maxRid) maxRid = r;
    }
    *minRidOut = minRid;
    return (int64_t)maxRid - (int64_t)minRid + 1;
}

/* Stable counting scatter of a[0, n) by rid: out[] receives make(i) for each
 * input i, grouped by ascending rid, keeping input order within a rid. cnt
 * must hold `range` entries; on return cnt[r] is the END offset of bucket r
 * (its start is cnt[r - 1], or 0 for r == 0). cnt is zeroed over [0, range)
 * here -- only the used prefix, not the whole capacity. */
template <typename T, typename Make>
static void smem_rid_scatter(const SMEM *a, int64_t n, uint32_t minRid, int64_t range,
                             int64_t *cnt, T *out, Make make)
{
    memset(cnt, 0, (size_t)range * sizeof(int64_t));
    for (int64_t i = 0; i < n; i++)
        cnt[a[i].rid - minRid]++;
    int64_t sum = 0;
    for (int64_t r = 0; r < range; r++) { int64_t c = cnt[r]; cnt[r] = sum; sum += c; }
    for (int64_t i = 0; i < n; i++) {
        int64_t b = (int64_t)a[i].rid - minRid;
        out[cnt[b]++] = make(i);
    }
}

void FMI_search::sortSMEMs(SMEM *matchArray,
        int64_t numTotalSmem[],
        int32_t numReads,
        int32_t readlength,
        int nthreads,
        SmemSortScratch &scratch)
{
    /* Stable counting sort of each per-thread block by rid: every SMEM of a
     * read (rid) becomes contiguous, reads in ascending rid order, and SMEMs of
     * one rid keep their incoming relative order. No (m, n) order is applied.
     *
     * The aligner no longer calls this (mem_collect_smem uses smem_sort_rid_mn,
     * which also orders each read by (m, n)); it is kept for the integration
     * tests, which diff its output, so its order must not change. The
     * smem_sort_rid_mn unit test pins it to a stable sort by rid.
     *
     * rid is a dense index in [0, numReads), so a counting sort by rid replaces
     * the old O(n log n) function-pointer qsort with an O(n + rid_range) pass.
     * The count/offset array (`cnt`) and the scatter buffer (`tmp`) live in
     * caller-owned scratch that is allocated once and reused across batches
     * (audit SEED-15). FMI_search is shared across worker threads, so no member
     * scratch may be used here; the scratch comes from the caller and is never
     * shared. */
    int tid;
    int32_t perThreadQuota = (numReads + (nthreads - 1)) / nthreads;
    for(tid = 0; tid < nthreads; tid++)
    {
        int64_t smem_count = numTotalSmem[tid];   /* scalar SMEM count; not the class `count[5]` array */
        if (smem_count <= 1) continue;   /* 0 or 1 element is already sorted */

        int32_t first = tid * perThreadQuota;
        SMEM *myMatchArray = matchArray + (int64_t)first * readlength;

        uint32_t minRid;
        int64_t range = smem_rid_range(myMatchArray, smem_count, &minRid);
        int64_t *cnt = smem_scratch_reserve(scratch.cnt, scratch.cntCap, range);
        SMEM    *tmp = smem_scratch_reserve(scratch.tmp, scratch.tmpCap, smem_count);

        smem_rid_scatter(myMatchArray, smem_count, minRid, range, cnt, tmp,
                         [myMatchArray](int64_t i) { return myMatchArray[i]; });
        /* tmp is distinct scratch storage and cannot overlap myMatchArray */
        memcpy(myMatchArray, tmp, (size_t)smem_count * sizeof(SMEM));
    }
}

/* (sort key, source index) pair for smem_sort_rid_mn; forward-declared as
 * SmemKeyIdx in FMI_search.h so the layout stays private to this file. */
struct smem_key_idx
{
    uint64_t key;   /* (uint64_t)m << 32 | n */
    uint64_t idx;   /* position of the SMEM in the unsorted array */
};

/* Per-rid order of mem_collect_smem: the key (m << 32 | n) ascending. This
 * MUST order exactly as bwamem.cpp's intv_lt1 (KSORT_INIT(mem_intv1)), which
 * defined the order before this sort replaced it; the unit test
 * test/unit/test_smem_sort_rid_mn.cpp enforces it against that comparator. */
#define smem_key_lt(a, b) ((a).key < (b).key)
KSORT_INIT(smem_key_idx, SmemKeyIdx, smem_key_lt)

void smem_sort_rid_mn(SMEM *a, int64_t n, SmemSortScratch &scratch)
{
    /* Same result as sortSMEMs (nthreads = 1) followed by a per-rid
     * ks_introsort(mem_intv1), by construction:
     *  - the stable counting scatter by rid is the same, except that it moves a
     *    16-byte (key, idx) pair instead of the 48-byte SMEM, so every rid
     *    bucket starts in the same incoming order;
     *  - each bucket is then sorted by the same ks_introsort algorithm with an
     *    identically-ordering comparator on the same keys, so it makes the same
     *    comparisons and applies the same swaps -- the same permutation, ties
     *    included;
     *  - each SMEM is then gathered into scratch and copied back (two 48-byte
     *    moves, against a scatter + copy-back + every sort shift before).
     * The rid buckets come straight from the counting offsets, so no rescan
     * for rid boundaries is needed either. */
    if (n <= 1) return;

    uint32_t minRid;
    int64_t range = smem_rid_range(a, n, &minRid);
    int64_t    *cnt = smem_scratch_reserve(scratch.cnt, scratch.cntCap, range);
    SMEM       *tmp = smem_scratch_reserve(scratch.tmp, scratch.tmpCap, n);
    SmemKeyIdx *p   = smem_scratch_reserve(scratch.pairs, scratch.pairsCap, n);

    smem_rid_scatter(a, n, minRid, range, cnt, p, [a](int64_t i) {
        SmemKeyIdx e;
        e.key = ((uint64_t)a[i].m << 32) | (uint64_t)a[i].n;
        e.idx = (uint64_t)i;
        return e;
    });

    int64_t start = 0;
    for (int64_t r = 0; r < range; r++) {
        int64_t end = cnt[r];
        if (end - start > 1)
            ks_introsort(smem_key_idx, (size_t)(end - start), p + start);
        start = end;
    }

    for (int64_t i = 0; i < n; i++)
        tmp[i] = a[p[i].idx];
    /* tmp is distinct scratch storage and cannot overlap a */
    memcpy(a, tmp, (size_t)n * sizeof(SMEM));
}


/* FMI_search::backwardExt is now defined inline in FMI_search.h
 * (issue #87) so the body inlines into all 9 hot callers, eliminating
 * the SysV-ABI struct-by-value pass and return-slot store that dominate
 * self-time on gcc 12+. */

/* Uncompressed SA accessors: these index sa_ms_byte[pos] / sa_ls_word[pos] by the
 * raw BWT row. Under SA_COMPRESSION (the shipped build) those arrays are sized to
 * reference_seq_len >> SA_COMPX, so a raw-row index reads ~SA_COMPX-fold out of
 * bounds. They are only correct for a full (uncompressed) SA, and every in-tree
 * caller is already guarded by `#if !SA_COMPRESSION` (see FMI_search.cpp's
 * sentinel scan and bwamem.cpp's per-read SA resolve), so fence the definitions to
 * match: under SA_COMPRESSION they must not be compiled or callable. The compressed
 * paths (get_sa_entry_compressed / get_sa_entries(..., tid) / get_sa_entries_prefetch)
 * below are the shipped equivalents. */
#if !SA_COMPRESSION
int64_t FMI_search::get_sa_entry(int64_t pos)
{
    int64_t sa_entry = sa_ms_byte[pos];
    sa_entry = sa_entry << 32;
    sa_entry = sa_entry + sa_ls_word[pos];
    return sa_entry;
}

void FMI_search::get_sa_entries(int64_t *posArray, int64_t *coordArray, uint32_t count, int32_t nthreads)
{
    uint32_t i;
// #pragma omp parallel for num_threads(nthreads)
    for(i = 0; i < count; i++)
    {
        /* Prefetch the SAL slot SAL_PFD iterations ahead. Both
         * sa_ms_byte and sa_ls_word are large random-access arrays
         * (separate cache lines), so issue two prefetches per slot. */
        if (i + SAL_PFD < count) {
            int64_t pf_pos = posArray[i + SAL_PFD];
            _mm_prefetch((const char *)(sa_ms_byte + pf_pos), _MM_HINT_T0);
            _mm_prefetch((const char *)(sa_ls_word + pf_pos), _MM_HINT_T0);
        }
        int64_t pos = posArray[i];
        int64_t sa_entry = sa_ms_byte[pos];
        sa_entry = sa_entry << 32;
        sa_entry = sa_entry + sa_ls_word[pos];
        coordArray[i] = sa_entry;
    }
}

void FMI_search::get_sa_entries(SMEM *smemArray, int64_t *coordArray, int32_t *coordCountArray, uint32_t count, int32_t max_occ)
{
    uint32_t i;
    int32_t totalCoordCount = 0;
    for(i = 0; i < count; i++)
    {
        int32_t c = 0;
        SMEM smem = smemArray[i];
        int64_t hi = smem.k + smem.s;
        int64_t step = (smem.s > max_occ) ? smem.s / max_occ : 1;
        int64_t j;
        for(j = smem.k; (j < hi) && (c < max_occ); j+=step, c++)
        {
            int64_t pos = j;
            /* Prefetch SAL_PFD iterations ahead. Both arrays sit on
             * separate cache lines, so issue both prefetches per slot. */
            int64_t pf_pos = pos + SAL_PFD * step;
            if (pf_pos < hi) {
                _mm_prefetch((const char *)(sa_ms_byte + pf_pos), _MM_HINT_T0);
                _mm_prefetch((const char *)(sa_ls_word + pf_pos), _MM_HINT_T0);
            }
            int64_t sa_entry = sa_ms_byte[pos];
            sa_entry = sa_entry << 32;
            sa_entry = sa_entry + sa_ls_word[pos];
            coordArray[totalCoordCount + c] = sa_entry;
        }
        coordCountArray[i] = c;
        totalCoordCount += c;
    }
}
#endif // !SA_COMPRESSION

// sa_compression
int64_t FMI_search::get_sa_entry_compressed(int64_t pos, int tid)
{
    if ((pos & sa_compx_mask) == 0) {
        
        #if  SA_COMPRESSION
        int64_t sa_entry = sa_ms_byte[pos >> sa_compx];
        #else
        int64_t sa_entry = sa_ms_byte[pos];     // simulation
        #endif
        
        sa_entry = sa_entry << 32;
        
        #if  SA_COMPRESSION
        sa_entry = sa_entry + sa_ls_word[pos >> sa_compx];
        #else
        sa_entry = sa_entry + sa_ls_word[pos];   // simulation
        #endif
        
        return sa_entry;        
    }
    else {
        // tprof[MEM_CHAIN][tid] ++;
        int64_t offset = 0; 
        int64_t sp = pos;
        while(true)
        {
            int64_t occ_id_pp_ = sp >> CP_SHIFT;
            int64_t y_pp_ = CP_BLOCK_SIZE - (sp & CP_MASK) - 1; 
            uint64_t *one_hot_bwt_str = cp_occ[occ_id_pp_].one_hot_bwt_str;
            uint8_t b;

            if((one_hot_bwt_str[0] >> y_pp_) & 1)
                b = 0;
            else if((one_hot_bwt_str[1] >> y_pp_) & 1)
                b = 1;
            else if((one_hot_bwt_str[2] >> y_pp_) & 1)
                b = 2;
            else if((one_hot_bwt_str[3] >> y_pp_) & 1)
                b = 3;
            else
                b = 4;

            if (b == 4) {
                return offset;
            }

            GET_OCC(sp, b, occ_id_sp, y_sp, occ_sp, one_hot_bwt_str_c_sp, match_mask_sp);

            sp = count[b] + occ_sp;
            
            offset ++;
            // tprof[ALIGN1][tid] ++;
            if ((sp & sa_compx_mask) == 0) break;
        }
        // assert((reference_seq_len >> sa_compx) - 1 >= (sp >> sa_compx));
        #if  SA_COMPRESSION
        int64_t sa_entry = sa_ms_byte[sp >> sa_compx];
        #else
        int64_t sa_entry = sa_ms_byte[sp];      // simultion
        #endif
        
        sa_entry = sa_entry << 32;

        #if  SA_COMPRESSION
        sa_entry = sa_entry + sa_ls_word[sp >> sa_compx];
        #else
        sa_entry = sa_entry + sa_ls_word[sp];      // simulation
        #endif
        
        sa_entry += offset;
        return sa_entry;
    }
}

void FMI_search::get_sa_entries(SMEM *smemArray, int64_t *coordArray, int32_t *coordCountArray, uint32_t count, int32_t max_occ, int tid)
{
    
    uint32_t i;
    int32_t totalCoordCount = 0;
    for(i = 0; i < count; i++)
    {
        int32_t c = 0;
        SMEM smem = smemArray[i];
        int64_t hi = smem.k + smem.s;
        int64_t step = (smem.s > max_occ) ? smem.s / max_occ : 1;
        int64_t j;
        for(j = smem.k; (j < hi) && (c < max_occ); j+=step, c++)
        {
            int64_t pos = j;
            int64_t sa_entry = get_sa_entry_compressed(pos, tid);
            coordArray[totalCoordCount + c] = sa_entry;
        }
        // coordCountArray[i] = c;
        *coordCountArray += c;
        totalCoordCount += c;
    }
}

// SA_COPMRESSION w/ PREFETCH
int64_t FMI_search::call_one_step(int64_t pos, int64_t &sa_entry, int64_t &offset,
                                  int drop_sentinel_offset)
{
    if ((pos & sa_compx_mask) == 0) {
        // A sampled row: either freshly staged (offset 0) or reached by an LF
        // step on this lane's previous visit (offset = steps walked), whose
        // sample lines the caller has since prefetched.
        sa_entry = sa_ms_byte[pos >> sa_compx];
        sa_entry = sa_entry << 32;
        sa_entry = sa_entry + sa_ls_word[pos >> sa_compx];
        sa_entry += offset;
        return 1;
    }
    else {
        // int64_t offset = 0; 
        int64_t sp = pos;

        int64_t occ_id_pp_ = sp >> CP_SHIFT;
        int64_t y_pp_ = CP_BLOCK_SIZE - (sp & CP_MASK) - 1; 
        uint64_t *one_hot_bwt_str = cp_occ[occ_id_pp_].one_hot_bwt_str;
        uint8_t b;

#if defined(__ARM_NEON) || defined(__aarch64__) || defined(APPLE_SILICON)
        /* The BWT symbol at sp is one-hot across the four bit-planes, and the
         * bases are close to uniformly distributed, so testing the planes one
         * branch at a time mispredicts on most steps. Test all four at once
         * and pick the symbol with a ternary over them; the empty case (no plane
         * set) is the sentinel row, kept as the one real branch since it is
         * nearly never taken. The occurrence count then reuses the k-only
         * backward step's NEON popcount, whose d-register loads keep the
         * masked word out of the general-purpose file. Symbol priority
         * (lowest plane wins) matches the plane-at-a-time test exactly. */
        const unsigned t0 = (unsigned)((one_hot_bwt_str[0] >> y_pp_) & 1);
        const unsigned t1 = (unsigned)((one_hot_bwt_str[1] >> y_pp_) & 1);
        const unsigned t2 = (unsigned)((one_hot_bwt_str[2] >> y_pp_) & 1);
        const unsigned t3 = (unsigned)((one_hot_bwt_str[3] >> y_pp_) & 1);
        if ((t0 | t1 | t2 | t3) == 0) {
            // Sentinel ($) row (no bit-plane set; the x86 arm's planes == 0):
            // its suffix-array value is 0 and we walked `offset` LF steps to
            // reach it, so the SA entry is 0 + offset. drop_sentinel_offset
            // (=> 0) reproduces bwa-mem2's dropped walk for --compat=bwa-mem2;
            // the default keeps the accumulated offset, matching the x86 arm
            // and the compressed sibling get_sa_entry_compressed.
            sa_entry = drop_sentinel_offset ? 0 : offset;
            return 1;
        }
        b = (uint8_t)(t0 ? 0 : t1 ? 1 : t2 ? 2 : 3);
        const int64_t occ_sp = cp_occ[occ_id_pp_].cp_count[b] +
            occ_popcount64(&one_hot_bwt_str[b], &one_hot_mask_array[sp & CP_MASK]);
#else
        /* As on NEON above: the symbol at sp is one-hot across the four
         * bit-planes and the bases are near uniform, so a plane-at-a-time
         * branch chain mispredicts on most LF steps, and in this
         * multi-lane pipelined resolver a mispredict also discards
         * the other lanes' in-flight work. Gather the four plane bits into
         * one 4-bit mask and take its lowest set bit, which matches the
         * plane-at-a-time test's priority (lowest plane wins). An empty
         * mask is the sentinel row, kept as the one real branch since it
         * is nearly never taken. */
        const unsigned planes =
              (unsigned)((one_hot_bwt_str[0] >> y_pp_) & 1)
            | ((unsigned)((one_hot_bwt_str[1] >> y_pp_) & 1) << 1)
            | ((unsigned)((one_hot_bwt_str[2] >> y_pp_) & 1) << 2)
            | ((unsigned)((one_hot_bwt_str[3] >> y_pp_) & 1) << 3);
        if (planes == 0) {
            // Sentinel ($) row: its suffix-array value is 0, and we have walked
            // `offset` LF steps to reach it, so the SA entry is 0 + offset. The
            // compressed sibling get_sa_entry_compressed returns `offset` here for
            // the same reason. bwa-mem2 set sa_entry = 0 here, dropping the walk
            // and reporting positions within the first `offset` bases of the
            // concatenated reference up to `offset` bases too far left; that is
            // kept only for --compat=bwa-mem2 (drop_sentinel_offset), whose
            // contract is to reproduce bwa-mem2's records.
            sa_entry = drop_sentinel_offset ? 0 : offset;
            return 1;
        }
        b = (uint8_t)__builtin_ctz(planes);

        GET_OCC(sp, b, occ_id_sp, y_sp, occ_sp, one_hot_bwt_str_c_sp, match_mask_sp);
#endif
        
        sp = count[b] + occ_sp;
        
        offset ++;
        // Do not read the sample here even when sp is a sampled row: the load
        // would be an un-prefetched miss into the multi-GB sample arrays. Hand
        // sp back as a continuation; the caller prefetches its sample lines and
        // the lane resolves it on its next visit (the sampled-row branch above).
        sa_entry = sp;
        return 0;
    } // else
}

/* Lane count of the pipelined compressed-SA resolver (get_sa_entries_prefetch);
 * see the BWA3_SA_LANES note there. The count only sets how many
 * checkpoint-block and SA-sample fetches are in flight; results are stored by
 * position and are byte-identical at every width.
 *
 * x86_64: 40. A 20-vs-40 A/B on a 5M-pair WGS slice (3 interleaved reps, user
 * CPU) measured -0.30% on AMD Zen3 (every rep faster) and no change on Intel
 * Sapphire Rapids. x86_64 has not been re-swept since call_one_step began
 * deferring the sampled-row read (below).
 *
 * Elsewhere: 48. Since call_one_step defers a sampled row's SA read to the
 * lane's next visit (prefetched), each visit keeps one more fetch in flight
 * and the walk gains from more lanes; the earlier 16..40 sweep, flat on
 * Graviton4, predates that. Re-swept on Graviton4 (clang 19, hg38, 3 rotated
 * reps per count, median MEM_SA vs 20 lanes): 5M-pair WGS slice at -t 16,
 * 12 lanes +35.0%, 28 -8.3%, 32 -9.4%, 40 -12.2%, 48 -12.8%, 64 -12.8%;
 * 200k pairs at -t 1, 40 -10.3%, 48 -11.2%, 64 -11.2%. The curve is flat from
 * 40 to 64 (48 was the lowest median at both thread counts, within the
 * rep-to-rep spread of 40 and 64). An A/B of 48 against 20 on the same slice,
 * 6 rotated pairs at -t 16: MEM_SA -13.3%, user CPU -0.59% and wall -0.84%
 * (user lower in 6/6). Only Graviton4 was measured; the value also applies to
 * every other non-x86 build, Apple Silicon included. */
#if defined(__x86_64__) || defined(_M_X64)
#define SA_RESOLVE_LANES_DEFAULT 40
#else
#define SA_RESOLVE_LANES_DEFAULT 48
#endif
#define SA_RESOLVE_LANES_MAX 64
/* The resolver's lane arrays (working_set/map_pos/offset) are sized to
 * SA_RESOLVE_LANES_MAX and indexed by sa_batch_size, which is clamped into
 * [1, SA_RESOLVE_LANES_MAX] for env overrides but defaults to
 * SA_RESOLVE_LANES_DEFAULT unclamped -- guard the default against silently
 * overrunning those stack arrays if it is ever raised past the array size. */
static_assert(SA_RESOLVE_LANES_DEFAULT >= 1 &&
              SA_RESOLVE_LANES_DEFAULT <= SA_RESOLVE_LANES_MAX,
              "SA_RESOLVE_LANES_DEFAULT must be within [1, SA_RESOLVE_LANES_MAX]");

/* Thread-local scratch for the pos_ar/map_ar staging buffers below. These were
 * two _mm_malloc/_mm_free per call, and this function runs once per read — so on
 * a WGS run that is tens of millions of aligned-allocation round-trips. The pair
 * is pure scratch (no state carried between calls), so a per-thread grow-only
 * holder reused across reads removes the churn. FMI_search is shared across
 * worker threads, hence thread_local (not a member). Freed on thread exit. */
namespace {
struct SaPrefetchScratch {
    // Only pos is staged now: the former map_ar array held map_ar[k] == k for
    // every entry (the staging loop writes map_ar[id] = totalCoordCount + c,
    // and totalCoordCount == id on entry to each SMEM while both advance in
    // lockstep with c), so the map index is just the buffer index and the
    // second allocation was dead. Downstream consumers use the index directly.
    int64_t *pos = nullptr;
    int64_t  cap = 0;
    void ensure(int64_t n) {
        if (n <= cap) return;
        _mm_free(pos);
        pos = (int64_t *) _mm_malloc((size_t)n * sizeof(int64_t), 64);
        // Fail loudly on OOM: without this, pos==NULL but cap=n reports capacity,
        // so the next call returns early and the staging loop writes through NULL.
        xassert(pos != NULL, "out of memory: SA prefetch staging buffer");
        cap = n;
    }
    ~SaPrefetchScratch() { _mm_free(pos); }
};
} // namespace

/* ---------------------------------------------------------------------------
 * Cross-read SMEM-interval (k,s) dedup (BWA3_KS_DEDUP).
 *
 * SA resolution is a pure function of (k, s, max_occ): an interval always
 * expands to the same coordinate sub-sequence (positions k, k+step, ... with
 * step=(s>max_occ)?s/max_occ:1, capped at min(s,max_occ)) and resolves to the
 * same coordinates. Different reads in the same SA-resolve chunk frequently
 * carry the SAME (k,s) (shared repeats / common k-mers; measured 17-28% dup).
 * So within a chunk we resolve each distinct (k,s) once (the "rep") and memcpy
 * its resolved slot-run to every duplicate's slot-run -- byte-identical, and
 * skipping the duplicate LF-walks is the whole saving.
 *
 * Modes: 'off' (never), 'on' (always, count>1), 'auto' (default). 'auto' is a
 * net-cycles adaptive controller (mirroring the shipped extension-DP job dedup
 * in read_memo.cpp): because the SA-resolve share is small, dedup bookkeeping
 * must cost less than the LF-walks it removes, so 'auto' MEASURES that trade
 * per chunk and latches ON/OFF from the observed net rather than a dup-rate
 * threshold. A plain on/off knob and optional stats remain available via
 * --ks-dedup / BWA3_KS_DEDUP_STATS. Alignment output is byte-identical in every
 * mode -- the flag only trades duplicate LF-walks for a per-chunk (k,s) dedup.
 * ------------------------------------------------------------------------- */
/* ---- mode + config (mirrors mem_dedup_configure / --dedup, bwamem.cpp) ---- */
enum { KS_DEDUP_OFF = 0, KS_DEDUP_ON = 1, KS_DEDUP_AUTO = 2 };

static int ks_parse_dedup_mode(const char *e)
{
    if (!strcmp(e, "0") || !strcmp(e, "off"))                 return KS_DEDUP_OFF;
    if (!strcmp(e, "1") || !strcmp(e, "on"))                  return KS_DEDUP_ON;
    if (!strcmp(e, "auto") || !strcmp(e, "2"))                return KS_DEDUP_AUTO;
    return -1;
}

/* mode -1 = unresolved (only valid before ks_dedup_configure() has run once);
 * every reader goes through ks_dedup_cfg(), which lazily resolves from env/
 * defaults if some entry point (unit binaries, library use) never called it. */
struct KsDedupCfg { int mode; double z; int64_t reprobe_positions; };
static KsDedupCfg g_ks_cfg = { -1, 2.0, 200000000 };   /* mode -1 = unresolved */

/* CLI(--ks-dedup) > env BWA3_KS_DEDUP > default 'auto'; fatal on a bad value.
 * Declared in FMI_search.h and called once from fastmap.cpp after getopt.
 * Mirrors mem_dedup_configure: a non-empty CLI value wins outright; otherwise
 * consult the env, distinguishing "unset" (default) from "set but empty" (an
 * explicit misconfiguration -> fatal). Expert knobs BWA3_KS_DEDUP_{Z,REPROBE}
 * are env-only, full-string validated. */
void ks_dedup_configure(const char *mode_arg)
{
    if (mode_arg && *mode_arg) {
        const int v = ks_parse_dedup_mode(mode_arg);
        if (v < 0) { fprintf(stderr, "ERROR: --ks-dedup: expected off|on|auto, got '%s'\n", mode_arg); exit(1); }
        g_ks_cfg.mode = v;
    } else {
        const char *m = getenv("BWA3_KS_DEDUP");
        if (m == NULL) g_ks_cfg.mode = KS_DEDUP_AUTO;    /* unset -> default */
        else {
            const int v = ks_parse_dedup_mode(m);         /* "" -> -1 -> fatal */
            if (v < 0) { fprintf(stderr, "ERROR: BWA3_KS_DEDUP: expected off|on|auto, got '%s'\n", m); exit(1); }
            g_ks_cfg.mode = v;
        }
    }
    const char *z = getenv("BWA3_KS_DEDUP_Z");           /* expert knob, env-only */
    if (z) {
        char *end = NULL; errno = 0;
        const double zv = strtod(z, &end);
        if (end == z || end == NULL || *end != '\0' || errno == ERANGE || !std::isfinite(zv) || zv <= 0.0) {
            fprintf(stderr, "ERROR: BWA3_KS_DEDUP_Z: must be a finite number > 0, got '%s'\n", z);
            exit(1);
        }
        g_ks_cfg.z = zv;
    }
    const char *r = getenv("BWA3_KS_DEDUP_REPROBE");     /* expert knob, env-only; positions, 0=off */
    if (r) {
        char *end = NULL; errno = 0;
        const long long rv = strtoll(r, &end, 10);
        if (end == r || end == NULL || *end != '\0' || errno == ERANGE || rv < 0) {
            fprintf(stderr, "ERROR: BWA3_KS_DEDUP_REPROBE: must be a non-negative integer number of positions (0 = off), got '%s'\n", r);
            exit(1);
        }
        g_ks_cfg.reprobe_positions = rv;
    }
}

static inline const KsDedupCfg &ks_dedup_cfg(void)
{
    static std::once_flag once;
    std::call_once(once, []() { if (g_ks_cfg.mode < 0) ks_dedup_configure(NULL); });
    return g_ks_cfg;
}
static inline int    ks_dedup_mode(void)          { return ks_dedup_cfg().mode; }
static inline double ks_net_z_crit(void)          { return ks_dedup_cfg().z; }
static inline int64_t ks_dedup_reprobe_positions(void) { return ks_dedup_cfg().reprobe_positions; }

/* ---------------------------------------------------------------------------
 * Net-cycles adaptive controller (--ks-dedup=auto). Decide ON/OFF from MEASURED
 * time, not a dup-rate threshold. get_sa_entries_prefetch can time BOTH sides
 * of the trade in one call, so this uses the single-observation net model (the
 * shipped extension-DP dedup, "#415"), not read_memo's cross-chunk A/B:
 *
 *   overhead_ns = dedup bookkeeping (build the (k,s) table + the scatter memcpy)
 *   c_pos       = resolve_kernel_ns / distinct_positions  (measured per-position
 *                 LF-walk cost on THIS host/index/tier)
 *   benefit_ns  = dup_positions * c_pos   (the walks we skipped)
 *   net_ns      = overhead_ns - benefit_ns          (< 0  =>  dedup won here)
 *
 * c_pos is measured, so the break-even self-calibrates -- no magic dup-rate
 * constant. While measuring we run dedup ON (real data is dup-rich, so ON is
 * the productive default and timing a chunk we already run is a few clock
 * reads), pool per-chunk net/position across threads (Welford), and latch when
 * a two-sided z-test clears |z| >= z_crit. Below NET_MIN_BATCHES we never latch
 * (CLT warmup); past NET_MAX_POS with no significance we keep ON (break-even).
 * ------------------------------------------------------------------------- */
enum { KS_NET_MEASURING = 0, KS_NET_LATCH_ON = 1, KS_NET_LATCH_OFF = 2 };
static std::atomic<int>      g_ks_net_state{KS_NET_MEASURING};
static const uint64_t        KS_NET_MIN_BATCHES = 24;         /* CLT floor before any latch */
static const uint64_t        KS_NET_MAX_POS     = 300000000;  /* give up measuring past this */
static const uint64_t        KS_NET_PROBE_MAX_POS = 150000000;/* re-probe cap -> keep incumbent */
static std::atomic<int>      g_ks_net_incumbent{0};   /* 0 = none (initial phase) */
static std::atomic<int64_t>  g_ks_net_since_latch{0}; /* positions since last latch */
static std::mutex   g_ks_net_mtx;
static uint64_t     g_ks_net_k    = 0;    /* instrumented chunks */
static uint64_t     g_ks_net_pos  = 0;    /* positions measured */
static double       g_ks_net_mean = 0.0;  /* Welford mean of net_ns/chunk */
static double       g_ks_net_m2   = 0.0;  /* Welford sum of squares */

/* Feed one measured chunk into the controller and latch when significant.
 * Incumbent-aware: the initial decision uses a plain z-test with the MAX-pos
 * break-even fallback; a re-probe confirms the incumbent at the ordinary z (or
 * keeps it on a probe-cap timeout), and only flips at a higher bar (z+1) so
 * noise near break-even cannot flap the decision. */
static inline void ks_net_observe(double net_ns, int64_t positions)
{
    if (positions <= 0) return;
    std::lock_guard<std::mutex> lk(g_ks_net_mtx);
    if (g_ks_net_state.load(std::memory_order_relaxed) != KS_NET_MEASURING) return;
    g_ks_net_pos += (uint64_t)positions;
    g_ks_net_k   += 1;
    const double d = net_ns - g_ks_net_mean;
    g_ks_net_mean += d / (double)g_ks_net_k;
    g_ks_net_m2   += d * (net_ns - g_ks_net_mean);
    const int incumbent = g_ks_net_incumbent.load(std::memory_order_relaxed);
    int decision = 0; const char *why = "";
    if (g_ks_net_k >= KS_NET_MIN_BATCHES) {
        const double var = g_ks_net_m2 / (double)(g_ks_net_k - 1);
        const double sd  = var > 0.0 ? sqrt(var) : 0.0;
        const double z   = sd > 0.0 ? g_ks_net_mean * sqrt((double)g_ks_net_k) / sd : 0.0;
        const int    fav = (g_ks_net_mean < 0.0) ? KS_NET_LATCH_ON : KS_NET_LATCH_OFF;
        const double zc  = ks_net_z_crit();
        if (incumbent == 0) {
            if      (fabs(z) >= zc)                 { decision = fav;             why = "measured"; }
            else if (g_ks_net_pos >= KS_NET_MAX_POS){ decision = KS_NET_LATCH_ON; why = "break-even, keeping default"; }
        } else if (fav == incumbent) {
            if      (fabs(z) >= zc)                       { decision = incumbent; why = "re-probe confirmed"; }
            else if (g_ks_net_pos >= KS_NET_PROBE_MAX_POS){ decision = incumbent; why = "re-probe inconclusive, kept"; }
        } else {
            if      (fabs(z) >= zc + 1.0)                 { decision = fav;       why = "re-probe FLIPPED"; }
            else if (g_ks_net_pos >= KS_NET_PROBE_MAX_POS){ decision = incumbent; why = "below flip margin, kept"; }
        }
    }
    if (decision) {
        g_ks_net_incumbent.store(decision, std::memory_order_relaxed);
        g_ks_net_since_latch.store(0, std::memory_order_relaxed);
        g_ks_net_state.store(decision, std::memory_order_relaxed);
        if (bwa_verbose >= 3)
            fprintf(stderr, "[ks-dedup-auto] %s: %llu pos / %llu chunks, mean net=%.1f us/chunk -> ks-dedup %s\n",
                    why, (unsigned long long)g_ks_net_pos, (unsigned long long)g_ks_net_k,
                    g_ks_net_mean / 1e3, decision == KS_NET_LATCH_ON ? "ON" : "OFF");
    }
}

/* Work-based re-probe: called from both latched paths with this chunk's
 * position count. When the cadence elapses, reset the accumulator and re-enter
 * MEASURING. Byte-identity is unaffected: a probe only changes whether the
 * dedup path or the plain path runs, both of which produce identical coords. */
static inline void ks_net_maybe_reprobe(int64_t positions)
{
    const int64_t period = ks_dedup_reprobe_positions();
    if (period <= 0) return;
    if (g_ks_net_since_latch.fetch_add(positions, std::memory_order_relaxed) + positions < period) return;
    /* Relaxed fast-path: once the counter passes `period` it STAYS past it until
     * a probe latches and resets it, so every eligible AUTO chunk would otherwise
     * take the global mutex only to hit the identical measuring-state check
     * below. Skip the lock while a probe is already running. The locked check is
     * retained because the state can change between here and the lock. */
    if (g_ks_net_state.load(std::memory_order_relaxed) == KS_NET_MEASURING) return;
    std::lock_guard<std::mutex> lk(g_ks_net_mtx);
    const int st = g_ks_net_state.load(std::memory_order_relaxed);
    if (st == KS_NET_MEASURING) return;                                   /* probe already running */
    if (g_ks_net_since_latch.load(std::memory_order_relaxed) < period) return;  /* lost the race */
    g_ks_net_k = 0; g_ks_net_pos = 0; g_ks_net_mean = 0.0; g_ks_net_m2 = 0.0;
    g_ks_net_since_latch.store(0, std::memory_order_relaxed);
    if (bwa_verbose >= 3)
        fprintf(stderr, "[ks-dedup-auto] re-probe (incumbent %s)\n", st == KS_NET_LATCH_ON ? "ON" : "OFF");
    g_ks_net_state.store(KS_NET_MEASURING, std::memory_order_relaxed);
}

static inline bool ks_dedup_stats_on(void)
{
    static const bool on = []() {
        const char *e = getenv("BWA3_KS_DEDUP_STATS");
        return e != NULL && *e != '\0' && !(e[0] == '0' && e[1] == '\0');
    }();
    return on;
}

/* Chunk-scoped (k,s) dedup table: open-addressing on a 64-bit mix of (k,s),
 * linear probe, holding a 1-based entry index into parallel key vectors so a
 * fresh (zeroed) slot reads as empty. Thread-local, grow-only, cleared per
 * call -- no allocation churn across chunks. */
namespace {
struct KsDedupTable {
    std::vector<int64_t>  k;        /* rep interval k, per entry            */
    std::vector<int64_t>  s;        /* rep interval s, per entry            */
    std::vector<int64_t>  base;     /* rep's coordArray slot-run start      */
    std::vector<uint32_t> slot;     /* hash slot -> 1-based entry idx, 0=empty */
    uint64_t              mask = 0; /* slot count - 1 (power of two)        */

    static inline uint64_t mix(int64_t kk, int64_t ss) {
        /* FNV-1a-style fold of the two 64-bit words; the table verifies exact
         * (k,s) on probe, so this only needs to spread well, not perfectly. */
        uint64_t h = 1469598103934665603ULL;
        h ^= (uint64_t)kk; h *= 1099511628211ULL;
        h ^= (uint64_t)ss; h *= 1099511628211ULL;
        return h;
    }

    void reset(int64_t n_hint) {
        k.clear(); s.clear(); base.clear();
        /* Load factor ~0.5: slots = next pow2 >= 2*n_hint (min 64). */
        uint64_t want = 64;
        while (want < (uint64_t)(n_hint * 2 + 1)) want <<= 1;
        if (want != mask + 1) { slot.assign(want, 0); mask = want - 1; }
        else                  { std::fill(slot.begin(), slot.end(), 0u); }
    }

    /* Return the rep entry index for (kk,ss); if absent, insert with rep_base
     * and return -1 (caller then stages this interval as the rep). */
    int64_t find_or_insert(int64_t kk, int64_t ss, int64_t rep_base) {
        uint64_t h = mix(kk, ss) & mask;
        for (;;) {
            uint32_t e = slot[h];
            if (e == 0) {                       /* empty -> insert new rep */
                k.push_back(kk); s.push_back(ss); base.push_back(rep_base);
                slot[h] = (uint32_t)k.size();   /* store 1-based */
                return -1;
            }
            if (k[e - 1] == kk && s[e - 1] == ss) return (int64_t)(e - 1);
            h = (h + 1) & mask;                 /* linear probe */
        }
    }
};
} // namespace

/* Stats: distinct-vs-total resolved POSITIONS (weighted by min(s,max_occ)),
 * summed across threads, printed once at exit. These count coordinates: an SMEM
 * whose coordinate comes from its SA tag (trust_s1_tag) counts as a distinct
 * position although nothing is walked for it. The AUTO controller's inputs
 * (ks_distinct_positions, ks_dup_positions) count LF walks instead, so the two
 * differ when tags are trusted; "saved" here is the share of coordinates that
 * dedup copies, not of LF walks. */
static std::atomic<uint64_t> g_ksdedup_total_pos{0};
static std::atomic<uint64_t> g_ksdedup_distinct_pos{0};

/* Test/observability hook (declared in FMI_search.h): read the process-global
 * (k,s)-dedup position counters. They only accumulate while BWA3_KS_DEDUP_STATS
 * is set (the fetch_adds are gated off the hot staging loop otherwise), so a
 * caller must enable that env before resolving to see nonzero counts. */
void ks_dedup_position_counts(uint64_t *total, uint64_t *distinct)
{
    if (total)    *total    = g_ksdedup_total_pos.load(std::memory_order_relaxed);
    if (distinct) *distinct = g_ksdedup_distinct_pos.load(std::memory_order_relaxed);
}

static struct KsDedupStatsDumper {
    ~KsDedupStatsDumper() {
        if (!ks_dedup_stats_on()) return;
        uint64_t t = g_ksdedup_total_pos.load(), d = g_ksdedup_distinct_pos.load();
        fprintf(stderr,
            "[ks-dedup-stats] total_pos=%llu distinct_pos=%llu saved=%.4f\n",
            (unsigned long long)t, (unsigned long long)d,
            t ? 1.0 - (double)d / (double)t : 0.0);
    }
} g_ksdedup_stats_dumper;

void FMI_search::get_sa_entries_prefetch(SMEM *smemArray, int64_t *coordArray,
                                         int64_t *coordCountArray, int64_t count,
                                         const int32_t max_occ, int tid, int64_t &id_,
                                         int drop_sentinel_offset, bool trust_s1_tag)
{
    // Internal invariant for direct callers. max_occ is a divisor below
    // (`step = s / max_occ`) and the dedup slot-run width `c = min(s, max_occ)`,
    // so a non-positive value is a divide-by-zero / meaningless copy length, not
    // a resolvable no-op. The `-c` option parser now rejects <= 0 fatally and
    // the public fmi_seed_sa_prefetch facade still no-ops on <= 0 for its own
    // documented contract (it returns before reaching here), so this guard only
    // catches a future direct caller passing a bad value -- fail loud rather than
    // silently produce no seeds. xassert survives an EXTRA_CXXFLAGS=-DNDEBUG opt-out build.
    xassert(max_occ > 0, "get_sa_entries_prefetch: max_occ must be > 0");

    // uint32_t i;
    // totalCoordCount and id (below) both count entries staged into the int64
    // coordArray/pos_ar/map_ar buffers and are paired with the int64 mem_lim,
    // so keep them int64 to avoid truncating the running offset.
    int64_t totalCoordCount = 0;
    // mem_lim is the exact number of entries the staging loop below writes:
    // each SMEM contributes min(smem.s, max_occ) (the inner loop stops at
    // c < max_occ). Summing min(s, max_occ) instead of the uncapped interval
    // size s both right-sizes the pos_ar/map_ar allocation -- repetitive seeds
    // have s up to the reference length but still write only max_occ entries --
    // and keeps every term bounded by max_occ so the int64 sum cannot run away.
    // id indexes those buffers and is paired with mem_lim, so it is int64 too.
    int64_t mem_lim = 0, id = 0;

    for(int i = 0; i < count; i++)
    {
        SMEM smem = smemArray[i];
        // SMEM.s is a public int64_t. A validated index never yields s < 0
        // (count[] monotonicity guarantees s = count[a+1] - count[a] >= 0), but
        // the field is caller-supplied, and a negative s would make this
        // mem_lim contribution negative (undersizing the scratch), the dedup
        // slot-run width `c` negative, and the scatter memcpy length wrap huge
        // as a size_t. Clamp it to 0 -- a negative/empty interval resolves
        // nothing -- so every downstream use of smem.s (hi, step, c, the (k,s)
        // key) is consistent. Done here and at each SMEM read below.
        if (smem.s < 0) smem.s = 0;
        mem_lim += (smem.s > max_occ) ? max_occ : smem.s;
    }

    /* Reuse the per-thread staging buffers, grown to the high-water mem_lim.
     * mem_lim == 0 leaves the pointers null, but the staging loop below writes
     * nothing in that case (id stays 0), so they are never dereferenced. */
    static thread_local SaPrefetchScratch t_sa;
    t_sa.ensure(mem_lim);
    int64_t *pos_ar = t_sa.pos;

    /* Decide dedup on/off for this chunk from the mode. OFF: never; ON: always
     * (count>1); AUTO: on while MEASURING or latched ON, off when latched OFF.
     * count<2 can't have a cross-read duplicate, so skip the machinery. */
    const int ks_mode = ks_dedup_mode();
    bool dedup;
    if (ks_mode == KS_DEDUP_OFF || count < 2) {
        dedup = false;
    } else if (ks_mode == KS_DEDUP_ON) {
        dedup = true;
    } else { /* AUTO */
        ks_net_maybe_reprobe(mem_lim);
        dedup = (g_ks_net_state.load(std::memory_order_relaxed) != KS_NET_LATCH_OFF);
    }
    /* While AUTO is measuring we run the dedup path AND time it, so the ON/OFF
     * decision comes from observed net (see ks_net_observe). Timers below span
     * only the dedup-specific work; the resolve pipeline is timed separately and
     * its per-position rate turns skipped duplicate positions into a benefit. */
    const bool ks_measure = (ks_mode == KS_DEDUP_AUTO && dedup
                        && g_ks_net_state.load(std::memory_order_relaxed) == KS_NET_MEASURING);
    std::chrono::steady_clock::time_point ks_t_stage0, ks_t_kern0, ks_t_scat0;
    double  ks_stage_ns = 0.0, ks_kernel_ns = 0.0, ks_scatter_ns = 0.0;
    int64_t ks_dup_positions = 0, ks_distinct_positions = 0;

    /* Copy jobs recorded during the dedup staging pass: after the resolve
     * pipeline fills every REP slot-run, each (dup_base, rep_base, n) is a
     * memcpy that reproduces the duplicate interval's coordinates byte-for-byte.
     * Thread-local + grow-only to avoid per-chunk allocation. */
    static thread_local std::vector<int64_t> ks_dup_base, ks_rep_base, ks_dup_cnt;
    static thread_local KsDedupTable ks_tab;
    /* Maps each staged pos_ar entry -> its destination coordArray slot. When
     * every SMEM is staged (the plain path with SA tags not trusted), staging
     * index == coordArray slot (1:1). The dedup path skips duplicate intervals
     * and, with SA tags trusted (use_tag below), both paths skip tagged SMEMs,
     * so pos_ar compacts and the two diverge: a staged position must land in
     * its OWN slot (a rep's slot-run base + k), not in the compacted pos_ar
     * index. The pipeline reads this to place each result. (Sized to mem_lim
     * like pos_ar; only used when dedup is on or SA tags are trusted.) */
    static thread_local std::vector<int64_t> ks_dst;
    /* A unique match whose text position the seeding path already resolved
     * (ls_s1_sa_tag, FMI_search.cpp's lockstep SMEM driver) gets its one
     * coordinate written here, at staging, and no LF walk. Only a caller whose
     * SMEMs all come from bwa-mem3's own seeding (which leaves info 0 on every
     * SMEM it does not tag) passes trust_s1_tag; the fmi_seed_sa_prefetch
     * facade takes SMEMs from external code and does not.
     * CONTRACT: trusting the tag is only correct if every SMEM producer on the
     * seeding path leaves info 0 on an SMEM it does not tag. Nothing checks it
     * here (a stale bit 63 on an s == 1 SMEM would be written as a wrong
     * coordinate), so a new or changed producer must zero info. The SMEM parity
     * tests check info == 0 on the untagged output of every producer
     * (getSMEMsOnePosOneThread and its lockstep driver, and both
     * bwtSeedStrategyAllPosOneThread drivers). The tag holds the
     * exact SA value, so it is used only at drop_sentinel_offset 0: under
     * --compat=bwa-mem2 a walk that reaches the sentinel row reports 0, and
     * how far a walk goes before it meets a sampled row is not bounded (the SA
     * is sampled by row). Skipping SMEMs breaks the plain path's staging index
     * == coordArray slot identity, so with the tag trusted every staged
     * position records its slot in ks_dst, as the dedup path does. */
    const bool use_tag = trust_s1_tag && drop_sentinel_offset == 0;
    const bool remap = dedup || use_tag;
    /* Grow-only (like the pos_ar scratch): only indices [0, id) are written
     * before they are read, so a plain resize every call would re-zero the tail
     * and thrash the allocator for no benefit. Keep the high-water capacity. */
    if (remap && ks_dst.size() < (size_t)mem_lim) ks_dst.resize((size_t)mem_lim);

    if (!dedup) {
        for(int i = 0; i < count; i++)
        {
            int32_t c = 0;
            SMEM smem = smemArray[i];
            if (smem.s < 0) smem.s = 0;   // clamp negative/empty interval (see mem_lim loop)
            if (use_tag) {
                if (ls_s1_tagged(smem)) {
                    coordArray[totalCoordCount] = ls_s1_sa_coord(smem);
                    *coordCountArray += 1;
                    totalCoordCount += 1;
                    continue;
                }
                int64_t hi = smem.k + smem.s;
                int64_t step = (smem.s > max_occ) ? smem.s / max_occ : 1;
                for(int64_t j = smem.k; (j < hi) && (c < max_occ); j+=step, c++)
                {
                    ks_dst[id] = totalCoordCount + c;
                    pos_ar[id++] = j;
                }
                *coordCountArray += c;
                totalCoordCount += c;
                continue;
            }
            int64_t hi = smem.k + smem.s;
            int64_t step = (smem.s > max_occ) ? smem.s / max_occ : 1;
            int64_t j;
            for(j = smem.k; (j < hi) && (c < max_occ); j+=step, c++)
            {
                int64_t pos = j;
                 pos_ar[id++]  = pos;
                // map_ar[k] == k (== id here), so the staging index is stored
                // implicitly; map_pos below reads the index directly.
                // int64_t sa_entry = get_sa_entry_compressed(pos, tid);
                // coordArray[totalCoordCount + c] = sa_entry;
            }
            //coordCountArray[i] = c;
            *coordCountArray += c;
            totalCoordCount += c;
        }
    } else {
        /* Dedup staging: each SMEM i owns the contiguous coordArray slot-run
         * [totalCoordCount, totalCoordCount + c). For a FIRST-seen (k,s) (a
         * "rep") we stage its positions into pos_ar exactly as the plain path
         * (its slots get resolved). For a DUPLICATE (k,s) we stage nothing --
         * skipping those LF-walks is the win -- and record a copy job that,
         * after resolution, replicates the rep's slot-run into this SMEM's. */
        if (ks_measure) ks_t_stage0 = std::chrono::steady_clock::now();
        ks_tab.reset(count);
        ks_dup_base.clear(); ks_rep_base.clear(); ks_dup_cnt.clear();
        uint64_t stat_total = 0, stat_distinct = 0, stat_dup_walk = 0;
        for(int i = 0; i < count; i++)
        {
            SMEM smem = smemArray[i];
            if (smem.s < 0) smem.s = 0;   // clamp negative/empty interval (see mem_lim loop)
            /* A tagged SMEM (see use_tag) keeps its rep/duplicate role: a rep
             * has its coordinate written here instead of staged, and a
             * duplicate copies its rep's slot as any duplicate does. */
            const bool tagged = use_tag && ls_s1_tagged(smem);
            int64_t hi = smem.k + smem.s;
            int64_t step = (smem.s > max_occ) ? smem.s / max_occ : 1;
            /* c is min(s, max_occ) -- the slot-run width. Compute it up front
             * (same rule the resolve/consumer loops use) so the rep bookkeeping
             * and the dup copy length agree exactly. s is clamped non-negative
             * above and max_occ > 0 is enforced at entry, so c is in [0, max_occ]
             * -- the scatter memcpy length can never be negative (the retained
             * internal defense against a wrapped size_t). */
            int32_t c = (smem.s > max_occ) ? max_occ : (int32_t)smem.s;
            const int64_t run_base = totalCoordCount;

            int64_t rep = ks_tab.find_or_insert(smem.k, smem.s, run_base);
            if (rep < 0) {
                /* rep: stage its positions, recording each staged entry's true
                 * coordArray destination (run_base + cc) so the pipeline writes
                 * results into this rep's own slot-run, not the compacted
                 * pos_ar index. */
                if (tagged) {
                    coordArray[run_base] = ls_s1_sa_coord(smem);
                } else {
                    int32_t cc = 0;
                    for(int64_t j = smem.k; (j < hi) && (cc < max_occ); j+=step, cc++) {
                        ks_dst[id] = run_base + cc;
                        pos_ar[id++] = j;
                    }
                }
                stat_distinct += (uint64_t)c;
            } else {
                /* duplicate: no staging; record the memcpy job. */
                ks_dup_base.push_back(run_base);
                ks_rep_base.push_back(ks_tab.base[rep]);
                ks_dup_cnt.push_back(c);
                if (!tagged) stat_dup_walk += (uint64_t)c;
            }
            *coordCountArray += c;
            totalCoordCount  += c;
            stat_total       += (uint64_t)c;
        }
        if (ks_measure) {
            ks_stage_ns = std::chrono::duration<double, std::nano>(
                              std::chrono::steady_clock::now() - ks_t_stage0).count();
            /* The controller prices LF walks: the positions walked (id; every
             * distinct position but a tagged rep's) and the duplicates whose
             * walk dedup skips (a tagged duplicate has none to skip). Without
             * tags these are stat_distinct and stat_total - stat_distinct. */
            ks_distinct_positions = id;
            ks_dup_positions      = (int64_t)stat_dup_walk;
        }
        if (ks_dedup_stats_on()) {
            g_ksdedup_total_pos.fetch_add(stat_total, std::memory_order_relaxed);
            g_ksdedup_distinct_pos.fetch_add(stat_distinct, std::memory_order_relaxed);
        }
    }
    
    /* Out-param is the resolved-COORDINATE count, not the compacted staging
     * count. These are equal only when every SMEM is staged (the plain path
     * with SA tags not trusted). The dedup path stages only reps into pos_ar
     * (local `id`), and with SA tags trusted both paths skip tagged SMEMs, so
     * `id` under-counts by the skipped duplicates and tagged SMEMs.
     * totalCoordCount advances by min(s,max_occ) for every SMEM in both paths,
     * so accumulating it keeps *id_ invariant to dedup and to tags -- the
     * fmi_seed_sa_prefetch facade documents *id as the coordinate count, and
     * external consumers must not see it shrink when 'auto'/'on' dedup engages.
     * (`id` stays the pipeline's compacted bound.) */
    id_ += totalCoordCount;

    /* Number of SA rows walked concurrently. Each lane's next checkpoint block
     * is prefetched when the lane is (re)staged, so the lane count sets how
     * many block fetches are in flight; it does not affect results, which are
     * stored by position. BWA3_SA_LANES overrides the default for tuning; out
     * of range values (below 1 or above SA_RESOLVE_LANES_MAX) are ignored. */
    static const int32_t sa_batch_size = []() {
        int32_t lanes = SA_RESOLVE_LANES_DEFAULT;
        const char *e = getenv("BWA3_SA_LANES");
        if (e != NULL) {
            // Strict parse: atoi accepts a numeric prefix ("1x" -> 1) and has
            // undefined behavior when the value overflows int. Match the
            // BWA3_LOAD_THREADS parse in index_load_threads -- reject empty,
            // trailing-garbage, and overflowed input, then keep values in
            // [1, SA_RESOLVE_LANES_MAX]. On anything else warn and keep the
            // default (ignore-and-continue, no hard exit).
            char *end = NULL;
            errno = 0;
            const long v = strtol(e, &end, 10);
            if (errno != 0 || end == e || *end != '\0' ||
                v < 1 || v > SA_RESOLVE_LANES_MAX) {
                fprintf(stderr, "[W::get_sa_entries_prefetch] ignoring BWA3_SA_LANES=%s (expected 1..%d)\n",
                        e, (int) SA_RESOLVE_LANES_MAX);
            } else {
                lanes = (int32_t) v;
            }
        }
        return lanes;
    }();
    int64_t working_set[SA_RESOLVE_LANES_MAX], map_pos[SA_RESOLVE_LANES_MAX];
    int64_t offset[SA_RESOLVE_LANES_MAX] = {-1};

    // Prefetch what a lane's next visit to row p reads: the SA sample lines
    // for a sampled row (including one an LF step just reached; call_one_step
    // defers that read to the next visit), else the row's checkpoint block.
    auto prefetch_row = [&](int64_t p) {
        if ((p & sa_compx_mask) == 0) {
            _mm_prefetch(&sa_ms_byte[p >> sa_compx], _MM_HINT_T0);
            _mm_prefetch(&sa_ls_word[p >> sa_compx], _MM_HINT_T0);
        }
        else
            _mm_prefetch(&cp_occ[p >> CP_SHIFT], _MM_HINT_T0);
    };

    if (ks_measure) ks_t_kern0 = std::chrono::steady_clock::now();
    int i = 0, j = 0;    
    while(i<id && j<sa_batch_size)
    {
        int64_t pos =  pos_ar[i];
        working_set[j] = pos;
        map_pos[j] = remap ? ks_dst[i] : i;   // dedup or a tag remaps to the slot
        offset[j] = 0;
        
        prefetch_row(pos);
        i++;
        j++;
    }
        
    // all_quit counts up to id (int64); lim stays int (bounded by sa_batch_size).
    int lim = j;
    int64_t all_quit = 0;
    while (all_quit < id)
    {
        
        for (int k=0; k<lim; k++)
        {
            int64_t sp = 0, pos = 0;
            bool quit;
            if (offset[k] >= 0) {
                quit = call_one_step(working_set[k], sp, offset[k], drop_sentinel_offset);
            }
            else
                continue;
            
            if (quit) {
                coordArray[map_pos[k]] = sp;
                all_quit ++;
                
                if (i < id)
                {
                    pos = pos_ar[i];
                    working_set[k] = pos;
                    map_pos[k] = remap ? ks_dst[i] : i;   // dedup or a tag remaps to the slot
                    i++;
                    offset[k] = 0;
                    
                    prefetch_row(pos);
                }
                else
                    offset[k] = -1;
            }
            else {
                working_set[k] = sp;
                prefetch_row(sp);
            }
        }
    }

    if (ks_measure)
        ks_kernel_ns = std::chrono::duration<double, std::nano>(
                           std::chrono::steady_clock::now() - ks_t_kern0).count();

    /* Dedup scatter pass: every REP slot-run is now fully resolved in
     * coordArray (reps are always staged, so their slots were written by the
     * pipeline above). Replicate each rep's coordinates into the duplicate
     * intervals' slot-runs. Byte-identical: the duplicate interval would have
     * resolved to exactly these coordinates in exactly this order. */
    if (dedup) {
        if (ks_measure) ks_t_scat0 = std::chrono::steady_clock::now();
        const size_t ndup = ks_dup_base.size();
        for (size_t d = 0; d < ndup; d++)
            memcpy(coordArray + ks_dup_base[d],
                   coordArray + ks_rep_base[d],
                   (size_t)ks_dup_cnt[d] * sizeof(int64_t));
        if (ks_measure)
            ks_scatter_ns = std::chrono::duration<double, std::nano>(
                                std::chrono::steady_clock::now() - ks_t_scat0).count();
    }

    /* Net-cycles decision (AUTO, measuring): overhead is the dedup-specific
     * bookkeeping (building the (k,s) table during staging + the scatter
     * memcpy); benefit is the LF-walk time the skipped duplicate positions
     * would have cost, at this chunk's MEASURED per-distinct-position resolve
     * rate. Feed the signed net (overhead - benefit) to the controller. Guard
     * on real signal (some distinct work resolved, some dup skipped, timers
     * nonzero) so a degenerate chunk doesn't feed noise into the z-test. */
    if (ks_measure && ks_distinct_positions > 0 && ks_dup_positions > 0
        && ks_kernel_ns > 0.0) {
        const double c_pos      = ks_kernel_ns / (double)ks_distinct_positions;
        const double benefit_ns = (double)ks_dup_positions * c_pos;
        const double overhead_ns = ks_stage_ns + ks_scatter_ns;
        ks_net_observe(overhead_ns - benefit_ns, ks_distinct_positions + ks_dup_positions);
    }
    /* pos_ar is the reused thread-local scratch — no free here. */
}
