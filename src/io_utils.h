#ifndef BWA_IO_UTILS_H
#define BWA_IO_UTILS_H

#include <cerrno>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <sys/types.h>
#include <unistd.h>

#include "utils.h"

// Single per-call byte cap shared by both the read and the write path. macOS
// rejects any pread()/pwrite() whose count exceeds INT_MAX (2GiB) with EINVAL
// outright -- it is not reported as a short transfer, so the retry loops never
// get a chance and the index build aborts. The meth doubled .pac for a
// human-scale genome is ~3.2GB, which trips this on the write side; the base
// doubled .pac (~1.6GB) stays under it, which is why only --meth failed. Linux
// glibc caps a single read/write at 0x7ffff000 but returns a short count the
// loops already handle. Clamping every request to 1GiB (well under INT_MAX)
// makes both platforms take the short-transfer path. This must NOT be
// "simplified" away. This is the ONE place the cap is defined; the read-side
// fmi_pread_request_size() delegates here too.
static const size_t IO_MAX_ONCE = (size_t)1 << 30;   /* 1GiB, well under INT_MAX */

// Bytes to hand a single pread()/pwrite() call: the smaller of what remains and
// `cap`. Shared by pwrite_request_size() and fmi_pread_request_size().
static inline size_t io_request_size(size_t remaining, size_t cap)
{
    return remaining > cap ? cap : remaining;
}

// Bytes to hand a single pwrite() call: the smaller of what remains and `cap`
// (the 1GiB IO_MAX_ONCE by default). Kept as a distinct, defaulted-cap entry
// point so the write-side tests can call it directly.
static inline size_t pwrite_request_size(size_t remaining, size_t cap = IO_MAX_ONCE)
{
    return io_request_size(remaining, cap);
}

// Test-only observability: counts the completed write chunks (one per positive
// pwrite(), not per EINTR retry) in the most
// recent pwrite_all()/pwrite_all_status() call on this thread, so a test can prove it
// actually split its buffer into chunks. A regression that drops the clamp
// would write the whole buffer in one pwrite() and the count would fall to 1.
// The single increment is negligible next to a 1GiB write, so production
// behaviour is unchanged. `inline` (not `static inline`) so every translation
// unit shares one counter -- the test reads the same one pwrite_all() bumps.
inline unsigned long& pwrite_all_chunk_counter()
{
    static thread_local unsigned long count = 0;
    return count;
}
inline void          pwrite_all_reset_chunk_count() { pwrite_all_chunk_counter() = 0; }
inline unsigned long pwrite_all_chunk_count()       { return pwrite_all_chunk_counter(); }

// Write the entire `len`-byte buffer at file offset `off`, retrying on EINTR
// and looping on short writes. Both are permitted by POSIX; treating either as
// a hard failure can turn a transient signal into a spurious index-build abort.
// Each request is clamped to `max_chunk` (see pwrite_request_size) so a buffer
// larger than 2GiB does not trip the macOS single-write EINVAL cap. `max_chunk`
// defaults to the real 1GiB cap and is a parameter only so a test can force the
// multi-chunk path with a small cap instead of a real >2GiB buffer.
//
// Returns 0 on success, -1 for a `pwrite` that returned 0 (an unexpected
// zero-byte short write), or a positive errno on a write error. This is the
// shared core: the fatal `pwrite_all` wrapper below aborts on a non-zero return,
// while callers inside an OpenMP region (where err_fatal()->exit() would tear
// down libomp with sibling threads mid-loop) call this directly, record the
// first non-zero return, leave the region, and err_fatal outside it.
static inline int pwrite_all_status(int fd, const void* buf, size_t len, off_t off,
                                    size_t max_chunk = IO_MAX_ONCE)
{
    const uint8_t* p = static_cast<const uint8_t*>(buf);
    size_t remaining = len;
    pwrite_all_reset_chunk_count();
    while (remaining > 0) {
        ssize_t w = pwrite(fd, p, pwrite_request_size(remaining, max_chunk), off);
        if (w < 0) {
            if (errno == EINTR) continue;
            return errno ? errno : EIO;
        }
        if (w == 0) return -1;   // pwrite returned 0 (distinct from an errno failure)
        ++pwrite_all_chunk_counter();   // one bump per completed write chunk (not per
                                        // EINTR retry, which does not advance the buffer)
        p         += (size_t)w;
        remaining -= (size_t)w;
        off       += (off_t)w;
    }
    return 0;
}

// Fatal wrapper: aborts via err_fatal on any failure. `what` names the buffer.
static inline void pwrite_all(int fd, const void* buf, size_t len, off_t off,
                              const char* what, size_t max_chunk = IO_MAX_ONCE)
{
    int rc = pwrite_all_status(fd, buf, len, off, max_chunk);
    if (rc < 0) err_fatal("pwrite_all", "pwrite(%s) returned 0", what);
    if (rc > 0) err_fatal("pwrite_all", "pwrite(%s) failed: %s", what, strerror(rc));
}

/* ------------------------------------------------------------------------- *
 * Parallel index-array read (the read-side counterpart of the pwrite family
 * above). Definitions live in FMI_search.cpp; declared here — the shared
 * low-level IO header — so index-array and index-element loaders both reach
 * them without depending on the derived FM-index header. FMI_PREAD_MIN_CHUNK
 * (the worker-count floor) stays with the definition in FMI_search.h.
 * ------------------------------------------------------------------------- */

/* Number of workers to split an `nbytes` index-array read across, given the
 * caller's requested `nthreads`. Never returns more than `nthreads` (when
 * positive) nor so many that a chunk would fall below FMI_PREAD_MIN_CHUNK,
 * except the unavoidable single-worker case where `nbytes` is itself below the
 * floor. Non-positive `nthreads` clamps UP to 1; the result is always >= 1.
 * Exposed so the chunk arithmetic is unit-testable without a real index. */
int fmi_pread_worker_count(size_t nbytes, int nthreads);

/* Bytes to request from a single pread() call, given how many remain in this
 * worker's chunk. macOS fails a pread() whose count exceeds INT_MAX with EINVAL,
 * so this clamps to IO_MAX_ONCE (delegating to io_request_size above). Exposed
 * so the clamp is unit-testable without materialising a multi-GB file. */
size_t fmi_pread_request_size(size_t remaining);

/* Read the next `nbytes` of `fp` into `dst` using up to `nthreads` pread
 * workers, then leave the stream positioned exactly past them so a following
 * sequential read still lands correctly. Aborts the process on a read error or
 * short file. */
void fmi_pread_from_stream(FILE *fp, void *dst, size_t nbytes, int nthreads);

/* Worker count for the index load: the caller's request clamped to [1, 8]
 * (bandwidth-bound past ~8), overridable via BWA3_LOAD_THREADS (clamped to 64).
 * A malformed override is warned about once and ignored. Always returns >= 1.
 * Exposed so the fail-closed env-parse is unit-testable without an index. */
int index_load_threads(int n_threads);

#endif
