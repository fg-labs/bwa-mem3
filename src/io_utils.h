#ifndef BWA_IO_UTILS_H
#define BWA_IO_UTILS_H

#include <cerrno>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <sys/types.h>
#include <unistd.h>

#include "utils.h"

// pwrite the entire `len`-byte buffer at file offset `off`, retrying on
// EINTR and looping on short writes. Both are permitted by POSIX; treating
// either as a hard failure can turn a transient signal into a spurious
// index-build abort.
static inline void pwrite_all(int fd, const void* buf, size_t len, off_t off,
                              const char* what)
{
    const uint8_t* p = static_cast<const uint8_t*>(buf);
    size_t remaining = len;
    while (remaining > 0) {
        ssize_t w = pwrite(fd, p, remaining, off);
        if (w < 0) {
            if (errno == EINTR) continue;
            err_fatal("pwrite_all", "pwrite(%s) failed: %s", what, strerror(errno));
        }
        if (w == 0) err_fatal("pwrite_all", "pwrite(%s) returned 0", what);
        p         += (size_t)w;
        remaining -= (size_t)w;
        off       += (off_t)w;
    }
}

#endif
