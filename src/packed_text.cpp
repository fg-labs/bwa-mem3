#include "packed_text.h"
#include "utils.h"
#include <cassert>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>

PackedText::PackedText(const std::string& pac_path, int64_t n_bases)
{
    fd_ = open(pac_path.c_str(), O_RDONLY | O_CLOEXEC);
    if (fd_ < 0)
        err_fatal(__func__, "failed to open '%s': %s", pac_path.c_str(), strerror(errno));
    struct stat st;
    if (fstat(fd_, &st) != 0)
        err_fatal(__func__, "failed to fstat '%s': %s", pac_path.c_str(), strerror(errno));
    map_bytes_ = (size_t)st.st_size;
    size_t expected_bytes = (size_t)((n_bases + 3) >> 2);
    if (map_bytes_ < expected_bytes)
        err_fatal(__func__, "pac file '%s' is too small: expected at least %zu bytes, got %zu",
                  pac_path.c_str(), expected_bytes, map_bytes_);
    if (map_bytes_ == 0) {
        // mmap() with length=0 is EINVAL on Linux; short-circuit so a
        // zero-base / zero-byte pac instantiates cleanly.
        data_ = nullptr;
    } else {
        void* p = mmap(nullptr, map_bytes_, PROT_READ, MAP_PRIVATE, fd_, 0);
        if (p == MAP_FAILED)
            err_fatal(__func__, "failed to mmap '%s': %s", pac_path.c_str(), strerror(errno));
        // Both documented hot consumers (emit_0123, compute_counts, write_doubled_pac)
        // stream the mapping linearly, so the kernel can prefetch aggressively.
        (void)posix_madvise(p, map_bytes_, POSIX_MADV_SEQUENTIAL);
        data_ = (const uint8_t*)p;
    }
    n_ = n_bases;
}

PackedText::~PackedText()
{
    if (data_) munmap((void*)data_, map_bytes_);
    if (fd_ >= 0) close(fd_);
}

uint64_t PackedText::get_kmer(int64_t pos, int k) const
{
    // Runtime guards (not just asserts): k > 32 would compute a negative
    // shift (62 - 2*i) and trigger UB; k <= 0 trivially returns 0. Both
    // also kept as asserts so debug builds catch the misuse loudly while
    // release builds (-DNDEBUG) still degrade gracefully.
    assert(k >= 0 && k <= 32);
    assert(pos >= 0);
    if (k <= 0) return 0;
    if (k > 32) k = 32;
    if (pos < 0) return 0;
    uint64_t v = 0;
    for (int i = 0; i < k; ++i) {
        int64_t p = pos + i;
        // Once we walk past the end, every later position is also past
        // the end and contributes a zero base, so stop calling get_base.
        if (p >= n_) break;
        v |= (uint64_t)get_base(p) << (62 - 2 * i);
    }
    return v;
}
