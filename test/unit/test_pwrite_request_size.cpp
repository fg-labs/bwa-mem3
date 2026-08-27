// test/unit/test_pwrite_request_size.cpp
//
// Per-call size clamp for pwrite_all, mirroring the read-side clamp exercised
// in test_fmi_pread_request_size.cpp.
//
// pwrite_all loops pwrite() until the whole buffer is written, retrying on
// EINTR and short writes. Handing the whole remaining buffer to a single
// pwrite() is fine on Linux but fails on macOS, where a count greater than
// INT_MAX returns EINVAL outright -- it is not a short write, so the loop never
// gets a chance and the build aborts. The meth doubled .pac for a human-scale
// genome is ~3.2GB, which trips this; the base doubled .pac (~1.6GB) stays
// under it, which is why only --meth failed. pwrite_request_size clamps every
// request so the requested count is always a legal pwrite() count on both
// platforms.

#include "doctest/doctest.h"
#include "../../src/io_utils.h"

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <climits>
#include <cstring>
#include <string>
#include <vector>
#include <fcntl.h>
#include <unistd.h>

TEST_CASE("pwrite_request_size never asks pwrite() for more than INT_MAX"
          * doctest::test_suite("unit/pwrite_request_size"))
{
    const size_t int_max = (size_t)INT_MAX;

    const size_t sizes[] = {
        0,
        1,
        (size_t)1 << 20,                 // 1MiB
        int_max - 1,
        int_max,
        int_max + 1,                     // first illegal count on macOS
        (size_t)3435973836u,             // ~3.2GB meth doubled .pac
        (size_t)20912755063u,            // a whole human-scale --meth index
    };
    for (size_t n : sizes)
        CHECK(pwrite_request_size(n) <= int_max);
}

TEST_CASE("pwrite_request_size clamps to exactly IO_MAX_ONCE, not just under INT_MAX"
          * doctest::test_suite("unit/pwrite_request_size"))
{
    // Pin the *value*, not merely `<= INT_MAX`: a regression that silently
    // changed the cap to any other sub-INT_MAX value (e.g. 512MiB, or the whole
    // INT_MAX) would still satisfy the boundary test above but ship a different
    // chunk size. IO_MAX_ONCE is the one shared 1GiB cap (io_utils.h).
    CHECK(IO_MAX_ONCE == (size_t)1 << 30);

    // Anything at or below the cap is returned whole.
    CHECK(pwrite_request_size(0) == 0);
    CHECK(pwrite_request_size(1) == 1);
    CHECK(pwrite_request_size(IO_MAX_ONCE - 1) == IO_MAX_ONCE - 1);
    CHECK(pwrite_request_size(IO_MAX_ONCE) == IO_MAX_ONCE);

    // Anything above the cap is clamped to exactly the cap.
    CHECK(pwrite_request_size(IO_MAX_ONCE + 1) == IO_MAX_ONCE);
    CHECK(pwrite_request_size((size_t)3435973836u) == IO_MAX_ONCE);   // ~3.2GB
    CHECK(pwrite_request_size((size_t)20912755063u) == IO_MAX_ONCE);  // whole --meth index
}

TEST_CASE("pwrite_request_size makes progress and never over-writes"
          * doctest::test_suite("unit/pwrite_request_size"))
{
    // Must never return more than remains, or the loop would advance past the
    // buffer and write past its end.
    CHECK(pwrite_request_size(0) == 0);
    CHECK(pwrite_request_size(1) == 1);
    CHECK(pwrite_request_size(4096) == 4096);

    // Must always make progress on a huge buffer, otherwise the loop spins
    // forever instead of completing -- and each step advances by exactly the cap.
    CHECK(pwrite_request_size((size_t)20912755063u) == IO_MAX_ONCE);

    // A buffer larger than the clamp is consumed in a finite number of calls,
    // each returning exactly the cap until less than a cap remains.
    size_t remaining = (size_t)3435973836u, calls = 0;
    while (remaining > 0) {
        size_t want = pwrite_request_size(remaining);
        REQUIRE(want > 0);
        REQUIRE(want <= remaining);
        CHECK(want == (remaining > IO_MAX_ONCE ? IO_MAX_ONCE : remaining));
        remaining -= want;
        if (++calls > 64) break;         // guard against a non-progressing clamp
    }
    CHECK(remaining == 0);
    CHECK(calls <= 64);
}

namespace {

// Open an anonymous temp fd: honour $TMPDIR (fall back to /tmp) like the
// sibling read-side test so CI and macOS use the per-user temp dir, and unlink
// immediately so nothing lingers even if a later REQUIRE aborts the case.
int open_anon_tmp_fd()
{
    const char* dir = getenv("TMPDIR");
    std::string tmpl = std::string(dir != NULL ? dir : "/tmp");
    if (!tmpl.empty() && tmpl[tmpl.size() - 1] != '/') tmpl += '/';
    tmpl += "bwa_pwrite_all_test_XXXXXX";
    std::vector<char> tmpl_buf(tmpl.begin(), tmpl.end());
    tmpl_buf.push_back('\0');
    int fd = mkstemp(tmpl_buf.data());
    if (fd >= 0) unlink(tmpl_buf.data());
    return fd;
}

// Write `src` at `off` via pwrite_all with the given per-call cap, then read it
// back and confirm the file length and every byte survive. A real I/O error in
// pwrite_all calls err_fatal, which exits the process -- so a failure here is a
// hard binary exit, not a localized CHECK failure; that is inherent to testing
// a function that fatals internally.
void check_pwrite_all_round_trip(const std::vector<uint8_t>& src, off_t off, size_t max_chunk)
{
    int fd = open_anon_tmp_fd();
    REQUIRE(fd >= 0);

    const size_t n = src.size();
    pwrite_all(fd, src.data(), n, off, "test buffer", max_chunk);

    // Observability: pwrite_all must have actually split the buffer into
    // ceil(n / max_chunk) chunks. Without this, a regression that dropped the
    // clamp and wrote the whole buffer in one pwrite() would still round-trip
    // green (the test buffers are < INT_MAX); the counter is what catches it.
    // Regular-file pwrite() writes the full requested count, so iterations
    // equal the chunk count exactly.
    const size_t expected_chunks = n == 0 ? 0 : (n + max_chunk - 1) / max_chunk;
    CHECK(pwrite_all_chunk_count() == expected_chunks);

    // Length: the file must be exactly off + n bytes.
    off_t end = lseek(fd, 0, SEEK_END);
    REQUIRE(end == off + (off_t)n);

    // Content: read the region back and compare byte-for-byte.
    std::vector<uint8_t> dst(n, 0);
    size_t done = 0;
    while (done < n) {
        ssize_t r = pread(fd, dst.data() + done, n - done, off + (off_t)done);
        REQUIRE(r > 0);
        done += (size_t)r;
    }
    close(fd);

    CHECK(done == n);
    CHECK(memcmp(src.data(), dst.data(), n) == 0);
}

std::vector<uint8_t> make_pattern(size_t n)
{
    std::vector<uint8_t> src(n);
    for (size_t i = 0; i < n; i++) src[i] = (uint8_t)((i * 2654435761u) >> 13);
    return src;
}

}  // namespace

TEST_CASE("pwrite_all round-trips a buffer byte-identically at a non-zero offset"
          * doctest::test_suite("unit/pwrite_request_size"))
{
    // Single-chunk functional coverage: a modest buffer, the real 1GiB cap
    // (so the helper's chunk-count check expects exactly one pwrite()).
    check_pwrite_all_round_trip(make_pattern((size_t)5 << 20), 1234567, IO_MAX_ONCE);
}

TEST_CASE("pwrite_all advances p/off/remaining across multiple clamped chunks"
          * doctest::test_suite("unit/pwrite_request_size"))
{
    // The behaviour the fix actually introduces is the multi-chunk loop: a
    // buffer larger than the per-call cap must be written in several pwrite()
    // calls with p, off, and remaining advancing in lockstep. Force that path
    // with a small cap over a small buffer instead of a real >2GiB write, so a
    // transposition bug in the advance (e.g. bumping off but not p, or clamping
    // the wrong quantity) is caught. cap*3 + 123 forces four iterations, the
    // last partial; the offset is deliberately non-aligned to the cap.
    const size_t cap = 4096;
    const size_t n   = cap * 3 + 123;
    check_pwrite_all_round_trip(make_pattern(n), (off_t)(cap * 5 + 7), cap);
}
