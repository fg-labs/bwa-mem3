// Exercises PackedText::get_base / get_kmer on a deterministic pac buffer
// built with the same _set_pac convention used by bwa-mem3 index.
#include "packed_text.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>

#define _set_pac(pac, l, c) ((pac)[(l)>>2] |= (c)<<((~(l)&3)<<1))

// Always-on check that doesn't get elided under -DNDEBUG (which would
// strip plain assert() and turn this test into a no-op).
#define CHECK(cond) do { \
    if (!(cond)) { \
        std::fprintf(stderr, "CHECK failed at %s:%d: %s\n", __FILE__, __LINE__, #cond); \
        std::fflush(stderr); \
        return 1; \
    } \
} while (0)

int main() {
    // Text: ACGTACGTACGT... length 64. With A=0,C=1,G=2,T=3 the 2-bit code
    // at position l is just l mod 4 — same identity used at the get_base
    // assertion below, but stated independently of PackedText's internals.
    const int N = 64;
    uint8_t pac[(N + 3) / 4] = {};
    for (int l = 0; l < N; ++l) {
        _set_pac(pac, l, l & 3);
    }
    // Spill to tmp file and mmap.
    char path[] = "/tmp/packed_text_test_XXXXXX";
    int fd = mkstemp(path);
    CHECK(fd >= 0);
    ssize_t nw = write(fd, pac, sizeof(pac));
    CHECK(nw == (ssize_t)sizeof(pac));
    close(fd);

    PackedText pt(path, N);
    for (int l = 0; l < N; ++l) {
        int expected = l & 3;
        CHECK(pt.get_base(l) == expected);
    }
    // k-mer at pos 0, k=8: ACGT ACGT -> 0,1,2,3,0,1,2,3 packed MSB-first in a uint64:
    // bits 62-60: 0, 60-58: 1, 58-56: 2, 56-54: 3, ...
    uint64_t expected_kmer = 0;
    for (int i = 0; i < 8; ++i)
        expected_kmer |= (uint64_t)(i & 3) << (62 - 2 * i);
    CHECK(pt.get_kmer(0, 8) == expected_kmer);

    // k-mer crossing a byte boundary (pos 1, k=8): CGTA CGTA
    uint64_t expected_kmer_1 = 0;
    for (int i = 0; i < 8; ++i)
        expected_kmer_1 |= (uint64_t)((i + 1) & 3) << (62 - 2 * i);
    CHECK(pt.get_kmer(1, 8) == expected_kmer_1);

    // End-of-text: k-mer starting 4 bases from end should pad with 0s.
    uint64_t near_end = pt.get_kmer(N - 4, 8);
    // First 4 bases are real; last 4 are padded 0. Real: bases 60..63 -> 0,1,2,3.
    uint64_t expected_end = ((uint64_t)0 << 62) | ((uint64_t)1 << 60) |
                            ((uint64_t)2 << 58) | ((uint64_t)3 << 56);
    CHECK(near_end == expected_end);

    // Past-end k-mer (pos == N): every position is past the end, so the
    // short-circuit must produce zero. Same for any pos > N.
    CHECK(pt.get_kmer(N, 8) == 0);
    CHECK(pt.get_kmer(N + 100, 4) == 0);

    unlink(path);

    // Zero-byte pac: n_bases == 0 means the file has zero bytes, and
    // mmap(length=0) is EINVAL on Linux. PackedText must short-circuit
    // and instantiate cleanly with data_ == nullptr / n_ == 0.
    char zpath[] = "/tmp/packed_text_test_zero_XXXXXX";
    int zfd = mkstemp(zpath);
    CHECK(zfd >= 0);
    close(zfd);
    {
        PackedText empty(zpath, 0);
        CHECK(empty.get_kmer(0, 4) == 0);
    }
    unlink(zpath);

    printf("packed_text_test OK\n");
    return 0;
}
