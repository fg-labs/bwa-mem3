// test/unit/test_neon_soa_pack.cpp
//
// neon_soa_pack16 byte-identity vs the scalar SoA fill it replaced (the NEON
// counterpart of test/x86_soa_pack_test.cpp). The 16-lane 8-bit kernels
// (banded extension, kswv mate rescue) read their sequences as SoA rows; the
// contract is "byte-for-byte what the scalar fill wrote", pad rows included,
// and a SAM comparison is blind to pad bytes. Every lane length from 0 through
// a few tiles is covered, with the padA run ending before, at, inside and past
// the boundary tile, sequences shorter than one tile (the two-load path), the
// ambiguity remap, pad bytes equal to the remapped code, and canary rows past
// nrows. The pack must never read outside a lane's sequence (a lane's bytes may
// end at the end of a mapping): every group is packed twice, once with each
// sequence flush against an inaccessible page after it and once flush against
// one before it, so an over- or under-read faults the test.

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <random>
#include <vector>

#include <sys/mman.h>
#include <unistd.h>

#include "doctest/doctest.h"

#if defined(__ARM_NEON) || defined(__aarch64__)

#ifndef AMBIG_
#define AMBIG_ 4
#endif
#ifndef AMBQ
#define AMBQ 8
#endif
#include "bandedSWA.h"   // SIMD_WIDTH8 / SIMD_WIDTH16
#include "neon_soa_pack.h"

namespace {

// A sequence copy placed flush against a PROT_NONE guard page: after its last
// byte (guardAfter) or before its first byte. Reading outside [data, data + len)
// faults.
class GuardedSeq {
  public:
    GuardedSeq(const std::vector<uint8_t> &bytes, bool guardAfter) : len_(bytes.size()) {
        page_ = (size_t)sysconf(_SC_PAGESIZE);
        const size_t body = ((len_ + page_ - 1) / page_) * page_ + page_;   // at least one page
        map_size_ = body + 2 * page_;
        void *m = mmap(nullptr, map_size_, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANON, -1, 0);
        REQUIRE(m != MAP_FAILED);
        map_ = (uint8_t *)m;
        REQUIRE(mprotect(map_, page_, PROT_NONE) == 0);
        REQUIRE(mprotect(map_ + page_ + body, page_, PROT_NONE) == 0);
        data_ = guardAfter ? map_ + page_ + body - len_ : map_ + page_;
        if (len_) memcpy(data_, bytes.data(), len_);
    }
    ~GuardedSeq() { munmap(map_, map_size_); }
    GuardedSeq(const GuardedSeq &) = delete;
    GuardedSeq &operator=(const GuardedSeq &) = delete;
    const uint8_t *data() const { return len_ ? data_ : nullptr; }

  private:
    size_t len_, page_ = 0, map_size_ = 0;
    uint8_t *map_ = nullptr, *data_ = nullptr;
};

// The per-lane contract the scalar loops implemented: real bases (4 -> 8 when
// remap4to8) for k < len, padA on [len, padStart), padB from padStart on.
void scalar_fill16(uint8_t *soa, const std::vector<std::vector<uint8_t>> &seq,
                   const int len[16], const int padStart[16], int nrows,
                   uint8_t padA, uint8_t padB, bool remap4to8) {
    for (int j = 0; j < 16; j++)
        for (int k = 0; k < nrows; k++) {
            uint8_t v;
            if (k < len[j]) {
                v = seq[j][k];
                if (remap4to8 && v == AMBIG_) v = AMBQ;
            } else {
                v = (k < padStart[j]) ? padA : padB;
            }
            soa[(size_t)k * 16 + j] = v;
        }
}

struct PackCase {
    uint8_t padA, padB;
    bool remap;
};

// The wrapper calls (bandedSWA.cpp / kswv.cpp) and adversarial pad values:
// a pad byte equal to the ambiguity code (pad bytes must not be remapped) and
// a pad byte equal to the remapped code.
const PackCase kCases[] = {
    {0xFF, 0xFF, false},   // rescue reference
    {16, 0xFF, true},      // rescue query (NEON_QPAD8 to the quantum, then 0xFF)
    {99, 99, false},       // extension reference
    {100, 100, true},      // extension query
    {AMBIG_, AMBIG_, true},
    {AMBQ, 7, true},
};

// Packs one group, with the sequences guarded after and then before, and
// compares the whole buffer (two canary rows included) with the scalar fill.
// Returns false (after a CHECK failure) on a mismatch.
bool pack_matches(const PackCase &c, const int len[16], const int padStart[16], int nrows,
                  std::mt19937 &rng) {
    std::vector<std::vector<uint8_t>> seq(16);
    for (int j = 0; j < 16; j++) {
        seq[j].resize(len[j]);
        for (int k = 0; k < len[j]; k++) {
            uint8_t b = (uint8_t)(rng() % 5);
            if (rng() % 9 == 0) b = c.padA;   // pad values as real bases
            if (rng() % 13 == 0) b = (uint8_t)(rng() % 256);
            seq[j][k] = b;
        }
    }
    const size_t bytes = (size_t)(nrows + 2) * 16;
    std::vector<uint8_t> want(bytes, 0xA5);
    scalar_fill16(want.data(), seq, len, padStart, nrows, c.padA, c.padB, c.remap);
    for (bool guardAfter : {true, false}) {
        std::vector<std::unique_ptr<GuardedSeq>> guarded;
        const uint8_t *seqp[16];
        for (int j = 0; j < 16; j++) {
            guarded.emplace_back(new GuardedSeq(seq[j], guardAfter));
            seqp[j] = guarded.back()->data();
        }
        void *gotp = nullptr;
        REQUIRE(posix_memalign(&gotp, 64, bytes) == 0);
        uint8_t *got = (uint8_t *)gotp;
        memset(got, 0xA5, bytes);
        neon_soa_pack16(got, seqp, len, padStart, nrows, c.padA, c.padB, c.remap);
        size_t i = 0;
        while (i < bytes && got[i] == want[i]) i++;
        const bool same = (i == bytes);
        if (!same) {
            const int row = (int)(i / 16), lane = (int)(i % 16);
            CHECK_MESSAGE(same, "nrows " << nrows << " first mismatch row " << row << " lane " << lane
                                         << " (len " << len[lane] << ", padStart " << padStart[lane]
                                         << ", padA " << (int)c.padA << ", padB " << (int)c.padB
                                         << ", remap " << c.remap << "): got " << (int)got[i]
                                         << " want " << (int)want[i]);
        }
        free(got);
        if (!same) return false;
    }
    return true;
}

}  // namespace

TEST_SUITE("unit/neon_soa_pack") {

TEST_CASE("neon_soa_pack16 matches the scalar SoA fill for every lane length and pad split") {
    std::mt19937 rng(0x5eed0016u);
    // Every length 0..56 (shorter than one tile, on and around the tile edges, a few tiles in),
    // with padStart at the length or 1, 7, 15, 16, 17 and 40 past it.
    const int padExtra[] = {0, 1, 7, 15, 16, 17, 40};
    for (const PackCase &c : kCases) {
        for (int L = 0; L <= 56; L++) {
            for (int pe : padExtra) {
                int len[16], padStart[16];
                for (int j = 0; j < 16; j++) {
                    // lane 0 has length L; the others spread around it
                    len[j] = j == 0 ? L : (int)(rng() % 60);
                    padStart[j] = len[j] + (j == 0 ? pe : padExtra[rng() % 7]);
                }
                int maxlen = 0;
                for (int j = 0; j < 16; j++) maxlen = len[j] > maxlen ? len[j] : maxlen;
                // nrows as the wrappers pass it (max length + 1), and a few rows past the pads
                for (int nrows : {maxlen + 1, maxlen + 1 + (int)(rng() % 40)})
                    if (!pack_matches(c, len, padStart, nrows, rng)) return;
            }
        }
    }
}

TEST_CASE("neon_soa_pack16 matches the scalar SoA fill on random lane groups") {
    std::mt19937 rng(0x5eed1016u);
    for (const PackCase &c : kCases) {
        for (int t = 0; t < 400; t++) {
            const int nrows = 1 + (int)(rng() % (t % 2 ? 20 : 300));
            int len[16], padStart[16];
            for (int j = 0; j < 16; j++) {
                len[j] = (int)(rng() % nrows);
                padStart[j] = len[j] + (int)(rng() % 40);   // may pass nrows
            }
            if (!pack_matches(c, len, padStart, nrows, rng)) return;
        }
    }
}

}  // TEST_SUITE

#else

TEST_CASE("neon_soa_pack16: skipped, not a NEON build") {
    MESSAGE("neon_soa_pack.h is NEON-only");
    CHECK(true);
}

#endif
