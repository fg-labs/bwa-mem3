// Unit tests for src/seq_reverse.h. The mate-rescue revseq and bwa_gen_cigar4's
// reverse-strand flips must leave exactly what the scalar swap loop they
// replaced leaves, so that loop is the oracle here. Both entry points are
// checked: seq_reverse_inplace (the 32-byte AVX2 blocks and the 16-byte
// NEON/SSSE3 blocks where the target has them, as every shipped target has the
// latter) and seq_reverse_inplace_u64 (the 8-byte bswap + single-byte fallback,
// which seq_reverse_inplace alone would never run on a host with a vector
// path). Lengths 0..300 cover the block loops, the overlapping block pairs and
// the single-byte tail, at every alignment within a 32-byte block, with guard
// bytes around.

#include "doctest/doctest.h"

#include "seq_reverse.h"

#include <climits>
#include <cstdint>
#include <random>
#include <string>
#include <vector>

namespace {

void oracle_revseq(int l, uint8_t *s)
{
    int i, t;
    for (i = 0; i < l >> 1; ++i)
        t = s[i], s[i] = s[l - 1 - i], s[l - 1 - i] = t;
}

/* Returns the number of (length, alignment, content) cases compared. */
template <typename F>
long check_against_oracle(F reverse)
{
    std::mt19937 rng(0x5eed);
    const int guard = 32;
    long cases = 0, mismatches = 0;
    for (int n = 0; n <= 300; ++n) {
        for (int off = 0; off < 32; ++off) {
            for (int rep = 0; rep < 2; ++rep) {
                std::vector<uint8_t> buf(guard + off + n + guard);
                for (auto &b : buf) b = (uint8_t)rng();
                if (rep == 1) /* base codes, as the rescue buffers hold */
                    for (auto &b : buf) b &= 3;
                std::vector<uint8_t> want = buf;
                oracle_revseq(n, want.data() + guard + off);
                reverse(buf.data() + guard + off, n);
                ++cases;
                if (buf != want && ++mismatches <= 5) {
                    INFO("n=" << n << " off=" << off << " rep=" << rep);
                    CHECK(buf == want);
                }
            }
        }
    }
    CHECK(mismatches == 0);
    return cases;
}

const long kCases = 301L * 32 * 2;

} // namespace

TEST_CASE("seq_reverse_inplace matches the scalar swap loop for lengths 0..300"
          * doctest::test_suite("unit/seq_reverse"))
{
#if defined(__AVX2__)
    MESSAGE(std::string("32-byte vector path: AVX2"));
#endif
#if defined(BWA3_SEQ_REV_VEC16)
    MESSAGE(std::string(BWA3_SEQ_REV_VEC16 == 1 ? "16-byte vector path: NEON" : "16-byte vector path: SSSE3"));
#else
    MESSAGE(std::string("no 16-byte vector path on this target: fallback only"));
#endif
    CHECK(check_against_oracle([](uint8_t *s, int n) { seq_reverse_inplace(s, n); }) == kCases);
}

TEST_CASE("seq_reverse_inplace_u64 (fallback) matches the scalar swap loop for lengths 0..300"
          * doctest::test_suite("unit/seq_reverse"))
{
    CHECK(check_against_oracle([](uint8_t *s, int n) { seq_reverse_inplace_u64(s, n); }) == kCases);
}

TEST_CASE("seq_reverse treats a length below 2 as a no-op, including INT_MIN"
          * doctest::test_suite("unit/seq_reverse"))
{
    const std::vector<uint8_t> orig = {1, 2, 3, 4, 5, 6, 7, 8};
    for (int n : {1, 0, -1, -17, INT_MIN}) {
        INFO("n=" << n);
        std::vector<uint8_t> a = orig;
        seq_reverse_inplace(a.data(), n);
        CHECK(a == orig);
        std::vector<uint8_t> b = orig;
        seq_reverse_inplace_u64(b.data(), n);
        CHECK(b == orig);
    }
}
