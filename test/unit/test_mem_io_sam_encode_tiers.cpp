// test/unit/test_mem_io_sam_encode_tiers.cpp
//
// The SEQ/QUAL encoders mem_aln2sam and the BAM writer call
// (sam_encode_seq_fwd, sam_encode_seq_rev, sam_encode_qual_rev) on every SIMD
// tier `mem` runs: x86's avx2 and avx512bw copies by their suffixed names
// (src/kernel_dispatch.h; libbwa.a carries every tier), the NEON build's
// unsuffixed ones on arm64. test_sam_encode.cpp checks the copy this binary's
// own baseline dispatches to; the per-tier copies are what `mem` actually
// runs, and a tier-specific rewrite (a wider vector loop, a different tail)
// could break one tier alone.
//
// Oracle: the encoders' documented model, byte by byte. A base code c (0-4 for
// A C G T N) prints as "ACGTN"[min(c, 4)]; the reverse encoder prints the
// complement "TGCAN"[min(c, 4)] of the input read backwards; the quality
// encoder reverses its input. Every length from 0 to 300 (each vector width's
// full blocks, overlapped tails and the short scalar path) plus a few long
// ones, all 256 byte values for the clamp, and guard bytes on both sides of
// the destination, so a write outside [dst, dst + n) is caught.

#include "doctest/doctest.h"
#include "sam_encode.h"

#include <cstdint>
#include <cstring>
#include <random>
#include <string>
#include <vector>

namespace {

typedef void (*SeqFn)(char *, const uint8_t *, int);
typedef void (*QualFn)(char *, const char *, int);

struct Tier {
    const char *name;
    SeqFn fwd, rev;
    QualFn qual_rev;
};

}  // namespace

#if defined(__x86_64__) || defined(__i386__)
extern "C" void sam_encode_seq_fwd_avx2(char *dst, const uint8_t *src, int n);
extern "C" void sam_encode_seq_rev_avx2(char *dst, const uint8_t *src, int n);
extern "C" void sam_encode_qual_rev_avx2(char *dst, const char *src, int n);
extern "C" void sam_encode_seq_fwd_avx512bw(char *dst, const uint8_t *src, int n);
extern "C" void sam_encode_seq_rev_avx512bw(char *dst, const uint8_t *src, int n);
extern "C" void sam_encode_qual_rev_avx512bw(char *dst, const char *src, int n);
#endif

namespace {

std::vector<Tier> runnable_tiers() {
    std::vector<Tier> ts;
#if defined(__x86_64__) || defined(__i386__)
    __builtin_cpu_init();
    if (__builtin_cpu_supports("avx2")) {
        Tier t = {"avx2", sam_encode_seq_fwd_avx2, sam_encode_seq_rev_avx2, sam_encode_qual_rev_avx2};
        ts.push_back(t);
    } else {
        MESSAGE("no AVX2 on this host; the avx2 tier is skipped");
    }
    if (__builtin_cpu_supports("avx512bw")) {
        Tier t = {"avx512bw", sam_encode_seq_fwd_avx512bw, sam_encode_seq_rev_avx512bw,
                  sam_encode_qual_rev_avx512bw};
        ts.push_back(t);
    } else {
        MESSAGE("no AVX-512BW on this host; the avx512bw tier is skipped");
    }
#else
    Tier t = {"native", sam_encode_seq_fwd, sam_encode_seq_rev, sam_encode_qual_rev};
    ts.push_back(t);
#endif
    return ts;
}

char fwd_char(uint8_t c) { return "ACGTN"[c <= 4 ? c : 4]; }
char rev_char(uint8_t c) { return "TGCAN"[c <= 4 ? c : 4]; }

const int kGuard = 64;
const char kCanary = (char)0xA5;

// Runs one encoder into a guarded buffer: the n output bytes, or a marker
// string when a guard byte changed.
template <class Fn, class T>
std::string run_guarded(Fn fn, const std::vector<T> &src, int n) {
    std::vector<char> buf((size_t)n + 2 * kGuard, kCanary);
    const std::vector<T> before = src;
    fn(buf.data() + kGuard, src.data(), n);
    for (int i = 0; i < kGuard; ++i)
        if (buf[(size_t)i] != kCanary || buf[(size_t)(kGuard + n + i)] != kCanary)
            return "<write outside [dst, dst + n)>";
    if (src != before) return "<source modified>";
    return std::string(buf.data() + kGuard, (size_t)n);
}

std::vector<int> lengths() {
    std::vector<int> ls;
    for (int n = 0; n <= 300; ++n) ls.push_back(n);
    const int longer[] = {511, 512, 513, 1023, 1024, 1025, 4097, 10007};
    for (int n : longer) ls.push_back(n);
    return ls;
}

}  // namespace

TEST_CASE("mem_io sam_encode: every runnable tier encodes SEQ and QUAL byte for byte"
          * doctest::test_suite("unit/sam_encode")) {
    const std::vector<Tier> tiers = runnable_tiers();
    REQUIRE(!tiers.empty());
    std::mt19937 rng(424242u);

    SUBCASE("A C G T N codes") {
        for (const Tier &t : tiers) {
            CAPTURE(t.name);
            for (int n : lengths()) {
                CAPTURE(n);
                std::vector<uint8_t> codes((size_t)n);
                std::vector<char> qual((size_t)n);
                for (int i = 0; i < n; ++i) {
                    codes[(size_t)i] = (uint8_t)(rng() % 5);
                    qual[(size_t)i] = (char)('!' + rng() % 94);
                }
                std::string fwd, rev, qrev;
                for (int i = 0; i < n; ++i) {
                    fwd += fwd_char(codes[(size_t)i]);
                    rev += rev_char(codes[(size_t)(n - 1 - i)]);
                    qrev += qual[(size_t)(n - 1 - i)];
                }
                CHECK(run_guarded(t.fwd, codes, n) == fwd);
                CHECK(run_guarded(t.rev, codes, n) == rev);
                CHECK(run_guarded(t.qual_rev, qual, n) == qrev);
            }
        }
    }
    SUBCASE("every byte value: codes above 4 print N, quality bytes pass through") {
        for (const Tier &t : tiers) {
            CAPTURE(t.name);
            // 0..255 in order, then shuffled by a fixed permutation of the
            // positions, at lengths that put each value in vector bodies and tails.
            const int ns[] = {256, 255, 17, 63, 65, 129};
            for (int n : ns) {
                CAPTURE(n);
                std::vector<uint8_t> codes((size_t)n);
                std::vector<char> qual((size_t)n);
                for (int i = 0; i < n; ++i) {
                    codes[(size_t)i] = (uint8_t)((i * 151 + 7) % 256);
                    qual[(size_t)i] = (char)codes[(size_t)i];
                }
                std::string fwd, rev, qrev;
                for (int i = 0; i < n; ++i) {
                    fwd += fwd_char(codes[(size_t)i]);
                    rev += rev_char(codes[(size_t)(n - 1 - i)]);
                    qrev += qual[(size_t)(n - 1 - i)];
                }
                CHECK(run_guarded(t.fwd, codes, n) == fwd);
                CHECK(run_guarded(t.rev, codes, n) == rev);
                CHECK(run_guarded(t.qual_rev, qual, n) == qrev);
            }
        }
    }
    SUBCASE("an unaligned source and destination") {
        for (const Tier &t : tiers) {
            CAPTURE(t.name);
            for (int off = 1; off < 64; off += 7) {
                CAPTURE(off);
                const int n = 150;
                std::vector<uint8_t> raw((size_t)(n + off));
                std::vector<char> rawq((size_t)(n + off));
                for (size_t i = 0; i < raw.size(); ++i) {
                    raw[i] = (uint8_t)(rng() % 5);
                    rawq[i] = (char)('!' + rng() % 94);
                }
                std::vector<char> out((size_t)(n + off + 1), kCanary);
                t.fwd(out.data() + off, raw.data() + off, n);
                std::string want;
                for (int i = 0; i < n; ++i) want += fwd_char(raw[(size_t)(off + i)]);
                CHECK(std::string(out.data() + off, (size_t)n) == want);
                t.rev(out.data() + off, raw.data() + off, n);
                want.clear();
                for (int i = 0; i < n; ++i) want += rev_char(raw[(size_t)(off + n - 1 - i)]);
                CHECK(std::string(out.data() + off, (size_t)n) == want);
                t.qual_rev(out.data() + off, rawq.data() + off, n);
                want.clear();
                for (int i = 0; i < n; ++i) want += rawq[(size_t)(off + n - 1 - i)];
                CHECK(std::string(out.data() + off, (size_t)n) == want);
                CHECK(out[(size_t)(off + n)] == kCanary);
            }
        }
    }
}
