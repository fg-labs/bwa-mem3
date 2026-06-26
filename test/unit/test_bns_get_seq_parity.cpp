// Parity test: bns_get_seq_into (caller-buffer write) vs bns_get_seq_v2
// (zero-copy pointer into pre-unpacked ref_string).
//
// Why this matters: bns_get_seq_into and bns_get_seq_v2 are the two ways
// in-tree callers fetch ref bases — the former does inline 2-bit decode
// into a caller buffer; the latter returns a pointer into a pre-decoded
// .0123 ref_string. Both must return byte-identical bases for the same
// input range, since callers pick one path based on whether ref_string
// is materialized in their context.
//
// Test strategy: build a small synthetic packed reference (.pac-style 2-bit
// packed), independently build the equivalent ref_string (forward + reverse
// complement, the format .0123 stores), and walk both functions over many
// (beg, end) ranges including:
//   - forward strand only (beg < l_pac, end <= l_pac)
//   - reverse strand only (beg >= l_pac)
//   - bridging forward/reverse boundary (returns *len = 0)
//   - boundary edges (beg = 0, end = 2*l_pac)
//   - swapped beg/end (legacy code does an XOR-swap on entry)

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <random>
#include <vector>

#include "doctest/doctest.h"

#include "bntseq.h"

// _set_pac is published by bntseq.h.
//
// Pac layout: byte (k>>2) holds 4 bases at positions 4k..4k+3 (top-first).
//   pac[i] bit 7:6 = base at position 4i+0
//   pac[i] bit 5:4 = base at position 4i+1
//   pac[i] bit 3:2 = base at position 4i+2
//   pac[i] bit 1:0 = base at position 4i+3

namespace {

// Build a packed-2bit reference of `l_pac` bases driven by `rng`, plus the
// ref_string format (forward bases 0..l_pac-1 followed by reverse-complemented
// bases l_pac..2*l_pac-1). Both representations describe the same forward
// sequence; ref_string also materializes the reverse strand.
void build_synthetic_ref(std::mt19937 &rng, int64_t l_pac, std::vector<uint8_t> &pac,
                         std::vector<uint8_t> &ref_string) {
    pac.assign((l_pac + 3) / 4, 0);
    ref_string.assign(static_cast<size_t>(2 * l_pac), 0);

    std::uniform_int_distribution<int> base(0, 3);
    for (int64_t i = 0; i < l_pac; ++i) {
        uint8_t c = static_cast<uint8_t>(base(rng));
        _set_pac(pac.data(), i, c);
        ref_string[i] = c;
    }
    // Reverse-complement: ref_string[2*l_pac - 1 - i] = 3 ^ ref_string[i].
    for (int64_t i = 0; i < l_pac; ++i) {
        ref_string[2 * l_pac - 1 - i] = static_cast<uint8_t>(3 ^ ref_string[i]);
    }
}

// Run both functions for one (beg, end) and assert byte-identical output.
// Skip cases where bns_get_seq_into reports *len_out = 0 (boundary-bridging).
void check_range(int64_t l_pac, const uint8_t *pac, const uint8_t *ref_string,
                 int64_t beg, int64_t end) {
    // bns_get_seq_into needs caller buffer >= |end - beg|. Allocate worst-case
    // (2 * l_pac) and zero so any unwritten bytes show up as a CHECK miss.
    std::vector<uint8_t> into_buf(2 * l_pac, 0xff);
    int64_t len_into = 0;
    bns_get_seq_into(l_pac, pac, beg, end, into_buf.data(), &len_into);

    // bns_get_seq_v2 expects a scratch buffer (`seqb`); the function ignores
    // it (returns a ref_string pointer directly), but we pass a real one to
    // match the production call shape.
    std::vector<uint8_t> scratch(2 * l_pac, 0xff);
    int64_t len_v2 = 0;
    uint8_t *v2 = bns_get_seq_v2(l_pac, pac, beg, end, &len_v2,
                                 const_cast<uint8_t *>(ref_string), scratch.data());

    CHECK(len_into == len_v2);
    if (len_into == 0) {
        // Boundary-bridging: both should report empty.
        return;
    }
    REQUIRE(v2 != nullptr);
    for (int64_t i = 0; i < len_into; ++i) {
        CAPTURE(beg);
        CAPTURE(end);
        CAPTURE(i);
        CHECK(static_cast<int>(into_buf[i]) == static_cast<int>(v2[i]));
    }
}

}  // namespace

TEST_CASE("bns_get_seq_into vs v2: forward strand exhaustive small ref") {
    std::mt19937 rng(0x5EEDu);
    constexpr int64_t L = 64;
    std::vector<uint8_t> pac, ref_string;
    build_synthetic_ref(rng, L, pac, ref_string);

    // All (beg, end) pairs with 0 <= beg <= end <= L (forward-only).
    for (int64_t beg = 0; beg <= L; ++beg) {
        for (int64_t end = beg; end <= L; ++end) {
            check_range(L, pac.data(), ref_string.data(), beg, end);
        }
    }
}

TEST_CASE("bns_get_seq_into vs v2: reverse strand exhaustive small ref") {
    std::mt19937 rng(0xC0FFEEu);
    constexpr int64_t L = 64;
    std::vector<uint8_t> pac, ref_string;
    build_synthetic_ref(rng, L, pac, ref_string);

    // All (beg, end) pairs entirely in the reverse strand: L <= beg <= end <= 2L.
    for (int64_t beg = L; beg <= 2 * L; ++beg) {
        for (int64_t end = beg; end <= 2 * L; ++end) {
            check_range(L, pac.data(), ref_string.data(), beg, end);
        }
    }
}

TEST_CASE("bns_get_seq_into vs v2: bridging forward/reverse boundary returns empty") {
    std::mt19937 rng(0xBEEFu);
    constexpr int64_t L = 64;
    std::vector<uint8_t> pac, ref_string;
    build_synthetic_ref(rng, L, pac, ref_string);

    // beg < l_pac, end > l_pac → should return *len = 0 in both.
    for (int64_t beg : {int64_t{0}, int64_t{1}, int64_t{L / 2}, int64_t{L - 1}}) {
        for (int64_t end : {int64_t{L + 1}, int64_t{L + L / 2}, int64_t{2 * L}}) {
            int64_t len_into = 0, len_v2 = 0;
            std::vector<uint8_t> into_buf(2 * L, 0xff);
            bns_get_seq_into(L, pac.data(), beg, end, into_buf.data(), &len_into);
            std::vector<uint8_t> scratch(2 * L, 0xff);
            (void)bns_get_seq_v2(L, pac.data(), beg, end, &len_v2,
                                 ref_string.data(), scratch.data());
            CHECK(len_into == 0);
            CHECK(len_v2 == 0);
        }
    }
}

TEST_CASE("bns_get_seq_into vs v2: swapped beg/end (XOR-swap path)") {
    std::mt19937 rng(0xDEADu);
    constexpr int64_t L = 64;
    std::vector<uint8_t> pac, ref_string;
    build_synthetic_ref(rng, L, pac, ref_string);

    // Swap path: bns_get_seq_into does `if (end < beg) end ^= beg, beg ^= end, end ^= beg`.
    // v2 does the same. Both should produce the same canonicalized [beg, end).
    check_range(L, pac.data(), ref_string.data(), /*beg=*/40, /*end=*/10);
    check_range(L, pac.data(), ref_string.data(), /*beg=*/L + 30, /*end=*/L + 5);
    check_range(L, pac.data(), ref_string.data(), /*beg=*/L, /*end=*/0);  // bridging
}

TEST_CASE("bns_get_seq_v2: NULL ref_string pac-fetches bases identical to .0123") {
    // pac-fetch contract: a NULL ref_string is no longer a "return NULL" guard —
    // it means "reconstruct this window from `pac` on demand" (the .0123 view was
    // not loaded). The bases returned MUST be byte-identical to what a non-NULL
    // ref_string (the unpacked .0123) view would yield for the same [beg, end),
    // including the reverse-complement half [l_pac, 2*l_pac). This is the
    // identity gate for dropping .0123 in favor of unpacking from .pac.
    std::mt19937 rng(0x5EEDu);
    constexpr int64_t L = 64;
    std::vector<uint8_t> pac, ref_string;
    build_synthetic_ref(rng, L, pac, ref_string);

    // Forward windows (beg,end < L) and reverse-complement windows (beg >= L);
    // bridging windows return len 0 and are covered by the boundary test above.
    for (int64_t beg : {int64_t{0}, int64_t{10}, int64_t{L - 20}, int64_t{L}, int64_t{L + 5}}) {
        for (int64_t span : {int64_t{1}, int64_t{17}, int64_t{40}}) {
            int64_t end = beg + span;
            if (end > 2 * L || (beg < L && end > L)) continue;  // skip bridging
            int64_t len_ref = 0, len_pf = 0;
            std::vector<uint8_t> scratch_ref(2 * L, 0xff);
            uint8_t *seq_ref = bns_get_seq_v2(L, pac.data(), beg, end, &len_ref,
                                              ref_string.data(), scratch_ref.data());
            // ref_string == NULL → pac-fetch path. Result points into a thread-local
            // buffer valid until the next pac-fetch call; compare immediately.
            std::vector<uint8_t> scratch_pf(2 * L, 0xff);
            uint8_t *seq_pf = bns_get_seq_v2(L, pac.data(), beg, end, &len_pf,
                                             /*ref_string=*/nullptr, scratch_pf.data());
            REQUIRE(len_pf == len_ref);
            REQUIRE(len_pf == end - beg);
            REQUIRE(seq_ref != nullptr);
            REQUIRE(seq_pf != nullptr);
            CHECK(std::memcmp(seq_pf, seq_ref, (size_t)len_pf) == 0);
        }
    }
}

TEST_CASE("bns_get_seq_v2: NULL ref_string bridging window returns empty without bases") {
    // pac-fetch bridging contract: a window crossing the forward/reverse boundary
    // (beg < l_pac < end) yields nothing, exactly like the .0123 path. The
    // implementation short-circuits this BEFORE growing its scratch buffer so a
    // bad bridge query can't balloon RSS toward the doubled reference. Observable
    // surface here is the empty result; the no-alloc guarantee is internal.
    std::mt19937 rng(0xB1D6Eu);
    constexpr int64_t L = 64;
    std::vector<uint8_t> pac, ref_string;
    build_synthetic_ref(rng, L, pac, ref_string);

    for (int64_t beg : {int64_t{0}, int64_t{1}, int64_t{L / 2}, int64_t{L - 1}}) {
        for (int64_t end : {int64_t{L + 1}, int64_t{L + L / 2}, int64_t{2 * L}}) {
            CAPTURE(beg);
            CAPTURE(end);
            int64_t len_pf = -1;
            std::vector<uint8_t> scratch(2 * L, 0xff);
            uint8_t *seq_pf = bns_get_seq_v2(L, pac.data(), beg, end, &len_pf,
                                             /*ref_string=*/nullptr, scratch.data());
            CHECK(len_pf == 0);
            CHECK(seq_pf == nullptr);
        }
    }
    // Swapped bridging (end < beg) canonicalizes to the same empty result.
    int64_t len_pf = -1;
    std::vector<uint8_t> scratch(2 * L, 0xff);
    uint8_t *seq_pf = bns_get_seq_v2(L, pac.data(), /*beg=*/L + 10, /*end=*/L - 10,
                                     &len_pf, /*ref_string=*/nullptr, scratch.data());
    CHECK(len_pf == 0);
    CHECK(seq_pf == nullptr);
}

TEST_CASE("bns_get_seq_into vs v2: large random ref, randomized ranges") {
    std::mt19937 rng(0x900DCAFEu);
    constexpr int64_t L = 16384;  // larger to exercise unaligned reads
    std::vector<uint8_t> pac, ref_string;
    build_synthetic_ref(rng, L, pac, ref_string);

    std::uniform_int_distribution<int64_t> pos(0, 2 * L);
    for (int trial = 0; trial < 1000; ++trial) {
        int64_t a = pos(rng);
        int64_t b = pos(rng);
        check_range(L, pac.data(), ref_string.data(), a, b);
    }
}
