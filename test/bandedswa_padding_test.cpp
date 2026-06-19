// Regression test for the BandedPairWiseSW::getScores8/getScores16
// padding-lane / prefetch contract (issue 146 item: "harden the getScores8/16
// padding/prefetch contract").
//
// The batched banded-SW kernels round numPairs up to a whole number of SIMD
// lanes (SIMD_WIDTH8 / SIMD_WIDTH16) and (a) WRITE the trailing padding lanes
// pairArray[numPairs .. roundNumPairs) and (b) PREFETCH-READ pairArray[i+j+PFD]
// in the SoA-packing loop. A caller that sizes pairArray to exactly
// roundup(numPairs, SIMD_WIDTH) is honoring the documented contract — yet
// before the prefetch was bounded, the last group's pairArray[i+j+PFD] read ran
// PFD SeqPairs past the end of that allocation. Harmless in production (the
// pipeline over-allocates) but a real heap-buffer-overflow that AddressSanitizer
// aborts on, and exactly the crash a direct caller hit on AVX2/AVX-512.
//
// This test sizes pairArray to the MINIMAL contract size (no slack) and calls
// both entry points. With the `(i + j + PFD) < roundNumPairs` guard it runs
// clean; without it, the ASan CI lane (`make ASAN=1 bandedswa_padding_test`)
// reports a heap-buffer-overflow READ inside smithWatermanBatchWrapper{8,16}.
//
// Arch coverage is implicit in the build: the binary is compiled -march=native
// and linked against a native-tier src/bandedSWA.native.o, so the host's
// widest wrappers (NEON 256/16, or x86 512/8 + 512/16 with their PFD8/PFD16
// prefetch sites) are the ones exercised.

#include <cstdio>
#include <cstdint>
#include <vector>

#include "bandedSWA.h"

static int32_t roundup(int32_t n, int32_t w) { return ((n + w - 1) / w) * w; }

// Build a batch of `numPairs` valid pairs, allocate EXACTLY roundup(numPairs,
// width) SeqPair slots (the minimal contract — the kernel fills the padding
// lanes itself), and score it. Every pair points at the same modest ref/query
// prefix so the DP runs and the SoA-packing loop (which holds the prefetch)
// executes for every lane.
static void run(BandedPairWiseSW &bsw, int32_t numPairs, int32_t width,
                bool eightBit) {
    const int32_t segLen = 48;     // modest, < int8 range, exercises real DP
    const int32_t w       = 100;   // band

    // Size the ref/query buffers past segLen so the kernel's lookahead
    // _mm_prefetch(seqBufRef + idr + 64) hint lands inside the allocation (a
    // prefetch never faults, but keeping it in-bounds avoids any false read
    // signal and mirrors production, where seqBufRef is large).
    std::vector<uint8_t> ref((size_t)segLen + 128, 0);
    std::vector<uint8_t> qer((size_t)segLen + 128, 0);

    const int32_t rounded = roundup(numPairs, width);
    std::vector<SeqPair> pairs(rounded);   // EXACTLY the contract size, no slack
    for (int32_t i = 0; i < numPairs; i++) {
        SeqPair sp = {};
        sp.id = i; sp.seqid = i; sp.regid = i;
        sp.idr = 0; sp.idq = 0;
        sp.len1 = segLen; sp.len2 = segLen;
        sp.h0 = 1;
        sp.score = sp.tle = sp.gtle = sp.qle = sp.gscore = sp.max_off = -1;
        pairs[i] = sp;
    }

    if (eightBit) bsw.getScores8(pairs.data(), ref.data(), qer.data(), numPairs, 1, w);
    else          bsw.getScores16(pairs.data(), ref.data(), qer.data(), numPairs, 1, w);
}

int main() {
    const int8_t a = 1, b = 4, ambig = -1;
    const int o = 6, e = 1, zdrop = 100, end_bonus = 5;
    int8_t mat[25];
    { int k = 0;
      for (int i = 0; i < 4; ++i) { for (int j = 0; j < 4; ++j) mat[k++] = (i == j) ? a : -b; mat[k++] = ambig; }
      for (int j = 0; j < 5; ++j) mat[k++] = ambig; }

    BandedPairWiseSW bsw(o, e, o, e, zdrop, end_bonus, mat, a, b, 1);

    // numPairs that are not multiples of the width (so padding lanes are
    // written) plus exact multiples (so the prefetch reaches the very last
    // group). 1 forces a single partly-padded group; the others span several.
    const int32_t cases[] = {1, 3, 7, 33, 64, 100, 129};
    for (int32_t n : cases) {
        fprintf(stderr, "[bandedswa-padding] 8-bit  numPairs=%d ...\n", n);
        run(bsw, n, SIMD_WIDTH8, /*eightBit=*/true);
        fprintf(stderr, "[bandedswa-padding] 16-bit numPairs=%d ...\n", n);
        run(bsw, n, SIMD_WIDTH16, /*eightBit=*/false);
    }

    fprintf(stderr, "bandedswa_padding_test: OK\n");
    return 0;
}
