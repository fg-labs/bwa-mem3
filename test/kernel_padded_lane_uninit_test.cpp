// Uninitialized-padded-lane gate for the batched banded-SW kernels
// (BandedPairWiseSW::getScores8/16).
//
// (A kswv mate-rescue counterpart is a follow-up: kswv's SoA packing needs a
// more careful DP setup than the banded kernels, and a naive raw-batch caller
// trips an unrelated SoA-buffer write before the padded-lane read is reached.
// The kswv padded-lane fix is guarded meanwhile by the whole-aligner
// SAM-identity comparison; see the sibling PR.)
//
// The kernels round numPairs up to a whole SIMD width and process the trailing
// padding lanes pairArray[numPairs .. roundup(numPairs, SIMD_WIDTH)) alongside
// the real ones. Their per-lane loops read the padded lanes' input fields back:
// idr/idq form seq1/seq2 pointers, and h0 seeds the banded DP and the kswv
// minsc/endsc masks (which feed the cross-lane all-lanes-done exit reduction).
// A correct kernel must therefore INITIALIZE every padded-lane field it reads;
// if it leaves one uninitialized, that is a use of an indeterminate value (UB),
// and because h0 feeds a cross-lane reduction it can in principle perturb a real
// lane's result.
//
// The in-tree correctness tests build the batch with `std::vector<SeqPair>`,
// which VALUE-INITIALIZES the padding lanes to 0 -- so a kernel that forgets to
// initialize, say, h0 reads the same 0 a correct kernel writes, and the bug is
// invisible. Production, by contrast, reuses a non-zeroed pair buffer, so the
// padding lanes hold nonzero stale data.
//
// This test allocates the pair array with a RAW malloc (never memset / value-
// initialized), fills ONLY the real pairs [0, numPairs), and leaves the padding
// lanes genuinely uninitialized. Run it under MemorySanitizer or Valgrind
// memcheck: if any tier reads a padded-lane field the kernel did not first
// write, the tool reports a use of an uninitialized value. It PASSES only when
// every wrapper initializes every padded-lane field it reads back.
//
// Arch coverage is implicit in the build: compiled -march=native against the
// native-tier kernel objects, so the host's widest wrappers are exercised.

#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "bandedSWA.h"

static int32_t roundup(int32_t n, int32_t w) { return ((n + w - 1) / w) * w; }

// Allocate roundup(numPairs, width) SeqPair with raw malloc (padding lanes stay
// uninitialized), fill only the real pairs, and return the buffer. The caller
// frees it after the kernel runs.
static SeqPair *make_uninit_batch(int32_t numPairs, int32_t width, int32_t segLen) {
    const int32_t rounded = roundup(numPairs, width);
    SeqPair *pairs = (SeqPair *)malloc((size_t)rounded * sizeof(SeqPair));
    for (int32_t i = 0; i < numPairs; i++) {
        // Fully initialize every field of every REAL pair so the only
        // indeterminate memory in the batch is the padding region.
        SeqPair sp;
        memset(&sp, 0, sizeof sp);
        sp.id = i; sp.seqid = i; sp.regid = i;
        sp.idr = 0; sp.idq = 0;
        sp.len1 = segLen; sp.len2 = segLen;
        sp.h0 = 1;                            // nonzero seed so the banded DP runs
        sp.score = sp.tle = sp.gtle = sp.qle = sp.gscore = sp.max_off = -1;
        pairs[i] = sp;
    }
    return pairs;   // [numPairs, rounded) deliberately left uninitialized
}

static void run_bsw(BandedPairWiseSW &bsw, int32_t numPairs, int32_t width,
                    bool eightBit, uint8_t *ref, uint8_t *qer, int32_t segLen) {
    SeqPair *pairs = make_uninit_batch(numPairs, width, segLen);
    if (eightBit) bsw.getScores8 (pairs, ref, qer, numPairs, 1, /*w=*/100);
    else          bsw.getScores16(pairs, ref, qer, numPairs, 1, /*w=*/100);
    free(pairs);
}

int main() {
    const int8_t a = 1, b = 4, ambig = -1;
    const int o = 6, e = 1, zdrop = 100, end_bonus = 5;
    const int32_t segLen = 48;   // < int8 range; a real DP runs for every lane

    int8_t mat[25];
    { int k = 0;
      for (int i = 0; i < 4; ++i) { for (int j = 0; j < 4; ++j) mat[k++] = (i == j) ? a : -b; mat[k++] = ambig; }
      for (int j = 0; j < 5; ++j) mat[k++] = ambig; }

    // Ref/query large enough that any in-bounds prefetch/read of a REAL lane
    // stays inside the allocation (a padded lane's idr/idq is what we're testing).
    std::vector<uint8_t> ref((size_t)segLen + 256, 0);
    std::vector<uint8_t> qer((size_t)segLen + 256, 0);

    BandedPairWiseSW bsw(o, e, o, e, zdrop, end_bonus, mat, a, b, 1);

    // Counts that are NOT multiples of the tier widths (16/8, 32/16, 64/32), so
    // every tier gets a partial final group -- the only case with padding lanes.
    const int32_t cases[] = {1, 3, 7, 23, 45, 100, 129};
    for (int32_t n : cases) {
        fprintf(stderr, "[padded-lane] bsw  8-bit numPairs=%d\n", n);
        run_bsw(bsw, n, SIMD_WIDTH8,  /*eightBit=*/true,  ref.data(), qer.data(), segLen);
        fprintf(stderr, "[padded-lane] bsw 16-bit numPairs=%d\n", n);
        run_bsw(bsw, n, SIMD_WIDTH16, /*eightBit=*/false, ref.data(), qer.data(), segLen);

    }

    fprintf(stderr, "kernel_padded_lane_uninit_test: OK (no uninitialized padded-lane read)\n");
    return 0;
}
