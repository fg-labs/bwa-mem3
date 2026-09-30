/* ungapped_ext.h -- the ungapped diagonal fast path in front of the banded-SW
 * extension ladder (bwamem.cpp): the mismatch scan, the diagonal walk that
 * reproduces the extension kernel's diagonal, and ungapped_analyze, which
 * decides HIT (skip the DP), TIGHT (carry a band bound) or FALLBACK.
 *
 * Header-only so test/unit/test_ungapped_fastpath.cpp can run it
 * differentially against the scalar extension kernels. Only bwamem.cpp
 * includes it in the aligner. */
#ifndef BWA_MEM3_UNGAPPED_EXT_H
#define BWA_MEM3_UNGAPPED_EXT_H

#include <stdint.h>
#include "macro.h"
#include "utils.h"

/* SIMD compatibility layer for ARM/x86 (same selection as bandedSWA.h). */
#if defined(__ARM_NEON) || defined(__aarch64__) || defined(APPLE_SILICON)
    #include "simd_compat.h"
#elif (__AVX512BW__ || __AVX2__)
    #include <immintrin.h>
#else
    #include <smmintrin.h>  // for SSE4.1
#endif

// ungapped diagonal-extension analyzer.
//
// Walks the N-step diagonal once (SIMD scan + bitmap + scalar trajectory).
// Emits three possible status codes:
//
//   FP_STATUS_HIT      — ungapped is provably optimal (≤ x_threshold
//                        mismatches, no ambig). Caller skips banded SW
//                        and uses sp.score/qle/gscore/gtle directly.
//
//   FP_STATUS_TIGHT    — fast-path fails (too many mismatches) but the
//                        ungapped score bounds the useful SW band:
//                          tight_band = ceil((min(len1,len2)·a - ungapped_score)
//                                            / (o_min + e_min))
//                        The certified probe rung (bsw_tb_probe_rung)
//                        uses it to pick a narrow probe width; the
//                        --adaptive-band ladder uses it to stop early
//                        (the bound is an upper bound on any gapped
//                        score). The exact ladders do NOT stop on it
//                        -- see ACCEPT_PAIR.
//
//   FP_STATUS_FALLBACK — ambig base or out-of-range length. Caller uses
//                        opt->w with the full retry loop.
//
// Derivation of tight_band: for any alignment with band offset B from
// diagonal, min B gaps are required; cost ≥ B · (o_min + e_min). Max
// alignment score (all matches) ≤ min(len1, len2) · a. For any gapped
// alignment to beat the observed ungapped max score S:
//   min_len·a - B·(o_min+e_min) > S
//   B < (min_len·a - S) / (o_min + e_min)
// So any band ≥ tight_band is sufficient; narrower bands suffice too
// when the gapped alternative score is < min_len·a. Starting SW at
// tight_band is strictly correct and avoids over-banded DP work.
#define FP_N_MAX 512
/* Length cap for emitting a tight_band at all (longer pairs HIT or FALLBACK).
 * tight_band feeds the certified probe rung (bsw_tb_probe_rung), which finalizes on
 * band_cert_ok at its probe width -- not on tb -- inside the
 * mem_band_cert_params_safe_w envelope, and the --adaptive-band ladder's narrowing-tier
 * stop. The exact ladders never stop on it (see ACCEPT_PAIR). Keep this <= the old
 * scanner cap so the set of pairs that carry a tight_band is unchanged. */
#define FP_TIGHT_MAX 128
/* Words in the per-pair mismatch bitmap. Sized to FP_N_MAX so every scanned
 * position 0..FP_N_MAX-1 has a bit; a shift by the in-word offset (< 64) can
 * never alias. (Must track FP_N_MAX: at 128 this is 2 words = the old
 * mis_lo/mis_hi pair.) */
#define FP_MIS_NWORDS ((FP_N_MAX + 63) / 64)
#define FP_STATUS_FALLBACK  0
#define FP_STATUS_HIT       1
#define FP_STATUS_TIGHT     2

/* Q3 helper: would-be ungapped extension score for arbitrary N. Mirrors the
 * HIT-path scalar walk in ungapped_analyze (cur with floor at 0; max_sc
 * tracker; ambig terminates the walk to match analyze's FALLBACK semantics).
 * Returned score is what an ungapped extension would produce on this pair
 * regardless of whether the fast-path actually triggered HIT. */
static inline int ungapped_walk_score(const uint8_t *qs, const uint8_t *rs,
                                       int N, int h0, int a, int b)
{
    int cur = h0, max_sc = h0;
    for (int j = 0; j < N; j++) {
        uint8_t qj = qs[j], rj = rs[j];
        if (qj >= 4 || rj >= 4) break;       /* ambig: terminate the walk */
        if (cur == 0) continue;
        if (qj == rj) cur += a;
        else {
            cur -= b;
            if (cur < 0) cur = 0;
        }
        if (cur >= max_sc) max_sc = cur;
    }
    return max_sc;
}

/* Mismatch-run form of the ungapped walk (byte-identical). ungapped_analyze has
 * already rejected every pair holding an ambiguous base in [0, N) and built the
 * mismatch bitmask `mis`, so the per-base walk it then runs (the HIT walk
 * below, and ungapped_walk_score on the TIGHT path) only ever sees match /
 * mismatch steps. With a > 0, b >= 0 and h0 > 0 those steps have a closed form
 * between mismatches:
 *  - cur > 0 on entry to a run of r matches stays > 0 and rises by a per step,
 *    so the run's per-step `cur >= max_sc` updates reduce to one test of the
 *    run's last value (strictly increasing: the last step is both the largest
 *    and the rightmost, which is what the >= tie-break selects); if the last
 *    value is below max_sc, so is every earlier one.
 *  - a mismatch step is applied exactly as in the loop (cur -= b, floor at 0,
 *    then the >= test).
 *  - once cur reaches 0 the loop's `if (cur == 0) continue` freezes every
 *    output for the rest of the walk, so the run form stops there.
 * So the walk costs one step per mismatch (<= x_threshold on a HIT, a handful
 * on TIGHT) instead of one per base. Outside that envelope (a <= 0, b < 0,
 * h0 <= 0) the original per-base loop runs. Outputs: max_sc, max_i (the
 * rightmost position reaching it), and the final cur. */
/* The per-base walk over the mismatch bitmask: match +a, mismatch -b floored at 0,
 * and frozen once the score reaches 0; max_sc and max_i (the rightmost position
 * reaching it, the >= tie-break) as the banded SW sees them. The reference form
 * the run-length walk below must equal, and its fallback outside that envelope. */
static inline void ungapped_walk_mis_perbase(const uint64_t *mis, int N, int h0, int a, int b,
                                             int *out_max_sc, int *out_max_i, int *out_cur)
{
    int cur = h0, max_sc = h0, max_i = 0;
    for (int j = 0; j < N; j++) {
        int is_mis = (int)((mis[j >> 6] >> (j & 63)) & 1ULL);
        if (cur == 0) continue;
        if (!is_mis) cur += a;
        else {
            cur -= b;
            if (cur < 0) cur = 0;
        }
        if (cur >= max_sc) { max_sc = cur; max_i = j + 1; }
    }
    *out_max_sc = max_sc; *out_max_i = max_i; *out_cur = cur;
}

/* Mismatch-run form of the ungapped walk (byte-identical to
 * ungapped_walk_mis_perbase). ungapped_analyze has already rejected every pair
 * holding an ambiguous base in [0, N) and built the mismatch bitmask `mis`, with
 * no bit set at or past N, so the per-base walk it then runs (the HIT walk below,
 * and the TIGHT path's walk score) only ever sees match / mismatch steps. With
 * a > 0, b >= 0 and h0 > 0 those steps have a closed form between mismatches:
 *  - cur > 0 on entry to a run of r matches stays > 0 and rises by a per step,
 *    so the run's per-step `cur >= max_sc` updates reduce to one test of the
 *    run's last value (strictly increasing: the last step is both the largest
 *    and the rightmost, which is what the >= tie-break selects); if the last
 *    value is below max_sc, so is every earlier one.
 *  - a mismatch step is applied exactly as in the loop (cur -= b, floor at 0,
 *    then the >= test).
 *  - once cur reaches 0 the per-base walk's `if (cur == 0) continue` freezes
 *    every output for the rest of the walk, so the run form stops there.
 * So the walk costs one step per mismatch (<= x_threshold on a HIT, a handful
 * on TIGHT) instead of one per base. Outside that envelope (a <= 0, b < 0,
 * h0 <= 0) the per-base walk runs. */
static inline void ungapped_walk_mis(const uint64_t *mis, int N, int h0, int a, int b,
                                     int *out_max_sc, int *out_max_i, int *out_cur)
{
    if (!(a > 0 && b >= 0 && h0 > 0)) {
        ungapped_walk_mis_perbase(mis, N, h0, a, b, out_max_sc, out_max_i, out_cur);
        return;
    }
    int cur = h0, max_sc = h0, max_i = 0;
    int j = 0;
    while (j < N) {
        /* m = first mismatch at or after j, or N if none (no bit is set past N) */
        int m = N;
        for (int w = j >> 6; w < FP_MIS_NWORDS && (w << 6) < N; ++w) {
            uint64_t bits = mis[w];
            if (w == (j >> 6)) bits &= ~0ULL << (j & 63);
            if (bits) { m = (w << 6) + __builtin_ctzll(bits); break; }
        }
        if (m > j) {
            cur += (m - j) * a;
            if (cur >= max_sc) { max_sc = cur; max_i = m; }
        }
        if (m >= N) break;
        cur -= b;
        if (cur < 0) cur = 0;
        if (cur >= max_sc) { max_sc = cur; max_i = m + 1; }
        if (cur == 0) break;
        j = m + 1;
    }
#ifdef BWA_MEM3_DEBUG_UNGAPPED_XCHECK
    {   /* cross-check the run form against the per-base walk it replaces */
        int rmax, rmi, rc;
        ungapped_walk_mis_perbase(mis, N, h0, a, b, &rmax, &rmi, &rc);
        xassert(rmax == max_sc && rmi == max_i && rc == cur,
                "ungapped_walk_mis disagrees with the per-base walk");
    }
#endif
    *out_max_sc = max_sc; *out_max_i = max_i; *out_cur = cur;
}

/* ungapped_walk_score (the byte walk above) over the mismatch bitmask: the same
 * floored walk score. PRECONDITION: no base in [0, N) is ambiguous -- the byte
 * walk then never takes its early break, and both walks see the identical
 * match/mismatch sequence. Named apart from the byte walk because of that
 * precondition. */
static inline int ungapped_walk_score_mis(const uint64_t *mis, int N, int h0, int a, int b)
{
    int max_sc, max_i, cur;
    ungapped_walk_mis(mis, N, h0, a, b, &max_sc, &max_i, &cur);
    return max_sc;
}

static inline int ungapped_analyze(const uint8_t *qs, const uint8_t *rs, int N,
                                    int h0, int a, int b,
                                    int o_min, int e_min,
                                    int x_threshold, int default_w,
                                    int *out_score, int *out_qle,
                                    int *out_gscore, int *out_gtle,
                                    int *out_tight_band)
{
    if (N <= 0 || N > FP_N_MAX || x_threshold < 0) return FP_STATUS_FALLBACK;

    const __m128i v3 = _mm_set1_epi8(3);
    uint64_t mis[FP_MIS_NWORDS] = {0};
    int total_mis = 0;
    int i = 0;
    for (; i + 16 <= N; i += 16) {
        __m128i qv = _mm_loadu_si128((const __m128i *)(qs + i));
        __m128i rv = _mm_loadu_si128((const __m128i *)(rs + i));
        // AMBIG (code >= 4) detected via max > 3. ACGT = 0..3.
        __m128i mxv = _mm_max_epu8(qv, rv);
        if ((unsigned)_mm_movemask_epi8(_mm_cmpgt_epi8(mxv, v3)))
            return FP_STATUS_FALLBACK;  // ambig: give up on both paths
        __m128i eqv = _mm_cmpeq_epi8(qv, rv);
        unsigned mism_mask = (~(unsigned)_mm_movemask_epi8(eqv)) & 0xFFFFu;
        total_mis += __builtin_popcount(mism_mask);
        /* i is 16-aligned, so (i & 63) in {0,16,32,48}: a 16-bit mask shifted
         * by at most 48 fills bits [i&63, (i&63)+15] within a single word and
         * never straddles a 64-bit boundary. */
        mis[i >> 6] |= ((uint64_t)mism_mask) << (i & 63);
    }
    for (; i < N; i++) {
        uint8_t qi = qs[i], ri = rs[i];
        if (qi >= 4 || ri >= 4) return FP_STATUS_FALLBACK;
        if (qi != ri) {
            total_mis++;
            mis[i >> 6] |= (1ULL << (i & 63));
        }
    }

    if (total_mis > x_threshold) {
        /* Only a pair within FP_TIGHT_MAX carries a tight_band (see the define);
         * a longer pair that is not a HIT falls back with tight_band = 0. */
        if (N > FP_TIGHT_MAX) return FP_STATUS_FALLBACK;
        // S in the band proof must be REALIZABLE by an actual offset-0 extension.
        // The tight_band derivation shows every out-of-band alignment (band
        // offset B >= tb) scores <= S, so a band >= tb is sufficient ONLY IF the
        // in-band (offset-0) run actually achieves S. ungapped_walk_score is the
        // floored ungapped score under ksw_extend local-truncation semantics
        // (once the running score hits 0 it stays 0) -- exactly the score the
        // rung-1 banded DP can reach on the diagonal, so it is a valid lower
        // bound on the in-band optimum.
        //
        // A no-floor score (which lets a negative prefix recover via later
        // matches) can EXCEED the floor-killed value the DP actually reaches;
        // feeding that larger S shrinks tb below what is sound, so the certified
        // probe rung or the --adaptive-band stop could skip the wider band a gapped
        // alignment in (default_w, wider] genuinely needs -- a CIGAR/coordinate
        // divergence. The walk is O(N) and runs only on the
        // TIGHT branch after the FP_TIGHT_MAX gate above (N <= 128), so its cost
        // is negligible.
        //
        // Band derivation (see comment above):
        //     B < (min_len·a − S − o_min) / e_min        with S = max_sc − h0
        // For LEFT extensions where the caller gates on len1 ≥ len2,
        // min_len = N. Substituting:
        //     numerator = N·a − (max_sc_proof − h0) − o_min
        //               = N·a + h0 − max_sc_proof − o_min
        /* Bitmask form of ungapped_walk_score(qs, rs, N, h0, a, b): no base in
         * [0, N) is ambiguous here (checked above), so the byte walk never
         * breaks early and both return the same floored score. */
        int max_sc_proof = ungapped_walk_score_mis(mis, N, h0, a, b);
        int64_t numerator = (int64_t)N * a + h0 - max_sc_proof - o_min;
        int band;
        if (numerator <= 0) {
            // Ungapped-optimal (same "tight_band = 0 sentinel" semantic as
            // walk-derived numerator ≤ 0: let SW run with fallback width).
            band = 0;
        } else if (e_min <= 0) {
            band = default_w;
        } else {
            band = (int)((numerator + e_min - 1) / e_min);
            // The band proof only excludes diagonal offsets >= band. When the
            // proven band exceeds default_w, clamping to default_w would falsely
            // certify that a width-default_w run is complete, letting the
            // certified probe rung or the --adaptive-band stop skip the wider band
            // a gapped alignment in the unproven window (default_w, band)
            // genuinely needs. Emit the tight_band = 0 fallback sentinel instead,
            // so no tight_band consumer acts on the unproven band.
            if (band > default_w) band = 0;
        }
        *out_tight_band = band;
        // Outputs unused on TIGHT; set sane values for any future caller.
        *out_score  = max_sc_proof;
        *out_qle    = 0;
        *out_gscore = 0;
        *out_gtle   = N;
        return FP_STATUS_TIGHT;
    }

    // Closed-form fast path: total_mis == 0 (perfect-match HIT).
    //
    // With no mismatches the scalar walk below is deterministic:
    //   cur starts at h0 and only increases (+a per iter), so it never
    //   touches the `cur == 0` early-skip after the first iteration.
    //   max_sc rises to h0 + N*a; max_i ends at N (the >= tie-break
    //   updates max_i on every step).
    //
    // The h0 > 0 guard preserves bit-identicality with the loop: when
    // h0 == 0 the loop's `if (cur == 0) continue` inhibits all updates,
    // leaving max_sc = 0 / max_i = 0; the closed form would instead
    // compute max_sc = N*a, breaking parity. h0 == 0 is pathological in
    // practice (caller passes a seed score) but the guard is a single
    // compare and free.
    if (total_mis == 0 && h0 > 0) {
        int score = h0 + N * a;
        *out_score      = score;
        *out_qle        = N;
        *out_gscore     = score;
        *out_gtle       = N;
        *out_tight_band = 0;   // unused on HIT (SW skipped)
        return FP_STATUS_HIT;
    }

    // HIT candidate: run the scalar walk for precise qle / gscore / max_sc.
    //
    // MAIN_CODE* local-SW semantics: once cur==0 in the ungapped path it
    // stays 0 (no e/f restart).
    //
    // Tie-break: SW's maxRS tracker (bandedSWA.cpp MAIN_CODE) updates the
    // position on BOTH strictly-greater and tied equal-to-current
    // comparisons — equivalent to "pick the rightmost position where the
    // max was achieved". We must mirror that (use >=) or qle/tle diverge
    // from SW on tied-score walks, breaking byte-identical SAM.
    int cur, max_sc, max_i;
    ungapped_walk_mis(mis, N, h0, a, b, &max_sc, &max_i, &cur);

    *out_score      = max_sc;
    *out_qle        = max_i;
    *out_gscore     = cur;
    *out_gtle       = N;
    *out_tight_band = 0;   // unused on HIT (SW skipped)
    return FP_STATUS_HIT;
}

#endif /* BWA_MEM3_UNGAPPED_EXT_H */
