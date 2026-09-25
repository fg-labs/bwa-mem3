/* BWA3_RESCUE_* environment toggles: the one list of them, and the helpers that read them.
 *
 * Every toggle is a developer lever for A/B-ing the mate-rescue path; none is a user option, and
 * every default is the production setting. Every alternative is exact (argued byte-identical
 * where it is defined), so no setting changes the output -- only the stats flag and the timing.
 * Flags are ON unless their value starts with '0' (rescue_env_flag); the stats flag is OFF unless
 * its value starts with '1' (rescue_env_opt_in); numbers are atoi'd (rescue_env_int).
 *
 *   kswv kernels (src/kswv.cpp; read once per batch, so a unit test can flip them in-process)
 *     BWA3_RESCUE_USQADD         1     u8: one saturating add per cell instead of the biased
 *                                      add / subtract pair
 *     BWA3_RESCUE_ROWPAIR        1     NEON: sweep two rows at a time (0: one row)
 *     BWA3_RESCUE_LAZYQE         1     NEON two-row sweep: recover the query end lazily
 *                                      (0: inline argmax)
 *     BWA3_RESCUE_FSCAN          1     every SIMD body: the 11-op G-based cell (any value not
 *                                      starting with '0' is on). The banded rescue kernel reads
 *                                      it once as a kernel number: 0 the original cell, 1 the
 *                                      fused cell, 2 (default) the fused cell on two rows per
 *                                      step with the direct qe scan (rb_dp_wave2)
 *   exact window pruning (src/bwamem_pair.cpp, rescue_prune.h; read once)
 *     BWA3_RESCUE_PRUNE          1     prune rescue windows by exact 5-mer bounds
 *     BWA3_RESCUE_PRUNE_MAX_HITS auto  skip the filter on windows with more 5-mer hits than
 *                                      this; auto is 1000 on aarch64 when banding is on and
 *                                      minsc == 19 (the NEON filter), else 400 (always on x86)
 *   banded rescue (src/rescue_band.{h,cpp}, NEON or AVX2 kernel; read once)
 *     BWA3_RESCUE_BAND           1     run pruned pass-0 jobs in diagonal bands (0: kswv on the
 *                                      hull, and MAX_HITS auto back to 400)
 *     BWA3_RESCUE_BAND_COST      85    band a pass-0 parent iff its band cells are < this % of
 *                                      the hull's
 *     BWA3_RESCUE_BAND_R2        1     run round 2 banded (0: kswv on the hull)
 *     BWA3_RESCUE_BAND_TIGHT     8     delta of the tight top band (0: off)
 *     BWA3_RESCUE_BAND_P1        2     pass-1 banding: 0 none, 1 banded parents only, 2 every
 *                                      eligible 8-bit default-scoring job
 *     BWA3_RESCUE_BAND_P1_COST   130   band a pass-1 job iff its per-row cells are < this % of
 *                                      kswv's
 *     BWA3_RESCUE_BAND_SHIFT     1     shift each narrower lane's spare diagonals below its
 *                                      band so the query offsets of a 16-lane group align
 *   mate-rescue dedup (src/bwamem_pair.cpp; read once)
 *     BWA3_RESCUE_DEDUP_SKIP     1     skip a post-rescue dedup that is provably a no-op (nothing
 *                                      added since a dedup that reported a fixed point)
 *   diagnostics
 *     BWA3_RESCUE_PRUNE_STATS    0     1: print the [RESCUE_PRUNE] and [RESCUE_BAND] counters
 *                                      and stage times at exit (both passes, not only pruning),
 *                                      plus the dedup_run / dedup_skip counts
 *
 * test/rescue_band_harness.cpp reads the same variables (plus its own RB_* knobs), so a harness
 * run and a whole-aligner run with one environment take the same paths. */
#ifndef BWA_MEM3_RESCUE_ENV_H
#define BWA_MEM3_RESCUE_ENV_H

#include <stdlib.h>

/* A default-ON flag: off only when the variable is set and starts with '0'. */
static inline bool rescue_env_flag(const char *var)
{
    const char *e = getenv(var);
    return !e || e[0] != '0';
}

/* A default-OFF flag: on only when the variable is set and starts with '1'. */
static inline bool rescue_env_opt_in(const char *var)
{
    const char *e = getenv(var);
    return e && e[0] == '1';
}

/* An integer knob: the variable's atoi value when set, else dflt. */
static inline int rescue_env_int(const char *var, int dflt)
{
    const char *e = getenv(var);
    return e ? atoi(e) : dflt;
}

#endif
