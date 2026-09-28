/* Shared readers for the BWA3_RESCUE_* knobs (kswv.cpp, rescue_prune.h, rescue_band.h).
 *
 * rescue_env_int: an integer tuning knob, read once per process into a function-local static by
 * its caller. It must be a non-negative decimal integer. An unset or empty knob silently takes the
 * default; anything else (a sign, leading or trailing characters, overflow) is reported to stderr
 * and the default used, so a typo cannot silently set a knob.
 *
 * rescue_env_on: an on/off toggle, default on; a value starting with '0' turns it off and
 * anything else leaves it on. Not reported and not cached: kswv reads its toggles on every call so
 * a unit test can flip them in-process, and a per-call report would repeat once per batch;
 * callers that want the value once wrap it in a function-local static.
 *
 * rescue_env_opt_in: a diagnostic switch, default off; only a value starting with '1' turns it on.
 *
 * The knobs, one per line between the markers: name, default, meaning (continuation lines are
 * indented past the name column). scripts/rescue_knobs.sh renders this list as
 * docs/_generated/rescue/knobs.md, and test/regression/rescue_docs_lint.sh checks it against the
 * knobs src/ reads, their defaults there, and the user-facing table in
 * docs/src/whats-different/performance.md. Design notes: docs/src/developer-guide/rescue.md.
 *
 * rescue-knobs:begin
 * BWA3_RESCUE_PRUNE           1     exact K-mer pruning of rescue windows (aarch64 and x86 AVX2 /
 *                                   AVX-512BW builds; scope: rescue-pruning.md); 0 also turns
 *                                   the banded passes off
 * BWA3_RESCUE_PRUNE_MAX_HITS  auto  keep the full window above this many K-mer hits; auto is 1000
 *                                   on aarch64 with banding on, unless a --meth run leaves its
 *                                   pruned windows unbanded (EM-seq by default); else 400
 * BWA3_RESCUE_PRUNE_KMAX      auto  the longest K-mer the filter may use where the scoring admits
 *                                   it (5 to 8); auto is 8 on aarch64 and the x86 AVX2 /
 *                                   AVX-512BW builds, else 5; --meth always 5
 * BWA3_RESCUE_PRUNE_REL       1     --meth genomic / neutral scoring, aarch64 (x86 --meth is not
 *                                   pruned): 1 filter TAPS under the exact relation (EM-seq and
 *                                   collapsed scoring match converted copies), 0 converted copies
 *                                   only (TAPS not pruned), 2 or more the relation for every
 *                                   genomic / neutral run
 * BWA3_RESCUE_PRUNE_STATS     0     1 prints the RESCUE_PRUNE and RESCUE_BAND counters and stage
 *                                   times at exit (counted per thread: no shared-counter
 *                                   contention in an instrumented run)
 * BWA3_RESCUE_DEDUP_SKIP      1     skip a post-rescue dedup proven to be a no-op, and add a
 *                                   single new region in one pass where that is provably exact
 * BWA3_RESCUE_REPEAT          1     a rescue job repeating one of the last eight filtered jobs byte
 *                                   for byte reads that job's result instead of being enqueued again
 * BWA3_RESCUE_FSCAN           1     the 11-op kswv rescue cell in every SIMD body; 0 also selects
 *                                   the banded DP's original cell
 * BWA3_RESCUE_USQADD          1     u8 kswv: one saturating add per cell instead of the biased
 *                                   add / subtract pair (NEON; the AVX2 and AVX-512BW FScan bodies
 *                                   do it in the signed H - 128 domain, at open-plus-extend sums
 *                                   up to 127)
 * BWA3_RESCUE_ROWPAIR         1     NEON kswv (8- and 16-bit): sweep two target rows per pass
 * BWA3_RESCUE_LAZYQE          1     NEON kswv two-row sweep: recover the query end after the row
 *                                   instead of inline
 * BWA3_RESCUE_BAND            1     run rescue jobs as diagonal bands (NEON or AVX2 kernel)
 * BWA3_RESCUE_BAND_COST       auto  band a pass-0 parent iff its band cells cost less than this %
 *                                   of the hull's; auto is 0 (no pass-0 banding) where kswv runs
 *                                   at the AVX-512BW tier, whose 64-lane hull undercuts the
 *                                   32-lane band kernel, else 85
 * BWA3_RESCUE_BAND_R2         1     run round 2 banded (0: kswv on the hull)
 * BWA3_RESCUE_BAND_TIGHT      8     threshold offset of the tight top band (0: off)
 * BWA3_RESCUE_BAND_P1         2     pass-1 banding: 0 none, 1 banded parents only, 2 every
 *                                   eligible 8-bit job
 * BWA3_RESCUE_BAND_P1_COST    130   band a pass-1 job iff its per-row cells cost less than this %
 *                                   of kswv's
 * BWA3_RESCUE_BAND_KERNEL     2     banded-DP kernel while BWA3_RESCUE_FSCAN is on: 0 original
 *                                   cell, 1 fused cell, 2 fused cell on two rows per step
 * BWA3_RESCUE_BAND_SHIFT      1     shift each lane's band so a lane group's query offsets align
 * BWA3_RESCUE_BAND_METH       1     --meth banding (both passes, the group's matrix): 1 except on top
 *                                   of converted-copy pruning (EM-seq), 0 never, 2 or more always
 * rescue-knobs:end */
#ifndef BWA_MEM3_RESCUE_ENV_H
#define BWA_MEM3_RESCUE_ENV_H

#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

static inline int rescue_env_int(const char *name, int dflt)
{
    const char *e = getenv(name);
    if (!e || !*e) return dflt;
    char *end = NULL;
    errno = 0;
    /* strtol skips leading whitespace and takes a sign, so require a digit first. */
    const long x = (e[0] >= '0' && e[0] <= '9') ? strtol(e, &end, 10) : -1;
    if (x < 0 || errno || *end || x > INT_MAX) {
        fprintf(stderr, "ERROR: %s=\"%s\" is not a non-negative integer; using %d.\n", name, e, dflt);
        return dflt;
    }
    return (int)x;
}

static inline bool rescue_env_on(const char *name)
{
    const char *e = getenv(name);
    return !e || e[0] != '0';
}

static inline bool rescue_env_opt_in(const char *name)
{
    const char *e = getenv(name);
    return e && e[0] == '1';
}

#endif
