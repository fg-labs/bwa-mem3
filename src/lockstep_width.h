#ifndef LOCKSTEP_WIDTH_H
#define LOCKSTEP_WIDTH_H

#include <stdint.h>
#include <stddef.h>  /* size_t */

/* Phase-2 SMEM lockstep width.
 *
 * SMEM_LOCKSTEP_N is three things at once: the compile-time default, the
 * regression floor for the runtime auto-tune (the startup probe may only
 * RAISE the width, never pick below this), and the enable guard
 * (`#if SMEM_LOCKSTEP_N > 1`) that selects the lockstep driver over the
 * scalar path at build time. Overridable at build time (a `-D` on the
 * command line wins over this default). See getSMEMsOnePosOneThread_lockstep
 * in FMI_search.cpp for the driver this width feeds. */
#ifndef SMEM_LOCKSTEP_N
#define SMEM_LOCKSTEP_N 16
#endif

/* Upper bound on the runtime width. Sizes the driver's on-stack slot array
 * (BatchSlot slots[SMEM_LOCKSTEP_N_MAX]) and caps the value the startup
 * memory-level-parallelism probe may select. Raise only in step with the
 * stack-array cost (~80 B per slot). */
#ifndef SMEM_LOCKSTEP_N_MAX
#define SMEM_LOCKSTEP_N_MAX 64
#endif

/* Both widths are `-D`-overridable, and the driver indexes an on-stack
 * BatchSlot slots[SMEM_LOCKSTEP_N_MAX] array with a runtime width in
 * [SMEM_LOCKSTEP_N, SMEM_LOCKSTEP_N_MAX]. An override that inverts the range
 * (SMEM_LOCKSTEP_N > SMEM_LOCKSTEP_N_MAX) or a non-positive bound would index
 * out of bounds, so reject it at compile time rather than miscompile.
 *
 * The upper bound is a STORAGE-safe cap, not merely INT32_MAX: both the
 * driver's BatchSlot slots[SMEM_LOCKSTEP_N_MAX] (~100 B/slot, FMI_search.cpp)
 * and the probe's uint64_t cur[SMEM_LOCKSTEP_N_MAX] (8 B/slot, lockstep_width.cpp)
 * are fixed-size STACK arrays sized by this bound. An INT32_MAX bound is
 * representable in int32_t but would demand a multi-gigabyte stack frame and
 * crash before the probe ever runs. Cap the override at a width whose
 * worst-case frame (~28 KB at 256) is safe on any thread stack, while still
 * leaving 4x headroom over the shipped default 64. */
#define SMEM_LOCKSTEP_N_MAX_CAP 256
#if !(SMEM_LOCKSTEP_N >= 1 && SMEM_LOCKSTEP_N <= SMEM_LOCKSTEP_N_MAX && SMEM_LOCKSTEP_N_MAX <= SMEM_LOCKSTEP_N_MAX_CAP)
#error "require 1 <= SMEM_LOCKSTEP_N <= SMEM_LOCKSTEP_N_MAX <= SMEM_LOCKSTEP_N_MAX_CAP (256; stack-array storage-safe)"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* Runtime phase-2 lockstep width, resolved once at startup. Initialized to
 * SMEM_LOCKSTEP_N so a binary that never runs the probe behaves exactly like
 * the compile-time constant. Read (not written) on the seeding hot path. */
extern int32_t g_smem_lockstep_n;

/* Clamp a raw memory-level-parallelism estimate to the usable width range
 * [SMEM_LOCKSTEP_N, SMEM_LOCKSTEP_N_MAX]. A non-positive estimate (a probe
 * that failed to measure) returns the floor. The floor is what guarantees the
 * auto-tune never selects a width below the shipped default. */
int32_t bwa3_lockstep_width_from_probe(int32_t raw_mlp);

/* Classify a BWA3_SMEM_LOCKSTEP_N override (env may be NULL/empty) WITHOUT
 * running the probe or touching global state. Returns:
 *   > 0  the width to pin, clamped up to SMEM_LOCKSTEP_N_MAX -- a valid positive
 *        integer override (it MAY sit below the floor, an explicit operator/gate
 *        escape hatch); the caller skips the probe and uses this;
 *     0  env is unset or empty -- no override; the caller keeps the compile-time
 *        default (and probes only if BWA3_SMEM_LOCKSTEP_PROBE is set);
 *    -1  env is set but malformed, non-positive, or overflowed -- an invalid
 *        override; the caller reports it and resolves it as if unset.
 * Pure and side-effect-free (emits no diagnostic), so it is unit-testable; the
 * ERROR message for the -1 case is the initializer's responsibility. */
int32_t bwa3_lockstep_width_parse_env(const char *env);

/* Gate the startup MLP probe opt-in on a BWA3_SMEM_LOCKSTEP_PROBE value (env may
 * be NULL/empty). Returns nonzero (enable the probe) only for a truthy value:
 * present, non-empty, and not "0"; unset, empty, and "0" return 0 (disabled).
 * Pure and side-effect-free, so the gating is unit-testable independently of the
 * initializer that consults it. */
int bwa3_lockstep_probe_enabled(const char *env);

/* Measure the core's memory-level parallelism: the number of independent
 * dependent-load chains it keeps outstanding before per-access latency stops
 * falling. The chase walks a value-dependent pseudo-random cycle THROUGH the
 * FM-index's own cp_occ checkpoint array -- gigabytes, already resident, the
 * same array the seeding walk hammers, with the same TLB behaviour -- rather
 * than a scratch buffer (a fresh buffer measures whatever cache level it fits
 * in and over-picks). `base`/`n_blocks`/`stride` describe cp_occ as an opaque
 * array of `n_blocks` fixed-size blocks; `word_off` is the byte offset of a
 * 64-bit word to read from each block. Returns the knee width (raw, pre-clamp),
 * or 0 if the array is absent or too small to measure (caller floors to the
 * default). Cost is host-dependent: ~10-20 ms on low-latency memory (e.g. Apple
 * Silicon), ~0.4 s on server-class DRAM (a ~100 ns dependent-load latency over
 * the sweep). This is the raw measurement primitive: it retains only its input
 * guard (returns 0 when base == NULL or n_blocks < 4096) and has no one-shot or
 * opt-in guard, so a direct caller with a large enough array runs the full sweep
 * on every call. The "at most once at startup, and only when opted in" contract
 * is enforced by the sole intended caller, bwa3_init_smem_lockstep_width
 * (idempotent, gated on BWA3_SMEM_LOCKSTEP_PROBE) -- see below. */
int32_t bwa3_measure_mlp(const void *base, int64_t n_blocks,
                         size_t stride, size_t word_off);

/* Resolve and install g_smem_lockstep_n once. Idempotent (subsequent calls are
 * no-ops). Resolution order:
 *   1. BWA3_SMEM_LOCKSTEP_N=<n>   -- explicit pin, taken alone.
 *   2. BWA3_SMEM_LOCKSTEP_PROBE=1 -- opt into the startup MLP probe, which
 *      chases the cp_occ array (base/n_blocks/stride/word_off) to self-calibrate
 *      the width for this host. Truthy opt-in (bwa3_lockstep_probe_enabled): any
 *      value other than unset, empty, or "0" enables it, so unset the variable
 *      (or set it to 0) to disable -- BWA3_SMEM_LOCKSTEP_PROBE=0 does NOT enable.
 *   3. neither set (the default)  -- keep the compile-time SMEM_LOCKSTEP_N.
 * The probe is OPT-IN because the width is a measured ~1% flat knob end-to-end
 * across architectures, core counts, and page sizes, so paying its startup cost
 * on every run buys no throughput; the default constant is what the probe
 * converges to on all benchmarked hardware. An invalid env value is reported
 * (ERROR to stderr) and ignored, resolving as if unset.
 * Call after the index is loaded and before the seeding workers spawn. */
void bwa3_init_smem_lockstep_width(const void *base, int64_t n_blocks,
                                   size_t stride, size_t word_off);

#ifdef __cplusplus
}
#endif

#endif /* LOCKSTEP_WIDTH_H */
