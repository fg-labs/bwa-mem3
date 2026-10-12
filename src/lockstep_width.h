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

/* Third-pass (bwtSeedStrategy) re-seeding lockstep WIDTH -- how many reads'
 * forward-extension walks the driver overlaps so their cp_occ misses issue
 * together (distinct from the on/off switch g_bwtseed_lockstep below, and from
 * the phase-2 SMEM width above). The optimum is the core's usable memory-level
 * parallelism; measured across five microarchitectures on a WGS slice (at
 * threads <= physical cores, holding the regime fixed across the sweep) it
 * converges on ~24:
 * x86 Sapphire Rapids (~-1.8% vs 8, flat plateau to 64) and Zen 4 (peak at 24,
 * rolls off after -- which caps a portable default at 24), ARM Graviton4
 * (~-0.8%), and Apple M3 Ultra (~-2.2%, still improving) and M2 Max (flat, no
 * regression). So 24 is a win or neutral everywhere and the portable default;
 * the prior compile-time 8 was under-tuned on every architecture. Hosts that
 * want more (SPR, M3 Ultra) are served by BWA3_BWTSEED_LOCKSTEP_N or a host-tuned
 * default (AMD Zen 5 gets 32; see "Host-tuned default widths" below). Byte-identical across widths (batching only; the bwtseed lockstep
 * parity harness pins the emission order). */
#ifndef BWTSEED_LOCKSTEP_N
#define BWTSEED_LOCKSTEP_N 24
#endif
/* Upper bound on the runtime width; sizes the driver's on-stack
 * BwtSeedSlot slots[BWTSEED_LOCKSTEP_N_MAX] array and caps any pin/calibration. */
#ifndef BWTSEED_LOCKSTEP_N_MAX
#define BWTSEED_LOCKSTEP_N_MAX 64
#endif
#if !(BWTSEED_LOCKSTEP_N >= 1 && BWTSEED_LOCKSTEP_N <= BWTSEED_LOCKSTEP_N_MAX && BWTSEED_LOCKSTEP_N_MAX <= SMEM_LOCKSTEP_N_MAX_CAP)
#error "require 1 <= BWTSEED_LOCKSTEP_N <= BWTSEED_LOCKSTEP_N_MAX <= 256 (stack-array storage-safe)"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* Runtime phase-2 lockstep width, resolved once at startup by
 * bwa3_init_smem_lockstep_width (env pin, opt-in probe, host-tuned default, or
 * compile-time default). Initialized to SMEM_LOCKSTEP_N so a binary that never
 * calls the initializer behaves exactly like the compile-time constant. Read
 * (not written) on the seeding hot path. */
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
 *     0  env is unset or empty -- no override; the caller probes if
 *        BWA3_SMEM_LOCKSTEP_PROBE is set, else uses the host default (the
 *        compile-time default, or the host-tuned value on AMD Zen 5);
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
 * or 0 if the array is absent or too small to measure (the initializer then
 * falls back to the host default). Cost is host-dependent: ~10-20 ms on low-latency memory (e.g. Apple
 * Silicon), ~0.4 s on server-class DRAM (a ~100 ns dependent-load latency over
 * the sweep). This is the raw measurement primitive: it retains only its input
 * guard (returns 0 when base == NULL or n_blocks < 4096) and has no one-shot or
 * opt-in guard, so a direct caller with a large enough array runs the full sweep
 * on every call. The "at most once at startup, and only when opted in" contract
 * is enforced by the sole intended caller, bwa3_init_smem_lockstep_width
 * (idempotent, gated on BWA3_SMEM_LOCKSTEP_PROBE) -- see below. */
int32_t bwa3_measure_mlp(const void *base, int64_t n_blocks,
                         size_t stride, size_t word_off);

/* ---- Host-tuned default widths ----------------------------------------------
 *
 * With no env pin and no opt-in probe, both lockstep widths fall back to a
 * per-host default. Every host gets the compile-time SMEM_LOCKSTEP_N /
 * BWTSEED_LOCKSTEP_N except AMD Zen 5, which gets 32 / 32. Zen 5 is matched by
 * vendor, family AND model: AMD family 1Ah also holds Zen 6, so the family alone
 * is not enough. The model ranges follow the Linux kernel's Zen generation table
 * (bsp_init_amd, arch/x86/kernel/cpu/amd.c): family 1Ah models 00h-2Fh, 40h-4Fh,
 * 60h-7Fh and D0h-D7h are Zen 5; 50h-5Fh, 80h-AFh, C0h-CFh and D8h-EFh are Zen 6.
 * Zen 6 and every unlisted family-1Ah model keep the compile-time defaults (not
 * measured). Zen 5c shares the Zen 5 ranges, so it gets 32 / 32 (also not
 * measured).
 *
 * Measured whole-aligner wall, 5M-pair WGS (HG00096) and WES (HG00100) slices,
 * hg38, -t 16 on 16 cores without SMT, clang 19, paired interleaved reps, the
 * 32/32 host default vs the 16/24 default build (negative = faster):
 *   AMD Zen 5   m8a.4xlarge  avx512bw  WGS -1.26 +/- 0.16% (12 reps)  WES -1.59 +/- 0.48% (12 reps)
 * An earlier sweep pinning the widths through the env vars measured the same
 * pair at WGS -1.44 +/- 0.29% and WES -1.22 +/- 0.28% (10 reps each), and found
 * 32/32 the most consistent pair across WGS and WES among SMEM {16, 24, 32} x
 * bwtSeed {24, 32, 48}. The other hosts measured do not want it:
 *   AMD Zen 4   m7a.4xlarge  avx512bw  SMEM 24: WGS +0.4%, SMEM 32: WGS +0.7% (6 reps, within noise)
 *   Graviton4   m8g.4xlarge  neon      WGS flat from SMEM 16 to 32; WES best at 16/24
 * so the default is keyed on the CPU, not on x86 or the SIMD tier. Intel was not
 * measured. Width is byte-identical by contract (batching only), so this is
 * purely a scheduling choice.
 *
 * Precedence, highest first: an explicit BWA3_SMEM_LOCKSTEP_N /
 * BWA3_BWTSEED_LOCKSTEP_N pin, then (SMEM only) the opt-in startup probe, then
 * this host default, then the compile-time default. The host default only ever
 * RAISES a width: it is clamped to [compile-time default, *_LOCKSTEP_N_MAX], so
 * a build that overrides either bound with -D never gets a width outside its
 * own range (a -DSMEM_LOCKSTEP_N_MAX=24 build gets 24 on Zen 5, a
 * -DSMEM_LOCKSTEP_N=48 build keeps 48).
 *
 * Memory: each worker thread keeps grow-only lockstep scratch proportional to
 * width x max read length (SMEM: two SMEM buffers per slot; bwtSeed: one), so
 * on Zen 5 that scratch is 2x (SMEM) and 4/3x (bwtSeed) the 16/24 default. It
 * is small for short reads; for very long reads, pin BWA3_SMEM_LOCKSTEP_N /
 * BWA3_BWTSEED_LOCKSTEP_N lower to cap it. */
#ifndef BWA3_ZEN5_SMEM_LOCKSTEP_N
#define BWA3_ZEN5_SMEM_LOCKSTEP_N 32
#endif
#ifndef BWA3_ZEN5_BWTSEED_LOCKSTEP_N
#define BWA3_ZEN5_BWTSEED_LOCKSTEP_N 32
#endif

/* Host-tuned default widths for a given CPU identity (vendor, display family,
 * display model; see bwa3_host_cpu_identity in simd_dispatch.h), clamped to
 * [compile-time default, *_LOCKSTEP_N_MAX]. Return the compile-time default for
 * every host without a tuned value. Pure, so the policy is unit-testable with
 * an injected identity. */
int32_t bwa3_smem_host_default_width(int is_amd, uint32_t family, uint32_t model);
int32_t bwa3_bwtseed_host_default_width(int is_amd, uint32_t family, uint32_t model);

/* Where a resolved lockstep width came from (for the -v 3 startup log). */
enum {
    BWA3_LOCKSTEP_WIDTH_SRC_DEFAULT = 0,  /* compile-time default */
    BWA3_LOCKSTEP_WIDTH_SRC_HOST    = 1,  /* host-tuned default (AMD Zen 5) */
    BWA3_LOCKSTEP_WIDTH_SRC_PROBE   = 2,  /* opt-in startup MLP probe (SMEM only) */
    BWA3_LOCKSTEP_WIDTH_SRC_ENV     = 3   /* explicit env pin */
};

/* Log label of a BWA3_LOCKSTEP_WIDTH_SRC_* value: "default", "host default: AMD
 * Zen 5", "startup probe", "env pin"; any other value returns "default". */
const char *bwa3_lockstep_width_source_name(int src);

/* Source of the installed g_smem_lockstep_n / g_bwtseed_lockstep_n, recorded by
 * the initializers below (BWA3_LOCKSTEP_WIDTH_SRC_DEFAULT until they run). */
int bwa3_smem_lockstep_n_source(void);
int bwa3_bwtseed_lockstep_n_source(void);

/* Resolve and install g_smem_lockstep_n once. Idempotent (subsequent calls are
 * no-ops). Resolution order:
 *   1. BWA3_SMEM_LOCKSTEP_N=<n>   -- explicit pin, taken alone.
 *   2. BWA3_SMEM_LOCKSTEP_PROBE=1 -- opt into the startup MLP probe, which
 *      chases the cp_occ array (base/n_blocks/stride/word_off) to self-calibrate
 *      the width for this host. Truthy opt-in (bwa3_lockstep_probe_enabled): any
 *      value other than unset, empty, or "0" enables it, so unset the variable
 *      (or set it to 0) to disable -- BWA3_SMEM_LOCKSTEP_PROBE=0 does NOT enable.
 *   3. the host default (bwa3_smem_host_default_width): the host-tuned value
 *      on AMD Zen 5 (32), else the compile-time SMEM_LOCKSTEP_N. Also the
 *      fallback when the opted-in probe cannot measure (returns 0).
 * The probe is OPT-IN: it costs a startup pass and picks the per-access latency
 * knee, which sits at or below the whole-aligner wall optimum (on Zen 5 it picked
 * 16, 16 and 32 on three startups against a wall optimum of 32), so the measured
 * host default is the better mechanism. An invalid env value is reported
 * (ERROR to stderr) and ignored, resolving as if unset.
 * Call after the index is loaded and before the seeding workers spawn. */
void bwa3_init_smem_lockstep_width(const void *base, int64_t n_blocks,
                                   size_t stride, size_t word_off);


/* ---- Third-pass bwtseed re-seeding lockstep: on/off policy -------------------
 *
 * FMI_search::bwtSeedStrategyAllPosOneThread_lockstep overlaps BWTSEED_LOCKSTEP_N
 * reads' forward-extension walks so their cp_occ cache misses issue together. It
 * is on by default on every platform and at every thread count;
 * BWA3_BWTSEED_LOCKSTEP=0|1 pins it (test/regression/bwtseed_lockstep_driver_identity.sh
 * pins the scalar driver in CI and diffs it against the default).
 *
 * x86 used to keep it off whenever -t exceeded the physical core count, on the
 * theory that a busy SMT sibling already hides the same latency. That was measured
 * at lockstep width 8 (wall bimodal on Zen 3 at -t = 2x cores). At the width-24
 * default it does not hold: with both siblings busy the sibling hides little, and
 * the scalar driver (bwtSeedStrategyAllPosOneThread) was ~13% of cycles. Measured
 * wall, lockstep off -> on, 5M-pair WGS / WES slices (hg38, clang 19, AVX2 build),
 * 5-7 interleaved reps per arm, rep spread <= 2%, no bimodality:
 *   AMD Zen 3              c6a.4xlarge  -t 16,  8 cores   WGS -8.2%   WES -7.2%
 *   Intel Sapphire Rapids  c7i.4xlarge  -t 16,  8 cores   WGS -9.6%   WES -8.2%
 *   AMD Zen 3              c6a.8xlarge  -t 32, 16 cores   WGS -8.2%   WES -7.2%
 *   Intel Sapphire Rapids  c7i.8xlarge  -t 32, 16 cores   WGS -8.9%   WES -8.3%
 * At -t = physical cores (one thread per core) it was already on, and wins by
 * more there (-13% WGS wall at -t 8 on both 8-core hosts). Not measured: SMT
 * hosts of other microarchitectures (e.g. Zen 4 with SMT) and oversubscription
 * beyond -t = 2x physical cores were not measured for speed; output identity
 * there rests on the driver's SMEM-order parity design (pinned by the parity
 * harness and the driver-identity regression), not on a measurement there.
 *
 * Output is byte-identical either way (same SMEM emission order; the lockstep
 * parity harness pins it) -- this is a scheduling choice, never a result one. */

/* Runtime switch read by mem_collect_smem: nonzero selects the lockstep driver,
 * zero the scalar one. Initialized to the compiled default (on) so a binary that
 * never runs bwa3_init_bwtseed_lockstep behaves the same. Read (not written) on the seeding hot
 * path. */
extern int32_t g_bwtseed_lockstep;

/* Runtime third-pass lockstep WIDTH, resolved once at startup by
 * bwa3_init_bwtseed_lockstep_width (env pin, host-tuned default, or compile-time
 * default). Initialized to the compile-time BWTSEED_LOCKSTEP_N so a binary that
 * never calls the initializer behaves exactly like the constant. Read (not written) on the seeding hot path by
 * bwtSeedStrategyAllPosOneThread_lockstep. */
extern int32_t g_bwtseed_lockstep_n;

/* Classify a BWA3_BWTSEED_LOCKSTEP_N override (env may be NULL/empty), WITHOUT
 * touching global state. Returns > 0 the width to pin (clamped up to
 * BWTSEED_LOCKSTEP_N_MAX), 0 when unset/empty (use the host default), or -1 when
 * set but malformed / non-positive / overflowed (caller reports it and uses the
 * host default).
 * Pure and side-effect-free, so it is unit-testable. */
int32_t bwa3_bwtseed_width_parse_env(const char *env);

/* Resolve and install g_bwtseed_lockstep_n once (idempotent). Resolution order:
 * BWA3_BWTSEED_LOCKSTEP_N pin -> host default (bwa3_bwtseed_host_default_width:
 * 32 on AMD Zen 5, else the compile-time BWTSEED_LOCKSTEP_N, 24). Call
 * before the seeding workers spawn; an invalid env value is reported to stderr and
 * resolved as if unset. */
void bwa3_init_bwtseed_lockstep_width(void);

/* Classify a BWA3_BWTSEED_LOCKSTEP override (env may be NULL/empty) without
 * touching global state:  1 = pin on, 0 = pin off, -1 = unset/empty (use the
 * compiled default), -2 = set but not "0"/"1" (invalid: the caller reports it and
 * uses the default). Pure, so it is unit-testable. */
int32_t bwa3_bwtseed_lockstep_parse_env(const char *env);

/* Resolve and install g_bwtseed_lockstep, re-resolving on every call so a library
 * caller that runs main_mem repeatedly gets each run's own pin, never a leaked
 * earlier one: a 0/1 pin is taken alone, unset and invalid both install the
 * compile-time default (BWA3_BWTSEED_LOCKSTEP_DEFAULT in lockstep_width.cpp, the
 * single place it is written down). An invalid env value is reported (ERROR to
 * stderr). Call on the main thread before the seeding workers spawn. */
void bwa3_init_bwtseed_lockstep(void);

#ifdef __cplusplus
}
#endif

#endif /* LOCKSTEP_WIDTH_H */
