# `BASELINE_ARCH=avx512bw` build flag

This page documents the empirical perf characterization of building
bwa-mem3 with `BASELINE_ARCH=avx512bw` and the
`-mprefer-vector-width=256` mitigation that ships as part of that
target.

## Background: BASELINE_ARCH

bwa-mem3 ships a single x86 binary with all five SIMD tiers
(sse41 / sse42 / avx / avx2 / avx512bw) compiled for the hand-tuned
kernel TUs (`KERNEL_SRCS` in the Makefile: `bandedSWA`, `kswv`, `ksw`,
`sam_encode`). The runtime dispatcher in `src/simd_dispatch.cpp` picks
the right tier per kernel call based on `__builtin_cpu_supports`.

Everything outside `KERNEL_SRCS` — `bwamem.cpp`, `bwamem_pair.cpp`,
`FMI_search.cpp`, `fastmap.cpp`, `bntseq.cpp`, etc. — is compiled
**once** at the tier set by the `BASELINE_ARCH` Makefile variable
(default: `avx2`). The compiler can auto-vectorize loops in those TUs
at up to that tier's width.

PR #84 raised the default from `sse41` to `avx2` after measuring
~10-15% wall-time gains on AVX2 hosts (c6a, etc.) when the
auto-vectorizer could finally widen hot non-kernel loops to 256-bit.

## The naive expectation: avx2 → avx512bw should give another tier

Following PR #84's logic, you might expect `BASELINE_ARCH=avx512bw` to
unlock another ~10-15% on AVX-512BW hosts (c7a, c7i, m7i) by widening
auto-vectorization to 512-bit. **It does not.** The avx2 → avx512bw
transition has fundamentally different hardware economics from the
sse41 → avx2 transition.

## The two AVX-512 perf hazards

### 1. AMD Zen 4 µop-split (c7a)

AMD's Zen 4 cores (c7a / Genoa, c8a / Bergamo) implement 512-bit AVX-512
operations by issuing **2× 256-bit µops per 512-bit op**. For
auto-vectorized loops:

- Iteration latency doubles.
- Iteration count only halves if the trip count is large enough.
  Short-trip loops eat the 2× latency without amortizing.
- 512-bit instruction encodings are larger → more I-cache pressure.

Net: loops that auto-vectorized productively at 256-bit AVX2 lose
performance when the compiler widens them to 512-bit.

### 2. Intel Sapphire Rapids transition + downclock (c7i / m7i)

Intel's Sapphire Rapids has native 512-bit execution units, so the
µop-split issue does not apply. But it pays:

- ~3-5% AVX-512 frequency downclock under sustained heavy 512-bit use.
- AVX-512 ↔ AVX2 transition penalties when non-kernel TUs running
  512-bit code call into the 256-bit hand-tuned kernel TUs (which
  always run at host tier via the dispatcher).

Net: small or zero gain from widening, often offset by the transition
costs.

## Mitigation: `-mprefer-vector-width=256`

The canonical mitigation (used by FFmpeg, libvpx, Intel ISPC) is to
**keep AVX-512BW *capabilities* available but cap auto-vectorization at
256-bit width**. The flag `-mprefer-vector-width=256` (gcc / clang) /
`-qopt-zmm-usage=low` (icpc) does exactly that:

- The compiler can still emit AVX-512BW instructions where it
  explicitly needs them (mask registers, byte/word lane permutes,
  gather/scatter, the 32-zmm register file).
- The auto-vectorizer's preferred SIMD width stays at 256-bit, dodging
  Zen 4's µop-split and Intel's downclock/transition costs.

The Makefile bakes this flag into `arch=avx512bw` directly. Hand-tuned
512-bit kernel intrinsics in `KERNEL_SRCS` are unaffected — the cap is
about auto-vec, not intrinsics.

## Empirical numbers

c7a.4xlarge (AMD EPYC 9R14, Zen 4) and c7i.4xlarge
(Intel Xeon Platinum 8488C, Sapphire Rapids) running the bench's
wgs-5M sample (1kg HG00096, 5M PE reads on hg38), shm-warmed via
`bwa-mem3 shm`, 3 reps median, timing via [tricord]
(`fg-labs/tricord`):

| host | avx2 | avx512bw | avx512bw + pvw256 (default) |
|------|------|----------|-----------------------------|
| c7a (Zen 4) | 105.70 s | 103.40 s (−2.2%) | **101.03 s (−4.4%)** |
| c7i (Sapphire Rapids) | 156.50 s | 155.47 s (−0.7%) | 155.41 s (−0.7%) |

[tricord]: https://github.com/fg-labs/tricord

The gain is real but small. Defaulting `BASELINE_ARCH=avx2` for x86
distribution is still correct: it's portable across every x86 host and
loses only ~2-4% to the host-locked avx512bw build on AVX-512 hosts.

## When to use `BASELINE_ARCH=avx512bw`

- **Production fleets pinned to AVX-512BW hosts (c7a / c7i / m7i):** ship
  a host-locked build for the small (~2-4%) extra gain. The Makefile's
  `arch=avx512bw` includes the `-mprefer-vector-width=256` cap by
  default, which is empirically the right choice for both Zen 4 and
  Sapphire Rapids. The binary will SIGILL on hosts below avx2; pair
  with explicit Batch queue / image plumbing.
- **Mixed fleets / generic x86 distribution:** stay on the default
  `BASELINE_ARCH=avx2`. The 2-4% gap is small enough that portability
  is worth it.

## Benchmarking it yourself

Build the `avx2` and `avx512bw` variants (the latter already includes
`-mprefer-vector-width=256` by default) and A/B them on your target host with
[bwa-mem3-bench](https://github.com/fg-labs/bwa-mem3-bench); GCUPS is
within-host only, so compare on the hardware you will deploy on. (An
earlier bench-infrastructure regression here was traced to a toolchain
difference, not avx512bw itself, and is resolved.)

## See also

- [SIMD dispatch architecture](../developer-guide/simd-dispatch.md)
  — how the runtime kernel dispatcher picks a tier
- [Build & infrastructure](build-infra.md) — the broader build layout
- [Architecture support](arch-support.md) — per-host SIMD coverage
