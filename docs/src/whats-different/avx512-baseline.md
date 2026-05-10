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

## Why the runtime warning was misleading

Earlier versions of `src/simd_dispatch.cpp` printed at startup:

```text
[W::bwamem3_simd_init_body] build baseline avx2 < host tier avx512bw;
non-kernel TUs are not auto-vectorized at the higher width (expect
10-15% slower hot paths). Rebuild with BASELINE_ARCH=avx512bw to recover.
```

The "10-15%" figure was the **sse41 → avx2** transition on AVX2-only
hosts (PR #84's measurement, before avx2 became the default). It did
not generalize to **avx2 → avx512bw** for the µop-split / downclock /
transition reasons above. The warning was demoted to BWAMEM3_DEBUG_SIMD
gating in a follow-up commit; the recommendation reflects the actual
measurements (typically <2% wall-time gain on AVX-512 hosts).

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

## Reproducer

The investigation harness lives at
`scripts/perf-diff-baseline-arch.sh`. It builds N variants of bwa-mem3
with different `BASELINE_ARCH` and `EXTRA_CXXFLAGS` settings, runs each
through `tricord` (or `perf record` for hot-function diffs), and emits
per-variant median tables. Example usage:

```bash
scripts/perf-diff-baseline-arch.sh \
    --ref hg38.fasta --r1 r1.fq.gz --r2 r2.fq.gz \
    --out out/ --reps 3 --threads 16 \
    --variants 'avx2:,avx512bw:,avx512bw-pvw256:-mprefer-vector-width=256'
```

Requires `tricord` on the PATH (`cargo install tricord`).
`/dev/shm` must be ≥18 GB to stage the hg38 FMI index — on a default
EC2 instance that means `mount -o remount,size=28g /dev/shm`.

## Bench-side caveats

The bench's Phase C report (May 2026) reported a +14% c7a regression
when comparing the bench's portable image vs the host-locked
`avx512bw` image inside AWS Batch. That delta does not reproduce on a
single-instance bare-metal measurement. A 4-way disambiguation test
(c7a.4xlarge, wgs-5M, shm-warmed, 3 reps each — see "Reproducer" above)
attributes the gap to two independent bench-side factors:

| variant | wall (s) | vs A |
|---|---|---|
| A: AL2023-built avx512bw, bare-metal | 99.50 | (baseline) |
| B: AL2023-built avx512bw, in bench Docker image | 102.30 | +2.8% |
| C: bench-image avx512bw binary, bare-metal | 117.03 | +17.6% |
| D: bench-image avx512bw binary, in bench Docker image | 118.44 | +19.0% |

The findings:

1. **Build-environment matters more than container.** The bench's
   `:316dba6-avx512bw` image binary, run **bare-metal** on the same
   c7a.4xlarge with the same input, is +17.6% slower than a fresh
   AL2023-built binary from the same SHA with the same `BASELINE_ARCH`.
   - **AL2023** ships `gcc 11.5.0` and defaults to `-no-pie`. Output is
     a non-PIE ELF.
   - **Debian bookworm** (the bench's Dockerfile base) ships `gcc 12.x`
     and defaults to `-pie`. Output is a PIE ELF.
   - PIE adds indirection through a GOT for every global reference and
     is well-known to cost 5–15% on tight CPU loops. Combined with
     gcc-11 vs gcc-12 codegen differences this comfortably accounts
     for the +17.6%.
2. **Docker container overhead is small (~3%).** A→B and C→D both show
   ~2–3% wall-time delta when wrapping the same binary in the bench's
   image. Consistent with the broader literature on cgroup-namespaced
   compute-bound workloads.
3. **The bench's portable `:316dba6` image is built with
   `BASELINE_ARCH=sse41`, not `avx2`.** Direct evidence from the
   binary's startup banner inside the container:
   `[W::bwamem3_simd_init_body] build baseline sse41 < host tier
   avx512bw`. The banner is generated from compile-time macros in
   `simd_dispatch.cpp` and unambiguously testifies to the build flags
   used. Why it's sse41 (rather than the post-PR-#84 default of avx2)
   is bench-side mystery — possibilities include an explicit
   `BASELINE_ARCH=sse41` build-arg, a Makefile-default change between
   the prior portable build and SHA 316dba6, or an environment quirk
   in the Docker build. The Phase C report's "42/42 saw avx2 warning"
   summary likely reflects a different prior run, not the current
   `:316dba6` portable image.

The Phase C report's headline "+14% c7a regression for avx512bw" is
therefore comparing an **sse41-built portable image at +17.6%
build-environment penalty** against a **`BASELINE_ARCH=avx512bw`
binary at the same +17.6% penalty plus Zen 4's µop-split cost** —
which roughly cancels for a small absolute delta in either direction.
None of it is bwa-mem3's BASELINE_ARCH knob's fault.

These are bench-side concerns. The bwa-mem3 fix
(`-mprefer-vector-width=256` for `arch=avx512bw`) stands on its own
bare-metal merit: −2.2% on Zen 4 vs avx2 vanilla, plus −2.2%
incremental from the cap; wash on Sapphire Rapids.

### Bench-side toolchain attribution

The +17.6% bare-metal delta between the bench's bookworm-built binary
and an AL2023-built binary was decomposed via a six-variant single-host
test (c7a.4xlarge, wgs-5M, shm-warmed, 5 reps each, tricord median):

| variant | gcc | PIE | CET (`endbr64`) | wall (s) | vs AL2023 |
|---|---|---|---|---|---|
| AL2023 default                                  | 11.5.0 | no  | 8 (libc init only)  | **98.28**  | (baseline) |
| bookworm + gcc-11                               | 11.4.0 | yes | 3 (libc init only)  | **100.16** | +1.9% |
| bookworm default                                | 12.2.0 | yes | 3 (libc init only)  | **117.21** | +19.3% |
| bookworm + `-no-pie`                            | 12.2.0 | no  | 3                   | **118.98** | +21.1% |
| bookworm + `-fcf-protection=none`               | 12.2.0 | yes | 3                   | **118.40** | +20.5% |
| bookworm + `-no-pie -fcf-protection=none`       | 12.2.0 | no  | 3                   | **116.84** | +18.9% |

Findings:

- **The +17.6% delta is gcc-11 vs gcc-12 codegen, not any hardening
  flag.** Switching gcc inside Debian bookworm (keeping PIE on, keeping
  whatever default CET there is, keeping every other Debian default)
  recovers the perf to within 2% of AL2023.
- **Disabling PIE in bookworm gcc-12 has no measurable effect.** Median
  118.98s with `-no-pie` vs 117.21s default — within rep-to-rep noise.
  The 5–15% PIE penalty cited in the literature doesn't manifest for
  bwa-mem3's intrinsic-heavy hot path; GOT indirection is rare on
  tight inner loops.
- **Disabling CET in bookworm gcc-12 has no measurable effect.** Both
  bookworm-default and `-fcf-protection=none` builds emit only ~3
  `endbr64` instructions in 7.6 MB of binary (libc init code).
  Something — probably bwa-mem3's `-mavx512f -mavx512bw` per-tier
  kernel flags, or a `#pragma GCC` somewhere — suppresses CET emission
  regardless of the flag. So disabling it removes nothing.
- **The combined `-no-pie -fcf-protection=none` build is statistically
  indistinguishable from the default.** Both PIE and CET are noise.

So the actionable bench-side fix isn't `-no-pie` or
`-fcf-protection=none` — it's "use gcc-11". The six-variant table
above is the only data we have; gcc-13 was not tested here, and the
postscript below shows that gcc-13 / gcc-14 do *not* recover the gap
without #88's source-side fix. A one-liner for the bench Dockerfile:

```dockerfile
RUN apt-get install -y gcc-11 g++-11
RUN cd fg-labs && CC=gcc-11 CXX=g++-11 make BASELINE_ARCH=... -j
```

…recovers ~17% wall-time **uniformly across every Batch worker, every
arch**. That's a much bigger lever than the bwa-mem3-side
`-mprefer-vector-width=256` mitigation here, which is ~2–4% on c7a and
wash on c7i. But it's bench infra, not bwa-mem3 source.

(The deeper question — *why* is gcc-12 codegen ~19% slower than gcc-11
on Zen 4 for this workload? — was followed up in [#88]; see the
postscript below.)

The bwa-mem3-side conclusion stands on its own bare-metal merit:
`-mprefer-vector-width=256` for `arch=avx512bw` is a real ~2–4% win on
c7a and a wash on c7i, independent of toolchain and container
concerns.

### Postscript: gcc-12 attribution closed by #88

The "use gcc-11" recommendation above was the actionable bench-side
fix at the time this page was written. [#88] (`perf(fmi): inline
backwardExt to recover gcc 12+ wall-clock regression`) has since
identified the underlying mechanism and closed the gap in source — so
on current main no compiler pin is required.

Profile attribution on a fresh c7a.8xlarge run with `perf record
--no-children` localized ~9 percentage points of wall-time to
`FMI_search::backwardExt`'s self-time (12.5% on gcc-11 vs 21.5% on
gcc-14). Disassembly histograms were nearly identical between
compilers (~110 instructions, 8 scalar `popcntq`, 25 mov each); IPC
fell from 1.77 to 1.60. `perf annotate` isolated a single instruction
at 42% of the function's samples on gcc-14: `vmovdqu %ymm0, (%r8)` —
the 32-byte AVX store of the `SMEM` return value through SysV's
hidden-pointer convention. The matching argument-load
(`mov 0x30(%rbp), %r10` for `smem.s` from the caller's stack push) was
next-hottest at 17%. Together those two instructions accounted for
~60% of the function's self-time.

The fix in #88 is a one-line attribute change: marking
`backwardExt` `__attribute__((always_inline))` removes the call
boundary at all 9 hot call sites in `getSMEMs*` and
`ls_advance_*`. Without a call boundary the `SMEM` struct stays in
caller registers across the would-be call site — no struct push, no
return-slot store, no `vzeroupper`.

Post-#88 numbers on c7a.8xlarge (hg38 wgs-5M, shm-warmed, `-t 32`, 5
reps mean, single-binary at `BASELINE_ARCH=avx2`):

| binary           | wall (s) | vs gcc-14 | vs gcc-11 | IPC  |
|------------------|----------|-----------|-----------|------|
| main + gcc-11    | 64.59    | −6.6%     | baseline  | 1.77 |
| main + gcc-14    | 69.16    | baseline  | +7.07%    | 1.60 |
| **#88 + gcc-14** | **61.94**| **−10.4%**| **−4.10%**| **1.83** |

So the gcc-11 vs gcc-12 attribution above was the surface symptom of
an ABI-level inefficiency that was costing cycles on gcc-11 too — the
always-inline fix beats gcc-11 baseline by 4.10%, not just gcc-14. The
empirical table and findings list in the previous section remain
accurate as an investigation snapshot at commit `316dba6` (the pre-#88
state); the bench-side `gcc-11` recommendation it produced is
obsolete.

[#88]: https://github.com/fg-labs/bwa-mem3/pull/88

## See also

- [SIMD dispatch architecture](../developer-guide/simd-dispatch.md)
  — how the runtime kernel dispatcher picks a tier
- [Build & infrastructure](build-infra.md) — the broader build layout
- [Architecture support](arch-support.md) — per-host SIMD coverage
