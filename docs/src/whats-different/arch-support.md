# Architecture Support

This page covers the architecture-specific build and runtime work carried in
bwa-mem3. The goal is a single codebase that builds cleanly on all supported
targets and runs the best available SIMD kernels on each.

For the full dispatch matrix and runtime selection logic, see
[Performance → SIMD dispatch matrix](../performance/simd-dispatch.md) and
[Developer Guide → SIMD dispatch architecture](../developer-guide/simd-dispatch.md).

## Linux ARM64 / aarch64 build (PR #1)

The Apple Silicon work that reached the fork in commit `ae73227` gated ARM
behavior on `$(UNAME_M) == arm64`. On macOS, `uname -m` returns `arm64`. On
Linux ARM64, it returns `aarch64`. The Makefile's `ifeq` check therefore fell
through to the x86 `multi` target on every Linux aarch64 host, failing with:

```text
g++: error: unrecognized command-line option '-msse'
```

PR #1 introduces an `IS_ARM` variable (`$(filter $(UNAME_M),arm64 aarch64)`)
that matches both names. All four architecture-conditional blocks in the
Makefile are rewritten to use `IS_ARM`: the NEON/sse2neon flag block, the x86
arch-specific block, the ARM64 single-binary build block, and the `multi`
target ARM64 short-circuit. The CI workflow is extended to trigger on pushes
to `fg-main` (the integration branch) and adds an `ubuntu-24.04-arm` matrix
row so the aarch64 path is exercised on every PR.

## `arch=avx512bw` explicit build target (PR #16)

The AVX-512 Smith-Waterman kernels in bwa-mem2 are guarded by the
`__AVX512BW__` preprocessor macro — not `__AVX512F__`. The only way to build
them before this PR was `arch=avx512`, but `make multi` emitted the dispatch
binary as `bwa-mem2.avx512bw`. The build selector (`avx512`), the preprocessor
guard (`__AVX512BW__`), and the dispatcher suffix (`.avx512bw`) disagreed.

PR #16 adds `arch=avx512bw` as an explicit Makefile target with flags
`-mavx512f -mavx512bw` and switches `make multi` to invoke `arch=avx512bw`
when emitting `bwa-mem3.avx512bw`. The legacy `arch=avx512` is preserved as
an alias with identical flags. No C++ is changed; the fix is 11 insertions and
2 deletions in the Makefile.

This is a pure build-correctness fix: before PR #16, `arch=avx512bw` and
`arch=multi` builds on AVX-512BW hardware silently compiled the wrong kernel
(see [Correctness → AVX-512BW dispatch guard](correctness.md) for the
downstream effect).

## NEON kswv mate-rescue (PR #18)

bwa-mem2 has a batched mate-rescue Smith-Waterman path (`BWAMEM_BATCHED_MATESW`)
that uses SIMD kswv kernels to score rescue candidates in parallel. On ARM64
the gate was `__AVX512BW__`, which is never true on NEON hardware. The NEON
`kswv::getScores8` kernel existed in the source but was unreachable in
production.

PR #18 enables this path on ARM64 by replacing the `__AVX512BW__` gate with a
new `BWAMEM_BATCHED_MATESW` macro that fires on NEON/Apple Silicon as well.
Along the way, four kernel bugs were found and fixed:

1. **te split** — the `te` (traceback end) value needed separate hi/lo tracking
   for 16-lane u8 batches.
2. **Freeze mask** — a `frozen_vec` mask now gates `gmax`/`te`/`qe` updates
   after `KSW_XSTOP` fires, preventing stale values from escaping to the score2
   scan.
3. **Per-lane score2 exclusion** — `len1`, `low`/`high`, and `qe` masks were
   not applied per-lane in Loop 1, allowing lanes without a valid primary to
   contribute spurious suboptimal scores.
4. **minsc filter on rowMax** — sub-`minsc` plateau scores were leaking into
   `score2` because the scalar `ksw_u8` gating condition (`imax >= minsc`) was
   not replicated.

Measured on an M-series Mac (8 threads, 500k PE 100 bp reads on chr17):
**1.42× speedup** (−29.4% wall time) with byte-identical sorted SAM output.

## AVX2 kswv mate-rescue (PR #20)

PR #18 enabled batched mate-rescue on ARM64. Most x86 production deployments
(AWS c6a, c6i, older Xeons) use AVX2 without AVX-512BW and were excluded from
the same gate. PR #20 extends the batched path to AVX2 by adding a 256-bit
`kswv256_u8` kernel and widening `BWAMEM_BATCHED_MATESW` to fire on `__AVX2__`.

The AVX2 kernel is a direct port of the corrected NEON kernel from PR #18,
with an additional fix for per-lane `te2` tracking (`_mm256_blendv_epi8` on a
sign-extended 8→16 bit mask). Verified byte-identical sorted SAM vs the
pre-`BWAMEM_BATCHED_MATESW` scalar control on EC2 m5.xlarge (Skylake-SP, 4
threads, 500k chr17 PE pairs).

Note: PR #20 introduced a score2 plateau regression in the AVX2 kernel that
was identified and fixed in the [correctness series](correctness.md)
(PRs #27, #28, #29).

---

## Changes catalog

| Item | bwa-mem3 PR | Upstream PR/issue | Status |
|------|-------------|-------------------|--------|
| Linux ARM64 / aarch64 build + CI | [#1](https://github.com/fg-labs/bwa-mem3/pull/1) | [bwa-mem2#288](https://github.com/bwa-mem2/bwa-mem2/pull/288) | fork-only (upstream PR open) |
| `arch=avx512bw` explicit target | [#16](https://github.com/fg-labs/bwa-mem3/pull/16) | — | fork-only |
| NEON kswv mate-rescue kernel | [#18](https://github.com/fg-labs/bwa-mem3/pull/18) | — | fork-only |
| AVX2 kswv mate-rescue kernel | [#20](https://github.com/fg-labs/bwa-mem3/pull/20) | — | fork-only |

---

**See also:**
[Performance → SIMD dispatch matrix](../performance/simd-dispatch.md) ·
[Developer Guide → SIMD dispatch architecture](../developer-guide/simd-dispatch.md) ·
[Developer Guide → Apple Silicon / NEON port](../developer-guide/neon-port.md) ·
[Correctness fixes](correctness.md) ·
[Performance → PGO build](../performance/pgo.md)
