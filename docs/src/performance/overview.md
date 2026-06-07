# Performance Overview

Performance claims in this section are benchmarked, not asserted. The canonical source of truth for benchmark methodology, hardware configurations, and current numbers is [bwa-mem3-bench](https://github.com/fg-labs/bwa-mem3-bench), a reproducible benchmarking harness that runs across AWS Batch architectures (x86 AVX2, AVX-512, ARM Graviton). Consult that repository before drawing conclusions from isolated anecdotal timings.

## What drives bwa-mem3's performance

bwa-mem3 inherits the SIMD-vectorized alignment kernels of bwa-mem2 and adds several improvements of its own. The headline gains relative to a stock bwa-mem2 build fall into four categories.

**Vectorized alignment kernels.** The Smith-Waterman and banded-SWA kernels (kswv, bandedSWA) are compiled against the widest SIMD ISA the current CPU supports — SSE4.1 through AVX-512BW on x86, or native NEON on ARM. On Apple Silicon, native NEON intrinsics replaced the sse2neon shim in the two hottest kernels, delivering roughly 10% additional throughput over the pure-translation baseline. See [SIMD dispatch matrix](simd-dispatch.md) for the full picture.

**libsais FM-index construction.** The indexing step uses the linear-time suffix-array/BWT construction library libsais in place of the original quadratic-time approach. This cuts `bwa-mem3 index` wall time substantially on large references. See [What's Different — Performance improvements](../whats-different/performance.md) for the corresponding PR details.

**mimalloc allocator.** bwa-mem3 vendors and statically links [mimalloc](https://github.com/microsoft/mimalloc), replacing the system `malloc`/`free` for all allocations. On Linux the library is injected via `--whole-archive`; on macOS it uses dyld interposition. The allocator shows consistent throughput gains on multi-threaded workloads because mimalloc avoids the lock contention in glibc's `ptmalloc` at high thread counts. See [User Guide — Memory allocator](../user-guide/allocator.md) for details.

**Profile-Guided Optimization (PGO).** The build system provides `make pgo-generate` and `make pgo-use` targets that compile an instrumented binary, gather branch-probability and call-frequency profiles from a representative workload, and then recompile with those profiles applied. On Apple Silicon the measured gain is approximately 3%; on x86 the gain depends on the workload mix. PGO is opt-in and is not applied to the default `make` output. See [PGO build](pgo.md) for the full workflow.

## Consolidated mapping speedups

PR [#58](https://github.com/fg-labs/bwa-mem3/pull/58) and the related lockstep SMEM-batching work ([#33](https://github.com/fg-labs/bwa-mem3/pull/33)) reduced per-read overhead in the main mapping loop beyond what upstream bwa-mem2 carries. The batch `-H` ingestion improvement ([#49](https://github.com/fg-labs/bwa-mem3/pull/49)) further reduces header-processing latency for large sample sets.

## Reference numbers across architectures

Wall-time medians from [bwa-mem3-bench](https://github.com/fg-labs/bwa-mem3-bench) at SHA `bffae5a` (2026-06-07), 5 reps per cell, t≈16, hg38, paired-end 150 bp:

| sample | c6a (AVX2, Zen3) | c7a (AVX-512, Zen4) | c7i (AVX-512, SPR) | c7g (NEON, Graviton3) | c8g (NEON, Graviton4) |
|---|---:|---:|---:|---:|---:|
| wgs-5M | 142.88 s | **97.94 s** | 147.83 s | 176.49 s | 150.76 s |
| wes-5M | 78.87 s | **59.81 s** | 79.25 s | 81.81 s | 67.47 s |
| panel-twist-5M | 157.12 s | **105.99 s** | 166.15 s | 192.53 s | 164.16 s |

Concordance vs upstream `bwa-mem2 v2.2.1` on these cells, measured over primary-alignment records: **wgs-5M 99.9893%, wes-5M 99.9996%, panel-twist-5M 99.9414%**. bwa-mem3 is intentionally **not** byte-identical to bwa-mem2 — the residual differences are additive SAM tags, per-architecture SIMD `score2`/`MAPQ` convergence, deterministic tie-breaks, and a small number of additional supplementary alignments; see [Equivalence with bwa-mem2](../whats-different/equivalence.md) for the full audited breakdown. NEON-vs-x86 cross-architecture concordance on the same builds remains **100.0000%** (the ARM and x86 fg-labs builds produce identical records). Spot-pool noise envelope (rep-to-rep CV): ~1% on c6a / c7a / c7g / c8g, ~8–9% on c7i. See the bench repo for the methodology, the full per-rep table, and noisier instance classes excluded from this summary.

## Benchmarking responsibly

Alignment throughput is sensitive to read length, error rate, reference size, thread count, CPU architecture, NUMA topology, and whether the index is cold (in-kernel page cache) or warm. The [bwa-mem3-bench](https://github.com/fg-labs/bwa-mem3-bench) harness controls for these variables by running standardized workloads on defined instance types. If you need numbers for a procurement or publication decision, run the harness against your target hardware.

---

**See also:**
[SIMD dispatch matrix](simd-dispatch.md) ·
[PGO build](pgo.md) ·
[Tuning checklist](tuning.md) ·
[What's Different — Performance improvements](../whats-different/performance.md) ·
[bwa-mem3-bench](../related-projects/bwa-mem3-bench.md)
