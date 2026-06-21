# Performance Overview

Performance claims in this section are benchmarked, not asserted. The canonical source of truth for benchmark methodology, hardware configurations, and current numbers is [bwa-mem3-bench](https://github.com/fg-labs/bwa-mem3-bench), a reproducible benchmarking harness that runs across AWS Batch architectures (x86 AVX2, AVX-512, ARM Graviton). Consult that repository before drawing conclusions from isolated anecdotal timings.

## What drives bwa-mem3's performance

There is no single "the speedup." bwa-mem3 inherits bwa-mem2's SIMD-vectorized
core and layers on a series of independent improvements, each targeting a
different bottleneck in the align pipeline. How much any one of them helps —
and therefore the total — depends heavily on the workload: read length, error
rate, reference size, thread count, CPU architecture (AVX2 vs AVX-512 vs NEON),
and whether the index is cold or warm in the page cache. A short-read whole-genome
run is dominated by seeding and FM-index walks; a long-read or high-error run
spends most of its cycles in the Smith–Waterman kernels; a many-sample run can be
bottlenecked on header ingestion or decompression. The drivers below group by
*what part of the machine or algorithm they fix*. For real, reproducible numbers
on specific hardware, always defer to [bwa-mem3-bench](https://github.com/fg-labs/bwa-mem3-bench)
rather than any single anecdote here.

For the full per-change list with PR links and status, see
[What's Different — Performance improvements](../whats-different/performance.md).

### 1. Getting more out of the machine (SIMD / microarchitecture)

The alignment kernels (`kswv` mate-rescue, `bandedSWA`) are compiled for the
widest SIMD ISA the host supports — SSE4.1 through AVX-512BW on x86, native NEON
on ARM — and selected at runtime by an in-process dispatcher (a single binary
that picks the right kernel per host, [#83](https://github.com/fg-labs/bwa-mem3/pull/83))
rather than the older multi-binary `execv` launcher.

- **Native NEON kernels** replaced the `sse2neon` translation shim in the two
  hottest kernels on ARM, worth roughly **10% additional throughput** on Apple
  Silicon over the pure-translation baseline.
- **AVX2 as the x86 baseline** ([#84](https://github.com/fg-labs/bwa-mem3/pull/84)):
  restoring the non-kernel translation units to `-mavx2` recovered ~**+15%** user
  time on `wgs-5M` / ~**+11%** on `wes-5M` that had been lost when the baseline
  briefly dropped to SSE4.1 — hot non-kernel paths (chain extension, FM-index BWT
  walks, mate scoring) auto-vectorize at 256-bit width again.
- **Capping AVX-512BW auto-vectorization at 256-bit** ([#86](https://github.com/fg-labs/bwa-mem3/pull/86))
  avoids the frequency downclock that wide 512-bit auto-vec code triggers on some
  x86 parts, where the wider vectors cost more in clock than they return.
- **Per-strip L1 prefetches** added to the 8-/16-bit `kswv` kernels that lacked
  them ([#70](https://github.com/fg-labs/bwa-mem3/pull/70), bringing them in line
  with `kswv512_16`) stop the inner SW loop from stalling on first-touch L1
  misses.
- **A recovered 8-bit banded SW path** for reads ≥128 bp
  ([#140](https://github.com/fg-labs/bwa-mem3/pull/140) and follow-ups) keeps
  long-read alignment in the cheaper 8-bit lane width where it is valid.
- **Per-ISA SW kernel hand-tuning.** A NEON pass
  ([#160](https://github.com/fg-labs/bwa-mem3/pull/160)) replaced the multi-instruction
  `sse2neon` "all lanes zero" tests in the hot band-narrowing scans with a single
  `vmaxvq` horizontal reduction; an AVX2 pass
  ([#161](https://github.com/fg-labs/bwa-mem3/pull/161)) relieved the port-5
  `vpblendvb`/`vpshufb` bottleneck that ceilings the inner loop on Zen3. Both are
  byte-identical to the prior output.
- **AVX2 16-bit mate rescue** ([#162](https://github.com/fg-labs/bwa-mem3/pull/162)):
  AVX2-only hosts (e.g. Zen3, which run the AVX2 default) previously fell back to
  scalar `ksw_align2` for 16-bit mate rescue because only NEON and AVX-512 had a
  batched 16-bit `kswv`; `kswv256_16` closes that gap.

See the [SIMD dispatch matrix](simd-dispatch.md) for the full ISA picture.

### 2. Fixing bad patterns in the hot path (algorithmic rewrites)

Several gains come simply from rewriting code that did more work than it needed
to — the "rewrites wind up fixing bad patterns" effect — without changing where
reads map.

- **Lockstep SMEM batching** ([#33](https://github.com/fg-labs/bwa-mem3/pull/33),
  widened from 8 to 16 reads in [#75](https://github.com/fg-labs/bwa-mem3/pull/75)):
  advances several reads' seed walks in interleaved order so the out-of-order
  engine overlaps the random `cp_occ` checkpoint-array cache misses of one read
  with the compute of another. Measured ~**−6% wall time** on 150 bp WGS, with the
  hot `cp_occ` load share dropping from 65.5% to 53.3% of seeding time.
- **O(n²) → O(n) header ingestion** ([#49](https://github.com/fg-labs/bwa-mem3/pull/49)):
  the `-H` path re-ran `strlen` + `realloc` per line; batching the read took a
  ~70 MB / 1.5 M-line header from **>10 minutes to under a second** before
  alignment even starts.
- **Closed-form ungapped scoring** when there are no mismatches
  ([#77](https://github.com/fg-labs/bwa-mem3/pull/77)) replaces a per-base
  Kadane-style walk with a direct store.
- **On-stack sort buffers for small arrays** ([#78](https://github.com/fg-labs/bwa-mem3/pull/78))
  removes a per-read `malloc`/`free` that was dominating the sort of the typical
  5–30-element alignment-region arrays.
- **Inlining `backwardExt`** ([#88](https://github.com/fg-labs/bwa-mem3/pull/88))
  eliminates a struct-by-value ABI pass that gcc 12+ could not optimize away,
  recovering (and beating) the older compiler's baseline.
- **`pdqsort` with stable tie-breaks** at the dedup-patch sort sites
  ([#123](https://github.com/fg-labs/bwa-mem3/pull/123)).
- The **consolidated mapping-speedup audit** ([#58](https://github.com/fg-labs/bwa-mem3/pull/58))
  bundles tuning across the ksw2 band loop, SMEM batching, suffix-array lookup
  prefetch, and SAM record building — historically the single largest wall-time
  step in `main`.

### 3. Indexing

`bwa-mem3 index` builds the FM-index with the linear-time
[libsais](https://github.com/IlyaGrebnov/libsais) library
([#57](https://github.com/fg-labs/bwa-mem3/pull/57)) instead of the older
sais-lite, cutting both wall time and peak memory while producing a
byte-identical index (existing indexes need no rebuild). Construction also skips
the wasted zero-initialization of unpack and suffix-array buffers
([#80](https://github.com/fg-labs/bwa-mem3/pull/80)), which on a doubled-human
input avoided tens of GiB of write-then-overwrite zero-fill, and right-sizes the
SA-entry staging buffers to the actual write count rather than the uncapped
SA-interval sum ([#157](https://github.com/fg-labs/bwa-mem3/pull/157)).

### 4. Memory allocation and I/O

- **mimalloc**, vendored and statically linked by default
  ([#19](https://github.com/fg-labs/bwa-mem3/pull/19)), replaces the system
  allocator and avoids glibc `ptmalloc`'s lock contention at high thread counts —
  a consistent multi-threaded throughput win. See
  [User Guide — Memory allocator](../user-guide/allocator.md).
- **Faster read ingestion**: a content-detecting FASTQ fast path over libdeflate
  BGZF ([#128](https://github.com/fg-labs/bwa-mem3/pull/128), merged) cuts the
  cost of decompressing and parsing input, which matters most when the aligner
  would otherwise be I/O- or parse-bound. A vendored zlib-ng inflate path with a
  third pipeline worker ([#153](https://github.com/fg-labs/bwa-mem3/pull/153),
  merged) extends the same idea: it cuts the read stage ~**2.2×** on gzip input
  and, because the serial read stage becomes a larger share of the wall as the
  compute stage parallelizes, improves end-to-end wall by up to **−7.8%** at 96
  cores. The accompanying `--profile` stage-timing mode
  ([#152](https://github.com/fg-labs/bwa-mem3/pull/152)) is what attributes wall
  time across the read‖proc‖write stages to surface exactly this effect.

### 5. Build-time optimization (PGO)

The build provides opt-in `make pgo-generate` / `make pgo-use` targets that
recompile with branch-probability and call-frequency profiles gathered from a
representative workload — ~**3%** on Apple Silicon, workload-dependent on x86.
PGO is not applied to the default `make` output. See [PGO build](pgo.md).

## Reference numbers across architectures

Wall-time medians from [bwa-mem3-bench](https://github.com/fg-labs/bwa-mem3-bench) at SHA `a02fcb4` (2026-06-20), 5 reps per cell, t≈16, hg38, paired-end 150 bp:

| sample | c6a (x86-64, AVX2, Zen3) | c7a (x86-64, AVX-512, Zen4) | c7i (x86-64, AVX-512, SPR) | c7g (arm64, NEON, Graviton3) | c8g (arm64, NEON, Graviton4) |
|---|---:|---:|---:|---:|---:|
| wgs-5M | 131.63 s | **92.66 s** | 112.98 s | 154.27 s | 127.93 s |
| wes-5M | 76.57 s | **56.53 s** | 67.66 s | 79.67 s | 65.26 s |
| panel-twist-5M | 150.84 s | **102.72 s** | 156.42 s | 181.37 s | 148.50 s |

Concordance vs upstream `bwa-mem2 v2.2.1` on these cells, measured over primary-alignment records: **wgs-5M 99.9893%, wes-5M 99.9996%, panel-twist-5M 99.9414%**. bwa-mem3 is intentionally **not** byte-identical to bwa-mem2 — the residual differences are additive SAM tags, per-architecture SIMD `score2`/`MAPQ` convergence, deterministic tie-breaks, and a small number of additional supplementary alignments; see [Equivalence with bwa-mem2](../whats-different/equivalence.md) for the full audited breakdown. NEON-vs-x86 cross-architecture concordance on the same builds remains **100.0000%** (the ARM and x86 fg-labs builds produce identical records). Spot-pool noise envelope (rep-to-rep CV) for this run: ~1–2% on c7a / c7g / c8g, ~1–5% on c6a, ~10–16% on c7i — the c7i medians in particular carry wide error bars and should be read as directional. See the bench repo for the methodology, the full per-rep table, and noisier instance classes (e.g. m7i) excluded from this summary.

Release-to-release speedups are deliberately uneven across this grid. A workload's gain scales with the share of its wall time spent in the Smith–Waterman kernels — highest on `wgs-5M` (~85% of cycles), lowest on the seed/IO-bound `wes-5M`, and intermediate on `panel-twist-5M`, whose deep target coverage produces many split alignments that each re-enter SW. That per-workload factor is multiplied by how heavily a given release retuned the host's per-ISA SW kernel: v0.3.0 concentrated on the NEON ([#160](https://github.com/fg-labs/bwa-mem3/pull/160), [#166](https://github.com/fg-labs/bwa-mem3/pull/166)) and AVX2 ([#161](https://github.com/fg-labs/bwa-mem3/pull/161), [#162](https://github.com/fg-labs/bwa-mem3/pull/162)) kernels, so Graviton and AVX2 hosts moved more than the already-tuned AVX-512BW path. See [SIMD dispatch](simd-dispatch.md) for how the per-host kernel is selected.

## Benchmarking responsibly

Alignment throughput is sensitive to read length, error rate, reference size, thread count, CPU architecture, NUMA topology, and whether the index is cold (in-kernel page cache) or warm. The [bwa-mem3-bench](https://github.com/fg-labs/bwa-mem3-bench) harness controls for these variables by running standardized workloads on defined instance types. If you need numbers for a procurement or publication decision, run the harness against your target hardware.

---

**See also:**
[SIMD dispatch matrix](simd-dispatch.md) ·
[PGO build](pgo.md) ·
[Tuning checklist](tuning.md) ·
[What's Different — Performance improvements](../whats-different/performance.md) ·
[bwa-mem3-bench](../related-projects/bwa-mem3-bench.md)
