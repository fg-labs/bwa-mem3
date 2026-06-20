# Performance Improvements

This page covers the performance work carried in bwa-mem3 on top of upstream
bwa-mem2. Almost every change listed here is a pure throughput or memory win
that preserves the aligner's output; the one exception is the deterministic
tie-break ordering in [#123](https://github.com/fg-labs/bwa-mem3/pull/123),
which can reorder equal-scoring alignments relative to upstream (see
[Equivalence with bwa-mem2](equivalence.md) for the full, audited list of where
bwa-mem3 output diverges).

For a reader-friendly grouping of *what drives the speedup* — by machine
architecture, hot-path rewrite, indexing, allocation/I/O, and build-time — see
[Performance → Overview](../performance/overview.md). For current benchmark
numbers across architectures and workloads, see
[bwa-mem3-bench](https://github.com/fg-labs/bwa-mem3-bench), the canonical
source of truth for benchmark methodology and results.

## Lockstep SMEM batching (PR #33)

Seeding in bwa-mem2 advances one read's SMEM walk at a time. Because each
forward/backward extension step issues a random access into the `cp_occ`
checkpoint array (~4 GB for human genome), the CPU stalls on cache misses
between steps. Lockstep batching advances `SMEM_LOCKSTEP_N` reads' SMEM walks
in slot-interleaved round-robin order so that the out-of-order engine can
overlap the `cp_occ` cache-miss loads for read _i+N_ with the compute-bound
walk of read _i_.

Each read slot (`BatchSlot`) carries its own `prev[]` walk buffer and
`match_buf[]` reorder buffer. A tight recycling loop assigns finished slots to
the next unprocessed read immediately. The match-emit cursor enforces
input-index order so output is byte-identical to scalar. `SMEM_LOCKSTEP_N` is
compile-time tunable; `N=1` dispatches to the unchanged scalar path for
bisection.

Measured improvement on 150 bp NovaSeq WGS (1M pairs, hg38, Graviton3 r7g.4xlarge,
8 threads): **−6.1% wall time** (82 s → 77 s). The `backwardExt` hot
`cp_occ` load share dropped from 65.5% to 53.3% of function time — direct
evidence that the OoO engine is overlapping cross-slot loads. On 300 bp MiSeq
reads the workload is SW-dominated (~85% of cycles in kswv kernels) and the
SMEM improvement is within noise; parity holds.

Supersedes PR #15 (cross-read `_mm_prefetch` shape), which regressed on
Graviton3.

## Batched `-H` header ingestion (PR #49, closes issue #37)

Passing a large header file via `-H <file>` re-ran `strlen` on the growing
header string and called `realloc` on every input line, making ingestion O(n²)
in the number of header lines. For a ~70 MB / ~1.5 M-line header (reported in
upstream [bwa-mem2#204](https://github.com/bwa-mem2/bwa-mem2/pull/204)) this
caused runtimes exceeding 10 minutes before alignment started.

The fix introduces `bwa_insert_header_file`, a batched helper that determines
the file size with `fseek`/`ftell`, allocates a single buffer, copies all
`@`-prefixed lines in one pass, and calls `bwa_insert_header` once. The fix
also addresses four correctness gaps in the upstream PR #204: the return-value
assignment was dropped (leaving `hdr_line` stale after realloc), `const FILE*`
caused compiler warnings, empty files were not guarded, and each `fgets` was
not bounded by remaining buffer. A regression test
(`test/header_insert_test.cpp`) diffs the batched path against the pre-patch
per-line baseline across eight edge cases.

## libsais FM-index construction (PR #57)

`bwa-mem3 index` now builds the FM-index using
[libsais v2.9.1](https://github.com/IlyaGrebnov/libsais) (Ilya Grebnov)
instead of the sais-lite (Yuta Mori saisxx) library that bwa-mem2 inherited.
libsais is actively maintained, supports OpenMP-parallel induced sorting, and
produces a byte-identical FM-index. No changes are required to existing
indexes — bwa-mem3 reads index files built by `bwa-mem2 index` without
re-indexing.

For a human reference (GRCh38 + decoys), libsais reduces indexing wall time
and peak memory vs sais-lite. Exact numbers depend on thread count and
available RAM; see the PR body for measurements on Graviton3.

## Consolidated mapping speedups (PR #58)

PR #58 is a multi-phase performance audit of bwa-mem2's hot path, squashed and
rebased onto `main`. It incorporates improvements across five subsystems:

- **ksw2 banded SW** — tuned the band extension loop to reduce redundant
  computation in the common case.
- **SMEM lockstep batching** — additional refinements on top of PR #33.
- **SAL prefetch** — prefetch hints for the suffix array lookup hot path.
- **SAM record building** — reduced per-record allocation in the text
  formatting path.
- **PGO build** — the opt-in profile-guided optimization target (see also
  [Performance → PGO build](../performance/pgo.md)) is included in this
  suite.

On the smoke-1M workload (1M PE 150 bp reads, hg38, Graviton3 r7g.4xlarge, 16
threads, warm page cache), this PR contributed the largest single-step wall
time reduction in the `main` branch's performance history. Benchmark details
are maintained at [bwa-mem3-bench](https://github.com/fg-labs/bwa-mem3-bench).

---

## Changes catalog

The mechanism column says, in one line, *why it is faster*. The stage column
groups changes by where in the pipeline they act: **seed** (SMEM/FM-index
walks), **sw** (Smith–Waterman / banded kernels), **index** (`bwa-mem3 index`),
**i/o** (read decompression + header ingestion), **mem** (allocation),
**dispatch** (SIMD/runtime selection), or **sort**.

### Merged

| Item | Stage | Mechanism | bwa-mem3 PR | Upstream |
|------|-------|-----------|-------------|----------|
| mimalloc by default | mem | Vendored + statically linked allocator; avoids glibc `ptmalloc` lock contention at high thread counts | [#19](https://github.com/fg-labs/bwa-mem3/pull/19) | — |
| Lockstep SMEM batching | seed | Interleaves several reads' seed walks so the OoO engine overlaps `cp_occ` cache misses across reads | [#33](https://github.com/fg-labs/bwa-mem3/pull/33) | — |
| Batched `-H` header ingestion | i/o | One-pass buffered read; fixes O(n²) `strlen`+`realloc` per line (>10 min → <1 s on large headers) | [#49](https://github.com/fg-labs/bwa-mem3/pull/49) | [bwa-mem2#204](https://github.com/bwa-mem2/bwa-mem2/pull/204) (open); closes [#37](https://github.com/fg-labs/bwa-mem3/issues/37) |
| libsais FM-index construction | index | Linear-time suffix-array/BWT build; less wall time + memory; byte-identical index | [#57](https://github.com/fg-labs/bwa-mem3/pull/57) | — |
| Consolidated mapping speedups | seed/sw | Multi-phase hot-path audit across ksw2 band loop, SMEM batching, SAL prefetch, SAM building | [#58](https://github.com/fg-labs/bwa-mem3/pull/58) | — |
| kswv per-strip L1 prefetches | sw | Adds the per-strip L1 prefetches that `kswv512_16` already had to the four `kswv` kernels that lacked them; stops first-touch L1 stalls | [#70](https://github.com/fg-labs/bwa-mem3/pull/70) | — |
| `SMEM_LOCKSTEP_N` 8 → 16 | seed | Wider lockstep batch exposes more memory-level parallelism on modern cores | [#75](https://github.com/fg-labs/bwa-mem3/pull/75) | — |
| Closed-form ungapped HIT | sw | Replaces the per-base Kadane walk with a direct store when `total_mis == 0` | [#77](https://github.com/fg-labs/bwa-mem3/pull/77) | — |
| On-stack ksort buffer | sort | Removes a per-call `malloc` for small arrays (typical n = 5–30 alignment regions) | [#78](https://github.com/fg-labs/bwa-mem3/pull/78) | — |
| Skip zero-init in libsais build | index | Drops value-initialization of unpack + SA buffers later fully overwritten (tens of GiB of zero-fill avoided) | [#80](https://github.com/fg-labs/bwa-mem3/pull/80) | — |
| Single-binary SIMD dispatch | dispatch | Replaces the multi-binary `execv` launcher with in-process per-host kernel selection | [#83](https://github.com/fg-labs/bwa-mem3/pull/83) | — |
| x86 baseline → avx2 | dispatch/sw | Restores 256-bit auto-vec on non-kernel TUs (~+15% wgs / +11% wes vs the sse41 baseline) | [#84](https://github.com/fg-labs/bwa-mem3/pull/84) | — |
| Cap avx512bw autovec at 256-bit | dispatch/sw | Avoids the 512-bit auto-vec frequency downclock on some x86 parts (the PR also adds an unrelated `bwa_shm` `/dev/shm` preflight) | [#86](https://github.com/fg-labs/bwa-mem3/pull/86) | — |
| Inline `backwardExt` | seed | Forces inlining at the hot SMEM call sites; kills a struct-by-value ABI pass gcc 12+ couldn't elide | [#88](https://github.com/fg-labs/bwa-mem3/pull/88) | — |
| SIMD host-floor precheck | dispatch | Fails fast with a clear message when a host lacks the build's required ISA (safe multi-arch deployment) | [#95](https://github.com/fg-labs/bwa-mem3/pull/95) | — |
| Stable tie-breaks + pdqsort | sort | Deterministic alnreg ordering and `pdqsort` at the dedup-patch sort sites (see equivalence note above) | [#123](https://github.com/fg-labs/bwa-mem3/pull/123) | — |
| FASTQ reader fast path | i/o | Content-detecting reader over libdeflate BGZF for faster input decompression/parse | [#128](https://github.com/fg-labs/bwa-mem3/pull/128) | — |
| Recover 8-bit banded SW (≥128 bp) | sw | Keeps long reads in the cheaper 8-bit lane width where valid | [#140](https://github.com/fg-labs/bwa-mem3/pull/140) | — |
| Gotoh gaps from H | sw | Derives extension gaps from H (standard Gotoh) rather than M in the recovered 8-bit path | [#141](https://github.com/fg-labs/bwa-mem3/pull/141) | — |
| Drop dead `qlen[]` param | sw | Removes an unread parameter from the three 8-bit kernels (cleanup; compile-output-identical) | [#143](https://github.com/fg-labs/bwa-mem3/pull/143) | — |
| Long-read kernel parity test | sw/test | Promotes the 8-bit/16-bit byte-identity harness into a CI doctest | [#144](https://github.com/fg-labs/bwa-mem3/pull/144) | — |
| Short-circuit re-baseline scan | sw | Skips the inert per-row re-baseline scan in the banded kernel | [#147](https://github.com/fg-labs/bwa-mem3/pull/147) | — |
| Remove dead SW code paths | sw | Deletes unreachable `SORT_PAIRS` / non-CORE / SSE2-polyfill code (no behavior change) | [#148](https://github.com/fg-labs/bwa-mem3/pull/148) | — |
| Vectorize epilogue side-channel | sw | Vectorizes the per-row epilogue side-channel loop | [#149](https://github.com/fg-labs/bwa-mem3/pull/149) | — |
| Bound getScores prefetch reads | sw | Bounds `getScores8/16` prefetch reads to the padding contract (hardening) | [#150](https://github.com/fg-labs/bwa-mem3/pull/150) | — |
| Unsigned 8-bit h0-prefix seed | sw | Widens the 8-bit h0-prefix seed to unsigned `[0,255]` | [#151](https://github.com/fg-labs/bwa-mem3/pull/151) | — |
| `--profile` stage timing | prof | Off-by-default read/proc/write + disk/decompress/parse breakdown for diagnosing pipeline scaling (no overhead when off) | [#152](https://github.com/fg-labs/bwa-mem3/pull/152) | — |
| zlib-ng inflate + 3rd worker | i/o | Vendored zlib-ng inflate path with chunk cap and an added pipeline worker (~2.2× faster read stage; up to −7.8% wall at 96 cores) | [#153](https://github.com/fg-labs/bwa-mem3/pull/153) | — |
| Right-size SA staging buffers | index | Sizes the `pos_ar`/`map_ar` staging buffers to the actual `min(s, max_occ)` write count instead of the uncapped SA-interval sum | [#157](https://github.com/fg-labs/bwa-mem3/pull/157) | — |
| `gtle` contract test | sw/test | Enforces `gtle` byte-identity when `gscore > 0`; documents the `gscore == 0` query-end tail divergence | [#158](https://github.com/fg-labs/bwa-mem3/pull/158) | — |
| NEON SW tuning | sw | Replaces `sse2neon` all-zero `movemask` tests with single-instruction `vmaxvq` horizontal reductions in the hot SW scans | [#160](https://github.com/fg-labs/bwa-mem3/pull/160) | — |
| AVX2 SW tuning | sw | Relieves the port-5 `vpblendvb`/`vpshufb` bottleneck in the two SW kernels (byte-identical, Zen3-verified) | [#161](https://github.com/fg-labs/bwa-mem3/pull/161) | — |
| AVX2 16-bit `kswv256_16` | sw | Adds the AVX2 16-bit mate-rescue kernel so AVX2 hosts no longer fall back to scalar `ksw_align2` for 16-bit rescue | [#162](https://github.com/fg-labs/bwa-mem3/pull/162) | — |
| NEON `movemask` parity test | sw/test | Scalar-vs-NEON unit test for the `movemask` helpers; guards the #160 rewrite against silent mask collapse | [#164](https://github.com/fg-labs/bwa-mem3/pull/164) | — |

### Open / in progress

Not yet merged to `main`:

| Item | Stage | Mechanism | bwa-mem3 PR |
|------|-------|-----------|-------------|
| AVX2 8-bit wrapper prefetch | sw | Adds the missing next-batch ref/query software-prefetch to `smithWatermanBatchWrapper8` — the lone SW wrapper that lacked it (follow-up to #161) | [#163](https://github.com/fg-labs/bwa-mem3/pull/163) |

Several correctness and crash fixes underpin the long-read SW, indexing, and
high-throughput work rather than adding speed themselves: SMEM read positions
widened `int16_t` → `int32_t` to stop a long-read `SIGSEGV`
([#142](https://github.com/fg-labs/bwa-mem3/pull/142), merged); a persistent
`kt_for` worker pool that fixes a multi-chunk `SIGSEGV` under mimalloc v3
([#154](https://github.com/fg-labs/bwa-mem3/pull/154), merged); and `mem_lim`
widened to `int64` to stop an SA-staging buffer overflow on highly repetitive
seeds ([#156](https://github.com/fg-labs/bwa-mem3/pull/156), merged).

---

**See also:**
[Performance → Overview](../performance/overview.md) ·
[Performance → PGO build](../performance/pgo.md) ·
[Correctness fixes](correctness.md) ·
[Build & infrastructure](build-infra.md) ·
[bwa-mem3-bench](https://github.com/fg-labs/bwa-mem3-bench)
