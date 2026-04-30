# Performance Improvements

This page covers the performance work carried in bwa-mem3 on top of upstream
bwa-mem2. Every change listed here preserves byte-identical SAM/BAM output vs
the upstream baseline it was benchmarked against.

For current benchmark numbers across architectures and workloads, see
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

| Item | bwa-mem3 PR | Upstream PR/issue | Status |
|------|-------------|-------------------|--------|
| Lockstep SMEM batching | [#33](https://github.com/fg-labs/bwa-mem3/pull/33) | — | fork-only |
| Batched `-H` header ingestion | [#49](https://github.com/fg-labs/bwa-mem3/pull/49) | [bwa-mem2#204](https://github.com/bwa-mem2/bwa-mem2/pull/204) | fork-only (upstream PR open) |
| Large header performance (issue) | — | [issue #37](https://github.com/fg-labs/bwa-mem3/issues/37) | closed by #49 |
| libsais FM-index construction | [#57](https://github.com/fg-labs/bwa-mem3/pull/57) | — | fork-only |
| Consolidated mapping speedups | [#58](https://github.com/fg-labs/bwa-mem3/pull/58) | — | fork-only |

---

**See also:**
[Performance → Overview](../performance/overview.md) ·
[Performance → PGO build](../performance/pgo.md) ·
[Correctness fixes](correctness.md) ·
[Build & infrastructure](build-infra.md) ·
[bwa-mem3-bench](https://github.com/fg-labs/bwa-mem3-bench)
