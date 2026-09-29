# Performance Improvements

This page covers the performance work carried in bwa-mem3 on top of upstream
bwa-mem2. Almost every change listed here is a throughput, memory, or
supporting (test/hardening/cleanup) change that preserves the aligner's output;
the one exception is the deterministic
tie-break ordering in [#123](https://github.com/fg-labs/bwa-mem3/pull/123),
which can reorder equal-scoring alignments relative to upstream (see
[Equivalence with bwa-mem2](equivalence.md) for the full, audited list of where
bwa-mem3 output diverges).

For a reader-friendly grouping of *what drives the speedup* — by machine
architecture, hot-path rewrite, indexing, allocation/I/O, and build-time — see
[Performance → Overview](../performance/overview.md). For current benchmark
numbers across architectures and workloads, see
[Benchmarks](../performance/benchmarks.md), the canonical
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

## Third-pass re-seeding lockstep width (PR #514)

The third seeding pass (`bwtSeedStrategy` occurrence-bounded re-seeding) has its
own lockstep driver, overlapping `BWTSEED_LOCKSTEP_N` reads' forward-extension
walks so their `cp_occ` cache misses issue together. Since
[#530](https://github.com/fg-labs/bwa-mem3/pull/530) it is on by default on every
platform and at every thread count (`BWA3_BWTSEED_LOCKSTEP=0`/`1` pins it). x86
used to turn it off whenever `-t` exceeded the physical core count, on the theory
that a busy SMT sibling hides the same latency; at the depth-24 default it does
not. Turning it on at `-t` = 2× physical cores cut whole-aligner wall on every SMT
host measured, with byte-identical alignment records at fixed `-t`/`-K` (headers
excluded). Wall change, lockstep off → on, on 5M-pair WGS and WES slices (150 bp
paired, hg38; AVX2 build; 5–7 interleaved reps per arm):

| Host | CPU | Threads / physical cores | WGS | WES |
|---|---|---|---|---|
| c6a.4xlarge | AMD Zen 3 | 16 / 8 | −8.2% | −7.2% |
| c7i.4xlarge | Intel Sapphire Rapids | 16 / 8 | −9.6% | −8.2% |
| c6a.8xlarge | AMD Zen 3 | 32 / 16 | −8.2% | −7.2% |
| c7i.8xlarge | Intel Sapphire Rapids | 32 / 16 | −8.9% | −8.3% |

Other SMT microarchitectures (for example Zen 4 with
SMT) and oversubscription beyond `-t` = 2× physical cores were not measured for
speed; output identity there rests on the driver's SMEM-order parity design, not
on a measurement. `test/regression/bwtseed_lockstep_driver_identity.sh`
pins the scalar driver in CI and checks it against the default.

Its depth is a distinct knob from the phase-2 SMEM width above. The default is
**24**, resolved at startup into a runtime value. Measured whole-aligner wall on a
5M-read WGS slice (150 bp paired NovaSeq reads, hg38), each host at threads ≤
physical cores (held fixed so the sweep compares depths, not thread regimes),
sweeping the depth via
`BWA3_BWTSEED_LOCKSTEP_N` on a single binary per host:

| Host | SIMD tier | Threads | Best depth | Δ vs depth 8 |
|---|---|---|---|---|
| Intel Sapphire Rapids (c7i.8xlarge, 16 cores) | avx512bw | 16 | ~24 (flat plateau to 64) | −1.8% |
| AMD Zen 4 (c7a.8xlarge, 16 cores) | avx512bw | 16 | 24 (peak, rolls off after) | −1.7% |
| AWS Graviton4 (c8g.4xlarge, 16 cores) | neon | 16 | ~16–24 | −0.8% |
| Apple M3 Ultra (20 perf cores) | neon | 8 | ≥24 | −2.2% |
| Apple M2 Max (8 perf cores) | neon | 8 | flat | ~0% (no regression) |

24 is a win or neutral on all of them; Zen 4 rolling off past 24 caps a portable
default there. At a fixed `-t`/`-K` (equal thread count and identical `-K` batch size,
so batch boundaries and `mem_pestat` are unchanged), the **alignment records** are
byte-identical across depths — the depth changes only how forward-extension walks are
batched within a run, not which seeds are emitted; the header block is not part of this
claim. Verified on the same 5M-read WGS parity workload by parity tests
(lockstep == scalar at depths 8 / 24 / 64) and a whole-genome alignment-record
md5 gate matching the prior depth-8 output on x86 (clang-19, avx512bw tier) and arm64
(NEON). So it is safe to tune per host.

| Variable | Effect |
|---|---|
| `BWA3_BWTSEED_LOCKSTEP_N=<n>` | Pin the third-pass lockstep depth to `<n>` (1–64). A value above 64 is clamped; a non-positive or malformed value is reported to stderr and ignored (the default is used). An unset or empty value (`BWA3_BWTSEED_LOCKSTEP_N=`) silently keeps the compiled default, with no message. |
| `BWA3_BWTSEED_LOCKSTEP=0`\|`1` | Force the third-pass lockstep driver off / on (default: on). |

## Batched `-H` header ingestion (PR #49, closes issue #37)

Passing a large header file via `-H <file>` re-ran `strlen` on the growing
header string and called `realloc` on every input line, making ingestion O(n²)
in the number of header lines. For a ~70 MB / ~1.5 M-line header (reported in
upstream [bwa-mem2#204](https://github.com/bwa-mem2/bwa-mem2/pull/204)) this
caused runtimes exceeding 10 minutes before alignment started.

The fix introduces `bwa_insert_header_file`, a batched helper that reads the
`@`-prefixed lines into a single growable buffer in one pass, escaping each
retained line in place as soon as it is complete (so a trailing backslash on
one line can never consume the separator before the next), then joins the
result onto the existing header with a single realloc rather than one per
retained line. It reads line by line rather than sizing the buffer with
`fseek`/`ftell`, so it also works on a non-seekable `-H` stream (a pipe,
`/dev/stdin`, or `-H <(...)` process substitution) — where `ftell` returns
`-1` and the earlier size-first approach silently discarded the whole header
file. A single header line longer than a 64 KiB budget is rejected with a
fatal error rather than silently truncated. A regression test
(`test/header_insert_test.cpp`) diffs the batched path against the pre-patch
per-line baseline across the supported edge cases, including the oversize-line
rejection.

## libsais FM-index construction (PR #57)

`bwa-mem3 index` now builds the FM-index using
[libsais v2.9.1](https://github.com/IlyaGrebnov/libsais) (Ilya Grebnov)
instead of the sais-lite (Yuta Mori saisxx) library that bwa-mem2 inherited.
libsais is actively maintained, supports OpenMP-parallel induced sorting, and
produces a byte-identical FM-index. An index built by `bwa-mem2 index` is read
without re-indexing; a `bwa` (v1) index uses a different format and must be
rebuilt with `bwa-mem3 index` (see
[Coming from bwa or bwa-mem2](../getting-started/migrating.md)).

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
are maintained under [Benchmarks](../performance/benchmarks.md).

## Chaining and extension setup

Byte-identical changes to the path between seeding and extension. Each keeps the
code it replaces reachable as a fallback and pins itself against it. Output was
verified identical (md5 of the non-`@PG` records) to the previous `main` on 5M-pair
WGS and WES slices (150 bp paired, hg38) at the default scoring, `-A 2`, `-B 6` and
`-O 8 -E 2`, on AMD Zen 3 (AVX2 build) and AWS Graviton 4 (NEON), and on three
`--meth` datasets on Graviton 4.

- **Run-length ungapped walk.** The ungapped fast path's score walk steps one
  mismatch at a time instead of one base at a time: between mismatches the score
  rises strictly, so a run's maximum is its last value. Outside the scoring
  envelope where that holds (`a > 0`, `b >= 0`, positive start score) the
  per-base walk runs. Cross-checked against the per-base walk under the opt-in
  `BWA_MEM3_DEBUG_UNGAPPED_XCHECK` build.
- **Pass-3 kept-set index.** The post-extension sweep that drops seeds contained
  in an earlier alignment answers each seed from an index over the kept seeds'
  alignments, instead of scanning every alignment of the read. If a kept
  alignment is ever purged the read falls back to the scan. Cross-checked against
  the scan under `BWA_MEM3_DEBUG_P3_XCHECK`.
- **`cal_max_gap` memo.** The per-seed reference-window derivation reads the
  maximum gap from a per-thread table over query lengths, keyed on the scoring
  options.
- **Flat chaining index.** Seed chaining finds each seed's closest chain in a
  sorted array of chain positions with a branchless binary search instead of the
  B-tree probe. Where the B-tree's answer depends on its node layout (two chains
  with the same position) or a read has more than `BWA3_CHAIN_FLAT_CAP` chains,
  the read is replayed through the unchanged B-tree path.

| Variable | Effect |
|---|---|
| `BWA3_CHAIN_STATS=1` | Print, once at exit, how the chaining and Pass-3 fast paths resolved (`[chain-stats] …`): reads indexed, reads that fell back, queries answered by bucket walks versus full scans, and reads the flat chaining index handed to the B-tree (equal positions, or over the cap). Measurement only; output is unchanged. |
| `BWA3_CHAIN_FLAT_CAP=<n>` | Largest number of chains a read may have on the flat chaining index before it is replayed through the B-tree (default 512, bounding the index's O(n) sorted insert). `0` sends every read to the B-tree. A malformed or negative value is reported to stderr and the default used. Output is identical at every value. |

## Mate rescue

Byte-identical changes to mate rescue, the batched Smith-Waterman of a read's
mate against the reference window its insert size allows. Output was verified
identical (md5 of the non-`@PG` records) to the previous `main` on 5M-pair WGS and
WES slices (150 bp paired, hg38) at the default scoring on AWS Graviton 4 (NEON),
both sides run at `-t 16 -K 160000000` (a pinned `-K` fixes the batch composition,
which the insert-size estimate, and so rescue and pairing, depend on).

- **Exact rescue pruning.** Before a rescue job is staged, a filter bounds the
  best local score from the exact 5-mer matches between the mate and the window.
  A job proven unable to reach the rescue threshold (`min_seed_len * a`) is not
  run and takes the ordinary failing-rescue path; otherwise only the rows that
  can hold an alignment at that threshold are computed, which reproduces every
  field the rescue consumes (score, positions, suboptimal score). Derived for the
  default scoring only (`-A 1 -B 4 -O 6 -E 1`); other scorings, `--meth`,
  `--rescue-kmer`, windows or mates with an N, and the 16-bit path keep the full
  window. It runs on aarch64 only, where a NEON filter carries it; elsewhere the
  full window is always computed.
- **11-op 8-bit rescue cell (NEON).** When the cost of opening a gap (open plus
  extend, `-O` + `-E`) is the same for insertions and deletions and fits a byte,
  and no gap cost is negative, the 8-bit rescue kernel builds each cell from the
  score before the in-row gap and opens both gaps from one saturating subtract:
  11 vector operations per cell instead of 13, and row i+1 no longer waits on row
  i's gap chain. Every score and position it emits is unchanged.

| Variable | Effect |
|---|---|
| `BWA3_RESCUE_PRUNE=0` | Turn off exact rescue pruning (the reference path for identity checks): no window is dropped or narrowed by it. `--rescue-kmer`, which narrows windows on its own, is unaffected. Default on where pruning runs. |
| `BWA3_RESCUE_PRUNE_MAX_HITS=<n>` | Keep the full window when the mate and window share more than `n` exact 5-mer hits, where the filter would cost more than it saves (default 400). A malformed or negative value, or one above 2147483647, is reported to stderr and the default used. Output is identical at every value. |
| `BWA3_RESCUE_PRUNE_STATS=1` | Print, once at exit, how the filter decided (`[RESCUE_PRUNE] jobs=… full=… b1=… b2=… rows_in=… rows_kept=…`): jobs filtered, and of them how many kept the full window, were proven to fail (`b1`) or were narrowed (`b2`), with the window rows before and after. Measurement only; output is unchanged. |
| `BWA3_RESCUE_FSCAN=0` | Use the original 13-op 8-bit rescue cell on NEON. Output is identical either way. |

---

## Full change list

The sections above narrate the load-bearing performance PRs. For the exhaustive
list, the [merged PR list](https://github.com/fg-labs/bwa-mem3/pulls?q=is%3Apr+is%3Amerged+base%3Amain)
(filter on `perf:`) and `git log master..main` are the source of truth; any
performance change that was proposed upstream is recorded on the
[Fork changes vs. upstream](../reference/pr-catalog.md) page.

## Open / in progress

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
[Benchmarks](../performance/benchmarks.md)
