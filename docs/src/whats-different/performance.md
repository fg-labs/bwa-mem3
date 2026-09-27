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
| `BWA3_CHAIN_STATS=1` | Print, once at exit, how the chaining, Pass-3 and contained-seed-skip fast paths resolved (`[chain-stats] …`): reads indexed, reads that fell back, queries answered by bucket walks versus full scans, reads the flat chaining index handed to the B-tree (equal positions, or over the cap), and how the contained-seed extension skip resolved its deferred seeds (`contained_deferred`, of which `contained_purged` had their banded-SW skipped and `contained_extended` ran in the second batch). Measurement only; output is unchanged. |
| `BWA3_CHAIN_FLAT_CAP=<n>` | Largest number of chains a read may have on the flat chaining index before it is replayed through the B-tree (default 512, bounding the index's O(n) sorted insert). `0` sends every read to the B-tree. A malformed or negative value is reported to stderr and the default used. Output is identical at every value. |

## Mate rescue

Byte-identical changes to mate rescue, the batched Smith-Waterman of a read's
mate against the reference window its insert size allows. Output was verified
identical (md5 of the non-`@PG` records) to the previous `main` on 5M-pair WGS and
WES slices (150 bp paired, hg38) at the default scoring on AWS Graviton 4 (NEON),
both sides run at `-t 16 -K 160000000` (a pinned `-K` fixes the batch composition,
which the insert-size estimate, and so rescue and pairing, depend on). The x86
port of pruning and banding was verified the same way against v0.13.0 on the same
slices, from multi-tier builds (`BASELINE_ARCH=avx2`, clang 19), on AMD EPYC 7R13
(Zen 3, kswv at the AVX2 tier) and AMD EPYC 9R45 (Zen 5, kswv at the AVX-512BW
tier); on both, the filter is the SSE4.1 / SSSE3 port and the band kernel the
AVX2 one. That real-data check predates the x86 hookup of the repeated-job
reuse below, whose x86 identity rests on CI only. CI's x86 checks run on a
generated fixture only, not on real data: the x86 rows (AVX2 build) run
`rescue_prune_identity.sh` under forced AVX2 and AVX-512BW kswv tiers
(AVX-512BW where the runner has it) and the filter and band harnesses on
generated jobs.

- **Exact rescue pruning.** Before a rescue job is staged, a filter bounds the
  best local score from the exact 5-mer matches between the mate and the window.
  A job proven unable to reach the rescue threshold (`min_seed_len * a`) is not
  run and takes the ordinary failing-rescue path; otherwise only the rows that
  can hold an alignment at that threshold are computed, which reproduces every
  field the rescue consumes (score, positions, suboptimal score). Derived for any
  scoring whose mismatch penalty (`-B`), deletion open alone and insertion open
  plus extend each cost at least four matches, the insertion strictly more, with
  a threshold above four matches and each gap type's open plus extend at most
  255 (the byte kswv's 8-bit kernels hold it in): the default
  `-A 1 -B 4 -O 6 -E 1`, and for example `-B 6` or `-O 8 -E 2`; the bound then
  counts exact 5-mers (`rescue_prune_params::from` states the conditions). Where
  every separator costs more, the bound counts longer K-mers, which cuts random
  hits 4x per step: K - 1 matches may cost at most the mismatch penalty, the
  deletion open and the insertion open plus extend, so `-B 6` counts 7-mers and
  `-x intractg` 8-mers, the longest any filter takes. On aarch64 only
  (`BWA3_RESCUE_PRUNE_KMAX`), and never under `--meth`, which stays at 5-mers. Other
  scorings, windows or mates with an N, mates over 1024 bases and windows over
  30000 rows keep the full window. Jobs under `--rescue-kmer` bypass exact
  pruning, but its own k-mer anchor path can still narrow the window. Jobs on
  kswv's 16-bit kernels (a mate too long for the 8-bit ones at the run's `-A`)
  are pruned too: the bound is on scores, not SIMD lanes. Under
  `--meth` it runs on aarch64 only. With EM-seq chemistry it matches
  C-to-T (or G-to-A) converted copies of the window and the mate, at
  `--meth -B 4` and the genomic and neutral scorings; with TAPS chemistry
  (mostly unconverted reads, where converted copies leave little to prune) it
  matches the exact relation of the genomic and neutral matrices instead (a
  read T may pair with a reference T or C on the top strand, a read A with a
  reference A or G on the bottom one, every other base only with itself), at
  its default neutral scoring and at genomic scoring. The collapsed scoring
  (the `--meth` default) keeps the full window under either chemistry. It runs where a
  SIMD filter carries it: on aarch64 (NEON) and on the x86 AVX2 and AVX-512BW
  builds (an SSE4.1 / SSSE3 port of the filter), at every seed length except
  from `-k 25` at the AVX-512BW kswv tier (at any `-A`), where the 64-lane kswv
  is cheap enough that the filter costs more than it saves (wgs-5M and wes-5M,
  prune on vs off at `-k 25` to 40, on an x86_64 host at the AVX-512BW
  tier). The SIMD filters cover every
  admitted scoring with `-A` of at most 16; the others run a scalar filter,
  which pays on aarch64 but not against x86's kswv, so x86 keeps the full window
  there (wgs-5M and wes-5M at `-O 8 -E 2` and `-x intractg`: a gain on
  an aarch64 host at the NEON tier, a loss on the x86_64 AVX-512BW host).
  Elsewhere the full window is always computed.
  ([#541](https://github.com/fg-labs/bwa-mem3/pull/541); x86:
  [#538](https://github.com/fg-labs/bwa-mem3/pull/538); 16-bit jobs and the
  AVX-512BW seed-length gate: [#542](https://github.com/fg-labs/bwa-mem3/pull/542);
  the TAPS relation: [#546](https://github.com/fg-labs/bwa-mem3/pull/546); K-mers longer than 5:
  [#548](https://github.com/fg-labs/bwa-mem3/pull/548))
- **Banded rescue DP (NEON, AVX2).** For a narrowed job, the rescue DP runs only
  inside the diagonal bands of the K-mer components that can reach the
  threshold, 16 bands per NEON vector or 32 per AVX2 vector, and the job's
  score, end positions and suboptimal score are reassembled from the bands'
  per-row maxima. A per-job
  cost model keeps the full hull when banding would not pay. Same scope as the
  pruning on 8-bit jobs, at every scoring it admits; under `--meth`, by
  default (`BWA3_RESCUE_BAND_METH`), only where it filters under the relation
  (TAPS), with each hypothesis's bisulfite matrix. The kernels take the run's
  match, mismatch and both gap types' costs.
  ([#541](https://github.com/fg-labs/bwa-mem3/pull/541); AVX2:
  [#538](https://github.com/fg-labs/bwa-mem3/pull/538); kept to 8-bit jobs
  when pruning took in 16-bit ones: [#542](https://github.com/fg-labs/bwa-mem3/pull/542);
  under `--meth` with TAPS pruning: [#546](https://github.com/fg-labs/bwa-mem3/pull/546); K-mer
  components: [#548](https://github.com/fg-labs/bwa-mem3/pull/548))
- **Banded start recovery (NEON, AVX2).** The second rescue pass, which finds where
  the best alignment starts, runs in a diagonal band derived from the first
  pass's score and end: an alignment of that score can hold only a bounded
  number of gapped bases, so the band holds it, and the first row and column
  reaching the score are the same as in the full pass. Used for every 8-bit job
  (banded or not in the first pass, at any seed length and any scoring the band
  kernels take, including ones the pruning refuses, such as `-B 3`) when the
  band is cheaper than the full pass; under `--meth`, with the bisulfite
  matrices, by default (`BWA3_RESCUE_BAND_METH`) wherever `--meth` does not
  prune on converted copies (on x86, at
  the default collapsed scoring, and with TAPS chemistry, pruned or not).
  ([#541](https://github.com/fg-labs/bwa-mem3/pull/541); AVX2:
  [#538](https://github.com/fg-labs/bwa-mem3/pull/538); TAPS with pruning: [#546](https://github.com/fg-labs/bwa-mem3/pull/546))
- **11-op rescue cell (NEON, AVX2, AVX-512BW; 8- and 16-bit).** When the
  open-plus-extend sums (`-O` + `-E`) of insertions and deletions are equal and
  fit the kernel's lane (a byte for the 8-bit kernels), and no gap cost is
  negative, the rescue kernels build each cell from the score before the in-row
  gap and open both gaps from one subtract: 11 vector operations per cell instead
  of 13 in the NEON 8-bit kernel (one max/subtract fewer per cell, and no
  reference-padding mask, in the others), and row i+1 no longer waits on row
  i's gap chain. Every score and position it emits is unchanged: output (md5 of
  the non-`@PG` records) matched `BWA3_RESCUE_FSCAN=0` on the same 5M-pair WGS
  and WES slices, at `-t 16 -K 160000000`, at the default scoring (8-bit bodies)
  and at `-A 2` (every rescue job on the 16-bit bodies), on AWS Graviton 4
  (NEON), AMD EPYC 7R13 (Zen 3, AVX2) and AMD EPYC 9R45 (Zen 5, AVX-512BW).
- **AVX2 and AVX-512BW rescue kernels: signed-domain 8-bit cell, deferred 16-bit query end.**
  x86 has no unsigned-plus-signed saturating add, so the AVX2 and AVX-512BW
  8-bit FScan bodies now keep their scores as H - 128 in signed bytes, where each cell scores with
  one signed saturating add instead of a biased add and a de-biasing subtract
  (as NEON's USQADD body does; `BWA3_RESCUE_USQADD=0` restores the biased form).
  Each is selected only when every gap constant fits a positive signed byte, i.e.
  the open-plus-extend sum is at most 127; above that the biased body runs even
  with `BWA3_RESCUE_USQADD` on.
  The AVX2 16-bit body finds each row's query end after the row, from per-block
  checkpoints of the row maximum, instead of tracking it in every cell. Every
  score and position is unchanged: output (md5 of the non-`@PG` records)
  matched v0.13.0 on the 5M-pair WGS slice at `-t 16 -K 160000000` (the same
  fixed batch size, and so the same batches, in both runs) on AMD EPYC 7R13
  (Zen 3, AVX2): the signed-domain cell at the default scoring (also with
  `BWA3_RESCUE_USQADD=0`), `-B 6`, `-O 8 -E 2`, `-x intractg`, `-k 25`,
  `-A 2` and `-A 3 -B 12 -O 18 -E 3`; the deferred query end at the default
  scoring, `-A 2` (every rescue job on the 16-bit body), `-A 2 -B 6`,
  `-A 2 -B 8 -O 12 -E 2` and `-A 3 -B 12 -O 18 -E 3`. The AVX-512BW 8-bit
  body takes the same signed domain by the same argument; its gate is the kswv
  unit tests on an AVX-512BW host. The NEON and SSE4.1 kernels are unchanged.
  ([#542](https://github.com/fg-labs/bwa-mem3/pull/542); AVX-512BW:
  [#546](https://github.com/fg-labs/bwa-mem3/pull/546))
- **Fused, two-row banded cell (NEON, AVX2).** The banded rescue DP uses the same
  fused cell (both gaps opened from one saturating subtract), steps two rows at
  a time, and reads the query end directly rather than from a per-row snapshot.
  Same scope as the banded DP; every value the rescue reads is unchanged.
- **Skipped or incremental post-rescue dedup.** After each rescue, the mate's
  region list was sorted and deduplicated again, including after rescues that
  added nothing. A dedup that finds every reference end distinct and drops
  nothing is a fixed point, so repeating it on the unchanged list is skipped;
  adding one new region to such a list is done in one linear pass when that
  region's reference end is distinct, it forms no redundant pair, and no two
  regions share a score, reference start and query start, and by the full
  dedup otherwise. The output of every dedup is unchanged. Runs on every
  architecture and at every scoring; it is covered by the Graviton 4 check
  above and, since v0.13.0 predates it, by the Zen 3 and Zen 5 x86 checks, as
  well as by the dedup unit tests and `rescue_prune_identity.sh` on a generated
  fixture in CI (AVX2, and AVX-512BW where the runner has it).
  ([#537](https://github.com/fg-labs/bwa-mem3/pull/537))
- **Band planning from the filter's own components (NEON, x86).** Band planning
  takes the diagonal components at the rescue threshold from the SIMD filter,
  which has just found them, instead of rescanning the filter's per-diagonal
  arrays; it skips components that cannot hold a higher threshold, and the
  filter returns its previous decision for a job repeating the previous one
  byte for byte. Plans and decisions are unchanged.
  ([#537](https://github.com/fg-labs/bwa-mem3/pull/537); x86:
  [#538](https://github.com/fg-labs/bwa-mem3/pull/538))
- **Repeated rescue jobs answered from an earlier result (NEON, x86).** A rescue
  job that repeats one of the thread's last eight filtered jobs byte for byte,
  the same mate against an identical window as anchors in identical repeat
  copies produce, reads that job's rescue-kernel result instead of being
  enqueued again. The filter's key of the inputs picks the candidates, and the
  staged window, mate, lengths, score gate and hull offset are compared byte
  for byte before a result is reused, so every value the rescue reads is
  unchanged. `BWA3_RESCUE_REPEAT=0` turns it off.
  ([#537](https://github.com/fg-labs/bwa-mem3/pull/537); x86:
  [#538](https://github.com/fg-labs/bwa-mem3/pull/538))
- **Window bounds only after the rescue.** The step that reads rescue results
  back computed each job's reference window again, bases included, though it
  needs only the clamped bounds and the contig; the bases are now fetched only
  for the rare job that falls back to the scalar aligner.
  ([#537](https://github.com/fg-labs/bwa-mem3/pull/537))
- **Rescue window prefetch.** The reference windows mate rescue will fetch for
  a pair are hinted to the cache two pairs ahead (every window of a read's top
  anchor lies within the largest usable insert-size bound of it). Pure hints;
  skipped when the insert-size model is too loose for the span to be small.
  ([#537](https://github.com/fg-labs/bwa-mem3/pull/537))

| Variable | Effect |
|---|---|
| `BWA3_RESCUE_PRUNE=0` | Turn off exact rescue pruning and the banded passes: every rescue window runs in full through the rescue kernel (the reference path for identity checks). `--rescue-kmer`, which narrows windows on its own, is unaffected. Default on where pruning runs. |
| `BWA3_RESCUE_PRUNE_MAX_HITS=<n>` | Keep the full window when the mate and window share more than `n` exact K-mer hits, where the filter would cost more than it saves (default 1000 on aarch64 when banding is on, unless a `--meth` run leaves its pruned windows unbanded, as EM-seq does by default; 400 otherwise, and always 400 on x86; [#541](https://github.com/fg-labs/bwa-mem3/pull/541), TAPS: [#546](https://github.com/fg-labs/bwa-mem3/pull/546), K-mers: [#548](https://github.com/fg-labs/bwa-mem3/pull/548)). It only chooses between exact paths, so output does not depend on its value by design. |
| `BWA3_RESCUE_PRUNE_KMAX=<n>` | The longest K-mer the filter may count where the scoring admits it, from 5 to 8 (values outside act as the nearer end): default 8 on aarch64, 5 elsewhere, and always 5 under `--meth`. The default scoring admits only 5. It only chooses between exact paths, so output does not depend on its value by design.([#548](https://github.com/fg-labs/bwa-mem3/pull/548)) |
| `BWA3_RESCUE_PRUNE_REL=<n>` | Under `--meth` at the genomic and neutral scorings, on aarch64 (x86 does not prune `--meth`), how the filter relates bases: `1` (default) the exact relation for TAPS chemistry and converted copies for EM-seq, `0` converted copies only (TAPS is then not pruned), `2` the relation for every such run; values above 2 act as 2. It only chooses between exact paths, so output does not depend on its value by design.([#546](https://github.com/fg-labs/bwa-mem3/pull/546)) |
| `BWA3_RESCUE_PRUNE_STATS=1` | Print, once at exit, how the filter decided (`[RESCUE_PRUNE] jobs=… full=… b1=… b2=… rows_in=… rows_kept=… jobs16=… b1_16=… b2_16=… kmer_jobs=… memo_hits=… reused=… filter_s=… kswv_pass0_s=… band_pass0_s=… kswv_pass1_s=… band_pass1_s=… dedup_run=… dedup_skip=… dedup_run_regs=… dedup_skip_regs=… dedup_insert1=… dedup_insert1_fast=… dedup_s=…`): jobs filtered, and of them how many kept the full window, were proven to fail (`b1`) or were narrowed (`b2`), with the window rows before and after, the number of 16-bit rescue jobs and how many of them were proven to fail or narrowed (`b1_16`, `b2_16`, included in `b1` and `b2`), the jobs the SIMD filter decided with K-mers longer than 5 (`kmer_jobs`), the jobs the SIMD filter answered from its repeat memo, the jobs answered from an identical recent job's result instead of being enqueued, and the thread-summed seconds of each rescue stage; how many post-rescue dedups ran, were skipped, took the one-region insert and of those were done in one pass; and how the banded DP resolved (`[RESCUE_BAND] banded_parents=… …`); the 16-bit decisions (`b1_16`, `b2_16`): [#542](https://github.com/fg-labs/bwa-mem3/pull/542), `kmer_jobs`: [#548](https://github.com/fg-labs/bwa-mem3/pull/548). Measurement only; output is unchanged. |
| `BWA3_RESCUE_DEDUP_SKIP=0` | Run every post-rescue dedup in full instead of skipping one proven to be a no-op or adding a single new region in one pass. Default on. It only chooses between exact paths, so output is the same either way by design. |
| `BWA3_RESCUE_REPEAT=0` | Enqueue every rescue job instead of answering one that repeats one of the thread's last eight filtered jobs byte for byte (the same mate against an identical window) from that job's result. Default on. It only chooses between exact paths, so output is the same either way by design. |
| `BWA3_RESCUE_BAND=0` | Run every narrowed job through the rescue kernel on its whole hull instead of banded (and the hit gate back to 400). Default on where banding runs. |
| `BWA3_RESCUE_BAND_COST=<pct>` | Band a narrowed job only when its band cells cost less than `pct` % of the hull's (default 85; 0, i.e. no first-pass banding, where the rescue kernel runs at the AVX-512BW tier, whose 64-lane sweep of the hull is cheaper than the 32-lane banded DP). It only chooses between exact paths, so output does not depend on its value by design. |
| `BWA3_RESCUE_BAND_R2=0` | Run the rare second round (a first round that cannot prove its result final) through the rescue kernel on the hull instead of banded. It only chooses between exact paths, so output is the same either way by design. |
| `BWA3_RESCUE_BAND_P1=<n>` | Which second-pass (start recovery) jobs run banded: `0` none, `1` only jobs banded in the first pass, `2` (default) every eligible 8-bit job; values above 2 act as 2. It only chooses between exact paths, so output does not depend on its value by design. |
| `BWA3_RESCUE_BAND_P1_COST=<pct>` | Band a second-pass job only when its band's per-row cells cost less than `pct` % of the full pass's (default 130). It only chooses between exact paths, so output does not depend on its value by design. |
| `BWA3_RESCUE_BAND_TIGHT=<n>` | Threshold offset of the first-round band for a lone near-perfect primary (default 8; 0 disables). It only chooses between exact paths, so output does not depend on its value by design. |
| `BWA3_RESCUE_FSCAN=0` | Use the original rescue cells instead of the 11-op cell, in the rescue kernels and in the banded DP. A value starting with `0` turns it off and anything else leaves it on; it is not reported. It only chooses between exact paths, so output is the same either way by design. |
| `BWA3_RESCUE_BAND_KERNEL=<n>` | Which banded-DP kernel runs while `BWA3_RESCUE_FSCAN` is on: `0` the original cell, `1` the fused cell one row at a time, `2` (default) the fused cell two rows at a time; values above 2 act as 2. It only chooses between exact paths, so output does not depend on its value by design. |
| `BWA3_RESCUE_USQADD=0` | In the NEON 8-bit rescue kernel and the AVX2 and AVX-512BW 8-bit FScan bodies, use the biased add / subtract pair per cell instead of one saturating add (on x86, of the signed H - 128 domain, which runs only when the open-plus-extend sum is at most 127; AVX2: [#542](https://github.com/fg-labs/bwa-mem3/pull/542), AVX-512BW: [#546](https://github.com/fg-labs/bwa-mem3/pull/546)). It only chooses between exact paths, so output is the same either way by design. |
| `BWA3_RESCUE_ROWPAIR=0` | In the NEON rescue kernels (8- and 16-bit), sweep one target row at a time instead of two. It only chooses between exact paths, so output is the same either way by design. |
| `BWA3_RESCUE_LAZYQE=0` | In the NEON two-row sweep, find each row's query end inline instead of after the row. It only chooses between exact paths, so output is the same either way by design. |
| `BWA3_RESCUE_BAND_SHIFT=0` | Turn off the per-lane band shift that aligns the query offsets of the 16 bands in a vector. It only chooses between exact paths, so output is the same either way by design. |
| `BWA3_RESCUE_BAND_METH=<n>` | Under `--meth`, where the banded DP runs (both passes, with each hypothesis's matrix): `1` (default) everywhere except on top of converted-copy pruning (EM-seq), `0` never, `2` always; values above 2 act as 2. It only chooses between exact paths, so output does not depend on its value by design. ([#546](https://github.com/fg-labs/bwa-mem3/pull/546)) |

Every integer knob above (`=<n>` or `=<pct>`) takes a non-negative decimal integer of at most 2147483647, with no sign or surrounding whitespace. An unset or empty knob silently takes the default; any other invalid value is reported to stderr and the default used. The on/off knobs (`=0`) are off for any value starting with `0` and on otherwise.

The byte-identity above was measured at the default knob values. Beyond the defaults, CI tests only these values, on generated inputs. On aarch64 (NEON): `BWA3_RESCUE_BAND=0` and `BWA3_RESCUE_BAND_COST=100000000` (the cost gate forced open) against the full-window reference (`BWA3_RESCUE_PRUNE=0 BWA3_RESCUE_BAND=0`) on a generated paired-end fixture at `-t 1` and `-t 4`, and `BWA3_RESCUE_BAND_COST=100000000` on generated rescue jobs checked field by field, alone and with each of `BWA3_RESCUE_PRUNE_MAX_HITS=100000`, `BWA3_RESCUE_BAND_SHIFT=0`, `BWA3_RESCUE_BAND_KERNEL=1` and `BWA3_RESCUE_FSCAN=0`. On aarch64 and x86 (AVX2): `BWA3_RESCUE_FSCAN=0` against the default in the rescue kernels' unit tests, field by field. Other values are exact by construction but untested.

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
