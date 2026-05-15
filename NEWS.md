Release notes for 0.3.0 and later live in [`CHANGELOG.md`](./CHANGELOG.md) and the [GitHub Releases](https://github.com/fg-labs/bwa-mem3/releases) page. The 0.2.0 entry below is preserved as historical context.

Release 0.2.0 (2026-05-13)
---------------------------------

Operational / packaging

* Single-binary SIMD dispatch on x86 (#83). The previous multi-binary
  build (`make multi` producing five `bwa-mem3.<tier>` ISA variants plus
  a `runsimd.cpp` launcher that `execv`'d the matching tier) is replaced
  by a single binary that contains compiled kernels for every supported
  tier (`sse41` / `sse42` / `avx` / `avx2` / `avx512bw`) and selects one
  in process at startup via `__builtin_cpu_supports`. Install size drops
  from ~120 MB to ~25 MB; per-call overhead is one indirect branch
  (~0.3 ns after BTB warm-up). No `.<tier>` companion files are produced
  or needed. See `docs/src/developer-guide/launcher.md`.
* `BWAMEM3_FORCE_TIER=<tier>` and `BWAMEM3_DEBUG_SIMD=1` env vars (#83).
  `BWAMEM3_FORCE_TIER` is downgrade-only and replaces the prior
  "exec the `bwa-mem3.sse41` binary" A/B-testing pattern; up-tier or
  unrecognized requests are rejected with a stderr warning.
* `BASELINE_ARCH=avx2` is the new default for non-kernel translation
  units on x86 (#84, supersedes the SSE4.1 floor that PR #83 originally
  shipped with). Override via `make BASELINE_ARCH=<tier>`. AVX-512BW
  hosts using `BASELINE_ARCH=avx512bw` see a small additional speedup
  on Zen 4 with `-mprefer-vector-width=256` (#86) and roughly flat
  results on Sapphire Rapids — see
  `docs/src/whats-different/avx512-baseline.md` for the characterization.
* Host-floor precheck (#95). `bwa-mem3 mem`, `bwa-mem3 index`, and
  `bwa-mem3 shm` refuse to run with exit code 2 and an `[E::bwamem3]`
  stderr message when the host CPU does not meet the build's compile-time
  SIMD floor, instead of SIGILL-ing deep in alignment. `bwa-mem3 version`,
  `--help`, and `-h` are exempt and always succeed.
* `bwa-mem3 version` now prints `SIMD floor:` (build's required minimum)
  and `SIMD runtime:` (resolved tier) lines on stdout, plus a
  `[W::bwa-mem3]` warning on stderr (exit 0) if the host is below the
  floor. See `docs/src/getting-started/host-requirements.md`.
* `bwa-mem3 shm` performs a `statvfs("/dev/shm")` capacity preflight
  (#86). When `/dev/shm` is too small for the index, the stage aborts
  with an `[E::bwa_shm_stage]` message naming `/dev/shm`, the required
  size, and a `mount -o remount,size=...` hint — replacing the prior
  `[fread] Bad address` failure mode. `statvfs` failures (no `/dev/shm`,
  restricted sandbox) are non-fatal and the stage proceeds.
* `bwa-mem3 shm` `/bwactl` registry RMW is now serialized via a POSIX
  named semaphore (#82, closes #66). Concurrent `shm stage` / `shm drop`
  invocations across processes no longer race when updating the
  registry; the prior best-effort flock was per-`open` and did not
  cover the read-modify-write window.

Methylation

* `mem --meth` emits Bismark-compatible auxiliary tags `XR:Z`
  (read conversion `CT`/`GA`), `XG:Z` (genome strand `CT`/`GA`),
  and `XM:Z` (per-base methylation call string) (#90). These replace
  the prior bwameth-style `YS:Z` / `YC:Z` / `YD:Z` on output (still
  used internally for SEQ restoration). The reference-annotation
  `XR:Z` from `-V` is suppressed under `--meth` to avoid colliding
  with the Bismark semantics. Downstream tools that previously read
  `YS:Z` / `YC:Z` / `YD:Z` must be pointed at the corresponding
  `XR:Z` / `XG:Z` and the per-base `XM:Z`. See
  `docs/src/methylation/tags.md`.

Correctness

* Fixed SIGSEGV in `mem_matesw` on shm-backed `ref_string` (#85).
  `ksw_align2` mutates its reference slice in place; when the slice
  pointed into a read-only shm segment, this faulted. Now copies the
  slice before passing it in.
* `FMI_search` sampled-SA prefetch: parenthesized `SA_COMPX_MASK`
  precedence so the masked offset is computed against the correct
  operand (#73). The unparenthesized form was silently producing
  wrong-but-harmless prefetch addresses; no alignment output was
  affected.
* `bntseq` `.alt` parser bounds the line buffer to prevent a
  stack-overflow on malicious or malformed `.alt` files (#74).
* `display_stats` clamps the per-thread bucket count to `LIM_C` so
  `--profile` with `-t` greater than the compiled-in limit no longer
  writes past the end of the stats array (#81).

Performance

* x86 wall-time improvements on the bench (vs the 0.1.0-pre baseline):
  AVX2 (c6a) −17 to −22%, AVX-512 AMD Zen4 (c7a) −16 to −24%, AVX-512
  Intel SPR (c7i) −28 to −30% across wgs / wes / panel-twist 5M-read
  samples. Concordance vs upstream `bwa-mem2 v2.2.1` remains
  100.0000% on all non-methylation cells. arm64 (c7g / c8g) is flat
  (within ±2%). The wins are attributable primarily to (a) capping
  AVX-512BW auto-vectorization at 256-bit on the `avx512bw` target
  (#86) and (b) inlining `FMI_search::backwardExt` to recover a
  gcc 12+ wall-clock regression (#88). See
  `docs/src/performance/overview.md` for the reference numbers across
  architectures.
* Smaller contributions in the release window: per-strip L1 prefetches
  across all `kswv` u8/u16 kernels (#70); `SMEM_LOCKSTEP_N` bumped from
  8 to 16 (#75); closed-form ungapped HIT path when `total_mis == 0`
  (#77); `ksort` switched to an on-stack buffer for small `n` to drop
  a per-call `malloc` (#78); `libsais_build` skips a wasted zero-init
  pass on its unpack and SA buffers, trimming index-build time (#80).

Release 0.1.0-pre (2026-04-28)
---------------------------------

* Project renamed from `bwa-mem2` to `bwa-mem3`. The new project tracks
  Fulcrum Genomics' performance and feature work on top of the upstream
  bwa-mem2 codebase.
* Default branch renamed from `fg-main` to `main`.
* Binary renamed from `bwa-mem2` to `bwa-mem3`. Arch-suffixed variants
  (`bwa-mem3.sse41`, `.sse42`, `.avx`, `.avx2`, `.avx512bw`, `.arm64`,
  `.pgo`, `.profile`, `.lto`) renamed to match.
* `@PG` SAM header tags now read `ID:bwa-mem3 PN:bwa-mem3` (and
  `bwa-mem3-meth` for `--meth` mode).
* Test binaries renamed: `bwa_mem2_tests_unit` →
  `bwa_mem3_tests_unit`, `bwa_mem2_tests_integration` →
  `bwa_mem3_tests_integration`.
* `.bwt.2bit.64` index file format unchanged — bwa-mem3 reads indexes
  built by `bwa-mem2 index` without re-indexing.

Release 2.2.1 (17 March 2021)
---------------------------------
Hotfix for v2.2: Fixed the bug mentioned in #135.


Release 2.2 (8 March 2021)
---------------------------------
Changes since the last release (2.1):

* Passed the validation test on ~88 billions reads (Credits: Keiran Raine, CASM division, Sanger Institute)
* Fixed bugs reported in #109 causing mismatch between bwa-mem and bwa-mem2
* Fixed the issue (# 112) causing crash due to corrupted thread id 
* Using all the SSE flags to create optimized SSE41 and SSE42 binaries


Release 2.1 (16 October 2020)
---------------------------------
Release 2.1 of BWA-MEM2.

Changes since the last release (2.0):
* *Smaller index*: the index size on disk is down by 8 times and in memory by 4 times due to moving to only one type of FM-index (2bit.64 instead of 2bit.64 and 8bit.32) and 8x compression of suffix array. For example, for human genome, index size on disk is down to ~10GB from ~80GB and memory footprint is down to ~10GB from ~40GB. There is a substantial decrease in index IO time due to the reduction and hardly any performance impact on read mapping.

* Added support for 2 more execution modes: sse4.2 and avx.

* Fixed multiple bugs including those reported in Issues #71, #80 and #85.

* Merged multiple pull requests.


Release 2.0 (9 July 2020)
---------------------------------
This is the first production release of BWA-MEM2.

Changes since the last release:
* Made the source code more secure with more than 300 changes all across it.

* Added support for memory re-allocations in case the pre-allocated fixed memory is insufficient.

* Added support for MC flag in the sam file and support for -5, -q flags in the command line.

* The output is now identical to the output of bwa-mem-0.7.17.

* Merged index building code with FMI_Search class.

* Added support for different ways to input read files, now, it is same as bwa-mem.

* Fixed a bug in AVX512 sam processing part, which was leading to incorrect output.


Release 2.0pre2 (4 February 2020)
---------------------------------

Miscellaneous changes:

 * Changed the license from GPL to MIT.

 * IMPORTANT: the index structure has changed since commit 6743183. Please
   rebuild the index if you are using a later commit or the new release.

 * Added charts in README.md comparing the performance of bwa-mem2 with bwa-mem.

Major code changes:

 * Fixed working for variable length reads.

 * Fixed a bug involving reads of length greater than 250bp.

 * Added support for allocation of more memory in small chunks if large
   pre-allocated fixed memory is insufficient. This is needed very rarely
   (thus, having no impact on performance) but prevents asserts from failing
   (code from crashing) in that scenario. 

 * Fixed a memory leak due to not releasing the memory allocated for seeds
   after smem.

 * Fixed a segfault due to non-alignment of small allocated memory in the
   optimized banded Smith-Waterman.

 * Enabled working with genomes larger than 7-8 billion nucleotides (e.g. Wheat
   genome).

 * Fixed a segfault occuring (with gcc compiler) while reading the index.
