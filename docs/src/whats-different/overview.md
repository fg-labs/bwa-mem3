# What's Different from bwa-mem2

This section tracks every change that bwa-mem3 carries on top of upstream
[bwa-mem2/bwa-mem2](https://github.com/bwa-mem2/bwa-mem2)'s `master` branch,
explains why each change was made, and records its upstream disposition.

## How this section is organized

Each deep-dive page covers one category of change:

- [Correctness fixes](correctness.md) — bugs in upstream bwa-mem2 that are
  fixed in bwa-mem3, including the kswv SIMD score2 plateau series, the
  proper-pair flag regression, the zero-init crash, the SMEM buffer overflow,
  and the `@PG` tab-escape issue.
- [Performance improvements](performance.md) — lockstep SMEM batching, batched
  `-H` header ingestion, libsais FM-index construction, and the consolidated
  mapping speedup suite.
- [Features](features.md) — `--meth` bisulfite mode, mimalloc allocator,
  `--supp-rep-hard-cap`, `bwa-mem3 shm`, `shm --meth`, the `HN:i` tag, and
  the `--bam=LEVEL` output flag.
- [Architecture support](arch-support.md) — the Linux ARM64/aarch64 build,
  the `arch=avx512bw` Makefile target, the NEON kswv mate-rescue kernel, and
  the AVX2 kswv mate-rescue kernel.
- [Build & infrastructure](build-infra.md) — the doctest framework, Codecov
  integration, `PACKAGE_VERSION` from `git describe`, PGO target
  parameterization, `CXXFLAGS`/`CPPFLAGS`/`LDFLAGS` forwarding, the unit-test
  harness, and the CI matrix expansion.
- [Upstream PR status](upstream-prs.md) — a single table cross-referencing
  every fork-carried change to its corresponding upstream PR or issue, with
  current upstream disposition.

## Carried on top of upstream

<!-- FG-MAIN-TABLE:start -->

Auto-generated from `git log --reverse --no-merges master..main` and the conventional-commits scope on each PR-merge title; do not edit by hand. For per-PR upstream disposition (bwa-mem2 PR / issue refs and status), see [Upstream PR status](upstream-prs.md).

| Commit    | Topic                                                                                                      | PR |
|-----------|------------------------------------------------------------------------------------------------------------|----|
| `ae73227` | Add Apple Silicon (ARM64/NEON) support with native optimizations | — |
| `744a9e7` | feat: add CI workflow with cross-platform build and end-to-end test | — |
| `490502b` | fix: drop unused global `stat` that shadows libc | — |
| `9364cfc` | ci: pin GitHub Actions to full-length commit SHAs | [#4](https://github.com/fg-labs/bwa-mem3/pull/4) |
| `b6eaba1` | chore: configure CodeRabbit to review PRs against fg-main | [#2](https://github.com/fg-labs/bwa-mem3/pull/2) |
| `db5086a` | docs: add FG-MAIN.md documenting the fork's relationship to upstream | [#3](https://github.com/fg-labs/bwa-mem3/pull/3) |
| `5132582` | feat(arm64): make Linux aarch64 build + CI-test on every fg-main push | [#1](https://github.com/fg-labs/bwa-mem3/pull/1) |
| `96016a5` | ci: pin dwgsim seed (-z 42) to stop parity-test flakiness | [#10](https://github.com/fg-labs/bwa-mem3/pull/10) |
| `246b528` | fix(hdr): align bwamem.h declarations with bwamem.cpp definitions | [#5](https://github.com/fg-labs/bwa-mem3/pull/5) |
| `b27f374` | feat(hdr): export mem_infer_dir for external consumers | [#6](https://github.com/fg-labs/bwa-mem3/pull/6) |
| `62700b1` | chore: move profiling globals out of main.cpp | [#7](https://github.com/fg-labs/bwa-mem3/pull/7) |
| `6b76c7b` | feat: expose worker_alloc/worker_free, the core worker_t pre-allocation helpers | [#8](https://github.com/fg-labs/bwa-mem3/pull/8) |
| `e80765b` | feat: split mem_sam_pe into mem_pair_resolve + thin emission wrapper | [#9](https://github.com/fg-labs/bwa-mem3/pull/9) |
| `84defc3` | feat: --bam[=LEVEL] output flag for direct BAM emission | [#12](https://github.com/fg-labs/bwa-mem3/pull/12) |
| `73907d7` | feat: vendor mimalloc v3.3.0 and link by default | [#19](https://github.com/fg-labs/bwa-mem3/pull/19) |
| `7641ebf` | feat(meth): --meth + `index --meth` — bwameth.py-equivalent bisulfite mode | [#13](https://github.com/fg-labs/bwa-mem3/pull/13) |
| `0165b6c` | fix: zero bseq1_t in kseq2bseq1 so realloc'd entries don't carry garbage | [#22](https://github.com/fg-labs/bwa-mem3/pull/22) |
| `e7cb763` | [proto] NEON kswv mate-rescue — correctness + perf harness | [#18](https://github.com/fg-labs/bwa-mem3/pull/18) |
| `a5aab04` | test(ci): add unit-test harness, fixtures, and ARM build support | [#23](https://github.com/fg-labs/bwa-mem3/pull/23) |
| `2fddafd` | [proto] AVX2 kswv mate-rescue — stacked on PR 18 | [#20](https://github.com/fg-labs/bwa-mem3/pull/20) |
| `8944028` | fix: compute no_pairing 0x2 flag from the emitted alignment | [#17](https://github.com/fg-labs/bwa-mem3/pull/17) |
| `2fd0e96` | fix(kswv): apply NEON score2-scan fixes to AVX-512BW kernel | [#21](https://github.com/fg-labs/bwa-mem3/pull/21) |
| `68adecd` | ci: expand workflow matrix + add canonical deep-test row | [#24](https://github.com/fg-labs/bwa-mem3/pull/24) |
| `690914f` | build(make): add explicit arch=avx512bw target | [#16](https://github.com/fg-labs/bwa-mem3/pull/16) |
| `0bb9402` | fix(kswv): gate AVX2 arch dispatch on !__AVX512BW__ | [#26](https://github.com/fg-labs/bwa-mem3/pull/26) |
| `43457e8` | fix(kswv): consolidate score2 plateaus per-lane to match scalar ksw_align2 | [#28](https://github.com/fg-labs/bwa-mem3/pull/28) |
| `2311f11` | fix(kswv): port score2 plateau consolidation to NEON + AVX-512BW | [#29](https://github.com/fg-labs/bwa-mem3/pull/29) |
| `75c709a` | fix(kswv): apply score2 plateau fix + missing filters to kswv_512_16 | [#30](https://github.com/fg-labs/bwa-mem3/pull/30) |
| `61813ef` | fix(kswv): rewrite kswv_neon_16 — real SIMD kernel with correct table + score2 | [#31](https://github.com/fg-labs/bwa-mem3/pull/31) |
| `1f76655` | perf(seed): lockstep SMEM batching across N reads | [#33](https://github.com/fg-labs/bwa-mem3/pull/33) |
| `93a79ec` | feat(mem): emit HN:i tag with total hit count per primary | [#42](https://github.com/fg-labs/bwa-mem3/pull/42) |
| `dd3a82c` | chore: port four nh13 lh3/bwa PRs into bwa-mem2 (-z, -u/XB, MQ, @HD order) | [#35](https://github.com/fg-labs/bwa-mem3/pull/35) |
| `98ba6ab` | build(make): forward user CXXFLAGS/CPPFLAGS/LDFLAGS to final link steps | [#50](https://github.com/fg-labs/bwa-mem3/pull/50) |
| `e9302a1` | fix(kswv): guard post-loop rowMax store on nrow==0 batches | [#51](https://github.com/fg-labs/bwa-mem3/pull/51) |
| `9b702ca` | fix(sam): sanitize whitespace in -R when embedding into @PG CL: field | [#54](https://github.com/fg-labs/bwa-mem3/pull/54) |
| `ed63fad` | perf(header): batch -H ingestion to fix O(n^2) header read (closes #37) | [#49](https://github.com/fg-labs/bwa-mem3/pull/49) |
| `595d8e5` | feat(mapq): add --supp-rep-hard-cap opt-in supp MAPQ rescoring | [#56](https://github.com/fg-labs/bwa-mem3/pull/56) |
| `79628c3` | chore(version): stamp PACKAGE_VERSION from git describe at build time | [#52](https://github.com/fg-labs/bwa-mem3/pull/52) |
| `e22dade` | fix(smem): size SMEM buffers from observed max read length (closes #44) | [#55](https://github.com/fg-labs/bwa-mem3/pull/55) |
| `03688a0` | chore: normalize CRLF line endings to LF (#43) | [#53](https://github.com/fg-labs/bwa-mem3/pull/53) |
| `57e21bd` | feat(makefile): parameterize PGO targets by arch + profile dir | [#59](https://github.com/fg-labs/bwa-mem3/pull/59) |
| `79b90ce` | feat(index): libsais-based memory-bounded FM-index construction | [#57](https://github.com/fg-labs/bwa-mem3/pull/57) |
| `d8d4a6d` | feat(cli): wire up --help across commands; add -h to top-level and index | [#60](https://github.com/fg-labs/bwa-mem3/pull/60) |
| `7301762` | perf: consolidated mapping speedups (ksw2, SMEM, SAL, SAM) | [#58](https://github.com/fg-labs/bwa-mem3/pull/58) |
| `eaf4ed6` | test: doctest-based test framework scaffolding + Codecov | [#34](https://github.com/fg-labs/bwa-mem3/pull/34) |
| `bbbecd3` | ci(proto-neon-kswv): split into fan-out/fan-in jobs with caching | [#63](https://github.com/fg-labs/bwa-mem3/pull/63) |
| `20f77e9` | feat(shm): port `bwa shm` from bwa-mem v1 | [#65](https://github.com/fg-labs/bwa-mem3/pull/65) |
| `c20f61c` | feat(shm): add `bwa-mem2 shm --meth` for symmetric meth UX | [#67](https://github.com/fg-labs/bwa-mem3/pull/67) |
| `ee18a3b` | refactor: rename bwa_mem2idx to bwa_mem3idx | — |
| `bb919f2` | feat: rename PG header to bwa-mem3 (ID, PN, usage strings) | — |
| `7a56f9a` | feat: rename meth PG header to bwa-mem3-meth and drive VN from PACKAGE_VERSION | — |
| `95c673d` | build: rename binary to bwa-mem3, update version guard and fallback | — |
| `2d5ad10` | test: rename test binaries from bwa_mem2_tests_* to bwa_mem3_tests_* | — |
| `31f214a` | chore: sweep bwa-mem2 -> bwa-mem3 in source comments and log messages | — |
| `ff40f96` | chore: rename BWAMEM2_* header guards to BWAMEM3_* | — |
| `c2c786a` | test: sweep bwa-mem2 -> bwa-mem3 in test, bench, and scripts | — |
| `148e431` | ci: update workflows for bwa-mem3 rename and main branch | — |
| `dddd8dd` | docs: rewrite README for bwa-mem3 (lineage attribution, drop upstream-only sections) | — |
| `34c8ea3` | docs: rename FG-MAIN.md to docs/whats-different.md | — |
| `5719617` | docs: update whats-different.md for bwa-mem3 and main branch | — |
| `bdd67f3` | docs: drop README-ori.md (lineage preserved in README + git history) | — |
| `924a70f` | docs: add 0.1.0-pre release notes and update status.md | — |
| `85d3b3b` | ci: drop master from branch filter (master branch removed from remote) | — |
| `2ea69db` | fix(test/meth): alias bwa-mem2 -> bwa-mem3 on PATH for bwameth.py oracle | [#72](https://github.com/fg-labs/bwa-mem3/pull/72) |
| `4f805e6` | chore: rename shell vars BWAMEM2/BWA_MEM2[_*] to BWAMEM3/BWA_MEM3[_*] | [#68](https://github.com/fg-labs/bwa-mem3/pull/68) |
| `8137740` | perf(kswv): add per-strip L1 prefetches to all u8/16 kernels | [#70](https://github.com/fg-labs/bwa-mem3/pull/70) |
| `41e1f3c` | docs: add comprehensive mdbook on Read the Docs | [#71](https://github.com/fg-labs/bwa-mem3/pull/71) |
| `442de25` | fix(fmi): parenthesize SA_COMPX_MASK precedence in sampled-SA prefetch | [#73](https://github.com/fg-labs/bwa-mem3/pull/73) |
| `000c0fd` | perf(fmi): bump SMEM_LOCKSTEP_N from 8 to 16 | [#75](https://github.com/fg-labs/bwa-mem3/pull/75) |
| `b3a665e` | fix(bntseq): bound .alt parse buffer to prevent stack overflow | [#74](https://github.com/fg-labs/bwa-mem3/pull/74) |
| `af33cdd` | feat(bns): convert mem_matesw_batch_{pre,post} to bns_fetch_seq_v2 | [#76](https://github.com/fg-labs/bwa-mem3/pull/76) |
| `9bb277a` | Update index.md | [#79](https://github.com/fg-labs/bwa-mem3/pull/79) |
| `fdb244d` | perf(libsais_build): skip wasted zero-init on unpack + SA buffers | [#80](https://github.com/fg-labs/bwa-mem3/pull/80) |
| `ff95a4f` | perf(ksort): replace per-call malloc with on-stack buffer for small n | [#78](https://github.com/fg-labs/bwa-mem3/pull/78) |
| `7caf77c` | perf(ungapped): closed-form HIT for total_mis == 0 | [#77](https://github.com/fg-labs/bwa-mem3/pull/77) |
| `e65ceb2` | fix(profiling): clamp display_stats nthreads to LIM_C | [#81](https://github.com/fg-labs/bwa-mem3/pull/81) |
| `ddfb0da` | feat(shm): serialize /bwactl RMW with a POSIX named semaphore | [#82](https://github.com/fg-labs/bwa-mem3/pull/82) |
| `b9e0b66` | feat(simd): replace multi-binary execv launcher with single-binary in-process dispatch | [#83](https://github.com/fg-labs/bwa-mem3/pull/83) |
| `7d27f23` | perf(build): default x86 single-binary baseline to avx2 (was sse41) | [#84](https://github.com/fg-labs/bwa-mem3/pull/84) |
| `316dba6` | fix(matesw): copy ref slice before ksw_align2 to avoid SIGSEGV on shm-backed ref_string | [#85](https://github.com/fg-labs/bwa-mem3/pull/85) |
| `427c81c` | perf(fmi): inline backwardExt to recover gcc 12+ wall-clock regression | [#88](https://github.com/fg-labs/bwa-mem3/pull/88) |
| `c96d31a` | perf(x86): cap avx512bw autovec at 256-bit; bwa_shm /dev/shm preflight | [#86](https://github.com/fg-labs/bwa-mem3/pull/86) |
| `23f528d` | ci: migrate parity tests from dwgsim/phiX174 to holodeck/chr22 | [#89](https://github.com/fg-labs/bwa-mem3/pull/89) |
| `ec67b09` | feat(meth): emit Bismark-compatible XR/XG/XM auxiliary tags | [#90](https://github.com/fg-labs/bwa-mem3/pull/90) |
| `652ce0f` | docs(install): list autoconf/automake/libomp/zlib system prereqs | [#93](https://github.com/fg-labs/bwa-mem3/pull/93) |
| `296b1b9` | docs(install): fix RHEL/Fedora package name pkgconfig → pkgconf-pkg-config | [#94](https://github.com/fg-labs/bwa-mem3/pull/94) |
| `dc7fcfe` | feat(simd): add SIMD host-floor precheck for multi-arch deployment | [#95](https://github.com/fg-labs/bwa-mem3/pull/95) |
| `3bc64b0` | docs: pre-release documentation pass for v0.2.0-pre | [#96](https://github.com/fg-labs/bwa-mem3/pull/96) |
| `27a60c9` | chore(release): prep v0.2.0 release notes and metadata | [#97](https://github.com/fg-labs/bwa-mem3/pull/97) |

<!-- FG-MAIN-TABLE:end -->

### Additional fork-level changes

- **Vendored mimalloc allocator**: `ext/mimalloc` is pinned at `v3.3.0` and
  linked into every binary by default (`USE_MIMALLOC=1`). Linux uses
  `--whole-archive` static linkage; macOS uses dyld-interposed shared linkage.
  `USE_MIMALLOC=1` is the supported and recommended default on all
  platforms; `USE_MIMALLOC=0` is provided as a best-effort opt-out and is
  CI-gated on Linux x86 only. See [Features](features.md) for details.

- **`--supp-rep-hard-cap INT`** (opt-in, default disabled): forces MAPQ=0 on
  supplementary alignments whose chain contains a seed with `>=INT` genome
  occurrences. Addresses the long-standing bwa/bwa-mem2 issue where a supp
  fragment that maps to many places standalone (e.g. a short read in a CCATCC
  repeat) inherits a high MAPQ from its primary because the supp's competing
  repetitive chains get filtered out during the full-read pipeline and
  therefore never contribute to its `sub`/`sub_n`. See
  [upstream #260](https://github.com/bwa-mem2/bwa-mem2/issues/260) for the
  reporter case. Primary MAPQ is unaffected; default output is byte-identical
  to stock bwa-mem2. Typical values are 5–20 (lower = more aggressive); the
  upstream #260 repro drops from MAPQ=60 to MAPQ=0 at `--supp-rep-hard-cap 18`.


## Version stamping

`PACKAGE_VERSION` (the value reported by `bwa-mem3 version` and written to
the `@PG VN:` SAM header field) is generated at build time by the Makefile
from `git describe --tags --dirty`, e.g. `v2.3-30-g61813ef` for a tree 30
commits past upstream tag `v2.3` at commit `61813ef`.

- No manual bumping required: cut a fresh release by tagging the commit
  (`git tag -a vX.Y-fg-labs.N -m ...`) and the next build picks it up.
- Builds where `git describe --tags` fails (source-tarball extractions, or
  shallow clones / checkouts with no tag reachable from `HEAD` — including
  CI's default `actions/checkout` fetch-depth of 1) fall back to the static
  `FG_LABS_VERSION_FALLBACK` in `Makefile`. Bump that when cutting a
  release that will be consumed as a tarball, or in CI artifacts.
- `src/version.h` is generated and `.gitignore`d; `make clean` removes it.

## Branching and update policy

- `master` tracks upstream unchanged.
- `main` is `upstream/master` plus the commits above. Rebased onto upstream roughly quarterly, or sooner when an upstream release we care about lands.
- Contributions go via PR targeting `main`. CI and CodeRabbit gate merges.
- Any PR that adds or removes a fork-carried commit must update the table above in the same PR.

## Consuming

Clone this repo and check out `main`:

```bash
git clone https://github.com/fg-labs/bwa-mem3.git
cd bwa-mem3
git checkout main
```

Or vendor the branch into a downstream repo by pinning to a specific
commit (not the branch tip) so your build is reproducible.

## Relationship to upstream

We submit the generally-useful fixes and features carried here as PRs against
[bwa-mem2/bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) when the upstream
maintainers are actively merging; while they are not, fixes land here first
and we drop them from `main` once they appear upstream.

---

**See also:**
[Correctness fixes](correctness.md) ·
[Performance improvements](performance.md) ·
[Features](features.md) ·
[Upstream PR status](upstream-prs.md) ·
[Developer Guide → Contributing](../developer-guide/contributing.md)
