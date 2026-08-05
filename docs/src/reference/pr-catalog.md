# PR catalog

The single, consolidated record of every change carried in bwa-mem3 `main` on
top of upstream bwa-mem2 — its bwa-mem3 PR, a one-line description, its class,
and its upstream bwa-mem2 disposition (PR/issue and status). The narrative
*What's Different* pages (correctness, performance, features, architecture
support, build & infrastructure) explain the *why*; this page is the flat,
scannable *what*.

> **Note — generated table; do not hand-edit between the markers**
>
> The table below is consolidated from the per-area change history
> (`git log --reverse --no-merges master..main` for the change list, merged with
> the upstream-disposition tracking). Every PR that adds a fork-carried commit
> must add its row here (the `FG-MAIN-TABLE` rule — see
> [Contributing](../developer-guide/contributing.md#the-fg-main-table-rule)).
> "Upstream status" of *fork-only* means no upstream PR exists; *open* means an
> upstream PR/issue existed at implementation time but was unmerged.

## Fork-carried changes

<!-- FG-MAIN-TABLE:start -->
| PR | Change | Class | Upstream PR / issue | Upstream status |
|----|--------|-------|---------------------|-----------------|
| [#1](https://github.com/fg-labs/bwa-mem3/pull/1) | feat(arm64): make Linux aarch64 build + CI-test on every fg-main push | Architecture support | [bwa-mem2#288](https://github.com/bwa-mem2/bwa-mem2/pull/288) | open PR |
| [#2](https://github.com/fg-labs/bwa-mem3/pull/2) | chore: configure CodeRabbit to review PRs against fg-main | Build & infrastructure | — | fork-only |
| [#3](https://github.com/fg-labs/bwa-mem3/pull/3) | docs: add FG-MAIN.md documenting the fork's relationship to upstream | Build & infrastructure | — | fork-only |
| [#4](https://github.com/fg-labs/bwa-mem3/pull/4) | ci: pin GitHub Actions to full-length commit SHAs | Build & infrastructure | — | fork-only |
| [#5](https://github.com/fg-labs/bwa-mem3/pull/5) | fix(hdr): align bwamem.h declarations with bwamem.cpp definitions | Correctness | — | fork-only |
| [#6](https://github.com/fg-labs/bwa-mem3/pull/6) | feat(hdr): export mem_infer_dir for external consumers | Features | — | fork-only |
| [#7](https://github.com/fg-labs/bwa-mem3/pull/7) | chore: move profiling globals out of main.cpp | Build & infrastructure | — | fork-only |
| [#8](https://github.com/fg-labs/bwa-mem3/pull/8) | feat: expose worker_alloc/worker_free, the core worker_t pre-allocation helpers | Features | — | fork-only |
| [#9](https://github.com/fg-labs/bwa-mem3/pull/9) | feat: split mem_sam_pe into mem_pair_resolve + thin emission wrapper | Features | — | fork-only |
| [#10](https://github.com/fg-labs/bwa-mem3/pull/10) | ci: pin dwgsim seed (-z 42) to stop parity-test flakiness | Build & infrastructure | — | fork-only |
| [#12](https://github.com/fg-labs/bwa-mem3/pull/12) | feat: --bam[=LEVEL] output flag for direct BAM emission | Features | — | fork-only |
| [#13](https://github.com/fg-labs/bwa-mem3/pull/13) | feat(meth): --meth + `index --meth` — bwameth.py-equivalent bisulfite mode | Features | — | fork-only |
| [#16](https://github.com/fg-labs/bwa-mem3/pull/16) | build(make): add explicit arch=avx512bw target | Architecture support | — | fork-only |
| [#17](https://github.com/fg-labs/bwa-mem3/pull/17) | fix: compute no_pairing 0x2 flag from the emitted alignment | Correctness | — | fork-only |
| [#18](https://github.com/fg-labs/bwa-mem3/pull/18) | [proto] NEON kswv mate-rescue — correctness + perf harness | Architecture support | — | fork-only |
| [#19](https://github.com/fg-labs/bwa-mem3/pull/19) | feat: vendor mimalloc v3.3.0 and link by default | Features | — | fork-only |
| [#20](https://github.com/fg-labs/bwa-mem3/pull/20) | [proto] AVX2 kswv mate-rescue — stacked on PR 18 | Architecture support | — | fork-only |
| [#21](https://github.com/fg-labs/bwa-mem3/pull/21) | fix(kswv): apply NEON score2-scan fixes to AVX-512BW kernel | Correctness | — | fork-only |
| [#22](https://github.com/fg-labs/bwa-mem3/pull/22) | fix: zero bseq1_t in kseq2bseq1 so realloc'd entries don't carry garbage | Correctness | — | fork-only |
| [#23](https://github.com/fg-labs/bwa-mem3/pull/23) | test(ci): add unit-test harness, fixtures, and ARM build support | Build & infrastructure | — | fork-only |
| [#24](https://github.com/fg-labs/bwa-mem3/pull/24) | ci: expand workflow matrix + add canonical deep-test row | Build & infrastructure | — | fork-only |
| [#26](https://github.com/fg-labs/bwa-mem3/pull/26) | fix(kswv): gate AVX2 arch dispatch on !__AVX512BW__ | Correctness | — | fork-only |
| [#28](https://github.com/fg-labs/bwa-mem3/pull/28) | fix(kswv): consolidate score2 plateaus per-lane to match scalar ksw_align2 | Correctness | — | fork-only |
| [#29](https://github.com/fg-labs/bwa-mem3/pull/29) | fix(kswv): port score2 plateau consolidation to NEON + AVX-512BW | Correctness | — | fork-only |
| [#30](https://github.com/fg-labs/bwa-mem3/pull/30) | fix(kswv): apply score2 plateau fix + missing filters to kswv_512_16 | Correctness | — | fork-only |
| [#31](https://github.com/fg-labs/bwa-mem3/pull/31) | fix(kswv): rewrite kswv_neon_16 — real SIMD kernel with correct table + score2 | Correctness | — | fork-only |
| [#33](https://github.com/fg-labs/bwa-mem3/pull/33) | perf(seed): lockstep SMEM batching across N reads | Performance | — | fork-only |
| [#34](https://github.com/fg-labs/bwa-mem3/pull/34) | test: doctest-based test framework scaffolding + Codecov | Build & infrastructure | — | fork-only |
| [#35](https://github.com/fg-labs/bwa-mem3/pull/35) | chore: port four nh13 lh3/bwa PRs into bwa-mem2 (-z, -u/XB, MQ, `@HD` order) | Features | [lh3/bwa#330](https://github.com/lh3/bwa/pull/330) | merged (lh3 only) |
| [#42](https://github.com/fg-labs/bwa-mem3/pull/42) | feat(mem): emit HN:i tag with total hit count per primary | Features | [lh3/bwa#438](https://github.com/lh3/bwa/pull/438) | analogous to bwa aln; no direct upstream port |
| [#49](https://github.com/fg-labs/bwa-mem3/pull/49) | perf(header): batch -H ingestion to fix O(n^2) header read (closes #37) | Performance | [bwa-mem2#204](https://github.com/bwa-mem2/bwa-mem2/pull/204) | open PR |
| [#50](https://github.com/fg-labs/bwa-mem3/pull/50) | build(make): forward user CXXFLAGS/CPPFLAGS/LDFLAGS to final link steps | Build & infrastructure | [bwa-mem2#290](https://github.com/bwa-mem2/bwa-mem2/pull/290) | open upstream PR |
| [#51](https://github.com/fg-labs/bwa-mem3/pull/51) | fix(kswv): guard post-loop rowMax store on nrow==0 batches | Correctness | [bwa-mem2#289](https://github.com/bwa-mem2/bwa-mem2/pull/289) | open PR (upstream covers AVX-512BW only) |
| [#52](https://github.com/fg-labs/bwa-mem3/pull/52) | chore(version): stamp PACKAGE_VERSION from git describe at build time | Build & infrastructure | [bwa-mem2#283](https://github.com/bwa-mem2/bwa-mem2/issues/283), [bwa-mem2#284](https://github.com/bwa-mem2/bwa-mem2/pull/284) | open issue + open PR |
| [#53](https://github.com/fg-labs/bwa-mem3/pull/53) | chore: normalize CRLF line endings to LF (#43) | Build & infrastructure | — | fork-only |
| [#54](https://github.com/fg-labs/bwa-mem3/pull/54) | fix(sam): sanitize whitespace in -R when embedding into `@PG` CL: field | Correctness | [bwa-mem2#293](https://github.com/bwa-mem2/bwa-mem2/issues/293) | open issue |
| [#55](https://github.com/fg-labs/bwa-mem3/pull/55) | fix(smem): size SMEM buffers from observed max read length (closes #44) | Correctness | [bwa-mem2#238](https://github.com/bwa-mem2/bwa-mem2/pull/238), [bwa-mem2#210](https://github.com/bwa-mem2/bwa-mem2/issues/210) | PR closed without merge; issue open |
| [#56](https://github.com/fg-labs/bwa-mem3/pull/56) | feat(mapq): add --supp-rep-hard-cap opt-in supp MAPQ rescoring | Features | [bwa-mem2#260](https://github.com/bwa-mem2/bwa-mem2/issues/260) | open issue |
| [#57](https://github.com/fg-labs/bwa-mem3/pull/57) | feat(index): libsais-based memory-bounded FM-index construction | Performance | — | fork-only |
| [#58](https://github.com/fg-labs/bwa-mem3/pull/58) | perf: consolidated mapping speedups (ksw2, SMEM, SAL, SAM) | Performance | — | fork-only |
| [#59](https://github.com/fg-labs/bwa-mem3/pull/59) | feat(makefile): parameterize PGO targets by arch + profile dir | Build & infrastructure | — | fork-only |
| [#60](https://github.com/fg-labs/bwa-mem3/pull/60) | feat(cli): wire up --help across commands; add -h to top-level and index | Features | — | fork-only |
| [#63](https://github.com/fg-labs/bwa-mem3/pull/63) | ci(proto-neon-kswv): split into fan-out/fan-in jobs with caching | Build & infrastructure | — | fork-only |
| [#65](https://github.com/fg-labs/bwa-mem3/pull/65) | feat(shm): port `bwa shm` from bwa-mem v1 | Features | — | fork-only (v1 feature port) |
| [#67](https://github.com/fg-labs/bwa-mem3/pull/67) | feat(shm): add `bwa-mem2 shm --meth` for symmetric meth UX | Features | — | fork-only |
| [#68](https://github.com/fg-labs/bwa-mem3/pull/68) | chore: rename shell vars BWAMEM2/BWA_MEM2[_*] to BWAMEM3/BWA_MEM3[_*] | Build & infrastructure | — | fork-only |
| [#70](https://github.com/fg-labs/bwa-mem3/pull/70) | perf(kswv): add per-strip L1 prefetches to all u8/16 kernels | Performance | — | fork-only |
| [#71](https://github.com/fg-labs/bwa-mem3/pull/71) | docs: add comprehensive mdbook on Read the Docs | Build & infrastructure | — | fork-only |
| [#72](https://github.com/fg-labs/bwa-mem3/pull/72) | fix(test/meth): alias bwa-mem2 -> bwa-mem3 on PATH for bwameth.py oracle | Build & infrastructure | — | fork-only |
| [#73](https://github.com/fg-labs/bwa-mem3/pull/73) | fix(fmi): parenthesize SA_COMPX_MASK precedence in sampled-SA prefetch | Correctness | — | fork-only |
| [#74](https://github.com/fg-labs/bwa-mem3/pull/74) | fix(bntseq): bound .alt parse buffer to prevent stack overflow | Correctness | — | fork-only |
| [#75](https://github.com/fg-labs/bwa-mem3/pull/75) | perf(fmi): bump SMEM_LOCKSTEP_N from 8 to 16 | Performance | — | fork-only |
| [#76](https://github.com/fg-labs/bwa-mem3/pull/76) | feat(bns): convert mem_matesw_batch_{pre,post} to bns_fetch_seq_v2 | Architecture support | — | fork-only |
| [#77](https://github.com/fg-labs/bwa-mem3/pull/77) | perf(ungapped): closed-form HIT for total_mis == 0 | Performance | — | fork-only |
| [#78](https://github.com/fg-labs/bwa-mem3/pull/78) | perf(ksort): replace per-call malloc with on-stack buffer for small n | Performance | — | fork-only |
| [#79](https://github.com/fg-labs/bwa-mem3/pull/79) | Update index.md | Build & infrastructure | — | fork-only |
| [#80](https://github.com/fg-labs/bwa-mem3/pull/80) | perf(libsais_build): skip wasted zero-init on unpack + SA buffers | Performance | — | fork-only |
| [#81](https://github.com/fg-labs/bwa-mem3/pull/81) | fix(profiling): clamp display_stats nthreads to LIM_C | Correctness | — | fork-only |
| [#82](https://github.com/fg-labs/bwa-mem3/pull/82) | feat(shm): serialize /bwactl RMW with a POSIX named semaphore | Features | — | fork-only |
| [#83](https://github.com/fg-labs/bwa-mem3/pull/83) | feat(simd): replace multi-binary execv launcher with single-binary in-process dispatch | Architecture support | — | fork-only |
| [#84](https://github.com/fg-labs/bwa-mem3/pull/84) | perf(build): default x86 single-binary baseline to avx2 (was sse41) | Architecture support | — | fork-only |
| [#85](https://github.com/fg-labs/bwa-mem3/pull/85) | fix(matesw): copy ref slice before ksw_align2 to avoid SIGSEGV on shm-backed ref_string | Correctness | — | fork-only |
| [#86](https://github.com/fg-labs/bwa-mem3/pull/86) | perf(x86): cap avx512bw autovec at 256-bit; bwa_shm /dev/shm preflight | Features | — | fork-only |
| [#88](https://github.com/fg-labs/bwa-mem3/pull/88) | perf(fmi): inline backwardExt to recover gcc 12+ wall-clock regression | Performance | — | fork-only |
| [#89](https://github.com/fg-labs/bwa-mem3/pull/89) | ci: migrate parity tests from dwgsim/phiX174 to holodeck/chr22 | Build & infrastructure | — | fork-only |
| [#90](https://github.com/fg-labs/bwa-mem3/pull/90) | feat(meth): emit Bismark-compatible XR/XG/XM auxiliary tags | Features | — | fork-only |
| [#93](https://github.com/fg-labs/bwa-mem3/pull/93) | docs(install): list autoconf/automake/libomp/zlib system prereqs | Build & infrastructure | — | fork-only |
| [#94](https://github.com/fg-labs/bwa-mem3/pull/94) | docs(install): fix RHEL/Fedora package name pkgconfig → pkgconf-pkg-config | Build & infrastructure | — | fork-only |
| [#95](https://github.com/fg-labs/bwa-mem3/pull/95) | feat(simd): add SIMD host-floor precheck for multi-arch deployment | Features | — | fork-only |
| [#96](https://github.com/fg-labs/bwa-mem3/pull/96) | docs: pre-release documentation pass for v0.2.0-pre | Build & infrastructure | — | fork-only |
| [#97](https://github.com/fg-labs/bwa-mem3/pull/97) | chore(release): prep v0.2.0 release notes and metadata | Build & infrastructure | — | fork-only |
| [#123](https://github.com/fg-labs/bwa-mem3/pull/123) | Stable tie-breaks + pdqsort | Performance | — | fork-only |
| [#128](https://github.com/fg-labs/bwa-mem3/pull/128) | FASTQ reader fast path | Performance | — | fork-only |
| [#140](https://github.com/fg-labs/bwa-mem3/pull/140) | Recover 8-bit banded SW (≥128 bp) | Performance | — | fork-only |
| [#141](https://github.com/fg-labs/bwa-mem3/pull/141) | Gotoh gaps from H | Performance | — | fork-only |
| [#143](https://github.com/fg-labs/bwa-mem3/pull/143) | Drop dead `qlen[]` param | Performance | — | fork-only |
| [#144](https://github.com/fg-labs/bwa-mem3/pull/144) | Long-read kernel parity test | Performance | — | fork-only |
| [#147](https://github.com/fg-labs/bwa-mem3/pull/147) | Short-circuit re-baseline scan | Performance | — | fork-only |
| [#148](https://github.com/fg-labs/bwa-mem3/pull/148) | Remove dead SW code paths | Performance | — | fork-only |
| [#149](https://github.com/fg-labs/bwa-mem3/pull/149) | Vectorize epilogue side-channel | Performance | — | fork-only |
| [#150](https://github.com/fg-labs/bwa-mem3/pull/150) | Bound getScores prefetch reads | Performance | — | fork-only |
| [#151](https://github.com/fg-labs/bwa-mem3/pull/151) | Unsigned 8-bit h0-prefix seed | Performance | — | fork-only |
| [#152](https://github.com/fg-labs/bwa-mem3/pull/152) | `--profile` stage timing | Performance | — | fork-only |
| [#153](https://github.com/fg-labs/bwa-mem3/pull/153) | zlib-ng inflate + 3rd worker | Performance | — | fork-only |
| [#157](https://github.com/fg-labs/bwa-mem3/pull/157) | Right-size SA staging buffers | Performance | — | fork-only |
| [#158](https://github.com/fg-labs/bwa-mem3/pull/158) | `gtle` contract test | Performance | — | fork-only |
| [#160](https://github.com/fg-labs/bwa-mem3/pull/160) | NEON SW tuning | Performance | — | fork-only |
| [#161](https://github.com/fg-labs/bwa-mem3/pull/161) | AVX2 SW tuning | Performance | — | fork-only |
| [#162](https://github.com/fg-labs/bwa-mem3/pull/162) | AVX2 16-bit `kswv256_16` | Performance | — | fork-only |
| [#164](https://github.com/fg-labs/bwa-mem3/pull/164) | NEON `movemask` parity test | Performance | — | fork-only |
| [#363](https://github.com/fg-labs/bwa-mem3/pull/363) | fix(pair): derive FLAG `0x2` from `a[0]` by default, matching both upstreams; [#17](https://github.com/fg-labs/bwa-mem3/pull/17)'s emitted-alignment derivation moves behind `--proper-pair-from-emitted` ([#362](https://github.com/fg-labs/bwa-mem3/issues/362)) | Correctness | — | fork-only |
<!-- FG-MAIN-TABLE:end -->

## Upstream issues tracked but not yet fixed

These upstream issues are tracked in the bwa-mem3 issue list but do not yet have
a corresponding fix in `main`:

| Issue | Upstream reference | Notes |
|-------|--------------------|-------|
| Split-alignment evidence loss vs bwa 0.7.17 | [bwa-mem2#273](https://github.com/bwa-mem2/bwa-mem2/issues/273) | [issue #47](https://github.com/fg-labs/bwa-mem3/issues/47) — under investigation |
| MAPQ/coordinate parity vs bwa mem 0.7.18 | [bwa-mem2#262](https://github.com/bwa-mem2/bwa-mem2/issues/262), [bwa-mem2#246](https://github.com/bwa-mem2/bwa-mem2/issues/246), [bwa-mem2#239](https://github.com/bwa-mem2/bwa-mem2/issues/239) | [issue #48](https://github.com/fg-labs/bwa-mem3/issues/48) — tracking only |

---

**See also:**
[What's Different — Overview](../whats-different/overview.md) ·
[Correctness fixes](../whats-different/correctness.md) ·
[Performance improvements](../whats-different/performance.md) ·
[Features](../whats-different/features.md) ·
[Architecture support](../whats-different/arch-support.md) ·
[Build & infrastructure](../whats-different/build-infra.md)
