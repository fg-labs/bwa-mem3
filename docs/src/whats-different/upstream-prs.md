# Upstream PR Status

This table cross-references every change carried in bwa-mem3 `main` to its
corresponding upstream bwa-mem2 PR or issue. "Fork-only" means no upstream PR
exists; the change may be submitted upstream in the future or may be
fork-specific by design. "Open" means the upstream PR or issue existed at the
time of bwa-mem3's implementation but had not been merged. Upstream status is
current as of the bwa-mem3 0.2.0-pre release.

For prose descriptions of each change, follow the links in the "bwa-mem3 PR"
column to the relevant deep-dive page section.

## Full cross-reference table

| Topic | bwa-mem3 PR | Upstream PR / Issue | Upstream status |
|-------|-------------|---------------------|-----------------|
| **Correctness** | | | |
| `@PG CL:` tab escaping | [#54](https://github.com/fg-labs/bwa-mem3/pull/54) | [bwa-mem2#293](https://github.com/bwa-mem2/bwa-mem2/issues/293) | open issue |
| SMEM buffer overflow on >151 bp reads | [#55](https://github.com/fg-labs/bwa-mem3/pull/55) | [bwa-mem2#238](https://github.com/bwa-mem2/bwa-mem2/pull/238), [bwa-mem2#210](https://github.com/bwa-mem2/bwa-mem2/issues/210) | PR closed without merge; issue open |
| kswv nrow==0 guard (all 5 kernels) | [#51](https://github.com/fg-labs/bwa-mem3/pull/51) | [bwa-mem2#289](https://github.com/bwa-mem2/bwa-mem2/pull/289) | open PR (upstream covers AVX-512BW only) |
| AVX-512BW dispatch guard (`!__AVX512BW__`) | [#26](https://github.com/fg-labs/bwa-mem3/pull/26) | — | fork-only |
| AVX2 score2 plateau consolidation | [#28](https://github.com/fg-labs/bwa-mem3/pull/28) | — | fork-only |
| NEON + AVX-512BW 8-bit score2 fix | [#29](https://github.com/fg-labs/bwa-mem3/pull/29) | — | fork-only |
| AVX-512BW 16-bit score2 fix | [#30](https://github.com/fg-labs/bwa-mem3/pull/30) | — | fork-only |
| NEON 16-bit kernel rewrite | [#31](https://github.com/fg-labs/bwa-mem3/pull/31) | — | fork-only |
| kseq2bseq1 zero-initialization | [#22](https://github.com/fg-labs/bwa-mem3/pull/22) | — | fork-only |
| Proper-pair flag from emitted alignment | [#17](https://github.com/fg-labs/bwa-mem3/pull/17) | — | fork-only |
| `@HD` emitted before `@SQ` per SAM spec | [#35](https://github.com/fg-labs/bwa-mem3/pull/35) | [lh3/bwa#345](https://github.com/lh3/bwa/pull/345) | closed (lh3 only) |
| `mem_matesw` SIGSEGV on shm-backed `ref_string` | [#85](https://github.com/fg-labs/bwa-mem3/pull/85) | — | fork-only |
| **Performance** | | | |
| Lockstep SMEM batching | [#33](https://github.com/fg-labs/bwa-mem3/pull/33) | — | fork-only |
| Batched `-H` header ingestion (O(n) fix) | [#49](https://github.com/fg-labs/bwa-mem3/pull/49) | [bwa-mem2#204](https://github.com/bwa-mem2/bwa-mem2/pull/204) | open PR |
| libsais FM-index construction | [#57](https://github.com/fg-labs/bwa-mem3/pull/57) | — | fork-only |
| Consolidated mapping speedups | [#58](https://github.com/fg-labs/bwa-mem3/pull/58) | — | fork-only |
| **Features** | | | |
| `--bam=LEVEL` direct BAM output | [#12](https://github.com/fg-labs/bwa-mem3/pull/12) | — | fork-only |
| `--meth` bisulfite alignment mode | [#13](https://github.com/fg-labs/bwa-mem3/pull/13) | — | fork-only |
| Vendored mimalloc allocator | [#19](https://github.com/fg-labs/bwa-mem3/pull/19) | — | fork-only |
| `HN:i` hit count tag | [#42](https://github.com/fg-labs/bwa-mem3/pull/42) | [lh3/bwa#438](https://github.com/lh3/bwa/pull/438) | analogous to bwa aln; no direct upstream port |
| `--supp-rep-hard-cap` MAPQ rescoring | [#56](https://github.com/fg-labs/bwa-mem3/pull/56) | [bwa-mem2#260](https://github.com/bwa-mem2/bwa-mem2/issues/260) | open issue |
| `bwa-mem3 shm` shared-memory index | [#65](https://github.com/fg-labs/bwa-mem3/pull/65) | — | fork-only (v1 feature port) |
| `shm --meth` symmetry | [#67](https://github.com/fg-labs/bwa-mem3/pull/67) | — | fork-only |
| `-z FLOAT` (XA_drop_ratio CLI knob) | [#35](https://github.com/fg-labs/bwa-mem3/pull/35) | [lh3/bwa#294](https://github.com/lh3/bwa/pull/294) | merged (lh3 only) |
| `-u` flag — widen `XA:Z` records with `,score,mapq` | [#35](https://github.com/fg-labs/bwa-mem3/pull/35) | [lh3/bwa#293](https://github.com/lh3/bwa/pull/293) | merged (lh3 only) |
| `MQ:i` mate mapping quality tag | [#35](https://github.com/fg-labs/bwa-mem3/pull/35) | [lh3/bwa#330](https://github.com/lh3/bwa/pull/330) | merged (lh3 only) |
| Bismark-compatible `XR:Z` / `XG:Z` / `XM:Z` tags | [#90](https://github.com/fg-labs/bwa-mem3/pull/90) | — | fork-only |
| **Architecture support** | | | |
| Linux ARM64 / aarch64 build + CI | [#1](https://github.com/fg-labs/bwa-mem3/pull/1) | [bwa-mem2#288](https://github.com/bwa-mem2/bwa-mem2/pull/288) | open PR |
| `arch=avx512bw` explicit Makefile target | [#16](https://github.com/fg-labs/bwa-mem3/pull/16) | — | fork-only |
| NEON kswv mate-rescue kernel | [#18](https://github.com/fg-labs/bwa-mem3/pull/18) | — | fork-only |
| AVX2 kswv mate-rescue kernel | [#20](https://github.com/fg-labs/bwa-mem3/pull/20) | — | fork-only |
| **Build & infrastructure** | | | |
| doctest framework + Codecov | [#34](https://github.com/fg-labs/bwa-mem3/pull/34) | — | fork-only |
| `PACKAGE_VERSION` from `git describe` | [#52](https://github.com/fg-labs/bwa-mem3/pull/52) | [bwa-mem2#283](https://github.com/bwa-mem2/bwa-mem2/issues/283), [bwa-mem2#284](https://github.com/bwa-mem2/bwa-mem2/pull/284) | open issue + open PR |
| PGO target parameterization | [#59](https://github.com/fg-labs/bwa-mem3/pull/59) | — | fork-only |
| `CXXFLAGS`/`CPPFLAGS`/`LDFLAGS` forwarding | [#50](https://github.com/fg-labs/bwa-mem3/pull/50) | [bwa-mem2#290](https://github.com/bwa-mem2/bwa-mem2/pull/290) | open upstream PR |
| Unit-test harness + ARM CI | [#23](https://github.com/fg-labs/bwa-mem3/pull/23) | — | fork-only |
| CI matrix expansion (7 rows) | [#24](https://github.com/fg-labs/bwa-mem3/pull/24) | — | fork-only |

## Upstream issues tracked but not yet fixed in bwa-mem3

The following upstream issues are tracked in the bwa-mem3 issue list but do
not yet have corresponding fixes in `main`:

| Issue | Upstream reference | Notes |
|-------|--------------------|-------|
| Split-alignment evidence loss vs bwa 0.7.17 | [bwa-mem2#273](https://github.com/bwa-mem2/bwa-mem2/issues/273) | [issue #47](https://github.com/fg-labs/bwa-mem3/issues/47) — under investigation |
| MAPQ/coordinate parity vs bwa mem 0.7.18 | [bwa-mem2#262](https://github.com/bwa-mem2/bwa-mem2/issues/262), [bwa-mem2#246](https://github.com/bwa-mem2/bwa-mem2/issues/246), [bwa-mem2#239](https://github.com/bwa-mem2/bwa-mem2/issues/239) | [issue #48](https://github.com/fg-labs/bwa-mem3/issues/48) — tracking only |

---

**See also:**
[Correctness fixes](correctness.md) ·
[Performance improvements](performance.md) ·
[Features](features.md) ·
[Architecture support](arch-support.md) ·
[Build & infrastructure](build-infra.md)
