# Equivalence with bwa-mem2 (bit-identity)

bwa-mem3 is **not** byte-identical to bwa-mem2, and that is intentional. Upstream bwa-mem2 advertises exact equivalence with the original `bwa` — "produces alignment identical to bwa" and "exact same output as bwa-mem(2)" — and that guarantee was the right bar for a project whose sole charter was to reproduce `bwa mem` faster. bwa-mem3 has a broader charter: it adds informative SAM tags, fixes crashes and undefined behavior, corrects SIMD scoring kernels, and makes tie resolution deterministic. Several of those changes necessarily move output away from a byte-for-byte match. We have consciously stepped below the bit-identity bar, and this page records exactly where, and gives an auditable trail back to every merged pull request so a reader can decide for themselves whether each divergence matters for their workflow.

The short version: on the data we have tested, the **core alignment** — where each read maps and how — is preserved on essentially every read. The **SAM byte stream** is not, primarily because bwa-mem3 emits additional auxiliary tags that upstream never wrote, and secondarily because it adds a handful of supplementary alignments and shifts MAPQ/CIGAR on a small per-architecture fraction of reads. Beyond the additive tags, the remaining divergences are latent, opt-in, or per-architecture; they are described under "What differs" below and in more depth on the [Correctness fixes](correctness.md), [Performance improvements](performance.md), and [Features](features.md) pages.

## What is preserved

We ran an empirical concordance check with [bwa-mem3-bench](https://github.com/fg-labs/bwa-mem3-bench) at commit `bffae5a` (current `main`), comparing bwa-mem3 against upstream `bwa-mem2` v2.2.1 on x86 hosts across whole-genome, whole-exome, and panel workloads. **Primary-alignment concordance** — reference name, position, CIGAR, MAPQ, and placement flags compared per read end — is:

| sample | primary concordance | primary records |
|---|---:|---:|
| wes-5M | 99.9996% | 10,051,170 |
| wgs-5M | 99.9893% | 9,980,872 |
| panel-twist-5M | 99.9414% | 7,913,324 |

Where each read maps is preserved on essentially every read. The well-under-0.1% of primary records that differ do so in `MAPQ`, `CIGAR`, or position, and are accounted for by the per-architecture SIMD `score2`/`MAPQ` convergence and the deterministic tie-break change described under "What differs" below and on the [Correctness fixes](correctness.md) page. (On the 1M-read `smoke-1M` cell the figure is 99.946%; the larger exome and genome cells above are more representative.) Cross-architecture, the NEON (ARM) and x86 builds are byte-identical to each other — 100.0000% concordance over all records, supplementary alignments included.

## What differs

### Additive SAM tags

The most pervasive difference is two **additive** auxiliary tags that bwa-mem3 emits and upstream does not:

- `MQ:i` — mate mapping quality, present on ~100% of bwa-mem3 records and absent from upstream output.
- `HN:i` — total hit count per primary, present on 54,188 of the 64,763 bwa-mem3 records and on 0 upstream records.

A representative record (`SRR34589119.1`) makes the shape of the difference concrete. bwa-mem3 emits:

```text
… MC:Z:53S12M1D6M2D22M58S  MQ:i:19  AS:i:25  XS:i:23  HN:i:9
```

while upstream emits:

```text
… MC:Z:53S12M1D6M2D22M58S  AS:i:25  XS:i:23
```

Same alignment, same scores — two extra tags. Because these tags are inserted into the optional-field area of the record, the line is no longer byte-identical to upstream even though the alignment it describes is. `MQ:i` is one of the lh3/bwa tags ported in [#35](https://github.com/fg-labs/bwa-mem3/pull/35); `HN:i` is added in [#42](https://github.com/fg-labs/bwa-mem3/pull/42). See [Features → `HN:i` hit count tag](features.md) for the full semantics.

Separately, the `@PG` header line reports `ID:bwa-mem3` / `PN:bwa-mem3` rather than `bwa-mem2`, which is also a byte-level header difference by design.

### Additional supplementary alignments

On the default build, with no special flags, bwa-mem3 emits a small number of **additional supplementary (chimeric/split) alignments** that upstream `bwa-mem2` v2.2.1 does not. On `wes-5M` (5,025,585 read pairs) bwa-mem3 emits 5,123 supplementary records versus upstream's 5,118 — five extra split alignments, on five templates (0.0001%). The **primary** alignment of every affected pair is unchanged; only an extra supplementary record is added. This is measured by bwa-mem3-bench's `compare-bams`, which reports per-template supplementary count-mismatch and position-unmatched rates alongside primary concordance. These additions are default-on behavior — they occur with no special flags — and are *not* a product of the opt-in `--supp-rep-hard-cap` rescoring, which only lowers MAPQ and never adds records. Pinning them to a specific upstream-divergence PR is tracked as follow-up.

### Divergences that are latent, opt-in, or per-architecture

The following changes can move alignments, scores, or MAPQ relative to upstream, but did **not** surface as primary-alignment differences on the measured cells because they are gated, latent, or only active on other inputs or architectures:

- **Proper-pair `FLAG` recompute ([#17](https://github.com/fg-labs/bwa-mem3/pull/17), default-on).** bwa-mem3 computes the `0x2` bit from the alignment actually emitted rather than from the below-threshold primary. This only changes the flag in the rare case where the primary's score is under `opt->T` but an ALT hit clears it; on `smoke-1M` no record hit that path, so the full `FLAG` matched upstream exactly. See [Correctness fixes → Proper-pair flag](correctness.md).
- **SIMD scoring-kernel fixes ([#21](https://github.com/fg-labs/bwa-mem3/pull/21), [#26](https://github.com/fg-labs/bwa-mem3/pull/26), [#28](https://github.com/fg-labs/bwa-mem3/pull/28), [#29](https://github.com/fg-labs/bwa-mem3/pull/29), [#30](https://github.com/fg-labs/bwa-mem3/pull/30), [#31](https://github.com/fg-labs/bwa-mem3/pull/31)).** These correct the batched mate-rescue `kswv` kernels so the suboptimal score (`score2` → `XS:i`/`MAPQ`) converges toward the scalar `ksw_align2` reference. They move `XS`/`MAPQ` on the minority of reads where the SIMD kernel previously diverged, and the affected reads differ by architecture (AVX2 vs NEON vs AVX-512BW). See [Correctness fixes → kswv score2 plateau series](correctness.md).
- **Seeding correctness fixes ([#55](https://github.com/fg-labs/bwa-mem3/pull/55), [#73](https://github.com/fg-labs/bwa-mem3/pull/73), [#100](https://github.com/fg-labs/bwa-mem3/pull/100)).** These fix buffer sizing and a prefetch-mask precedence bug. They change alignments only where the old bug actually triggered (e.g. reads longer than 151 bp for [#55](https://github.com/fg-labs/bwa-mem3/pull/55); [#73](https://github.com/fg-labs/bwa-mem3/pull/73) is a prefetch hint with no semantic change).
- **Opt-in MAPQ rescoring ([#56](https://github.com/fg-labs/bwa-mem3/pull/56), [#101](https://github.com/fg-labs/bwa-mem3/pull/101), [#118](https://github.com/fg-labs/bwa-mem3/pull/118), default-off).** `--supp-rep-hard-cap INT` forces MAPQ=0 on supplementary alignments anchored in repetitive seeds. With no flag the output is unchanged; [#101](https://github.com/fg-labs/bwa-mem3/pull/101) makes the flag actually take effect (it shipped as a silent no-op before), and [#118](https://github.com/fg-labs/bwa-mem3/pull/118) is its regression test. See [Features → `--supp-rep-hard-cap`](features.md).
- **Tie-break determinism ([#123](https://github.com/fg-labs/bwa-mem3/pull/123)).** Makes secondary-alignment ordering deterministic across runs; can reorder equal-scoring ties relative to upstream's order.

## Auditable PR list

Every merged pull request on `main` is linked below, grouped by whether and how it can affect output, output-affecting groups first. The build / CI / docs / chore group at the end has no output impact and is collapsed but still links each PR. PR link format is `https://github.com/fg-labs/bwa-mem3/pull/<N>`. For the commit-level table and upstream disposition see the [Overview](overview.md) and [Upstream PR status](upstream-prs.md) pages.

### Additive SAM tags (not present in bwa-mem2)

These add optional fields and therefore break byte-identity by construction, without changing the alignment.

- [#42](https://github.com/fg-labs/bwa-mem3/pull/42) — `HN:i` total-hit-count tag per primary.
- [#35](https://github.com/fg-labs/bwa-mem3/pull/35) — ports four lh3/bwa PRs, including the `MQ` mate-MAPQ tag, the `-u`/`XB` tag, and `@HD` ordering.
- [#90](https://github.com/fg-labs/bwa-mem3/pull/90) — Bismark-compatible `XR`/`XG`/`XM` methylation tags (meth mode only).

### Proper-pair FLAG semantics (default-on, latent on smoke-1M)

- [#17](https://github.com/fg-labs/bwa-mem3/pull/17) — compute the `0x2` proper-pair flag from the emitted alignment.

### SIMD scoring-kernel fixes (per-arch `score2` → `XS`/MAPQ; converge toward scalar)

- [#21](https://github.com/fg-labs/bwa-mem3/pull/21) — apply NEON score2-scan fixes to the AVX-512BW kernel.
- [#26](https://github.com/fg-labs/bwa-mem3/pull/26) — gate AVX2 arch dispatch on `!__AVX512BW__`.
- [#28](https://github.com/fg-labs/bwa-mem3/pull/28) — consolidate AVX2 score2 plateaus to match scalar `ksw_align2`.
- [#29](https://github.com/fg-labs/bwa-mem3/pull/29) — port score2 plateau consolidation to NEON + AVX-512BW.
- [#30](https://github.com/fg-labs/bwa-mem3/pull/30) — port the score2 plateau fix to `kswv_512_16` (AVX-512BW 16-bit).
- [#31](https://github.com/fg-labs/bwa-mem3/pull/31) — rewrite `kswv_neon_16` as a real SIMD kernel with correct table + score2.

### Seeding correctness fixes (can change alignments where the old bug triggered)

- [#55](https://github.com/fg-labs/bwa-mem3/pull/55) — size SMEM buffers from observed max read length (fixes >151 bp corruption).
- [#100](https://github.com/fg-labs/bwa-mem3/pull/100) — track `enc_qdb` byte capacity separately from `wsize_mem`.
- [#73](https://github.com/fg-labs/bwa-mem3/pull/73) — parenthesize `SA_COMPX_MASK` precedence in the sampled-SA prefetch (prefetch hint; no semantic change).

### Opt-in MAPQ rescoring (default OFF)

- [#56](https://github.com/fg-labs/bwa-mem3/pull/56) — `--supp-rep-hard-cap` opt-in supplementary MAPQ rescoring.
- [#101](https://github.com/fg-labs/bwa-mem3/pull/101) — propagate SMEM SA-count to seed `n_hits` so `--supp-rep-hard-cap` actually fires.
- [#118](https://github.com/fg-labs/bwa-mem3/pull/118) — regression test for `--supp-rep-hard-cap` on a repetitive workload.

### Tie-break determinism

- [#123](https://github.com/fg-labs/bwa-mem3/pull/123) — deterministic tie resolution (open at the time this page shipped; this docs page stacks on it).

### New modes / output formats (separate features)

- [#13](https://github.com/fg-labs/bwa-mem3/pull/13) — `--meth` bisulfite alignment mode.
- [#90](https://github.com/fg-labs/bwa-mem3/pull/90) — Bismark-compatible `XR`/`XG`/`XM` meth tags (also listed under additive tags).
- [#12](https://github.com/fg-labs/bwa-mem3/pull/12) — `--bam[=LEVEL]` direct BAM output.
- [#65](https://github.com/fg-labs/bwa-mem3/pull/65) — `bwa-mem3 shm` shared-memory index.
- [#67](https://github.com/fg-labs/bwa-mem3/pull/67) — `shm --meth` symmetry.
- [#82](https://github.com/fg-labs/bwa-mem3/pull/82) — serialize `/bwactl` RMW with a POSIX named semaphore.
- [#83](https://github.com/fg-labs/bwa-mem3/pull/83) — single-binary in-process SIMD dispatch (replaces the multi-binary `execv` launcher).
- [#57](https://github.com/fg-labs/bwa-mem3/pull/57) — libsais-based memory-bounded FM-index construction.

### Core-kernel performance (intended output-neutral; smoke evidence consistent)

- [#33](https://github.com/fg-labs/bwa-mem3/pull/33) — lockstep SMEM batching across N reads.
- [#49](https://github.com/fg-labs/bwa-mem3/pull/49) — batch `-H` header ingestion (fixes O(n²) header read).
- [#58](https://github.com/fg-labs/bwa-mem3/pull/58) — consolidated mapping speedups (ksw2, SMEM, SAL, SAM).
- [#70](https://github.com/fg-labs/bwa-mem3/pull/70) — per-strip L1 prefetches across all u8/16 kswv kernels.
- [#75](https://github.com/fg-labs/bwa-mem3/pull/75) — bump `SMEM_LOCKSTEP_N` from 8 to 16.
- [#76](https://github.com/fg-labs/bwa-mem3/pull/76) — convert `mem_matesw_batch_{pre,post}` to `bns_fetch_seq_v2`.
- [#77](https://github.com/fg-labs/bwa-mem3/pull/77) — closed-form HIT for `total_mis == 0` in the ungapped path.
- [#78](https://github.com/fg-labs/bwa-mem3/pull/78) — replace per-call malloc with an on-stack buffer for small `n` in ksort.
- [#80](https://github.com/fg-labs/bwa-mem3/pull/80) — skip wasted zero-init on libsais unpack + SA buffers.
- [#86](https://github.com/fg-labs/bwa-mem3/pull/86) — cap AVX-512BW autovec at 256-bit; `bwa_shm` `/dev/shm` preflight.
- [#88](https://github.com/fg-labs/bwa-mem3/pull/88) — inline `backwardExt` to recover a gcc 12+ wall-clock regression.

### Crash / UB / safety fixes (output-neutral except where UB previously corrupted results)

- [#22](https://github.com/fg-labs/bwa-mem3/pull/22) — zero `bseq1_t` in `kseq2bseq1` (crash on realloc growth).
- [#51](https://github.com/fg-labs/bwa-mem3/pull/51) — guard the post-loop `rowMax` store on `nrow==0` batches.
- [#74](https://github.com/fg-labs/bwa-mem3/pull/74) — bound the `.alt` parse buffer to prevent stack overflow.
- [#85](https://github.com/fg-labs/bwa-mem3/pull/85) — copy the ref slice before `ksw_align2` to avoid SIGSEGV on a shm-backed `ref_string`.
- [#117](https://github.com/fg-labs/bwa-mem3/pull/117) — free lockstep SMEM caches at thread exit.

### Build / CI / docs / chore (no output impact)

These do not change alignment output. Each is linked for completeness.

- Build / Make / link: [#16](https://github.com/fg-labs/bwa-mem3/pull/16), [#50](https://github.com/fg-labs/bwa-mem3/pull/50), [#52](https://github.com/fg-labs/bwa-mem3/pull/52), [#53](https://github.com/fg-labs/bwa-mem3/pull/53), [#59](https://github.com/fg-labs/bwa-mem3/pull/59), [#84](https://github.com/fg-labs/bwa-mem3/pull/84), [#105](https://github.com/fg-labs/bwa-mem3/pull/105), [#108](https://github.com/fg-labs/bwa-mem3/pull/108), [#122](https://github.com/fg-labs/bwa-mem3/pull/122).
- CLI / API plumbing (no default-output change): [#5](https://github.com/fg-labs/bwa-mem3/pull/5), [#6](https://github.com/fg-labs/bwa-mem3/pull/6), [#7](https://github.com/fg-labs/bwa-mem3/pull/7), [#8](https://github.com/fg-labs/bwa-mem3/pull/8), [#9](https://github.com/fg-labs/bwa-mem3/pull/9), [#54](https://github.com/fg-labs/bwa-mem3/pull/54), [#60](https://github.com/fg-labs/bwa-mem3/pull/60), [#81](https://github.com/fg-labs/bwa-mem3/pull/81).
- Allocator: [#19](https://github.com/fg-labs/bwa-mem3/pull/19) — vendored mimalloc (output-neutral; can change crash/no-crash behavior on the buggy paths above).
- CI / test harness / fixtures: [#1](https://github.com/fg-labs/bwa-mem3/pull/1), [#2](https://github.com/fg-labs/bwa-mem3/pull/2), [#4](https://github.com/fg-labs/bwa-mem3/pull/4), [#10](https://github.com/fg-labs/bwa-mem3/pull/10), [#18](https://github.com/fg-labs/bwa-mem3/pull/18), [#20](https://github.com/fg-labs/bwa-mem3/pull/20), [#23](https://github.com/fg-labs/bwa-mem3/pull/23), [#24](https://github.com/fg-labs/bwa-mem3/pull/24), [#34](https://github.com/fg-labs/bwa-mem3/pull/34), [#63](https://github.com/fg-labs/bwa-mem3/pull/63), [#68](https://github.com/fg-labs/bwa-mem3/pull/68), [#72](https://github.com/fg-labs/bwa-mem3/pull/72), [#89](https://github.com/fg-labs/bwa-mem3/pull/89), [#102](https://github.com/fg-labs/bwa-mem3/pull/102), [#111](https://github.com/fg-labs/bwa-mem3/pull/111), [#112](https://github.com/fg-labs/bwa-mem3/pull/112), [#113](https://github.com/fg-labs/bwa-mem3/pull/113), [#114](https://github.com/fg-labs/bwa-mem3/pull/114), [#115](https://github.com/fg-labs/bwa-mem3/pull/115).
- Arch support / SIMD floor precheck: [#1](https://github.com/fg-labs/bwa-mem3/pull/1) (arm64 build, also above), [#95](https://github.com/fg-labs/bwa-mem3/pull/95).
- Docs / release / metadata: [#3](https://github.com/fg-labs/bwa-mem3/pull/3), [#71](https://github.com/fg-labs/bwa-mem3/pull/71), [#79](https://github.com/fg-labs/bwa-mem3/pull/79), [#93](https://github.com/fg-labs/bwa-mem3/pull/93), [#94](https://github.com/fg-labs/bwa-mem3/pull/94), [#96](https://github.com/fg-labs/bwa-mem3/pull/96), [#97](https://github.com/fg-labs/bwa-mem3/pull/97), [#98](https://github.com/fg-labs/bwa-mem3/pull/98), [#106](https://github.com/fg-labs/bwa-mem3/pull/106), [#107](https://github.com/fg-labs/bwa-mem3/pull/107).

---

**See also:**
[Overview](overview.md) ·
[Correctness fixes](correctness.md) ·
[Performance improvements](performance.md) ·
[Features](features.md) ·
[Upstream PR status](upstream-prs.md)
