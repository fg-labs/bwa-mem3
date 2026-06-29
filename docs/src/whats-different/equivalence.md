# Equivalence with bwa-mem2 (bit-identity)

bwa-mem3 is **not** byte-identical to bwa-mem2, and that is intentional. Upstream bwa-mem2 advertises exact equivalence with the original `bwa` — "produces alignment identical to bwa" and "exact same output as bwa-mem(2)" — and that guarantee was the right bar for a project whose sole charter was to reproduce `bwa mem` faster. bwa-mem3 has a broader charter: it adds informative SAM tags, fixes crashes and undefined behavior, corrects SIMD scoring kernels, and makes tie resolution deterministic. Several of those changes necessarily move output away from a byte-for-byte match. We have consciously stepped below the bit-identity bar, and this page records exactly where, and gives an auditable trail back to every merged pull request so a reader can decide for themselves whether each divergence matters for their workflow.

The short version: on the data we have tested, the **core alignment** — where each read maps and how — is preserved on essentially every read. The **SAM byte stream** is not, primarily because bwa-mem3 emits additional auxiliary tags that upstream never wrote, and secondarily because it adds a handful of supplementary alignments and shifts MAPQ/CIGAR on a small per-architecture fraction of reads. Beyond the additive tags, the remaining divergences are latent, opt-in, or per-architecture; they are described under "What differs" below and in more depth on the [Correctness fixes](correctness.md), [Performance improvements](performance.md), and [Features](features.md) pages.

## What is preserved

We ran an empirical concordance check with [bwa-mem3-bench](https://github.com/fg-labs/bwa-mem3-bench) at commit `a02fcb4` (current `main`), comparing bwa-mem3 against upstream `bwa-mem2` v2.2.1 on x86 hosts across whole-genome, whole-exome, and panel workloads. **Primary-alignment concordance** — reference name, position, CIGAR, MAPQ, and placement flags compared per read end — is:

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

On the default build, with no special flags, bwa-mem3 emits a small number of **additional supplementary (chimeric/split) alignments** that upstream `bwa-mem2` v2.2.1 does not. On `wes-5M` (5,025,585 read pairs) bwa-mem3 emits 5,123 supplementary records versus upstream's 5,118 — five extra split alignments, on five templates (0.0001%). The **primary** alignment of every affected pair is unchanged; only an extra supplementary record is added. This is measured by bwa-mem3-bench's `compare-bams`, which reports per-template supplementary count-mismatch and position-unmatched rates alongside primary concordance. These additions are default-on behavior — they occur with no special flags — and are *not* a product of the opt-in `--supp-rep-hard-cap` rescoring, which only lowers MAPQ and never adds records. Pinning them to a specific upstream-divergence PR is tracked as follow-up ([#127](https://github.com/fg-labs/bwa-mem3/issues/127)).

### Divergences that are latent, opt-in, or per-architecture

The following changes can move alignments, scores, or MAPQ relative to upstream, but did **not** surface as primary-alignment differences on the measured cells because they are gated, latent, or only active on other inputs or architectures:

- **Proper-pair `FLAG` recompute ([#17](https://github.com/fg-labs/bwa-mem3/pull/17), default-on).** bwa-mem3 computes the `0x2` bit from the alignment actually emitted rather than from the below-threshold primary. This only changes the flag in the rare case where the primary's score is under `opt->T` but an ALT hit clears it; on `smoke-1M` no record hit that path, so the full `FLAG` matched upstream exactly. See [Correctness fixes → Proper-pair flag](correctness.md).
- **SIMD scoring-kernel fixes ([#21](https://github.com/fg-labs/bwa-mem3/pull/21), [#26](https://github.com/fg-labs/bwa-mem3/pull/26), [#28](https://github.com/fg-labs/bwa-mem3/pull/28), [#29](https://github.com/fg-labs/bwa-mem3/pull/29), [#30](https://github.com/fg-labs/bwa-mem3/pull/30), [#31](https://github.com/fg-labs/bwa-mem3/pull/31)).** These correct the batched mate-rescue `kswv` kernels so the suboptimal score (`score2` → `XS:i`/`MAPQ`) converges toward the scalar `ksw_align2` reference. They move `XS`/`MAPQ` on the minority of reads where the SIMD kernel previously diverged, and the affected reads differ by architecture (AVX2 vs NEON vs AVX-512BW). See [Correctness fixes → kswv score2 plateau series](correctness.md).
- **Seeding correctness fixes ([#55](https://github.com/fg-labs/bwa-mem3/pull/55), [#73](https://github.com/fg-labs/bwa-mem3/pull/73), [#100](https://github.com/fg-labs/bwa-mem3/pull/100)).** These fix buffer sizing and a prefetch-mask precedence bug. They change alignments only where the old bug actually triggered (e.g. reads longer than 151 bp for [#55](https://github.com/fg-labs/bwa-mem3/pull/55); [#73](https://github.com/fg-labs/bwa-mem3/pull/73) is a prefetch hint with no semantic change).
- **Opt-in MAPQ rescoring ([#56](https://github.com/fg-labs/bwa-mem3/pull/56), [#101](https://github.com/fg-labs/bwa-mem3/pull/101), [#118](https://github.com/fg-labs/bwa-mem3/pull/118), default-off).** `--supp-rep-hard-cap INT` forces MAPQ=0 on supplementary alignments anchored in repetitive seeds. With no flag the output is unchanged; [#101](https://github.com/fg-labs/bwa-mem3/pull/101) makes the flag actually take effect (it shipped as a silent no-op before), and [#118](https://github.com/fg-labs/bwa-mem3/pull/118) is its regression test. See [Features → `--supp-rep-hard-cap`](features.md).
- **Tie-break determinism ([#123](https://github.com/fg-labs/bwa-mem3/pull/123)).** Makes secondary-alignment ordering deterministic across runs; can reorder equal-scoring ties relative to upstream's order.

(The resolve→order→chain seeding refactor that backs `--seed-order` is byte-identical in its default `off` mode; it is described in the dedicated section below rather than listed here, since only its non-`off` modes are divergent.)

### Seed ordering (`--seed-order`, opt-in)

`--seed-order off` (the default) preserves byte-identical output. The internal resolve→order→chain refactor that enables it was verified to be bit-for-bit identical to the previous code path across single-end, paired-end, and threaded runs.

`--seed-order local-longest` (and the unadvertised advanced modes `global-longest`, `absorb-count`, `most-absorb`) are **opt-in** and are **not byte-identical**. They reorder each read's SA-resolved seeds before chaining so that the longest seeds anchor chains first; contained shorter seeds are then absorbed rather than extended. The output shift is non-trivial: secondary alignments, `XA:Z:`, `XS:i`, and `HN:i` can all change, and a small number of primary alignments may shift as well.

Accuracy on an easy simulated profile (holodeck, ~94.4 % F1) is flat relative to `off`. Hard-data F1 validation on divergent/indel-rich reads and GIAB benchmarks is **not yet complete**. Because the accuracy gate is not fully passed, all non-`off` modes remain opt-in only; the default stays `off`. See [Optimization checklist → Reorder seeds longest-first](../best-practices/optimization-checklist.md#6-reorder-seeds-longest-first---seed-order-local-longest) for usage guidance.

## Declared divergence catalog

The divergences described above are tracked as a structured registry in [bwa-mem3-bench](https://github.com/fg-labs/bwa-mem3-bench) (`docs/expected-divergences.yaml`). Each carries a per-sample concordance-drift budget that the benchmark gates against on every run — a new bwa-mem3 build that drifts beyond its budget fails CI rather than silently shipping a regression. The table below is generated from that registry; do not edit it by hand (see [bwa-mem3-bench → Per-release concordance history](../related-projects/bwa-mem3-bench.md) for how it is regenerated).

<!-- FG-DIVERGENCE-CATALOG:start -->
| id | pr | affected | samples | budget_% | summary |
| --- | --- | --- | --- | --- | --- |
| FG-PRIMARY-DRIFT | fg-labs/bwa-mem3#123 | primary_alignment | wgs-5M, wes-5M, panel-twist-5M, smoke-1M | 0.1000 | Per-architecture SIMD score2/MAPQ convergence (#21, #26, #28-#31) and deterministic tie-break ordering (#123) shift MAPQ, CIGAR, or position on a small fraction of primary alignments relative to bwa-mem2 v2.2.1. Where each read maps is preserved; the affected reads differ in placement detail, and the set varies by SIMD architecture. |
| FG-METH-DIVERGENCE | fg-labs/bwa-mem3#90 | meth_alignment | meth-twist-emseq-5M, smoke-meth | 1.5000 | Bisulfite (--meth) mode against the bwameth.py baseline diverges beyond the ignored YD/XM/XG tag set (Bismark-compatible XR/XG/XM tags and C->T/G->A conversion handling), giving a larger but still-bounded concordance drift on methylation workloads. |
| FG-SUPP-ADDITIONS | TBD | supplementary_alignment | all | 0.0000 | bwa-mem3 emits a small number of additional supplementary (split/chimeric) alignments vs bwa-mem2 v2.2.1 on the default build (e.g. wes-5M: 5123 vs 5118). Primary alignments are unchanged. Tracked as a supplementary-count metric (compare-bams supp_count_mismatch / supp_unmatched); it does not affect the primary-concordance drift budget, hence 0.0 here. |
<!-- FG-DIVERGENCE-CATALOG:end -->

## Per-PR audit trail

Every fork-carried change is listed — with its class and upstream bwa-mem2
disposition — in the [PR catalog](../reference/pr-catalog.md). The declared
divergence catalog above calls out the entries that actually affect output.

---

**See also:**
[Overview](overview.md) ·
[Correctness fixes](correctness.md) ·
[Performance improvements](performance.md) ·
[Features](features.md) ·
[PR catalog](../reference/pr-catalog.md)
