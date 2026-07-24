# Equivalence with bwa-mem2 (bit-identity)

<div class="warning" style="border:3px solid #2e7d32;border-radius:8px;padding:1rem 1.25rem;margin:1.25rem 0;background:#e8f5e9;color:#1b1b1b;">
<strong style="font-size:1.35em;display:block;margin-bottom:0.4rem;">✅ No alignment record bwa-mem2 emits is changed or removed — and on the cell re-measured for 0.7.1 its primary alignments are reproduced exactly.</strong>
On the drop-in profile (plain <code>bwa-mem3 mem</code>, no flags), release 0.7.1 gives the <strong>same primary alignment as bwa-mem2</strong> for <strong>every</strong> read of the slice re-measured after the parity work landed (1,066,668 records, HG00096 WGS, hg38, x86): identical <code>FLAG</code>, <code>RNAME</code>, <code>POS</code>, <code>MAPQ</code>, <code>CIGAR</code>, <code>AS</code>, and <code>XS</code>. On those records the SAM byte stream differs by <strong>three things — two additive tags plus a renamed header line</strong>:
<ul style="margin:0.5rem 0 0 0;">
<li><code>MQ:i</code> — mate mapping quality (an extra tag bwa-mem2 never wrote)</li>
<li><code>HN:i</code> — hit count (an extra tag bwa-mem2 never wrote)</li>
<li>the <code>@PG</code> line names <code>bwa-mem3</code> instead of <code>bwa-mem2</code></li>
</ul>
Strip those two tags and normalize <code>@PG</code> and those records are byte-for-byte identical.
<br><br>
<strong>Two things are not yet closed out.</strong> A residual <strong>+4 supplementary alignments</strong> on the 5 M WGS cell — also purely additive, but a real parity gap: it was last measured after <a href="https://github.com/fg-labs/bwa-mem3/pull/257">#257</a> and <em>not</em> re-measured on a build carrying <a href="https://github.com/fg-labs/bwa-mem3/pull/268">#268</a>. And the genome-wide, multi-sample, cross-architecture re-run that would extend all of the above past the cells named on this page. Treat every claim here as scoped to the evidence it cites, not as a guarantee across all inputs and architectures. (The opt-in <code>--fast</code> speed levers are a separate, deliberate deviation — see below.)
</div>

As of **release 0.7.1**, the bwa-mem3 **drop-in profile** (plain `bwa-mem3 mem`, no
extra flags) produces the **same primary alignment** as bwa-mem2 — same `FLAG`,
`RNAME`, `POS`, `MAPQ`, `CIGAR`, `AS`, `XS` — on the data we have re-measured since
the restoration. The **SAM byte stream** is still not a literal byte-for-byte
match, for two deliberate reasons only: bwa-mem3 emits extra auxiliary tags
(`MQ:i`, `HN:i`) that upstream never wrote, and its `@PG` line names `bwa-mem3`. Strip
those tags and normalize `@PG` and those records are identical. This is a change in
posture from earlier releases, which carried three changes that *did* move primary
alignments; all three were reverted for 0.7.1 (see [Restored bwa-mem2 parity
(0.7.1)](#restored-bwa-mem2-parity-071) below).

Two things scope that statement, and both are tracked below rather than assumed away.
The **supplementary** set is not yet confirmed to match: the extra split alignments
were cut from +239 to +4 by [#257](https://github.com/fg-labs/bwa-mem3/pull/257), and
that +4 has not been re-measured on a build carrying
[#268](https://github.com/fg-labs/bwa-mem3/pull/268) (see [Additional supplementary
alignments](#additional-supplementary-alignments-reduced-in-071-re-measurement-pending)).
And the post-restoration measurements cover one sample on one architecture; the
genome-wide, multi-sample, cross-architecture concordance re-run is still pending. Read
every figure on this page as scoped to the cell it names.

Upstream bwa-mem2 advertises exact equivalence with the original `bwa` — "produces
alignment identical to bwa" and "exact same output as bwa-mem(2)" — and that guarantee
was the right bar for a project whose sole charter was to reproduce `bwa mem` faster.
bwa-mem3 has a broader charter: it adds informative SAM tags, fixes crashes and
undefined behavior, corrects SIMD scoring kernels, and makes tie resolution
deterministic. Where those changes would have moved *which alignment* a read gets, we
have chosen bwa-mem2 output compatibility on the default path; where they are purely
additive (the extra tags), confined to the `@PG` header naming, or strictly opt-in (the
`--fast` speed levers), they remain. This page records exactly where output does and does not match, and gives an auditable
trail back to every merged pull request — both the ones that *introduced* a divergence
and the ones that *restored* parity — so a reader can decide for themselves whether each
difference matters for their workflow.

The short version: on the drop-in profile the **core alignment** — where each read maps
and how, plus its scores — is preserved on essentially every read across every cell we
have measured, and on the one cell re-measured after the 0.7.1 restoration it is
preserved exactly. On those records the **SAM byte stream** differs only by the additive `MQ:i`/`HN:i` tags and the
`@PG` line. The **supplementary set** is the one piece not yet re-confirmed (a residual
+4, pending re-measurement). The **opt-in speed levers** (`--fast` and its constituent flags) are
separately, deliberately *not* byte-identical; they are listed under "What differs" and
described in more depth on the [Correctness fixes](correctness.md), [Performance
improvements](performance.md), and [Features](features.md) pages.

## Restored bwa-mem2 parity (0.7.1)

Three fork-carried changes shipped between 0.5.0 and 0.7.0 each moved **primary**
alignments away from bwa-mem2 on the default path. Each was root-caused to a specific
mechanism and reverted (or repacked) for 0.7.1, pairing every divergence-introducing PR
with the PR that restored parity:

| divergence introduced by | mechanism | restored by |
|---|---|---|
| [#141](https://github.com/fg-labs/bwa-mem3/pull/141) — banded-SW extension opened gaps from the cell max `H` (standard Gotoh) instead of the diagonal `M` | `gscore` inflates upward-only, flipping the `gscore <= score - pen_clip3` local-vs-extend test in `mem_chain2aln` | [#256](https://github.com/fg-labs/bwa-mem3/pull/256) — derive extension gaps from `M` again (bwa-mem2's `bandedSWA` convention) |
| [#123](https://github.com/fg-labs/bwa-mem3/pull/123) — replaced bwa-mem2's partial-order dedup comparator with a strict total order + pdqsort | `mem_sort_dedup_patch`'s merge loop is order-sensitive; a total order picks different survivors among equal-`re` ties, moving `XS`/`MAPQ` and adding split alignments | [#257](https://github.com/fg-labs/bwa-mem3/pull/257) — restore bwa-mem2's `re`-only comparator + introsort outcome on the default path (the pdqsort/total-order sort moves behind `--fast`), plus [#261](https://github.com/fg-labs/bwa-mem3/pull/261) — run pdqsort only when the comparator sees no tie, recovering its speed *with* byte-identity |
| [#174](https://github.com/fg-labs/bwa-mem3/pull/174) — the `--meth` feature appended `int8_t meth_hypothesis` to `mem_chain_t`, growing it 48 → 56 B | klib's kbtree derives node fan-out from `sizeof(key_t)`; the changed fan-out reorders `.pos`-tied chains, regrouping seeds and flipping co-optimal primaries — even on **non-`--meth`** runs | [#268](https://github.com/fg-labs/bwa-mem3/pull/268) — repack `meth_hypothesis` into the existing bitfield word so `sizeof(mem_chain_t)` returns to 48 B (guarded by `static_assert`) |

The reductions were measured per-PR against a real bwa-mem2 v2.2.1 golden (5 M HG00096
WGS pairs, hg38, `-t 16` on both sides so batch boundaries and `mem_pestat` match):
[#256](https://github.com/fg-labs/bwa-mem3/pull/256) removed ~14 divergent primaries
(all `XS`); [#257](https://github.com/fg-labs/bwa-mem3/pull/257) took the residual from
1,501 (0.0150%) to 508 (0.0051%) and the extra supplementary alignments from +239 to +4;
and [#268](https://github.com/fg-labs/bwa-mem3/pull/268), measured on a build that
already carried `#256` and `#257`, took the remaining 54
divergent primaries on a 1,066,668-record batch-1 slice to **0** — every changed record
matching bwa-mem2 exactly, with zero regressions. A combined genome-wide, multi-sample,
cross-architecture concordance re-run to republish the aggregate figures below and
regenerate the [declared divergence catalog](#declared-divergence-catalog) is pending.

> **This is a compatibility choice, not a correctness claim.** Gap-from-`H` is textbook
> Gotoh and what ksw2/minimap2 do; gap-from-`M` is a deliberate bwa idiosyncrasy. `#256`
> chooses bwa-mem2 output compatibility over the convention because bwa-mem3 is a drop-in
> replacement. `IPNP-BIPN/bwa-mem4` made the same call (`990df3a`) for the same reason.

## What is preserved

The figures in this section were measured **before** the 0.7.1 parity restoration and so
capture the residual drift of the three now-reverted changes ([#141](https://github.com/fg-labs/bwa-mem3/pull/141),
[#123](https://github.com/fg-labs/bwa-mem3/pull/123), [#174](https://github.com/fg-labs/bwa-mem3/pull/174);
see [Restored bwa-mem2 parity (0.7.1)](#restored-bwa-mem2-parity-071)). They are
**pre-restoration measurements**, not a bound on current behavior: the 0.7.1 drop-in
profile measured tighter where it has been re-run (0 diverging primaries on the batch-1
WGS slice), and it has not been re-run on the other cells at all. We ran an empirical
concordance check with
[bwa-mem3-bench](https://github.com/fg-labs/bwa-mem3-bench) at commit `a02fcb4`, comparing
bwa-mem3 against upstream `bwa-mem2` v2.2.1 on x86 hosts across whole-genome, whole-exome,
and panel workloads. **Primary-alignment concordance** — reference name, position, CIGAR,
MAPQ, and placement flags compared per read end — was:

| sample | primary concordance (pre-0.7.1) | primary records |
|---|---:|---:|
| wes-5M | 99.9996% | 10,051,170 |
| wgs-5M | 99.9893% | 9,980,872 |
| panel-twist-5M | 99.9414% | 7,913,324 |

Even before restoration, where each read maps was preserved on essentially every read; the
well-under-0.1% of primary records that differed did so in `MAPQ`, `CIGAR`, or position, and
were accounted for by the deterministic tie-break change ([#123](https://github.com/fg-labs/bwa-mem3/pull/123))
and gap-from-`H`/`sizeof` changes now reverted for 0.7.1. (On the 1M-read `smoke-1M` cell the
pre-restoration figure was 99.946%.) Cross-architecture, the NEON (ARM) and x86 builds are
byte-identical to each other — 100.0000% concordance over all records, supplementary
alignments included — and that has held throughout.

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

### Additional supplementary alignments (reduced in 0.7.1, re-measurement pending)

Earlier releases emitted a small number of **additional supplementary (chimeric/split) alignments** that upstream `bwa-mem2` v2.2.1 did not (e.g. on `wes-5M`, 5,123 vs 5,118). This was the previously-unattributed [#127](https://github.com/fg-labs/bwa-mem3/issues/127) finding; it has since been root-caused to the [#123](https://github.com/fg-labs/bwa-mem3/pull/123) total-order dedup sort and **largely reduced by [#257](https://github.com/fg-labs/bwa-mem3/pull/257)**, which cut the extra split alignments from +239 to **+4** on the 5 M WGS cell by restoring bwa-mem2's dedup outcome on the default path. The primary alignment of every affected pair was unchanged throughout; only an extra supplementary record was ever added.

**That +4 is the current published figure, not zero.** It was measured on a build carrying [#256](https://github.com/fg-labs/bwa-mem3/pull/256) and [#257](https://github.com/fg-labs/bwa-mem3/pull/257) but *not* [#268](https://github.com/fg-labs/bwa-mem3/pull/268), and no supplementary-alignment count has been published for a build that carries `#268`. `#268`'s own measurement covered *primary* divergence on a 1,066,668-record batch-1 slice, which is a different quantity. Until a 0.7.1 supplementary count lands, treat the residual +4 as outstanding and [`FG-SUPP-ADDITIONS`](#declared-divergence-catalog) as unresolved.

(These additions were never a product of the opt-in `--supp-rep-hard-cap` rescoring, which only lowers MAPQ and never adds records.) The `--fast` profile, which keeps the [#123](https://github.com/fg-labs/bwa-mem3/pull/123) sort, still emits them by design.

### Divergences that are latent, opt-in, or per-architecture

The following changes can move alignments, scores, or MAPQ relative to upstream, but did **not** surface as primary-alignment differences on the measured cells because they are gated, latent, or only active on other inputs or architectures:

- **Proper-pair `FLAG` recompute ([#17](https://github.com/fg-labs/bwa-mem3/pull/17), default-on).** bwa-mem3 computes the `0x2` bit from the alignment actually emitted rather than from the below-threshold primary. This only changes the flag in the rare case where the primary's score is under `opt->T` but an ALT hit clears it; on `smoke-1M` no record hit that path, so the full `FLAG` matched upstream exactly. See [Correctness fixes → Proper-pair flag](correctness.md).
- **SIMD scoring-kernel fixes ([#21](https://github.com/fg-labs/bwa-mem3/pull/21), [#26](https://github.com/fg-labs/bwa-mem3/pull/26), [#28](https://github.com/fg-labs/bwa-mem3/pull/28), [#29](https://github.com/fg-labs/bwa-mem3/pull/29), [#30](https://github.com/fg-labs/bwa-mem3/pull/30), [#31](https://github.com/fg-labs/bwa-mem3/pull/31)).** These correct the batched mate-rescue `kswv` kernels so the suboptimal score (`score2` → `XS:i`/`MAPQ`) converges toward the scalar `ksw_align2` reference. They move `XS`/`MAPQ` on the minority of reads where the SIMD kernel previously diverged, and the affected reads differ by architecture (AVX2 vs NEON vs AVX-512BW). See [Correctness fixes → kswv score2 plateau series](correctness.md).
- **Seeding correctness fixes ([#55](https://github.com/fg-labs/bwa-mem3/pull/55), [#73](https://github.com/fg-labs/bwa-mem3/pull/73), [#100](https://github.com/fg-labs/bwa-mem3/pull/100)).** These fix buffer sizing and a prefetch-mask precedence bug. They change alignments only where the old bug actually triggered (e.g. reads longer than 151 bp for [#55](https://github.com/fg-labs/bwa-mem3/pull/55); [#73](https://github.com/fg-labs/bwa-mem3/pull/73) is a prefetch hint with no semantic change).
- **Opt-in MAPQ rescoring ([#56](https://github.com/fg-labs/bwa-mem3/pull/56), [#101](https://github.com/fg-labs/bwa-mem3/pull/101), [#118](https://github.com/fg-labs/bwa-mem3/pull/118), default-off).** `--supp-rep-hard-cap INT` forces MAPQ=0 on supplementary alignments anchored in repetitive seeds. With no flag the output is unchanged; [#101](https://github.com/fg-labs/bwa-mem3/pull/101) makes the flag actually take effect (it shipped as a silent no-op before), and [#118](https://github.com/fg-labs/bwa-mem3/pull/118) is its regression test. See [Features → `--supp-rep-hard-cap`](features.md).
- **Tie-break determinism ([#123](https://github.com/fg-labs/bwa-mem3/pull/123), reverted on the default path in 0.7.1).** `#123`'s strict total order + pdqsort at the dedup-patch sort sites reordered equal-scoring ties relative to bwa-mem2's partial-order outcome. As of 0.7.1 the default path restores bwa-mem2's comparator + introsort outcome ([#257](https://github.com/fg-labs/bwa-mem3/pull/257)) and recovers pdqsort's speed only on tie-free inputs ([#261](https://github.com/fg-labs/bwa-mem3/pull/261)); `#123`'s behavior now lives behind `--fast`, where it remains **not** byte-identical.

(The resolve→order→chain seeding refactor that backs `--seed-order` is byte-identical in its default `off` mode; it is described in the dedicated section below rather than listed here, since only its non-`off` modes are divergent.)

### Seed ordering (`--seed-order`, opt-in)

`--seed-order off` (the default) preserves byte-identical output. The internal resolve→order→chain refactor that enables it was verified to be bit-for-bit identical to the previous code path across single-end, paired-end, and threaded runs.

`--seed-order local-longest` (and the unadvertised advanced modes `global-longest`, `absorb-count`, `most-absorb`) are **opt-in** and are **not byte-identical**. They reorder each read's SA-resolved seeds before chaining so that the longest seeds anchor chains first; contained shorter seeds are then absorbed rather than extended. The output shift is non-trivial: secondary alignments, `XA:Z:`, `XS:i`, and `HN:i` can all change, and a small number of primary alignments may shift as well.

Accuracy on an easy simulated profile (holodeck, ~94.4 % F1) is flat relative to `off`. Hard-data F1 validation on divergent/indel-rich reads and GIAB benchmarks is **not yet complete**. Because the accuracy gate is not fully passed, all non-`off` modes remain opt-in only; the default stays `off`. See [Optimization checklist → Reorder seeds longest-first](../best-practices/optimization-checklist.md#6-reorder-seeds-longest-first---seed-order-local-longest) for usage guidance.

## Declared divergence catalog

The divergences described above are tracked as a structured registry in [bwa-mem3-bench](https://github.com/fg-labs/bwa-mem3-bench) (`docs/expected-divergences.yaml`). Each carries a per-sample concordance-drift budget that the benchmark gates against on every run — a new bwa-mem3 build that drifts beyond its budget fails CI rather than silently shipping a regression. The table below is generated from that registry; do not edit it by hand (see [bwa-mem3-bench → Per-release concordance history](../related-projects/bwa-mem3-bench.md) for how it is regenerated).

> **⚠️ This catalog predates the 0.7.1 parity restoration and is pending regeneration.**
> Its `FG-PRIMARY-DRIFT` row still attributes primary drift to [#123](https://github.com/fg-labs/bwa-mem3/pull/123)
> (and to per-architecture SIMD `score2`/MAPQ convergence — an attribution the
> measurements so far do not confirm, since the NEON and x86 builds came out
> byte-identical to each other, and which remains unchecked on AVX-512BW pending the
> cross-architecture re-run); on the default path that drift is **resolved in 0.7.1**
> ([#256](https://github.com/fg-labs/bwa-mem3/pull/256), [#257](https://github.com/fg-labs/bwa-mem3/pull/257),
> [#261](https://github.com/fg-labs/bwa-mem3/pull/261), [#268](https://github.com/fg-labs/bwa-mem3/pull/268))
> on the slice re-measured so far. `FG-SUPP-ADDITIONS` tracks the [#127](https://github.com/fg-labs/bwa-mem3/issues/127)
> extra supplementaries and **remains unresolved**: [#257](https://github.com/fg-labs/bwa-mem3/pull/257)
> reduced them from +239 to +4 and no count has been published for a build carrying
> [#268](https://github.com/fg-labs/bwa-mem3/pull/268).
> The registry in `bwa-mem3-bench` must be re-run against a 0.7.1 build and its budgets
> tightened; until then, read the table as the pre-0.7.1 state, not current behavior.

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
