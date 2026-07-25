# Equivalence with bwa-mem2 (bit-identity)

<div class="warning" style="border:3px solid #2e7d32;border-radius:8px;padding:1rem 1.25rem;margin:1.25rem 0;background:#e8f5e9;color:#1b1b1b;">
<strong style="font-size:1.35em;display:block;margin-bottom:0.4rem;">✅ No alignment record bwa-mem2 emits is changed or removed — and on the cell re-measured for 0.7.1 its primary alignments are reproduced exactly.</strong>
On the drop-in profile (plain <code>bwa-mem3 mem</code>, no flags), release 0.7.1 gives the <strong>same primary alignment as bwa-mem2</strong> for <strong>every</strong> read of the slice re-measured after the parity work landed (1,066,668 records, HG00096 WGS, hg38, x86): identical <code>FLAG</code>, <code>RNAME</code>, <code>POS</code>, <code>MAPQ</code>, <code>CIGAR</code>, <code>AS</code>, and <code>XS</code>. Those record bytes differ by <strong>two additive tags</strong> and nothing else:
<ul style="margin:0.5rem 0 0 0;">
<li><code>MQ:i</code> — mate mapping quality (an extra tag bwa-mem2 never wrote)</li>
<li><code>HN:i</code> — hit count (an extra tag bwa-mem2 never wrote)</li>
</ul>
Strip those two tags and those records are byte-for-byte identical. The <strong>header</strong> is a separate question that these measurements do not cover: the <code>@PG</code> line names <code>bwa-mem3</code>, and bwa-mem3 emits a default <code>@HD</code> line that bwa-mem2 has no code to write at all (<a href="https://github.com/fg-labs/bwa-mem3/issues/288">#288</a>). Every byte-identity claim on this page is about <em>alignment records</em>, not the header block.
<br><br>
<strong>The supplementary gap is now closed.</strong> The residual <strong>+4 supplementary alignments</strong> has been re-measured on <code>e722ed0</code> (a build carrying <a href="https://github.com/fg-labs/bwa-mem3/pull/268">#268</a>) and is <strong>+0</strong>: on <code>wgs-5M</code>, <code>wes-5M</code> and <code>hic-1M</code> the full <em>alignment-record</em> stream is byte-identical to bwa-mem2 v2.2.1 once <code>MQ:i</code>/<code>HN:i</code> are stripped — 22.5 M records including 436 k supplementary alignments. That comparison covers alignment records only, not the header block. What remains open is the genome-wide, multi-sample re-run that would extend this past the cells named on this page. Treat every claim here as scoped to the evidence it cites, not as a guarantee across all inputs and architectures. (The opt-in <code>--fast</code> speed levers are a separate, deliberate deviation — see below.)
</div>

As of **release 0.7.1**, the bwa-mem3 **drop-in profile** (plain `bwa-mem3 mem`, no
extra flags) produces the **same primary alignment** as bwa-mem2 — same `FLAG`,
`RNAME`, `POS`, `MAPQ`, `CIGAR`, `AS`, `XS` — on the data we have re-measured since
the restoration. The **SAM byte stream** is still not a literal byte-for-byte
match: bwa-mem3 emits extra auxiliary tags (`MQ:i`, `HN:i`) that upstream never wrote.
Strip those two tags and the *records* are identical. The header is not covered by that
statement — the `@PG` line names `bwa-mem3`, and the default `@HD` bwa-mem3 writes has no
bwa-mem2 counterpart ([#288](https://github.com/fg-labs/bwa-mem3/issues/288)). This is a change in
posture from earlier releases, which carried three changes that *did* move primary
alignments; all three were reverted for 0.7.1 (see [Restored bwa-mem2 parity
(0.7.1)](#restored-bwa-mem2-parity-071) below).

Two things scope that statement, and both are tracked below rather than assumed away.
The **supplementary** set is now confirmed to match: the extra split alignments were cut
from +239 to +4 by [#257](https://github.com/fg-labs/bwa-mem3/pull/257), and re-measurement on
`e722ed0` (which carries [#268](https://github.com/fg-labs/bwa-mem3/pull/268)) puts it at **+0**
(see [Additional supplementary
alignments](#additional-supplementary-alignments-resolved)).
And the post-restoration evidence base is still narrow, in two different ways: the
**0-divergent-primaries** count comes from a single cell (the 1,066,668-record batch-1 WGS
slice, x86), while the **byte-identity** result spans three cells (`wgs-5M`, `wes-5M`,
`hic-1M`, x86) with a `c6a`/`c8g` cross-architecture check on top. The genome-wide,
multi-sample concordance re-run that would extend either past those cells is still pending.
Read every figure on this page as scoped to the cell it names.

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
have measured, and on the cells re-measured after the 0.7.1 restoration it is
preserved exactly. On those records the **SAM byte stream** differs only by the additive
`MQ:i`/`HN:i` tags; the header block is a separate matter (the `@PG` line, and the default
`@HD` — [#288](https://github.com/fg-labs/bwa-mem3/issues/288)). The **supplementary set** is re-confirmed at +0 on `e722ed0`. The **opt-in speed levers** (`--fast` and its constituent flags) are
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
**pre-restoration measurements**, not a bound on current behavior. The 0.7.1 drop-in profile
measured tighter everywhere it has been re-run: 0 diverging primaries on the batch-1 WGS
slice, and a byte-identical complete alignment-record stream (`MQ:i`/`HN:i` stripped) on `wgs-5M`,
`wes-5M` and `hic-1M` (see
[Additional supplementary alignments](#additional-supplementary-alignments-resolved)) — which
supersedes the `wes-5M` and `wgs-5M` rows below. `panel-twist-5M` has been re-run only for the
cross-architecture check, and `smoke-1M` has not been re-run at all, so those two figures still
stand as written. We ran an empirical concordance check with
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
pre-restoration figure was 99.946%.) Cross-architecture, the NEON (ARM) and x86 builds came
out byte-identical to each other on these cells — 100.0000% concordance over all records,
supplementary alignments included. Read that as scoped to the cells and builds it was measured
on rather than a standing property: a later `panel-twist-5M` run found exactly one record
differing between `c6a` and `c8g` (see [the cross-architecture
exception](#additional-supplementary-alignments-resolved) below).

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

Separately, the `@PG` header line reports `ID:bwa-mem3` / `PN:bwa-mem3` rather than `bwa-mem2`, which is also a byte-level header difference by design. It is not the only one: bwa-mem3 emits a default `@HD` line, which bwa-mem2 v2.2.1 has no code to write at all, so the header block differs by more than `@PG` alone ([#288](https://github.com/fg-labs/bwa-mem3/issues/288)). None of the record-level byte-identity results above depend on the header — every one of them compares alignment records only.

### Additional supplementary alignments (resolved)

Earlier releases emitted a small number of **additional supplementary (chimeric/split) alignments** that upstream `bwa-mem2` v2.2.1 did not (e.g. on `wes-5M`, 5,123 vs 5,118). This was the previously-unattributed [#127](https://github.com/fg-labs/bwa-mem3/issues/127) finding; it has since been root-caused to the [#123](https://github.com/fg-labs/bwa-mem3/pull/123) total-order dedup sort and **largely reduced by [#257](https://github.com/fg-labs/bwa-mem3/pull/257)**, which cut the extra split alignments from +239 to **+4** on the 5 M WGS cell by restoring bwa-mem2's dedup outcome on the default path. The primary alignment of every affected pair was unchanged throughout; only an extra supplementary record was ever added.

**That +4 is now measured at +0 on a `#268`-carrying build.** Re-measured 2026-07-24 against
bwa-mem2 v2.2.1 on `e722ed0`, comparing the *complete* record stream (all 11 mandatory fields and
every optional tag, byte-for-byte) with only the additive `MQ:i`/`HN:i` removed:

| cell | records | supplementary (bwa-mem2 / bwa-mem3) | result |
|---|---|---|---|
| `wgs-5M` | 10,030,558 | 49,686 / 49,686 | byte-identical |
| `wes-5M` | 10,056,288 | 5,118 / **5,118** | byte-identical |
| `hic-1M` | 2,381,418 | 381,418 / 381,418 | byte-identical |

`hic-1M` is the strongest of the three for this question — it is the cell chosen in
`config/samples.yaml` precisely because chimeric Hi-C ligation products generate heavy
supplementary load, and it shows zero disagreement across 381,418 supplementary records. The
earlier `wes-5M` figure of 5,123 vs 5,118 no longer reproduces; both are 5,118.

Two independent corroborations landed with it. **Cross-architecture:** bwa-mem3 `c6a` (AVX2) and
`c8g` (NEON) are byte-identical on all three cells with *nothing* stripped — `MQ`/`HN` included.
**A third aligner:** `bwa` v0.7.19 was run on `hic-1M` and `wes-5M` (`c6a`, `-t 16`, matching
`mem_flags`) and agrees with both bwa-mem2 and bwa-mem3 over 12.4 M records. That agreement is
byte-identical only *after* the same normalization applied above — alignment records only, since
the `@PG` line names a different program in each of the three outputs, and with the additive
`MQ:i`/`HN:i` stripped from the bwa-mem3 side. The raw SAM streams are **not** byte-identical;
the normalized records are, so the agreement is not an artifact of the shared bwa-mem2 lineage.

`FG-SUPP-ADDITIONS` can be marked resolved for these cells. Scope: three cells on x86, plus the
cross-arch check; it is not a genome-wide guarantee.

**One cross-architecture exception, tracked separately.** On `panel-twist-5M`, `c6a` and `c8g`
differ on exactly **1 of 8,100,270 records** — a phantom `XS` from unzeroed query padding in the
8-bit `kswv` mate-rescue kernels on NEON and AVX-512BW, root-caused in
[#290](https://github.com/fg-labs/bwa-mem3/pull/290), where the fix is still open — so 0.7.1
still carries it. Those two kernels are the entire divergent
set: the AVX2 8-bit kernel and all three 16-bit kernels kept the correct per-cell padding test,
and the SSE tiers have no batched `kswv` kernel at all (mate rescue there runs the scalar
`ksw_align2` path). The split therefore falls between AVX2 and AVX-512BW/NEON, not between x86
and ARM.

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
> and to per-architecture SIMD `score2`/MAPQ convergence. That second mechanism did not surface
> on the three cells re-measured for 0.7.1 — NEON and x86 came out byte-identical there — but it
> is not hypothetical: `panel-twist-5M` turned up one record where the AVX-512BW and NEON 8-bit
> `kswv` kernels produced a phantom `XS` that AVX2 did not
> ([#290](https://github.com/fg-labs/bwa-mem3/pull/290), fix still open). On the default path the
> `#123` drift is **resolved in 0.7.1**
> ([#256](https://github.com/fg-labs/bwa-mem3/pull/256), [#257](https://github.com/fg-labs/bwa-mem3/pull/257),
> [#261](https://github.com/fg-labs/bwa-mem3/pull/261), [#268](https://github.com/fg-labs/bwa-mem3/pull/268))
> on the slice re-measured so far. `FG-SUPP-ADDITIONS` tracks the [#127](https://github.com/fg-labs/bwa-mem3/issues/127)
> extra supplementaries; its row below has been **hand-updated ahead of regeneration** to carry
> the [+0 result](#additional-supplementary-alignments-resolved) measured on a
> [#268](https://github.com/fg-labs/bwa-mem3/pull/268)-carrying build, which is why it reads
> `RESOLVED` where the registry still says `TBD`. That resolution is scoped to the three
> measured cells (`wgs-5M`, `wes-5M`, `hic-1M`, on x86); the remaining samples have not been
> re-measured. Its `samples` and `budget_%` columns are inherited unchanged from the
> un-regenerated registry, where `samples = all` and `budget_% = 0.0000` mean this class
> contributes nothing to the *primary*-concordance drift budget — they are not a validated
> zero-tolerance gate across every sample.
> The registry in `bwa-mem3-bench` (`docs/expected-divergences.yaml`, which still carries
> `pr: TBD`, `samples: []`, and `expected_drift_pct: 0.0` for this entry) must be updated and
> re-run against a 0.7.1 build and its budgets tightened. Until that lands, only the
> `FG-SUPP-ADDITIONS` row's `pr` and `summary` cells reflect current behavior — they are the
> hand-update. Read the rest of the table as the pre-0.7.1 state: the other two rows in full,
> and that row's `samples` and `budget_%` columns.

<!-- FG-DIVERGENCE-CATALOG:start -->
| id | pr | affected | samples | budget_% | summary |
| --- | --- | --- | --- | --- | --- |
| FG-PRIMARY-DRIFT | fg-labs/bwa-mem3#123 | primary_alignment | wgs-5M, wes-5M, panel-twist-5M, smoke-1M | 0.1000 | Per-architecture SIMD score2/MAPQ convergence (#21, #26, #28-#31) and deterministic tie-break ordering (#123) shift MAPQ, CIGAR, or position on a small fraction of primary alignments relative to bwa-mem2 v2.2.1. Where each read maps is preserved; the affected reads differ in placement detail, and the set varies by SIMD architecture. |
| FG-METH-DIVERGENCE | fg-labs/bwa-mem3#90 | meth_alignment | meth-twist-emseq-5M, smoke-meth | 1.5000 | Bisulfite (--meth) mode against the bwameth.py baseline diverges beyond the ignored YD/XM/XG tag set (Bismark-compatible XR/XG/XM tags and C->T/G->A conversion handling), giving a larger but still-bounded concordance drift on methylation workloads. |
| FG-SUPP-ADDITIONS | RESOLVED | supplementary_alignment | all | 0.0000 | Earlier builds emitted additional supplementary (split/chimeric) alignments vs bwa-mem2 v2.2.1 (e.g. wes-5M: 5123 vs 5118). Re-measured 2026-07-24 on `e722ed0`: **+0** on wgs-5M, wes-5M and hic-1M (x86; the other samples have not been re-measured), with the complete alignment-record stream byte-identical (MQ/HN stripped; header lines not compared) across 22.5 M records incl. 436 k supplementary. Primary alignments were unchanged throughout. |
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
