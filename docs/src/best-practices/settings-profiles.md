# Settings profiles: bwa drop-in vs recommended

bwa-mem3's **defaults track bwa-mem / bwa-mem2 behavior** — it is a faster drop-in for an existing
bwa-mem2 pipeline. Some defaults the bwa lineage inherited are, however, more conservative than
necessary. This page defines two profiles:

- a **drop-in profile** that keeps the bwa-mem/bwa-mem2 defaults, for migrations and parity checks;
- a **recommended profile** that deviates from those defaults where we have benchmarked the change
  to be near-neutral on accuracy and clearly faster.

The defaults ship as the drop-in profile. The recommended profile is **opt-in** — you turn it on
with explicit flags — so upgrading bwa-mem3 never silently changes your alignments.

> For the higher-level picture — how the drop-in default, `--compat`, and the `--fast` preset
> differ in the alignments they produce — see
> [Alignment modes: plain, `--compat`, `--fast`](../whats-different/modes.md). This page is the
> detail on the individual recommended levers.

**Which profile?**

| Your situation | Profile | Invocation |
|---|---|---|
| Migrating a bwa-mem2 pipeline, or validating against bwa/bwa-mem2 | **Drop-in** | `bwa-mem3 mem` (no extra flags) |
| New pipeline, or migration already validated | **Recommended** | `bwa-mem3 mem -m 10 -y 0` (add `-s 2` under `--meth`) |

> **As of 0.7.1 the drop-in profile produces the same primary alignment as bwa-mem2 on every
> cell re-measured since the parity restoration** — identical
> `FLAG`/`RNAME`/`POS`/`MAPQ`/`CIGAR`/`AS`/`XS` on all 1,066,668 records of the HG00096 WGS
> slice re-run on x86. Those record bytes differ by **two additive tags** and nothing else: the
> extra `MQ:i` and `HN:i`. No alignment record bwa-mem2 emits is changed or
> removed; strip those two tags and the records are byte-for-byte identical. The **header** is
> not covered by that statement — the `@PG` line names `bwa-mem3`, and bwa-mem3 writes a default
> `@HD` that bwa-mem2 never emits ([#288](https://github.com/fg-labs/bwa-mem3/issues/288)).
> The three earlier changes that *did* move primaries were reverted for 0.7.1
> ([#256](https://github.com/fg-labs/bwa-mem3/pull/256), [#257](https://github.com/fg-labs/bwa-mem3/pull/257)/[#261](https://github.com/fg-labs/bwa-mem3/pull/261),
> [#268](https://github.com/fg-labs/bwa-mem3/pull/268)).
>
> **Scope:** the earlier residual of **+4 supplementary alignments** has since been
> re-measured on a build carrying [#268](https://github.com/fg-labs/bwa-mem3/pull/268) and is
> **+0** on `wgs-5M`, `wes-5M` and `hic-1M` (x86), with the complete *alignment-record* stream
> byte-identical once `MQ:i`/`HN:i` are stripped — records only; the header block was not part
> of that comparison. The remaining samples have not been re-measured, and the genome-wide,
> multi-sample, cross-architecture re-run is still pending —
> so validate your own workload before treating parity as a guarantee. See
> [Equivalence with bwa-mem2](../whats-different/equivalence.md) for the full audit trail. The
> **recommended** profile (and `--fast`) is a further, deliberate, documented deviation — it is
> *not* byte-identical.

## Drop-in profile (default)

```bash
bwa-mem3 mem -t <N> ref.fa R1.fq R2.fq > out.sam
```

No extra flags. Use this when you are:

- reproducing or validating against bwa-mem / bwa-mem2 output, or
- migrating an existing bwa-mem2 pipeline and want to change one variable (the aligner binary) at a
  time before tuning anything else.

> **Add `-K` if you are comparing against bwa/bwa-mem2, or need reproducible output.**
> The default batch size is `chunk_size × -t`, so it moves with the thread count, and batch
> boundaries determine each batch's `mem_pestat` insert-size cohort — which feeds pairing, mate
> rescue and MAPQ. `bwa` and `bwa-mem2` use the identical formula, so this is inherited behaviour,
> not a bwa-mem3 difference; but it does mean a drop-in comparison is only apples-to-apples when
> both sides see the same batching. Pass the **same `-K` to both binaries** (or match the `-t`):
>
> ```bash
> bwa-mem3 mem -t <N> -K 100000000 ref.fa R1.fq R2.fq > mem3.sam
> bwa-mem2 mem -t <M> -K 100000000 ref.fa R1.fq R2.fq > mem2.sam   # any -t; -K pins the batching
> ```
>
> See [Aligning → `-K`](../user-guide/aligning.md) for the full explanation.

## Recommended profile

```bash
bwa-mem3 mem -t <N> -m 10 -y 0 ref.fa R1.fq R2.fq > out.sam
```

> **Shorthand:** `bwa-mem3 mem --fast` applies `-m 10 -y 0 --min-ext-len 30
> --smem-dedup --skip-contained-ext --max-extend-chains 20 --adaptive-band
> --extend-mate-concordant` (plus `-s 2` and a lower `--max-extend-chains 10`
> under `--meth`) in one flag. Explicit flags
> still override individual levers where applicable;
> [`--smem-dedup`](https://github.com/fg-labs/bwa-mem3/pull/187) and
> [`--skip-contained-ext`](https://github.com/fg-labs/bwa-mem3/pull/192) are forced on with no opt-out, and
> [`--adaptive-band`](https://github.com/fg-labs/bwa-mem3/pull/194) is
> forced on but can be opted back out with `--no-adaptive-band` (which restores the exact
> full-width extension step — byte-identical to a run without `--adaptive-band` — while
> keeping the rest of `--fast`, whose other levers can still change output)
> (`--adaptive-band` is a no-op on short reads, a ~25% speedup on medium-length runs
> such as SBX ~240 bp; kilobase-scale HiFi/ONT do not run at default settings — see
> [Situational: `--adaptive-band`](#situational---adaptive-band-for-medium-length-reads)).
> `--skip-contained-ext` no-ops under `--meth` (its own internal gate disables it
> there), so on a `--meth` run the effective levers are
> `-m 10 -y 0 --min-ext-len 30 --smem-dedup --max-extend-chains 10 --adaptive-band -s 2 --extend-mate-concordant`.
> See [`mem` → `--fast`](../cli/mem.md#--fast--speed-preset-opt-in-not-byte-identical).

Use this for new pipelines, or once a drop-in migration is validated and you want bwa-mem3's best
speed/accuracy trade-off. Current recommended deviations:

| flag | default (drop-in) | recommended | effect |
|---|---|---|---|
| `-m` (mate-rescue depth) | 50 | **10** | ~11–22% less alignment CPU; near-neutral accuracy (see below) |
| `-y` (3rd-round seeding occurrence) | 20 | **0** | ~11–30% less alignment CPU; F1 near-neutral across regimes (within ±0.02; better on divergent/repeat — see below) |
| `-s` (Pass-2 re-seed width), **`--meth` only** | 10 | **2** | light re-seed: ~same speed as `-s 0` but recovers the MAPQ/placement `-s 0` lost (see below) |
| `--min-ext-len` (skip short-seed extension), **standard-error reads** | 0 | **30** | ~10–20% less alignment CPU; accuracy change confined to the already-low-confidence tail (see below) |
| `--skip-contained-ext` (skip contained-seed extension), **short/medium non-meth reads** | off | **on** | ~10% less alignment CPU; **byte-identical** — zero accuracy change (see below) |

This table will grow as we benchmark additional tunings; each entry is gated on the same
"measurably faster, near-neutral accuracy" bar.

For a bisulfite (`--meth`) pipeline the recommended invocation is therefore:

```bash
bwa-mem3 mem -t <N> --meth -m 10 -s 2 -y 0 ref.fa R1.fq R2.fq > out.sam
```

## Mate-rescue depth: `-m 10`

`-m` caps how many near-best candidate loci per read get a mate-rescue Smith–Waterman pass. bwa-mem
and bwa-mem2 both default to **50**; we measured the marginal value of that depth directly —
end-to-end and at the read level against simulated golden truth.

**The depth beyond ~10 is a measured wash.** Mate rescue earns its keep on the first few candidates,
but candidates 11–50 only ever engage on reads with many near-best loci — repetitive and
cross-contig placements — and on truth data deep rescue places about **as many of those reads
correctly as it mis-places**. Lowering `-m` 50 → 10 changes which low-confidence reads win, but the
**net true-recall cost is ≈ 0.003–0.004%** (≈200 reads in 5.4 M), on both methylation and
non-methylation data, and the accuracy-vs-depth curve is flat (no depth between 10 and 50 is
meaningfully better). Disabling rescue entirely (`-S`, or `-m 0`) *is* a real regression — so rescue
matters, just not 50-deep.

**In exchange you get a real speed-up**, largest exactly where deep rescue binds (high-depth
amplicon/panel and repetitive regions):

| dataset | alignment CPU | wall |
|---|---:|---:|
| panel-twist (high depth) | **−22%** | −21% |
| wgs | −19% | −19% |
| wes | −18% | −14% |
| meth (em-seq) | −11% | −9% |

Aggregate `samtools flagstat` impact is a uniform, directional sub-0.2 pp reduction in
mapped%/proper-pair% (≤0.02 pp on WGS/WES; ≤0.12–0.19 pp on amplicon/meth), entirely in the
low-MAPQ deep-rescue tail: the mapped count only ever *decreases* (it never newly maps a
previously-unmapped read), and the handful of reads whose placement does change net to ≈0 recall
against golden truth.

One caveat: this golden-truth metric scores a single best placement per read, so it does not capture downstream tools that *aggregate over repeats* — depth-based CNV, SV/mobile-element calling in repetitive regions, or repeat-region methylation. If your pipeline relies on which repeat copy is reported (rather than just the confident, MAPQ-high tier), validate `-m 10` against your own analysis rather than assuming neutrality.

Numbers are from [bwa-mem3-bench](../related-projects/bwa-mem3-bench.md) on 5 M-read real datasets
plus a multi-contig holodeck golden-truth ablation; consult the bench for methodology and current
figures.

## Third-round seeding: `-y 0`

bwa-mem seeds in up to three rounds: (1) SMEMs, (2) reseeding long SMEMs (`-r`), and (3) an
occurrence-bounded "seed strategy" round (`-y`) that, for every read position, grows an exact match
until it occurs fewer than `-y` (default 20) times in the genome and emits it. Round 3 is a
repeat-region safety net. `-y 0` disables it entirely.

**It is net-useless-to-harmful, even in the repeats it was built for.** We swept `-y 0` (and the
more-aggressive `-y 0 -r 10`, see below) against the stock default across four golden-truth regimes
(holodeck, read-name truth) and on real WES + WGS:

| regime | ΔF1 (`-y 0` − default) | ΔF1 (`-y 0 -r 10` − default) | ΔF1 in the MAPQ-0 (repeat) bin (`-y 0` / `-y 0 -r 10`) |
|---|---|---|---|
| easy (150 bp, default error) | −0.010 | −0.001 | −0.14 / −0.04 |
| substitution-divergent (2–15%) | **+0.008** | −0.074 | **+0.13** / −0.07 |
| indel-rich (37 % indels) | −0.016 | −0.017 | −0.02 / −0.03 |
| **repeat-enriched** (reads simulated from RepeatMasker, 49 % of hg38) | **+0.007** | +0.004 | **+0.16** / −0.02 |

The tiny easy/indel deltas are ≤0.016 in absolute F1; on the divergent and *repeat-enriched* regimes —
where round 3 is supposed to earn its keep — removing it makes F1 **better**. Mechanism: its
low-occurrence (<20-hit) seeds add spurious near-tied repeat candidates that slightly degrade
placement; dropping them lets the correct unique anchor win. (The regime sweep was on standard,
non-`--meth` reads; `-y 0` is a generic seeding change and is recommended for `--meth` on the same
basis, but was not separately F1-measured there.)

**On real data the confident core is untouched.** On 20 M real WES+WGS reads, **zero** confidently and
uniquely mapped reads (MAPQ ≥ 60, NM ≤ 1, no soft-clip) became unmapped under `-y 0`, and 3–4 per 10 M
lost any alignment score. Every regression lands on the already-low-confidence tail (base MAPQ ≈ 0,
heavily soft-clipped or high-NM), and round 3 was so often producing *worse* primaries that removing it
*improves* a large share of as many reads as it perturbs (on WES, 3,194 reads improve vs 3,784 that
worsen).

**In exchange, alignment CPU drops ~11–30 % single-thread**, larger on data with more repeat content
(WGS > WES). The win is confirmed cross-architecture (Graviton4/Linux: realistic +29.9 %, HG002 WGS
+19.8 %, matching macOS) — it cuts cold-memory seeding work, not arithmetic, so it ports.

**More aggressive (clean-data only): also drop round 2 with `-y 0 -r 10`.** That roughly halves
alignment CPU (~50–63 %), but round 2 *is* genuine split-read/divergence sensitivity — it costs
−0.074 F1 on divergent data (the `-y 0 -r 10` column above, concentrated in the repeat/multimapper
bin). Use `-y 0 -r 10` only on known-clean, low-divergence libraries; `-y 0` alone is the
broadly-safe recommendation.

## Pass-2 re-seeding under `--meth`: `-s 2`

> **`--fast --meth` uses `-s 2` (light Pass-2 re-seed), not `-s 0`.** Earlier releases set `-s 0`
> (no re-seed). Read-level analysis on holodeck sim reads showed `-s 0` **inflates MAPQ**: the affected
> reads map on occurrence-1 SMEMs that hide an *interior* repeat only Pass-2 would surface, so without
> it they look uniquely placed and MAPQ is pushed toward 60 — the calibration caveat noted below, now
> quantified. `-s 2` re-seeds exactly those occurrence-1 SMEMs (the cheap subset), recovering MAPQ
> **and** placement at ≈ the speed of `-s 0` (3.6× vs 3.25× on twist-em-seq 5M; MAPQ matches the `-s 10`
> default on 85/86 changed reads; placement 97.6% = default). This is the "cheap middle ground" the
> read-confidence fallback at the end of this section did not find (it re-seeds by SMEM *occurrence*,
> not by read confidence). **Validation status:** the `-s 0` recall/speed figures below are the
> genome-wide em-seq bench; the `-s 2` MAPQ/placement recovery is read-level + chr1 sim, with a
> genome-wide em-seq re-run pending.

`-s` controls bwa-mem's Pass-2 "re-seeding" — after the first seeding pass, long super-maximal exact
matches whose occurrence count is `≤ -s` are re-seeded from their midpoint to recover shorter, more
frequent matches that a long seed would otherwise swallow. Under `--meth`, reads are projected into a
3-letter alphabet (C→T and G→A), which collapses sequence complexity: in that space Pass-2 mostly
mines low-complexity sub-seeds that inflate the candidate set without changing placement. Setting
`-s 0` disables Pass-2 entirely (an exact seed's interval size is always ≥ 1, so the re-seed gate
never fires). This is **methylation-specific** — in the full 4-letter alphabet Pass-2 still earns its
keep, so the recommendation is scoped to `--meth`.

**Placement is a measured wash, even where re-seeding should matter most.** We aligned ~33 M
holodeck em-seq golden-truth reads with `-s 0` vs the `-s 10` default across two references —
single-contig chr22 and a deliberately hard multi-contig build (chr19–22 + their ALT scaffolds +
decoys + HLA, 2,977 contigs) that stresses cross-contig multi-mapping. In both regimes the net
true-recall difference is **within noise**:

| reference (10 seeds) | records | net recall Δ (`-s 0` − `-s 10`) | significance |
|---|---:|---:|---|
| single-contig chr22 | 16.9 M | −22 (−0.0001%) | McNemar p = 0.87 |
| multi-contig (ALT/decoy/HLA) | 16.4 M | −453 (−0.0028%) | McNemar p = 0.23 |

The direction (a sub-0.003% loss concentrated in the low-MAPQ cross-contig tail) mirrors `-m 10`'s,
and `-s 0` is *not* the same as disabling all seeding — it only removes the redundant second pass.

**The one caveat is mapping-quality calibration.** Removing Pass-2 shrinks the candidate set, so MAPQ
trends upward (≈95% of changed reads). On easy references this is benign re-distribution; in the hard
multi-contig regime the dominant high-confidence bin (MAPQ 50–60, ~84% of reads) stays *identically*
calibrated (0.020% mis-placement either way), but the small mid-MAPQ bins (~3% of reads) are modestly
to ~2× more error-prone at a given reported MAPQ under `-s 0`. This is negligible for methylation
quantification but worth knowing for callers that hard-filter on MAPQ in repetitive regions.

**In exchange, seeding is ~20% cheaper.** On chr22 (single thread, best-of-5): alignment user-CPU
**−19.7%** (12.96 s vs 16.14 s) and peak RSS **−6.5%**, with the gain largest on high-depth and
repetitive inputs where the redundant Pass-2 candidate set is largest.

A selective fallback — skip Pass-2 globally but re-seed only low-confidence reads — was evaluated and
**does not help**: in the multi-contig regime it would re-seed ~13% of reads yet recovers about as
many placements as it sacrifices (net within noise), so there is no cheap middle ground.

## Short-seed extension: `--min-ext-len 30`

`--min-ext-len INT` skips banded Smith–Waterman extension of seeds shorter than `INT` bp **in chains
that still hold a longer anchor seed** — those short seeds are collinear with the anchor, whose
extension already covers them, so dropping them is near output-neutral and pure speed. A chain whose
seeds are *all* short is left untouched, so the filter never empties a chain and never drops a read;
such all-short chains (common on low-mappability / repetitive loci) extend exactly as the default
does. Off by default (`0` → byte-identical to baseline). Extension is ~60% of `mem` CPU and almost
all of it is spent on short seeds — seeds ≤40 bp hold ~90% of all banded-SW *cells* yet are ~99%
wasted, because long seeds already resolve via the ungapped fast-path. Skipping the *redundant* ones
thins the extension stage without touching seeding or chaining.

> **Recall-safe (non-emptying filter).** Earlier releases dropped *every* short seed unconditionally,
> which emptied all-short chains and silently unmapped reads whose only evidence was short —
> negligible on clean WGS but catastrophic on low-mappability short-read data (a 151 bp
> low-mappability sample lost 63% of its mappings). The current filter only drops a short seed when a
> longer anchor survives in the same chain, making it a *strict recall improvement*: it can reduce
> extension work but can never lose a read.

**The accuracy change is small, confined to low-confidence reads, and costs no recall.** On real
HG002 1M PE WGS at `--min-ext-len 30` (non-emptying filter), the mapped count is unchanged from
default (99.75% both), and only ~0.10% of reads change locus — down from ~0.40% under the old
emptying filter, because the reads that used to vanish now map identically to default. The locus
changes concentrate in the low-confidence tail (repeat/paralog churn); ~0.005% of reads (≈100 of
2 M) change at MAPQ ≥ 60.

**In exchange, alignment CPU drops ~10% single-thread at `--min-ext-len 30`** (measured −9.5%
`main_mem` on HG002 WGS-1M), largest on data carrying many short seeds. Because the speedup thins
extension rather than speeding seeding — and seeding dominates the wall — the wall-clock gain is
smaller than the ~90% banded-SW *cell* reduction would suggest.

**Higher thresholds no longer cliff.** Because all-short chains are now protected, raising
`--min-ext-len` to `40`–`50` no longer drops reads: on HG002 WGS-1M the mapped count holds at 99.75%
and divergence stays ~0.10% across `30`–`50`. The previously-documented high-error F1 cliff (F1
−1.4% at `40`, a cliff at `50` on 2–15% simulated substitution error) was a property of the
*emptying* behavior — its mechanism, correct alignments carried *solely* on short seeds being
dropped, is exactly what the non-emptying filter now protects — so it is expected to be largely
removed, **pending re-validation** of the simulated-error sweep under the new filter. **Indels and
structural variants were never a contraindication** — an indel splits a read into two still-long
exact segments the fast-path handles, so indel-rich data is free.

`30` remains the recommended value. **Validation status:** the non-emptying filter changes the
observable output of `--min-ext-len` (and therefore `--fast`), so the cross-architecture speed
figures and the golden-truth F1 sweep warrant a fresh
[bwa-mem3-bench](../related-projects/bwa-mem3-bench.md) run (multi-thread, all regimes) to confirm
these `--fast` accuracy figures genome-wide. `--min-ext-len` and `--fast` stay opt-in; the drop-in
defaults are unchanged.

## Contained-seed extension skip: `--skip-contained-ext`

`--skip-contained-ext` skips banded Smith–Waterman extension of a seed that is **strictly contained**
(same diagonal, query subinterval) in a longer seed of the same chain. The longest seed on a diagonal
is always extended, and its alignment covers the contained seed's region, so the contained seed's own
alignment is purged after extension anyway (the post-extension containment purge). Detecting that
before extension and skipping the redundant banded-SW pair thins the extension stage without changing
the seed set (the seed stays in the chain, so seed coverage and MAPQ are untouched). Off by default.

**For non-meth short- and medium-length reads, this recommended lever has zero accuracy cost: it is
byte-identical.** The skip set is a proven subset of the post-extension purge set, so on those inputs
the output is bit-for-bit identical to the default — verified (records + header + count) on WGS 5M
PE, HiC 1M PE, and an ALT/HLA-enriched set. It no-ops under `--meth` (its own internal gate disables
it, since bisulfite's asymmetric scoring breaks the dominance argument), so the flag neither helps
nor hurts a `--meth` run.

**In exchange, alignment CPU drops ~10%** (measured −10.6% single-thread CPU on 200k HG002 WGS
2×150 single-end reads, md5-identical output; consistent with the −9.7% seen on SBX-SE in the
original characterization). Because the win is pure — no accuracy trade-off to weigh on non-meth
short- and medium-length reads — this is the safest entry in the recommended table.

> **Not byte-identical on kilobase-scale long reads.** On PacBio HiFi / ONT the "longest contained
> seed dominates" assumption breaks (extensions terminate early, so the container no longer covers
> the contained seed), and output diverges — dropped supplementary alignments, occasionally a
> degraded primary. Kilobase-scale long reads are outside the scope of this recommendation, which
> covers non-meth short- and medium-length reads only; such reads do not run at default settings
> today in any case (see
> [`--adaptive-band`](#situational---adaptive-band-for-medium-length-reads) and
> [fg-labs/bwa-mem3#238](https://github.com/fg-labs/bwa-mem3/issues/238)).

Already bundled into `--fast`. Recommended value: **on** for standard-error non-meth short- and
medium-read pipelines.

## Chain extension cap: `--max-extend-chains`

`--max-extend-chains INT` caps the number of chains that reach banded Smith–Waterman extension to the
top-`INT` by chain weight (applied after `mem_chain_flt`); the dropped chains are the lowest-weight
secondaries. It is the only lever that reduces the *number of chains extended* per read — the other
`--fast` levers (`-y 0`, `--min-ext-len`, `--smem-dedup`) cut seeds and SW-per-chain but leave chains
extended nearly unchanged — so it is orthogonal and adds a real marginal speedup on top of them. Off
by default (`0` → byte-identical to baseline). As a safety fallback the cap is a no-op for pathological
reads with more than 4096 chains (`MAX_EXTEND_CHAINS_CAP`): those reads extend all of their chains as
usual, so the option has no effect on them.

**Not byte-identical.** Dropping candidate chains removes low-weight secondaries, so `XS`, secondary
alignments, and `MAPQ` can move on multi-mapping reads. High-confidence (uniquely-placed) reads are
unaffected; the default path (`0`) is verified byte-identical to the base branch.

**The accuracy/speed tradeoff is a smooth monotonic curve with a knee at `N = 4–5`.** On holodeck
truth (`sim-wgs-place`, 10.7 M reads), standalone `N = 5` is **−23% total alignment CPU** at +20
high-confidence (MAPQ ≥ 60) mismaps out of 9.5 M (1529 → 1549). Stacked on top of `--fast` it is
**−15% marginal CPU** at +21 high-confidence mismaps (1559 → 1580, +0.0002% absolute; overall
−0.045 pp). `N = 2` is too aggressive (+13.5% high-confidence mismaps). `--fast` sets `20` and
pairs it with `--extend-mate-concordant` (below): the standalone MAPQ ≥ 1 tail rises 3.8× at cap 5,
but mate-concordant retention closes most of that at cap 20 for ~−20% aligner CPU (fg-labs/bwa-mem3#202).
`--max-extend-chains` and `--fast` stay opt-in; the drop-in defaults are unchanged.

## Mate-concordant chain retention: `--extend-mate-concordant`

The chain cap interacts badly with paired-end pairing, most acutely under `--meth`. Bisulfite reads are projected into a
3-letter alphabet (C→T, G→A), which collapses sequence complexity and *flattens chain weights*:
a read carries many similarly-weighted chains instead of one clear winner. Capping to the top-5 by
weight then frequently drops a read's *true* low-weight chain — not because it was chain-filtered
away (in 89% of the regressions instrumented on a 50k-pair slice the true chain is still a candidate
at cap-5), but because capping starves PE pairing and mate rescue of the secondary anchors that let
the true concordant pair win. Both mates then flip together to a wrong concordant locus (99% flagged
proper-pair). On a 1M-pair `sim-meth-place` slice this dropped correct placement 98.10% → 97.63% and
raised confident (MAPQ ≥ 30) mismaps 4.4× (0.036% → 0.160%).

`--extend-mate-concordant` fixes this at the chain-cap stage: when `--max-extend-chains` would cap a
paired-end read, it additionally **retains any chain concordant with one of the mate's chains** —
same contig, FR ("innie") orientation, within a window — even if it ranks below the cap. This keeps
the true pair's low-weight anchor while still dropping the far/redundant chains the cap targets. The
mate scan is bounded (`MATE_SCAN_MAX`, 256) to keep it O(n) in practice.

The window is sized to the aligner's own proper-pair insert bound. `--extend-mate-concordant` (bare)
is **auto**: it uses the estimated `pes[FR].high` (inferred from the data each chunk and read by the
next chunk's cap, which runs before pairing), falling back to a built-in default until the insert size
is known. `--extend-mate-concordant=INT` pins a fixed bp window; `=0` disables. Matching the window to
the insert bound matters because retained chains are then **extended** — a wide window admits far and
spurious concordant chains, adding alignment CPU on chain-rich reads, so auto keeps only genuine pair
anchors.

**Validation.** On the canonical bench (holodeck eval, 5 reps on m7i), the option recovers **part** of
the `--fast --meth` regression at **~1% alignment CPU** with the auto window:

| dataset (`--fast --meth`) | placement correct% | MAPQ≥30 mismap% | align CPU vs 0.5.0 |
|---|---:|---:|---:|
| default `--meth` (reference) | sim-meth-place 94.27 / sim-meth-vars 95.11 | 5.73 / 4.89 | — |
| 0.5.0 cap-5 (regressed) | 93.14 / 94.22 | 6.86 / 5.78 | baseline |
| + `--extend-mate-concordant` (auto) | 93.71 / 94.46 | 6.29 / 5.54 | **+1%** |

So placement recovers roughly a quarter to a half of the gap (not to parity), confident mismaps recover
similarly, RSS is unchanged, and the CPU cost is ~1%. The window sizing is what makes it cheap: a fixed
2000 bp window recovered the same accuracy but cost **+16–21%** CPU on chain-rich simulated reads,
because it retained and extended far/spurious concordant chains; the auto `pes[FR].high` window admits
only genuine pair anchors. (Turning the chain cap off entirely under `--meth` recovers similar accuracy
but at a larger, uniform CPU cost.) Full per-dataset figures:
[fg-labs/bwa-mem3#195](https://github.com/fg-labs/bwa-mem3/pull/195).

**`--fast` enables it automatically** for both non-meth and `--meth` runs. The
non-meth case (fg-labs/bwa-mem3#202): the top-`INT` cap inflates the confident (MAPQ ≥ 1) mis-placement
tail 3.8× on `sim-wgs-place` (3,626 → 13,921 vs uncapped) by the same mechanism — the true chain is
low-weight but mate-concordant in ~98% of cap-dropped reads. Paired with `--max-extend-chains 20`,
mate-concordant retention closes most of that tail (→ ~4,450, verified on a rebuilt `--fast`) at
~−20% aligner CPU. It is a no-op unless a chain cap is actually in effect,
and like `--max-extend-chains` it is **not** byte-identical when it retains a chain. `--extend-mate-concordant`
and `--fast` stay opt-in; the drop-in defaults are unchanged.

## Speed: drop-in and recommended

When comparing bwa-mem3 to a stock bwa-mem2 baseline, report both layers:

1. **Drop-in speed-up** — bwa-mem3 vs bwa-mem2 at *identical* settings, from its vectorized kernels,
   lockstep SMEM batching, libsais indexing, and mimalloc. See
   [Performance Overview](../performance/overview.md) and
   [What's Different — Performance](../whats-different/performance.md).
2. **Recommended-profile speed-up** — the additional ~11–22% alignment-CPU reduction from `-m 10`,
   *on top of* the drop-in gain.

A benchmark that quotes only default-vs-default settings understates the throughput available to a
caller who has adopted the recommended profile; quote both so readers can pick the comparison that
matches their deployment.

## Situational: `--supp-rep-hard-cap` for SV-aware pipelines

This is **not** part of the recommended profile — it is a situational knob, not a free speed-up, and
the default (`0`, off) is correct for plain alignment.

If you run structural-variant or breakpoint detection downstream (e.g. `fgsv SvPileup`),
`--supp-rep-hard-cap 20` suppresses repeat-induced *spurious* supplementary breakpoints — split
alignments anchored in repeats whose mapping quality is overestimated — at no measured cost to real
SV breakpoints. On GIAB HG002 CMRG (GRCh38, 2×250) it stripped repeat artifacts while preserving
every credible real-SV breakpoint, including ones in moderately-repetitive regions.

Do **not** use more aggressive values: `--supp-rep-hard-cap` ≤ 10 began suppressing *real* GIAB SV
breakpoints (a verified insertion at chr3:45890256 and deletion at chr18:79739776) whose supporting
split reads happen to land in repetitive regions.

Caveats: this was measured on a single repeat-enriched truth set at the breakpoint level (not
end-to-end SV calls), so treat `20` as a sensible starting point for SV workflows rather than a
universal recommendation. It has no effect on primary-alignment MAPQ or on non-SV pipelines.

## Situational: `--adaptive-band` for medium-length reads

This is **not** part of the recommended (short-read) profile — it is a
**medium-length-read** lever, and the default (off) is correct for standard Illumina WGS/WES.

> **Not for kilobase-scale long reads.** `bwa-mem3 mem` does not currently run on PacBio HiFi
> (~15–20 kb) or ONT reads at default settings — the short-read-sized extension buffers cannot
> hold a kilobase-scale extension and the run aborts. `--adaptive-band` does not change that.
> Use minimap2 for true long-read alignment. Tracking issue: fg-labs/bwa-mem3#238.

`--adaptive-band` starts banded Smith-Waterman tight and expands each extension only to the band its
chain's seed geometry implies, instead of the fixed `-w` band (100) for every extension. The band
only constrains the DP matrix when the extension's reference window exceeds it — which begins around
the **~200 bp** mark — so:

- **Use it for medium-length reads:** SBX (HG002, ~240 bp), where it cut alignment CPU by **~25 %**.
  It nominally covers longer reads, but kilobase-scale HiFi/ONT do not run at default settings (see
  the note above), so its usable range today is medium reads.
- **No-op on short reads:** WGS (~150 bp) and WES (~76 bp) extensions are already smaller than the
  band, so there is nothing to trim; those reads run on the 8-bit kernel, which the option leaves
  untouched. Enabling it on a short-read run neither helps nor hurts.

Accuracy is preserved: placement is identical to default on holodeck `sim-wgs-place` (MAPQ-60+ mismaps
unchanged), and indel representation matches the `-w 100` default (indels up to the chaining limit
still emit a single `D`/`I` CIGAR, so small/mid-size indel callers are unaffected). Like `--fast`, it
is **not** byte-identical when enabled (it shifts a small number of borderline secondary alignments),
which is why it is an opt-in flag rather than a default.

---

**See also:**
[Optimization checklist](optimization-checklist.md) ·
[Performance Overview](../performance/overview.md) ·
[Equivalence with bwa-mem2](../whats-different/equivalence.md) ·
[CLI Reference — `mem`](../cli/mem.md)
