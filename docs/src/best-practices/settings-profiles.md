# Settings profiles: bwa drop-in vs recommended

bwa-mem3's **defaults track bwa-mem / bwa-mem2 behavior** — it is a faster drop-in for an existing
bwa-mem2 pipeline. Some defaults the bwa lineage inherited are, however, more conservative than
necessary. This page defines two profiles:

- a **drop-in profile** that keeps the bwa-mem/bwa-mem2 defaults, for migrations and parity checks;
- a **recommended profile** that deviates from those defaults where we have benchmarked the change
  to be near-neutral on accuracy and clearly faster.

The defaults ship as the drop-in profile. The recommended profile is **opt-in** — you turn it on
with explicit flags — so upgrading bwa-mem3 never silently changes your alignments.

**Which profile?**

| Your situation | Profile | Invocation |
|---|---|---|
| Migrating a bwa-mem2 pipeline, or validating against bwa/bwa-mem2 | **Drop-in** | `bwa-mem3 mem` (no extra flags) |
| New pipeline, or migration already validated | **Recommended** | `bwa-mem3 mem -m 10 -y 0` (add `-s 2` under `--meth`) |

> bwa-mem3 is already, by design, *not* byte-identical to bwa-mem2 even at default settings
> (additive SAM tags, per-architecture SIMD `score2`/MAPQ convergence, deterministic tie-breaks, a
> few extra supplementary alignments) — it is 99.94%–99.9996% concordant on primary records
> (panel-twist to WES). See
> [Equivalence with bwa-mem2](../whats-different/equivalence.md). "Drop-in" therefore means *as close
> to bwa-mem2 as bwa-mem3 gets*; the recommended profile is a further, documented deviation.

## Drop-in profile (default)

```bash
bwa-mem3 mem -t <N> ref.fa R1.fq R2.fq > out.sam
```

No extra flags. Use this when you are:

- reproducing or validating against bwa-mem / bwa-mem2 output, or
- migrating an existing bwa-mem2 pipeline and want to change one variable (the aligner binary) at a
  time before tuning anything else.

## Recommended profile

```bash
bwa-mem3 mem -t <N> -m 10 -y 0 ref.fa R1.fq R2.fq > out.sam
```

> **Shorthand:** `bwa-mem3 mem --fast` applies `-m 10 -y 0 --min-ext-len 30
> --smem-dedup --skip-contained-ext` (and `-s 2` under `--meth`) in one flag.
> Explicit flags still override individual levers where applicable; `--smem-dedup`
> and `--skip-contained-ext` are forced on with no opt-out. `--skip-contained-ext`
> no-ops under `--meth` (its own internal gate disables it there), so on a `--meth`
> run the effective levers are `-m 10 -y 0 --min-ext-len 30 --smem-dedup -s 2`.
> See [`mem` → `--fast`](../cli/mem.md#--fast--speed-preset-opt-in-not-byte-identical).

Use this for new pipelines, or once a drop-in migration is validated and you want bwa-mem3's best
speed/accuracy trade-off. Current recommended deviations:

| flag | default (drop-in) | recommended | effect |
|---|---|---|---|
| `-m` (mate-rescue depth) | 50 | **10** | ~11–22% less alignment CPU; near-neutral accuracy (see below) |
| `-y` (3rd-round seeding occurrence) | 20 | **0** | ~11–30% less alignment CPU; F1 near-neutral across regimes (within ±0.02; better on divergent/repeat — see below) |
| `-s` (Pass-2 re-seed width), **`--meth` only** | 10 | **2** | light re-seed: ~same speed as `-s 0` but recovers the MAPQ/placement `-s 0` lost (see below) |
| `--min-ext-len` (skip short-seed extension), **standard-error reads** | 0 | **30** | ~10–20% less alignment CPU; accuracy change confined to the already-low-confidence tail (see below) |

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

---

**See also:**
[Optimization checklist](optimization-checklist.md) ·
[Performance Overview](../performance/overview.md) ·
[Equivalence with bwa-mem2](../whats-different/equivalence.md) ·
[CLI Reference — `mem`](../cli/mem.md)
