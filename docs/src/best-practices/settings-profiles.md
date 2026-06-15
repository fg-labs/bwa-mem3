# Settings profiles: bwa drop-in vs recommended

bwa-mem3's **defaults are chosen to track bwa-mem / bwa-mem2 behavior** — it is meant to be a
faster drop-in for an existing bwa-mem2 pipeline. Some defaults that the bwa lineage inherited are,
however, more conservative than necessary. This page defines two profiles:

- a **drop-in profile** that keeps the bwa-mem/bwa-mem2 defaults, for migrations and parity checks;
- a **recommended profile** that deviates from those defaults where we have benchmarked the change
  to be near-neutral on accuracy and clearly faster.

The defaults ship as the drop-in profile. The recommended profile is **opt-in** — you turn it on
with explicit flags — so upgrading bwa-mem3 never silently changes your alignments.

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
bwa-mem3 mem -t <N> -m 10 ref.fa R1.fq R2.fq > out.sam
```

Use this for new pipelines, or once a drop-in migration is validated and you want bwa-mem3's best
speed/accuracy trade-off. Current recommended deviations:

| flag | default (drop-in) | recommended | effect |
|---|---|---|---|
| `-m` (mate-rescue depth) | 50 | **10** | ~11–22% less alignment CPU; near-neutral accuracy (see below) |

This table will grow as we benchmark additional tunings; each entry is gated on the same
"measurably faster, near-neutral accuracy" bar.

## Why we recommend `-m 10`

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

Numbers are from [bwa-mem3-bench](../related-projects/bwa-mem3-bench.md) on 5 M-read real datasets
plus a multi-contig holodeck golden-truth ablation; consult the bench for methodology and current
figures.

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
[Tuning checklist](../performance/tuning.md) ·
[Performance Overview](../performance/overview.md) ·
[Equivalence with bwa-mem2](../whats-different/equivalence.md) ·
[CLI Reference — `mem`](../cli/mem.md)
