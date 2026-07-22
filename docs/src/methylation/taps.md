# TAPS libraries

TAPS (TET-assisted pyridine borane sequencing) is supported via
`--meth=taps`. This page is the short version; the mechanics live on the
existing pages and are linked from here rather than repeated.

## TL;DR

```bash
# Index once — the SAME index serves both chemistries.
bwa-mem3 index --meth ref.fa

# Align TAPS reads. Note the '=' — see the warning below.
bwa-mem3 mem --meth=taps -t 16 ref.fa R1.fq.gz R2.fq.gz | samtools sort -o out.bam
```

> **Do not align TAPS data with a bare `--meth`.** It will succeed, produce
> well-placed reads, and emit `XM:Z` methylation calls that are **completely
> inverted** — percent methylation comes out as `1 − truth`, with no warning.
> Measured against simulated ground truth: Pearson **r = −0.956** with `--meth`
> versus **+0.956** with `--meth=taps`.
>
> **Use `--meth=taps`, not `--meth taps`.** The argument is optional, and
> `getopt_long` only accepts optional arguments in the `=` form. A separate word
> is parsed as a positional — i.e. as your first FASTQ.

## Why one flag is enough

Bisulfite/EM-seq converts the **un**methylated cytosines; TAPS converts the
**methylated** ones (5mC/5hmC). But both produce the *same* base change as far
as the aligner can see — ref `C` → read `T` on the top strand, ref `G` → read
`A` on the bottom. So seeding, the converted index, and the scoring matrices
are shared, and TAPS reads align correctly with the existing machinery.

The chemistries differ only in *which* cytosines converted, which inverts the
meaning of the observed base when **calling** methylation. `--meth=taps`
therefore changes the `XM:Z` call polarity and nothing else. CpG/CHG/CHH
context classification is a property of the reference, so it is identical
either way. See [SAM tags](tags.md#xmz--methylation-call-string).

## What `--meth=taps` changes

| | `--meth` (em-seq) | `--meth=taps` |
|---|---|---|
| `XM:Z` polarity | retained `C` = methylated | retained `C` = **un**methylated |
| `--meth-scoring` default | `collapsed` | **`neutral`** |
| Index, seeding, `XR`/`XG` | identical | identical |

The scoring default differs because TAPS conversions are far sparser — 3.2 % of
cytosines versus 94.6 % for em-seq on a matched simulation — so the collapsed
3-letter alphabet costs specificity it can no longer repay. `neutral` goes one
step further than `genomic`: it scores the conversion cell as `0` (tolerated but
not rewarded) rather than a full match, so a sparse TAPS conversion no longer
over-credits a spurious C→T alignment. Measured placement on simulated TAPS
against full hg38, across three methylation loads:

| arm | placement |
|---|---|
| `--meth=taps` (i.e. `neutral`) | **95.95–96.01 %** |
| `--meth-scoring genomic` | 95.68–95.73 % |
| plain (no `--meth`) | ~95.53 % |
| `--meth-scoring collapsed` | 95.37–95.43 % |
| *unconverted baseline* | *96.02 %* |

`neutral` reaches essentially the unconverted ceiling — a robust +0.24–0.28 pp
over `genomic` — with identical MAPQ calibration (60+ bin 99.99 % either way) and
`NM` (mean 1.966 vs 1.969), so the gain is not bought with over-confidence or
inflated edit distance. An explicit `--meth-scoring` always overrides the default.

## Downstream methylation calling

**Most extractors assume bisulfite polarity.** MethylDackel reports the
fraction of reads showing `C` at a reference C, which under TAPS is the
**un**methylated fraction — so TAPS methylation is `100 − MethylDackel_rate`.
Either apply that inversion explicitly, or use a TAPS-aware caller such as
[Rastair](https://bitbucket.org/bsblabludwig/rastair) (what nf-core/methylseq
uses for its TAPS path).

`XM:Z` from `--meth=taps` is already in the correct polarity, so tools that
consume `XM:Z` directly need no adjustment.

## Library directionality

`--meth` assumes a **directional** library (R1 → `CT`, R2 → `GA`). Standard TAPS
protocols are directional; verified on real data (`SRR25526768`), where the
`XG:Z` split is 49.8 % `CT` / 49.9 % `GA`. Undirectional libraries are not
supported.

## Validation status

- Simulated TAPS with ground truth: `r = +0.956`, RMSE 0.117 (em-seq control:
  `+0.954`, 0.119) over ~600 k CpG sites.
- Real human TAPS (`SRR25526768`, HAP1): `--meth=taps` and `--meth` produce
  **exact complement** call sets on identical alignments.
- **Not yet validated:** absolute call accuracy on real data (no truth set),
  x86 SIMD tiers (results are arm64/NEON), a methylation-rate sweep, and
  undirectional libraries.

Full experiment, including the controls and what was not tested:
`reports/2026-07-20-taps-alignment-experiment-results.md`.

---

**See also:**
[Flags](flags.md#--methemseqtaps) ·
[SAM tags: XR, XG, XM](tags.md) ·
[Conversion details](conversion.md) ·
[Overview](overview.md)
