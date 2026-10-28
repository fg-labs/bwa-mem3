# Illumina 5-base libraries

Illumina 5-base libraries are supported via `--meth=5base` (alias `5-base`;
[#650](https://github.com/fg-labs/bwa-mem3/pull/650)). There is deliberately no
`illumina` alias: Illumina also sells bisulfite methylation kits, whose calls
`--meth=5base` would invert. The aligner treats 5-base exactly as
[TAPS](taps.md): the option parser maps `5base` to the TAPS chemistry, so
alignment, `XM:Z` polarity and the default scoring are those of `--meth=taps`.
Every alignment record is therefore identical to `--meth=taps` by construction,
a property of the code rather than of any one host or SIMD tier.
`test/regression/meth_5base_alias.sh` checks it record for record in CI,
and the whole SAM output of the simulation
[below](#measured-on-simulated-5-base-reads) (Apple M-series, arm64 NEON)
matched: its 400,000 primary records and every secondary and supplementary one.
Only the `@PG` header's `CL:` differs, since it records the spelling used; what
a `Z` call measures differs too (see 5hmC below), but not how reads align.

## TL;DR

```bash
# Index once — the same index serves every chemistry.
bwa-mem3 index --meth ref.fa

# Align 5-base reads. Note the '=' (see TAPS libraries for why).
bwa-mem3 mem --meth=5base -t 16 ref.fa R1.fq.gz R2.fq.gz | samtools sort -o out.bam
```

> **Do not align 5-base data with a bare `--meth`.** As with TAPS, the reads
> place well but every `XM:Z` call is inverted, with no warning.

## Why 5-base is TAPS to the aligner

The 5-base chemistry converts **methylated** cytosine (5mC) to T and leaves
unmethylated C alone, the reverse of bisulfite and EM-seq. Adapters are
ligated before the DNA is denatured and converted, so the library is
directional: read 1 is the converted original strand, read 2 its reverse
complement. That gives the same aligner-visible change as the other
chemistries (ref `C` → read `T` on the top strand, ref `G` → read `A` on the
bottom) and the TAPS meaning of the observed base:

| Read base at a reference C | `--meth` (em-seq) | `--meth=taps` / `--meth=5base` |
|---|---|---|
| `C` (retained) | methylated (`Z`) | **un**methylated (`z`) |
| `T` (converted) | **un**methylated (`z`) | methylated (`Z`) |

Conversions are also sparse, since only methylated C converts, so `neutral`
scoring (the TAPS default) applies for the same reason it does for TAPS. See
[What `--meth=taps` changes](taps.md#what---methtaps-changes).

The two chemistries differ in 5-hydroxymethylcytosine: TAPS converts 5hmC along
with 5mC, while the 5-base enzyme has minimal activity on 5hmC. That changes
what a `Z` call measures (5mC only, versus 5mC + 5hmC), not how reads align.

## Measured on simulated 5-base reads

200,000 directional 2×150 bp pairs simulated from a 5 Mb chr22 slice of hg38:
methylation only at CpG (strand-symmetric per-site levels), only methylated C
converted, 0.2% sequencing error. Measured on Apple M-series (arm64, NEON
tier), clang 21.

| configuration | placement | CpG `XM:Z` vs truth (Pearson r) |
|---|---|---|
| `--meth=5base` (`neutral`) | **98.27 %** | **+0.963** |
| `--meth=5base --meth-scoring genomic` | 98.08 % | +0.963 |
| `--meth=5base --meth-scoring collapsed` | 97.82 % | +0.961 |
| bare `--meth` (em-seq polarity) | 97.82 % | **−0.961** |
| no `--meth` | 98.15 % | — |

Methylated CpG calls were unbiased on all four strand × read combinations
(OT/OB × R1/R2), and 0.07 % of non-CpG calls read as methylated, which is the
simulated sequencing-error floor.

## Choosing `--meth-context cpg`

5-base converts only methylated cytosines, and in mammalian somatic tissue
those are almost all at CpG. Under the default `--meth-context all` the aligner
still tolerates a conversion at every `C`, so every candidate locus that
differs from the read only by C/T scores close to the true one and `MAPQ`
drops. `--meth-context cpg` tolerates a conversion only at a CpG
([model](conversion.md#cpg-context---meth-context-cpg)); a C→T at a non-CpG `C`
(off-target conversion, damage, a variant) costs a mismatch.

```bash
bwa-mem3 mem --meth=5base --meth-context cpg -t 16 ref.fa R1.fq.gz R2.fq.gz
```

Measured on simulated directional 5-base reads from hg38 chr22, with
chemistry noise modelled on real 5-base libraries (1.2 % off-target conversion
of unmethylated C in every context, 3 % failed conversion of methylated CpG).
Correct placement / wrong placement, primary records:

| set | `MAPQ` | `all` | `cpg` | CpG-only upper bound |
|---|---|---|---|---|
| 2×150 bp, 338,790 reads | ≥ 20 | 97.43 / 0.05 % | **97.59 / 0.06 %** | 97.84 / 0.06 % |
| 2×150 bp | any | 98.73 / 1.27 % | 98.74 / 1.26 % | 98.81 / 1.19 % |
| 20–40 bp SE, `-k 14 -T 18`, 228,683 reads | ≥ 20 | 74.26 / 0.13 % | **76.24 / 0.18 %** | 75.52 / 0.08 % |
| 20–40 bp SE | any | 84.74 / 13.54 % | 84.21 / 13.12 % | 85.80 / 12.46 % |

The upper bound is the same reads with every true conversion undone and the
off-target noise kept, aligned without `--meth`. At `MAPQ` ≥ 20 `cpg` gains
at every conversion load and on reads carrying a real SNV, in both sets. On
2×150 bp reads the gain is confident placements (any-`MAPQ` placement is
unchanged), and it closes about 40 % of the gap to the bound. On 20–40 bp
reads it adds 2 points at `MAPQ` ≥ 20, with a slightly higher wrong rate, but
places 0.5 points fewer reads overall: when a 20–29 bp read
carries an off-target conversion inside its seed, that conversion now costs a
mismatch, which can drop the read below `-T`. `XM:Z` is the same in both modes
(Z-count equal to the read's true CpG conversions for 99.97 % of correctly
placed 2×150 reads).

`cpg` costs run time: on the 2×150 bp set, +5–7 % wall on x86 (AVX-512BW,
clang 19.1.7, 16 threads) and +7–9 % on arm64 (Apple M2 Max, NEON, clang 21,
8 threads), split between the pairing stage (mate rescue, CIGAR and score
regeneration) and banded extension.

**The default stays `all`.** These measurements cover chr22 only. A later
release may make `cpg` the 5-base default once the same comparison holds on the
full genome, overall and on the hypermethylated and SNV-carrying reads.
Until then, choose `cpg` for mammalian 5-base libraries when confident
placement matters, and keep `all` where non-CpG methylation is real biology
(see [TAPS](taps.md#--meth-context-stays-all-for-taps)) or where recall of very
short fragments at any `MAPQ` matters most.

## Before aligning

- **Trim the UMI cycles.** 5-base reads start with a 7-base UMI and one skipped
  base. Remove them at demultiplexing (the demultiplexer's cycle-override
  setting, `OverrideCycles` `U7N1Y<n>;I10;I10;U7N1Y<n>`, as in Illumina's
  [5-base pipeline guide](https://help.dragen.illumina.com/dragen-v4.5/product-guides/dragen-v4.5/dragen-5-base/dragen-5base-pipeline)).
  Left in, those eight bases are not genomic: they clip or mismatch at the read
  start, and any reference C they cover can be called.
- **Controls.** The kits spike in unmethylated lambda and CpG-methylated pUC19.
  To measure conversion from them, append both sequences to the reference
  before `index --meth`.

## Downstream methylation calling

Extractors that assume bisulfite polarity report the **un**methylated
fraction on 5-base data, exactly as on TAPS: use `100 − rate`, or a TAPS-aware
caller. `XM:Z` from `--meth=5base` is already in
the correct polarity. See [TAPS: downstream methylation calling](taps.md#downstream-methylation-calling).

A C→T SNV at a reference C looks like a conversion within a single read. The
5-base analysis model separates them using the opposite strand (a methylated C
pairs with G, a C→T variant with A), which is a caller's job rather than the
aligner's.
