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
