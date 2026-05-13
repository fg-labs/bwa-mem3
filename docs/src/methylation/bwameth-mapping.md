# bwameth.py Drop-In Mapping

`bwa-mem3 --meth` is designed to produce output that is equivalent to the
bwameth.py pipeline for the standard paired-end case. This page explains what
changes between the two approaches and what stays the same.

## Command comparison

### bwameth.py pipeline (multi-step)

```sh
# Step 1: build a doubled reference with bwameth.py
bwameth.py index ref.fa                # writes ref.fa.bwameth.c2t + bwa-mem2 FMI

# Step 2: align (bwameth.py converts reads, calls bwa-mem2, post-processes)
bwameth.py map --bwa-mem2 -t 16 ref.fa R1.fq.gz R2.fq.gz \
  | samtools sort -o out.bam
samtools index out.bam
```

### bwa-mem3 --meth (single binary)

```sh
# Step 1: build the doubled reference with bwa-mem3
bwa-mem3 index --meth ref.fa           # same ref.fa.bwameth.c2t layout as bwameth.py

# Step 2: align (inline c2t conversion + post-processing, no Python)
bwa-mem3 mem --meth -t 16 ref.fa R1.fq.gz R2.fq.gz \
  | samtools sort -o out.bam
samtools index out.bam
```

The index files produced by `bwa-mem3 index --meth` and `bwameth.py index` are
identical in layout: the same `ref.fa.bwameth.c2t` doubled-reference FASTA
followed by the bwa-mem2 FM-index files (`.bwt.2bit.64`, `.0123`, `.pac`,
`.amb`, `.ann`).

## What is gained

**No Python or bwameth.py dependency.** The entire pipeline — read conversion,
alignment, and BAM post-processing — runs inside a single `bwa-mem3` process.
This simplifies deployment: one binary, no virtual environment, no version
pinning of bwameth.py.

**No intermediate files.** bwameth.py writes a converted FASTQ (or pipes it)
before handing off to the aligner. `bwa-mem3 --meth` performs the C→T / G→A
conversion in-memory on each read batch before passing it to the alignment
kernel. No temporary FASTQ is written and no extra pipe stage is needed.

**Inline BAM post-processing.** Header rewriting, Bismark `XR:Z` / `XG:Z` /
`XM:Z` tag emission, opt-in chimera QC (`--chimera-qc`), and QC-fail
propagation all happen inside the same process and the same pass over the
alignments. There is no separate post-processing step. Output is
written as uncompressed BAM (`wb0`) — a near-zero-cost format that downstream
`samtools sort` reads natively.

**Same flag defaults.** `--meth` applies `-B 2 -L 10 -U 100 -T 40 -CM`
automatically, matching bwameth.py's default scoring. All parameters can be
overridden.

## What stays the same

The output BAM is field-compatible with bwameth.py output for the standard
methylation tag set, flags, and SEQ representation (the `@PG` provenance line
intentionally differs — see below):

| Field | bwameth.py | bwa-mem3 --meth |
|-------|-----------|-----------------|
| `@SQ` headers | One per real chromosome | One per real chromosome |
| Methylation aux tags | `YS:Z`, `YC:Z`, `YD:Z` (bwameth) | `XR:Z`, `XG:Z`, `XM:Z` (Bismark-compatible) |
| `@PG` | `ID:bwameth` | `ID:bwa-mem3-meth` |
| Chimera QC threshold | Longest M < 44% of read | Same (44%), opt-in via `--chimera-qc` |
| Chimera QC flags | `0x200`, clear `0x2`, MAPQ ≤ 1 | Same |
| SEQ field | Pre-c2t bases (RC-flipped when `is_rev`) | Same |

The `@PG ID:` is intentionally different so provenance is unambiguous.
bwa-mem3 `--meth` emits the **Bismark-compatible** `XR:Z` / `XG:Z` /
`XM:Z` tag set rather than the bwameth-style `YS:Z` / `YC:Z` / `YD:Z`
set, which means the output is directly consumable by
`bismark_methylation_extractor`, methylKit, methtuple, DMRfinder, and
epialleleR in addition to MethylDackel and biscuit. Downstream tools
that read `YS:Z` / `YC:Z` / `YD:Z` will not find those tags and must be
pointed at the corresponding `XR:Z` / `XG:Z` (and the per-base `XM:Z`
methylation call string) instead.

> **Info — End-to-end regression coverage**
>
> PR #13 includes a three-layer regression test that verifies 100% chrom+pos
> match, 100% CIGAR match, and byte-identical SEQ across 92,684 paired-end
> records compared to a bwameth.py reference run.

## When to prefer bwameth.py

If your workflow requires bwameth.py-specific features (e.g. `bwameth.py
markduplicates` or non-standard bwameth.py post-processors), continue using
bwameth.py. `bwa-mem3 --meth` targets the indexing + alignment + standard
post-processing path only.

---

**See also:**
[Overview](overview.md) ·
[Conversion details](conversion.md) ·
[SAM tags: XR, XG, XM](tags.md) ·
[Chimera QC and header rewriting](post-processing.md) ·
[Related Projects: bwameth.py](../related-projects/bwameth.md)
