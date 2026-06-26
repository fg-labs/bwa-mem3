# Quick start: methylation alignment

bwa-mem3 supports bisulfite-converted (WGBS/RRBS/EM-seq) read alignment through a single `--meth`
flag on both `index` and `mem`. No Python interpreter, no piped preprocessor, and no separate
postprocessing step are required.

> **Note — bwameth-compatible by default, variant-aware on request**
>
> By default (`--meth-scoring collapsed`) bwa-mem3 closely tracks bwameth.py's read **placement** and
> emits the standard Bismark tags methylation callers expect (MethylDackel, Bismark, PileOMeth,
> etc.). It is a placement drop-in — *not* a byte-for-byte reproduction of bwameth: a small fraction
> of records (~1% on typical WGBS/EM-seq) still differ in `POS`/`CIGAR`/`MAPQ`, so re-validate if you
> are pinned to a specific bwameth release (see
> [bwameth.py drop-in mapping](../methylation/bwameth-mapping.md)). Add
> `--meth-scoring genomic` to opt into variant-aware scoring (truthful `NM`/`MD`; one BAM for both
> methylation and variant calling).

## Index the reference for methylation

Build the methylation index once:

```bash
bwa-mem3 index --meth ref.fa
```

This builds the normal 4-letter index at the bare prefix **plus** a converted seed index:

| File(s) | Description |
|---------|-------------|
| `ref.fa.amb`, `ref.fa.ann`, `ref.fa.bwt.2bit.64`, `ref.fa.pac` | Normal index over the original reference (used for scoring/extension; bases are pac-fetched from `.pac`, so no `.0123` is built). |
| `ref.fa.meth.fa` | Per-strand C→T / G→A converted FASTA (`f`/`r` doubled contigs). |
| `ref.fa.meth.*` | FM-index over the converted FASTA (used only for **seeding**). |

Both indexes live next to `ref.fa`; `bwa-mem3 index ref.fa` (no `--meth`) builds only the first set.

## Align bisulfite-converted reads

```bash
# Default: collapsed (bwameth-compatible placement)
bwa-mem3 mem --meth -t 16 ref.fa R1.fq.gz R2.fq.gz \
  | samtools sort -o out.bam
samtools index out.bam

# Opt into variant-aware scoring
bwa-mem3 mem --meth --meth-scoring genomic -t 16 ref.fa R1.fq.gz R2.fq.gz \
  | samtools sort -o out.bam

# Single-end (RRBS etc.): pass one FASTQ. The C→T (R1) projection is used for all reads.
bwa-mem3 mem --meth -t 16 ref.fa R1.fq.gz \
  | samtools sort -o out.bam
```

Pass the original (unconverted) reference path. bwa-mem3 finds the `ref.fa.meth.*` seed index
automatically when `--meth` is active. The paired-end-only flags (`-U 100`, `-p`) do not apply
to single-end input.

## What `--meth` does

`--meth` activates a pipeline of in-process steps that would otherwise require external tools:

1. **Seed in 3-letter space, score in 4-letter space.** Each read is projected (R1 `C→T`, R2 `G→A`)
   to find seeds in the `ref.fa.meth.*` index, then extended and **scored against the original
   4-letter reference** with a per-strand asymmetric matrix. The original bases are restored into
   the BAM `SEQ` field on emit, and the conversion direction is reported per record in the Bismark
   `XR:Z` tag (`CT` for R1/SE, `GA` for R2).

2. **Scoring defaults.** `--meth` sets `-L 10 -U 100 -T 40 -M -C` in both modes, plus the
   mode-dependent mismatch penalty: `-B 2` for `collapsed`, `-B 4` for `genomic`. These mirror
   bwameth's `bwa mem -T 40 -B 2 -L 10 -CM` (with `-U 100` for paired-end). Any value can be
   overridden on the command line after `--meth`.

3. **Inline BAM post-processing.** After alignment, bwa-mem3 rewrites the stream in-process:
   - `@SQ` headers with `f`/`r` prefixes (e.g. `fchr1`, `rchr1`) are consolidated back to one entry
     per real chromosome (`chr1`); record `RNAME`/`RNEXT` fields are rewritten to match.
   - Each mapped record gains Bismark `XG:Z` (genome strand) and `XM:Z` (per-base methylation call
     string).
   - Optional chimera QC (`--chimera-qc`, **off by default** to match Bismark): reads whose longest
     `M`/`=`/`X` run is under 44% of the read are flagged `0x200`, have `0x2` cleared, and MAPQ
     capped at 1. QC-fail flags then propagate across the read group.
   - A `@PG ID:bwa-mem3-meth` program record is appended to the header.

4. **Uncompressed BAM output.** The stream is written as uncompressed BAM (`wb0`) rather than SAM
   text, so downstream `samtools sort` reads it natively. It is still readable by any htslib tool.

For the scoring modes, each tag, the optional chimera QC heuristic, and the `--set-as-failed` /
`--chimera-qc` flags, see the [Methylation Reference](../methylation/overview.md).

---

**See also:**
[Methylation Reference — Overview](../methylation/overview.md) ·
[Methylation Reference — Flags](../methylation/flags.md) ·
[Methylation Reference — SAM tags](../methylation/tags.md) ·
[Best Practices — Methylation defaults](../best-practices/methylation.md) ·
[CLI Reference — mem](../cli/mem.md)
