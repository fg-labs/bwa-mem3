# Quick start: methylation alignment

bwa-mem3 supports bisulfite-converted (WGBS/RRBS/EM-seq) read alignment through a single `--meth`
flag on both `index` and `mem`. No Python interpreter, no piped preprocessor, and no separate
postprocessing step are required.

> **Note — Drop-in replacement for bwameth.py**
>
> bwa-mem3 with `--meth` is a single-binary drop-in replacement for the `bwameth.py` pipeline.
> The output BAM is byte-compatible for the standard tags used by methylation callers (Bismark,
> MethylDackel, PileOMeth, etc.).

## Index the reference for methylation

Build the c2t doubled reference once:

```bash
bwa-mem3 index --meth ref.fa
```

This writes two additional files next to the standard index:

| File | Description |
|------|-------------|
| `ref.fa.bwameth.c2t` | C→T converted reference (forward strand) with G→A reverse complement interleaved |
| `ref.fa.bwameth.c2t.*` | FM-index files for the c2t reference |

The c2t index is separate from the standard index produced by `bwa-mem3 index ref.fa`. You need
both if you intend to run standard and methylation alignments against the same reference.

## Align bisulfite-converted reads

```bash
bwa-mem3 mem --meth -t 16 ref.fa R1.fq.gz R2.fq.gz \
  | samtools sort -o out.bam
samtools index out.bam
```

Pass the original (unconverted) reference path, not the `.bwameth.c2t` file. bwa-mem3
auto-appends `.bwameth.c2t` to the reference path when `--meth` is active.

## What `--meth` does

`--meth` activates a pipeline of in-process transformations that would otherwise require
external tools:

1. **Inline c2t read conversion.** R1 reads have every `C` converted to `T` before alignment;
   R2 reads have every `G` converted to `A`. The original unconverted sequence is preserved in
   the `YS:Z:` SAM tag. The conversion direction for each read is recorded in `YC:Z:` (value
   `CT` or `GA`), matching the bwameth.py convention.

2. **bwameth.py-equivalent scoring defaults.** `--meth` sets `-B 2 -L 10 -U 100 -T 40 -CM`
   automatically. These match the defaults used by bwameth.py and are optimized for
   bisulfite-converted reads where C→T mismatches carry no penalty. Any of these values can be
   overridden on the command line.

3. **Inline BAM post-processing.** After alignment, bwa-mem3 rewrites the SAM stream in-process:
   - `@SQ` headers with `f`/`r` prefixes (e.g. `fchr1`, `rchr1`) are collapsed back to one
     entry per real chromosome (`chr1`). Read-level `RNAME` fields are rewritten to match.
   - Each mapped record gains a `YD:Z:` tag (`f` for forward-strand, `r` for reverse-strand)
     indicating which converted strand the read aligned to.
   - Chimera QC: reads whose longest `M`/`=`/`X` run is less than 44% of the read length are
     flagged `0x200` (QC-fail), have flag `0x2` (proper pair) cleared, and have MAPQ capped at 1.
   - Pair-level QC-fail propagation: if one mate is QC-failed, the other mate is also flagged.
   - A `@PG ID:bwa-mem3-meth` program record is appended to the header.

4. **Uncompressed BAM output.** The post-processed stream is written as uncompressed BAM
   (`wb0`) rather than SAM text. This eliminates text serialization overhead and allows
   downstream `samtools sort` to read BAM natively. The stream is still fully readable by any
   htslib-based tool.

For full details on each tag, the chimera QC heuristic, and the `--set-as-failed` and
`--do-not-penalize-chimeras` flags, see the [Methylation Reference](../methylation/overview.md).

---

**See also:**
[Methylation Reference — Overview](../methylation/overview.md) ·
[Methylation Reference — SAM tags](../methylation/tags.md) ·
[Best Practices — Methylation defaults](../best-practices/methylation.md) ·
[CLI Reference — mem](../cli/mem.md)
