# Methylation Defaults

`bwa-mem3 mem --meth` ships with a set of scoring and filtering defaults that
match the bwameth.py reference implementation. This page describes what those
defaults are, when to keep them, and when to override them.

## What `--meth` sets

When `--meth` is passed, the following flags are applied automatically in
addition to enabling inline c2t conversion and BAM post-processing:

| Flag | Value | Purpose |
|------|-------|---------|
| `-B` | `2` | Mismatch penalty. Reduced from the bwa-mem2 default of 4. Bisulfite-treated reads carry C→T and G→A mismatches at converted positions; a lower penalty prevents these from causing spurious soft-clipping or unmapped reads. |
| `-L` | `10` | Clipping penalty. Increased from the bwa-mem2 default of 5 to discourage clipping of read ends that carry converted bases at positions that look like mismatches. |
| `-U` | `100` | Unpaired read penalty. Higher than default; methylation libraries typically have well-defined insert sizes and anomalous pairing usually reflects a mapping artifact. |
| `-T` | `40` | Minimum alignment score threshold. Higher than default; raises the bar to report an alignment, reducing spurious low-quality hits against the doubled reference. |
| `-CM` | — | Treats soft-clipped bases as matches in CIGAR output. Required for correct behavior of downstream methylation callers (e.g. Bismark, MethylDackel) that count clipped bases. |

These defaults can all be overridden on the command line. The `--meth` flag
sets them first; any explicit flag that follows overrides the `--meth`-set
value.

## When to keep the defaults

For standard whole-genome bisulfite sequencing (WGBS) workflows, the defaults
are appropriate as-is. They were derived from the bwameth.py codebase and are
expected by most downstream methylation calling tools. Unless you have a
specific reason to deviate, use:

```bash
bwa-mem3 mem --meth --bam=0 -t 16 ref.fa R1.fq.gz R2.fq.gz \
  | samtools sort -@ 4 -o out.bam -
samtools index out.bam
```

## When to override

**Low-coverage or targeted bisulfite sequencing.** If your library covers a
small target region and insert sizes are more variable, consider lowering `-T`
(e.g. `-T 20`) to recover short or soft-clipped alignments in the target.

**Amplicon bisulfite sequencing.** Amplicon reads have uniform insert sizes;
the default `-U 100` is appropriate. However, if your amplicons are short
(< 100 bp), consider lowering `-L` further to reduce clipping at read ends.

**Non-standard conversion chemistry.** Some library preparations use only one
strand conversion (C→T only, not G→A). In such cases, `--set-as-failed r`
suppresses alignments to the reverse-complement strand, which reduces noise
from strand-ambiguous alignments:

```bash
bwa-mem3 mem --meth --set-as-failed r --bam=0 -t 16 ref.fa R1.fq.gz R2.fq.gz \
  | samtools sort -@ 4 -o out.bam -
```

**Chimeric reads from long-read-length protocols.** By default, `--meth`
applies a chimera QC heuristic: if the longest matching run (CIGAR `M`/`=`/`X`
operations) is less than 44% of the read length, the alignment is flagged
`0x200` (QC fail), the paired flag `0x2` is cleared, and MAPQ is capped at 1.
If your protocol produces legitimate long reads where this heuristic
over-aggressively flags alignments, pass `--do-not-penalize-chimeras`:

```bash
bwa-mem3 mem --meth --do-not-penalize-chimeras --bam=0 -t 16 ref.fa R1.fq.gz R2.fq.gz \
  | samtools sort -@ 4 -o out.bam -
```

> **Note — Overrides are positional**
>
> Flags supplied after `--meth` on the command line override the defaults set by
> `--meth`. For example, `bwa-mem3 mem --meth -B 4 ...` uses `-B 4` (not 2).
> Flags supplied before `--meth` are silently overwritten by `--meth`'s defaults,
> so always place overrides after `--meth`.

## Downstream tool compatibility

The `--meth` output BAM is designed to be a drop-in replacement for the output
of the `bwameth.py` pipeline. The following downstream tools have been used
successfully with `bwa-mem3 --meth` output:

- **MethylDackel** — extracts methylation calls from the `YD:Z:` strand tag.
- **Bismark** — accepts the bwameth-convention `YD:Z:` strand annotation.
- **PileOMeth** — reads the standard bisulfite BAM format.

If a tool requires the `XB:Z:` tag convention used by Bismark's own aligner
rather than the `YD:Z:` convention, a conversion step is needed before
methylation calling.

---

**See also:**
[Methylation Reference: Overview](../methylation/overview.md) ·
[SAM tags: YS, YC, YD](../methylation/tags.md) ·
[Flags: --set-as-failed, --do-not-penalize-chimeras](../methylation/flags.md) ·
[Quick start: methylation alignment](../getting-started/quick-meth.md) ·
[Output format](output-format.md)
