# Methylation Defaults

`bwa-mem3 mem --meth` ships with scoring and filtering defaults aligned with the
bwameth.py reference implementation (in the default `collapsed` scoring mode).
This page is prescriptive: when to keep those defaults and when to override them.
For *what* the flags are and *how* the scoring modes work, see the reference:

- **Scoring modes (`collapsed` vs `genomic`) and the full `--meth` flag set** —
  [Methylation Reference → Flags](../methylation/flags.md#--meth-scoring-collapsedgenomic).
- **Placement compatibility with bwameth.py** (the default `collapsed` mode is a
  placement drop-in, *not* byte-identical — ~1% of records differ in
  `POS`/`CIGAR`/`MAPQ`, so **re-validate if you are pinned to a specific bwameth
  release**) — [bwameth.py drop-in mapping](../methylation/bwameth-mapping.md).

## When to keep the defaults

For standard whole-genome bisulfite sequencing (WGBS) workflows, the defaults
(`collapsed` scoring plus `-B 2 -L 10 -U 100 -T 40 -M -C`) are appropriate as-is.
They were derived from the bwameth.py codebase and are expected by most
downstream methylation callers. Unless you have a specific reason to deviate, use:

```bash
bwa-mem3 mem --meth -t 16 ref.fa R1.fq.gz R2.fq.gz \
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
bwa-mem3 mem --meth --set-as-failed r -t 16 ref.fa R1.fq.gz R2.fq.gz \
  | samtools sort -@ 4 -o out.bam -
```

**Chimera QC is opt-in** (off by default, matching Bismark). Leave it off for
standard directional EM-seq / WGBS. Turn it on for PBAT / scBS-Seq libraries
(where intra-fragment chimerism is common) or when you want bwameth.py-equivalent
flagging — pass `--chimera-qc` (see [Flags](../methylation/flags.md#--chimera-qc)
for the exact heuristic):

```bash
bwa-mem3 mem --meth --chimera-qc -t 16 ref.fa R1.fq.gz R2.fq.gz \
  | samtools sort -@ 4 -o out.bam -
```

> **Note — Overrides are positional**
>
> Flags supplied after `--meth` on the command line override the defaults set by
> `--meth`. For example, `bwa-mem3 mem --meth -B 4 ...` uses `-B 4` (not 2).
> Flags supplied before `--meth` are silently overwritten by `--meth`'s defaults,
> so always place overrides after `--meth`.

## Downstream tool compatibility

The `--meth` output BAM carries the Bismark `XR:Z` / `XG:Z` / `XM:Z` tag set, so
it feeds Bismark-aware callers directly — `bismark_methylation_extractor`,
methylKit, methtuple, DMRfinder, epialleleR, MethylDackel, and biscuit per-read
tools have all been used successfully. See [SAM tags](../methylation/tags.md) for
the tag definitions and which tools read which tag.

---

**See also:**
[Methylation Reference: Overview](../methylation/overview.md) ·
[SAM tags: XR, XG, XM](../methylation/tags.md) ·
[Flags: --set-as-failed, --chimera-qc](../methylation/flags.md) ·
[Quick start: methylation alignment](../getting-started/quick-meth.md) ·
[Output format](output-format.md)
