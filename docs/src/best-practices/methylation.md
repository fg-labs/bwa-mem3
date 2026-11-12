# Methylation Defaults

`bwa-mem3 mem --meth` ships with filtering defaults aligned with the bwameth.py
reference implementation (`-L 10 -U 100 -T 40 -M -C`) and, for EM-seq/bisulfite,
the variant-aware `genomic` scoring mode.
This page is prescriptive: when to keep those defaults and when to override them.
For *what* the flags are and *how* the scoring modes work, see the reference:

- **Scoring modes (`collapsed`, `genomic`, `neutral`) and the full `--meth` flag
  set** —
  [Methylation Reference → Flags](../methylation/flags.md#--meth-scoring-collapsedgenomicneutral).
- **Placement compatibility with bwameth.py** (`--meth-scoring collapsed` is a
  placement drop-in, *not* byte-identical — ~1% of records differ in
  `POS`/`CIGAR`/`MAPQ`, so **re-validate if you are pinned to a specific bwameth
  release**) — [bwameth.py drop-in mapping](../methylation/bwameth-mapping.md).

## When to keep the defaults

First verify the [directional library requirement](../methylation/overview.md#library-orientation):
R1/single-end uses C-to-T and R2 uses G-to-A as input conversion hypotheses,
not fixed genomic mapping strands. Neither chemistry selection nor QC flags
infer or change library orientation.

For standard directional whole-genome bisulfite and EM-seq workflows, the defaults
(`genomic` scoring plus `-B 4 -L 10 -U 100 -T 40 -M -C`) are appropriate as-is.
The filtering values come from the bwameth.py codebase. `genomic` keeps the
mirror change (reference T read as C) a mismatch, which gives better-calibrated
MAPQ than `collapsed`: on simulated EM-seq (chr20–22 at 2×, 2×150 bp, full hg38,
Apple M2 Max, arm64 NEON) it puts 1.3 pp more reads at MAPQ ≥ 20 in the right place and makes
about 20× fewer wrong placements at MAPQ 60
([measurements](../methylation/overview.md#three-scoring-modes---meth-scoring),
[#756](https://github.com/fg-labs/bwa-mem3/pull/756)).
Unless you have a specific reason to deviate, use:

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

**Strand-specific QC.** If a supported directional protocol requires
alignments under one genome conversion hypothesis (`XG:Z`) to be marked
QC-failed, use `--set-as-failed r` (`XG:Z:GA`) or `f` (`XG:Z:CT`). This marks
matching records with `0x200`; it does not remove them, change the input
conversion hypotheses, or enable an unsupported library orientation:

```bash
bwa-mem3 mem --meth --set-as-failed r -t 16 ref.fa R1.fq.gz R2.fq.gz \
  | samtools sort -@ 4 -o out.bam -
```

**Chimera QC is opt-in** (off by default, matching Bismark). Leave it off for
standard directional EM-seq / WGBS. Consider it for chimera-prone supported
directional data, or when you want bwameth.py-style heuristic flagging.
PBAT / scBS-Seq protocols must meet the
[library orientation requirement](../methylation/overview.md#library-orientation).
Pass `--chimera-qc` to enable it (see
[Flags](../methylation/flags.md#--chimera-qc) for the exact heuristic):

```bash
bwa-mem3 mem --meth --chimera-qc -t 16 ref.fa R1.fq.gz R2.fq.gz \
  | samtools sort -@ 4 -o out.bam -
```

> **Note — Overriding scoring defaults**
>
> Scoring flags supplied on the command line override the defaults set by
> `--meth`. For example, `bwa-mem3 mem --meth -B 2 ...` uses `-B 2`, not the
> default `-B 4`. Order does not matter: the defaults are applied after the
> whole command line is parsed, so `-B 2 --meth` and `--meth -B 2` behave
> identically.
>
> This applies to `-B`, `-L`, `-U`, and `-T`. It does **not** apply to `-M` and
> `-C`: bwa has no option that unsets either one, so `--meth` sets them
> unconditionally and they cannot be turned off from the command line. See
> [Flags → `-M` and split alignments](../methylation/flags.md#-m-and-split-alignments).

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
