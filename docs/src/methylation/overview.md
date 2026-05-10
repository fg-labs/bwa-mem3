# Methylation Reference Overview

`bwa-mem3 --meth` is a single-binary, single-command drop-in replacement for
the [bwameth.py](https://github.com/brentp/bwa-meth) bisulfite-sequencing
alignment pipeline. No Python installation, no piped preprocessing step, and no
separate post-processing script — one `bwa-mem3 index --meth` builds the
reference, and one `bwa-mem3 mem --meth` aligns and post-processes reads from
raw FASTQ to sorted-ready BAM.

The output BAM is structurally equivalent to what the bwameth.py pipeline
produces: consolidated `@SQ` headers (one entry per real chromosome rather
than one per doubled-reference contig), Bismark-compatible `XR:Z` (read
conversion `CT`/`GA`), `XG:Z` (genome strand `CT`/`GA`), and `XM:Z`
(per-base methylation call string) auxiliary tags, optional chimera QC flags
(`--chimera-qc`, off by default to match Bismark), and a
`@PG ID:bwa-mem3-meth` provenance entry. Every Bismark-native tool
(`bismark_methylation_extractor`, methylKit, methtuple, DMRfinder,
epialleleR), MethylDackel, and biscuit's per-read methylation tools read
the BAM directly without conversion.

## Pipeline at a glance

The diagram below shows the internal flow when `bwa-mem3 mem --meth` runs.
Every step executes inside the single process; no external programs or temporary
files are required.

```mermaid
flowchart LR
    A[Raw FASTQ\nR1 / R2] -->|inline C→T / G→A| B[c2t-converted reads\n+ internal YS/YC carrier]
    B -->|bwa mem core| C[mem_aln_t\nalignments vs doubled ref]
    C -->|chrom map\nf/r → real chr| D[header rewrite\n@SQ consolidated]
    D -->|XR/XG/XM Bismark tags\noptional --chimera-qc\nQC-fail propagation| E[BAM output\nwb0 uncompressed]
```

Steps:

1. **FASTQ ingest with inline c2t conversion.** R1 bases have every `C`
   replaced with `T`; R2 bases have every `G` replaced with `A`. The
   original bases and conversion direction are kept on an internal carrier
   on each read (in `bseq1_t.comment`); they are never emitted to BAM as
   tags themselves but feed the BAM-write step (SEQ restoration, `XR:Z`
   derivation). This conversion happens in-memory — the FASTQ is never
   written to disk in converted form.

2. **Alignment against the doubled reference.** The converted reads are aligned
   against the `ref.fa.bwameth.c2t` reference, which contains both a forward
   C→T projection (`f`-prefixed contigs) and a reverse G→A projection
   (`r`-prefixed contigs) of each chromosome.

3. **Header rewriting and chrom consolidation.** The `f`/`r`-prefixed contig
   names used internally are collapsed: every pair `fchr1` / `rchr1` becomes a
   single `@SQ SN:chr1` entry in the output BAM header. RNAME and RNEXT fields
   in each record are rewritten to the consolidated name.

4. **Tag emission and QC.** Each aligned record receives Bismark-compatible
   `XR:Z` (read conversion direction), `XG:Z` (genome strand), and `XM:Z`
   (per-base methylation call string) auxiliary tags. With opt-in
   `--chimera-qc` (off by default — matches Bismark), records whose longest
   M/=/X CIGAR run covers less than 44 % of the read length are flagged
   `0x200`; QC-fail flags then propagate across all records in a read group.
   The original pre-c2t sequence is copied back into the BAM SEQ field so
   methylation callers see real cytosines rather than the converted
   sequence.

5. **BAM output.** Records are written as uncompressed BAM (`wb0` mode via
   htslib). The `@PG ID:bwa-mem3-meth` line records the exact command line.
   The caller pipes directly to `samtools sort`.

## Quick-start commands

```sh
# Index the reference once (builds ref.fa.bwameth.c2t + FMI)
bwa-mem3 index --meth ref.fa

# Align paired-end FASTQs
bwa-mem3 mem --meth -t 16 ref.fa R1.fq.gz R2.fq.gz \
  | samtools sort -o out.bam
samtools index out.bam
```

> **Note — bwameth.py compatibility**
>
> The default scoring parameters applied by `--meth` (`-B 2 -L 10 -U 100 -T 40 -CM`)
> match those used by bwameth.py so outputs are comparable. Any parameter can be
> overridden on the command line.

---

**See also:**
[bwameth.py drop-in mapping](bwameth-mapping.md) ·
[Conversion details](conversion.md) ·
[SAM tags: XR, XG, XM](tags.md) ·
[Chimera QC and header rewriting](post-processing.md) ·
[Quick start: methylation alignment](../getting-started/quick-meth.md)
