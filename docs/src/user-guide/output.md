# Output: SAM/BAM, headers, tags

bwa-mem3 writes output in either SAM (default) or BAM (`--bam`) format.
This page covers the header structure and every non-standard SAM tag emitted
by bwa-mem3.

## Output format

By default, `bwa-mem3 mem` writes SAM to stdout. Pass `--bam` (or `--bam=N`
for a specific compression level) to write BAM. Level 0 (uncompressed) is the
default when `--bam` is given without an argument, which is optimal when piping
to a downstream `samtools sort`.

```sh
# SAM (default)
bwa-mem3 mem -t 16 ref.fa R1.fq.gz R2.fq.gz > out.sam

# Uncompressed BAM — best for piping
bwa-mem3 mem --bam -t 16 ref.fa R1.fq.gz R2.fq.gz | samtools sort -@ 8 -o out.bam -

# Compressed BAM — useful when the output is the final file
bwa-mem3 mem --bam=6 -t 16 ref.fa R1.fq.gz R2.fq.gz > out.bam
```

## SAM header

### `@HD`

A default `@HD VN:1.6 SO:unsorted` line is emitted unless the user supplies
one via `-H`. The sort order is `unsorted` because bwa-mem3 writes records in
input read order; downstream sorting is always a separate step.

### `@SQ`

One `@SQ` line is written per reference sequence, with the sequence name
(`SN:`) and length (`LN:`) derived from the FM-index. If the index was built
with a `.dict` or `.hdr` file that supplies `@SQ` records, those records are
used instead of the auto-generated ones.

In methylation mode (`--meth`), the doubled reference contains sequences with
an `f` or `r` prefix in their names. The inline BAM post-processor collapses
these back to canonical chromosome names so that the output `@SQ` lines match
a standard non-methylation alignment. See
[Chimera QC and header rewriting](../methylation/post-processing.md).

### `@PG`

One `@PG` entry is written in standard mode:

| ID | Description |
|----|-------------|
| `bwa-mem3` | The alignment step. `VN:` is the bwa-mem3 version string; `CL:` is the full command line. |

In methylation mode (`--meth`), a second `@PG` entry is appended:

| ID | Description |
|----|-------------|
| `bwa-mem3-meth` | The inline post-processor. `VN:` carries the version with `-meth` suffix; `CL:` is the full command line. |

The `bwa-mem3-meth` entry follows immediately after the `bwa-mem3` entry and
records the post-processing step as a distinct pipeline node, matching the
convention of separate-tool pipelines.

## Tags emitted by bwa-mem3

### Standard tags

bwa-mem3 emits the same standard tags as bwa-mem2 (`NM:i`, `MD:Z`, `AS:i`,
`XS:i`, `SA:Z`, `RG:Z`, `XA:Z`, `MC:Z`, etc.). These are documented in the
SAM specification and are not described further here.

bwa-mem3 additionally emits `MQ:i` on paired-end records — the mate's
mapping quality, set alongside `MC:Z` (the mate's CIGAR) so callers that
key off the mate's MAPQ don't need to look at the mate record. Both SAM
and `--bam` output paths emit it. Backported from `lh3/bwa` PR #330 in
fg-labs PR #35.

The `XA:Z` field set widens from `chr,pos,CIGAR,NM` to
`chr,pos,CIGAR,NM,score,mapq` when `-u` (a.k.a. the upstream "XB" toggle)
is passed; the tag name itself remains `XA:Z` for downstream
compatibility. Tools that parse `XA:Z` need to be aware of the two
possible field widths.

### `HN:i` — total alignment hit count

```text
HN:i:<count>
```

The total number of primary alignments (both reported and suppressed) that
the aligner found for this read, before the `-h` supplementary cap is applied.
Useful for distinguishing "uniquely mapped" from "multi-mapped" reads without
relying solely on MAPQ.

`HN:i` is emitted on the primary alignment record only.

### Methylation-only tags

The following Bismark-compatible tags are emitted only when `--meth` is
active. See [SAM tags: XR, XG, XM](../methylation/tags.md) for the full
per-tag reference, including the `XM:Z` character alphabet and the
`XG:Z` strand-pick semantics.

| Tag | Type | Description |
|-----|------|-------------|
| `XR:Z` | string | Read conversion direction: `CT` (R1 / SE) or `GA` (R2) |
| `XG:Z` | string | Genome strand of the alignment: `CT` (OT) or `GA` (OB) |
| `XM:Z` | string | Per-base methylation call string (length = `SEQ`) |

The bwameth-style `YS:Z` / `YC:Z` tags exist only as an internal carrier
on `bseq1_t.comment` for SEQ restoration and `XR:Z` derivation; they
are suppressed at BAM emit and never appear in output. The bwameth
`YD:Z` strand tag has been replaced by Bismark `XG:Z` and is not
emitted.

## MAPQ semantics

MAPQ semantics are inherited from bwa-mem2 and follow the same scoring model.
In methylation mode, alignments identified as chimeras (longest `M`/`=`/`X`
run covering less than 44% of the read length) have their MAPQ capped at 1 and
the `0x200` (QC fail) flag set. See
[Chimera QC and header rewriting](../methylation/post-processing.md).

---

**See also:**
[Aligning short reads (mem)](aligning.md) ·
[Methylation Reference: SAM tags](../methylation/tags.md) ·
[Methylation Reference: post-processing](../methylation/post-processing.md) ·
[CLI Reference: mem](../cli/mem.md) ·
[Best Practices: output format](../best-practices/output-format.md)
