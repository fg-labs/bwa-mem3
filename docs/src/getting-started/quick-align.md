# Quick start: align paired-end FASTQs

This page walks through the two-command workflow: index the reference once, then align reads.

> **Switching from `bwa` or `bwa-mem2`?** The command line is unchanged, but
> rebuild the index if you're coming from `bwa` v1 (a `bwa-mem2` index is reused
> as-is), and expect output that is functionally equivalent but **not**
> byte-identical. See [Coming from bwa or bwa-mem2](migrating.md).

## Index the reference

```bash
bwa-mem3 index ref.fa
```

This produces four index files alongside `ref.fa`:

| File | Description |
|------|-------------|
| `ref.fa.bwt.2bit.64` | FM-index in 2-bit packed format |
| `ref.fa.amb` | Ambiguous base positions |
| `ref.fa.ann` | Sequence name and length annotations |
| `ref.fa.pac` | 2-bit packed reference sequence |

(No `ref.fa.0123` is written: `mem` reconstructs reference bases from `.pac` on
demand. Pass `index --emit-unpacked-ref` if a tool such as bwa-mem2 needs it.)

Indexing hg38 takes a few minutes on a multi-core host (construction is
multi-threaded and memory-bounded — see [Indexing the reference](../user-guide/indexing.md));
the final index is roughly 11 GB on disk. The index is read once per `mem` invocation; for workloads that align many
samples, load it into shared memory first (see [Quick start: shared-memory index](quick-shm.md)).

## Align paired-end reads

```bash
bwa-mem3 mem -t 16 ref.fa r1.fq.gz r2.fq.gz > out.sam
```

`-t 16` sets the thread count to 16. bwa-mem3 scales well up to the number of physical CPU
cores; hyperthreading provides diminishing returns above that point. See
[User Guide — Threading and resource use](../user-guide/threading.md) for recommendations at
different core counts.

The default output is SAM on stdout. To write BAM directly, add `--bam`
(uncompressed by default — best for piping to `samtools sort`; use `--bam=6` for
BGZF-compressed BAM):

```bash
bwa-mem3 mem --bam -t 16 ref.fa r1.fq.gz r2.fq.gz \
  | samtools sort -@ 8 -o out.bam -
samtools index out.bam
```

> **Tip — Prefer BAM output in production**
>
> Piping BAM (`--bam`) to `samtools sort` avoids the text formatting and parsing overhead of SAM
> on both sides of the pipe. For large cohorts this yields a measurable wall-clock reduction.
> See [Best Practices — Output format](../best-practices/output-format.md) for the recommended
> pipeline and a discussion of when SAM is still useful.

## Read group tagging

For downstream tools that require a `@RG` header (most variant callers), pass `-R`:

```bash
bwa-mem3 mem -t 16 \
  -R '@RG\tID:sample1\tSM:sample1\tPL:ILLUMINA\tLB:lib1' \
  ref.fa r1.fq.gz r2.fq.gz > out.sam
```

The value is a tab-delimited string following BWA conventions. Every aligned record receives an
`RG:Z:` tag matching the `ID` field of the read-group header.

bwa-mem3 emits the standard SAM tags plus the fork's `HN:i` tag; for the full
list see [Output: SAM/BAM, headers, tags](../user-guide/output.md) (and
[Methylation SAM tags](../methylation/tags.md) under `--meth`).

---

**See also:**
[Quick start: methylation alignment](quick-meth.md) ·
[Quick start: shared-memory index](quick-shm.md) ·
[User Guide — Aligning short reads](../user-guide/aligning.md) ·
[CLI Reference — mem](../cli/mem.md)
