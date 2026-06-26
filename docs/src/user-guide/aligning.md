# Aligning short reads (mem)

`bwa-mem3 mem` aligns one or two FASTQ files against an indexed reference and
writes SAM (default) or BAM (`--bam`) to stdout. It is a drop-in replacement
for `bwa-mem2 mem` and supports all standard bwa-mem flags.

## Basic usage

**Paired-end:**

```bash
bwa-mem3 mem -t 16 ref.fa R1.fq.gz R2.fq.gz > out.sam
```

**Single-end:**

```bash
bwa-mem3 mem -t 16 ref.fa reads.fq.gz > out.sam
```

**Pipe directly to samtools:**

```bash
bwa-mem3 mem --bam=0 -t 16 ref.fa R1.fq.gz R2.fq.gz \
  | samtools sort -@ 8 -o out.bam -
samtools index out.bam
```

Using `--bam=0` (uncompressed BAM) avoids SAM text formatting on the write side
and SAM parsing on the samtools side, and skips the wasted compression that
`samtools sort` would immediately decompress; the BAM bytes flow between
processes in the pipe.

## Key flags

### Threading: `-t`

```text
-t INT   number of threads [1]
```

Performance scales well through 8–16 threads on most machines. Beyond 32
threads, returns diminish on typical workloads because inter-thread locking
and IO become the bottleneck. See [Threading and resource use](threading.md)
for detailed guidance.

### Read-group header: `-R`

```text
-R STR   read group header line, e.g. '@RG\tID:sample1\tSM:sample1\tLB:lib1\tPL:ILLUMINA'
```

Every production alignment should include a `@RG` header. The ID in the `-R`
string is embedded as an `RG:Z:` tag on every output record.

> **Tip — Escape the tab correctly**
>
> Pass `-R` with a literal `\t` between fields. Most shells require single quotes
> or `$'...'` quoting to prevent interpretation of the backslash:
>
> ```bash
> bwa-mem3 mem -R $'@RG\tID:s1\tSM:sample1' -t 16 ref.fa R1.fq.gz R2.fq.gz
> ```

### Chunk size: `-K`

```text
-K INT   process INT input bases in each batch [10000000]
```

Larger `-K` values increase memory use but can improve throughput on very deep
or very wide batches. The default is appropriate for most workloads.

### SAM output control: `-S`, `-P`

```text
-S    skip mate rescue
-P    skip pairing; mate rescue performed unless -S also in use
```

These flags are primarily useful for debugging or non-standard workflows.
Normal paired-end alignments should leave both at their defaults.

## Output modes

### SAM (default)

```bash
bwa-mem3 mem -t 16 ref.fa R1.fq.gz R2.fq.gz > out.sam
```

Plain-text SAM. Suitable for inspection, compatibility testing, and piping to
tools that consume SAM.

### BAM (`--bam=0`)

```bash
bwa-mem3 mem --bam=0 -t 16 ref.fa R1.fq.gz R2.fq.gz > out.bam
```

Writes BAM directly. `--bam=0` is uncompressed BAM, which avoids
double-compression when piping into a downstream sorter and is roughly
10–15% faster end-to-end. Pass `--bam=6` to write a fully compressed BAM if
the output is the final product.

> **Note — --bam=0 is the recommended output mode**
>
> For production pipelines, always use `--bam=0` and pipe to `samtools sort`.
> See [Best Practices: output format](../best-practices/output-format.md) for
> the canonical pipeline.

## Methylation alignment (`--meth`)

Pass `--meth` for bisulfite/RRBS samples. This activates inline C-to-T read
conversion, bwameth.py-compatible flag defaults, and inline BAM post-processing.
See [Quick start: methylation alignment](../getting-started/quick-meth.md) for
the two-command workflow and the [Methylation Reference](../methylation/overview.md)
for full detail.

## Shared-memory index auto-attach

When `bwa-mem3 shm` has staged the index into shared memory, `bwa-mem3 mem`
attaches automatically — no extra flag is required. The shared-memory path
is transparent to users.

## Cross-references

The full flag list is in the [CLI Reference: mem](../cli/mem.md) page.

---

**See also:**
[Output: SAM/BAM, headers, tags](output.md) ·
[Threading and resource use](threading.md) ·
[Best Practices: output format](../best-practices/output-format.md) ·
[CLI Reference: mem](../cli/mem.md) ·
[Methylation Reference: overview](../methylation/overview.md)
