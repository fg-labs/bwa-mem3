# Threading and resource use

## The `-t` flag

```text
-t INT   number of threads [1]
```

bwa-mem3 parallelizes alignment by dividing the input into fixed-size batches
(controlled by `-K`) and processing batches concurrently. Threads share the
in-memory FM-index; there is no per-thread copy.

## How threads interact with performance

### Where threads help

- Seed finding (SMEM enumeration) is fully parallel across reads in a batch.
- Extension (banded Smith-Waterman) is fully parallel.
- Pair rescue is parallel.
- BAM encoding (when `--bam` is active) is parallel.

### Where threads stop helping

Thread count and wall-clock alignment time scale well to approximately 16–32
threads on a modern CPU. Beyond that, several effects conspire to flatten the
curve:

1. **FM-index bandwidth.** The resident index for hg38 is ~15 GB and does not fit
   in the L3 cache of any current server. At high thread counts, threads contend
   for memory bandwidth accessing the BWT.
2. **IO contention.** On spinning disk or a shared network filesystem,
   concurrent reads of the same large index file saturate IO bandwidth before
   the CPU is saturated.
3. **Output serialization.** SAM output is serialized per-record to stdout.
   BAM output with `--bam` reduces this bottleneck but does not eliminate it
   entirely.

### Recommended thread counts

| Machine | Recommended `-t` | Notes |
|---------|-----------------|-------|
| 16-core workstation | 12–14 | Leave 2 cores for `samtools sort` |
| 32-core server | 24–28 | Leave cores for downstream and OS overhead |
| 64-core server | 40–48 | Marginal returns above 48; test with your workload |
| Multiple parallel runs | divide evenly | See below |

These are starting points. Profile with your specific data and storage
configuration to find the practical optimum.

## Running multiple parallel alignments

When running multiple `bwa-mem3 mem` processes on the same machine, divide
threads so that the total does not exceed the physical core count. For example,
on a 32-core machine running four concurrent samples:

```bash
# Four parallel runs, 8 threads each
for sample in a b c d; do
  bwa-mem3 mem --bam -t 8 ref.fa ${sample}_R1.fq.gz ${sample}_R2.fq.gz \
    | samtools sort -@ 2 -o ${sample}.bam - &
done
wait
```

Using shared memory (`bwa-mem3 shm`) amortizes the index read-in cost across
all four runs. See [Quick start: shared-memory index](../getting-started/quick-shm.md)
and [Best Practices: multi-sample workflows](../best-practices/multi-sample.md).

## Memory use

> **`-t` changes the output, not just the speed.** The default batch size is
> `chunk_size × -t`, and each batch's `mem_pestat` insert-size estimate is drawn
> from the reads in that batch — so the same input aligned at different thread
> counts can differ on a small number of records. `bwa` and `bwa-mem2` behave
> identically here. Pass **`-K INT`** to pin the batch size and make output
> independent of `-t`; do this for regression tests, release gating, and any
> comparison against `bwa`/`bwa-mem2`. See [Aligning → `-K`](aligning.md).

Peak RAM is the resident index (~15 GB for hg38, ~22 GB under `--meth`) plus a
per-batch working set that scales with the *effective* batch size
(`chunk_size × n_threads`), and is fixed with respect to `-t`. The per-batch term
is what tips memory-constrained or wide-window (e.g. Hi-C) runs into OOM, and
`bwa-mem3 shm` lets concurrent processes share one physical copy of the index.
For the full budgeting model and the `-K`/`-t` interaction, see
[Memory budgeting and data-type tuning](memory-and-data-types.md).

## IO recommendations

- **Use local NVMe storage** for the index files when possible. The ~11 GB index
  read is the dominant IO event at the start of each `mem` run.
- **Write BAM (`--bam`) to a fast local disk** or pipe directly to
  `samtools sort`. Avoid writing uncompressed SAM to a network filesystem.
- **Separate read and write paths** if your storage topology allows it:
  read the index from one volume and write sorted BAM to another.

---

**See also:**
[Aligning short reads (mem)](aligning.md) ·
[Memory budgeting and data-type tuning](memory-and-data-types.md) ·
[Memory allocator (mimalloc)](allocator.md) ·
[Quick start: shared-memory index](../getting-started/quick-shm.md) ·
[Best Practices: multi-sample workflows](../best-practices/multi-sample.md) ·
[Best Practices: optimization checklist](../best-practices/optimization-checklist.md)
