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

1. **FM-index bandwidth.** The index for hg38 is ~28 GB and does not fit in
   the L3 cache of any current server. At high thread counts, threads contend
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

```sh
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

Peak RAM during alignment is dominated by the in-memory FM-index. For hg38,
expect roughly 28 GB of resident memory per `bwa-mem3 mem` process. Additional
memory is used per batch (`-K` reads × read length × a small constant).

With `bwa-mem3 shm`, the index is mapped from a shared-memory segment, so
multiple concurrent `mem` processes share the same physical pages. The OS
deduplicates the pages; total RAM use is approximately one index, not one per
process.

> **Tip — Use shm for repeated runs on the same machine**
>
> If you run more than a few samples on the same machine without rebooting,
> `bwa-mem3 shm` pays off immediately. The index is read from disk once and
> stays in RAM for all subsequent `mem` invocations.

## IO recommendations

- **Use local NVMe storage** for the index files when possible. The ~28 GB BWT
  read is the dominant IO event at the start of each `mem` run.
- **Write BAM (`--bam`) to a fast local disk** or pipe directly to
  `samtools sort`. Avoid writing uncompressed SAM to a network filesystem.
- **Separate read and write paths** if your storage topology allows it:
  read the index from one volume and write sorted BAM to another.

---

**See also:**
[Aligning short reads (mem)](aligning.md) ·
[Memory allocator (mimalloc)](allocator.md) ·
[Quick start: shared-memory index](../getting-started/quick-shm.md) ·
[Best Practices: multi-sample workflows](../best-practices/multi-sample.md) ·
[Performance: tuning checklist](../performance/tuning.md)
