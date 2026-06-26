# Multi-Sample Workflows

When you need to align many samples back-to-back against the same reference on
a single machine, loading the FM-index into shared memory once — and keeping it
resident across all alignment jobs — eliminates the index I/O cost for every
sample after the first.

## The problem: repeated index loads

The bwa-mem3 index for hg38 is roughly 11 GB on disk (~15 GB resident once
loaded). Without shared memory, `bwa-mem3 mem` reads the entire index from disk
on every invocation.
On a fast NVMe drive this takes 30–60 seconds; on a network-attached or
spinning-disk filesystem it can take several minutes. For a batch of 100
samples, that adds hours of pure I/O overhead.

## Staging the index once with `bwa-mem3 shm`

```bash
# Stage the index into shared memory (one-time cost, ~15 GB for hg38).
bwa-mem3 shm ref.fa

# Align each sample. bwa-mem3 mem attaches automatically — no extra flag.
bwa-mem3 mem --bam=0 -t 16 ref.fa sample1_R1.fq.gz sample1_R2.fq.gz \
  | samtools sort -@ 4 -o sample1.bam -
bwa-mem3 mem --bam=0 -t 16 ref.fa sample2_R1.fq.gz sample2_R2.fq.gz \
  | samtools sort -@ 4 -o sample2.bam -
# ...

# When done, release the segment.
bwa-mem3 shm -d
```

For methylation workflows, stage the c2t index instead:

```bash
bwa-mem3 shm --meth ref.fa
bwa-mem3 mem --meth --bam=0 -t 16 ref.fa sample1_R1.fq.gz sample1_R2.fq.gz \
  | samtools sort -@ 4 -o sample1.bam -
bwa-mem3 shm -d
```

## Confirming the index is staged

```bash
bwa-mem3 shm -l
# Prints the basename and memory usage of each staged segment.
```

If the listing is empty, the index is not staged and `bwa-mem3 mem` will fall
back to loading from disk.

## Thread layout for parallel alignment

Running multiple `bwa-mem3 mem` instances in parallel is efficient when the
samples are independent and the machine has enough cores. The shared-memory
index eliminates disk contention, so the bottleneck becomes CPU and memory
bandwidth.

Guidelines for `N`-core machines:

- **N = 32:** Two instances at `-t 14` each, with `-@ 4` for `samtools sort`.
  Keeps 4 cores reserved for OS and I/O.
- **N = 64:** Two to four instances at `-t 14` to `-t 16`, each with `-@ 4`
  for `samtools sort`.
- **N = 128:** Four to eight instances; keep at least 8–16 cores free for
  `samtools sort` threads and OS scheduling.

> **Tip — Memory bandwidth limit**
>
> The FM-index lookup is memory-bandwidth bound. On machines with NUMA topology
> (multi-socket or multi-chiplet), binding each bwa-mem3 instance to a NUMA
> node with `numactl --cpunodebind=N --membind=N` can improve throughput by
> reducing cross-node memory traffic.

## Scripting a batch with a loop

```bash
bwa-mem3 shm ref.fa

for sample in sample1 sample2 sample3; do
  bwa-mem3 mem --bam=0 -t 16 ref.fa "${sample}_R1.fq.gz" "${sample}_R2.fq.gz" \
    | samtools sort -@ 4 -o "${sample}.bam" -
  samtools index "${sample}.bam"
done

bwa-mem3 shm -d
```

For parallel execution, replace the `for` loop body with a background job (or
use a workflow manager such as Snakemake or Nextflow) and limit the degree of
parallelism to match available cores.

> **Warning — Stale segment footgun**
>
> If you need to re-index the reference (e.g. after updating it), always run
> `bwa-mem3 shm -d` before `bwa-mem3 index`. There is no automatic staleness
> check. See [Anti-patterns](anti-patterns.md) for details.

---

**See also:**
[Quick start: shared-memory index](../getting-started/quick-shm.md) ·
[CLI Reference: shm](../cli/shm.md) ·
[Output format](output-format.md) ·
[Threading and resource use](../user-guide/threading.md) ·
[Anti-patterns](anti-patterns.md)
