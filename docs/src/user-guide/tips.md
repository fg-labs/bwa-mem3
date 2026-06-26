# Tips and best practices

This page collects the most commonly useful operational tips for running
bwa-mem3. Each tip is a short actionable point; the linked pages provide the
full rationale.

## Index once, align many times

Build the FM-index once per reference version. The on-disk index format is
stable across bwa-mem3 releases and across every SIMD tier inside the
single binary — the AVX2 and AVX-512BW kernel paths read the same files.
You do not need to re-index when upgrading bwa-mem3 unless the release
notes say otherwise.

```sh
# Build once
bwa-mem3 index ref.fa

# Align many samples
for sample in a b c d; do
  bwa-mem3 mem --bam -t 16 ref.fa ${sample}_R1.fq.gz ${sample}_R2.fq.gz \
    | samtools sort -@ 4 -o ${sample}.bam -
done
```

## Pipe to `samtools sort -@`

Never write an intermediate unsorted BAM to disk and then sort it in a second
step. bwa-mem3's `--bam` mode + `samtools sort` in a single pipeline avoids the
extra write/read cycle and is significantly faster:

```sh
bwa-mem3 mem --bam -t 16 ref.fa R1.fq.gz R2.fq.gz \
  | samtools sort -@ 8 -o out.bam -
samtools index out.bam
```

Allocate roughly `2/3` of available threads to `bwa-mem3 mem` and `1/3` to
`samtools sort`. On a 24-core machine, `-t 16` for bwa-mem3 and `-@ 8` for
samtools is a good starting point.

## Stage the index in shared memory for batch workloads

When aligning more than a few samples on the same machine, reading the ~28 GB
hg38 index from disk on every `mem` invocation is the dominant wall-clock cost.
Stage it once:

```sh
bwa-mem3 shm ref.fa
```

Subsequent `bwa-mem3 mem` invocations attach automatically. The shared-memory
segment persists until explicitly dropped (`bwa-mem3 shm -d`) or the machine
reboots.

> **Warning — Always drop the segment before re-indexing**
>
> There is no staleness check. If you rebuild the index without first dropping
> the shared-memory segment, `bwa-mem3 mem` will attach to the stale segment and
> produce incorrect alignments without any warning. Always run `bwa-mem3 shm -d`
> before `bwa-mem3 index`.

## Pin threads when running concurrent jobs

When running multiple `bwa-mem3 mem` processes in parallel, divide threads
explicitly so that the total does not exceed the physical core count. Avoid
relying on the scheduler to balance over-subscribed threads — each process
will spin waiting for CPU time, and total throughput drops.

```sh
# Good: 4 jobs × 6 threads = 24 cores, on a 24-core machine
for sample in a b c d; do
  bwa-mem3 mem --bam -t 6 ref.fa ${sample}_R1.fq.gz ${sample}_R2.fq.gz \
    | samtools sort -@ 2 -o ${sample}.bam - &
done
wait
```

See [Threading and resource use](threading.md) for per-machine thread count
recommendations.

## Confirm the binary's SIMD tier matches your CPU

bwa-mem3 ships **one binary per platform** that contains every supported
x86 SIMD tier (or the single NEON path on arm64) and picks the right
tier in process at startup. There are no per-tier companion binaries to
copy or call directly.

| CPU generation | Resolved tier |
|---------------|---------------|
| Modern Intel/AMD (2018+) | `avx512bw` or `avx2` |
| Older x86 | `sse42` or `sse41` |
| Apple Silicon / AWS Graviton | `neon` |

Verify the resolved tier with `bwa-mem3 version` (prints `SIMD floor:` and
`SIMD runtime:` lines on stdout) or set `BWAMEM3_DEBUG_SIMD=1` to get a
startup banner from `bwa-mem3 mem`. If you need to force a lower tier
for A/B regression testing, set `BWAMEM3_FORCE_TIER=<tier>` — upgrade
requests above the host's capability are rejected.

See [Performance: SIMD dispatch matrix](../performance/simd-dispatch.md).

## Include a read-group header

Always pass `-R` with at minimum `ID:` and `SM:` fields. Many downstream tools
(GATK, fgbio, Picard) require a `@RG` header and will fail or warn without one.

```sh
bwa-mem3 mem \
  -R $'@RG\tID:run1\tSM:sample1\tLB:lib1\tPL:ILLUMINA' \
  -t 16 ref.fa R1.fq.gz R2.fq.gz
```

## Further reading

The Best Practices section covers these topics in depth:

- [Best Practices: build](../best-practices/build.md) — PGO builds, arch selection
- [Best Practices: output format](../best-practices/output-format.md) — the canonical pipeline
- [Best Practices: multi-sample workflows](../best-practices/multi-sample.md) — shared-memory batch jobs
- [Best Practices: anti-patterns](../best-practices/anti-patterns.md) — common mistakes and how to avoid them

---

**See also:**
[Aligning short reads (mem)](aligning.md) ·
[Threading and resource use](threading.md) ·
[Memory allocator (mimalloc)](allocator.md) ·
[Best Practices: optimization checklist](../best-practices/optimization-checklist.md) ·
[Best Practices: anti-patterns](../best-practices/anti-patterns.md)
