# Tuning Checklist

The items below are ordered by expected impact for most workloads. Work through them in sequence; there is little point optimizing output format before confirming you are running the right binary for your CPU.

## 1. Confirm the resolved SIMD tier matches your CPU

The default `make` produces a single binary that contains every supported
x86 SIMD tier and selects one in process at startup. Verify which tier
is running:

```bash
bwa-mem3 version
# expect: SIMD floor: <build_floor>; SIMD runtime: <resolved_tier>
```

If the runtime tier is below what your CPU supports, double-check
whether you accidentally built with a lower `BASELINE_ARCH=` or set
`BWAMEM3_FORCE_TIER` in the environment. Set `BWAMEM3_DEBUG_SIMD=1` to
get a startup banner on stderr at the start of a `mem` run.

On ARM / Apple Silicon, the binary has one NEON tier; `bwa-mem3 version`
reports `SIMD runtime: neon`.

See [SIMD dispatch matrix](simd-dispatch.md) for the full dispatch
logic and the minimum CPU requirements for each tier.

> **Tip — Single-arch deployments**
>
> On a cluster where every node has the same CPU, build with
> `make arch=avx2` (or the appropriate ISA). The runtime dispatch
> overhead is negligible, but a single-arch build trims the binary
> and removes any chance of `BWAMEM3_FORCE_TIER` accidentally
> downgrading throughput in production.

## 2. Build with PGO if you will run repeatedly

For production pipeline nodes that will process many samples against the same reference, a PGO build provides an additional 2–5% throughput at the cost of one extra build pass and a training run:

```bash
make pgo-generate PGO_ARCH=avx2
./bwa-mem3.pgo-instr.avx2 mem -t 16 ref.fa R1.fq.gz R2.fq.gz > /dev/null
make pgo-use PGO_ARCH=avx2
```

See [PGO build](pgo.md) for the full workflow, including multi-arch and profile portability notes.

## 3. Use shared memory for many small samples

When aligning many samples on one machine against the same reference, loading the index into POSIX shared memory once and reusing it across all `mem` invocations eliminates redundant I/O and reduces per-sample startup time significantly. The benefit grows with the number of samples and the size of the reference.

```bash
# Load the index into shared memory once
bwa-mem3 shm ref.fa

# Align each sample against the in-memory index
bwa-mem3 mem -t 16 ref.fa R1.fq.gz R2.fq.gz | samtools sort -@ 4 -o sample.bam -

# When finished with all samples, drop the shared segment
bwa-mem3 shm -d
```

> **Warning — No staleness check**
>
> `bwa-mem3 shm` does not detect whether the on-disk index has changed after the segment was loaded. Always run `bwa-mem3 shm -d` before re-indexing a reference and re-loading with `bwa-mem3 shm`. Failing to do so results in alignments against a stale index.

See [Getting Started — Shared-memory index](../getting-started/quick-shm.md) and [Best Practices — Multi-sample workflows](../best-practices/multi-sample.md) for complete workflows.

## 4. Emit BAM directly

Use `--bam` (or `--bam=0` for uncompressed BAM) to emit BAM instead of SAM. Uncompressed BAM avoids the text-formatting cost on the aligner side and the text-parsing cost on the downstream side. `samtools sort` reads BAM natively and is fastest when the input is uncompressed:

```bash
bwa-mem3 mem --bam=0 -t 16 ref.fa R1.fq.gz R2.fq.gz \
  | samtools sort -@ 8 -o out.bam -
samtools index out.bam
```

The `--bam` flag (without `=0`) produces BGZF-compressed BAM. This is useful when writing directly to disk without a downstream piped tool.

See [Best Practices — Output format](../best-practices/output-format.md) for guidance on when SAM is still appropriate.

## 5. Pipe to a multi-threaded sorter

Sorting is typically the bottleneck after alignment. Keep a separate thread budget for `samtools sort`:

```bash
bwa-mem3 mem --bam=0 -t 12 ref.fa R1.fq.gz R2.fq.gz \
  | samtools sort -@ 8 -m 2G -o out.bam -
```

On a 16-core machine, allocating 12 threads to `mem` and 8 to `samtools sort` (with overlap via the pipe) is a common starting point. The aligner is generally CPU-bound; the sorter is I/O-bound during merge. Profile both stages to find the right split for your hardware.

> **Tip — Thread count tuning**
>
> `bwa-mem3 mem` scales well to 16–32 threads on most workloads. Beyond 32 threads the per-thread work unit becomes small enough that synchronization overhead starts to erode gains. See [User Guide — Threading and resource use](../user-guide/threading.md) for thread-scaling data.

## Summary table

| Item | Action | Reference |
|---|---|---|
| Right SIMD tier for CPU | `bwa-mem3 version`; verify `SIMD runtime:` | [SIMD dispatch matrix](simd-dispatch.md) |
| PGO for production | `pgo-generate` → train → `pgo-use` | [PGO build](pgo.md) |
| Shared-memory index | `bwa-mem3 shm ref.fa` before batch runs | [Quick start: shm](../getting-started/quick-shm.md) |
| Emit uncompressed BAM | `--bam=0` | [Best Practices — Output format](../best-practices/output-format.md) |
| Multi-threaded sort | `samtools sort -@` with appropriate thread split | [User Guide — Threading](../user-guide/threading.md) |

---

**See also:**
[Performance overview](overview.md) ·
[SIMD dispatch matrix](simd-dispatch.md) ·
[PGO build](pgo.md) ·
[Best Practices — Build](../best-practices/build.md) ·
[User Guide — Threading and resource use](../user-guide/threading.md)
