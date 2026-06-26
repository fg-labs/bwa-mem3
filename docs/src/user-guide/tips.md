# Tips and best practices

Quick operational pointers for running bwa-mem3. Each links to the page with the
full rationale.

- **Index once, align many.** The on-disk index is stable across bwa-mem3
  releases and across every SIMD tier in the binary; you do not re-index when
  upgrading unless the release notes say so. See
  [Indexing the reference](indexing.md).

- **Pipe `--bam` straight into `samtools sort`** — never write an intermediate
  unsorted BAM. See [Best Practices → Output format](../best-practices/output-format.md).

- **Stage the index in shared memory for batch workloads** (`bwa-mem3 shm`), and
  always `shm -d` before re-indexing — there is no staleness check. See
  [Quick start: shared-memory index](../getting-started/quick-shm.md) and
  [Anti-patterns](../best-practices/anti-patterns.md).

- **Divide threads explicitly across concurrent jobs** so the total stays at or
  below the physical core count. See [Threading and resource use](threading.md).

- **Confirm the resolved SIMD tier** with `bwa-mem3 version` (or
  `BWAMEM3_DEBUG_SIMD=1`); one binary carries every tier and selects at startup.
  See [Performance → SIMD dispatch matrix](../performance/simd-dispatch.md).

- **Always pass a read-group** (`-R` with at least `ID:` and `SM:`) — GATK,
  fgbio, and Picard require an `@RG` header. See [Aligning short reads](aligning.md).

For an impact-ordered tuning walkthrough, see
[Best Practices → Optimization checklist](../best-practices/optimization-checklist.md).

---

**See also:**
[Aligning short reads (mem)](aligning.md) ·
[Threading and resource use](threading.md) ·
[Best Practices: optimization checklist](../best-practices/optimization-checklist.md) ·
[Best Practices: anti-patterns](../best-practices/anti-patterns.md)
