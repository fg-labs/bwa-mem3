# Quick start: shared-memory index

The bwa-mem3 hg38 index is roughly 11 GB on disk (~15 GB once loaded into RAM). By default, every
`bwa-mem3 mem` invocation reads the index from disk, which can take 30–60 seconds on a spinning
disk and several seconds even on fast NVMe storage. For workloads that align many small samples
in sequence on the same machine, this per-invocation overhead accumulates.

`bwa-mem3 shm` stages the index once into POSIX shared memory. Subsequent `mem` invocations
attach to the in-memory segment instead of reading from disk, reducing per-sample startup time
to near zero.

## Stage the index

```bash
bwa-mem3 shm ref.fa
```

This reads the index files from disk and copies them into a POSIX shared-memory segment. The
command returns when staging is complete. The index stays in memory until it is explicitly
dropped or the system is rebooted.

To stage a methylation (`--meth`) index:

```bash
bwa-mem3 shm --meth ref.fa
```

A standard and a methylation index for the same reference can be staged simultaneously; they
occupy separate named segments.

## Align using the staged index

No extra flag is needed. When `bwa-mem3 mem` starts, it checks whether a matching shared-memory
segment exists. If one does, it attaches automatically:

```bash
bwa-mem3 mem -t 16 ref.fa R1.fq.gz R2.fq.gz > out.sam
```

## Inspect and drop staged segments

List all currently staged indices:

```bash
bwa-mem3 shm -l
```

Drop all staged segments:

```bash
bwa-mem3 shm -d
```

## When to use shared-memory indexing

Shared-memory indexing is most beneficial when:

- Aligning tens to hundreds of small samples (e.g. amplicon panels, targeted sequencing) where
  per-sample read time dominates the per-sample alignment time.
- Running a batch pipeline on a single large machine where the index fits comfortably in RAM
  (approximately 15 GB for hg38 with the standard index).
- The same reference is used for all samples in the batch; a new `shm` invocation is required
  for each distinct reference.

It provides little benefit when:

- Aligning a small number of large samples (WGS), where alignment time far exceeds index
  load time.
- The available RAM is insufficient to hold the index alongside the operating system and
  alignment worker processes.

> **Warning — No staleness check — always drop before re-indexing**
>
> `bwa-mem3 shm` does not detect whether the on-disk index files have changed after staging. If
> you run `bwa-mem3 index ref.fa` again (e.g. to rebuild after a reference update), the
> shared-memory segment is not invalidated. Subsequent `mem` invocations will attach to the stale
> segment and produce silently incorrect alignments.
>
> Always drop the segment before re-indexing:
>
> ```bash
> bwa-mem3 shm -d
> bwa-mem3 index ref.fa
> bwa-mem3 shm ref.fa
> ```

---

**See also:**
[CLI Reference — shm](../cli/shm.md) ·
[Best Practices — Multi-sample workflows](../best-practices/multi-sample.md) ·
[Best Practices — Anti-patterns](../best-practices/anti-patterns.md) ·
[User Guide — Threading and resource use](../user-guide/threading.md)
