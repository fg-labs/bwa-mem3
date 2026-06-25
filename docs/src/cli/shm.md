# shm

`bwa-mem3 shm` stages an FM-index into POSIX shared memory so that subsequent
`bwa-mem3 mem` invocations on the same machine attach to the in-memory segment
instead of re-reading the index files from disk. For workloads that align many
small samples back-to-back against the same reference — such as clinical
panels or amplicon sequencing — this removes the dominant I/O bottleneck.
`shm` also lists and destroys staged segments.

## Synopsis

```text
{{#include ../../_generated/cli/shm.txt}}
```

## Common usage

Stage a standard index, align two samples, then release the segment:

```sh
bwa-mem3 shm ref.fa
bwa-mem3 mem -t 16 ref.fa sample1_R1.fq sample1_R2.fq > sample1.sam
bwa-mem3 mem -t 16 ref.fa sample2_R1.fq sample2_R2.fq > sample2.sam
bwa-mem3 shm -d
```

Stage a methylation index and align:

```sh
bwa-mem3 shm --meth ref.fa
bwa-mem3 mem --meth -t 16 ref.fa R1.fq R2.fq | samtools sort -o out.bam -
bwa-mem3 shm -d
```

List all currently staged segments:

```sh
bwa-mem3 shm -l
```

## Flag reference

### (no flags) `<idxbase>` — stage an index

Loads all index files for `<idxbase>` into a POSIX shared-memory segment.
After staging, any `bwa-mem3 mem <idxbase> ...` on the same machine
auto-attaches and reads from memory rather than disk.

### `-d` — destroy all segments

Removes every bwa-mem3 shared-memory segment on the machine. This is the
correct clean-up command after a batch job and the required step before
re-building the index (see the footgun warning below).

### `-l` — list staged indices

Prints the names of all currently staged segments. Useful to confirm that
staging succeeded before launching alignment jobs.

### `--meth` — stage a methylation index

Stages the `.meth` seed index (`<idxbase>.meth.*`) into shared memory, mirroring
the behavior of `bwa-mem3 index --meth` and `bwa-mem3 mem --meth`. Pass the same
plain `<idxbase>` to all three commands; the `.meth` suffix is handled
transparently.

The staged seed segment is **seed-only**: it holds the seed FM-index and contig
metadata (BNS) but omits the packed reference (PAC), because `mem --meth` extends
against the original reference and never reads the seed's bases. This trims
~1.6 GB from the staged segment on hg38. (The unpacked `.0123` reference is never
staged for **any** index — plain or seed — because `mem` pac-fetches reference
bases from `.pac` on demand; that saves the seed's ~13 GB and the original's
~6.4 GB versus staging `.0123`.)

## Notes / Gotchas

> **Warning — No staleness check — always destroy before re-indexing**
>
> There is no staleness check. If you re-run `bwa-mem3 index ref.fa` after
> staging, the on-disk index files will not match the in-memory segment, but
> `bwa-mem3 mem` will still attach to the stale segment and silently produce
> incorrect alignments. Always run `bwa-mem3 shm -d` before re-indexing.
>
> **Note — Platform limits**
>
> **macOS:** POSIX shared memory has implementation-defined per-segment size
> caps. Staging a full hg38 index (~18 GB; ~21 GB for a `--meth` seed segment)
> may fail silently or with a
> cryptic error. Prefer Linux for production use with large references.
>
> **Linux containers:** `/dev/shm` typically defaults to ~50% of physical RAM
> on bare metal but is often much smaller inside Docker containers or
> Kubernetes pods. Raise the limit with `--shm-size` (Docker) or an
> `emptyDir` tmpfs volume with an explicit size (Kubernetes) before attempting
> to stage a large index.
>
> **Note — `/dev/shm` capacity preflight (PR #86)**
>
> Before opening the segment, `bwa-mem3 shm` calls `statvfs("/dev/shm")` and
> compares the available bytes against the index's `total_size`. If `/dev/shm`
> is too small the stage aborts cleanly with an `[E::bwa_shm_stage]` message
> that names `/dev/shm`, the required size, and a `mount -o remount,size=...`
> hint. This replaces the previous failure mode where `ftruncate` succeeded
> lazily and `pack_into` later surfaced ENOSPC as `[fread] Bad address` with
> no indication that `/dev/shm` was the cause. The preflight is best-effort:
> a `statvfs` failure (no `/dev/shm`, restricted sandbox, ENOSYS) is
> non-fatal and the stage proceeds. As a rough sizing guide, hg38 stages
> ~17 GB; AWS instances default to RAM/2 (so c7a.4xlarge / c7i.4xlarge at
> 32 GB get ~16 GB of `/dev/shm`, which is just under the index size — a
> `remount,size=28g` is the documented fix).
>
> **Note — Stuck-lock recovery**
>
> Concurrent `bwa-mem3 shm <prefix>` invocations are serialized by a named
> POSIX semaphore (`/bwactl_lock`) so the registry stays consistent. POSIX
> semaphores have no `SEM_UNDO` equivalent: if a stager segfaults or is
> `kill -9`'d while holding the lock, every subsequent stage will block in
> `sem_wait` forever. Run `bwa-mem3 shm -d` to recover — it unlinks the
> semaphore alongside the registry, freeing the next stager.

---

**See also:**
[Getting Started — Quick start: shared-memory index](../getting-started/quick-shm.md) ·
[CLI Reference — index](index-cmd.md) ·
[CLI Reference — mem](mem.md) ·
[Best Practices — Multi-sample workflows](../best-practices/multi-sample.md) ·
[Best Practices — Anti-patterns](../best-practices/anti-patterns.md)
