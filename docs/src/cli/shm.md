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

```bash
bwa-mem3 shm ref.fa
bwa-mem3 mem -t 16 ref.fa sample1_R1.fq sample1_R2.fq > sample1.sam
bwa-mem3 mem -t 16 ref.fa sample2_R1.fq sample2_R2.fq > sample2.sam
bwa-mem3 shm -d
```

Stage a methylation index and align:

```bash
bwa-mem3 shm --meth ref.fa
bwa-mem3 mem --meth -t 16 ref.fa R1.fq R2.fq | samtools sort -o out.bam -
bwa-mem3 shm -d
```

List all currently staged segments:

```bash
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

### `-u INT` — densify the staged SA sample table

Stages a **denser** suffix-array (SA) sample table than the one on disk, without
rebuilding the index. SA resolution walks back (LF-mapping) to the nearest stored
sample; one row in `1 << u` is sampled (the on-disk default is `u = 3`, one in 8).
A denser table (smaller `INT`) means fewer LF steps per resolve, so faster
alignment, in exchange for a larger staged segment. For example, `-u 2` stages a
stride-4 table from a stock stride-8 index.

The added samples are synthesised at stage time by the same LF-walk the resolver
uses, so the staged table is **byte-identical** to one built on disk at that rate
with [`index -u`](index-cmd.md#-u-int--sa-sample-rate) — alignments are unchanged.
This byte-identity is a deterministic property of the construction (the LF-walk
recovers exactly the SA value the resolver would compute for each added row), not
a benchmarked measurement, so it holds independent of host, architecture, or SIMD
tier; it is enforced by the shm densify-parity test (`test/shm_pack_round_trip_test.sh`),
which checks the staged table against an independent `index -u` build.
`INT` must be in `[0, 6]`; a value that is not strictly denser than the on-disk
rate is ignored with a warning (the disk rate is staged unchanged). The extra
memory is real: on hg38, `-u 2` adds ~4 GB to the staged segment.

To make the denser table permanent on disk (so every `mem` and every future
`shm` stage picks it up with no flag, at the cost of a larger index file), use
[`re-sa`](re-sa.md) instead — it is the on-disk counterpart of this flag.

### `-t, --threads INT` — densify worker threads

Number of worker threads for the `-u` densify pass, default `1`. The per-sample
LF-walks are independent and write disjoint slots, so densification scales
near-linearly with thread count (on hg38, Graviton4, arm64/NEON tier, the
densify pass runs ~11 min at `-t 1` and ~40 s at `-t 16`). Meaningful only with
`-u`; a value `> 1` without `-u` warns and is ignored.

## Performance

Denser SA sampling lowers alignment CPU by shortening the LF-walk each SA
resolve performs. Representative figures on hg38, a 5M-read-pair WGS slice
(Graviton4, arm64/NEON tier, 16 threads), measured same-boot and interleaved —
these are relative user-CPU deltas versus a plain stride-8 index, not absolute
times:

| `-u` | stride | user-CPU Δ vs stride-8 | added shared memory (hg38) |
|---|---|---|---|
| _(none)_ | 8 | — | — |
| `-u 2` | 4 | **≈ −2.5 %** | +4 GB |
| `-u 1` | 2 | **≈ −4.2 %** | +12 GB |

**`-u 2` is the best speed-per-GB operating point and the recommended setting.**
(Staging without `-u` does not select it — no-flag staging keeps the on-disk
default rate, `u = 3` / stride 8; `-u 2` is an explicit opt-in.)
The second halving (`-u 1`) is not sharply diminishing — it buys roughly another
−1.7 % for +8 GB more — so it is worthwhile where the RAM is available (peak
staging RSS is ~26 GB on hg38), but `-u 2` is the sweet spot. `-u 0` (stride-1)
would need a host with well over 32 GB just to stage hg38 and is not recommended.
The resolved coordinates are byte-identical either way — a deterministic property
of the construction (denser sampling changes only LF-walk step counts, never the
resolved coordinate), independent of host, architecture, or SIMD tier and enforced
by `test/shm_pack_round_trip_test.sh`, not a benchmarked result. The *speed-up*
above is a measurement (see the scope stated with the table) and scales with how
resolve-heavy the workload is and how large the reference is, so treat the
percentages as indicative rather than a guarantee for your data.

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
[CLI Reference — re-sa](re-sa.md) ·
[CLI Reference — mem](mem.md) ·
[Best Practices — Multi-sample workflows](../best-practices/multi-sample.md) ·
[Best Practices — Anti-patterns](../best-practices/anti-patterns.md)
