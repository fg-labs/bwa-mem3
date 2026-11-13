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

Removes every registered bwa-mem3 shared-memory segment on the machine. This is
the correct clean-up command after a batch job and the required step before
re-building the index (see the footgun warning below).

A stage that is killed (Ctrl-C, `kill`, a crash) or that exits on a fatal read
error leaves an incomplete segment: `shm -l` lists it as `incomplete`, and `mem`
loads that index from disk instead of attaching. A stage that fails with an
ordinary error removes its own segment and registry entry. Rerun
`bwa-mem3 shm <idxbase>` to reclaim an incomplete segment: the new stage sees
that the old stager is gone, removes the segment and stages again.

While a stage is running, a second `shm <idxbase>` for the same basename prints
`'<name>' is being staged by PID N; waiting` and waits, holding no lock, then
returns as soon as the first stage finishes (or stages the index itself if the
first stage died or failed). After an hour it gives up with
`gave up after 3600 s waiting for PID N to stage '<name>'; if PID N is not a
running bwa-mem3 shm, run bwa-mem3 shm -d`, for example when PID N was reused by
an unrelated process. `shm -d` during a stage unlinks the segment and that stage
then fails; the stager's memory is freed only when it exits.

A complete segment with no registry entry is never adopted or overwritten: a
stage with the same basename fails and names it. `shm -d` cannot find such an
unregistered segment. Remove it with `shm_unlink(3)` (on Linux, unlink its file
under `/dev/shm`, for example `/dev/shm/bwaidx-ref.fa`), then stage again.

For an index staged with [`--huge-pages`](#--huge-pagesdir--stage-the-index-on-1-gb-huge-pages-linux-opt-in),
`shm -d` also removes its file on the hugetlbfs mount, which returns its 1 GB
pages to the pool, and so does the reclaim of an interrupted 1 GB page stage. If
that file cannot be removed, `shm -d` says so and exits non-zero: its pages stay
reserved until it is removed.

### `-l` — list staged indices

Prints the name and size in bytes of each registered segment, tab-separated,
one per line. A third column flags an entry `mem` will not attach to:

- `incomplete`: its stage is still running or was interrupted;
- `missing`: the entry has no segment (a stage that lost its lock; the next
  stage drops the entry), or it is a `--huge-pages` redirect whose file is gone
  (removed by hand, or its hugetlbfs mount is gone; the next stage reclaims it);
- `unusable`: the segment cannot be opened, is not a bwa-mem3 index of this
  version, or is shorter than its header declares; or it is a `--huge-pages`
  redirect whose file does not check out (see below).

Useful to confirm that staging succeeded before launching alignment jobs.

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

### `--huge-pages[=DIR]` — stage the index on 1 GB huge pages (Linux, opt-in)

Writes the packed index to a file on a `hugetlbfs` mount with 1 GB pages
instead of into `/dev/shm`, so every `mem` that attaches maps it on 1 GB pages.
This is the shared-memory counterpart of
[`mem --huge-pages`](mem.md#--huge-pages--put-the-index-on-1-gb-huge-pages-linux-opt-in),
which has no effect on an attached index: it cuts data-TLB page walks in
seeding on hosts whose 1 GB pages reach the TLB. `/dev/shm` is a tmpfs, which
can back a segment with 2 MB pages at most (see
[Huge pages](#huge-pages-mount-devshm-with-hugeadvise-or-hugewithin_size) below).

With no `DIR`, `shm` uses the first `hugetlbfs` mount with `pagesize=1G` in
`/proc/mounts` that it can create files in. The file is `DIR/bwaidx-<name>`,
and the whole index goes on 1 GB pages, rounded up to whole pages: 11 pages for
hg38, about 4 more with `-u 2`. The pages are reserved when the file is mapped
and faulted in over several threads before the index is written. The registry
entry and the name are the same as for a `/dev/shm` stage, so `mem`, `shm -l`
and `shm -d` need no flag: `/dev/shm` holds a small redirect record that points
at the file, and `shm -d` removes the file too, which returns its pages to the
pool.

Set up once per boot (as root), then stage as the user who will run `mem`:

```bash
echo 12 | sudo tee /sys/kernel/mm/hugepages/hugepages-1048576kB/nr_hugepages
sudo mkdir -p /mnt/huge-1g
sudo mount -t hugetlbfs -o pagesize=1G,uid=$(id -u),gid=$(id -g),mode=0755 none /mnt/huge-1g
bwa-mem3 shm --huge-pages ref.fa        # or --huge-pages=/mnt/huge-1g
bwa-mem3 mem -t 16 ref.fa R1.fq R2.fq > out.sam
bwa-mem3 shm -d                         # removes /mnt/huge-1g/bwaidx-ref.fa too
```

**Safe by default:** when there is no such mount, the pool has too few free
pages, the mount's `size=` limit is reached, the pages cannot be faulted in
(this needs Linux 5.14 or later), or the index is 512 MiB or less, `shm` stages
into `/dev/shm` as usual and a one-line `[M::bwa_shm_stage]` note says why. The
staged bytes are the same either way, so **alignments are byte-identical**.

An index already staged stays where it is: `shm --huge-pages` of an index
already in `/dev/shm` (or `shm` without the flag of one already on 1 GB pages)
returns at once with a note. Run `shm -d` first to restage it the other way.

The redirect record is written before the file is created and removed after
it, so a `--huge-pages` stage that is killed or fails leaves a redirect that
`shm -l` lists as `incomplete`; the next stage of that name (with or without
`--huge-pages`) or `shm -d` removes the file and frees its pages. A redirect
is followed only to `bwaidx-<name>` for its own name, as a regular file (not a
symlink) owned by the redirect's owner, in a directory on a `hugetlbfs` mount
with 1 GB pages, of the size it records and holding the index it records;
otherwise `shm -l` lists it as `unusable`, `mem` loads from disk, and `shm -d`
removes the redirect but leaves the file.

Notes:

- An existing `bwaidx-<name>` on the mount that no registered redirect names
  is never adopted or removed: the stage fails, names it, and gives the `rm`
  to run. It may be left by a stage whose redirect was removed (by an older
  `bwa-mem3`, or because `/dev/shm` was cleared, for example by systemd's
  `RemoveIPC` at logout), or belong to another registry sharing the mount (a
  container with its own `/dev/shm`). Such a file pins its pages until it is
  removed by hand.
- Do not mix with an older `bwa-mem3` on the same host: one that predates this
  option sees a complete redirect as `unusable` (and loads from disk), and its
  `shm -d` or reclaim removes the redirect but not the file.
- The pages stay reserved, pinned and out of the page cache until `shm -d`.
  `--huge-pages` stages are not counted against `/dev/shm` by the
  [capacity preflight](#notes--gotchas) of other stages.
- hugetlbfs pages are not counted in a process's RSS. To check placement, a
  running `mem`'s `/proc/<pid>/smaps` entry for the file shows
  `KernelPageSize: 1048576 kB` and `Shared_Hugetlb` (or `Private_Hugetlb`).
- As with `mem --huge-pages`, the gain depends on the host: it measured
  about 0.7 % faster for a disk-loaded index on Intel Sapphire Rapids, and
  within noise on AMD Zen 4 instances whose hypervisor backs guest memory with
  2 MB pages. Measure on your own host.
- NUMA placement of the pages is not controlled, and `DIR` must be given as
  `--huge-pages=DIR`: `--huge-pages DIR idxbase` is rejected as two positional
  arguments.

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

### Huge pages: mount `/dev/shm` with `huge=advise` or `huge=within_size`

On Linux, an index loaded from disk requests transparent huge pages; actual
backing is best-effort and depends on the kernel's THP policy and available
memory. (`mem --huge-pages` instead asks for 1 GB pages for the occurrence
table, also best-effort; see [`--huge-pages`](mem.md#--huge-pages--put-the-index-on-1-gb-huge-pages-linux-opt-in).)
A staged segment lives in the `/dev/shm` tmpfs, and its page size also depends
on that mount's `huge=` option, which defaults to `never`. In the measurements
below, disk-loaded and huge-page-backed shm indices used 2 MB pages. With 4 KB
shm pages, FM-index and SA lookups missed the TLB more often and alignment used
more CPU than a disk load:

| hg38, 5M read pairs, `mem -t 16`, clang 19 | disk load | shm, 2 MB pages | shm, 4 KB pages (`huge=never`) |
|---|---|---|---|
| Zen 4 (m7a.4xlarge), AVX-512BW: CPU s (wall s) | 497 (32.6) | 495 (32.4) | 661 (43.1), +33 % CPU |
| Graviton4 (m8g.4xlarge), NEON: CPU s (wall s) | 529 (34.5) | 527 (34.3) | 665 (43.2), +26 % CPU |

Output is byte-identical in every column. `bwa-mem3 shm` and `mem` ask for huge
pages on the segment (`madvise(MADV_HUGEPAGE)`), and with either of these mounts
a kernel with transparent huge pages backs it with 2 MB pages, which gave
disk-load speed in the benchmark above. The advice is best-effort: when the
kernel can't supply huge pages the segment stays on 4 KB pages, so confirm the
backing as shown below. Remount, then re-stage, since a segment keeps the page
size it was written with:

```bash
sudo mount -o remount,huge=advise /dev/shm   # or huge=within_size
bwa-mem3 shm -d
bwa-mem3 shm ref.fa
grep ShmemHugePages /proc/meminfo            # host-wide: non-zero once anything is on shmem huge pages
```

`ShmemHugePages` totals every shmem huge page on the host, so it can be non-zero
because of another tmpfs. To check this index's segment, look at its entry in
`/proc/<pid>/smaps` of a running `mem`: `ShmemPmdMapped` is non-zero when the
segment is mapped on 2 MB pages.

On a `huge=never` mount the request is ignored. Staging still saves the index
load, but each run then pays the CPU cost above, which on a large sample can
exceed the load time saved.

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
> After creating the (still empty) index segment, and before sizing or
> filling it, `bwa-mem3 shm` calls `statvfs("/dev/shm")` and compares the
> available bytes, less the part of each other stage still in progress that it
> has not yet written (tmpfs allocates pages only as they are written), against
> the index's `total_size` (a `--huge-pages` stage in progress takes its pages
> from the hugetlb pool, not `/dev/shm`, so it is not counted). Creating the segment
> first means an orphaned segment is reported by name before this check can
> fail. If `/dev/shm` is too small, the stage unlinks the empty segment and
> aborts cleanly with an `[E::bwa_shm_stage]` message that names `/dev/shm`,
> the required size, and a `mount -o remount,size=...` hint. This replaces the
> previous failure mode where `ftruncate` succeeded lazily and `pack_into`
> later surfaced ENOSPC as `[fread] Bad address` with no indication that
> `/dev/shm` was the cause. The preflight is best-effort:
> a `statvfs` failure (no `/dev/shm`, restricted sandbox, ENOSYS) is
> non-fatal and the stage proceeds. As a rough sizing guide, hg38 stages
> ~17 GB; AWS instances default to RAM/2 (so c7a.4xlarge / c7i.4xlarge at
> 32 GB get ~16 GB of `/dev/shm`, which is just under the index size — a
> `remount,size=28g` is the documented fix).
>
> **Note — Staging lock**
>
> Concurrent `bwa-mem3 shm <prefix>` invocations take a named POSIX semaphore
> (`/bwactl_lock`) while they update the registry, which takes milliseconds.
> The index is packed without the lock, so interrupting a stage during the
> pack leaves the lock free. POSIX semaphores have no `SEM_UNDO` equivalent:
> only a stager killed inside the short registry update leaves the lock held,
> and later stages then block in `sem_wait`. Run `bwa-mem3 shm -d` to recover:
> it takes the lock for its walk of the registry, so no stage rewrites the
> registry under it. A lock still held after 5 s is taken to be one left by a
> killed stager: `shm -d` warns and proceeds without it, and unlinks the
> semaphore alongside the registry. Unlinking does not wake a stager already
> blocked in `sem_wait`, so stop any such stager first, then start a new stage.

---

**See also:**
[Getting Started — Quick start: shared-memory index](../getting-started/quick-shm.md) ·
[CLI Reference — index](index-cmd.md) ·
[CLI Reference — re-sa](re-sa.md) ·
[CLI Reference — mem](mem.md) ·
[Best Practices — Multi-sample workflows](../best-practices/multi-sample.md) ·
[Best Practices — Anti-patterns](../best-practices/anti-patterns.md)
