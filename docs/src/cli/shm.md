# shm

`bwa-mem3 shm` stages an FM-index into shared memory (POSIX shared memory, or
with `--huge-pages` a System V segment on 1 GB pages) so that subsequent
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
The one exception is a `--huge-pages` segment that cannot be removed (another
user's, or another IPC namespace's): it and its entry are kept and `shm -d`
exits non-zero; see [`--huge-pages`](#--huge-pages--stage-the-index-on-1-gb-huge-pages-linux-opt-in).

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

### `-l` — list staged indices

Prints the name and size in bytes of each registered segment, tab-separated,
one per line. The size is that of the packed index; for a `--huge-pages` stage
the pages it pins are that size rounded up to whole 1 GB pages. A third column flags an entry `mem` will not attach to:

- `incomplete`: its stage is still running or was interrupted, or it was staged
  with `--huge-pages` and its 1 GB segment has since been removed;
- `missing`: the entry has no segment (a stage that lost its lock; the next
  stage drops the entry);
- `unusable`: the segment cannot be opened, is not a bwa-mem3 index of this
  version, or is shorter than its header declares; or it was staged with
  `--huge-pages` and its 1 GB segment is another program's or belongs to another
  IPC namespace.

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

### `--huge-pages` — stage the index on 1 GB huge pages (Linux, opt-in)

Stages the index on **1 GB huge pages** from the host's hugetlb pool, to cut the
data-TLB page walks of seeding below what the 2 MB pages of `/dev/shm` allow
(see [Performance](#1-gb-pages-shm---huge-pages)). `mem`, `shm -l` and `shm -d`
are used exactly as for any other staged index; `mem` needs no flag.

`/dev/shm` cannot hold 1 GB pages, so the index goes into a System V shared
memory segment created with `SHM_HUGETLB | SHM_HUGE_1GB`, and `/dev/shm` holds
only a 64-byte stub that names it. The index is rounded up to whole pages: a
default hg38 stage (11,261,255,872 bytes) takes 11, and `-u` adds the pages of
its denser table. An index of any size takes at least one whole page.

Which hosts: it is a gain on Intel x86 and, by less, on Graviton4 (arm64);
those and AMD Zen 4 are the CPUs it was measured on. On
an AMD x86 host the stage is **refused and stages nothing**, because 1 GB
staging measured no gain there (AMD Zen 4 ran 0.4 to 0.6 % slower than with the
index in `/dev/shm` on 2 MB pages; see the table under Performance). The
message says so. To measure another AMD CPU anyway, set
`BWAMEM3_SHM_HUGE_PAGES_AMD=1` in the environment of `bwa-mem3 shm`. The stage
is refused in the same way on a host with `kernel.shm_rmid_forced=1`, where the
kernel removes a segment as soon as no process is attached to it, which for a
stage is when the staging command exits, and where that setting cannot be read
(`/proc/sys` masked, as in some containers), since it may then be 1.

Requirements:

- Enough free, unreserved 1 GB pages in the host's pool: 11 for a default hg38
  stage. `nr_hugepages` sets the pool's total size; it does not add pages, and
  other workloads on the host may be using or have reserved some of them. So
  check what is available first, and grow the pool only by the shortfall, never
  shrinking it:

  ```bash
  pool=/sys/kernel/mm/hugepages/hugepages-1048576kB
  need=11
  avail=$(( $(cat $pool/free_hugepages) - $(cat $pool/resv_hugepages) ))
  if [ "$avail" -lt "$need" ]; then
    echo $(( $(cat $pool/nr_hugepages) + need - avail )) | sudo tee $pool/nr_hugepages
  fi
  ```

  Read `free_hugepages` again afterwards: the kernel adds fewer pages than asked
  when memory is fragmented. Reserving at boot, before memory fragments
  (`hugepagesz=1G hugepages=N` on the kernel command line, with N counting the
  other workloads' pages too), avoids that. The pool is the host's, shared by
  every container on it. A staged index does not survive a reboot: reserve the
  pool again and stage again. Lowering `nr_hugepages` does not take pages from a
  staged index; they return to the pool only after `shm -d`.
- Permission to create a hugetlb segment: `CAP_IPC_LOCK`, or membership of the
  group named by `vm.hugetlb_shm_group`
  (`sudo sysctl vm.hugetlb_shm_group=$(id -g)`; the setting is host-wide and
  lasts until reboot unless it is also put in `/etc/sysctl.d/`). A stage run in
  a container needs one of the two there. Attaching with `mem` needs neither;
  the segment is world-readable, like a segment in `/dev/shm`.
- One IPC namespace for the stage, every `mem` and `shm -d`. System V segments,
  unlike `/dev/shm`, cannot be bind-mounted: a container sees the index only if
  it shares the host's IPC namespace (`docker run --ipc=host`, or
  `--ipc=container:<id>`), not just `/dev/shm`. From another namespace a
  completed stage is listed as `unusable` (an interrupted one as `incomplete`),
  `mem` loads from disk, and `shm -d` and a restage leave it alone and fail,
  naming the cause. The stage itself needs `/proc` mounted, to record its IPC
  and PID namespaces.

  To remove an index staged in another IPC namespace, run `shm -d` in that
  namespace. Only once that IPC namespace no longer exists has the kernel
  removed the segment with it, leaving only the stub; then delete the stub
  (`rm /dev/shm/bwaidx-<name>`) and run `shm -d`. A container's exit does not
  show that: the namespace lives on while any process is in it (another
  container started with `--ipc=container:<id>`, for one) or anything holds it
  open, and its segment then stays allocated. Check first that
  `sudo lsns -t ipc` no longer lists the namespace that the `shm -d` message
  names (`ipc:[N]`), and that nothing bind-mounts it; deleting the stub of a
  live segment leaves its pages pinned with nothing that names them.

  A stage interrupted before it recorded its segment is identified by the
  segment's creator PID, which the kernel reports per PID namespace, so
  `shm -d` from another PID namespace (a container run with `--ipc=host` but
  its own PID namespace, for one) keeps it and fails, naming the cause. Run
  `shm -d` in the stager's PID namespace. The end of a PID namespace frees
  nothing, so if that one is gone the segment is still allocated: find it with
  `ipcs -m` by the key the message gives, and only if it lists the size the
  message gives and `ipcs -m -i <id>` shows no attached process (nattch 0),
  remove it with `ipcrm -m <id>`, then delete the stub and run `shm -d`. (The
  creator PID that `ipcs` shows is the one seen from your PID namespace, not the
  stager's, so it is no check here.)

When the segment cannot be created (too few free pages, no permission, no 1 GB
page size on this kernel or CPU, a `kernel.shmmax` or `kernel.shmall` limit),
the stage **fails and stages nothing**; the message gives the system error and
what it means here. Rerun without `--huge-pages` to stage in `/dev/shm`.
There is no silent fallback, because a staged index on the wrong page size is
slower on every later run and nothing else would say so.

The pages are pinned and unswappable until `shm -d`. An interrupted stage keeps
its pages too, exactly as it keeps an incomplete segment in `/dev/shm`:
`shm -l` lists it as `incomplete`, and a restage or `shm -d` returns the pages
to the pool. A stage that stops on an unreadable index file removes its segment
as it exits, and leaves only the `incomplete` entry. `shm -d` during a running
`mem` returns the pages when that `mem` exits. A stage that completes prints its
segment's key; `ipcs -m` shows every segment with its key, owner and size, and
`ipcrm -m <id>` removes one by hand if the stub in `/dev/shm` was deleted some
other way. That includes `shm -d` run with a bwa-mem3 older than this option,
and its restage of an interrupted `--huge-pages` stage: it removes the stub and
does not know about the segment. Before removing a segment by hand, check that
it is that stage's: `ipcs -m` lists the key the stage printed and a size of
whole 1 GB pages, and `ipcs -m -i <id>` shows no attached process (nattch 0).
A key alone does not identify a segment, since another program's can hold it. A segment
that bwa-mem3 did not create is never removed: the kernel records each
segment's creator, and only the segment made by the stage that a stub records
is touched. If a segment cannot be removed (it belongs to another user),
`shm -d` says so, exits non-zero, and keeps that index's stub and entry so that
the owner's `shm -d` still finds it; everything else is removed.

The pool reservation covers the whole segment, but a NUMA memory policy
(`numactl --membind`) or a hugetlb cgroup limit narrower than the pool can
still make a page unavailable when it is first written. The stage writes every
page before it reads the index, so it then fails and stages nothing (on a Linux
older than 5.14 the kernel kills it with SIGBUS instead, leaving an
`incomplete` entry). Stage without such a restriction.

An index that is already staged is left as it is, whatever kind of pages it is
on. A `--huge-pages` stage of an index that is staged in `/dev/shm` fails, since
it did not get what it asked for; a plain stage of an index that is on 1 GB
pages succeeds and says where it is. To change the kind of pages, run `shm -d`
(which removes every staged index) and stage again.

Two host settings remove System V segments on their own.
`kernel.shm_rmid_forced=1` removes a segment as soon as no process is attached,
so the stage is refused there (it does not affect `/dev/shm`). If it is turned
on later, or a segment is removed with `ipcrm`, `mem` loads from disk, `shm -l`
lists the index as `incomplete`, and a restage reclaims the entry. systemd's
`RemoveIPC=yes` removes the segment when the staging user's last session ends,
together with that user's objects in `/dev/shm` (the stub, and the registry if
that user created it), so the index is then listed as `missing` or not at all;
stage it again.

Alignment records are byte-identical to every other way of loading the index:
the segment holds the same packed index bytes as a stage in `/dev/shm`. That is
a property of the format, independent of host, architecture and SIMD tier, and
`shm_stage_test` checks the attached bytes against the packed index.

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

### Huge pages: mount `/dev/shm` with `huge=advise` or `huge=within_size`

On Linux, an index loaded from disk requests transparent huge pages; actual
backing is best-effort and depends on the kernel's THP policy and available
memory. (`mem --huge-pages` instead asks for 1 GB pages for the occurrence
table, also best-effort; see [`--huge-pages`](mem.md#--huge-pages--put-the-index-on-1-gb-huge-pages-linux-opt-in).
`shm --huge-pages` stages the whole index on 1 GB pages, outside `/dev/shm`; see
[below](#1-gb-pages-shm---huge-pages).)
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

### 1 GB pages: `shm --huge-pages`

`/dev/shm` backs a staged index with 2 MB pages at best. `shm --huge-pages`
puts it on 1 GB pages instead (see
[`--huge-pages`](#--huge-pages--stage-the-index-on-1-gb-huge-pages-linux-opt-in)),
which removes nearly every data-TLB page walk on the index. Whether that saves
time depends on the CPU: it does on Intel and Graviton4, and it does not on AMD
Zen 4, which is why the flag is refused on AMD x86. Measured on hg38 with 5M
read pairs of 150 bp WGS, `mem -t 16`, clang 19, Linux 6.1; five runs of each
arm, interleaved on one host, against the same index staged in `/dev/shm` with
every page of the segment on 2 MB pages (checked in `smaps`):

| `shm --huge-pages` vs shm on 2 MB pages | wall | CPU (user + sys) | user-mode cycles | page walks on the index |
|---|---|---|---|---|
| Intel Sapphire Rapids (Xeon 8488C), AVX-512BW | **−1.7 %** | −1.3 % | −1.4 % | cycles spent in a walk: 7.2 % of all cycles → 0.02 % |
| Graviton4, NEON | −0.7 % | −0.5 % | −0.1 % | walks: 1.57 G → 1.1 M, with no change in cycles |
| AMD Zen 4 (EPYC 9R14), AVX-512BW, two stagings | +0.6 %, +0.4 % | +0.5 %, +0.5 % | +0.3 %, +0.6 % | not measured |

The figures are medians of the per-run ratios. On the Intel host every 1 GB run
was faster than its 2 MB partner, and so was every run on Graviton4, where the
walks go away but cost little to begin with. The Zen 4 host gained nothing:
both stagings measured slightly slower, within that host's run-to-run spread
(about 1 %). Output is byte-identical in every arm.

1 GB pages also do not depend on the kernel finding free 2 MB pages at staging
time. On a host whose memory is fragmented or full of page cache, a segment in
`/dev/shm` can end up partly or mostly on 4 KB pages (the cost is in the table
[above](#huge-pages-mount-devshm-with-hugeadvise-or-hugewithin_size)), while a
1 GB stage either gets its reserved pages or fails.

`mem --huge-pages` is the counterpart for an index loaded from disk: it puts
the occurrence table on 1 GB pages, which removes most of the same page walks
(its measured effect is in
[Allocator → Measured effect](../user-guide/allocator.md#measured-effect)), but
every run still pays for loading the index (about 0.5 s on the Intel host
above).

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
> the index's `total_size`. Creating the segment
> first means an orphaned segment is reported by name before this check can
> fail. If `/dev/shm` is too small, the stage unlinks the empty segment and
> aborts cleanly with an `[E::bwa_shm_stage]` message that names `/dev/shm`,
> the required size, and a `mount -o remount,size=...` hint. This replaces the
> previous failure mode where `ftruncate` succeeded lazily and `pack_into`
> later surfaced ENOSPC as `[fread] Bad address` with no indication that
> `/dev/shm` was the cause. The preflight is best-effort:
> a `statvfs` failure (no `/dev/shm`, restricted sandbox, ENOSYS) is
> non-fatal and the stage proceeds. A `--huge-pages` stage needs room in
> `/dev/shm` for its 64-byte stub only. As a rough sizing guide, a default hg38
> stage is ~11.3 GB and `/dev/shm` defaults to half of RAM, so a host with
> 16 GB of RAM (8 GB of `/dev/shm`) needs a remount (`remount,size=12g`) and one
> with 32 GB does not.
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
