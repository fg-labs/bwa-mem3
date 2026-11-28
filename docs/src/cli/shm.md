# shm

`bwa-mem3 shm` stages an FM-index into shared memory (POSIX shared memory, or
with `--huge-pages` a System V segment on 1 GB pages) so that subsequent
`bwa-mem3 mem` invocations on the same machine attach to the in-memory segment
instead of re-reading the index files from disk. For workloads that align many
small samples back-to-back against the same reference — such as clinical
panels or amplicon sequencing — this removes the dominant I/O bottleneck.
`shm` also lists and destroys staged segments.

> **Important — on Linux, mount `/dev/shm` with `huge=advise` before staging**
>
> On a stock `/dev/shm` mount (`huge=never`) a staged index sits on 4 KB pages, and
> alignment against it uses more CPU than an index loaded from disk, by a quarter to a
> third on a human genome
> ([measurements below](#huge-pages-mount-devshm-with-hugeadvise-or-hugewithin_size)).
> When that happens, `shm` (after staging) and `mem` (when it attaches) print a one-line
> `[W::bwa_shm_stage]` or `[W::bwa_shm_attach]` warning; exit status and output do not
> change. The warning is skipped where the segment cannot be on huge pages: a segment smaller than one
> huge page (2 MB; 512 MB on an arm64 kernel with 64 KB pages), a kernel whose
> `/sys/kernel/mm/transparent_hugepage/shmem_enabled` is `[deny]`, which keeps every tmpfs
> mount on 4 KB pages, and an index staged with `--huge-pages`; and it is not printed
> with `-v 1` or lower, or when `/proc/self/smaps` cannot be read. To fix it, remount
> and restage (but see `[force]` below). The page size is fixed when the segment is written, so an index staged
> before the remount stays on 4 KB pages until restaged:
>
> ```bash
> sudo mount -o remount,huge=advise /dev/shm
> bwa-mem3 shm -d
> bwa-mem3 shm ref.fa
> ```
>
> To keep it across reboots, add `huge=advise` to the `/dev/shm` options in
> `/etc/fstab`, e.g. `tmpfs /dev/shm tmpfs defaults,huge=advise 0 0`.
>
> Where `shmem_enabled` is `[force]`, every tmpfs mount already gets huge pages whatever
> its `huge=` option, so a segment left on 4 KB pages means the kernel found no free huge
> page when it was written (memory too low or too fragmented). The warning then says so
> and advises freeing memory or compacting it
> (`echo 1 | sudo tee /proc/sys/vm/compact_memory`, with
> `/sys/kernel/mm/transparent_hugepage/defrag` not `[never]`) before restaging, not a
> remount.

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
The exceptions are segments it cannot safely remove: one whose unlink fails
(another user's, for example), one whose header cannot be read, and a
`--huge-pages` segment that cannot be removed or identified (another user's,
another IPC namespace's, or a stage interrupted in another PID namespace). Each
is kept with its entry, and `shm -d` exits non-zero; see
[`--huge-pages`](#--huge-pages--stage-the-index-on-1-gb-huge-pages-linux-opt-in).

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
The exception is the 64-byte stub of a `--huge-pages` stage, which names a
System V segment holding pinned 1 GB pages: unlinking the stub alone leaves
those pages allocated with nothing that names them. For it the stage's message
names that segment's id and the checks to make before removing the segment
with `ipcrm -m <id>`, and only then the stub (see
[`--huge-pages`](#--huge-pages--stage-the-index-on-1-gb-huge-pages-linux-opt-in)).

### `-l` — list staged indices

Prints the name and size in bytes of each registered segment, tab-separated,
one per line. The size is that of the packed index; for a `--huge-pages` stage
the pages it pins are that size rounded up to whole 1 GB pages. A third column gives an entry's
state when it is not a plain complete index: `compact` for one `mem` attaches (below), and
otherwise why `mem` will not attach to it:

- `incomplete`: its stage is still running or was interrupted, or it was staged
  with `--huge-pages` and its 1 GB segment has since been removed;
- `missing`: the entry has no segment (a stage that lost its lock, or a
  `shm -d` that removed the segment without the stager lock and so left the
  registry as it was; the next stage drops the entry);
- `unusable`: the segment cannot be opened, is not a bwa-mem3 index of this
  version, or is shorter than its header declares; or it was staged with
  `--huge-pages` and its 1 GB segment is another program's, belongs to another
  IPC namespace, is not the size its stub records, or does not hold the complete
  index its stage packed (its header, and every section a load reads, of the
  size that load requires).

- `compact`: a complete segment staged with
  [`--compact-occ`](#--compact-occ--stage-the-compact-checkpoint-layout); informational, `mem`
  attaches it.

Useful to confirm that staging succeeded before launching alignment jobs.

After an upgrade that changes the staged layout, an index staged by the older
bwa-mem3 is listed as `unusable`: `mem` loads from disk instead, and restaging
that index fails. It stays staged, holding its memory (or its 1 GB pages), until
`shm -d` removes it; then stage it again.

### `--meth` — stage a methylation index

Stages the `.meth` seed index (`<idxbase>.meth.*`) into shared memory, mirroring
the behavior of `bwa-mem3 index --meth` and `bwa-mem3 mem --meth`. Pass the same
plain `<idxbase>` to all three commands; the `.meth` suffix is handled
transparently.

The staged seed segment is **seed-only**: it holds the seed FM-index and contig
metadata (BNS) but omits the packed reference (PAC), because `mem --meth` extends
against the original reference and never reads the seed's bases. This trims
~1.6 GB from the staged segment on hg38. A load that does read the packed
reference (a plain `mem` given the seed prefix itself) declines such a segment
and loads from disk, saying so. (The unpacked `.0123` reference is never
staged for **any** index — plain or seed — because `mem` pac-fetches reference
bases from `.pac` on demand; that saves the seed's ~13 GB and the original's
~6.4 GB versus staging `.0123`.)

### `--huge-pages` — stage the index on 1 GB huge pages (Linux, opt-in)

Stages the index on **1 GB huge pages** from the host's hugetlb pool, to cut the
data-TLB page walks of seeding below what the 2 MB pages of `/dev/shm` allow
(see [Performance](#1-gb-pages-shm---huge-pages)). `mem`, `shm -l` and `shm -d`
are used exactly as for any other staged index; `mem` needs no flag. A 1 GB stage
holds the default checkpoint layout: it is refused with
[`--compact-occ`](#--compact-occ--stage-the-compact-checkpoint-layout).

`/dev/shm` cannot hold 1 GB pages, so the index goes into a System V shared
memory segment created with `SHM_HUGETLB | SHM_HUGE_1GB`, and `/dev/shm` holds
only a 64-byte stub that names it. The index is rounded up to whole pages: a
default hg38 stage (10,557,461,248 bytes, 9.83 GiB) takes 10, and `-u` adds the pages of
its denser table. An index of any size takes at least one whole page.

Which hosts: it was measured on three CPUs only, on hg38 with 5M read pairs
(`mem -t 16`, clang 19): it was faster on an Intel Sapphire Rapids host
(AVX-512BW tier, −1.7 % wall) and, by less, on Graviton4 (arm64, NEON tier,
−0.7 % wall), and not on AMD Zen 4 (see
[Performance](#1-gb-pages-shm---huge-pages)). Other CPUs and workloads are not
measured. On
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

- Enough free, unreserved 1 GB pages in the host's pool: 10 for a default hg38
  stage. `nr_hugepages` sets the pool's total size; it does not add pages, and
  other workloads on the host may be using or have reserved some of them. So
  check what is available first, and grow the pool only by the shortfall, never
  shrinking it:

  ```bash
  pool=/sys/kernel/mm/hugepages/hugepages-1048576kB
  need=10
  avail=$(( $(cat $pool/free_hugepages) - $(cat $pool/resv_hugepages) ))
  if [ "$avail" -lt "$need" ]; then
    echo $(( $(cat $pool/nr_hugepages) + need - avail )) | sudo tee $pool/nr_hugepages
    avail=$(( $(cat $pool/free_hugepages) - $(cat $pool/resv_hugepages) ))
    [ "$avail" -ge "$need" ] || echo "only $avail of $need 1 GB pages available" >&2
  fi
  ```

  The second check matters: the kernel adds fewer pages than asked when memory
  is fragmented, so recount the unreserved pages (`free_hugepages` minus
  `resv_hugepages`) after growing the pool. Reserving at boot, before memory
  fragments (`hugepagesz=1G hugepages=N` on the kernel command line, with N
  counting the other workloads' pages too), avoids that. The pool is the host's, shared by
  every container on it. A staged index does not survive a reboot: reserve the
  pool again and stage again. Lowering `nr_hugepages` does not take pages from a
  staged index; they return to the pool only after `shm -d`.
- Permission to create a hugetlb segment: `CAP_IPC_LOCK`, or membership of the
  group named by `vm.hugetlb_shm_group`
  (`sudo sysctl vm.hugetlb_shm_group=$(id -g)`; the setting is host-wide and
  lasts until reboot unless it is also put in `/etc/sysctl.d/`. Check the
  current value first with `sysctl -n vm.hugetlb_shm_group`: changing it takes
  the permission from the group it named). A stage run in
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
  namespace, against the `/dev/shm` it was staged in: the segment is in the
  IPC namespace, but the stub and the registry are files in `/dev/shm`, which
  belongs to the mount namespace. Run it inside the container
  (`docker exec <id> bwa-mem3 shm -d`), or
  `sudo nsenter -t <pid> --ipc --mount bwa-mem3 shm -d` with the PID of a
  process there (the container must have `bwa-mem3`); either removes both the
  segment and the stub. `--ipc` without `--mount` is enough only when this
  `/dev/shm` is the one the stage used (as when `shm -d` run here names the
  index's IPC namespace): with a container's private `/dev/shm` it reads the
  host's registry and leaves the container's index staged. For a stage
  interrupted before it recorded its segment's id, enter the stager's PID
  namespace as well (add `--pid`, with a process in all of them): such a
  segment is identified by its creator PID, which only that PID namespace
  can check (see below), and the `shm -d` message names both namespaces.
  The kernel frees the segment without that only once the IPC namespace no
  longer exists: no process is in it and nothing holds it open (an open file
  descriptor or a bind mount of its
  `/proc/<pid>/ns/ipc`). A container's exit does not show that, since another
  container started with `--ipc=container:<id>` keeps it alive, and `lsns` may
  not list a namespace that only an open file holds. Only when you are certain
  the namespace that the `shm -d` message names (`ipc:[N]`) is gone, delete the
  stub (`rm /dev/shm/bwaidx-<name>`) and run `shm -d`. If unsure, leave the
  stub: it costs one page of `/dev/shm`, while deleting it in front of a live
  segment leaves that segment's pages pinned with nothing that names them.

  A stage interrupted before it recorded its segment is identified by the
  segment's creator PID, which the kernel reports per PID namespace, so
  `shm -d` from another PID namespace (a container run with `--ipc=host` but
  its own PID namespace, for one) keeps it and fails, naming the cause. Run
  `shm -d` in the stager's PID namespace, where the creator can be checked. A
  restage from another PID namespace likewise refuses to reclaim, or wait for,
  such a stage, since whether its stager is alive cannot be told from there.
  The end of a PID namespace frees nothing: if the stager's is gone, nothing
  proves which segment on that key, if any, the stage created (a key, a size
  and an attach count can all match another program's segment), so bwa-mem3
  leaves the segment and the stub as they are. Until that segment is gone,
  `shm -d` keeps that index's entry and exits non-zero (it still removes every
  other index), and a restage of that basename fails. The segment goes when its
  IPC namespace ends, which for the host's namespace is at reboot.

When the segment cannot be created (too few free pages, no permission, no 1 GB
page size on this kernel or CPU, a `kernel.shmmax` or `kernel.shmall` limit),
the stage **fails and stages nothing**; the message gives the system error and
what it means here. Rerun without `--huge-pages` to stage in `/dev/shm`.
There is no silent fallback: the stage would not give the page size it was
asked for, and nothing later would say so.

The pages are pinned and unswappable until `shm -d`. An interrupted stage keeps
its pages too, exactly as it keeps an incomplete segment in `/dev/shm`:
`shm -l` lists it as `incomplete`, and a restage or `shm -d` returns the pages
to the pool. A stage that stops on an unreadable index file removes its segment
as it exits, and leaves only the `incomplete` entry. `shm -d` during a running
`mem` returns the pages when that `mem` exits. A stage that completes prints its
segment's id and key. If the stub in `/dev/shm` is deleted some other way (for
example by `shm -d` run with a bwa-mem3 older than this option, or its restage
of an interrupted `--huge-pages` stage, which remove the stub and do not know
about the segment), bwa-mem3 can no longer find the segment, and the id that
the stage printed is the proof of which one it is: remove that id with
`ipcrm -m <id>` once `ipcs -m` still lists it with that key and a size of the
stage's whole 1 GB pages, and `ipcs -m -i <id>` shows no attached process
(nattch 0). The key and size check that the id was not reused by another
segment after the stage's own was removed. Without that record, leave the
segment: a key, a size and an attach count do not identify one, since another
program's can match all three. A segment
that bwa-mem3 did not create is never removed: the kernel records each
segment's creator, and only the segment made by the stage that a stub records
is touched. If a segment cannot be removed (it belongs to another user),
`shm -d` says so, exits non-zero, and keeps that index's stub and entry so that
the owner's `shm -d` still finds it; everything else is removed.

The pool reservation covers the whole segment, but a NUMA memory policy
(`numactl --membind`) or a hugetlb cgroup limit narrower than the pool can
still make a page unavailable when it is first written. The stage writes every
page before it reads the index, so it then fails and stages nothing. On a Linux
older than 5.14, which cannot write the pages ahead, the refused page kills the
stage with SIGBUS while it packs; the stage then removes its segment as it
dies, so no pages stay pinned, and leaves an `incomplete` entry that a restage
or `shm -d` clears. Stage without such a restriction.

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
`RemoveIPC=yes` (in `logind.conf`) removes the segment when the staging user's
last session ends, together with that user's objects in `/dev/shm` (the stub,
and the registry if that user created it), so the index is then listed as
`missing` or not at all; stage it again. It applies to regular users only:
root and other system users are exempt, so a stage run as one of them keeps its
pages after logout, until `shm -d`.

The segment holds the same packed index bytes as a stage in `/dev/shm`:
`shm_stage_test` checks the attached bytes against the packed index on every
platform it runs on. Alignment records were compared in the measured runs
under Performance (hg38, 5M read pairs, `mem -t 16`, clang 19, on Intel
Sapphire Rapids at the AVX-512BW tier, Graviton4 at the NEON tier and AMD Zen 4
at the AVX-512BW tier), where they were byte-identical to the same index
staged in `/dev/shm` in every arm.

### `--compact-occ` — stage the compact checkpoint layout

Stages the FM-index occurrence checkpoints in the 32-byte compact layout that
[`mem --compact-occ`](mem.md#--compact-occ--hold-the-occurrence-checkpoints-in-half-the-memory-opt-in)
uses, converted while staging: about 3 GB less shared memory on hg38. The staged
index wins: `mem` attaches a staged segment in the layout it holds, with or without
`--compact-occ`. The design contract is byte-identical alignment records for the same
reads and run settings: the compact layout is a lossless representation, not a change
in alignment policy. The [`mem --compact-occ` validation scope](mem.md#--compact-occ--hold-the-occurrence-checkpoints-in-half-the-memory-opt-in)
describes the fixture-based correctness evidence, not a guarantee across all inputs,
hosts or SIMD tiers. `shm_compact_occ_test.sh` checks staged alignment parity on
phiX and a generated seed-only fixture; `shm_pack_round_trip_test` validates compact
segment sections and row resolution. These tests do not establish exhaustive
staged-layout validation. A run that requested the other
layout (`--compact-occ` against a default segment, or the reverse
with `BWA3_COMPACT_OCC=0`) prints a `[W::load_index]` line saying the staged layout is
used and how to restage.

An index has one staged segment, in one layout: to switch, `bwa-mem3 shm -d` (which
drops every staged index) and stage again; staging the other layout over a staged
index is an error. Only this flag selects the layout: `shm` does not read
`BWA3_COMPACT_OCC`, which applies to `mem`. The conversion uses
[`-t`](#-t---threads-int--worker-threads) workers. Not combinable with `-u` (any value),
`-l` or `-d`, nor with
[`--huge-pages`](#--huge-pages--stage-the-index-on-1-gb-huge-pages-linux-opt-in): a 1 GB
stage holds the default checkpoint layout only, so `shm --huge-pages --compact-occ` is
refused, staging nothing, before any host check. `BWA3_COMPACT_OCC` does not change
`shm --huge-pages`.

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
memory is real: on hg38, `-u 2` adds ~3.3 GB to the staged segment (4.125 bytes
per added sample: a 4-byte low word and a 1-bit high part).

To make the denser table permanent on disk (so every `mem` and every future
`shm` stage picks it up with no flag, at the cost of a larger index file), use
[`re-sa`](re-sa.md) instead — it is the on-disk counterpart of this flag.

### `-t, --threads INT` — worker threads

Number of worker threads for the `-u` densify pass and the `--compact-occ`
conversion, default `1`. The per-sample
LF-walks are independent and write disjoint slots, so densification scales
near-linearly with thread count (on hg38, Graviton4, arm64/NEON tier, the
densify pass runs ~11 min at `-t 1` and ~40 s at `-t 16`). Meaningful only with
`-u` or `--compact-occ`; a value `> 1` without either warns and is ignored.

## Performance

Denser SA sampling lowers alignment CPU by shortening the LF-walk each SA
resolve performs. Representative figures on hg38, a 5M-read-pair WGS slice
(Graviton4, arm64/NEON tier, 16 threads), measured same-boot and interleaved —
these are relative user-CPU deltas versus a plain stride-8 index, not absolute
times:

| `-u` | stride | user-CPU Δ vs stride-8 | added shared memory (hg38) |
|---|---|---|---|
| _(none)_ | 8 | — | — |
| `-u 2` | 4 | **≈ −2.5 %** | +3.3 GB |
| `-u 1` | 2 | **≈ −4.2 %** | +10 GB |

**`-u 2` is the best speed-per-GB operating point and the recommended setting.**
(Staging without `-u` does not select it — no-flag staging keeps the on-disk
default rate, `u = 3` / stride 8; `-u 2` is an explicit opt-in.)
The second halving (`-u 1`) is not sharply diminishing — it buys roughly another
−1.7 % for +6.7 GB more — so it is worthwhile where the RAM is available (peak
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

On a `huge=never` mount the request is ignored, and `shm` and `mem` say so
with a `[W::]` line on stderr. Staging still saves the index
load, but each run then pays the CPU cost above, which on a large sample can
exceed the load time saved.

### 1 GB pages: `shm --huge-pages`

`/dev/shm` backs a staged index with 2 MB pages at best. `shm --huge-pages`
puts it on 1 GB pages instead (see
[`--huge-pages`](#--huge-pages--stage-the-index-on-1-gb-huge-pages-linux-opt-in)),
which removes nearly every data-TLB page walk on the index. Whether that saves
time depends on the CPU: it did on the Intel and Graviton4 hosts below, and it
did not on the AMD Zen 4 host, which is why the flag is refused on AMD x86. Measured on hg38 with 5M
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
> stage is ~10.6 GB and `/dev/shm` defaults to half of RAM, so a host with
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
> and later stages then wait for it, retrying every 10 ms. Run
> `bwa-mem3 shm -d` to recover: it takes the lock for its walk of the
> registry, so no stage rewrites the registry under it. A lock still held after
> 5 s is taken to be one left by a killed stager: `shm -d` warns and proceeds
> without it, and unlinks the semaphore alongside the registry. Without the
> lock it does not rewrite the registry, so the entries of the segments it
> removed stay, listed as `missing`, until a stage drops them. It keeps the
> stub of a `--huge-pages` stage that may still be creating its segment (its
> stager is alive, or in another PID namespace), and then leaves the semaphore
> in place too, so that no later `shm -d` can take a new one and remove that
> stub; run `shm -d` again once that stage has finished or its stager has
> exited. Until then new stages wait on that lock. If the stager is dead but
> cannot be seen to be (it was in another PID namespace, which is gone), remove
> the lock by hand (`rm /dev/shm/sem.bwactl_lock` on Linux) and run `shm -d`
> again.
>
> A waiting stage, and a waiting `shm -d`, open the lock by name on each try,
> so once the semaphore is unlinked they move to the new one. When `shm -d`
> removes everything, it leaves the semaphore it unlinked taken. A stage that
> was waiting on that semaphore therefore cannot take it while a new stage
> holds the new one.
>
> Run the same release for every process that shares the registry. A stage from
> a release without this fix waits on one semaphore: if it is already waiting
> when `shm -d` unlinks the lock, it never wakes, so stop it and run it again.
> A `shm -d` from such a release releases the semaphore it unlinked, so a stage
> waiting on that semaphore and a stage that starts afterwards can both hold
> the lock for a moment.
>
> Known limitation: the 5 s timeout cannot tell a killed stager from a live
> one that is merely slow to finish its registry update. A `shm -d` that times
> out removes segments without the lock, and can then race that stager: the
> stager may reclaim and restage an entry between `shm -d` reading its header
> and removing it, and a `--huge-pages` stub removed that way leaves its 1 GB
> segment pinned with nothing that names it. Run `shm -d` only when no stage
> is in progress.

---

**See also:**
[Getting Started — Quick start: shared-memory index](../getting-started/quick-shm.md) ·
[CLI Reference — index](index-cmd.md) ·
[CLI Reference — re-sa](re-sa.md) ·
[CLI Reference — mem](mem.md) ·
[Best Practices — Multi-sample workflows](../best-practices/multi-sample.md) ·
[Best Practices — Anti-patterns](../best-practices/anti-patterns.md)
