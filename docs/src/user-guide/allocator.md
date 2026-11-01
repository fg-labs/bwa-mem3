# Memory allocator (mimalloc)

bwa-mem3 vendors and links [mimalloc](https://github.com/microsoft/mimalloc),
Microsoft's high-performance memory allocator, into every binary by default.
On multi-threaded alignment workloads, mimalloc reduces wall-clock time by
replacing the system allocator with one optimized for many small, short-lived
allocations — exactly the access pattern produced by the inner alignment loops.

## What mimalloc replaces

The system allocator (`glibc` malloc on Linux, `libSystem` malloc on macOS) is
a general-purpose allocator with a global lock. Under heavy multi-threaded
allocation pressure — 16+ threads each issuing thousands of short-lived
allocations per batch — the lock becomes a measurable bottleneck. mimalloc uses
per-thread free lists and a segment-based heap to eliminate most of this
contention.

## Platform-specific linkage

The linkage strategy differs by OS:

| Platform | Mechanism |
|----------|-----------|
| Linux | Static linkage with `--whole-archive`. The entire mimalloc static library is embedded into the `bwa-mem3` binary; its `malloc`/`free` symbols take precedence over `glibc`'s at link time. |
| macOS | Dynamic linkage via dyld interposing. `libmimalloc.dylib` is built alongside the binary; dyld's `DYLD_INSERT_LIBRARIES` interposing mechanism replaces `malloc`/`free` at load time. The dylib ships next to the binary. |

> **Warning — macOS: keep libmimalloc.dylib next to the binary**
>
> On macOS, `libmimalloc.dylib` must remain in the same directory as the
> `bwa-mem3` binary (or be reachable via the embedded rpath). The binary carries
> a hard dynamic dependency on `libmimalloc.dylib`; if you move `bwa-mem3`
> without it, dyld can no longer resolve the library and the binary fails to
> launch (`dyld: Library not loaded`) rather than silently running on the system
> allocator. Whenever mimalloc *is* linked and loaded, `bwa-mem3 version` reports
> its status explicitly — see
> [Verifying that mimalloc is active](#verifying-that-mimalloc-is-active) below.

## Verifying that mimalloc is active

Run:

```bash
./bwa-mem3 version
```

The `version` output always carries a mimalloc line whenever the library is
linked, with a status suffix reporting whether it is actually intercepting the
standard allocator:

```text
mimalloc 3.x.x (active)
```

`(active)` means a standard `malloc` allocation was routed to mimalloc — the
allocator is doing its job. If instead you see:

```text
mimalloc 3.x.x (linked but NOT overriding malloc)
```

then a libmimalloc that exports only the `mi_*` API is linked — e.g. a
distro/conda `libmimalloc` built without the malloc override — and every real
`malloc`/`free` is still going to the system allocator. If the mimalloc line is
absent entirely, the binary was built with `USE_MIMALLOC=0`.

The status is probed at runtime by allocating through the standard `malloc` and
asking mimalloc whether the pointer lives in one of its heap regions, so it
reflects the allocator that is genuinely in effect — not merely what was linked.

## Opting out

Pass `USE_MIMALLOC=0` at build time to produce a binary linked against the
system allocator:

```bash
make USE_MIMALLOC=0
```

Reasons to opt out:

- **AddressSanitizer (ASAN) builds.** The Makefile automatically sets
  `USE_MIMALLOC=0` when `ASAN_FLAGS` is detected, because ASAN and mimalloc's
  malloc interposing cannot coexist cleanly.
- **Container environments** where distributing a dylib alongside the binary
  is inconvenient.
- **Reproducibility testing** to isolate whether a behavioral difference is
  allocator-related.

> **Note — Default is on**
>
> `USE_MIMALLOC=1` is the default. Opt-out is not recommended for production
> workloads — mimalloc measurably reduces wall time on multi-threaded runs.

## Page purging (`MIMALLOC_PURGE_DELAY`)

mimalloc returns freed memory to the operating system once it has sat unused
for its *purge delay* (1000 ms by default in mimalloc v3). `bwa-mem3 mem` frees
and reallocates the same large per-batch buffers on every batch, so with the
stock delay those pages are handed back and then page-faulted straight back in
by the next batch.

`bwa-mem3 mem` therefore sets the purge delay to **-1 (never purge)**, and
reports the value in effect on stderr:

```text
[M::main] mimalloc purge delay: -1 ms (-1 = never purge; set MIMALLOC_PURGE_DELAY to override)
```

Setting `MIMALLOC_PURGE_DELAY` yourself always wins, for example to let RSS
shrink between batches on a memory-constrained host. mimalloc's legacy name for
the same option, `MIMALLOC_RESET_DELAY`, is honoured too. mimalloc matches both
names case-insensitively, so `mimalloc_purge_delay` works as well, and gets the
same warning below.

```bash
MIMALLOC_PURGE_DELAY=1000 bwa-mem3 mem ref.fa r1.fq r2.fq > out.sam  # mimalloc's stock behaviour
```

Give the delay as a whole number of milliseconds (or `-1` for never). mimalloc
reads the value itself, and two kinds of input do not mean what they look like:

- **Boolean words and empty values are accepted as 0 or 1.** `off`, `no` and
  `false` mean 0 (purge immediately, the opposite of the default); `on`, `yes`,
  `true` and an empty value mean 1 ms.
- **Anything else that is not an integer is ignored**, for example `1s` or
  `250ms`. `mem` then keeps its -1 rather than mimalloc's own 1000 ms, and warns:

  ```text
  [W::main] ignoring MIMALLOC_PURGE_DELAY='1s': mimalloc could not parse it; give whole milliseconds, or -1 to never purge
  ```

Only `mem` changes the default; `index` and the other subcommands run with
mimalloc's own. A `USE_MIMALLOC=0` build has no purge delay and prints no line,
and a build whose mimalloc is linked but not overriding `malloc` says the delay
has no effect.

Measured on a 5 M-read-pair short-read paired-end WGS slice (HG00096, hg38);
host Graviton4, 32 vCPU; SIMD tier neon; built with clang-19; index warm in the
page cache; 2 reps per arm, alternated, wall-time rep-to-rep spread ≤ 0.14 %.
Never purging against the stock 1000 ms delay:

| Threads | Wall time | CPU time | System CPU | Peak RSS |
|--------:|----------:|---------:|-----------:|---------:|
| 32 | −0.6 % | ≈ flat (−0.2 %) | −28 % | +0.5 GB (15.4 → 15.9 GB) |
| 8  | −0.5 % | ≈ flat (−0.4 %) | −41 % | +0.3 GB (14.0 → 14.3 GB) |

The CPU-time deltas are too small to resolve with 2 reps; the wall-time and
system-CPU deltas are well outside the spread.

On the workload above (same slice, Graviton4, neon tier, clang-19 build), the
SAM output was **byte-identical** between never purging and the stock 1000 ms
delay, and the CI regression (phix-derived reads, x86-64 Linux, avx2 build)
checks the same between never purging and purging immediately. Other hosts and
modes are not measured, but none is expected to differ: the purge delay only
decides when free memory goes back to the operating system, never what is
computed. The cost is the higher peak RSS in the table, since freed pages stay
resident for reuse. That cost has been measured only for the workload above:
x86 hosts, `--meth`, the long-read presets (`-x pacbio`, `-x ont2d`) and
`--hic` have not been measured, so set `MIMALLOC_PURGE_DELAY=1000` if peak RSS
matters more than the wall-time gain there.

## Large pages for the index (Linux deployment lever)

Seeding is memory-latency-bound: the FM-index and suffix-array lookups are
effectively random accesses over multi-gigabyte arrays, so the data-TLB and the
page-table walk sit on the critical path. By default bwa-mem3 already hints these
arrays for transparent huge pages (`madvise(MADV_HUGEPAGE)`), which makes them
eligible for huge pages (2 MB on x86-64 and on arm64 with 4 KB base pages) and
helps but does not fully cover the working set. The kernel still falls back to
base pages when it cannot assemble a huge page, for example on fragmented memory.
The hint covers the page-aligned interior of each array, so it also applies in a
`USE_MIMALLOC=0` build, whose system-allocator arrays do not start on a page
boundary. With the kernel's THP mode set to `madvise`, only those advised arrays
are eligible for huge pages in such a build.

On **Linux** you can put the FM-index occurrence table (`cp_occ`)
on **explicit 1 GB huge pages** with `--huge-pages`, with
byte-identical alignment output (see [Measured effect](#measured-effect) below).
The occurrence table's lookups are where nearly all of seeding's page walks
are; the suffix-array samples stay on the default (THP-advised) allocation,
since on a host where 1 GB pages reach the TLB, mapping them as well removed
few further walks. The table is mapped directly from the host's 1 GB hugetlb pool (`mmap`
with `MAP_HUGETLB`), not through the allocator, so this works in any Linux build,
`USE_MIMALLOC=0` included.

### The `--huge-pages` flag

Reserve 1 GB hugepages on the host, then pass `--huge-pages`:

```bash
# Reserve 1 GB hugepages (root). Works at runtime on a freshly booted machine
# with enough free contiguous RAM — no reboot needed. hg38 needs 6 pages for
# its occurrence table (5.99 GiB, independent of `index -u`). A --meth index,
# whose seed reference holds two converted copies of the genome, has a table
# about twice that size, so it needs about 12: the table's size divided by
# 1 GiB, rounded up. Set `want` below accordingly.
# Reserving more than needed is harmless.
# 1 GB HugeTLB support is architecture-dependent, so the pool directory only
# exists on hosts that provide it — reserve only when it does (the flag itself
# falls back safely, but the reservation write would error otherwise).
#
# nr_hugepages is the HOST-WIDE TOTAL, not a delta: `echo 6` sets the whole 1 GB
# pool to 6 pages (it does not add 6). On a shared host, writing an absolute
# value BELOW the current pool shrinks it and can starve other processes already
# using it — read the current count first and pick a non-decreasing target:
# the block below does exactly that.
pool=/sys/kernel/mm/hugepages/hugepages-1048576kB
if [ -d "$pool" ]; then
    cur=$(cat "$pool/nr_hugepages")
    want=6   # hg38; about 12 for a --meth index
    echo $(( cur > want ? cur : want )) | sudo tee "$pool/nr_hugepages"
fi

bwa-mem3 mem --huge-pages ref.fa r1.fq r2.fq > out.sam
```

When the index loads, `--huge-pages` maps the occurrence table on whole 1 GB
pages, faults the pages in over several threads, and reads the table into them.
The pages come out of the pool when the table is mapped and go back when the
process exits. It is **safe by default**: when the pool is missing or has too few
free pages, a process limit (`RLIMIT_AS`, `vm.max_map_count`) refuses the
mapping, or its pages cannot be faulted in, the table uses the default allocation
and a one-line `[M::]` note says so — it does not fail the alignment. A table of
512 MiB or less (the index of a small genome or an amplicon panel) is not put on
1 GB pages at all: rounding it up to a whole page would more than double it,
pinning pool memory other processes share. A larger table is also rounded up to
whole 1 GB pages, so one just over a page boundary reserves nearly a whole extra
page of the pool (hg38's 5.99 GiB table uses six pages almost fully). The default
allocation only asks for transparent huge pages (`madvise(MADV_HUGEPAGE)`); the
kernel may still back it with 4 KB pages when THP is off or memory is
fragmented. Faulting the pages in uses `MADV_POPULATE_WRITE`, so the 1 GB
mapping needs Linux 5.14 or later; on an older kernel the table falls back:

```text
[M::load_index] --huge-pages: occurrence table (5.99 GiB) on 6 x 1 GB pages
# or, when the pool is too small:
[M::load_index] --huge-pages: could not map the occurrence table (5.99 GiB) on 6 x 1 GB pages (Cannot allocate memory: the host's 1 GB pool is missing or too small, see /sys/kernel/mm/hugepages/hugepages-1048576kB; or a hugetlb cgroup reservation limit, or a process limit such as RLIMIT_AS or vm.max_map_count, was reached); using default pages for it
```

The flag applies to an index loaded from disk. An index attached from shared
memory (`bwa-mem3 shm`) lives in the shared segment, so `--huge-pages` does not
move it, and an `[M::]` note says so.

Pages in the hugetlb pool are not counted in a process's resident set size, so
with `--huge-pages` the peak RSS that `/usr/bin/time` reports drops by about the
size of the occurrence table (about 6 GB on hg38). The memory is still in use; it is
accounted under `HugetlbPages` in `/proc/<pid>/status` instead.

`MIMALLOC_RESERVE_HUGE_OS_PAGES` (mimalloc's own setting) does not help here:
mimalloc serves allocations over 2 GiB, which the occurrence table is on hg38,
from the OS rather than from its reserved pages. Use `--huge-pages`.

### Verifying

Confirm the table landed on 1 GB pages from the `[M::]` line above, or from the
process's mappings while a run is in flight:

```bash
# resolve exactly one 'bwa-mem3 mem' PID first, so a stray second match (e.g.
# this pipeline, another run) can't skew the count.
pid=$(pgrep -f 'bwa-mem3 mem')
[ "$(printf '%s\n' "$pid" | grep -c .)" -eq 1 ] \
  || { echo "expected exactly one 'bwa-mem3 mem' PID, got: ${pid:-none}"; }

# Count mappings actually backed by 1 GB HugeTLB pages. KernelPageSize alone
# does not prove HugeTLB backing, so require all three per mapping: HugeTLB flag
# (VmFlags: ht), 1 GB page size (1048576 kB), and nonzero Private_Hugetlb.
awk '
  /^[0-9a-f]+-[0-9a-f]+ /  { kps=0; ph=0 }
  $1 == "KernelPageSize:"  { kps = ($2 == 1048576) }
  $1 == "Private_Hugetlb:" { ph  = ($2 > 0) }
  $1 == "VmFlags:"         { if (kps && ph && $0 ~ /(^| )ht( |$)/) n++ }
  END { print n + 0 }
' /proc/"$pid"/smaps
# 1 = the occurrence table is on 1 GB pages
```

### Measured effect

Measured on AWS c7i.8xlarge (Intel Sapphire Rapids, 32 vCPUs) and c7a.8xlarge
(AMD Zen 4, EPYC 9R14, 32 vCPUs), Amazon Linux 2023 (kernel 6.1), built with
clang 19.1.7 and LTO, 1 GB pages reserved at runtime. The workload was 5 M pairs
of HG00096 WGS aligned to hg38 at `-t 16` with the index warm in the page cache,
compared with the same binary without the flag in eight rotated repetitions per
SIMD tier (means). Zen 4 was measured on three separate instances:

| host | tier | wall time | user + sys CPU | runs faster with the flag |
|---|---|---|---|---|
| Sapphire Rapids | avx512bw | 29.76 s → 29.56 s (−0.7 %) | 475.4 s → 472.6 s (−0.6 %) | 8 of 8 |
| Sapphire Rapids | avx2 | 31.15 s → 30.92 s (−0.7 %) | 497.0 s → 494.2 s (−0.6 %) | 8 of 8 |
| Zen 4, instance 1 | avx512bw | 24.61 s → 24.88 s (+1.1 %) | 395.0 s → 399.8 s (+1.2 %) | 1 of 8 |
| Zen 4, instance 1 | avx2 | 25.13 s → 25.29 s (+0.6 %) | 403.2 s → 406.2 s (+0.7 %) | 3 of 8 |
| Zen 4, instance 2 | avx512bw | 24.52 s → 24.55 s (+0.1 %) | 393.6 s → 394.4 s (+0.2 %) | 3 of 8 |
| Zen 4, instance 2 | avx2 | 25.05 s → 25.00 s (−0.2 %) | 402.0 s → 401.6 s (−0.1 %) | 5 of 8 |
| Zen 4, instance 3 | avx512bw | 24.75 s → 24.71 s (−0.2 %) | 397.2 s → 397.0 s (−0.0 %) | 4 of 8 |
| Zen 4, instance 3 | avx2 | 25.15 s → 25.21 s (+0.2 %) | 403.5 s → 404.9 s (+0.3 %) | 4 of 8 |

System CPU rises by about 0.7 s on both hosts, the cost of zeroing the six 1 GB
pages when they are faulted in. On Sapphire Rapids user CPU falls by about
3.5 s: on the same workload, the cycles with a data-TLB page walk active fell from
7.4 % to 0.8 % with the occurrence table on 1 GB pages. The CPU saving is much
smaller than that share because most of the walk latency overlapped other work.
On Zen 4 the data-TLB counters did not change with the flag: the hypervisor
backs guest memory with 2 MB pages there, and under nested paging a TLB entry is
no larger than the smaller of the guest and host pages, so the guest's 1 GB
mappings give no extra TLB reach. A random-access microbenchmark over 6 GiB
confirmed it, running at the same speed on 1 GB and 2 MB pages. The flag
therefore only adds the zeroing cost there, and the wall-time differences on
instances 2 and 3 are within run-to-run noise: two builds with the same hot code
differed by up to 0.5 % at eight repetitions on these hosts. Instance 1 measured
0.6–1.1 % slower, with user CPU 2 to 4 s higher; nothing in the code or the
counters accounts for that, and it did not reproduce on instances 2 and 3.
The effect on other hosts, thread counts and references is not characterized
here. That includes arm64 Linux, where the flag has not been measured at all
(1 GB huge pages there depend on the kernel's base page size), and multi-socket
hosts: the
pages come from the pool of whichever NUMA node faults them in, with no
interleaving, and placement on a host with more than one node is unmeasured.
Measure on your own hosts before relying on it.

The **alignment records are byte-identical** to the default path — page size does
not change alignments; only the `@PG` header line differs, since it records the
`--huge-pages` flag on the command line.

> **Note — why this is opt-in, not a default**
>
> 1 GB hugepages must be pre-reserved (privileged, host-global, pinned
> non-swappable RAM), the number needed depends on the reference size, and with
> no pages reserved there is nothing to use. `--huge-pages` therefore stays off
> by default and falls back when the host is not configured for it. It is also
> **Linux-only**, and transparent huge pages (2 MB) alone are not equivalent to
> it.

## Build internals

The mimalloc source lives in `ext/mimalloc/` as a git submodule. The Makefile
target builds it via CMake before linking `bwa-mem3`. The relevant Makefile
variables are `MIMALLOC_SRC`, `MIMALLOC_BUILD`, and `MIMALLOC_LIB`.

The feature was introduced in bwa-mem3 as part of the performance improvement
work. See [Features](../whats-different/features.md) and
[Build & infrastructure](../whats-different/build-infra.md) for the PR history.

---

**See also:**
[Threading and resource use](threading.md) ·
[Features: mimalloc](../whats-different/features.md) ·
[Getting Started: installation](../getting-started/installation.md) ·
[Developer Guide: building from source](../developer-guide/building.md) ·
[Best Practices: optimization checklist](../best-practices/optimization-checklist.md)
