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

## Large pages for the index (Linux deployment lever)

Seeding is memory-latency-bound: the FM-index and suffix-array lookups are
effectively random accesses over multi-gigabyte arrays, so the data-TLB and the
page-table walk sit on the critical path. By default bwa-mem3 already hints these
arrays for transparent huge pages (`madvise(MADV_HUGEPAGE)` → 2 MB pages), which
helps but does not fully cover the working set.

On **Linux** you can back the index arrays with **explicit 1 GB huge pages** for
a further speedup, with byte-identical alignment output (see
[Measured effect](#measured-effect) below). bwa-mem3's allocator is mimalloc
(above), so the reservation goes through it — either automatically via the
`--huge-pages` flag, or manually via mimalloc's own environment variable. Both
paths reserve through the **shipping mimalloc build**, where the executable
whole-archives libmimalloc so its `malloc`/`free` are the ones the index arrays
actually allocate through — the build `bwa-mem3 version` reports as
`mimalloc … (active)`. A `USE_MIMALLOC=0` build has no mimalloc to reserve
through, and a libmimalloc linked without the malloc override
(`linked but NOT overriding malloc`) would reserve pages the index never
allocates from; the reservation code checks only that mimalloc's reservation
entry point is present, so confirm `(active)` before relying on the win.

### The `--huge-pages` flag (recommended)

Reserve 1 GB hugepages on the host, then pass `--huge-pages`:

```bash
# Reserve 1 GB hugepages (root). Works at runtime on a freshly booted machine
# with enough free contiguous RAM — no reboot needed. hg38 needs ~13 pages;
# rounding up is harmless. 1 GB HugeTLB support is architecture-dependent, so the
# pool directory only exists on hosts that provide it — reserve only when it does
# (the flag itself no-ops safely, but the reservation write would error otherwise).
#
# nr_hugepages is the HOST-WIDE TOTAL, not a delta: `echo 16` sets the whole 1 GB
# pool to 16 pages (it does not add 16). On a shared host, writing an absolute
# value BELOW the current pool shrinks it and can starve other processes already
# using it — read the current count first and pick a non-decreasing target:
#   cur=$(cat /sys/kernel/mm/hugepages/hugepages-1048576kB/nr_hugepages)
#   echo $(( cur > 16 ? cur : 16 )) | sudo tee .../nr_hugepages
if [ -d /sys/kernel/mm/hugepages/hugepages-1048576kB ]; then
    echo 16 | sudo tee /sys/kernel/mm/hugepages/hugepages-1048576kB/nr_hugepages
fi

bwa-mem3 mem --huge-pages ref.fa r1.fq r2.fq > out.sam
```

`--huge-pages` sizes the reservation from the index footprint automatically and
reserves the pages before the index loads. It is **safe by default**: if the host
has no 1 GB pool, or too few free pages, it prints a one-line `[M::]` note and
runs on the default page size — it does not fail the alignment. The reservation
itself makes only a bounded attempt (a finite timeout, at least 5 s and up to
~1 s per requested page) before falling back to default pages, so a heavily
fragmented host may see a brief startup delay, never a failed run:

```text
[M::bwamem_reserve_huge_pages] --huge-pages: reserved 13 x 1 GB huge pages for the index via mimalloc
# or, when unavailable:
[M::bwamem_reserve_huge_pages] --huge-pages: the index needs ~13 free 1 GB pages but only 0 are available; not reserving (raise nr_hugepages). Running on default pages
```

### Manual: `MIMALLOC_RESERVE_HUGE_OS_PAGES`

Equivalently, without the flag, reserve the pages and tell mimalloc directly. `N`
is the 1 GB page count, sized to the resident index (hg38 ≈ 11 GB → 11 pages + a
2-page margin = 13, the same count `--huge-pages` auto-sizes):

```bash
N=13   # hg38 example; size to your reference
# Fail closed: reserve and run only with an active mimalloc override AND a 1 GB
# HugeTLB pool. A USE_MIMALLOC=0 build silently ignores MIMALLOC_RESERVE_HUGE_OS_PAGES.
# As above, nr_hugepages is the HOST-WIDE TOTAL: `echo "$N"` sets the whole pool
# to N (it does not add N). On a shared host, don't write a value below the
# current pool or you shrink it under other processes — use a non-decreasing
# target: cur=$(cat .../nr_hugepages); N=$(( cur > N ? cur : N )).
if bwa-mem3 version | grep -q 'mimalloc.*(active)' \
   && [ -d /sys/kernel/mm/hugepages/hugepages-1048576kB ]; then
    echo "$N" | sudo tee /sys/kernel/mm/hugepages/hugepages-1048576kB/nr_hugepages
    MIMALLOC_RESERVE_HUGE_OS_PAGES="$N" bwa-mem3 mem ref.fa r1.fq r2.fq > out.sam
else
    echo "1 GB huge pages unavailable (need active mimalloc + 1 GB HugeTLB); running on default pages"
fi
```

### Verifying

Confirm the index actually landed on 1 GB pages, rather than silently falling
back to the default page size:

```bash
# while a run is in flight — resolve exactly one 'bwa-mem3 mem' PID first, so a
# stray second match (e.g. this pipeline, another run) can't skew the count.
pid=$(pgrep -f 'bwa-mem3 mem')
[ "$(printf '%s\n' "$pid" | grep -c .)" -eq 1 ] \
  || { echo "expected exactly one 'bwa-mem3 mem' PID, got: ${pid:-none}"; }

# Count index VMAs actually backed by 1 GB HugeTLB pages. KernelPageSize alone
# does not prove HugeTLB backing, so require all three per mapping: HugeTLB flag
# (VmFlags: ht), 1 GB page size (1048576 kB), and nonzero Private_Hugetlb.
awk '
  /^[0-9a-f]+-[0-9a-f]+ /  { kps=0; ph=0 }
  $1 == "KernelPageSize:"  { kps = ($2 == 1048576) }
  $1 == "Private_Hugetlb:" { ph  = ($2 > 0) }
  $1 == "VmFlags:"         { if (kps && ph && $0 ~ /(^| )ht( |$)/) n++ }
  END { print n + 0 }
' /proc/"$pid"/smaps
# a nonzero count (roughly one mapping per large index array) confirms 1 GB backing
```

### Measured effect

Measured on a 5 M-read WGS slice (5 M paired reads — the total workload, not the
batch size — HG00096 aligned to hg38); host AMD Zen3, 16 physical cores / 32
threads; SIMD tier avx2; built with clang-19; index warm in the page cache;
compared against the default path (transparent huge pages, 2 MB) on the same host
and input ([PR #405](https://github.com/fg-labs/bwa-mem3/pull/405); the
reproducible, multi-architecture throughput methodology this follows is
[bwa-mem3-bench](../related-projects/bwa-mem3-bench.md)). Both runs used matching
`-t 32` and no explicit `-K`, i.e. the default per-batch size (`chunk_size` =
10 Mbp per thread), and differed only by `--huge-pages`, so the batch boundaries
were identical:

- **wall time:** ~1.5 % lower
- **user CPU:** ~0.9 % lower
- **system CPU:** ~25 % lower, from the eliminated page-table walks

The **alignment records are byte-identical** to the default path — page size does
not change alignments; only the `@PG` header line differs, since it records the
`--huge-pages` flag on the command line. This is a single-host, single-config
measurement; the effect at other thread counts, architectures, or references is
not characterized here.

> **Note — why this is opt-in, not a default**
>
> 1 GB hugepages must be pre-reserved (privileged, host-global, pinned
> non-swappable RAM), the number needed depends on the reference size, and with
> no pages reserved there is nothing to use. `--huge-pages` therefore stays off
> by default and no-ops when the host is not configured for it. It is also
> **Linux-only** (and needs the bundled mimalloc), and transparent huge pages
> (2 MB) alone do not reproduce the win.
>
> Both paths go through **mimalloc**, never an `LD_PRELOAD` `posix_memalign`
> shim: mimalloc overrides the allocator, so a libc-level preload never sees the
> index allocations and quietly does nothing.

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
