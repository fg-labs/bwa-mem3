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
a further speedup; the output stays byte-identical (see
[Measured effect](#measured-effect) below). Because the allocator is mimalloc
(above), the lever is one of mimalloc's own options — no rebuild and no code
change. It requires an **active mimalloc build** (`bwa-mem3 version` reports
`mimalloc … (active)`); a `USE_MIMALLOC=0` build ignores the reservation:

```bash
# N = number of 1 GB pages; must cover the resident index. hg38 (≈ 11 GB) → 14.
# Size this to *your* reference; the value below is the hg38 example.
N=14

# Preflight (fail closed): the lever works only when BOTH hold, so gate the
# reservation AND the run on both, and do nothing otherwise:
#   * bwa-mem3 must report an ACTIVE mimalloc — a USE_MIMALLOC=0 build (no
#     mimalloc line) or a linked-but-not-overriding build silently ignores
#     MIMALLOC_RESERVE_HUGE_OS_PAGES;
#   * the 1 GB HugeTLB pool must exist (architecture- and kernel-dependent).
if bwa-mem3 version | grep -q 'mimalloc.*(active)' \
   && [ -d /sys/kernel/mm/hugepages/hugepages-1048576kB ]; then
    # Reserve N x 1 GB hugepages on the host (root). Works at runtime on a freshly
    # booted machine with enough free contiguous RAM — no reboot needed.
    echo "$N" | sudo tee /sys/kernel/mm/hugepages/hugepages-1048576kB/nr_hugepages
    # Run with mimalloc reserving them at startup.
    MIMALLOC_RESERVE_HUGE_OS_PAGES="$N" bwa-mem3 mem ref.fa r1.fq r2.fq > out.sam
else
    echo "1 GB huge pages unavailable (need active mimalloc + 1 GB HugeTLB); running on default pages"
fi
```

Verify the index actually landed on 1 GB pages, rather than silently falling back
to the default page size:

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

Measured on a 5 M-read WGS slice (HG00096, hg38), AMD Zen3 (AVX2 tier), 32
threads, with the index warm in the page cache and everything else held fixed —
compared against the default path (transparent huge pages, 2 MB) on the same
host and input:

- **wall time:** ~1.5 % lower
- **user CPU:** ~0.9 % lower
- **system CPU:** ~25 % lower, from the eliminated page-table walks

Output was **byte-identical** to the default path — page size does not change
alignments, so this is purely a performance lever, not an output change. This is
a single-configuration measurement; the effect at other thread counts, hosts, or
references is not characterized here.

> **Note — why this is a deployment lever, not a default**
>
> 1 GB hugepages must be pre-reserved (privileged, host-global, pinned
> non-swappable RAM), the right `N` depends on the reference size, and with no
> pages reserved the reservation fails and falls back. It is therefore opt-in
> host configuration, not something the binary enables by default. It is also
> **Linux-only**, and transparent huge pages (2 MB) alone do not reproduce it.
>
> The reservation must go through **mimalloc**
> (`MIMALLOC_RESERVE_HUGE_OS_PAGES`), not an `LD_PRELOAD` `posix_memalign` shim:
> mimalloc overrides the allocator, so a libc-level preload never sees the index
> allocations and quietly does nothing. A future opt-in `--huge-pages` flag that
> auto-sizes `N` from the index and no-ops when pages are unavailable is tracked
> in [issue #402](https://github.com/fg-labs/bwa-mem3/issues/402).

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
