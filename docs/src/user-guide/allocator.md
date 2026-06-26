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
> `bwa-mem3` binary (or be reachable via the embedded rpath). If you move
> `bwa-mem3` without also moving `libmimalloc.dylib`, the binary will fall back
> to the system allocator silently — `bwa-mem3 version` will not print a mimalloc
> line, which is the indicator that the allocator is active.

## Verifying that mimalloc is active

Run:

```bash
./bwa-mem3 version
```

When mimalloc is linked and loaded, the output includes a line like:

```text
mimalloc 3.x.x
```

If that line is absent, mimalloc is not active.

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
