# Build

This page describes the recommended build configuration for production use of bwa-mem3.

## Choose the right arch target

The default `make` invocation builds a single multi-tier binary on x86
(or a single NEON binary on arm64). For production clusters where the
CPU family is uniform, you can trim further by building one tier only —
the binary drops the per-tier dispatch table and ships a single kernel
path:

```bash
# Most modern x86-64 servers (Haswell or later):
make arch=avx2

# Intel Cascade Lake / Sapphire Rapids, AWS c7i/m7i:
make arch=avx512bw

# Apple Silicon / AWS Graviton:
make arch=arm64
```

Omit `arch=` if the deployment target is heterogeneous or unknown; the
default `make` produces a single binary that includes every supported
x86 tier and dispatches at runtime via `__builtin_cpu_supports`. Tune
the non-kernel TU compile baseline with `BASELINE_ARCH=` (default
`avx2`) — see
[Single-binary SIMD dispatch (x86)](../developer-guide/launcher.md).

See [SIMD dispatch matrix](../performance/simd-dispatch.md) for the full list
of targets and which kernels each vectorizes.

## Build with clang, not g++ (strongly recommended)

**clang produces a materially faster bwa-mem3 than g++ on every platform we
measure — build with clang unless you have a specific reason not to.** Pass
`CXX=clang++ CC=clang` to `make` (the default `make` uses g++ and emits a
warning nudging you here). Use the newest clang available; compiler *version*
matters nearly as much as the vendor.

### x86-64

The bigger, and more surprising, win is on x86 — clang's optimizer handles
bwa-mem3's C++ notably better than g++ here. Measured on AWS `c6a.8xlarge`
(AVX2, `-t 16`), hg38, 5M read pairs, `make arch=avx2`, median-of-3:

| Compiler | CPU-seconds | wall | vs g++ |
|---|---:|---:|---:|
| g++ 14.2 | 1486 | 1:39 | — |
| clang 19.1 | 1246 | 1:23 | **~16% faster** |

### ARM / aarch64

The aarch64 build runs its SIMD through the
[sse2neon](../developer-guide/neon-port.md) translation layer rather than
hand-written intrinsics, so codegen quality depends heavily on the compiler
**and its version**. Measured on AWS Graviton4 (`c8g.4xlarge`, 16 cores), hg38,
5M read pairs, `make arm64`, best-of-3 CPU-seconds:

| Compiler | CPU-seconds | vs gcc 15.2 |
|---|---:|---:|
| gcc 15.2 | 1779 | — |
| clang 22.1 | 1679 | ~6% faster |

Takeaways:

- **clang beats g++ on both ISAs** — ~16% on x86 (AVX2) and ~6% on ARM here.
  The x86 gain is the larger one and applies to the mainstream deployment
  target, so clang is the recommended production compiler across the board.
- **Compiler *version* matters too.** On ARM a larger ~18% clang-over-gcc gap
  has been reported against an older gcc (~13); against a modern gcc (15.2) it
  narrows to ~6% because recent gcc closed much of the NEON gap. Prefer clang,
  but if you must use gcc, use the newest one.

If you build with clang, note the OpenMP runtime changes from `libgomp` to
`libomp` (`llvm-openmp`) — see
[Multi-architecture deployment](multi-arch-deployment.md). Confirm which
compiler a given binary was built with via `bwa-mem3 version`, which now prints
a `Compiler:` line.

## Profile-Guided Optimization (PGO)

PGO adds 3–5% throughput on real workloads and is recommended for any
installation that runs many alignment jobs against the same reference. It is
opt-in — the default `make` does not use it. The generate → train → use
workflow, the `PGO_ARCH=` selector, `PGO_PROFILE_DIR=`, and the training-data
caveats are all in [Performance → PGO build](../performance/pgo.md); the Summary
below shows the production recipe.

## mimalloc

mimalloc is compiled in by default (`USE_MIMALLOC=1`). The allocator
improves multi-threaded throughput by reducing lock contention on `malloc`
and `free` hot paths. Run `bwa-mem3 version` to confirm it is active:

```text
bwa-mem3 version
# Expected output includes a line like:
#   mimalloc 3.x.x
```

To build without mimalloc (for example, when using AddressSanitizer or on a
system with a known-incompatible allocator):

```bash
make USE_MIMALLOC=0
```

## Summary

For a production installation on a known x86 server with AVX2 — build with
clang and apply PGO on top:

```bash
make pgo-generate PGO_ARCH=avx2 CXX=clang++ CC=clang
./bwa-mem3.pgo-instr.avx2 mem -t 16 ref.fa R1.fq.gz R2.fq.gz > /dev/null
make pgo-use PGO_ARCH=avx2 CXX=clang++ CC=clang
# Deploy: bwa-mem3.pgo.avx2
```

On ARM/aarch64 (Apple Silicon, AWS Graviton), likewise build with a recent
`clang` and apply PGO on top:

```bash
make pgo-generate PGO_ARCH=arm64 CXX=clang++ CC=clang
./bwa-mem3.pgo-instr mem -t 16 ref.fa R1.fq.gz R2.fq.gz > /dev/null
make pgo-use PGO_ARCH=arm64 CXX=clang++ CC=clang
# Deploy: bwa-mem3.pgo
```

---

**See also:**
[SIMD dispatch matrix](../performance/simd-dispatch.md) ·
[PGO build](../performance/pgo.md) ·
[Memory allocator (mimalloc)](../user-guide/allocator.md) ·
[Building from source](../developer-guide/building.md) ·
[Anti-patterns](anti-patterns.md)
