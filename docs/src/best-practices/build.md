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

## Profile-Guided Optimization (PGO)

PGO typically yields 3–5% throughput improvement on real workloads.  It is opt-in — the standard
`make` target does not use it — but is recommended for any installation that
will run many alignment jobs against the same reference.

The workflow is three steps:

```bash
# Step 1: Build an instrumented binary (produces bwa-mem3.pgo-instr).
make pgo-generate

# Step 2: Run a representative training workload.
#   Use reads and a reference that reflect actual production input.
#   About 10–30 million read pairs is sufficient.
./bwa-mem3.pgo-instr mem -t 16 ref.fa R1.fq.gz R2.fq.gz > /dev/null

# Step 3: Build the PGO-optimized binary (produces bwa-mem3.pgo).
make pgo-use
```

To target a specific SIMD level, pass `PGO_ARCH=`:

```bash
make pgo-generate PGO_ARCH=avx2
./bwa-mem3.pgo-instr.avx2 mem -t 16 ref.fa R1.fq.gz R2.fq.gz > /dev/null
make pgo-use PGO_ARCH=avx2
# Produces: bwa-mem3.pgo.avx2
```

Profile data is written to `pgo_profiles/` by default. Pass
`PGO_PROFILE_DIR=<path>` to change the location.

> **Tip — Training data matters**
>
> The training workload should resemble production input in read length, base
> quality distribution, and reference composition. A read set that is too short,
> too long, or too easy (low mismatch rate) will bias the branch predictions and
> may produce a build that is slower than the non-PGO baseline on real data.

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

For a production installation on a known x86 server with AVX2:

```bash
make pgo-generate PGO_ARCH=avx2
./bwa-mem3.pgo-instr.avx2 mem -t 16 ref.fa R1.fq.gz R2.fq.gz > /dev/null
make pgo-use PGO_ARCH=avx2
# Deploy: bwa-mem3.pgo.avx2
```

---

**See also:**
[SIMD dispatch matrix](../performance/simd-dispatch.md) ·
[PGO build](../performance/pgo.md) ·
[Memory allocator (mimalloc)](../user-guide/allocator.md) ·
[Building from source](../developer-guide/building.md) ·
[Anti-patterns](anti-patterns.md)
