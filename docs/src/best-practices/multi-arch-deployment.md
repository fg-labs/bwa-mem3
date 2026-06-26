# Multi-architecture deployment

This page covers running bwa-mem3 in heterogeneous compute environments — AWS Batch with mixed instance families, GCP Batch with mixed CPU platforms, on-prem Slurm with mixed nodes, Kubernetes clusters with mixed node pools.

## Within x86_64: one binary, dynamic dispatch

bwa-mem3 ships a single x86_64 binary that contains five SIMD kernel tiers (`sse41`, `sse42`, `avx`, `avx2`, `avx512bw`) and selects the best one at runtime via `__builtin_cpu_supports`. See `src/simd_dispatch.cpp` for the dispatcher and `src/kernel_dispatch.h` for the per-tier symbol mangling.

Build once at the `BASELINE_ARCH` floor that matches your fleet's oldest x86 host. The default `BASELINE_ARCH=avx2` covers Intel Haswell (2013) and AMD Zen (2017) onward. Within that floor, every host transparently uses its best available tier for the hot kernel paths.

## Across x86_64 and arm64

A single ELF binary cannot span CPU families. You must build two binaries — one for x86_64, one for arm64 — and package them so the right one runs on each host.

The recommended approach is a **Docker manifest-list container of your own making**, with one layer per architecture under a single tag. Example:

```dockerfile
FROM ubuntu:24.04 AS build
RUN apt-get update && apt-get install -y \
    build-essential git cmake pkg-config \
    autoconf automake autoconf-archive libtool \
    zlib1g-dev libdeflate-dev
WORKDIR /src
RUN git clone --recursive https://github.com/fg-labs/bwa-mem3 .
RUN make -j

FROM ubuntu:24.04
COPY --from=build /src/bwa-mem3 /usr/local/bin/bwa-mem3
RUN apt-get update && apt-get install -y libgomp1 zlib1g libdeflate0 \
 && rm -rf /var/lib/apt/lists/*
ENTRYPOINT ["bwa-mem3"]
```

Build for both architectures with one command:

```bash
docker buildx build --platform linux/amd64,linux/arm64 \
  -t <registry>/<image>:<tag> --push .
```

> **Building the arm64 layer with clang?** A recent `clang` produces faster
> NEON code than `gcc` on aarch64 (see
> [Best Practices → Build](build.md#use-a-recent-compiler-especially-on-arm)).
> If you switch the arm64 build to clang, swap the OpenMP runtime: clang links
> `libomp` (LLVM) rather than `libgomp`, so the runtime stage needs `libomp5`
> (or `llvm-openmp`) in place of `libgomp1`.

AWS Batch, GCP Batch, Kubernetes, and containerd all read the manifest list and pull the correct layer based on the host's architecture. The submitter references one tag; the runtime picks the right binary automatically.

## Verifying at runtime

`bwa-mem3 version` reports the build's floor, the kernels compiled in, and the resolved runtime tier. Use this in CI or in your Batch job's startup script to confirm the right layer was pulled:

```text
$ bwa-mem3 version
v0.2.0-12-gabcdef1
SIMD floor: avx2 (x86-64-v3, Haswell 2013+); kernels: sse41 sse42 avx avx2 avx512bw
SIMD runtime: avx512bw (BWAMEM3_FORCE_TIER unset)
```

Grep for `SIMD runtime:` to record the tier each job ran at — useful for post-mortem diagnosis of perf regressions.

## Pre-Haswell hosts

If your fleet really must include pre-Haswell x86 (c4, m4, pre-Skylake Xeons), rebuild with a lower floor:

```bash
make BASELINE_ARCH=sse41
```

Expect roughly 10-15% slower wall time on AVX2 hosts in the same container compared to a default `BASELINE_ARCH=avx2` build. This is the trade-off for broader host coverage; only do it if you actually need pre-Haswell support.

The default `BASELINE_ARCH=avx2` covers virtually every modern compute environment. AWS, GCP, and Azure all default to Haswell-or-newer instance types in current-generation compute environments.

## What the host-floor precheck does

If a job is scheduled onto a host that doesn't meet the build's floor (e.g. an `avx2`-baseline binary lands on a pre-Haswell host), `bwa-mem3 mem` refuses to run with exit code 2 and a clear stderr message:

```text
[E::bwamem3] this binary was compiled for SIMD floor avx2 and emits avx2
instructions in non-kernel translation units. The host CPU does not support
avx2 (detected: sse42). Running would SIGILL on the first avx2 instruction.

To run on this host, rebuild bwa-mem3 with BASELINE_ARCH=sse42 (or lower),
or use a binary built for a lower SIMD floor.
```

This is a clean failure: the job exits before any billable alignment work starts. Compare to the alternative without the precheck (SIGILL deep inside an alignment job, opaque process death, wasted compute).

**Defence-in-depth recommendation:** configure your AWS Batch compute environment (or equivalent) to exclude instance families older than your binary's floor. The precheck protects against accidental scheduling; an allowlist at the orchestrator level prevents the scheduling decision in the first place.

---

**See also:**
[Getting Started → Host requirements](../getting-started/host-requirements.md) ·
[Getting Started → Installation](../getting-started/installation.md) ·
[Build](build.md) ·
[Performance → SIMD dispatch matrix](../performance/simd-dispatch.md)
