# Host requirements

bwa-mem3 runs on the hosts in the table below. Verify your host with `bwa-mem3 version` — the SIMD floor and runtime lines tell you what the binary needs and what your host provides.

| Platform | Default build floor | Earliest supported CPU | Notes |
|---|---|---|---|
| Linux x86_64 | AVX2 (`BASELINE_ARCH=avx2`) | Intel Haswell (2013); AMD Zen / Naples (2017) | AVX2 is the minimum floor; kernels self-dispatch to the best of `avx2 / avx512bw` at runtime |
| Linux arm64 | NEON (aarch64 ABI baseline) | Any aarch64 host | Single tier; NEON is mandatory in the aarch64 ABI |
| macOS arm64 | NEON | Apple M1 (2020) | Apple Silicon only; macOS x86_64 is unsupported |

## How to verify

```text
$ bwa-mem3 version
v0.2.0-12-gabcdef1
Compiler: clang X.Y.Z
SIMD floor: avx2 (x86-64-v3, Haswell 2013+); kernels: sse41 sse42 avx avx2 avx512bw
SIMD runtime: avx512bw (BWAMEM3_FORCE_TIER unset)
mimalloc 3.x.x (active)
```

- The **`SIMD floor:`** line tells you what host features the binary requires.
- The **`SIMD runtime:`** line tells you what kernel tier was selected at startup.
- On a host below the floor, `bwa-mem3 version` writes a `[W::bwa-mem3]` warning line to **stderr** (not stdout) and still exits 0, so the diagnostic command stays usable even on hosts that cannot run alignment. The floor + runtime lines remain on stdout, so `bwa-mem3 version | grep '^SIMD'` works in CI scripts even on too-old hosts.

## Failure mode on too-old hosts

If you run `bwa-mem3 mem` (or another alignment subcommand) on a host below the floor, the binary refuses with exit code 2 and a stderr message identifying the gap:

```text
[E::bwamem3] this binary was compiled for SIMD floor avx2 and emits avx2
instructions in non-kernel translation units. The host CPU does not support
avx2 (detected: sse42). Running would SIGILL on the first avx2 instruction.

bwa-mem3's lowest supported floor is BASELINE_ARCH=avx2; this host is below
avx2 and cannot run bwa-mem3.
```

AVX2 is the minimum SIMD floor: it cannot be lowered, so a host below AVX2 (Intel pre-Haswell, AMD pre-Zen) is unsupported. On an AVX2 host running an AVX-512 build, rebuilding with `BASELINE_ARCH=avx2` produces a runnable binary.

The `version` subcommand stays exit-0 so introspection still works on the same host.

## Mixed-architecture fleets

For AWS Batch and other heterogeneous compute environments where the same job may schedule onto x86_64 *or* arm64 hosts, see [Best Practices → Multi-architecture deployment](../best-practices/multi-arch-deployment.md).

---

**See also:**
[Installation](installation.md) ·
[Coming from bwa or bwa-mem2](migrating.md) ·
[Quick start: align paired-end FASTQs](quick-align.md) ·
[Performance → SIMD dispatch matrix](../performance/simd-dispatch.md) ·
[Best Practices → Multi-architecture deployment](../best-practices/multi-arch-deployment.md)
