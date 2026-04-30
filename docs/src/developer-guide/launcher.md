# Multi-binary launcher (x86)

On x86 Linux and x86 macOS, `bwa-mem3` is a thin launcher binary rather than the aligner itself. Its sole job is to detect the host CPU's capabilities at startup and exec the best-matching ISA-specific binary in the same directory.

ARM / Apple Silicon does not use this mechanism. The `make arm64` target creates a symlink `bwa-mem3 -> bwa-mem3.arm64`; there is only one NEON instruction-set level on all current ARM64 CPUs.

## What `make multi` produces

```bash
make multi
```

Produces six files in the repo root:

| File | ISA | ARCH_FLAGS |
|---|---|---|
| `bwa-mem3` | launcher (no aligner code) | compiled from `src/runsimd.cpp` |
| `bwa-mem3.sse41` | SSE4.1 | `-msse -msse2 -msse3 -mssse3 -msse4.1` |
| `bwa-mem3.sse42` | SSE4.2 | adds `-msse4.2` |
| `bwa-mem3.avx` | AVX | `-mavx` |
| `bwa-mem3.avx2` | AVX2 | `-mavx2` |
| `bwa-mem3.avx512bw` | AVX-512BW | `-mavx512f -mavx512bw` |

The six binaries must reside in the same directory for the launcher to find them.

## How the launcher selects a binary

`src/runsimd.cpp` calls `cpuid` (via the `__cpuid` intrinsic or a hand-rolled `CPUID` wrapper) to read the CPU's feature flags and picks the highest ISA level supported by the CPU:

1. Check `CPUID` leaf 7 for AVX-512BW → exec `bwa-mem3.avx512bw`
2. Check `CPUID` leaf 7 for AVX2 → exec `bwa-mem3.avx2`
3. Check `CPUID` leaf 1 for AVX → exec `bwa-mem3.avx`
4. Check `CPUID` leaf 1 for SSE4.2 → exec `bwa-mem3.sse42`
5. Fallback → exec `bwa-mem3.sse41`

The launcher calls `execv` with the same `argv` that was passed to it. The selected binary's `main()` therefore receives the original arguments unchanged. The `@PG CL:` tag in the output SAM/BAM records the original invocation, not the ISA-suffixed binary name.

> **Note — exec replaces the process**
>
> The launcher does not fork. It calls execv(), which replaces the launcher process
> image with the ISA-specific binary. There is no wrapper process resident in memory
> during alignment.

## Using a specific ISA binary directly

You can bypass the launcher and invoke a specific binary directly:

```bash
./bwa-mem3.avx2 mem -t 16 ref.fa r1.fq.gz r2.fq.gz
```

This is useful when benchmarking a particular ISA level, testing a regression, or deploying in an environment where only one binary is installed. The ISA-specific binary behaves identically to the launcher output for that ISA — there is no functional difference.

## Distribution layout

When packaging or deploying on x86, include all five ISA binaries plus the launcher in the same directory:

```text
bin/
  bwa-mem3           ← launcher
  bwa-mem3.sse41
  bwa-mem3.sse42
  bwa-mem3.avx
  bwa-mem3.avx2
  bwa-mem3.avx512bw
```

On ARM, only `bwa-mem3` (the symlink) and `bwa-mem3.arm64` are needed.

## The `mem` SIMD banner

After selecting and executing a binary, the `mem` subcommand prints a single-line banner to stderr before alignment begins:

```text
-----------------------------
Executing in AVX2 mode!!
-----------------------------
```

The banner text is set at compile time via `#if __AVX512BW__` / `#elif __AVX2__` / … preprocessor guards in `src/main.cpp`. This confirms at runtime which ISA path is active.

---

**See also:**
[SIMD dispatch architecture](simd-dispatch.md) ·
[Apple Silicon / NEON port](neon-port.md) ·
[Building from source](building.md) ·
[Performance → SIMD dispatch matrix](../performance/simd-dispatch.md)
