# Installation

## Bioconda (coming soon)

A Bioconda package for bwa-mem3 is in preparation. Once published, installation will be:

```bash
conda install -c bioconda bwa-mem3
```

This will be the recommended path for most users. Check back here or watch the
[fg-labs/bwa-mem3](https://github.com/fg-labs/bwa-mem3) repository for the announcement.

## Build from source

Until the Bioconda package is available, build from source using the steps below.

### Prerequisites

- A C++17 compiler (GCC 8+ or Clang 7+)
- GNU make
- Rust toolchain (for `cargo install` of mdbook tools, not required for the aligner itself)
- Git (for submodule checkout)

### Clone and build

```bash
git clone --recursive https://github.com/fg-labs/bwa-mem3
cd bwa-mem3
make
```

The `--recursive` flag is required. bwa-mem3 vendors several libraries (mimalloc, sse2neon, and
others) as git submodules. A shallow or non-recursive clone will fail to compile.

> **Warning — Shallow clone submodule pitfall**
>
> If you cloned without `--recursive`, initialize the submodules before running `make`:
>
> ```bash
> git submodule update --init --recursive
> ```
>
> Forgetting this step is the most common source of build failures.

### Target architecture

By default, `make` builds a general-purpose binary that runs on any supported CPU. For maximum
performance, specify the architecture that matches your deployment target:

| Flag | Requires | Notes |
|------|----------|-------|
| `make` | SSE4.1 or better (x86), any (ARM) | Default; selects best dispatch at runtime on x86 |
| `make arch=avx2` | AVX2 (e.g. Haswell, Zen 2) | Recommended for modern x86 servers |
| `make arch=avx512bw` | AVX-512BW (e.g. Skylake-X, Ice Lake, Sapphire Rapids) | Maximum x86 performance |
| `make arch=arm64` | Apple Silicon / AWS Graviton | NEON-vectorized build |

See [Performance — SIMD dispatch matrix](../performance/simd-dispatch.md) for the full matrix
of which kernels are vectorized under each target.

### Memory allocator (mimalloc)

bwa-mem3 bundles [mimalloc](https://github.com/microsoft/mimalloc) and links it into every binary
by default. mimalloc reduces allocator contention under high thread counts and lowers wall-clock
time on multi-threaded alignment runs.

To build without mimalloc, pass `USE_MIMALLOC=0`:

```bash
make USE_MIMALLOC=0
```

See [User Guide — Memory allocator](../user-guide/allocator.md) for details on how mimalloc is
linked on Linux versus macOS and when opting out is appropriate.

## Smoke test

After building, run the smoke test to confirm the binary works and report which allocator is
active:

```bash
./bwa-mem3 version
```

Expected output (with mimalloc):

```text
bwa-mem3-0.1.0-pre-12-gabcdef1
mimalloc 3.x.x
```

If the `mimalloc` line is absent, the build linked the system allocator (expected when
`USE_MIMALLOC=0` was passed or when the vendor submodule was not initialized).

---

**See also:**
[Quick start: align paired-end FASTQs](quick-align.md) ·
[User Guide — Memory allocator](../user-guide/allocator.md) ·
[Developer Guide — Building from source](../developer-guide/building.md)
