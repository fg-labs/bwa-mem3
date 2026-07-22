# Installation

Already using `bwa` or `bwa-mem2`? After installing, see
[Coming from bwa or bwa-mem2](migrating.md) — the command line is unchanged, but
the index and output have a couple of migration notes.

## Bioconda

bwa-mem3 is published on [Bioconda](https://bioconda.github.io/) for linux-64,
linux-aarch64, osx-64, and osx-arm64. This is the recommended path for most
users:

```bash
conda install -c bioconda -c conda-forge bwa-mem3
```

Build from source (below) if you need a compiler/arch-tuned build — in
particular a clang-built, arch-specific, or PGO binary for maximum throughput
(see [Best Practices → Build](../best-practices/build.md)).

## Build from source

### Prerequisites

bwa-mem3 vendors several libraries as git submodules. Building from source requires the toolchain to compile bwa-mem3 itself plus the bootstrap tools each vendored library needs.

| Tool | Why it's needed | Minimum version |
|------|-----------------|-----------------|
| C++14 compiler (**Clang recommended**, or GCC) | bwa-mem3 itself — clang builds run ~5–10% faster on x86, see [Build](../best-practices/build.md) | Linux: Clang 7+ / GCC 8+ · macOS: Clang 15+ (Xcode) |
| GNU make | top-level build | 3.81+ |
| Git | submodule checkout (with `--recursive`) | any recent |
| autoconf, automake, autoconf-archive, libtool | `ext/htslib` runs `autoreconf -i && ./configure` during build | any recent |
| pkg-config | htslib's `configure` uses it to locate zlib | any recent |
| zlib development headers | htslib links against zlib | any recent |
| libdeflate development headers | `src/fast_reader.c` uses libdeflate for BGZF block decode | any recent |
| OpenMP runtime | `ext/libsais` uses OpenMP for parallel suffix-array construction | see notes below |
| CMake 3.12+ | building bundled mimalloc (default; skip if you pass `USE_MIMALLOC=0`) | 3.12+ |

> **OpenMP notes.**
> - On **Linux with GCC**, libgomp ships with the compiler — no extra package needed.
> - On **Linux with Clang**, install `libomp-dev` (Debian/Ubuntu) or `libomp-devel` (RHEL/Fedora).
> - On **macOS**, install Homebrew's `libomp` (`brew install libomp`). The Makefile auto-detects the Homebrew prefix; set `LIBOMP_PREFIX=/path/to/libomp` if you installed it elsewhere.

#### Install prerequisites by platform

**Debian / Ubuntu:**
```bash
sudo apt-get install \
    build-essential git cmake pkg-config \
    autoconf automake autoconf-archive libtool \
    zlib1g-dev libdeflate-dev \
    libomp-dev          # only needed if building with Clang
```

**RHEL / Fedora / Amazon Linux:**
```bash
sudo dnf install \
    gcc gcc-c++ make git cmake pkgconf-pkg-config \
    autoconf automake autoconf-archive libtool \
    zlib-devel libdeflate-devel \
    libomp-devel        # only needed if building with Clang
```

> **Amazon Linux 2023 has no `libdeflate-devel`.** The `dnf install` above will fail on that one
> package. Build *and install* libdeflate from source instead (e.g. v1.22): configure with
> `cmake -B build -DCMAKE_INSTALL_PREFIX=/usr/local -DLIBDEFLATE_BUILD_SHARED_LIB=OFF`, then
> `cmake --build build && sudo cmake --install build` to put the headers and static library under
> `/usr/local`. CMake installs the library to `lib` or `lib64` depending on the distro, so set
> `LIBRARY_PATH=/usr/local/lib64:/usr/local/lib` (covering both) before running bwa-mem3's `make`.
> Other RHEL/Fedora releases have the package.

**macOS (Homebrew):**
```bash
xcode-select --install   # Apple Clang + git + make
brew install \
    cmake pkg-config \
    autoconf automake autoconf-archive libtool \
    libdeflate libomp
```

> **What happens if a prereq is missing.** The Makefile fails fast with an actionable error: a missing `libomp` on macOS, a missing `autoreconf`, or a missing `cmake` each produce a one-line hint pointing at the install command above. There is no need to install everything optimistically — install only what the error message asks for if you prefer.

### Clone and build

```bash
git clone --recursive https://github.com/fg-labs/bwa-mem3
cd bwa-mem3
make CXX=clang++ CC=clang
```

> **Build with clang.** bwa-mem3 is ~5–10% faster on x86 (AVX2) and ~6% faster on
> ARM when built with clang instead of g++ — see
> [Best Practices → Build](../best-practices/build.md). A bare `make` uses g++
> and prints a warning nudging you to clang. Pass `CXX=clang++ CC=clang` as
> shown. (On Linux with clang, install `libomp-dev` / `libomp-devel` for
> OpenMP; see the prerequisites above.)

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
v0.2.0-12-gabcdef1
Compiler: clang X.Y.Z
SIMD floor: avx2 (x86-64-v3, Haswell 2013+); kernels: sse41 sse42 avx avx2 avx512bw
SIMD runtime: avx2 (BWAMEM3_FORCE_TIER unset)
mimalloc 3.x.x (active)
```

The `Compiler:` line reports which toolchain built the binary and its version —
`clang`, `gcc`, `icpx`, or `icpc`. The exact version varies with your toolchain;
what matters is that a production build says `clang` (see
[Best Practices → Build](../best-practices/build.md)).

If the `mimalloc` line is absent, the build linked the system allocator (expected when
`USE_MIMALLOC=0` was passed or when the vendor submodule was not initialized). If it is
present but reads `(linked but NOT overriding malloc)` instead of `(active)`, mimalloc is
linked but not intercepting allocations — see
[User Guide → Memory allocator (mimalloc)](../user-guide/allocator.md).

## Next: host requirements

If you're planning to deploy bwa-mem3 across a heterogeneous fleet (AWS Batch, mixed compute clusters), read [Host requirements](host-requirements.md) for the supported CPU floor and [Best Practices → Multi-architecture deployment](../best-practices/multi-arch-deployment.md) for the deployment recipe.

---

**See also:**
[Quick start: align paired-end FASTQs](quick-align.md) ·
[User Guide — Memory allocator](../user-guide/allocator.md) ·
[Developer Guide — Building from source](../developer-guide/building.md)
