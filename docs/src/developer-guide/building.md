# Building from source

This page documents every build target available in the Makefile and what each produces. For the recommended production build workflow see [Best Practices → Build](../best-practices/build.md).

## Prerequisites

- A C++14-capable compiler: GCC 7+ or Clang 6+ on Linux; Clang 15+ (Xcode) on macOS.
- GNU make 3.81+.
- CMake 3.12+ (required only when `USE_MIMALLOC=1`, which is the default).
- libomp (macOS only): `brew install libomp`. libsais uses OpenMP for parallel suffix-array construction.
- Git submodules initialised: `git submodule update --init --recursive`.

> **Warning — Submodules must be present**
>
> The build will fail with a clear error message if any of the required submodules
> (`ext/libsais`, `ext/htslib`, `ext/safestringlib`, `ext/mimalloc`, `ext/sse2neon`) are missing.
> Always clone with `--recursive` or run `git submodule update --init --recursive` before `make`.

## Standard builds

### Default build (host-native)

```bash
make
```

On x86 hosts this is equivalent to `make multi` (see below). On Apple Silicon and other aarch64 hosts the Makefile detects the architecture and builds a single ARM64 binary instead.

The resulting binary is `bwa-mem3` in the repo root.

### Single-arch x86 builds

Pass `arch=<target>` to compile a single binary with a specific ISA level:

| Command | SIMD level | `ARCH_FLAGS` |
|---|---|---|
| `make arch=sse41` | SSE4.1 | `-msse … -msse4.1` |
| `make arch=sse42` | SSE4.2 | `-msse … -msse4.2` |
| `make arch=avx` | AVX | `-mavx` |
| `make arch=avx2` | AVX2 | `-mavx2` |
| `make arch=avx512bw` | AVX-512BW | `-mavx512f -mavx512bw` |
| `make arch=native` | host CPU features | `-march=native` |

For Intel compiler (`icpc` / `icpx`) the flags differ slightly; see the Makefile for the `ifeq ($(CXX), icpc)` branches.

### Multi-binary x86 build (default on x86)

```bash
make multi
```

Builds five ISA-specific binaries (`bwa-mem3.sse41`, `bwa-mem3.sse42`, `bwa-mem3.avx`, `bwa-mem3.avx2`, `bwa-mem3.avx512bw`) plus the thin launcher `bwa-mem3` that execs the best-matching binary at runtime. See [Multi-binary launcher](launcher.md) for details.

### ARM64 / Apple Silicon build

```bash
make arch=arm64
```

Compiles a single binary `bwa-mem3.arm64` and creates a symlink `bwa-mem3 -> bwa-mem3.arm64`. See [Apple Silicon / NEON port](neon-port.md) for background.

## Tuned builds

### Profile-Guided Optimization (PGO)

PGO produces the best single-binary performance. The workflow is two-phase:

```bash
# Phase 1: instrument binary
make pgo-generate                              # builds bwa-mem3.pgo-instr (arm64 default)
make pgo-generate PGO_ARCH=avx2               # or a specific x86 target

# Run your training workload with the instrumented binary
./bwa-mem3.pgo-instr mem -t 16 ref.fa r1.fq.gz r2.fq.gz > /dev/null

# Phase 2: optimised binary
make pgo-use                                   # builds bwa-mem3.pgo
make pgo-use PGO_ARCH=avx2                     # matching arch
```

`PGO_ARCH` accepts the same values as `arch=`. `PGO_PROFILE_DIR` defaults to `pgo_profiles/` but can be overridden. Output binaries are named `bwa-mem3.pgo` (default arch) or `bwa-mem3.pgo.<arch>` when a non-default arch is specified, so multiple arch builds coexist.

Clean up instrumented objects and profile data:

```bash
make pgo-clean
```

### Link-Time Optimization (LTO)

```bash
make lto-build                                 # builds bwa-mem3.lto (native arch)
make lto-build LTO_ARCH=avx2                   # explicit arch
```

LTO compiles bwa-mem3's own translation units with `-flto` (thin LTO on Clang, full LTO on GCC) plus `-fno-semantic-interposition` on GCC. Third-party libraries (`htslib`, `mimalloc`, `safestringlib`) are linked without LTO. Clean:

```bash
make lto-clean
```

### Compute-only profile binary

Used when profiling CPU hotspots without I/O noise. The `-DDISABLE_OUTPUT` flag short-circuits all BAM/SAM write paths and the file-open / header-emit step, so only alignment work contributes to wall time.

```bash
make profile-build                             # builds bwa-mem3.profile (native)
make profile-build PROFILE_ARCH=avx2          # explicit arch
./bwa-mem3.profile mem -t 16 ref.fa r1.fq.gz r2.fq.gz

make profile-clean
```

## Build knobs

| Variable | Default | Effect |
|---|---|---|
| `USE_MIMALLOC` | `1` | Include mimalloc; set `0` to use the system allocator |
| `ASAN` | _(unset)_ | Set to any non-empty value to enable AddressSanitizer (forces `USE_MIMALLOC=0`) |
| `COVERAGE` | _(unset)_ | Set to enable `--coverage` + `-O0` for gcov line-level coverage |
| `EXTRA_CXXFLAGS` | _(empty)_ | Appended to `CXXFLAGS`; forwarded through PGO / LTO targets |
| `DISABLE_BATCHED_MATESW` | _(unset)_ | Set to `1` to disable the batched mate-rescue SW path on ARM |
| `CXX` | `c++` | Compiler. Paired `CC` is auto-derived from `CXX` for libsais. |

## Cleaning

```bash
make clean
```

Removes object files, `libbwa.a`, all binaries, test binaries, libsais objects, safestringlib, htslib, and the mimalloc build tree.

```bash
make docs-clean
```

Removes only the mdbook build output (`docs/book/`). Covered in [Developer Guide → Building](building.md) context; see the Makefile docs targets for the full list.

## Documentation targets

| Target | Action |
|---|---|
| `make docs` | Build the mdbook into `docs/book/` |
| `make docs-serve` | Live-preview at `http://localhost:3000` |
| `make docs-cli` | Capture `--help` output for each subcommand into `docs/_generated/cli/` |
| `make docs-clean` | Remove `docs/book/` |
| `make docs-install-tools` | `cargo install` mdbook + three plugins |

---

**See also:**
[SIMD dispatch architecture](simd-dispatch.md) ·
[Multi-binary launcher](launcher.md) ·
[Best Practices → Build](../best-practices/build.md) ·
[Performance → PGO build](../performance/pgo.md) ·
[Apple Silicon / NEON port](neon-port.md)
