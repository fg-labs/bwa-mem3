# Building from source

This page documents every build target available in the Makefile and what each produces. For the recommended production build workflow see [Best Practices → Build](../best-practices/build.md).

## Prerequisites

- A C++14-capable compiler. **Clang is recommended** — bwa-mem3 runs ~5–10% faster on x86 (AVX2) and ~6% faster on ARM when built with clang than with g++ (see [Best Practices → Build](../best-practices/build.md)). Clang 7+ or GCC 8+ on Linux; Clang 15+ (Xcode) on macOS. A bare `make` defaults to g++ and prints a warning suggesting clang; pass `CXX=clang++ CC=clang` to build with clang.
- GNU make 3.81+.
- CMake 3.12+ (required only when `USE_MIMALLOC=1`, which is the default).
- autoconf, automake, autoconf-archive, libtool, pkg-config — `ext/htslib`'s build runs `autoreconf -i && ./configure` and locates zlib via `pkg-config`.
- zlib development headers — htslib links against zlib.
- libdeflate development headers — `src/fast_reader.c` uses libdeflate for BGZF block decode (htslib also links it transitively). Debian/Ubuntu: `libdeflate-dev`; RHEL/Fedora: `libdeflate-devel`; macOS: `brew install libdeflate` (the Makefile auto-detects the Homebrew prefix or honours `LIBDEFLATE_PREFIX`). **Amazon Linux 2023 ships no `libdeflate-devel`** — build *and install* it from source (e.g. libdeflate v1.22): `cmake -B build -DCMAKE_INSTALL_PREFIX=/usr/local -DLIBDEFLATE_BUILD_SHARED_LIB=OFF`, then `cmake --build build && sudo cmake --install build`, and set `LIBRARY_PATH=/usr/local/lib64:/usr/local/lib` (CMake installs to `lib` or `lib64` depending on the distro) before `make`.
- OpenMP runtime — libsais uses OpenMP for parallel suffix-array construction. Linux + GCC: libgomp ships with the compiler, nothing extra to install. Linux + Clang: `libomp-dev` (Debian) / `libomp-devel` (RHEL). macOS: `brew install libomp`; the Makefile auto-detects the Homebrew prefix or honours `LIBOMP_PREFIX`.
- Git submodules initialised: `git submodule update --init --recursive`.

See [Getting Started → Installation](../getting-started/installation.md) for the full per-platform install commands.

> **Warning — Submodules must be present**
>
> The build will fail with a clear error message if any of the required submodules
> (`ext/libsais`, `ext/htslib`, `ext/mimalloc`, `ext/sse2neon`) are missing.
> Always clone with `--recursive` or run `git submodule update --init --recursive` before `make`.

## Standard builds

### Default build (host-native)

```bash
make                       # uses g++ (warns, suggesting clang)
make CXX=clang++ CC=clang  # recommended: ~5–10% faster on x86, ~6% on ARM
```

On x86 hosts this is equivalent to `make single` (see below): one binary
containing all five SIMD tiers, dispatched in process at startup. On
Apple Silicon and other aarch64 hosts the Makefile detects the
architecture and builds a single ARM64 binary with one NEON kernel TU.

The resulting binary is `bwa-mem3` in the repo root.

### Single multi-tier x86 build (default on x86)

```bash
make single                       # alias of the default `make`
make BASELINE_ARCH=avx512bw       # raise non-kernel TU compile baseline
make BASELINE_ARCH=sse41          # lower it for pre-Haswell hosts
```

Builds one `bwa-mem3` binary. The four hand-tuned kernel TUs in
`KERNEL_SRCS` (`bandedSWA.cpp`, `kswv.cpp`, `ksw.cpp`,
`sam_encode.cpp`) are compiled five times each — once per supported
tier (`sse41` / `sse42` / `avx` / `avx2` / `avx512bw`) — and dispatched
at runtime via `__builtin_cpu_supports`. Non-kernel TUs compile once at
`BASELINE_ARCH` (default `avx2` since PR #84). See
[Single-binary SIMD dispatch (x86)](launcher.md) for the full design.

### Single-tier x86 builds

Pass `arch=<target>` to compile a single binary with kernels for **one
tier only** (no runtime dispatch table — useful on clusters with uniform
hardware):

| Command | SIMD level | `ARCH_FLAGS` |
|---|---|---|
| `make arch=sse41` | SSE4.1 | `-msse … -msse4.1` |
| `make arch=sse42` | SSE4.2 | `-msse … -msse4.2` |
| `make arch=avx` | AVX | `-mavx` |
| `make arch=avx2` | AVX2 | `-mavx2` |
| `make arch=avx512bw` | AVX-512BW | `-mavx512f -mavx512bw -mprefer-vector-width=256` |
| `make arch=native` | host CPU features | `-march=native` |

For Intel compiler (`icpc` / `icpx`) the flags differ slightly; see the
Makefile for the `ifeq ($(CXX), icpc)` branches. The `avx512bw` target
keeps the `-mprefer-vector-width=256` cap from PR #86 — see
[`BASELINE_ARCH=avx512bw` build flag](../whats-different/avx512-baseline.md)
for the empirical perf characterization.

### ARM64 / Apple Silicon build

```bash
make arch=arm64
```

Compiles a single binary `bwa-mem3` with one NEON kernel TU. See
[Apple Silicon / NEON port](neon-port.md) for background.

## Tuned builds

### Profile-Guided Optimization (PGO)

PGO produces the best single-binary performance. The workflow is two-phase:

```bash
# Phase 1: instrument binary
make pgo-generate                              # builds bwa-mem3.pgo-instr (arm64 default)
make pgo-generate PGO_ARCH=avx2               # or a specific x86 target

# Run your training workload with the instrumented binary
./bwa-mem3.pgo-instr mem -t 16 ref.fa R1.fq.gz R2.fq.gz > /dev/null

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

LTO compiles bwa-mem3's own translation units with `-flto` (thin LTO on Clang, full LTO on GCC) plus `-fno-semantic-interposition` on GCC. Third-party libraries (`htslib`, `mimalloc`) are linked without LTO. Clean:

```bash
make lto-clean
```

### Compute-only profile binary

Used when profiling CPU hotspots without I/O noise. The `-DDISABLE_OUTPUT` flag short-circuits all BAM/SAM write paths and the file-open / header-emit step, so only alignment work contributes to wall time.

```bash
make profile-build                             # builds bwa-mem3.profile (native)
make profile-build PROFILE_ARCH=avx2          # explicit arch
./bwa-mem3.profile mem -t 16 ref.fa R1.fq.gz R2.fq.gz

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

## Incremental builds and header dependencies

Every compile emits a sidecar `.d` file next to its object listing the headers that
translation unit included (`-MMD -MP`), and the Makefile includes them, so editing a header
rebuilds exactly the objects that read it. There is no `make depend` step to remember.

Objects additionally depend on `.build-flags`, a stamp holding the expanded compile flags.
It is rewritten only when that text changes, so switching flags (`make arch=avx2` after
`arch=sse41`, adding `ASAN=1`, changing `CXX`) rebuilds the objects built with the old ones,
while editing an unrelated part of the Makefile does not.

Both are build artifacts: `make clean` removes them and `.gitignore` covers them.
`test/regression/make_header_deps.sh` guards the mechanism — it enumerates the Makefile's
object lists and checks that every one of them that the tree has actually *built* has a
generated `.d`, so a compile rule added without `$(DEPFLAGS)` fails the check. Objects absent
from the build tree are skipped (a run that skipped everything fails), so run it against a
tree built the way you want covered. CI runs it on the canonical row.

The Makefile also states its default goal (`.DEFAULT_GOAL`) rather than letting make fall back
to the first target of the first rule. `myall` has a rule only for a bare `make`; once `arch=`
is set — which CI's build step and `single:`'s own recursion both do — the goal is `all`, and
without the explicit setting any rule added above `all:` would take that slot instead. An
object landing there makes a build compile one file, exit 0, and produce no binary.
`test/regression/make_default_goal.sh` guards it, also on the canonical row.

## Cleaning

```bash
make clean
```

Removes object files, dependency files, `libbwa.a`, all binaries, test binaries, libsais objects, htslib, and the mimalloc build tree.

```bash
make docs-clean
```

Removes only the mdbook build output (`docs/book/`). See the [Documentation targets](#documentation-targets) below for the full list.

## Documentation targets

| Target | Action |
|---|---|
| `make docs` | Build the mdbook into `docs/book/html/` |
| `make docs-serve` | Live-preview at `http://localhost:3000` |
| `make docs-cli` | Capture `--help` output for each subcommand into `docs/_generated/cli/` |
| `make docs-clean` | Remove `docs/book/` |
| `make docs-install-tools` | `cargo install` mdbook, mdbook-mermaid, and mdbook-linkcheck2 |

The build runs the [mdbook-linkcheck2](https://github.com/marxin/mdbook-linkcheck2) backend, which **fails the build on a dead internal link** (a link to a page that does not exist). This guards against broken cross-references reaching the published site — mdBook on its own only warns. External (web) links are not checked, and bracketed literal text in the captured CLI snippets (e.g. `[P]`) is reported only as a non-fatal warning. Because a second output backend is configured, the HTML site is written to `docs/book/html/` rather than `docs/book/`.

---

**See also:**
[SIMD dispatch architecture](simd-dispatch.md) ·
[Single-binary SIMD dispatch (x86)](launcher.md) ·
[Best Practices → Build](../best-practices/build.md) ·
[Performance → PGO build](../performance/pgo.md) ·
[Apple Silicon / NEON port](neon-port.md)
