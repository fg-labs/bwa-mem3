# Building from source

This page documents every build target available in the Makefile and what each produces. For the recommended production build workflow see [Best Practices → Build](../best-practices/build.md).

## Prerequisites

- A C++14-capable compiler, **modern by default**. **Clang is recommended** — across the benchmark sweep (hg38, 5M read pairs, `-t 16`) bwa-mem3 built with clang runs ~5–10% faster on x86 (AVX2, on AWS `c6a`; the AVX-512 `c7a` gain is smaller — see the per-arch record) and ~6% faster on ARM (NEON, on AWS Graviton4 `c8g`) than the same source built with g++, because older compilers emit redundant SIMD ops a modern clang already elides (full per-arch record — workload, host, architecture, and SIMD tier — in [Best Practices → Build](../best-practices/build.md)). To keep an accidentally-slow build from being mistaken for bwa-mem3 being slow, the build **fails by default on a compiler below the supported floor**: **clang ≥ 19** (recommended — what bioconda and the performance benchmarks ship), **GCC ≥ 15** (acceptable, ~5% slower on ARM), or **Apple clang ≥ 15** on macOS. A bare `make` defaults to g++: it builds (with a warning suggesting clang) only when g++ meets the **GCC ≥ 15** floor, and otherwise fails at the floor check as above; pass `CXX=clang++ CC=clang` to build with clang. To build anyway on an older toolchain (packagers on older distros, compatibility testing, deliberate A/B), pass `ALLOW_UNSUPPORTED_COMPILER=1`. The floor versions are overridable — e.g. `GCC_MIN=13` — for a distro that sets its own policy.
- GNU make 3.81+.
- CMake 3.12+ (required only when `USE_MIMALLOC=1`, which is the default).
- autoconf, automake, autoconf-archive, libtool, pkg-config — `ext/htslib`'s build runs `autoreconf -i && ./configure` and locates zlib via `pkg-config`.
- zlib development headers — htslib links against zlib.
- libdeflate development headers — `src/fast_reader.c` uses libdeflate for BGZF block decode (htslib also links it transitively). Debian/Ubuntu: `libdeflate-dev`; RHEL/Fedora: `libdeflate-devel`; macOS: `brew install libdeflate` (the Makefile auto-detects the Homebrew prefix or honours `LIBDEFLATE_PREFIX`). **Amazon Linux 2023 ships no `libdeflate-devel`** — build *and install* it from source (e.g. libdeflate v1.22): `cmake -B build -DCMAKE_INSTALL_PREFIX=/usr/local -DLIBDEFLATE_BUILD_SHARED_LIB=OFF`, then `cmake --build build && sudo cmake --install build`, and set `LIBRARY_PATH=/usr/local/lib64:/usr/local/lib` (CMake installs to `lib` or `lib64` depending on the distro) before `make`.
- OpenMP runtime — libsais uses OpenMP for parallel suffix-array construction. Linux + GCC: libgomp ships with the compiler, nothing extra to install. Linux + Clang: `libomp-dev` (Debian) / `libomp-devel` (RHEL). macOS: `brew install libomp`; the Makefile auto-detects the Homebrew prefix or honours `LIBOMP_PREFIX`.
- `ld.lld` (optional, recommended with clang on Linux) — the default build links `bwa-mem3` with LTO through lld; without it the build still succeeds, without LTO. Debian/Ubuntu: `lld` (or `lld-19` from apt.llvm.org); RHEL/Fedora: `lld`; the LLVM release tarballs include it. See [Link-Time Optimization](#link-time-optimization-lto).
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
make                       # g++ ≥ 15: builds + warns; older g++ fails (see Prerequisites)
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
| `make arch=avx2` | AVX2 | `-mavx2` |
| `make arch=avx512bw` | AVX-512BW | `-mavx512f -mavx512bw -mprefer-vector-width=256` |
| `make arch=native` | host CPU features | `-march=native` |

`arch=sse41`, `arch=sse42`, and `arch=avx` are refused at configure time: the batched mate-rescue kernels require AVX2+ and the pre-AVX2 scalar fallback was removed, so AVX2 is the minimum x86 tier.

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

The default build links the `bwa-mem3` binary with clang link-time optimization where the toolchain supports it, so a plain `make CXX=clang++ CC=clang` on Linux with `ld.lld` installed already produces an LTO binary. `LTO=` selects the mode:

| `LTO=` | Effect |
|---|---|
| `auto` _(default)_ | Full LTO when `CXX` and `CC` are clang on Linux and a test compile + link with `-flto=full -ffat-lto-objects -fuse-ld=lld` against the system libraries the build links succeeds. Otherwise a non-LTO build, with a `NOTE` on stderr when clang is found but cannot be used (the test link fails, e.g. because `ld.lld` is not installed, or `CC` is not clang). |
| `full` | Full (monolithic) LTO. Fails at parse time when it cannot be had: the test link fails, `CXX` or `CC` is not clang, or the host is macOS. |
| `thin` | ThinLTO; same requirements. Measured 0.3-0.9% slower wall time than `full`, with a faster, parallel link. |
| `off` | No LTO. |

`LTO` resolves only for goals that build (clean-only goals leave it off), and `make -s print-BWA_LTO_MODE` prints the mode a build would use.

What it needs: clang with `-ffat-lto-objects` (clang 18 and later; the build's floor is clang 19) for both `CXX` and `CC`, and `ld.lld` (Debian/Ubuntu `lld` or `lld-19`; it is included in the LLVM release tarballs). GCC builds, macOS builds (Apple clang has no `-ffat-lto-objects` for Mach-O) and `ASAN=1` / `MSAN=1` / `COVERAGE=1` builds link without LTO. Unlike `ld.bfd`, lld has no built-in library search path, so a library installed under `/usr/local` (the Amazon Linux libdeflate route above) needs the `LIBRARY_PATH` setting described there; without it the test link fails and the build falls back to no LTO. A packaged build (for example a conda recipe) gets LTO only if its build environment provides `ld.lld`; otherwise it builds exactly as before.

Only the objects that go into the binary are compiled for LTO, and only the link that produces `bwa-mem3` is an LTO link. The objects are *fat* LTO objects: each `.o` carries native code as well as the bitcode the LTO link reads. `libbwa.a` therefore stays an ordinary static archive for anything else that links it (the standalone tests, the doctest suites in `test/`, out-of-tree harnesses): those links use the native code, GNU `ar` indexes it without an LTO plugin, and no test binary pays for an LTO link. Third-party libraries (`htslib`, `mimalloc`, `zlib-ng`) are linked without LTO; libsais is LTO'd with the rest.

Measured on hg38 with clang 19.1.7 (HG00096 5M-pair WGS slice at `-t 16`, its first 200k pairs at `-t 1`), 6 rotated rounds per build, medians of the per-round change against `LTO=off`. Every LTO build was faster than `LTO=off` in all 6 rounds, and the SAM output (without `@PG`) was byte-identical across `off`, `thin` and `full` in every run, including `--compat=bwa-mem2` and, on x86, with `BWAMEM3_FORCE_TIER=avx2` and `avx512bw`.

| Host | Run | `LTO=off` wall | `thin` wall / CPU | `full` wall / CPU |
|---|---|---|---|---|
| AWS Graviton4 `c8g.4xlarge` (NEON) | 5M, `-t 16` | 42.6 s | -1.2% / -1.2% | **-1.7% / -1.7%** |
| AWS Graviton4 `c8g.4xlarge` (NEON) | 200k, `-t 1` | 29.5 s | -1.0% / -1.1% | **-1.7% / -1.7%** |
| AMD Zen 4 `c7a.4xlarge` (AVX-512BW) | 5M, `-t 16` | 42.3 s | -2.8% / -3.1% | **-3.1% / -3.1%** |
| AMD Zen 4 `c7a.4xlarge` (AVX-512BW) | 200k, `-t 1` | 27.4 s | -2.6% / -2.9% | **-3.5% / -3.8%** |

The binary shrinks from 8.7 MB to 7.9 MB on arm64 and from 9.1 MB to 8.1 MB on x86 with full LTO.

The gain comes from cross-file inlining of the glue around the SIMD kernels (seeding, chaining, the extension driver, paired-end post-processing); the kernels themselves are single translation units and do not change. On x86 the per-tier kernel objects keep their own ISA: the inliner does not move a function into a caller compiled for a narrower ISA, and a scan of the linked binaries for AVX-512 instructions finds them in the same functions with and without LTO. The cost is build time: the full-LTO link of `bwa-mem3` is single-threaded and adds about 15-20 s to a clean `-j16` build (ThinLTO: 2.5-6 s). For a faster edit-compile loop, build with `LTO=thin` or `LTO=off`.

The older `lto-build` target is still available for GCC LTO and single-ISA cross-builds; it turns the default LTO off for its own build:

```bash
make lto-build                                 # builds bwa-mem3.lto (native arch)
make lto-build LTO_ARCH=avx2                   # explicit arch
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
| `COVERAGE` | _(unset)_ | Set to enable `--coverage` at `COVERAGE_OPT` for gcov line-level coverage |
| `COVERAGE_OPT` | `-O0` | Optimization level of a `COVERAGE=1` build. `-O0` maps gcov lines 1:1 to source; CI's coverage job uses `-O1`, which runs the tests faster but reports some inlined lines and small functions under their caller |
| `LTO` | `auto` | Link-time optimization: `auto`, `thin`, `full` or `off` (see [Link-Time Optimization](#link-time-optimization-lto)) |
| `EXTRA_CXXFLAGS` | _(empty)_ | Appended to `CXXFLAGS`; forwarded through the PGO targets (`pgo-generate`, `pgo-use`), `profile-build` and `lto-build`, ahead of the flags each of those adds |
| `EXTRA_CFLAGS` | _(empty)_ | Appended last to the FASTQ reader's C flags (`-g -O3` plus the arch flags), e.g. `EXTRA_CFLAGS=-O0` to debug it; forwarded through the PGO targets and `lto-build` |
| `CPPFLAGS` | _(environment)_ | Preprocessor flags from the environment, e.g. a packager's `-D_FORTIFY_SOURCE=2`; the Makefile appends its own defines to them. Any `-O` level in an environment `CPPFLAGS` (conda-forge's compiler activation exports `-O2` there) is stripped, so the Makefile's per-target levels win: `-O3` for the aligner, `-O1` under `ASAN`/`MSAN`, and `COVERAGE_OPT` under `COVERAGE=1`. Set the optimization level with `EXTRA_CXXFLAGS` / `EXTRA_CFLAGS` instead. A `-DNDEBUG` in `CPPFLAGS` (conda's has one) is stripped too (the single-token spellings `-DNDEBUG` and `-DNDEBUG=<v>` only), so plain `assert()` calls stay live in packaged builds as in every other build of this repository without an explicit opt-out; to build bwa-mem3's own sources without them on purpose, pass `EXTRA_CXXFLAGS=-DNDEBUG` (plus `EXTRA_CFLAGS=-DNDEBUG` for the C reader), which is not filtered (`xassert()` stays live either way). The opt-out does not reach the bundled htslib, which is configured with the filtered `CPPFLAGS` and so keeps its own asserts; a package that links a system htslib (Bioconda uses conda-forge's) is unaffected. Pass it in the environment, not as `make CPPFLAGS=...`: a command-line value replaces every define the Makefile adds itself (`-DV17=1`, `-DMATE_SORT=0`, `-DENABLE_PREFETCH`, `-DUSE_MIMALLOC=1`, the sse2neon `-D__SSE*__` shims, `-DCACHE_LINE_BYTES`, ...), and reaches `test/Makefile` unfiltered |
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
