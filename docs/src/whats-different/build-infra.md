# Build & Infrastructure

This page covers the build-system, testing, and CI infrastructure changes
carried in bwa-mem3 on top of upstream bwa-mem2.

## doctest framework and Codecov (PR #34)

PR #34 establishes the long-term test infrastructure for bwa-mem3:

- **doctest 2.4.11** is vendored as a single-header under `ext/doctest/`, with
  the SHA256 recorded in `ext/doctest/VERSION`.
- A new `test/framework/` static library provides shared helpers: scoring
  matrices, deterministic sequence-pair generators, kswv-style batch packers,
  scalar and SIMD runners, kswr comparators, a JUnit reporter hook, and a
  shared `main`.
- Two test binaries are produced: `bwa_mem3_tests_unit` (runs on every CI
  matrix row) and `bwa_mem3_tests_integration` (runs on a subset of rows).
- The existing `kswv_selftest` is ported to
  `test/unit/test_kswv_correctness.cpp` — 30,049 assertions against scalar
  `ksw_align2` on 10k random plus curated edge pairs.
- Five legacy integration sources are moved to `test/integration/` via
  `git mv`; their binaries still emit at `test/<name>` so existing scripts
  keep working.
- Five inline CI bash regression blocks are extracted to
  `test/regression/*.sh` (phix_parity, chr22_parity, thread_determinism,
  bam_roundtrip, meth_oracle).
- A `coverage` CI job builds `libbwa.a` and both test binaries with
  `COVERAGE=1` (`-O0 --coverage`), runs both test binaries, collects Cobertura
  XML via `gcovr`, and uploads to Codecov via `codecov/codecov-action`.

## `PACKAGE_VERSION` from `git describe` (PR #52)

Before PR #52, `src/main.cpp` hardcoded `PACKAGE_VERSION "2.2.1"`. This string
appeared in `bwa-mem3 version` output and in the `@PG VN:` SAM header field
but was never updated, causing every build to report an outdated version.

The Makefile now generates `src/version.h` from `git describe --tags --dirty`,
falling back to a static `FG_LABS_VERSION_FALLBACK` when `git describe` cannot
reach a tag (source-tarball extractions, shallow clones — e.g. CI with the
default `fetch-depth: 1`). A write-if-changed mechanism (`cmp -s` + `mv`)
regenerates the file on every invocation but only bumps its mtime when the
stamped string changes, so only `main.o` is rebuilt when the version changes,
not the entire tree. `src/version.h` is `.gitignore`d and removed by
`make clean`. Fixes
[issue #40](https://github.com/fg-labs/bwa-mem3/issues/40). Related upstream:
[bwa-mem2#283](https://github.com/bwa-mem2/bwa-mem2/issues/283),
[bwa-mem2#284](https://github.com/bwa-mem2/bwa-mem2/pull/284).

## PGO target parameterization (PR #59)

The original `pgo-generate` and `pgo-use` Makefile targets hardcoded
`arch=arm64` and a single shared `pgo_profiles/` directory. PR #59 generalizes
both:

- `PGO_ARCH` (default: `arm64` on ARM hosts, `native` otherwise) passes
  through to the recursive `make` invocation as `arch=$(PGO_ARCH)`. Accepts
  the same values as the rest of the Makefile: `arm64`, `sse41`, `avx2`,
  `avx512bw`, `native`, etc.
- `PGO_PROFILE_DIR` is now overridable (`?=` instead of `=`). Each
  `(arch × training-regime)` combination can capture into its own directory.
- When `PGO_ARCH != arm64`, the output binaries are named
  `bwa-mem3.pgo-instr.<arch>` and `bwa-mem3.pgo.<arch>` so multiple per-arch
  PGO builds coexist. The default arm64 names are unchanged for backward
  compatibility.
- `pgo-clean` now removes arch-suffixed PGO binaries in addition to bare names.

This enables the benchmarking workflow at
[bwa-mem3-bench](https://github.com/fg-labs/bwa-mem3-bench), which requires
per-arch × per-regime profile capture. See also
[Performance → PGO build](../performance/pgo.md).

## `CXXFLAGS`/`CPPFLAGS`/`LDFLAGS` forwarding (PR #50)

At the time of PR #50, the Makefile's `multi:` rule compiled
`runsimd.cpp` (the x86 multi-binary launcher) without honoring
`CXXFLAGS`, `CPPFLAGS`, or `LDFLAGS`. The `$(EXE)` link honored
`CXXFLAGS` and `LDFLAGS` but not `CPPFLAGS`. PR #83 has since replaced
the multi-binary scheme with a single binary that builds via the
`single:` target (the default), and that target inherits the same
flag-forwarding behavior.

PR #50 mirrored upstream [bwa-mem2#290](https://github.com/bwa-mem2/bwa-mem2/pull/290):
the compile rules now honor all three variables, and `$(EXE)` link adds
`$(CPPFLAGS)`. This allows downstream packagers (Debian, Bioconda) and
reproducible-build systems to inject hardening flags (`-D_FORTIFY_SOURCE=2`,
`-fstack-protector-strong`, `-Wl,-z,relro`) through the environment without
patching the Makefile. No functional change unless the env vars are set.
Closes [issue #39](https://github.com/fg-labs/bwa-mem3/issues/39).

## Unit-test harness and ARM CI (PR #23)

Historically, PR #23 added a local bash harness (`test/run_unit_tests.sh`) that
built and ran the five C++ unit binaries under `test/` against committed fixtures
in `test/fixtures/`, asserting exit 0 and non-empty diff-able output (those
binaries have since been consolidated into the doctest harness — see the section
above). It also
fixes several pre-existing issues blocking the harness:

- `test/Makefile` defaulted to `icpc` (Intel compiler, not available on GitHub
  runners); changed to `g++` on Linux x86.
- ARM flags are mirrored from the parent Makefile so `cd test && make` builds
  on macOS arm64 and Linux aarch64.
- Three test sources (`smem2_test`, `bwt_seed_strategy_test`, `sa2ref_test`)
  were missing the `fmiSearch->load_index()` call that `fmi_test.cpp` has,
  causing immediate segfaults on run.
- `test/main_banded.cpp` opened `fksw.txt` but never wrote to it; output is
  now written and `main()` returns 0 on success.
- Fixtures are added under `test/fixtures/` covering phiX174, 50 bp test
  reads, BWT seed strategy inputs, SA pairs, and SW pairs.

## CI matrix expansion (PR #24)

PR #24 stacks on PR #23 and expands the GitHub workflow `.github/workflows/ci.yml` from 5 matrix
rows to 7:

| Row | Runner | Arch | Role |
|-----|--------|------|------|
| 1 | ubuntu-latest | sse41 | smoke + unit tests |
| 2 | ubuntu-latest | avx2 | canonical deep tests |
| 3 | ubuntu-latest | avx2 (no mimalloc) | unchanged |
| 4 | ubuntu-24.04-arm | arm64 | unchanged |
| 5 | macos-latest | arm64 | unchanged |
| 6 (new) | ubuntu-latest | multi | runsimd dispatcher smoke |
| 7 (new) | ubuntu-latest | avx2 clang++ | Linux Clang smoke |

The canonical row (row 2) adds: `--bam=6` roundtrip record-count parity,
thread-determinism (`-t1` vs `-t4` sorted diff), unit-test harness, chr22
pipeline parity vs bwa, SE smoke, interleaved smoke, and `--meth` Layers 1–3.

---

## Changes catalog

| Item | bwa-mem3 PR | Upstream PR/issue | Status |
|------|-------------|-------------------|--------|
| doctest framework + Codecov | [#34](https://github.com/fg-labs/bwa-mem3/pull/34) | — | fork-only |
| `PACKAGE_VERSION` from `git describe` | [#52](https://github.com/fg-labs/bwa-mem3/pull/52) | [bwa-mem2#283](https://github.com/bwa-mem2/bwa-mem2/issues/283), [bwa-mem2#284](https://github.com/bwa-mem2/bwa-mem2/pull/284) | fork-only (upstream issue + PR open) |
| PGO target parameterization | [#59](https://github.com/fg-labs/bwa-mem3/pull/59) | — | fork-only |
| `CXXFLAGS`/`CPPFLAGS`/`LDFLAGS` forwarding | [#50](https://github.com/fg-labs/bwa-mem3/pull/50) | [bwa-mem2#290](https://github.com/bwa-mem2/bwa-mem2/pull/290) | fork-only (mirrors open upstream PR) |
| Unit-test harness + ARM CI | [#23](https://github.com/fg-labs/bwa-mem3/pull/23) | — | fork-only |
| CI matrix expansion | [#24](https://github.com/fg-labs/bwa-mem3/pull/24) | — | fork-only |

---

**See also:**
[Developer Guide → Regression test framework](../developer-guide/regression-tests.md) ·
[Developer Guide → Release process](../developer-guide/release.md) ·
[Performance → PGO build](../performance/pgo.md) ·
[Performance improvements](performance.md) ·
[Upstream PR status](upstream-prs.md)
