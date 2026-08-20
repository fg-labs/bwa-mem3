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
  bam_roundtrip, meth_oracle). `phix_parity.sh` no longer exists:
  [PR #89](https://github.com/fg-labs/bwa-mem3/pull/89) migrated the parity
  checks from dwgsim/phiX174 to holodeck/chr22 and folded it into
  `chr22_parity.sh`.
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

## Testing and CI (PR #23, #24)

The test harness and the GitHub Actions CI matrix (multi-arch builds, the
canonical deep-test row, `--bam` roundtrip and thread-determinism checks, chr22
parity vs bwa, and the `--meth` regressions) are contributor-facing. See the
[Developer Guide → Regression test framework](../developer-guide/regression-tests.md)
for what runs and how to run it locally.

---

**See also:**
[Developer Guide → Regression test framework](../developer-guide/regression-tests.md) ·
[Developer Guide → Release process](../developer-guide/release.md) ·
[Performance → PGO build](../performance/pgo.md) ·
[Performance improvements](performance.md) ·
[Fork changes vs. upstream](../reference/pr-catalog.md)
