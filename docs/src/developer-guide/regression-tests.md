# Regression test framework

bwa-mem3 has three categories of tests — **unit**, **integration**, and **regression** — plus a separate **benchmark harness** in `bench/`. Understanding the distinction helps you choose where to add a new test and what to expect from CI.

## Test categories

| Category | Binary / runner | Fixtures | CI scope |
|---|---|---|---|
| unit | `test/bwa_mem3_tests_unit` | None; all inputs synthetic | Every matrix row |
| integration | `test/bwa_mem3_tests_integration` | Small committed FASTAs / FMI in `test/fixtures/` | SSE4.1, AVX2, ARM64 Linux, macOS ARM |
| regression | `test/regression/*.sh` | Downloaded references (phiX, chr22) + bwa + dwgsim | Canonical AVX2 row only |

**Unit tests** must use only synthetic inputs generated programmatically and complete in under 100 ms each. They exercise individual kernels in isolation: kswv scoring, banded Smith-Waterman, KSW, FM-index operations, SMEM extraction, BAM encoding, and pair handling.

**Integration tests** may load small committed fixtures from `test/fixtures/` and have a per-test budget of 10 seconds. They exercise cross-component paths: index loading, SMEM-to-alignment pipelines, and output format validation.

**Regression tests** are standalone bash scripts that shell out to the `bwa-mem3` binary, may diff against third-party tool output (bwa, bwa-meth, samtools), and require fixtures that are either committed to the fixtures directory or downloaded by CI at run time.

## Running tests locally

```bash
# Build the aligner and test binaries
make
make -C test -j$(nproc 2>/dev/null || sysctl -n hw.ncpu)

# Run all unit tests
./test/bwa_mem3_tests_unit

# Run all integration tests
./test/bwa_mem3_tests_integration

# Run a specific test case or suite
./test/bwa_mem3_tests_unit --test-case="*kswv*"
./test/bwa_mem3_tests_unit --test-suite="unit/kswv"
./test/bwa_mem3_tests_unit --test-suite-exclude=slow

# Verbose output (also print passing assertions)
./test/bwa_mem3_tests_unit --success
```

The `make test` target is a convenience shortcut that builds and runs the unit and integration binaries plus the two legacy standalone regression tests (`kswv_nrow_zero_test` and `shm_section_find_test`):

```bash
make test
```

### Running a regression test locally

Regression scripts expect certain environment variables to point at fixtures. The phiX parity test requires `dwgsim`:

```bash
mkdir -p /tmp/ci-test && cd /tmp/ci-test
curl -sL "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/819/615/GCF_000819615.1_ViralProj14015/GCF_000819615.1_ViralProj14015_genomic.fna.gz" | gunzip > phix174.fa
dwgsim -z 42 -N 500 -1 150 -2 150 -r 0.001 -S 2 phix174.fa reads
cd -
BWA_MEM2="$(pwd)/bwa-mem3" CI_TEST_DIR=/tmp/ci-test bash test/regression/phix_parity.sh
```

## Test framework

The unit and integration binaries are built on [doctest](https://github.com/doctest/doctest), a single-header C++ test framework. Tests are discovered by file glob: any `test/unit/test_*.cpp` file is compiled into the unit binary; any `test/integration/test_*.cpp` file is compiled into the integration binary. No Makefile edit is needed when adding a new `test_*.cpp`.

### Test organisation

Tag each `TEST_CASE` with `doctest::test_suite("category/module")`:

```cpp
TEST_CASE("nrow==0 batch does not store out of bounds"
          * doctest::test_suite("unit/kswv")) {
    // ...
}
```

The `test_suite` decorator is overriding (not additive). Encode the category (`unit` or `integration`) and module (`kswv`, `bandedsw`, `ksw`, `fmindex`, `smem`, `bam`, `pair`, `cigar`, `util`) as a single slash-separated string.

### Framework helpers

The `test/framework/` directory provides helpers shared across test files:

| Header | Provides |
|---|---|
| `scoring.h` | `ScoringMatrix`, `build_scoring_matrix`, `default_scoring_matrix` |
| `seqpair.h` | `TestPair` struct |
| `seqpair_gen.h` | Deterministic pair generators: random, exact-match, all-mismatch, homopolymer, sub-cluster, N-bases |
| `seqpair_batch.h` | `BatchBuffers` — flat-layout packer for kswv batch input |
| `ksw_runner.h` | `run_scalar_ksw`, default gap/extra parameters |
| `kswv_runner.h` | Two-pass `run_kswv_batch` |
| `kswr_cmp.h` | Score / coordinate / score2 comparators |
| `junit_reporter.h` | CI matrix-row banner and JUnit XML output |

### Debugging a failing test

```bash
# Break into debugger at the first failing assertion
./test/bwa_mem3_tests_unit --test-case="*kswv*" --break

# Run a single SUBCASE
./test/bwa_mem3_tests_unit --test-case="*foo*" --subcase="bar"

# Enable per-phase diagnostics for kswv tests
BWA_TESTS_DEBUG_PHASE0=1 BWA_TESTS_DEBUG_PHASE1=1 \
  ./test/bwa_mem3_tests_unit --test-suite="unit/kswv"
```

JUnit artifacts are uploaded per CI matrix row (`unit-results-<name>.xml`, `integration-results-<name>.xml`) and available on the Actions run page.

> **Tip — Use ASAN for memory bugs**
>
> Build with `make ASAN=1 test` to catch out-of-bounds writes in vectorised kernels.
> The kswv_nrow_zero_test specifically exercises the nrow==0 path that triggered a
> pre-allocation store bug; ASAN reports this immediately rather than at a later
> allocator operation.

## Standalone regression tests

Three standalone regression tests live outside the doctest harness because they predated it. The two binaries are built and run by `make test`; the third is script-driven:

- `kswv_nrow_zero_test` — binary; exercises the all-len1==0 batch path in every SIMD kswv variant. Catches the nrow==0 rowMax store overrun from issue #38 / upstream bwa-mem2 PR #289.
- `shm_section_find_test` — binary; exercises the shared-memory index section-find logic.
- `shm_pack_round_trip_test` — script-driven, invoked via `test/shm_pack_round_trip_test.sh`, which builds the phiX index first.

Additional integration shell scripts in `test/`:

| Script | What it tests |
|---|---|
| `pg_cl_escape_test.sh` | `@PG CL:` tab/newline escape in SAM headers |
| `mimalloc_loaded_test.sh` | mimalloc override is active when `USE_MIMALLOC=1` |
| `shm_round_trip_test.sh` | `bwa-mem3 shm` load / list / drop cycle |
| `shm_meth_test.sh` | `--meth` index compatibility with `shm` |
| `help_prescan_test.sh` | `--help` prints without running alignment |
| `proc_freq_calibration_test.sh` | startup does not sleep to calibrate `proc_freq`, and the calibrated tick rate is plausible |
| `bam_compress_warn_test.sh` | the `--bam=N` single-writer-thread warning fires once, for the resolved level only |
| `libsais_*.sh` | libsais index correctness vs. BWA / determinism |

## Benchmark harness (`bench/`)

`bench/` is a separate performance measurement harness used during development to gate performance PRs. It is not part of the CI test suite.

```bash
cp bench/config.env.example bench/config.env
# Edit config.env to point at your index, reads, and binary paths
bench/run.sh baseline         # N trials; appends to bench/results.csv
bench/run.sh candidate        # N trials on the candidate binary
bench/compare.sh baseline candidate  # wall-clock / RSS / md5 delta report
```

Each run records: tag, host, architecture, binary path, thread count, trial index, wall-clock seconds, max RSS (KB), and a golden md5 (single-threaded, `@PG`-stripped SAM). The md5 verifies byte-identical output across builds; wall-clock is the primary performance metric.

---

**See also:**
[Building from source](building.md) ·
[SIMD dispatch architecture](simd-dispatch.md) ·
[Contributing](contributing.md) ·
[What's Different → Correctness fixes](../whats-different/correctness.md) ·
[What's Different → Build & infrastructure](../whats-different/build-infra.md)
