# Testing bwa-mem3

bwa-mem3 tests are organized into three categories — **unit**, **integration**, and **regression** — each with its own build target, CI scope, and runtime cost.

| Category       | Binary / runner                    | Fixtures                                     | CI scope                             |
|----------------|------------------------------------|----------------------------------------------|--------------------------------------|
| **unit**       | `test/bwa_mem3_tests_unit`         | None. All inputs synthetic (scratch files a case generates itself are allowed). | Every matrix row                     |
| **integration**| `test/bwa_mem3_tests_integration`  | Small committed FASTAs / FMI under `test/fixtures/` | SSE4.1, AVX2, ARM64 Linux, macOS ARM |
| **regression** | `test/regression/*.sh`             | Committed phiX in `test/fixtures/`; chr22 downloaded by CI, reads simulated by holodeck; bwa for the parity diffs | Canonical AVX2 row, except `chr22_parity.sh` (every row), `mem_cli_basics.sh`, `mem_io_cli.sh` and `mem_seed_cli.sh` (also Linux ARM64) and the source-only lints (own jobs) |

Design rationale: see the internal test framework design spec (not committed).

## Running tests locally

All commands run from the repo root unless noted.

### Build everything

```bash
make                                      # builds bwa-mem3 binary
make -C test -j$(nproc 2>/dev/null || sysctl -n hw.ncpu)   # builds framework + both doctest binaries + legacy integration binaries
```

The doctest binaries link `libbwa.a` with the decompression libraries the aligner links: htslib's static dependencies, and the two the default FASTQ reader (`src/fast_reader.c`) uses, zlib-ng to inflate plain gzip (the vendored `ext/zlib-ng/build/libz-ng.a`, or the host's `-lz-ng` with `ZLIBNG_USE_VENDORED=0`) and libdeflate to decode BGZF blocks. A unit test can therefore call the reader directly (`bseq_read_fast`, `fr_fastq_next`).

### Run unit tests

```bash
./test/bwa_mem3_tests_unit                # run all unit tests
./test/bwa_mem3_tests_unit --list-test-cases
./test/bwa_mem3_tests_unit --test-case="*kswv*"
./test/bwa_mem3_tests_unit --test-suite="unit/kswv"
./test/bwa_mem3_tests_unit --test-suite-exclude=slow
./test/bwa_mem3_tests_unit --success      # also print passing assertions
./test/bwa_mem3_tests_unit --reporters=junit --out=unit-results.xml
./test/bwa_mem3_tests_unit --help         # full doctest flag list
bash test/run_unit_shards.sh ./test/bwa_mem3_tests_unit   # run as one concurrent shard per CPU (what CI does)
```

### Run unit tests under AddressSanitizer

`ASAN=1` builds `libbwa.a` with AddressSanitizer, and the doctest binaries are then compiled and linked with `-fsanitize=address` too:

```bash
make ASAN=1 test-binaries                 # add the row's arch= (e.g. arch=avx2) to match a CI build
ASAN_OPTIONS=allocator_may_return_null=1:detect_leaks=0 ./test/bwa_mem3_tests_unit
```

A run of the whole unit binary needs both options for now, because of two existing tests: `test_bam_rec_scratch.cpp` requests a `SIZE_MAX` allocation on purpose to check the overflow guard, which AddressSanitizer aborts on unless `allocator_may_return_null=1`; and LeakSanitizer reports 25 KB that `read_memo_prepass` (`src/read_memo.cpp`) allocates for the `read_memo_state` buffers in `test_read_memo.cpp`, whose cases do not release them with `read_memo_state_free`. A run filtered to other files (`-sf`, `-ts`) needs neither option.

### Run integration tests

```bash
./test/bwa_mem3_tests_integration         # same doctest CLI as unit binary
```

Plus the legacy integration binaries (being migrated to the doctest integration binary one-by-one):

```bash
bash test/run_unit_tests.sh               # runs fmi_test, smem2_test, sa2ref_test, bwt_seed_strategy_test, xeonbsw
```

### Run regression checks

Most regression scripts expect the caller to have staged their inputs, named by environment variables — see the comment block at the top of each script, and `test/regression/README.md` for the contract. The source-only lints take neither a binary nor fixtures, so they are the cheapest thing to run:

```bash
# No build, no fixtures: these read src/ and ci.yml directly.
bash test/regression/ndebug_gate_lint.sh
bash test/regression/debug_macro_flag_lint.sh
```

For a script that does need the binary, the phiX-backed ones use a committed fixture, so they need no download and no third-party tool:

```bash
make                                             # builds ./bwa-mem3
BWA_MEM3="$(pwd)/bwa-mem3" \
COMPAT_PHIX_FA="$(pwd)/test/fixtures/phix.fa" \
COMPAT_WORK_DIR=/tmp/compat-bytes \
bash test/regression/compat_byte_identical.sh
```

The chr22 scripts are the expensive end: they want a chr22 reference indexed for both `bwa` and `bwa-mem3` plus a holodeck read simulation, which is why CI caches them. `ci.yml` is the reference for how each script's inputs get staged.

## Running tests in CI

- **ci.yml** runs unit tests on every matrix row, integration tests on the four widened canonical rows, and regression tests on the canonical AVX2 configuration — split between the canonical matrix row (the scripts that need its built tree, such as the Makefile checks) and the `canonical-regressions` job, which builds the same binary alongside the matrix and runs the ones that need only the binary and fixtures — except `chr22_parity.sh`, which runs on every row; the checks that need a differently-configured build of that row (debug macros, `TESTING_BUILD`, ASAN), which run in the `rebuild-variants` job alongside the matrix; and the source-only lints, which run in their own jobs alongside the build matrix (they need no binary, so they do not wait on it).
- **CI depth.** Every change lands through the merge queue, and its `merge_group` run checks the exact commit that lands at full depth. A `pull_request` run, repeated on every push to the branch, is light unless the PR's diff touches a kernel, harness or build path: the SIMD kernels and the headers they include, the rescue harnesses and the kernels' unit tests, `Makefile`, `test/Makefile`, the clang-19 install action, `ci.yml` or the script itself. The list and the diff live in `.github/scripts/ci-depth.sh`, which exports `CI_DEPTH=full|light` to the steps and falls back to `full` whenever it cannot compute the diff. A light run still runs every job and every check except three deep safety nets for the kernels: the rescue kernel sweeps (`rescue_prune_eq`, `rescue_band_harness`) run at about a fifth of their jobs on every row instead of full counts on the canonical, ARM64 and gcc rows; the canonical row reports `SKIP` for the AVX-512BW kswv unit tests; and the MSan padded-lane gate (not a required check) does not run. A light PR that goes green can therefore still fail one of those in the queue, though only through a change outside the paths that guard them.
- JUnit artifacts are uploaded per row (artifacts `unit-results-<row>` and `integration-results-<row>`; the unit artifact holds one `unit-results-<job-index>-shard<k>.xml` per unit-test shard). Download them from a failed run's Actions page to see fine-grained assertion output.
- **proto-neon-kswv.yml** runs `./test/bwa_mem3_tests_unit --test-suite="unit/kswv"` on proto branches.

## Adding a unit test

Unit tests are doctest `TEST_CASE`s in `test/unit/test_*.cpp`. They must:

1. Use only synthetic inputs: read no committed fixture (a test that does is an integration test). A case may write scratch files it generates itself, inside a `bwa_tests::TempDir` (`temp_dir.h`), so long as they are removed when it ends: a corrupt index header, captured stderr, or a real index of a few kilobases of synthetic contigs built and loaded with `bwa_tests::TempIndex` (`temp_index.h`). An index a fixture shares across cases is loaded once and its files deleted at once with `drop_files()`, so the case that builds it leaves none behind.
2. Complete in <100 ms each.
3. Be tagged with `doctest::test_suite("unit/<module>")` (or wrapped in a `TEST_SUITE("unit/<module>") { ... }` block), where `<module>` names the code under test. Reuse an existing module where one fits; add a new one when none does. The modules in use in `test/unit/` are: `alnreg_sort_dedup`, `bam_aux_append_propagation`, `bam_pa_tag_parity`, `bandedswa`, `bandedswa-asym`, `bandedswa-signed8`, `bandedswa-zdrop`, `bandedswa_cpg`, `bns_get_seq`, `bns_pos2rid`, `bns_restore`, `bwa_thp`, `chain`, `chain_cap`, `chain_flat_key`, `chain_flt`, `chain_weight_clamp`, `cigar`, `ext-tier-partition`, `extend_csub`, `extension`, `flt_chained_seeds`, `gen_cigar`, `huge_pages`, `index_load_threads`, `ksort_permutation`, `ksw_global2`, `kswv`, `kswv_cpg`, `kt_for`, `madvise`, `mem_io`, `mem_seed`, `meth-asym-directed`, `meth_model`, `meth_nm_md`, `meth_scoring`, `min_ext_len`, `neon_movemask`, `neon_soa_pack`, `pac_slurp_and_close`, `pair`, `pair64_sort`, `proc_freq`, `pwrite_request_size`, `read_memo`, `rescue_band`, `rescue_rb_index`, `sam_encode`, `seq_reverse`, `smem`, `smem_dedup`, `test_framework`, `ungapped_fastpath`, `util`, `xa_regen` (regenerate from the repo root with `grep -rhoE '(test_suite|TEST_SUITE)\("unit/[^"]+"' test/unit/*.cpp | sed -E 's/.*unit\///; s/"$//' | LC_ALL=C sort -u`). Some older files in `test/unit/` carry no suite tag; doctest runs them in its default suite, so `--test-suite` filters skip them and they are selected by name (`-tc`) instead. Tag new tests. doctest's `test_suite` decorator is overriding (not additive) — chaining `* doctest::test_suite("a") * doctest::test_suite("b")` keeps only the last one — so encode the category and module as a single slash-separated string.
4. Stand alone. CI runs the binary as concurrent shards (`test/run_unit_shards.sh`), each a separate process running a contiguous slice of the cases, so a case cannot rely on an earlier case having run in the same process, and any scratch file must get a unique name (`mkstemp(3)`, or the pid in the name) rather than a fixed path.

Template:

```cpp
// test/unit/test_<module>.cpp
#include "doctest/doctest.h"
#include "<framework-header>.h"

TEST_CASE("short description of behavior"
          * doctest::test_suite("unit/<module>")) {
    // Arrange
    auto input = /* construct synthetic input */;
    // Act
    auto result = function_under_test(input);
    // Assert
    CHECK(result == expected);
}
```

Use `SUBCASE` for parameterized variants:

```cpp
TEST_CASE("behavior across inputs" * doctest::test_suite("unit/<module>")) {
    SUBCASE("small input")   { /* ... */ }
    SUBCASE("empty input")   { /* ... */ }
    SUBCASE("huge input")    { /* ... */ }
}
```

To filter by category or module, use a glob: `--test-suite="unit/*"` selects everything under `unit/`, `--test-suite="*/kswv"` selects every category's `kswv` tests.

Framework helpers available in `test/framework/`:
- `bwa_tests::build_scoring_matrix`, `default_scoring_matrix` (scoring.h)
- `bwa_tests::gen_{random,exact_match,all_mismatch,homopolymer,sub_cluster,with_n_bases}_pair` (seqpair_gen.h)
- `bwa_tests::BatchBuffers` (seqpair_batch.h)
- `bwa_tests::run_scalar_ksw` (ksw_runner.h)
- `bwa_tests::run_kswv_batch` (kswv_runner.h)
- `bwa_tests::kswr_{score,coords,score2}_eq` (kswr_cmp.h)
- `bwa_tests::TempDir` (temp_dir.h)
- `bwa_tests::StderrCapture`, `QuietStderr` (stderr_capture.h)
- `bwa_tests::TempIndex`, `Contig`, `nt4`, `revcomp` (temp_index.h)

The file is picked up automatically by the wildcard rule in `test/Makefile` — no Makefile edits needed when adding a new `test_*.cpp`.

## Adding an integration test

Same as unit tests, but:
- Lives in `test/integration/test_*.cpp` (note the `test_` prefix — plain `.cpp` files in that dir are reserved for the legacy binaries in flight).
- Tagged with `doctest::test_suite("integration/<module>")` (same overriding-decorator caveat as unit).
- May load fixtures from `test/fixtures/`. Fixtures must be committed, ≤100 KB, and come with a documented regeneration procedure.
- Runtime budget per `TEST_CASE`: <10 s.

## Adding a regression test

Write a standalone bash script in `test/regression/<name>.sh`:
- `set -euo pipefail` at the top.
- Document required env-var inputs in the comment block.
- Emit `PASS:` / `FAIL:` lines; exit nonzero on failure.
- Add a row to `test/regression/README.md`.
- Wire it into `.github/workflows/ci.yml`: a script that needs only the built `bwa-mem3` and fixtures goes in the `canonical-regressions` job; one that needs the canonical row's built tree goes in `build-and-test` with the `if: matrix.canonical == true` gate.

A script that runs the binary under each SIMD tier with a batched kswv kernel that the host has sources `test/framework/simd_tiers.sh`, whose `simd_sweep_tiers` lists them (avx2 and, where the runner has it, avx512bw on x86; neon on arm64), as `all_tiers_parity.sh`, `band_cert_identity.sh`, `dash_read_as_n.sh` and `meth_context_identity.sh` do; a script that needs only the tier the binary dispatches to calls its `simd_detect_tier` (`rescue_prune_identity.sh`). `all_tiers_parity.sh` and `band_cert_identity.sh` run only in `canonical-regressions`, so on x86 alone; the `mem` CLI scripts (`mem_cli_basics.sh` and the `mem_*_cli.sh` scripts built on the same helpers) also run on the Linux ARM64 row, gated on `matrix.os == 'ubuntu-24.04-arm'` in `build-and-test`, to cover neon.

Scripts that check `bwa-mem3 mem` output field by field can source the shared helpers in `test/framework/mem_cli.sh` (see `test/regression/mem_cli_basics.sh` for a complete example): `mem_cli_setup_synth` copies and indexes the committed synthetic reference `test/fixtures/mem_synth.fa` (the header of its generator, `test/tools/mem_synth_ref.c`, lists the contigs and planted repeats, so expected coordinates can be worked out by hand), `mem_cli_mem` runs `mem` under one tier (`BWAMEM3_FORCE_TIER`) for each of `mem_cli_tiers` (`mem_cli_mem_status` for a run expected to fail, and `mem_cli_expect_fail` to check its exit and error message in one go; `mem_cli_ref_name` names a reference a script indexes itself), `mem_cli_primary`, `mem_cli_col`, `mem_cli_tag` and `mem_cli_tags` read records back, `mem_cli_count` counts the records an awk condition selects, and `mem_cli_verdict` decides the verdict the script prints (a run with no check, a failed check, or a `mem_cli_seq` slice that failed inside a command substitution is a failure). `test/regression/mem_cli_selftest.sh` checks the helpers themselves. The helpers are written for bash 3.2 and POSIX awk so the scripts also run on a macOS checkout, but CI runs them only on Linux (`canonical-regressions` and the ARM64 row), so that portability is not checked there.

Choose regression over integration when:
- The test shells out to bwa-mem3 (or any other binary) as a process.
- The test diffs against a third-party tool's output (bwa, bwa-meth, samtools).
- The test needs fixtures >100 KB or downloaded at CI time.

## Debugging a failing test

```bash
./test/bwa_mem3_tests_unit --test-case="*kswv*" --success  # verbose
./test/bwa_mem3_tests_unit --test-case="*kswv*" --break    # break into debugger on failure
./test/bwa_mem3_tests_unit --test-case="*foo*" --subcase="bar"  # run single SUBCASE

# Re-run kswv tests with the original phase-0/phase-1 diagnostics:
BWA_TESTS_DEBUG_PHASE0=1 BWA_TESTS_DEBUG_PHASE1=1 \
  ./test/bwa_mem3_tests_unit --test-suite="unit/kswv"
```

For a test using `std::mt19937(seed)`, reproduce locally by editing the seed and running the full binary — there is no persistent "last-failing-seed" record.

## Upgrading doctest

1. Pick a new upstream tag from <https://github.com/doctest/doctest/releases>.
2. Fetch the header:
   ```bash
   curl -sSL -o ext/doctest/doctest.h \
     https://raw.githubusercontent.com/doctest/doctest/<new-tag>/doctest/doctest.h
   ```
3. Compute and record the new SHA256 in `ext/doctest/VERSION`:
   ```bash
   shasum -a 256 ext/doctest/doctest.h
   ```
4. Run the full test suite locally on x86 and ARM.
5. Push to a branch and verify CI passes on every matrix row before merging.
6. Commit with message `build(ext): upgrade doctest to <new-tag>`.

## Framework reference

| Header                           | What                                           |
|----------------------------------|------------------------------------------------|
| `framework/scoring.h`            | `ScoringMatrix` type, `build_scoring_matrix`   |
| `framework/seqpair.h`            | `TestPair` struct                              |
| `framework/seqpair_gen.h`        | Deterministic TestPair generators              |
| `framework/seqpair_batch.h`      | `BatchBuffers` — flat-layout packer            |
| `framework/ksw_runner.h`         | Scalar `run_scalar_ksw`, default gap/xtra      |
| `framework/kswv_runner.h`        | Two-pass `run_kswv_batch`                      |
| `framework/kswr_cmp.h`           | Score / coord / score2 comparators             |
| `framework/ext_result.h`         | `ExtResult`: one extension's six result fields |
| `framework/junit_reporter.h`     | CI matrix-row banner                           |
| `framework/test_main.cpp`        | Shared `main()` (single doctest IMPL TU)       |
| `framework/temp_dir.h`           | `TempDir`: scratch directory, removed on destruction or earlier with `remove()`; `anon_temp_fd`: an already-unlinked scratch file |
| `framework/temp_path.h`          | C99 scratch files and directories under `$TMPDIR` (`bwa_test_mkstemp`, `bwa_test_mkdtemp`, `bwa_test_temp_template`), for the standalone C tests and any test that wants a plain path in a buffer; `temp_dir.h` builds on it |
| `framework/stderr_capture.h`     | `StderrCapture`, `QuietStderr`: scoped stderr redirection |
| `framework/temp_index.h`         | `TempIndex`: a real index of test contigs, built and loaded in a `TempDir` (`drop_files()` deletes the files and keeps the loaded index); `nt4`, `revcomp` |
| `framework/mem_record_fixture.h` | Hand-built inputs for the record writers (`mem_aln2sam`, `mem_aln_to_bam`): a two-contig reference, a read and alignments, with no index |
| `framework/simd_tiers.sh`        | `simd_sweep_tiers`: the SIMD tiers a regression script sweeps; `simd_detect_tier`, `simd_tier_rank`: the tier a binary dispatches to, and tier order |
| `framework/dna.sh`               | `dna_revcomp`: the reverse complement regression scripts build reads with |
| `framework/mem_cli.sh`           | Source-only shell helpers for the `mem` CLI regression scripts |
