# Regression scripts

End-to-end parity and invariant checks. `chr22_parity.sh` and
`version_banner.sh` run on every matrix row; the rest run on the canonical
AVX2 row only, except for the ones wired to a job of their own:
`profile_slice_cpu.sh` runs from `profiling-build` (see below), and the
`ndebug_gate_lint*`, `debug_macro_flag_lint*`, `shell_lint*`,
`regression_coverage_lint*` and `readme_contract_lint*` pairs need no binary at
all and run from `ndebug-gate-lint`, `debug-macro-flag-lint`, `shell-lint`,
`regression-coverage-lint` and `readme-contract-lint` respectively. Each script:

- is self-contained (set -euo pipefail; explicit input contract — env vars for
  every script but the source-only lints, which instead take an optional
  positional directory so their self-tests can aim the same checks at a fixture
  tree, and `meth_oracle.sh`, which reads no input at all)
- emits `PASS:` on success and `FAIL:` on failure
- returns nonzero on failure

| Script                       | What it checks                                                        | Origin in ci.yml                                    |
|------------------------------|-----------------------------------------------------------------------|-----------------------------------------------------|
| `chr22_parity.sh`            | bwa vs bwa-mem3 SAM parity on ~50k PE holodeck reads (chr22)          | "Compare chr22 bwa vs bwa-mem3 (parity)"            |
| `thread_determinism.sh`      | `-t 1` == `-t 4` output after sort on chr22                           | "Thread-determinism smoke (chr22, -t 1 vs -t 4)"    |
| `bam_roundtrip.sh`           | `--bam=6` BAM decodes and has same record count as SAM (chr22)        | "--bam=6 roundtrip smoke (chr22)"                   |
| `short_read_smoke.sh`        | ASAN SE on 25-50 bp variable-length dense-chr22 reads (PR #100 fix)   | "Short-read SE smoke (chr22, ASAN, dense+variable-length)" |
| `supp_rep_hard_cap.sh`       | `--supp-rep-hard-cap` forces MAPQ=0 on repetitive-seed supps (#101)   | "--supp-rep-hard-cap repetitive-seed regression"    |
| `compat_byte_identical.sh`   | `--compat=bwa-mem2` suppresses only MQ/HN/@HD; rest byte-identical    | "--compat byte-identical regression (phix)"         |
| `header_parity.sh`           | `AH:*` on generated @SQ (#281); `--compat` skips the .hdr/.dict sidecar | "header parity regression (AH:*, --compat @HD/@SQ)" |
| `default_hd_parity.sh`       | default `@HD` byte-identical across SAM/`--bam`/`--meth` (#288)        | "default @HD parity across output paths"            |
| `meth_sam_output.sh`         | `--meth` emits SAM text by default and BAM under `--bam`, same records | "--meth output container follows --bam"             |
| `meth_oracle.sh`             | `--meth` Layer 1 (valid BAM emission) via the harness under `test/meth/`; Layers 2–3 retired in D3 | "Run --meth Layer 1"                                |
| `cohort_ramp_validation.sh`  | `--cohort-ramp-first`/`-ratio` reject malformed values; env warns and falls back | "Cohort ramp values are validated (flag and environment alike)" |
| `rescue_kmer_options.sh`     | `--rescue-kmer`/`--rescue-band` reject malformed and out-of-range values; the `=` form is required | "--rescue-kmer/--rescue-band reject malformed values" |
| `profile_slice_cpu.sh`       | `--profile` accounts for a partial cohort slice's compute CPU (needs `STAGE_PROF=1`) | "Partial cohort slices report their compute CPU"    |
| `ndebug_gate_lint.sh`        | no `#if`/`#ifdef`/`#ifndef`/`#elif` NDEBUG gates in `src/` — nothing here defines NDEBUG, so they never compile out | "No NDEBUG preprocessor gates in src/"              |
| `ndebug_gate_lint_selftest.sh` | the lint above still flags real gates, so its `PASS` means something | "NDEBUG gate lint still detects gates"              |
| `debug_macro_flag_lint.sh`   | the opt-in macro build's `-D` list and the `BWA_MEM3_DEBUG_*` macros in `src/` still name each other | "Opt-in macro -D list matches the macros in src/"   |
| `debug_macro_flag_lint_selftest.sh` | the lint above still detects a drifted list, so its `PASS` means something | "Macro list lint still detects drift"               |
| `shell_lint.sh`              | every tracked `*.sh` is shellcheck-clean and shfmt-formatted           | "Tracked shell scripts are shellcheck-clean and shfmt-formatted" |
| `shell_lint_selftest.sh`     | the lint above still rejects bad scripts, so its `PASS` means something | "Shell lint still detects bad scripts"             |
| `regression_coverage_lint.sh` | every script in this directory is named by a CI workflow, or by a Makefile target CI invokes — not just by `make test` | "Every regression script is run by CI"              |
| `regression_coverage_lint_selftest.sh` | the lint above still detects an unrun script, so its `PASS` means something | "Coverage lint still detects an unrun script"       |
| `host_floor_enforce.sh`      | below-floor hosts get exit 2 + a readable error, not a SIGILL (needs `TESTING_BUILD=1`) | "SIMD floor enforcement (TESTING_BUILD, injected below-floor tier)" |
| `version_banner.sh`          | `bwa-mem3 version` prints the SIMD floor and runtime tier lines        | "Version banner regression"                         |
| `readme_contract_lint.sh`    | this README names no script that was deleted, and its source-only-lint block lists exactly the scripts that read no environment | "README still describes the regression scripts"     |
| `readme_contract_lint_selftest.sh` | the lint above still detects a stale README, so its `PASS` means something | "README lint still detects drift"                  |

The table is a reading guide, not an inventory — `ls test/regression/*.sh` is
the authoritative list, and `ci.yml` is where each one is actually wired up.

Every script but the ones named in the block below reads its inputs from
environment variables — see the comment block at the top of each file.
`.github/workflows/ci.yml` sets the required vars and invokes the scripts.

<!-- source-only lints: begin -->
These scripts read no environment at all.

The source-only lints each take the directory to work on as an optional
positional argument, so that its self-test can aim the same checks at a fixture
tree; each self-test takes no input, since it builds the trees it aims its lint
at.

| Lint | Positional argument | Self-test |
|------|---------------------|-----------|
| `ndebug_gate_lint.sh`         | directory to scan (default `src/`)                  | `ndebug_gate_lint_selftest.sh`         |
| `debug_macro_flag_lint.sh`    | repository root to check (default: this repository) | `debug_macro_flag_lint_selftest.sh`    |
| `regression_coverage_lint.sh` | repository root to check (default: this repository) | `regression_coverage_lint_selftest.sh` |
| `readme_contract_lint.sh`     | repository root to check (default: this repository) | `readme_contract_lint_selftest.sh`     |

`meth_oracle.sh` is the one env-free script that is not a lint, and it takes no
argument either. It wraps the `--meth` harness under `test/meth/`, whose inputs
are the gitignored fixtures CI copies in there, plus an optional `SAMTOOLS`
naming a samtools off `PATH`; the wrapper's whole job is to invoke that harness
and mark the result, so it passes the environment through rather than reading
or pinning any of it. "Reads no environment" is a claim about each script
itself — which is also how `readme_contract_lint.sh` decides it — not about
everything it may go on to run.
<!-- source-only lints: end -->

`readme_contract_lint.sh` checks that block against the scripts themselves, so
a script that stops reading the environment without gaining a mention here
fails CI rather than going unnoticed.

One exception to "any binary will do": `profile_slice_cpu.sh` asserts on
`--profile` output, which a default build compiles out entirely, so it needs a
binary from `make STAGE_PROF=1`. It fails loudly rather than skipping if handed
one without `--profile`, and CI runs it from the `profiling-build` job.
