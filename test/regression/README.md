# Regression scripts

End-to-end parity and invariant checks. `chr22_parity.sh` and
`version_banner.sh` run on every matrix row; the rest run on the canonical
`Linux x86_64 AVX2 (mimalloc)` row only, except for the ones wired to a job of
their own: `profile_slice_cpu.sh` runs from `profiling-build` (see below), and
the `ndebug_gate_lint*`, `debug_macro_flag_lint*`, `shell_lint*`,
`regression_coverage_lint*` and `readme_contract_lint*` pairs need no binary at
all and run from `ndebug-gate-lint`, `debug-macro-flag-lint`, `shell-lint`,
`regression-coverage-lint` and `readme-contract-lint` respectively. Each script:

- is self-contained (set -euo pipefail; explicit input contract — env vars for
  every script but the source-only lints, which instead take an optional
  positional directory so their self-tests can aim the same checks at a fixture
  tree, and `meth_oracle.sh`, which reads no input at all)
- emits `PASS:` on success and `FAIL:` on failure, and `SKIP:` where a check
  cannot run (a missing tool, a host with nothing to compare) — a skipped check
  is reported as skipped, never as a pass
- returns nonzero on failure

Where a script runs is also what backs its claim: unless a row says otherwise,
read that row's byte-identity or parity claim as holding on the canonical row —
Linux x86_64, AVX2 tier, mimalloc. Two rows say otherwise, and neither is one of
the job-of-its-own cases above. `chr22_parity.sh` runs on every matrix row, so
its parity claim is not AVX2-specific. `all_tiers_parity.sh` does run on the
canonical row, but it drives the binary through each usable tier with
`BWAMEM3_FORCE_TIER`, so the tier is the one thing its claim deliberately does
not hold fixed. (`version_banner.sh` runs everywhere too, but its row asserts no
parity or byte-identity, so nothing here scopes it.)

| Script                       | What it checks                                                        | Origin in ci.yml                                    |
|------------------------------|-----------------------------------------------------------------------|-----------------------------------------------------|
| `chr22_parity.sh`            | bwa vs bwa-mem3 SAM parity on ~50k PE holodeck reads (chr22)          | "Compare chr22 bwa vs bwa-mem3 (parity)"            |
| `thread_determinism.sh`      | `-t 1` == `-t 4` output after sort on chr22                           | "Thread-determinism smoke (chr22, -t 1 vs -t 4)"    |
| `all_tiers_parity.sh`        | one binary, one host: byte-identical SAM on the staged `PARITY_FA`/`PARITY_R1`/`PARITY_R2` workload (the same chr22 holodeck PE reads as `chr22_parity.sh`) whichever kernel `BWAMEM3_FORCE_TIER` selects. Real coverage needs an x86 host where both avx2 and avx512bw are sweepable; anywhere else — arm64, or x86 without AVX-512BW — there is one tier and it reports `SKIP:` | "SIMD tier parity (avx2 vs avx512bw, chr22)"        |
| `bam_roundtrip.sh`           | `--bam=6` BAM decodes and has same record count as SAM (chr22)        | "--bam=6 roundtrip smoke (chr22)"                   |
| `short_read_smoke.sh`        | ASAN SE on 25-50 bp variable-length dense-chr22 reads (PR #100 fix)   | "Short-read SE smoke (chr22, ASAN, dense+variable-length)" |
| `supp_rep_hard_cap.sh`       | `--supp-rep-hard-cap` forces MAPQ=0 on repetitive-seed supps (#101)   | "--supp-rep-hard-cap repetitive-seed regression"    |
| `compat_byte_identical.sh`   | one binary, one host: `--compat=bwa-mem2` suppresses only `MQ`/`HN`/`@HD`, and the rest of the SAM is byte-identical to the same binary's default run on PE reads sliced from phix. bwa-mem2 itself is never run — this is a default-vs-`--compat` invariant, not a cross-aligner one. `@PG` is excluded from both sides: its `CL:` records the argv, which necessarily differs between the two runs | "--compat byte-identical regression (phix)"         |
| `hic_alias_identical.sh`     | one binary, one host: `--hic` is byte-identical to `-5SP` on PE reads sliced from phix (#368). Two guards keep it from passing vacuously — `--hic` output must differ from a plain run, and the whole-file `FLAG` `0x2` counts must show the plain run with at least one proper-pair record while the `--hic` run has exactly zero (a total count over all records, not restricted to records that overlap the plain run). `@PG` is excluded from both sides: its `CL:` records the argv, which differs by construction. Env: `BWA_MEM3`, `HIC_PHIX_FA` (must point at the committed `test/fixtures/phix.fa`), `HIC_WORK_DIR`; no staged inputs — the FASTQs are generated from phix | "--hic alias regression (byte-identical to -5SP)"    |
| `header_parity.sh`           | `AH:*` on generated `@SQ` (#281); `--compat` skips the .hdr/.dict sidecar | "header parity regression (AH:*, --compat @HD/@SQ)" |
| `proper_pair_alt.sh`         | `--proper-pair-from-emitted` derives `FLAG` `0x2` from the emitted `a[which]` rather than the default `a[0]` (#17, #362), on a generated 2-contig ALT fixture whose one divergent pair emits on the ALT contig over a sub-`T` primary hit. Asserts the bit flips, that `FLAG` is the only field that differs and differs by exactly bit `0x2`, and that the option is inert on the same reference without a `.alt` sidecar. Generates its own fixture — no staged inputs, but `PROPER_PAIR_ALT_WORK_DIR` must name a directory of its own, since the script deletes the fixture files it owns there on entry | "proper-pair ALT derivation regression (--proper-pair-from-emitted)" |
| `default_hd_parity.sh`       | one binary, one host: default `@HD` byte-identical across SAM/`--bam`/`--meth` on PE reads sliced from phix (#288) | "default @HD parity across output paths"            |
| `alt_pa_parity.sh`           | one binary, one host: `pa:f:` decodes to the same float32 in SAM text and `--bam`, and carries no more than three decimals, on a generated ALT-aware fixture whose reads come from the ALT contig — the only non-`--meth` fixture in this directory that makes `alt_sc` positive, so the only one that emits `pa` at all (#365). Default path only (`mem_aln2sam` vs `mem_aln_to_bam`); the `--meth` writer is covered by `meth_alt_pa_parity.sh`. Compares the decoded float, not the text `samtools view` prints, which renders an `f` field with `%g` | "pa:f: SAM vs --bam parity (ALT-aware fixture)"     |
| `meth_alt_pa_parity.sh`      | one binary, one host: `--meth` renders `pa:f:` from the same shared definition the SAM-text writer uses — every value is the float32 of a three-decimal rendering, secondary records carry none, and `SA:Z` precedes it (#365). Needed separately because `meth_mem_aln_to_bam` is a third, independent writer and its fixture has to be bisulfite-converted as well as ALT-aware | "pa:f: --meth parity (ALT-aware bisulfite fixture)" |
| `meth_sam_output.sh`         | `--meth` emits SAM text by default and BAM under `--bam`, same records | "--meth output container follows --bam"             |
| `meth_collapsed_scoring.sh`  | `--meth-scoring collapsed` also frees the conversion mirror cell, so a ref-T→read-C substitution scores `a+b` above `genomic` and reports NM=0 against its NM=1 | "--meth whole-aligner regressions (D3)"             |
| `meth_oracle.sh`             | `--meth` Layer 1 (valid BAM emission) via the harness under `test/meth/`; Layers 2–3 retired in D3 | "Run --meth Layer 1"                                |
| `cohort_ramp_validation.sh`  | `--cohort-ramp-first`/`-ratio` reject malformed values; env warns and falls back | "Cohort ramp values are validated (flag and environment alike)" |
| `rescue_kmer_options.sh`     | `--rescue-kmer`/`--rescue-band` reject malformed and out-of-range values; the `=` form is required | "--rescue-kmer/--rescue-band reject malformed values" |
| `rescue_skip_options.sh`     | `--rescue-skip` requires `--rescue-kmer` (order-independent), takes no argument, and is NOT enabled by `--fast` | "--rescue-skip requires --rescue-kmer and composes with --fast" |
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
| `readme_contract_lint.sh`    | this README names no script that was deleted, its source-only-lint block lists exactly the scripts that read no environment, every row's `Origin in ci.yml` names a step a workflow defines, and every script can emit the `PASS:`/`FAIL:` markers above | "README still describes the regression scripts"     |
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
