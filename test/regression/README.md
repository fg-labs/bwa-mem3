# Regression scripts

End-to-end parity and invariant checks. `chr22_parity.sh` runs on every
matrix row; the rest run on the canonical AVX2 row only, except
`profile_slice_cpu.sh`, which runs from the `profiling-build` job (see below),
and the two `ndebug_gate_lint*` scripts, which need no binary at all and run
from the `ndebug-gate-lint` job. Each script:

- is self-contained (set -euo pipefail; explicit env-var contract)
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
| `meth_oracle.sh`             | `--meth` Layers 1–3 match bwa-meth oracle                             | "Run --meth Layers 1-3"                             |
| `cohort_ramp_validation.sh`  | `--cohort-ramp-first`/`-ratio` reject malformed values; env warns and falls back | "Cohort ramp values are validated (flag and environment alike)" |
| `rescue_kmer_options.sh`     | `--rescue-kmer`/`--rescue-band` reject malformed and out-of-range values; the `=` form is required | "--rescue-kmer/--rescue-band reject malformed values" |
| `profile_slice_cpu.sh`       | `--profile` accounts for a partial cohort slice's compute CPU (needs `STAGE_PROF=1`) | "Partial cohort slices report their compute CPU"    |
| `ndebug_gate_lint.sh`        | no `#if`/`#ifdef`/`#ifndef`/`#elif` NDEBUG gates in `src/` — nothing here defines NDEBUG, so they never compile out | "No NDEBUG preprocessor gates in src/"              |
| `ndebug_gate_lint_selftest.sh` | the lint above still flags real gates, so its `PASS` means something | "NDEBUG gate lint still detects gates"              |

The table is a reading guide, not an inventory — `ls test/regression/*.sh` is
the authoritative list, and `ci.yml` is where each one is actually wired up.

Each script reads its inputs from environment variables — see the comment
block at the top of each file. `.github/workflows/ci.yml` sets those vars
and invokes the scripts.

One exception to "any binary will do": `profile_slice_cpu.sh` asserts on
`--profile` output, which a default build compiles out entirely, so it needs a
binary from `make STAGE_PROF=1`. It fails loudly rather than skipping if handed
one without `--profile`, and CI runs it from the `profiling-build` job.
