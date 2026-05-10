# Regression scripts

End-to-end parity and invariant checks. `chr22_parity.sh` runs on every
matrix row; the rest run on the canonical AVX2 row only. Each script:

- is self-contained (set -euo pipefail; explicit env-var contract)
- emits `PASS:` on success and `FAIL:` on failure
- returns nonzero on failure

| Script                       | What it checks                                                        | Origin in ci.yml                                    |
|------------------------------|-----------------------------------------------------------------------|-----------------------------------------------------|
| `chr22_parity.sh`            | bwa vs bwa-mem3 SAM parity on ~50k PE holodeck reads (chr22)          | "Compare chr22 bwa vs bwa-mem3 (parity)"            |
| `thread_determinism.sh`      | `-t 1` == `-t 4` output after sort on chr22                           | "Thread-determinism smoke (chr22, -t 1 vs -t 4)"    |
| `bam_roundtrip.sh`           | `--bam=6` BAM decodes and has same record count as SAM (chr22)        | "--bam=6 roundtrip smoke (chr22)"                   |
| `meth_oracle.sh`             | `--meth` Layers 1–3 match bwa-meth oracle                             | "Run --meth Layers 1-3"                             |

Each script reads its inputs from environment variables — see the comment
block at the top of each file. `.github/workflows/ci.yml` sets those vars
and invokes the scripts.
