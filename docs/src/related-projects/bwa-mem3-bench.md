# bwa-mem3-bench

bwa-mem3-bench is a benchmarking suite that measures the alignment performance
of bwa-mem3 against the upstream bwa-mem2 v2.2.1 baseline. It runs on AWS Batch
spot instances across four dataset types — whole-genome sequencing (WGS), whole-exome
sequencing (WES), panel, and bisulfite-sequencing (methylation) — all aligned against
the hg38 reference. The suite covers three CPU microarchitectures: ARM Neon, x86 AVX2,
and x86 AVX-512. Results are collected into a SQLite database for local analysis and
reporting. The project is implemented in Python (orchestration, reporting, and CLI),
Rust (BAM comparison tool), Snakemake (alignment workflow), and AWS CDK (cloud
infrastructure).

## When you'd use it

Use bwa-mem3-bench when you need reproducible, multi-architecture throughput numbers
before committing a bwa-mem3 change to production or before deciding whether to adopt
bwa-mem3 in place of bwa-mem2. It provides a structured "bless baseline, then compare"
workflow: an upstream bwa-mem2 run is blessed once per upstream tag and stored in S3;
subsequent bwa-mem3 runs are measured against that fixed baseline. Running a full
benchmark fires a Snakemake coordinator job on AWS Batch and costs roughly $10 in spot
capacity.

## How it relates to bwa-mem3

bwa-mem3-bench is the authoritative source of benchmark evidence for every performance
claim made in the bwa-mem3 documentation and changelog. When the [Performance
Overview](../performance/overview.md) cites speedup numbers, those numbers come from
bwa-mem3-bench runs collected after the relevant PR was merged. The suite also
validates that bwa-mem3 does not regress relative to bwa-mem2 on any supported
architecture before a new release is tagged.

## Links

- GitHub: <https://github.com/fg-labs/bwa-mem3-bench>
- License: MIT

---

**See also:**
[Performance Overview](../performance/overview.md) ·
[SIMD dispatch matrix](../performance/simd-dispatch.md) ·
[bwa-mem2 (upstream)](bwa-mem2.md) ·
[Release process](../developer-guide/release.md)
