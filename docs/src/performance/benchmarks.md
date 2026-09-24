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

## Per-release concordance history

Per-(release, sample) primary-alignment concordance against upstream bwa-mem2 v2.2.1, with supplementary-alignment counts, across released bwa-mem3 versions. Concordance is the minimum vs-baseline value over reps and x86 architectures (deterministic per sample); `supp_query`/`supp_baseline` are total supplementary records emitted by bwa-mem3 and bwa-mem2, and `count_mismatch` is the number of templates whose supplementary count differs. The [divergence catalog](../whats-different/equivalence.md#declared-divergence-catalog) explains what each kind of drift is and its budget.

This table and the divergence catalog are both generated from the benchmark database — do not edit them by hand. Regenerate after a new release is collected with `pixi run python -m bwa_mem3_bench.cli bench docs --releases v0.2.0=<sha>,v0.2.1=<sha>,...` (in the bwa-mem3-bench repo), then replace the content between the `FG-DIVERGENCE-CATALOG` / `FG-RELEASE-TABLE` markers with the emitted `runs/docs/{divergence-catalog.md,release-table.md}` (the `inject_between_markers` helper in `bwa_mem3_bench.report.docs` does exactly this splice).

<!-- FG-RELEASE-TABLE:start -->
| release | sample | concordance_% | supp_query | supp_baseline | count_mismatch |
| --- | --- | --- | --- | --- | --- |
| v0.2.0 | meth-twist-emseq-5M | 98.8852 | 0 | 0 | 0 |
| v0.2.0 | panel-twist-5M | 100.0000 | 186946 | 186946 | 0 |
| v0.2.0 | smoke-1M | 100.0000 | 1455 | 1455 | 0 |
| v0.2.0 | smoke-meth | 98.8573 | 0 | 0 | 0 |
| v0.2.0 | wes-5M | 100.0000 | 5118 | 5118 | 0 |
| v0.2.0 | wgs-5M | 100.0000 | 49686 | 49686 | 0 |
| v0.2.1 | meth-twist-emseq-5M | 98.8852 | 0 | 0 | 0 |
| v0.2.1 | panel-twist-5M | 100.0000 | 186946 | 186946 | 0 |
| v0.2.1 | smoke-1M | 100.0000 | 1455 | 1455 | 0 |
| v0.2.1 | smoke-meth | 98.8573 | 0 | 0 | 0 |
| v0.2.1 | wes-5M | 100.0000 | 5118 | 5118 | 0 |
| v0.2.1 | wgs-5M | 100.0000 | 49686 | 49686 | 0 |
| v0.2.2 | meth-twist-emseq-5M | 98.8773 | 0 | 0 | 0 |
| v0.2.2 | panel-twist-5M | 99.9414 | 187039 | 186946 | 199 |
| v0.2.2 | smoke-1M | 99.9460 | 0 | 0 | 0 |
| v0.2.2 | smoke-meth | 98.8429 | 0 | 0 | 0 |
| v0.2.2 | wes-5M | 99.9996 | 5123 | 5118 | 5 |
| v0.2.2 | wgs-5M | 99.9893 | 49926 | 49686 | 256 |
<!-- FG-RELEASE-TABLE:end -->

## Links

- GitHub: <https://github.com/fg-labs/bwa-mem3-bench>
- License: MIT

---

**See also:**
[Performance Overview](../performance/overview.md) ·
[SIMD dispatch matrix](../performance/simd-dispatch.md) ·
[bwa-mem2 (upstream)](bwa-mem2.md) ·
[Release process](../developer-guide/release.md)
