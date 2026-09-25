# Benchmarks

Every performance and equivalence figure in this documentation comes from
[bwa-mem3-bench](https://github.com/fg-labs/bwa-mem3-bench), a public benchmarking suite that runs bwa-mem3 alongside
bwa, bwa-mem2 v2.2.1, minibwa and bwameth on AWS instances across the x86 AVX2,
x86 AVX-512 and Arm NEON SIMD tiers. Each release candidate is run through its
full matrix before being blessed. Other pages cite this one instead of linking
to the suite directly.

## Published results

- **[Results for each release](https://github.com/fg-labs/bwa-mem3-bench/blob/main/results/README.md)**: the
  latest release's same-host speed against bwa, bwa-mem2 and minibwa, then per
  release:
  - speed, memory and CPU efficiency on each dataset across six instance types;
  - agreement with bwa-mem2 (and bwameth for methylation), with the previous
    release, between Arm and x86, and between `--fast` and the default preset;
  - accuracy against simulated reads with known origins;
  - thread scaling from 1 to 64 threads.
- **[Methodology](https://github.com/fg-labs/bwa-mem3-bench/blob/main/results/methodology.md)**: how the numbers are
  produced, which host each measurement ran on, and which tags each agreement
  check ignores.

The headline speed table in the
[README](https://github.com/fg-labs/bwa-mem3#performance) comes from the same
same-host runs.

## Datasets

| dataset | description |
| --- | --- |
| `wgs-5M` | Whole genome, 1000 Genomes HG00096, 5M read pairs |
| `wes-5M` | Whole exome, 1000 Genomes HG00100, 5M read pairs |
| `panel-agilent-qxt-5M` | Agilent SureSelect QXT hereditary-cancer panel, 5M read pairs (earlier releases used a Twist panel, `panel-twist-5M`) |
| `hic-1M` | HG002 Hi-C, 1M read pairs |
| `sbx-1M` | HG002 Roche SBX, ~1M single-end reads |
| `meth-twist-emseq-5M` | Twist EM-seq methylation, 5M read pairs |
| `sim-wgs-*`, `sim-meth-*` | Reads simulated from GRCh38 with known origins, for accuracy |

All are aligned to GRCh38 (the hg38 analysis set with decoys and ALTs) with
16 threads. Where to obtain each input is documented in the suite's
[data-setup guide](https://github.com/fg-labs/bwa-mem3-bench/blob/main/docs/data-setup.md).

## Release validation

Each release candidate runs the full matrix and must pass three gates before it
is blessed as the new reference ("golden") for the next release:

1. **vs the baseline aligner**: drift stays within each dataset's budget in the
   [expected-divergence registry](https://github.com/fg-labs/bwa-mem3-bench/blob/main/docs/expected-divergences.yaml)
   (see [Equivalence with bwa-mem2](../whats-different/equivalence.md)). DNA
   datasets are scored on concordance against bwa-mem2 v2.2.1; methylation
   datasets on confident relocation against bwameth (the share of primaries
   either aligner maps at MAPQ ≥ 20 that the two place at a different locus).
2. **vs the previous release**: at least 99.999% concordance on every cell that
   has a previous-release result. An intended alignment change is recorded in
   the release's entry in the
   [release ledger](https://github.com/fg-labs/bwa-mem3-bench/blob/main/docs/release-allowances.yaml),
   which is required to promote the release; the entry does not relax this
   check.
3. **Thread scaling**: pipeline efficiency does not regress.

The suite's [release runbook](https://github.com/fg-labs/bwa-mem3-bench/blob/main/docs/RELEASE.md) is the full
procedure; the [release process](../developer-guide/release.md) page covers the
bwa-mem3 side.

## Earlier one-off measurements

- The ALT-aware `--compat` byte-identity figures on
  [Equivalence with bwa-mem2](../whats-different/equivalence.md) have their
  per-cell breakdown in [bwa-mem3-bench#47](https://github.com/fg-labs/bwa-mem3-bench/pull/47).

## Per-release concordance history

Per-(release, sample) primary-alignment concordance against upstream bwa-mem2 v2.2.1 (bwameth.py for the methylation samples), with supplementary-alignment counts, for the bwa-mem3 releases shown below. Concordance is the minimum vs-baseline value over reps and x86 architectures (deterministic per sample); `supp_query`/`supp_baseline` are total supplementary records emitted by bwa-mem3 and the baseline, and `count_mismatch` is the number of templates whose supplementary count differs. The [divergence catalog](../whats-different/equivalence.md#declared-divergence-catalog) explains what each kind of drift is and its budget. The suite's [methodology](https://github.com/fg-labs/bwa-mem3-bench/blob/main/results/methodology.md#instance-types) lists the instance types behind each architecture and the SIMD path each one runs.

This table and the divergence catalog are generated from the benchmark database; do not edit them by hand. After a new release is collected, regenerate both and replace the content between the `FG-DIVERGENCE-CATALOG` / `FG-RELEASE-TABLE` markers.

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
[Performance Overview](overview.md) ·
[SIMD dispatch matrix](simd-dispatch.md) ·
[Equivalence with bwa-mem2](../whats-different/equivalence.md) ·
[Release process](../developer-guide/release.md)
