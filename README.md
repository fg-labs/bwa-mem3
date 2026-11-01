# bwa-mem3

[![CI](https://github.com/fg-labs/bwa-mem3/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/fg-labs/bwa-mem3/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/fg-labs/bwa-mem3/branch/main/graph/badge.svg)](https://codecov.io/gh/fg-labs/bwa-mem3)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/bwa-mem3.svg?label=bioconda)](https://anaconda.org/bioconda/bwa-mem3)
[![Documentation](https://img.shields.io/readthedocs/bwa-mem3?label=docs)](https://bwa-mem3.readthedocs.io)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/fg-labs/bwa-mem3/blob/main/LICENSE)

bwa-mem3 is a short-read aligner derived from [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2),
carrying correctness fixes, performance improvements, and new features (methylation alignment,
shared-memory index, mimalloc allocator) maintained by [Fulcrum Genomics](https://fulcrumgenomics.com).

## Performance

Wall-clock speedup of the current release (v0.15.1\*) against `bwa` 0.7.19, `bwa-mem2` v2.2.1, and `minibwa`, on the `wgs-5M` sample. Cells are `stock / --fast`.

| arch | wall_s | vs bwa | vs bwa-mem2 | vs minibwa |
|---|---:|---:|---:|---:|
| ARM | 32.38 / 15.07 | 7.22x / 15.52x | — | 1.26x / 2.70x |
| x86 | 19.62 / 11.56 | 9.97x / 16.92x | 4.38x / 7.43x | 1.55x / 2.63x |

\* v0.15.1, like v0.13.0 through v0.15.0, is measured with a **denser suffix-array index** — stride 2 instead of the stock stride 8 — built with `bwa-mem3 re-sa -u 1` (or staged with `bwa-mem3 shm -u 1`), both new in v0.13.0. Output is byte-identical by construction (each added suffix-array sample is computed from the on-disk ones), and measured so: a full-content `fgumi compare bams` of `wgs-5M` (10,030,558 records) aligned on the stride-2 index against the stock stride-8 index — v0.13.0 against the output-identical v0.12.0 — found 0 differences on both benchmark hosts, Graviton4 m8g (NEON) and AMD m8a (AVX-512), and the change that added it ([#510](https://github.com/fg-labs/bwa-mem3/pull/510)) found 0 differences on 5M-pair hg38 WGS and WES slices. It costs ~12 GB more memory for the human genome. v0.13.0 through v0.15.1 are all measured on it in the tables below, so v0.15.1's `vs prev` is a like-for-like comparison. Every older release and every comparator uses the stock index.

> [!TIP]
> **📈 Full release-history table** — every bwa-mem3 release since v0.2.1, full methodology, and version pins.
>
> <details>
> <summary><strong>Click to expand</strong></summary>
>
> **Graviton4 (m8g, arm64/NEON)**
>
> | release | wall_s | vs bwa | vs bwa-mem2 | vs minibwa | vs prev |
> |---|---:|---:|---:|---:|---:|
> | bwa | 233.94 | 1.00x | — | 0.17x | — |
> | bwa-mem2 | — | — | — | — | — |
> | minibwa | 40.70 | 5.75x | — | 1.00x | — |
> | v0.2.1 | 121.08 | 1.93x | — | 0.34x | — |
> | v0.2.2 | 118.12 | 1.98x | — | 0.34x | 1.025x |
> | v0.3.0 | 103.20 | 2.27x | — | 0.39x | 1.145x |
> | v0.4.0 | 104.25 | 2.24x | — | 0.39x | 0.990x |
> | v0.5.0 | 104.45 / 37.16 | 2.24x / 6.30x | — | 0.39x / 1.10x | 0.998x |
> | v0.6.0 | 97.21 / 36.49 | 2.41x / 6.41x | — | 0.42x / 1.12x | 1.074x / 1.018x |
> | v0.7.0 | 91.55 / 39.69 | 2.56x / 5.89x | — | 0.44x / 1.03x | 1.062x / 0.919x |
> | v0.8.0 | 76.32 / 28.46 | 3.07x / 8.22x | — | 0.53x / 1.43x | 1.200x / 1.395x |
> | v0.9.0 | 76.43 / 28.47 | 3.06x / 8.22x | — | 0.53x / 1.43x | 0.999x / 1.000x |
> | v0.10.0 | 71.05 / 27.90 | 3.29x / 8.38x | — | 0.57x / 1.46x | 1.076x / 1.020x |
> | v0.11.0 | 63.93 / 26.31 | 3.66x / 8.89x | — | 0.64x / 1.55x | 1.111x / 1.060x |
> | v0.12.0 | 58.22 / 19.24 | 4.02x / 12.16x | — | 0.70x / 2.12x | 1.098x / 1.367x |
> | v0.13.0\* | 55.08 / 18.72 | 4.25x / 12.50x | — | 0.74x / 2.17x | 1.057x / 1.028x |
> | v0.14.0\* | 41.72 / 17.60 | 5.61x / 13.29x | — | 0.98x / 2.31x | 1.320x / 1.064x |
> | v0.14.1\* | 35.55 / 16.59 | 6.58x / 14.10x | — | 1.14x / 2.45x | 1.174x / 1.061x |
> | v0.15.0\* | 34.85 / 16.09 | 6.71x / 14.54x | — | 1.17x / 2.53x | 1.020x / 1.031x |
> | **v0.15.1**\* | **32.38 / 15.07** | **7.22x / 15.52x** | **—** | **1.26x / 2.70x** | **1.076x / 1.068x** |
>
> **AMD (m8a, x86-64/AVX-512)**
>
> | release | wall_s | vs bwa | vs bwa-mem2 | vs minibwa | vs prev |
> |---|---:|---:|---:|---:|---:|
> | bwa | 195.57 | 1.00x | 0.44x | 0.16x | — |
> | bwa-mem2 | 85.92 | 2.28x | 1.00x | 0.35x | — |
> | minibwa | 30.42 | 6.43x | 2.82x | 1.00x | — |
> | v0.2.1 | 56.39 | 3.47x | 1.52x | 0.54x | 1.524x |
> | v0.2.2 | 56.78 | 3.44x | 1.51x | 0.54x | 0.993x |
> | v0.3.0 | 55.35 | 3.53x | 1.55x | 0.55x | 1.026x |
> | v0.4.0 | 52.73 | 3.71x | 1.63x | 0.58x | 1.050x |
> | v0.5.0 | 52.55 / 22.81 | 3.72x / 8.57x | 1.64x / 3.77x | 0.58x / 1.33x | 1.003x |
> | v0.6.0 | 50.37 / 21.84 | 3.88x / 8.95x | 1.71x / 3.93x | 0.60x / 1.39x | 1.043x / 1.044x |
> | v0.7.0 | 48.90 / 22.29 | 4.00x / 8.77x | 1.76x / 3.85x | 0.62x / 1.36x | 1.030x / 0.980x |
> | v0.8.0 | 43.43 / 18.71 | 4.50x / 10.45x | 1.98x / 4.59x | 0.70x / 1.63x | 1.126x / 1.191x |
> | v0.9.0 | 43.42 / 18.21 | 4.50x / 10.74x | 1.98x / 4.72x | 0.70x / 1.67x | 1.000x / 1.027x |
> | v0.10.0 | 43.42 / 18.70 | 4.50x / 10.46x | 1.98x / 4.59x | 0.70x / 1.63x | 1.000x / 0.974x |
> | v0.11.0 | 41.83 / 17.67 | 4.68x / 11.07x | 2.05x / 4.86x | 0.73x / 1.72x | 1.038x / 1.058x |
> | v0.12.0 | 33.18 / 14.65 | 5.89x / 13.35x | 2.59x / 5.86x | 0.92x / 2.08x | 1.261x / 1.206x |
> | v0.13.0\* | 30.05 / 14.18 | 6.51x / 13.79x | 2.86x / 6.06x | 1.01x / 2.15x | 1.104x / 1.033x |
> | v0.14.0\* | 24.15 / 13.07 | 8.10x / 14.96x | 3.56x / 6.57x | 1.26x / 2.33x | 1.244x / 1.085x |
> | v0.14.1\* | 24.04 / 13.07 | 8.14x / 14.96x | 3.57x / 6.57x | 1.27x / 2.33x | 1.005x / 1.000x |
> | v0.15.0\* | 24.34 / 12.07 | 8.03x / 16.20x | 3.53x / 7.12x | 1.25x / 2.52x | 0.988x‡ / 1.083x |
> | **v0.15.1**\* | **19.62 / 11.56** | **9.97x / 16.92x** | **4.38x / 7.43x** | **1.55x / 2.63x** | **1.241x / 1.044x** |
>
> `vs prev` is the release-over-release speedup (`prev_wall / this_wall`, `>1` = faster) vs the previous release on this same host, `stock / --fast`. ‡ On x86 this table run put v0.15.0 at 0.988x of v0.14.1 (24.34 vs 24.04 s); v0.15.0's own release run measured 1.056x (22.83 vs 24.11 s). v0.15.1's speedup over v0.15.0 on the same host holds on alignment-only (`PROCESS()`) time, which excludes index loading and so is unaffected by a cold first rep: 1.23x (22.81 vs 18.48 s). The first release's predecessor is upstream `bwa-mem2` — bwa-mem3 is its successor — so v0.2.1's `vs prev` is its speedup over bwa-mem2 (blank on ARM, where upstream has no build).
>
> Version pins: `bwa` 0.7.19 · `bwa-mem2` v2.2.1 · `minibwa` commit [`d6d9f87d`](https://github.com/lh3/minibwa) (`minibwa-0.7`). "ARM" = Graviton4 m8g (arm64/NEON, no SMT); "x86" = AMD m8a (x86-64/AVX-512, no SMT). Both are the general-purpose siblings of the c8g/c8a hosts used through v0.12.0 — same CPU family, same core count, no SMT, but 4 GiB/vCPU instead of 2 so every historical arm fits in memory — so absolute times are not comparable with earlier versions of this table, only ratios within it. (The x86 arm replaced an earlier Intel c7i arm, which ran 16 vCPUs over 8 physical cores under 2-way SMT and so wasn't a real core-for-core match for Graviton's 16 real cores.) No ARM `bwa-mem2` build exists, hence the blank cells there. Every arm for a given arch ran interleaved on one fixed on-demand host — 3 reps each, median wall-clock shown — so these are same-host comparisons, not medians pooled across separate runs. `—` means the release predates the comparator or predates `--fast`. \* v0.13.0 through v0.15.1 rows use the stride-2 suffix-array index described above; all other rows use the stock stride-8 index. Regenerated at each release; see [Benchmarks](https://bwa-mem3.readthedocs.io/en/latest/performance/benchmarks.html).
>
> </details>

> [!WARNING]
> `--fast` is **not alignment-identical** to the default preset — it trades some sensitivity/specificity at the extremes (repetitive/multi-mapping regions, low-`MAPQ` reads) for the speedup above. See "Three ways to run it" below before switching a production pipeline to it.

## Three ways to run it — plain, `--compat`, `--fast`

bwa-mem3 has three alignment modes that differ in *what alignments come out*, not just in speed:

| mode | where reads align | when to use |
|---|---|---|
| **plain** (default) | bwa-mem2's alignments **plus bonafide correctness fixes**, with two extra tags (`MQ:i`, `HN:i`) and an enriched header. On the cells re-measured for release 0.7.1, the complete alignment-record stream (tags stripped) is byte-identical to bwa-mem2 v2.2.1 on `wgs-5M`/`wes-5M`/`hic-1M` (x86 `c6a` AVX2, with a `c6a`/`c8g` cross-arch check confirming the Arm `c8g` NEON build matches) — differing only by those additive tags and the header. Separately, a 1.07M-record HG00096 WGS slice shows zero diverging **primary** alignments (x86, primary-only; not part of the cross-arch or complete-stream checks). | Migrating a pipeline, validating against bwa/bwa-mem2, or any new pipeline. |
| **`--compat=bwa-mem2` / `--compat=bwa-mem`** | Byte-for-byte identical **alignment records** to a **specific** upstream (bwa-mem2 v2.2.1 or bwa 0.7.19), `@PG` excluded and `-t`/`-K` matched. The two targets are **not** interchangeable. | Diff-clean validation against an existing bwa/bwa-mem2 golden. |
| **`--fast`** | Faster, and **not** record-compatible with the default: it reshuffles the low-confidence tail (~85% of the reads it re-places had `MAPQ 0`; the confident `MAPQ`-60 core moves on ≤0.5%, 0.011% on `wgs-5M`) while staying accuracy-neutral against golden truth (≤0.02 pp across the WGS and methylation sims). Figures from the [benchmark](https://bwa-mem3.readthedocs.io/en/latest/performance/benchmarks.html) release-validation cells (`wgs-5M`/`wes-5M`/`panel-twist-5M` at 5 M reads, `hic-1M`/`sbx-1M` at 1 M) across every SIMD tier (AVX2 `c6a`, AVX-512 `c7a`/`c7i`, NEON `c7g`/`c8g`; meth on `m7i`), each a `.4xlarge` host at `-t 16`, `-K 160000000`. | High-throughput pipelines where you care about the confident, uniquely-mapped calls. |

`--compat` is mutually exclusive with `--fast` (and with `--meth` and `--proper-pair-from-emitted`). See [Alignment modes](https://bwa-mem3.readthedocs.io/en/latest/whats-different/modes.html) for the full side-by-side and [Equivalence with bwa-mem2](https://bwa-mem3.readthedocs.io/en/latest/whats-different/equivalence.html) for the field-by-field audit.

By default bwa-mem3 keeps bwa-mem2's command-line defaults, so it drops into an existing pipeline unchanged. For the fastest configuration — and what each recommended deviation from the bwa defaults trades for speed — see [Settings profiles: bwa drop-in vs recommended](https://bwa-mem3.readthedocs.io/en/latest/best-practices/settings-profiles.html).

**Full documentation:** <https://bwa-mem3.readthedocs.io>

## Install

The recommended way to install bwa-mem3 is via [bioconda](https://bioconda.github.io):

```sh
mamba install -c bioconda bwa-mem3
bwa-mem3 version
```

Prebuilt packages are available for `linux-64`, `linux-aarch64`, and `osx-arm64`.

### Build from source

```sh
git clone --recursive https://github.com/fg-labs/bwa-mem3.git
cd bwa-mem3
make
./bwa-mem3 version
```

See the [installation guide](https://bwa-mem3.readthedocs.io/en/latest/getting-started/installation.html) for prerequisites and architecture-specific notes.

## Other subcommands

`bwa-mem3 aln` aligns ultrashort, ancient or degraded DNA with bwa's backtracking algorithm. It fuses `bwa aln` and `bwa samse`/`bwa sampe` into one command and writes SAM that is byte-identical to bwa 0.7.19's on the generated reads the parity checks use, on Linux x86_64, Linux arm64 and macOS arm64 (see [aln](https://bwa-mem3.readthedocs.io/en/latest/cli/aln.html)). `index`, `mem`, `shm`, `re-sa` and `version` are the others.

## Quick links

- [Benchmarks](https://bwa-mem3.readthedocs.io/en/latest/performance/benchmarks.html) — published results for every release, and how they are measured
- [bwa-mem3-rs](https://github.com/fg-labs/bwa-mem3-rs) — Rust bindings for bwa-mem3
- [bioconda recipe](https://github.com/bioconda/bioconda-recipes/tree/master/recipes/bwa-mem3) — conda package on bioconda
- [fgumi](https://github.com/fulcrumgenomics/fgumi) — UMI-aware consensus and deduplication
- [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) — upstream project

## License

MIT. See the [License page](https://bwa-mem3.readthedocs.io/en/latest/reference/license.html) in the docs.

## Citation

Please cite the bwa-mem2 paper (Vasimuddin Md et al., IPDPS 2019). See the [Citation page](https://bwa-mem3.readthedocs.io/en/latest/reference/citation.html) for BibTeX.

## Issues / contributing

File [issues](https://github.com/fg-labs/bwa-mem3/issues) and [pull requests](https://github.com/fg-labs/bwa-mem3/pulls) on [fg-labs/bwa-mem3](https://github.com/fg-labs/bwa-mem3).
