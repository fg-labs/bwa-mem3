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

Wall-clock speedup of the current release (v0.15.0\*) against `bwa` 0.7.19, `bwa-mem2` v2.2.1, and `minibwa`, on the `wgs-5M` sample. Cells are `stock / --fast`.

| arch | wall_s | vs bwa | vs bwa-mem2 | vs minibwa |
|---|---:|---:|---:|---:|
| ARM | 36.85 / 17.08 | 6.93x / 14.96x | — | 1.13x / 2.44x |
| x86 | 22.83 / 12.07 | 8.65x / 16.36x | 3.70x / 7.00x | 1.33x / 2.52x |

\* v0.15.0, like v0.13.0, v0.14.0 and v0.14.1, is measured with a **denser suffix-array index** — stride 2 instead of the stock stride 8 — built with `bwa-mem3 re-sa -u 1` (or staged with `bwa-mem3 shm -u 1`), both new in v0.13.0. Output is byte-identical by construction (each added suffix-array sample is computed from the on-disk ones), and measured so: a full-content `fgumi compare bams` of `wgs-5M` (10,030,558 records) aligned on the stride-2 index against the stock stride-8 index — v0.13.0 against the output-identical v0.12.0 — found 0 differences on both benchmark hosts, Graviton4 m8g (NEON) and AMD m8a (AVX-512), and the change that added it ([#510](https://github.com/fg-labs/bwa-mem3/pull/510)) found 0 differences on 5M-pair hg38 WGS and WES slices. It costs ~12 GB more memory for the human genome. v0.13.0 through v0.15.0 are all measured on it in the tables below, so v0.15.0's `vs prev` is a like-for-like comparison. Every older release and every comparator uses the stock index.

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
> | bwa | 255.55 | 1.00x | — | 0.16x | — |
> | bwa-mem2 | — | — | — | — | — |
> | minibwa | 41.68 | 6.13x | — | 1.00x | — |
> | v0.2.1 | 129.89 | 1.97x | — | 0.32x | — |
> | v0.2.2 | 120.02 | 2.13x | — | 0.35x | 1.082x |
> | v0.3.0 | 106.00 | 2.41x | — | 0.39x | 1.132x |
> | v0.4.0 | 105.08 | 2.43x | — | 0.40x | 1.009x |
> | v0.5.0 | 105.40 / 37.99 | 2.42x / 6.73x | — | 0.40x / 1.10x | 0.997x |
> | v0.6.0 | 100.43 / 37.20 | 2.54x / 6.87x | — | 0.42x / 1.12x | 1.049x / 1.021x |
> | v0.7.0 | 93.35 / 40.90 | 2.74x / 6.25x | — | 0.45x / 1.02x | 1.076x / 0.910x |
> | v0.8.0 | 76.88 / 28.38 | 3.32x / 9.00x | — | 0.54x / 1.47x | 1.214x / 1.441x |
> | v0.9.0 | 78.03 / 28.93 | 3.28x / 8.83x | — | 0.53x / 1.44x | 0.985x / 0.981x |
> | v0.10.0 | 72.60 / 27.88 | 3.52x / 9.17x | — | 0.57x / 1.49x | 1.075x / 1.038x |
> | v0.11.0 | 64.74 / 27.35 | 3.95x / 9.34x | — | 0.64x / 1.52x | 1.121x / 1.019x |
> | v0.12.0 | 59.02 / 19.82 | 4.33x / 12.89x | — | 0.71x / 2.10x | 1.097x / 1.380x |
> | v0.13.0\* | 56.24 / 19.80 | 4.54x / 12.91x | — | 0.74x / 2.11x | 1.049x / 1.001x |
> | v0.14.0\* | 44.06 / 18.61 | 5.80x / 13.73x | — | 0.95x / 2.24x | 1.276x / 1.064x |
> | v0.14.1\* | 36.04 / 17.60 | 7.09x / 14.52x | — | 1.16x / 2.37x | 1.223x / 1.057x |
> | **v0.15.0**\* | **36.85 / 17.08** | **6.93x / 14.96x** | **—** | **1.13x / 2.44x** | **1.015x† / 1.030x** |
>
> **AMD (m8a, x86-64/AVX-512)**
>
> | release | wall_s | vs bwa | vs bwa-mem2 | vs minibwa | vs prev |
> |---|---:|---:|---:|---:|---:|
> | bwa | 197.47 | 1.00x | 0.43x | 0.15x | — |
> | bwa-mem2 | 84.53 | 2.34x | 1.00x | 0.36x | — |
> | minibwa | 30.44 | 6.49x | 2.78x | 1.00x | — |
> | v0.2.1 | 60.31 | 3.27x | 1.40x | 0.50x | 1.402x |
> | v0.2.2 | 57.24 | 3.45x | 1.48x | 0.53x | 1.054x |
> | v0.3.0 | 53.47 | 3.69x | 1.58x | 0.57x | 1.071x |
> | v0.4.0 | 53.20 | 3.71x | 1.59x | 0.57x | 1.005x |
> | v0.5.0 | 52.33 / 22.85 | 3.77x / 8.64x | 1.62x / 3.70x | 0.58x / 1.33x | 1.017x |
> | v0.6.0 | 52.66 / 22.35 | 3.75x / 8.84x | 1.61x / 3.78x | 0.58x / 1.36x | 0.994x / 1.022x |
> | v0.7.0 | 48.98 / 22.36 | 4.03x / 8.83x | 1.73x / 3.78x | 0.62x / 1.36x | 1.075x / 1.000x |
> | v0.8.0 | 45.01 / 18.67 | 4.39x / 10.58x | 1.88x / 4.53x | 0.68x / 1.63x | 1.088x / 1.198x |
> | v0.9.0 | 44.70 / 18.19 | 4.42x / 10.86x | 1.89x / 4.65x | 0.68x / 1.67x | 1.007x / 1.026x |
> | v0.10.0 | 43.47 / 18.68 | 4.54x / 10.57x | 1.94x / 4.53x | 0.70x / 1.63x | 1.028x / 0.974x |
> | v0.11.0 | 42.25 / 17.66 | 4.67x / 11.18x | 2.00x / 4.79x | 0.72x / 1.72x | 1.029x / 1.058x |
> | v0.12.0 | 33.18 / 14.65 | 5.95x / 13.48x | 2.55x / 5.77x | 0.92x / 2.08x | 1.273x / 1.205x |
> | v0.13.0\* | 31.10 / 14.16 | 6.35x / 13.95x | 2.72x / 5.97x | 0.98x / 2.15x | 1.067x / 1.035x |
> | v0.14.0\* | 25.13 / 13.58 | 7.86x / 14.54x | 3.36x / 6.22x | 1.21x / 2.24x | 1.238x / 1.043x |
> | v0.14.1\* | 24.11 / 13.08 | 8.19x / 15.10x | 3.51x / 6.46x | 1.26x / 2.33x | 1.042x / 1.038x |
> | **v0.15.0**\* | **22.83 / 12.07** | **8.65x / 16.36x** | **3.70x / 7.00x** | **1.33x / 2.52x** | **1.056x / 1.084x** |
>
> `vs prev` is the release-over-release speedup (`prev_wall / this_wall`, `>1` = faster) vs the previous release on this same host, `stock / --fast`. † On ARM the three interleaved reps above put v0.15.0 at 0.978x of v0.14.1 (36.85 vs 36.04 s); a dedicated 7-rep interleaved A/B on `wgs-5M` (fresh m8g.4xlarge, arm64/NEON, `-t 16`, one discarded warmup, prewarmed page cache) measured 1.015x (36.17 vs 36.72 s), faster in all 7 pairs, so the A/B figure is shown. The first release's predecessor is upstream `bwa-mem2` — bwa-mem3 is its successor — so v0.2.1's `vs prev` is its speedup over bwa-mem2 (blank on ARM, where upstream has no build).
>
> Version pins: `bwa` 0.7.19 · `bwa-mem2` v2.2.1 · `minibwa` commit [`d6d9f87d`](https://github.com/lh3/minibwa) (`minibwa-0.7`). "ARM" = Graviton4 m8g (arm64/NEON, no SMT); "x86" = AMD m8a (x86-64/AVX-512, no SMT). Both are the general-purpose siblings of the c8g/c8a hosts used through v0.12.0 — same CPU family, same core count, no SMT, but 4 GiB/vCPU instead of 2 so every historical arm fits in memory — so absolute times are not comparable with earlier versions of this table, only ratios within it. (The x86 arm replaced an earlier Intel c7i arm, which ran 16 vCPUs over 8 physical cores under 2-way SMT and so wasn't a real core-for-core match for Graviton's 16 real cores.) No ARM `bwa-mem2` build exists, hence the blank cells there. Every arm for a given arch ran interleaved on one fixed on-demand host — 3 reps each, median wall-clock shown — so these are same-host comparisons, not medians pooled across separate runs. `—` means the release predates the comparator or predates `--fast`. \* v0.13.0 through v0.15.0 rows use the stride-2 suffix-array index described above; all other rows use the stock stride-8 index. Regenerated at each release; see [Benchmarks](https://bwa-mem3.readthedocs.io/en/latest/performance/benchmarks.html).
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
