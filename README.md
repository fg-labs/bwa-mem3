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

Wall-clock speedup of the current release (v0.14.0\*) against `bwa` 0.7.19, `bwa-mem2` v2.2.1, and `minibwa`, on the `wgs-5M` sample. Cells are `stock / --fast`.

| arch | wall_s | vs bwa | vs bwa-mem2 | vs minibwa |
|---|---:|---:|---:|---:|
| ARM | 42.23 / 18.08 | 5.69x / 13.30x | — | 0.96x / 2.24x |
| x86 | 24.15 / 13.07 | 8.07x / 14.91x | 3.52x / 6.51x | 1.26x / 2.32x |

\* v0.14.0, like v0.13.0, is measured with a **denser suffix-array index** — stride 2 instead of the stock stride 8 — built with `bwa-mem3 re-sa -u 1` (or staged with `bwa-mem3 shm -u 1`), both new in v0.13.0. Output is byte-identical by construction (each added suffix-array sample is computed from the on-disk ones), and measured so: a full-content `fgumi compare bams` of `wgs-5M` (10,030,558 records) aligned on the stride-2 index against the stock stride-8 index — v0.13.0 against the output-identical v0.12.0 — found 0 differences on both benchmark hosts, Graviton4 m8g (NEON) and AMD m8a (AVX-512), and the change that added it ([#510](https://github.com/fg-labs/bwa-mem3/pull/510)) found 0 differences on 5M-pair hg38 WGS and WES slices. It costs ~12 GB more memory for the human genome. Both v0.13.0 and v0.14.0 are measured on it in the tables below, so v0.14.0's `vs prev` is a like-for-like comparison. Every older release and every comparator uses the stock index.

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
> | bwa | 240.38 | 1.00x | — | 0.17x | — |
> | bwa-mem2 | — | — | — | — | — |
> | minibwa | 40.59 | 5.92x | — | 1.00x | — |
> | v0.2.1 | 120.54 | 1.99x | — | 0.34x | — |
> | v0.2.2 | 122.73 | 1.96x | — | 0.33x | 0.982x |
> | v0.3.0 | 102.11 | 2.35x | — | 0.40x | 1.202x |
> | v0.4.0 | 104.14 | 2.31x | — | 0.39x | 0.980x |
> | v0.5.0 | 104.47 / 37.23 | 2.30x / 6.46x | — | 0.39x / 1.09x | 0.997x |
> | v0.6.0 | 97.69 / 36.59 | 2.46x / 6.57x | — | 0.42x / 1.11x | 1.069x / 1.017x |
> | v0.7.0 | 92.06 / 39.78 | 2.61x / 6.04x | — | 0.44x / 1.02x | 1.061x / 0.920x |
> | v0.8.0 | 76.34 / 28.33 | 3.15x / 8.48x | — | 0.53x / 1.43x | 1.206x / 1.404x |
> | v0.9.0 | 76.65 / 28.34 | 3.14x / 8.48x | — | 0.53x / 1.43x | 0.996x / 1.000x |
> | v0.10.0 | 71.11 / 27.80 | 3.38x / 8.65x | — | 0.57x / 1.46x | 1.078x / 1.020x |
> | v0.11.0 | 64.32 / 26.30 | 3.74x / 9.14x | — | 0.63x / 1.54x | 1.106x / 1.057x |
> | v0.12.0 | 58.68 / 19.22 | 4.10x / 12.51x | — | 0.69x / 2.11x | 1.096x / 1.369x |
> | v0.13.0\* | 55.02 / 18.73 | 4.37x / 12.84x | — | 0.74x / 2.17x | 1.066x / 1.026x |
> | **v0.14.0**\* | **42.23 / 18.08** | **5.69x / 13.30x** | **—** | **0.96x / 2.24x** | **1.303x / 1.036x** |
>
> **AMD (m8a, x86-64/AVX-512)**
>
> | release | wall_s | vs bwa | vs bwa-mem2 | vs minibwa | vs prev |
> |---|---:|---:|---:|---:|---:|
> | bwa | 194.94 | 1.00x | 0.44x | 0.16x | — |
> | bwa-mem2 | 85.04 | 2.29x | 1.00x | 0.36x | — |
> | minibwa | 30.36 | 6.42x | 2.80x | 1.00x | — |
> | v0.2.1 | 56.63 | 3.44x | 1.50x | 0.54x | 1.502x |
> | v0.2.2 | 53.06 | 3.67x | 1.60x | 0.57x | 1.067x |
> | v0.3.0 | 51.30 | 3.80x | 1.66x | 0.59x | 1.034x |
> | v0.4.0 | 51.49 | 3.79x | 1.65x | 0.59x | 0.996x |
> | v0.5.0 | 51.44 / 22.84 | 3.79x / 8.53x | 1.65x / 3.72x | 0.59x / 1.33x | 1.001x |
> | v0.6.0 | 50.03 / 22.34 | 3.90x / 8.73x | 1.70x / 3.81x | 0.61x / 1.36x | 1.028x / 1.023x |
> | v0.7.0 | 48.34 / 22.31 | 4.03x / 8.74x | 1.76x / 3.81x | 0.63x / 1.36x | 1.035x / 1.001x |
> | v0.8.0 | 43.37 / 18.69 | 4.49x / 10.43x | 1.96x / 4.55x | 0.70x / 1.62x | 1.115x / 1.194x |
> | v0.9.0 | 43.29 / 18.17 | 4.50x / 10.73x | 1.96x / 4.68x | 0.70x / 1.67x | 1.002x / 1.028x |
> | v0.10.0 | 42.82 / 18.16 | 4.55x / 10.73x | 1.99x / 4.68x | 0.71x / 1.67x | 1.011x / 1.001x |
> | v0.11.0 | 41.77 / 17.64 | 4.67x / 11.05x | 2.04x / 4.82x | 0.73x / 1.72x | 1.025x / 1.030x |
> | v0.12.0 | 32.59 / 14.64 | 5.98x / 13.31x | 2.61x / 5.81x | 0.93x / 2.07x | 1.282x / 1.205x |
> | v0.13.0\* | 29.83 / 13.65 | 6.54x / 14.28x | 2.85x / 6.23x | 1.02x / 2.22x | 1.093x / 1.073x |
> | **v0.14.0**\* | **24.15 / 13.07** | **8.07x / 14.91x** | **3.52x / 6.51x** | **1.26x / 2.32x** | **1.235x / 1.045x** |
>
> `vs prev` is the release-over-release speedup (`prev_wall / this_wall`, `>1` = faster) vs the previous release on this same host, `stock / --fast`. The first release's predecessor is upstream `bwa-mem2` — bwa-mem3 is its successor — so v0.2.1's `vs prev` is its speedup over bwa-mem2 (blank on ARM, where upstream has no build).
>
> Version pins: `bwa` 0.7.19 · `bwa-mem2` v2.2.1 · `minibwa` commit [`d6d9f87d`](https://github.com/lh3/minibwa) (`minibwa-0.7`). "ARM" = Graviton4 m8g (arm64/NEON, no SMT); "x86" = AMD m8a (x86-64/AVX-512, no SMT). Both are the general-purpose siblings of the c8g/c8a hosts used through v0.12.0 — same CPU family, same core count, no SMT, but 4 GiB/vCPU instead of 2 so every historical arm fits in memory — so absolute times are not comparable with earlier versions of this table, only ratios within it. (The x86 arm replaced an earlier Intel c7i arm, which ran 16 vCPUs over 8 physical cores under 2-way SMT and so wasn't a real core-for-core match for Graviton's 16 real cores.) No ARM `bwa-mem2` build exists, hence the blank cells there. Every arm for a given arch ran interleaved on one fixed on-demand host — 3 reps each, median wall-clock shown — so these are same-host comparisons, not medians pooled across separate runs. `—` means the release predates the comparator or predates `--fast`. \* v0.13.0 and v0.14.0 rows use the stride-2 suffix-array index described above; all other rows use the stock stride-8 index. Regenerated at each release; see [Benchmarks](https://bwa-mem3.readthedocs.io/en/latest/performance/benchmarks.html).
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
