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

Wall-clock speedup of the current release (v0.12.0) against `bwa` 0.7.19, `bwa-mem2` v2.2.1, and `minibwa`, on the `wgs-5M` sample. Cells are `stock / --fast`.

| arch | wall_s | vs bwa | vs bwa-mem2 | vs minibwa |
|---|---:|---:|---:|---:|
| ARM | 59.85 / 20.86 | 4.11x / 11.79x | — | 0.70x / 2.00x |
| x86 | 36.99 / 15.77 | 5.46x / 12.80x | 2.70x / 6.34x | 0.87x / 2.04x |

> [!TIP]
> **📈 Full release-history table** — every bwa-mem3 release since v0.2.1, full methodology, and version pins.
>
> <details>
> <summary><strong>Click to expand</strong></summary>
>
> **Graviton4 (c8g, arm64/NEON)**
>
> | release | wall_s | vs bwa | vs bwa-mem2 | vs minibwa | vs prev |
> |---|---:|---:|---:|---:|---:|
> | bwa | 245.88 | 1.00x | — | 0.17x | — |
> | bwa-mem2 | — | — | — | — | — |
> | minibwa | 41.65 | 5.90x | — | 1.00x | — |
> | v0.2.1 | 143.88 | 1.71x | — | 0.29x | — |
> | v0.2.2 | 143.22 | 1.72x | — | 0.29x | 1.005x |
> | v0.3.0 | 128.74 | 1.91x | — | 0.32x | 1.112x |
> | v0.4.0 | 106.87 | 2.30x | — | 0.39x | 1.205x |
> | v0.5.0 | 110.18 / 40.10 | 2.23x / 6.13x | — | 0.38x / 1.04x | 0.970x |
> | v0.6.0 | 100.06 / 38.56 | 2.46x / 6.38x | — | 0.42x / 1.08x | 1.101x / 1.040x |
> | v0.7.0 | 96.07 / 41.32 | 2.56x / 5.95x | — | 0.43x / 1.01x | 1.042x / 0.933x |
> | v0.8.0 | 77.72 / 29.06 | 3.16x / 8.46x | — | 0.54x / 1.43x | 1.236x / 1.422x |
> | v0.9.0 | 78.12 / 29.30 | 3.15x / 8.39x | — | 0.53x / 1.42x | 0.995x / 0.992x |
> | v0.10.0 | 72.89 / 30.45 | 3.37x / 8.07x | — | 0.57x / 1.37x | 1.072x / 0.962x |
> | v0.11.0 | 66.38 / 27.23 | 3.70x / 9.03x | — | 0.63x / 1.53x | 1.098x / 1.118x |
> | **v0.12.0** | **59.85 / 20.86** | **4.11x / 11.79x** | **—** | **0.70x / 2.00x** | **1.109x / 1.305x** |
>
> **AMD (c8a, x86-64/AVX-512)**
>
> | release | wall_s | vs bwa | vs bwa-mem2 | vs minibwa | vs prev |
> |---|---:|---:|---:|---:|---:|
> | bwa | 201.79 | 1.00x | 0.50x | 0.16x | — |
> | bwa-mem2 | 99.92 | 2.02x | 1.00x | 0.32x | — |
> | minibwa | 32.15 | 6.28x | 3.11x | 1.00x | — |
> | v0.2.1 | 76.49 | 2.64x | 1.31x | 0.42x | 1.306x |
> | v0.2.2 | 77.15 | 2.62x | 1.30x | 0.42x | 0.991x |
> | v0.3.0 | 71.93 | 2.81x | 1.39x | 0.45x | 1.073x |
> | v0.4.0 | 57.29 | 3.52x | 1.74x | 0.56x | 1.256x |
> | v0.5.0 | 55.74 / 23.44 | 3.62x / 8.61x | 1.79x / 4.26x | 0.58x / 1.37x | 1.028x |
> | v0.6.0 | 55.08 / 23.17 | 3.66x / 8.71x | 1.81x / 4.31x | 0.58x / 1.39x | 1.012x / 1.012x |
> | v0.7.0 | 52.67 / 23.85 | 3.83x / 8.46x | 1.90x / 4.19x | 0.61x / 1.35x | 1.046x / 0.971x |
> | v0.8.0 | 47.94 / 19.22 | 4.21x / 10.50x | 2.08x / 5.20x | 0.67x / 1.67x | 1.099x / 1.241x |
> | v0.9.0 | 47.22 / 19.93 | 4.27x / 10.12x | 2.12x / 5.01x | 0.68x / 1.61x | 1.015x / 0.964x |
> | v0.10.0 | 47.73 / 20.21 | 4.23x / 9.98x | 2.09x / 4.94x | 0.67x / 1.59x | 0.989x / 0.986x |
> | v0.11.0 | 46.29 / 18.21 | 4.36x / 11.08x | 2.16x / 5.49x | 0.69x / 1.77x | 1.031x / 1.110x |
> | **v0.12.0** | **36.99 / 15.77** | **5.46x / 12.80x** | **2.70x / 6.34x** | **0.87x / 2.04x** | **1.251x / 1.155x** |
>
> `vs prev` is the release-over-release speedup (`prev_wall / this_wall`, `>1` = faster) vs the previous release on this same host, `stock / --fast`. The first release's predecessor is upstream `bwa-mem2` — bwa-mem3 is its successor — so v0.2.1's `vs prev` is its speedup over bwa-mem2 (blank on ARM, where upstream has no build).
>
> Version pins: `bwa` 0.7.19 · `bwa-mem2` v2.2.1 · `minibwa` commit [`d6d9f87d`](https://github.com/lh3/minibwa) (`minibwa-0.7`). "ARM" = Graviton4 c8g (arm64/NEON, no SMT); "x86" = AMD c8a (x86-64/AVX-512, no SMT — replaces an earlier Intel c7i arm, which ran 16 vCPUs over 8 physical cores under 2-way SMT and so wasn't a real core-for-core match for Graviton's 16 real cores); no ARM `bwa-mem2` build exists, hence the blank cells there. Every arm for a given arch ran interleaved on one fixed on-demand host — 3 reps each, median wall-clock shown — so these are same-host comparisons, not medians pooled across separate runs. `—` means the release predates the comparator or predates `--fast`. Regenerate via `bench release-speedup` in [bwa-mem3-bench](https://github.com/fg-labs/bwa-mem3-bench).
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
| **`--fast`** | Faster, and **not** record-compatible with the default: it reshuffles the low-confidence tail (~85% of the reads it re-places had `MAPQ 0`; the confident `MAPQ`-60 core moves on ≤0.5%, 0.011% on `wgs-5M`) while staying accuracy-neutral against golden truth (≤0.02 pp across the WGS and methylation sims). Figures from the [bwa-mem3-bench](https://github.com/fg-labs/bwa-mem3-bench) release-validation cells (`wgs-5M`/`wes-5M`/`panel-twist-5M` at 5 M reads, `hic-1M`/`sbx-1M` at 1 M) across every SIMD tier (AVX2 `c6a`, AVX-512 `c7a`/`c7i`, NEON `c7g`/`c8g`; meth on `m7i`), each a `.4xlarge` host at `-t 16`, `-K 160000000`. | High-throughput pipelines where you care about the confident, uniquely-mapped calls. |

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

- [bwa-mem3-bench](https://github.com/fg-labs/bwa-mem3-bench) — benchmarking harness across CPU architectures
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
