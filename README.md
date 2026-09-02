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

Wall-clock speedup of the current release (v0.11.0) against `bwa` 0.7.19, `bwa-mem2` v2.2.1, and `minibwa`, on the `wgs-5M` sample. Cells are `stock / --fast`.

| arch | wall_s | vs bwa | vs bwa-mem2 | vs minibwa |
|---|---:|---:|---:|---:|
| ARM | 67.54 / 26.37 | 3.58x / 9.18x | — | 0.60x / 1.55x |
| x86 | 43.97 / 19.06 | 4.55x / 10.48x | 2.31x / 5.32x | 0.71x / 1.65x |

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
> | bwa | 242.00 | 1.00x | — | 0.17x | — |
> | bwa-mem2 | — | — | — | — | — |
> | minibwa | 40.77 | 5.94x | — | 1.00x | — |
> | v0.2.1 | 143.57 | 1.69x | — | 0.28x | — |
> | v0.2.2 | 142.50 | 1.70x | — | 0.29x | 1.007x |
> | v0.3.0 | 124.70 | 1.94x | — | 0.33x | 1.143x |
> | v0.4.0 | 106.24 | 2.28x | — | 0.38x | 1.174x |
> | v0.5.0 | 108.76 / 38.02 | 2.23x / 6.37x | — | 0.37x / 1.07x | 0.977x |
> | v0.6.0 | 100.08 / 37.47 | 2.42x / 6.46x | — | 0.41x / 1.09x | 1.087x / 1.015x |
> | v0.7.0 | 94.74 / 41.12 | 2.55x / 5.89x | — | 0.43x / 0.99x | 1.056x / 0.911x |
> | v0.8.0 | 76.96 / 28.70 | 3.14x / 8.43x | — | 0.53x / 1.42x | 1.231x / 1.433x |
> | v0.9.0 | 77.58 / 29.16 | 3.12x / 8.30x | — | 0.53x / 1.40x | 0.992x / 0.984x |
> | v0.10.0 | 75.03 / 27.85 | 3.23x / 8.69x | — | 0.54x / 1.46x | 1.034x / 1.047x |
> | **v0.11.0** | **67.54 / 26.37** | **3.58x / 9.18x** | **—** | **0.60x / 1.55x** | **1.111x / 1.056x** |
>
> **AMD (c8a, x86)**
>
> | release | wall_s | vs bwa | vs bwa-mem2 | vs minibwa | vs prev |
> |---|---:|---:|---:|---:|---:|
> | bwa | 199.88 | 1.00x | 0.51x | 0.16x | — |
> | bwa-mem2 | 101.39 | 1.97x | 1.00x | 0.31x | — |
> | minibwa | 31.39 | 6.37x | 3.23x | 1.00x | — |
> | v0.2.1 | 78.42 | 2.55x | 1.29x | 0.40x | 1.293x |
> | v0.2.2 | 71.83 | 2.78x | 1.41x | 0.44x | 1.092x |
> | v0.3.0 | 70.02 | 2.85x | 1.45x | 0.45x | 1.026x |
> | v0.4.0 | 56.51 | 3.54x | 1.79x | 0.56x | 1.239x |
> | v0.5.0 | 56.21 / 23.99 | 3.56x / 8.33x | 1.80x / 4.23x | 0.56x / 1.31x | 1.005x |
> | v0.6.0 | 53.23 / 23.56 | 3.76x / 8.48x | 1.90x / 4.30x | 0.59x / 1.33x | 1.056x / 1.018x |
> | v0.7.0 | 52.07 / 24.33 | 3.84x / 8.22x | 1.95x / 4.17x | 0.60x / 1.29x | 1.022x / 0.968x |
> | v0.8.0 | 44.78 / 19.38 | 4.46x / 10.31x | 2.26x / 5.23x | 0.70x / 1.62x | 1.163x / 1.255x |
> | v0.9.0 | 45.55 / 19.35 | 4.39x / 10.33x | 2.23x / 5.24x | 0.69x / 1.62x | 0.983x / 1.001x |
> | v0.10.0 | 44.35 / 18.68 | 4.51x / 10.70x | 2.29x / 5.43x | 0.71x / 1.68x | 1.027x / 1.036x |
> | **v0.11.0** | **43.97 / 19.06** | **4.55x / 10.48x** | **2.31x / 5.32x** | **0.71x / 1.65x** | **1.009x / 0.980x** |
>
> `vs prev` is the release-over-release speedup (`prev_wall / this_wall`, `>1` = faster) vs the previous release on this same host, `stock / --fast`. The first release's predecessor is upstream `bwa-mem2` — bwa-mem3 is its successor — so v0.2.1's `vs prev` is its speedup over bwa-mem2 (blank on ARM, where upstream has no build).
>
> Version pins: `bwa` 0.7.19 · `bwa-mem2` v2.2.1 · `minibwa` commit [`d6d9f87d`](https://github.com/lh3/minibwa) (`minibwa-0.7`). "ARM" = Graviton4 c8g (arm64/NEON, no SMT); "x86" = AMD c8a (no SMT — replaces an earlier Intel c7i arm, which ran 16 vCPUs over 8 physical cores under 2-way SMT and so wasn't a real core-for-core match for Graviton's 16 real cores); no ARM `bwa-mem2` build exists, hence the blank cells there. Every arm for a given arch ran interleaved on one fixed on-demand host — 3 reps each, median wall-clock shown — so these are same-host comparisons, not medians pooled across separate runs. `—` means the release predates the comparator or predates `--fast`. Regenerate via `bench release-speedup` in [bwa-mem3-bench](https://github.com/fg-labs/bwa-mem3-bench).
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
