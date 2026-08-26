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

Wall-clock speedup of the current release (v0.10.0) against `bwa` 0.7.19, `bwa-mem2` v2.2.1, and `minibwa`, on the `wgs-5M` sample. Cells are `stock / --fast`.

| arch | wall_s | vs bwa | vs bwa-mem2 | vs minibwa |
|---|---:|---:|---:|---:|
| ARM | 72.23 / 28.26 | 3.32x / 8.47x | — | 0.57x / 1.46x |
| x86 | 44.70 / 19.40 | 4.44x / 10.24x | 2.22x / 5.12x | 0.69x / 1.59x |

> [!TIP]
> **📈 Full release-history table** — every bwa-mem3 release since v0.2.1, full methodology, and version pins.
>
> <details>
> <summary><strong>Click to expand</strong></summary>
>
> **Graviton4 (c8g, arm64/NEON)**
>
> | release | wall_s | vs bwa | vs bwa-mem2 | vs minibwa |
> |---|---:|---:|---:|---:|
> | bwa | 239.49 | 1.00x | — | 0.17x |
> | bwa-mem2 | — | — | — | — |
> | minibwa | 41.35 | 5.79x | — | 1.00x |
> | v0.2.1 | 140.74 | 1.70x | — | 0.29x |
> | v0.2.2 | 141.57 | 1.69x | — | 0.29x |
> | v0.3.0 | 123.35 | 1.94x | — | 0.34x |
> | v0.4.0 | 108.09 | 2.22x | — | 0.38x |
> | v0.5.0 | 106.18 / 38.03 | 2.26x / 6.30x | — | 0.39x / 1.09x |
> | v0.6.0 | 99.01 / 37.06 | 2.42x / 6.46x | — | 0.42x / 1.12x |
> | v0.7.0 | 96.57 / 40.38 | 2.48x / 5.93x | — | 0.43x / 1.02x |
> | v0.8.0 | 77.21 / 28.95 | 3.10x / 8.27x | — | 0.54x / 1.43x |
> | v0.9.0 | 77.49 / 29.00 | 3.09x / 8.26x | — | 0.53x / 1.43x |
> | **v0.10.0** | **72.23 / 28.26** | **3.32x / 8.47x** | — | **0.57x / 1.46x** |
>
> **AMD (c8a, x86)**
>
> | release | wall_s | vs bwa | vs bwa-mem2 | vs minibwa |
> |---|---:|---:|---:|---:|
> | bwa | 198.63 | 1.00x | 0.50x | 0.16x |
> | bwa-mem2 | 99.33 | 2.00x | 1.00x | 0.31x |
> | minibwa | 30.90 | 6.43x | 3.21x | 1.00x |
> | v0.2.1 | 80.04 | 2.48x | 1.24x | 0.39x |
> | v0.2.2 | 72.85 | 2.73x | 1.36x | 0.42x |
> | v0.3.0 | 64.43 | 3.08x | 1.54x | 0.48x |
> | v0.4.0 | 54.57 | 3.64x | 1.82x | 0.57x |
> | v0.5.0 | 55.02 / 23.87 | 3.61x / 8.32x | 1.81x / 4.16x | 0.56x / 1.29x |
> | v0.6.0 | 55.87 / 23.23 | 3.56x / 8.55x | 1.78x / 4.28x | 0.55x / 1.33x |
> | v0.7.0 | 51.64 / 23.78 | 3.85x / 8.35x | 1.92x / 4.18x | 0.60x / 1.30x |
> | v0.8.0 | 47.05 / 18.82 | 4.22x / 10.55x | 2.11x / 5.28x | 0.66x / 1.64x |
> | v0.9.0 | 45.96 / 18.74 | 4.32x / 10.60x | 2.16x / 5.30x | 0.67x / 1.65x |
> | **v0.10.0** | **44.70 / 19.40** | **4.44x / 10.24x** | **2.22x / 5.12x** | **0.69x / 1.59x** |
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
