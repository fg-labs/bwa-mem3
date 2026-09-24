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

Wall-clock speedup of the current release (v0.13.0\*) against `bwa` 0.7.19, `bwa-mem2` v2.2.1, and `minibwa`, on the `wgs-5M` sample. Cells are `stock / --fast`.

| arch | wall_s | vs bwa | vs bwa-mem2 | vs minibwa |
|---|---:|---:|---:|---:|
| ARM | 55.26 / 19.26 | 4.40x / 12.63x | — | 0.74x / 2.14x |
| x86 | 30.57 / 14.16 | 6.47x / 13.97x | 2.75x / 5.94x | 0.99x / 2.15x |

\* v0.13.0 is measured with a **denser suffix-array index** — stride 2 instead of the stock stride 8 — built with `bwa-mem3 re-sa -u 1` (or staged with `bwa-mem3 shm -u 1`), both new in v0.13.0. Output is byte-identical; it costs ~12 GB more memory for the human genome and accounts for roughly 4% of v0.13.0's speedup over v0.12.0. Every other release and every comparator uses the stock index.

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
> | bwa | 243.30 | 1.00x | — | 0.17x | — |
> | bwa-mem2 | — | — | — | — | — |
> | minibwa | 41.14 | 5.91x | — | 1.00x | — |
> | v0.2.1 | 121.94 | 2.00x | — | 0.34x | — |
> | v0.2.2 | 119.32 | 2.04x | — | 0.34x | 1.022x |
> | v0.3.0 | 103.42 | 2.35x | — | 0.40x | 1.154x |
> | v0.4.0 | 105.35 | 2.31x | — | 0.39x | 0.982x |
> | v0.5.0 | 105.53 / 38.36 | 2.31x / 6.34x | — | 0.39x / 1.07x | 0.998x |
> | v0.6.0 | 98.74 / 37.75 | 2.46x / 6.45x | — | 0.42x / 1.09x | 1.069x / 1.016x |
> | v0.7.0 | 92.31 / 40.76 | 2.64x / 5.97x | — | 0.45x / 1.01x | 1.070x / 0.926x |
> | v0.8.0 | 77.19 / 28.91 | 3.15x / 8.41x | — | 0.53x / 1.42x | 1.196x / 1.410x |
> | v0.9.0 | 77.87 / 28.89 | 3.12x / 8.42x | — | 0.53x / 1.42x | 0.991x / 1.001x |
> | v0.10.0 | 72.14 / 28.37 | 3.37x / 8.57x | — | 0.57x / 1.45x | 1.079x / 1.018x |
> | v0.11.0 | 65.11 / 26.82 | 3.74x / 9.07x | — | 0.63x / 1.53x | 1.108x / 1.058x |
> | v0.12.0 | 58.97 / 20.24 | 4.13x / 12.02x | — | 0.70x / 2.03x | 1.104x / 1.325x |
> | **v0.13.0**\* | **55.26 / 19.26** | **4.40x / 12.63x** | **—** | **0.74x / 2.14x** | **1.067x / 1.051x** |
>
> **AMD (m8a, x86-64/AVX-512)**
>
> | release | wall_s | vs bwa | vs bwa-mem2 | vs minibwa | vs prev |
> |---|---:|---:|---:|---:|---:|
> | bwa | 197.78 | 1.00x | 0.43x | 0.15x | — |
> | bwa-mem2 | 84.09 | 2.35x | 1.00x | 0.36x | — |
> | minibwa | 30.38 | 6.51x | 2.77x | 1.00x | — |
> | v0.2.1 | 55.49 | 3.56x | 1.52x | 0.55x | 1.516x |
> | v0.2.2 | 54.01 | 3.66x | 1.56x | 0.56x | 1.027x |
> | v0.3.0 | 50.11 | 3.95x | 1.68x | 0.61x | 1.078x |
> | v0.4.0 | 51.39 | 3.85x | 1.64x | 0.59x | 0.975x |
> | v0.5.0 | 51.91 / 22.31 | 3.81x / 8.86x | 1.62x / 3.77x | 0.59x / 1.36x | 0.990x |
> | v0.6.0 | 51.00 / 21.85 | 3.88x / 9.05x | 1.65x / 3.85x | 0.60x / 1.39x | 1.018x / 1.021x |
> | v0.7.0 | 48.70 / 22.30 | 4.06x / 8.87x | 1.73x / 3.77x | 0.62x / 1.36x | 1.047x / 0.980x |
> | v0.8.0 | 43.79 / 18.68 | 4.52x / 10.59x | 1.92x / 4.50x | 0.69x / 1.63x | 1.112x / 1.194x |
> | v0.9.0 | 43.78 / 18.18 | 4.52x / 10.88x | 1.92x / 4.62x | 0.69x / 1.67x | 1.000x / 1.027x |
> | v0.10.0 | 43.29 / 18.66 | 4.57x / 10.60x | 1.94x / 4.51x | 0.70x / 1.63x | 1.011x / 0.974x |
> | v0.11.0 | 41.73 / 17.66 | 4.74x / 11.20x | 2.01x / 4.76x | 0.73x / 1.72x | 1.037x / 1.057x |
> | v0.12.0 | 33.02 / 15.14 | 5.99x / 13.06x | 2.55x / 5.55x | 0.92x / 2.01x | 1.264x / 1.166x |
> | **v0.13.0**\* | **30.57 / 14.16** | **6.47x / 13.97x** | **2.75x / 5.94x** | **0.99x / 2.15x** | **1.080x / 1.069x** |
>
> `vs prev` is the release-over-release speedup (`prev_wall / this_wall`, `>1` = faster) vs the previous release on this same host, `stock / --fast`. The first release's predecessor is upstream `bwa-mem2` — bwa-mem3 is its successor — so v0.2.1's `vs prev` is its speedup over bwa-mem2 (blank on ARM, where upstream has no build).
>
> Version pins: `bwa` 0.7.19 · `bwa-mem2` v2.2.1 · `minibwa` commit [`d6d9f87d`](https://github.com/lh3/minibwa) (`minibwa-0.7`). "ARM" = Graviton4 m8g (arm64/NEON, no SMT); "x86" = AMD m8a (x86-64/AVX-512, no SMT). Both are the general-purpose siblings of the c8g/c8a hosts used through v0.12.0 — same CPU family, same core count, no SMT, but 4 GiB/vCPU instead of 2 so every historical arm fits in memory — so absolute times are not comparable with earlier versions of this table, only ratios within it. (The x86 arm replaced an earlier Intel c7i arm, which ran 16 vCPUs over 8 physical cores under 2-way SMT and so wasn't a real core-for-core match for Graviton's 16 real cores.) No ARM `bwa-mem2` build exists, hence the blank cells there. Every arm for a given arch ran interleaved on one fixed on-demand host — 3 reps each, median wall-clock shown — so these are same-host comparisons, not medians pooled across separate runs. `—` means the release predates the comparator or predates `--fast`. \* v0.13.0 rows use the stride-2 suffix-array index described above; all other rows use the stock stride-8 index. Regenerated at each release; see [Benchmarks](https://bwa-mem3.readthedocs.io/en/latest/performance/benchmarks.html).
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
