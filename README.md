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

Wall-clock speedup of the current release (v0.16.0\*) against `bwa` 0.7.19, `bwa-mem2` v2.2.1, and `minibwa`, on the `wgs-5M` sample. Cells are `stock / --fast`.

| arch | wall_s | vs bwa | vs bwa-mem2 | vs minibwa |
|---|---:|---:|---:|---:|
| ARM | 28.72 / 12.46 | 8.16x / 18.81x | — | 1.40x / 3.24x |
| x86 | 16.36 / 9.05 | 12.15x / 21.96x | 5.41x / 9.79x | 1.84x / 3.33x |

\* v0.16.0, like v0.13.0 through v0.15.1, is measured with a **denser suffix-array index** — stride 2 instead of the stock stride 8 — built with `bwa-mem3 re-sa -u 1` (or staged with `bwa-mem3 shm -u 1`), both new in v0.13.0. Output is byte-identical by construction (each added suffix-array sample is computed from the on-disk ones), and measured so: a full-content `fgumi compare bams` of `wgs-5M` (10,030,558 records) aligned on the stride-2 index against the stock stride-8 index — v0.13.0 against the output-identical v0.12.0 — found 0 differences on both benchmark hosts, Graviton4 m8g (NEON) and AMD m8a (AVX-512), and the change that added it ([#510](https://github.com/fg-labs/bwa-mem3/pull/510)) found 0 differences on 5M-pair hg38 WGS and WES slices. It costs ~12 GB more memory for the human genome. v0.13.0 through v0.16.0 are all measured on it in the tables below, so v0.16.0's `vs prev` is a like-for-like comparison. Every older release and every comparator uses the stock index.

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
> | bwa | 234.36 | 1.00x | — | 0.17x | — |
> | bwa-mem2 | — | — | — | — | — |
> | minibwa | 40.35 | 5.81x | — | 1.00x | — |
> | v0.2.1 | 120.02 | 1.95x | — | 0.34x | — |
> | v0.2.2 | 117.51 | 1.99x | — | 0.34x | 1.021x |
> | v0.3.0 | 101.31 | 2.31x | — | 0.40x | 1.160x |
> | v0.4.0 | 103.46 | 2.27x | — | 0.39x | 0.979x |
> | v0.5.0 | 103.58 / 36.50 | 2.26x / 6.42x | — | 0.39x / 1.11x | 0.999x |
> | v0.6.0 | 96.61 / 35.67 | 2.43x / 6.57x | — | 0.42x / 1.13x | 1.072x / 1.023x |
> | v0.7.0 | 91.00 / 39.04 | 2.58x / 6.00x | — | 0.44x / 1.03x | 1.062x / 0.914x |
> | v0.8.0 | 75.74 / 27.51 | 3.09x / 8.52x | — | 0.53x / 1.47x | 1.201x / 1.419x |
> | v0.9.0 | 75.81 / 27.63 | 3.09x / 8.48x | — | 0.53x / 1.46x | 0.999x / 0.996x |
> | v0.10.0 | 70.30 / 26.89 | 3.33x / 8.72x | — | 0.57x / 1.50x | 1.078x / 1.028x |
> | v0.11.0 | 63.38 / 25.36 | 3.70x / 9.24x | — | 0.64x / 1.59x | 1.109x / 1.060x |
> | v0.12.0 | 57.63 / 18.55 | 4.07x / 12.63x | — | 0.70x / 2.18x | 1.100x / 1.367x |
> | v0.13.0\* | 53.91 / 17.86 | 4.35x / 13.12x | — | 0.75x / 2.26x | 1.069x / 1.039x |
> | v0.14.0\* | 40.99 / 16.94 | 5.72x / 13.83x | — | 0.98x / 2.38x | 1.315x / 1.054x |
> | v0.14.1\* | 34.56 / 15.72 | 6.78x / 14.91x | — | 1.17x / 2.57x | 1.186x / 1.078x |
> | v0.15.0\* | 33.82 / 15.11 | 6.93x / 15.51x | — | 1.19x / 2.67x | 1.022x / 1.040x |
> | v0.15.1\* | 31.29 / 14.08 | 7.49x / 16.64x | — | 1.29x / 2.87x | 1.081x / 1.073x |
> | **v0.16.0**\* | **28.72 / 12.46** | **8.16x / 18.81x** | **—** | **1.40x / 3.24x** | **1.089x / 1.130x** |
>
> **AMD (m8a, x86-64/AVX-512)**
>
> | release | wall_s | vs bwa | vs bwa-mem2 | vs minibwa | vs prev |
> |---|---:|---:|---:|---:|---:|
> | bwa | 198.77 | 1.00x | 0.45x | 0.15x | — |
> | bwa-mem2 | 88.57 | 2.24x | 1.00x | 0.34x | — |
> | minibwa | 30.10 | 6.60x | 2.94x | 1.00x | — |
> | v0.2.1 | 54.58 | 3.64x | 1.62x | 0.55x | 1.623x |
> | v0.2.2 | 52.41 | 3.79x | 1.69x | 0.57x | 1.041x |
> | v0.3.0 | 49.07 | 4.05x | 1.80x | 0.61x | 1.068x |
> | v0.4.0 | 51.14 | 3.89x | 1.73x | 0.59x | 0.960x |
> | v0.5.0 | 51.07 / 21.97 | 3.89x / 9.05x | 1.73x / 4.03x | 0.59x / 1.37x | 1.001x |
> | v0.6.0 | 49.97 / 21.53 | 3.98x / 9.23x | 1.77x / 4.11x | 0.60x / 1.40x | 1.022x / 1.020x |
> | v0.7.0 | 47.96 / 21.87 | 4.14x / 9.09x | 1.85x / 4.05x | 0.63x / 1.38x | 1.042x / 0.984x |
> | v0.8.0 | 43.05 / 17.81 | 4.62x / 11.16x | 2.06x / 4.97x | 0.70x / 1.69x | 1.114x / 1.228x |
> | v0.9.0 | 42.73 / 17.77 | 4.65x / 11.19x | 2.07x / 4.98x | 0.70x / 1.69x | 1.007x / 1.002x |
> | v0.10.0 | 42.68 / 17.88 | 4.66x / 11.12x | 2.08x / 4.95x | 0.71x / 1.68x | 1.001x / 0.994x |
> | v0.11.0 | 41.38 / 17.10 | 4.80x / 11.62x | 2.14x / 5.18x | 0.73x / 1.76x | 1.031x / 1.046x |
> | v0.12.0 | 32.61 / 14.13 | 6.10x / 14.07x | 2.72x / 6.27x | 0.92x / 2.13x | 1.269x / 1.210x |
> | v0.13.0\* | 29.56 / 13.47 | 6.72x / 14.76x | 3.00x / 6.58x | 1.02x / 2.23x | 1.103x / 1.049x |
> | v0.14.0\* | 23.46 / 12.43 | 8.47x / 15.99x | 3.78x / 7.13x | 1.28x / 2.42x | 1.260x / 1.084x |
> | v0.14.1\* | 22.04 / 12.14 | 9.02x / 16.37x | 4.02x / 7.30x | 1.37x / 2.48x | 1.064x / 1.024x |
> | v0.15.0\* | 20.68 / 11.41 | 9.61x / 17.42x | 4.28x / 7.76x | 1.46x / 2.64x | 1.066x / 1.064x |
> | v0.15.1\* | 18.65 / 10.60 | 10.66x / 18.75x | 4.75x / 8.36x | 1.61x / 2.84x | 1.109x / 1.076x |
> | **v0.16.0**\* | **16.36 / 9.05** | **12.15x / 21.96x** | **5.41x / 9.79x** | **1.84x / 3.33x** | **1.140x / 1.171x** |
>
> `vs prev` is the release-over-release speedup (`prev_wall / this_wall`, `>1` = faster) vs the previous release on this same host, `stock / --fast`. Since v0.16.0 every bwa-mem3 release in these tables, old ones included, writes BAM natively (`--bam=0`) while `bwa`, `bwa-mem2` and `minibwa` pipe SAM through `samtools view -u`; earlier tables sent every arm through SAM, so `vs prev` stays like-for-like but the comparator ratios are not comparable across that change. The first release's predecessor is upstream `bwa-mem2` — bwa-mem3 is its successor — so v0.2.1's `vs prev` is its speedup over bwa-mem2 (blank on ARM, where upstream has no build).
>
> Version pins: `bwa` 0.7.19 · `bwa-mem2` v2.2.1 · `minibwa` commit [`d6d9f87d`](https://github.com/lh3/minibwa) (`minibwa-0.7`). "ARM" = Graviton4 m8g (arm64/NEON, no SMT); "x86" = AMD m8a (x86-64/AVX-512, no SMT). Both are the general-purpose siblings of the c8g/c8a hosts used through v0.12.0 — same CPU family, same core count, no SMT, but 4 GiB/vCPU instead of 2 so every historical arm fits in memory — so absolute times are not comparable with earlier versions of this table, only ratios within it. (The x86 arm replaced an earlier Intel c7i arm, which ran 16 vCPUs over 8 physical cores under 2-way SMT and so wasn't a real core-for-core match for Graviton's 16 real cores.) No ARM `bwa-mem2` build exists, hence the blank cells there. Every arm for a given arch ran interleaved on one fixed on-demand host — 5 rounds for v0.16.0 and v0.15.1, 2 for the rest, median wall-clock shown — so these are same-host comparisons, not medians pooled across separate runs. `—` means the release predates the comparator or predates `--fast`. \* v0.13.0 through v0.16.0 rows use the stride-2 suffix-array index described above; all other rows use the stock stride-8 index. Regenerated at each release; see [Benchmarks](https://bwa-mem3.readthedocs.io/en/latest/performance/benchmarks.html).
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

`bwa-mem3 aln` aligns ultrashort, ancient or degraded DNA with bwa's backtracking algorithm. It fuses `bwa aln` and `bwa samse`/`bwa sampe` into one command and writes SAM that is byte-identical to bwa 0.7.19's, apart from the `@PG` line, on the generated reads the parity checks use, on Linux x86_64, Linux arm64 and macOS arm64 (see [aln](https://bwa-mem3.readthedocs.io/en/latest/cli/aln.html)). `index`, `mem`, `shm`, `re-sa` and `version` are the others.

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
