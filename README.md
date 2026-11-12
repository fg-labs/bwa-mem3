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

Wall-clock speedup of the current release (v0.17.0\*) against `bwa` 0.7.19, `bwa-mem2` v2.2.1, and `minibwa`, on the `wgs-5M` sample. Cells are `stock / --fast`.

| arch | wall_s | vs bwa | vs bwa-mem2 | vs minibwa |
|---|---:|---:|---:|---:|
| ARM | 28.19 / 12.33 | 8.36x / 19.12x | — | 1.43x / 3.28x |
| x86 | 16.19 / 9.69 | 12.86x / 21.48x | 5.75x / 9.60x | 1.99x / 3.33x |

\* v0.17.0, like v0.13.0 through v0.16.0, is measured with a **denser suffix-array index** — stride 2 instead of the stock stride 8 — built with `bwa-mem3 re-sa -u 1` (or staged with `bwa-mem3 shm -u 1`), both new in v0.13.0. Alignment records are byte-identical by construction (each added suffix-array sample is computed from the on-disk ones), and measured so; only the `@PG` header line differs, since its command line names the index path: a full-content `fgumi compare bams` of `wgs-5M` (10,030,558 records) aligned on the stride-2 index against the stock stride-8 index — v0.13.0 against the output-identical v0.12.0 — found 0 differences on both benchmark hosts, Graviton4 m8g (NEON) and AMD m8a (AVX-512), and the change that added it ([#510](https://github.com/fg-labs/bwa-mem3/pull/510)) found 0 differences in a record-by-record comparison of core SAM fields and tags on 5M-pair hg38 WGS and WES slices (that PR does not record the host or SIMD tier). It costs ~12 GB more memory for the human genome. v0.13.0 through v0.17.0 are all measured on it in the tables below, so v0.17.0's `vs prev` is a like-for-like comparison. Every older release and every comparator uses the stock index.

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
> | bwa | 235.69 | 1.00x | — | 0.17x | — |
> | bwa-mem2 | — | — | — | — | — |
> | minibwa | 40.43 | 5.83x | — | 1.00x | — |
> | v0.2.1 | 120.01 | 1.96x | — | 0.34x | — |
> | v0.2.2 | 117.71 | 2.00x | — | 0.34x | 1.020x |
> | v0.3.0 | 101.37 | 2.33x | — | 0.40x | 1.161x |
> | v0.4.0 | 103.54 | 2.28x | — | 0.39x | 0.979x |
> | v0.5.0 | 103.74 / 36.54 | 2.27x / 6.45x | — | 0.39x / 1.11x | 0.998x |
> | v0.6.0 | 96.67 / 35.65 | 2.44x / 6.61x | — | 0.42x / 1.13x | 1.073x / 1.025x |
> | v0.7.0 | 91.10 / 39.10 | 2.59x / 6.03x | — | 0.44x / 1.03x | 1.061x / 0.912x |
> | v0.8.0 | 75.82 / 27.50 | 3.11x / 8.57x | — | 0.53x / 1.47x | 1.202x / 1.422x |
> | v0.9.0 | 75.78 / 27.60 | 3.11x / 8.54x | — | 0.53x / 1.46x | 1.001x / 0.996x |
> | v0.10.0 | 70.32 / 26.99 | 3.35x / 8.73x | — | 0.57x / 1.50x | 1.078x / 1.023x |
> | v0.11.0 | 63.50 / 25.40 | 3.71x / 9.28x | — | 0.64x / 1.59x | 1.107x / 1.063x |
> | v0.12.0 | 57.66 / 18.56 | 4.09x / 12.70x | — | 0.70x / 2.18x | 1.101x / 1.369x |
> | v0.13.0\* | 53.98 / 17.92 | 4.37x / 13.15x | — | 0.75x / 2.26x | 1.068x / 1.036x |
> | v0.14.0\* | 41.05 / 16.96 | 5.74x / 13.90x | — | 0.98x / 2.38x | 1.315x / 1.057x |
> | v0.14.1\* | 34.64 / 15.72 | 6.80x / 14.99x | — | 1.17x / 2.57x | 1.185x / 1.079x |
> | v0.15.0\* | 33.93 / 15.15 | 6.95x / 15.56x | — | 1.19x / 2.67x | 1.021x / 1.038x |
> | v0.15.1\* | 31.09 / 14.07 | 7.58x / 16.75x | — | 1.30x / 2.87x | 1.091x / 1.077x |
> | v0.16.0\* | 28.53 / 12.45 | 8.26x / 18.93x | — | 1.42x / 3.25x | 1.090x / 1.130x |
> | **v0.17.0**\* | **28.19 / 12.33** | **8.36x / 19.12x** | **—** | **1.43x / 3.28x** | **1.012x / 1.010x** |
>
> **AMD (m8a, x86-64/AVX-512)**
>
> | release | wall_s | vs bwa | vs bwa-mem2 | vs minibwa | vs prev |
> |---|---:|---:|---:|---:|---:|
> | bwa | 208.17 | 1.00x | 0.45x | 0.15x | — |
> | bwa-mem2 | 93.06 | 2.24x | 1.00x | 0.35x | — |
> | minibwa | 32.24 | 6.46x | 2.89x | 1.00x | — |
> | v0.2.1 | 58.05 | 3.59x | 1.60x | 0.56x | 1.603x |
> | v0.2.2 | 55.77 | 3.73x | 1.67x | 0.58x | 1.041x |
> | v0.3.0 | 52.80 | 3.94x | 1.76x | 0.61x | 1.056x |
> | v0.4.0 | 54.20 | 3.84x | 1.72x | 0.59x | 0.974x |
> | v0.5.0 | 54.02 / 23.43 | 3.85x / 8.88x | 1.72x / 3.97x | 0.60x / 1.38x | 1.003x |
> | v0.6.0 | 52.60 / 23.14 | 3.96x / 9.00x | 1.77x / 4.02x | 0.61x / 1.39x | 1.027x / 1.013x |
> | v0.7.0 | 50.80 / 23.39 | 4.10x / 8.90x | 1.83x / 3.98x | 0.63x / 1.38x | 1.035x / 0.989x |
> | v0.8.0 | 46.06 / 19.14 | 4.52x / 10.88x | 2.02x / 4.86x | 0.70x / 1.68x | 1.103x / 1.222x |
> | v0.9.0 | 46.03 / 19.03 | 4.52x / 10.94x | 2.02x / 4.89x | 0.70x / 1.69x | 1.001x / 1.006x |
> | v0.10.0 | 45.16 / 19.10 | 4.61x / 10.90x | 2.06x / 4.87x | 0.71x / 1.69x | 1.019x / 0.996x |
> | v0.11.0 | 43.78 / 18.16 | 4.75x / 11.46x | 2.13x / 5.12x | 0.74x / 1.78x | 1.032x / 1.052x |
> | v0.12.0 | 34.95 / 14.98 | 5.96x / 13.90x | 2.66x / 6.21x | 0.92x / 2.15x | 1.253x / 1.212x |
> | v0.13.0\* | 30.84 / 13.94 | 6.75x / 14.93x | 3.02x / 6.68x | 1.05x / 2.31x | 1.133x / 1.075x |
> | v0.14.0\* | 24.71 / 12.86 | 8.42x / 16.19x | 3.77x / 7.24x | 1.30x / 2.51x | 1.248x / 1.084x |
> | v0.14.1\* | 22.90 / 12.68 | 9.09x / 16.42x | 4.06x / 7.34x | 1.41x / 2.54x | 1.079x / 1.014x |
> | v0.15.0\* | 21.43 / 11.62 | 9.71x / 17.91x | 4.34x / 8.01x | 1.50x / 2.77x | 1.069x / 1.091x |
> | v0.15.1\* | 19.70 / 11.01 | 10.57x / 18.91x | 4.72x / 8.45x | 1.64x / 2.93x | 1.088x / 1.055x |
> | v0.16.0\* | 17.21 / 9.78 | 12.10x / 21.29x | 5.41x / 9.52x | 1.87x / 3.30x | 1.145x / 1.126x |
> | **v0.17.0**\* | **16.19 / 9.69** | **12.86x / 21.48x** | **5.75x / 9.60x** | **1.99x / 3.33x** | **1.063x / 1.009x** |
>
> `vs prev` is the release-over-release speedup (`prev_wall / this_wall`, `>1` = faster) vs the previous release on this same host, `stock / --fast`. Since v0.16.0 every bwa-mem3 release in these tables, old ones included, writes BAM natively (`--bam=0`) while `bwa`, `bwa-mem2` and `minibwa` pipe SAM through `samtools view -u`; earlier tables sent every arm through SAM, so `vs prev` stays like-for-like but the comparator ratios are not comparable across that change. Since v0.17.0 every bwa-mem3 release whose Makefile supports link-time optimization (LTO), v0.15.1 onward, is built with it, old ones included; releases before v0.15.1 predate it and are built without. Every row was re-measured, so neither absolute times nor ratios are comparable with earlier versions of this table. The first release's predecessor is upstream `bwa-mem2` — bwa-mem3 is its successor — so v0.2.1's `vs prev` is its speedup over bwa-mem2 (blank on ARM, where upstream has no build).
>
> Version pins: `bwa` 0.7.19 · `bwa-mem2` v2.2.1 · `minibwa` commit [`d6d9f87d`](https://github.com/lh3/minibwa) (`minibwa-0.7`). "ARM" = Graviton4 m8g (arm64/NEON, no SMT); "x86" = AMD m8a (x86-64/AVX-512, no SMT). Both are the general-purpose siblings of the c8g/c8a hosts used through v0.12.0 — same CPU family, same core count, no SMT, but 4 GiB/vCPU instead of 2 so every historical arm fits in memory — so absolute times are not comparable with earlier versions of this table, only ratios within it. (The x86 arm replaced an earlier Intel c7i arm, which ran 16 vCPUs over 8 physical cores under 2-way SMT and so wasn't a real core-for-core match for Graviton's 16 real cores.) No ARM `bwa-mem2` build exists, hence the blank cells there. Every arm for a given arch ran interleaved on one fixed on-demand host — 5 rounds for v0.13.0 through v0.17.0, 2 for the rest, median wall-clock shown — so these are same-host comparisons, not medians pooled across separate runs. `—` means the release predates the comparator or predates `--fast`. \* v0.13.0 through v0.17.0 rows use the stride-2 suffix-array index described above; all other rows use the stock stride-8 index. Regenerated at each release; see [Benchmarks](https://bwa-mem3.readthedocs.io/en/latest/performance/benchmarks.html).
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

`bwa-mem3 aln` aligns ultrashort, ancient or degraded DNA with bwa's backtracking algorithm. It fuses `bwa aln` and `bwa samse`/`bwa sampe` into one command and writes SAM that is byte-identical to bwa 0.7.19's, apart from the `@PG` line, on the generated reads the parity checks use, in CI on Linux x86_64 (AVX2 build) and Linux arm64 (NEON), against a stock bwa built without floating-point contraction (a default stock build on arm64 can differ at rare rounding boundaries; see [aln](https://bwa-mem3.readthedocs.io/en/latest/cli/aln.html)). `index`, `mem`, `shm`, `re-sa` and `version` are the others.

## Quick links

- [Benchmarks](https://bwa-mem3.readthedocs.io/en/latest/performance/benchmarks.html) — published results for every release, and how they are measured
- [bwa-mem3-rs](https://github.com/fg-labs/bwa-mem3-rs) — Rust bindings for bwa-mem3
- [bioconda recipe](https://github.com/bioconda/bioconda-recipes/tree/master/recipes/bwa-mem3) — conda package on bioconda
- [fgumi](https://github.com/fulcrumgenomics/fgumi) — UMI-aware consensus and deduplication
- [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) — upstream project

## License

bwa-mem3 is primarily MIT-licensed; the `aln` code derived from bwa is licensed under Apache-2.0. Binary distributors must include both the root `LICENSE` and `src/aln/LICENSE` with their binaries. See the [License page](https://bwa-mem3.readthedocs.io/en/latest/reference/license.html) in the docs.

## Citation

Please cite the bwa-mem2 paper (Vasimuddin Md et al., IPDPS 2019). See the [Citation page](https://bwa-mem3.readthedocs.io/en/latest/reference/citation.html) for BibTeX.

## Issues / contributing

File [issues](https://github.com/fg-labs/bwa-mem3/issues) and [pull requests](https://github.com/fg-labs/bwa-mem3/pulls) on [fg-labs/bwa-mem3](https://github.com/fg-labs/bwa-mem3).
