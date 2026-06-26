# bwa-mem2 (upstream)

bwa-mem2 is the direct predecessor of bwa-mem3 and the project from which the
bwa-mem3 fork is derived. It was created at Intel's Parallel Computing Lab by
Vasimuddin Md and Sanchit Misra to accelerate the alignment algorithm originally
written by Heng Li in [bwa](https://github.com/lh3/bwa). bwa-mem2 achieves a
1.3–3.1x throughput improvement over the original bwa-mem by replacing key inner
loops with vectorised implementations (SSE4.1, SSE4.2, AVX2, and AVX-512) and by
switching to a more compact FM-index encoding. Its output is identical to bwa-mem
at the alignment level, and it is distributed under the MIT license.

## Lineage

The bwa alignment family has evolved through three generations, each building on
the last:

1. **bwa** — Written by Heng Li. Established the BWA-MEM algorithm, the SAM output
   format conventions, and the `.bwt` / `.pac` / `.ann` / `.amb` index layout.
2. **bwa-mem2** (Vasimuddin et al., Intel) — Replaced scalar inner loops with SIMD
   kernels; introduced the compact `.bwt.2bit.64` and `.0123` index formats;
   retained full output compatibility with bwa-mem. (bwa-mem3 reads the same
   formats but no longer **builds** the unpacked `.0123` by default — it
   pac-fetches reference bases from `.pac`. If you want to share a bwa-mem3 index
   with bwa-mem2, which loads `.0123` directly, build it with
   `bwa-mem3 index --emit-unpacked-ref`.)
3. **bwa-mem3** (Fulcrum Genomics fork) — Carries correctness fixes, performance
   improvements, new features (bisulfite alignment, mimalloc, ARM Neon), and
   expanded architecture support on top of the bwa-mem2 codebase. See
   [What's Different from bwa-mem2](../whats-different/overview.md) for the full
   change catalog.

## When you'd use it

Use bwa-mem2 directly when you need a stable, widely validated aligner with
precompiled binaries available via Bioconda and the project's GitHub releases page,
and when you do not require the features or fixes that bwa-mem3 adds. bwa-mem2 is
also the right choice when you are working in an environment where the bwa-mem3 fork
has not yet been validated against your specific reference or sequencing library type.

## How it relates to bwa-mem3

bwa-mem3 tracks bwa-mem2's `master` branch and periodically rebases fork-carried
commits on top of upstream changes. The [What's Different](../whats-different/overview.md)
section documents every divergence between the two projects, and the
[Upstream PR status](../whats-different/upstream-prs.md) page tracks which bwa-mem3
changes have been proposed back to bwa-mem2. The goal is to keep the fork divergence
minimal and to upstream as many fixes as practical.

## Links

- GitHub: <https://github.com/bwa-mem2/bwa-mem2>
- Citation: Vasimuddin Md, Sanchit Misra, Heng Li, Srinivas Aluru. "Efficient
  Architecture-Aware Acceleration of BWA-MEM for Multicore Systems." *IEEE IPDPS 2019.*
- License: MIT (with third-party components under their respective licenses)

---

**See also:**
[What's Different from bwa-mem2](../whats-different/overview.md) ·
[Upstream PR status](../whats-different/upstream-prs.md) ·
[bwa-mem3-bench](bwa-mem3-bench.md) ·
[Citation](../reference/citation.md)
