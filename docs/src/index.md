<h1 class="bm3-title"><span class="bm3-base">bwa-mem</span><span class="bm3-num">3</span></h1>

**A faster, more correct, drop-in replacement for `bwa mem` and `bwa-mem2`.**

If you align short reads with bwa or bwa-mem2 today, bwa-mem3 will give you
the same answers — only quicker, with fewer rough edges, and with first-class
support for things you used to need a wrapper script for.

## Why bwa-mem3

- **Drop in, go faster.** Same algorithm, same outputs, same flags as
  bwa-mem2 — but consolidated mapping speedups, a memory-bounded index
  builder, batched header ingestion, and a tuned allocator add up to
  measurable wall-clock wins on real workloads.
- **Methylation in one binary.** A `--meth` flag turns bwa-mem3 into a
  drop-in replacement for the entire `bwameth.py` pipeline. No Python, no
  inline conversion script, no separate post-processing step. One
  `bwa-mem3 index --meth ref.fa`, one `bwa-mem3 mem --meth ref.fa
  R1.fq R2.fq`, done — header collapsed, tags emitted, chimeras flagged.
- **Stage the index once, align many.** A `bwa-mem3 shm` subcommand
  pins the FM-index in shared memory so back-to-back runs on the same host
  skip the 28 GB read every time.
- **Correctness fixes upstream haven't merged yet.** Tabs in `-R`,
  151+ bp reads, AVX-512 mate-rescue, kswv `score2` plateau across
  NEON/AVX2/AVX-512BW, mem_sam_pe proper-pair flag — every fix tracked
  back to the upstream PR or issue that found it.
- **Architecture-aware out of the box.** SSE4.1, SSE4.2, AVX, AVX2,
  AVX-512BW, and ARM64/NEON. One binary per platform; the dispatcher
  picks the right tier for your CPU in process at startup.

## Get started in 30 seconds

```bash
git clone --recursive https://github.com/fg-labs/bwa-mem3
cd bwa-mem3 && make
./bwa-mem3 index ref.fa
./bwa-mem3 mem -t 16 ref.fa R1.fq.gz R2.fq.gz \
  | samtools sort -@ 8 -o out.bam
```

> **Tip — Emit BAM directly**
>
> For production pipelines, add `--bam=0` to skip the SAM text round-trip
> entirely. See [Best Practices: Output format](best-practices/output-format.md).

## Where to start

- **[Installation](getting-started/installation.md)** — Build from source
  (Bioconda is on the way).
- **[Quick start: align paired-end FASTQs](getting-started/quick-align.md)**
  — Two commands to your first alignment.
- **[Quick start: methylation](getting-started/quick-meth.md)** — The
  one-binary `bwameth.py` replacement, in two commands.
- **[Best Practices](best-practices/build.md)** — The five things that
  actually move the needle for production runs.
- **[What's different from bwa-mem2](whats-different/overview.md)** —
  Every fix and feature, with upstream cross-references.

## What's in this book

- **[Getting Started](getting-started/installation.md)** — Install and run
  your first alignment.
- **[User Guide](user-guide/indexing.md)** — Indexing, alignment, output,
  threading, allocator notes.
- **[Performance](performance/overview.md)** — Where the speed comes from
  and how to get more.
- **[Best Practices](best-practices/build.md)** — Build, run, and deploy
  recommendations.
- **[CLI Reference](cli/overview.md)** — Every flag, auto-captured from
  `--help`.
- **[Methylation Reference](methylation/overview.md)** — `--meth` mode in
  full.
- **[What's Different from bwa-mem2](whats-different/overview.md)** — The
  full changelog, by category.
- **[Developer Guide](developer-guide/building.md)** — Build matrix, SIMD
  dispatch, regression tests, contributing.
- **[Related Projects](related-projects/bwa-mem3-bench.md)** —
  bwa-mem3-bench, bwa-mem3-rs, fgumi, bwa-mem2 upstream.
- **[Reference](reference/glossary.md)** — Glossary, citation, license,
  changelog.

---

bwa-mem3 is a derivative of [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2)
maintained by [Fulcrum Genomics](https://www.fulcrumgenomics.com). MIT
licensed. See [License](reference/license.md) and
[Citation](reference/citation.md).
