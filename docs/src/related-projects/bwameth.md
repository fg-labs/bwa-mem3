# bwameth.py

bwameth.py is a Python script written by Brent Pedersen that implements bisulfite
sequencing (BS-Seq) alignment using the in-silico three-letter genome approach. It
converts all cytosines to thymines in both the reference and the reads (C-to-T on the
forward strand, G-to-A on the reverse), aligns the converted sequences with bwa-mem (or
optionally bwa-mem2), and then recovers the original read sequence from the aligner's
tag output to tabulate methylation. bwameth.py supports single-end and paired-end reads
from the directional bisulfite protocol and is published at
<https://arxiv.org/abs/1401.1129>.

## When you'd use it

Use bwameth.py when you need a battle-tested, community-supported bisulfite aligner that
runs on top of the standard bwa-mem or bwa-mem2 you have already installed, and when you
prefer a Python wrapper over a self-contained binary. It also remains the reference for
downstream tabulation tools such as [MethylDackel](https://github.com/dpryan79/MethylDackel)
and SNP callers such as [biscuit](https://github.com/huishenlab/biscuit) that expect the
bwameth.py output format. For the actual methylation tabulation and variant calling steps,
bwameth.py's author recommends those dedicated tools rather than the tabulation utilities
bundled with the original script.

## How it relates to bwa-mem3

`bwa-mem3 mem --meth` is a single-binary drop-in replacement for the bwameth.py
alignment pipeline. It inlines the C-to-T and G-to-A conversion, runs the bwa-mem3
alignment engine (with all of its correctness fixes and SIMD speedups), rewrites the
`@SQ` headers to collapse the per-strand contig pairs back to canonical chromosome names,
emits Bismark-compatible `XR:Z` / `XG:Z` / `XM:Z` auxiliary tags, and writes a
`@PG ID:bwa-mem3-meth` header. The bwameth.py-style chimera QC heuristic is
available via `--chimera-qc` (off by default — Bismark behavior).
The [Methylation Reference](../methylation/overview.md) section documents the full
implementation in detail, including the Bismark `XR:Z` / `XG:Z` / `XM:Z` tags and
the `--set-as-failed` / `--chimera-qc` flags.

> **Tip — Interop with the bwameth.py c2t step**
>
> If your pipeline already performs its own C-to-T conversion before alignment, see
> [Interop with external bwameth.py c2t](../methylation/external-c2t.md) for how to pass
> pre-converted reads to `bwa-mem3 mem --meth` without double-conversion.

## Links

- GitHub: <https://github.com/brentp/bwa-meth>
- Paper: <https://arxiv.org/abs/1401.1129>
- License: MIT

---

**See also:**
[Methylation Reference: Overview](../methylation/overview.md) ·
[Quick start: methylation alignment](../getting-started/quick-meth.md) ·
[Best Practices — Methylation defaults](../best-practices/methylation.md) ·
[Interop with external bwameth.py c2t](../methylation/external-c2t.md)
