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

`bwa-mem3 mem --meth` is a single-binary alignment pipeline that, in its default
`collapsed` mode, reproduces bwameth.py's read *placement* (it is a placement
drop-in, not a byte-for-byte clone). The key difference is *where it scores*:
where bwameth.py converts both reads and reference to 3-letter space and aligns
there, bwa-mem3 uses the 3-letter projection only to **find seeds**, then extends
and scores against the **original 4-letter reference**. That enables two further
variant-aware modes — `--meth-scoring genomic`, which scores the conversion cell as
a full match, and `--meth-scoring neutral` (the `--meth=taps` default), which scores
it `0` so a sparse conversion is tolerated but not rewarded. Both keep a real C/T or
G/A variant in the direction *opposite* the conversion as a mismatch (visible in
`NM`/`MD`), something a collapsed-space aligner cannot do; a variant in the conversion
direction itself is indistinguishable from a conversion under either design
([why](../methylation/overview.md#which-real-variants-stay-visible)).

It writes `@SQ` headers directly from the original reference, so output carries
the original reference contig names, emits Bismark-compatible `XR:Z` / `XG:Z` / `XM:Z`
auxiliary tags, and writes a `@PG ID:bwa-mem3-meth` header. The bwameth.py-style
chimera QC heuristic is available via `--chimera-qc` (off by default — Bismark
behavior). The [Methylation Reference](../methylation/overview.md) documents the
full implementation, including the three
[`--meth-scoring` modes](../methylation/overview.md#three-scoring-modes---meth-scoring),
the Bismark tags, and the `--set-as-failed` / `--chimera-qc` flags.

> **Note — external c2t interop was removed**
>
> Because scoring runs on the original bases, `bwa-mem3 mem --meth` can no longer
> consume pre-converted reads or a `.bwameth.c2t` reference. Pass raw FASTQ and the
> original `ref.fa` prefix; see
> [Migrating from bwameth.py c2t](../methylation/external-c2t.md).

## Links

- GitHub: <https://github.com/brentp/bwa-meth>
- Paper: <https://arxiv.org/abs/1401.1129>
- License: MIT

---

**See also:**
[Methylation Reference: Overview](../methylation/overview.md) ·
[Quick start: methylation alignment](../getting-started/quick-meth.md) ·
[Best Practices — Methylation defaults](../best-practices/methylation.md) ·
[Migrating from bwameth.py c2t](../methylation/external-c2t.md)
