# bwameth.py Drop-In Mapping

`bwa-mem3 --meth` is designed so that, in its **default `collapsed` mode**, read
*placement* closely tracks the bwameth.py pipeline for the standard case (with a
small, bounded divergence — see the callout below), while emitting the Bismark
tag set methylation callers expect. This page explains what is the same, what
differs, and where the two approaches diverge by design.

> **Important — placement drop-in, *not* byte-identical, and *not* drift-free**
>
> `collapsed` **approximates** bwameth's placement (where reads map and their
> primary/MAPQ behavior) because both treat C/T and G/A as interchangeable — but
> it is **neither** a byte-for-byte reproduction **nor** drift-free. bwa-mem3
> scores against the original 4-letter reference rather than in collapsed space,
> so a small but nonzero fraction of records differ from bwameth.py in **`POS`**,
> **`CIGAR`**, or **`MAPQ`** — on the order of ~1% of records on typical
> WGBS/EM-seq, with true mapped-position (`POS`) changes affecting a smaller
> subset (a few tenths of a percent). The 4-letter scoring path **widens** this
> versus a pure collapsed-space aligner — and versus the older 3-letter `--meth`
> releases — and the variant-aware modes (`genomic`, `neutral`) diverge further
> on purpose (they penalize real variants). The tag schema is also Bismark
> (`XR`/`XG`/`XM`), not bwameth (`YS`/`YC`/`YD`).
>
> **If you are pinned to a specific bwameth release** (e.g. a clinical pipeline
> validated against a bwameth version), treat `collapsed` as a *new aligner* and
> **re-validate against your own bwameth output** — do not assume placement
> equivalence. The divergence is small and bounded, but it is real and it is
> larger than the older 3-letter `--meth` path.

## Command comparison

### bwameth.py pipeline (multi-step)

```bash
# Step 1: build a single doubled (c2t) reference
bwameth.py index ref.fa                # writes ref.fa.bwameth.c2t + FMI

# Step 2: align (bwameth.py converts reads, calls bwa/bwa-mem2, post-processes)
bwameth.py map --bwa-mem2 -t 16 ref.fa R1.fq.gz R2.fq.gz \
  | samtools sort -o out.bam
samtools index out.bam
```

### bwa-mem3 --meth (single binary)

```bash
# Step 1: build a dual index (original index + .meth seed index)
bwa-mem3 index --meth ref.fa           # writes ref.fa.* AND ref.fa.meth.*

# Step 2: align (inline seed projection + 4-letter scoring + post-processing)
bwa-mem3 mem --meth -t 16 ref.fa R1.fq.gz R2.fq.gz \
  | samtools sort -o out.bam
samtools index out.bam
```

The index layouts **differ**. bwameth.py builds one collapsed doubled reference
(`ref.fa.bwameth.c2t` + FMI) and aligns entirely against it. bwa-mem3 builds two
indexes: the normal 4-letter index at the bare prefix (for scoring/extension) and
a converted **seed** index `ref.fa.meth.*` (for seeding only). A legacy bwameth
`.bwameth.c2t` index is not used directly — rebuild with `index --meth` (see
[Migrating from bwameth.py c2t](external-c2t.md)).

## What is gained

**No Python or bwameth.py dependency.** Read seeding, 4-letter scoring, and BAM
post-processing all run inside a single `bwa-mem3` process. One binary, no virtual
environment, no bwameth.py version pinning.

**No intermediate files.** No converted FASTQ is written; the C→T / G→A projection
is applied in-memory to the seeding copy of each read.

**Variant-aware options.** `--meth-scoring genomic` and `--meth-scoring neutral`
(the `--meth=taps` default) both score a real C/T or G/A variant on the mirror cell
as a mismatch, so a single BAM supports both methylation calling and variant
calling — something a collapsed-space aligner cannot produce. (A variant in the
conversion direction itself remains hidden; see
[which real variants stay visible](overview.md#which-real-variants-stay-visible).)
They differ only in what the *conversion* cell scores: a full match under
`genomic`, `0` under `neutral`.

**Inline post-processing.** `@SQ` built from the original reference, Bismark
`XR`/`XG`/`XM` tags, opt-in chimera QC (`--chimera-qc`), and QC-fail propagation
happen in the same pass. Output is SAM text by default and uncompressed BAM
(`wb0`) under `--bam`; `samtools sort` reads either natively.

**bwameth-aligned defaults (collapsed).** `--meth-scoring collapsed` applies
`-B 2 -L 10 -U 100 -T 40 -M -C`, mirroring bwameth's `bwa mem -T 40 -B 2 -L 10
-CM` (plus `-U 100` for paired-end). `genomic` and `neutral` use the same set but
keep `-B 4`. The scoring parameters (`-B`, `-L`, `-U`, `-T`) can be overridden on the
command line, in any position relative to `--meth`. `-M` and `-C` cannot: bwa
has no option that unsets them, so `--meth` applies them unconditionally.

These constants are quoted at bwa's default match score (`-A 1`, what bwameth
runs). Like every other score-derived default, they scale with `-A`: under
`-A 2` the effective values are `-L 20 -U 200 -T 80` and `-B 4` (collapsed) /
`-B 8` (genomic and neutral).

## What stays the same (collapsed mode)

The output BAM carries the standard methylation tag set, flags, and SEQ
representation, and read placement **closely tracks** bwameth at the standard
case — but "stays the same" here means *functionally equivalent*, not
*identical*: as the callout above notes, ~1% of records still differ in
`POS`/`CIGAR`/`MAPQ`. The `@PG` provenance line and the tag schema intentionally
differ:

| Field | bwameth.py | bwa-mem3 --meth |
|-------|-----------|-----------------|
| `@SQ` headers | One per original reference sequence | One per original reference sequence |
| Read placement (collapsed) | reference | Closely tracks at the standard case; ~1% of records differ in `POS`/`CIGAR`/`MAPQ` on typical WGBS/EM-seq, `POS` alone in a few tenths of a percent (see the callout at the top of this page) — re-validate if pinned to a bwameth release |
| Methylation aux tags | `YS:Z`, `YC:Z`, `YD:Z` | `XR:Z`, `XG:Z`, `XM:Z` (Bismark) |
| `@PG` | `ID:bwameth` | `ID:bwa-mem3-meth` |
| Chimera QC threshold | Longest M < 44% of read | Same (44%), opt-in via `--chimera-qc` |
| Chimera QC flags | `0x200`, clear `0x2`, MAPQ ≤ 1 | Same |
| SEQ field | Pre-conversion bases (RC-flipped when `is_rev`) | Same |
| `NM`/`MD` | Collapsed (conversions and real variants both hidden) | Conversions hidden; real variants hidden in `collapsed`, **shown in `genomic` and `neutral`**[^nmmd] |

[^nmmd]: `genomic` and `neutral` both free only the conversion direction — they
    differ solely in what that cell scores (a full match vs. `0`) and both leave
    the *mirror* cell penalised — so both keep a real variant in the *opposite*
    direction (ref `T` × read `C` on OT) visible. A real variant in the
    conversion direction itself (a genuine C→T SNP at a reference C) is
    indistinguishable from a conversion in a single read under any mode — that
    aliasing is resolved downstream at the pileup, not by the aligner. See
    [which real variants stay visible](overview.md#which-real-variants-stay-visible).

bwa-mem3 emits the **Bismark-compatible** `XR:Z` / `XG:Z` / `XM:Z` tag set rather
than bwameth's `YS:Z` / `YC:Z` / `YD:Z`, so output is directly consumable by
`bismark_methylation_extractor`, methylKit, methtuple, DMRfinder, and epialleleR
in addition to MethylDackel and biscuit. Tools that expect `YS`/`YC`/`YD` must be
pointed at the corresponding `XR`/`XG` (and per-base `XM`) tags.

## When to prefer bwameth.py

If your workflow requires bwameth.py-specific features (e.g. `bwameth.py
markduplicates` or non-standard post-processors), or strict byte-for-byte
reproduction of a bwameth release, continue using bwameth.py. `bwa-mem3 --meth`
targets the indexing + alignment + standard post-processing path, with
bwameth-compatible placement (`collapsed`) or variant-aware scoring (`genomic`,
`neutral`).

---

**See also:**
[Overview](overview.md) ·
[Conversion details](conversion.md) ·
[SAM tags: XR, XG, XM](tags.md) ·
[Chimera QC and header construction](post-processing.md) ·
[Related Projects: bwameth.py](../related-projects/bwameth.md)
