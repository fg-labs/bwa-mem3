# bwameth.py Drop-In Mapping

`bwa-mem3 --meth` is designed so that, in its **default `collapsed` mode**, read
*placement* matches the bwameth.py pipeline for the standard case, while emitting
the Bismark tag set methylation callers expect. This page explains what is the
same, what differs, and where the two approaches diverge by design.

> **Important — placement drop-in, not byte-identical**
>
> `collapsed` reproduces bwameth's **placement** (where reads map and their
> primary/MAPQ behavior) because both treat C/T and G/A as interchangeable. It is
> **not** a byte-for-byte reproduction: bwa-mem3 scores against the original
> 4-letter reference rather than in collapsed space, so scores and some CIGARs can
> differ at the margins, and the tag schema is Bismark (`XR`/`XG`/`XM`), not
> bwameth (`YS`/`YC`/`YD`). The opt-in `genomic` mode diverges from bwameth on
> purpose (it penalizes real variants).

## Command comparison

### bwameth.py pipeline (multi-step)

```sh
# Step 1: build a single doubled (c2t) reference
bwameth.py index ref.fa                # writes ref.fa.bwameth.c2t + FMI

# Step 2: align (bwameth.py converts reads, calls bwa/bwa-mem2, post-processes)
bwameth.py map --bwa-mem2 -t 16 ref.fa R1.fq.gz R2.fq.gz \
  | samtools sort -o out.bam
samtools index out.bam
```

### bwa-mem3 --meth (single binary)

```sh
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

**Variant-aware option.** `--meth-scoring genomic` scores real C/T and G/A
variants as mismatches, so a single BAM supports both methylation calling and
variant calling — something a collapsed-space aligner cannot produce.

**Inline BAM post-processing.** Header rewriting, Bismark `XR`/`XG`/`XM` tags,
opt-in chimera QC (`--chimera-qc`), and QC-fail propagation happen in the same
pass. Output is uncompressed BAM (`wb0`) that `samtools sort` reads natively.

**bwameth-aligned defaults (collapsed).** `--meth-scoring collapsed` applies
`-B 2 -L 10 -U 100 -T 40 -M -C`, mirroring bwameth's `bwa mem -T 40 -B 2 -L 10
-CM` (plus `-U 100` for paired-end). `genomic` uses the same set but keeps
`-B 4`. All parameters can be overridden.

## What stays the same (collapsed mode)

The output BAM carries the standard methylation tag set, flags, and SEQ
representation. The `@PG` provenance line and the tag schema intentionally differ:

| Field | bwameth.py | bwa-mem3 --meth |
|-------|-----------|-----------------|
| `@SQ` headers | One per real chromosome | One per real chromosome |
| Read placement (collapsed) | reference | Matches at the standard case |
| Methylation aux tags | `YS:Z`, `YC:Z`, `YD:Z` | `XR:Z`, `XG:Z`, `XM:Z` (Bismark) |
| `@PG` | `ID:bwameth` | `ID:bwa-mem3-meth` |
| Chimera QC threshold | Longest M < 44% of read | Same (44%), opt-in via `--chimera-qc` |
| Chimera QC flags | `0x200`, clear `0x2`, MAPQ ≤ 1 | Same |
| SEQ field | Pre-conversion bases (RC-flipped when `is_rev`) | Same |
| `NM`/`MD` | Collapsed (conversions and real variants both hidden) | Conversions hidden; real variants hidden in `collapsed`, **shown in `genomic`** |

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
bwameth-compatible placement (`collapsed`) or variant-aware scoring (`genomic`).

---

**See also:**
[Overview](overview.md) ·
[Conversion details](conversion.md) ·
[SAM tags: XR, XG, XM](tags.md) ·
[Chimera QC and header rewriting](post-processing.md) ·
[Related Projects: bwameth.py](../related-projects/bwameth.md)
