# Conversion Details (C→T, G→A)

Bisulfite sequencing converts unmethylated cytosines to uracil (read as thymine
after PCR). `bwa-mem3 --meth` uses a C→T / G→A projection to **find seeds**, but
— unlike a classic 3-letter aligner — it then scores and reports against the
original 4-letter reference. This page describes what gets projected, where, and
how the original bases come back for scoring and output.

## What gets projected (for seeding only)

Paired-end bisulfite reads follow a strand convention:

- **R1** (read 1): every `C` is replaced with `T` (models the OT / CTOB strands).
- **R2** (read 2): every `G` is replaced with `A` (models the OB / CTOT strands).

Single-end mode uses the R1 (`C→T`) rule for all reads.

This projection is applied **only to the copy of the read used for seeding**. The
`.meth` doubled seed index built by `bwa-mem3 index --meth` holds two projections
of each chromosome:

- `f`-prefixed contigs (e.g. `fchr1`): the chromosome with every `C` → `T`.
- `r`-prefixed contigs (e.g. `rchr1`): the reverse-complement strand with every
  `G` → `A`.

Projected R1 reads seed against `f`-prefixed contigs and projected R2 reads
against `r`-prefixed contigs. The contig prefix records the strand hypothesis
(`f` → OT, `r` → OB), which both selects the per-strand scoring matrix and feeds
the Bismark `XG:Z` tag (`CT` for OT, `GA` for OB).

> **Key difference from bwameth / 3-letter aligners**
>
> A classic 3-letter aligner *also* converts the reference and the read, then
> scores in collapsed space — so a real C/T or G/A variant is invisible. Here the
> projection is used only to locate seeds; **extension and scoring run on the
> original bases against the original reference** (next sections). That is what
> lets `--meth-scoring genomic` tell a real variant apart from a conversion.

## Where projection happens

Seed projection runs inside `src/fastmap.cpp` in the `meth_mode` ingest block,
right after sequence parsing. It writes the projected bases into the in-memory
seeding buffer; the original FASTQ is never rewritten.

Before projecting, the original sequence is preserved on a first-class per-read
field, `bseq1_t.meth_orig_seq`. (For interoperability with an external c2t
converter that does not populate that field, the original bases can also be
carried on the read comment as `YS:Z:<l_seq bases>\tYC:Z:<direction>`, where
`<direction>` is `CT` for R1 and `GA` for R2.) These carriers are internal: they
are **not** emitted to the output BAM.

## Scoring against the original reference

After seeds are remapped to original coordinates with their OT/OB hypothesis,
`bwa-mem3` extends and scores the **original** read against the **original**
4-letter reference window using the per-strand asymmetric matrix:

- **OT** frees ref-`C` × read-`T` (the expected unmethylated C→T conversion).
- **OB** frees ref-`G` × read-`A` (the expected bottom-strand G→A conversion).

Under `--meth-scoring collapsed` (default) the mirror cell is freed too (ref-`T` ×
read-`C`, ref-`A` × read-`G`), so C/T and G/A are interchangeable and placement
matches bwameth. Under `--meth-scoring genomic` only the conversion direction is
freed, so a real variant stays a mismatch. See
[Overview → `--meth-scoring`](overview.md#two-scoring-modes---meth-scoring).

The seed's own ungapped score is recomputed in the same matrix (not assumed to be
a perfect `len × match`), so under `--meth-scoring genomic` a seed-internal C/T or
G/A variant correctly lowers the alignment score, `AS`, and `MAPQ`. Under
`--meth-scoring collapsed` the mirror cell is freed, so such a variant is scored
as a conversion and does not penalize placement.

## Sequence restoration in the BAM SEQ field

Methylation callers (MethylDackel, Bismark tools) read the BAM `SEQ` field to see
real `C`/`T` bases, not the projected `T`/`A`. `meth_mem_aln_to_bam` (in
`src/meth_bam.cpp`) restores the original bases before writing each record:

1. The original bases come from `bseq1_t.meth_orig_seq` (the first-class field),
   falling back to the `YS:Z` comment carrier only when that field is absent.
2. For forward-aligned records (`!p.is_rev`), the original bases are copied
   directly into `SEQ`.
3. For reverse-aligned records (`p.is_rev`), they are reverse-complemented with
   the standard `TGCAN` table.
4. If neither carrier is available (e.g. an external c2t converter that emits
   neither), the code falls back to the seeding buffer in `s->seq`, with the same
   RC flip.

> **Warning — Soft-clip and supplementary trimming**
>
> When computing the SEQ range for supplementary alignments, the `qb`/`qe`
> boundaries account for soft-clip / hard-clip operations at the CIGAR ends. The
> restoration applies over the same trimmed range, so SEQ length always matches
> the emitted CIGAR.

## QUAL field handling

The QUAL field is taken directly from the original FASTQ (`bseq1_t.qual`) over the
same `[qb, qe)` range and is never modified. Quality scores correspond to the
original base calls.

## Relationship to the reference index

`bwa-mem3 index --meth ref.fa` writes **two** indexes:

- the normal 4-letter index at the bare prefix (`ref.fa.0123`, `.amb`, `.ann`,
  `.bwt.2bit.64`, `.pac`) — used for scoring/extension against the original
  reference, and
- the converted **seed** index `ref.fa.meth.*` (built over a per-strand-converted
  FASTA `ref.fa.meth.fa` with the `f`/`r` doubled contigs) — used only for
  seeding.

This dual-index layout differs from bwameth.py, which builds a single
`ref.fa.bwameth.c2t` doubled reference and aligns entirely against it. A legacy
bwameth `.bwameth.c2t` index cannot be reused directly — rebuild with
`index --meth` (see [Migrating from bwameth.py c2t](external-c2t.md)).

---

**See also:**
[Overview](overview.md) ·
[SAM tags: XR, XG, XM](tags.md) ·
[Migrating from bwameth.py c2t](external-c2t.md) ·
[User Guide → Indexing the reference](../user-guide/indexing.md) ·
[Best Practices → Methylation defaults](../best-practices/methylation.md)
