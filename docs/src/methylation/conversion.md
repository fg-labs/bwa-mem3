# Conversion Details (C→T, G→A)

Bisulfite sequencing relies on chemical conversion of unmethylated cytosines to
uracil (read as thymine after PCR). `bwa-mem3 --meth` models this with an
in-memory read transformation applied to every read before the alignment kernel
sees the bases.

## What gets converted

Paired-end bisulfite reads follow a strand convention:

- **R1** (read 1): every `C` in the base sequence is replaced with `T`.
  This models the OT (original top) and CTOB (complementary to original bottom)
  strands as they appear after bisulfite treatment and PCR.

- **R2** (read 2): every `G` in the base sequence is replaced with `A`.
  This models the OB (original bottom) and CTOT strands.

Single-end mode uses the R1 (`C→T`) rule for all reads.

The doubled reference built by `bwa-mem3 index --meth` (or bwameth.py) contains
two projections of each chromosome:

- `f`-prefixed contigs (e.g. `fchr1`): the chromosome with every `C` replaced
  by `T`.
- `r`-prefixed contigs (e.g. `rchr1`): the reverse complement of the chromosome
  with every `G` replaced by `A`.

Converted R1 reads are therefore alignable to `f`-prefixed contigs and
converted R2 reads to `r`-prefixed contigs. The contig prefix records which
strand hypothesis was used and feeds the Bismark `XG:Z` tag directly
(`CT` for `f`-prefixed / OT, `GA` for `r`-prefixed / OB).

## Where conversion happens

Read conversion runs inside `src/fastmap.cpp` in the `meth_mode` ingest block,
immediately after sequence parsing and before any alignment work. The
transformation is applied to the in-memory `bseq1_t.seq` buffer; the original
FASTQ file is never rewritten.

Before the bases are modified, the original sequence is recorded in the read's
comment buffer as:

```text
YS:Z:<l_seq bases>\tYC:Z:<direction>
```

where `<direction>` is `CT` for R1 (C→T) and `GA` for R2 (G→A). These
fields pass through the alignment kernel untouched and serve two
internal purposes at BAM-write time: `YS:Z` is the source for SEQ
restoration (see next section), and `YC:Z` is the source for the
emitted `XR:Z:` Bismark tag. They are **not** emitted to the output
BAM — `bam_writer.cpp` suppresses them under `--meth`. See
[SAM tags: XR, XG, XM](tags.md) for the per-tag output reference.

## Sequence restoration in the BAM SEQ field

Methylation callers such as MethylDackel identify methylated cytosines by
examining the BAM SEQ field at each CpG site. They need to see real `C`/`T`
bases — not the uniformly-converted `T`/`A` bases that were used for alignment.

`meth_mem_aln_to_bam` (in `src/meth_bam.cpp`) restores the original
sequence from the internal `YS:Z` carrier on `bseq1_t.comment` before
writing the BAM record (the carrier itself is suppressed at BAM emit
by `bam_writer.cpp` under `--meth`, so it never reaches the output
file):

1. The `YS:Z:` payload is located at the start of the `bseq1_t.comment` field
   (offset `+5` past the `YS:Z:` header bytes).
2. For forward-aligned records (`!p.is_rev`), the pre-c2t bases are copied
   directly into the BAM SEQ buffer.
3. For reverse-aligned records (`p.is_rev`), the bases are reverse-complemented
   using the standard `TGCAN` table before being placed in SEQ.
4. If `YS:Z:` is absent (e.g. when running with an external c2t converter that
   does not emit it), the code falls back to the converted sequence in `s->seq`,
   with the same RC flip logic.

> **Warning — Soft-clip and supplementary trimming**
>
> When computing the SEQ range for supplementary alignments, the `qb`/`qe`
> boundaries account for soft-clip or hard-clip operations at the CIGAR ends.
> The YS:Z: restoration applies over the same trimmed range so SEQ length
> always matches the emitted CIGAR.

## QUAL field handling

The QUAL field is taken directly from the original FASTQ (`bseq1_t.qual`) over
the same `[qb, qe)` range and is never modified by the c2t process. Quality
scores correspond to the original base calls, not the converted ones.

## Relationship to the reference index

`bwa-mem3 index --meth ref.fa` writes `ref.fa.bwameth.c2t`, which applies the
same C→T / G→A projection to the reference sequence. The resulting file is
compatible with what `bwameth.py index` produces, so the same doubled-reference
FASTA can be used interchangeably with either tool across tested versions.

---

**See also:**
[Overview](overview.md) ·
[SAM tags: XR, XG, XM](tags.md) ·
[Interop with external bwameth.py c2t](external-c2t.md) ·
[User Guide → Indexing the reference](../user-guide/indexing.md) ·
[Best Practices → Methylation defaults](../best-practices/methylation.md)
