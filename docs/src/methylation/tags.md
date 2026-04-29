# SAM Tags: YS, YC, YD

`bwa-mem3 --meth` emits three methylation-specific auxiliary tags that carry
the information downstream methylation callers need. Two of these (`YS:Z:` and
`YC:Z:`) are set during FASTQ ingest and pass through the alignment kernel
unchanged. The third (`YD:Z:`) is set during BAM post-processing based on the
contig name of the alignment.

## Tag reference

### `YS:Z:` — original (pre-conversion) sequence

| Property | Value |
|----------|-------|
| Type | `Z` (NUL-terminated string) |
| Length | Equal to `l_seq` (full read length) |
| Set by | FASTQ ingest (`src/fastmap.cpp` meth_mode block) |
| Emitted on | All records (mapped and unmapped) |

`YS:Z:` holds the original base sequence of the read **before** the C→T or
G→A conversion. The value is the ASCII string of bases as read from the FASTQ,
in read order (not reverse-complemented).

This tag serves two purposes:

1. **SEQ restoration.** `meth_mem_aln_to_bam` copies the `YS:Z:` payload back
   into the BAM SEQ field (with reverse-complement when `is_rev` is set) so
   that methylation callers see real cytosines. Without this restoration the SEQ
   field would show only `T`s where `C`s existed in the original read.

2. **Downstream inspection.** Tools that need to examine the unconverted sequence
   independently of the BAM SEQ field can read `YS:Z:` directly.

> **Note — Format inside the comment buffer**
>
> Internally, the ingest code stores the value as
> `YS:Z:<bases>\tYC:Z:<dir>` starting at offset 0 of `bseq1_t.comment`.
> `meth_mem_aln_to_bam` locates the payload at `comment + 5` (past the
> `YS:Z:` prefix). The two tags are always co-emitted in this order.

### `YC:Z:` — conversion direction

| Property | Value |
|----------|-------|
| Type | `Z` (NUL-terminated string) |
| Values | `CT` (R1, C→T) or `GA` (R2, G→A) |
| Set by | FASTQ ingest (`src/fastmap.cpp` meth_mode block) |
| Emitted on | All records |

`YC:Z:` records which conversion was applied to the read:

- `CT` — C→T conversion applied; this is an R1 read (or a single-end read).
- `GA` — G→A conversion applied; this is an R2 read.

bwameth.py uses `YC:Z:` for the same purpose and with the same values. Tools
such as MethylDackel use `YC:Z:` to determine which cytosines to call as
methylated. `YC:Z:CT` records are candidates for CpG methylation on the top
strand; `YC:Z:GA` records are candidates on the bottom strand.

### `YD:Z:` — strand hypothesis

| Property | Value |
|----------|-------|
| Type | `Z` (NUL-terminated string) |
| Values | `f` (forward / top strand) or `r` (reverse / bottom strand) |
| Set by | `meth_mem_aln_to_bam` (`src/meth_bam.cpp`) |
| Emitted on | Mapped records only (not unmapped) |

`YD:Z:` records which strand of the doubled reference the read aligned to. The
value is derived from the `f`/`r` prefix on the internal contig name via the
`meth_chrom_map_t.direction` array. Unmapped reads do not receive `YD:Z:`.

- `f` — the read aligned to an `f`-prefixed contig (the C→T projection of the
  top strand).
- `r` — the read aligned to an `r`-prefixed contig (the G→A projection of the
  bottom strand).

This tag is used by `--set-as-failed` (see [Flags](flags.md)) and is also
consumed by downstream methylation callers to confirm which strand each
alignment supports.

## Tag emission summary

| Tag | Records | Source |
|-----|---------|--------|
| `YS:Z:` | All | FASTQ ingest (comment buffer) |
| `YC:Z:` | All | FASTQ ingest (comment buffer) |
| `YD:Z:` | Mapped only | `meth_mem_aln_to_bam` from chrom map |

> **Tip — Checking tags with samtools**
>
> To inspect these tags on a BAM file:
>
>     samtools view out.bam | cut -f12- | grep -oP 'Y[SCD]:Z:[^\t]+'
>
> Or use `samtools view -H` to confirm the `@PG ID:bwa-mem3-meth` entry is
> present and the `@SQ` lines are consolidated (no `f`/`r` prefixes).

---

**See also:**
[Overview](overview.md) ·
[Conversion details](conversion.md) ·
[Chimera QC and header rewriting](post-processing.md) ·
[Flags: --set-as-failed, --do-not-penalize-chimeras](flags.md) ·
[User Guide → Output: SAM/BAM, headers, tags](../user-guide/output.md)
