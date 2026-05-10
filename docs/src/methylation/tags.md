# SAM Tags: XR, XG, XM (Bismark-compatible)

`bwa-mem3 mem --meth` emits three Bismark-compatible auxiliary tags on
each output record: `XR:Z`, `XG:Z`, and `XM:Z`. These tags are read by
`bismark_methylation_extractor`, `deduplicate_bismark`, methylKit
`processBismarkAln`, methtuple, DMRfinder, epialleleR, MethylDackel, and
biscuit's per-read methylation tools.

## Tag reference

### `XR:Z` — read conversion direction

| Property | Value |
|----------|-------|
| Type | `Z` (NUL-terminated string) |
| Values | `CT` (R1 / SE) or `GA` (R2) |
| Set by | `meth_mem_aln_to_bam` from FASTQ-ingest carrier (`s->comment`'s YC payload) |
| Emitted on | All records (mapped and unmapped) |

`XR:Z` records which conversion was applied to the read at FASTQ ingest:

- `CT` — C→T conversion applied; this is an R1 read or single-end read.
- `GA` — G→A conversion applied; this is an R2 read.

### `XG:Z` — genome strand of the alignment

| Property | Value |
|----------|-------|
| Type | `Z` (NUL-terminated string) |
| Values | `CT` (aligned to original top, `f-`-prefixed contig) or `GA` (aligned to original bottom, `r-`-prefixed contig) |
| Set by | `meth_mem_aln_to_bam` from `meth_chrom_map_t.direction` |
| Emitted on | Mapped records only |

`XG:Z` indicates which doubled-reference strand the read aligned to:

- `CT` — read aligned to the C→T-projected forward strand (OT).
- `GA` — read aligned to the G→A-projected forward strand (OB).

For properly paired directional reads, R1 and R2 of a fragment naturally
share `XG:Z`. Discordant pairs (already flagged with `0x200` by the
chimera-QC heuristic) may see `XG:Z` diverge between mates.

### `XM:Z` — methylation call string

| Property | Value |
|----------|-------|
| Type | `Z` (NUL-terminated string) |
| Length | Equal to `SEQ` length |
| Set by | `meth_build_xm` (`src/meth_xm.cpp`) walking SEQ-orientation read against un-converted ref |
| Emitted on | Mapped records only |

Per-base methylation call. Each character corresponds to one SEQ base:

| char | meaning |
|------|---------|
| `z` / `Z` | unmethylated / methylated C in CpG context |
| `x` / `X` | unmethylated / methylated C in CHG context |
| `h` / `H` | unmethylated / methylated C in CHH context |
| `u` / `U` | unmethylated / methylated C in unknown context (N within 1 or 2 bp downstream of the C, on the read's source strand) |
| `.`       | non-C reference at this position, sequencing mismatch (read base ≠ C/T at a ref C), insertion, soft clip, or N at the C position itself |

The string is in SEQ orientation (matches the BAM SEQ field): for reads
with the `0x10` flag set, both SEQ and `XM:Z` are reverse-complemented
relative to FASTQ-original orientation.

## Computation

Under `--meth`, the doubled c2t reference (`<prefix>.bwameth.c2t.*`) is
folded once at startup into an in-memory un-converted pac (the
`meth_orig_ref` module — `src/meth_orig_ref.cpp`). The fold uses
`(f, r) → original` recovery on every position via a 5-row table:

| f[P] | r[P] | original[P] |
|------|------|-------------|
| T    | T    | T           |
| T    | C    | C           |
| G    | A    | G           |
| A    | A    | A           |
| N    | N    | N (via `bns->ambs`) |

Per mapped record, `meth_build_xm` slices the un-converted forward-strand
window at the read's footprint plus 2 bp of context on either side, then
walks the BAM CIGAR jointly over the restored SEQ and the ref window.
The classifier matches Bismark's `methylation_call`:

```
match position with ref[t] == 'C' (top strand) or 'G' (bottom strand):
    determine context from ref[t±1], ref[t±2]
        N in either context base    -> u/U (unknown context)
        ref[t±1] == G/C (per strand) -> z/Z (CpG)
        ref[t±2] == G/C (per strand) -> x/X (CHG)
        otherwise                    -> h/H (CHH)
    determine methylation:
        read base == C/G (per strand) -> uppercase (methylated)
        read base == T/A (per strand) -> lowercase (unmethylated)
        otherwise                      -> '.'
insertion / soft clip                  -> '.' per consumed read base
deletion / N op                         -> no XM emit
hard clip / pad                         -> no XM emit
```

The `top vs bottom strand` choice is driven by `XG:Z` (= cmap
direction), **not** by the SAM 0x10 (RC) flag. CTOT reads (R2 mapped
forward to a top-strand contig with 0x10 set) and OB reads (R1 mapped
RC to a bottom-strand contig) are both handled by reading the rule
table from the strand encoded in `XG`. The walk runs in SEQ orientation
throughout — no RC of the ref slice or the read.

For bottom-strand methylation, the C of interest at forward position P
is encoded as a G on the forward strand (complement of bottom-strand C).
The downstream context on the bottom strand corresponds to *upstream*
positions on the forward strand; the classifier indexes `ref[t-1]` and
`ref[t-2]` instead of `ref[t+1]` and `ref[t+2]`, and looks for a `C`
(forward) instead of a `G` to flag CpG.

## Inspecting tags with samtools

```
samtools view out.bam | head -1 | tr '\t' '\n' | grep -E '^X[RGM]:'
```

Expected output looks like:

```
XR:Z:CT
XG:Z:CT
XM:Z:..z..h..Z..x..h.....Z..
```

## See also

- [Overview](overview.md)
- [Conversion details](conversion.md)
- [Chimera QC and header rewriting](post-processing.md)
- [Flags: --set-as-failed, --chimera-qc](flags.md)
- [User Guide → Output: SAM/BAM, headers, tags](../user-guide/output.md)
