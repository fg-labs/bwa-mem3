# Chimera QC and Header Rewriting

After the alignment kernel produces `mem_aln_t` records, `bwa-mem3 --meth`
applies a set of post-processing steps before writing BAM output. These steps
are implemented in `src/meth_bam.cpp` and run in the same process, in the same
pass over the aligned records.

## `@SQ` header consolidation

The doubled reference (`ref.fa.bwameth.c2t`) contains two contigs for each
chromosome:

- `fchr1`, `fchr2`, … — C→T projections of each chromosome.
- `rchr1`, `rchr2`, … — G→A projections of each chromosome.

If the raw alignment header were written directly, every downstream tool would
see twice as many sequences as there are real chromosomes, with unfamiliar
`f`/`r`-prefixed names. `meth_bam_writer_open` instead builds a consolidated
header using the `meth_chrom_map_t`:

1. `meth_chrom_map_build_from_bns` iterates over `bns->anns` and strips the
   leading `f`/`r` from each contig name.
2. The first contig with a given stripped name registers that name in the output
   list; subsequent contigs with the same stripped name map to the same output
   index.
3. The BAM `@SQ` lines are written from the consolidated list — one `SN:` per
   real chromosome.

RNAME, RNEXT, and SA/XA tag contig references in every record are rewritten
through `cmap->out_tid` and `cmap->output_names` so they reference the
consolidated names. The mapping from internal (doubled-ref) contig index to
output contig index is `cmap->out_tid[p.rid]`.

> **Note — TLEN computation uses consolidated TIDs**
>
> Template length (TLEN) is computed using the consolidated output TIDs, not the
> internal `p.rid` values. Two mates that rescue onto `fchr1` and `rchr1`
> respectively both map to output `chr1`, so TLEN is reported as a non-zero
> distance rather than zero (which would happen if the mismatched internal TIDs
> were used).

## Chimera QC heuristic

bwameth.py applies a heuristic to flag reads that look like chimeric fragments:
if the longest contiguous alignment run (sum of M/=/X CIGAR operations) covers
less than 44% of the read length, the read is considered a potential chimera.

`bwa-mem3 --meth` applies the same heuristic inside `meth_mem_aln_to_bam`:

```text
if (100 * longest_M_run < 44 * l_seq):
    flag  |=  0x200   # set QC fail
    flag  &= ~0x2     # clear proper pair
    mapq   =  min(mapq, 1)
```

The threshold constant is `MIN_LONGEST_M_PCT = 44` (defined at the top of
`src/meth_bam.cpp`). The longest run is computed by `cigar_longest_m_mem`
from `src/cigar_util.cpp`, which counts M, `=`, and X operations.

The chimera heuristic is only applied to mapped records (`!(flag & 0x4) &&
direction != 0`). Unmapped records are not touched.

To disable this heuristic, pass `--do-not-penalize-chimeras`. See
[Flags](flags.md) for details.

## `--set-as-failed` strand filtering

Before the chimera check, `meth_mem_aln_to_bam` checks whether
`opt->meth_set_as_failed` is set and matches the record's strand direction:

```text
if (meth_set_as_failed != 0 && meth_set_as_failed == direction):
    flag |= 0x200
```

This unconditionally marks all alignments to the specified strand (`f` or `r`)
as QC-failed before chimera logic runs. The chimera check then applies on top
of the already-set fail flag.

## Pair-level QC-fail propagation

Once per read group (all records sharing the same query name), after individual
records have been processed:

```text
meth_bam_group_propagate_qcfail(group, n)
```

This function scans all records in the group. If any record has `0x200` set, it
propagates that flag to every other record in the group and clears `0x2`
(proper pair) on all of them. This ensures that a chimeric or strand-filtered
primary alignment also marks its supplementary alignments and the mate as
QC-failed, preventing inconsistent flag states in the output BAM.

## `@PG ID:bwa-mem3-meth` insertion

`meth_bam_writer_open` appends a `@PG` line to the header after the original
bwa-mem3 `@PG` entry:

```text
@PG  ID:bwa-mem3-meth  PN:bwa-mem3-meth  VN:<version>-meth  CL:<command line>
```

The `<command line>` field is the full `bwa-mem3 mem --meth ...` invocation with
embedded tab characters replaced by spaces (htslib does not permit literal tabs
in `@PG CL:` fields). This records the exact parameters used for provenance
and reproducibility.

> **Tip — Verifying the header**
>
> After alignment, confirm consolidation and provenance with:
>
> ```sh
> samtools view -H out.bam | grep -E '^@SQ|^@PG'
> ```
>
> You should see one @SQ line per chromosome (no f/r prefixes) and both
> `@PG ID:bwa-mem3` and `@PG ID:bwa-mem3-meth` entries.

---

**See also:**
[Overview](overview.md) ·
[SAM tags: YS, YC, YD](tags.md) ·
[Flags: --set-as-failed, --do-not-penalize-chimeras](flags.md) ·
[Conversion details](conversion.md) ·
[User Guide → Output: SAM/BAM, headers, tags](../user-guide/output.md)
