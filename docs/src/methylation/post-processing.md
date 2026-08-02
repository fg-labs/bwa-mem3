# Chimera QC and header construction

After the alignment kernel produces `mem_aln_t` records, `bwa-mem3 --meth`
applies a set of post-processing steps before writing its output. These steps
are implemented in `src/meth_bam.cpp` and run in the same process, in the same
pass over the aligned records. They are independent of the output container —
`--meth` writes SAM text by default and BAM under `--bam`, and the records are
identical either way.

## `@SQ` headers come from the original reference

The `.meth` **seed** reference (`ref.fa.meth.fa`) contains two contigs for each
chromosome:

- `fchr1`, `fchr2`, … — C→T projections of each chromosome.
- `rchr1`, `rchr2`, … — G→A projections of each chromosome.

Those `f`/`r` contigs exist only for seeding, and they never reach the output.
Since D3 (#174), `--meth` emits native original-alphabet output: seeds are
remapped to original coordinates plus an OT/OB hypothesis at the end of the
seeding phase (`orig_tid = seed_rid / 2`, `hypothesis = seed_rid & 1`), so every
mapped `mem_aln_t` that reaches the writer already carries an **original** rid.
Unmapped records keep the `rid = -1` sentinel `mem_reg2aln` gives them.

`meth_bam_writer_open` therefore builds `@SQ` straight from the original
(un-converted) reference's `bns->anns` — one `SN:` per reference sequence, in the
original reference's order — and enriches each line with the `M5`/`UR`/`AS`/`SP`
identity tags from that reference's `.hdr`/`.dict` sidecar when one is present.

Per record, `meth_mem_aln_to_bam` sets RNAME and RNEXT from `p.rid` and
`mp->rid` directly. `SA:Z` names are looked up in `bns->anns` and `XA:Z` is
emitted verbatim. There is no name rewriting and no contig-index translation
anywhere on the path.

> **Historical note**
>
> Before D3, output carried the doubled reference's coordinates and a
> `meth_chrom_map_t` collapsed the `f`/`r` names at emit time, translating every
> RNAME/RNEXT/SA/XA reference through an `out_tid` table. That map and its
> consolidation step were removed when the coordinate cutover landed; the
> original-alphabet rids make it unnecessary. Documentation describing
> `meth_chrom_map_*`, `out_tid`, or `output_names` refers to the pre-D3 design.

TLEN is computed from `p.pos`/`mp->pos` when both mates carry a CIGAR and land
on the same contig (`tid == mtid`). Because rids are already original, that test
is a plain same-chromosome check — the pre-D3 subtlety, where two mates rescued
onto `fchr1` and `rchr1` had different internal rids for the same real
chromosome, cannot arise.

## Chimera QC heuristic (opt-in)

bwameth.py applies a heuristic to flag reads that look like chimeric
fragments: if the longest contiguous alignment run (sum of M/=/X CIGAR
operations) covers less than 44 % of the read length, the read is
considered a potential chimera. **Bismark does not apply this kind of
heuristic.**

`bwa-mem3 --meth` makes this opt-in via `--chimera-qc` (default off, so
the runtime posture matches Bismark). When enabled, the check inside
`meth_mem_aln_to_bam` does:

```text
if (100 * longest_M_run < 44 * l_seq):
    flag  |=  0x200   # set QC fail
    flag  &= ~0x2     # clear proper pair
    mapq   =  min(mapq, 1)
```

The threshold constant is `MIN_LONGEST_M_PCT = 44` (defined at the top
of `src/meth_bam.cpp`). The longest run is computed by
`cigar_longest_m_mem` from `src/cigar_util.cpp`, which counts M, `=`,
and X operations.

The chimera heuristic is only applied to mapped records (`!(flag & 0x4)
&& direction != 0`). Unmapped records are not touched.

See [Flags](flags.md) for when to use `--chimera-qc` (PBAT / scBS-Seq;
bwameth.py-equivalence runs).

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
primary alignment also marks its split hits and the mate as QC-failed,
preventing inconsistent flag states in the output BAM. (Under `--meth` those
split hits carry `0x100`, not `0x800` — see
[`-M` and split alignments](flags.md#-m-and-split-alignments).)

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
> After alignment, confirm the header and provenance with (samtools autodetects
> both containers, so this works on `--meth` and `--meth --bam` output alike):
>
> ```bash
> samtools view -H out.bam | grep -E '^@SQ|^@PG'
> ```
>
> You should see one @SQ line per reference sequence (no f/r prefixes) and both
> `@PG ID:bwa-mem3` and `@PG ID:bwa-mem3-meth` entries.

---

**See also:**
[Overview](overview.md) ·
[bwameth.py drop-in mapping](bwameth-mapping.md) ·
[SAM tags: XR, XG, XM](tags.md) ·
[Flags: --set-as-failed, --chimera-qc](flags.md) ·
[Conversion details](conversion.md) ·
[User Guide → Output: SAM/BAM, headers, tags](../user-guide/output.md)
