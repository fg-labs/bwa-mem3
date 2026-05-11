# Flags: --set-as-failed, --chimera-qc

`bwa-mem3 --meth` adds two flags that control QC behavior during BAM
post-processing. Both flags affect the chimera QC and strand-filtering logic
inside `meth_mem_aln_to_bam` (`src/meth_bam.cpp`).

## `--set-as-failed {f|r}`

Marks every alignment to the specified strand as QC-failed (`0x200`) regardless
of alignment quality or CIGAR structure.

**Accepted values:**

- `f` — flag all alignments to `f`-prefixed contigs (C→T top-strand projection).
- `r` — flag all alignments to `r`-prefixed contigs (G→A bottom-strand
  projection).

**Effect on records:**

When `--set-as-failed f` (or `r`) is set and a mapped record's strand matches
the specified value, the record's SAM flag has `0x200` set. If `--chimera-qc`
is also active, the chimera heuristic runs on top, possibly clearing `0x2` and
capping MAPQ. QC-fail propagation then spreads the flag to all records in the
read group.

**When to use it:**

Some experimental designs produce reads that are expected to align exclusively
to one strand. Flagging the other strand as QC-failed before downstream
analysis prevents spurious methylation calls from mis-strand alignments. It is
also useful for diagnosing library preparation issues: run once with
`--set-as-failed r` and once without to compare yield on each strand.

> **Warning — All records on the strand are flagged**
>
> `--set-as-failed` is a blunt instrument. It marks every alignment to the chosen
> strand, including correctly aligned reads that simply happened to land on the
> complementary strand due to library structure. Use this flag only when your
> library is expected to be strand-specific.

## `--chimera-qc`

Enables the bwameth.py-style longest-M chimera heuristic. **Off by default**;
this is the Bismark-equivalent posture, since Bismark itself does not apply
this kind of QC heuristic.

When `--chimera-qc` is set, any mapped record whose longest M/=/X CIGAR run
covers less than 44 % of the read length receives:

- `0x200` (QC fail) set.
- `0x2` (proper pair) cleared.
- MAPQ capped at 1.

QC-fail propagation across the read group also applies.

**When to use it:**

The 44 % threshold was calibrated by bwameth.py for standard mammalian whole-
genome bisulfite-sequencing (WGBS) libraries with typical read lengths and is
helpful on PBAT / scBS-Seq libraries where intra-fragment chimeras are common.
For Bismark-equivalent output (and most directional EM-seq / WGBS workflows),
leave it off.

It is also useful when benchmarking: comparing `bwa-mem3 --meth` output
against bwameth.py output is cleaner with `--chimera-qc` enabled, since
bwameth.py's chimera logic always runs.

> **Note — Pair-level propagation still applies**
>
> `--chimera-qc` controls only whether the heuristic itself runs.
> `--set-as-failed` is independent: when active, those flags are still set,
> and `meth_bam_group_propagate_qcfail` propagates any `0x200` flags across
> the read group regardless of `--chimera-qc`.

## Flag interaction summary

| Condition | `0x200` set? | `0x2` cleared? | MAPQ capped? |
|-----------|-------------|----------------|--------------|
| Normal aligned record (default, no flags) | No | No | No |
| `--chimera-qc` triggers (longest M/=/X < 44%) | Yes | Yes | Yes (≤1) |
| `--set-as-failed` strand matches | Yes | No | No |
| Both `--chimera-qc` + `--set-as-failed` active | Yes | Yes | Yes (≤1) |

## `-V` reference annotation `XR:Z` is suppressed under `--meth`

`bwa-mem3 mem -V` normally emits the contig annotation as an `XR:Z`
auxiliary field. Under `--meth`, `XR:Z` carries the Bismark
read-conversion direction (`CT`/`GA`) instead. The reference-annotation
`XR:Z` is silently suppressed when `--meth` is active so the two uses
don't collide. There is no flag to override this — `-V` is a no-op for
`XR:Z` under `--meth`. See [tags.md](tags.md).

---

**See also:**
[Overview](overview.md) ·
[Chimera QC and header rewriting](post-processing.md) ·
[SAM tags: XR, XG, XM](tags.md) ·
[Best Practices → Methylation defaults](../best-practices/methylation.md) ·
[CLI Reference → mem](../cli/mem.md)
