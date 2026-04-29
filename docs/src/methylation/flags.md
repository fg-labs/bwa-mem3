# Flags: --set-as-failed, --do-not-penalize-chimeras

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

When `--set-as-failed f` (or `r`) is set and a mapped record's `YD:Z:` strand
matches the specified value, the record's SAM flag has `0x200` set. The chimera
heuristic then runs on top of the already-set flag, but since `0x200` is already
present, it can only enforce additional constraints (clearing `0x2`, capping
MAPQ). QC-fail propagation then spreads the flag to all records in the read
group.

**When to use it:**

Some experimental designs produce reads that are expected to align exclusively to
one strand. Flagging the other strand as QC-failed before downstream analysis
prevents spurious methylation calls from mis-strand alignments. It is also
useful for diagnosing library preparation issues: run once with `--set-as-failed r`
and once without to compare yield on each strand.

> **Warning — All records on the strand are flagged**
>
> `--set-as-failed` is a blunt instrument. It marks every alignment to the chosen
> strand, including correctly aligned reads that simply happened to land on the
> complementary strand due to library structure. Use this flag only when your
> library is expected to be strand-specific.

## `--do-not-penalize-chimeras`

Disables the chimera QC heuristic entirely.

Without this flag, any mapped record whose longest M/=/X CIGAR run covers less
than 44% of the read length receives:

- `0x200` (QC fail) set.
- `0x2` (proper pair) cleared.
- MAPQ capped at 1.

With `--do-not-penalize-chimeras`, none of these penalties are applied. Records
are written with the MAPQ and flags as determined by the alignment kernel. The
chimera check in `meth_mem_aln_to_bam` is skipped entirely.

**When to use it:**

The 44% threshold was calibrated for standard mammalian whole-genome bisulfite
sequencing (WGBS) libraries with typical read lengths. For short reads (< 50 bp),
reads with large insertions, or amplicon designs where short alignments are
expected, the heuristic can produce excessive false-positive flagging. In those
cases, disable it and apply a custom chimera filter downstream.

It is also useful when benchmarking: comparing `bwa-mem3 --meth` output against
bwameth.py output on a specific dataset is cleaner when chimera filtering is
disabled, since bwameth.py's chimera logic may differ in edge cases.

> **Note — Pair-level propagation still applies**
>
> `--do-not-penalize-chimeras` only suppresses the chimera heuristic.
> If `--set-as-failed` is also active, those flags are still set, and
> `meth_bam_group_propagate_qcfail` still propagates any `0x200` flags
> across the read group.

## Flag interaction summary

| Condition | `0x200` set? | `0x2` cleared? | MAPQ capped? |
|-----------|-------------|----------------|--------------|
| Normal aligned record | No | No | No |
| Chimera heuristic triggers (default) | Yes | Yes | Yes (≤1) |
| `--set-as-failed` strand matches | Yes | No | No |
| Both chimera + `--set-as-failed` active | Yes | Yes | Yes (≤1) |
| `--do-not-penalize-chimeras` only | No | No | No |

---

**See also:**
[Overview](overview.md) ·
[Chimera QC and header rewriting](post-processing.md) ·
[SAM tags: YS, YC, YD](tags.md) ·
[Best Practices → Methylation defaults](../best-practices/methylation.md) ·
[CLI Reference → mem](../cli/mem.md)
