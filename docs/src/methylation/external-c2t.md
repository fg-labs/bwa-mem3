# Migrating from external bwameth.py c2t

Earlier (D1) `bwa-mem3 --meth` releases mirrored bwameth.py: they aligned
pre-converted reads against a single `ref.fa.bwameth.c2t` doubled reference, so a
bwameth-style external c2t workflow could be wired up directly. **The current
(D3) design does not support that pattern**, because it seeds in 3-letter space
but scores against the *original* 4-letter reference. This page explains what
changed and how to migrate.

> **Important — external c2t interop is no longer supported**
>
> `bwa-mem3 mem --meth` no longer aligns against a `.bwameth.c2t` reference, and
> it cannot consume pre-converted reads. It must do the C→T / G→A projection
> itself (for seeding) so it can keep the original bases for scoring against the
> original reference. Pass **raw FASTQ** and the **original** `ref.fa` prefix.

## What changed

| | D1 (old) | D3 (current) |
|---|----------|--------------|
| Reference used for alignment | `ref.fa.bwameth.c2t` (converted, doubled) | `ref.fa` original 4-letter index + `ref.fa.meth.*` seed index |
| Reads | pre-converted, or inline-converted | raw FASTQ; projected internally for seeding only |
| External `bwameth.py c2t` reads piped in | supported | **not supported** |
| Passing a `.bwameth.c2t` reference path | used as-is | **errors**: "found a legacy '.bwameth.c2t' index … Re-run: bwa-mem3 index --meth" |

Because scoring now runs on the original bases, feeding pre-converted reads would
make the converted bases *look* like the truth and corrupt scoring, `NM`/`MD`,
and the `XM` methylation string. That is why the external-c2t path was removed
rather than adapted.

## How to migrate

1. **Rebuild the index.** A legacy `ref.fa.bwameth.c2t` index is not usable.
   Build the dual index from the original FASTA:

   ```sh
   bwa-mem3 index --meth ref.fa     # writes ref.fa.* and ref.fa.meth.*
   ```

2. **Pass raw FASTQ and the original prefix.** Drop any `bwameth.py c2t`
   preprocessing step and the `-p /dev/stdin` plumbing:

   ```sh
   bwa-mem3 mem --meth -t 16 ref.fa R1.fq.gz R2.fq.gz \
     | samtools sort -o out.bam
   ```

   bwa-mem3 finds `ref.fa.meth.*` automatically and does the seed-time projection
   internally. (You can pass the `ref.fa.meth` seed-index path directly if you
   prefer, but the original-reference handles must sit alongside it.)

If you specifically need bwameth.py's collapsed-space alignment or its own c2t
tooling, continue to use bwameth.py itself — see
[bwameth.py drop-in mapping](bwameth-mapping.md) for how `--meth-scoring collapsed`
reproduces its placement instead.

---

**See also:**
[Overview](overview.md) ·
[Conversion details](conversion.md) ·
[SAM tags: XR, XG, XM](tags.md) ·
[bwameth.py drop-in mapping](bwameth-mapping.md) ·
[Related Projects: bwameth.py](../related-projects/bwameth.md)
