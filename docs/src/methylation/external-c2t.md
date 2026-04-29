# Interop with External bwameth.py c2t

Some workflows use bwameth.py's `c2t` subcommand to convert reads before
passing them to an aligner. `bwa-mem3 --meth` supports this pattern by
detecting whether the caller has already provided a pre-converted FASTQ and
whether the reference path already points to the doubled-reference FASTA.

## Auto-detect logic for the reference path

When `--meth` is active, `bwa-mem3 mem` ordinarily appends `.bwameth.c2t` to
the reference path so the user can pass the original FASTA prefix:

```sh
bwa-mem3 mem --meth -t 16 ref.fa R1.fq.gz R2.fq.gz
# internally uses ref.fa.bwameth.c2t as the reference
```

If the reference path already ends with `.bwameth.c2t`, the auto-append is
skipped:

```sh
bwa-mem3 mem --meth -t 16 ref.fa.bwameth.c2t R1.fq.gz R2.fq.gz
# no suffix appended; ref.fa.bwameth.c2t is used as-is
```

This detection is a simple suffix check on the path string. It allows callers
that manage the doubled-reference path explicitly to pass it without triggering
double-append.

## Using bwameth.py c2t as the read preprocessor

If your pipeline already runs `bwameth.py c2t` to convert reads (for example,
because it needs to reuse converted reads across multiple aligners), you can
pipe the output directly to `bwa-mem3 mem --meth`:

```sh
bwameth.py c2t R1.fq.gz R2.fq.gz \
  | bwa-mem3 mem --meth -p -t 16 ref.fa.bwameth.c2t /dev/stdin \
  | samtools sort -o out.bam
```

Key points for this pattern:

- Pass the `.bwameth.c2t` reference path explicitly so the auto-append is
  suppressed.
- Use `-p` to tell `bwa-mem3 mem` that the input contains interleaved
  paired-end reads (bwameth.py c2t emits interleaved output to stdout).
- Use `/dev/stdin` as the reads argument to read from the pipe.
- The `bwa-mem3 --meth` inline c2t conversion is **not** applied when the
  reads arrive pre-converted; however, `YS:Z:` and `YC:Z:` tags are only
  written by the inline conversion path. If you need those tags in the output,
  you must either use the integrated mode (no external c2t step) or ensure
  your external preprocessor emits compatible comments in the FASTQ.

> **Warning — YS:Z: and YC:Z: require integrated mode**
>
> When reads are pre-converted by an external tool and piped in, the inline
> c2t step in `src/fastmap.cpp` is bypassed. `YS:Z:` and `YC:Z:` tags will
> not be present in the output BAM unless the external converter writes them
> as FASTQ comment fields in the format `YS:Z:<seq>\tYC:Z:<dir>` and those
> comments are passed through. MethylDackel and similar callers use `YS:Z:`
> to restore original bases for methylation calling; if the tag is absent,
> they fall back to reading SEQ directly, which may affect accuracy.

## Header rewriting and BAM post-processing with external c2t

Whether reads are converted inline or externally, all BAM post-processing steps
apply identically when `--meth` is active:

- `@SQ` header consolidation (f/r contigs → one entry per chromosome).
- `YD:Z:` tag emission from the contig name prefix.
- Chimera QC heuristic (unless `--do-not-penalize-chimeras` is set).
- Pair-level QC-fail propagation.
- `@PG ID:bwa-mem3-meth` insertion.

The post-processing pipeline depends only on the reference contig names (to
determine `YD:Z:`) and the alignment flags — not on whether reads were
converted inline or externally.

## Summary of path variants

| Reference arg | Read source | Auto-append? | Inline c2t? | `YS:Z:` emitted? |
|---------------|-------------|-------------|------------|-----------------|
| `ref.fa` | Raw FASTQ | Yes (→ `ref.fa.bwameth.c2t`) | Yes | Yes |
| `ref.fa.bwameth.c2t` | Raw FASTQ | No | Yes | Yes |
| `ref.fa.bwameth.c2t` | Pre-converted (pipe) | No | No | No (unless pre-emitted) |

---

**See also:**
[Overview](overview.md) ·
[Conversion details](conversion.md) ·
[SAM tags: YS, YC, YD](tags.md) ·
[bwameth.py drop-in mapping](bwameth-mapping.md) ·
[Related Projects: bwameth.py](../related-projects/bwameth.md)
