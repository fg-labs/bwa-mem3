# Quick start: methylation alignment

bwa-mem3 supports bisulfite-converted (WGBS/RRBS/EM-seq) read alignment through a single `--meth`
flag on both `index` and `mem`. No Python interpreter, no piped preprocessor, and no separate
postprocessing step are required.

> **Note — variant-aware by default, bwameth-compatible on request**
>
> By default (`--meth-scoring genomic`) bwa-mem3 frees only the conversion direction and keeps
> other C/T and G/A differences as mismatches, so variants outside the conversion direction count
> in `NM` (a genuine C→T at a reference `C` is indistinguishable from a conversion and is not
> counted, see
> [which real variants stay visible](../methylation/overview.md#which-real-variants-stay-visible)).
> It emits the standard Bismark tags methylation callers expect (MethylDackel, Bismark, PileOMeth,
> etc.). Add `--meth-scoring collapsed` to closely track bwameth.py's read **placement** instead:
> a placement drop-in, *not* a byte-for-byte reproduction (~1% of records differ in
> `POS`/`CIGAR`/`MAPQ`; see [bwameth.py drop-in mapping](../methylation/bwameth-mapping.md)).
> `--meth-scoring neutral`, the `--meth=taps` / `--meth=5base` default, is variant-aware too but
> scores the conversion `0` instead of a match.
>
> **Trim inline UMIs first.** If your reads carry an inline UMI or spacer (unconverted synthetic
> bases at the start of each read), extract or trim it before `--meth`. Otherwise those bases are
> sometimes soft-clipped, but can remain aligned and receive methylation calls under either mode
> (`collapsed` aligns them more often, since it treats C/T and G/A as interchangeable).

## Index the reference for methylation

Build the methylation index once:

```bash
bwa-mem3 index --meth ref.fa
```

This builds the normal 4-letter index at the bare prefix **plus** a converted seed index:

| File(s) | Description |
|---------|-------------|
| `ref.fa.amb`, `ref.fa.ann`, `ref.fa.bwt.2bit.64`, `ref.fa.pac` | Normal index over the original reference (used for scoring/extension; bases are pac-fetched from `.pac`, so no `.0123` is built). |
| `ref.fa.meth.fa` (only with `--keep-meth-fasta`) | Per-strand C→T / G→A converted FASTA (`f`/`r` doubled contigs). |
| `ref.fa.meth.*` | FM-index over the converted FASTA (used only for **seeding**). |

Both indexes live next to `ref.fa`; `bwa-mem3 index ref.fa` (no `--meth`) builds only the first set.

To write elsewhere, add `-p /idx/ref` and align against `/idx/ref`. For methylation-only use, add `--meth-only`: it omits the original `.bwt.2bit.64`, leaving seven index files and any source sidecars. Ordinary alignment then requires rebuilding without `--meth-only`. Existing outputs require `--force`; the converted FASTA is temporary unless explicitly retained. See [indexing](../user-guide/indexing.md#minimal-methylation-output-and-safe-rebuilds) for safety and scratch details.

## Align bisulfite-converted reads

```bash
# Default: genomic (variant-aware)
bwa-mem3 mem --meth -t 16 ref.fa R1.fq.gz R2.fq.gz \
  | samtools sort -o out.bam
samtools index out.bam

# Opt into bwameth-compatible placement
bwa-mem3 mem --meth --meth-scoring collapsed -t 16 ref.fa R1.fq.gz R2.fq.gz \
  | samtools sort -o out.bam

# Single-end (RRBS etc.): pass one FASTQ. The C→T (R1) projection is used for all reads.
bwa-mem3 mem --meth -t 16 ref.fa R1.fq.gz \
  | samtools sort -o out.bam
```

Pass the original-reference output prefix (the FASTA path only when `-p` was omitted). bwa-mem3 finds the `<prefix>.meth.*` seed index
automatically when `--meth` is active. The paired-end-only flags (`-U 100`, `-p`) do not apply
to single-end input.

## What `--meth` does

In one process, `--meth` replaces what bwameth.py does with external tools:

- **Seeds in 3-letter space, scores in 4-letter space** — projects each read
  (R1 `C→T`, R2 `G→A`) to seed against `ref.fa.meth.*`, then extends and scores
  against the *original* reference, restoring the original bases in `SEQ`.
- **Applies bwameth-aligned defaults** — `-L 10 -U 100 -T 40 -M -C` plus
  mode-dependent `-B` (`4` for `genomic`, the default; `2` for `collapsed`).
- **Post-processes records inline** — writes `@SQ` from the original reference
  (the `f`/`r` seed contigs never reach the output), emits Bismark `XR`/`XG`/`XM`
  tags, and runs optional `--chimera-qc`. Output is SAM text by default and
  uncompressed BAM under `--bam`, either one ready for `samtools sort`.

See the [Methylation Reference — Overview](../methylation/overview.md) for the
mechanism in detail, [Flags](../methylation/flags.md) for the scoring modes and
QC flags, and [SAM tags](../methylation/tags.md) for the tag definitions.

---

**See also:**
[Methylation Reference — Overview](../methylation/overview.md) ·
[Methylation Reference — Flags](../methylation/flags.md) ·
[Methylation Reference — SAM tags](../methylation/tags.md) ·
[Best Practices — Methylation defaults](../best-practices/methylation.md) ·
[CLI Reference — mem](../cli/mem.md)
