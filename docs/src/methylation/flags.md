# Flags: --meth chemistry, --meth-scoring, --meth-tags, --set-as-failed, --chimera-qc

`--meth` itself takes an optional chemistry argument. Alongside it,
`--meth-scoring` selects the scoring model and `--meth-tags` selects which
Bismark tags are emitted; `--set-as-failed` and `--chimera-qc`
control QC behavior during BAM post-processing (both affect the chimera QC and
strand-filtering logic inside `meth_mem_aln_to_bam`, `src/meth_bam.cpp`).

## `--meth[={emseq|taps}]`

Selects the **methylation chemistry**. Defaults to `emseq`, so a bare `--meth`
behaves exactly as it always has.

| Value | Chemistry | Which cytosines convert to T |
|-------|-----------|------------------------------|
| `emseq` (default) | bisulfite / EM-seq | **un**methylated C (aliases: `em-seq`, `bisulfite`) |
| `taps` | TET-assisted pyridine borane | **methylated** 5mC/5hmC |

> **Use `--meth=taps`, not `--meth taps`.** The argument is optional, and
> `getopt_long` only accepts optional arguments in the `=` form; a separate word
> is parsed as a positional (i.e. your first FASTQ).

**What the chemistry changes — and what it does not.** Both chemistries produce
the *same* base change as far as the aligner can see (ref C → read T on the top
strand, ref G → read A on the bottom). Seeding, the converted index, and the
scoring matrices are therefore **shared**, which is why TAPS reads align correctly
even without this flag. The chemistries differ only in *which* cytosines convert,
and that inverts the meaning of the observed base when **calling** methylation:

| Read base at a reference C | `emseq` | `taps` |
|---|---|---|
| `C` (retained) | methylated (`Z`) | **un**methylated (`z`) |
| `T` (converted) | **un**methylated (`z`) | methylated (`Z`) |

So the flag changes the `XM:Z` call polarity, and nothing else. CpG/CHG/CHH
context classification is a property of the reference, not the chemistry, and is
identical either way.

> **Using `--meth` on TAPS data inverts every methylation call.** Measured on
> 4.07 M simulated TAPS reads against full hg38, aligner `XM:Z` correlated
> **r = −0.956** with truth over ~600 k CpG sites; with `--meth=taps` the same
> reads give **r = +0.956** (RMSE 0.777 → 0.117). Percent methylation comes out
> as `1 − truth`, silently. Always pass `=taps` for TAPS libraries, or suppress
> `XM:Z` and call with a TAPS-aware extractor such as Rastair.

**Scoring default.** `--meth=taps` defaults `--meth-scoring` to `neutral`
(bare `--meth` still defaults to `collapsed`). TAPS converts only the methylated
cytosines, so conversions are far sparser than under EM-seq — 3.2 % of cytosines
versus 94.6 % on a matched chr22 simulation — and the collapsed 3-letter alphabet
costs specificity it no longer buys back. `neutral` goes one step past `genomic`:
it scores the conversion cell as `0` rather than a full match, so a sparse TAPS
conversion is *tolerated* without *rewarding* spurious C→T alignments. Measured on
4.07 M simulated TAPS reads against full hg38, across three methylation loads:
`neutral` placed 95.95–96.01 % (essentially the 96.02 % unconverted ceiling),
`genomic` 95.68–95.73 %, `collapsed` 95.37–95.43 % — a robust +0.24–0.28 pp over
`genomic` — with **identical** MAPQ calibration (60+ bin 99.99 % either way) and
`NM` (mean 1.966 vs 1.969), so the gain is not bought with over-confidence or
inflated edit distance. An explicit `--meth-scoring` always wins.

## `--meth-scoring {collapsed|genomic|neutral}`

Selects how the 4-letter scoring matrix treats bisulfite/TAPS-converted bases.
`bwa-mem3 --meth` scores against the original reference, so it can collapse the
conversion (bwameth-style), keep it variant-aware, or merely tolerate it.

**Accepted values:**

- `collapsed` (**default for `--meth`/`--meth=emseq`**) — free **both** conversion
  directions: C↔T *and* G↔A are interchangeable (a two-cell matrix). Reproduces
  bwameth's collapsed-space **placement** and sets the mismatch penalty to `-B 2`.
  Use this when you need bwameth-compatible read placement (the drop-in default).
- `genomic` — free **only** the conversion direction (a one-cell matrix), scored
  as a full match (`+A`), so the mirror cell stays a real mismatch. A genuine C/T
  or G/A variant on that mirror cell (ref `T` × read `C` on OT, ref `A` × read
  `G` on OB) is penalized and stays visible in `NM`/`MD`, making the BAM usable
  for variant calling; a variant in the conversion direction itself shares the
  freed cell with a conversion and stays hidden
  ([why](overview.md#which-real-variants-stay-visible)). Keeps bwa's default `-B 4`.
- `neutral` (**default for `--meth=taps`**) — free only the conversion direction,
  but score it `0` rather than `+A`: tolerated but not rewarded. Best for TAPS,
  whose conversions are sparse, where a full-match reward over-credits spurious
  C→T alignments (see the measurement above). Keeps `-B 4`; the mirror cell stays
  a mismatch, so real variants remain visible in `NM`/`MD` on the same terms as
  `genomic`.

> **Important — `collapsed` is a placement drop-in, not byte-identical to bwameth**
>
> `collapsed` *closely tracks* bwameth's placement but scores against the original
> 4-letter reference, so ~1% of records differ from bwameth.py in `POS`/`CIGAR`/`MAPQ`.
> **If you are pinned to a specific bwameth release, re-validate against your own
> bwameth output** — see [bwameth.py drop-in mapping](bwameth-mapping.md) for the
> full caveat.

**Effect on output:**

The mode changes alignment score, `MAPQ`, `NM`, `MD`, and occasionally placement
and CIGAR. On a real C/T (or G/A) variant under a seed, `genomic` lowers the
score by `-A + -B` (the match score plus the mismatch penalty) relative to
`collapsed` — the freed match becomes a mismatch — which can break paralog ties
in `genomic`'s favor and avoid spurious indels. `neutral` scores the conversion
at `0`, i.e. `-A` relative to `genomic`, so a converted base neither helps nor
hurts a seed the way a match or mismatch would; it is the strictest of the three.
The Bismark `XR`/`XG`/`XM` tags and the `SEQ` field are identical in all modes.

**When to use it:**

Keep `collapsed` for methylation-only workflows that must match bwameth placement
(e.g. clinical pipelines validated against a bwameth release). Choose `genomic`
when you want one BAM that serves both methylation *and* variant calling, or want
the aligner to distinguish real variants from conversions. Prefer `neutral` for
TAPS (its default), where conversions are sparse and rewarding them as matches
costs placement specificity — it keeps variants visible for variant calling on the
same terms as `genomic` while placing sparse-conversion reads more accurately.

> **Note — `-B` follows the mode, but you can override it**
>
> `collapsed` sets `-B 2`; `genomic` and `neutral` keep `-B 4` by default. An
> explicit `-B` overrides the mode default and still reaches the per-strand
> matrices, whether it appears before or after `--meth`. The other `--meth`
> defaults (`-L 10 -U 100 -T 40 -M -C`) are the same across modes.

## `--meth-tags SPEC`

Selects which of the Bismark tags — [`XR:Z`, `XG:Z`, `XM:Z`](tags.md) — are
emitted. **Default: `all`.** A tag that is not selected is never computed, so
deselecting `XM` skips its per-read pass as well as its bytes.

| Spec | Meaning |
|------|---------|
| `all` | emit `XR`, `XG`, `XM` (the default) |
| `none` | emit none of the three |
| `XR,XG` | inclusion list — emit exactly these |
| `^XM` | exclusion — everything except `XM` |

A spec is either an inclusion list or a set of `^` exclusions; mixing the two
(`XR,^XM`) is rejected rather than given an arbitrary reading. Tag names are
case-insensitive.

```bash
# Drop the methylation call string; keep the strand labels.
bwa-mem3 mem --meth --meth-tags '^XM' ref.fa r1.fq r2.fq > out.bam

# Identical, spelled as an inclusion list (no quoting needed).
bwa-mem3 mem --meth --meth-tags XR,XG ref.fa r1.fq r2.fq > out.bam
```

> **Quote `^`-prefixed specs.** In `zsh` with `EXTENDED_GLOB` enabled — the
> default under oh-my-zsh and prezto — a bare `^XM` is a *negated glob* and
> expands to every file in the working directory except `XM`. The aligner then
> receives a filename as the tag spec and the remaining filenames as positional
> arguments. `bash` and stock `zsh` leave `^` alone, so this bites only some
> users; quoting it (`'^XM'`) is always safe, and `XR,XG` sidesteps it entirely.

> **The list is comma-separated, not space-separated.** `--meth-tags XR XG`
> passes only `XR` to the flag; `XG` becomes the next positional argument, i.e.
> the reference. The failure surfaces as a confusing
> `--meth seed index 'XG.meth.*' not found` rather than a flag error. This is the
> same shape as the `--meth taps` trap above, though the cause differs:
> `--meth` takes an *optional* argument and so requires `=`, whereas
> `--meth-tags` takes a *required* argument and happily consumes the first word
> but not the second.

**Why you would drop `XM`.** It is a read-length string, so it dominates the
aux payload: on a 151 bp EM-seq-style library it is **~33 % of the BAM**
(~155 B/read uncompressed, ~42 B/read at BGZF level 6), while `XR` and `XG`
together cost ~1 B/read. The alignments themselves are unaffected either way —
same records, positions, CIGARs, and scores — but dropping `XM` removes both its
bytes and the per-read methylation-call pass that produces them.

**Why the default keeps it.** `XM` is what makes the BAM directly consumable by
the Bismark family — `bismark_methylation_extractor`, methylKit, methtuple,
DMRfinder, epialleleR. Those tools cannot reconstruct it. Callers that work
from the reference FASTA instead — **MethylDackel** and **biscuit**, the most
common pairing for a bwameth-style BAM — never read `XM`, so `^XM` is free for
them.

The asymmetry is why `all` remains the default: a pipeline that needed `XM` and
silently lost it produces empty or wrong methylation calls, whereas a pipeline
that did not need it merely carries a larger BAM and can opt out with one flag.

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

## `-M` and split alignments

`--meth` sets `-M` and `-C` unconditionally (bwa has no option that unsets
either), so unlike the scoring defaults these two cannot be overridden.

`-M` changes how a split (chimeric) hit is **flagged**, not whether it is
emitted: supplementary (`0x800`) becomes secondary (`0x100`). The record is
otherwise identical — `SA:Z`, SEQ, QUAL, `NM`/`MD` and the methylation tags are
all retained. So `samtools flagstat` on a `--meth` BAM always reports
`0 supplementary`, and anything selecting split-read evidence by the `0x800` bit
will miss it. Methylation calling is unaffected.

> **Note — on short reads, `-T 40` filters far more split hits than `-M` does**
>
> A split arm covers only part of the read, so its score scales with arm length;
> below `T` it is dropped outright, before `-M` is consulted. On 1M Twist EM-seq
> pairs (75 bp) vs hg38, `--meth` emitted 94 split records against 4039 at
> `-T 30` — and of those 4039, 3945 scored `AS < 40` versus 94 at `AS >= 40`,
> exactly the default-threshold survivors. (`-L 10` accounts for ~18 more.)
>
> The effect shrinks with read length — a 150 bp read splits into ~75 bp arms
> scoring well clear of 40 — though that has not been measured. Passing `-T 30`
> (optionally `-L 5`) alongside `--meth` recovers those dropped split records,
> at the cost of divergence from bwameth placement.
>
> **This is not a substitute for supplementary alignments.** `-M` is still
> unconditional, so the recovered records are flagged secondary (`0x100`), not
> supplementary (`0x800`) — an SV caller that selects split-read evidence by the
> `0x800` bit will ignore them no matter how low `-T` goes. Lowering `-T` is
> useful for manual inspection, or for callers that accept `0x100` records and
> read `SA:Z` directly.

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
