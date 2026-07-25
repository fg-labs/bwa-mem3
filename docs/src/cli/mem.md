# mem

`bwa-mem3 mem` aligns short DNA reads against an indexed reference genome
using the BWA-MEM algorithm. It accepts one or two FASTQ files (single-end or
paired-end) and writes alignments to stdout in SAM or BAM format. It is the
primary alignment subcommand; nearly all bwa-mem3 usage flows through it.

## Synopsis

```text
{{#include ../../_generated/cli/mem.txt}}
```

## Common usage

Paired-end alignment, 16 threads, SAM to stdout:

```bash
bwa-mem3 mem -t 16 ref.fa R1.fq.gz R2.fq.gz > out.sam
```

Paired-end alignment, emit uncompressed BAM, pipe directly to `samtools sort`:

```bash
bwa-mem3 mem --bam -t 16 ref.fa R1.fq.gz R2.fq.gz \
  | samtools sort -@ 8 -o out.bam -
samtools index out.bam
```

Paired-end methylation alignment with a read group header:

```bash
bwa-mem3 mem --meth -t 16 \
  -R '@RG\tID:lib1\tSM:sample1\tPL:ILLUMINA' \
  ref.fa R1.fq.gz R2.fq.gz \
  | samtools sort -o out.bam -
```

## Flag reference

### Input / output

#### `-o STR` — output file

Write output to `STR` instead of stdout. Honored for both SAM and `--bam`
output; the path is opened lazily so BAM mode can hand it to htslib instead of
truncating it as a SAM-text file. Stdout redirection (`>`) remains an
alternative.

#### `--bam[=N]` — emit BAM

Emit BAM instead of SAM. `N` controls BGZF compression: `0` (default when
`--bam` is used without `=`) writes uncompressed BAM, which costs almost no
CPU and is the recommended mode for piping to `samtools sort`. Values `1`–`9`
select increasing BGZF deflate levels; use `--bam=6` or `--bam=9` only when
writing directly to final storage without a downstream sort step.

> **Tip — Prefer --bam for production pipelines**
>
> Uncompressed BAM (`--bam` or `--bam=0`) eliminates the text-formatting cost on
> the aligner side and the text-parse cost on the `samtools sort` side. For any
> pipeline that immediately sorts or processes the output, this is faster than
> SAM at no quality cost.

#### `--compat=TARGET` — byte-identical output for another aligner

Shapes bwa-mem3's output so it is byte-for-byte identical to another aligner's.
`--compat` takes a required target:

| target | meaning |
|---|---|
| `bwa-mem2` (alias `mem2`) | match bwa-mem2 v2.2.1 |
| `off` | bwa-mem3's native output — identical to omitting the flag |

```bash
bwa-mem3 mem --compat=bwa-mem2 -t 16 ref.fa R1.fq R2.fq > out.sam
bwa-mem3 mem --compat bwa-mem2 -t 16 ref.fa R1.fq R2.fq > out.sam   # same
```

Under `--compat=bwa-mem2`, bwa-mem3 suppresses exactly the differences that are
bwa-mem3-side additions, and nothing else:

- **`MQ:i`** (mate MAPQ) and **`HN:i`** (hit count) are not emitted. Note only
  `HN` is genuinely a bwa-mem3 invention — bwa emits `MQ` too
  ([lh3/bwa#330](https://github.com/lh3/bwa/pull/330), merged 2022-03-06);
  bwa-mem2 lacks it only because it forked at 0.7.17, before that landed.
- **The default `@HD` line** is not emitted. bwa-mem2 v2.2.1 writes no `@HD` at
  all; bwa gained one in 0.7.18 (`6b18630`), again after the fork. A user `@HD`
  from `-H` is still honored.
- **The `<prefix>.hdr` / `<baseprefix>.dict` sidecar is not read**, so `@SQ` is
  generated from the index as bare `SN`/`LN` (plus `AH:*` on ALT contigs, which
  bwa-mem2 does emit). The sidecar is a bwa-mem3-only feature — a port of
  [lh3/bwa#348](https://github.com/lh3/bwa/pull/348), which upstream closed
  unmerged — so its `M5`/`AS`/`UR`/`SP` have no counterpart to match.

`--compat` is **output-shaping only** — it changes no alignment, no score, no
flag, and no tag's value. It exists so pipelines validating bwa-mem3 against a
bwa-mem2 golden can diff raw output without a post-processing step.

Verify it with:

```bash
# records
diff <(samtools view out.compat.bam) <(samtools view out.mem2.bam)
# header, @PG excluded on both sides (see below)
diff <(samtools view -H --no-PG out.compat.bam | grep -v '^@PG') \
     <(samtools view -H --no-PG out.mem2.bam   | grep -v '^@PG')
```

Notes and caveats:

- **`@PG` is still emitted, by design.** bwa-mem2 writes its own `@PG`, so
  suppressing ours would turn a *changed* line into a *missing* one — the diff is
  non-empty either way, and suppression would only cost the record of what
  actually produced the file. `@PG` is unmatchable by construction regardless:
  `CL:` embeds the invocation and its paths, so even two bwa-mem2 runs from
  different directories differ. Exclude `@PG` on both sides when comparing.
- **Mutually exclusive with `--fast`.** Passing both is a hard error
  (`--compat and --fast are mutually exclusive`), because `--fast` deliberately
  *changes alignments* while `--compat` only shapes output — so `--fast
  --compat=bwa-mem2` would produce a diff-clean-*looking* stream over genuinely
  different alignments, defeating the parity-validation purpose. `--compat` is
  for the drop-in profile. (`--compat=off --fast` is fine: `off` selects no
  target.)
- **Non-`--meth` only; combining them is a hard error**
  (`--compat is not supported with --meth`). bwa-mem2 has no bisulfite mode, so
  byte-identity is undefined under `--meth`, which also emits
  methylation-specific tags that no target models.
- **`--compat=bwa-mem` is recognized but rejected.** bwa-mem3 and bwa still
  differ on 224 of 63,583 records (53 above MAPQ 30) in MAPQ, CIGAR, POS and
  TLEN, so byte-identity is not achievable no matter what output is suppressed.
  The target's output shaping is already specified in `src/compat_target.cpp`;
  it becomes selectable when the alignment work lands.

See [Equivalence with bwa-mem2 → Byte-identical output
(`--compat`)](../whats-different/equivalence.md#byte-identical-output---compat).

#### `-R STR` — read group header

Injects a `@RG` header line and tags every alignment with `RG:Z:<ID>`. The
value is a tab-separated `@RG` line with literal `\t` escapes, for example:

```bash
-R '@RG\tID:run1\tSM:HG001\tPL:ILLUMINA\tLB:lib1'
```

bwa-mem3 escapes any literal tab characters inside `-R` values before writing
them to the `@PG CL:` field, preventing header corruption (fix for issue #45).

#### `-H STR/FILE` — extra header lines

If `STR` begins with `@`, it is injected verbatim as a header line. Otherwise
`STR` is treated as a path and every line in the file is injected. Useful for
adding `@CO` comments or custom `@RG` / `@PG` entries.

#### `-p` — smart pairing

Reads interleaved paired-end data from a single FASTQ file (`in1.fq`) rather
than two separate files. The second positional argument (`in2.fq`) is ignored.

#### `-5` — leftmost-coordinate primary

For split alignments, designates the alignment with the smallest genomic
coordinate as primary, rather than the longest alignment. Useful for some
downstream tools that expect the leftmost alignment to be primary.

#### `-q` — preserve supplementary MAPQ

By default, bwa-mem3 may downgrade the MAPQ of supplementary alignments.
`-q` suppresses that adjustment.

#### `-K INT` — fixed batch size

Forces each thread batch to process exactly `INT` input bases regardless of
the number of threads. Useful when you need bit-for-bit reproducible output
across runs with different `-t` values: fix `-K` to the same value and the
output is deterministic.

#### `-v INT` — verbosity

Controls stderr diagnostic output: `1` = errors only, `2` = warnings,
`3` = informational messages (default), `4+` = debugging.

#### `-a` — all alignments

Output all alignments for single-end or unpaired paired-end reads, including
secondary alignments. Equivalent to enabling secondary-alignment reporting.

#### `-C` — append FASTA/FASTQ comment

Appends the comment field from the FASTA/FASTQ header to the SAM output as
an additional column. Useful when the comment carries barcodes or UMIs.

#### `-V` — reference header in XR tag

Emits the reference FASTA header line for each alignment position as an `XR`
SAM tag.

Under `--meth`, `XR:Z` instead carries the Bismark read-conversion direction
(`CT`/`GA`) and this reference-annotation use of `XR` is suppressed — see
[Methylation Reference → Flags](../methylation/flags.md).

#### `-Y` — soft-clip supplementary alignments

Uses soft clipping instead of hard clipping for supplementary alignments.
Some downstream tools require this.

#### `-M` — mark shorter split hits as secondary

Marks the shorter alignment in a split read as secondary (sets `0x100` flag)
rather than supplementary. Required for compatibility with tools that do not
handle supplementary alignments (e.g. Picard's duplicate-marking before
certain versions).

#### `-j` — treat ALT contigs as primary

Treats ALT contigs as part of the primary assembly by ignoring the
`<idxbase>.alt` file. Use when your workflow does not include ALT-aware
postprocessing.

### Scoring

All scoring flags accept integer values. Changing `-A` (match score) scales
the penalty flags that default to multiples of `-A`; explicit overrides of
individual flags are unaffected.

| Flag | Default | Meaning |
|------|---------|---------|
| `-A INT` | 1 | Score for a sequence match. Scales `-T`, `-d`, `-B`, `-O`, `-E`, `-L`, `-U` unless overridden. |
| `-B INT` | 4 | Mismatch penalty. |
| `-O INT[,INT]` | 6,6 | Gap open penalty for deletions and insertions respectively. |
| `-E INT[,INT]` | 1,1 | Gap extension penalty per base. A gap of length k costs `-O + -E * k`. |
| `-L INT[,INT]` | 5,5 | Clipping penalty for 5' and 3' ends. |
| `-U INT` | 17 | Penalty for an unpaired read pair (affects mate-rescue scoring). |
| `-T INT` | 30 | Minimum alignment score to output. Alignments below this threshold are not reported. |

> **Note — --meth overrides scoring defaults**
>
> When `--meth` is active, bwa-mem3 applies `-L 10 -U 100 -T 40 -M -C` plus a
> mode-dependent mismatch penalty: `-B 2` for `--meth-scoring collapsed` (the
> `--meth`/`--meth=emseq` default, bwameth-compatible) and `-B 4` for
> `--meth-scoring genomic` and `neutral` (`neutral` is the `--meth=taps` default).
> This mirrors bwameth's `bwa mem -T 40 -B 2 -L 10 -CM` (with `-U 100` for
> paired-end). The scoring values (`-B`, `-L`, `-U`, `-T`) can still be overridden
> by passing the flag explicitly, in any position relative to `--meth`; `-M` and
> `-C` cannot, since bwa has no option that unsets them.
>
> Those constants are quoted at bwameth's match score (`-A 1`) and scale with
> `-A` like every other score-derived default above: under `-A 2` the effective
> values are `-L 20 -U 200 -T 80` and `-B 4` (collapsed) / `-B 8` (genomic and
> neutral).

### Paired-end

#### `-I FLOAT[,FLOAT[,INT[,INT]]]` — insert size distribution

Specifies the mean, standard deviation (default: 10% of mean), maximum
(default: 4 sigma above mean), and minimum of the insert size distribution for
FR-orientation paired-end reads. By default bwa-mem3 infers these parameters
from the first batch of reads. Provide them explicitly for speed or when the
reference is short and inference may be inaccurate.

#### `-m INT` — mate rescue rounds

Maximum number of mate-rescue attempts per read. Reduce to speed up alignment
on data where the default (50) wastes time on unrescuable pairs. See
[Settings profiles](../best-practices/settings-profiles.md) for the benchmarked
`-m 10` recommendation.

#### `-S` — skip mate rescue

Disables mate rescue entirely. Faster but may reduce sensitivity for
discordant pairs.

#### `-P` — skip pairing

Skips the pairing step; mate rescue still runs unless `-S` is also given.

### Filtering

#### `-c INT` — skip repetitive seeds

Seeds with more than `INT` occurrences in the reference are skipped. Lowering
this (e.g. to 50) speeds up alignment of highly repetitive reads but may
reduce sensitivity. Raising it increases sensitivity in repeat-heavy regions
at a cost in runtime.

#### `-D FLOAT` — chain length fraction

Drops chains shorter than `FLOAT` times the longest overlapping chain. The
default (0.50) discards chains that are less than half the length of the best
chain.

#### `-W INT` — minimum seeded bases

Discards chains with fewer than `INT` seeded bases. Raising this filters out
very short, low-confidence chains.

#### `--min-ext-len INT` — skip Smith-Waterman extension of short seeds

Off by default (`0`) → output byte-identical to baseline. When `INT > 0`, a short
seed (< `INT` bp) is dropped before banded Smith-Waterman **only if its chain
still has a longer anchor seed** — its extension is then redundant (the anchor
already covers it), so skipping it is near output-neutral (~10 % less alignment
CPU at `30`). A chain whose seeds are *all* short is left untouched, so the filter
never empties a chain or drops a read: it is recall-safe by construction. `30` is
the recommended value. For the benchmarks, behavior details, and validation
status, see
[Settings profiles → `--min-ext-len 30`](../best-practices/settings-profiles.md#short-seed-extension---min-ext-len-30).

#### `--max-extend-chains INT` — cap chains extended per read

Off by default (`0`) → output byte-identical to baseline. When `INT > 0`, only the
top-`INT` chains by weight (after chain filtering) reach banded Smith-Waterman
extension; the remaining lower-weight chains are dropped before extension. This is
the only lever that reduces the *number of chains extended* per read, so it is
orthogonal to the seed- and SW-per-chain levers and adds a real marginal speedup on
top of them (~15 % marginal alignment CPU on top of `--fast`, ~23 % standalone, at
`5`). It is **not** byte-identical: dropping candidate chains removes low-weight
secondaries, so `XS`, secondary alignments, and `MAPQ` can shift on multi-mapping
reads. High-confidence (uniquely-placed) reads are unaffected. The cap is a safety
no-op for pathological reads with more than 4096 chains (`MAX_EXTEND_CHAINS_CAP`):
those reads extend all of their chains as usual, so `--max-extend-chains` has no
effect on them. `--fast` sets `20`. For the accuracy/speed curve and validation status,
see
[Settings profiles → `--max-extend-chains`](../best-practices/settings-profiles.md#chain-extension-cap---max-extend-chains).

#### `--adaptive-band` — adaptive banded Smith-Waterman for medium-length reads

> **Kilobase-scale long reads are not currently supported at default settings.**
> `bwa-mem3 mem` on PacBio HiFi (~15–20 kb) or ONT reads exhausts the
> short-read-sized extension buffers and aborts (or, before the guard, OOM'd /
> segfaulted) — the buffer model cannot represent a kilobase-scale extension.
> This lever helps the **medium-length** end of its range (e.g. SBX at ~240 bp,
> where it is validated); it does **not** make HiFi/ONT usable. For a true
> long-read workload use minimap2. Tracking issue: bounding long-read extension
> memory (fg-labs/bwa-mem3#238).

Off by default → output byte-identical to baseline. When set, banded extension
starts at a tight band and expands each pair only to the band its chain's seed
geometry actually needs (the inter-seed indel), rather than the fixed `-w` band
(100) for every extension.

**When to use it: medium-length reads.** The band only constrains the DP matrix
when the extension's reference window exceeds it (`ref_window > 2·w+1`), which
happens once reads are roughly ≥ 200 bp. So this is a lever for the
**medium-length** range — SBX (HG002, ~240 bp), where it cuts alignment CPU by
**~25 %**. It nominally applies to PacBio HiFi and ONT too, but kilobase-scale
reads do not run at default settings (see the caveat above), so in practice its
usable range today is medium reads, not true long reads. On short-read data (WGS ~150 bp, WES ~76 bp) the extension matrix is
already smaller than the band, so there is nothing to trim: those reads run on the
8-bit kernel, which this option deliberately leaves untouched, making it a **no-op
on short reads** (enabling it on a WGS/WES run neither helps nor hurts).

**Accuracy:** placement is unchanged (holodeck `sim-wgs-place`: MAPQ-60+ mismaps
identical to default) and indel representation is preserved — indels up to the
chaining limit still emit a single `D`/`I` CIGAR, matching the `-w 100` default,
so small/mid-size indel callability is unaffected.

**Not byte-identical when on.** Like `--fast`, enabling it shifts a small number of
borderline secondary alignments (starting tight and expanding can change which of
several near-tied placements wins). It is therefore an opt-in flag, not a default.

#### `--extend-mate-concordant` — retain mate-concordant chains under a chain cap

Takes an optional window: `--extend-mate-concordant` (bare) = **auto**, sizing the
window to the estimated proper-pair insert bound (`pes[FR].high`, inferred from the
data during the run); `--extend-mate-concordant=INT` pins a fixed window in bp;
`--extend-mate-concordant=0` disables it. Off by default → no effect. When on (and
`--max-extend-chains` is capping a paired-end read), a chain that would be dropped by
the cap is instead **retained if it is concordant with one of the mate's chains** —
same contig, FR ("innie") orientation, within the window. It only does anything when
a chain cap is in effect, so it is a strict no-op without `--max-extend-chains`.

The window matters: too wide and it retains — and then extends — far/spurious
concordant chains, adding alignment CPU on chain-rich reads; sizing it to the
aligner's own proper-pair insert bound (the auto default) admits only genuine pair
anchors. Before the insert size is estimated (the first chunk), auto falls back to a
built-in default.

**When to use it: `--meth`.** Bisulfite's collapsed 3-letter alphabet flattens chain
weights, so under `--max-extend-chains` the cap often drops a read's true low-weight
chain and starves PE pairing of the anchor that lets the true concordant pair win —
flipping both mates to a wrong concordant locus. This option recovers that anchor.
`--fast` enables it (auto) automatically **under `--meth` only**; on non-meth data the
cap does not regress placement, so `--fast` leaves it off to preserve the speedup.
The recovery is **partial** (it narrows, but does not fully close, the placement gap
to default), and with the auto window the alignment-CPU cost is **~1%** — sizing the
window to the insert bound is what keeps it there (a wide fixed window instead retains
and extends far/spurious concordant chains, costing 15–20% on chain-rich reads). See
the benchmarked per-dataset figures in
[fg-labs/bwa-mem3#195](https://github.com/fg-labs/bwa-mem3/pull/195).

**Not byte-identical when it retains a chain.** Like `--max-extend-chains`, keeping an
extra candidate can move `XS`, secondaries, and `MAPQ` on multi-mapping reads;
high-confidence placement is unaffected. For the placement/mismap validation, see
[Settings profiles → `--extend-mate-concordant`](../best-practices/settings-profiles.md#mate-concordant-chain-retention---extend-mate-concordant).

#### `-h INT[,INT]` — secondary alignment reporting

If there are fewer than `INT` hits with score exceeding `FLOAT` (see `-z`)
times the maximum score, all of them are output in the `XA` auxiliary tag.
The second integer is a hard cap on the number of XA entries. Defaults: 5, 200.

#### `-z FLOAT` — secondary score fraction

Fraction of the maximum alignment score used as the threshold for secondary
hit reporting with `-h`. Default: 0.80.

#### `-u` — emit XB instead of XA

Outputs `XB` in place of `XA`. `XB` is an extension of `XA` that also carries
the alignment score and mapping quality for each secondary hit.

### Speed preset

#### `--fast` — speed preset (opt-in, not byte-identical)

`--fast` is a one-flag shorthand for the characterized speed levers:

```text
bwa-mem3 mem --fast  ≡  -m 10 -y 0 --min-ext-len 30 --smem-dedup --skip-contained-ext --max-extend-chains 20 --adaptive-band --extend-mate-concordant
```

`--skip-contained-ext` is byte-identical to the default on non-meth short- and
medium-length single- and paired-end reads (validated on WGS/WES/HiC and SBX, and on an
ALT/HLA-enriched set) and no-ops under `--meth` (via its own internal gate). It is **not**
byte-identical on kilobase-scale long reads (PacBio HiFi, ONT): the "longest contained seed
dominates" assumption its skip relies on breaks when extensions terminate early on such
reads, dropping supplementary alignments and occasionally degrading a primary. Those inputs
are in any case not currently usable at default settings (see the long-read caveat under
[`--adaptive-band`](#--adaptive-band--adaptive-banded-smith-waterman-for-medium-length-reads)).

`--adaptive-band` (see above) is included because it is a strict no-op on short reads
(the reads `--fast` primarily targets) and a ~25% alignment-CPU speedup on medium-length
runs (e.g. SBX ~240 bp), so bundling it only helps. Note kilobase-scale HiFi/ONT reads do
not run at default settings regardless of this flag (see
[`--adaptive-band`](#--adaptive-band--adaptive-banded-smith-waterman-for-medium-length-reads)).

`--extend-mate-concordant` repairs the chain-cap pairing regression — the true, low-weight but
mate-concordant chain the cap would otherwise drop — and is included for both non-meth and `--meth`
`--fast` (see
[Settings profiles → `--extend-mate-concordant`](../best-practices/settings-profiles.md#mate-concordant-chain-retention---extend-mate-concordant)).

`--fast` also flips one lever that has **no flag of its own**: the sort that orders alignment
regions by reference end position before the order-sensitive dedup/patch pass. The default
build uses bwa-mem2's comparator — a partial order on the end position alone — and reproduces the
permutation bwa-mem2's unstable introsort leaves, and so its dedup outcome, exactly. `--fast`
switches to a strict total order (end position, then start, score and query start) under pdqsort,
which can resolve equal-end-position ties differently. The two orderings only diverge when the
strict order's deterministic result differs from the permutation that unstable sort leaves;
because bwa-mem2's output is defined by that permutation, a strict total order and exact
bwa-mem2 parity cannot both hold in general. The
effect is confined to reads carrying several regions that end
at the same reference base, and it surfaces as different `XS`/MAPQ and occasionally a different
primary or supplementary record. Runs report the lever as `alnreg-sort=fast` on the audit line
below.

The lever is no longer worth much speed. The default path only needs the unstable introsort when
the comparator actually sees a tie — with no tie the ordering is unique, so pdqsort's result is
what introsort would have produced — and ties occur on under 1% of calls. Detecting them costs a
save-copy and a linear scan, which leaves the bwa-mem2-compatible default within ~0.7% of the
strict-total-order path on 5M HG00096 WGS pairs (arm64/NEON, `-t 16`, 5 interleaved reps), at or
below this benchmark's run-to-run spread.

Under `--meth` it additionally sets `-s 2` (light Pass-2 re-seeding) and lowers the chain cap to `10`.
Earlier releases used `-s 0` (no re-seed), which inflated
MAPQ on bisulfite reads; `-s 2` recovers the MAPQ/placement at nearly the same speed
(see [Settings profiles → Pass-2 re-seeding](../best-practices/settings-profiles.md#pass-2-re-seeding-under---meth--s-2)).

Each lever is applied only if you did not set it explicitly, so explicit flags
win where applicable (`--fast -m 30` keeps `-m 30`; `--fast --max-extend-chains 8`
keeps `8`); `--smem-dedup` and `--skip-contained-ext` are always enabled and
cannot be opted back out of once
`--fast` is set. Output is **not**
byte-identical to the default; the accuracy cost of each lever is characterized in
[Settings profiles](../best-practices/settings-profiles.md) and is confined to
the already-low-confidence tail. `bwa-mem3 mem` prints the resolved preset to
stderr (`[M::main_mem] --fast: ...`) so runs are self-documenting.

### Methylation (`--meth`)

#### `--meth` — enable bisulfite alignment mode

Activates bisulfite alignment: each read is projected (R1 `C→T`, R2 `G→A`) to find
seeds in the converted `.meth` seed index, then extended and scored against the
**original** 4-letter reference, with inline BAM post-processing and forced
`--bam` output. The reference must have been indexed with `bwa-mem3 index --meth`.

Pass the original FASTA prefix as `<idxbase>` (e.g. `ref.fa`); the `ref.fa.meth.*`
seed index alongside it is found automatically. A legacy bwameth `.bwameth.c2t`
index is not used directly — rebuild with `index --meth` (see
[Migrating from bwameth.py c2t](../methylation/external-c2t.md)).

See [Methylation Reference](../methylation/overview.md) for the full treatment.

#### `--meth-scoring {collapsed|genomic|neutral}` — bisulfite/TAPS scoring model

Selects how the 4-letter matrix treats converted bases. `collapsed` (the
`--meth`/`--meth=emseq` default) frees C↔T and G↔A both ways (bwameth-compatible
placement, sets `-B 2`); `genomic` frees only the conversion direction scored as
a full match, keeping real variants as mismatches (variant-aware, truthful
`NM`/`MD`, keeps `-B 4`); `neutral` (the `--meth=taps` default) frees only the
conversion direction but scores it `0` — tolerated, not rewarded — best for the
sparse conversions of TAPS (keeps `-B 4`, truthful `NM`/`MD`). Only meaningful
with `--meth`.
See [Flags → --meth-scoring](../methylation/flags.md#--meth-scoring-collapsedgenomicneutral).

#### `--set-as-failed {f|r}` — strand QC-fail flag

Forces the QC-fail bit (`0x200`) on all alignments to the forward (`f`) or
reverse (`r`) bisulfite strand. Used when one strand is known to be
unreliable for a given library preparation.

#### `--chimera-qc` — opt in to bwameth.py-style chimera heuristic

Off by default (matches Bismark, which has no equivalent heuristic).
When set, mapped records whose longest M/=/X CIGAR run is less than 44 % of
the read length get `0x200` set, `0x2` cleared, and MAPQ capped at 1. Useful
for PBAT / scBS-Seq libraries where intra-fragment chimerism is common, or
when reproducing bwameth.py output bit-for-bit.

### Threading

#### `-t INT` — number of threads

Number of worker threads. Defaults to 1. Set to the number of physical cores
available to this job. Scaling is workload- and hardware-dependent: on typical
machines the curve flattens around 16–32 threads (FM-index bandwidth and I/O
contention dominate); on high-memory / fast-I/O servers the aligner can keep
scaling toward ~64 threads on hg38 before saturating. See the threading guide
for measured guidance and per-machine recommendations.

See [User Guide — Threading and resource use](../user-guide/threading.md) for
guidance on thread counts at various machine sizes.

### Supplementary MAPQ rescoring

#### `--supp-rep-hard-cap INT` — cap MAPQ for repetitive supplementary alignments

Forces MAPQ=0 for supplementary alignments whose chain contains any seed with
at least `INT` occurrences in the genome. This targets supplementary
alignments anchored in repetitive regions that upstream MAPQ scoring may
overestimate. `0` disables the cap (default). Typical values are 5–20; lower
values are more aggressive. Primary alignment MAPQ is unaffected.

### Debug

#### `-k INT` — minimum seed length

Minimum exact-match seed length. Shorter seeds increase sensitivity but raise
runtime. The default (19) is calibrated for 100–150 bp Illumina reads.

#### `-w INT` — band width

Band width for the banded Smith-Waterman extension. Wider bands can recover
alignments with long indels at greater CPU cost.

#### `-d INT` — X-dropoff

Off-diagonal X-dropoff for the Z-drop heuristic. Controls how far an alignment
extension continues after a score drop.

#### `-r FLOAT` — re-seeding factor

Seeds longer than `-k * FLOAT` are re-seeded internally to find sub-seeds (bwa-mem's
second seeding round). Lowering produces more seeds / higher sensitivity at greater
cost; raising (e.g. `-r 10`) suppresses the round. Round 2 is genuine
split-read/divergence sensitivity, so only suppress it on known-clean data — see
[Settings profiles → `-y 0`](../best-practices/settings-profiles.md#third-round-seeding--y-0).

#### `-y INT` — third-round seed occurrence threshold

bwa-mem's third seeding round: for each read position, grow an exact match until it
occurs fewer than `INT` times in the genome (default 20), then emit it as a seed — a
repeat-region safety net. **`-y 0` disables the round**, cutting ~11–30 % of alignment
CPU with F1-near-neutral accuracy; it is part of the recommended profile. For the
regime sweep and rationale, see
[Settings profiles → `-y 0`](../best-practices/settings-profiles.md#third-round-seeding--y-0).

#### `--legacy-reader` — use the legacy input reader

Read input with the legacy `gzFile`/`kseq` reader instead of the default
content-detecting fast reader. An escape hatch for A/B baselining or working
around an input the fast reader mishandles; not needed in normal use.

## Notes / Gotchas

> **Warning — --meth requires a --meth index**
>
> Running `bwa-mem3 mem --meth` against a standard (non-c2t) index produces
> incorrect alignments without an error. Confirm that the index was built with
> `bwa-mem3 index --meth` before aligning bisulfite data.
>
> **Note — SIMD variant printed to stderr at startup**
>
> When mem starts it prints a banner (`Executing in AVX512 mode!!` etc.) to
> stderr. This is informational and does not affect stdout output.

---

**See also:**
[User Guide — Aligning short reads](../user-guide/aligning.md) ·
[User Guide — Output: SAM/BAM, headers, tags](../user-guide/output.md) ·
[CLI Reference — index](index-cmd.md) ·
[Methylation Reference — Overview](../methylation/overview.md) ·
[Best Practices — Output format](../best-practices/output-format.md)
