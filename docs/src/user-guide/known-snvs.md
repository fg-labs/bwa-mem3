# Counting alleles at known SNVs

`bwa-mem3 aln --known-snvs FILE` places short reads as if the genome held both alleles of every
SNV in FILE, and tags each read that covers one of those sites with `XV:Z`. This page explains
what the tag records and which reads to count when you want allele counts without reference
bias. The option itself is described in [aln](../cli/aln.md#known-snvs---known-snvs).

## Why short reads need it

A read that carries the ALT allele of an SNV has one more mismatch against the reference than the
same read carrying REF. For a 100-base read that hardly matters. For a read of 20-30 bases,
`aln` allows only about two differences, so the extra mismatch can:

- push the read over its mismatch budget, leaving it unmapped;
- make it tie with a paralog that matches it equally well, giving MAPQ 0;
- drop its MAPQ below the usual cut-off of 20.

REF reads at the same site are not affected, so a plain count of REF and ALT reads at MAPQ 20
under-counts ALT. On a synthetic set of 457k read pairs on chr1, where each read carried REF or ALT
at a known SNV, stock `aln` misplaced 92,416 ALT mates at MAPQ 20 or above (compared to aligning
each read to a genome that carries its true allele). With `--known-snvs` that fell to 185, and the
REF reads were unchanged.

Placing ALT reads correctly is half the job. The other half is choosing which reads to count, and
that is what `XV:Z` is for.

## Reading `XV:Z`

`XV:Z` has one entry for each site in FILE that the read's alignment covers, separated by `;`:

```text
contig:pos:REF>ALT,allele,away,max_diff,twin_same,twin_mapq
```

| Field | Meaning |
|---|---|
| `contig:pos:REF>ALT` | The site, with a 1-based position. |
| `allele` | The read's base at the site: `R` (REF), `A` (ALT) or `O` (another base). Under `--meth`, `H` means the conversion could turn either allele into the base the read shows, so the read cannot tell them apart. |
| `away` | The read's differences (mismatches, gap opens and gap extensions) outside this site. |
| `max_diff` | The read's mismatch budget, set by its length and `-n`. |
| `twin_same` | 1 if the read's twin lands on the same strand within 5 bases of the read, else 0. |
| `twin_mapq` | The twin's MAPQ. |

The **twin** is a copy of the read with REF and ALT swapped at the site. An ALT read's twin
carries REF and a REF read's twin carries ALT. `aln` aligns the twin the same way as the read and
reports where it lands. For paired input the twin is aligned as a pair, with a copy of the mate
(with the mate's base swapped too if it covers the same site). The twin only adds these two fields:
it never changes the read's own record.

The twin fields are `NA,NA` for:

- every entry after the first `R` or `A` entry, because each read is twinned at one site only;
- `O` and `H` entries;
- reads that keep their stock result because they have more than 256 candidate loci or are longer
  than about 160 bases.

CIGAR, MD and NM are still computed against the reference, so an ALT read shows its mismatch
there. `XM:i` gives the read's mismatch count with the ALT treated as a match.

### Examples

Three 24-base reads over the same site. Each has a mismatch budget of 2:

```text
XV:Z:chr7:1048576:C>T,A,0,2,1,37
XV:Z:chr7:1048576:C>T,R,1,2,1,37
XV:Z:chr7:1048576:C>T,A,2,2,1,25
```

- The first read carries ALT, has no other differences, and its REF twin lands in the same place at
  MAPQ 37. Count it as ALT.
- The second carries REF with one difference elsewhere, and its ALT twin places too. Count it as
  REF.
- The third carries ALT, but it already has two differences outside the site, its whole budget.
  Do not count it (see the third rule below).

A read over two sites:

```text
XV:Z:chr1:200100:G>A,R,0,2,1,30;chr1:200112:A>G,A,0,2,NA,NA
```

It can be counted at `chr1:200100` only. It has no twin at `chr1:200112`, so it is not counted
there.

## Which reads to count

Count a read at a site when all of the following hold:

1. **Allele `R` or `A`.** `O` is neither allele, and `H` cannot tell them apart.
2. **The read has MAPQ at or above your cut-off** (20 is usual).
3. **Its twin also passes:** `twin_same` is 1 and `twin_mapq` is at or above the same cut-off.
4. **`away` is at most `max_diff - 1`.**

Each rule removes one source of bias.

- **The twin rule** makes placement independent of the allele. A read is counted only if it would
  have placed at the same locus, just as confidently, carrying the other allele. This keeps out
  sites where one allele is unique but the other matches a paralog, and it also removes the
  reference bias left in stock placements: REF reads whose ALT twin would not have placed are
  dropped too.
- **The budget rule** stops sequencing errors from becoming ALT counts. A read with `max_diff`
  differences elsewhere is already at its budget. Stock `aln` drops it whatever base sits at the
  site. With `--known-snvs` it places only when the base at the site is the ALT, and at that depth
  a base that matches the ALT is often a sequencing error. Its REF twin has the same differences
  elsewhere and places as well, so the twin rule cannot catch it. Requiring one spare difference
  means the base at the site never decides whether the read places.

All four rules treat a REF read and an ALT read the same way, so they do not favour either allele.
The cost is depth: the budget rule removes a few percent of the reads that would otherwise be
counted. Sites in repeats, where neither allele places
uniquely, end up with few or no counted reads instead of a biased count.

## Paired-end reads

Each mate has its own `XV:Z`. Count fragments, not mates:

- If only one mate covers the site, the fragment's count is that mate's.
- If both cover it, count the fragment once, and only when both mates pass the rules above. If
  the mates disagree at the site, the fragment is `O`.

The twin of a paired read is a pair, so `twin_mapq` is the twin pair's paired MAPQ.

## What it does not handle

- **Nearby sites.** Each window carries one site's ALT at a time, so a read that carries the ALT at
  two sites is credited for one of them only, and may still be misplaced.
- **Long reads.** Reads longer than about 160 bases keep their stock result. They rarely need
  crediting.
- **Paralogous sites.** If a site's ALT allele matches a paralog somewhere else in the genome, reads
  from the paralog carry the ALT. Crediting breaks their tie toward the site, and their REF twin is
  unique there, so no read-level rule removes them. This leaves a small excess of ALT at such sites.
- **Indels and MNVs.** Only single-base substitutions in FILE are used.

At the end of a run `aln` prints how many reads were credited, how many were at a site with no
change, and how many kept their stock result:

```text
[M::aln] --known-snvs: 1203 reads credited, 51877 at a site with no change, 12 kept stock (too many hits or too long)
```

---

**See also:**
[aln](../cli/aln.md#known-snvs---known-snvs) ·
[Output: SAM/BAM, headers, tags](output.md) ·
[Methylation SAM tags](../methylation/tags.md)
