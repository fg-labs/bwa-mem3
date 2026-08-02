# Glossary

Terms used throughout this book, listed alphabetically.

---

**@HD header**
The first line of a SAM file header. Specifies the SAM format version (`VN`) and sort order (`SO`). Required when any other header lines are present. See [Output: SAM/BAM, headers, tags](../user-guide/output.md).

**@PG header**
A SAM header line recording a program that processed the file, including `ID`, `PN`, `VN`, and `CL` fields. bwa-mem3 inserts `ID:bwa-mem3` (or `ID:bwa-mem3-meth` in methylation mode). See [Output: SAM/BAM, headers, tags](../user-guide/output.md).

**@SQ header**
A SAM header line describing a reference sequence (chromosome). Contains the sequence name (`SN`) and length (`LN`). In methylation mode, `@SQ` is written directly from the original (un-converted) reference, so the `f`/`r`-prefixed contigs of the `.meth` seed index never appear in output. See [Chimera QC and header construction](../methylation/post-processing.md).

**BAM**
Binary Alignment Map — a compressed, binary encoding of SAM. Produced by bwa-mem3 when the `--bam` flag is given or when output is piped through samtools. See [Output: SAM/BAM, headers, tags](../user-guide/output.md).

**Banded Smith-Waterman (banded SWA)**
A heuristic variant of the Smith-Waterman alignment algorithm that restricts the dynamic programming to a band of width *w* around the main diagonal. bwa-mem3 uses banded SWA for extension alignment; bwa-mem2 kernels are SIMD-vectorized and bwa-mem3 adds NEON implementations for Apple Silicon. See [SIMD dispatch architecture](../developer-guide/simd-dispatch.md).

**c2t**
Cytosine-to-thymine in-silico conversion applied to reads (or reference) before methylation alignment. In `--meth` mode, bwa-mem3 converts R1 reads C→T and R2 reads G→A inline, without writing intermediate FASTQ files. See [Conversion details (C->T, G->A)](../methylation/conversion.md).

**Chimera**
A read alignment where the aligned portion is short relative to the read length, often indicating a mapping artefact or a true chimeric molecule. In methylation mode, bwa-mem3 applies a chimera QC heuristic: if the longest contiguous M/=/X CIGAR run is less than 44% of the read length, the alignment is flagged 0x200, the proper-pair bit is cleared, and MAPQ is capped at 1. See [Chimera QC and header construction](../methylation/post-processing.md).

**FASTQ**
A text format for raw sequencing reads. Each record contains a sequence identifier, the nucleotide sequence, a separator, and per-base quality scores in ASCII-encoded Phred format. bwa-mem3 accepts gzip-compressed FASTQ as input. See [Quick start: align paired-end FASTQs](../getting-started/quick-align.md).

**FM-index**
Ferragina-Manzini index — a full-text index over the Burrows-Wheeler Transform of a sequence. bwa-mem3 uses the compressed `.bwt.2bit.64` FM-index for seed finding (SMEM lookup). See [Indexing the reference](../user-guide/indexing.md).

**Hard clip**
A CIGAR operation (`H`) indicating that bases at the read end are absent from the SEQ field of the alignment record. Hard clipping is used in supplementary alignments to avoid duplicating the read sequence. See [Output: SAM/BAM, headers, tags](../user-guide/output.md).

**kswv**
The SIMD-vectorized kernel implementing the inner loop of the Smith-Waterman extension alignment in bwa-mem2/bwa-mem3. bwa-mem3 carries correctness fixes for the score-saturation edge case across all SIMD width variants (NEON, AVX2, AVX-512BW). See [Correctness fixes](../whats-different/correctness.md).

**libsais**
A library implementing the suffix-array induced sorting (SAIS) algorithm. bwa-mem3 optionally uses libsais for FM-index construction, reducing indexing time compared to the default suffix-array builder. See [Performance improvements](../whats-different/performance.md).

**LTO**
Link-Time Optimization — a compiler mode that defers optimization to link time, enabling cross-compilation-unit inlining. Activated via `make lto-build`. See [Building from source](../developer-guide/building.md).

**MAPQ**
Mapping quality — a Phred-scaled probability that a read alignment is incorrectly mapped. Reported in SAM field 5. bwa-mem3 follows bwa-mem2 MAPQ semantics; chimera QC in methylation mode caps MAPQ at 1 for chimeric alignments. See [Output: SAM/BAM, headers, tags](../user-guide/output.md).

**Mate rescue**
A step in paired-end alignment where, if one mate lacks a confident seed, bwa-mem3 attempts to find it by performing Smith-Waterman alignment in the region near the mapped mate. bwa-mem3 adds NEON and AVX2 implementations of the mate-rescue kernel. See [Architecture support](../whats-different/arch-support.md).

**mimalloc**
A high-performance memory allocator from Microsoft. bwa-mem3 vendors mimalloc and links it into every binary by default. To disable, build with `USE_MIMALLOC=0`. See [Memory allocator (mimalloc)](../user-guide/allocator.md).

**Single-binary SIMD dispatch**
On x86, bwa-mem3 ships one binary that contains compiled kernels for every supported SIMD tier (`sse41` / `sse42` / `avx` / `avx2` / `avx512bw`) and selects one in process at startup via `__builtin_cpu_supports`. There are no per-tier companion binaries. On ARM64 the binary contains a single NEON kernel TU. Replaces the prior multi-binary `execv` launcher (PR #83). See [Single-binary SIMD dispatch (x86)](../developer-guide/launcher.md).

**PGO**
Profile-Guided Optimization — a two-pass build where the first pass instruments the binary, a representative workload is run to collect profiles, and the second pass uses those profiles to guide inlining and branch layout. Activated via `make pgo-generate` then `make pgo-use`. See [PGO build](../performance/pgo.md).

**Primary alignment**
The alignment record for a read that represents the aligner's best placement. A read has exactly one primary alignment (or is reported as unmapped). All other alignments for the same read are marked supplementary (chimeric split read) or secondary (alternative mapping). See [Output: SAM/BAM, headers, tags](../user-guide/output.md).

**Proper-pair flag (0x2)**
SAM flag bit indicating that both mates of a pair are mapped in the expected orientation and insert-size range. In bwa-mem3, the `mem_sam_pe` function sets this flag; a correctness fix (PR #17) ensures it is propagated correctly under all conditions. See [Correctness fixes](../whats-different/correctness.md).

**SAM**
Sequence Alignment Map — a tab-delimited text format for read alignments. Each record contains mandatory fields (QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL) plus optional tags. See [Output: SAM/BAM, headers, tags](../user-guide/output.md).

**SIMD dispatch**
Runtime selection of the fastest available SIMD instruction set (SSE4.1, SSE4.2, AVX, AVX2, AVX-512BW, NEON) for hot alignment kernels. On x86 this is implemented in process by `src/simd_dispatch.cpp` via `__builtin_cpu_supports`; on ARM64 a single NEON tier covers every supported CPU. See [SIMD dispatch matrix](../performance/simd-dispatch.md).

**SMEM**
Super-Maximal Exact Match — a seed found by extending a read's position in the FM-index as far as possible in both directions. SMEMs form the initial seeds for chaining and extension in the BWA-MEM algorithm. See [Performance improvements](../whats-different/performance.md).

**Soft clip**
A CIGAR operation (`S`) indicating that bases at the read end were not part of the alignment, but are still present in the SEQ field. Soft clipping commonly appears at adapter-containing or low-quality read ends. See [Output: SAM/BAM, headers, tags](../user-guide/output.md).

**Supplementary alignment**
A SAM record (FLAG bit 0x800 set) representing a chimeric read split across two or more genomic loci. The segment with the longest aligned span is typically designated primary; remaining segments are supplementary. Hard clipping is used to avoid duplicating the SEQ field. See [Output: SAM/BAM, headers, tags](../user-guide/output.md).

---

**See also:**
[Citation](citation.md) ·
[License](license.md) ·
[Changelog](changelog.md) ·
[Output: SAM/BAM, headers, tags](../user-guide/output.md) ·
[What's Different — Overview](../whats-different/overview.md)
