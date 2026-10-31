# src/aln: bwa `aln`/`samse`/`sampe` port

This directory holds the code behind `bwa-mem3 aln`, a fused port of bwa's
short-read `aln` + `samse`/`sampe` pipeline whose output is byte-identical to
stock bwa 0.7.19.

## Provenance

| Item | Value |
| --- | --- |
| Upstream repository | https://github.com/lh3/bwa |
| Branch | `Apache2` |
| Commit | `a08c7c0` (a08c7c05517ab534b1618514507e95cfeaf8cdc7) |
| License | Apache License 2.0 (`LICENSE` in this directory); relicensing statement at https://github.com/lh3/bwa/issues/389 |

Code is copied only from the `Apache2` branch, never from `master` or a release
tag. Each file keeps its upstream copyright block.

## Files

| File | Origin |
| --- | --- |
| `bwtgap.h`, `bwtgap.cpp` | `bwtgap.h`/`bwtgap.c`: the backtracking search (`bwt_match_gap` and its priority stack) |
| `bwtaln.cpp` | the search half of `bwtaln.c`: `gap_init_opt`, `bwa_cal_maxdiff`, `bwt_cal_width`, `bwa_cal_sa_reg_gap` and the thread fan-out |
| `bwase.h`, `bwase.cpp` | `bwase.h`/`bwase.c`: the SAM half of samse (`bwa_aln2seq_core`, `bwa_sa2pos`, `bwa_cal_pac_pos_core`, `bwa_refine_gapped`, `bwa_print_sam1`), plus `bns_get_seq` and `bns_cnt_ambi` from `bntseq.c` |
| `bwape.h`, `bwape.cpp` | `bwape.c`: the pairing half of sampe (`bwa_init_pe_opt`, `infer_isize`, `pairing`, `bwa_cal_pac_pos_pe`, `bwa_sw_core`, `bwa_paired_sw`) |
| `aln_out.h`, `aln_out.cpp`, `aln_main.h`, `aln_main.cpp` | new code: the SAM/BAM sink and the `aln` subcommand driver (`bwa_aln` and `bwa_sai2sam_se` option handling, fused) |
| `aln_bwt.h`, `aln_bwt.cpp` | bwa's FM-index API (`bwt_occ`, `bwt_2occ`, `bwt_2occ4`, `bwt_sa`) re-implemented over bwa-mem3's own index (`FMI_search`); `bwt_match_exact_alt` is copied from upstream `bwt.c` |

## Modifications

Each task that changes or adds upstream-derived code appends a line here.

- `aln_bwt.cpp`: `bwt_match_exact_alt` is upstream's body verbatim, moved into
  `namespace bwa_aln`. The occ and SA functions are new code with bwa's
  semantics, computed from bwa-mem3's checkpointed occ table instead of bwa's
  interleaved 2-bit `.bwt` layout.
- `bwtaln.h`: upstream's types and constants unchanged, moved into
  `namespace bwa_aln` (the `bwt_t`/`bwtint_t` types come from `aln_bwt.h`). Dropped:
  `SAI_MAGIC`, the `extern "C"` function declarations and `bwa_bam_open`.
- `bwaseqio.cpp`/`bwaseqio.h`: `bwa_seq_open`, `bwa_seq_close`, `seq_reverse`,
  `bwa_trim_read`, `bwa_read_seq` and `bwa_free_read_seq` are upstream's bodies, in
  `namespace bwa_aln`. Dropped: BAM input (`bwa_bam_open`, `bwa_read_bam`, the `is_bam`
  branches and `bamlite.h`). Input is opened with `xzopen` (zlib `gzopen`, or
  `gzdopen` on stdin for `-`), as bwa does. `nst_nt4_table` is upstream's table, kept
  local because bwa-mem3's differs at `'-'` (5 upstream, 4 here). Parsing uses
  bwa-mem3's `kseq.h`, whose record handling matches upstream's; its differences are
  an OOM guard, a wider `kseq_read` return type and extra sequence-buffer slack.
- `bwtgap.cpp`/`bwtgap.h`: `gap_init_stack*`, `gap_destroy_stack`, `bwt_match_gap` and the
  stack helpers are upstream's bodies, in `namespace bwa_aln`. Dropped: the `extern "C"`
  block, the malloc wrappers and the `bwa_aln2seq` declaration (ported with the SAM
  code). Added to `bwtgap.h`: declarations for `gap_init_opt`, `bwa_cal_maxdiff` and
  `bwt_cal_width` (upstream declares them in `bwtaln.h`), plus `aln_search_batch` and
  `aln_batch_max_gapo`, which are new.
- `bwtaln.cpp`: `gap_init_opt`, `bwa_cal_maxdiff` and `bwt_cal_width` are upstream's
  bodies. `bwa_cal_sa_reg_gap` differs in three ways. (1) The read is complemented into
  a per-thread scratch buffer instead of in place, and `name`, `seq`, `rseq` and `qual`
  are no longer freed, so the SAM stage can use the same reads. (2) The batch-wide
  `max_len`, the `max_diff` derived from it and the `max_gapo` cap are computed once in
  `aln_search_batch` before the threads start (upstream repeats the same computation in
  every thread) and passed down; the per-read `max_diff` is recomputed per read as
  upstream does. (3) The thread count is a parameter rather than `opt->n_threads`, and
  `i % n_threads == tid` is applied for one thread too, so a count of 0 cannot divide by
  zero. Dropped: `bwa_aln_core` (the `.sai` writer and read loop), `bwa_open_reads`,
  `bwa_aln` (option parsing and usage) and the `HAVE_PTHREAD` conditionals (pthreads
  are always used). The code is built with `-ffp-contract=off` (Makefile) so
  `bwa_cal_maxdiff` rounds as stock's does.
- `bwase.cpp`/`bwase.h`: `bwa_aln2seq_core`, `bwa_aln2seq`, `bwa_approx_mapQ`, `bwa_sa2pos`,
  `bwa_cal_pac_pos_core`, `bwa_refine_gapped_core`, `bwa_cal_md1`, `bwa_correct_trimmed`,
  `pos_end`, `bwa_print_sam1` and `bwase_initialize` are upstream's bodies, in
  `namespace bwa_aln`, with these changes: `bwa_print_sam1` and the header writer append to
  a `kstring_t` instead of writing to stdout (`bwa_print_seq` likewise); `bwa_cal_pac_pos`
  takes the `bwt_t` as a parameter instead of reloading `.bwt` and `.sa` from disk, and
  `bwa_refine_gapped` always receives the packed reference (upstream's `NULL` branch that
  reads `.pac` is dropped); the `.sai` reading, `bwa_sai2sam_se_core` and the samse command
  line are not ported (alignments come from the search's `aln` arrays, and the driver is
  `aln_main.cpp`); and `bwa_print_sam1` / `bwa_refine_gapped_core` have C++ casts. Output
  goes to a sink that writes per batch, in read order.
- `bwase.cpp`: `bns_get_seq` and `bns_cnt_ambi` are upstream `bntseq.c`'s bodies, kept as
  file-local `aln_bns_get_seq` / `aln_bns_cnt_ambi`. bwa-mem3's own `bns_cnt_ambi` counts every
  overlapping N run (upstream stops at the first), which would change the `XN` tag.
- `aln_print_sam_hdr` follows bwa **0.7.19's** header output, not the Apache2 branch's
  `bwa_print_sam_hdr`. The branch prints `@SQ` before `@HD` and finds `@HD` anywhere in the
  `-r` text; the release prints `@HD` first (unless a line of the `-r` text starts `@HD<tab>`),
  then `@SQ`, then the `-r` text. The release is the parity target, so its order is used. The
  `@PG` line is bwa-mem3's own (same builder as `mem`); parity checks drop `@PG` lines.
- `bwa-mem3 aln` differs from running `bwa aln` then `bwa samse` in its command line: it takes
  `aln`'s options, `-f` and `-r` (samse's), and `--max-xa-hits` for samse's `-n`. `-b`, `-0`, `-1`
  and `-2` (BAM input) are not accepted. `-r` with a literal tab is rejected as stock does
  (bwa-mem3's `bwa_set_rg` alone would accept it). `--bam` and `--bam-level` are new.
- `srand48(bns->seed)` is called once per run, before the first batch, and `bwa_aln2seq_core` runs
  serially in read order, so the random choice among equally good hits is the same as stock's
  and does not depend on `-t`.
- `bwape.cpp`/`bwape.h`: `bwa_init_pe_opt`, `infer_isize`, `pairing`, `bwa_sw_core` and `bwa_paired_sw` are
  upstream's bodies in `namespace bwa_aln`; `bwa_cal_pac_pos_pe` is upstream's with these changes. Each read's
  alignments come from its own `aln` / `n_aln` (upstream reads them from the `.sai` and keeps a backup copy per
  read, which is not needed here); the BWT is a parameter; `infer_isize` is exported (the unit tests call it).
  Dropped: the dead `if (0 && ...)` hash-table branch and `g_hash` (so `-P` has no effect), the `.sai` reading
  and the sampe command line, the per-step timing lines on stderr, and the branch that loads the reference when
  the caller gives none (it is always the in-memory packed copy). `g_log_n` (bwase.cpp) is now shared with it. The
  mate alignment in `bwa_sw_core` calls bwa-mem3's `ksw_align` and `ksw_global`, which match stock's up to cosmetics;
  `aln_parity.sh` is the check. `infer_isize`'s `erfc`/`sqrt` loop, the pairing score's `log`/`erfc` and the
  `int += double` arithmetic in `bwa_paired_sw` are kept verbatim; the objects are built with `-ffp-contract=off`.
- Paired-end driver (`aln_run_pe`, new code). `bwa-mem3 aln [options] <prefix> <R1.fq> <R2.fq>` is `aln` on each
  file followed by `sampe`. Per batch of 0x40000 reads from each file: R1's reads are read, then R2's (stock
  order); each file is searched as its own batch, because stock runs `aln` on each file separately and the
  batch-wide gap-open cap depends on the longest read of the batch; pairing runs serially in read order (it
  draws from `drand48`), with the previous batch's insert-size estimate carried over; the mate Smith-Waterman, the
  gapped refinement and the output follow, each pair printed as read 1 then read 2. sampe's options are `-a`, `-s`,
  `-A`, `-c` and `-P` (as stock), `--max-occ` (sampe `-o`), `--max-xa-hits` (sampe `-n`, shared with the
  single-end XA limit) and `--max-xa-hits-discordant` (sampe `-N`). Where stock sampe takes the options saved in
  R2's `.sai` for both ends, one option set serves both files here. When the two files give different read
  counts in a batch, sampe's rule is kept (the count of R2) and a warning is written to stderr; if R2 has more
  reads than R1 the count of R1 is used, where sampe would read past R1's reads. A pair whose names differ after
  the `/1` and `/2` trim stops the run with an error after that pair is written, as sampe does.
- The batch-wide `max_gapo` cap does not change any output, so searching the two files of a pair separately or
  together gives the same hits and the same `-m` behaviour. A gap open is tried only from an entry with at least one
  difference still allowed, so its gap opens stay below `min(-o, differences allowed for the file's longest read)`,
  and a larger cap never lets a further entry in. The driver still searches each file as its own batch, as stock does.
- `bwase.cpp`: `aln_print_sam_hdr` is new code that follows the 0.7.19 release's header
  order; its line-start counting is a small helper (`count_line_start`).
- `aln_main.cpp`: new code. Inputs and `-f` are checked before the index loads, so a failed
  run leaves no output; help text lists each option with its stock default and source command.
