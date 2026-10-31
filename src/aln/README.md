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
