# bwa-mem2 performance porting harness

Minimal reproducible benchmark for the performance porting effort tracked in `PERF_PORTING_PLAN.md`. Each PR in that plan uses this harness to (a) gate byte-identical SAM output via a golden md5 and (b) measure wall-clock + max-RSS deltas.

## Quick start

```sh
cp bench/config.env.example bench/config.env
$EDITOR bench/config.env                 # point at your index + reads + binary
bench/run.sh baseline                    # runs N trials, appends rows to bench/results.csv
bench/run.sh candidate                   # runs N trials on the candidate binary
bench/compare.sh baseline candidate      # prints wall-clock / RSS / md5 deltas
```

## What it measures

For each run:
- **Golden md5** — single-threaded, `@PG`-stripped SAM md5. Stable across builds when output is byte-identical.
- **Wall clock** — N multi-threaded trials, median reported.
- **Max RSS** — from `/usr/bin/time` (`-l` on macOS, `-v` on Linux). Median reported.

Results append to `bench/results.csv` as rows:
```
tag,host,arch,binary,threads,trial,wall_s,max_rss_kb,md5
```

The md5 is populated only on the single-thread row (marked `trial=golden`).

## Conventions

- **Golden diff** uses `-t 1`. Multi-threaded output ordering is not stable across runs — only wall-clock is comparable.
- `bench/config.env` is `.gitignore`'d. `config.env.example` is the template.
- Inputs must stay fixed across PRs for comparability. The default config points at `/Volumes/scratch-00001/bwa-mem3-bench/data/smoke/1M` (1M read pairs) + hg38.
- For representative multi-arch ranking, rerun the harness on an AVX-512 host (see `CLAUDE.md` AWS section).

## Files

- `run.sh` — main driver. Runs 1 golden + N perf trials against one binary.
- `compare.sh` — diff two tags from `results.csv`.
- `normalize_sam.sh` — `grep -v '^@PG' | md5sum` helper. Used by `run.sh`.
- `config.env.example` — path template.
- `results.csv` — append-only measurement log.
