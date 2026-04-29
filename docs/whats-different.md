# What's different from bwa-mem2

This document tracks the changes carried by `bwa-mem3` on top of upstream
[`bwa-mem2/bwa-mem2`](https://github.com/bwa-mem2/bwa-mem2)'s `master` branch.

## Carried on top of upstream

<!-- FG-MAIN-TABLE:start -->
| Commit    | Topic                                                           | Upstream status |
|-----------|-----------------------------------------------------------------|-----------------|
| `ae73227` | Apple Silicon / ARM64 NEON support (PR #288 work)               | [PR #288][pr288] open |
| `744a9e7` | `ci`: cross-platform build + dwgsim phiX end-to-end test        | fork-only       |
| `490502b` | `fix`: drop unused global `stat` that shadows libc              | fork-only       |
<!-- FG-MAIN-TABLE:end -->

### Additional fork-level changes

- **Vendored mimalloc allocator**: `ext/mimalloc` is pinned at `v3.3.0` and
  linked into every binary by default (`USE_MIMALLOC=1`). Linux uses
  `--whole-archive` static linkage; macOS uses dyld-interposed shared linkage.
  `USE_MIMALLOC=1` is the supported and recommended default on all
  platforms; `USE_MIMALLOC=0` is provided as a best-effort opt-out and is
  CI-gated on Linux x86 only. See `README.md` → "Memory allocator" for
  details.

- **`--supp-rep-hard-cap INT`** (opt-in, default disabled): forces MAPQ=0 on
  supplementary alignments whose chain contains a seed with `>=INT` genome
  occurrences. Addresses the long-standing bwa/bwa-mem2 issue where a supp
  fragment that maps to many places standalone (e.g. a short read in a CCATCC
  repeat) inherits a high MAPQ from its primary because the supp's competing
  repetitive chains get filtered out during the full-read pipeline and
  therefore never contribute to its `sub`/`sub_n`. See
  [upstream #260](https://github.com/bwa-mem2/bwa-mem2/issues/260) for the
  reporter case. Primary MAPQ is unaffected; default output is byte-identical
  to stock bwa-mem2. Typical values are 5–20 (lower = more aggressive); the
  upstream #260 repro drops from MAPQ=60 to MAPQ=0 at `--supp-rep-hard-cap 18`.

[pr288]: https://github.com/bwa-mem2/bwa-mem2/pull/288

## Version stamping

`PACKAGE_VERSION` (the value reported by `bwa-mem3 version` and written to
the `@PG VN:` SAM header field) is generated at build time by the Makefile
from `git describe --tags --dirty`, e.g. `v2.3-30-g61813ef` for a tree 30
commits past upstream tag `v2.3` at commit `61813ef`.

- No manual bumping required: cut a fresh release by tagging the commit
  (`git tag -a vX.Y-fg-labs.N -m …`) and the next build picks it up.
- Builds where `git describe --tags` fails (source-tarball extractions, or
  shallow clones / checkouts with no tag reachable from `HEAD` — including
  CI's default `actions/checkout` fetch-depth of 1) fall back to the static
  `FG_LABS_VERSION_FALLBACK` in `Makefile`. Bump that when cutting a
  release that will be consumed as a tarball, or in CI artifacts.
- `src/version.h` is generated and `.gitignore`d; `make clean` removes it.

## Branching and update policy

- `master` tracks upstream unchanged.
- `main` is `upstream/master` plus the commits above. Rebased onto upstream roughly quarterly, or sooner when an upstream release we care about lands.
- Contributions go via PR targeting `main`. CI and CodeRabbit gate merges.
- Any PR that adds or removes a fork-carried commit must update the table above in the same PR.

## Consuming

Clone this repo and check out `main`:

```bash
git clone https://github.com/fg-labs/bwa-mem3.git
cd bwa-mem3
git checkout main
```

Or vendor the branch into a downstream repo by pinning to a specific
commit (not the branch tip) so your build is reproducible.

## Relationship to upstream

We submit the generally-useful fixes and features carried here as PRs against [`bwa-mem2/bwa-mem2`](https://github.com/bwa-mem2/bwa-mem2) when the upstream maintainers are actively merging; while they are not, fixes land here first and we drop them from `main` once they appear upstream.
