# Release process

bwa-mem3 follows [semantic versioning](https://semver.org). Releases are driven by git tags. The version string is derived automatically from `git describe` and embedded in every binary at compile time.

## Version stamping

The Makefile computes the version string at parse time:

```makefile
FG_LABS_VERSION_FALLBACK := 0.1.0-pre
VERSION_STRING := $(shell git describe --tags --dirty 2>/dev/null || echo $(FG_LABS_VERSION_FALLBACK))
```

`git describe` produces a string such as `v0.1.0` (on a tag), `v0.1.0-3-gabcdef1` (three commits past the tag), or `v0.1.0-dirty` (uncommitted changes). If `git describe` fails — for example in a source tarball or a shallow clone without tag history — the build falls back to `FG_LABS_VERSION_FALLBACK`.

The string is written into `src/version.h` by the `src/version.h: FORCE` rule, which runs on every `make` invocation but only touches the file when the string changes. This minimises unnecessary recompilation of `src/main.o`.

`PACKAGE_VERSION` from `src/version.h` appears in:

- `bwa-mem3 version` output (stdout).
- The `@PG VN:` field in every SAM/BAM file produced by `bwa-mem3 mem`.

### Verifying the version

```bash
./bwa-mem3 version
# Example output on a tagged commit:
# v0.1.0
# mimalloc 3.3.0        ← if USE_MIMALLOC=1
```

On an untagged commit the string includes the commit distance and short SHA:

```text
v0.1.0-12-g3f7ab2e
```

## Semver policy

bwa-mem3 follows semver, interpreted for an alignment tool as follows.

**MAJOR (`X.0.0`)** — bump when the change would break a downstream
consumer that pinned the previous version without checking release
notes. Concretely:

- An on-disk index file format change (a re-index is required to use
  the new version).
- Removal or rename of a CLI flag or subcommand.
- A SAM/BAM tag is removed, renamed, or its type/value space changes
  incompatibly (a column-fixed downstream parser would break). Adding
  a new tag is **not** a major change.
- A change to the resolved primary alignment that is intentional and
  affects more than a negligible fraction of reads (e.g. a MAPQ
  recalibration applied unconditionally). Concordance regressions
  attributable to bug fixes are not major changes — call them out in
  the release notes under "Correctness" instead.
- Dropping support for a previously supported host class
  (e.g. raising the build's compiled-in `BASELINE_ARCH` floor in a way
  that excludes hosts the previous release ran on).

**MINOR (`0.X.0`)** — bump for any user-visible new functionality that
does not break consumers pinned to the previous minor. Examples:

- A new CLI flag or subcommand.
- A new SAM aux tag emitted on output (e.g. `HN:i` in 0.1.0, the
  Bismark `XR:Z` / `XG:Z` / `XM:Z` set in 0.2.0).
- A new operational feature (e.g. `bwa-mem3 shm`, in-process SIMD
  dispatch).
- A user-facing default change that is documented in release notes but
  does not require any consumer action (e.g. `BASELINE_ARCH=avx2` as
  the build default).
- New performance characteristics that change wall-time meaningfully.

**PATCH (`0.0.X`)** — bump for bug fixes, doc-only changes, build
fixes, and internal refactors that have no user-visible behavioral
delta. Pre-existing-bug fixes that incidentally shift output for a
small fraction of reads are patch-level when called out in the
release notes; widespread output shifts (>0.1% of reads on a typical
WGS bench cell) deserve MINOR or MAJOR depending on the source.

While the project is pre-1.0, the leading `0.` is treated literally —
`0.2.0` may make breaking changes vs `0.1.0` if called out clearly in
the release notes. After 1.0, MAJOR bumps are reserved for genuinely
breaking changes.

## Release-readiness checklist

Run through this list on the commit you intend to tag, before any
git-tag command. **Every item must pass.**

### Build and test

- [ ] `make clean && make` succeeds at the default `BASELINE_ARCH`
      (`avx2`) on a Linux x86_64 host.
- [ ] `make clean && make BASELINE_ARCH=sse41` succeeds on the same
      host — confirms the portability floor still compiles.
- [ ] `make clean && make` succeeds on an arm64 host (Apple Silicon
      or aarch64 Linux).
- [ ] `make test` passes on both x86_64 and arm64.
- [ ] `test/regression/all_tiers_parity.sh` produces byte-identical
      SAM across `BWAMEM3_FORCE_TIER=sse41 → sse42 → avx → avx2 → avx512bw`
      on an AVX-512BW host. Failures here indicate a per-tier kernel
      or dispatcher-wiring regression — fix before tagging.

### Bench

- [ ] [bwa-mem3-bench](https://github.com/fg-labs/bwa-mem3-bench) run
      submitted on the candidate SHA via
      `bwa_mem3_bench.cli submit --fg-labs-sha <sha>` (or the local
      smoke path for a fast sanity check).
- [ ] `bench regression --prev <previous-tag-sha>` reports gate `PASS`
      — concordance ≥ 99.999% on every `vs-baseline.json` cell except
      methylation (which is expected to drift vs the bwameth
      baseline; see the methylation carve-out below) and no cell
      labeled `REGRESSION`.
- [ ] Methylation cells reviewed for expected-drift consistency: the
      `meth-twist-emseq-5M` concordance vs the bwameth baseline should
      sit at ~98.9% post-PR-#90, with the per-class breakdown matching
      the entry in
      [`bwa-mem3-bench/docs/expected-divergences.yaml`](https://github.com/fg-labs/bwa-mem3-bench/blob/main/docs/expected-divergences.yaml)
      (or the entry added in this release — the file is in the bench
      repo, not in this repo).

### Docs

- [ ] `make docs` builds cleanly with no mdbook warnings.
- [ ] `NEWS.md` has a top-section entry for the new version with
      Operational / packaging, Correctness, Performance, and
      Methylation subheadings as applicable; every user-visible PR in
      the release window is listed with its number.
- [ ] `docs/src/whats-different/overview.md` `FG-MAIN-TABLE` is
      regenerated to cover the new PRs (see
      [Contributing](contributing.md) for the regeneration command).
- [ ] `docs/src/whats-different/upstream-prs.md` has rows for every
      user-visible PR landed since the previous tag.
- [ ] `docs/src/reference/changelog.md` and `docs/src/cli/version.md`
      examples reference the new release string.
- [ ] Spot-check the bwa-mem3-bench reference numbers in
      `docs/src/performance/overview.md` against the bench's
      `regression.md` for the tagging SHA.

## Tagging the release

Run these in order; each command depends on the previous one.

1. Pre-flight (confirms the readiness checklist):

   ```bash
   make clean && make
   make test
   make docs
   ```

2. Confirm `NEWS.md` is current. The top entry header line must match
   the tag you are about to create (e.g. `Release 0.2.0 (YYYY-MM-DD)`).

3. Tag the release commit. Prefer a signed tag (`-s`); fall back to an
   annotated tag (`-a`) only when signing is unavailable:

   ```bash
   git tag -s v0.X.Y -m "Release v0.X.Y"
   ```

4. Push the tag to the `fg-labs` remote:

   ```bash
   git push fg-labs v0.X.Y
   ```

   Read the Docs activates a versioned build at `/v0.X.Y/` automatically
   when the tag appears on the remote.

5. Create a GitHub release from the tag via `gh`. The body should be
   the matching `NEWS.md` section, no preamble:

   ```bash
   gh release create v0.X.Y --repo fg-labs/bwa-mem3 \
     --title "bwa-mem3 v0.X.Y" \
     --notes-file <(awk '/^Release / {p = ($0 ~ /^Release 0\.X\.Y /)} p' NEWS.md)
   ```

   Substitute `0.X.Y` with the exact version literal you are tagging,
   including any pre-release suffix — e.g. for `v0.2.0-pre` use
   `/^Release 0\.2\.0-pre /`. The trailing space in the inner pattern
   anchors the match to a complete version token so that `0.2.0` does
   not also match `Release 0.2.0-pre (...)`. The awk script prints
   lines while `p` is true: it flips on at the matching `Release` line
   and back off at the next `Release` line, which gives a clean
   section without needing a trailing `sed '$d'`.

> **Note — Tarball builds**
>
> Source tarballs created by GitHub (or `git archive`) do not include git history,
> so `git describe` fails and the version falls back to `FG_LABS_VERSION_FALLBACK`.
> For reproducible tarball builds, set `VERSION_STRING` explicitly on the command line:
> `make VERSION_STRING=v0.X.Y`.

## Post-release verification

After the tag is pushed and the GitHub release is published:

- Wait ~5 minutes for Read the Docs to build the new version, then
  open `https://bwa-mem3.readthedocs.io/en/v0.X.Y/` and confirm:
  - The version selector lists `v0.X.Y`.
  - The home page renders with no missing-page errors.
  - `developer-guide/launcher.md`, `performance/overview.md`, and
    `methylation/tags.md` all render with their mermaid diagrams and
    tables intact (these are the most diagram-heavy pages).
- Pull the tag in a clean clone and verify `bwa-mem3 version`
  reports the bare tag string (no `-N-gSHA` distance suffix):

  ```bash
  git clone -b v0.X.Y --depth 1 https://github.com/fg-labs/bwa-mem3.git
  cd bwa-mem3 && make
  ./bwa-mem3 version | head -1
  # expect: v0.X.Y
  ```

- If the docs build failed on RTD or the version string is wrong, do
  not delete or move the tag. Tags are immutable in practice — open a
  follow-up `v0.X.(Y+1)` patch release with the fix instead.

## Branch and tag conventions

- All release tags are on the `main` branch, which carries both upstream bwa-mem2 commits and fork-carried changes. See [Branch and worktree conventions](branches.md) for the full branching model.
- Tags are prefixed with `v`: `v0.1.0`, `v0.2.0`, etc.
- Pre-release tags use a `-pre` suffix: `v0.1.0-pre`.
- Patch releases increment the third component: `v0.1.1`.

## What's Different table update

When a release bundles new fork-carried commits that were not previously documented, update the `FG-MAIN-TABLE` in `docs/src/whats-different/overview.md` in the same PR before tagging. See [Contributing](contributing.md) for the rule.

---

**See also:**
[Branch and worktree conventions](branches.md) ·
[What's Different → Overview](../whats-different/overview.md) ·
[Reference → Changelog](../reference/changelog.md) ·
[Building from source](building.md)
