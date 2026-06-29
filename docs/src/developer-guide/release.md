# Release process

bwa-mem3 follows [semantic versioning](https://semver.org). Releases are automated with [release-please](https://github.com/googleapis/release-please): every push to `main` updates a standing "release PR" that bumps the version and regenerates the changelog from the [Conventional Commits](contributing.md#commit-message-conventions) history; merging that PR tags the release and publishes a GitHub Release. The version string is derived from `version.txt` and embedded in every binary at compile time.

## Version stamping

`version.txt` at the repo root is the **single source of truth** for the version. release-please rewrites it (and `.release-please-manifest.json`) in the release PR; nothing else should edit it by hand.

The Makefile computes the build's version string from it at parse time:

```makefile
# version.txt is the single source of truth; scripts/version.sh reads it
# and appends an informational git-describe-style dev suffix.
VERSION_STRING := $(shell scripts/version.sh)
```

[`scripts/version.sh`](https://github.com/fg-labs/bwa-mem3/blob/main/scripts/version.sh) reads the base version from `version.txt` and, when git is available, appends a dev suffix that surfaces how far the working tree is from the matching tag:

| Working-tree state | Example string |
|---|---|
| Source tarball / shallow clone (no git) | `0.4.0` |
| Clean and exactly at tag `v0.4.0` | `0.4.0` |
| Clean, not at the tag | `0.4.0-3f7ab2e` |
| Uncommitted changes at the tag | `0.4.0-dirty` |
| Uncommitted changes, not at the tag | `0.4.0-3f7ab2e-dirty` |

"At the tag" means HEAD is precisely the commit pointed to by tag `v<base>` (or `<base>`). Manifest drift — HEAD at a tag that disagrees with `version.txt` — is treated as "not at tag", so the `-<sha>` suffix surfaces the drift visibly rather than silently printing a wrong bare version.

The string is written into `src/version.h` by the `src/version.h: FORCE` rule, which runs on every `make` invocation but only touches the file when the string changes. This minimises unnecessary recompilation of `src/main.o`.

`PACKAGE_VERSION` from `src/version.h` appears in:

- `bwa-mem3 version` output (stdout).
- The `@PG VN:` field in every SAM/BAM file produced by `bwa-mem3 mem`.

Note the version string has **no leading `v`** (it mirrors `version.txt`, e.g. `0.4.0`), even though the git tags do (`v0.4.0`).

### Verifying the version

```bash
./bwa-mem3 version
# Example output on a release tag:
# 0.4.0
# mimalloc 3.3.0        ← if USE_MIMALLOC=1
```

On a commit past the tag the string carries the short SHA suffix:

```text
0.4.0-3f7ab2e
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

Run through this list on the candidate `main` commit **before merging the
release PR** (merging it is what tags and publishes the release — see
[Cutting a release](#cutting-a-release)). **Every item must pass.**

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
- [ ] The release notes are generated automatically by release-please
      from the conventional-commit history, so the real check is
      upstream: every user-visible PR in the release window has a
      correct conventional-commit type, and any breaking change carries
      a `!` / `BREAKING CHANGE:` marker so it lands in the `⚠ BREAKING
      CHANGES` section (see
      [Flagging breaking changes](contributing.md#flagging-breaking-changes)).
      `NEWS.md` is **not** updated — it is frozen at 0.2.0; 0.3.0 and
      later live only in `CHANGELOG.md` and the GitHub Releases page.
- [ ] `docs/src/reference/pr-catalog.md` `FG-MAIN-TABLE` block has a row for
      every fork-carried PR landed since the previous tag, with its upstream
      disposition (see [Contributing](contributing.md#the-fg-main-table-rule)).
- [ ] `docs/src/reference/changelog.md` and `docs/src/cli/version.md`
      examples reference the new release string.
- [ ] Spot-check the bwa-mem3-bench reference numbers in
      `docs/src/performance/overview.md` against the bench's
      `regression.md` for the tagging SHA.

## Cutting a release

Releases are not tagged by hand. The [`.github/workflows/release.yml`](https://github.com/fg-labs/bwa-mem3/blob/main/.github/workflows/release.yml) workflow runs release-please on every push to `main` and a tarball job on each published release.

1. **release-please maintains a standing release PR.** After PRs land on
   `main`, release-please opens (or updates) a PR titled
   `chore(main): release X.Y.Z`. It computes the next version from the
   conventional-commit history since the last tag, bumps `version.txt`
   and `.release-please-manifest.json`, and prepends the generated
   section to `CHANGELOG.md`. The bump level is driven entirely by the
   commit types: `feat:` → minor, `fix:`/`perf:`/etc. → patch, and a
   `!` / `BREAKING CHANGE:` marker forces the breaking bump (a minor
   pre-1.0, since `bump-minor-pre-major` is set). If the proposed version
   is wrong, the fix is upstream — correct the offending commit's type or
   add a breaking marker (see
   [Flagging breaking changes](contributing.md#flagging-breaking-changes)),
   not the release PR.

2. **Review the release PR.** Confirm the proposed version matches the
   change set per the [semver policy](#semver-policy), and that the
   generated `CHANGELOG.md` reads correctly — in particular that any
   breaking change appears under `⚠ BREAKING CHANGES`. Run the
   [release-readiness checklist](#release-readiness-checklist) against the
   PR's base commit.

3. **Merge the release PR.** This is the action that ships the release.
   release-please creates the `vX.Y.Z` tag and a GitHub Release whose body
   is the generated changelog section. Read the Docs activates a versioned
   build at `/vX.Y.Z/` automatically once the tag appears.

4. **The tarball job runs automatically** once the release is created. It
   checks out the tag with submodules, verifies `version.txt` matches the
   tag, builds the vendored `Source_code_including_submodules.tar.gz`
   (all submodules bundled, no `.git/`), smoke-tests that it compiles and
   reports the right version, uploads the asset plus its `.sha256`, and
   appends a "For packagers" block (with the asset URL and sha256) to the
   release body. Bioconda recipes pin against this asset.

> **Note — Manual rebuild of a tarball asset**
>
> The tarball job can be re-run for an existing tag via the workflow's
> `workflow_dispatch` input (e.g. to repair a missing asset), but `v0.2.0`
> is rejected — its asset is pinned by sha256 in an open bioconda PR and
> must not change.
>
> <hr>
>
> **Note — Tarball builds and the version string**
>
> A source tarball has no git history, so `scripts/version.sh` cannot append a
> dev suffix — but it still reads the git-tracked `version.txt`, so a `make`
> from the tarball prints the bare base version (e.g. `0.4.0`) with no further
> action needed.

## Post-release verification

After the release PR is merged and the GitHub release is published:

- Wait ~5 minutes for Read the Docs to build the new version, then
  open `https://bwa-mem3.readthedocs.io/en/v0.X.Y/` and confirm:
  - The version selector lists `v0.X.Y`.
  - The home page renders with no missing-page errors.
  - `developer-guide/launcher.md`, `performance/overview.md`, and
    `methylation/tags.md` all render with their mermaid diagrams and
    tables intact (these are the most diagram-heavy pages).
- Pull the tag in a clean clone and verify `bwa-mem3 version`
  reports the bare version string (no `-<sha>` dev suffix):

  ```bash
  git clone -b v0.X.Y --depth 1 https://github.com/fg-labs/bwa-mem3.git
  cd bwa-mem3 && make
  ./bwa-mem3 version | head -1
  # expect: 0.X.Y   (no leading 'v'; mirrors version.txt)
  ```

- If the docs build failed on RTD or the version string is wrong, do
  not delete or move the tag. Tags are immutable in practice — let
  release-please open the next release PR with the fix (a follow-up
  `0.X.(Y+1)` patch) instead.

## Branch and tag conventions

- All release tags are on the `main` branch, which carries both upstream bwa-mem2 commits and fork-carried changes. See [Branch and worktree conventions](branches.md) for the full branching model.
- Tags are prefixed with `v`: `v0.1.0`, `v0.2.0`, etc.
- Pre-release tags use a `-pre` suffix: `v0.1.0-pre`.
- Patch releases increment the third component: `v0.1.1`.

## What's Different table update

When a release bundles new fork-carried commits that were not previously documented, update the `FG-MAIN-TABLE` in `docs/src/reference/pr-catalog.md` in the same PR before tagging. See [Contributing](contributing.md#the-fg-main-table-rule) for the rule.

---

**See also:**
[Branch and worktree conventions](branches.md) ·
[What's Different → Overview](../whats-different/overview.md) ·
[Reference → Changelog](../reference/changelog.md) ·
[Building from source](building.md)
