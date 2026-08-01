# Contributing

This page covers the mechanics of submitting changes to bwa-mem3: commit conventions, PR workflow, CI requirements, and the rule for keeping the fork-lineage table current.

## Before you start

1. Check the [open issues](https://github.com/fg-labs/bwa-mem3/issues) and [existing PRs](https://github.com/fg-labs/bwa-mem3/pulls) to avoid duplicate work.
2. For substantial changes, open an issue first to discuss scope and approach.
3. Fork or branch from `fg-labs/bwa-mem3 main`. See [Branch and worktree conventions](branches.md) for the branching model.

## Commit message conventions

bwa-mem3 follows [Conventional Commits](https://www.conventionalcommits.org/) (`v1.0.0`). Every commit message must start with a type prefix:

| Prefix | Use |
|---|---|
| `feat:` | New feature or capability |
| `fix:` | Bug fix |
| `perf:` | Performance improvement |
| `test:` | Test additions or changes |
| `docs:` | Documentation only |
| `ci:` | CI / build-system changes |
| `refactor:` | Restructuring without behaviour change |
| `chore:` | Maintenance (dependency bumps, version pins) |

The subject line is lowercase after the prefix, imperative mood, no trailing period. Keep it under 72 characters. Body lines wrap at 100 characters.

**Good:**

```text
fix: kswv nrow==0 batch skips rowMax store when i==0

Exercises the all-len1==0 path across SSE4.1, AVX2, AVX-512BW, and ARM NEON.
Without the `if (i > 0)` guard, the store writes SIMD_WIDTH* bytes before the
allocation.

Closes #38.
```

**Not acceptable:**

```text
Fixed stuff
Updated kswv
WIP
```

### Flagging breaking changes

Releases are generated automatically by [release-please](https://github.com/googleapis/release-please) from the conventional-commit history (see [Release process](release.md)). A change that breaks a downstream consumer — see the [semver policy](release.md#semver-policy) for what qualifies — **must** be marked so the breaking-change notice lands in `CHANGELOG.md` and the GitHub release body. Describing the break only in the commit body is not enough: release-please does not read prose, so an unmarked break is silently filed under its plain type (e.g. a `perf:` commit lands under "Performance" with no warning) and downstream users never see it.

Mark a break either way:

- Append a `!` after the type/scope: `perf(index)!: stop building .0123 by default`, **and/or**
- Add a `BREAKING CHANGE:` footer (the colon is required) describing the break and the migration path.

```text
perf(index)!: stop building the unpacked .0123 reference by default

BREAKING CHANGE: `bwa-mem3 index` no longer writes `.0123`; `mem` reconstructs
reference bases from `.pac` on demand. External tools that read `.0123` directly
(e.g. sharing an index with bwa-mem2) must re-run with `index --emit-unpacked-ref`.
```

release-please then emits a `⚠ BREAKING CHANGES` section automatically. While the project is pre-1.0 this still bumps only the minor version (`bump-minor-pre-major` is set), matching the [semver policy](release.md#semver-policy)'s rule that a pre-1.0 break is allowed when it is called out clearly — the footer is how it gets called out.

> **Note — Squash merges**
>
> The project [squash-merges single-commit PRs and rebase-merges multi-commit
> PRs with a clean history](branches.md). For a **squash-merged** PR,
> release-please reads only the squash-merge commit subject (which defaults to
> the PR title), **not** the individual commit bodies — so put the `!` /
> `BREAKING CHANGE:` marker on the PR title or the squash commit message itself,
> because a marker buried in a sub-commit body is discarded at squash time. A
> **rebase-merged** PR keeps its individual commits, so a marker in any sub-commit
> message is preserved and detected; even so, prefer putting it on the PR title
> as well so the convention is uniform regardless of how the PR is merged.

## Pull request workflow

1. Push your branch to `fg-labs/bwa-mem3` (or your fork) and open a PR targeting `fg-labs/bwa-mem3 main`.
2. The PR description should explain the motivation, summarise the change, and note any benchmarks or test results.
3. All CI jobs must pass before merge. See [CI matrix](#ci-matrix) below.
4. [CodeRabbitAI](https://www.coderabbit.ai) reviews every PR automatically. Address all comments, including inline suggestions, summary comments, and nitpicks. Do not dismiss comments without a reply explaining why the suggestion was not adopted.
5. A project maintainer will review and merge once CI is green and all comments are resolved.

> **Note — Draft PRs first**
>
> Open PRs as drafts while CI is running or while you are actively revising. Convert
> to ready-for-review only when the branch is stable, CI is green, and you have
> self-reviewed the diff.

## The FG-MAIN-TABLE rule

Every PR that introduces a new fork-carried commit — a commit that is on `main` but not on `master` (the upstream bwa-mem2 mirror) — **must** update the `FG-MAIN-TABLE` block in `docs/src/reference/pr-catalog.md` in the same PR.

The table records each fork-carried change, its bwa-mem3 PR number, the corresponding upstream bwa-mem2 PR or issue (if any), and its upstream status. Keeping this table current is the primary mechanism by which the project maintains transparency about its relationship to upstream.

> **Warning — Do not skip the table update**
>
> A PR that adds a fork-carried commit but omits the table update will be sent back
> for revision. The table is reviewed as part of the standard PR checklist.

### What counts as a fork-carried commit

A commit is fork-carried if:

- It adds new behaviour, fixes a bug, or changes build infrastructure in a way that diverges from upstream bwa-mem2 `master`.
- It is present on `fg-labs/bwa-mem3 main` but not (yet) merged upstream.

Pure documentation commits, CI-only changes, and upstream-rebase bookkeeping commits do not need a table entry.

## CI matrix

CI runs on every PR and on push to `main`. The matrix covers:

| Row | Architecture | ISA | Platform |
|---|---|---|---|
| `sse41` | x86_64 | SSE4.1 | Ubuntu |
| `avx2` | x86_64 | AVX2 | Ubuntu (canonical) |
| `avx512bw` | x86_64 | AVX-512BW | Ubuntu |
| `arm64-linux` | aarch64 | NEON | Ubuntu ARM |
| `arm64-macos` | arm64 | NEON | macOS |

The canonical row (`avx2`) is the only one that runs regression tests (shell scripts in `test/regression/`). Unit tests run on every row. Integration tests run on the four widened canonical rows (SSE4.1, AVX2, ARM64 Linux, macOS ARM).

Alongside the matrix, three source-only lint jobs run without a toolchain and report in seconds: the NDEBUG gate lint, the opt-in debug macro list, and [shell lint](#shell-scripts).

A PR must pass all rows before merge.

## Code style

- C++14, `gnu++14` dialect.
- Match the style of the surrounding code. The codebase inherits the upstream bwa-mem2 style, which is C-ish C++ with minimal STL use in hot paths.
- For new test code, follow the [doctest patterns](regression-tests.md) documented in the test framework.
- New SIMD code must include `src/simd_compat.h` rather than platform-specific headers directly. See [SIMD dispatch architecture](simd-dispatch.md).

### Shell scripts

Shell is a load-bearing part of the test surface — `test/regression/` alone is dozens of scripts — so tracked `*.sh` files are gated in CI on two tools:

| Tool | Checks | Version |
|---|---|---|
| [shellcheck](https://www.shellcheck.net/) | correctness, at default severity | 0.11.0 |
| [shfmt](https://github.com/mvdan/sh) | layout, `-i 4 -ci -bn -sr` | v3.13.1 |

Both versions are pinned, in `test/regression/shell_lint.sh` and in the `shell-lint` job in `.github/workflows/ci.yml`. Bump them together; a new upstream release that adds a check would otherwise fail a PR that never touched shell.

Run the gate, or let it fix the formatting for you:

```bash
make shell-lint    # check (this is what CI runs)
make shell-fix     # apply shfmt, then report anything shellcheck still wants
```

`shell-fix` only fixes layout. A shellcheck finding needs a human: either correct the code, or add a scoped disable *with a comment saying why the diagnostic does not apply*.

```bash
# shellcheck disable=SC2016  # $1/$NF here are awk fields, not shell vars
```

Both targets skip with a visible notice if shellcheck and shfmt are not installed, so `make test` works without them. Install with `brew install shellcheck shfmt` or your distribution's equivalent.

CI does *not* get that leniency: the `shell-lint` job sets `SHELL_LINT_REQUIRED=1`, which turns a missing tool into a hard error rather than a skip — otherwise a broken install step would leave the job permanently, silently green. To reproduce that locally, run `SHELL_LINT_REQUIRED=1 make shell-lint`.

Two things worth knowing before you write a new script:

- New tracked `*.sh` files are picked up automatically — discovery is `git ls-files`, so there is no list to add yourself to.
- A comment line whose first word is `shellcheck` is parsed as a directive, not prose. Reword it, or the file will fail to lint with SC1072/SC1073.

## Adding a test for your change

- **Bug fix** → add a unit test or integration test that fails without the fix and passes with it.
- **New feature** → add unit tests for the core logic and, if the feature is end-to-end testable with a shell invocation, a regression test in `test/regression/`.
- **Performance change** → run the benchmark harness (`bench/`) to confirm the improvement and include median wall-clock numbers in the PR description.

See [Regression test framework](regression-tests.md) for the full guide on where to add tests and how to organise them.

---

**See also:**
[Branch and worktree conventions](branches.md) ·
[Regression test framework](regression-tests.md) ·
[Release process](release.md) ·
[What's Different → Overview](../whats-different/overview.md) ·
[Building from source](building.md)
