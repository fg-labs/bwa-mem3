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

A PR must pass all rows before merge.

## Code style

- C++14, `gnu++14` dialect.
- Match the style of the surrounding code. The codebase inherits the upstream bwa-mem2 style, which is C-ish C++ with minimal STL use in hot paths.
- For new test code, follow the [doctest patterns](regression-tests.md) documented in the test framework.
- New SIMD code must include `src/simd_compat.h` rather than platform-specific headers directly. See [SIMD dispatch architecture](simd-dispatch.md).

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
