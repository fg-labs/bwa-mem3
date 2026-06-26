# Branch and worktree conventions

This page describes how the bwa-mem3 repository branches relate to upstream bwa-mem2, the policy for where PRs land, and the conventions for local worktrees when working on multiple branches simultaneously.

## Branch model

### `master` — upstream mirror

`master` tracks the upstream [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) `master` branch verbatim. No fork-carried changes are applied here. When upstream bwa-mem2 merges new commits, `master` is fast-forwarded to match.

`master` is the starting point for upstream rebase operations. It is never the target of fork PRs.

### `main` — fork integration branch

`main` carries all fork-carried commits on top of a rebased upstream baseline. This is the branch that:

- All new feature, fix, and improvement PRs target.
- All git tags (`v0.X.Y`) are placed on.
- Read the Docs `/latest/` follows.

When upstream bwa-mem2 makes significant changes, `master` is fast-forwarded and then `main` is rebased onto the new `master` tip. The rebase is verified by running the full test suite before the result is pushed.

### Feature and fix branches

All development work happens on short-lived branches that are merged into `main` via pull request. Branch name conventions:

| Prefix | Use |
|---|---|
| `feat/` | New features or capabilities |
| `fix/` | Bug fixes |
| `perf/` | Performance improvements |
| `test/` | Test additions or improvements |
| `docs/` | Documentation changes |
| `ci/` | CI / build system changes |
| `refactor/` | Code restructuring without behaviour change |

Branch names use kebab-case after the prefix: `fix/kswv-nrow-zero`, `perf/libsais-fm-index`, `test/regression-tests`.

## Upstream rebase cadence

`main` is rebased onto `master` (i.e., onto upstream bwa-mem2) periodically — not on every upstream commit, but when upstream merges a batch of changes worth incorporating. The process is:

1. Fast-forward `master` to the new upstream tip.
2. Rebase `main` onto `master`, resolving any conflicts.
3. Run `make && make test` to confirm the rebase is correct.
4. Push `master` and `main` to the `fg-labs` remote.

> **Warning — Do not merge upstream into main**
>
> Always rebase rather than merge when incorporating upstream changes. Merge commits
> obscure the fork-carried commit history and make the What's Different table harder
> to maintain.

## Worktrees for parallel branches

When working on multiple branches simultaneously, use git worktrees instead of stashing or switching branches. Each worktree is a sibling directory of the main clone.

### Creating a worktree for a PR branch

```bash
# Fetch the PR's head branch from the fg-labs remote
git fetch fg-labs <head-branch-name>

# Create a worktree with a local branch tracking the remote branch
git worktree add ../pr-<N> -b pr-<N> --track fg-labs/<head-branch-name>
```

The local branch name and directory name match the PR number (`pr-N`).

### Creating a worktree for a new issue branch

```bash
# Fetch the latest main from fg-labs
git fetch fg-labs main

# Create a new feature branch off fg-labs/main
git worktree add ../issue-<N> -b <prefix>/issue-<N>-<short-slug> fg-labs/main

# Unset the upstream so the branch is untracked until first push
git -C ../issue-<N> branch --unset-upstream
```

On first push, push to `fg-labs` so the head branch is in the same organisation as the PR base:

```bash
git push -u fg-labs HEAD
```

### Worktree naming conventions

| Directory name | Branch type |
|---|---|
| `main/` | Primary checkout; tracks `fg-labs/main` |
| `pr-<N>/` | PR review; local branch `pr-N` tracks `fg-labs/<head-branch>` |
| `issue-<N>/` | Issue work; local branch `<prefix>/issue-N-<slug>` |
| Descriptive name | Feature work not yet tied to a PR or issue |

### Listing and removing worktrees

```bash
# List all worktrees
git worktree list

# Remove a worktree after the PR is merged
git worktree remove ../pr-<N>
git branch -D pr-<N>

# Remove an issue worktree
git worktree remove ../issue-<N>
git branch -D <prefix>/issue-<N>-<slug>
```

> **Note — Worktree directories are siblings, not nested**
>
> All worktree directories sit next to the main clone at the same directory level,
> not inside it. This avoids confusing `git` commands that walk parent directories
> looking for `.git`.

## PR policy

- All PRs target `main`.
- PRs from fork contributors should be opened against `fg-labs/bwa-mem3 main`.
- Every PR that adds a fork-carried commit must update the `FG-MAIN-TABLE` in `docs/src/reference/pr-catalog.md` in the same PR. See [Contributing](contributing.md#the-fg-main-table-rule).
- Merge policy: squash-merge for single-commit changes; rebase-merge for multi-commit PRs with a clean commit history.

---

**See also:**
[Contributing](contributing.md) ·
[Release process](release.md) ·
[What's Different → Overview](../whats-different/overview.md) ·
[Building from source](building.md)
