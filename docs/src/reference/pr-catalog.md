# PR catalog

> **Retired.** This page used to hold a hand-maintained, per-PR table
> (`FG-MAIN-TABLE`) of every fork-carried change. It was enforced only by a
> review checklist, not by CI, so it drifted far behind `main` and became
> misleading — an authoritative-looking record that omitted most of the fork.
> It has been retired rather than backfilled by hand.

## Where the information lives now

- **What changed on top of upstream** — the git history is the source of truth:

  ```sh
  # every fork-carried, non-merge commit on main but not on the upstream mirror
  git log --reverse --no-merges master..main
  ```

  and the [merged pull requests](https://github.com/fg-labs/bwa-mem3/pulls?q=is%3Apr+is%3Amerged+base%3Amain)
  on GitHub, filterable by label, author, and date.

- **Which changes move alignment output vs. upstream** — the part that actually
  needs a curated, reviewed record — lives in the
  [`FG-DIVERGENCE-CATALOG`](../whats-different/equivalence.md), which is
  **auto-generated** from `bwa-mem3-bench`'s `docs/expected-divergences.yaml`
  and so cannot silently drift.

- **The *why* behind each class of change** — the narrative *What's Different*
  pages ([correctness](../whats-different/correctness.md),
  [performance](../whats-different/performance.md),
  [features](../whats-different/features.md),
  [architecture support](../whats-different/arch-support.md),
  [build & infrastructure](../whats-different/build-infra.md)) explain the
  rationale without trying to enumerate every PR.

If a curated per-PR catalog is wanted again, generate it from the git history
above (classified by conventional-commit type) with a CI check that fails when
the checked-in copy is stale — never hand-maintain it.
