# Fork changes vs. upstream

bwa-mem3 is a fork of [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2). This page
records the fork's relationship to upstream: **which fork changes were offered to
an upstream project, and where each stands today.** It is deliberately *not* an
exhaustive per-PR catalog — that record already exists, authoritatively and
without hand-maintenance, in git history and the GitHub PR list (see
[Where the full change list lives](#where-the-full-change-list-lives)). The
narrative *What's Different* pages explain the *why* behind each class of change.

> **Why this page is no longer a full table**
>
> Through early 2026 this page carried a hand-maintained row for every
> fork-carried PR (the "FG-MAIN-TABLE"). That model drifted badly — hundreds of
> merged PRs went unrecorded — because git history and the GitHub PR list are the
> real source of truth and a parallel hand-typed copy adds no information for the
> ~90 % of changes that are fork-only. What a copy *can* add — the upstream
> disposition of the handful of changes we tried to send upstream — is preserved
> below and is small enough to keep current by hand. If a full per-PR catalog is
> ever wanted back, *generate* it from git history (classified by
> conventional-commit type) with a CI check that fails when the checked-in copy
> is stale — never hand-maintain it again.

## Changes offered upstream

Ten early fork changes were submitted to an upstream project. Their current
disposition (re-verified 2026-08-19):

| Fork PR | Change | Class | Upstream | Status |
|---------|--------|-------|----------|--------|
| [#35](https://github.com/fg-labs/bwa-mem3/pull/35) | MQ tag + `-z` / `-u`/XB / `@HD` order port | Features | [lh3/bwa#330](https://github.com/lh3/bwa/pull/330) | **merged** (lh3/bwa only) |
| [#49](https://github.com/fg-labs/bwa-mem3/pull/49) | Batch `-H` header ingestion (O(n²) fix) | Performance | [bwa-mem2#204](https://github.com/bwa-mem2/bwa-mem2/pull/204) | open PR |
| [#50](https://github.com/fg-labs/bwa-mem3/pull/50) | Forward `CXXFLAGS`/`CPPFLAGS`/`LDFLAGS` to link | Build & infrastructure | [bwa-mem2#290](https://github.com/bwa-mem2/bwa-mem2/pull/290) | open PR |
| [#51](https://github.com/fg-labs/bwa-mem3/pull/51) | Guard kswv `rowMax` store on `nrow==0` | Correctness | [bwa-mem2#289](https://github.com/bwa-mem2/bwa-mem2/pull/289) | open PR (upstream covers AVX-512BW only) |
| [#52](https://github.com/fg-labs/bwa-mem3/pull/52) | Stamp `PACKAGE_VERSION` from `git describe` | Build & infrastructure | [bwa-mem2#284](https://github.com/bwa-mem2/bwa-mem2/pull/284) + [#283](https://github.com/bwa-mem2/bwa-mem2/issues/283) | open PR + open issue |
| [#54](https://github.com/fg-labs/bwa-mem3/pull/54) | Sanitize whitespace in `-R` for `@PG` `CL:` | Correctness | [bwa-mem2#293](https://github.com/bwa-mem2/bwa-mem2/issues/293) | open issue |
| [#55](https://github.com/fg-labs/bwa-mem3/pull/55) | Size SMEM buffers from observed max read length | Correctness | [bwa-mem2#238](https://github.com/bwa-mem2/bwa-mem2/pull/238) + [#210](https://github.com/bwa-mem2/bwa-mem2/issues/210) | PR closed unmerged; issue open |
| [#56](https://github.com/fg-labs/bwa-mem3/pull/56) | `--supp-rep-hard-cap` supp-MAPQ rescoring | Features | [bwa-mem2#260](https://github.com/bwa-mem2/bwa-mem2/issues/260) | open issue |
| [#1](https://github.com/fg-labs/bwa-mem3/pull/1) | Linux aarch64 (ARM64/NEON) build + CI | Architecture support | [bwa-mem2#288](https://github.com/bwa-mem2/bwa-mem2/pull/288) | **closed unmerged** (2026-07-10) |
| [#42](https://github.com/fg-labs/bwa-mem3/pull/42) | Emit `HN:i` total-hit-count tag | Features | [lh3/bwa#438](https://github.com/lh3/bwa/pull/438) | **closed unmerged** (2025-03) |

In short: one change merged upstream (#35, into Heng Li's `lh3/bwa`); three were
closed without merge, so the fork carries them permanently — arm64 support (#1),
the `HN:i` tag (#42), and the SMEM allocation fix (#55, whose companion upstream
*issue* remains open); the remaining six are still open bwa-mem2 PRs or issues.
Every fork change not listed here is fork-only — originated in bwa-mem3 with no
upstream counterpart.

## Where the full change list lives

The exhaustive "what changed on top of upstream" record is not hand-maintained.
Use, in decreasing order of rawness:

- **Git history** — `git log --reverse --no-merges master..main` lists every
  fork-carried commit on top of the upstream bwa-mem2 mirror (`master`).
- **The GitHub PR list** —
  [merged PRs against `main`](https://github.com/fg-labs/bwa-mem3/pulls?q=is%3Apr+is%3Amerged+base%3Amain),
  filterable and searchable.
- **The narrative [*What's Different*](../whats-different/overview.md) pages** —
  Correctness, Performance, Features, Architecture support, and Build &
  infrastructure — for the curated *why* behind each class of change.
- **The auto-generated [equivalence / divergence catalog](../whats-different/equivalence.md)** —
  for the reachable output divergences vs bwa-mem2, which is the part of the
  "is it still equivalent?" story that matters most.

## Upstream issues tracked but not yet fixed

These upstream issues are tracked in the bwa-mem3 issue list but do not yet have
a corresponding fix in `main`:

| Issue | Upstream reference | Notes |
|-------|--------------------|-------|
| Split-alignment evidence loss vs bwa 0.7.17 | [bwa-mem2#273](https://github.com/bwa-mem2/bwa-mem2/issues/273) | [issue #47](https://github.com/fg-labs/bwa-mem3/issues/47) — under investigation |
| MAPQ/coordinate parity vs bwa mem 0.7.18 | [bwa-mem2#262](https://github.com/bwa-mem2/bwa-mem2/issues/262), [bwa-mem2#246](https://github.com/bwa-mem2/bwa-mem2/issues/246), [bwa-mem2#239](https://github.com/bwa-mem2/bwa-mem2/issues/239) | [issue #48](https://github.com/fg-labs/bwa-mem3/issues/48) — tracking only |

---

**See also:**
[What's Different — Overview](../whats-different/overview.md) ·
[Correctness fixes](../whats-different/correctness.md) ·
[Performance improvements](../whats-different/performance.md) ·
[Features](../whats-different/features.md) ·
[Architecture support](../whats-different/arch-support.md) ·
[Build & infrastructure](../whats-different/build-infra.md)
