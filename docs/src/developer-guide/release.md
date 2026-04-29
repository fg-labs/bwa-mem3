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

## Tagging a release

1. Confirm `main` builds cleanly and all tests pass:

   ```bash
   make clean && make
   make test
   make docs
   ```

2. Update `NEWS.md` with the release notes for the new version.

3. Tag the release commit:

   ```bash
   git tag -s v0.X.Y -m "Release v0.X.Y"
   ```

   Use a signed tag (`-s`) when your GPG key is available. Annotated tags (`-a`) are acceptable when signing is not possible.

4. Push the tag to the `fg-labs` remote:

   ```bash
   git push fg-labs v0.X.Y
   ```

   Read the Docs activates a versioned build at `/v0.X.Y/` automatically when the tag appears on the remote.

5. Create a GitHub release from the tag via `gh`:

   ```bash
   gh release create v0.X.Y --repo fg-labs/bwa-mem3 \
     --title "bwa-mem3 v0.X.Y" \
     --notes-file <(sed -n '/^## \[0\.X\.Y\]/,/^## /p' NEWS.md | head -n -1)
   ```

> **Note — Tarball builds**
>
> Source tarballs created by GitHub (or `git archive`) do not include git history,
> so `git describe` fails and the version falls back to `FG_LABS_VERSION_FALLBACK`.
> For reproducible tarball builds, set `VERSION_STRING` explicitly on the command line:
> `make VERSION_STRING=v0.X.Y`.

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
