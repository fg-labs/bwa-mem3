# version

`bwa-mem3 version` prints the release version of the binary and, when
mimalloc is compiled in, the mimalloc version. It is a quick way to confirm
which build is on `PATH` and whether the allocator is active.

## Synopsis

```text
{{#include ../../_generated/cli/version.txt}}
```

## Common usage

Confirm the installed version and allocator:

```sh
./bwa-mem3 version
```

A typical run prints two lines (the mimalloc line goes to stderr and the
version string to stdout, so the order in a merged stream is not guaranteed):

```text
mimalloc 3.3.0
v0.1.0-pre-7-g613a4dc
```

The first line is the mimalloc version (present only when `USE_MIMALLOC=1`,
the default). The second line is the bwa-mem3 version string, derived from
`git describe` at build time and stored as `PACKAGE_VERSION` in the binary.
When building from a tarball without git history, the fallback value is set
via `FG_LABS_VERSION_FALLBACK` at compile time.

## Notes / Gotchas

> **Note — mimalloc line goes to stderr, version to stdout**
>
> The mimalloc line is printed to stderr and the version string to stdout.
> Scripts that capture the version should redirect stderr appropriately or
> use `bwa-mem3 version 2>/dev/null`.
>
> **Tip — No mimalloc line means USE_MIMALLOC=0**
>
> If only the version string appears and no mimalloc line, the binary was built
> without the bundled allocator (`make USE_MIMALLOC=0`). See
> [User Guide — Memory allocator (mimalloc)](../user-guide/allocator.md) for
> when this is appropriate.

---

**See also:**
[User Guide — Memory allocator (mimalloc)](../user-guide/allocator.md) ·
[Developer Guide — Release process](../developer-guide/release.md) ·
[Getting Started — Installation](../getting-started/installation.md) ·
[What's Different — Build & infrastructure](../whats-different/build-infra.md)
