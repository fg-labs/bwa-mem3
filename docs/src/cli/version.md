# version

`bwa-mem3 version` prints the release version, the build's compiled-in
SIMD floor, the SIMD tier resolved at runtime, and (when mimalloc is
compiled in) the mimalloc version. It is the canonical way to confirm
which build is on `PATH`, what host class it requires, and what kernel
path it will dispatch to.

`bwa-mem3 version` always exits 0 — even on a host below the build's
SIMD floor — so operators can introspect a binary on a host that cannot
actually run alignment. `bwa-mem3 <subcommand> --help` and `-h` share
the same property.

## Synopsis

```text
{{#include ../../_generated/cli/version.txt}}
```

## Common usage

Confirm the installed version, SIMD floor, and resolved tier:

```sh
./bwa-mem3 version
```

A typical run on an AVX-512BW host with the default `BASELINE_ARCH=avx2`
build prints (mimalloc line on stderr, the rest on stdout — order in a
merged stream is not guaranteed):

```text
v0.3.0-12-gabcdef1
SIMD floor: avx2 (x86-64-v3, Haswell 2013+); kernels: sse41 sse42 avx avx2 avx512bw
SIMD runtime: avx512bw (BWAMEM3_FORCE_TIER unset)
mimalloc 3.3.0
```

- **version line** — bwa-mem3's release string, derived from
  `git describe` at build time and stored as `PACKAGE_VERSION` in the
  binary. When building from a tarball without git history, the
  fallback value is set via `FG_LABS_VERSION_FALLBACK` at compile time.
- **`SIMD floor:`** — the compile-time minimum the binary requires.
  Set by `BASELINE_ARCH` (default `avx2`) and listed alongside the
  per-tier kernel set the binary carries.
- **`SIMD runtime:`** — the tier resolved at startup by
  `__builtin_cpu_supports` (or `BWAMEM3_FORCE_TIER` if set in the
  environment). On arm64 this is always `neon`.
- **mimalloc line** — present only when `USE_MIMALLOC=1` (the default).

On a host below the SIMD floor, `version` also writes a
`[W::bwa-mem3]` warning on **stderr** identifying the gap (and the
alignment subcommands will refuse to run with exit code 2 — see
[Host requirements](../getting-started/host-requirements.md) for the
exit-2 message format and rebuild instructions).

## Notes / Gotchas

> **Tip — `version | grep` is safe in CI**
>
> The version, `SIMD floor:`, and `SIMD runtime:` lines all go to stdout;
> the mimalloc line and any host-below-floor warning go to stderr. So
> `bwa-mem3 version | grep '^SIMD'` works in CI scripts even on hosts
> that cannot run alignment. Use `2>/dev/null` to suppress the mimalloc
> and warning lines if you want stdout only.
>
> **Tip — No mimalloc line means USE_MIMALLOC=0**
>
> If no mimalloc line appears, the binary was built without the bundled
> allocator (`make USE_MIMALLOC=0`). See
> [User Guide — Memory allocator (mimalloc)](../user-guide/allocator.md)
> for when this is appropriate.

---

**See also:**
[User Guide — Memory allocator (mimalloc)](../user-guide/allocator.md) ·
[Developer Guide — Release process](../developer-guide/release.md) ·
[Getting Started — Installation](../getting-started/installation.md) ·
[What's Different — Build & infrastructure](../whats-different/build-infra.md)
