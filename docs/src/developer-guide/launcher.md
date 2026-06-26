# Single-binary SIMD dispatch (x86)

On x86 Linux and x86 macOS, `bwa-mem3` is a **single binary** that
contains compiled kernels for every supported SIMD tier
(`sse41` / `sse42` / `avx` / `avx2` / `avx512bw`). At startup the binary
detects the host CPU's capabilities and selects the matching tier in
process, without `fork` or `exec`. There is no separate launcher binary
and no `bwa-mem3.<tier>` variant files on disk.

ARM / Apple Silicon does not need tier dispatch at all: there is only
one NEON instruction-set level across current ARM64 CPUs, so the arm64
build is a single binary with one kernel TU. The dispatch machinery
described below is only meaningful on x86.

This design replaces the multi-binary `execv` launcher inherited from
bwa-mem2. The motivation, validation, and trade-offs are tracked in
[PR #83](https://github.com/fg-labs/bwa-mem3/pull/83); the AVX-512
auto-vectorization cap that ships alongside it is documented in
[`BASELINE_ARCH=avx512bw` build flag](../whats-different/avx512-baseline.md).

## What the build produces

```bash
make            # default: single multi-tier binary, BASELINE_ARCH=avx2
make single     # explicit alias of the default target
```

Produces **one** file in the repo root:

| File | Contains | Non-kernel TU compile flags |
|---|---|---|
| `bwa-mem3` | All 5 x86 tier kernels + dispatcher + non-kernel TUs | `BASELINE_ARCH` (default `avx2`) |

The five kernel translation units listed in the Makefile's `KERNEL_SRCS`
(`bandedSWA.cpp`, `kswv.cpp`, `ksw.cpp`, `sam_encode.cpp`) are compiled
five times each, once per tier, with tier-specific `-m...` flags. Every
non-kernel TU is compiled once at the `BASELINE_ARCH` tier.
`BASELINE_ARCH` defaults to `avx2` (PR #84) and can be set on the
`make` line:

```bash
make BASELINE_ARCH=avx512bw       # for an AVX-512BW-only fleet
make BASELINE_ARCH=sse41          # for pre-Haswell hosts (~10–15% slower on AVX2)
```

Lowering `BASELINE_ARCH` reduces the supported host floor and is the
documented escape hatch for vintage hardware. Raising it locks the
binary to that host class and disables the host-floor precheck for
lower tiers. The `bwa-mem3 version` banner prints the resulting
`SIMD floor:` line so operators can confirm the build matches the
intended deployment surface — see
[Host requirements](../getting-started/host-requirements.md) and
[`BASELINE_ARCH=avx512bw` build flag](../whats-different/avx512-baseline.md).

For ARM, `make arm64` produces a single binary with a single NEON
kernel TU; no dispatch table is generated.

## Runtime tier selection

`src/simd_dispatch.cpp` provides three pieces:

- `bwamem3_simd_init()` — idempotent initializer called from
  `main.cpp`. Caches the host's raw capability into a file-scope
  `g_host_capability` and the effective dispatch tier into a separate
  `g_tier` (the two differ when `BWAMEM3_FORCE_TIER` is set).
- An enum of supported tiers (`sse41` → `sse42` → `avx` → `avx2` → `avx512bw`,
  plus `neon` on arm64) and `bwamem3_simd_tier_name()` for stderr
  reporting.
- Per-kernel factory functions (`make_bsw_kernel_<tier>`,
  `make_kswv_kernel_<tier>`) and free-function dispatch wrappers
  (`ksw_extend2`, `ksw_global2`, `ksw_extend`, `ksw_global`,
  `ksw_align2`, `ksw_align`, `sam_encode_*`) that switch on `g_tier`
  and call into the matching mangled per-tier symbol.

x86 detection uses `__builtin_cpu_supports` directly; arm64 reports
`neon` unconditionally. The selection happens once at startup and the
result is cached in a TU-level global — subsequent kernel calls pay a
single indirect-call overhead through a vtable (for the
`BandedPairWiseSW` / `kswv` factories) or an `extern "C"` wrapper (for
the `ksw_*` free functions). Per PR #83 measurement, the indirect call
costs ~0.3 ns after BTB warm-up, so a 1M-read alignment with ~100M
kernel calls adds roughly 30 ms — well below run-to-run noise on every
tested host.

## Symbol mangling per tier

The compile-time machinery — the `KERNEL_VARIANT` symbol-rename scheme, the
`IBandedPairWiseSW` / `Ikswv` interface split that keeps the dispatcher TU free
of class-layout knowledge, and the ODR-collision avoidance — is documented once
in [SIMD dispatch architecture → Per-tier compilation and symbol mangling](simd-dispatch.md#per-tier-compilation-and-symbol-mangling).
This page covers the *runtime* side: tier selection, host-floor enforcement, and
distribution.

## Environment overrides

Two environment variables exposed at runtime:

| Variable | Behavior |
|---|---|
| `BWAMEM3_FORCE_TIER=<tier>` | Force the dispatcher to use `<tier>` (one of `sse41` `sse42` `avx` `avx2` `avx512bw`). **Downgrade-only:** requests above the detected host tier (which would SIGILL on the first wider instruction) and unrecognized names are rejected with a stderr warning and the dispatcher falls back to the detected tier. Replaces the prior "exec the `bwa-mem3.sse41` binary" pattern for A/B regression testing on AVX-512 hosts. |
| `BWAMEM3_DEBUG_SIMD=1` | Print a one-line `[I::bwamem3_simd_init_body]` banner at startup naming the build baseline (`g_build_tier`), the detected host capability, and the resolved dispatch tier. Also enables the build-baseline-vs-host gap warning that PR #84 originally emitted unconditionally and PR #86 demoted to debug-only. |

Both are read once during `bwamem3_simd_init()` and ignored after that
call returns.

## Host-floor enforcement

`bwa-mem3 mem`, `bwa-mem3 index`, and `bwa-mem3 shm` all call
`bwamem3_enforce_host_floor()` early in `main()` (PR #95). The check
compares `g_host_capability` against the compile-time `g_build_tier`
(derived from compiler predefined macros, reflecting whichever
`BASELINE_ARCH` was set at build time) and exits with code `2` and an
`[E::bwamem3]` message naming the gap if the host cannot execute the
binary's compiled-in instructions. This converts what would otherwise
be an unhelpful SIGILL deep in alignment into a clean abort at startup.

Diagnostic invocations opt out: `bwa-mem3 version`,
`bwa-mem3 <subcommand> --help`, and `bwa-mem3 <subcommand> -h` always
succeed regardless of host capability, so operators can introspect a
binary on a host that cannot run alignment. The `version` command
prints `SIMD floor:` (the build's required minimum) and
`SIMD runtime:` (the resolved tier) on stdout; on a too-old host it
also emits a `[W::bwa-mem3]` warning on stderr.

The `simd_dispatch.cpp` translation unit is compiled at `BASELINE_ARCH`
like every other non-kernel TU; an earlier draft forced it to
`-march=x86-64` to keep the precheck SIGILL-safe, but that broke
`g_build_tier` (a `static constexpr` derived from `__AVX2__` /
`__SSE4_1__` / etc., which are only defined when the matching `-m`
flag is in scope) — every binary reported its floor as `scalar` and
the precheck became a no-op. In practice the precheck path is
scalar-only (`std::call_once`, integer comparisons, `getenv`,
`snprintf`, `fputs`, `exit`) with no array loops the compiler could
autovectorize, so it stays SIGILL-safe even when `BASELINE_ARCH=avx2`
(or higher) for the rest of the binary.

## Per-tier parity validation

`test/regression/all_tiers_parity.sh` runs `bwa-mem3 mem` with
`BWAMEM3_FORCE_TIER` walking the full ladder
(`sse41 → sse42 → avx → avx2 → avx512bw`) on the same input and
diff's the BAM output. The expected result is byte-identical SAM
across every tier; any divergence is a bug in either a kernel TU or
the per-kernel factory wiring. CI runs this script on the x86 matrix
row.

## Trade-offs vs the prior multi-binary launcher

| Property | Pre-PR-#83 (multi-binary `execv`) | Current (single binary, in-process dispatch) |
|---|---|---|
| Install size | ~120 MB (5 ISA binaries + launcher) | ~25 MB (one binary) |
| Build cost | 5 sequential clean rebuilds + launcher | One parallel build |
| Process model | `bwa-mem3` (launcher) → `execv` → `bwa-mem3.<tier>` | One process, one `main()` |
| Per-call overhead | Direct call (tier fixed at launch via separate binary) | Indirect call through factory vtable or `extern "C"` wrapper (~0.3 ns / call) |
| Non-kernel auto-vectorization | At each binary's compile tier | At `BASELINE_ARCH` (default `avx2`); raise via `BASELINE_ARCH=` |
| Tier override | Run the `.<tier>` binary directly | `BWAMEM3_FORCE_TIER=<tier>` (downgrade-only) |
| `runsimd.cpp` (220-line launcher) | Required | Removed |

The ~0.3 ns indirect-call cost is amortized across alignment work and
has not been measurable in any bench cell. The non-kernel
auto-vectorization at `BASELINE_ARCH` is what closes the gap PR #84
identified after PR #83 originally regressed by silently hardcoding
the non-kernel compile to `sse41`.

## Distribution layout

For deployment on any x86 host meeting the build's floor:

```text
bin/
  bwa-mem3       ← single binary, dispatches in-process
```

For ARM:

```text
bin/
  bwa-mem3       ← single binary, NEON kernels only
```

No `.<tier>`-suffixed companion files are produced or needed. When
shipping a Docker image intended for a mixed-microarch fleet, build at
the lowest expected tier (e.g. `BASELINE_ARCH=avx2` for "AVX2 and
newer") — the runtime dispatcher will still pick AVX-512BW kernels on
AVX-512 hosts via the per-tier factory tables. See
[Multi-architecture deployment](../best-practices/multi-arch-deployment.md)
for the `docker buildx` manifest-list recipe.

## The `mem` SIMD banner

The legacy `Executing in AVX2 mode!!` banner is gone. Use either:

- `bwa-mem3 version` — prints `SIMD floor:` and `SIMD runtime:` lines
  on stdout (always available, no alignment required).
- `BWAMEM3_DEBUG_SIMD=1 bwa-mem3 mem …` — prints a one-line
  `[I::bwamem3_simd_init_body]` banner on stderr at the start of the
  run.

---

**See also:**
[SIMD dispatch architecture](simd-dispatch.md) ·
[Apple Silicon / NEON port](neon-port.md) ·
[Building from source](building.md) ·
[Performance → SIMD dispatch matrix](../performance/simd-dispatch.md) ·
[Host requirements](../getting-started/host-requirements.md) ·
[`BASELINE_ARCH=avx512bw` build flag](../whats-different/avx512-baseline.md) ·
[Multi-architecture deployment](../best-practices/multi-arch-deployment.md)
