# PGO Build

Profile-Guided Optimization (PGO) is a two-pass compiler technique. In the first pass (`pgo-generate`) the compiler inserts counters into every branch, call site, and loop back-edge. You run a representative training workload against the instrumented binary so those counters accumulate real branch-probability data. In the second pass (`pgo-use`) the compiler recompiles every translation unit using the collected profiles to make better inlining, branch-prediction, and code-layout decisions.

bwa-mem3's Makefile provides three targets that implement this workflow.

## Observed gains

On Apple Silicon (M-series), PGO delivered approximately **3%** throughput improvement over the native NEON build. The gain on x86 depends on the workload — short-read paired-end alignment on `avx2` or `avx512bw` hardware typically sees 2–5%. PGO is most useful when you will run the same binary on the same hardware against the same workload repeatedly (e.g. a production pipeline node). It is not worth the extra build time for one-off or exploratory runs.

## Workflow

### Step 1: Build the instrumented binary

```bash
make pgo-generate
```

By default `PGO_ARCH` is set to `arm64` on Apple Silicon / aarch64 hosts and `native` on x86 hosts. To target a specific ISA, pass `PGO_ARCH` explicitly:

```bash
make pgo-generate PGO_ARCH=avx2
```

This produces a binary named `bwa-mem3.pgo-instr` (or `bwa-mem3.pgo-instr.avx2` for non-default arch). Profiles are written to the directory `pgo_profiles/` by default. Override with `PGO_PROFILE_DIR`:

```bash
make pgo-generate PGO_ARCH=avx2 PGO_PROFILE_DIR=/scratch/pgo-profiles-avx2
```

### Step 2: Run the training workload

Run a workload that is representative of your production use. A single-end or paired-end alignment run against the same reference and similar read length is sufficient. A larger training run produces more stable profiles but 5–10 million read pairs is generally enough.

```bash
./bwa-mem3.pgo-instr mem -t 16 ref.fa R1.fq.gz R2.fq.gz > /dev/null
```

The run discards output so you are measuring the alignment work alone.

> **Tip — Training workload size**
>
> Aim for a training run that exercises the same code paths as your production workload. If you align 150 bp paired-end reads in production, train on 150 bp reads. If you use `--meth`, include a methylation alignment run in training. A few million read pairs is sufficient; a full WGS run provides diminishing returns.

### Step 3: Build the optimized binary

```bash
make pgo-use
```

Or with matching arch and profile dir:

```bash
make pgo-use PGO_ARCH=avx2 PGO_PROFILE_DIR=/scratch/pgo-profiles-avx2
```

This produces `bwa-mem3.pgo` (or `bwa-mem3.pgo.avx2`). The binary is ready to use in production.

### Step 4: Clean up instrumentation artifacts

```bash
make pgo-clean
```

This removes the profile directory and all `bwa-mem3.pgo-instr*` and `bwa-mem3.pgo*` files.

## Multi-arch builds with PGO

Each architecture requires its own profile because the instrumentation counters are embedded in arch-specific code. Run the full three-step workflow once per arch and keep the profiles in separate directories:

```bash
# AVX2 profile
make pgo-generate PGO_ARCH=avx2 PGO_PROFILE_DIR=pgo_profiles_avx2
./bwa-mem3.pgo-instr.avx2 mem -t 16 ref.fa R1.fq.gz R2.fq.gz > /dev/null
make pgo-use PGO_ARCH=avx2 PGO_PROFILE_DIR=pgo_profiles_avx2

# AVX-512BW profile (separate host or same host with matching CPU)
make pgo-generate PGO_ARCH=avx512bw PGO_PROFILE_DIR=pgo_profiles_avx512bw
./bwa-mem3.pgo-instr.avx512bw mem -t 16 ref.fa R1.fq.gz R2.fq.gz > /dev/null
make pgo-use PGO_ARCH=avx512bw PGO_PROFILE_DIR=pgo_profiles_avx512bw
```

> **Warning — Profile portability**
>
> Profile data collected on one microarchitecture is not portable to a different one. An AVX2 profile collected on a Haswell CPU will not improve — and may pessimize — an AVX-512BW build run on a Sapphire Rapids CPU. Always collect profiles on the same hardware class where the optimized binary will run.

## PGO and the single-binary multi-tier build

The PGO targets produce one optimized binary for a single `arch=` target. They do not yet rebuild the default `make single` multi-tier binary's per-tier kernel TUs. If you need PGO across more than one host class, build and profile each `arch=` variant separately and deploy whichever matches the target fleet — `bwa-mem3 version` will report the resolved tier so you can confirm. PGO for the in-process multi-tier dispatch path is tracked as a future enhancement.

## Relationship to LTO

`make lto-build` produces a Link-Time Optimization binary; `make pgo-use` produces a PGO-optimized binary. Both are independent opt-in targets. You can combine them by passing `-flto` (or `-flto=thin` for clang) as part of `EXTRA_CXXFLAGS` during the `pgo-use` step, but the combination has not been systematically benchmarked. In practice, LTO and PGO each provide modest single-digit gains; their interaction is compiler-specific.

---

**See also:**
[Performance overview](overview.md) ·
[SIMD dispatch matrix](simd-dispatch.md) ·
[Optimization checklist](../best-practices/optimization-checklist.md) ·
[Best Practices — Build](../best-practices/build.md) ·
[Developer Guide — Building from source](../developer-guide/building.md)
