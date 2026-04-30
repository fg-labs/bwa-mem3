# bwa-mem3-rs

bwa-mem3-rs is a Rust crate that provides idiomatic bindings to the bwa-mem family of
short-read aligners — bwa (original), bwa-mem2, and bwa-mem3. It exposes a safe Rust
API over the underlying C++ alignment engine, allowing Rust programs to index a
reference, configure alignment parameters, and align reads without shelling out to an
external process. The bindings link statically against the chosen backend, so a binary
built with bwa-mem3-rs carries the aligner and its SIMD kernels as a self-contained
artifact.

## When you'd use it

Use bwa-mem3-rs when you are building a Rust bioinformatics tool or pipeline that needs
short-read alignment as an in-process library call rather than a subprocess invocation.
It is especially useful when latency between reads arriving and alignments being available
matters (no process-startup overhead), or when you want tight integration between the
aligner's output and downstream Rust code such as UMI grouping, consensus calling, or
duplicate marking.

## How it relates to bwa-mem3

bwa-mem3-rs targets bwa-mem3 as its primary high-performance backend. It is the intended
integration path for [fgumi](fgumi.md) and other Fulcrum Genomics tools that need
alignment as a library dependency. Changes to bwa-mem3's public API, flag semantics,
or output format are coordinated with bwa-mem3-rs to keep the bindings current.

## Links

- GitHub: <https://github.com/fg-labs/bwa-mem3-rs>
- License: MIT

---

**See also:**
[fgumi](fgumi.md) ·
[bwa-mem3-bench](bwa-mem3-bench.md) ·
[Aligning short reads (mem)](../user-guide/aligning.md) ·
[Developer Guide — Contributing](../developer-guide/contributing.md)
