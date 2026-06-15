# Summary

[Home](index.md)

# Getting Started

- [Installation](getting-started/installation.md)
- [Host requirements](getting-started/host-requirements.md)
- [Quick start: align paired-end FASTQs](getting-started/quick-align.md)
- [Quick start: methylation alignment](getting-started/quick-meth.md)
- [Quick start: shared-memory index](getting-started/quick-shm.md)

# User Guide

- [Indexing the reference](user-guide/indexing.md)
- [Aligning short reads (mem)](user-guide/aligning.md)
- [Output: SAM/BAM, headers, tags](user-guide/output.md)
- [Threading and resource use](user-guide/threading.md)
- [Memory allocator (mimalloc)](user-guide/allocator.md)
- [Tips and best practices](user-guide/tips.md)

# Performance

- [Overview](performance/overview.md)
- [SIMD dispatch matrix](performance/simd-dispatch.md)
- [PGO build](performance/pgo.md)
- [Tuning checklist](performance/tuning.md)

# Best Practices

- [Settings profiles: bwa drop-in vs recommended](best-practices/settings-profiles.md)
- [Build](best-practices/build.md)
- [Output format](best-practices/output-format.md)
- [Multi-sample workflows](best-practices/multi-sample.md)
- [Methylation defaults](best-practices/methylation.md)
- [Multi-architecture deployment](best-practices/multi-arch-deployment.md)
- [Anti-patterns](best-practices/anti-patterns.md)

# CLI Reference

- [Overview](cli/overview.md)
- [index](cli/index-cmd.md)
- [mem](cli/mem.md)
- [shm](cli/shm.md)
- [version](cli/version.md)

# Methylation Reference

- [Overview](methylation/overview.md)
- [bwameth.py drop-in mapping](methylation/bwameth-mapping.md)
- [Conversion details (C->T, G->A)](methylation/conversion.md)
- [SAM tags: XR, XG, XM](methylation/tags.md)
- [Chimera QC and header rewriting](methylation/post-processing.md)
- [Flags: --set-as-failed, --chimera-qc](methylation/flags.md)
- [Interop with external bwameth.py c2t](methylation/external-c2t.md)

# What's Different from bwa-mem2

- [Overview](whats-different/overview.md)
- [Equivalence with bwa-mem2](whats-different/equivalence.md)
- [Correctness fixes](whats-different/correctness.md)
- [Performance improvements](whats-different/performance.md)
- [Features](whats-different/features.md)
- [Architecture support](whats-different/arch-support.md)
- [Build & infrastructure](whats-different/build-infra.md)
- [`BASELINE_ARCH=avx512bw` build flag](whats-different/avx512-baseline.md)
- [Upstream PR status](whats-different/upstream-prs.md)

# Developer Guide

- [Building from source](developer-guide/building.md)
- [SIMD dispatch architecture](developer-guide/simd-dispatch.md)
- [Single-binary SIMD dispatch (x86)](developer-guide/launcher.md)
- [Apple Silicon / NEON port](developer-guide/neon-port.md)
- [Regression test framework](developer-guide/regression-tests.md)
- [Release process](developer-guide/release.md)
- [Branch and worktree conventions](developer-guide/branches.md)
- [Contributing](developer-guide/contributing.md)

# Related Projects

- [bwa-mem3-bench](related-projects/bwa-mem3-bench.md)
- [bwa-mem3-rs](related-projects/bwa-mem3-rs.md)
- [bwa-mem2 (upstream)](related-projects/bwa-mem2.md)
- [fgumi](related-projects/fgumi.md)
- [bwameth.py](related-projects/bwameth.md)

# Reference

- [Glossary](reference/glossary.md)
- [Citation](reference/citation.md)
- [License](reference/license.md)
- [Changelog](reference/changelog.md)
