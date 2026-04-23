#!/usr/bin/env bash
# Normalize a SAM stream for golden md5 comparison across builds.
# Strips @PG (contains version/cmdline) and md5sums the rest.
# Usage: bwa-mem2 mem ... | bench/normalize_sam.sh
set -euo pipefail
grep -v '^@PG' | md5sum | awk '{print $1}'
