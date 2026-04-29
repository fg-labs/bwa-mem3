#!/usr/bin/env bash
# Normalize a SAM stream for golden md5 comparison across builds.
# Strips @PG (contains version/cmdline) and md5sums the rest.
# Usage: bwa-mem3 mem ... | bench/normalize_sam.sh
set -euo pipefail
if command -v md5sum >/dev/null 2>&1; then
  grep -v '^@PG' | md5sum | awk '{print $1}'
else
  # macOS: BSD `md5 -q` emits just the hash.
  grep -v '^@PG' | md5 -q
fi
