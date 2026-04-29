#!/usr/bin/env bash
# test/mimalloc_loaded_test.sh
#
# Asserts that a bwa-mem3 binary built with USE_MIMALLOC=1 reports a
# mimalloc version line from the `version` subcommand. Used in CI to prove
# the allocator isn't silently absent on a given build.
#
# Usage: test/mimalloc_loaded_test.sh <path-to-bwa-mem3>

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 <bwa-mem3-binary>" >&2
    exit 2
fi

bin="$1"
out="$("$bin" version 2>&1 || true)"

if grep -Eq '^mimalloc [0-9]+\.[0-9]+\.[0-9]+$' <<<"$out"; then
    echo "PASS: $bin reports mimalloc"
    grep '^mimalloc' <<<"$out"
    exit 0
else
    echo "FAIL: $bin does not report mimalloc in 'version' output"
    echo "---- actual output ----"
    echo "$out"
    echo "-----------------------"
    exit 1
fi
