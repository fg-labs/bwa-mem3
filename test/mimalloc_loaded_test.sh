#!/usr/bin/env bash
# test/mimalloc_loaded_test.sh
#
# Asserts that a bwa-mem3 binary built with USE_MIMALLOC=1 reports mimalloc as
# the *active* allocator from the `version` subcommand. The version line now
# always carries a status suffix — "(active)" when standard malloc is routed to
# mimalloc, or "(linked but NOT overriding malloc)" when the library is linked
# but not intercepting malloc. Used in CI to prove the allocator isn't silently
# absent or merely linked-but-inert on a given build.
#
# Usage: test/mimalloc_loaded_test.sh <path-to-bwa-mem3>

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "usage: $0 <bwa-mem3-binary>" >&2
    exit 2
fi

bin="$1"
out="$("$bin" version 2>&1 || true)"

if grep -Eq '^mimalloc [0-9]+\.[0-9]+\.[0-9]+ \(active\)$' <<<"$out"; then
    echo "PASS: $bin reports mimalloc as the active allocator"
    grep '^mimalloc' <<<"$out"
    exit 0
elif grep -Eq '^mimalloc [0-9]+\.[0-9]+\.[0-9]+ ' <<<"$out"; then
    echo "FAIL: $bin links mimalloc but it is NOT the active allocator"
    grep '^mimalloc' <<<"$out"
    exit 1
else
    echo "FAIL: $bin does not report mimalloc in 'version' output"
    echo "---- actual output ----"
    echo "$out"
    echo "-----------------------"
    exit 1
fi
