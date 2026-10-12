#!/usr/bin/env bash
# test/regression/thp_check_smoke.sh
#
# Smoke test for the transparent-huge-page (THP) fallback check in `mem`
# (src/bwa_thp.cpp), which runs from the pipeline's write step on Linux.
#
# The unit tests drive the check's state machine through an injected reader;
# this one runs the real binary, so it covers what they cannot: the write-step
# hook, the real /proc and /sys reads, and the `-v 4` debug line. On phiX:
#
#   * `mem -v 4` logs a `[D::main_mem] THP check batch 1:` line that reports a
#     real smaps walk (a VMA count above 0), or a `check disabled: <reason>`
#     line when this host is out of the check's scope (base page not 4 KB, THP
#     mode unreadable, another mTHP size enabled, smaps unreadable); the latter
#     is reported as SKIP for the walk check.
#   * no THP warning is printed: phiX is far below the 4 GiB floor.
#   * the alignment records are identical between `-v 3` and `-v 4` (the check
#     is stderr-only; only the @PG CL: field may differ).
#
# Linux only: off Linux the check is compiled out, so the script reports SKIP.
#
# Inputs:
#   BWA_MEM3 — path to bwa-mem3 binary
#   FIXTURES — directory containing phix.fa and reads.fa (default: test/fixtures)
set -euo pipefail

: "${BWA_MEM3:?BWA_MEM3 must be set}"
FIXTURES="${FIXTURES:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../fixtures" && pwd)}"

if [[ "$(uname -s)" != Linux ]]; then
    echo "SKIP: THP check is Linux-only (host is $(uname -s))"
    exit 0
fi

src_ref="$FIXTURES/phix.fa"
reads="$FIXTURES/reads.fa"
[[ -x "$BWA_MEM3" ]] || {
    echo "FAIL: binary not executable: $BWA_MEM3" >&2
    exit 1
}
[[ -s "$src_ref" ]] || {
    echo "FAIL: phix.fa missing: $src_ref" >&2
    exit 1
}
[[ -s "$reads" ]] || {
    echo "FAIL: reads.fa missing: $reads" >&2
    exit 1
}

# Index a private copy so the test never writes into the fixtures tree.
mdir="$(mktemp -d)"
trap 'rm -rf "$mdir"' EXIT
ref="$mdir/phix.fa"
cp "$src_ref" "$ref"
"$BWA_MEM3" index "$ref" > /dev/null 2>&1 || {
    echo "FAIL: index phix.fa" >&2
    exit 1
}

run_mem() { # $1 = verbosity, $2 = SAM out, $3 = stderr out
    "$BWA_MEM3" mem -v "$1" "$ref" "$reads" > "$2" 2> "$3" || {
        echo "FAIL: mem -v $1 exited nonzero" >&2
        cat "$3" >&2
        exit 1
    }
}
run_mem 3 "$mdir/v3.sam" "$mdir/v3.err"
run_mem 4 "$mdir/v4.sam" "$mdir/v4.err"

fails=0
skipped=0

# The check never touches output: records identical, @PG aside.
grep -v '^@PG' "$mdir/v3.sam" > "$mdir/v3.body"
grep -v '^@PG' "$mdir/v4.sam" > "$mdir/v4.body"
if [[ ! -s "$mdir/v3.body" ]] || ! grep -qv '^@' "$mdir/v3.body"; then
    echo "FAIL: mem -v 3 produced no alignment records" >&2
    fails=$((fails + 1))
elif ! cmp -s "$mdir/v3.body" "$mdir/v4.body"; then
    echo "FAIL: alignment records differ between -v 3 and -v 4" >&2
    diff "$mdir/v3.body" "$mdir/v4.body" | head -20 >&2 || true
    fails=$((fails + 1))
else
    echo "PASS: alignment records identical between -v 3 and -v 4"
fi

# phiX is far below the 4 GiB floor: no THP warning at either verbosity.
if grep -qE '^\[W::main_mem\] .*(transparent huge pages|THP-eligible)' "$mdir/v3.err" "$mdir/v4.err"; then
    echo "FAIL: THP warning printed on a phiX run" >&2
    grep -hE '^\[W::main_mem\]' "$mdir/v3.err" "$mdir/v4.err" >&2
    fails=$((fails + 1))
else
    echo "PASS: no THP warning on a phiX run"
fi

# The -v 3 run logs no THP debug line.
if grep -q '^\[D::main_mem\] THP check' "$mdir/v3.err"; then
    echo "FAIL: THP debug line printed at -v 3" >&2
    fails=$((fails + 1))
fi

line="$(grep -m1 '^\[D::main_mem\] THP check batch 1:' "$mdir/v4.err" || true)"
if [[ -z $line ]]; then
    echo "FAIL: no '[D::main_mem] THP check batch 1:' line at -v 4" >&2
    cat "$mdir/v4.err" >&2
    fails=$((fails + 1))
elif [[ $line == *"check disabled: "* ]]; then
    echo "SKIP: THP check out of scope on this host (${line#*check disabled: })"
    skipped=1
else
    vmas="$(sed -nE 's/.*, ([0-9]+) VMAs,.*/\1/p' <<< "$line")"
    if [[ -z $vmas ]] || ((vmas <= 0)); then
        echo "FAIL: THP sample line reports no VMAs: $line" >&2
        fails=$((fails + 1))
    else
        echo "PASS: THP check sampled batch 1 at -v 4 ($vmas VMAs)"
    fi
fi

if ((fails > 0)); then
    echo "FAIL: $fails THP check smoke assertion(s) failed" >&2
    exit 1
fi
if ((skipped)); then
    echo "SKIP: THP sample check skipped; the output and no-warning checks passed"
    exit 0
fi
echo "PASS: THP check smoke"
