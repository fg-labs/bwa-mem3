#!/usr/bin/env bash
# test/regression/meth_oracle.sh
#
# Regression: bwa-mem3 --meth Layer 1 -- the test/meth/test.sh harness.
#
# Was: the "Run --meth Layers 1-3" step inline in ci.yml, from when that harness
# still asserted bwameth structural/byte equivalence. Layers 2-3 were retired in
# D3; this is a thin wrapper that re-invokes what is left of it.
#
# Inputs:
#   The test/meth fixtures (ref.fa, t_R1.fastq.gz), copied into test/meth/ by
#   the caller -- see the "Copy bwa-meth fixtures into test/meth/" step in
#   ci.yml. test/meth/test.sh checks for them and says where they come from.
#
#   Nothing is read from the environment HERE, which is what puts this script in
#   the source-only block of test/regression/README.md -- that contract is about
#   what a script itself reads, and readme_contract_lint.sh decides it by reading
#   the file. The child honors ${SAMTOOLS:-samtools}, and the environment is
#   passed through to it untouched on purpose: pinning SAMTOOLS here would take
#   away the one knob a caller has for aiming the harness at a samtools that is
#   not on PATH, and would buy nothing -- no caller sets it today.
#
# This used to demand BWAMETH_DIR. Nothing reads it any more: test/meth/test.sh
# dropped its last use when Layers 2-3 (bwameth structural/byte equivalence)
# were retired, so the bwa-meth checkout is now a source of FIXTURES rather
# than of an oracle binary. Keeping the guard meant this wrapper refused to run
# over a variable no code consumed, while the precondition that can actually be
# missing -- the fixtures -- went unchecked.
set -euo pipefail

# Run from the repository root, the way the rest of this directory does. The
# harness path below is repo-relative, so resolving it against the caller's
# working directory made every invocation from anywhere else fail with a bare
# 127 -- and this script's whole point is that its failures say what is wrong.
cd "$(dirname "${BASH_SOURCE[0]}")/../.."

# Not `exec`: that forwards the child's exit status but emits neither of the
# PASS:/FAIL: markers every script in this directory owes its caller (see
# test/regression/README.md). The status is captured explicitly rather than
# left to `set -e`, which would abort before the FAIL: line could print. The
# child's own OK:/FAIL: lines still reach our streams untouched; this only adds
# the regression-level marker on top, and re-raises the child's status
# unchanged so a failure is never downgraded to a pass.
status=0
bash test/meth/test.sh || status=$?
if [[ $status -ne 0 ]]; then
    echo "FAIL: test/meth/test.sh exited $status"
    exit "$status"
fi

# Deliberately not "Layers 1-3": test/meth/test.sh retired Layers 2-3 (bwameth
# structural/byte equivalence) in D3 and now asserts Layer 1 only. Claiming the
# retired layers passed would be the same vacuous green this script exists to
# rule out. The ci.yml step name and the README table row that quotes it are
# renamed to match, so all three now say the same thing.
echo "PASS: --meth Layer 1 (test/meth/test.sh; Layers 2-3 retired in D3)"
