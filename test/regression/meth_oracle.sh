#!/usr/bin/env bash
# test/regression/meth_oracle.sh
#
# Regression: bwa-mem3 --meth Layers 1-3 match the bwa-meth oracle.
#
# Was: the "Run --meth Layers 1-3" step inline in ci.yml. This is a thin
# wrapper that re-invokes the existing test/meth/test.sh harness.
#
# Inputs:
#   The test/meth fixtures (ref.fa, t_R1.fastq.gz), copied into test/meth/ by
#   the caller -- see the "Copy bwa-meth fixtures into test/meth/" step in
#   ci.yml. test/meth/test.sh checks for them and says where they come from.
#
# This used to demand BWAMETH_DIR. Nothing reads it any more: test/meth/test.sh
# dropped its last use when Layers 2-3 (bwameth structural/byte equivalence)
# were retired, so the bwa-meth checkout is now a source of FIXTURES rather
# than of an oracle binary. Keeping the guard meant this wrapper refused to run
# over a variable no code consumed, while the precondition that can actually be
# missing -- the fixtures -- went unchecked.
set -euo pipefail

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
# rule out. The ci.yml step name and the README table row still say "Layers
# 1-3"; the README table maps each script to its step name, so those two have
# to be renamed together and that is a separate change from this one.
echo "PASS: --meth Layer 1 (test/meth/test.sh; Layers 2-3 retired in D3)"
