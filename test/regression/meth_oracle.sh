#!/usr/bin/env bash
# test/regression/meth_oracle.sh
#
# Regression: bwa-mem3 --meth Layers 1-3 match the bwa-meth oracle.
#
# Was: the "Run --meth Layers 1-3" step inline in ci.yml. This is a thin
# wrapper that re-invokes the existing test/meth/test.sh harness.
#
# Inputs:
#   BWAMETH_DIR — path to the checked-out bwa-meth oracle repo
#   (and test/meth fixtures copied into test/meth/ by the caller —
#    see ci.yml "Copy bwa-meth fixtures into test/meth/" step)
set -euo pipefail

: "${BWAMETH_DIR:?BWAMETH_DIR must be set}"

exec bash test/meth/test.sh
