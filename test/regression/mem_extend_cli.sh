#!/usr/bin/env bash
# test/regression/mem_extend_cli.sh
#
# End-to-end checks of `bwa-mem3 mem` chaining and extension: chaining and
# chain filtering (-w, -D, -W, -m), seed extension and the banded
# Smith-Waterman kernels, the scoring options (-A, -B, -O, -E, -L, -d, -T,
# -h, -a) and CIGAR generation. Each check runs under every SIMD tier this
# host can run (test/framework/mem_cli.sh: mem_cli_tiers), on the synthetic
# reference (test/fixtures/mem_synth_ref.awk, whose header lists its contigs
# and planted repeats) or the committed phiX, with every expected value worked
# out from the reference by hand.
#
# Checks:
#   -B    a 100 bp slice of ctg1 (1001-1100) with one substitution at read base
#         50 scores AS 99 - 6 = 93 under -B 6, with 100M and NM 1.
#
# Inputs (env vars):
#   BWA_MEM3       path to the bwa-mem3 binary under test
#   MEM_CLI_TIERS  (optional) tiers to run instead of every runnable one
#   MEM_CLI_KEEP   (optional) 1 keeps the work directory
# No staged inputs: the reference is generated, the reads are sliced from it.

set -euo pipefail

# shellcheck source=test/framework/mem_cli.sh
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../framework" && pwd)/mem_cli.sh"

mem_cli_init mem_extend_cli
mem_cli_setup_synth
ref="$MEM_CLI_SYNTH"
work="$MEM_CLI_WORK"

read1="$(mem_cli_seq "$ref" ctg1 1001 100)"
case "${read1:49:1}" in
    A) sub=T ;;
    *) sub=A ;;
esac
mem_cli_fq mm "$(mem_cli_subst "$read1" 50 "$sub")" > "$work/reads.fq"

for tier in $(mem_cli_tiers); do
    sam="$work/b6.$tier.sam"
    mem_cli_mem "$tier" "$sam" -B 6 "$ref" "$work/reads.fq"
    rec="$(mem_cli_primary "$sam" mm)"
    mem_cli_expect "$tier: -B 6 mm AS CIGAR NM" "93 100M 1" \
        "$(mem_cli_tag "$rec" AS) $(mem_cli_col "$rec" 6) $(mem_cli_tag "$rec" NM)"
done

if ((MEM_CLI_FAILS > 0)); then
    echo "FAIL: mem_extend_cli: $MEM_CLI_FAILS of $MEM_CLI_CHECKS checks failed"
    exit 1
fi
echo "PASS: mem_extend_cli: $MEM_CLI_CHECKS checks, $MEM_CLI_SKIPS skipped"
