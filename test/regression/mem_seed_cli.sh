#!/usr/bin/env bash
# test/regression/mem_seed_cli.sh
#
# End-to-end checks of `bwa-mem3 mem` seeding: SMEM search over the FM index,
# re-seeding (-r, -y, -c), the minimum seed length (-k), SA lookup, and index
# variants given as input (SA sampling stride, a re-sa'd index). Each check
# runs under every SIMD tier this host can run (test/framework/mem_cli.sh:
# mem_cli_tiers), on the synthetic reference (test/fixtures/mem_synth_ref.awk,
# whose header lists its contigs and planted repeats) or the committed phiX,
# with every expected value worked out from the reference by hand.
#
# Checks:
#   -k    a 50 bp exact slice of ctg1 (501-550) maps at POS 501 with 50M under
#         the default -k 19, and comes out unmapped (FLAG 4) under -k 51: no
#         seed can be longer than the read.
#
# Inputs (env vars):
#   BWA_MEM3       path to the bwa-mem3 binary under test
#   MEM_CLI_TIERS  (optional) tiers to run instead of every runnable one
#   MEM_CLI_KEEP   (optional) 1 keeps the work directory
# No staged inputs: the reference is generated, the reads are sliced from it.

set -euo pipefail

# shellcheck source=test/framework/mem_cli.sh
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../framework" && pwd)/mem_cli.sh"

mem_cli_init mem_seed_cli
mem_cli_setup_synth
ref="$MEM_CLI_SYNTH"
work="$MEM_CLI_WORK"

mem_cli_fq r50 "$(mem_cli_seq "$ref" ctg1 501 50)" > "$work/reads.fq"

for tier in $(mem_cli_tiers); do
    sam="$work/k19.$tier.sam"
    mem_cli_mem "$tier" "$sam" "$ref" "$work/reads.fq"
    rec="$(mem_cli_primary "$sam" r50)"
    mem_cli_expect "$tier: -k 19 r50 FLAG RNAME POS CIGAR" "0 ctg1 501 50M" \
        "$(mem_cli_col "$rec" 2) $(mem_cli_col "$rec" 3) $(mem_cli_col "$rec" 4) $(mem_cli_col "$rec" 6)"

    sam="$work/k51.$tier.sam"
    mem_cli_mem "$tier" "$sam" -k 51 "$ref" "$work/reads.fq"
    mem_cli_expect "$tier: -k 51 r50 FLAG" 4 "$(mem_cli_col "$(mem_cli_primary "$sam" r50)" 2)"
done

if ((MEM_CLI_FAILS > 0)); then
    echo "FAIL: mem_seed_cli: $MEM_CLI_FAILS of $MEM_CLI_CHECKS checks failed"
    exit 1
fi
echo "PASS: mem_seed_cli: $MEM_CLI_CHECKS checks, $MEM_CLI_SKIPS skipped"
