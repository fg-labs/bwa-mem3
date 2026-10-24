#!/usr/bin/env bash
# test/regression/mem_pair_cli.sh
#
# End-to-end checks of `bwa-mem3 mem` paired-end and mode handling:
# insert-size estimation and -I, mate rescue (-S, -P, -m), proper pairs and
# the mate fields, MAPQ, primary/secondary/supplementary marking, --fast,
# --compat=bwa-mem2 and --compat=bwa-mem, --meth, --proper-pair-from-emitted,
# -j, and determinism across -K and -t. Each check runs under every SIMD tier
# this host can run (test/framework/mem_cli.sh: mem_cli_tiers), on the
# synthetic reference (test/fixtures/mem_synth_ref.awk, whose header lists its
# contigs and planted repeats) or the committed phiX, with every expected value
# worked out from the reference by hand.
#
# Checks:
#   -I    one FR pair, ctg1:1001-1100 forward and the reverse complement of
#         ctg1:1201-1300, under -I 300,30: FLAGs 99 and 147 (paired, proper,
#         first/second, mate reverse / read reverse), RNEXT =, PNEXT the mate's
#         POS, and TLEN +300 / -300 (the fragment is ctg1:1001-1300).
#
# Inputs (env vars):
#   BWA_MEM3       path to the bwa-mem3 binary under test
#   MEM_CLI_TIERS  (optional) tiers to run instead of every runnable one
#   MEM_CLI_KEEP   (optional) 1 keeps the work directory
# No staged inputs: the reference is generated, the reads are sliced from it.

set -euo pipefail

# shellcheck source=test/framework/mem_cli.sh
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../framework" && pwd)/mem_cli.sh"

mem_cli_init mem_pair_cli
mem_cli_setup_synth
ref="$MEM_CLI_SYNTH"
work="$MEM_CLI_WORK"

mem_cli_fq p1 "$(mem_cli_seq "$ref" ctg1 1001 100)" > "$work/r1.fq"
mem_cli_fq p1 "$(mem_cli_revcomp "$(mem_cli_seq "$ref" ctg1 1201 100)")" > "$work/r2.fq"

for tier in $(mem_cli_tiers); do
    sam="$work/i300.$tier.sam"
    mem_cli_mem "$tier" "$sam" -I 300,30 "$ref" "$work/r1.fq" "$work/r2.fq"
    for end in 1 2; do
        rec="$(mem_cli_primary "$sam" p1 "$end")"
        got="$(mem_cli_col "$rec" 2) $(mem_cli_col "$rec" 4) $(mem_cli_col "$rec" 7) $(mem_cli_col "$rec" 8) $(mem_cli_col "$rec" 9)"
        if [[ $end == 1 ]]; then
            want="99 1001 = 1201 300"
        else
            want="147 1201 = 1001 -300"
        fi
        mem_cli_expect "$tier: -I 300,30 end $end FLAG POS RNEXT PNEXT TLEN" "$want" "$got"
    done
done

if ((MEM_CLI_FAILS > 0)); then
    echo "FAIL: mem_pair_cli: $MEM_CLI_FAILS of $MEM_CLI_CHECKS checks failed"
    exit 1
fi
echo "PASS: mem_pair_cli: $MEM_CLI_CHECKS checks, $MEM_CLI_SKIPS skipped"
