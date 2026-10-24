#!/usr/bin/env bash
# test/regression/mem_io_cli.sh
#
# End-to-end checks of `bwa-mem3 mem` input and output: FASTQ/FASTA readers
# (plain and gzip, -p interleaved), -C, -R/-H headers, -o, --bam output, the
# SAM record fields and tags (FLAG, CIGAR, NM, MD, AS, XS, XA, SA, MQ, HN),
# -Y, -M, -V, -5, clipping and unmapped records. Each check runs under every
# SIMD tier this host can run (test/framework/mem_cli.sh: mem_cli_tiers), on
# the synthetic reference (test/fixtures/mem_synth_ref.awk, whose header lists
# its contigs and planted repeats) or the committed phiX, with every expected
# value worked out from the reference by hand.
#
# Checks:
#   @SQ   one line per contig, in reference order, with the generator's
#         lengths: ctg1 5000, ctg2 4000, ctg3 3000, ctg4 600.
#
# Inputs (env vars):
#   BWA_MEM3       path to the bwa-mem3 binary under test
#   MEM_CLI_TIERS  (optional) tiers to run instead of every runnable one
#   MEM_CLI_KEEP   (optional) 1 keeps the work directory
# No staged inputs: the reference is generated, the reads are sliced from it.

set -euo pipefail

# shellcheck source=test/framework/mem_cli.sh
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../framework" && pwd)/mem_cli.sh"

mem_cli_init mem_io_cli
mem_cli_setup_synth
ref="$MEM_CLI_SYNTH"
work="$MEM_CLI_WORK"

mem_cli_fq exact "$(mem_cli_seq "$ref" ctg1 1001 100)" > "$work/reads.fq"

for tier in $(mem_cli_tiers); do
    sam="$work/default.$tier.sam"
    mem_cli_mem "$tier" "$sam" "$ref" "$work/reads.fq"
    mem_cli_expect "$tier: @SQ lines" \
        "SN:ctg1 LN:5000;SN:ctg2 LN:4000;SN:ctg3 LN:3000;SN:ctg4 LN:600;" \
        "$(awk -F '\t' '$1 == "@SQ" { printf "%s %s;", $2, $3 }' "$sam")"
done

if ((MEM_CLI_FAILS > 0)); then
    echo "FAIL: mem_io_cli: $MEM_CLI_FAILS of $MEM_CLI_CHECKS checks failed"
    exit 1
fi
echo "PASS: mem_io_cli: $MEM_CLI_CHECKS checks, $MEM_CLI_SKIPS skipped"
