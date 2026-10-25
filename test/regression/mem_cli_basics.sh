#!/usr/bin/env bash
# test/regression/mem_cli_basics.sh
#
# End-to-end `bwa-mem3 mem` records whose every field is worked out by hand,
# on the synthetic reference (test/fixtures/mem_synth_ref.awk; contigs and
# planted repeats are listed in its header), under each SIMD tier this host can
# run (test/framework/mem_cli.sh: mem_cli_tiers).
#
# The reads are 100 bp slices of the reference, so every expected value
# follows from the default scoring (-A 1 -B 4 -O 6 -E 1) without running an
# aligner:
#
#   exact     ctg1:1001-1100 forward: FLAG 0, POS 1001, MAPQ 60, 100M, NM 0,
#             MD 100, AS 100 (100 matches), XS 0 (no second hit).
#   rev       reverse complement of ctg1:1201-1300: FLAG 16, POS 1201, 100M,
#             and SEQ printed on the forward strand, i.e. the reference bases.
#   mm        exact with read base 50 changed: NM 1, MD <49><ref base><50>,
#             AS 99 - 4 = 95.
#   del2      ctg1:1401-1450 + ctg1:1453-1502 (two reference bases skipped):
#             50M2D50M, NM 2, MD 50^<the two bases>50, AS 100 - (6 + 2*1) = 92.
#   inR       ctg1:2101-2200, inside block R, which occurs three times
#             (ctg1:+2001, ctg2:+1001, ctg4:-101 as its reverse complement):
#             MAPQ 0, AS 100, XS 100, and the primary hit plus the XA list are
#             exactly {ctg1 +2101, ctg2 +1101, ctg4 -201}, each 100M with NM 0.
#   inS       ctg1:3101-3200, inside block S, whose ctg2 copy carries a
#             substitution at block offsets 10, 30, ..., 490: POS 3101 on
#             ctg1, AS 100, XS 95 - 5*4 = 75 (the copy has 5 mismatches in
#             offsets 101-200, at 110, 130, 150, 170 and 190).
#   unmapped  ctg1:3701-3800 reversed (not complemented), which shares no
#             19-mer with the reference: FLAG 4, RNAME *, POS 0, MAPQ 0,
#             CIGAR *, SEQ as given.
#   -A 2      scales -B -O -E -L -T -d -U with it: exact AS 200, mm
#             AS 2*99 - 8 = 190, del2 AS 200 - (12 + 2*2) = 184.
#
# and the SAM text (@PG aside) must be byte-identical across tiers.
#
# The helpers it uses, shared with the other mem_*_cli.sh scripts beside it,
# are in test/framework/mem_cli.sh.
#
# Inputs (env vars):
#   BWA_MEM3       path to the bwa-mem3 binary under test
#   MEM_CLI_TIERS  (optional) tiers to run instead of every runnable one
#   MEM_CLI_KEEP   (optional) 1 keeps the work directory
# No staged inputs: the reference is generated, the reads are sliced from it.

set -euo pipefail

# shellcheck source=test/framework/mem_cli.sh
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../framework" && pwd)/mem_cli.sh"

mem_cli_init mem_cli_basics
mem_cli_setup_synth
ref="$MEM_CLI_SYNTH"
work="$MEM_CLI_WORK"

exact="$(mem_cli_seq "$ref" ctg1 1001 100)"
rev_fwd="$(mem_cli_seq "$ref" ctg1 1201 100)"
mm_ref_base="${exact:49:1}"
case "$mm_ref_base" in
    A) mm_read_base=T ;;
    *) mm_read_base=A ;;
esac
mm="$(mem_cli_subst "$exact" 50 "$mm_read_base")"
del_skipped="$(mem_cli_seq "$ref" ctg1 1451 2)"
del2="$(mem_cli_seq "$ref" ctg1 1401 50)$(mem_cli_seq "$ref" ctg1 1453 50)"
unmapped="$(mem_cli_seq "$ref" ctg1 3701 100 | rev)"

{
    mem_cli_fq exact "$exact"
    mem_cli_fq rev "$(mem_cli_revcomp "$rev_fwd")"
    mem_cli_fq mm "$mm"
    mem_cli_fq del2 "$del2"
    mem_cli_fq inR "$(mem_cli_seq "$ref" ctg1 2101 100)"
    mem_cli_fq inS "$(mem_cli_seq "$ref" ctg1 3101 100)"
    mem_cli_fq unmapped "$unmapped"
} > "$work/reads.fq"

# check_record <label> <sam> <qname> <FLAG> <RNAME> <POS> <MAPQ> <CIGAR>
check_record() {
    local rec
    rec="$(mem_cli_primary "$2" "$3")"
    mem_cli_expect "$1: $3 FLAG" "$4" "$(mem_cli_col "$rec" 2)"
    mem_cli_expect "$1: $3 RNAME" "$5" "$(mem_cli_col "$rec" 3)"
    mem_cli_expect "$1: $3 POS" "$6" "$(mem_cli_col "$rec" 4)"
    mem_cli_expect "$1: $3 MAPQ" "$7" "$(mem_cli_col "$rec" 5)"
    mem_cli_expect "$1: $3 CIGAR" "$8" "$(mem_cli_col "$rec" 6)"
}

# check_tags <label> <sam> <qname> <NM> <MD> <AS> <XS>
check_tags() {
    local rec
    rec="$(mem_cli_primary "$2" "$3")"
    mem_cli_expect "$1: $3 NM" "$4" "$(mem_cli_tag "$rec" NM)"
    mem_cli_expect "$1: $3 MD" "$5" "$(mem_cli_tag "$rec" MD)"
    mem_cli_expect "$1: $3 AS" "$6" "$(mem_cli_tag "$rec" AS)"
    mem_cli_expect "$1: $3 XS" "$7" "$(mem_cli_tag "$rec" XS)"
}

first_tier=""
for tier in $(mem_cli_tiers); do
    sam="$work/default.$tier.sam"
    mem_cli_mem "$tier" "$sam" "$ref" "$work/reads.fq"

    check_record "$tier" "$sam" exact 0 ctg1 1001 60 100M
    check_tags "$tier" "$sam" exact 0 100 100 0

    check_record "$tier" "$sam" rev 16 ctg1 1201 60 100M
    check_tags "$tier" "$sam" rev 0 100 100 0
    mem_cli_expect "$tier: rev SEQ is the forward-strand reference" "$rev_fwd" \
        "$(mem_cli_col "$(mem_cli_primary "$sam" rev)" 10)"

    check_record "$tier" "$sam" mm 0 ctg1 1001 60 100M
    check_tags "$tier" "$sam" mm 1 "49${mm_ref_base}50" 95 0

    check_record "$tier" "$sam" del2 0 ctg1 1401 60 50M2D50M
    check_tags "$tier" "$sam" del2 2 "50^${del_skipped}50" 92 0

    rec="$(mem_cli_primary "$sam" inR)"
    mem_cli_expect "$tier: inR MAPQ" 0 "$(mem_cli_col "$rec" 5)"
    mem_cli_expect "$tier: inR AS" 100 "$(mem_cli_tag "$rec" AS)"
    mem_cli_expect "$tier: inR XS" 100 "$(mem_cli_tag "$rec" XS)"
    # The primary hit and the XA entries, as "contig,strand+pos,CIGAR,NM", sorted.
    strand=+
    if (($(mem_cli_col "$rec" 2) / 16 % 2 == 1)); then strand=-; fi
    hits="$(
        {
            printf '%s,%s%s,%s,%s;' "$(mem_cli_col "$rec" 3)" "$strand" \
                "$(mem_cli_col "$rec" 4)" "$(mem_cli_col "$rec" 6)" "$(mem_cli_tag "$rec" NM)"
            mem_cli_tag "$rec" XA
        } | tr ';' '\n' | awk 'NF' | LC_ALL=C sort | tr '\n' ' '
    )"
    mem_cli_expect "$tier: inR primary + XA hits" \
        "ctg1,+2101,100M,0 ctg2,+1101,100M,0 ctg4,-201,100M,0 " "$hits"

    check_record "$tier" "$sam" inS 0 ctg1 3101 60 100M
    check_tags "$tier" "$sam" inS 0 100 100 75

    rec="$(mem_cli_primary "$sam" unmapped)"
    mem_cli_expect "$tier: unmapped FLAG RNAME POS MAPQ CIGAR" "4 * 0 0 *" \
        "$(mem_cli_col "$rec" 2) $(mem_cli_col "$rec" 3) $(mem_cli_col "$rec" 4) $(mem_cli_col "$rec" 5) $(mem_cli_col "$rec" 6)"
    mem_cli_expect "$tier: unmapped SEQ" "$unmapped" "$(mem_cli_col "$rec" 10)"

    a2="$work/a2.$tier.sam"
    mem_cli_mem "$tier" "$a2" -A 2 "$ref" "$work/reads.fq"
    mem_cli_expect "$tier: -A 2 exact AS" 200 "$(mem_cli_tag "$(mem_cli_primary "$a2" exact)" AS)"
    mem_cli_expect "$tier: -A 2 mm AS" 190 "$(mem_cli_tag "$(mem_cli_primary "$a2" mm)" AS)"
    mem_cli_expect "$tier: -A 2 del2 AS" 184 "$(mem_cli_tag "$(mem_cli_primary "$a2" del2)" AS)"

    if [[ -z $first_tier ]]; then
        first_tier="$tier"
    else
        mem_cli_same_file "SAM identical, $first_tier vs $tier" "$work/default.$first_tier.sam" "$sam"
        mem_cli_same_file "-A 2 SAM identical, $first_tier vs $tier" "$work/a2.$first_tier.sam" "$a2"
    fi
done

if ((MEM_CLI_FAILS > 0)); then
    echo "FAIL: mem_cli_basics: $MEM_CLI_FAILS of $MEM_CLI_CHECKS checks failed"
    exit 1
fi
echo "PASS: mem_cli_basics: $MEM_CLI_CHECKS checks on tiers: $(mem_cli_tiers | tr '\n' ' ')"
