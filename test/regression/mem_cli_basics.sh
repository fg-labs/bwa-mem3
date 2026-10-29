#!/usr/bin/env bash
# test/regression/mem_cli_basics.sh
#
# End-to-end `bwa-mem3 mem` records whose every field is worked out by hand,
# on the synthetic reference (test/fixtures/mem_synth.fa; contigs and planted
# repeats are listed in the header of test/tools/mem_synth_ref.c, which
# generates and checks it), under each SIMD tier this host can run
# (test/framework/mem_cli.sh: mem_cli_tiers).
#
# The reads are 100 bp slices of the reference, so every expected value
# follows from the default scoring (-A 1 -B 4 -O 6 -E 1) without running an
# aligner:
#
#   exact     ctg1:1001-1100 forward: FLAG 0, POS 1001, MAPQ 60, 100M, NM 0,
#             MD 100, AS 100 (100 matches), XS 0 (no second hit).
#   rev       reverse complement of ctg1:1201-1300, with a quality string
#             that is not a palindrome: FLAG 16, POS 1201, 100M, SEQ printed
#             on the forward strand (the reference bases) and QUAL reversed
#             with it.
#   mm        exact with read base 50 changed: NM 1, MD <49><ref base><50>,
#             AS 99 - 4 = 95.
#   del2      ctg1:1401-1450 + ctg1:1453-1502 (two reference bases skipped):
#             50M2D50M, NM 2, MD 50^<the two bases>50, AS 100 - (6 + 2*1) = 92.
#             The deletion has one placement only because the bases at its
#             ends differ from their neighbours (1450 vs 1452, 1451 vs 1453);
#             the script checks that before relying on it.
#   inR       ctg1:2101-2200, inside block R, which occurs three times
#             (ctg1:+2001, ctg2:+1001, ctg4:-101 as its reverse complement):
#             MAPQ 0, AS 100, XS 100, and the primary hit plus the XA list are
#             exactly {ctg1 +2101, ctg2 +1101, ctg4 -201}, each 100M with NM 0.
#   inS       ctg1:3101-3200, inside block S, whose ctg2 copy carries a
#             substitution at block offsets 10, 30, ..., 490: POS 3101 on
#             ctg1, AS 100, XS 95 - 5*4 = 75 (the copy has 5 mismatches in
#             offsets 101-200, at 110, 130, 150, 170 and 190).
#   unmapped  ctg1:3701-3800 reversed (not complemented), which is expected to
#             share no 19-mer with the reference and so to leave no seed; the
#             record itself is what is checked: FLAG 4, RNAME *, POS 0,
#             MAPQ 0, CIGAR *, SEQ as given.
#   -A 2      scales -B -O -E -L -T -d -U with it: exact AS 200, mm
#             AS 2*99 - 8 = 190, del2 AS 200 - (12 + 2*2) = 184, each at the
#             default run's POS, CIGAR and NM; rev AS 200.
#
# Each read must give exactly one record (XA carries the alternative hits),
# and the SAM text (@PG aside) must be byte-identical across tiers; on a host
# with one tier, that comparison is reported as SKIP.
#
# The helpers it uses are in test/framework/mem_cli.sh.
#
# Inputs (env vars):
#   BWA_MEM3       path to the bwa-mem3 binary under test
#   MEM_CLI_TIERS  (optional) tiers to run instead of every runnable one
#   MEM_CLI_KEEP   (optional) 1 keeps the work directory
# No staged inputs: the reference is the committed test/fixtures/mem_synth.fa,
# and the reads are sliced from it.

set -euo pipefail

# shellcheck source=test/framework/mem_cli.sh
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../framework" && pwd)/mem_cli.sh"

mem_cli_init mem_cli_basics
mem_cli_setup_synth
ref="$MEM_CLI_SYNTH"
work="$MEM_CLI_WORK"

# Each slice on its own line: a failed mem_cli_seq ends the script only when
# its substitution is the whole assignment.
exact="$(mem_cli_seq "$ref" ctg1 1001 100)"
rev_fwd="$(mem_cli_seq "$ref" ctg1 1201 100)"
mm_ref_base="${exact:49:1}"
case "$mm_ref_base" in
    A) mm_read_base=T ;;
    *) mm_read_base=A ;;
esac
mm="$(mem_cli_subst "$exact" 50 "$mm_read_base")"
del_left="$(mem_cli_seq "$ref" ctg1 1401 50)"
del_right="$(mem_cli_seq "$ref" ctg1 1453 50)"
del2="$del_left$del_right"
del_skipped="$(mem_cli_seq "$ref" ctg1 1451 2)"
del_ends="$(mem_cli_seq "$ref" ctg1 1450 4)"
in_r="$(mem_cli_seq "$ref" ctg1 2101 100)"
in_s="$(mem_cli_seq "$ref" ctg1 3101 100)"
unmapped_fwd="$(mem_cli_seq "$ref" ctg1 3701 100)"
unmapped="$(printf '%s\n' "$unmapped_fwd" | rev)"
# 60 x Q40 then 40 x Q20: reversed, it reads 40 x Q20 then 60 x Q40.
rev_qual="$(printf 'I%.0s' $(seq 60))$(printf '5%.0s' $(seq 40))"
rev_qual_sam="$(printf '%s\n' "$rev_qual" | rev)"

{
    mem_cli_fq exact "$exact"
    mem_cli_fq rev "$(mem_cli_revcomp "$rev_fwd")" "$rev_qual"
    mem_cli_fq mm "$mm"
    mem_cli_fq del2 "$del2"
    mem_cli_fq inR "$in_r"
    mem_cli_fq inS "$in_s"
    mem_cli_fq unmapped "$unmapped"
} > "$work/reads.fq"

# The 2 bp deletion slides (and its CIGAR and MD move) if the base before it
# equals its last base, or its first base equals the base after it.
del_unique=yes
if [[ ${del_ends:0:1} == "${del_ends:2:1}" || ${del_ends:1:1} == "${del_ends:3:1}" ]]; then
    del_unique="no (ctg1:1450-1453 is $del_ends)"
fi
mem_cli_expect "del2: the deletion has one placement in the reference" yes "$del_unique"

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

# hits <record>: the primary hit and the XA entries, as
# "contig,strand+pos,CIGAR,NM", sorted, space-separated.
hits() {
    local strand=+
    if (($(mem_cli_col "$1" 2) / 16 % 2 == 1)); then strand=-; fi
    {
        printf '%s,%s%s,%s,%s;' "$(mem_cli_col "$1" 3)" "$strand" \
            "$(mem_cli_col "$1" 4)" "$(mem_cli_col "$1" 6)" "$(mem_cli_tag "$1" NM)"
        mem_cli_tag "$1" XA
    } | tr ';' '\n' | awk 'NF' | LC_ALL=C sort | tr '\n' ' '
}

# placement <record>: its FLAG, RNAME, POS and CIGAR.
placement() {
    printf '%s %s %s %s\n' "$(mem_cli_col "$1" 2)" "$(mem_cli_col "$1" 3)" "$(mem_cli_col "$1" 4)" "$(mem_cli_col "$1" 6)"
}

first_tier=""
for tier in $(mem_cli_tiers); do
    sam="$work/default.$tier.sam"
    mem_cli_mem "$tier" "$sam" "$ref" "$work/reads.fq"

    mem_cli_expect "$tier: one record per read" "del2 exact inR inS mm rev unmapped " \
        "$(mem_cli_records "$sam" | awk -F '\t' '{ print $1 }' | LC_ALL=C sort | tr '\n' ' ')"

    check_record "$tier" "$sam" exact 0 ctg1 1001 60 100M
    mem_cli_expect "$tier: exact SEQ is the forward-strand reference" "$exact" \
        "$(mem_cli_col "$(mem_cli_primary "$sam" exact)" 10)"
    check_tags "$tier" "$sam" exact 0 100 100 0

    check_record "$tier" "$sam" rev 16 ctg1 1201 60 100M
    check_tags "$tier" "$sam" rev 0 100 100 0
    rec="$(mem_cli_primary "$sam" rev)"
    mem_cli_expect "$tier: rev SEQ is the forward-strand reference" "$rev_fwd" "$(mem_cli_col "$rec" 10)"
    mem_cli_expect "$tier: rev QUAL is reversed with SEQ" "$rev_qual_sam" "$(mem_cli_col "$rec" 11)"

    check_record "$tier" "$sam" mm 0 ctg1 1001 60 100M
    check_tags "$tier" "$sam" mm 1 "49${mm_ref_base}50" 95 0
    mem_cli_expect "$tier: mm SEQ is the read with base 50 changed" "$mm" \
        "$(mem_cli_col "$(mem_cli_primary "$sam" mm)" 10)"

    check_record "$tier" "$sam" del2 0 ctg1 1401 60 50M2D50M
    check_tags "$tier" "$sam" del2 2 "50^${del_skipped}50" 92 0

    rec="$(mem_cli_primary "$sam" inR)"
    mem_cli_expect "$tier: inR MAPQ" 0 "$(mem_cli_col "$rec" 5)"
    mem_cli_expect "$tier: inR AS" 100 "$(mem_cli_tag "$rec" AS)"
    mem_cli_expect "$tier: inR XS" 100 "$(mem_cli_tag "$rec" XS)"
    mem_cli_expect "$tier: inR primary + XA hits" \
        "ctg1,+2101,100M,0 ctg2,+1101,100M,0 ctg4,-201,100M,0 " "$(hits "$rec")"

    check_record "$tier" "$sam" inS 0 ctg1 3101 60 100M
    check_tags "$tier" "$sam" inS 0 100 100 75

    rec="$(mem_cli_primary "$sam" unmapped)"
    mem_cli_expect "$tier: unmapped FLAG RNAME POS MAPQ CIGAR" "4 * 0 0 *" \
        "$(mem_cli_col "$rec" 2) $(mem_cli_col "$rec" 3) $(mem_cli_col "$rec" 4) $(mem_cli_col "$rec" 5) $(mem_cli_col "$rec" 6)"
    mem_cli_expect "$tier: unmapped SEQ" "$unmapped" "$(mem_cli_col "$rec" 10)"

    a2="$work/a2.$tier.sam"
    mem_cli_mem "$tier" "$a2" -A 2 "$ref" "$work/reads.fq"
    mem_cli_expect "$tier: -A 2 one record per read" "del2 exact inR inS mm rev unmapped " \
        "$(mem_cli_records "$a2" | awk -F '\t' '{ print $1 }' | LC_ALL=C sort | tr '\n' ' ')"
    # Scaling the scores moves none of the alignments: POS, CIGAR and NM as
    # under the default scoring, AS scaled.
    for e in "exact 1001 100M 0 200" "mm 1001 100M 1 190" "del2 1401 50M2D50M 2 184"; do
        read -r q want_rest <<< "$e"
        rec="$(mem_cli_primary "$a2" "$q")"
        mem_cli_expect "$tier: -A 2 $q POS CIGAR NM AS" "$want_rest" \
            "$(mem_cli_col "$rec" 4) $(mem_cli_col "$rec" 6) $(mem_cli_tag "$rec" NM) $(mem_cli_tag "$rec" AS)"
    done
    # Every read keeps the default run's FLAG, RNAME, POS and CIGAR (and inR its
    # primary + XA set).
    for q in exact mm del2 rev inR inS unmapped; do
        mem_cli_expect "$tier: -A 2 $q FLAG RNAME POS CIGAR as by default" \
            "$(placement "$(mem_cli_primary "$sam" "$q")")" "$(placement "$(mem_cli_primary "$a2" "$q")")"
    done
    mem_cli_expect "$tier: -A 2 rev AS" 200 "$(mem_cli_tag "$(mem_cli_primary "$a2" rev)" AS)"
    # -A 2 scales the penalties with the match score, so every score doubles:
    # inR's three exact copies score 200 each (AS = XS), and inS's 5%-diverged
    # copy, 95 matches and 5 mismatches, scores 190 - 5 * 8 = 150.
    for e in "inR 200 200" "inS 200 150"; do
        read -r q want_rest <<< "$e"
        rec="$(mem_cli_primary "$a2" "$q")"
        mem_cli_expect "$tier: -A 2 $q AS XS" "$want_rest" "$(mem_cli_tag "$rec" AS) $(mem_cli_tag "$rec" XS)"
    done
    mem_cli_expect "$tier: -A 2 inR primary + XA hits as by default" \
        "$(hits "$(mem_cli_primary "$sam" inR)")" "$(hits "$(mem_cli_primary "$a2" inR)")"

    if [[ -z $first_tier ]]; then
        first_tier="$tier"
    else
        mem_cli_same_file "SAM identical, $first_tier vs $tier" "$work/default.$first_tier.sam" "$sam"
        mem_cli_same_file "-A 2 SAM identical, $first_tier vs $tier" "$work/a2.$first_tier.sam" "$a2"
    fi
done

if [[ $(mem_cli_tiers | wc -l) -lt 2 ]]; then
    mem_cli_skip "SAM identical across tiers" "only one tier ($first_tier) runs on this host"
fi

if mem_cli_verdict; then
    echo "PASS: $MEM_CLI_VERDICT"
else
    echo "FAIL: $MEM_CLI_VERDICT"
    exit 1
fi
