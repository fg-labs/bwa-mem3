#!/usr/bin/env bash
# test/regression/mem_extend_cli.sh
#
# End-to-end checks of the `bwa-mem3 mem` chaining and extension options that
# no other script reaches: chaining and the band (-w), z-drop (-d) and region
# patching, clipping (-L), the scoring options (-A, -B, -O, -E, -x), chain
# filtering (-N) and the XA options (-h, -z, -a). Each check runs under every
# SIMD tier this host can run (test/framework/mem_cli.sh: mem_cli_tiers), on the
# synthetic reference (test/fixtures/mem_synth.fa; the header of
# test/tools/mem_synth_ref.c lists its contigs and planted repeats), with every
# expected value worked out from the reference by hand. The kernels themselves
# are checked against oracles by unit tests (test/unit/:
# test_mem_extend_bsw_oracle.cpp, test_mem_extend_local_oracle.cpp,
# test_ksw_global2_oracle.cpp); this script checks that the options reach them
# and that the aligner commits what they decide. The byte-identity switches
# (BWAMEM3_DISABLE_BSW8, --keep-contained-ext, --no-band-cert, --dedup, the
# flat chaining index) and -t have scripts of their own (see README.md).
#
# Reads (positions 1-based on the synthetic reference; "mm at k" is a
# substitution at read base k by a base differing from the read there and from
# its two neighbours as already substituted, so no shifted placement scores
# better):
#   mm     ctg1:1001-1100, mm at 50            AS 99 - 4 = 95
#   c5f    ctg1:401-500, mm at 3               seed 4-100 (97); to the 5' end
#                                              97 - 4 + 2 = 95
#   c3f    ctg1:601-700, mm at 98              the same at the 3' end
#   c5r    reverse complement of ctg1:801-900 with mm at forward base 98: the
#          read's 5' end, at the right of the forward-strand CIGAR
#   del2   ctg1:1401-1450 + ctg1:1453-1502     50M2D50M (mem_cli_basics.sh's read)
#   ins2   ctg1:1601-1650 + 2 bases + ctg1:1651-1700, the bases differing from
#          the reference on either side        50M2I50M
#   zd     ctg1:2501-2600, mm at 51 and 52     100M, AS 98 - 8 = 90
#   zd3    ctg1:2701-2800, mm at 51, 52, 53    100M, AS 97 - 12 = 85
#   w10    ctg1:1301-1360 + ctg1:1371-1410     60M10D40M, AS 100 - (6 + 10) = 84
#          (bases 1360/1370 and 1361/1371 differ: one placement)
#   r2     ctg1:4001-4050 + ctg1:4059-4088 + ctg1:4097-4136: two 8-base
#          deletions, 50M8D28M8D42M, AS 120 - 2*(6 + 8) = 92. Each deletion
#          is checked at its leftmost equal-scoring placement: the second
#          moves 2 bases left (4087-4088 equal 4095-4096, 4086 does not equal
#          4094); the first cannot (4050 differs from 4058), though it could
#          equally sit one base right (4051 equals 4059)
#   inR    ctg1:2101-2200, inside block R (three exact copies)
#   inS    ctg1:3101-3200, inside block S (ctg2 copy with 5 mismatches here:
#          AS 95 - 20 = 75 at ctg2:2101)
#   l250   ctg1:3901-4150, mm at 2, 5, 246, 249  seed 6-245 (240), each end
#                                              240 - 8 + 3 = 235 <= 240 - 5
#
# Checks:
#   -B        -B 6: mm AS 99 - 6 = 93, 100M, NM 1.
#   -L        default: l250 5S240M5S POS 3906 AS 240 (235 <= 235: clipped at
#             -L 5 exactly). -L 2,3: c5f 3S97M POS 404, c3f 100M, c5r 97M3S;
#             -L 3,2: c5f 100M POS 401, c3f 97M3S, c5r 100M (clip when the
#             end-to-end score <= best - penalty; -L is 5',3' of the read).
#   -O -E     -O 3,9 -E 2,1: del2 AS 100 - (3 + 2*2) = 93, ins2 AS
#             100 - (9 + 2*1) = 89, CIGARs unchanged.
#   -A -B     -A 2 -B 4 (an explicit -B is not scaled): mm AS 198 - 4 = 194.
#   -d        zd3 at -d 3 stops both extensions at the mismatches: primary
#             50M50S AS 50 and supplementary 53H47M AS 47. zd at -d 3 splits
#             the same way, but the two regions are patched back into one
#             record, 100M AS 90 (global score 90 / predicted 100 = 0.90, the
#             patch threshold; zd3's 85 / 100 is under it).
#   -w        -w 9: w10's seeds are 10 diagonals apart, out of one chain and out
#             of the band: primary 60M40S AS 60 + supplementary 60H40M AS 40;
#             -w 10: 60M10D40M AS 84. r2 at -w 8: one chain (diagonals 0, 8,
#             16, each step 8), but the second deletion is 16 off the first
#             seed's diagonal, outside band 8: only the band-doubling retry
#             reaches it, giving the full 50M8D28M8D42M AS 92 with MD
#             50^<4051-4058>28^<4087-4094>42 and NM 16.
#   -N        -N 1 drops every chain marked kept 1 or 2 (one that overlaps a
#             heavier chain): inR MAPQ 60 XS 0 no XA (its hit one of the three
#             equal copies); -N 2 keeps the first such chain: the primary and
#             the one XA hit are two different copies of R.
#   -h -z -a  inR's two other copies are not in XA at -h 1 (HN 2); inS's ctg2
#             hit (75) is in XA at -z 0.75, not at -z 0.76; -a writes it as a
#             secondary record, ctg2:2101 100M NM 5 AS 75, and inR's two other
#             copies as secondaries (primary and secondaries the three copies).
#   -x        intractg (-B 9 -O 16): mm AS 90, c5f 3S97M (97 - 9 + 2 = 90 <=
#             92); pacbio (-B 1 -O 1 -E 1 -L 0): mm AS 98, c5f 100M AS 98.
#   tiers     every SAM above is byte-identical across tiers (@PG aside).
#
# Inputs (env vars):
#   BWA_MEM3       path to the bwa-mem3 binary under test
#   MEM_CLI_TIERS  (optional) tiers to run instead of every runnable one
#   MEM_CLI_KEEP   (optional) 1 keeps the work directory
# No staged inputs: the reference is committed, the reads are sliced from it.

set -euo pipefail

# shellcheck source=test/framework/mem_cli.sh
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../framework" && pwd)/mem_cli.sh"

mem_cli_init mem_extend_cli
mem_cli_setup_synth
ref="$MEM_CLI_SYNTH"
work="$MEM_CLI_WORK"

# other_base <base>...: the first of A C G T that is none of the arguments.
other_base() {
    local b
    for b in A C G T; do
        case " $* " in
            *" $b "*) ;;
            *)
                echo "$b"
                return 0
                ;;
        esac
    done
}

# mutate <seq> <pos>...: substitute each 1-based position by a base differing
# from the read there and from its two neighbours.
mutate() {
    local s="$1" p
    shift
    for p in "$@"; do
        s="$(mem_cli_subst "$s" "$p" "$(other_base "${s:$((p - 1)):1}" "${s:$((p - 2)):1}" "${s:$p:1}")")"
    done
    printf '%s\n' "$s"
}

# record <sam> <qname> <flag-bit>: the one record of a read with that FLAG bit
# set (256 secondary, 2048 supplementary); empty when not exactly one.
record() {
    awk -F '\t' -v q="$2" -v bit="$3" '
        /^@/ || $1 != q { next }
        int($2 / bit) % 2 == 1 { rec = $0; n++ }
        END { if (n == 1) print rec }' "$1"
}

# fields <record>: "FLAG RNAME POS CIGAR AS" of a record.
fields() {
    printf '%s %s %s %s %s' "$(mem_cli_col "$1" 2)" "$(mem_cli_col "$1" 3)" "$(mem_cli_col "$1" 4)" \
        "$(mem_cli_col "$1" 6)" "$(mem_cli_tag "$1" AS)"
}

# locus <record>: "RNAME:<strand>POS", the form of an XA hit's locus.
locus() {
    local strand=+
    if (($(mem_cli_col "$1" 2) / 16 % 2)); then strand=-; fi
    printf '%s:%s%s' "$(mem_cli_col "$1" 3)" "$strand" "$(mem_cli_col "$1" 4)"
}

# loci_with_xa <record>: the record's locus and each XA hit's, sorted, one
# space apart.
loci_with_xa() {
    {
        locus "$1"
        echo
        mem_cli_tag "$1" XA | tr ';' '\n' | awk -F , 'NF { print $1 ":" $2 }'
    } | LC_ALL=C sort | paste -sd ' ' -
}

# read_loci <sam> <qname> <awk-condition>: the loci of the read's records for
# which the condition holds, sorted, one space apart.
read_loci() {
    awk -F '\t' -v q="$2" '
        !/^@/ && $1 == q && ('"$3"') { print $3 ":" (int($2 / 16) % 2 ? "-" : "+") $4 }' "$1" \
        | LC_ALL=C sort | paste -sd ' ' -
}

# The three exact copies of block R, as loci_with_xa and read_loci list them.
r_copies="ctg1:+2101 ctg2:+1101 ctg4:-201"

# of_r <loci>: "<n> of R", n the number of distinct copies of R among the
# loci, followed by any locus that is not a copy of R and any copy listed again,
# so a duplicate locus shows instead of passing for a second copy.
of_r() {
    tr ' ' '\n' <<< "$1" | LC_ALL=C sort | awk -v r="$r_copies" '
        BEGIN { split(r, c, " "); for (i in c) isr[c[i]] = 1 }
        NF && ($0 in isr) && !seen[$0]++ { n++; next }
        NF { other = other " " $0 }
        END { print (n + 0) " of R" other }'
}

# expect_primary <label> <sam> <qname> <"FLAG RNAME POS CIGAR AS">
expect_primary() {
    mem_cli_expect "$1: $3 FLAG RNAME POS CIGAR AS" "$4" "$(fields "$(mem_cli_primary "$2" "$3")")"
}

seq_at() { mem_cli_seq "$ref" "$@"; }

# The fixed reads (see the header).
r="$(seq_at ctg1 1001 100)"
mm="$(mutate "$r" 50)"
c5f="$(mutate "$(seq_at ctg1 401 100)" 3)"
c3f="$(mutate "$(seq_at ctg1 601 100)" 98)"
c5r="$(dna_revcomp "$(mutate "$(seq_at ctg1 801 100)" 98)")"
del2="$(seq_at ctg1 1401 50)$(seq_at ctg1 1453 50)"
ins_left="$(seq_at ctg1 1650 1)"
ins_right="$(seq_at ctg1 1651 1)"
ins2="$(seq_at ctg1 1601 50)$(other_base "$ins_right")$(other_base "$ins_left")$(seq_at ctg1 1651 50)"
zd="$(mutate "$(seq_at ctg1 2501 100)" 51 52)"
zd3="$(mutate "$(seq_at ctg1 2701 100)" 51 52 53)"
w10="$(seq_at ctg1 1301 60)$(seq_at ctg1 1371 40)"
l250="$(mutate "$(seq_at ctg1 3901 250)" 2 5 246 249)"

r2="$(seq_at ctg1 4001 50)$(seq_at ctg1 4059 30)$(seq_at ctg1 4097 40)"
r2_md="50^$(seq_at ctg1 4051 8)28^$(seq_at ctg1 4087 8)42"
placed=yes
[[ "$(seq_at ctg1 4050 1)" != "$(seq_at ctg1 4058 1)" ]] || placed=no
[[ "$(seq_at ctg1 4087 2)" == "$(seq_at ctg1 4095 2)" ]] || placed=no
[[ "$(seq_at ctg1 4086 1)" != "$(seq_at ctg1 4094 1)" ]] || placed=no
[[ "$(seq_at ctg1 4051 1)" == "$(seq_at ctg1 4059 1)" ]] || placed=no
if [[ $placed != yes ]]; then
    echo "FAIL: mem_extend_cli: r2's deletions are not placed as the header says on this reference" >&2
    exit 1
fi

# The 10-base deletion has one placement only if it cannot shift by one.
if [[ "$(seq_at ctg1 1360 1)" == "$(seq_at ctg1 1370 1)" || "$(seq_at ctg1 1361 1)" == "$(seq_at ctg1 1371 1)" ]]; then
    echo "FAIL: mem_extend_cli: the w10 deletion is not uniquely placed on this reference" >&2
    exit 1
fi

{
    mem_cli_fq mm "$mm"
    mem_cli_fq c5f "$c5f"
    mem_cli_fq c3f "$c3f"
    mem_cli_fq c5r "$c5r"
    mem_cli_fq del2 "$del2"
    mem_cli_fq ins2 "$ins2"
    mem_cli_fq zd "$zd"
    mem_cli_fq zd3 "$zd3"
    mem_cli_fq w10 "$w10"
    mem_cli_fq r2 "$r2"
    mem_cli_fq inR "$(seq_at ctg1 2101 100)"
    mem_cli_fq inS "$(seq_at ctg1 3101 100)"
    mem_cli_fq l250 "$l250"
} > "$work/reads.fq"

# Every run: a name and its options, run on reads.fq under each tier.
runs=(
    "default|"
    "b6|-B 6"
    "l23|-L 2,3"
    "l32|-L 3,2"
    "oe|-O 3,9 -E 2,1"
    "a2b4|-A 2 -B 4"
    "d3|-d 3"
    "w8|-w 8"
    "w9|-w 9"
    "w10|-w 10"
    "n1|-N 1"
    "n2|-N 2"
    "h1|-h 1"
    "z75|-z 0.75"
    "z76|-z 0.76"
    "xa|-a"
    "intractg|-x intractg"
    "pacbio|-x pacbio"
)

# s <run-name>: that run's SAM under the current tier.
s() { printf '%s' "$work/$1.$tier.sam"; }

first_tier=""
for tier in $(mem_cli_tiers); do
    for run in "${runs[@]}"; do
        name="${run%%|*}"
        # shellcheck disable=SC2086 # the options are word-split on purpose
        mem_cli_mem "$tier" "$work/$name.$tier.sam" ${run#*|} "$ref" "$work/reads.fq"
    done

    # -B with the default -A.
    expect_primary "$tier: -B 6" "$(s b6)" mm "0 ctg1 1001 100M 93"
    mem_cli_expect "$tier: -B 6 mm NM" 1 "$(mem_cli_tag "$(mem_cli_primary "$(s b6)" mm)" NM)"

    # -L: the 5' and 3' penalties, on both strands, at the boundary.
    expect_primary "$tier: default" "$(s default)" l250 "0 ctg1 3906 5S240M5S 240"
    expect_primary "$tier: -L 2,3" "$(s l23)" c5f "0 ctg1 404 3S97M 97"
    expect_primary "$tier: -L 2,3" "$(s l23)" c3f "0 ctg1 601 100M 97"
    expect_primary "$tier: -L 2,3" "$(s l23)" c5r "16 ctg1 801 97M3S 97"
    expect_primary "$tier: -L 3,2" "$(s l32)" c5f "0 ctg1 401 100M 97"
    expect_primary "$tier: -L 3,2" "$(s l32)" c3f "0 ctg1 601 97M3S 97"
    expect_primary "$tier: -L 3,2" "$(s l32)" c5r "16 ctg1 801 100M 97"

    # -O -E with deletion and insertion priced apart; -A with an explicit -B.
    expect_primary "$tier: -O 3,9 -E 2,1" "$(s oe)" del2 "0 ctg1 1401 50M2D50M 93"
    expect_primary "$tier: -O 3,9 -E 2,1" "$(s oe)" ins2 "0 ctg1 1601 50M2I50M 89"
    expect_primary "$tier: -A 2 -B 4" "$(s a2b4)" mm "0 ctg1 1001 100M 194"

    # -d: the split z-drop causes, and the patch that undoes it.
    expect_primary "$tier: -d 3" "$(s d3)" zd3 "0 ctg1 2701 50M50S 50"
    mem_cli_expect "$tier: -d 3 zd3 supplementary" "2048 ctg1 2754 53H47M 47" "$(fields "$(record "$(s d3)" zd3 2048)")"
    expect_primary "$tier: -d 3 (patched)" "$(s d3)" zd "0 ctg1 2501 100M 90"
    # shellcheck disable=SC2016 # an awk condition: its $n are awk fields
    mem_cli_expect "$tier: -d 3 zd records" 1 "$(mem_cli_count "$(s d3)" '$1 == "zd"')"

    # -w: chaining, the band, and the band retry.
    expect_primary "$tier: -w 10" "$(s w10)" w10 "0 ctg1 1301 60M10D40M 84"
    expect_primary "$tier: -w 9" "$(s w9)" w10 "0 ctg1 1301 60M40S 60"
    mem_cli_expect "$tier: -w 9 w10 supplementary" "2048 ctg1 1371 60H40M 40" "$(fields "$(record "$(s w9)" w10 2048)")"
    expect_primary "$tier: -w 8 (band retry)" "$(s w8)" r2 "0 ctg1 4001 50M8D28M8D42M 92"
    rec="$(mem_cli_primary "$(s w8)" r2)"
    mem_cli_expect "$tier: -w 8 r2 MD NM" "$r2_md 16" "$(mem_cli_tag "$rec" MD) $(mem_cli_tag "$rec" NM)"

    # -N: chains extended.
    rec="$(mem_cli_primary "$(s n1)" inR)"
    hit="$(mem_cli_col "$rec" 3):$(mem_cli_col "$rec" 4):$(mem_cli_col "$rec" 2)"
    case "$hit" in
        ctg1:2101:0 | ctg2:1101:0 | ctg4:201:16) hit=one-of-R ;;
    esac
    mem_cli_expect "$tier: -N 1 inR hit MAPQ AS XS XA" "one-of-R 60 100 0 " \
        "$hit $(mem_cli_col "$rec" 5) $(mem_cli_tag "$rec" AS) $(mem_cli_tag "$rec" XS) $(mem_cli_tag "$rec" XA)"
    mem_cli_expect "$tier: -N 2 inR primary and XA loci" "2 of R" \
        "$(of_r "$(loci_with_xa "$(mem_cli_primary "$(s n2)" inR)")")"

    # -h, -z, -a.
    rec="$(mem_cli_primary "$(s h1)" inR)"
    mem_cli_expect "$tier: -h 1 inR XA HN" " 2" "$(mem_cli_tag "$rec" XA) $(mem_cli_tag "$rec" HN)"
    mem_cli_expect "$tier: -z 0.75 inS XA" "ctg2,+2101,100M,5;" "$(mem_cli_tag "$(mem_cli_primary "$(s z75)" inS)" XA)"
    # An absent primary also gives an empty XA, so require the one primary first.
    # shellcheck disable=SC2016 # an awk condition: its $n are awk fields
    mem_cli_expect "$tier: -z 0.76 inS exactly one primary" 1 \
        "$(mem_cli_count "$(s z76)" '$1 == "inS" && int($2 / 256) % 2 == 0 && int($2 / 2048) % 2 == 0')"
    mem_cli_expect "$tier: -z 0.76 inS XA" "" "$(mem_cli_tag "$(mem_cli_primary "$(s z76)" inS)" XA)"
    rec="$(record "$(s xa)" inS 256)"
    mem_cli_expect "$tier: -a inS secondary" "256 ctg2 2101 100M 75 5" "$(fields "$rec") $(mem_cli_tag "$rec" NM)"
    # shellcheck disable=SC2016 # an awk condition: its $n are awk fields
    mem_cli_expect "$tier: -a inR primary and secondary loci" "$r_copies" \
        "$(read_loci "$(s xa)" inR 'int($2 / 2048) % 2 == 0')"
    # shellcheck disable=SC2016 # an awk condition: its $n are awk fields
    mem_cli_expect "$tier: -a inR exactly one primary, at a copy of R" "1 of R" \
        "$(of_r "$(read_loci "$(s xa)" inR 'int($2 / 256) % 2 == 0 && int($2 / 2048) % 2 == 0')")"

    # -x presets.
    expect_primary "$tier: -x intractg" "$(s intractg)" mm "0 ctg1 1001 100M 90"
    expect_primary "$tier: -x intractg" "$(s intractg)" c5f "0 ctg1 404 3S97M 97"
    expect_primary "$tier: -x pacbio" "$(s pacbio)" mm "0 ctg1 1001 100M 98"
    expect_primary "$tier: -x pacbio" "$(s pacbio)" c5f "0 ctg1 401 100M 98"

    if [[ -z $first_tier ]]; then
        first_tier="$tier"
    else
        for run in "${runs[@]}"; do
            name="${run%%|*}"
            mem_cli_same_file "$name SAM identical, $first_tier vs $tier" \
                "$work/$name.$first_tier.sam" "$work/$name.$tier.sam"
        done
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
