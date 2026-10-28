#!/usr/bin/env bash
# test/regression/mem_extend_cli.sh
#
# End-to-end checks of `bwa-mem3 mem` chaining and extension: chaining and
# chain filtering (-w, -W, -N), seed extension through the banded
# Smith-Waterman kernels and the ungapped fast path, z-drop (-d), clipping
# (-L), region patching, the scoring options (-A, -B, -O, -E, -T, -x), the XA
# options (-h, -z, -a) and CIGAR generation. Each check runs under every SIMD
# tier this host can run (test/framework/mem_cli.sh: mem_cli_tiers), on the
# synthetic reference (test/fixtures/mem_synth_ref.awk, whose header lists its
# contigs and planted repeats), with every expected value worked out from the
# reference by hand. Unit tests (test/unit/test_mem_extend_*.cpp) check the
# kernels and filters themselves against oracles; this script checks that the
# options reach them and that the aligner commits what they decide.
#
# Reads (positions 1-based on the synthetic reference; "mm at k" is a
# substitution at read base k by a base differing from the reference there and
# at its neighbours, so no shifted placement scores better):
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
#          deletions, 50M8D28M8D42M, AS 120 - 2*(6 + 8) = 92 (the second
#          deletion moves 2 bases left: 4087-4088 equal 4095-4096, 4086 does
#          not equal 4094; the first cannot: 4050 differs from 4058)
#   r3     ctg1:4301-4350 + ctg1:4357-4380 + ctg1:4387-4400: two 6-base
#          deletions, 50M6D24M6D14M, AS 88 - 2*(6 + 6) = 64 (neither moves:
#          4350 differs from 4356 and 4380 from 4386); the 14-base tail is
#          too short to seed
#   inR    ctg1:2101-2200, inside block R (three exact copies)
#   inS    ctg1:3101-3200, inside block S (ctg2 copy with 5 mismatches here:
#          AS 95 - 20 = 75 at ctg2:2101)
#   l150   ctg2:101-250, mm at 3, 6, 145, 148  seed 7-144 (138), each end
#                                              138 - 8 + 4 = 134 > 138 - 5
#   l250   ctg1:3901-4150, mm at 2, 5, 246, 249  seed 6-245 (240), each end
#                                              240 - 8 + 3 = 235 <= 240 - 5
#   l300   ctg1:1501-1800, mm at 3, 6, 295, 298  seed 7-294 (288), each end 284
# With two mismatches in an extension the ungapped fast path declines it, so
# l150 and l250 (seed + extension < 254) run the 8-bit banded SW kernel and
# l300 (seed 288) the 16-bit one.
#
# Checks:
#   default   c5f 100M AS 97 (95 > 97 - 5: extended, AS stays the best local
#             score); l150 150M AS 138 NM 4; l300 300M AS 288 NM 4; l250
#             5S240M5S POS 3906 AS 240 (235 <= 235: clipped at -L 5 exactly);
#             zd 100M AS 90; zd3 100M AS 85; w10 60M10D40M AS 84.
#   -L        -L 2,3: c5f 3S97M POS 404, c3f 100M, c5r 97M3S;
#             -L 3,2: c5f 100M POS 401, c3f 97M3S, c5r 100M (clip when the
#             end-to-end score <= best - penalty; -L is 5',3' of the read).
#   -O -E     -O 3,9 -E 2,1: del2 AS 100 - (3 + 2*2) = 93, ins2 AS
#             100 - (9 + 2*1) = 89, CIGARs unchanged.
#   -A -B     -A 2 -B 4 (an explicit -B is not scaled): mm AS 198 - 4 = 194,
#             del2 AS 200 - (12 + 2*2) = 184.
#   -T        -T 95 keeps mm (AS 95); -T 96 reports it unmapped (FLAG 4).
#   -d        -d 4: c5f 100M (the 5' extension's drop is exactly 4); -d 3:
#             c5f 3S97M. zd3 at -d 3 stops both extensions at the mismatches:
#             primary 50M50S AS 50 and supplementary 53H47M AS 47. zd at -d 3
#             splits the same way, but the two regions are patched back into
#             100M AS 90 (global score 90 / predicted 100 = 0.90, the patch
#             threshold; zd3's 85 / 100 is under it).
#   -w -W     -w 9: w10's seeds are 10 diagonals apart, out of one chain and out
#             of the band: primary 60M40S AS 60 + supplementary 60H40M AS 40;
#             -w 10: 60M10D40M AS 84. With -W 70, -w 10 still gives 60M10D40M
#             (one chain, weight 100) and -w 9 leaves two chains of 60 and 40,
#             both dropped: unmapped. r2 at -w 8: one chain (diagonals 0, 8,
#             16, each step 8), but the second deletion is 16 off the first
#             seed's diagonal, outside band 8: only the band-doubling retry
#             reaches it, giving the full 50M8D28M8D42M AS 92 with MD
#             50^<4051-4058>28^<4087-4094>42 and NM 16; at -w 7 the three
#             seeds are three chains, so the primary is the first seed's
#             51M69S (read base 51 matches 4051 by chance) AS 51, with two
#             supplementary records. r3 at -w 8: one chain of two seeds
#             (diagonals 0 and 6); the first band reaches the first deletion
#             with max_off 6, not below 8/2 + 8/4 = 6, so the band is doubled,
#             which alone reaches the second deletion (diagonal 12) and the
#             unseeded tail: 50M6D24M6D14M AS 64, MD
#             50^<4351-4356>24^<4381-4386>14, NM 12, the same as at -w 100.
#   -N        -N 1 drops every chain marked kept 1 or 2 (one that overlaps a
#             heavier chain), so each of these reads keeps one: inR MAPQ 60
#             XS 0 no XA (its hit one of the three equal copies), inS XS 0;
#             -N 2 keeps the first such chain: inR XS 100, one XA hit.
#   -h -z -a  inR's two other copies are in XA at -h 2, not at -h 1 (HN 2
#             both); inS's ctg2 hit (75) is in XA at -z 0.75, not at -z 0.76;
#             -a writes it as a secondary record, ctg2:2101 100M NM 5 AS 75,
#             and inR's two other copies as secondaries.
#   -x        intractg (-B 9 -O 16): mm AS 90, del2 AS 100 - 18 = 82, c5f
#             3S97M (97 - 9 + 2 = 90 <= 92); pacbio (-B 1 -O 1 -E 1 -L 0): mm
#             AS 98, del2 AS 97, c5f 100M AS 98.
#   kernels   BWAMEM3_DISABLE_BSW8=1 (every pair on the 16-bit kernel) gives
#             the default SAM byte for byte, on these reads and on 300
#             generated ones (substitutions and indels, both strands).
#             Not a check: a "note:" line prints the banded-SW kernel widths
#             the BWAMEM3_DUMP_PAIRS diagnostic lists for these reads (8 and 16
#             today), to show that the comparison covers both kernels. Which
#             kernel a pair is routed to, and the diagnostic itself, may change
#             without changing the output, so the note never fails the script.
#   threads   -t 1 and -t 4 at -K 10000 give the same SAM on the generated
#             reads.
#   identity  on the generated reads, each path documented as byte-identical
#             gives the default SAM: --keep-contained-ext, --skip-contained-ext,
#             --no-band-cert, --no-adaptive-band, --dedup off, --dedup on,
#             BWA3_CHAIN_FLAT_CAP=0 and =1, BWAMEM3_FORCE_GENMAT=1.
#   tiers     every SAM above is byte-identical across tiers (@PG aside).
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

# xa_hits <record>: the number of hits in its XA tag.
xa_hits() {
    mem_cli_tag "$1" XA | tr ';' '\n' | awk 'NF' | wc -l | tr -d ' '
}

# widths_used <pair-dump>: "8 16" when BWAMEM3_DUMP_PAIRS rows name both
# banded-SW kernel widths ("-" for one that does not occur).
widths_used() {
    awk -F '\t' '
        NR > 1 && ($5 == 8 || $5 == 16) { seen[$5] = 1 }
        END { printf "%s %s", (8 in seen ? "8" : "-"), (16 in seen ? "16" : "-") }' "$1"
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
c5r="$(mem_cli_revcomp "$(mutate "$(seq_at ctg1 801 100)" 98)")"
del2="$(seq_at ctg1 1401 50)$(seq_at ctg1 1453 50)"
ins_left="$(seq_at ctg1 1650 1)"
ins_right="$(seq_at ctg1 1651 1)"
ins2="$(seq_at ctg1 1601 50)$(other_base "$ins_right")$(other_base "$ins_left")$(seq_at ctg1 1651 50)"
zd="$(mutate "$(seq_at ctg1 2501 100)" 51 52)"
zd3="$(mutate "$(seq_at ctg1 2701 100)" 51 52 53)"
w10="$(seq_at ctg1 1301 60)$(seq_at ctg1 1371 40)"
l150="$(mutate "$(seq_at ctg2 101 150)" 3 6 145 148)"
l250="$(mutate "$(seq_at ctg1 3901 250)" 2 5 246 249)"
l300="$(mutate "$(seq_at ctg1 1501 300)" 3 6 295 298)"

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

r3="$(seq_at ctg1 4301 50)$(seq_at ctg1 4357 24)$(seq_at ctg1 4387 14)"
r3_md="50^$(seq_at ctg1 4351 6)24^$(seq_at ctg1 4381 6)14"
if [[ "$(seq_at ctg1 4350 1)" == "$(seq_at ctg1 4356 1)" || "$(seq_at ctg1 4380 1)" == "$(seq_at ctg1 4386 1)" ]]; then
    echo "FAIL: mem_extend_cli: r3's deletions are not placed as the header says on this reference" >&2
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
    mem_cli_fq r3 "$r3"
    mem_cli_fq inR "$(seq_at ctg1 2101 100)"
    mem_cli_fq inS "$(seq_at ctg1 3101 100)"
    mem_cli_fq l150 "$l150"
    mem_cli_fq l250 "$l250"
    mem_cli_fq l300 "$l300"
} > "$work/reads.fq"

# 300 generated reads: 60-400 bases from any contig, either strand, with up to
# four substitutions and, in one read of three, a 1-6 base insertion or
# deletion; every tenth is written twice (g<n> and g<n>d), so identical
# extension jobs occur within a batch. Positions from the Park-Miller generator
# (exact in awk doubles).
awk -v n=300 '
    function rnd(k) { x = (x * 48271) % 2147483647; return x % k }
    function rc(s,   i, o, c) {
        o = ""
        for (i = length(s); i >= 1; i--) {
            c = substr(s, i, 1)
            o = o (c == "A" ? "T" : c == "C" ? "G" : c == "G" ? "C" : c == "T" ? "A" : "N")
        }
        return o
    }
    /^>/ { name[++nc] = substr($1, 2); next }
    { seq[nc] = seq[nc] $0 }
    END {
        x = 4242
        for (r = 1; r <= n; r++) {
            c = 1 + rnd(nc)
            len = 60 + rnd(341)
            if (len > length(seq[c]) - 20) len = length(seq[c]) - 20
            p = 1 + rnd(length(seq[c]) - len - 10)
            s = substr(seq[c], p, len + 10)
            if (rnd(3) == 0) {
                k = 10 + rnd(len - 20)
                g = 1 + rnd(6)
                if (rnd(2)) s = substr(s, 1, k) substr(s, k + g + 1)
                else s = substr(s, 1, k) substr("ACGTTGCAAC", 1 + rnd(4), g) substr(s, k + 1)
            }
            s = substr(s, 1, len)
            for (m = rnd(5); m > 0; m--) {
                k = 1 + rnd(len)
                b = substr("ACGT", 1 + rnd(4), 1)
                s = substr(s, 1, k - 1) b substr(s, k + 1)
            }
            if (rnd(2)) s = rc(s)
            q = s
            gsub(/./, "I", q)
            printf "@g%d\n%s\n+\n%s\n", r, s, q
            if (r % 10 == 0) printf "@g%dd\n%s\n+\n%s\n", r, s, q
        }
    }' "$ref" > "$work/generated.fq"
cat "$work/reads.fq" "$work/generated.fq" > "$work/all.fq"

# Every run: a name and its options, run on reads.fq under each tier.
runs=(
    "default|"
    "b6|-B 6"
    "l23|-L 2,3"
    "l32|-L 3,2"
    "oe|-O 3,9 -E 2,1"
    "a2b4|-A 2 -B 4"
    "t95|-T 95"
    "t96|-T 96"
    "d3|-d 3"
    "d4|-d 4"
    "w7|-w 7"
    "w8|-w 8"
    "w9|-w 9"
    "w10|-w 10"
    "w9W70|-w 9 -W 70"
    "w10W70|-w 10 -W 70"
    "n1|-N 1"
    "n2|-N 2"
    "h1|-h 1"
    "h2|-h 2"
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

    # -B (the placeholder this script started with).
    rec="$(mem_cli_primary "$(s b6)" mm)"
    mem_cli_expect "$tier: -B 6 mm AS CIGAR NM" "93 100M 1" \
        "$(mem_cli_tag "$rec" AS) $(mem_cli_col "$rec" 6) $(mem_cli_tag "$rec" NM)"

    # Defaults.
    expect_primary "$tier: default" "$(s default)" c5f "0 ctg1 401 100M 97"
    expect_primary "$tier: default" "$(s default)" l150 "0 ctg2 101 150M 138"
    expect_primary "$tier: default" "$(s default)" l250 "0 ctg1 3906 5S240M5S 240"
    expect_primary "$tier: default" "$(s default)" l300 "0 ctg1 1501 300M 288"
    nm150="$(mem_cli_tag "$(mem_cli_primary "$(s default)" l150)" NM)"
    nm300="$(mem_cli_tag "$(mem_cli_primary "$(s default)" l300)" NM)"
    mem_cli_expect "$tier: default l150 l300 NM" "4 4" "$nm150 $nm300"
    expect_primary "$tier: default" "$(s default)" zd "0 ctg1 2501 100M 90"
    expect_primary "$tier: default" "$(s default)" zd3 "0 ctg1 2701 100M 85"
    expect_primary "$tier: default" "$(s default)" w10 "0 ctg1 1301 60M10D40M 84"

    # -L: the 5' and 3' penalties, on both strands, at the boundary.
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
    expect_primary "$tier: -A 2 -B 4" "$(s a2b4)" del2 "0 ctg1 1401 50M2D50M 184"

    # -T at the score of the read.
    expect_primary "$tier: -T 95" "$(s t95)" mm "0 ctg1 1001 100M 95"
    mem_cli_expect "$tier: -T 96 mm FLAG RNAME" "4 *" \
        "$(mem_cli_col "$(mem_cli_primary "$(s t96)" mm)" 2) $(mem_cli_col "$(mem_cli_primary "$(s t96)" mm)" 3)"

    # -d: the z-drop boundary, the split it causes, and the patch that undoes it.
    expect_primary "$tier: -d 4" "$(s d4)" c5f "0 ctg1 401 100M 97"
    expect_primary "$tier: -d 3" "$(s d3)" c5f "0 ctg1 404 3S97M 97"
    expect_primary "$tier: -d 3" "$(s d3)" zd3 "0 ctg1 2701 50M50S 50"
    mem_cli_expect "$tier: -d 3 zd3 supplementary" "2048 ctg1 2754 53H47M 47" "$(fields "$(record "$(s d3)" zd3 2048)")"
    expect_primary "$tier: -d 3 (patched)" "$(s d3)" zd "0 ctg1 2501 100M 90"
    # shellcheck disable=SC2016 # an awk condition: its $n are awk fields
    mem_cli_expect "$tier: -d 3 zd records" 1 "$(mem_cli_count "$(s d3)" '$1 == "zd"')"

    # -w and -W: chaining and the band.
    expect_primary "$tier: -w 10" "$(s w10)" w10 "0 ctg1 1301 60M10D40M 84"
    expect_primary "$tier: -w 9" "$(s w9)" w10 "0 ctg1 1301 60M40S 60"
    mem_cli_expect "$tier: -w 9 w10 supplementary" "2048 ctg1 1371 60H40M 40" "$(fields "$(record "$(s w9)" w10 2048)")"
    expect_primary "$tier: -w 10 -W 70" "$(s w10W70)" w10 "0 ctg1 1301 60M10D40M 84"
    mem_cli_expect "$tier: -w 9 -W 70 w10 FLAG RNAME" "4 *" \
        "$(mem_cli_col "$(mem_cli_primary "$(s w9W70)" w10)" 2) $(mem_cli_col "$(mem_cli_primary "$(s w9W70)" w10)" 3)"

    expect_primary "$tier: default" "$(s default)" r2 "0 ctg1 4001 50M8D28M8D42M 92"
    expect_primary "$tier: -w 8 (band retry)" "$(s w8)" r2 "0 ctg1 4001 50M8D28M8D42M 92"
    rec="$(mem_cli_primary "$(s w8)" r2)"
    mem_cli_expect "$tier: -w 8 r2 MD NM" "$r2_md 16" "$(mem_cli_tag "$rec" MD) $(mem_cli_tag "$rec" NM)"
    expect_primary "$tier: default" "$(s default)" r3 "0 ctg1 4301 50M6D24M6D14M 64"
    expect_primary "$tier: -w 8 (band retry)" "$(s w8)" r3 "0 ctg1 4301 50M6D24M6D14M 64"
    rec="$(mem_cli_primary "$(s w8)" r3)"
    mem_cli_expect "$tier: -w 8 r3 MD NM" "$r3_md 12" "$(mem_cli_tag "$rec" MD) $(mem_cli_tag "$rec" NM)"
    expect_primary "$tier: -w 7" "$(s w7)" r2 "0 ctg1 4001 51M69S 51"
    # shellcheck disable=SC2016 # an awk condition: its $n are awk fields
    mem_cli_expect "$tier: -w 7 r2 supplementary records" 2 \
        "$(mem_cli_count "$(s w7)" '$1 == "r2" && int($2 / 2048) % 2 == 1')"

    # -N: chains extended.
    rec="$(mem_cli_primary "$(s n1)" inR)"
    hit="$(mem_cli_col "$rec" 3):$(mem_cli_col "$rec" 4):$(mem_cli_col "$rec" 2)"
    case "$hit" in
        ctg1:2101:0 | ctg2:1101:0 | ctg4:201:16) hit=one-of-R ;;
    esac
    mem_cli_expect "$tier: -N 1 inR hit MAPQ AS XS XA" "one-of-R 60 100 0 " \
        "$hit $(mem_cli_col "$rec" 5) $(mem_cli_tag "$rec" AS) $(mem_cli_tag "$rec" XS) $(mem_cli_tag "$rec" XA)"
    mem_cli_expect "$tier: -N 1 inS XS" 0 "$(mem_cli_tag "$(mem_cli_primary "$(s n1)" inS)" XS)"
    rec="$(mem_cli_primary "$(s n2)" inR)"
    mem_cli_expect "$tier: -N 2 inR MAPQ XS XA-hits" "0 100 1" \
        "$(mem_cli_col "$rec" 5) $(mem_cli_tag "$rec" XS) $(xa_hits "$rec")"
    mem_cli_expect "$tier: default inS XS" 75 "$(mem_cli_tag "$(mem_cli_primary "$(s default)" inS)" XS)"

    # -h, -z, -a.
    rec="$(mem_cli_primary "$(s h2)" inR)"
    mem_cli_expect "$tier: -h 2 inR XA-hits HN" "2 2" \
        "$(xa_hits "$rec") $(mem_cli_tag "$rec" HN)"
    rec="$(mem_cli_primary "$(s h1)" inR)"
    mem_cli_expect "$tier: -h 1 inR XA HN" " 2" "$(mem_cli_tag "$rec" XA) $(mem_cli_tag "$rec" HN)"
    mem_cli_expect "$tier: -z 0.75 inS XA" "ctg2,+2101,100M,5;" "$(mem_cli_tag "$(mem_cli_primary "$(s z75)" inS)" XA)"
    mem_cli_expect "$tier: -z 0.76 inS XA" "" "$(mem_cli_tag "$(mem_cli_primary "$(s z76)" inS)" XA)"
    rec="$(record "$(s xa)" inS 256)"
    mem_cli_expect "$tier: -a inS secondary" "256 ctg2 2101 100M 75 5" "$(fields "$rec") $(mem_cli_tag "$rec" NM)"
    # shellcheck disable=SC2016 # an awk condition: its $n are awk fields
    mem_cli_expect "$tier: -a inR secondaries" 2 "$(mem_cli_count "$(s xa)" '$1 == "inR" && int($2 / 256) % 2 == 1')"
    # shellcheck disable=SC2016 # an awk condition: its $n are awk fields
    mem_cli_expect "$tier: default inS secondaries" 0 \
        "$(mem_cli_count "$(s default)" '$1 == "inS" && int($2 / 256) % 2 == 1')"

    # -x presets.
    expect_primary "$tier: -x intractg" "$(s intractg)" mm "0 ctg1 1001 100M 90"
    expect_primary "$tier: -x intractg" "$(s intractg)" del2 "0 ctg1 1401 50M2D50M 82"
    expect_primary "$tier: -x intractg" "$(s intractg)" c5f "0 ctg1 404 3S97M 97"
    expect_primary "$tier: -x pacbio" "$(s pacbio)" mm "0 ctg1 1001 100M 98"
    expect_primary "$tier: -x pacbio" "$(s pacbio)" del2 "0 ctg1 1401 50M2D50M 97"
    expect_primary "$tier: -x pacbio" "$(s pacbio)" c5f "0 ctg1 401 100M 98"

    # The 8-bit kernel against the 16-bit one, and threads at a fixed -K. The
    # pair dump (BWAMEM3_DUMP_PAIRS, one row per banded-SW pair with the
    # kernel width it was routed to) is reported as a note, never checked: it
    # shows that both kernels ran, but routing is not part of the output.
    mem_cli_mem "$tier" "$work/all.$tier.sam" "$ref" "$work/all.fq"
    widths=unknown
    if BWAMEM3_DUMP_PAIRS="$work/pairs.$tier.tsv" BWAMEM3_FORCE_TIER="$tier" "$BWA_MEM3" mem "$ref" \
        "$work/all.fq" > "$work/dump.$tier.sam" 2> "$work/dump.$tier.sam.err" \
        && [[ -s "$work/pairs.$tier.tsv" ]]; then
        widths="$(widths_used "$work/pairs.$tier.tsv")"
    fi
    echo "note: $tier: banded-SW kernel widths in the pair dump: $widths (8 16 while both kernels take these reads)"
    if ! BWAMEM3_DISABLE_BSW8=1 BWAMEM3_FORCE_TIER="$tier" "$BWA_MEM3" mem "$ref" "$work/all.fq" \
        > "$work/nobsw8.$tier.sam" 2> "$work/nobsw8.$tier.sam.err"; then
        echo "FAIL: mem_extend_cli: BWAMEM3_DISABLE_BSW8=1 run failed (tier $tier)" >&2
        exit 1
    fi
    mem_cli_same_file "$tier: BWAMEM3_DISABLE_BSW8=1 SAM identical to default" \
        "$work/all.$tier.sam" "$work/nobsw8.$tier.sam"
    mem_cli_mem "$tier" "$work/t1.$tier.sam" -t 1 -K 10000 "$ref" "$work/all.fq"
    mem_cli_mem "$tier" "$work/t4.$tier.sam" -t 4 -K 10000 "$ref" "$work/all.fq"
    mem_cli_same_file "$tier: -t 1 and -t 4 at -K 10000 identical" "$work/t1.$tier.sam" "$work/t4.$tier.sam"
    # The paths the code documents as byte-identical to the default, each run
    # with the default's fast path turned off or forced: the contained-seed
    # extension skip, the certified band, extension-job dedup, the flat
    # chaining index (0: every read through the B-tree; 1: a replay mid-read),
    # and the generic-matrix prepass of the banded SW.
    for alt in "keep|--keep-contained-ext" "skipdep|--skip-contained-ext" "nocert|--no-band-cert" \
        "noadapt|--no-adaptive-band" "dedupoff|--dedup off" "dedupon|--dedup on"; do
        # shellcheck disable=SC2086 # the options are word-split on purpose
        mem_cli_mem "$tier" "$work/${alt%%|*}.$tier.sam" ${alt#*|} "$ref" "$work/all.fq"
        mem_cli_same_file "$tier: ${alt#*|} SAM identical to default" "$work/all.$tier.sam" "$work/${alt%%|*}.$tier.sam"
    done
    for env in BWA3_CHAIN_FLAT_CAP=0 BWA3_CHAIN_FLAT_CAP=1 BWAMEM3_FORCE_GENMAT=1; do
        if ! env "$env" BWAMEM3_FORCE_TIER="$tier" "$BWA_MEM3" mem "$ref" "$work/all.fq" \
            > "$work/env.$tier.sam" 2> "$work/env.$tier.sam.err"; then
            echo "FAIL: mem_extend_cli: $env run failed (tier $tier)" >&2
            exit 1
        fi
        mem_cli_same_file "$tier: $env SAM identical to default" "$work/all.$tier.sam" "$work/env.$tier.sam"
    done

    # Most generated reads map, so the comparisons above are not of empty output.
    # shellcheck disable=SC2016 # an awk condition: its $n are awk fields
    mapped="$(mem_cli_count "$work/t1.$tier.sam" \
        '$1 ~ /^g/ && int($2 / 4) % 2 == 0 && int($2 / 256) % 2 == 0 && int($2 / 2048) % 2 == 0')"
    mem_cli_expect "$tier: at least 300 of the 330 generated reads mapped" yes \
        "$( ((mapped >= 300)) && echo yes || echo "no ($mapped)")"

    if [[ -z $first_tier ]]; then
        first_tier="$tier"
    else
        for run in "${runs[@]}" "all|" "t4|"; do
            name="${run%%|*}"
            mem_cli_same_file "$name SAM identical, $first_tier vs $tier" \
                "$work/$name.$first_tier.sam" "$work/$name.$tier.sam"
        done
    fi
done

if ((MEM_CLI_FAILS > 0)); then
    echo "FAIL: mem_extend_cli: $MEM_CLI_FAILS of $MEM_CLI_CHECKS checks failed"
    exit 1
fi
echo "PASS: mem_extend_cli: $MEM_CLI_CHECKS checks, $MEM_CLI_SKIPS skipped"
