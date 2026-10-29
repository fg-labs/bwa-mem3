#!/usr/bin/env bash
# test/regression/meth_context_cpg.sh
#
# Regression: what `--meth-context cpg` does to the alignment. Under cpg a
# conversion (ref C / read T on the top strand, ref G / read A on the bottom) is
# free only at a CpG; at any other C or G it is an ordinary mismatch, in the
# extension score (AS), in NM and in MD alike.
#
# Fixture: a deterministic two-contig reference (2 kb each, built below by a
# fixed LCG), indexed with `index --meth`, and hand-built 2x100 bp pairs whose
# conversions sit at positions the script picks from the reference itself.
# Every read is scored under --meth=5base, whose default scoring is neutral: a
# freed conversion scores 0 and a mismatch -4, so with 100 bases
# AS = (matches) - 4 * (mismatches).
#
# Asserted here, for read 1 of each pair:
#  1. OT pair, 3 C->T at CpG Cs: AS 97, NM 0, MD lists the three reference Cs,
#     XM has 3 Z (5-base polarity: a converted CpG C is methylated).
#  2. The same pair plus one C->T at a non-CpG C: under cpg AS 92 and NM 1, and
#     MD shows the reference C there; under `--meth-context all` the same read
#     is AS 96 and NM 0 (the documented difference).
#  3. OB pair, 3 G->A (reference frame) at CpG Gs, read 1 on the reverse strand:
#     AS 97, NM 0.
#  4. The reverse-strand twin: pairs 1 and 2 aligned to a reference whose first
#     contig is reverse-complemented put read 1 on the reverse strand. AS, NM,
#     the count of MD mismatch bases and the count of XM Z calls equal the
#     forward alignment's.
#  5. A CpG whose C is read 1's last aligned base (its G is the next reference
#     base, outside the alignment) converted: freed (100M, AS 99, NM 0).
#  6. Conversions the banded extension scores (cases 1-5 are one seed end to
#     end, since conversions are invisible to the converted seed index): an OT
#     read 1 whose first 16 bases carry two CpG and one non-CpG C->T, cut from
#     its seed by a transversion at offset 17 (the one case the SIMD routing
#     under cpg is visible in). Under cpg AS = 100 - 4 - 8 = 88
#     and NM 2, under all AS 92 and NM 1; forward and reverse-strand twin agree.
#  7. Conversions scored by the mate-rescue SW: a mate that cannot seed (an
#     A <-> T transversion every 15 bases leaves no 19-mer) is rescued from its
#     mate's anchor, so its AS is the rescue SW's score. Six transversions, so
#     with c CpG conversions free and p non-CpG ones: AS = (100 - 6 - c - p)
#     - 4 * (6 + p) under cpg, and p is free under all. Read 2 rescued from
#     read 1 (matrix OT: read 2's chemistry flipped by the rescue strand) with
#     two CpG conversions: AS 68, NM 6; with one non-CpG conversion more: cpg
#     AS 63 / NM 7, all AS 67 / NM 6. Read 1 rescued from read 2 (matrix OB),
#     two CpG and one non-CpG conversion: the same 63 / 7 and 67 / 6. All
#     three equal on the reverse-strand twin. Between the two references each
#     hypothesis is rescued once from a forward-strand anchor (window in the
#     forward half) and once from a reverse-strand one (window in the reverse
#     half); an FR pair always rescues its mate reverse-complemented
#     (is_rev = 1). Pairs use -I 300,20 so rescue needs no insert-size
#     estimate.
#  7b. The same three rescues with --rescue-kmer, on both references: the
#     batched rescue narrows the window to a k-mer anchor's diagonal before the
#     staged copy is labelled. Same AS / NM as 7.
#  8. All of the above reads under `--meth=5base --meth-context all` are
#     record-for-record identical to `--meth=5base` alone.
#
# Inputs:
#   BWA_MEM3 — the bwa-mem3 binary under test: a bare command name resolved via
#              PATH, or a path to the binary (a relative one is made absolute
#              before the `cd "$WORK"` below)
set -euo pipefail
# shellcheck source=test/framework/dna.sh
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../framework" && pwd)/dna.sh"
: "${BWA_MEM3:?BWA_MEM3 must be set}"
case "$BWA_MEM3" in
    */*) BWA_MEM3="$(cd "$(dirname "$BWA_MEM3")" && pwd)/$(basename "$BWA_MEM3")" ;;
esac
command -v samtools > /dev/null 2>&1 || {
    echo "SKIP: samtools not on PATH (--meth emits BAM)"
    exit 0
}
command -v mawk > /dev/null 2>&1 || {
    echo "SKIP: mawk not on PATH"
    exit 0
}

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
cd "$WORK"
fail() {
    echo "FAIL: meth_context_cpg: $*" >&2
    exit 1
}

# A fixed 31-bit LCG, so the reference is the same on every host.
lcg_seq() { # lcg_seq SEED LEN
    local x=$1 n=$2 i seq="" acgt=(A C G T)
    for ((i = 0; i < n; i++)); do
        x=$(((x * 1103515245 + 12345) % 2147483648))
        seq+=${acgt[$(((x >> 16) & 3))]}
    done
    printf '%s' "$seq"
}
# set_base STR OFFSET BASE: STR with the base at 0-based OFFSET replaced.
set_base() { printf '%s' "${1:0:$2}$3${1:$(($2 + 1))}"; }

C1=$(lcg_seq 12345 2000)
C2=$(lcg_seq 67890 2000)
printf '>chr1\n%s\n>chr2\n%s\n' "$C1" "$C2" > ref.fa
# The twin reference: chr1 reverse-complemented, chr2 unchanged.
printf '>chr1rc\n%s\n>chr2\n%s\n' "$(dna_revcomp "$C1")" "$C2" > ref_rc.fa

# Offsets (within win[0,100)) of the Cs followed by G inside the window, of
# the Cs not followed by G, and of the Gs preceded by C. The neighbour is read
# from the contig, so the window's last base sees the base after it.
cpg_c() { # cpg_c CONTIG START
    local s=$1 b=$2 k
    for ((k = 0; k < 100; k++)); do
        [ "${s:$((b + k)):2}" = CG ] && printf '%d ' "$k"
    done
}
plain_c() {
    local s=$1 b=$2 k
    for ((k = 0; k < 100; k++)); do
        [ "${s:$((b + k)):1}" = C ] && [ "${s:$((b + k + 1)):1}" != G ] && printf '%d ' "$k"
    done
}
cpg_g() {
    local s=$1 b=$2 k
    for ((k = 0; k < 100; k++)); do
        [ "${s:$((b + k - 1)):2}" = CG ] && printf '%d ' "$k"
    done
}
# The first N offsets of a space-separated list that lie in [lo, hi].
pick() { # pick N LO HI "OFFSETS"
    local n=$1 lo=$2 hi=$3 o offs picked=()
    read -r -a offs <<< "$4"
    for o in "${offs[@]}"; do
        if [ "$o" -ge "$lo" ] && [ "$o" -le "$hi" ] && [ "${#picked[@]}" -lt "$n" ]; then picked+=("$o"); fi
    done
    [ "${#picked[@]}" -eq "$n" ] || fail "fixture: wanted $n offsets in [$lo,$hi], have: $4"
    printf '%s\n' "${picked[*]}"
}

Q=$(printf 'I%.0s' $(seq 1 100))
: > r1.fq
: > r2.fq
add_pair() { # add_pair NAME R1 R2
    printf '@%s\n%s\n+\n%s\n' "$1" "$2" "$Q" >> r1.fq
    printf '@%s\n%s\n+\n%s\n' "$1" "$3" "$Q" >> r2.fq
}

# Pairs 1 and 2: an OT fragment chr1[300, 600); read 1 = chr1[300, 400) with
# C->T at three CpG Cs (pair 2 adds one at a non-CpG C), read 2 = the
# reverse complement of chr1[500, 600), unconverted.
A=300
W=${C1:$A:100}
read -r -a CPG_A <<< "$(pick 3 10 89 "$(cpg_c "$C1" $A)")"
R1A=$W
for o in "${CPG_A[@]}"; do R1A=$(set_base "$R1A" "$o" T); done
NONCPG=""
for o in $(plain_c "$C1" $A); do
    [ "$o" -lt 30 ] || [ "$o" -gt 70 ] && continue
    far=1
    for c in "${CPG_A[@]}"; do [ $((o - c)) -le 2 ] && [ $((c - o)) -le 2 ] && far=0; done
    [ "$far" = 1 ] && {
        NONCPG=$o
        break
    }
done
[ -n "$NONCPG" ] || fail "fixture: no non-CpG C in chr1[$((A + 30)), $((A + 70))]"
R1B=$(set_base "$R1A" "$NONCPG" T)
R2A=$(dna_revcomp "${C1:$((A + 200)):100}")
add_pair ot3 "$R1A" "$R2A"
add_pair ot3plain "$R1B" "$R2A"

# Pair 3: an OB fragment chr2[200, 500); read 1 = reverse complement of
# chr2[400, 500) with G->A (reference frame) at three CpG Gs, read 2 =
# chr2[200, 300), unconverted.
B=400
read -r -a CPG_B <<< "$(pick 3 10 89 "$(cpg_g "$C2" $B)")"
WB=${C2:$B:100}
for o in "${CPG_B[@]}"; do WB=$(set_base "$WB" "$o" A); done
add_pair ob3 "$(dna_revcomp "$WB")" "${C2:200:100}"

# Pair 4: a CpG at read 1's last base: chr2[e - 99, e + 1) with C at e and
# G at e + 1, so the G lies just past the alignment. Read 2 = the reverse
# complement of chr2[e + 101, e + 201).
E=""
for ((k = 1000; k < 1600; k++)); do
    [ "${C2:$k:2}" = CG ] && {
        E=$k
        break
    }
done
[ -n "$E" ] || fail "fixture: no CpG in chr2[1000, 1600)"
R1E=$(set_base "${C2:$((E - 99)):100}" 99 T)
add_pair edge "$R1E" "$(dna_revcomp "${C2:$((E + 101)):100}")"

# Pair 5: an OT fragment chr1[x, x + 300) whose read 1 has two CpG Cs in
# offsets [1, 15], a plain C in [8, 15] and an A or T at offset 17. Read 1
# converts the three Cs and swaps the base at 17 (A <-> T: a mismatch under
# either hypothesis), so the converted seed starts at 18 and the left
# extension scores [0, 18). The plain C sits at 8 or later so that at least
# four matches follow it on the way to the read start: the extension's best
# score is then the end-to-end one, and AS is that score (bwa reports the
# local maximum otherwise). Read 2 = the reverse complement of
# chr1[x + 200, x + 300).
X=""
for ((x = 700; x < 1600; x++)); do
    [ "${C1:$((x + 17)):1}" = A ] || [ "${C1:$((x + 17)):1}" = T ] || continue
    read -r -a cg <<< "$(cpg_c "$C1" $x)"
    n=0
    for o in "${cg[@]}"; do [ "$o" -ge 1 ] && [ "$o" -le 15 ] && n=$((n + 1)); done
    [ "$n" -ge 2 ] || continue
    pc=""
    for o in $(plain_c "$C1" $x); do
        [ "$o" -ge 8 ] && [ "$o" -le 15 ] && {
            pc=$o
            break
        }
    done
    [ -n "$pc" ] && {
        X=$x
        PLAIN_X=$pc
        break
    }
done
[ -n "$X" ] || fail "fixture: no extension window in chr1[700, 1600)"
read -r -a CPG_X <<< "$(pick 2 1 15 "$(cpg_c "$C1" $X)")"
R1X=${C1:$X:100}
for o in "${CPG_X[@]}" "$PLAIN_X"; do R1X=$(set_base "$R1X" "$o" T); done
[ "${R1X:17:1}" = A ] && R1X=$(set_base "$R1X" 17 T) || R1X=$(set_base "$R1X" 17 A)
add_pair ext "$R1X" "$(dna_revcomp "${C1:$((X + 200)):100}")"

# Rescue pairs. mutate_unseedable CONTIG START NCPG NPLAIN: the forward 100-mer
# CONTIG[START, START + 100) with C->T at NCPG CpG Cs in [5, 95] and NPLAIN
# non-CpG Cs in [30, 70], and an A <-> T swap at the A or T nearest each of
# offsets 10, 25, ..., 85 (never next to a conversion), so no 19 bases in a
# row match. Prints the read.
mutate_unseedable() {
    local s=$1 b=$2 ncpg=$3 nplain=$4 w o t d k conv=() cg pl far
    w=${s:$b:100}
    read -r -a cg <<< "$(pick "$ncpg" 5 95 "$(cpg_c "$s" "$b")")"
    conv+=("${cg[@]}")
    if [ "$nplain" -gt 0 ]; then
        pl=""
        for o in $(plain_c "$s" "$b"); do
            [ "$o" -ge 30 ] && [ "$o" -le 70 ] || continue
            far=1
            for k in "${conv[@]}"; do [ $((o - k)) -le 2 ] && [ $((k - o)) -le 2 ] && far=0; done
            [ "$far" = 1 ] && {
                pl=$o
                break
            }
        done
        [ -n "$pl" ] || fail "fixture: no non-CpG C in [30, 70] of $b"
        conv+=("$pl")
    fi
    for o in "${conv[@]}"; do w=$(set_base "$w" "$o" T); done
    for t in 10 25 40 55 70 85; do
        local hit=""
        for d in 0 1 -1 2 -2; do
            o=$((t + d))
            case "${s:$((b + o)):1}" in A | T) ;; *) continue ;; esac
            far=1
            for k in "${conv[@]}"; do [ $((o - k)) -le 1 ] && [ $((k - o)) -le 1 ] && far=0; done
            [ "$far" = 1 ] && {
                hit=$o
                break
            }
        done
        [ -n "$hit" ] || fail "fixture: no A/T near offset $t of $b"
        [ "${s:$((b + hit)):1}" = A ] && w=$(set_base "$w" "$hit" T) || w=$(set_base "$w" "$hit" A)
    done
    printf '%s\n' "$w"
}
# Pair 6 (a, b): an OT fragment chr1[1200, 1500); read 1 = chr1[1200, 1300),
# read 2 = the reverse complement of the unseedable chr1[1400, 1500).
add_pair resc2cpg "${C1:1200:100}" "$(dna_revcomp "$(mutate_unseedable "$C1" 1400 2 0)")"
add_pair resc2plain "${C1:1200:100}" "$(dna_revcomp "$(mutate_unseedable "$C1" 1400 2 1)")"
# Pair 7: an OT fragment chr1[1650, 1950); read 1 = the unseedable
# chr1[1650, 1750), read 2 = the reverse complement of chr1[1850, 1950).
add_pair resc1plain "$(mutate_unseedable "$C1" 1650 2 1)" "$(dna_revcomp "${C1:1850:100}")"

"$BWA_MEM3" index --meth ref.fa > /dev/null 2>&1 || fail "index --meth ref.fa nonzero exit"
"$BWA_MEM3" index --meth ref_rc.fa > /dev/null 2>&1 || fail "index --meth ref_rc.fa nonzero exit"

run() { # run BAM REF ARGS...
    local bam=$1 ref=$2
    shift 2
    "$BWA_MEM3" mem -t 1 -I 300,20 "$@" "$ref" r1.fq r2.fq > "$bam" 2> "$bam.log" \
        || {
            cat "$bam.log" >&2
            fail "mem $* $ref nonzero exit"
        }
    samtools quickcheck "$bam" || fail "mem $* $ref produced an invalid BAM"
}
run cpg.bam ref.fa --meth=5base --meth-context cpg
run all.bam ref.fa --meth=5base --meth-context all
run base.bam ref.fa --meth=5base
run cpg_rc.bam ref_rc.fa --meth=5base --meth-context cpg
run cpg_kmer.bam ref.fa --meth=5base --meth-context cpg --rescue-kmer
run cpg_rc_kmer.bam ref_rc.fa --meth=5base --meth-context cpg --rescue-kmer

# Read 1 of NAME (primary record): rname pos rev cigar AS NM MD #MD-bases #XM-Z.
# READ=128 selects read 2 instead.
r1() { # r1 BAM NAME
    samtools view -F 0x900 "$1" | mawk -v name="$2" -v bit="${READ:-64}" '
    $1 == name && (int($2/bit) % 2) == 1 {
        as = nm = md = xm = "NA"
        for (i = 12; i <= NF; i++) {
            if ($i ~ /^AS:i:/) as = substr($i, 6)
            if ($i ~ /^NM:i:/) nm = substr($i, 6)
            if ($i ~ /^MD:Z:/) md = substr($i, 6)
            if ($i ~ /^XM:Z:/) xm = substr($i, 6)
        }
        t = md; nb = gsub(/[ACGTN]/, "", t)
        t = xm; nz = gsub(/Z/, "", t)
        print $3, $4, (int($2/16) % 2), $6, as, nm, md, nb, nz
        exit
    }'
}
# The MD string of a 100M read 1 whose reference C at each listed offset is
# read as T and every other base matches.
md_of() {
    local prev=-1 o s=""
    for o in "$@"; do
        s+="$((o - prev - 1))C"
        prev=$o
    done
    printf '%s%d' "$s" $((99 - prev))
}

# 1. OT, three CpG conversions.
read -r rn _ rev cig as nm md nb nz < <(r1 cpg.bam ot3)
[ "$rn/$rev/$cig" = "chr1/0/100M" ] || fail "ot3: placed $rn rev=$rev $cig, want chr1 forward 100M"
[ "$as" = 97 ] || fail "ot3: AS $as, want 97 (three CpG conversions free under cpg)"
[ "$nm" = 0 ] || fail "ot3: NM $nm, want 0 (CpG conversions are not edits)"
want_md=$(md_of "${CPG_A[@]}")
[ "$md" = "$want_md" ] || fail "ot3: MD $md, want $want_md (literal: the three reference Cs)"
[ "$nz" = 3 ] || fail "ot3: $nz XM Z calls, want 3"
OT3_CPG="$as $nm $nb $nz"

# 2. The same read plus a non-CpG conversion: a mismatch under cpg only.
read -r rn _ rev cig as nm md nb nz < <(r1 cpg.bam ot3plain)
[ "$rn/$rev/$cig" = "chr1/0/100M" ] || fail "ot3plain: placed $rn rev=$rev $cig, want chr1 forward 100M"
[ "$as" = 92 ] || fail "ot3plain: AS $as, want 92 (the non-CpG C->T costs a mismatch under cpg)"
[ "$nm" = 1 ] || fail "ot3plain: NM $nm, want 1 (the non-CpG C->T is an edit under cpg)"
mapfile -t SORTED_B < <(printf '%s\n' "${CPG_A[@]}" "$NONCPG" | sort -n)
want_md=$(md_of "${SORTED_B[@]}")
[ "$md" = "$want_md" ] || fail "ot3plain: MD $md, want $want_md (the reference C at offset $NONCPG)"
[ "$nz" = 3 ] || fail "ot3plain: $nz XM Z calls, want 3"
OT3PLAIN_CPG="$as $nm $nb $nz"
read -r _ _ _ cig as nm md _ _ < <(r1 all.bam ot3plain)
[ "$cig" = 100M ] || fail "ot3plain (all): CIGAR $cig, want 100M"
[ "$as" = 96 ] || fail "ot3plain (all): AS $as, want 96 (every conversion free under all)"
[ "$nm" = 0 ] || fail "ot3plain (all): NM $nm, want 0 (every conversion free under all)"
[ "$md" = "$want_md" ] || fail "ot3plain (all): MD $md, want $want_md"

# 3. OB, three CpG conversions, read 1 on the reverse strand.
read -r rn _ rev cig as nm _ _ _ < <(r1 cpg.bam ob3)
[ "$rn/$rev/$cig" = "chr2/1/100M" ] || fail "ob3: placed $rn rev=$rev $cig, want chr2 reverse 100M"
[ "$as" = 97 ] || fail "ob3: AS $as, want 97 (three CpG G->A free under cpg)"
[ "$nm" = 0 ] || fail "ob3: NM $nm, want 0"

# 4. The reverse-strand twin of pairs 1 and 2.
read -r rn _ rev cig as nm _ nb nz < <(r1 cpg_rc.bam ot3)
[ "$rn/$rev/$cig" = "chr1rc/1/100M" ] || fail "ot3 (twin): placed $rn rev=$rev $cig, want chr1rc reverse 100M"
[ "$as $nm $nb $nz" = "$OT3_CPG" ] || fail "ot3 (twin): AS NM #MD #Z = $as $nm $nb $nz, want the forward $OT3_CPG"
read -r rn _ rev cig as nm _ nb nz < <(r1 cpg_rc.bam ot3plain)
[ "$rn/$rev/$cig" = "chr1rc/1/100M" ] || fail "ot3plain (twin): placed $rn rev=$rev $cig, want chr1rc reverse 100M"
[ "$as $nm $nb $nz" = "$OT3PLAIN_CPG" ] || fail "ot3plain (twin): AS NM #MD #Z = $as $nm $nb $nz, want the forward $OT3PLAIN_CPG"

# 5. A CpG straddling read 1's last aligned base.
read -r rn pos rev cig as nm _ _ _ < <(r1 cpg.bam edge)
[ "$rn/$pos/$rev/$cig" = "chr2/$((E - 98))/0/100M" ] \
    || fail "edge: placed $rn:$pos rev=$rev $cig, want chr2:$((E - 98)) forward 100M"
[ "$as" = 99 ] || fail "edge: AS $as, want 99 (the last base's C->T is at a CpG)"
[ "$nm" = 0 ] || fail "edge: NM $nm, want 0 (the CpG's G is outside the alignment but still makes it a CpG)"

# 6. Conversions scored by the banded extension.
read -r rn _ rev cig as nm _ nb nz < <(r1 cpg.bam ext)
[ "$rn/$rev/$cig" = "chr1/0/100M" ] || fail "ext: placed $rn rev=$rev $cig, want chr1 forward 100M"
[ "$as" = 88 ] || fail "ext: AS $as, want 88 (extension: two CpG conversions free, the non-CpG one and the transversion -4 each)"
[ "$nm" = 2 ] || fail "ext: NM $nm, want 2 (the non-CpG conversion and the transversion)"
EXT_CPG="$as $nm $nb $nz"
read -r rn _ rev cig as nm _ nb nz < <(r1 cpg_rc.bam ext)
[ "$rn/$rev/$cig" = "chr1rc/1/100M" ] || fail "ext (twin): placed $rn rev=$rev $cig, want chr1rc reverse 100M"
[ "$as $nm $nb $nz" = "$EXT_CPG" ] || fail "ext (twin): AS NM #MD #Z = $as $nm $nb $nz, want the forward $EXT_CPG"
read -r _ _ _ cig as nm _ _ _ < <(r1 all.bam ext)
[ "$cig/$as/$nm" = "100M/92/1" ] || fail "ext (all): $cig AS $as NM $nm, want 100M AS 92 NM 1 (every conversion free under all)"

# 7. Conversions scored by the mate-rescue SW, forward and on the twin.
rescued() { # rescued NAME READ-BIT WANT-CPG WANT-ALL
    local name=$1 bit=$2 want_cpg=$3 want_all=$4 rn rev cig as nm
    read -r rn _ rev cig as nm _ _ _ < <(READ=$bit r1 cpg.bam "$name")
    [ "$cig" = 100M ] || fail "$name: read $((bit / 64)) $rn rev=$rev $cig, want 100M (rescued)"
    [ "$as $nm" = "$want_cpg" ] || fail "$name: read $((bit / 64)) AS NM = $as $nm, want $want_cpg under cpg (rescue SW on the labelled window)"
    read -r rn _ rev cig as nm _ _ _ < <(READ=$bit r1 cpg_rc.bam "$name")
    [ "$cig" = 100M ] || fail "$name (twin): read $((bit / 64)) $rn rev=$rev $cig, want 100M"
    [ "$as $nm" = "$want_cpg" ] || fail "$name (twin): read $((bit / 64)) AS NM = $as $nm, want $want_cpg under cpg"
    read -r _ _ _ cig as nm _ _ _ < <(READ=$bit r1 all.bam "$name")
    [ "$cig $as $nm" = "100M $want_all" ] || fail "$name (all): read $((bit / 64)) $cig AS NM = $as $nm, want 100M $want_all"
}
rescued resc2cpg 128 "68 6" "68 6"
rescued resc2plain 128 "63 7" "67 6"
rescued resc1plain 64 "63 7" "67 6"

# 7b. The same rescues with --rescue-kmer: the batched rescue scans the plain
# window for a k-mer anchor and narrows it to the anchor's diagonal before the
# staged copy is CpG-labelled, so the labels at the narrowed edges must come
# from their real neighbours. Each rescued mate is ungapped and the anchor lies
# on its diagonal, so the narrowed SW finds the same alignment and score.
rescued_kmer() { # rescued_kmer NAME READ-BIT WANT-CPG
    local name=$1 bit=$2 want=$3 bam rn rev cig as nm
    for bam in cpg_kmer.bam cpg_rc_kmer.bam; do
        read -r rn _ rev cig as nm _ _ _ < <(READ=$bit r1 "$bam" "$name")
        [ "$cig $as $nm" = "100M $want" ] \
            || fail "$name ($bam): read $((bit / 64)) $rn rev=$rev $cig AS NM = $as $nm, want 100M $want under cpg with --rescue-kmer"
    done
}
rescued_kmer resc2cpg 128 "68 6"
rescued_kmer resc2plain 128 "63 7"
rescued_kmer resc1plain 64 "63 7"

# 8. `all` is the default, record for record.
cmp -s <(samtools view all.bam) <(samtools view base.bam) \
    || fail "--meth-context all differs from --meth=5base alone"

echo "PASS: meth_context_cpg"
