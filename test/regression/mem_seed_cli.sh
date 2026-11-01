#!/usr/bin/env bash
# test/regression/mem_seed_cli.sh
#
# End-to-end checks of `bwa-mem3 mem` seeding options, worked out by hand on
# the synthetic reference (test/fixtures/mem_synth.fa; the header of
# test/tools/mem_synth_ref.c lists its contigs and planted repeats), and of the
# seeding knobs that promise identical output reaching the code they tune. The
# seeding functions themselves are checked against a brute-force suffix array
# in test/unit/test_mem_seed_oracle.cpp; --ks-dedup, BWA3_KS_DEDUP and
# BWA3_BWTSEED_LOCKSTEP are covered by ks_dedup_byte_identity.sh and
# bwtseed_lockstep_driver_identity.sh, re-sa by test/resa_byte_identity_test.sh,
# and -t by thread_determinism.sh.
#
# Checks under each SIMD tier this host can run (mem_cli_tiers), whose SAM of
# the hand reads must also be byte-identical across tiers:
#   -k    k19, ctg1:1101-1159 with read bases 20 and 40 changed: its only exact
#         matches of 19 bp or more are 1-19, 21-39 and 41-59: mapped at POS
#         1101, 59M, NM 2, AS 57 - 2*4 = 49 under the default -k 19 (a seed of
#         exactly k counts), unmapped under -k 20.
#   -r    inS60, ctg1:3001-3060 (block S, offsets 1-60): its ctg2 copy differs
#         at offsets 10, 30 and 50, so the copy is seeded only by re-seeding
#         inside the 60 bp SMEM, at its middle (read base 31, offset 31), where
#         the 19 bp window 31-49 occurs twice. Default -r 1.5 (60 >= 19*1.5):
#         XS 57 - 3*4 = 45; -r 10 (60 < 190): no re-seed, XS 0. AS 60, POS 3001
#         both times.
#   -y    the third round alone finds that copy under -k 15 -r 10: its walk
#         over inS60 emits 16-mers at read bases 1, 17 and 33, and only the
#         third (offsets 33-48) lies inside a window, occurring twice. A seed is
#         emitted once its occurrence count drops below -y: -y 3 gives XS 45,
#         -y 2 (count 2 is not < 2) and -y 0 (round off) give XS 0.
#   -c    Rend, ctg1:2251-2300 (the last 50 bp of block R), occurs three times:
#         ctg1 +2251, ctg2 +1251 and ctg4 -101. -c N keeps the first N rows of
#         its SA interval, i.e. the copies whose following text sorts first;
#         the script sorts the 30 bases after each copy (on its strand) to
#         derive them. -c 1: the one primary at the first copy, XS 0, no XA.
#         -c 2: the primary plus XA are exactly the first two copies.
#
# Checks of index and knob wiring, run once (the tier does not matter):
#   SA stride: indexes built with `index -u 0` and `index -u 6` (SA sampled
#         every row and every 64 rows) give the stride-8 index's SAM for a
#         300-read set under -a (every alignment, so every seed's SA lookup
#         shows), and ctg1:1-60, whose suffix sits at text position 0, maps at
#         POS 1.
#   Knobs documented as output-neutral, against the default run on that set
#         behind 2100 copies of (AC)x20 (each copy occurs 41 times, so the
#         chunked SA resolve crosses a 16384-position chunk boundary):
#         BWA3_SMEM_LOCKSTEP_N=7, BWA3_BWTSEED_LOCKSTEP_N=5, BWA3_SA_LANES=64,
#         BWA3_SA_PER_READ=1 and --huge-pages leave the SAM unchanged. The two
#         lockstep widths are also checked on stderr, where each driver reports
#         the width it took, so the check cannot pass for a variable that is
#         never read; the SA knobs report nothing, so for them only the SAM is
#         checked.
#   MLP probe: on a 150 kb random contig (4688 checkpoint blocks, enough for
#         the probe to measure), BWA3_SMEM_LOCKSTEP_PROBE=1 with
#         BWA3_MLP_DEBUG=1 leaves the SAM unchanged, and stderr shows the probe
#         ran and set the width.
#
# Inputs (env vars):
#   BWA_MEM3       path to the bwa-mem3 binary under test
#   MEM_CLI_TIERS  (optional) tiers to run instead of every runnable one
#   MEM_CLI_KEEP   (optional) 1 keeps the work directory
# No staged inputs: the reference is committed, the reads are sliced from it.

set -euo pipefail

# shellcheck source=test/framework/mem_cli.sh
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../framework" && pwd)/mem_cli.sh"

mem_cli_init mem_seed_cli
mem_cli_setup_synth
ref="$MEM_CLI_SYNTH"
work="$MEM_CLI_WORK"

# next_base <base>: the next base in ACGT order (A>C>G>T>A), as the synthetic
# reference's diverged copy uses.
next_base() {
    case "$1" in
        A) echo C ;;
        C) echo G ;;
        G) echo T ;;
        *) echo A ;;
    esac
}

# changed <seq> <pos>...: <seq> with the base at each 1-based <pos> changed.
changed() {
    local s="$1" p
    shift
    for p in "$@"; do
        s="$(mem_cli_subst "$s" "$p" "$(next_base "${s:$((p - 1)):1}")")"
    done
    printf '%s\n' "$s"
}

# seed_reads <ref.fa> <n> <seed>: n FASTQ reads sliced from the reference with
# a Park-Miller stream (exact in awk's doubles): a contig chosen by length,
# 30-150 bp, 0-3 substitutions, an N in one read of ten, half reverse-
# complemented.
seed_reads() {
    awk -v n="$2" -v x="$3" '
        function rnd(m) {
            x = (x * 48271) % 2147483647
            return x % m
        }
        function rc(s, r, i, c) {
            r = ""
            for (i = length(s); i > 0; i--) {
                c = substr(s, i, 1)
                r = r (c == "A" ? "T" : c == "C" ? "G" : c == "G" ? "C" : c == "T" ? "A" : "N")
            }
            return r
        }
        /^>/ { nc++; next }
        { seq[nc] = seq[nc] $0 }
        END {
            for (c = 1; c <= nc; c++) total += length(seq[c])
            for (i = 1; i <= n; i++) {
                p = rnd(total)
                for (c = 1; p >= length(seq[c]); c++) p -= length(seq[c])
                len = 30 + rnd(121)
                if (len > length(seq[c])) len = length(seq[c])
                s = substr(seq[c], 1 + rnd(length(seq[c]) - len + 1), len)
                nsub = rnd(4)
                for (k = 0; k < nsub; k++) {
                    p = 1 + rnd(len)
                    b = substr(s, p, 1)
                    b = b == "A" ? "C" : b == "C" ? "G" : b == "G" ? "T" : "A"
                    s = substr(s, 1, p - 1) b substr(s, p + 1)
                }
                if (rnd(10) == 0) {
                    p = 1 + rnd(len)
                    s = substr(s, 1, p - 1) "N" substr(s, p + 1)
                }
                if (rnd(2) == 1) s = rc(s)
                q = s
                gsub(/./, "I", q)
                printf "@g%d\n%s\n+\n%s\n", i, s, q
            }
        }' "$1"
}

# mem_env <tier> <out> <"NAME=value ..."> [mem args...]: mem_cli_mem with the
# given environment variables exported (in a subshell, so they do not leak).
mem_env() {
    local tier="$1" out="$2" vars="$3"
    shift 3
    (
        for v in $vars; do export "${v?}"; done
        mem_cli_mem "$tier" "$out" "$@"
    )
}

# --- Reads with hand-derived seeds ----------------------------------------

k19="$(changed "$(mem_cli_seq "$ref" ctg1 1101 59)" 20 40)"
{
    mem_cli_fq k19 "$k19"
    mem_cli_fq inS60 "$(mem_cli_seq "$ref" ctg1 3001 60)"
    mem_cli_fq Rend "$(mem_cli_seq "$ref" ctg1 2251 50)"
} > "$work/hand.fq"

# The order of Rend's three copies in its SA interval: by the text after each
# copy on its own strand (ASCII order is A < C < G < T, the index's order).
rend_order="$(
    {
        printf 'ctg1,+2251 %s\n' "$(mem_cli_seq "$ref" ctg1 2301 30)"
        printf 'ctg2,+1251 %s\n' "$(mem_cli_seq "$ref" ctg2 1301 30)"
        printf 'ctg4,-101 %s\n' "$(dna_revcomp "$(mem_cli_seq "$ref" ctg4 71 30)")"
    } | LC_ALL=C sort -k 2,2 | awk '{ print $1 }'
)"
rend_first="$(printf '%s\n' "$rend_order" | sed -n 1p)"
rend_first_two="$(printf '%s\n' "$rend_order" | sed -n 1,2p | LC_ALL=C sort | tr '\n' ' ')"

# hit <record>: "contig,strand+pos" of a record.
hit() {
    local strand=+
    if (($(mem_cli_col "$1" 2) / 16 % 2 == 1)); then strand=-; fi
    printf '%s,%s%s\n' "$(mem_cli_col "$1" 3)" "$strand" "$(mem_cli_col "$1" 4)"
}

# hits <record>: the primary hit and the XA hits, sorted, space-separated.
hits() {
    {
        hit "$1"
        mem_cli_tag "$1" XA | tr ';' '\n' | awk -F, 'NF { print $1 "," $2 }'
    } | LC_ALL=C sort | tr '\n' ' '
}

# fields <sam> <qname>: FLAG RNAME POS CIGAR AS XS of the primary record.
fields() {
    local rec
    rec="$(mem_cli_primary "$1" "$2")"
    printf '%s %s %s %s %s %s\n' "$(mem_cli_col "$rec" 2)" "$(mem_cli_col "$rec" 3)" \
        "$(mem_cli_col "$rec" 4)" "$(mem_cli_col "$rec" 6)" "$(mem_cli_tag "$rec" AS)" \
        "$(mem_cli_tag "$rec" XS)"
}

# --- The identity read sets -------------------------------------------------

seed_reads "$ref" 300 20261024 > "$work/gen.fq"
{
    mem_cli_fq junction "$(mem_cli_seq "$ref" ctg1 4971 30)$(mem_cli_seq "$ref" ctg2 1 30)"
    mem_cli_fq tail4 "$(mem_cli_seq "$ref" ctg4 541 60)"
    mem_cli_fq start "$(mem_cli_seq "$ref" ctg1 1 60)"
} >> "$work/gen.fq"
awk -v s="$(mem_cli_seq "$ref" ctg3 501 40)" 'BEGIN {
    q = s
    gsub(/./, "I", q)
    for (i = 1; i <= 2100; i++) printf "@ac%d\n%s\n+\n%s\n", i, s, q
}' > "$work/ident.fq"
cat "$work/gen.fq" >> "$work/ident.fq"

# --- Indexes at the finest and coarsest SA stride ----------------------------

# idx_fa <variant>: the reference of index variant <variant> (u0, u6), under a
# name unique to this run (mem_cli_ref_name).
idx_fa() {
    printf '%s/%s/%s\n' "$work" "$1" "$(mem_cli_ref_name "$1")"
}
for u in 0 6; do
    mkdir -p "$work/u$u"
    cp "$ref" "$(idx_fa "u$u")"
    mem_cli_index "$(idx_fa "u$u")" -u "$u"
done

# A 150 kb random contig: 300001 rows, 4688 checkpoint blocks of 64 rows.
mkdir -p "$work/big"
big="$work/big/$(mem_cli_ref_name big)"
awk 'BEGIN {
    x = 777
    print ">big"
    s = ""
    for (i = 0; i < 150000; i++) {
        x = (x * 48271) % 2147483647
        s = s substr("ACGT", int(x / 536870912) + 1, 1)
        if (length(s) == 60) {
            print s
            s = ""
        }
    }
}' > "$big"
mem_cli_index "$big"
for p in 1001 50001 90001; do
    mem_cli_fq "b$p" "$(mem_cli_seq "$big" big "$p" 100)"
done > "$work/big.fq"

first_tier=""
for tier in $(mem_cli_tiers); do
    hand="$work/hand.$tier.sam"
    mem_cli_mem "$tier" "$hand" "$ref" "$work/hand.fq"

    # --- -k ---------------------------------------------------------------
    mem_cli_expect "$tier: default k19 FLAG RNAME POS CIGAR AS XS" "0 ctg1 1101 59M 49 0" "$(fields "$hand" k19)"
    mem_cli_expect "$tier: default k19 NM" 2 "$(mem_cli_tag "$(mem_cli_primary "$hand" k19)" NM)"
    sam="$work/hk20.$tier.sam"
    mem_cli_mem "$tier" "$sam" -k 20 "$ref" "$work/hand.fq"
    mem_cli_expect "$tier: -k 20 k19 FLAG" 4 "$(mem_cli_col "$(mem_cli_primary "$sam" k19)" 2)"

    # --- -r ---------------------------------------------------------------
    mem_cli_expect "$tier: default inS60 FLAG RNAME POS CIGAR AS XS" "0 ctg1 3001 60M 60 45" \
        "$(fields "$hand" inS60)"
    sam="$work/r10.$tier.sam"
    mem_cli_mem "$tier" "$sam" -r 10 "$ref" "$work/hand.fq"
    mem_cli_expect "$tier: -r 10 inS60 FLAG RNAME POS CIGAR AS XS" "0 ctg1 3001 60M 60 0" \
        "$(fields "$sam" inS60)"

    # --- -y ---------------------------------------------------------------
    for y in 0:0 2:0 3:45; do
        sam="$work/y${y%%:*}.$tier.sam"
        mem_cli_mem "$tier" "$sam" -k 15 -r 10 -y "${y%%:*}" "$ref" "$work/hand.fq"
        mem_cli_expect "$tier: -k 15 -r 10 -y ${y%%:*} inS60 POS AS XS" "3001 60 ${y##*:}" \
            "$(fields "$sam" inS60 | awk '{ print $3, $5, $6 }')"
    done

    # --- -c ---------------------------------------------------------------
    sam="$work/c1.$tier.sam"
    mem_cli_mem "$tier" "$sam" -c 1 "$ref" "$work/hand.fq"
    rec="$(mem_cli_primary "$sam" Rend)"
    mem_cli_expect "$tier: -c 1 Rend hit" "$rend_first" "$(hit "$rec")"
    mem_cli_expect "$tier: -c 1 Rend XS XA" "0 " "$(mem_cli_tag "$rec" XS) $(mem_cli_tag "$rec" XA)"
    # shellcheck disable=SC2016 # an awk condition; $1 is awk's field, not the shell's
    mem_cli_expect "$tier: -c 1 Rend records" 1 "$(mem_cli_count "$sam" '$1 == "Rend"')"
    sam="$work/c2.$tier.sam"
    mem_cli_mem "$tier" "$sam" -c 2 "$ref" "$work/hand.fq"
    mem_cli_expect "$tier: -c 2 Rend primary + XA hits" "$rend_first_two" \
        "$(hits "$(mem_cli_primary "$sam" Rend)")"

    if [[ -z $first_tier ]]; then
        first_tier="$tier"
    else
        mem_cli_same_file "hand SAM identical, $first_tier vs $tier" "$work/hand.$first_tier.sam" "$hand"
    fi
done

if [[ $(mem_cli_tiers | wc -l) -lt 2 ]]; then
    mem_cli_skip "SAM identical across tiers" "only one tier ($first_tier) runs on this host"
fi

# ---------------------------------------------------------------------------
# Index and knob wiring, once (under the first tier)
# ---------------------------------------------------------------------------

tier="$first_tier"

# --- SA stride ------------------------------------------------------------
base="$work/gen.sam"
mem_cli_mem "$tier" "$base" -a "$ref" "$work/gen.fq"
mem_cli_expect "stride 8 start FLAG RNAME POS CIGAR" "0 ctg1 1 60M" \
    "$(fields "$base" start | awk '{ print $1, $2, $3, $4 }')"
for idx in u0 u6; do
    sam="$work/gen.$idx.sam"
    mem_cli_mem "$tier" "$sam" -a "$(idx_fa "$idx")" "$work/gen.fq"
    mem_cli_same_file "-a SAM, stride-8 index vs $idx" "$base" "$sam"
done

# --- Knobs that promise identical output ------------------------------------
ident="$work/ident.sam"
mem_cli_mem "$tier" "$ident" "$ref" "$work/ident.fq"
# Each entry: the variables, then (after |) the stderr line that shows the
# setting reached the code it tunes, or "" where the code reports nothing.
for spec in \
    "BWA3_SMEM_LOCKSTEP_N=7|phase-2 SMEM lockstep width: 7 (width source: env pin)" \
    "BWA3_BWTSEED_LOCKSTEP_N=5|width 5 (width source: env pin)" \
    "BWA3_SA_LANES=64|" \
    "BWA3_SA_PER_READ=1|"; do
    e="${spec%%|*}"
    line="${spec#*|}"
    sam="$work/env.sam"
    mem_env "$tier" "$sam" "$e" "$ref" "$work/ident.fq"
    mem_cli_same_file "$e leaves the SAM unchanged" "$ident" "$sam"
    if [[ -n $line ]]; then
        mem_cli_expect "$e reaches the code (stderr)" yes \
            "$(grep -qF -- "$line" "$sam.err" && echo yes || echo no)"
    fi
done
mem_cli_mem "$tier" "$work/hp.sam" --huge-pages "$ref" "$work/ident.fq"
mem_cli_same_file "--huge-pages leaves the SAM unchanged" "$ident" "$work/hp.sam"

# --- MLP probe ---------------------------------------------------------------
# The startup MLP probe needs an index of at least 4096 checkpoint blocks; it
# may pick any width, and the output must not change.
sam="$work/big.sam"
mem_cli_mem "$tier" "$sam" "$big" "$work/big.fq"
mem_env "$tier" "$sam.probe" "BWA3_SMEM_LOCKSTEP_PROBE=1 BWA3_MLP_DEBUG=1" "$big" "$work/big.fq"
mem_cli_same_file "BWA3_SMEM_LOCKSTEP_PROBE=1 leaves the SAM unchanged (150 kb index)" "$sam" "$sam.probe"
mem_cli_expect "the probe ran and set the width (stderr)" yes \
    "$(grep -q '^\[mlp\] ns/access chasing cp_occ:' "$sam.probe.err" && grep -qF '(width source: startup probe)' "$sam.probe.err" && echo yes || echo no)"

if mem_cli_verdict; then
    echo "PASS: $MEM_CLI_VERDICT"
else
    echo "FAIL: $MEM_CLI_VERDICT"
    exit 1
fi
