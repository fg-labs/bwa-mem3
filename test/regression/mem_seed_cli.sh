#!/usr/bin/env bash
# test/regression/mem_seed_cli.sh
#
# End-to-end checks of `bwa-mem3 mem` seeding: SMEM search over the FM index,
# re-seeding (-r, -y), the minimum seed length (-k), SA lookup (-c), and index
# variants given as input (SA sampling stride, a re-sa'd index). Each check
# runs under every SIMD tier this host can run (test/framework/mem_cli.sh:
# mem_cli_tiers), on the synthetic reference (test/fixtures/mem_synth_ref.awk,
# whose header lists its contigs and planted repeats). Expected values are
# worked out from the reference by hand (the arithmetic is below), or are the
# documented promise that an option or environment variable leaves the output
# byte-identical.
#
# Checks:
#   -k    a 50 bp exact slice of ctg1 (501-550) maps at POS 501 with 50M under
#         the default -k 19, and comes out unmapped (FLAG 4) under -k 51: no
#         seed can be longer than the read.
#   -k    k19, ctg1:1101-1159 with read bases 20 and 40 changed: its only exact
#         matches of 19 bp or more are 1-19, 21-39 and 41-59: mapped at POS
#         1101, 59M, NM 2, AS 57 - 2*4 = 49 under the default -k 19 (a seed of
#         exactly k counts), unmapped under -k 20. k20, ctg1:601-660 with bases
#         21 and 41 changed (matches of 20, 19 and 19 bp): mapped at POS 601,
#         60M, AS 58 - 8 = 50 under -k 20, unmapped under -k 21.
#   -r    inS60, ctg1:3001-3060 (block S, offsets 1-60): its ctg2 copy differs
#         at offsets 10, 30 and 50, so the copy is seeded only by re-seeding
#         inside the 60 bp SMEM, at its middle (read base 31, offset 31), where
#         the 19 bp window 31-49 occurs twice. Default -r 1.5 (60 >= 19*1.5):
#         XS 57 - 3*4 = 45; -r 10 (60 < 190): no re-seed, XS 0. AS 60, POS 3001
#         both times.
#   -y    the third round alone finds that copy under -k 15 -r 10: its walk
#         over inS60 emits 16-mers at read bases 1, 17 and 33, and only the
#         third (offsets 33-48) lies inside a window, occurring twice. A seed is
#         emitted once its occurrence count drops below -y: -y 3 and -y 20 give
#         XS 45, -y 2 (count 2 is not < 2) and -y 0 (round off) give XS 0.
#   -c    Rend, ctg1:2251-2300 (the last 50 bp of block R), occurs three times:
#         ctg1 +2251, ctg2 +1251 and ctg4 -101. -c N keeps the first N rows of
#         its SA interval, i.e. the copies whose following text sorts first;
#         the script sorts the 30 bases after each copy (on its strand) to
#         derive them. -c 1: the one primary at the first copy, XS 0, no XA.
#         -c 2: the primary plus XA are exactly the first two copies. Both
#         under the default --ks-dedup and under off (two staging paths).
#   SA stride: indexes built with `index -u 0..6` (SA sampled every 1 to 64
#         rows) and two re-sa'd copies of the stride-8 index (`re-sa -u 1`,
#         densified; `re-sa -u 6`, coarsened) give SAM identical to the
#         stride-8 index for a 300-read set under -a (every alignment, so every
#         seed's SA lookup shows). ctg1:1-60, whose suffix sits at text
#         position 0 (its SA walk can reach the sentinel row), maps at POS 1.
#   Byte identity, documented as output-neutral, against the default run on
#         that set behind 2100 copies of (AC)x20. Each copy occurs 41 times,
#         and a pass of more than two work items is seeded in full items of
#         512 reads (1024 on arm64), so each seeding call over the copies
#         resolves more SA positions than one 16384-position chunk and the
#         chunked resolve crosses a chunk boundary. Checked:
#         BWA3_SMEM_LOCKSTEP_N 1/7/64, BWA3_SMEM_LOCKSTEP_PROBE=1 with
#         BWA3_MLP_DEBUG=1,
#         BWA3_BWTSEED_LOCKSTEP 0/1, BWA3_BWTSEED_LOCKSTEP_N 1/5/64,
#         BWA3_SA_LANES 1/64, BWA3_SA_PER_READ 1/0, BWA3_LOAD_THREADS 1/4,
#         BWA3_KS_DEDUP_Z=0.5 with BWA3_KS_DEDUP_REPROBE=1, --ks-dedup
#         off/on/auto, --huge-pages, and BWA3_KS_DEDUP off/on with
#         BWA3_KS_DEDUP_STATS=1 (under on the stats line shows distinct_pos <
#         total_pos, since the copies share their SA intervals; under off
#         nothing is collapsed). Each width pin and switch is reported on
#         stderr as taken from the environment. Invalid values
#         (BWA3_SMEM_LOCKSTEP_N=abc, BWA3_BWTSEED_LOCKSTEP=2, BWA3_SA_LANES=0)
#         are diagnosed and resolved as unset, so the SAM is unchanged too.
#   Fatal values: --ks-dedup maybe, BWA3_KS_DEDUP=maybe, BWA3_KS_DEDUP_Z=-1
#         and BWA3_KS_DEDUP_REPROBE=x each exit 1.
#   MLP probe: on a 150 kb random contig (4688 checkpoint blocks, enough for
#         the probe to measure), BWA3_SMEM_LOCKSTEP_PROBE=1 with
#         BWA3_MLP_DEBUG=1 prints the sweep, takes the probe's width, and
#         leaves the SAM unchanged; three 100 bp slices map where they were
#         taken (POS 1001, 50001, 90001).
#   Threads: -t 1 and -t 4 at -K 20000 (several batches) give identical SAM.
#   --smem-dedup, --seed-order global-longest and local-longest leave the
#         records of k19, inS60 and the stride read as derived above (their
#         seeds all fall in one chain, so neither dedup nor order matters).
#   Tiers: the SAM of the hand-read, -a and identity-set runs is identical
#         across tiers.
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

mem_cli_fq r50 "$(mem_cli_seq "$ref" ctg1 501 50)" > "$work/reads.fq"
k19="$(changed "$(mem_cli_seq "$ref" ctg1 1101 59)" 20 40)"
k20="$(changed "$(mem_cli_seq "$ref" ctg1 601 60)" 21 41)"
{
    mem_cli_fq k19 "$k19"
    mem_cli_fq k20 "$k20"
    mem_cli_fq inS60 "$(mem_cli_seq "$ref" ctg1 3001 60)"
    mem_cli_fq Rend "$(mem_cli_seq "$ref" ctg1 2251 50)"
    mem_cli_fq start "$(mem_cli_seq "$ref" ctg1 1 60)"
} > "$work/hand.fq"

# The order of Rend's three copies in its SA interval: by the text after each
# copy on its own strand (ASCII order is A < C < G < T, the index's order).
rend_order="$(
    {
        printf 'ctg1,+2251 %s\n' "$(mem_cli_seq "$ref" ctg1 2301 30)"
        printf 'ctg2,+1251 %s\n' "$(mem_cli_seq "$ref" ctg2 1301 30)"
        printf 'ctg4,-101 %s\n' "$(mem_cli_revcomp "$(mem_cli_seq "$ref" ctg4 71 30)")"
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

# --- Indexes at every SA stride, and re-sa'd copies --------------------------

copy_index() {
    mkdir -p "$2"
    local ext
    for ext in "" .bwt.2bit.64 .ann .amb .pac; do cp "$1$ext" "$2/ref.fa$ext"; done
}
for u in 0 1 2 4 5 6; do
    mkdir -p "$work/u$u"
    cp "$ref" "$work/u$u/ref.fa"
    if ! "$BWA_MEM3" index -u "$u" "$work/u$u/ref.fa" > "$work/u$u/index.log" 2>&1; then
        echo "FAIL: mem_seed_cli: bwa-mem3 index -u $u failed:" >&2
        tail -n 20 "$work/u$u/index.log" >&2
        exit 1
    fi
done
for u in 1 6; do
    copy_index "$ref" "$work/resa$u"
    if ! "$BWA_MEM3" re-sa -u "$u" "$work/resa$u/ref.fa" > "$work/resa$u/resa.log" 2>&1; then
        echo "FAIL: mem_seed_cli: bwa-mem3 re-sa -u $u failed:" >&2
        tail -n 20 "$work/resa$u/resa.log" >&2
        exit 1
    fi
done

# A 150 kb random contig: 300001 rows, 4688 checkpoint blocks of 64 rows.
mkdir -p "$work/big"
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
}' > "$work/big/big.fa"
mem_cli_index "$work/big/big.fa"
for p in 1001 50001 90001; do
    mem_cli_fq "b$p" "$(mem_cli_seq "$work/big/big.fa" big "$p" 100)"
done > "$work/big.fq"

first_tier=""
for tier in $(mem_cli_tiers); do
    # --- -k ---------------------------------------------------------------
    sam="$work/k19.$tier.sam"
    mem_cli_mem "$tier" "$sam" "$ref" "$work/reads.fq"
    rec="$(mem_cli_primary "$sam" r50)"
    mem_cli_expect "$tier: -k 19 r50 FLAG RNAME POS CIGAR" "0 ctg1 501 50M" \
        "$(mem_cli_col "$rec" 2) $(mem_cli_col "$rec" 3) $(mem_cli_col "$rec" 4) $(mem_cli_col "$rec" 6)"

    sam="$work/k51.$tier.sam"
    mem_cli_mem "$tier" "$sam" -k 51 "$ref" "$work/reads.fq"
    mem_cli_expect "$tier: -k 51 r50 FLAG" 4 "$(mem_cli_col "$(mem_cli_primary "$sam" r50)" 2)"

    hand="$work/hand.$tier.sam"
    mem_cli_mem "$tier" "$hand" "$ref" "$work/hand.fq"
    mem_cli_expect "$tier: default k19 FLAG RNAME POS CIGAR AS XS" "0 ctg1 1101 59M 49 0" "$(fields "$hand" k19)"
    mem_cli_expect "$tier: default k19 NM" 2 "$(mem_cli_tag "$(mem_cli_primary "$hand" k19)" NM)"
    sam="$work/hk20.$tier.sam"
    mem_cli_mem "$tier" "$sam" -k 20 "$ref" "$work/hand.fq"
    mem_cli_expect "$tier: -k 20 k19 FLAG" 4 "$(mem_cli_col "$(mem_cli_primary "$sam" k19)" 2)"
    mem_cli_expect "$tier: -k 20 k20 FLAG RNAME POS CIGAR AS XS" "0 ctg1 601 60M 50 0" "$(fields "$sam" k20)"
    sam="$work/hk21.$tier.sam"
    mem_cli_mem "$tier" "$sam" -k 21 "$ref" "$work/hand.fq"
    mem_cli_expect "$tier: -k 21 k20 FLAG" 4 "$(mem_cli_col "$(mem_cli_primary "$sam" k20)" 2)"

    # --- -r ---------------------------------------------------------------
    mem_cli_expect "$tier: default inS60 FLAG RNAME POS CIGAR AS XS" "0 ctg1 3001 60M 60 45" \
        "$(fields "$hand" inS60)"
    sam="$work/r10.$tier.sam"
    mem_cli_mem "$tier" "$sam" -r 10 "$ref" "$work/hand.fq"
    mem_cli_expect "$tier: -r 10 inS60 FLAG RNAME POS CIGAR AS XS" "0 ctg1 3001 60M 60 0" \
        "$(fields "$sam" inS60)"

    # --- -y ---------------------------------------------------------------
    for y in 0:0 2:0 3:45 20:45; do
        sam="$work/y${y%%:*}.$tier.sam"
        mem_cli_mem "$tier" "$sam" -k 15 -r 10 -y "${y%%:*}" "$ref" "$work/hand.fq"
        mem_cli_expect "$tier: -k 15 -r 10 -y ${y%%:*} inS60 POS AS XS" "3001 60 ${y##*:}" \
            "$(fields "$sam" inS60 | awk '{ print $3, $5, $6 }')"
    done

    # --- -c ---------------------------------------------------------------
    # Under the default --ks-dedup (auto) and under off, which stage the SA
    # positions on separate paths.
    for kd in auto off; do
        sam="$work/c1.$kd.$tier.sam"
        mem_cli_mem "$tier" "$sam" -c 1 --ks-dedup "$kd" "$ref" "$work/hand.fq"
        rec="$(mem_cli_primary "$sam" Rend)"
        mem_cli_expect "$tier: -c 1 --ks-dedup $kd Rend hit" "$rend_first" "$(hit "$rec")"
        mem_cli_expect "$tier: -c 1 --ks-dedup $kd Rend XS XA" "0 " "$(mem_cli_tag "$rec" XS) $(mem_cli_tag "$rec" XA)"
        # shellcheck disable=SC2016 # an awk condition; $1 is awk's field, not the shell's
        mem_cli_expect "$tier: -c 1 --ks-dedup $kd Rend records" 1 "$(mem_cli_count "$sam" '$1 == "Rend"')"
        sam="$work/c2.$kd.$tier.sam"
        mem_cli_mem "$tier" "$sam" -c 2 --ks-dedup "$kd" "$ref" "$work/hand.fq"
        mem_cli_expect "$tier: -c 2 --ks-dedup $kd Rend primary + XA hits" "$rend_first_two" \
            "$(hits "$(mem_cli_primary "$sam" Rend)")"
    done

    # --- --smem-dedup and --seed-order ------------------------------------
    for o in --smem-dedup "--seed-order global-longest" "--seed-order local-longest"; do
        sam="$work/opt.$tier.sam"
        # shellcheck disable=SC2086 # $o is one option, or an option and its value
        mem_cli_mem "$tier" "$sam" $o "$ref" "$work/hand.fq"
        mem_cli_expect "$tier: $o k19" "0 ctg1 1101 59M 49 0" "$(fields "$sam" k19)"
        mem_cli_expect "$tier: $o inS60" "0 ctg1 3001 60M 60 45" "$(fields "$sam" inS60)"
        mem_cli_expect "$tier: $o start" "0 ctg1 1 60M 60 0" "$(fields "$sam" start)"
    done

    # --- SA stride and re-sa'd indexes ------------------------------------
    base="$work/gen.$tier.sam"
    mem_cli_mem "$tier" "$base" -a "$ref" "$work/gen.fq"
    mem_cli_expect "$tier: stride 8 start FLAG RNAME POS CIGAR" "0 ctg1 1 60M" \
        "$(fields "$base" start | awk '{ print $1, $2, $3, $4 }')"
    for idx in u0 u1 u2 u4 u5 u6 resa1 resa6; do
        sam="$work/gen.$idx.$tier.sam"
        mem_cli_mem "$tier" "$sam" -a "$work/$idx/ref.fa" "$work/gen.fq"
        mem_cli_same_file "$tier: -a SAM, stride-8 index vs $idx" "$base" "$sam"
        mem_cli_expect "$tier: $idx start FLAG RNAME POS CIGAR" "0 ctg1 1 60M" \
            "$(fields "$sam" start | awk '{ print $1, $2, $3, $4 }')"
    done

    # --- Environment variables and options that promise identical output --
    ident="$work/ident.$tier.sam"
    mem_cli_mem "$tier" "$ident" "$ref" "$work/ident.fq"
    # Each entry: the variables, then (after |) a line the run's stderr must
    # hold, showing the setting reached the code it tunes ("" when the code
    # reports nothing). An invalid value is diagnosed and resolved as unset.
    for spec in \
        "BWA3_SMEM_LOCKSTEP_N=1|phase-2 SMEM lockstep width: 1 (width source: env pin)" \
        "BWA3_SMEM_LOCKSTEP_N=7|phase-2 SMEM lockstep width: 7 (width source: env pin)" \
        "BWA3_SMEM_LOCKSTEP_N=64|phase-2 SMEM lockstep width: 64 (width source: env pin)" \
        "BWA3_SMEM_LOCKSTEP_N=abc|ERROR: BWA3_SMEM_LOCKSTEP_N=\"abc\" is not a positive integer" \
        "BWA3_SMEM_LOCKSTEP_PROBE=1 BWA3_MLP_DEBUG=1|" \
        "BWA3_BWTSEED_LOCKSTEP=0|third-pass bwtseed lockstep: off (BWA3_BWTSEED_LOCKSTEP)" \
        "BWA3_BWTSEED_LOCKSTEP=1|third-pass bwtseed lockstep: on (BWA3_BWTSEED_LOCKSTEP)" \
        "BWA3_BWTSEED_LOCKSTEP=2|ERROR: BWA3_BWTSEED_LOCKSTEP=\"2\" is not 0 or 1" \
        "BWA3_BWTSEED_LOCKSTEP_N=1|width 1 (width source: env pin)" \
        "BWA3_BWTSEED_LOCKSTEP_N=5|width 5 (width source: env pin)" \
        "BWA3_BWTSEED_LOCKSTEP_N=64|width 64 (width source: env pin)" \
        "BWA3_SA_LANES=1|" \
        "BWA3_SA_LANES=64|" \
        "BWA3_SA_LANES=0|ignoring BWA3_SA_LANES=0 (expected 1..64)" \
        "BWA3_SA_PER_READ=1|" \
        "BWA3_SA_PER_READ=0|" \
        "BWA3_LOAD_THREADS=1|" \
        "BWA3_LOAD_THREADS=4|" \
        "BWA3_KS_DEDUP=auto BWA3_KS_DEDUP_Z=0.5 BWA3_KS_DEDUP_REPROBE=1|"; do
        e="${spec%%|*}"
        line="${spec#*|}"
        sam="$work/env.$tier.sam"
        mem_env "$tier" "$sam" "$e" "$ref" "$work/ident.fq"
        mem_cli_same_file "$tier: $e leaves the SAM unchanged" "$ident" "$sam"
        if [[ -n $line ]]; then
            mem_cli_expect "$tier: $e reported" yes "$(grep -qF "$line" "$sam.err" && echo yes || echo no)"
        fi
    done
    for o in "--ks-dedup off" "--ks-dedup on" "--ks-dedup auto" --huge-pages; do
        sam="$work/opt.$tier.sam"
        # shellcheck disable=SC2086 # $o is one option, or an option and its value
        mem_cli_mem "$tier" "$sam" $o "$ref" "$work/ident.fq"
        mem_cli_same_file "$tier: $o leaves the SAM unchanged" "$ident" "$sam"
    done
    for m in on off; do
        sam="$work/ks$m.$tier.sam"
        mem_env "$tier" "$sam" "BWA3_KS_DEDUP=$m BWA3_KS_DEDUP_STATS=1" "$ref" "$work/ident.fq"
        mem_cli_same_file "$tier: BWA3_KS_DEDUP=$m leaves the SAM unchanged" "$ident" "$sam"
        stats="$(sed -n 's/^\[ks-dedup-stats\] total_pos=\([0-9]*\) distinct_pos=\([0-9]*\).*/\1 \2/p' "$sam.err")"
        if [[ $m == on ]]; then
            mem_cli_expect "$tier: BWA3_KS_DEDUP=on collapsed duplicate intervals (distinct < total)" yes \
                "$(printf '%s\n' "$stats" | awk 'NF == 2 && $2 < $1 { print "yes"; exit } { print "no: " $0 }')"
        else
            mem_cli_expect "$tier: BWA3_KS_DEDUP=off collapsed nothing (distinct == total)" yes \
                "$(printf '%s\n' "$stats" | awk 'NF == 2 && $2 == $1 { print "yes"; exit } { print "no: " $0 }')"
        fi
    done

    # A bad --ks-dedup value is fatal (exit 1), from the option or the environment.
    for spec in "--ks-dedup maybe|" "|BWA3_KS_DEDUP=maybe" "|BWA3_KS_DEDUP_Z=-1" "|BWA3_KS_DEDUP_REPROBE=x"; do
        o="${spec%%|*}"
        e="${spec#*|}"
        sam="$work/bad.$tier.sam"
        (
            if [[ -n $e ]]; then export "${e?}"; fi
            # shellcheck disable=SC2086 # $o is empty, or an option and its value
            mem_cli_mem_status "$tier" "$sam" $o "$ref" "$work/hand.fq"
            echo "$MEM_CLI_STATUS" > "$sam.status"
        )
        mem_cli_expect "$tier: ${o:-$e} exits 1" 1 "$(cat "$sam.status")"
    done

    # The startup MLP probe needs an index of at least 4096 checkpoint blocks;
    # it may pick any width, and the output must not change.
    sam="$work/big.$tier.sam"
    mem_cli_mem "$tier" "$sam" "$work/big/big.fa" "$work/big.fq"
    mem_env "$tier" "$sam.probe" "BWA3_SMEM_LOCKSTEP_PROBE=1 BWA3_MLP_DEBUG=1" "$work/big/big.fa" "$work/big.fq"
    mem_cli_same_file "$tier: BWA3_SMEM_LOCKSTEP_PROBE=1 leaves the SAM unchanged (150 kb index)" "$sam" "$sam.probe"
    mem_cli_expect "$tier: the probe ran and BWA3_MLP_DEBUG printed its sweep" "yes yes" \
        "$(grep -q '^\[mlp\] ns/access chasing cp_occ:' "$sam.probe.err" && echo yes || echo no) $(grep -qF '(width source: startup probe)' "$sam.probe.err" && echo yes || echo no)"
    mem_cli_expect "$tier: big-index reads map where they were sliced" "1001 50001 90001" \
        "$(for q in b1001 b50001 b90001; do mem_cli_col "$(mem_cli_primary "$sam" "$q")" 4; done | tr '\n' ' ' | sed 's/ $//')"

    # --- Threads at a fixed -K --------------------------------------------
    t1="$work/t1.$tier.sam"
    t4="$work/t4.$tier.sam"
    mem_cli_mem "$tier" "$t1" -t 1 -K 20000 "$ref" "$work/ident.fq"
    mem_cli_mem "$tier" "$t4" -t 4 -K 20000 "$ref" "$work/ident.fq"
    mem_cli_same_file "$tier: -t 1 vs -t 4 at -K 20000" "$t1" "$t4"

    if [[ -z $first_tier ]]; then
        first_tier="$tier"
    else
        mem_cli_same_file "hand SAM identical, $first_tier vs $tier" "$work/hand.$first_tier.sam" "$hand"
        mem_cli_same_file "-a SAM identical, $first_tier vs $tier" "$work/gen.$first_tier.sam" "$base"
        mem_cli_same_file "identity-set SAM identical, $first_tier vs $tier" "$work/ident.$first_tier.sam" "$ident"
    fi
done

if ((MEM_CLI_FAILS > 0)); then
    echo "FAIL: mem_seed_cli: $MEM_CLI_FAILS of $MEM_CLI_CHECKS checks failed"
    exit 1
fi
echo "PASS: mem_seed_cli: $MEM_CLI_CHECKS checks, $MEM_CLI_SKIPS skipped"
