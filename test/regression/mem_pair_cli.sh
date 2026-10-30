#!/usr/bin/env bash
# test/regression/mem_pair_cli.sh
#
# End-to-end checks of the `bwa-mem3 mem` paired-end options that no other
# script reaches: mate rescue (-S, -m, -P) in both directions and at the edge of
# the rescue window, an FF library, -I with max and min and with the std left
# to its default, the insert-size estimate feeding the proper-pair flag, -K's
# per-batch estimation, and -j. On the synthetic reference
# (test/fixtures/mem_synth.fa; the header of test/tools/mem_synth_ref.c lists
# its contigs and planted repeats), with every expected value worked out from
# the reference by hand. The pairing and MAPQ arithmetic is unit-tested
# (test/unit/test_mem_pair_*.cpp, test/unit/test_mapq_primary.cpp);
# determinism across -t, --dedup-reads, cohort slicing and streaming,
# --compat, --meth pairs and the rescue shortcuts have scripts of their own
# (see README.md).
#
# Reads are 100 bp slices of the reference. A pair "contig start tlen" is read
# 1 = the 100 bases at start (forward) and read 2 = the reverse complement of
# the 100 bases ending at start + tlen - 1, so the fragment is tlen long and an
# FR pair has TLEN +tlen / -tlen. The insert size bwa estimates and bounds is
# the distance from read 1's first base to read 2's last, tlen - 1. A read made
# "unseedable" carries substitutions at read offsets 10, 25, 40, 55, 70 and 85
# (each base -> the next in ACGT order): no exact match reaches the 19 bp
# minimum seed, so it can only be placed by mate rescue, where it scores
# 94 - 6 * 4 = 70 with NM 6.
#
# Checks under each SIMD tier this host can run (mem_cli_tiers; the rescue
# kernels differ per tier), whose SAM must also be byte-identical across tiers:
#   rescue  under -I 300,30 (window [180, 420]): pair a (ctg1:1001, tlen 300,
#         read 2 unseedable) and pair b (ctg2:501, tlen 300, read 1
#         unseedable) are rescued in both directions: 99/147 at the fragment
#         ends, MAPQ 60, 100M, TLEN +-300. Pair e (ctg1:2401, tlen 420, read 2
#         unseedable) puts the mate's last base 419 bases past read 1's first,
#         the last position inside the rescue window: rescued the same way,
#         99/147, TLEN +-420. Pair c: read 1 inside block R (ctg2:1101, three
#         exact copies), read 2 unseedable next to the ctg2 copy: only that copy
#         pairs, so read 1 lands at ctg2:1101 with MAPQ min(q_pe, 0 + 40) = 40;
#         read 2 MAPQ 60.
#   FF    12 pairs with read 2 forward 300 bp after read 1 (insert 300, a
#         400 bp fragment) fit FF to [300, 300]; a 13th pair's unseedable read
#         2 (ctg3:1001) is rescued on read 1's strand. FLAGs 67/131 (no strand
#         bits), MAPQ 60, TLEN +-301 (bwa measures TLEN between 5' ends:
#         301 = 1001 - 701 + 1).
#
# Checks of the remaining options, run once (under the first tier):
#   -S -m 0  no rescue, so the unseedable ends are unmapped (a: 73/133 at
#         ctg1:1001, b: 117/185 at ctg2:701, the unmapped end taking its mate's
#         position and strand) and c's read 1 has MAPQ 0.
#   -P    rescue still runs but nothing is paired: a, b and e are rescued at
#         their fragment ends as 97/145 (inside the window on one contig, so
#         without -P the unpaired rule would set 0x2), c's read 2 is rescued at
#         ctg2:1301 (AS 70) without 0x2, and c's read 1 has MAPQ 0.
#   -I max,min  -I 300,30,350,250 sets the FR window to [250, 350] directly:
#         pairs at inserts 249, 250, 350 and 351 are 97/145, 99/147, 99/147 and
#         97/145 (an edge pair is 1.67 std out, pair score 199 > 183).
#   -I mean  read 1 = ctg1:1201-1270 + ctg3:2601-2630 (two primary hits, so not
#         paired: 0x2 comes from the unpaired rule), read 2 reverse with its
#         last base at ctg1:1201 + d. Under -I 300, which takes the std as 10%
#         of the mean (the help text), the window is [180, 420]: d = 180 gives
#         99/147 and d = 421 gives 97/145, the supplementary ctg3 record
#         (70H30M) carrying the same bits; TLEN +-(d + 1), MAPQ 60.
#   pestat  32 FR pairs on unique sequence, inserts (tlen - 1) 4 x 240,
#         10 x 300, 10 x 330, 4 x 390, and probes at 146, 147, 483 and 484. In
#         bwa's estimator the 32 sorted inserts give p25 = 300 (index 8) and
#         p75 = 330 (index 24); [240, 390] (2 IQR) keeps the 28 main inserts:
#         avg 315, std sqrt(49500 / 28) = 42.046; the 3 IQR window [210, 420] is
#         widened to avg -/+ 4 std, [147, 483]. A probe at an edge is 3.996 std
#         out, pair score (int)(200 - 6.43 + .499) = 194 > 200 - 17, so every
#         pair is proper (99/147) except the probes at 146 and 484 (97/145);
#         MAPQ 60.
#   -K    15 pairs at tlen 301, then 16 at tlen 601 (200 bases per pair). The
#         insert size is estimated per batch, and a batch ends with the pair
#         that brings it to -K bases. -K 3000: batches of 15, 15 and 1 pairs;
#         the first two fit their own inserts exactly ([300, 300], [600, 600]),
#         the last fits nothing: every pair proper but the 31st. -K 3001:
#         batches of 16 and 15; the first fits [300, 300] from its 15 inserts
#         of 300, so its 16th pair (tlen 601) alone is not proper.
#   -j    a reference of ctgA = ctg1:1-2000 and ctgB = ctg1:501-1500 with one
#         substitution (at ctg1:1051): a read at ctg1:1001-1100 scores 100 on
#         ctgA and 95 on ctgB (XS 95). With ctgB in a .alt sidecar the ALT hit
#         does not compete: MAPQ 60, XS 0. -j ignores the sidecar: the output
#         equals the plain index's.
#
# Inputs (env vars):
#   BWA_MEM3       path to the bwa-mem3 binary under test
#   MEM_CLI_TIERS  (optional) tiers to run instead of every runnable one
#   MEM_CLI_KEEP   (optional) 1 keeps the work directory
# No staged inputs: the reference is committed, and the reads and the -j
# reference are built from it.

set -euo pipefail

# shellcheck source=test/framework/mem_cli.sh
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../framework" && pwd)/mem_cli.sh"

mem_cli_init mem_pair_cli
mem_cli_setup_synth
ref="$MEM_CLI_SYNTH"
work="$MEM_CLI_WORK"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# pair_reads <spec> <out-prefix> [ff]
#   One pair per spec line, "name contig start tlen mut1 mut2 ...": read 1 is
#   the 100 bases at start, read 2 the reverse complement of the 100 bases
#   ending at start + tlen - 1 (with "ff": those bases as they are, an FF
#   pair). mutK "s" makes read K unseedable (substitutions at offsets 10, 25,
#   ..., 85 of the forward-strand slice), "0" keeps it. Writes <prefix>.1.fq
#   and <prefix>.2.fq.
pair_reads() {
    awk -v out1="$2.1.fq" -v out2="$2.2.fq" -v ff="${3:-}" '
        function nextb(b) { return b == "A" ? "C" : b == "C" ? "G" : b == "G" ? "T" : "A" }
        function spaced(s, i, p) {
            for (i = 0; i < 6; i++) {
                p = 10 + 15 * i
                s = substr(s, 1, p - 1) nextb(substr(s, p, 1)) substr(s, p + 1)
            }
            return s
        }
        function rc(s, r, i, b) {
            r = ""
            for (i = length(s); i > 0; i--) {
                b = substr(s, i, 1)
                r = r (b == "A" ? "T" : b == "C" ? "G" : b == "G" ? "C" : b == "T" ? "A" : "N")
            }
            return r
        }
        function mut(s, m) { return m == "s" ? spaced(s) : s }
        BEGIN { for (i = 0; i < 100; i++) q = q "I" }
        FNR == NR {
            if (/^>/) c = substr($1, 2)
            else seq[c] = seq[c] $0
            next
        }
        {
            r1 = mut(substr(seq[$2], $3, 100), $5)
            r2 = mut(substr(seq[$2], $3 + $4 - 100, 100), $6)
            printf "@%s\n%s\n+\n%s\n", $1, r1, q > out1
            printf "@%s\n%s\n+\n%s\n", $1, (ff == "ff" ? r2 : rc(r2)), q > out2
        }' "$ref" "$1"
}

# pair_want <spec>
#   The expected primary records of an all-FR spec whose columns 7-9 are
#   "proper mapq1 mapq2", one line per end: "name end FLAG RNAME POS MAPQ CIGAR
#   RNEXT PNEXT TLEN". FLAG is 99/147 for a proper pair, 97/145 otherwise.
pair_want() {
    awk '{
        p2 = $3 + $4 - 100
        print $1, 1, ($7 ? 99 : 97), $2, $3, $8, "100M", "=", p2, $4
        print $1, 2, ($7 ? 147 : 145), $2, p2, $9, "100M", "=", $3, -$4
    }' "$1"
}

# pair_got <sam>
#   The same fields of every primary record (neither 0x100 nor 0x800), in
#   output order; "end" is 1 for 0x40 and 2 otherwise.
pair_got() {
    awk -F '\t' '!/^@/ {
        f = $2 + 0
        if (int(f / 256) % 2 || int(f / 2048) % 2) next
        print $1, (int(f / 64) % 2 ? 1 : 2), $2, $3, $4, $5, $6, $7, $8, $9
    }' "$1"
}

# expect_records <label> <want-file> <sam> [<qname>...]
#   One check: the pair_got lines of <sam> (only the named reads, when given)
#   equal <want-file>.
expect_records() {
    local label="$1" want="$2" sam="$3"
    shift 3
    pair_got "$sam" | awk -v names="$*" '
        BEGIN { n = split(names, a, " "); for (i = 1; i <= n; i++) keep[a[i]] = 1 }
        n == 0 || ($1 in keep)' > "$sam.got"
    mem_cli_same_file "$label" "$want" "$sam.got"
}

# primary_fields <sam> <qname> <end> <column-or-TAG>...
#   The listed columns (numbers, mem_cli_col) and optional-field values (tags,
#   mem_cli_tag) of the read's primary record (mem_cli_primary; <end> "" for a
#   single-end read), space-separated; a missing one is left empty.
primary_fields() {
    local rec f out=""
    rec="$(mem_cli_primary "$1" "$2" "$3")"
    shift 3
    for f in "$@"; do
        if [[ $f =~ ^[0-9]+$ ]]; then
            out="$out $(mem_cli_col "$rec" "$f")"
        else
            out="$out $(mem_cli_tag "$rec" "$f")"
        fi
    done
    printf '%s\n' "${out# }"
}

# primary_bit <sam> <qname> <end> <bit>
#   1 or 0 for FLAG bit <bit> (a power of two) of the read's primary record,
#   "none" when there is no single primary record.
primary_bit() {
    primary_fields "$1" "$2" "$3" 2 | awk -v b="$4" '{ print ($1 == "" ? "none" : int($1 / b) % 2) }'
}

# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

# Rescue pairs a, b, c, e (see the header). Columns 7-9 are the default-run
# expectations: proper, MAPQ of read 1, MAPQ of read 2.
cat > "$work/abc.spec" << 'EOF'
a ctg1 1001 300 0 s 1 60 60
b ctg2 501 300 s 0 1 60 60
c ctg2 1101 300 0 s 1 40 60
e ctg1 2401 420 0 s 1 60 60
EOF
pair_reads "$work/abc.spec" "$work/abc"
pair_want "$work/abc.spec" > "$work/abc.want"

# pair_spec <name-prefix> <seed> <proper-lo> <proper-hi> <tlen>...
#   A spec of one FR pair per <tlen>, each fragment inside one of six unique
#   stretches of ctg1/ctg2 (clear of blocks R and S and their copies) in turn,
#   at a Park-Miller offset; proper iff the insert tlen - 1 is in
#   [<proper-lo>, <proper-hi>], MAPQ 60 for both ends.
pair_spec() {
    local prefix="$1" seed="$2" lo="$3" hi="$4"
    shift 4
    awk -v prefix="$prefix" -v x="$seed" -v lo="$lo" -v hi="$hi" -v tlens="$*" 'BEGIN {
        split("ctg1 ctg1 ctg1 ctg2 ctg2 ctg2", wc, " ")
        split("1 2301 3501 1 1301 2501", wlo, " ")
        split("2000 3000 5000 1000 2000 4000", whi, " ")
        n = split(tlens, tl, " ")
        for (k = 0; k < n; k++) {
            t = tl[k + 1]
            w = k % 6 + 1
            x = (x * 48271) % 2147483647
            st = wlo[w] + x % (whi[w] - wlo[w] + 2 - t)
            printf "%s%02d %s %d %d 0 0 %d 60 60\n", prefix, k, wc[w], st, t, (t - 1 >= lo && t - 1 <= hi)
        }
    }'
}

# -I with max and min: pairs at inserts 249, 250, 350 and 351.
pair_spec m 2024 250 350 250 251 351 352 > "$work/imm.spec"
pair_reads "$work/imm.spec" "$work/imm"
pair_want "$work/imm.spec" > "$work/imm.want"

# Insert-size estimation (see the header): window [147, 483].
pair_spec e 4242 147 483 241 241 241 241 301 301 301 301 301 301 301 301 301 301 \
    331 331 331 331 331 331 331 331 331 331 391 391 391 391 147 148 484 485 > "$work/est.spec"
pair_reads "$work/est.spec" "$work/est"
pair_want "$work/est.spec" > "$work/est.want"

# Batch boundary: 15 pairs at tlen 301, then 16 at tlen 601 (all proper by
# their own batch's window, except where the header says).
pair_spec k 777 1 100000 301 301 301 301 301 301 301 301 301 301 301 301 301 301 301 \
    601 601 601 601 601 601 601 601 601 601 601 601 601 601 601 601 > "$work/kb.spec"
pair_reads "$work/kb.spec" "$work/kb"
awk '$1 == "k30" { $7 = 0 } { print }' "$work/kb.spec" > "$work/kb3000.spec"
pair_want "$work/kb3000.spec" > "$work/kb3000.want"
awk '$1 == "k15" { $7 = 0 } { print }' "$work/kb.spec" > "$work/kb3001.spec"
pair_want "$work/kb3001.spec" > "$work/kb3001.want"

# Split read 1 at distances 180 and 421 (see the header). The bases across the
# junction differ (ctg1:1271 A vs ctg3:2601 C; ctg3:2600 C vs ctg1:1270 T), so
# neither part extends into the other.
split_part2="$(mem_cli_seq "$ref" ctg3 2601 30)"
: > "$work/split.1.fq"
: > "$work/split.2.fq"
: > "$work/split.want"
for d in 180 421; do
    mem_cli_fq "s$d" "$(mem_cli_seq "$ref" ctg1 1201 70)$split_part2" >> "$work/split.1.fq"
    mem_cli_fq "s$d" "$(dna_revcomp "$(mem_cli_seq "$ref" ctg1 $((1201 + d - 99)) 100)")" >> "$work/split.2.fq"
    if ((d >= 180 && d <= 420)); then pp=2; else pp=0; fi
    printf '%s\n' \
        "s$d $((97 + pp)) ctg1 1201 60 70M30S = $((1201 + d - 99)) $((d + 1))" \
        "s$d $((2145 + pp)) ctg3 2601 60 70H30M ctg1 $((1201 + d - 99)) 0" \
        "s$d $((145 + pp)) ctg1 $((1201 + d - 99)) 60 100M = 1201 -$((d + 1))" >> "$work/split.want"
done

# An FF library: 12 FF pairs at insert 300 (tlen 400) and one whose read 2 is
# unseedable, rescued on the anchor's strand.
pair_spec f 31337 1 100000 400 400 400 400 400 400 400 400 400 400 400 400 > "$work/ff.spec"
echo 'fr ctg3 701 400 0 s 1 60 60' >> "$work/ff.spec"
pair_reads "$work/ff.spec" "$work/ff" ff
awk '{
    p2 = $3 + $4 - 100
    print $1, 1, 67, $2, $3, 60, "100M", "=", p2, $4 - 99
    print $1, 2, 131, $2, p2, 60, "100M", "=", $3, -($4 - 99)
}' "$work/ff.spec" > "$work/ff.want"

# -j: ctgA = ctg1:1-2000, ctgB = ctg1:501-1500 with its base 551 (ctg1:1051)
# moved to the next base in ACGT order. Indexed once; the ALT copy of the index
# adds ctgB.alt as the sidecar.
mkdir -p "$work/jplain" "$work/jalt"
jplain="$work/jplain/$(mem_cli_ref_name jplain)"
jalt="$work/jalt/$(mem_cli_ref_name jalt)"
ctgb="$(mem_cli_seq "$ref" ctg1 501 1000)"
ctgb_sub="$(tr ACGT CGTA <<< "${ctgb:550:1}")"
{
    echo ">ctgA"
    mem_cli_seq "$ref" ctg1 1 2000 | fold -w 60
    echo ">ctgB"
    mem_cli_subst "$ctgb" 551 "$ctgb_sub" | fold -w 60
} > "$jplain"
mem_cli_index "$jplain"
for f in "$jplain"*; do cp "$f" "$jalt${f#"$jplain"}"; done
printf 'ctgB\n' > "$jalt.alt"
mem_cli_fq j1 "$(mem_cli_seq "$ref" ctg1 1001 100)" > "$work/j.fq"

# ---------------------------------------------------------------------------
# Rescue, under each tier
# ---------------------------------------------------------------------------

first_tier=""
for tier in $(mem_cli_tiers); do
    first_tier="${first_tier:-$tier}"

    # Mate rescue in both directions, at the window edge and through a repeat
    # anchor.
    mem_cli_mem "$tier" "$work/abc.$tier.sam" -I 300,30 "$ref" "$work/abc.1.fq" "$work/abc.2.fq"
    expect_records "$tier: rescue: records of pairs a, b, c, e" "$work/abc.want" "$work/abc.$tier.sam"

    # An FF library: the FF window is fitted and its rescue runs forward.
    mem_cli_mem "$tier" "$work/ff.$tier.sam" "$ref" "$work/ff.1.fq" "$work/ff.2.fq"
    expect_records "$tier: FF library: 67/131, TLEN +-301" "$work/ff.want" "$work/ff.$tier.sam"

    if [[ $tier != "$first_tier" ]]; then
        for s in abc ff; do
            mem_cli_same_file "$tier: $s SAM equals $first_tier's" "$work/$s.$first_tier.sam" "$work/$s.$tier.sam"
        done
    fi
done

if [[ $(mem_cli_tiers | wc -l) -lt 2 ]]; then
    mem_cli_skip "SAM identical across tiers" "only one tier ($first_tier) runs on this host"
fi

# ---------------------------------------------------------------------------
# The remaining options, once (under the first tier)
# ---------------------------------------------------------------------------

tier="$first_tier"

# -S and -m 0: no rescue.
printf '%s\n' \
    'a 1 73 ctg1 1001 60 100M = 1001 0' 'a 2 133 ctg1 1001 0 * = 1001 0' \
    'b 1 117 ctg2 701 0 * = 701 0' 'b 2 185 ctg2 701 60 100M = 701 0' > "$work/norescue.want"
for opt in -S "-m 0"; do
    sam="$work/norescue.${opt// /}.sam"
    # shellcheck disable=SC2086 # "-m 0" is two arguments on purpose
    mem_cli_mem "$tier" "$sam" $opt -I 300,30 "$ref" "$work/abc.1.fq" "$work/abc.2.fq"
    expect_records "$opt: pairs a and b unrescued" "$work/norescue.want" "$sam" a b
    got="$(primary_fields "$sam" c 1 5) $(primary_bit "$sam" c 2 4)"
    mem_cli_expect "$opt: pair c read 1 MAPQ 0, read 2 unmapped" "0 1" "$got"
done

# -P: rescue runs, pairing does not. Pairs a, b and e are rescued on one contig
# inside the window, so only -P keeps the unpaired rule from setting 0x2.
awk '$1 != "c" { $7 = 0; print }' "$work/abc.spec" > "$work/nopair.spec"
pair_want "$work/nopair.spec" > "$work/nopair.want"
sam="$work/nopair.sam"
mem_cli_mem "$tier" "$sam" -P -I 300,30 "$ref" "$work/abc.1.fq" "$work/abc.2.fq"
expect_records "-P: pairs a, b and e rescued, no 0x2" "$work/nopair.want" "$sam" a b e
got="$(primary_fields "$sam" c 2 3 4 AS) $(primary_bit "$sam" c 2 2) $(primary_fields "$sam" c 1 5)"
mem_cli_expect "-P: pair c read 2 rescued at ctg2:1301, no 0x2, read 1 MAPQ 0" "ctg2 1301 70 0 0" "$got"

# -I mean,std,max,min: the window is [min, max].
sam="$work/imm.sam"
mem_cli_mem "$tier" "$sam" -I 300,30,350,250 "$ref" "$work/imm.1.fq" "$work/imm.2.fq"
expect_records "-I 300,30,350,250: proper exactly in [250, 350]" "$work/imm.want" "$sam"

# -I 300 alone implies std 30: the unpaired rule's window is [180, 420].
sam="$work/split.sam"
mem_cli_mem "$tier" "$sam" -I 300 "$ref" "$work/split.1.fq" "$work/split.2.fq"
awk -F '\t' '!/^@/ { print $1, $2, $3, $4, $5, $6, $7, $8, $9 }' "$sam" > "$sam.got"
mem_cli_same_file "-I 300: split read 1, unpaired 0x2 window [180, 420]" "$work/split.want" "$sam.got"

# Insert-size estimation from the data feeds the proper-pair flag.
sam="$work/est.sam"
mem_cli_mem "$tier" "$sam" "$ref" "$work/est.1.fq" "$work/est.2.fq"
expect_records "estimated window [147, 483]: proper exactly inside" "$work/est.want" "$sam"

# The batch boundary and per-batch estimation.
for k in 3000 3001; do
    sam="$work/kb.K$k.sam"
    mem_cli_mem "$tier" "$sam" -K "$k" "$ref" "$work/kb.1.fq" "$work/kb.2.fq"
    expect_records "-K $k: per-batch windows" "$work/kb$k.want" "$sam"
done

# -j and the ALT sidecar.
mem_cli_mem "$tier" "$work/jplain.sam" "$jplain" "$work/j.fq"
mem_cli_mem "$tier" "$work/jalt.sam" "$jalt" "$work/j.fq"
mem_cli_mem "$tier" "$work/jalt-j.sam" -j "$jalt" "$work/j.fq"
for e in "jalt 60 0" "jalt-j 10 95"; do
    read -r idx mapq xs <<< "$e"
    got="$(primary_fields "$work/$idx.sam" j1 "" 3 4 5 XS)"
    mem_cli_expect "$idx: j1 RNAME POS MAPQ XS" "ctgA 1001 $mapq $xs" "$got"
done
mem_cli_same_file "-j equals the index without a sidecar" "$work/jplain.sam" "$work/jalt-j.sam"

if mem_cli_verdict; then
    echo "PASS: $MEM_CLI_VERDICT"
else
    echo "FAIL: $MEM_CLI_VERDICT"
    exit 1
fi
