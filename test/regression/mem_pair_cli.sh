#!/usr/bin/env bash
# test/regression/mem_pair_cli.sh
#
# End-to-end checks of `bwa-mem3 mem` paired-end and mode handling:
# insert-size estimation and -I, mate rescue (-S, -P, -m), proper pairs and
# the mate fields, MAPQ, primary/secondary/supplementary marking, --fast,
# --compat=bwa-mem2 and --compat=bwa-mem, --meth, -j, and determinism across
# -K and -t. Each check runs under every SIMD tier this host can run
# (test/framework/mem_cli.sh: mem_cli_tiers), on the synthetic reference
# (test/fixtures/mem_synth_ref.awk, whose header lists its contigs and planted
# repeats), with every expected value worked out from the reference by hand.
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
# Checks:
#   -I    one FR pair, ctg1:1001-1100 forward and the reverse complement of
#         ctg1:1201-1300, under -I 300,30: FLAGs 99 and 147 (paired, proper,
#         first/second, mate reverse / read reverse), RNEXT =, PNEXT the mate's
#         POS, and TLEN +300 / -300 (the fragment is ctg1:1001-1300).
#   -I max,min  -I 300,30,350,250 sets the FR window to [250, 350] directly:
#         pairs at inserts 249, 250, 350 and 351 are 97/145, 99/147, 99/147 and
#         97/145 (an edge pair is 1.67 std out, pair score 199 > 183).
#   rescue  under -I 300,30 (window [180, 420]): pair a (ctg1:1001, tlen 300,
#         read 2 unseedable) and pair b (ctg2:501, tlen 300, read 1
#         unseedable) are rescued in both directions: 99/147 at the fragment
#         ends, MAPQ 60, 100M, TLEN +-300, and the rescued end has AS 70, NM 6
#         and an MD naming the six reference bases. Pair c: read 1 inside block
#         R (ctg2:1101, three exact copies), read 2 unseedable next to the ctg2
#         copy: only that copy pairs, so read 1 lands at ctg2:1101 with MAPQ
#         min(q_pe, 0 + 40), q_pe = 6.02 * (170 - (170 - U)) = 6.02 U capped at
#         60: 40 at the default -U 17, 24 at -U 4 and 12 at -U 2; read 2 MAPQ 60.
#         -S and -m 0: no rescue, so the unseedable ends are unmapped (a: 73/133
#         at ctg1:1001, b: 117/185 at ctg2:701, the unmapped end taking its
#         mate's position and strand) and c's read 1 has MAPQ 0. -P: rescue
#         still runs but nothing is paired: 97/145, no 0x2, TLEN still +-300.
#         With every BWA3_RESCUE_* shortcut off the output is byte-identical.
#         --fast places all three pairs the same way.
#   FF    12 pairs with read 2 forward 300 bp after read 1 (insert 300, a
#         400 bp fragment) fit FF to [300, 300]; a 13th pair's unseedable read
#         2 (ctg3:1001) is rescued on read 1's strand. FLAGs 67/131 (no strand
#         bits), MAPQ 60, TLEN +-301 (bwa measures TLEN between 5' ends:
#         301 = 1001 - 701 + 1), the rescued end AS 70, NM 6.
#   split  read 1 = ctg1:1201-1270 + ctg3:2601-2630 (two primary hits, 70 and
#         30 >= -T), read 2 reverse with its last base at ctg1:1201 + d, under
#         -I 300,30. An end with two primary hits is not paired, so 0x2 comes
#         from the unpaired rule: read 1's top hit and read 2 at distance d
#         inside [180, 420]. d = 179 and 421: 97/145; d = 180 and 420: 99/147,
#         and the supplementary ctg3 record (70H30M, MAPQ 60, TLEN 0) carries
#         the same bits (2145/2147). TLEN +-(d + 1), MAPQ 60.
#   --compat  bwa-mem2: the records equal the default's without MQ:i and HN:i,
#         and there is no @HD; bwa-mem: without HN:i only, and @HD is kept.
#   pestat  32 FR pairs on unique sequence, inserts (tlen - 1) 4 x 240,
#         10 x 300, 10 x 330, 4 x 390, and probes at 146, 147, 483 and 484. In
#         bwa's estimator the 32 sorted inserts give p25 = 300 (index 8) and
#         p75 = 330 (index 24); [240, 390] (2 IQR) keeps the 28 main inserts:
#         avg 315, std sqrt(49500 / 28) = 42.046; the 3 IQR window [210, 420] is
#         widened to avg -/+ 4 std, [147, 483]. A probe at an edge is 3.996 std
#         out, pair score (int)(200 - 6.43 + .499) = 194 > 200 - 17, so every
#         pair is proper (99/147) except the probes at 146 and 484 (97/145);
#         MAPQ 60. Nine pairs at tlen 301 fit nothing (fewer than 10): no 0x2;
#         ten fit [300, 300]: all proper.
#   -K    15 pairs at tlen 301, then 16 at tlen 601 (200 bases per pair). The
#         insert size is estimated per batch, and a batch ends with the pair
#         that brings it to -K bases. -K 3000: batches of 15, 15 and 1 pairs;
#         the first two fit their own inserts exactly ([300, 300], [600, 600]),
#         the last fits nothing: every pair proper but the 31st, at -t 1 and
#         -t 4. -K 3001: batches of 16 and 15; the first fits [300, 300] from
#         its 15 inserts of 300, so its 16th pair (tlen 601) alone is not
#         proper.
#   -t    101 pairs (the sets above, 4 unrelated pairs, 4 pairs inside block R
#         placed by the read id's hash, 10 pairs each followed by a copy) at
#         -K 2000: byte-identical at -t 1, 2 and 4, with -1, with --dedup-reads
#         off and on (the stats must report memoized duplicates; also with
#         BWAMEM3_DEDUP_READS_VERIFY=1 and _CHECK_BASES=1, which align every
#         duplicate anyway and abort on a difference), with --cohort-slices 0,
#         with every batch read in slices
#         (BWA_MEM3_COHORT_SLICE_ALL=1, whose slices must show on stderr) and
#         with output streaming off (BWA_MEM3_STREAM_OUTPUT=0). Under -I 300,30,
#         -K 1 (one pair per batch) equals the default batch size byte for
#         byte, and every pair comes out once, in input order. The -K 2000
#         output is identical across tiers.
#   supplementary  read = ctg1:2101-2170 (in R) + ctg3:2501-2550 (unique): the
#         primary is the 70 bp R part (MAPQ 0, three copies), the ctg3 part is
#         supplementary (FLAG 2048, POS 2501, 70H50M, AS 50) with its MAPQ 60
#         lowered to the primary's 0; -q keeps it at 60.
#   -j    a reference of ctgA = ctg1:1-2000 and ctgB = ctg1:501-1500 with one
#         substitution (at ctg1:1051): a read at ctg1:1001-1100 scores 100 on
#         ctgA and 95 on ctgB, MAPQ (int)(6.02 * 5 * (3/ln 100)^2 + .499) - 3 =
#         10, XS 95. With ctgB in a .alt sidecar the ALT hit does not compete:
#         MAPQ 60, XS 0. -j ignores the sidecar: the output equals the plain
#         index's.
#   --meth  on an `index --meth` copy, a directional OT pair (fragment
#         ctg1:1001-1300, top strand converted C->T except at CpG) and an OB pair
#         (fragment ctg1:3601-3900, bottom strand converted, i.e. G->A on the top
#         except at CpG): OT 99/147, OB read 1 83 at 3801 and read 2 163 at 3601,
#         TLEN +-300, MAPQ 60, AS 100, NM 0 (conversions are free), XR CT for
#         read 1 and GA for read 2, XG CT / GA, and XM computed from the
#         reference by Bismark's rule (z/x/h for a C in CpG/CHG/CHH context on
#         the converted strand, upper case when it was left unconverted).
#
# Inputs (env vars):
#   BWA_MEM3       path to the bwa-mem3 binary under test
#   MEM_CLI_TIERS  (optional) tiers to run instead of every runnable one
#   MEM_CLI_KEEP   (optional) 1 keeps the work directory
# No staged inputs: the references are generated, the reads are sliced from
# them.

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
#   ..., 85 of the forward-strand slice), "u" replaces it with unrelated bases
#   (a Park-Miller stream), "0" keeps it. Writes <prefix>.1.fq and
#   <prefix>.2.fq.
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
        function junk(n, s, i) {
            s = ""
            for (i = 0; i < n; i++) {
                x = (x * 48271) % 2147483647
                s = s substr("ACGT", int(x / 536870912) + 1, 1)
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
        function mut(s, m) { return m == "s" ? spaced(s) : m == "u" ? junk(length(s)) : s }
        BEGIN {
            x = 20261024
            for (i = 0; i < 100; i++) q = q "I"
        }
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

# spaced_md <contig> <start>
#   The MD of an unseedable read cut from <contig> at <start>, aligned in full:
#   9 matches, the reference base at offset 10, 14 matches, ..., 15 matches.
spaced_md() {
    local s
    s="$(mem_cli_seq "$ref" "$1" "$2" 100)"
    printf '9%s14%s14%s14%s14%s14%s15\n' "${s:9:1}" "${s:24:1}" "${s:39:1}" \
        "${s:54:1}" "${s:69:1}" "${s:84:1}"
}

# strip_tags <sam> <tag-regex>
#   The records of <sam> without the optional fields whose "TAG:TYPE" prefix
#   matches the regex (e.g. '^(MQ|HN):i:').
strip_tags() {
    awk -F '\t' -v re="$2" '!/^@/ {
        out = $1
        for (i = 2; i <= NF; i++) if (i < 12 || $i !~ re) out = out "\t" $i
        print out
    }' "$1"
}

# mem_env <tier> <out> <VAR=value>... -- <mem args>...
#   mem_cli_mem with extra environment variables for the one run.
mem_env() {
    local tier="$1" out="$2" envs=()
    shift 2
    while [[ $1 != -- ]]; do
        envs+=("$1")
        shift
    done
    shift
    if ! env ${envs[@]+"${envs[@]}"} BWAMEM3_FORCE_TIER="$tier" "$BWA_MEM3" mem "$@" > "$out" 2> "$out.err"; then
        echo "FAIL: $MEM_CLI_NAME: ${envs[*]} bwa-mem3 mem $* (tier $tier) exited nonzero:" >&2
        tail -n 20 "$out.err" >&2
        exit 1
    fi
    if ! grep -q "^Executing in $tier mode" "$out.err"; then
        echo "FAIL: $MEM_CLI_NAME: BWAMEM3_FORCE_TIER=$tier was not honoured:" >&2
        head -n 5 "$out.err" >&2
        exit 1
    fi
}

# record_fields <record> <column-or-TAG>...
#   The listed columns (numbers) and optional-field values (two-letter tags) of
#   one SAM record, space-separated; a missing one is left empty.
record_fields() {
    local rec="$1"
    shift
    printf '%s\n' "$rec" | awk -F '\t' -v want="$*" '{
        n = split(want, w, " ")
        out = ""
        for (k = 1; k <= n; k++) {
            v = ""
            if (w[k] ~ /^[0-9]+$/) v = $(w[k])
            else for (i = 12; i <= NF; i++) if (substr($i, 1, 3) == w[k] ":") { v = substr($i, 6); break }
            out = out (k > 1 ? " " : "") v
        }
        print out
    }'
}

# primary_fields <sam> <qname> <end> <column-or-TAG>...
#   record_fields of the read's primary record (mem_cli_primary; <end> "" for a
#   single-end read).
primary_fields() {
    local rec
    rec="$(mem_cli_primary "$1" "$2" "$3")"
    shift 3
    record_fields "$rec" "$@"
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

mem_cli_fq p1 "$(mem_cli_seq "$ref" ctg1 1001 100)" > "$work/r1.fq"
mem_cli_fq p1 "$(mem_cli_revcomp "$(mem_cli_seq "$ref" ctg1 1201 100)")" > "$work/r2.fq"

# Rescue pairs a, b, c (see the header). Columns 7-9 are the default-run
# expectations: proper, MAPQ of read 1, MAPQ of read 2.
cat > "$work/abc.spec" << 'EOF'
a ctg1 1001 300 0 s 1 60 60
b ctg2 501 300 s 0 1 60 60
c ctg2 1101 300 0 s 1 40 60
EOF
pair_reads "$work/abc.spec" "$work/abc"
md_a2="$(spaced_md ctg1 1201)"
md_b1="$(spaced_md ctg2 501)"
md_c2="$(spaced_md ctg2 1301)"

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
# Nine and ten pairs at tlen 301: proper only with ten.
pair_spec n 99 300 300 301 301 301 301 301 301 301 301 301 301 > "$work/cnt10.spec"
head -n 9 "$work/cnt10.spec" | awk '{ $7 = 0; print }' > "$work/cnt9.spec"
for n in 9 10; do
    pair_reads "$work/cnt$n.spec" "$work/cnt$n"
    pair_want "$work/cnt$n.spec" > "$work/cnt$n.want"
done

# Batch boundary: 15 pairs at tlen 301, then 16 at tlen 601 (all proper by
# their own batch's window, except where the header says).
pair_spec k 777 1 100000 301 301 301 301 301 301 301 301 301 301 301 301 301 301 301 \
    601 601 601 601 601 601 601 601 601 601 601 601 601 601 601 601 > "$work/kb.spec"
pair_reads "$work/kb.spec" "$work/kb"
awk '$1 == "k30" { $7 = 0 } { print }' "$work/kb.spec" > "$work/kb3000.spec"
pair_want "$work/kb3000.spec" > "$work/kb3000.want"
awk '$1 == "k15" { $7 = 0 } { print }' "$work/kb.spec" > "$work/kb3001.spec"
pair_want "$work/kb3001.spec" > "$work/kb3001.want"

# Split read 1 at distances 179, 180, 420 and 421 (see the header). The bases
# across the junction differ (ctg1:1271 A vs ctg3:2601 C; ctg3:2600 C vs
# ctg1:1270 T), so neither part extends into the other.
split_part2="$(mem_cli_seq "$ref" ctg3 2601 30)"
: > "$work/split.1.fq"
: > "$work/split.2.fq"
: > "$work/split.want"
for d in 179 180 420 421; do
    mem_cli_fq "s$d" "$(mem_cli_seq "$ref" ctg1 1201 70)$split_part2" >> "$work/split.1.fq"
    mem_cli_fq "s$d" "$(mem_cli_revcomp "$(mem_cli_seq "$ref" ctg1 $((1201 + d - 99)) 100)")" >> "$work/split.2.fq"
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
md_fr2="$(spaced_md ctg3 1001)"

# Determinism: everything above, with each of the first 10 estimation pairs
# followed by a copy under a new name (duplicate pairs in one batch, which
# --dedup-reads aligns once), 4 pairs of unrelated bases, 4 pairs wholly inside
# block R (three equal pairs each, so the pair is chosen by the read id's
# hash), the split pairs and the FF library.
{
    cat "$work/abc.spec"
    awk '{ print } NR <= 10 { $1 = "d" $1; print }' "$work/est.spec"
    cat "$work/kb.spec"
    printf 'u%d ctg1 1001 300 u u\n' 1 2 3 4
    printf '%s\n' 'r1 ctg1 2001 300 0 0' 'r2 ctg1 2011 280 0 0' 'r3 ctg1 2021 260 0 0' 'r4 ctg1 2031 250 0 0'
} > "$work/all.spec"
pair_reads "$work/all.spec" "$work/all"
cat "$work/split.1.fq" "$work/ff.1.fq" >> "$work/all.1.fq"
cat "$work/split.2.fq" "$work/ff.2.fq" >> "$work/all.2.fq"
{
    awk '{ print $1 }' "$work/all.spec"
    printf 's%d\n' 179 180 420 421
    awk '{ print $1 }' "$work/ff.spec"
} > "$work/all.names"
n_all="$(awk 'END { print NR }' "$work/all.names")"

# Supplementary: 70 bp of block R, then 50 bp of unique ctg3. The bases across
# the junction differ (ctg1:2171 A vs ctg3:2501 G; ctg3:2500 C vs ctg1:2170 G),
# so neither part extends into the other.
mem_cli_fq chim "$(mem_cli_seq "$ref" ctg1 2101 70)$(mem_cli_seq "$ref" ctg3 2501 50)" > "$work/chim.fq"

# -j: ctgA = ctg1:1-2000, ctgB = ctg1:501-1500 with its base 551 (ctg1:1051)
# moved to the next base in ACGT order. Indexed once; the ALT copy of the index
# adds ctgB.alt as the sidecar.
mkdir -p "$work/jplain" "$work/jalt"
ctgb="$(mem_cli_seq "$ref" ctg1 501 1000)"
case "${ctgb:550:1}" in
    A) ctgb_sub=C ;;
    C) ctgb_sub=G ;;
    G) ctgb_sub=T ;;
    *) ctgb_sub=A ;;
esac
{
    echo ">ctgA"
    mem_cli_seq "$ref" ctg1 1 2000 | fold -w 60
    echo ">ctgB"
    mem_cli_subst "$ctgb" 551 "$ctgb_sub" | fold -w 60
} > "$work/jplain/ref.fa"
mem_cli_index "$work/jplain/ref.fa"
cp "$work"/jplain/ref.fa* "$work/jalt/"
printf 'ctgB\n' > "$work/jalt/ref.fa.alt"
mem_cli_fq j1 "$(mem_cli_seq "$ref" ctg1 1001 100)" > "$work/j.fq"

# --meth: an `index --meth` copy of the reference, an OT and an OB pair.
mkdir -p "$work/meth"
cp "$ref" "$work/meth/ref.fa"
if ! "$BWA_MEM3" index --meth "$work/meth/ref.fa" > "$work/meth/index.log" 2>&1; then
    echo "FAIL: mem_pair_cli: bwa-mem3 index --meth failed:" >&2
    tail -n 20 "$work/meth/index.log" >&2
    exit 1
fi
# meth_pair <strand> <start> <out-prefix>: the directional pair of the 300 bp
# fragment at ctg1:<start>, converted on the top (CT) or bottom (GA) strand,
# CpG left unconverted. CT: read 1 = converted top [start, +100), read 2 =
# reverse complement of converted top [start+200, +100). GA: read 1 = reverse
# complement of converted [start+200, +100), read 2 = converted [start, +100).
meth_pair() {
    awk -v strand="$1" -v st="$2" -v out1="$3.1.fq" -v out2="$3.2.fq" '
        function rc(s, r, i, b) {
            r = ""
            for (i = length(s); i > 0; i--) {
                b = substr(s, i, 1)
                r = r (b == "A" ? "T" : b == "C" ? "G" : b == "G" ? "C" : "A")
            }
            return r
        }
        /^>/ { c = substr($1, 2); next }
        c == "ctg1" { g = g $0 }
        END {
            f = ""
            for (p = st; p < st + 300; p++) {
                b = substr(g, p, 1)
                if (strand == "CT" && b == "C" && substr(g, p + 1, 1) != "G") b = "T"
                if (strand == "GA" && b == "G" && substr(g, p - 1, 1) != "C") b = "A"
                f = f b
            }
            for (i = 0; i < 100; i++) q = q "I"
            a = substr(f, 1, 100)
            z = substr(f, 201, 100)
            if (strand == "CT") { r1 = a; r2 = rc(z) } else { r1 = rc(z); r2 = a }
            printf "@m%s\n%s\n+\n%s\n", strand, r1, q > out1
            printf "@m%s\n%s\n+\n%s\n", strand, r2, q > out2
        }' "$ref"
}
# meth_xm <strand> <pos>: Bismark's XM over ctg1:[pos, pos+100) for a read
# from the converted <strand>, every CpG unconverted (upper case Z) and every
# other C converted (lower case x or h).
meth_xm() {
    awk -v strand="$1" -v pos="$2" '
        /^>/ { c = substr($1, 2); next }
        c == "ctg1" { g = g $0 }
        END {
            s = ""
            for (p = pos; p < pos + 100; p++) {
                b = substr(g, p, 1)
                if (strand == "CT") {
                    if (b != "C") { s = s "."; continue }
                    n1 = substr(g, p + 1, 1); n2 = substr(g, p + 2, 1)
                } else {
                    if (b != "G") { s = s "."; continue }
                    n1 = substr(g, p - 1, 1); n2 = substr(g, p - 2, 1)
                    n1 = n1 == "C" ? "G" : "A"; n2 = n2 == "C" ? "G" : "A"
                }
                s = s (n1 == "G" ? "Z" : n2 == "G" ? "x" : "h")
            }
            print s
        }' "$ref"
}
meth_pair CT 1001 "$work/mct"
meth_pair GA 3601 "$work/mga"
cat "$work/mct.1.fq" "$work/mga.1.fq" > "$work/meth.1.fq"
cat "$work/mct.2.fq" "$work/mga.2.fq" > "$work/meth.2.fq"
xm_ct1="$(meth_xm CT 1001)"
xm_ct2="$(meth_xm CT 1201)"
xm_ga1="$(meth_xm GA 3801)"
xm_ga2="$(meth_xm GA 3601)"

# Every BWA3_RESCUE_* shortcut off (test/regression/rescue_prune_identity.sh's
# reference leg).
rescue_ref_env=(BWA3_RESCUE_PRUNE=0 BWA3_RESCUE_BAND=0 BWA3_RESCUE_DEDUP_SKIP=0
    BWA3_RESCUE_SKIP_REUSE=0 BWA3_RESCUE_REPEAT=0 BWA3_RESCUE_FSCAN=0
    BWA3_RESCUE_USQADD=0 BWA3_RESCUE_ROWPAIR=0 BWA3_RESCUE_LAZYQE=0
    BWA3_RESCUE_SCORE2_BITS=0 BWA3_RESCUE_PADCOLS=0 BWA3_RESCUE_QE_REVSCAN=0)

# ---------------------------------------------------------------------------
# Checks, under each tier
# ---------------------------------------------------------------------------

first_tier=""
for tier in $(mem_cli_tiers); do
    first_tier="${first_tier:-$tier}"

    # -I: one proper FR pair.
    sam="$work/i300.$tier.sam"
    mem_cli_mem "$tier" "$sam" -I 300,30 "$ref" "$work/r1.fq" "$work/r2.fq"
    for end in 1 2; do
        got="$(primary_fields "$sam" p1 "$end" 2 4 7 8 9)"
        if [[ $end == 1 ]]; then
            want="99 1001 = 1201 300"
        else
            want="147 1201 = 1001 -300"
        fi
        mem_cli_expect "$tier: -I 300,30 end $end FLAG POS RNEXT PNEXT TLEN" "$want" "$got"
    done

    # -I mean,std,max,min: the window is [min, max].
    sam="$work/imm.$tier.sam"
    mem_cli_mem "$tier" "$sam" -I 300,30,350,250 "$ref" "$work/imm.1.fq" "$work/imm.2.fq"
    expect_records "$tier: -I 300,30,350,250: proper exactly in [250, 350]" "$work/imm.want" "$sam"

    # Mate rescue in both directions, and through a repeat anchor.
    abc="$work/abc.$tier.sam"
    mem_cli_mem "$tier" "$abc" -I 300,30 "$ref" "$work/abc.1.fq" "$work/abc.2.fq"
    pair_want "$work/abc.spec" > "$work/abc.want"
    expect_records "$tier: rescue: records of pairs a, b, c" "$work/abc.want" "$abc"
    for e in "a 2 $md_a2" "b 1 $md_b1" "c 2 $md_c2"; do
        read -r qn end md <<< "$e"
        got="$(primary_fields "$abc" "$qn" "$end" AS NM MD)"
        mem_cli_expect "$tier: rescue: $qn end $end (rescued) AS NM MD" "70 6 $md" "$got"
    done

    # -U sets q_pe for pair c: read 1 MAPQ min(6.02 U, 40); nothing else moves.
    for u in 4:24 2:12; do
        sam="$work/abc.U${u%%:*}.$tier.sam"
        mem_cli_mem "$tier" "$sam" -U "${u%%:*}" -I 300,30 "$ref" "$work/abc.1.fq" "$work/abc.2.fq"
        awk -v q="${u##*:}" '$1 == "c" { $8 = q } { print }' "$work/abc.spec" > "$work/abc.U.spec"
        pair_want "$work/abc.U.spec" > "$work/abc.U.want"
        expect_records "$tier: -U ${u%%:*}: pair c read 1 MAPQ ${u##*:}" "$work/abc.U.want" "$sam"
    done

    # -S and -m 0: no rescue.
    printf '%s\n' \
        'a 1 73 ctg1 1001 60 100M = 1001 0' 'a 2 133 ctg1 1001 0 * = 1001 0' \
        'b 1 117 ctg2 701 0 * = 701 0' 'b 2 185 ctg2 701 60 100M = 701 0' > "$work/norescue.want"
    for opt in -S "-m 0"; do
        sam="$work/norescue.${opt// /}.$tier.sam"
        # shellcheck disable=SC2086 # "-m 0" is two arguments on purpose
        mem_cli_mem "$tier" "$sam" $opt -I 300,30 "$ref" "$work/abc.1.fq" "$work/abc.2.fq"
        expect_records "$tier: $opt: pairs a and b unrescued" "$work/norescue.want" "$sam" a b
        got="$(primary_fields "$sam" c 1 5) $(primary_bit "$sam" c 2 4)"
        mem_cli_expect "$tier: $opt: pair c read 1 MAPQ 0, read 2 unmapped" "0 1" "$got"
    done

    # -P: rescue runs, pairing does not.
    sam="$work/nopair.$tier.sam"
    mem_cli_mem "$tier" "$sam" -P -I 300,30 "$ref" "$work/abc.1.fq" "$work/abc.2.fq"
    awk '$1 != "c" { $7 = 0; print }' "$work/abc.spec" > "$work/nopair.spec"
    pair_want "$work/nopair.spec" > "$work/nopair.want"
    expect_records "$tier: -P: pairs a and b rescued, unpaired" "$work/nopair.want" "$sam" a b
    got="$(primary_fields "$sam" c 2 3 4 5 AS) $(primary_bit "$sam" c 2 2)"
    got="$got $(primary_bit "$sam" c 1 2) $(primary_fields "$sam" c 1 5)"
    mem_cli_expect "$tier: -P: pair c read 2 rescued at ctg2:1301, no 0x2, read 1 MAPQ 0" \
        "ctg2 1301 60 70 0 0 0" "$got"

    # The rescue shortcuts are output-neutral; --fast places these pairs alike.
    sam="$work/abc.refenv.$tier.sam"
    mem_env "$tier" "$sam" "${rescue_ref_env[@]}" -- -I 300,30 "$ref" "$work/abc.1.fq" "$work/abc.2.fq"
    mem_cli_same_file "$tier: rescue with every BWA3_RESCUE_* shortcut off" "$abc" "$sam"
    sam="$work/abc.fast.$tier.sam"
    mem_cli_mem "$tier" "$sam" --fast -I 300,30 "$ref" "$work/abc.1.fq" "$work/abc.2.fq"
    expect_records "$tier: --fast: records of pairs a, b, c" "$work/abc.want" "$sam"

    # An FF library: the FF window is fitted and its rescue runs forward.
    sam="$work/ff.$tier.sam"
    mem_cli_mem "$tier" "$sam" "$ref" "$work/ff.1.fq" "$work/ff.2.fq"
    expect_records "$tier: FF library: 67/131, TLEN +-301" "$work/ff.want" "$sam"
    got="$(primary_fields "$sam" fr 2 AS NM MD)"
    mem_cli_expect "$tier: FF library: read 2 of fr rescued forward at ctg3:1001, AS NM MD" "70 6 $md_fr2" "$got"

    # A split end: 0x2 from the unpaired rule, inclusive at both window ends.
    sam="$work/split.$tier.sam"
    mem_cli_mem "$tier" "$sam" -I 300,30 "$ref" "$work/split.1.fq" "$work/split.2.fq"
    awk -F '\t' '!/^@/ { print $1, $2, $3, $4, $5, $6, $7, $8, $9 }' "$sam" > "$sam.got"
    mem_cli_same_file "$tier: split read 1: unpaired 0x2 window [180, 420]" "$work/split.want" "$sam.got"

    # --compat shapes tags and @HD only.
    for target in bwa-mem2 bwa-mem; do
        sam="$work/abc.$target.$tier.sam"
        mem_cli_mem "$tier" "$sam" --compat="$target" -I 300,30 "$ref" "$work/abc.1.fq" "$work/abc.2.fq"
        if [[ $target == bwa-mem2 ]]; then
            re='^(MQ|HN):i:'
            want_hd=0
        else
            re='^HN:i:'
            want_hd=1
        fi
        strip_tags "$abc" "$re" > "$sam.want"
        mem_cli_records "$sam" > "$sam.got"
        mem_cli_same_file "$tier: --compat=$target: default records less $re" "$sam.want" "$sam.got"
        mem_cli_expect "$tier: --compat=$target: @HD lines" "$want_hd" "$(grep -c '^@HD' "$sam" || true)"
    done

    # Insert-size estimation from the data, and its 10-pair minimum.
    sam="$work/est.$tier.sam"
    mem_cli_mem "$tier" "$sam" "$ref" "$work/est.1.fq" "$work/est.2.fq"
    expect_records "$tier: estimated window [147, 483]: proper exactly inside" "$work/est.want" "$sam"
    for n in 9 10; do
        sam="$work/cnt$n.$tier.sam"
        mem_cli_mem "$tier" "$sam" "$ref" "$work/cnt$n.1.fq" "$work/cnt$n.2.fq"
        expect_records "$tier: $n pairs: proper only from 10" "$work/cnt$n.want" "$sam"
    done

    # The batch boundary, per-batch estimation, and -t at a fixed -K.
    for kt in 3000:1 3000:4 3001:1; do
        sam="$work/kb.K${kt%%:*}.t${kt##*:}.$tier.sam"
        mem_cli_mem "$tier" "$sam" -K "${kt%%:*}" -t "${kt##*:}" "$ref" "$work/kb.1.fq" "$work/kb.2.fq"
        expect_records "$tier: -K ${kt%%:*} -t ${kt##*:}: per-batch windows" "$work/kb${kt%%:*}.want" "$sam"
    done

    # Determinism across -t and --dedup-reads at a fixed -K; -K under -I.
    det="$work/all.K2000.t1.$tier.sam"
    mem_cli_mem "$tier" "$det" -K 2000 -t 1 "$ref" "$work/all.1.fq" "$work/all.2.fq"
    for args in "-t 2" "-t 4" "-t 4 -1" "-t 4 --dedup-reads off"; do
        sam="$work/all.K2000.${args// /}.$tier.sam"
        # shellcheck disable=SC2086 # $args is a list of options on purpose
        mem_cli_mem "$tier" "$sam" -K 2000 $args "$ref" "$work/all.1.fq" "$work/all.2.fq"
        mem_cli_same_file "$tier: -K 2000 $args equals -t 1" "$det" "$sam"
    done
    # Read memoization: on (its counters must show memoized duplicates), and on
    # with every duplicate aligned anyway and checked against its representative.
    sam="$work/all.K2000.memo.$tier.sam"
    mem_env "$tier" "$sam" BWAMEM3_DEDUP_READS_STATS=1 -- -K 2000 -t 4 --dedup-reads on \
        "$ref" "$work/all.1.fq" "$work/all.2.fq"
    mem_cli_same_file "$tier: -K 2000 -t 4 --dedup-reads on equals -t 1" "$det" "$sam"
    got="$(sed -n 's/^\[dedup-reads-stats\].* dup=\([0-9]*\).*/\1/p' "$sam.err")"
    if [[ -n $got ]] && ((got > 0)); then got=yes; else got="no ('$got')"; fi
    mem_cli_expect "$tier: --dedup-reads on: duplicate pairs were memoized" "yes" "$got"
    sam="$work/all.K2000.memoverify.$tier.sam"
    mem_env "$tier" "$sam" BWAMEM3_DEDUP_READS_VERIFY=1 BWAMEM3_DEDUP_READS_CHECK_BASES=1 -- \
        -K 2000 -t 4 --dedup-reads on "$ref" "$work/all.1.fq" "$work/all.2.fq"
    mem_cli_same_file "$tier: --dedup-reads on, duplicates verified, equals -t 1" "$det" "$sam"
    # Cohort slicing and output streaming are output-neutral: no slices, every
    # batch sliced (the stress knob; its slices must show on stderr), and
    # streaming off.
    sam="$work/all.K2000.noslice.$tier.sam"
    mem_cli_mem "$tier" "$sam" -K 2000 -t 4 --cohort-slices 0 "$ref" "$work/all.1.fq" "$work/all.2.fq"
    mem_cli_same_file "$tier: -K 2000 --cohort-slices 0 equals -t 1" "$det" "$sam"
    sam="$work/all.K2000.sliceall.$tier.sam"
    mem_env "$tier" "$sam" BWA_MEM3_COHORT_SLICE_ALL=1 -- -K 2000 -t 4 "$ref" "$work/all.1.fq" "$work/all.2.fq"
    mem_cli_same_file "$tier: -K 2000 with every batch sliced equals -t 1" "$det" "$sam"
    got="$(grep -c '(cohort slice)' "$sam.err" || true)"
    if ((got > 0)); then got=yes; else got="no ($got)"; fi
    mem_cli_expect "$tier: every batch sliced: slices were read" "yes" "$got"
    sam="$work/all.K2000.nostream.$tier.sam"
    mem_env "$tier" "$sam" BWA_MEM3_STREAM_OUTPUT=0 -- -K 2000 -t 4 "$ref" "$work/all.1.fq" "$work/all.2.fq"
    mem_cli_same_file "$tier: -K 2000 with output streaming off equals -t 1" "$det" "$sam"
    sam1="$work/all.I.K1.$tier.sam"
    sam2="$work/all.I.$tier.sam"
    mem_cli_mem "$tier" "$sam1" -K 1 -t 4 -I 300,30 "$ref" "$work/all.1.fq" "$work/all.2.fq"
    mem_cli_mem "$tier" "$sam2" -t 4 -I 300,30 "$ref" "$work/all.1.fq" "$work/all.2.fq"
    mem_cli_same_file "$tier: -I 300,30: -K 1 equals the default batch size" "$sam2" "$sam1"
    if pair_got "$sam1" | awk 'NR % 2 == 1 { print $1 }' | cmp -s - "$work/all.names"; then
        got=ok
    else
        got="names differ"
    fi
    mem_cli_expect "$tier: -K 1: each of the $n_all pairs once, in input order" "ok" "$got"
    got="$(pair_got "$sam1" | awk 'END { print NR }')"
    mem_cli_expect "$tier: -K 1: two primary records per pair" "$((2 * n_all))" "$got"

    # Supplementary marking and its MAPQ; -q keeps it.
    for q in default -q; do
        sam="$work/chim.$q.$tier.sam"
        if [[ $q == -q ]]; then
            mem_cli_mem "$tier" "$sam" -q "$ref" "$work/chim.fq"
            supp_mapq=60
        else
            mem_cli_mem "$tier" "$sam" "$ref" "$work/chim.fq"
            supp_mapq=0
        fi
        supp="$(awk -F '\t' '!/^@/ && int($2 / 2048) % 2' "$sam")"
        got="$(record_fields "$supp" 2 3 4 5 6 AS)"
        mem_cli_expect "$tier: $q: supplementary ctg3 part" "2048 ctg3 2501 $supp_mapq 70H50M 50" "$got"
        got="$(primary_fields "$sam" chim "" 5 AS) $(mem_cli_count "$sam" '1')"
        mem_cli_expect "$tier: $q: primary R part MAPQ 0 AS 70, two records" "0 70 2" "$got"
    done

    # -j and the ALT sidecar.
    for idx in jplain jalt jalt-j; do
        sam="$work/$idx.$tier.sam"
        if [[ $idx == jalt-j ]]; then
            mem_cli_mem "$tier" "$sam" -j "$work/jalt/ref.fa" "$work/j.fq"
        else
            mem_cli_mem "$tier" "$sam" "$work/$idx/ref.fa" "$work/j.fq"
        fi
    done
    for e in "jplain 10 95" "jalt 60 0" "jalt-j 10 95"; do
        read -r idx mapq xs <<< "$e"
        sam="$work/$idx.$tier.sam"
        got="$(primary_fields "$sam" j1 "" 3 4 5 XS)"
        mem_cli_expect "$tier: $idx: j1 RNAME POS MAPQ XS" "ctgA 1001 $mapq $xs" "$got"
    done
    mem_cli_same_file "$tier: -j equals the index without a sidecar" \
        "$work/jplain.$tier.sam" "$work/jalt-j.$tier.sam"

    # --meth: directional OT and OB pairs.
    sam="$work/meth.$tier.sam"
    mem_cli_mem "$tier" "$sam" --meth -I 300,30 "$work/meth/ref.fa" "$work/meth.1.fq" "$work/meth.2.fq"
    printf '%s\n' \
        'mCT 1 99 ctg1 1001 60 100M = 1201 300' 'mCT 2 147 ctg1 1201 60 100M = 1001 -300' \
        'mGA 1 83 ctg1 3801 60 100M = 3601 -300' 'mGA 2 163 ctg1 3601 60 100M = 3801 300' > "$work/meth.want"
    expect_records "$tier: --meth: OT and OB pair records" "$work/meth.want" "$sam"
    for e in "mCT 1 CT CT $xm_ct1" "mCT 2 GA CT $xm_ct2" "mGA 1 CT GA $xm_ga1" "mGA 2 GA GA $xm_ga2"; do
        read -r qn end xr xg xm <<< "$e"
        got="$(primary_fields "$sam" "$qn" "$end" AS NM XR XG XM)"
        mem_cli_expect "$tier: --meth: $qn end $end AS NM XR XG XM" "100 0 $xr $xg $xm" "$got"
    done

    # Every SAM above is identical across tiers.
    if [[ $tier != "$first_tier" ]]; then
        for s in i300 abc ff split est kb.K3001.t1 all.K2000.t1 all.I chim.default jalt meth; do
            mem_cli_same_file "$tier: $s SAM equals $first_tier's" \
                "$work/$s.$first_tier.sam" "$work/$s.$tier.sam"
        done
    fi
done

if ((MEM_CLI_FAILS > 0)); then
    echo "FAIL: mem_pair_cli: $MEM_CLI_FAILS of $MEM_CLI_CHECKS checks failed"
    exit 1
fi
echo "PASS: mem_pair_cli: $MEM_CLI_CHECKS checks, $MEM_CLI_SKIPS skipped"
