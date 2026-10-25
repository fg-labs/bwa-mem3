#!/usr/bin/env bash
# test/regression/mem_io_cli.sh
#
# End-to-end checks of `bwa-mem3 mem` input and output: FASTQ/FASTA readers
# (plain and gzip, -p interleaved), -C, -R/-H headers, -o, --bam output, the
# SAM record fields and tags (FLAG, CIGAR, NM, MD, AS, XS, XA, SA, MQ, HN),
# -Y, -M, -V, -5, clipping and unmapped records. Each check runs under every
# SIMD tier this host can run (test/framework/mem_cli.sh: mem_cli_tiers), on
# the synthetic reference (test/fixtures/mem_synth_ref.awk, whose header lists
# its contigs and planted repeats), with every expected value worked out from
# the reference by hand.
#
# The reads are slices of the synthetic reference, so each field follows from
# the slice and the default scoring (-A 1 -B 4 -O 6 -E 1 -L 5). A unique exact
# 100 bp hit has MAPQ 60 (bwa's mem_approx_mapq_se: 6.02 * (100 - 19) *
# (3 / ln 100)^2 = 207, with bwa's integer coefficient (int)ln 50 = 3, capped
# at 60), NM 0, MD 100, AS 100 and XS 0. Read
# set "fmt" (one record each):
#
#   exact   ctg1:1001-1100, name "exact/1", comment "BC:Z:ACGT<TAB>RX:Z:AAC".
#   rev     reverse complement of ctg1:1201-1300: FLAG 16, SEQ the reference
#           bases, QUAL reversed.
#   lower   ctg1:1401-1500 in lower case, name "lower/2": SEQ upper case.
#   iupac   ctg1:1601-1700 with read base 50 an R: SEQ shows N there, NM 1,
#           MD 49<ref>50, AS 99 - 1 = 98 (N scores -1).
#   revn    reverse complement of ctg1:2601-2700 with read base 30 an N:
#           FLAG 16, SEQ the reference with N at base 71 (100 - 30 + 1), NM 1,
#           MD 70<ref>29, AS 98.
#   name/12 ctg1:1801-1900; "/12" is not a /<digit> suffix, so the name is kept.
#   crlf    ctg1:2401-2500.
#   multi   ctg1:3601-3700.
#   qat     ctg1:3801-3900 with a quality string that starts with '@'.
#   long    ctg2:2601-3600 (1000 bp): 1000M, AS 1000.
#   longname ctg2:101-200, with a 200-character name of SAM-legal punctuation.
#   unmapped ctg1:3701-3800 reversed (not complemented), which shares no 19-mer
#           with the reference (its longest exact match is 11 bp), comment
#           "CO:Z:nohit": FLAG 4, RNAME *, POS 0, MAPQ 0, CIGAR *, SEQ and QUAL
#           as given, AS:i:0 XS:i:0 and nothing else.
#
# Checks:
#   @SQ     one line per contig, in reference order, with the generator's
#           lengths: ctg1 5000, ctg2 4000, ctg3 3000, ctg4 600.
#   fields  every field of every "fmt" record, as above, and records in input
#           order.
#   reader  the same reads written as messy FASTQ (text before the first
#           record, blank lines, CRLF, multi-line sequence and quality, a '+'
#           line repeating the name), as gzip, as gzip of several members, as
#           BGZF (when bgzip is installed), from a pipe (/dev/stdin, "-" and a
#           "<command" input; a command that exits nonzero fails the run),
#           with the reader's buffer forced down to 1, 7 and 61 bytes
#           (FR_FASTQ_TEST_BUFSZ, so every record straddles refills), and read
#           by --legacy-reader: the SAM is byte-identical to the clean FASTQ's.
#           As multi-line FASTA: the same records with QUAL *.
#   bad input  a quality string shorter than its sequence, and a gzip file cut
#           short: mem exits nonzero.
#   -C      appends the comment verbatim to the exact record and to the unmapped
#           record (after AS and XS); every other byte is unchanged.
#   -R -H   header: @HD, the four @SQ, the -H string, the @-lines of the -H file
#           (others dropped), the @RG line, @PG, in that order; RG:Z:grp1 on
#           every record. A -R line not starting with @RG, or without ID, exits
#           nonzero.
#   -o -f   the file holds what stdout would, and stdout is empty; -f is -o.
#   -v      -v 1 and -v 4 change only stderr: the SAM is the default one.
#   --bam   --bam, --bam=0, --bam=6 and --bam=6 --bam-threads 2 are
#           BGZF-wrapped BAM ("BAM\1" after gunzip); decoded with samtools
#           (when installed) they equal the SAM text, and carry bwa-mem3's
#           @PG. A malformed --bam-threads value exits nonzero.
#   split   a read made of ctg1:101-160 and ctg3:2601-2640 (forward, and with
#           the second part reverse-complemented): primary 60M40S on ctg1:101,
#           supplementary FLAG 2048 (2064 reverse) on ctg3:2601 as 60H40M
#           (40M60H reverse) with SEQ/QUAL of the 40 aligned bases only (reverse
#           complemented and reversed on the reverse strand), SA:Z naming the
#           other hit with S clips; MAPQ 60 on both.
#   -Y      the supplementary keeps S clips and the whole SEQ/QUAL.
#   -M      the SAM equals the default one with FLAG 2048 -> 256 (2064 -> 272).
#   -5      a read made of ctg3:2601-2640 then ctg1:4101-4160: by default the
#           longer ctg1 hit is primary; with -5 the ctg3 hit (the read's 5'
#           end).
#   -V      against a copy of the reference whose ctg1 header carries the
#           description "first contig<TAB>desc": XR:Z:first contig desc on the
#           ctg1 hit, none on a ctg2 hit, and none without -V.
#   XA HN   a read inside block R (ctg1:2101-2200), and its reverse complement:
#           the primary hit plus the XA entries are exactly the three copies
#           (ctg1 +2101, ctg2 +1101, ctg4 -201, strands flipped for the reverse
#           complement), each 100M NM 0, and HN:i:2; -u appends the score
#           (100) and MAPQ (0) to each XA entry. Only the entries are checked,
#           read from XA or XB: the -u help text and docs/src/cli/mem.md name
#           the tag XB, docs/src/user-guide/output.md says it stays XA, and
#           the tag name is left to that documentation to settle.
#   -p      interleaved pairs (with a single-end read between them) at
#           -I 300,30: a FR pair FLAG 99/147, TLEN +-300, MC 100M, MQ 60; a FR
#           pair whose read 1 is the reverse one, 83/163; a pair whose second
#           end is unmappable, 73/133, with the unmapped end placed at its mate
#           (POS 1001, CIGAR *); the single read FLAG 0 in its input position.
#           Without the single read, -p on one interleaved file equals the
#           two-file run byte for byte.
#   threads -t 1 and -t 4 at a fixed -K give byte-identical SAM, single-end and
#           -p, over several batches (one of which reaches -K at an odd record
#           count); records come out in input order, and every renamed copy of a
#           read or pair gets its original's records. The same for 3300 short
#           reads in 1200-read batches at -t 4 (where the process step streams
#           finished work items to the output; a "note:" line reports how many
#           reads went out early, a timing-dependent count that is not
#           checked), with streaming off (BWA_MEM3_STREAM_OUTPUT=0), and at -t 1.
#   tiers   every SAM above is byte-identical across tiers (@PG aside).
#
# Inputs (env vars):
#   BWA_MEM3       path to the bwa-mem3 binary under test
#   MEM_CLI_TIERS  (optional) tiers to run instead of every runnable one
#   MEM_CLI_KEEP   (optional) 1 keeps the work directory
# No staged inputs: the reference is generated, the reads are sliced from it.
# Optional tools: samtools (the --bam decode check) and bgzip (the BGZF input
# check); each check is reported as SKIP when its tool is missing.

# Many arguments are awk programs or awk conditions (mem_cli_count,
# same_after), so `$1`, `$NF` and the like in single quotes are awk's fields,
# never the shell's.
# shellcheck disable=SC2016
set -euo pipefail

# shellcheck source=test/framework/mem_cli.sh
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../framework" && pwd)/mem_cli.sh"

mem_cli_init mem_io_cli
mem_cli_setup_synth
ref="$MEM_CLI_SYNTH"
work="$MEM_CLI_WORK"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# qual_string <len> <salt>
#   A quality string of <len> printable characters ('#'..'J'), varied so that
#   an offset, trimming or reversal error moves a visible character.
qual_string() {
    awk -v n="$1" -v s="$2" 'BEGIN {
        for (i = 0; i < n; i++) printf "%c", 35 + (i * 7 + s) % 40
        printf "\n"
    }'
}

# reverse_string <s>
reverse_string() {
    printf '%s\n' "$1" | rev
}

# upper <s>
upper() {
    printf '%s\n' "$1" | tr 'acgtn' 'ACGTN'
}

# fq_rec <name> <comment> <seq> <qual>: one FASTQ record (comment may be "").
fq_rec() {
    if [[ -n $2 ]]; then
        printf '@%s %s\n%s\n+\n%s\n' "$1" "$2" "$3" "$4"
    else
        printf '@%s\n%s\n+\n%s\n' "$1" "$3" "$4"
    fi
}

# split_lines <s> <width>: <s> broken into lines of at most <width> characters.
split_lines() {
    printf '%s\n' "$1" | awk -v w="$2" '{
        for (i = 1; i <= length($0); i += w) print substr($0, i, w)
    }'
}

# check_fields <label> <sam> <qname> <FLAG RNAME POS MAPQ CIGAR> <SEQ> <QUAL>
#              <tags>
#   The one primary record of <qname>: columns 2-6, SEQ, QUAL and the optional
#   fields (TAB-joined) as expected.
check_fields() {
    local rec
    rec="$(mem_cli_primary "$2" "$3")"
    mem_cli_expect "$1: $3 FLAG RNAME POS MAPQ CIGAR" "$4" \
        "$(printf '%s\n' "$rec" | awk -F '\t' '{ print $2, $3, $4, $5, $6 }')"
    mem_cli_expect "$1: $3 SEQ" "$5" "$(mem_cli_col "$rec" 10)"
    mem_cli_expect "$1: $3 QUAL" "$6" "$(mem_cli_col "$rec" 11)"
    mem_cli_expect "$1: $3 tags" "$7" \
        "$(printf '%s\n' "$rec" | awk -F '\t' '{ s = $12; for (i = 13; i <= NF; i++) s = s "\t" $i; print s }')"
}

# record_of <sam> <qname> <flag>: the record of <qname> with exactly <flag>.
record_of() {
    awk -F '\t' -v q="$2" -v f="$3" '!/^@/ && $1 == q && $2 == f' "$1"
}

# overhang_clean <ref-bases> <read-bases>
#   Prints "yes" when extending an alignment base by base into <ref-bases>
#   with <read-bases> never reaches a positive running score (+1 match, -4
#   mismatch): the local extension then stops exactly at the junction, which
#   is what the hand-derived CIGARs of the split reads assume.
overhang_clean() {
    awk -v r="$1" -v q="$2" 'BEGIN {
        s = 0; ok = "yes"
        for (i = 1; i <= length(q) && i <= length(r); i++) {
            s += (substr(r, i, 1) == substr(q, i, 1)) ? 1 : -4
            if (s >= 0) ok = "no"
        }
        print ok
    }'
}

tab="$(printf '\t')"

# ---------------------------------------------------------------------------
# Reads
# ---------------------------------------------------------------------------

s_exact="$(mem_cli_seq "$ref" ctg1 1001 100)"
s_rev_fwd="$(mem_cli_seq "$ref" ctg1 1201 100)"
s_rev="$(mem_cli_revcomp "$s_rev_fwd")"
s_lower_ref="$(mem_cli_seq "$ref" ctg1 1401 100)"
s_lower="$(printf '%s\n' "$s_lower_ref" | tr 'ACGT' 'acgt')"
s_iupac_ref="$(mem_cli_seq "$ref" ctg1 1601 100)"
s_iupac="$(mem_cli_subst "$s_iupac_ref" 50 R)"
s_iupac_sam="$(mem_cli_subst "$s_iupac_ref" 50 N)"
iupac_ref_base="${s_iupac_ref:49:1}"
s_revn_ref="$(mem_cli_seq "$ref" ctg1 2601 100)"
s_revn="$(mem_cli_subst "$(mem_cli_revcomp "$s_revn_ref")" 30 N)"
s_revn_sam="$(mem_cli_subst "$s_revn_ref" 71 N)"
revn_ref_base="${s_revn_ref:70:1}"
s_slash="$(mem_cli_seq "$ref" ctg1 1801 100)"
s_crlf="$(mem_cli_seq "$ref" ctg1 2401 100)"
s_multi="$(mem_cli_seq "$ref" ctg1 3601 100)"
s_qat="$(mem_cli_seq "$ref" ctg1 3801 100)"
s_long="$(mem_cli_seq "$ref" ctg2 2601 1000)"
s_unm="$(reverse_string "$(mem_cli_seq "$ref" ctg1 3701 100)")"

q_exact="$(qual_string 100 1)"
q_rev="$(qual_string 100 2)"
q_lower="$(qual_string 100 3)"
q_iupac="$(qual_string 100 4)"
q_slash="$(qual_string 100 5)"
q_crlf="$(qual_string 100 6)"
q_multi="$(qual_string 100 7)"
q_qat="@$(qual_string 99 8)"
q_long="$(qual_string 1000 9)"
q_unm="$(qual_string 100 10)"
s_lname="$(mem_cli_seq "$ref" ctg2 101 100)"
q_lname="$(qual_string 100 13)"
q_revn="$(qual_string 100 15)"
lname="rd:1#ACGT+_=|$(printf '%188s' '' | tr ' ' 'x')"
c_exact="BC:Z:ACGT${tab}RX:Z:AAC"
c_unm="CO:Z:nohit"

# The canonical single-line FASTQ.
{
    fq_rec exact/1 "$c_exact" "$s_exact" "$q_exact"
    fq_rec rev "" "$s_rev" "$q_rev"
    fq_rec lower/2 "" "$s_lower" "$q_lower"
    fq_rec iupac "" "$s_iupac" "$q_iupac"
    fq_rec revn "" "$s_revn" "$q_revn"
    fq_rec name/12 "" "$s_slash" "$q_slash"
    fq_rec crlf "" "$s_crlf" "$q_crlf"
    fq_rec multi "" "$s_multi" "$q_multi"
    fq_rec qat "" "$s_qat" "$q_qat"
    fq_rec long "" "$s_long" "$q_long"
    fq_rec unmapped "$c_unm" "$s_unm" "$q_unm"
    fq_rec "$lname" "" "$s_lname" "$q_lname"
} > "$work/fmt.fq"
n_fmt=12

# The same records, written messily: text before the first record, blank
# lines, CRLF line ends, multi-line sequence and quality, a '+' line that
# repeats the name.
cr="$(printf '\r')"
{
    printf 'text before the first record\n\n'
    printf '@exact/1 %s%s\n%s%s\n+%s\n%s%s\n' "$c_exact" "$cr" "$s_exact" "$cr" "$cr" "$q_exact" "$cr"
    printf '\n@rev\n'
    split_lines "$s_rev" 37
    printf '+rev\n'
    split_lines "$q_rev" 50
    fq_rec lower/2 "" "$s_lower" "$q_lower"
    printf '\n\n'
    fq_rec iupac "" "$s_iupac" "$q_iupac"
    fq_rec revn "" "$s_revn" "$q_revn"
    fq_rec name/12 "" "$s_slash" "$q_slash"
    printf '@crlf%s\n%s%s\n+%s\n%s%s\n' "$cr" "$s_crlf" "$cr" "$cr" "$q_crlf" "$cr"
    printf '@multi\n'
    split_lines "$s_multi" 30
    printf '+\n'
    split_lines "$q_multi" 30
    fq_rec qat "" "$s_qat" "$q_qat"
    printf '@long\n'
    split_lines "$s_long" 80
    printf '+long\n'
    split_lines "$q_long" 333
    fq_rec unmapped "$c_unm" "$s_unm" "$q_unm"
    fq_rec "$lname" "" "$s_lname" "$q_lname"
} > "$work/messy.fq"

# Multi-line FASTA of the same reads (no quality).
awk 'NR % 4 == 1 { sub(/^@/, ">"); print; next }
     NR % 4 == 2 { for (i = 1; i <= length($0); i += 60) print substr($0, i, 60) }' \
    "$work/fmt.fq" > "$work/fmt.fa"

gzip -c "$work/fmt.fq" > "$work/fmt.fq.gz"
# Several gzip members: the messy file cut into three, each part compressed.
awk 'NR <= 9' "$work/messy.fq" | gzip -c > "$work/messy.fq.gz"
awk 'NR > 9 && NR <= 40' "$work/messy.fq" | gzip -c >> "$work/messy.fq.gz"
awk 'NR > 40' "$work/messy.fq" | gzip -c >> "$work/messy.fq.gz"
have_bgzip=0
if command -v bgzip > /dev/null 2>&1; then
    bgzip -c "$work/messy.fq" > "$work/messy.fq.bgz"
    have_bgzip=1
fi

# Malformed input: a quality string shorter than its sequence; a gzip file cut
# in half.
printf '@a\nACGTACGTAC\n+\nIIIII\n' > "$work/short_qual.fq"
gz_bytes="$(wc -c < "$work/fmt.fq.gz" | tr -d ' ')"
dd if="$work/fmt.fq.gz" of="$work/cut.fq.gz" bs=1 count="$((gz_bytes / 2))" 2> /dev/null

# Split reads. Part A ctg1:101-160, part B ctg3:2601-2640, part A5
# ctg1:4101-4160.
s_a="$(mem_cli_seq "$ref" ctg1 101 60)"
s_b="$(mem_cli_seq "$ref" ctg3 2601 40)"
s_a5="$(mem_cli_seq "$ref" ctg1 4101 60)"
q_split="$(qual_string 100 11)"
# The hand-derived CIGARs assume each part's alignment ends exactly at the
# junction; check that against the reference bases beyond each part.
# Each entry: label|reference bases beyond the part|read bases beyond it, both
# read outwards from the junction.
for pre in \
    "chim: ctg1 hit into B|$(mem_cli_seq "$ref" ctg1 161 40)|$s_b" \
    "chim: ctg3 hit into A|$(reverse_string "$(mem_cli_seq "$ref" ctg3 2541 60)")|$(reverse_string "$s_a")" \
    "chimrev: ctg1 hit into rc(B)|$(mem_cli_seq "$ref" ctg1 161 40)|$(mem_cli_revcomp "$s_b")" \
    "chimrev: ctg3 hit into rc(A)|$(mem_cli_seq "$ref" ctg3 2641 60)|$(mem_cli_revcomp "$s_a")" \
    "chim5: ctg3 hit into A5|$(mem_cli_seq "$ref" ctg3 2641 60)|$s_a5" \
    "chim5: ctg1 hit into B|$(reverse_string "$(mem_cli_seq "$ref" ctg1 4061 40)")|$(reverse_string "$s_b")"; do
    label="${pre%%|*}"
    rest="${pre#*|}"
    if [[ $(overhang_clean "${rest%%|*}" "${rest#*|}") != yes ]]; then
        echo "FAIL: mem_io_cli: split-read precondition '$label' does not hold for this reference" >&2
        exit 1
    fi
done
{
    fq_rec chim "" "$s_a$s_b" "$q_split"
    fq_rec chimrev "" "$s_a$(mem_cli_revcomp "$s_b")" "$q_split"
    fq_rec chim5 "" "$s_b$s_a5" "$q_split"
} > "$work/split.fq"

# Reads for the XA/HN checks: inside block R (three copies), and unique.
s_inR="$(mem_cli_seq "$ref" ctg1 2101 100)"
{
    fq_rec inR "" "$s_inR" "$(qual_string 100 12)"
    fq_rec inRrc "" "$(mem_cli_revcomp "$s_inR")" "$(qual_string 100 14)"
    fq_rec uniq "" "$s_exact" "$q_exact"
} > "$work/xa.fq"

# hits <sam> <qname>: the primary hit and its XA entries as
# "contig,strand+POS,CIGAR,NM", sorted, space-separated.
hits() {
    local rec strand=+
    rec="$(mem_cli_primary "$1" "$2")"
    if (($(mem_cli_col "$rec" 2) / 16 % 2 == 1)); then strand=-; fi
    {
        printf '%s,%s%s,%s,%s;' "$(mem_cli_col "$rec" 3)" "$strand" \
            "$(mem_cli_col "$rec" 4)" "$(mem_cli_col "$rec" 6)" "$(mem_cli_tag "$rec" NM)"
        mem_cli_tag "$rec" XA
    } | tr ';' '\n' | awk 'NF' | LC_ALL=C sort | tr '\n' ' '
}

# -V: a copy of the reference with descriptions in two contig headers.
vref="$work/vref/vref.fa"
mkdir -p "$work/vref"
awk -v t="$tab" '/^>ctg1$/ { print ">ctg1 first contig" t "desc"; next }
                 /^>ctg3$/ { print ">ctg3 third"; next }
                 { print }' "$ref" > "$vref"
mem_cli_index "$vref"
{
    fq_rec on1 "" "$(mem_cli_seq "$ref" ctg1 1001 100)" "$q_exact"
    fq_rec on2 "" "$(mem_cli_seq "$ref" ctg2 101 100)" "$q_exact"
} > "$work/v.fq"

# -H file: two header lines with a non-header line between them.
printf '@CO\tfrom file 1\nnot a header line\n@CO\tfrom file 2\n' > "$work/hdr.txt"

# Pairs. pp: FR, read 1 forward at ctg1:1001, read 2 the reverse complement of
# ctg1:1201-1300. qq: read 1 the reverse complement of ctg2:3201-3300, read 2
# forward at ctg2:3001. pu: read 1 at ctg1:1001, read 2 unmappable.
s_pp2="$(mem_cli_revcomp "$s_rev_fwd")"
s_qq1="$(mem_cli_revcomp "$(mem_cli_seq "$ref" ctg2 3201 100)")"
s_qq2="$(mem_cli_seq "$ref" ctg2 3001 100)"
s_single="$(mem_cli_seq "$ref" ctg1 4001 100)"
{
    fq_rec pp/1 "" "$s_exact" "$q_exact"
    fq_rec pp/2 "" "$s_pp2" "$q_rev"
    fq_rec single "" "$s_single" "$q_lower"
    fq_rec qq/1 "" "$s_qq1" "$q_iupac"
    fq_rec qq/2 "" "$s_qq2" "$q_slash"
    fq_rec pu/1 "" "$s_exact" "$q_crlf"
    fq_rec pu/2 "" "$s_unm" "$q_unm"
} > "$work/pairs.fq"
awk 'NR < 9 || NR > 12' "$work/pairs.fq" > "$work/pairs_nosingle.fq"
awk 'int((NR - 1) / 4) % 2 == 0' "$work/pairs_nosingle.fq" > "$work/pairs_r1.fq"
awk 'int((NR - 1) / 4) % 2 == 1' "$work/pairs_nosingle.fq" > "$work/pairs_r2.fq"

# Many reads, for the thread and batching checks: 30 renamed copies of the
# fmt reads (single-end) and of the pairs.
: > "$work/many.fq"
: > "$work/many_pairs.fq"
i=0
while ((i < 30)); do
    awk -v k="$i" 'NR % 4 == 1 { n = $1; sub(/^@/, "", n); sub(/\/[0-9]+$/, "", n); print "@" n "_" k; next } { print }' \
        "$work/fmt.fq" >> "$work/many.fq"
    awk -v k="$i" 'NR % 4 == 1 { n = $1; sub(/^@/, "", n); e = ""; if (n ~ /\/[12]$/) { e = substr(n, length(n) - 1); n = substr(n, 1, length(n) - 2) } print "@" n "_" k e; next } { print }' \
        "$work/pairs_nosingle.fq" >> "$work/many_pairs.fq"
    i=$((i + 1))
done

# Many short reads, for output streaming: 300 renamed copies of the eleven
# 100 bp fmt reads (all but "long"). At -K 120000 a batch holds 1200 reads,
# more than one of the work items (512 reads; 1024 on Apple silicon) whose
# records the process step may write as soon as they and every earlier item
# are done (SAM text, single-end, -t > 1).
awk 'NR % 4 == 1 { keep = ($1 != "@long") } keep' "$work/fmt.fq" > "$work/short.fq"
awk -v n=300 '{ line[NR] = $0 }
    END {
        for (k = 0; k < n; k++)
            for (i = 1; i <= NR; i++) {
                if (i % 4 != 1) { print line[i]; continue }
                nm = line[i]
                sub(/[ \t].*/, "", nm)
                sub(/\/[0-9]+$/, "", nm)
                print nm "_" k
            }
    }' "$work/short.fq" > "$work/stream.fq"

# ---------------------------------------------------------------------------
# Checks, under each tier
# ---------------------------------------------------------------------------

have_samtools=0
if command -v samtools > /dev/null 2>&1; then have_samtools=1; fi

# Every SAM named here is compared across tiers at the end.
sams="default messy fasta gz messy_gz legacy legacy_gz stdin dash kpipe c_opt rh o_opt f_opt v1 v4 split
      y_opt m_opt p5 v_opt nov xa u_opt pairs pairs_nosingle pairs_2file many_t1 many_t4 manyp_t1
      manyp_t4 stream_t4"
if ((have_bgzip)); then sams="$sams messy_bgz"; fi

# out <name>: the SAM path of run <name> under the current tier.
out() {
    printf '%s/%s.%s.sam' "$work" "$1" "$tier"
}

# same_after <label> <file-a> <awk-a> <file-b> <awk-b>
#   The records (no header) of two SAM files, each passed through its own awk
#   program, must be byte-identical; a failure shows the diff.
same_after() {
    mem_cli_records "$2" | awk -F '\t' -v OFS='\t' "$3" > "$work/same_after.a"
    mem_cli_records "$4" | awk -F '\t' -v OFS='\t' "$5" > "$work/same_after.b"
    mem_cli_same_file "$1" "$work/same_after.a" "$work/same_after.b"
}

first_tier=""
for tier in $(mem_cli_tiers); do

    # --- @SQ and the record fields -------------------------------------------
    sam="$(out default)"
    mem_cli_mem "$tier" "$sam" "$ref" "$work/fmt.fq"
    mem_cli_expect "$tier: @SQ lines" \
        "SN:ctg1 LN:5000;SN:ctg2 LN:4000;SN:ctg3 LN:3000;SN:ctg4 LN:600;" \
        "$(awk -F '\t' '$1 == "@SQ" { printf "%s %s;", $2, $3 }' "$sam")"
    t100="NM:i:0${tab}MD:Z:100${tab}AS:i:100${tab}XS:i:0${tab}HN:i:0"
    check_fields "$tier" "$sam" exact "0 ctg1 1001 60 100M" "$s_exact" "$q_exact" "$t100"
    check_fields "$tier" "$sam" rev "16 ctg1 1201 60 100M" "$s_rev_fwd" "$(reverse_string "$q_rev")" "$t100"
    check_fields "$tier" "$sam" lower "0 ctg1 1401 60 100M" "$s_lower_ref" "$q_lower" "$t100"
    check_fields "$tier" "$sam" iupac "0 ctg1 1601 60 100M" "$s_iupac_sam" "$q_iupac" \
        "NM:i:1${tab}MD:Z:49${iupac_ref_base}50${tab}AS:i:98${tab}XS:i:0${tab}HN:i:0"
    check_fields "$tier" "$sam" revn "16 ctg1 2601 60 100M" "$s_revn_sam" "$(reverse_string "$q_revn")" \
        "NM:i:1${tab}MD:Z:70${revn_ref_base}29${tab}AS:i:98${tab}XS:i:0${tab}HN:i:0"
    check_fields "$tier" "$sam" name/12 "0 ctg1 1801 60 100M" "$s_slash" "$q_slash" "$t100"
    check_fields "$tier" "$sam" crlf "0 ctg1 2401 60 100M" "$s_crlf" "$q_crlf" "$t100"
    check_fields "$tier" "$sam" multi "0 ctg1 3601 60 100M" "$s_multi" "$q_multi" "$t100"
    check_fields "$tier" "$sam" qat "0 ctg1 3801 60 100M" "$s_qat" "$q_qat" "$t100"
    check_fields "$tier" "$sam" long "0 ctg2 2601 60 1000M" "$s_long" "$q_long" \
        "NM:i:0${tab}MD:Z:1000${tab}AS:i:1000${tab}XS:i:0${tab}HN:i:0"
    check_fields "$tier" "$sam" unmapped "4 * 0 0 *" "$s_unm" "$q_unm" "AS:i:0${tab}XS:i:0"
    check_fields "$tier" "$sam" "$lname" "0 ctg2 101 60 100M" "$s_lname" "$q_lname" "$t100"
    mem_cli_expect "$tier: one record per read, in input order" \
        "exact rev lower iupac revn name/12 crlf multi qat long unmapped $lname " \
        "$(mem_cli_records "$sam" | awk -F '\t' '{ printf "%s ", $1 }')"
    mem_cli_expect "$tier: RNEXT PNEXT TLEN of single-end records" "$n_fmt" \
        "$(mem_cli_count "$sam" '$7 == "*" && $8 == 0 && $9 == 0')"

    # --- Reader equivalences ---------------------------------------------------
    mem_cli_mem "$tier" "$(out messy)" "$ref" "$work/messy.fq"
    mem_cli_same_file "$tier: messy FASTQ gives the clean FASTQ's SAM" "$sam" "$(out messy)"
    mem_cli_mem "$tier" "$(out gz)" "$ref" "$work/fmt.fq.gz"
    mem_cli_same_file "$tier: gzip FASTQ gives the plain FASTQ's SAM" "$sam" "$(out gz)"
    mem_cli_mem "$tier" "$(out messy_gz)" "$ref" "$work/messy.fq.gz"
    mem_cli_same_file "$tier: multi-member gzip gives the plain FASTQ's SAM" "$sam" "$(out messy_gz)"
    if ((have_bgzip)); then
        mem_cli_mem "$tier" "$(out messy_bgz)" "$ref" "$work/messy.fq.bgz"
        mem_cli_same_file "$tier: BGZF gives the plain FASTQ's SAM" "$sam" "$(out messy_bgz)"
    else
        mem_cli_skip "$tier: BGZF input" "bgzip not installed"
    fi
    mem_cli_mem "$tier" "$(out legacy)" --legacy-reader "$ref" "$work/messy.fq"
    mem_cli_same_file "$tier: --legacy-reader gives the default reader's SAM" "$sam" "$(out legacy)"
    mem_cli_mem "$tier" "$(out legacy_gz)" --legacy-reader "$ref" "$work/messy.fq.gz"
    mem_cli_same_file "$tier: --legacy-reader on multi-member gzip" "$sam" "$(out legacy_gz)"
    mem_cli_mem "$tier" "$(out fasta)" "$ref" "$work/fmt.fa"
    same_after "$tier: FASTA gives the FASTQ records with QUAL *" \
        "$sam" '{ $11 = "*"; print }' "$(out fasta)" '{ print }'
    mem_cli_expect "$tier: FASTA run has every record" "$n_fmt" "$(mem_cli_count "$(out fasta)" '$11 == "*"')"
    mem_cli_mem "$tier" "$(out stdin)" "$ref" /dev/stdin < "$work/messy.fq.gz"
    mem_cli_same_file "$tier: gzip from a pipe (/dev/stdin)" "$sam" "$(out stdin)"
    mem_cli_mem "$tier" "$(out dash)" "$ref" - < "$work/messy.fq"
    mem_cli_same_file "$tier: FASTQ from a pipe (-)" "$sam" "$(out dash)"
    mem_cli_mem "$tier" "$(out kpipe)" "$ref" "<gzip -dc '$work/messy.fq.gz'"
    mem_cli_same_file "$tier: FASTQ from a '<command' input" "$sam" "$(out kpipe)"
    mem_cli_mem_status "$tier" "$work/kpipe_fail.$tier.sam" "$ref" "<cat '$work/fmt.fq'; exit 3"
    mem_cli_expect "$tier: a '<command' input that exits nonzero fails the run" nonzero \
        "$(if ((MEM_CLI_STATUS != 0)); then echo nonzero; else echo zero; fi)"
    for bufsz in 1 7 61; do
        for f in messy.fq messy.fq.gz; do
            # Exported and unset around the call rather than prefixed to it:
            # some shells keep a prefix assignment after a function call.
            export FR_FASTQ_TEST_BUFSZ="$bufsz"
            mem_cli_mem "$tier" "$work/bufsz.$tier.sam" "$ref" "$work/$f"
            unset FR_FASTQ_TEST_BUFSZ
            mem_cli_same_file "$tier: $f with a $bufsz-byte reader buffer" "$sam" "$work/bufsz.$tier.sam"
        done
    done

    for bad in short_qual.fq cut.fq.gz; do
        mem_cli_mem_status "$tier" "$work/bad.$tier.sam" "$ref" "$work/$bad"
        mem_cli_expect "$tier: $bad makes mem exit nonzero" nonzero \
            "$(if ((MEM_CLI_STATUS != 0)); then echo nonzero; else echo zero; fi)"
    done

    # --- -C ------------------------------------------------------------------
    mem_cli_mem "$tier" "$(out c_opt)" -C "$ref" "$work/fmt.fq"
    mem_cli_expect "$tier: -C appends the comment to the exact record" \
        "${t100}${tab}${c_exact}" \
        "$(mem_cli_primary "$(out c_opt)" exact | awk -F '\t' '{ s = $12; for (i = 13; i <= NF; i++) s = s "\t" $i; print s }')"
    mem_cli_expect "$tier: -C appends the comment to the unmapped record" \
        "AS:i:0${tab}XS:i:0${tab}${c_unm}" \
        "$(mem_cli_primary "$(out c_opt)" unmapped | awk -F '\t' '{ s = $12; for (i = 13; i <= NF; i++) s = s "\t" $i; print s }')"
    same_after "$tier: -C changes nothing else" "$sam" '{ print }' "$(out c_opt)" \
        '$1 == "exact" || $1 == "unmapped" { n = ($1 == "exact") ? 2 : 1; s = $1; for (i = 2; i <= NF - n; i++) s = s OFS $i; print s; next } { print }'

    # --- -R, -H and the header -------------------------------------------------
    mem_cli_mem "$tier" "$(out rh)" -R '@RG\tID:grp1\tSM:s1' -H '@CO\tfrom string' -H "$work/hdr.txt" \
        "$ref" "$work/fmt.fq"
    mem_cli_expect "$tier: header record order" \
        "@HD@SQ@SQ@SQ@SQ@CO@CO@CO@RG@PG" \
        "$(awk -F '\t' '/^@/ { printf "%s", $1 }' "$(out rh)")"
    mem_cli_expect "$tier: -H lines and the @RG line" \
        "@CO from string|@CO from file 1|@CO from file 2|@RG ID:grp1 SM:s1|" \
        "$(awk -F '\t' '$1 == "@CO" || $1 == "@RG" { $1 = $1; printf "%s|", $0 }' "$(out rh)")"
    mem_cli_expect "$tier: @HD is the default one" "@HD${tab}VN:1.5${tab}SO:unsorted${tab}GO:query" \
        "$(head -n 1 "$(out rh)")"
    mem_cli_expect "$tier: RG:Z:grp1 after XS on every record" "$n_fmt" \
        "$(mem_cli_count "$(out rh)" '$0 ~ /\tXS:i:[0-9]+\tRG:Z:grp1(\t|$)/')"
    same_after "$tier: -R changes nothing else in the records" "$sam" '{ print }' "$(out rh)" \
        '{ s = $1; for (i = 2; i <= NF; i++) if ($i != "RG:Z:grp1") s = s OFS $i; print s }'
    for bad_rg in 'RG\tID:x' '@RG\tSM:x'; do
        mem_cli_mem_status "$tier" "$work/bad_rg.$tier.sam" -R "$bad_rg" "$ref" "$work/fmt.fq"
        mem_cli_expect "$tier: -R '$bad_rg' is refused" nonzero \
            "$(if ((MEM_CLI_STATUS != 0)); then echo nonzero; else echo zero; fi)"
    done

    # --- -o ------------------------------------------------------------------
    mem_cli_mem "$tier" "$work/o_stdout.$tier.txt" -o "$(out o_opt)" "$ref" "$work/fmt.fq"
    mem_cli_same_file "$tier: -o writes what stdout would" "$sam" "$(out o_opt)"
    mem_cli_expect "$tier: -o leaves stdout empty" 0 "$(wc -c < "$work/o_stdout.$tier.txt" | tr -d ' ')"
    mem_cli_mem "$tier" "$work/f_stdout.$tier.txt" -f "$(out f_opt)" "$ref" "$work/fmt.fq"
    mem_cli_same_file "$tier: -f is -o" "$sam" "$(out f_opt)"

    # --- -v ------------------------------------------------------------------
    for v in 1 4; do
        mem_cli_mem "$tier" "$(out "v$v")" -v "$v" "$ref" "$work/fmt.fq"
        mem_cli_same_file "$tier: -v $v changes only stderr" "$sam" "$(out "v$v")"
    done

    # --- --bam -----------------------------------------------------------------
    mem_cli_mem "$tier" "$work/rhc.$tier.sam" -C -R '@RG\tID:grp1\tSM:s1' "$ref" "$work/fmt.fq"
    for bam_opt in --bam --bam=0 --bam=6 "--bam=6 --bam-threads 2"; do
        bam="$work/out.$tier.bam"
        # Word splitting of bam_opt is intended: the last entry is two options.
        # shellcheck disable=SC2086
        mem_cli_mem "$tier" "$bam" $bam_opt -C -R '@RG\tID:grp1\tSM:s1' "$ref" "$work/fmt.fq"
        gzip -dc < "$bam" > "$bam.raw"
        mem_cli_expect "$tier: $bam_opt writes BGZF-wrapped BAM" "BAM 001" \
            "$(od -An -c -N 4 "$bam.raw" | awk '{ print $1 $2 $3, $4 }')"
        if ((have_samtools)); then
            samtools view -h "$bam" > "$work/bam.$tier.sam"
            mem_cli_same_file "$tier: $bam_opt decodes to the SAM text" "$work/rhc.$tier.sam" "$work/bam.$tier.sam"
            # same_file sets @PG aside (its CL: is the argv); the BAM header
            # must still carry bwa-mem3's one @PG (samtools adds its own).
            mem_cli_expect "$tier: $bam_opt header has bwa-mem3's @PG" "1" \
                "$(awk -F '\t' '$1 == "@PG" && $2 == "ID:bwa-mem3" && $3 == "PN:bwa-mem3" { n++ } END { print n + 0 }' "$work/bam.$tier.sam")"
        else
            mem_cli_skip "$tier: $bam_opt decodes to the SAM text" "samtools not installed"
        fi
    done
    for bad_bt in abc -1; do
        mem_cli_mem_status "$tier" "$work/bad_bt.$tier.bam" --bam=6 --bam-threads "$bad_bt" "$ref" "$work/fmt.fq"
        mem_cli_expect "$tier: --bam-threads $bad_bt is refused" nonzero \
            "$(if ((MEM_CLI_STATUS != 0)); then echo nonzero; else echo zero; fi)"
    done

    # --- Split reads: SA, hard and soft clips, -Y, -M, -5 ----------------------
    rc_q="$(reverse_string "$q_split")"
    sp="$(out split)"
    mem_cli_mem "$tier" "$sp" "$ref" "$work/split.fq"
    mem_cli_expect "$tier: chim primary" \
        "chim 0 ctg1 101 60 60M40S $s_a$s_b $q_split SA:Z:ctg3,2601,+,60S40M,60,0;" \
        "$(record_of "$sp" chim 0 | awk -F '\t' '{ print $1, $2, $3, $4, $5, $6, $10, $11, $16 }')"
    mem_cli_expect "$tier: chim supplementary (hard-clipped SEQ/QUAL)" \
        "chim 2048 ctg3 2601 60 60H40M $s_b ${q_split:60} NM:i:0 MD:Z:40 AS:i:40 SA:Z:ctg1,101,+,60M40S,60,0;" \
        "$(record_of "$sp" chim 2048 | awk -F '\t' '{ print $1, $2, $3, $4, $5, $6, $10, $11, $12, $13, $14, $16 }')"
    mem_cli_expect "$tier: chimrev primary" \
        "chimrev 0 ctg1 101 60 60M40S SA:Z:ctg3,2601,-,40M60S,60,0;" \
        "$(record_of "$sp" chimrev 0 | awk -F '\t' '{ print $1, $2, $3, $4, $5, $6, $16 }')"
    mem_cli_expect "$tier: chimrev reverse supplementary" \
        "chimrev 2064 ctg3 2601 60 40M60H $s_b $(reverse_string "${q_split:60}") SA:Z:ctg1,101,+,60M40S,60,0;" \
        "$(record_of "$sp" chimrev 2064 | awk -F '\t' '{ print $1, $2, $3, $4, $5, $6, $10, $11, $16 }')"
    mem_cli_expect "$tier: chim5 primary is the longer hit" \
        "chim5 0 ctg1 4101 60 40S60M|chim5 2048 ctg3 2601 60 40M60H $s_b ${q_split:0:40}" \
        "$(record_of "$sp" chim5 0 | awk -F '\t' '{ printf "%s %s %s %s %s %s|", $1, $2, $3, $4, $5, $6 }')$(record_of "$sp" chim5 2048 | awk -F '\t' '{ print $1, $2, $3, $4, $5, $6, $10, $11 }')"

    mem_cli_mem "$tier" "$(out y_opt)" -Y "$ref" "$work/split.fq"
    mem_cli_expect "$tier: -Y forward supplementary keeps S clips and the whole read" \
        "chim 2048 ctg3 2601 60S40M $s_a$s_b $q_split" \
        "$(record_of "$(out y_opt)" chim 2048 | awk -F '\t' '{ print $1, $2, $3, $4, $6, $10, $11 }')"
    mem_cli_expect "$tier: -Y reverse supplementary" \
        "chimrev 2064 ctg3 2601 40M60S $(mem_cli_revcomp "$s_a$(mem_cli_revcomp "$s_b")") $rc_q" \
        "$(record_of "$(out y_opt)" chimrev 2064 | awk -F '\t' '{ print $1, $2, $3, $4, $6, $10, $11 }')"
    same_after "$tier: -Y changes only the supplementary CIGAR, SEQ and QUAL" \
        "$sp" '$2 >= 2048 { $6 = $10 = $11 = "" } { print }' \
        "$(out y_opt)" '$2 >= 2048 { $6 = $10 = $11 = "" } { print }'

    mem_cli_mem "$tier" "$(out m_opt)" -M "$ref" "$work/split.fq"
    same_after "$tier: -M marks the shorter hits 256 and changes nothing else" \
        "$sp" '$2 >= 2048 { $2 -= 1792 } { print }' "$(out m_opt)" '{ print }'
    mem_cli_expect "$tier: -M output has the 256 records" 3 "$(mem_cli_count "$(out m_opt)" 'int($2 / 256) % 2 == 1')"

    mem_cli_mem "$tier" "$(out p5)" -5 "$ref" "$work/split.fq"
    mem_cli_expect "$tier: -5 makes the 5' hit primary" \
        "chim5 0 ctg3 2601 60 40M60S|chim5 2048 ctg1 4101 60 40H60M $s_a5 ${q_split:40}" \
        "$(record_of "$(out p5)" chim5 0 | awk -F '\t' '{ printf "%s %s %s %s %s %s|", $1, $2, $3, $4, $5, $6 }')$(record_of "$(out p5)" chim5 2048 | awk -F '\t' '{ print $1, $2, $3, $4, $5, $6, $10, $11 }')"
    same_after "$tier: -5 leaves reads whose 5' hit is primary alone" \
        "$sp" '$1 != "chim5"' "$(out p5)" '$1 != "chim5"'

    # --- -V --------------------------------------------------------------------
    mem_cli_mem "$tier" "$(out v_opt)" -V "$vref" "$work/v.fq"
    mem_cli_expect "$tier: -V XR from the ctg1 description, TAB as space" "first contig desc" \
        "$(mem_cli_tag "$(mem_cli_primary "$(out v_opt)" on1)" XR)"
    mem_cli_expect "$tier: -V ctg2 has no description, no XR" 0 \
        "$(mem_cli_count "$(out v_opt)" '$1 == "on2" && $0 ~ /\tXR:Z:/')"
    mem_cli_expect "$tier: XR is the last field" "XR:Z:first contig desc" \
        "$(mem_cli_primary "$(out v_opt)" on1 | awk -F '\t' '{ print $NF }')"
    mem_cli_mem "$tier" "$(out nov)" "$vref" "$work/v.fq"
    mem_cli_expect "$tier: without -V no XR" 0 "$(mem_cli_count "$(out nov)" '$0 ~ /\tXR:Z:/')"
    same_after "$tier: -V adds only XR" "$(out nov)" '{ print }' "$(out v_opt)" \
        '{ s = $1; for (i = 2; i <= NF; i++) if (substr($i, 1, 5) != "XR:Z:") s = s OFS $i; print s }'

    # --- XA, HN, -u ------------------------------------------------------------
    mem_cli_mem "$tier" "$(out xa)" "$ref" "$work/xa.fq"
    # R is at ctg1:2001 and ctg2:1001 (+) and reverse-complemented at ctg4:101;
    # the read's R offsets 100-199 sit at ctg4:201-300 on the reverse strand.
    mem_cli_expect "$tier: inR primary + XA hits" \
        "ctg1,+2101,100M,0 ctg2,+1101,100M,0 ctg4,-201,100M,0 " "$(hits "$(out xa)" inR)"
    mem_cli_expect "$tier: inRrc primary + XA hits (strands flipped)" \
        "ctg1,-2101,100M,0 ctg2,-1101,100M,0 ctg4,+201,100M,0 " "$(hits "$(out xa)" inRrc)"
    mem_cli_expect "$tier: HN counts the inR read's two other hits" 2 \
        "$(mem_cli_tag "$(mem_cli_primary "$(out xa)" inR)" HN)"
    mem_cli_expect "$tier: a unique read has HN 0 and no XA" "0|" \
        "$(mem_cli_tag "$(mem_cli_primary "$(out xa)" uniq)" HN)|$(mem_cli_tag "$(mem_cli_primary "$(out xa)" uniq)" XA)"
    mem_cli_mem "$tier" "$(out u_opt)" -u "$ref" "$work/xa.fq"
    xa_default="$(mem_cli_tag "$(mem_cli_primary "$(out xa)" inR)" XA)"
    u_rec="$(mem_cli_primary "$(out u_opt)" inR)"
    # The entries under whichever of XA or XB the record carries (see the header).
    xa_u="$(mem_cli_tag "$u_rec" XA)$(mem_cli_tag "$u_rec" XB)"
    mem_cli_expect "$tier: -u entries are the XA entries plus score 100 and MAPQ 0" \
        "$(printf '%s\n' "$xa_default" | sed 's/;/,100,0;/g')" "$xa_u"
    mem_cli_expect "$tier: -u leaves the record otherwise unchanged" \
        "$(mem_cli_primary "$(out xa)" inR | awk -F '\t' '{ print $1, $2, $3, $4, $5, $6, $10, $11 }')" \
        "$(printf '%s\n' "$u_rec" | awk -F '\t' '{ print $1, $2, $3, $4, $5, $6, $10, $11 }')"

    # --- -p and pairs ------------------------------------------------------------
    pr="$(out pairs)"
    mem_cli_mem "$tier" "$pr" -p -I 300,30 "$ref" "$work/pairs.fq"
    mem_cli_expect "$tier: -p records in input order" \
        "pp 99|pp 147|single 0|qq 83|qq 163|pu 73|pu 133|" \
        "$(mem_cli_records "$pr" | awk -F '\t' '{ printf "%s %s|", $1, $2 }')"
    mem_cli_expect "$tier: pp end 1 mate fields and tags" \
        "ctg1 1001 60 100M = 1201 300 MC:Z:100M MQ:i:60" \
        "$(record_of "$pr" pp 99 | awk -F '\t' '{ print $3, $4, $5, $6, $7, $8, $9, $14, $15 }')"
    mem_cli_expect "$tier: pp end 2 mate fields, SEQ and QUAL" \
        "ctg1 1201 60 100M = 1001 -300 $s_rev_fwd $(reverse_string "$q_rev")" \
        "$(record_of "$pr" pp 147 | awk -F '\t' '{ print $3, $4, $5, $6, $7, $8, $9, $10, $11 }')"
    mem_cli_expect "$tier: qq (read 1 reverse) mate fields" \
        "ctg2 3201 = 3001 -300|ctg2 3001 = 3201 300" \
        "$(record_of "$pr" qq 83 | awk -F '\t' '{ printf "%s %s %s %s %s|", $3, $4, $7, $8, $9 }')$(record_of "$pr" qq 163 | awk -F '\t' '{ print $3, $4, $7, $8, $9 }')"
    mem_cli_expect "$tier: the single read is single-end" \
        "ctg1 4001 60 100M * 0 0 $s_single" \
        "$(record_of "$pr" single 0 | awk -F '\t' '{ print $3, $4, $5, $6, $7, $8, $9, $10 }')"
    mem_cli_expect "$tier: pu mapped end with an unmapped mate" \
        "ctg1 1001 60 100M = 1001 0 NM:i:0 MD:Z:100 MQ:i:0 AS:i:100 XS:i:0 HN:i:0" \
        "$(record_of "$pr" pu 73 | awk -F '\t' '{ print $3, $4, $5, $6, $7, $8, $9, $12, $13, $14, $15, $16, $17 }')"
    mem_cli_expect "$tier: pu unmapped end placed at its mate" \
        "ctg1 1001 0 * = 1001 0 $s_unm $q_unm MC:Z:100M MQ:i:60 AS:i:0 XS:i:0" \
        "$(record_of "$pr" pu 133 | awk -F '\t' '{ print $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15 }')"
    mem_cli_mem "$tier" "$(out pairs_nosingle)" -p -I 300,30 "$ref" "$work/pairs_nosingle.fq"
    mem_cli_mem "$tier" "$(out pairs_2file)" -I 300,30 "$ref" "$work/pairs_r1.fq" "$work/pairs_r2.fq"
    mem_cli_same_file "$tier: -p interleaved equals two files" "$(out pairs_2file)" "$(out pairs_nosingle)"
    same_after "$tier: -p without the single read keeps the pairs' records" \
        "$pr" '$1 != "single"' "$(out pairs_nosingle)" '{ print }'

    # --- Threads and batches -------------------------------------------------------
    mem_cli_mem "$tier" "$(out many_t1)" -t 1 -K 10050 "$ref" "$work/many.fq"
    mem_cli_mem "$tier" "$(out many_t4)" -t 4 -K 10050 "$ref" "$work/many.fq"
    mem_cli_same_file "$tier: -t 1 and -t 4 at -K 10050 (single-end)" "$(out many_t1)" "$(out many_t4)"
    awk 'NR % 4 == 1 { sub(/^@/, "", $1); print $1 }' "$work/many.fq" > "$work/many.names"
    mem_cli_records "$(out many_t4)" | awk -F '\t' '{ print $1 }' > "$work/many_t4.names"
    mem_cli_same_file "$tier: many reads come out in input order" "$work/many.names" "$work/many_t4.names"
    # Copy k of read r is named r_k; with the name stripped of _k each copy's
    # record is its original's (the copies carry no comment, as the run has no -C).
    i=0
    : > "$work/many.expected"
    while ((i < 30)); do
        mem_cli_records "$sam" | awk -F '\t' -v OFS='\t' -v k="$i" '{ sub(/\/[0-9]+$/, "", $1); $1 = $1 "_" k; print }' >> "$work/many.expected"
        i=$((i + 1))
    done
    mem_cli_records "$(out many_t4)" > "$work/many_t4.records"
    mem_cli_same_file "$tier: every copy maps as its original" "$work/many.expected" "$work/many_t4.records"
    mem_cli_mem "$tier" "$(out manyp_t1)" -t 1 -K 10050 -p -I 300,30 "$ref" "$work/many_pairs.fq"
    mem_cli_mem "$tier" "$(out manyp_t4)" -t 4 -K 10050 -p -I 300,30 "$ref" "$work/many_pairs.fq"
    mem_cli_same_file "$tier: -t 1 and -t 4 at -K 10050 (-p)" "$(out manyp_t1)" "$(out manyp_t4)"
    awk 'NR % 4 == 1 { sub(/^@/, "", $1); sub(/\/[12]$/, "", $1); print $1 }' "$work/many_pairs.fq" \
        > "$work/manyp.names"
    mem_cli_records "$(out manyp_t4)" | awk -F '\t' '{ print $1 }' > "$work/manyp_t4.names"
    mem_cli_same_file "$tier: -p many pairs in input order" "$work/manyp.names" "$work/manyp_t4.names"
    # 100 bp reads and -K 10050: the first batch reaches the chunk size at an
    # odd record count, so it must take one more record to keep pair 51 whole.
    i=0
    : > "$work/manyp.expected"
    while ((i < 30)); do
        mem_cli_records "$(out pairs_nosingle)" | awk -F '\t' -v OFS='\t' -v k="$i" '{ $1 = $1 "_" k; print }' >> "$work/manyp.expected"
        i=$((i + 1))
    done
    mem_cli_records "$(out manyp_t4)" > "$work/manyp_t4.records"
    mem_cli_same_file "$tier: every pair copy maps as its original" "$work/manyp.expected" "$work/manyp_t4.records"

    # Streaming: -t 4 (records may be written while the batch is still being
    # processed), the same with streaming turned off, and -t 1 all give each
    # copy its original's records, in input order.
    mem_cli_records "$sam" | awk -F '\t' -v n=300 '$1 != "long" {
            nm = $1
            sub(/\/[0-9]+$/, "", nm)
            name[++m] = nm
            rest[m] = substr($0, index($0, "\t"))
        }
        END { for (k = 0; k < n; k++) for (i = 1; i <= m; i++) print name[i] "_" k rest[i] }' \
        > "$work/stream.expected"
    # BWA_MEM3_STREAM_STATS=1 reports on stderr how many reads of each batch
    # were written before the batch finished. That count depends on thread
    # timing, so it is printed as a note (it shows the run took the streaming
    # path), never checked.
    export BWA_MEM3_STREAM_STATS=1
    mem_cli_mem "$tier" "$(out stream_t4)" -t 4 -K 120000 "$ref" "$work/stream.fq"
    unset BWA_MEM3_STREAM_STATS
    echo "note: $tier: stream_t4 wrote $(awk '/stream_output:/ { n += $3; t += $5 } END { print n + 0 " of " t + 0 }' "$(out stream_t4).err") reads before their batch finished"
    export BWA_MEM3_STREAM_OUTPUT=0
    mem_cli_mem "$tier" "$(out stream_off)" -t 4 -K 120000 "$ref" "$work/stream.fq"
    unset BWA_MEM3_STREAM_OUTPUT
    mem_cli_mem "$tier" "$(out stream_t1)" -t 1 -K 120000 "$ref" "$work/stream.fq"
    for s in stream_t4 stream_off stream_t1; do
        mem_cli_records "$(out "$s")" > "$work/$s.records"
        mem_cli_same_file "$tier: $s: 3300 reads each get their original's records, in order" \
            "$work/stream.expected" "$work/$s.records"
    done

    # --- Across tiers ----------------------------------------------------------------
    if [[ -z $first_tier ]]; then
        first_tier="$tier"
    else
        for s in $sams; do
            mem_cli_same_file "$s SAM identical, $first_tier vs $tier" \
                "$work/$s.$first_tier.sam" "$(out "$s")"
        done
    fi
done

if ((MEM_CLI_FAILS > 0)); then
    echo "FAIL: mem_io_cli: $MEM_CLI_FAILS of $MEM_CLI_CHECKS checks failed"
    exit 1
fi
echo "PASS: mem_io_cli: $MEM_CLI_CHECKS checks, $MEM_CLI_SKIPS skipped"
