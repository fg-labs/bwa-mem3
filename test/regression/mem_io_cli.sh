#!/usr/bin/env bash
# test/regression/mem_io_cli.sh
#
# End-to-end checks of `bwa-mem3 mem` input and output that no other test
# makes, on the synthetic reference (test/fixtures/mem_synth.fa; the header of
# test/tools/mem_synth_ref.c lists its contigs and planted repeats), with every
# expected value worked out from the reference by hand. What the suite already
# covers is left out: the common records (mem_cli_basics.sh), the readers and
# codecs (test/unit/test_mem_io_reader.cpp, fr_fastq_diff_test.c,
# fast_reader_selftest.c, input_fail_fast.sh), the record and header writers
# (test/unit/test_mem_io_*.cpp), BAM output (bam_*.sh) and threading
# (thread_determinism.sh, stream_output_order.sh).
#
# The reads are slices of the synthetic reference, so each field follows from
# the slice and the default scoring (-A 1 -B 4 -O 6 -E 1 -L 5). A unique exact
# 100 bp hit has MAPQ 60, NM 0, MD 100, AS 100 and XS 0.
#
# Checks under each SIMD tier this host can run (mem_cli_tiers), whose SAM must
# also be byte-identical across tiers (@PG aside):
#   bases   ctg1:1401-1500 in lower case: SEQ upper case. ctg1:1601-1700 with
#           read base 50 an R: SEQ shows N there, NM 1, MD 49<ref>50, AS 98 (N
#           scores -1). The reverse complement of ctg1:2601-2700 with read base
#           30 an N: FLAG 16, SEQ the reference with N at base 71, NM 1,
#           MD 70<ref>29, AS 98. Records come out in input order.
#   split   a read made of ctg1:101-160 and ctg3:2601-2640: primary 60M40S on
#           ctg1:101, supplementary FLAG 2048 on ctg3:2601 as 60H40M with
#           SEQ/QUAL of the 40 aligned bases only, SA:Z naming the other hit
#           with S clips; MAPQ 60 on both. -Y keeps S clips and the whole
#           SEQ/QUAL on the supplementary; -M flags it 256 instead.
#   -5      a read made of ctg3:2601-2640 then ctg1:4101-4160: by default the
#           longer ctg1 hit is primary; with -5 the ctg3 hit (the read's 5'
#           end).
#   XA HN   a read inside block R (ctg1:2101-2200): the primary hit plus the XA
#           entries are exactly the three copies (ctg1 +2101, ctg2 +1101,
#           ctg4 -201), each 100M NM 0, and HN:i:2; -u appends the score (100)
#           and MAPQ (0) to each entry. Only the entries are checked, read from
#           XA or XB: the -u help text and docs/src/cli/mem.md name the tag XB,
#           docs/src/user-guide/output.md says it stays XA, and the tag name is
#           left to that documentation to settle.
#
# Checks of option and input wiring, run once (the tier does not matter):
#   reader  --legacy-reader, and gzip arriving on a pipe (/dev/stdin), give the
#           default run's SAM. A gzip file cut short after a whole member makes
#           mem exit nonzero with --legacy-reader, which reads the cut as end of
#           file, so the failure comes from main_mem closing the input.
#   -C      appends the comment verbatim after the tags.
#   -R -H   the @-lines of an -H file (a non-header line dropped) and the -R
#           line reach the header.
#   --bam   --bam=6 with -C and -R, decoded with samtools, equals the SAM text.
#   -V      against a copy of the reference whose ctg1 header carries the
#           description "first contig<TAB>desc": XR:Z:first contig desc.
#   -p      interleaved pairs with a single-end read between them: the single
#           read stays in its input position as a single-end record; without
#           it, -p on one interleaved file equals the two-file run byte for
#           byte.
#
# Inputs (env vars):
#   BWA_MEM3       path to the bwa-mem3 binary under test
#   MEM_CLI_TIERS  (optional) tiers to run instead of every runnable one
#   MEM_CLI_KEEP   (optional) 1 keeps the work directory
#   MEM_IO_CLI_REQUIRE_TOOLS (optional) 1 makes a missing samtools a failed
#                  check instead of a SKIP
# No staged inputs: the reference is committed, the reads are sliced from it.
# Optional tool: samtools (the --bam decode check), reported as SKIP when
# missing.

# Many arguments are awk programs or awk conditions (mem_cli_count), so `$1`,
# `$NF` and the like in single quotes are awk's fields, never the shell's.
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
        "$(mem_cli_tags "$rec")"
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
s_lower_ref="$(mem_cli_seq "$ref" ctg1 1401 100)"
s_lower="$(printf '%s\n' "$s_lower_ref" | tr 'ACGT' 'acgt')"
s_iupac_ref="$(mem_cli_seq "$ref" ctg1 1601 100)"
s_iupac="$(mem_cli_subst "$s_iupac_ref" 50 R)"
s_iupac_sam="$(mem_cli_subst "$s_iupac_ref" 50 N)"
iupac_ref_base="${s_iupac_ref:49:1}"
s_revn_ref="$(mem_cli_seq "$ref" ctg1 2601 100)"
s_revn="$(mem_cli_subst "$(dna_revcomp "$s_revn_ref")" 30 N)"
s_revn_sam="$(mem_cli_subst "$s_revn_ref" 71 N)"
revn_ref_base="${s_revn_ref:70:1}"

q_exact="$(qual_string 100 1)"
q_lower="$(qual_string 100 3)"
q_iupac="$(qual_string 100 4)"
q_revn="$(qual_string 100 15)"
c_exact="BC:Z:ACGT${tab}RX:Z:AAC"

{
    mem_cli_fq "exact/1 $c_exact" "$s_exact" "$q_exact"
    mem_cli_fq lower/2 "$s_lower" "$q_lower"
    mem_cli_fq iupac "$s_iupac" "$q_iupac"
    mem_cli_fq revn "$s_revn" "$q_revn"
} > "$work/fmt.fq"

# A gzip file cut short right after a whole member, in the 10-byte header of a
# second one, so what does decode is whole records and only the cut itself can
# fail the run.
gzip -c "$work/fmt.fq" > "$work/cut.fq.gz"
dd if="$work/cut.fq.gz" bs=10 count=1 2> /dev/null > "$work/cut.hdr"
cat "$work/cut.hdr" >> "$work/cut.fq.gz"

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
    mem_cli_fq chim "$s_a$s_b" "$q_split"
    mem_cli_fq chim5 "$s_b$s_a5" "$q_split"
} > "$work/split.fq"

# A read inside block R (three copies), for the XA/HN checks.
mem_cli_fq inR "$(mem_cli_seq "$ref" ctg1 2101 100)" "$(qual_string 100 12)" > "$work/xa.fq"

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

# -V: a copy of the reference with a description in the ctg1 header.
vref="$work/vref/$(mem_cli_ref_name vref)"
mkdir -p "$work/vref"
awk -v t="$tab" '/^>ctg1$/ { print ">ctg1 first contig" t "desc"; next } { print }' "$ref" > "$vref"
mem_cli_index "$vref"

# -H file: two header lines with a non-header line between them.
printf '@CO\tfrom file 1\nnot a header line\n@CO\tfrom file 2\n' > "$work/hdr.txt"

# Pairs. pp: FR, read 1 forward at ctg1:1001, read 2 the reverse complement of
# ctg1:1201-1300. qq: read 1 the reverse complement of ctg2:3201-3300, read 2
# forward at ctg2:3001. A single-end read between them.
{
    mem_cli_fq pp/1 "$s_exact" "$q_exact"
    mem_cli_fq pp/2 "$(dna_revcomp "$(mem_cli_seq "$ref" ctg1 1201 100)")" "$(qual_string 100 2)"
    mem_cli_fq single "$(mem_cli_seq "$ref" ctg1 4001 100)" "$q_lower"
    mem_cli_fq qq/1 "$(dna_revcomp "$(mem_cli_seq "$ref" ctg2 3201 100)")" "$q_iupac"
    mem_cli_fq qq/2 "$(mem_cli_seq "$ref" ctg2 3001 100)" "$(qual_string 100 5)"
} > "$work/pairs.fq"
awk 'NR < 9 || NR > 12' "$work/pairs.fq" > "$work/pairs_nosingle.fq"
awk 'int((NR - 1) / 4) % 2 == 0' "$work/pairs_nosingle.fq" > "$work/pairs_r1.fq"
awk 'int((NR - 1) / 4) % 2 == 1' "$work/pairs_nosingle.fq" > "$work/pairs_r2.fq"

# out <name>: the SAM path of run <name> under the current tier.
out() {
    printf '%s/%s.%s.sam' "$work" "$1" "$tier"
}

# ---------------------------------------------------------------------------
# Alignment checks, under each tier
# ---------------------------------------------------------------------------

# Every SAM named here is compared across tiers at the end.
sams="default split y_opt m_opt p5 xa u_opt"

first_tier=""
for tier in $(mem_cli_tiers); do

    # --- Read bases in SEQ, NM, MD and AS --------------------------------------
    sam="$(out default)"
    mem_cli_mem "$tier" "$sam" "$ref" "$work/fmt.fq"
    t100="NM:i:0${tab}MD:Z:100${tab}AS:i:100${tab}XS:i:0${tab}HN:i:0"
    check_fields "$tier" "$sam" lower "0 ctg1 1401 60 100M" "$s_lower_ref" "$q_lower" "$t100"
    check_fields "$tier" "$sam" iupac "0 ctg1 1601 60 100M" "$s_iupac_sam" "$q_iupac" \
        "NM:i:1${tab}MD:Z:49${iupac_ref_base}50${tab}AS:i:98${tab}XS:i:0${tab}HN:i:0"
    check_fields "$tier" "$sam" revn "16 ctg1 2601 60 100M" "$s_revn_sam" "$(reverse_string "$q_revn")" \
        "NM:i:1${tab}MD:Z:70${revn_ref_base}29${tab}AS:i:98${tab}XS:i:0${tab}HN:i:0"
    mem_cli_expect "$tier: one record per read, in input order" "exact lower iupac revn " \
        "$(mem_cli_records "$sam" | awk -F '\t' '{ printf "%s ", $1 }')"

    # --- Split reads: SA, hard and soft clips, -Y, -M, -5 ----------------------
    sp="$(out split)"
    mem_cli_mem "$tier" "$sp" "$ref" "$work/split.fq"
    mem_cli_expect "$tier: chim primary" \
        "chim 0 ctg1 101 60 60M40S $s_a$s_b $q_split SA:Z:ctg3,2601,+,60S40M,60,0;" \
        "$(record_of "$sp" chim 0 | awk -F '\t' '{ print $1, $2, $3, $4, $5, $6, $10, $11, $16 }')"
    mem_cli_expect "$tier: chim supplementary (hard-clipped SEQ/QUAL)" \
        "chim 2048 ctg3 2601 60 60H40M $s_b ${q_split:60} NM:i:0 MD:Z:40 AS:i:40 SA:Z:ctg1,101,+,60M40S,60,0;" \
        "$(record_of "$sp" chim 2048 | awk -F '\t' '{ print $1, $2, $3, $4, $5, $6, $10, $11, $12, $13, $14, $16 }')"
    mem_cli_expect "$tier: chim5 primary is the longer hit" \
        "chim5 0 ctg1 4101 60 40S60M|chim5 2048 ctg3 2601 60 40M60H $s_b ${q_split:0:40}" \
        "$(record_of "$sp" chim5 0 | awk -F '\t' '{ printf "%s %s %s %s %s %s|", $1, $2, $3, $4, $5, $6 }')$(record_of "$sp" chim5 2048 | awk -F '\t' '{ print $1, $2, $3, $4, $5, $6, $10, $11 }')"

    mem_cli_mem "$tier" "$(out y_opt)" -Y "$ref" "$work/split.fq"
    mem_cli_expect "$tier: -Y supplementary keeps S clips and the whole read" \
        "chim 2048 ctg3 2601 60S40M $s_a$s_b $q_split" \
        "$(record_of "$(out y_opt)" chim 2048 | awk -F '\t' '{ print $1, $2, $3, $4, $6, $10, $11 }')"

    mem_cli_mem "$tier" "$(out m_opt)" -M "$ref" "$work/split.fq"
    mem_cli_expect "$tier: -M flags the shorter hits 256, not 2048" "chim 256 ctg3 2601|chim5 256 ctg3 2601|" \
        "$(mem_cli_records "$(out m_opt)" | awk -F '\t' '$2 >= 256 { printf "%s %s %s %s|", $1, $2, $3, $4 }')"

    mem_cli_mem "$tier" "$(out p5)" -5 "$ref" "$work/split.fq"
    mem_cli_expect "$tier: -5 makes the 5' hit primary" \
        "chim5 0 ctg3 2601 60 40M60S|chim5 2048 ctg1 4101 60 40H60M $s_a5 ${q_split:40}" \
        "$(record_of "$(out p5)" chim5 0 | awk -F '\t' '{ printf "%s %s %s %s %s %s|", $1, $2, $3, $4, $5, $6 }')$(record_of "$(out p5)" chim5 2048 | awk -F '\t' '{ print $1, $2, $3, $4, $5, $6, $10, $11 }')"

    # --- XA, HN, -u ------------------------------------------------------------
    mem_cli_mem "$tier" "$(out xa)" "$ref" "$work/xa.fq"
    # R is at ctg1:2001 and ctg2:1001 (+) and reverse-complemented at ctg4:101;
    # the read's R offsets 100-199 sit at ctg4:201-300 on the reverse strand.
    mem_cli_expect "$tier: inR primary + XA hits" \
        "ctg1,+2101,100M,0 ctg2,+1101,100M,0 ctg4,-201,100M,0 " "$(hits "$(out xa)" inR)"
    mem_cli_expect "$tier: HN counts the inR read's two other hits" 2 \
        "$(mem_cli_tag "$(mem_cli_primary "$(out xa)" inR)" HN)"
    mem_cli_mem "$tier" "$(out u_opt)" -u "$ref" "$work/xa.fq"
    u_rec="$(mem_cli_primary "$(out u_opt)" inR)"
    # The entries under whichever of XA or XB the record carries (see the header).
    mem_cli_expect "$tier: -u entries are the XA entries plus score 100 and MAPQ 0" \
        "$(mem_cli_tag "$(mem_cli_primary "$(out xa)" inR)" XA | sed 's/;/,100,0;/g')" \
        "$(mem_cli_tag "$u_rec" XA)$(mem_cli_tag "$u_rec" XB)"

    if [[ -z $first_tier ]]; then
        first_tier="$tier"
    else
        for s in $sams; do
            mem_cli_same_file "$s SAM identical, $first_tier vs $tier" \
                "$work/$s.$first_tier.sam" "$(out "$s")"
        done
    fi
done

if [[ $(mem_cli_tiers | wc -l) -lt 2 ]]; then
    mem_cli_skip "SAM identical across tiers" "only one tier ($first_tier) runs on this host"
fi

# ---------------------------------------------------------------------------
# Option and input wiring, once (under the first tier)
# ---------------------------------------------------------------------------

tier="$first_tier"
sam="$(out default)"

# --- Readers -------------------------------------------------------------------
mem_cli_mem "$tier" "$(out legacy)" --legacy-reader "$ref" "$work/fmt.fq"
mem_cli_same_file "--legacy-reader gives the default reader's SAM" "$sam" "$(out legacy)"
gzip -c "$work/fmt.fq" | mem_cli_mem "$tier" "$(out stdin)" "$ref" /dev/stdin
mem_cli_same_file "gzip from a pipe (/dev/stdin)" "$sam" "$(out stdin)"
mem_cli_expect_fail "a gzip file cut short makes mem exit nonzero (--legacy-reader)" '\[gzclose\]' \
    "$tier" "$work/bad.sam" --legacy-reader "$ref" "$work/cut.fq.gz"

# --- -C, -R, -H ----------------------------------------------------------------
mem_cli_mem "$tier" "$(out c_opt)" -C "$ref" "$work/fmt.fq"
mem_cli_expect "-C appends the comment after the tags" "${t100}${tab}${c_exact}" \
    "$(mem_cli_tags "$(mem_cli_primary "$(out c_opt)" exact)")"
mem_cli_mem "$tier" "$(out rh)" -R '@RG\tID:grp1\tSM:s1' -H "$work/hdr.txt" "$ref" "$work/fmt.fq"
mem_cli_expect "the -H file's @-lines and the -R line reach the header" \
    "@CO from file 1|@CO from file 2|@RG ID:grp1 SM:s1|" \
    "$(awk -F '\t' '$1 == "@CO" || $1 == "@RG" { $1 = $1; printf "%s|", $0 }' "$(out rh)")"

# --- --bam -----------------------------------------------------------------------
if command -v samtools > /dev/null 2>&1; then
    mem_cli_mem "$tier" "$(out rhc)" -C -R '@RG\tID:grp1\tSM:s1' "$ref" "$work/fmt.fq"
    mem_cli_mem "$tier" "$work/out.bam" --bam=6 -C -R '@RG\tID:grp1\tSM:s1' "$ref" "$work/fmt.fq"
    if ! samtools view -h "$work/out.bam" > "$work/bam.sam" 2> "$work/bam.err"; then
        : > "$work/bam.sam" # mem_cli_same_file fails on an empty file
    fi
    mem_cli_same_file "--bam=6 decodes to the SAM text" "$(out rhc)" "$work/bam.sam"
elif [[ ${MEM_IO_CLI_REQUIRE_TOOLS:-0} == 1 ]]; then
    mem_cli_expect "--bam=6 decodes to the SAM text" "samtools installed" \
        "samtools not found (MEM_IO_CLI_REQUIRE_TOOLS=1)"
else
    mem_cli_skip "--bam=6 decodes to the SAM text" "samtools not installed"
fi

# --- -V --------------------------------------------------------------------------
mem_cli_mem "$tier" "$(out v_opt)" -V "$vref" "$work/fmt.fq"
mem_cli_expect "-V XR from the ctg1 description, TAB as space" "first contig desc" \
    "$(mem_cli_tag "$(mem_cli_primary "$(out v_opt)" exact)" XR)"

# --- -p ----------------------------------------------------------------------------
pr="$(out pairs)"
mem_cli_mem "$tier" "$pr" -p -I 300,30 "$ref" "$work/pairs.fq"
mem_cli_expect "-p keeps the single read in place, single-end" "pp pp single qq qq |0 * 0 0" \
    "$(mem_cli_records "$pr" | awk -F '\t' '{ printf "%s ", $1 }')|$(record_of "$pr" single 0 | awk -F '\t' '{ print $2, $7, $8, $9 }')"
mem_cli_mem "$tier" "$(out pairs_nosingle)" -p -I 300,30 "$ref" "$work/pairs_nosingle.fq"
mem_cli_mem "$tier" "$(out pairs_2file)" -I 300,30 "$ref" "$work/pairs_r1.fq" "$work/pairs_r2.fq"
mem_cli_same_file "-p interleaved equals two files" "$(out pairs_2file)" "$(out pairs_nosingle)"

if mem_cli_verdict; then
    echo "PASS: $MEM_CLI_VERDICT"
else
    echo "FAIL: $MEM_CLI_VERDICT"
    exit 1
fi
