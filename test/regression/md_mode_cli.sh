#!/usr/bin/env bash
# test/regression/md_mode_cli.sh
#
# Regression: the --md command line (full|short|none), end to end through the
# SAM writer and both BAM writers (the plain and the --meth one).
#
# Asserted here:
#  1. An unknown value is rejected with its exact message and a non-zero exit;
#     --compat refuses --md none, accepts --md full and --md short unchanged,
#     and --compat-allow-divergent overrides the refusal (with its warning).
#  2. --md full is record-for-record identical to the default, with and
#     without --meth, single- and paired-end, in SAM and in BAM.
#  3. --md none differs from full only by the MD:Z field: every record loses
#     it and nothing else changes, in the same eight combinations.
#  4. --md short without --meth is identical to full.
#  5. --md short under --meth: a fully converted read with one real variant
#     lists only that variant (MD:Z:<n>A<m>), where full is the literal MD of
#     every converted C plus the variant; the rest of the record matches
#     full; BAM equals SAM.
#  6. Non-primary records: a mapped secondary (under -a) and, without --meth, a
#     supplementary record are present and carry NM and MD by default; full
#     equals the default, none drops only MD, short changes only MD (to the
#     bare match length under --meth, where the reads carry no real variant;
#     not at all without it); BAM equals SAM.
#
# Inputs:
#   BWA_MEM3 — the bwa-mem3 binary under test: a bare command name resolved via
#              PATH, or a path to the binary (a relative one is made absolute
#              before the `cd "$WORK"` below)
set -euo pipefail
if [[ -z ${BWA_MEM3:-} ]]; then
    echo "FAIL: BWA_MEM3 must be set" >&2
    exit 1
fi
command -v samtools > /dev/null 2>&1 || {
    echo "SKIP: samtools not on PATH (needed to read BAM back as SAM)"
    exit 0
}
case "$BWA_MEM3" in
    */*) BWA_MEM3="$(cd "$(dirname "$BWA_MEM3")" && pwd)/$(basename "$BWA_MEM3")" ;;
esac

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
cd "$WORK"
fail() {
    echo "FAIL: $*" >&2
    exit 1
}

# Deterministic 1500 bp reference (same as meth_context_cli).
REF=TCATTGGCTATCCTAACCCGACCCTAGGAGCGGTTGGCGTGTATGCCGTGAATTTTCTCATTTCCGCTAGACATAATCGTTCTGCCTATATCTGGACAACATCCCGGCGACTTAGGCGACCCACAGAATCGTCCCTTCTAACGTAGTTCGCATAGTTCCCGTCCGTAGCCGGACTATTCGAACACCCAGTATTCGATTAACTCGGGCTTGACGTATTAGAGGCGTTAGTGTGCCAGGTAAGATACGCCAACGGAATTAACCTCTGTGACACTCCGCGGAGCCTTCGGACATATAAGTGATCGGGTCTACGTTTGTTAGACTTGAGACGTCTGTTAAGAGTTGGGTCTAATAAATCGCCTACACGTGGAGTCTAACGGGGAAGCGTCGAATCCTGATACATCATATAATGGAGCGTGTTATGAAAAAAGAGCATTCCATTGTACGAGCCGTGCCAGAAACGGCTTGACTACGTGAGCGTAGTGTTAGATAAACAGGAAACACTGACGCGGTTAGAAGGCGGATTGCCGGTAGGTTTTGGAAACATAAATACACACGGTATCATGTTGGGTCACGATTCCTATCACCGCACAGGGCCAACCATAGAAGAACTGAAAGAACTAATCTGGCGGCGGGCTCGGTGCTTATATTTTCCACCCAACATCGTGCACATTAGGCTCACCGCGCCCTACGGGCGAAGGGTGCGTACGGTGTTTATAAGGCGTGACGGCCCCAAGTAGAGGGTAATTCTGTGAAAGAATCTCAGGACGGTGGCATGAATTCAATTCCTTTTAAACCTATCGTTCCGACCTTATGCAATCCTTCAATGAAGATCGTCAACGACCATCGTTCTTCTGCTTTAAGTGTGAGTTCTCTCTTACAAGCTAATACACCCCAGCGTTCTCCGTACTCTTCACTGCCCAAGCGAGGCTAACCTTTTGAAATGTCACAGTCGAAGCATATCTCCCGTACATCTTTTTCGGAGATCGCAGCTCGCGGAGCTATAAGCGACTTAAGCCCTTGTGTCGGTGATCCCAAGGGTCTGACTCCTGTACCAGGGTTACTGTTTCGCTTTACGGAGTAGCCTGTGAGGTGAACTGAAAGGAGCATATTTGAGATCTAAGATAGGGTCCTCCTCTGCGTCTACGTTCTCTCCGTTACGTACGGCTTCGCACCGGAGTGCATCTTGGCCCCGAAACGCACTGTGTGTGCTGATACAGCGTCCCTGGCCGGCCATGGGTTCAGAACTCCCGGGAACGCTTTTCAACTTAGAGGAACCCCGTCATGGAAGTAGATCGCGTCGAATGAGGGAGTTAGTCCTCGTTCCAGCTGGTAATTGTTTTACCGCTTGGGACCACTATAGGCCGCGGGTAGAAGTTGCTGGGTGTTGATTCCAACCCTCGAACCACGATACGACCTGCCATTTATGGCACAGTAAGGTTCAAACAGCATAATGAATACAGTTATAGTAACTTCCTCACGTACGATTAGGACGCAGCCTTG
printf '>chrA\n%s\n' "$REF" > ref.fa
"$BWA_MEM3" index ref.fa > /dev/null 2>&1 || fail "index nonzero exit"
"$BWA_MEM3" index --meth ref.fa > /dev/null 2>&1 || fail "index --meth nonzero exit"

qual=$(printf 'I%.0s' $(seq 1 100))
# r: chrA:201-300 with every C converted to T (OT) and its first A changed to
# G (a real variant). p: the same locus unconverted with that variant, and a
# second read with a 2 bp deletion, for the non-meth runs.
region=${REF:200:100}
head_a=${region%%A*}
va=${#head_a} # offset of the first A
variant="${region:0:$va}G${region:$((va + 1))}"
converted=$(printf '%s' "$variant" | tr C T)
printf '@r\n%s\n+\n%s\n' "$converted" "$qual" > m.fq
del=${REF:700:50}${REF:752:50}
printf '@p\n%s\n+\n%s\n@d\n%s\n+\n%s\n' "$variant" "$qual" "$del" "$qual" > p.fq

# Every compared side is a file checked for at least one record, so a run or a
# samtools view that produced nothing fails here rather than diffing equal.
# Each is assigned in its own statement: set -e stops at a failing $(...) only
# when it is the last substitution of the statement.
recs() { # $1 = SAM or BAM; writes its records to $1.recs and echoes that path
    local out="$1.recs"
    samtools view "$1" > "$out" || fail "samtools view $1 failed"
    [ -s "$out" ] || fail "$1: no records"
    echo "$out"
}
nomd() { # $1 = SAM or BAM; its records with MD:Z removed, in $1.nomd
    local out="$1.nomd"
    samtools view -x MD "$1" > "$out" || fail "samtools view -x MD $1 failed"
    [ -s "$out" ] || fail "$1: no records"
    echo "$out"
}
run() { # $1 = output, rest = mem argv (before ref/reads)
    local out="$1"
    shift
    "$BWA_MEM3" mem -t 1 "$@" > "$out" 2> "$out.err" || fail "mem $* nonzero exit: $(head -2 "$out.err" | tr '\n' ' ')"
}
same() { # $1, $2 = record files, $3 = message
    cmp -s "$1" "$2" || fail "$3: $(diff "$1" "$2" | head -4 | cut -c1-200 | tr '\n' ' ')"
}
# The literal MD of an ungapped read against `$1` that differs at every column
# whose reference base is in the set $2, plus column $3.
literal_md() {
    local seq="$1" set="$2" extra="$3" i b n=0 md=""
    for ((i = 0; i < ${#seq}; i++)); do
        b=${seq:i:1}
        if [[ $set == *"$b"* || $i -eq $extra ]]; then
            md+="$n$b"
            n=0
        else
            n=$((n + 1))
        fi
    done
    echo "MD:Z:$md$n"
}

# --- 1. rejections and --compat ----------------------------------------------
if "$BWA_MEM3" mem --md bogus ref.fa p.fq > /dev/null 2> err.txt; then fail "--md bogus exited 0"; fi
grep -Fq -- "unknown --md 'bogus' (full|short|none)" err.txt || fail "--md bogus: message missing: $(head -2 err.txt | tr '\n' ' ')"
if "$BWA_MEM3" mem --compat=bwa-mem --md none ref.fa p.fq > /dev/null 2> err.txt; then fail "--compat --md none exited 0"; fi
# The refusal itself, not the usage text (which also names --md none).
grep -Fq -- "mutually exclusive: --md none" err.txt || fail "--compat --md none: refusal does not name it: $(head -2 err.txt | tr '\n' ' ')"
run compat.def.sam --compat=bwa-mem ref.fa p.fq
compat_def=$(recs compat.def.sam)
for m in full short; do
    run "compat.$m.sam" --compat=bwa-mem --md "$m" ref.fa p.fq
    compat_m=$(recs "compat.$m.sam")
    same "$compat_def" "$compat_m" "--compat --md $m differs from --compat"
done
run compat.none.sam --compat=bwa-mem --md none --compat-allow-divergent ref.fa p.fq
grep -Fq -- "--compat-allow-divergent: proceeding with --md none" compat.none.sam.err \
    || fail "--compat-allow-divergent --md none: no warning: $(head -2 compat.none.sam.err | tr '\n' ' ')"
compat_def_nomd=$(nomd compat.def.sam)
compat_none=$(recs compat.none.sam)
same "$compat_def_nomd" "$compat_none" "--compat-allow-divergent --md none differs from --compat beyond MD"

# --- 2-3. full == default; none drops only MD (SE and PE, SAM and BAM) -------
for mode in plain meth; do
    if [ "$mode" = meth ]; then
        fq=m.fq
        extra=(--meth)
    else
        fq=p.fq
        extra=()
    fi
    for ends in se pe; do
        reads=("$fq")
        [ "$ends" = se ] || reads=("$fq" "$fq")
        for out in def full none; do
            md=()
            [ "$out" = def ] || md=(--md "$out")
            for fmt in sam bam; do
                bam=()
                [ "$fmt" = sam ] || bam=(--bam)
                run "$mode.$ends.$out.$fmt" "${extra[@]+"${extra[@]}"}" "${bam[@]+"${bam[@]}"}" \
                    "${md[@]+"${md[@]}"}" ref.fa "${reads[@]}"
            done
        done
        for fmt in sam bam; do
            tag="$mode $ends $fmt"
            def=$(recs "$mode.$ends.def.$fmt")
            full=$(recs "$mode.$ends.full.$fmt")
            none=$(recs "$mode.$ends.none.$fmt")
            def_nomd=$(nomd "$mode.$ends.def.$fmt")
            same "$def" "$full" "$tag: --md full differs from the default"
            grep -q 'MD:Z:' "$def" || fail "$tag: default records carry no MD"
            if grep -q 'MD:Z:' "$none"; then fail "$tag: --md none still has MD"; fi
            same "$def_nomd" "$none" "$tag: --md none differs from full beyond MD"
        done
    done
done

# --- 4. short without --meth == full -----------------------------------------
run plain.short.sam --md short ref.fa p.fq
plain_def=$(recs plain.se.def.sam)
plain_short=$(recs plain.short.sam)
same "$plain_def" "$plain_short" "--md short without --meth differs from full"

# --- 5. short under --meth ---------------------------------------------------
run meth.short.sam --meth --md short ref.fa m.fq
run meth.short.bam --meth --bam --md short ref.fa m.fq
short=$(recs meth.short.sam)
short_bam=$(recs meth.short.bam)
def=$(recs meth.se.def.sam)
want="MD:Z:${va}A$((100 - va - 1))"
got=$(grep -o 'MD:Z:[^[:space:]]*' "$short" || true)
[ "$got" = "$want" ] || fail "--meth --md short: want $want, got $got"
# full lists every converted reference C and the variant A, literally.
want_full=$(literal_md "$region" C "$va")
got_full=$(grep -o 'MD:Z:[^[:space:]]*' "$def" || true)
[ "$got_full" = "$want_full" ] || fail "--meth --md full: want $want_full, got $got_full"
def_nomd=$(nomd meth.se.def.sam)
short_nomd=$(nomd meth.short.sam)
same "$def_nomd" "$short_nomd" "--meth --md short differs from full beyond MD"
same "$short" "$short_bam" "--meth --md short: BAM differs from SAM"

# --- 6. non-primary records (secondary and supplementary) --------------------
# A second contig copies chrA:1101-1400, so a read from inside it has a
# secondary hit under -a; a read joining chrA:501-560 and chrA:1001-1060 is
# split. Without --meth the split is a supplementary (0x800) record; --meth
# implies -M, which writes it as a secondary (0x100) one instead.
printf '>chrA\n%s\n>chrB\n%s\n' "$REF" "${REF:1100:300}" > np.fa
"$BWA_MEM3" index np.fa > /dev/null 2>&1 || fail "index np.fa nonzero exit"
"$BWA_MEM3" index --meth np.fa > /dev/null 2>&1 || fail "index --meth np.fa nonzero exit"
dup=${REF:1150:100}
chim=${REF:500:60}${REF:1000:60}
q120=$(printf 'I%.0s' $(seq 1 120))
printf '@s\n%s\n+\n%s\n@c\n%s\n+\n%s\n' "$dup" "$qual" "$chim" "$q120" > np.fq
printf '@s\n%s\n+\n%s\n@c\n%s\n+\n%s\n' "$(printf '%s' "$dup" | tr C T)" "$qual" \
    "$(printf '%s' "$chim" | tr C T)" "$q120" > npm.fq
# Succeeds when $1 (a record file) holds a mapped record whose flag has bit $2.
has_flag() {
    awk -F'\t' -v b="$2" 'int($2 / 4) % 2 == 0 && int($2 / b) % 2 == 1 { n++ } END { exit !n }' "$1"
}
for mode in plain meth; do
    if [ "$mode" = meth ]; then
        fq=npm.fq
        extra=(--meth)
    else
        fq=np.fq
        extra=()
    fi
    for out in def full short none; do
        md=()
        [ "$out" = def ] || md=(--md "$out")
        for fmt in sam bam; do
            bam=()
            [ "$fmt" = sam ] || bam=(--bam)
            run "np.$mode.$out.$fmt" -a "${extra[@]+"${extra[@]}"}" "${bam[@]+"${bam[@]}"}" \
                "${md[@]+"${md[@]}"}" np.fa "$fq"
        done
    done
    for fmt in sam bam; do
        tag="non-primary $mode $fmt"
        def=$(recs "np.$mode.def.$fmt")
        full=$(recs "np.$mode.full.$fmt")
        short=$(recs "np.$mode.short.$fmt")
        none=$(recs "np.$mode.none.$fmt")
        has_flag "$def" 256 || fail "$tag: no mapped secondary record"
        if [ "$mode" = plain ]; then
            has_flag "$def" 2048 || fail "$tag: no mapped supplementary record"
        else
            # Under -M the split segment is a secondary record that keeps SA:Z.
            awk -F'\t' 'int($2 / 4) % 2 == 0 && int($2 / 256) % 2 == 1 && /\tSA:Z:/ { n++ } END { exit !n }' "$def" \
                || fail "$tag: no mapped secondary split record (SA:Z)"
        fi
        # Every mapped record, primary or not, carries NM and MD by default.
        if awk -F'\t' '$3 != "*" && !(/\tNM:i:/ && /\tMD:Z:/)' "$def" | grep -q .; then
            fail "$tag: a mapped record lacks NM or MD"
        fi
        same "$def" "$full" "$tag: --md full differs from the default"
        if grep -q 'MD:Z:' "$none"; then fail "$tag: --md none still has MD"; fi
        def_nomd=$(nomd "np.$mode.def.$fmt")
        short_nomd=$(nomd "np.$mode.short.$fmt")
        same "$def_nomd" "$none" "$tag: --md none differs from full beyond MD"
        same "$def_nomd" "$short_nomd" "$tag: --md short differs from full beyond MD"
        if [ "$mode" = plain ]; then
            same "$def" "$short" "$tag: --md short without --meth differs from full"
        else
            # The reads carry no real variant, so short leaves out every
            # conversion: each MD is a bare match length.
            if grep -o 'MD:Z:[^[:space:]]*' "$short" | grep -qv '^MD:Z:[0-9]*$'; then
                fail "$tag: --md short lists a mismatch: $(grep -o 'MD:Z:[^[:space:]]*' "$short" | tr '\n' ' ')"
            fi
            if ! grep -o 'MD:Z:[^[:space:]]*' "$def" | grep -qv '^MD:Z:[0-9]*$'; then
                fail "$tag: --md full lists no conversion"
            fi
        fi
    done
    for out in def short none; do
        sam_recs=$(recs "np.$mode.$out.sam")
        bam_recs=$(recs "np.$mode.$out.bam")
        same "$sam_recs" "$bam_recs" "non-primary $mode --md $out: BAM differs from SAM"
    done
done

echo "PASS: md_mode_cli (unknown value and --compat --md none rejected, --compat accepts full/short, --compat-allow-divergent overrides none; full == default and literal; none drops only MD, SE and PE; short == full without --meth; short lists only the real variant under --meth; secondary and supplementary records follow the same rules; SAM and BAM agree)"
