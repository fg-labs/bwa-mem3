#!/usr/bin/env bash
# test/regression/aln_known_snvs.sh
#
# `bwa-mem3 aln --known-snvs` from the command line, on a seeded synthetic reference: a VCF with no
# SNV leaves single-end, paired-end and --meth output unchanged; with sites, records without XV:Z
# are unchanged, short reads carrying an ALT plus errors (over the mismatch budget, so stock leaves
# them unmapped) are placed at their true locus with an A vote in XV:Z, single- and paired-end;
# the output does not depend on the thread count; --bam carries XV:Z; bad VCFs are refused before
# any output.
#
# Env:
#   BWA_MEM3  the binary (default: the repository's ./bwa-mem3)
set -euo pipefail
trap 'echo "FAIL: line $LINENO: $BASH_COMMAND" >&2' ERR
here=$(cd "$(dirname "$0")" && pwd)
bin="${BWA_MEM3:-$here/../../bwa-mem3}"
[[ $bin = /* ]] || bin=$PWD/$bin
[[ -x $bin ]] || {
    echo "FAIL: BWA_MEM3 ('$bin') is not an executable binary" >&2
    exit 1
}
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
export TMPDIR=$tmp

"$bin" aln --help > "$tmp/help"
grep -q -- '--known-snvs FILE' "$tmp/help"
echo 'PASS: help'

# Reference: two random 6,000 bp contigs. Sites: every 150 bp of c1 from 600 to 5,400 (1-based),
# ALT another base. Reads (true position in the name):
#   ref*/alt*  25 bp over a site, REF or ALT there, with one error elsewhere (both place in stock)
#   over*      25 bp, ALT plus two errors: three differences, over the budget of 2 (stock: unmapped)
#   bg*        25 bp from c2, no site
# Pairs: R1 is one of the reads above (fragment start), R2 the reverse complement of the
# fragment's last 25 bases (fragments 60-90 bp), never over a site.
cat > "$tmp/gen.awk" << 'AWK'
function rseq(n,   s, i) { s = ""; for (i = 0; i < n; i++) s = s substr("ACGT", ri(4) + 1, 1); return s }
function mut(s, k) { return substr(s, 1, k - 1) other(substr(s, k, 1)) substr(s, k + 1) }
function put(name, s, f) { printf "@%s\n%s\n+\n%s\n", name, s, substr(Q, 1, length(s)) > f }
function pair(name, c, beg, r1,   flen, r2) {
    flen = 60 + ri(31)
    r2 = rc(substr(c, beg + flen - 25, 25))
    put(name, r1, dir "/r1.fq"); put(name, r2, dir "/r2.fq"); put(name, r1, dir "/se.fq")
}
BEGIN {
    seed_rng(20261109)
    Q = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
    c1 = rseq(6000); c2 = rseq(6000)
    printf ">c1\n%s\n>c2\n%s\n", c1, c2 > (dir "/ref.fa")
    print "##fileformat=VCFv4.2\n##contig=<ID=c1>\n##contig=<ID=c2>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" > (dir "/sites.vcf")
    n = 0
    for (p = 600; p <= 5400; p += 150) {
        r = substr(c1, p, 1); a = other(r)
        printf "c1\t%d\t.\t%s\t%s\t.\tPASS\t.\n", p, r, a > (dir "/sites.vcf")
        off = 5 + ri(15); beg = p - off      # the site is base off + 1 of the read
        base = substr(c1, beg, 25)
        e1 = (off + 6) % 25 + 1; e2 = (off + 14) % 25 + 1
        alt = substr(base, 1, off) a substr(base, off + 2)
        pair("ref" n "|c1|" beg, c1, beg, mut(base, e1))
        pair("alt" n "|c1|" beg, c1, beg, mut(alt, e1))
        pair("over" n "|c1|" beg, c1, beg, mut(mut(alt, e1), e2))
        g = 200 + ri(5500)
        pair("bg" n "|c2|" g, c2, g, substr(c2, g, 25))
        n++
    }
    # an empty set: one indel, no SNV
    print "##fileformat=VCFv4.2\n##contig=<ID=c1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" > (dir "/none.vcf")
    printf "c1\t600\t.\t%sA\t%s\t.\tPASS\t.\n", substr(c1, 600, 1), substr(c1, 600, 1) > (dir "/none.vcf")
    printf "c1\t%d\t.\t%s\t%s\t.\tPASS\t.\n", 700, other(substr(c1, 700, 1)), "A" > (dir "/badref.vcf")
    printf "chr9\t5\t.\tA\tC\t.\tPASS\t.\n" > (dir "/badctg.vcf")
}
AWK
awk -v dir="$tmp" -f "$here/aln_rng.awk" -f "$tmp/gen.awk"
for f in badref badctg; do
    {
        printf '##fileformat=VCFv4.2\n##contig=<ID=c1>\n##contig=<ID=chr9>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        cat "$tmp/$f.vcf"
    } > "$tmp/$f.tmp"
    mv "$tmp/$f.tmp" "$tmp/$f.vcf"
done
"$bin" index "$tmp/ref.fa" > "$tmp/index.log" 2>&1
"$bin" index --meth "$tmp/ref.fa" >> "$tmp/index.log" 2>&1

body() { grep -v '^@PG' "$1"; }
# bit(flag, b): whether SAM flag has bit b (a power of two); POSIX awk has no and()
bit='function bit(f, b) { return int(f / b) % 2 } '

recs() { grep -v '^@' "$1"; }
run() { # run NAME ARGS...: aln with ARGS (options, index, reads) into $tmp/NAME.sam
    local name=$1
    shift
    "$bin" aln -v 1 -f "$tmp/$name.sam" "$@" > /dev/null 2> "$tmp/$name.err" || {
        echo "FAIL: $name: aln exited non-zero" >&2
        cat "$tmp/$name.err" >&2
        exit 1
    }
}

# 1. No SNV in the VCF: byte-identical output.
for mode in plain meth; do
    m=()
    [ $mode = meth ] && m=(--meth)
    run "stock.se.$mode" "${m[@]}" "$tmp/ref.fa" "$tmp/se.fq"
    run "none.se.$mode" "${m[@]}" --known-snvs "$tmp/none.vcf" "$tmp/ref.fa" "$tmp/se.fq"
    run "stock.pe.$mode" "${m[@]}" "$tmp/ref.fa" "$tmp/r1.fq" "$tmp/r2.fq"
    run "none.pe.$mode" "${m[@]}" --known-snvs "$tmp/none.vcf" "$tmp/ref.fa" "$tmp/r1.fq" "$tmp/r2.fq"
    for e in se pe; do
        cmp <(body "$tmp/stock.$e.$mode.sam") <(body "$tmp/none.$e.$mode.sam") || {
            echo "FAIL: no-SNV VCF changed the $e $mode output" >&2
            exit 1
        }
    done
done
echo 'PASS: no-snv-identical (single, paired; plain, --meth)'

# 2. With sites (plain): records without XV are unchanged; over-budget ALT reads are placed.
run ks.se --known-snvs "$tmp/sites.vcf" "$tmp/ref.fa" "$tmp/se.fq"
run ks.pe --known-snvs "$tmp/sites.vcf" "$tmp/ref.fa" "$tmp/r1.fq" "$tmp/r2.fq"
for e in se pe; do
    # a fragment (a read, or a pair) with no XV on any of its records is unchanged
    paste -d $'\001' <(recs "$tmp/stock.$e.plain.sam") <(recs "$tmp/ks.$e.sam") | awk -F $'\001' -v pe="$([ "$e" = pe ] && echo 1 || echo 0)" '
        function done_frag() { if (!xv && !same) { print "differs: " name > "/dev/stderr"; bad = 1 } }
        { split($1, f, "\t") }
        !pe || NR % 2 == 1 { if (NR > 1) done_frag(); name = f[1]; xv = 0; same = 1 }
        { if ($2 ~ /\tXV:Z:/) xv = 1; if ($1 != $2) same = 0 }
        END { if (NR) done_frag(); exit bad }' || {
        echo "FAIL: $e: a read or pair without XV differs from stock" >&2
        exit 1
    }
    # stock leaves every over-budget read unmapped; crediting places each at its locus, with an A vote
    awk -F'\t' "$bit"'!/^@/ && $1 ~ /^over/ && (!bit($2, 1) || bit($2, 64)) { if (!bit($2, 4)) bad = 1; n++ }
        END { exit !(n > 0 && !bad) }' "$tmp/stock.$e.plain.sam" || {
        echo "FAIL: $e: stock placed an over-budget ALT read (the test reads are wrong)" >&2
        exit 1
    }
    awk -F'\t' "$bit"'!/^@/ && $1 ~ /^over/ && (!bit($2, 1) || bit($2, 64)) {
            n++; split($1, t, "|")
            if (!bit($2, 4) && $3 == t[2] && $4 == t[3] && $5 >= 20 && $0 ~ /\tXV:Z:c1:[0-9]+:[ACGT]>[ACGT],A,2,2,/) ok++
        }
        END { if (ok < n * 0.9) { printf "%d of %d placed\n", ok, n > "/dev/stderr"; exit 1 } }' "$tmp/ks.$e.sam" || {
        echo "FAIL: $e: over-budget ALT reads not placed with an A vote" >&2
        exit 1
    }
    # REF and ALT reads both carry their vote
    for a in ref:R alt:A; do
        awk -F'\t' -v p="${a%%:*}" -v v="${a##*:}" "$bit"'!/^@/ && index($1, p) == 1 && (!bit($2, 1) || bit($2, 64)) {
                n++; if ($0 ~ ("\tXV:Z:c1:[0-9]+:[ACGT]>[ACGT]," v ",1,2,1,")) ok++ }
            END { exit !(n > 0 && ok >= n * 0.9) }' "$tmp/ks.$e.sam" || {
            echo "FAIL: $e: ${a%%:*} reads lack an XV ${a##*:} vote with a same-locus twin" >&2
            exit 1
        }
    done
done
echo 'PASS: sites-plain (untouched records unchanged; over-budget ALT reads placed; votes and twins)'

# 3. --meth with sites: runs, single and paired, and tags reads over the sites.
run ks.se.meth --meth --known-snvs "$tmp/sites.vcf" "$tmp/ref.fa" "$tmp/se.fq"
run ks.pe.meth --meth --known-snvs "$tmp/sites.vcf" "$tmp/ref.fa" "$tmp/r1.fq" "$tmp/r2.fq"
for e in se pe; do
    [ "$(grep -c $'\tXV:Z:' "$tmp/ks.$e.meth.sam")" -gt 0 ]
done
echo 'PASS: sites-meth'

# 4. Thread count does not change the output.
for e in se pe se.meth pe.meth; do
    args=()
    [[ $e == *meth ]] && args=(--meth)
    in=("$tmp/se.fq")
    [[ $e == pe* ]] && in=("$tmp/r1.fq" "$tmp/r2.fq")
    run "ks.$e.t4" "${args[@]}" -t 4 --known-snvs "$tmp/sites.vcf" "$tmp/ref.fa" "${in[@]}"
    cmp <(body "$tmp/ks.$e.sam") <(body "$tmp/ks.$e.t4.sam") || {
        echo "FAIL: $e: -t 4 differs from -t 1" >&2
        exit 1
    }
done
echo 'PASS: threads'

# 4b. The ALT-window index is built in memory: a read-only TMPDIR is not touched or needed.
ro=$tmp/ro
mkdir "$ro"
chmod 500 "$ro"
for e in se se.meth; do
    args=()
    [[ $e == *meth ]] && args=(--meth)
    TMPDIR=$ro run "ks.$e.ro" "${args[@]}" --known-snvs "$tmp/sites.vcf" "$tmp/ref.fa" "$tmp/se.fq"
    cmp <(body "$tmp/ks.$e.sam") <(body "$tmp/ks.$e.ro.sam") || {
        echo "FAIL: $e: a read-only TMPDIR changed the output" >&2
        exit 1
    }
done
chmod 700 "$ro"
[ -z "$(ls -A "$ro")" ] || {
    echo "FAIL: --known-snvs wrote into TMPDIR" >&2
    exit 1
}
# chmod 500 does not stop root, so also point TMPDIR at a path no build can mkdtemp in.
for e in se se.meth; do
    args=()
    [[ $e == *meth ]] && args=(--meth)
    TMPDIR="$tmp/does-not-exist" run "ks.$e.nodir" "${args[@]}" --known-snvs "$tmp/sites.vcf" "$tmp/ref.fa" "$tmp/se.fq"
    cmp <(body "$tmp/ks.$e.sam") <(body "$tmp/ks.$e.nodir.sam") || {
        echo "FAIL: $e: a non-existent TMPDIR changed the output" >&2
        exit 1
    }
done
[ ! -e "$tmp/does-not-exist" ] || {
    echo "FAIL: --known-snvs created TMPDIR" >&2
    exit 1
}
echo 'PASS: no-temp-files (read-only and non-existent TMPDIR, empty afterwards)'

# 5. --bam carries XV:Z.
if command -v samtools > /dev/null; then
    "$bin" aln -v 1 --bam --known-snvs "$tmp/sites.vcf" "$tmp/ref.fa" "$tmp/r1.fq" "$tmp/r2.fq" > "$tmp/ks.bam" 2> /dev/null
    cmp <(samtools view "$tmp/ks.bam") <(recs "$tmp/ks.pe.sam") || {
        echo "FAIL: --bam records differ from SAM" >&2
        exit 1
    }
    echo 'PASS: bam'
else
    echo 'SKIP: bam (no samtools)'
fi

# 6. Bad VCFs are refused before any output.
for c in "badref:does not match the reference base" "badctg:is not in the reference" "missing:cannot open"; do
    name=${c%%:*} msg=${c#*:} rc=0
    "$bin" aln --known-snvs "$tmp/$name.vcf" "$tmp/ref.fa" "$tmp/se.fq" > "$tmp/$name.out" 2> "$tmp/$name.err" || rc=$?
    [ "$rc" = 1 ] && grep -qF -- "$msg" "$tmp/$name.err" && [ ! -s "$tmp/$name.out" ] || {
        echo "FAIL: $name: exit $rc, or no '$msg', or output written" >&2
        cat "$tmp/$name.err" >&2
        exit 1
    }
done
echo 'PASS: bad-vcf'
echo 'PASS: aln_known_snvs'
