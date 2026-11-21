#!/usr/bin/env bash
# test/regression/aln_known_snvs.sh
#
# `bwa-mem3 aln --known-snvs` from the command line, on a seeded synthetic reference: a VCF with no
# SNV leaves single-end, paired-end and --meth output unchanged but for @PG; with sites, records without XV:Z
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
# A unique basename (the temporary directory's own name): aln attaches an index staged in shared
# memory by basename, so a fixed name such as ref.fa could load another run's index.
reference="$tmp/${tmp##*/}.fa"

"$bin" aln --help > "$tmp/help"
grep -q -- '--known-snvs FILE' "$tmp/help"
echo 'PASS: help'

# Reference: two random 6,000 bp contigs. Sites: every 150 bp of c1 from 600 to 5,400 (1-based),
# ALT another base. Reads (true position in the name):
#   ref*/alt*  25 bp over a site, REF or ALT there, with one error elsewhere (both place in stock)
#   over*      25 bp, ALT plus two errors: three differences, over the budget of 2 (stock: unmapped)
#   bg*        25 bp from c2, no site
#   many*      25 bp over a site, exactly REF or exactly ALT (16 per site, for the thread check)
# Pairs: R1 is one of the reads above (fragment start), R2 the reverse complement of the
# fragment's last 25 bases (fragments 60-90 bp), never over a site.
cat > "$tmp/gen.awk" << 'AWK'
function rseq(n,   s, i) { s = ""; for (i = 0; i < n; i++) s = s substr("ACGT", ri(4) + 1, 1); return s }
function mut(s, k) { return substr(s, 1, k - 1) other(substr(s, k, 1)) substr(s, k + 1) }
function put(name, s, f) { printf "@%s\n%s\n+\n%s\n", name, s, substr(Q, 1, length(s)) > f }
function pair(name, c, beg, r1) { pair_at(name, c, beg, r1, 60 + ri(31)) }
function pair_at(name, c, beg, r1, flen,   r2) {
    r2 = rc(substr(c, beg + flen - 25, 25))
    put(name, r1, dir "/r1.fq"); put(name, r2, dir "/r2.fq"); put(name, r1, dir "/se.fq")
}
BEGIN {
    seed_rng(20261109)
    Q = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
    c1 = rseq(6000); c2 = rseq(6000)
    printf ">c1\n%s\n>c2\n%s\n", c1, c2 > reference
    print "##fileformat=VCFv4.2\n##contig=<ID=c1>\n##contig=<ID=c2>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" > (dir "/sites.vcf")
    n = 0
    for (p = 600; p <= 5400; p += 150) {
        r = substr(c1, p, 1); a = other(r)
        printf "c1\t%d\t.\t%s\t%s\t.\tPASS\t.\n", p, r, a > (dir "/sites.vcf")
        off = 5 + ri(15); beg = p - off      # the site is base off + 1 of the read
        base = substr(c1, beg, 25)
        e1 = (off + 6) % 25 + 1; e2 = (off + 14) % 25 + 1
        alt = substr(base, 1, off) a substr(base, off + 2)
        printf "%d\tc1:%d:%s>%s\n", n, p, r, a > (dir "/expected-sites.tsv")
        pair("ref" n "|c1|" beg, c1, beg, mut(base, e1))
        pair("alt" n "|c1|" beg, c1, beg, mut(alt, e1))
        pair("over" n "|c1|" beg, c1, beg, mut(mut(alt, e1), e2))
        g = 200 + ri(5500)
        pair("bg" n "|c2|" g, c2, g, substr(c2, g, 25))
        site_pos[n] = p; site_alt[n] = a
        n++
    }
    # More site-covering pairs, drawing no random numbers (the reads above stay as they were), so
    # that the paired twin pass holds more twins than one thread work item (128 reads on x86, 256 on
    # arm64): 16 per site, exact REF or exact ALT, each starting at its own offset.
    for (i = 0; i < n; i++)
        for (k = 0; k < 16; k++) {
            p = site_pos[i]; off = 2 + k; beg = p - off
            base = substr(c1, beg, 25)
            if (k % 2) base = substr(base, 1, off) site_alt[i] substr(base, off + 2)
            # many<letter k><site>: check_reads takes the site from the digits after the letters
            pair_at("many" substr("abcdefghijklmnop", k + 1, 1) i "|c1|" beg, c1, beg, base, 60 + (7 * k) % 31)
        }
    # an empty set: one indel, no SNV
    print "##fileformat=VCFv4.2\n##contig=<ID=c1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" > (dir "/none.vcf")
    printf "c1\t600\t.\t%sA\t%s\t.\tPASS\t.\n", substr(c1, 600, 1), substr(c1, 600, 1) > (dir "/none.vcf")
    printf "c1\t%d\t.\t%s\t%s\t.\tPASS\t.\n", 700, other(substr(c1, 700, 1)), "A" > (dir "/badref.vcf")
    printf "chr9\t5\t.\tA\tC\t.\tPASS\t.\n" > (dir "/badctg.vcf")
}
AWK
awk -v dir="$tmp" -v reference="$reference" -f "$here/aln_rng.awk" -f "$tmp/gen.awk"
for f in badref badctg; do
    {
        printf '##fileformat=VCFv4.2\n##contig=<ID=c1>\n##contig=<ID=chr9>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        cat "$tmp/$f.vcf"
    } > "$tmp/$f.tmp"
    mv "$tmp/$f.tmp" "$tmp/$f.vcf"
done
"$bin" index "$reference" > "$tmp/index.log" 2>&1
"$bin" index --meth "$reference" >> "$tmp/index.log" 2>&1

# body/recs FILE OUT: FILE without its @PG lines / without its header into OUT. Each is its own
# checked command (never a process substitution, whose exit status nothing sees), so a failed
# filter stops the test instead of leaving a comparison of two empty streams.
body() { grep -v '^@PG' "$1" > "$2"; }
# bit(flag, b): whether SAM flag has bit b (a power of two); POSIX awk has no and()
bit='function bit(f, b) { return int(f / b) % 2 } '

recs() { grep -v '^@' "$1" > "$2"; }
run() { # run NAME ARGS...: aln with ARGS (options, index, reads) into $tmp/NAME.sam
    local name=$1
    shift
    "$bin" aln -v 1 -f "$tmp/$name.sam" "$@" > /dev/null 2> "$tmp/$name.err" || {
        echo "FAIL: $name: aln exited non-zero" >&2
        cat "$tmp/$name.err" >&2
        exit 1
    }
}

# Require each generated read of PREFIX exactly once (R1 for pairs). Without VOTE, require
# unmapped stock reads; otherwise require the true locus, MAPQ >= 20 and an XV entry for the
# read's own site (from expected-sites.tsv) whose tail after the site, from its first comma, matches
# VOTE, an anchored regex such as '^,A,1,2,1,[0-9]+$'.
check_reads() {
    local sam=$1 prefix=$2 vote=${3:-}
    awk -F'\t' -v p="$prefix" -v vote="$vote" "$bit"'
        FILENAME == ARGV[1] { site[$1] = $2; next }
        FILENAME == ARGV[2] {
            if (FNR % 4 == 1 && index(substr($0, 2), p) == 1) {
                expected[substr($0, 2)] = 1; n++
            }
            next
        }
        !/^@/ && index($1, p) == 1 && (!bit($2, 1) || bit($2, 64)) {
            if (!($1 in expected) || ++seen[$1] != 1) bad = 1
            split($1, t, "|")
            if (vote == "") valid = bit($2, 4)
            else {
                valid = !bit($2, 4) && $3 == t[2] && $4 == t[3] && $5 >= 20
                id = t[1]; sub(/^[a-z]+/, "", id)
                matched = 0
                for (i = 12; i <= NF; i++) if (index($i, "XV:Z:") == 1) {
                    nv = split(substr($i, 6), entries, ";")
                    for (j = 1; j <= nv; j++) {
                        split(entries[j], fields, ",")
                        if (fields[1] == site[id] && substr(entries[j], length(fields[1]) + 1) ~ vote) matched = 1
                    }
                }
                valid = valid && matched
            }
            if (valid) ok++
            else bad = 1
        }
        END {
            for (name in expected) if (seen[name] != 1) bad = 1
            if (!n || bad || ok != n) {
                printf "%s: %d valid records for %d expected reads (each required exactly once)\n", p, ok, n > "/dev/stderr"
                exit 1
            }
        }' "$tmp/expected-sites.tsv" "$tmp/se.fq" "$sam"
}

# Exact QNAME/mate multiset from the generated FASTQs, including background records; every
# record must have all 11 mandatory fields.
check_all_records() {
    local sam=$1 paired=$2
    awk -F'\t' -v paired="$paired" "$bit"'
        FNR == NR { if (FNR % 4 == 1) { name = substr($0, 2); expected[name, paired ? 64 : 0]++; if (paired) expected[name, 128]++ }; next }
        !/^@/ {
            if (NF < 11) bad = 1 # every record, mapped or not, has the 11 mandatory fields
            mate = paired ? (bit($2, 64) ? 64 : 128) : 0
            if (paired && (!bit($2, 1) || bit($2, 64) == bit($2, 128))) bad = 1
            if (!paired && (bit($2, 1) || bit($2, 64) || bit($2, 128))) bad = 1
            if (!(($1, mate) in expected) || ++seen[$1, mate] != expected[$1, mate]) bad = 1
        }
        END { for (key in expected) if (seen[key] != expected[key]) bad = 1; exit bad }
    ' "$tmp/se.fq" "$sam"
}

# 1. No SNV in the VCF: every expected record, and output identical but for @PG (which records the
# command line).
for mode in plain meth; do
    m=()
    [ $mode = meth ] && m=(--meth)
    run "stock.se.$mode" "${m[@]}" "$reference" "$tmp/se.fq"
    run "none.se.$mode" "${m[@]}" --known-snvs "$tmp/none.vcf" "$reference" "$tmp/se.fq"
    run "stock.pe.$mode" "${m[@]}" "$reference" "$tmp/r1.fq" "$tmp/r2.fq"
    run "none.pe.$mode" "${m[@]}" --known-snvs "$tmp/none.vcf" "$reference" "$tmp/r1.fq" "$tmp/r2.fq"
    for e in se pe; do
        paired=0
        [ "$e" = pe ] && paired=1
        for v in stock none; do
            check_all_records "$tmp/$v.$e.$mode.sam" "$paired" || {
                echo "FAIL: $v $e $mode: missing, duplicate, or unexpected QNAME/mate" >&2
                exit 1
            }
            body "$tmp/$v.$e.$mode.sam" "$tmp/$v.$e.$mode.body"
        done
        cmp "$tmp/stock.$e.$mode.body" "$tmp/none.$e.$mode.body" || {
            echo "FAIL: no-SNV VCF changed the $e $mode output" >&2
            exit 1
        }
    done
done
echo 'PASS: no-snv-identical (single, paired; plain, --meth)'

# 2. With sites (plain): records without XV are unchanged; over-budget ALT reads are placed.
run ks.se --known-snvs "$tmp/sites.vcf" "$reference" "$tmp/se.fq"
run ks.pe --known-snvs "$tmp/sites.vcf" "$reference" "$tmp/r1.fq" "$tmp/r2.fq"
for e in se pe; do
    # a fragment (a read, or a pair) with no XV on any of its records is unchanged
    recs "$tmp/stock.$e.plain.sam" "$tmp/stock.$e.plain.recs"
    recs "$tmp/ks.$e.sam" "$tmp/ks.$e.recs"
    paste -d $'\001' "$tmp/stock.$e.plain.recs" "$tmp/ks.$e.recs" | awk -F $'\001' -v pe="$([ "$e" = pe ] && echo 1 || echo 0)" '
        function done_frag() { if (!xv && !same) { print "differs: " name > "/dev/stderr"; bad = 1 } }
        { split($1, f, "\t") }
        !pe || NR % 2 == 1 { if (NR > 1) done_frag(); name = f[1]; xv = 0; same = 1 }
        { if ($2 ~ /\tXV:Z:/) xv = 1; if ($1 != $2) same = 0 }
        END { if (NR) done_frag(); exit bad }' || {
        echo "FAIL: $e: a read or pair without XV differs from stock" >&2
        exit 1
    }
    # stock leaves every over-budget read unmapped; crediting places each at its locus, with an A vote
    check_reads "$tmp/stock.$e.plain.sam" over || {
        echo "FAIL: $e: stock must report each over-budget ALT read exactly once, unmapped" >&2
        exit 1
    }
    check_reads "$tmp/ks.$e.sam" over '^,A,2,2,(0|1|NA),([0-9]+|NA)$' || {
        echo "FAIL: $e: over-budget ALT reads not placed with an A vote" >&2
        exit 1
    }
    # every extra exact read carries its own allele's vote and a same-locus twin: the generator
    # gives letters at even k (a, c, ...) exact REF and at odd k (b, d, ...) exact ALT
    letters=abcdefghijklmnop
    for ((k = 0; k < 16; k++)); do
        allele=R
        if ((k % 2)); then allele=A; fi
        check_reads "$tmp/ks.$e.sam" "many${letters:k:1}" "^,$allele,0,2,1,[0-9]+$" || {
            echo "FAIL: $e: exact site reads many${letters:k:1}* not placed with an $allele vote and a same-locus twin" >&2
            exit 1
        }
    done
    # REF and ALT reads both carry their vote
    for a in ref:R alt:A; do
        check_reads "$tmp/ks.$e.sam" "${a%%:*}" "^,${a##*:},1,2,1,[0-9]+$" || {
            echo "FAIL: $e: ${a%%:*} reads not placed with an XV ${a##*:} vote and a same-locus twin" >&2
            exit 1
        }
    done
done
echo 'PASS: sites-plain (untouched records unchanged; over-budget ALT reads placed; votes and twins)'

# 3. --meth with sites: complete records and chemistry-specific placement/vote assertions.
run ks.se.meth --meth --known-snvs "$tmp/sites.vcf" "$reference" "$tmp/se.fq"
run ks.pe.meth --meth --known-snvs "$tmp/sites.vcf" "$reference" "$tmp/r1.fq" "$tmp/r2.fq"
for e in se pe; do
    paired=0
    [ "$e" = pe ] && paired=1
    check_all_records "$tmp/ks.$e.meth.sam" "$paired" || {
        echo "FAIL: $e meth: missing, duplicate, or unexpected QNAME/mate" >&2
        exit 1
    }
    awk -F'\t' "$bit"'
        BEGIN {
            vote["ref1|c1|742"] = "c1:750:A>C,R,1,2,1,"
            vote["alt1|c1|742"] = "c1:750:A>C,A,1,2,1,"
            vote["over0|c1|592"] = "c1:600:T>C,A,2,2,1,"
        }
        $1 in vote && (!bit($2, 1) || bit($2, 64)) {
            split($1, locus, "|"); found++
            if (bit($2, 4) || $3 != locus[2] || $4 != locus[3] || $5 < 20 || $6 != "25M") bad = 1
            ok = 0
            for (i = 12; i <= NF; i++) if ($i == "XV:Z:" vote[$1] $5) ok = 1
            if (!ok) bad = 1
        }
        END { exit bad || found != 3 }
    ' "$tmp/ks.$e.meth.sam" || {
        echo "FAIL: $e meth: selected placement or XV differs" >&2
        exit 1
    }
done
echo 'PASS: sites-meth'

# 4. Thread count does not change the output.
for e in se pe se.meth pe.meth; do
    args=()
    [[ $e == *meth ]] && args=(--meth)
    in=("$tmp/se.fq")
    [[ $e == pe* ]] && in=("$tmp/r1.fq" "$tmp/r2.fq")
    run "ks.$e.t4" "${args[@]}" -t 4 --known-snvs "$tmp/sites.vcf" "$reference" "${in[@]}"
    paired=0
    [[ $e == pe* ]] && paired=1
    check_all_records "$tmp/ks.$e.t4.sam" "$paired" || {
        echo "FAIL: $e -t 4: missing, duplicate, or unexpected QNAME/mate" >&2
        exit 1
    }
    body "$tmp/ks.$e.sam" "$tmp/ks.$e.body"
    body "$tmp/ks.$e.t4.sam" "$tmp/ks.$e.t4.body"
    cmp "$tmp/ks.$e.body" "$tmp/ks.$e.t4.body" || {
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
    TMPDIR=$ro run "ks.$e.ro" "${args[@]}" --known-snvs "$tmp/sites.vcf" "$reference" "$tmp/se.fq"
    body "$tmp/ks.$e.sam" "$tmp/ks.$e.body"
    body "$tmp/ks.$e.ro.sam" "$tmp/ks.$e.ro.body"
    cmp "$tmp/ks.$e.body" "$tmp/ks.$e.ro.body" || {
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
    TMPDIR="$tmp/does-not-exist" run "ks.$e.nodir" "${args[@]}" --known-snvs "$tmp/sites.vcf" "$reference" "$tmp/se.fq"
    body "$tmp/ks.$e.sam" "$tmp/ks.$e.body"
    body "$tmp/ks.$e.nodir.sam" "$tmp/ks.$e.nodir.body"
    cmp "$tmp/ks.$e.body" "$tmp/ks.$e.nodir.body" || {
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
    "$bin" aln -v 1 --bam --known-snvs "$tmp/sites.vcf" "$reference" "$tmp/r1.fq" "$tmp/r2.fq" > "$tmp/ks.bam" 2> /dev/null
    samtools view "$tmp/ks.bam" > "$tmp/ks.bam.sam" || {
        echo "FAIL: --bam could not be decoded by samtools" >&2
        exit 1
    }
    cmp "$tmp/ks.bam.sam" "$tmp/ks.pe.recs" || {
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
    "$bin" aln --known-snvs "$tmp/$name.vcf" "$reference" "$tmp/se.fq" > "$tmp/$name.out" 2> "$tmp/$name.err" || rc=$?
    [ "$rc" = 1 ] && grep -qF -- "$msg" "$tmp/$name.err" && [ ! -s "$tmp/$name.out" ] || {
        echo "FAIL: $name: exit $rc, or no '$msg', or output written" >&2
        cat "$tmp/$name.err" >&2
        exit 1
    }
done
echo 'PASS: bad-vcf'
echo 'PASS: aln_known_snvs'
