#!/usr/bin/env bash
# test/regression/u_xb_tag.sh
#
# Regression: `mem -u` adds each alternative hit's score and MAPQ to the hit
# list (chr,pos,CIGAR,NM,score,mapq) and names the list XB:Z, as bwa does and as
# the -u help text says. It used to keep the XA:Z name, so a reader parsing XA
# with the four-field layout got six fields under -u. The tag is written by
# three writers (the SAM text path, --bam, and --meth for both containers), so
# each is checked:
#
#   default            XA:Z, four fields per hit, no XB
#   -u                 XB:Z, six fields per hit, no XA; the same hits as the
#                      default's XA in the same order, plus score and MAPQ
#   -u --bam           the same records as -u's SAM text
#   --meth -u          XB:Z, six fields per hit, no XA
#   --meth -u --bam    the same records as --meth -u's SAM text
#
# The fixture is a reference whose second contig is a 2%-diverged copy of a
# window of the first, so reads from that window carry an alternative hit.
#
# Inputs (env vars):
#   BWA_MEM3        path to the bwa-mem3 binary under test
#   U_XB_WORK_DIR   directory for the fixture and outputs
# Needs samtools (to read the --bam output) and python3.

set -euo pipefail

trap 'echo "FAIL: u_xb_tag: line $LINENO: $BASH_COMMAND" >&2' ERR

: "${BWA_MEM3:?BWA_MEM3 must be set}"
: "${U_XB_WORK_DIR:?U_XB_WORK_DIR must be set}"

command -v samtools > /dev/null 2>&1 || {
    echo "FAIL: samtools not on PATH; the --bam records cannot be read" >&2
    exit 1
}

mkdir -p "$U_XB_WORK_DIR"
cd "$U_XB_WORK_DIR"
rm -f ref.fa ref.fa.* reads.fq reads_ct.fq ./*.sam ./*.bam ./*.log ./*.err

python3 - << 'PY'
import random

random.seed(2026)
B = "ACGT"
primary = "".join(random.choice(B) for _ in range(12000))
window = list(primary[3000:9000])
for i, b in enumerate(window):
    if random.random() < 0.02:
        window[i] = random.choice([x for x in B if x != b])
copy = "".join(window)


def wrap(s, w=60):
    return "\n".join(s[i:i + w] for i in range(0, len(s), w))


with open("ref.fa", "w") as fh:
    fh.write(">chr1\n%s\n>chr2\n%s\n" % (wrap(primary), wrap(copy)))
with open("reads.fq", "w") as fq, open("reads_ct.fq", "w") as ct:
    for i in range(400):
        st = random.randrange(0, len(copy) - 150)
        r = copy[st:st + 150]
        fq.write("@r%d\n%s\n+\n%s\n" % (i, r, "I" * 150))
        ct.write("@r%d\n%s\n+\n%s\n" % (i, r.replace("C", "T"), "I" * 150))
PY

"$BWA_MEM3" index ref.fa > index.log 2>&1
"$BWA_MEM3" index --meth ref.fa > index-meth.log 2>&1 || {
    echo "FAIL: bwa-mem3 index --meth failed" >&2
    tail -5 index-meth.log >&2
    exit 1
}

run() { # run <out> <mem args...>
    local out="$1"
    shift
    "$BWA_MEM3" mem -t 1 "$@" > "$out" 2> "$out.err" || {
        echo "FAIL: mem $* exited non-zero" >&2
        tail -5 "$out.err" >&2
        exit 1
    }
}
run default.sam ref.fa reads.fq
run u.sam -u ref.fa reads.fq
run u.bam -u --bam=0 ref.fa reads.fq
run meth_u.sam --meth -u ref.fa reads_ct.fq
run meth_u.bam --meth -u --bam=0 ref.fa reads_ct.fq
samtools view u.bam > u.bam.sam
samtools view meth_u.bam > meth_u.bam.sam

python3 - << 'PY' || exit 1
import sys

fails = []


def records(path):
    out = []
    for line in open(path):
        if line.startswith("@"):
            continue
        f = line.rstrip("\n").split("\t")
        tags = {t[:2]: t[5:] for t in f[11:]}
        out.append((tuple(f[:11]), tags))
    return out


def hits(v):
    return [h.split(",") for h in v.rstrip(";").split(";")]


def check_lists(name, recs, tag, other, width, minimum):
    before = len(fails)
    n = 0
    for key, tags in recs:
        if other in tags:
            fails.append("%s: %s carries %s:Z" % (name, key[0], other))
        if tag in tags:
            n += 1
            bad = [h for h in hits(tags[tag]) if len(h) != width]
            if bad:
                fails.append("%s: %s %s:Z hit with %d fields, expected %d: %s"
                             % (name, key[0], tag, len(bad[0]), width, ",".join(bad[0])))
    if n < minimum:
        fails.append("%s: only %d records carry %s:Z (expected >= %d)" % (name, n, tag, minimum))
    if len(fails) == before:
        print("ok: %s: %d records carry %s:Z, %d fields per hit, no %s:Z" % (name, n, tag, width, other))


default, u = records("default.sam"), records("u.sam")
check_lists("default", default, "XA", "XB", 4, 100)
check_lists("-u", u, "XB", "XA", 6, 100)

# -u changes the tag and widens each hit, nothing else.
if len(default) != len(u):
    fails.append("-u: %d records, default %d" % (len(u), len(default)))
same, before = 0, len(fails)
for (dk, dt), (uk, ut) in zip(default, u):
    if dk != uk or ("XA" in dt) != ("XB" in ut):
        fails.append("-u: record %s differs from the default's beyond the tag" % dk[0])
        break
    if "XA" in dt:
        if [h[:4] for h in hits(ut["XB"])] != hits(dt["XA"]):
            fails.append("-u: %s XB hits are not the default's XA hits" % dk[0])
            break
        same += 1
if len(fails) == before:
    print("ok: -u: %d XB lists are the default's XA hits with score and MAPQ added" % same)

meth = records("meth_u.sam")
check_lists("--meth -u", meth, "XB", "XA", 6, 100)

for name, a, b in (("-u --bam", u, records("u.bam.sam")),
                   ("--meth -u --bam", meth, records("meth_u.bam.sam"))):
    if a != b:
        fails.append("%s: records differ from the SAM text's" % name)
    else:
        print("ok: %s: %d records equal the SAM text's" % (name, len(a)))

for f in fails:
    print("FAIL: " + f, file=sys.stderr)
sys.exit(1 if fails else 0)
PY

echo "PASS: u_xb_tag (XB:Z under -u from the SAM, --bam and --meth writers; XA:Z otherwise)"
