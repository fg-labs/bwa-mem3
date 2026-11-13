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
#                      default's XA in the same order, plus score and MAPQ;
#                      every other field and tag (name, type and value) as the
#                      default's
#   -u vs -a           each XB hit's score and MAPQ equal the AS:i and MAPQ of
#                      the secondary record -a writes for the same hit (-a writes
#                      no XA/XB, so the values come from the record path)
#   -u --bam           a BAM container, with the same records as -u's SAM text
#   --meth -u          XB:Z, six fields per hit, no XA; the same hits as
#                      --meth's XA in the same order, plus score and MAPQ, and
#                      the same check against --meth -a
#   --meth -u --bam    a BAM container, with the same records as --meth -u's SAM
#                      text
#   --compat=bwa-mem -u  XB:Z as -u writes it (bwa 0.7.18+ has -u), no HN:i
#   --compat=bwa-mem2 -u refused, naming -u (bwa-mem2 has no -u); with
#                      --compat-allow-divergent it warns and writes -u's XB:Z
#
# The fixture is a reference whose second contig is a 2%-diverged copy of a
# window of the first, so reads from that window carry an alternative hit.
#
# Inputs (env vars):
#   BWA_MEM3        the bwa-mem3 binary under test: a bare command name
#                   resolved via PATH, or a path to the binary (either is made
#                   absolute before the `cd "$U_XB_WORK_DIR"` below)
#   U_XB_WORK_DIR   directory for the fixture and outputs; only the files this
#                   script writes are removed from it
# Needs samtools (to read the --bam output) and python3.

set -euo pipefail

trap 'echo "FAIL: u_xb_tag: line $LINENO: $BASH_COMMAND" >&2' ERR

if [[ -z ${BWA_MEM3:-} ]]; then
    echo "FAIL: u_xb_tag: BWA_MEM3 must be set" >&2
    exit 1
fi
# A bare name is looked up on PATH now, before the cd: PATH may hold relative
# entries (./build), which the cd would break. Any relative result is then made
# absolute.
if [[ $BWA_MEM3 != */* ]]; then
    BWA_MEM3="$(command -v "$BWA_MEM3")" || {
        echo "FAIL: u_xb_tag: BWA_MEM3 is not on PATH" >&2
        exit 1
    }
fi
case "$BWA_MEM3" in
    /*) ;;
    *) BWA_MEM3="$(CDPATH='' cd -- "$(dirname -- "$BWA_MEM3")" && pwd)/$(basename -- "$BWA_MEM3")" ;;
esac
if [[ -z ${U_XB_WORK_DIR:-} ]]; then
    echo "FAIL: u_xb_tag: U_XB_WORK_DIR must be set" >&2
    exit 1
fi

command -v samtools > /dev/null 2>&1 || {
    echo "FAIL: samtools not on PATH; the --bam records cannot be read" >&2
    exit 1
}
command -v python3 > /dev/null 2>&1 || {
    echo "FAIL: python3 not on PATH; the fixture cannot be generated or checked" >&2
    exit 1
}

mkdir -p "$U_XB_WORK_DIR"
cd "$U_XB_WORK_DIR"
OUTS=(default.sam u.sam u.bam all.sam meth_default.sam meth_u.sam meth_u.bam
    meth_all.sam compat_mem_u.sam compat_mem2_u_div.sam)
for out in "${OUTS[@]}"; do
    rm -f "$out" "$out.err"
done
rm -f ref.fa ref.fa.amb ref.fa.ann ref.fa.bwt.2bit.64 ref.fa.pac \
    ref.fa.meth.fa ref.fa.meth.amb ref.fa.meth.ann ref.fa.meth.bwt.2bit.64 ref.fa.meth.pac \
    reads.fq reads_ct.fq index.log index-meth.log \
    compat_mem2_u.sam compat_mem2_u.err u.bam.sam meth_u.bam.sam

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

"$BWA_MEM3" index ref.fa > index.log 2>&1 || {
    echo "FAIL: bwa-mem3 index failed" >&2
    tail -5 index.log >&2
    exit 1
}
"$BWA_MEM3" index --meth --force ref.fa > index-meth.log 2>&1 || {
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
run all.sam -a ref.fa reads.fq
run meth_default.sam --meth ref.fa reads_ct.fq
run meth_u.sam --meth -u ref.fa reads_ct.fq
run meth_u.bam --meth -u --bam=0 ref.fa reads_ct.fq
run meth_all.sam --meth -a ref.fa reads_ct.fq
run compat_mem_u.sam --compat=bwa-mem -u ref.fa reads.fq
run compat_mem2_u_div.sam --compat=bwa-mem2 -u --compat-allow-divergent ref.fa reads.fq
grep -q 'proceeding with -u under --compat=bwa-mem2' compat_mem2_u_div.sam.err || {
    echo "FAIL: --compat=bwa-mem2 -u --compat-allow-divergent did not warn about -u" >&2
    exit 1
}
if "$BWA_MEM3" mem -t 1 --compat=bwa-mem2 -u ref.fa reads.fq > compat_mem2_u.sam 2> compat_mem2_u.err; then
    echo "FAIL: --compat=bwa-mem2 -u was accepted; bwa-mem2 has no -u" >&2
    exit 1
fi
grep -q 'mutually exclusive: -u' compat_mem2_u.err || {
    echo "FAIL: --compat=bwa-mem2 -u was refused without naming -u" >&2
    cat compat_mem2_u.err >&2
    exit 1
}
echo "ok: --compat=bwa-mem2 -u: refused, naming -u"
samtools view u.bam > u.bam.sam
samtools view meth_u.bam > meth_u.bam.sam

# The ERR trap does not fire for a command left of ||, so name the failure here:
# the validator prints a FAIL: line per check it fails, but an uncaught
# exception (malformed output it did not anticipate) prints only a traceback.
validator_rc=0
python3 - << 'PY' || validator_rc=$?
import gzip
import sys

fails = []


def records(path):
    out = []
    for line in open(path):
        if line.startswith("@"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 11:
            fails.append("%s: a record with %d fields, expected >= 11: %r" % (path, len(f), line))
            continue
        bad = [t for t in f[11:] if len(t) < 5 or t[2] != ":" or t[4] != ":"]
        if bad:
            fails.append("%s: %s: malformed tag %r" % (path, f[0], bad[0]))
            continue
        # A repeated tag would collapse to its last copy in the dict below, so an
        # extra XA or XB ahead of the right one would go unseen.
        names = [t[:2] for t in f[11:]]
        dup = sorted({n for n in names if names.count(n) > 1})
        if dup:
            fails.append("%s: %s: repeated tag %s" % (path, f[0], ",".join(dup)))
            continue
        # tag -> (type, value), so a type change is caught as well as a value one
        tags = {t[:2]: (t[3], t[5:]) for t in f[11:]}
        for t in f[11:]:
            if t[:2] in ("XA", "XB") and t[2:5] != ":Z:":
                fails.append("%s: %s: %s has type %s, expected :Z:" % (path, f[0], t[:2], t[2:5]))
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
            bad = [h for h in hits(tags[tag][1]) if len(h) != width]
            if bad:
                fails.append("%s: %s %s:Z hit with %d fields, expected %d: %s"
                             % (name, key[0], tag, len(bad[0]), width, ",".join(bad[0])))
    if n < minimum:
        fails.append("%s: only %d records carry %s:Z (expected >= %d)" % (name, n, tag, minimum))
    if len(fails) == before:
        print("ok: %s: %d records carry %s:Z, %d fields per hit, no %s:Z" % (name, n, tag, width, other))


def check_widened(name, base, default, u):
    # -u changes the tag and widens each hit, nothing else: every other field
    # and tag (type and value) is the same, record for record.
    before = len(fails)
    if len(default) != len(u):
        fails.append("%s: %d records, %s %d" % (name, len(u), base, len(default)))
    same = 0
    for (dk, dt), (uk, ut) in zip(default, u):
        d_rest = {k: v for k, v in dt.items() if k != "XA"}
        u_rest = {k: v for k, v in ut.items() if k != "XB"}
        if dk != uk or ("XA" in dt) != ("XB" in ut) or d_rest != u_rest:
            fails.append("%s: record %s differs from %s's beyond the tag" % (name, dk[0], base))
            break
        if "XA" in dt:
            if [h[:4] for h in hits(ut["XB"][1])] != hits(dt["XA"][1]):
                fails.append("%s: %s XB hits are not %s's XA hits" % (name, dk[0], base))
                break
            same += 1
    if len(fails) == before:
        print("ok: %s: %d XB lists are %s's XA hits with score and MAPQ added"
              % (name, same, base))


def check_scores(name, u, all_path):
    # -a writes each alternative hit as a secondary record (and no XA/XB, which
    # is checked here too), so its AS:i and MAPQ give each XB hit's score and
    # MAPQ independently of the list writer. Secondary records hard-clip where
    # the list soft-clips.
    before = len(fails)
    sec = {}
    for key, tags in records(all_path):
        for t in ("XA", "XB"):
            if t in tags:
                fails.append("%s: %s: %s record carries %s:Z, expected none under -a"
                             % (name, all_path, key[0], t))
        flag = int(key[1])
        if flag & 0x100:
            hit = (key[0], key[2], "-+"[not flag & 16] + key[3], key[5].replace("H", "S"))
            for t in ("AS", "NM"):
                if tags.get(t, ("", ""))[0] != "i":
                    fails.append("%s: %s: %s secondary record %s is not an integer tag (type %r)"
                                 % (name, all_path, key[0], t, tags.get(t, ("", ""))[0]))
            sec[hit] = (tags.get("AS", ("", "missing"))[1], key[4],
                        tags.get("NM", ("", "missing"))[1])
    n = 0
    for key, tags in u:
        if "XB" not in tags:
            continue
        for h in hits(tags["XB"][1]):
            hit = (key[0], h[0], h[1], h[2])
            if not all(x.lstrip("-").isdigit() for x in h[3:]):
                fails.append("%s: %s XB hit with a non-integer NM, score or MAPQ: %s"
                             % (name, key[0], ",".join(h)))
            elif hit not in sec:
                fails.append("%s: %s XB hit %s has no secondary record under -a"
                             % (name, key[0], ",".join(h)))
            elif sec[hit] != (h[4], h[5], h[3]):
                fails.append("%s: %s XB hit %s: -a has AS %s, MAPQ %s, NM %s"
                             % ((name, key[0], ",".join(h)) + sec[hit]))
            else:
                n += 1
            if len(fails) > before:
                break
        if len(fails) > before:
            break
    if n < 100 and len(fails) == before:
        fails.append("%s: only %d XB hits checked against -a (expected >= 100)" % (name, n))
    if len(fails) == before:
        print("ok: %s: %d XB hits' score and MAPQ equal -a's AS:i and MAPQ" % (name, n))


default, u = records("default.sam"), records("u.sam")
check_lists("default", default, "XA", "XB", 4, 100)
check_lists("-u", u, "XB", "XA", 6, 100)
check_widened("-u", "the default", default, u)
check_scores("-u", u, "all.sam")

meth_default, meth = records("meth_default.sam"), records("meth_u.sam")
check_lists("--meth", meth_default, "XA", "XB", 4, 100)
check_lists("--meth -u", meth, "XB", "XA", 6, 100)
check_widened("--meth -u", "--meth", meth_default, meth)
check_scores("--meth -u", meth, "meth_all.sam")

# --compat only drops HN:i here (no MQ:i on single-end reads, no sidecar), so
# its records are -u's minus HN, with the same XB lists.
def no_hn(recs):
    return [(k, {t: v for t, v in tags.items() if t != "HN"}) for k, tags in recs]


for name, path in (("--compat=bwa-mem -u", "compat_mem_u.sam"),
                   ("--compat=bwa-mem2 -u --compat-allow-divergent", "compat_mem2_u_div.sam")):
    c = records(path)
    if any("HN" in tags for _, tags in c):
        fails.append("%s: HN:i present" % name)
    elif c != no_hn(u):
        fails.append("%s: records differ from -u's minus HN:i" % name)
    else:
        print("ok: %s: %d records equal -u's minus HN:i (XB:Z kept)" % (name, len(c)))

for name, a, b, path in (("-u --bam", u, records("u.bam.sam"), "u.bam"),
                         ("--meth -u --bam", meth, records("meth_u.bam.sam"), "meth_u.bam")):
    # samtools view also reads SAM, so check the container is BAM: BGZF is
    # gzip, and the decompressed stream starts with the BAM magic.
    try:
        with gzip.open(path, "rb") as fh:
            magic = fh.read(4)
    except OSError:
        magic = b""
    if magic != b"BAM\x01":
        fails.append("%s: %s is not a BAM file (decompressed magic %r)" % (name, path, magic))
    elif a != b:
        fails.append("%s: records differ from the SAM text's" % name)
    else:
        print("ok: %s: a BAM file whose %d records equal the SAM text's" % (name, len(a)))

for f in fails:
    print("FAIL: " + f, file=sys.stderr)
sys.exit(1 if fails else 0)
PY
if [[ $validator_rc -ne 0 ]]; then
    echo "FAIL: u_xb_tag: the output checks failed (see above)" >&2
    exit 1
fi

echo "PASS: u_xb_tag (XB:Z under -u from the SAM, --bam and --meth writers and --compat=bwa-mem; XA:Z otherwise; -u refused under --compat=bwa-mem2)"
