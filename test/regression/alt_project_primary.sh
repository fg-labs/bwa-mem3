#!/usr/bin/env bash
# test/regression/alt_project_primary.sh
#
# End-to-end behavior of `--alt-project primary`: the output names no ALT contig.
#
# Fixture (generated): chrA and chrB, a decoy-like chrD that sorts after the
# ALT contigs, chrA_alt1 (an ALT of a chrA segment with a private SNP, which
# lifts), and chrB_alt2 (600 bp of ALT-only sequence, a 600S clip in its .alt
# line, then 1000 bp of chrB). Read pairs, -I fixed:
#   lift_*: from chrA_alt1, so they lift onto chrA;
#   both_*: both mates inside the ALT-only 600 bp, so no hit outside the ALT;
#   half_*: R1 inside the ALT-only sequence, R2 on chrB's part of chrB_alt2;
#   dec_*:  from chrD;
#   uniq_*: from chrB past chrB_alt2's footprint, so no ALT hit at all.
# Two single-end fixtures (their own generators, as in alt_project.sh): tan, a
# tandem repeat U U U' on chrA whose ALT is U U' U' (--alt-project-min-gap 2
# vetoes the lift of a read across copies 1 and 2), and bal40, an ALT with a
# 40I...40D pair, past the re-alignment band, so its reads' ALT hits do not lift.
#
# Asserts:
#   * the mode is required: a bare --alt-project, an unknown mode, a missing
#     value and a repeated option are refused before the index loads, naming
#     keep-alt and primary;
#   * keep-alt (the non-vacuity control) writes the both_/half_ reads on
#     chrB_alt2 and lists the ALT @SQ;
#   * primary, SAM and --bam: the @SQ block is exactly chrA, chrB, chrD with no
#     AH; no RNAME, RNEXT, SA or XA names an ALT contig (XL:Z may);
#   * a read with no hit outside the ALTs is unmapped (0x4, MAPQ 0, CIGAR *,
#     TLEN 0, no NM/MD/SA/XA) with OA:Z naming its chrB_alt2 alignment, equal
#     (RNAME, POS, strand, CIGAR, MAPQ, NM) to the record keep-alt writes for
#     that end, paired and single-end; a
#     half_ R1 sits at its mate, and that mate has 0x8, no 0x2, RNEXT '=',
#     PNEXT its own POS and TLEN 0; lift_ reads are on chrA as proper pairs
#     with MAPQ > 0 and XL:Z; dec_ reads are on chrD;
#   * every primary record passes an awk mate-consistency check (RNEXT/PNEXT
#     resolve to the mate, 0x8/0x20 agree with the mate's 0x4/0x10, TLEN is
#     antisymmetric, 0x2 only when both are mapped);
#   * the stderr count of read ends written unmapped equals the OA:Z records;
#   * pairs with no ALT hit (uniq_, dec_) have records byte-identical to
#     keep-alt's under -I;
#   * tan under --alt-project-min-gap 2 (the veto shown to fire in the trace):
#     no ALT anywhere, and each read is on chrA at the copy one unit on (AS 150,
#     XS 145 from its own locus), not at its own locus through the vetoed lift;
#   * bal40: no ALT anywhere, and each read keeps keep-alt's primary record on
#     chrA (mapped, no OA:Z), the ALT supplementary keep-alt also writes gone;
#   * direct BAM passes samtools quickcheck and decodes to the SAM records, with
#     chrD's records remapped past the dropped ALT (needs samtools; OPEN without
#     it, a failure with ALT_PROJECT_REQUIRE_SAMTOOLS=1);
#   * ALT @SQ supplied with -H or in the index's .dict sidecar are left out with
#     a note, the other records kept;
#   * output is byte-identical at -t 1, 2 and 4, under --fast at -t 1 and 4, and
#     with --dedup-reads on (memo shown to engage, VERIFY passing) and off;
#   * with --alt-project-exclude-bed over chrA_alt1's footprint, the output still
#     names no ALT and the lift_ reads are on chrA through their own chrA hits;
#   * single-end input follows the same rules.
#
# Inputs:
#   BWA_MEM3 — path to the bwa-mem3 binary under test
#   ALT_PROJECT_REQUIRE_SAMTOOLS (optional) — 1 makes a missing samtools a
#     failure instead of an OPEN BAM check
set -euo pipefail

if [[ -z ${BWA_MEM3:-} ]]; then
    echo "FAIL: BWA_MEM3 must be set" >&2
    exit 1
fi
command -v python3 > /dev/null 2>&1 || {
    echo "FAIL: python3 not on PATH (it generates the fixture)" >&2
    exit 1
}
case "$BWA_MEM3" in
    /*) BIN="$BWA_MEM3" ;;
    *) BIN="$PWD/$BWA_MEM3" ;;
esac
[[ -x "$BIN" ]] || {
    echo "FAIL: BWA_MEM3 ($BIN) is not executable" >&2
    exit 1
}

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
cd "$WORK"
fail() {
    echo "FAIL: $*" >&2
    exit 1
}
nopg() { grep -v '^@PG' "$1"; }

python3 - << 'EOF'
import random
random.seed(23)
B = "ACGT"
def rnd(n): return "".join(random.choice(B) for _ in range(n))
def mut(c): return random.choice([x for x in B if x != c])
RC = str.maketrans("ACGT", "TGCA")
def rc(s): return s.translate(RC)[::-1]
X = rnd(2000)
chrA = rnd(3000) + X + rnd(3000)
chrB = rnd(5000)
chrD = rnd(3000)
alt1 = list(X)
for p in (300, 700, 1000, 1300, 1700): alt1[p] = mut(X[p])
alt1 = "".join(alt1)
Q = rnd(600)
alt2 = Q + chrB[1000:2000]
open("ref.fa", "w").write(f">chrA\n{chrA}\n>chrB\n{chrB}\n>chrA_alt1\n{alt1}\n>chrB_alt2\n{alt2}\n>chrD\n{chrD}\n")
open("ref.fa.alt", "w").write("chrA_alt1\t0\tchrA\t3001\t60\t2000M\t*\t0\t0\t*\t*\n"
                              "chrB_alt2\t0\tchrB\t1001\t60\t600S1000M\t*\t0\t0\t*\t*\n")
open("excl.bed", "w").write("chrA\t3000\t5000\n")
q = "I" * 150
def pair(o1, o2, name, frag):
    o1.write(f"@{name}\n{frag[:150]}\n+\n{q}\n")
    o2.write(f"@{name}\n{rc(frag[-150:])}\n+\n{q}\n")
with open("r1.fq", "w") as r1, open("r2.fq", "w") as r2:
    for i, s in enumerate(range(100, 1550, 50)):
        pair(r1, r2, f"lift_{s}", alt1[s:s + 400])
    for s in range(0, 201, 25):
        pair(r1, r2, f"both_{s}", alt2[s:s + 400])
    for s in range(350, 451, 25):
        pair(r1, r2, f"half_{s}", alt2[s:s + 400])
    for s in range(200, 2400, 300):
        pair(r1, r2, f"dec_{s}", chrD[s:s + 400])
    for s in range(2500, 4500, 400):
        pair(r1, r2, f"uniq_{s}", chrB[s:s + 400])
with open("se.fq", "w") as se:
    for s in (200, 900):
        se.write(f"@slift_{s}\n{alt1[s:s + 150]}\n+\n{q}\n")
    for s in (0, 200, 400):
        se.write(f"@sboth_{s}\n{alt2[s:s + 150]}\n+\n{q}\n")
    se.write(f"@sdec_500\n{chrD[500:650]}\n+\n{q}\n")
# Duplicated pairs (new names) for the --dedup-reads check.
with open("r1d.fq", "w") as o1, open("r2d.fq", "w") as o2:
    for src, dst in (("r1.fq", o1), ("r2.fq", o2)):
        lines = open(src).read().splitlines()
        for i in range(0, len(lines), 4):
            dst.write("\n".join(lines[i:i + 4]) + "\n")
            dst.write("\n".join([lines[i] + "_dup"] + lines[i + 1:i + 4]) + "\n")
# -H text whose @SQ block lists every contig, ALT ones included.
# Tandem fixture: U U U' on chrA, U U' U' on the ALT (U' = U with one SNP). A
# read t_<s> across copies 1 and 2 matches the ALT at its own locus (lifting to
# chrA:3001+s) and chrA one unit on (chrA:3101+s) exactly, and its own chrA
# locus with one mismatch.
rt_ = random.Random(500)
rr = lambda n: "".join(rt_.choice(B) for _ in range(n))
U = rr(100)
U2 = U[:20] + rt_.choice([x for x in B if x != U[20]]) + U[21:]
lf, rgt = rr(500), rr(500)
alt3 = lf + U + U2 + U2 + rgt
open("tan.fa", "w").write(f">chrA\n{rr(3000) + lf + U + U + U2 + rgt + rr(3000)}\n>chrA_alt1\n{alt3}\n")
open("tan.fa.alt", "w").write("chrA_alt1\t0\tchrA\t3001\t60\t1300M\t*\t0\t0\t*\t*\n")
with open("tan.fq", "w") as fq:
    for s in (500, 510, 530, 550):
        fq.write(f"@t_{s}\n{alt3[s:s + 150]}\n+\n{q}\n")
# bal40: an ALT of Z at chrA:3001 with 40 bp inserted after Z[1000) and 40
# deleted 30 bases later, the same fixture alt_project.sh builds.
rb_ = random.Random(140)
rz = lambda n: "".join(rb_.choice(B) for _ in range(n))
Z = rz(2000)
a40 = Z[:1000] + rz(40) + Z[1000:1030] + Z[1070:]
open("bal40.fa", "w").write(f">chrA\n{rz(3000) + Z + rz(3000)}\n>chrA_alt1\n{a40}\n")
open("bal40.fa.alt", "w").write("chrA_alt1\t0\tchrA\t3001\t60\t1000M40I30M40D930M\t*\t0\t0\t*\t*\n")
with open("bal40.fq", "w") as fq:
    for s in (930, 950, 970):
        fq.write(f"@i_{s}_40_40\n{a40[s:s + 150]}\n+\n{q}\n")
open("allsq.hdr", "w").write("@SQ\tSN:chrA\tLN:8000\n@SQ\tSN:chrB\tLN:5000\n@SQ\tSN:chrA_alt1\tLN:2000\n"
                             "@SQ\tSN:chrB_alt2\tLN:1600\n@SQ\tSN:chrD\tLN:3000\n@CO\tuser\n")
open("check.py", "w").write(r'''
import re, sys
path, kind = sys.argv[1], sys.argv[2]
ALTS = {"chrA_alt1", "chrB_alt2"}
OA = re.compile(r"^chrB_alt2,\d+,[+-],(\d+[MIDS])+,\d+,\d+;$")
hdr, recs = [], []
for l in open(path):
    l = l.rstrip("\n")
    (hdr if l.startswith("@") else recs).append(l.split("\t"))
sq = [h for h in hdr if h[0] == "@SQ"]
assert [h[1] for h in sq] == ["SN:chrA", "SN:chrB", "SN:chrD"], f"@SQ block: {[h[1] for h in sq]}"
assert not any(f.startswith("AH:") for h in sq for f in h), "an @SQ carries AH"
def tags(t): return {x[:2]: x[5:] for x in t[11:]}
oa = 0
prim = {}
for t in recs:
    flag = int(t[1])
    assert t[2] not in ALTS and t[6] not in ALTS, f"{t[0]} names an ALT in RNAME/RNEXT"
    tg = tags(t)
    for k in ("SA", "XA"):
        if k in tg:
            names = [e.split(",")[0] for e in tg[k].split(";") if e]
            assert not ALTS & set(names), f"{t[0]} has an ALT in {k}"
    if "OA" in tg:
        oa += 1
        assert flag & 4, f"{t[0]}: OA:Z on a mapped record"
        assert OA.match(tg["OA"]), f"{t[0]}: OA:Z {tg['OA']}"
    if flag & 0x900: continue
    end = 1 if (flag & 0x40 or kind == "se") else 2
    prim[(t[0], end)] = t
mates = lambda n: (prim[(n, 1)], prim[(n, 2)])
xl = 0
names = sorted({n for n, _ in prim})
for n in names:
    if kind == "se":
        t = prim[(n, 1)]; flag = int(t[1]); tg = tags(t)
        if n.startswith("sboth_"):
            assert flag == 4 and t[2] == "*" and t[4] == "0" and t[5] == "*", f"{n}: {t[:9]}"
            assert "OA" in tg and not {"NM", "MD", "SA", "XA"} & tg.keys(), f"{n}: tags {sorted(tg)}"
        elif n.startswith("slift_"):
            assert t[2] == "chrA" and int(t[4]) > 0, f"{n}: {t[:6]}"
            if tg.get("XL") == "chrA_alt1": xl += 1
        else:
            assert t[2] == "chrD", f"{n}: {t[:6]}"
        continue
    a, b = mates(n)
    fa, fb = int(a[1]), int(b[1])
    ta, tb = tags(a), tags(b)
    if n.startswith("both_"):
        for t, f, tg in ((a, fa, ta), (b, fb, tb)):
            assert f & 0xC == 0xC and not f & 2, f"{n}: flag {f}"
            assert t[2] == "*" and t[4] == "0" and t[5] == "*" and t[8] == "0", f"{n}: {t[:9]}"
            assert "OA" in tg and not {"NM", "MD", "SA", "XA"} & tg.keys(), f"{n}: tags {sorted(tg)}"
    elif n.startswith("half_"):
        assert fa & 4 and not fa & 2 and a[4] == "0" and a[5] == "*" and a[8] == "0", f"{n} R1: {a[:9]}"
        assert "OA" in ta and not {"NM", "MD", "SA", "XA"} & ta.keys(), f"{n} R1 tags {sorted(ta)}"
        assert a[2] == b[2] == "chrB" and a[3] == b[3], f"{n}: R1 not placed at its mate"
        assert not fb & 4 and fb & 8 and not fb & 2, f"{n} R2: flag {fb}"
        assert b[6] == "=" and b[7] == b[3] and b[8] == "0", f"{n} R2 mate fields {b[6:9]}"
    elif n.startswith("lift_"):
        for t, f, tg in ((a, fa, ta), (b, fb, tb)):
            assert t[2] == "chrA" and f & 2 and int(t[4]) > 0, f"{n}: {t[:6]}"
            # Only an end covering one of the ALT's SNPs beats chrA and lifts.
            if tg.get("XL") == "chrA_alt1": xl += 1
    elif n.startswith("dec_"):
        assert a[2] == b[2] == "chrD" and fa & 2, f"{n}: {a[:6]}"
    elif n.startswith("uniq_"):
        assert a[2] == b[2] == "chrB" and fa & 2 and fb & 2, f"{n}: {a[:6]}"
assert xl > 0, "no lift_ read was placed through its lifted hit (no XL:Z)"
assert oa > 0, "no record carries OA:Z"
print(oa)
''')
EOF

"$BIN" index ref.fa > idx.log 2>&1 || fail "index ref.fa failed"

# The awk mate-consistency checker: every pair's two primary records must agree.
matecheck() {
    awk -F'\t' '
        /^@/ { next }
        { f = $2 + 0 }
        int(f / 256) % 2 == 1 || int(f / 2048) % 2 == 1 { next }
        {
            e = (int(f / 64) % 2 == 1) ? 1 : 2
            n[$1] = 1; fl[$1, e] = f; rn[$1, e] = $3; po[$1, e] = $4
            nx[$1, e] = $7; pn[$1, e] = $8; tl[$1, e] = $9
        }
        function bit(x, b) { return int(x / b) % 2 }
        END {
            bad = 0
            for (q in n) for (e = 1; e <= 2; ++e) {
                o = 3 - e
                if (!((q, o) in fl)) { print "no mate for " q; bad = 1; continue }
                r = nx[q, e]; if (r == "=") r = rn[q, e]
                if (r != rn[q, o] || pn[q, e] != po[q, o]) { print q " " e ": RNEXT/PNEXT do not point at the mate"; bad = 1 }
                if (bit(fl[q, e], 8) != bit(fl[q, o], 4)) { print q " " e ": 0x8 disagrees with the mate 0x4"; bad = 1 }
                if (bit(fl[q, e], 32) != bit(fl[q, o], 16)) { print q " " e ": 0x20 disagrees with the mate 0x10"; bad = 1 }
                if (tl[q, e] + tl[q, o] != 0) { print q " " e ": TLEN not antisymmetric"; bad = 1 }
                if (bit(fl[q, e], 2) && (bit(fl[q, e], 4) || bit(fl[q, o], 4))) { print q " " e ": 0x2 with an unmapped end"; bad = 1 }
            }
            exit bad
        }' "$1"
}

# --- The mode is required and checked before the index loads. ---
"$BIN" mem --alt-project ref.fa r1.fq r2.fq > /dev/null 2> m.bare.log && fail "a bare --alt-project was accepted"
grep -q "takes a mode, keep-alt or primary (got 'ref.fa'" m.bare.log || fail "a bare --alt-project gave the wrong error"
"$BIN" mem --alt-project bogus missing.fa r1.fq r2.fq > /dev/null 2> m.bad.log && fail "an unknown mode was accepted"
grep -q "takes a mode, keep-alt or primary (got 'bogus'" m.bad.log || fail "an unknown mode gave the wrong error"
grep -q "missing.fa" m.bad.log && fail "an unknown mode was reported only after looking at the index"
"$BIN" mem ref.fa r1.fq r2.fq --alt-project > /dev/null 2> m.none.log && fail "--alt-project without a value was accepted"
grep -q "alt-project" m.none.log || fail "--alt-project without a value gave an unrelated error"
"$BIN" mem --alt-project primary --alt-project keep-alt ref.fa r1.fq r2.fq > /dev/null 2> m.twice.log \
    && fail "a repeated --alt-project was accepted"
grep -q "may be supplied only once" m.twice.log || fail "a repeated --alt-project gave the wrong error"
"$BIN" mem 2> help.log || true
grep -q -e "--alt-project MODE" help.log || fail "help does not show --alt-project MODE"

run() { # out extra...
    local out=$1
    shift
    "$BIN" mem -t 2 -I 400,10 "$@" ref.fa r1.fq r2.fq > "$out" 2> "$out.log" || fail "mem $* failed (see $out.log)"
}

# --- keep-alt: the non-vacuity control. ---
run keep.sam --alt-project keep-alt
grep -q $'^@SQ\tSN:chrB_alt2\t' keep.sam || fail "keep-alt has no ALT @SQ (the fixture is vacuous)"
awk -F'\t' '!/^@/ && ($1 ~ /^(both|half)_/) && $3 == "chrB_alt2"' keep.sam | grep -q . \
    || fail "keep-alt wrote no both_/half_ record on chrB_alt2 (the fixture is vacuous)"
grep -q "alt-project primary" keep.sam.log && fail "keep-alt reported the primary count"

# --- primary, paired, SAM. ---
run pri.sam --alt-project primary
n_oa=$(python3 check.py pri.sam pe) || fail "primary output broke a rule (above)"
matecheck pri.sam > pri.mate.txt || fail "primary output failed the mate-consistency check: $(head -3 pri.mate.txt)"
[[ $n_oa -gt 0 ]] || fail "no record carries OA:Z"
grep -q "alt-project primary: $n_oa read ends had only unlifted ALT hits" pri.sam.log \
    || fail "the stderr count is not the $n_oa OA:Z records (see pri.sam.log)"

# --- single-end. ---
"$BIN" mem -t 2 --alt-project primary ref.fa se.fq > se.sam 2> se.log || fail "single-end primary failed"
python3 check.py se.sam se > /dev/null || fail "single-end primary output broke a rule (above)"
"$BIN" mem -t 2 --alt-project keep-alt ref.fa se.fq > kse.sam 2> kse.log || fail "single-end keep-alt failed"

# --- OA:Z is the alignment keep-alt writes for that end. A demoted end's only
# usable hits were on chrB_alt2, so keep-alt writes its primary record there;
# OA:Z must be that record's RNAME,POS,strand,CIGAR,MAPQ,NM, and every such
# keep-alt end must be demoted with OA:Z. ---
cat > oa_oracle.py << 'EOF'
import sys
def ends(path, keep):
    out = {}
    for l in open(path):
        if l[0] == "@": continue
        t = l.rstrip("\n").split("\t")
        f = int(t[1])
        tg = {x[:2]: x[5:] for x in t[11:]}
        key = (t[0], 1 if f & 0x40 else 2 if f & 0x80 else 0)
        if keep:
            if t[2] == "chrB_alt2" and not f & 0x900:
                out[key] = f"{t[2]},{t[3]},{'-' if f & 16 else '+'},{t[5]},{t[4]},{tg['NM']};"
        elif "OA" in tg:
            out[key] = tg["OA"]
    return out
keep, pri = ends(sys.argv[1], True), ends(sys.argv[2], False)
assert keep, "keep-alt wrote no chrB_alt2 record (vacuous)"
assert set(keep) == set(pri), f"ends demoted {sorted(set(pri) - set(keep))[:3]}, keep-alt on chrB_alt2 {sorted(set(keep) - set(pri))[:3]}"
bad = [(k, pri[k], keep[k]) for k in sorted(keep) if pri[k] != keep[k]]
assert not bad, f"{len(bad)} OA:Z differ from keep-alt, e.g. {bad[0]}"
print(len(keep))
EOF
n_pe=$(python3 oa_oracle.py keep.sam pri.sam) || fail "paired OA:Z is not keep-alt's alignment of the end (above)"
n_se=$(python3 oa_oracle.py kse.sam se.sam) || fail "single-end OA:Z is not keep-alt's alignment of the read (above)"
[[ $n_pe == "$n_oa" ]] || fail "$n_pe paired ends checked against keep-alt, want all $n_oa OA:Z records"
[[ $n_se -gt 0 ]] || fail "no single-end OA:Z checked"

# --- Pairs with no ALT hit are untouched: under -I (no insert-size estimate)
# their records are byte-identical between keep-alt and primary. ---
nalt() { awk -F'\t' '!/^@/ && $1 ~ /^(uniq|dec)_/' "$1"; }
nalt keep.sam > keep.nalt
nalt pri.sam > pri.nalt
[[ $(grep -c '^uniq_' keep.nalt) -ge 10 ]] || fail "keep-alt wrote too few uniq_ records (the check is vacuous)"
cmp -s keep.nalt pri.nalt || fail "primary changed records of pairs with no ALT hit: $(diff keep.nalt pri.nalt | head -2)"

# --- ALT @SQ from -H are left out, with a note. ---
run hdr.sam --alt-project primary -H allsq.hdr
python3 check.py hdr.sam pe > /dev/null || fail "-H run broke a rule (above)"
grep -q "left out 2 ALT @SQ records" hdr.sam.log || fail "no note about the ALT @SQ dropped from -H"
grep -q $'^@CO\tuser' hdr.sam || fail "-H's other records were not kept"

# --- ALT @SQ from the index's .dict sidecar are left out too, other fields kept. ---
for f in ref.fa ref.fa.*; do ln -s "$f" "sc${f#ref}"; done
printf '@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chrA\tLN:8000\tM5:a\n@SQ\tSN:chrB\tLN:5000\tM5:b\n@SQ\tSN:chrA_alt1\tLN:2000\tAH:*\n@SQ\tSN:chrB_alt2\tLN:1600\tAH:*\n@SQ\tSN:chrD\tLN:3000\tM5:d\n' > sc.dict
"$BIN" mem -t 2 -I 400,10 --alt-project primary sc.fa r1.fq r2.fq > sc.sam 2> sc.sam.log || fail "the sidecar run failed (see sc.sam.log)"
python3 check.py sc.sam pe > /dev/null || fail "the sidecar run broke a rule (above)"
grep -q $'^@SQ\tSN:chrD\tLN:3000\tM5:d$' sc.sam || fail "the sidecar's own @SQ records were not used"
grep -q "left out 2 ALT @SQ records" sc.sam.log || fail "no note about the ALT @SQ dropped from the sidecar"

# --- Determinism. ---
for t in 1 4; do
    "$BIN" mem -t "$t" -I 400,10 --alt-project primary ref.fa r1.fq r2.fq > "det.t$t.sam" 2> "det.t$t.log" || fail "mem -t $t failed"
    "$BIN" mem -t "$t" -I 400,10 --fast --alt-project primary ref.fa r1.fq r2.fq > "detf.t$t.sam" 2> "detf.t$t.log" \
        || fail "mem -t $t --fast failed"
    BWAMEM3_DEDUP_READS_STATS=1 "$BIN" mem -t "$t" -I 400,10 --dedup-reads on --alt-project primary ref.fa r1d.fq r2d.fq \
        > "detd.on.t$t.sam" 2> "detd.on.t$t.log" || fail "mem -t $t --dedup-reads on failed"
    "$BIN" mem -t "$t" -I 400,10 --dedup-reads off --alt-project primary ref.fa r1d.fq r2d.fq \
        > "detd.off.t$t.sam" 2> "detd.off.t$t.log" || fail "mem -t $t --dedup-reads off failed"
    cmp -s <(nopg "detd.on.t$t.sam") <(nopg "detd.off.t$t.sam") || fail "output differs with --dedup-reads on and off (-t $t)"
    grep -q -E '^\[dedup-reads-stats\] pairs=[0-9]+ dup=[1-9][0-9]* .*final=ON' "detd.on.t$t.log" \
        || fail "--dedup-reads on memoized no pair (-t $t, see detd.on.t$t.log)"
done
cmp -s <(nopg det.t1.sam) <(nopg det.t4.sam) || fail "output differs between -t 1 and -t 4"
cmp -s <(nopg det.t1.sam) <(nopg pri.sam) || fail "output differs between -t 1 and -t 2"
cmp -s <(nopg detf.t1.sam) <(nopg detf.t4.sam) || fail "--fast output differs between -t 1 and -t 4"
python3 check.py detf.t1.sam fast > /dev/null || fail "--fast output broke a rule (above)"
cmp -s <(nopg detd.on.t1.sam) <(nopg detd.on.t4.sam) || fail "--dedup-reads on output differs between -t 1 and -t 4"
BWAMEM3_DEDUP_READS_VERIFY=1 BWAMEM3_DEDUP_READS_STATS=1 "$BIN" mem -t 4 -I 400,10 --dedup-reads on --alt-project primary \
    ref.fa r1d.fq r2d.fq > detd.verify.sam 2> detd.verify.log || fail "--dedup-reads VERIFY failed (see detd.verify.log)"
cmp -s <(nopg detd.verify.sam) <(nopg detd.off.t4.sam) || fail "--dedup-reads VERIFY output differs from the off run's"

# --- With the exclusion BED over chrA_alt1's footprint. ---
run bed.sam --alt-project primary --alt-project-exclude-bed excl.bed
python3 - << 'EOF' || fail "the BED run broke a rule"
ALTS = {"chrA_alt1", "chrB_alt2"}
for l in open("bed.sam"):
    if l[0] == "@":
        assert not any(f"SN:{a}\t" in l for a in ALTS), l
        continue
    t = l.split("\t")
    assert t[2] not in ALTS and t[6] not in ALTS, t[0]
    if t[0].startswith("lift_") and not int(t[1]) & 0x900:
        assert t[2] == "chrA", t[:6]
        assert not any(x.startswith("XL:") for x in t[11:]), t[0]
EOF
grep -q "excluded=1 retained=1" bed.sam.log || fail "the BED run did not exclude chrA_alt1 (see bed.sam.log)"

# --- Single-end reads whose ALT hit is not lifted, beyond the .alt clip. ---
for r in tan bal40; do "$BIN" index "$r.fa" > "idx.$r.log" 2>&1 || fail "index $r.fa failed"; done
noalt() { # sam: no @SQ, RNAME, RNEXT, SA or XA names chrA_alt1
    python3 - "$1" << 'EOF'
import sys
for l in open(sys.argv[1]):
    t = l.rstrip("\n").split("\t")
    if t[0] == "@SQ": assert t[1] != "SN:chrA_alt1", "an @SQ names the ALT"
    if l[0] == "@": continue
    tg = {x[:2]: x[5:] for x in t[11:]}
    assert "chrA_alt1" not in (t[2], t[6]), f"{t[0]} names the ALT in RNAME/RNEXT"
    for k in ("SA", "XA"):
        assert "chrA_alt1" not in [e.split(",")[0] for e in tg.get(k, "").split(";")], f"{t[0]} names the ALT in {k}"
EOF
}

# The min-gap veto: the lift of t_<s> ties the chrA copy one unit on, so
# --alt-project-min-gap 2 refuses it and the ALT hit stays on the ALT, where
# primary drops it. The read is then at that copy (chrA:3101+s, AS 150), with
# its own locus (AS 145) as XS -- not at chrA:3001+s through the lift, which
# would tie (XS 150, MAPQ 0).
BWAMEM3_ALT_PROJECT_TRACE=t_500 "$BIN" mem -t 2 --alt-project primary --alt-project-min-gap 2 tan.fa tan.fq \
    > tan.sam 2> tan.log || fail "primary with the min-gap veto failed (see tan.log)"
grep -q $'^\\[alt-project-trace\\]\tt_500\t.*\tveto\t' tan.log || fail "the min-gap veto did not fire on t_500 (the check is vacuous)"
noalt tan.sam || fail "primary with the min-gap veto named the ALT (above)"
python3 - tan.sam << 'EOF' || fail "primary with the min-gap veto misplaced a read"
import sys
seen = set()
for l in open(sys.argv[1]):
    if l[0] == "@": continue
    t = l.rstrip("\n").split("\t")
    if int(t[1]) & 0x900: continue
    s = int(t[0].split("_")[1])
    tg = {x[:2]: x[5:] for x in t[11:]}
    got = (t[1], t[2], t[3], t[5], tg.get("AS"), tg.get("XS"), "OA" in tg)
    want = ("0", "chrA", str(3101 + s), "150M", "150", "145", False)
    assert got == want, f"{t[0]}: {got}, want {want}"
    seen.add(t[0])
assert len(seen) == 4, f"{len(seen)} reads, want 4"
EOF

# Band refusal: bal40's 40I...40D is past the re-alignment band, so its reads'
# ALT hits do not lift. Each read still has chrA hits above -T: keep-alt writes
# its primary record there and an ALT supplementary; primary keeps that primary
# record (mapped, no OA:Z) and writes nothing on the ALT.
"$BIN" mem -t 2 --alt-project keep-alt bal40.fa bal40.fq > bal40.keep.sam 2> bal40.keep.log || fail "keep-alt on bal40 failed"
"$BIN" mem -t 2 --alt-project primary bal40.fa bal40.fq > bal40.sam 2> bal40.log || fail "primary on bal40 failed (see bal40.log)"
noalt bal40.sam || fail "primary on bal40 named the ALT (above)"
python3 - bal40.keep.sam bal40.sam << 'EOF' || fail "primary on bal40 misplaced or unmapped a read"
import sys
def recs(path):
    out = {}
    for l in open(path):
        if l[0] != "@":
            t = l.rstrip("\n").split("\t")
            out.setdefault(t[0], []).append(t)
    return out
keep, pri = recs(sys.argv[1]), recs(sys.argv[2])
assert len(keep) == 3 and set(keep) == set(pri), f"reads {sorted(pri)}, want keep-alt's {sorted(keep)}"
for n, rs in keep.items():
    assert any(t[2] == "chrA_alt1" for t in rs), f"{n}: keep-alt has no ALT record (vacuous)"
    kp = [t for t in rs if not int(t[1]) & 0x900]
    pp = [t for t in pri[n] if not int(t[1]) & 0x900]
    assert len(kp) == 1 and len(pp) == 1, f"{n}: primary records {len(kp)}/{len(pp)}"
    k, p = kp[0], pp[0]
    assert not int(p[1]) & 4 and p[2] == "chrA" and p[:11] == k[:11], f"{n}: {p[1:9]}, want keep-alt's {k[1:9]}"
    assert not any(x.startswith("OA:Z:") for x in p[11:]), f"{n}: OA:Z on a mapped record"
EOF

# --- Direct BAM. ---
if command -v samtools > /dev/null 2>&1; then
    run pri.bam --alt-project primary --bam=0
    samtools quickcheck pri.bam || fail "the primary BAM failed quickcheck"
    samtools view -h --no-PG pri.bam > pri.bam.sam
    cmp -s <(grep '^@SQ' pri.sam) <(grep '^@SQ' pri.bam.sam) || fail "the BAM @SQ block differs from the SAM's"
    cmp -s <(grep -v '^@' pri.sam) <(grep -v '^@' pri.bam.sam) || fail "the BAM records differ from the SAM records"
    echo "PASS: BAM checked"
elif [[ ${ALT_PROJECT_REQUIRE_SAMTOOLS:-0} == 1 ]]; then
    fail "samtools is required (ALT_PROJECT_REQUIRE_SAMTOOLS=1) for the BAM check"
else
    echo "OPEN: samtools not on PATH; the primary BAM check did not run"
fi

echo "PASS: alt_project_primary ($n_oa paired and $n_se single-end read ends unmapped with OA:Z equal to keep-alt's alignment; header lists only primary contigs; no ALT in RNAME, RNEXT, SA or XA; mate fields consistent; pairs with no ALT hit unchanged; min-gap veto and band refusal placed on chrA; byte-identical across -t, --fast and --dedup-reads; BED and -H handled)"
